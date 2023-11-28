from typing import Callable, Union
from modules.FBA_S6803 import growth, get_flux, get_fluxes, get_fva
from modules.import_elife import DataImporter
import numpy as np
import numpy.typing
import matplotlib.pyplot as plt
import cobra
from modules.add_glycogen import update_glycogen, update_prot_glyc
from modules.fit_glycogen_funcs import fit as fit_glyc

plt.rc('svg', hashsalt="abc")     # hashsalt for reproducable svg

zavrel = DataImporter()

blue_light = 27.5

fitted_glyc = fit_glyc(
    zavrel.get("light") + blue_light,
    zavrel.get("glycogen")/zavrel.get("DW")
)


def proteincontent(input: float) -> float:
    '''
    Interpolate protein content based on light intensity.

    This function calculates the protein content in relation to the light
    intensity. It performs linear interpolation between known data points
    of light intensity and protein content.

    Parameters:
    input (float): The light intensity for which protein content is to be
        interpolated.

    Returns:
    float: The interpolated protein content corresponding to the given light
        intensity.

    Notes:
    - The function uses a linear interpolation technique to estimate protein
      content based on the provided light intensity value.
    - The data points for light intensity and protein content are obtained
      from external sources (zavrel.get("light") and related data).
    - If the provided light intensity is outside the range of available data,
      the function extrapolates the protein content using linear interpolation
      based on the nearest data points.
    - If the input is beyond the range of data points on both ends, it will be
      extrapolated beyond the range.
    '''

    lights_in = zavrel.get("light") + blue_light
    protein_cont = zavrel.get("protein")/zavrel.get("DW")*100

    def linear_function(point1, point2, x):
        # Unpack the points into their x and y components
        x1, y1 = point1
        x2, y2 = point2

        # Calculate the slope and y-intercept of the line
        slope = (y2 - y1) / (x2 - x1)
        y_intercept = y1 - slope * x1

        # Calculate the y value for the given x value
        y = slope * x + y_intercept

        return y

    if input > lights_in[-2]:
        return linear_function(
            (lights_in[-2], protein_cont[-2]),
            (lights_in[-1], protein_cont[-1]),
            input
        )
    elif input < lights_in[0]:
        return linear_function(
            (0, 51),
            (lights_in[0], protein_cont[0]),
            input
        )
    for i in range(1, len(lights_in)):
        if input < lights_in[i]:
            return linear_function(
                (lights_in[i-1], protein_cont[i-1]),
                (lights_in[i], protein_cont[i]),
                input
            )


def get_kl(model: cobra.Model, pmax: float) -> float:
    '''calculates KL for a given mumax without maintenance'''
    with model as m1:
        # zero all N influxes:
        for rxn in model.boundary:
            if 'N' in rxn.check_mass_balance():
                rxn.lower_bound = 0.

        model.reactions.get_by_id(
            "EX_C00244_ext_b"
        ).lower_bound = -1000  # nitrate
        model.reactions.get_by_id(
            "EX_C00011_ext_b"
        ).lower_bound = -1000
        m1.reactions.R_ATP_maint.knock_out()  # no maintenance
        for con in [
            "respiration"
            "mehlerlike"
            "PSI_mehler"
            "PSII_mehler"
            "rubisco"
        ]:
            if con in m1.constraints._dict.keys():
                m1.remove_cons_vars((con))
        m1.reactions.get_by_id("EX_E1_ext_b").lower_bound = -1*pmax
        mumax = m1.slim_optimize()  # maximal growth rate mu
        # print(mumax)
        # get the yield (BM per photon):
        m1.reactions.get_by_id("EX_E1_ext_b").lower_bound = -1
        BMyield = m1.slim_optimize()
    kl = mumax/BMyield
    return kl


def photodamage_helper(
    model: cobra.Model, rxn_in: str, rxn_out: str,
    input1: float, kl: float, damage: float,
    alpha: float = 1, alpha2: float = 0, alpha3: float = 0,
    add_glyc: tuple = (False, 0), maint: float = 1.3, fva: bool = False
) -> Union[float, tuple]:
    '''Implements a photodamge with maintenance damage percent of photon flux
    returns either the objective value or a flux specified in rxn_out'''
    # input = alpha * input1
    input = -plus_blue(-input1, alpha, alpha2, alpha3)
    vstar = -1*input*kl/(-1*input + kl)  # reduced photon usage
    with model as m:
        atp = m.reactions.get_by_id("R_ATP_maint")
        photodamage = model.problem.Constraint(
            damage/100 * input + atp.flux_expression - maint,
            lb=0,
            ub=0
        )
        m.add_cons_vars(photodamage)
        try:
            if add_glyc[0]:
                if add_glyc[1] == 1:
                    for i in range(10):
                        gr = growth(m, rxn_in, vstar)
                        update_glycogen(m, 3.4+gr*170)
                elif add_glyc[1] == 2:
                    # m = update_glycogen(
                    #     m, fitted_glyc(-input1)
                    # )
                    m = update_prot_glyc(
                        m, proteincontent(-input1), fitted_glyc(-input1)
                    )
                elif add_glyc[1] == 3:
                    update_glycogen(m, 100*(0.0023*-vstar+0.0424))
        except NameError:
            pass
        if rxn_out == "BM0009":
            return growth(m, rxn_in, -1*vstar)
        elif isinstance(rxn_out, list):
            return get_fluxes(m, rxn_in, rxn_out, (-1*vstar, 0))
        elif fva:
            return get_fva(m, rxn_in, rxn_out, -1*vstar)
        else:
            return get_flux(m, rxn_in, rxn_out, -1*vstar)


def pd_fitting(
    model: cobra.Model, rxn_in: str,
    add_glyc: tuple = (False, 0), maint: float = 1.3
) -> Callable:
    '''takes model and input reaction id as arguments
    and returns the function to fit'''
    def pd_fba(
        input: np.typing.ArrayLike,
        pmax: float,
        alpha: float,
        alpha2: float,
        kd: float
    ) -> list:
        '''function to fit to the growth rate data
        input: array-like input data
        pmax, alpha, kd: parameters to fit'''
        # zero all N influxes:
        for rxn in model.boundary:
            if 'N' in rxn.check_mass_balance():
                rxn.lower_bound = 0.

        model.reactions.get_by_id(
            "EX_C00244_ext_b"  # nitrate
        ).lower_bound = -1000
        model.reactions.get_by_id(
            "EX_C00011_ext_b"  # co2
        ).lower_bound = -1000
        kl = get_kl(model, pmax)

        output = []
        for i1 in input:
            # i = alpha * i1
            i = plus_blue(i1, alpha, alpha2, 0)
            vstar = -1*i*kl/(i + kl)
            with model as m:
                atp = m.reactions.get_by_id("R_ATP_maint")
                atp.upper_bound = 1000
                photodamage = model.problem.Constraint(
                    kd/-100 * i + atp.flux_expression - maint,
                    lb=0,
                    ub=0
                )
                m.add_cons_vars(photodamage)
                try:
                    if add_glyc[0]:
                        if add_glyc[1] == 1:
                            for i in range(10):
                                gr = growth(m, rxn_in, vstar)
                                update_glycogen(m, 3.4+gr*170)
                        elif add_glyc[1] == 2:
                            # m = update_glycogen(
                            #     m, fitted_glyc(i1)
                            # )
                            # ind = list(input).index(i1)
                            m = update_prot_glyc(
                                m, proteincontent(i1), fitted_glyc(i1)
                            )
                        elif add_glyc[1] == 3:
                            update_glycogen(m, 100*(0.0023*-vstar+0.0424))
                except NameError:
                    print("could not update biomass")
                # print(m.reactions.BM0009.metabolites[m.metabolites.B10_cyt])
                output.append(growth(m, rxn_in, vstar))
        return output
    return pd_fba


def pd_double_fit(model: cobra.Model, rxn_in: str, rxns_out: list) -> Callable:
    '''takes model, input reaction id, and list of
    output reaction ids as arguments
    and returns the function to fit'''
    nitrate_id = "EX_C00244_ext_b"
    co2_id = "EX_C00011_ext_b"

    def pd_fba(input, pmax, alpha, kd):
        '''function to fit to the growth rate data
        input: array-like input data
        pmax, alpha, kd, kd2: parameters to fit'''
        # zero all N influxes:
        for rxn in model.boundary:
            if 'N' in rxn.check_mass_balance():
                rxn.lower_bound = 0.

        model.reactions.get_by_id(nitrate_id).lower_bound = -1000  # nitrate
        model.reactions.get_by_id(co2_id).lower_bound = -1000
        kl = get_kl(model, pmax)

        output1 = []
        output2 = []
        for i in input[:8]:  # only take first half of input, as it repeats
            i *= alpha
            vstar = -1*i*kl/(i + kl)
            with model as m:
                atp = m.reactions.get_by_id("R_ATP_maint")
                atp.upper_bound = 1000
                photodamage = model.problem.Constraint(
                    kd/-100 * i + atp.flux_expression - 1.3,
                    lb=0,
                    ub=0
                )
                m.add_cons_vars(photodamage)
                # nadph = m.reactions.get_by_id("NADPH_maint")
                # photodamage1 = model.problem.Constraint(
                #     kd2/-100 * i + nadph.flux_expression,
                #     lb=0,
                #     ub=0
                # )
                # m.add_cons_vars(photodamage1)
                out = get_fluxes(m, rxn_in, rxns_out, (vstar, vstar))
                output1.append(out[0]*100)
                if round(i/alpha, 2) in input[8:]:
                    output2.append(out[1])
                output = np.hstack((output1, output2))
        return output
    return pd_fba


def dark_resp(
    model,
    input1,
    kl,
    kd,
    alpha=1, alpha2=0, alpha3=0,
    add_glyc=(False, 0),
    maint=1.3
):
    '''return water splitting flux minus o2 export and mehler-like
    take model, kl and damage values'''
    with model as m4:
        ofluxes = [
            a*photodamage_helper(
                m4, "EX_E1_ext_b", id,
                input1, kl, kd, alpha, alpha2, alpha3, add_glyc, maint
            )
            for a, id in zip([-1, 1, 1], [
                "PR0043",
                "EX_C00007_ext_b",
                "PR0033"
                # "PR0032",
                # "PR0034"
            ])
        ]
    return sum(ofluxes)


def pd_fitting_w_resp(
    model: cobra.Model, rxn_in: str, add_glyc: tuple = (False, 0)
) -> Callable:
    '''takes model and input reaction id as arguments
    and returns the function to fit'''
    def pd_fba(
        input: np.typing.ArrayLike,
        pmax: float,
        alpha: float,
        alpha2: float,
        kd: float,
        res: float,
        maint: float
    ) -> list:
        '''function to fit to the growth rate data
        input: array-like input data
        pmax, alpha, kd: parameters to fit'''
        # zero all N influxes:
        for rxn in model.boundary:
            if 'N' in rxn.check_mass_balance():
                rxn.lower_bound = 0.

        model.reactions.get_by_id(
            "EX_C00244_ext_b"  # nitrate
        ).lower_bound = -1000
        model.reactions.get_by_id(
            "EX_C00011_ext_b"  # co2
        ).lower_bound = -1000
        kl = get_kl(model, pmax)

        output = []
        respir = []
        inputs = []
        for i1 in input:
            # i = alpha * i1
            i = plus_blue(i1, alpha, alpha2, 0)
            vstar = -1*i*kl/(i + kl)
            with model as m:
                atp = m.reactions.get_by_id("R_ATP_maint")
                atp.upper_bound = 1000
                photodamage = m.problem.Constraint(
                    kd/-100 * i + atp.flux_expression - maint,
                    lb=0,
                    ub=0
                )
                resp_rxns = [
                    m.reactions.get_by_id(X) for X in [
                        "PR0010", "PR0010_2", "PR0011", "PR0011_2"
                    ]
                ]
                use_sum = sum([x.flux_expression for x in resp_rxns])
                some_flux = m.problem.Constraint(
                    use_sum - res/100 * m.reactions.PR0043.flux_expression,
                    lb=0,
                    ub=1000,
                    name="respiration"
                )
                # mehler_like = m.problem.Constraint(
                #     m.reactions.PR0033.flux_expression -
                #     res/100 * m.reactions.PR0043.flux_expression,
                #     lb=0,
                #     ub=1000,
                #     name="mehlerlike"
                # )
                m.add_cons_vars([photodamage, some_flux])  # , mehler_like])
                try:
                    if add_glyc[0]:
                        if add_glyc[1] == 1:
                            for i in range(10):
                                gr = growth(m, rxn_in, vstar)
                                update_glycogen(m, 3.4+gr*170)
                        elif add_glyc[1] == 2:
                            m = update_prot_glyc(
                                m, proteincontent(i1), fitted_glyc(i1)
                            )
                        elif add_glyc[1] == 3:
                            update_glycogen(m, 100*(0.0023*-vstar+0.0424))
                except NameError:
                    print("could not update biomass")
                if i1 not in inputs:
                    output.append(growth(m, rxn_in, vstar)*1000)
                else:
                    respir.append(
                        dark_resp(
                            m, -i1, kl, kd, alpha, alpha2, 0, add_glyc, maint
                        )
                    )
                inputs.append(i1)
        return output+respir
    return pd_fba


def plus_blue(input, alpha, alpha2, alpha3=False):
    if alpha3:
        I_star = (I0 := input-27.5) * alpha + alpha2 * (1 + alpha3 * I0) * 27.5
    else:
        I_star = (input-27.5 + alpha2 * 27.5) * alpha
    return I_star * 3.6


def fit_haldane(input, mu_star, ka, gamma):
    # define the haldane function:
    growthrate = mu_star * input/(ka + input + gamma * (input**2)/ka)

    return growthrate
