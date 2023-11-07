from typing import Callable, Union
from modules.FBA_S6803 import growth, get_flux, get_fluxes, get_fva
from modules.import_model import import_model, split_rubisco
from modules.import_elife import DataImporter
import numpy as np
import numpy.typing
import matplotlib.pyplot as plt
from scipy import optimize as op
import multiprocessing as mp
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

    lights_in = zavrel.get("light")
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
                        m, proteincontent(-input), fitted_glyc(-input)
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
        I_star = (input-27.5 + alpha2 * 27.5) * alpha * 3.6
    return I_star


def fit_haldane(input, mu_star, ka, gamma):
    # define the haldane function:
    growthrate = mu_star * input/(ka + input + gamma * (input**2)/ka)

    return growthrate


if __name__ == "__main__":

    new_syn = import_model(
            '../models_sbml/Synechocystis6803_SBML_COBRA_24_4_17.sbml'
        )
    # knock out cell wall atpase, to prevent "atp from nothing" production:
    new_syn.reactions.PR0045.knock_out()
    split_rubisco(new_syn, 300/97)
    add_glyc = (True, 2)
    new_syn = update_glycogen(new_syn, 3.41)
    # new_syn = update_prot_glyc(new_syn, 30, 3.41)

    # add nadph maintenance reaction:
    nadph = cobra.Reaction('NADPH_maint')
    nadph.name = 'NADPH maintenance reaction'
    nadph.subsystem = 'cell maintenance'
    nadph.lower_bound = 0
    nadph.upper_bound = 1000

    nadph.add_metabolites({
        new_syn.metabolites.get_by_id('C00005_cyt'): -1.0,
        # new_syn.metabolites.get_by_id('C00007_cyt'): -2.0,
        new_syn.metabolites.get_by_id('C00080_cyt'): 1.0,
        # new_syn.metabolites.get_by_id('C00704_cyt'): 2.0,
        new_syn.metabolites.get_by_id('C00006_cyt'): 1.0
    })

    # new_syn.add_reactions([nadph])
    # print(nadph.check_mass_balance())

    nitrate_id = "EX_C00244_ext_b"
    photon_id = "EX_E1_ext_b"
    o2_id = "EX_C00007_ext_b"
    ps2_id = "PR0043"
    co2_id = "EX_C00011_ext_b"
    bof_id = "BM0009"

    light = zavrel.get("light")  # [4:]
    light1 = light[:-1]  # + light[-2:-1]
    growth_rate = zavrel.get("growth")  # [4:])
    o2_data = zavrel.get("O2")
    chl_data = zavrel.get("chl")
    o2_data_calc = zavrel.conv_for_fba(o2_data)
    co2_data_calc = zavrel.conv_for_fba(zavrel.get("CO2"))
    dark_resp_data = zavrel.conv_for_fba(zavrel.get("dark resp"))
    net_o2 = zavrel.value("net O2 fba")

    light_calc_with_error = (
        (zavrel.get("light cap with error") + blue_light * zavrel.get("light cap with error")
         /
         zavrel.get("light with error") * 0.66)
        /
        (zavrel.get("DW with error") * 0.02) * (60*60) / 1000
    )
    light_calc = np.array([l.n for l in light_calc_with_error])
    # light_calc = zavrel.get("light cap")

    fig, ax = plt.subplots()
    ax.plot(light, light_calc, "o")
    ax.plot(light, 0.58*light)
    ax.set_ylabel("mmol Photons/gDW/h")
    ax.set_xlabel("µE/m²/s")
    fig.savefig("zavrel/photonsperdw.svg", bbox_inches="tight")

    (pmax, alpha, alpha2, kd), _ = op.curve_fit(
        pd_fitting(new_syn, photon_id, add_glyc), light+blue_light, growth_rate,
        p0=(150, .5, .66, 6), bounds=((20, 0, 0, 0.01), (300, 1, 1, 100))
    )
    # alpha2 = alpha
    alpha3 = 0

    print(pmax, alpha, alpha2, alpha3, kd)
    new_syn = update_prot_glyc(new_syn)

    # plot the linear and calculated alpha:
    alpha_with_error = (
        light_calc_with_error / zavrel.get("light with error")
    )

    fig, ax = plt.subplots()
    ax.plot(light, [alpha for i in light], "k-", label="linear fitted alpha")
    ax.plot(light, light_calc/light, "kX", label="calculated alpha")
    ax.set_xlabel("Incident light [µE/m²/s]")
    ax.set_ylabel("alpha")
    ax.errorbar(
        light, light_calc/light, fmt="none",
        yerr=[a.s for a in alpha_with_error],
        elinewidth=.5, capsize=2, ecolor="k"
    )
    ax.set_ylim(0)
    ax.set_xlim(0)
    ax.legend()
    fig.savefig("zavrel/alphas.svg", bbox_inches="tight")
    fig.savefig("zavrel/alphas.png", bbox_inches="tight")

    fig, ax = plt.subplots()
    ax.plot(light, plus_blue(light, alpha, alpha2), "k-", label="linear fitted photons")
    ax.plot(light, light_calc, "kX", label="calculated photons")
    ax.set_xlabel("Incident light [µE/m²/s]")
    ax.set_ylabel("Photons [mmol/gDW/h]")
    ax.errorbar(
        light, light_calc, fmt="none",
        yerr=[l.s for l in light_calc_with_error],
        elinewidth=.5, capsize=2, ecolor="k"
    )
    ax.set_ylim(0)
    ax.set_xlim(0)
    ax.legend()
    fig.savefig("zavrel/photons.svg", bbox_inches="tight")
    fig.savefig("zavrel/photons.png", bbox_inches="tight")

    photons = np.linspace(0, light[-1]*1.1, 1000) + blue_light  # *alpha
    fig, ax = plt.subplots()
    ax.plot(photons, [fitted_glyc(p) for p in photons], "r")
    ax.plot(light, zavrel.get("glycogen")/zavrel.get("DW") * 100, "ro")
    fig.savefig("zavrel/glycogen_fitted.svg", bbox_inches="tight")
    with new_syn as m3:
        atp = m3.reactions.get_by_id("R_ATP_maint")
        atp.upper_bound = 1000
        atp.lower_bound = 1.3
        with mp.Pool(mp.cpu_count()-2) as pool:
            kl = get_kl(m3, pmax)
            with m3 as m4:

                stoich_gro = pool.starmap(
                    growth,
                    (
                        (update_glycogen(m4, fitted_glyc(p)), photon_id, plus_blue(p, alpha, alpha2, alpha3)*-1)
                        for p in photons
                    )
                )
            gro = pool.starmap(
                photodamage_helper,
                (
                    (m3, photon_id, "BM0009", p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc)
                    for p in photons
                )
            )
            fluxes = pool.starmap(
                photodamage_helper,
                (
                    (m3, photon_id, [o2_id, co2_id, ps2_id], p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc)
                    for p in photons
                )
            )
            o2 = np.array([f[0] for f in fluxes])
            co2 = np.array([f[1] for f in fluxes])
            h2osplit = np.array([f[2] for f in fluxes])
            o2_points = []
            gro1 = []
            for i in range(len(phots := zavrel.get("light")+blue_light)):
                with m3 as model3:
                    model3 = update_prot_glyc(model3, zavrel.get('protein')[i]/zavrel.get('DW')[i]*100, zavrel.get('glycogen')[i]/zavrel.get('DW')[i]*100)
                    model3.reactions.EX_E1_ext_b.lower_bound = -phots[i]
                    o2_points.append(photodamage_helper(model3, photon_id, o2_id, -phots[i], kl, kd, alpha, alpha2, alpha3, maint=1.3))
                    gro1.append(photodamage_helper(model3, photon_id, "BM0009", -phots[i], kl, kd, alpha, alpha2, alpha3, maint=1.3))
            # co2 = pool.starmap(
            #     photodamage_helper,
            #     (
            #         (m3, photon_id, co2_id, p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc)
            #         for p in photons
            #     )
            # )
            # splitminmax = pool.starmap(
            #     photodamage_helper,
            #     (
            #         (m3, photon_id, ps2_id, p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc, 1.3, True)
            #         for p in photons
            #     )
            # )
            # o2minmax = pool.starmap(
            #     photodamage_helper,
            #     (
            #         (m3, photon_id, o2_id, p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc, 1.3, True)
            #         for p in photons
            #     )
            # )
            o2resp = pool.starmap(
                dark_resp,
                ((m3, p*-1, kl, kd, alpha, alpha2, alpha3, add_glyc) for p in photons)
            )
            # o2resp = pool.starmap(
            #     dark_met_o, ((m3, i, p, add_glyc) for i, p in zip(gro, photons))
            # )
    print("YIELD:", stoich_gro[50]/photons[50])

    # Plot Growth rate and O2:
    light_calc -= blue_light
    # photons = np.linspace(0, light[-1]*1.1, 1000)
    fig, ax = plt.subplots()
    p1 = ax.plot(plus_blue(photons, alpha, alpha2, alpha3), gro, 'k-', label="fitted growth curve")
    p2 = ax.plot(plus_blue(light+blue_light, alpha, alpha2, alpha3), gro1, "kX", label="fitted with protein")
    p2 = ax.plot(plus_blue(light+blue_light, alpha, alpha2, alpha3), growth_rate, "ko", label="zavrel data")
    ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), growth_rate, fmt="none",
        # xerr=zavrel.error("light calc"),
        yerr=zavrel.error("growth"),
        elinewidth=.5, capsize=2, ecolor="k"
    )
    ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    ax.set_ylabel("Growth rate [$h^{-1}$]")
    ax.set_xlim(0)
    ax.set_ylim(0)
    p5 = ax.plot(plus_blue(photons, alpha, alpha2, alpha3), stoich_gro, "k--", label="stoichiometric growth")
    ax.text(250, 0.03, (
        f'K$_L$ = {round(kl, 2)} $\\frac{{mmol\\ \\mathrm{{photons}}}}'
        f'{{gDW\\cdot h}}$\n'
        f'$\\alpha$ = {round(alpha, 2)} $\\frac{{m^2}}{{gDW}}$\n'
        f'$\\alpha_2$ = {round(alpha2, 2)} $\\frac{{m^2}}{{gDW}}$\n'
        f'$\\alpha_3$ = {round(alpha3, 2)} $\\frac{{m^2}}{{gDW}}$\n'
        f'K$_d$ = {round(kd, 2)} %\n'
        # f'K$_{{d,NADPH}}$ = {round(kd2, 2)} %'
    ))
    ax.legend()
    fig.savefig("zavrel/fit_growth.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_growth_error.png", bbox_inches="tight")
    ax.clear()
    ax.plot(
        plus_blue(photons, alpha, alpha2), gro/(plus_blue(photons, alpha, alpha2))*1000,
        'b-', label="yield from fitted growth"
    )
    ax.plot(
        plus_blue(light+blue_light, alpha, alpha2), growth_rate/(plus_blue(light+blue_light, alpha, alpha2))*1000,
        "bo", label="zavrel data"
    )
    yield_error = np.array(
        [
            y.s for y in
            zavrel.get("growth with error")
            /
            (plus_blue(zavrel.get("light with error")+blue_light, alpha, alpha2))*1000
        ]
    )
    ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2), growth_rate/(plus_blue(light+blue_light, alpha, alpha2))*1000, fmt="none",
        # xerr=zavrel.error("light calc"),
        yerr=yield_error,
        # yerr=zavrel.error("growth"),
        elinewidth=.5, capsize=2, ecolor="b"
    )
    ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    ax.set_ylabel("Yield [$\\frac{gDW}{mol_{Photons}}$]")
    ax.set_xlim(0)
    ax.plot(
        plus_blue(photons, alpha, alpha2), stoich_gro/(plus_blue(photons, alpha, alpha2))*1000,
        "b--", label="stoichiometric yield"
    )
    ax.set_ylim(0)
    ax.set_xlim(0)
    ax.legend(loc="center right")
    fig.savefig("zavrel/fit_gro_yield.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_yield_error.png", bbox_inches="tight")
    ax.clear()
    # ax2 = ax.twinx()
    fac = 1
    # calculate a factor to "scale" the y-axis:
    # fac = np.mean([
    #     o2_data[i] / o2[int(np.where(np.arange(0,1200, 0.1) == light[i])[0])]
    #     for i in range(len(light[:-1]))
    # ])
    p3 = ax.plot(
        plus_blue(photons, alpha, alpha2, alpha3), [o*fac for o in h2osplit],
        '-', color="limegreen", label="O$_2$ evolution at PS2"
    )
    ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW\\cdot h}$]")
    ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    p4 = ax.plot(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), o2_data_calc,
        "o", color='limegreen', label="zavrel data O$_2$"
    )
    p5 = ax.plot(
        plus_blue(photons, alpha, alpha2, alpha3), h2osplit,
        '--', color="limegreen", label="O$_2$ evolution at ps2"
    )
    # p5 = ax.plot(
    #     plus_blue(photons, alpha, alpha2, alpha3), [o[0] for o in splitminmax],
    #     '--', color="limegreen", label="O$_2$ evolution variability"
    # )
    # ax.plot(
    #     plus_blue(photons, alpha, alpha2, alpha3), [o[1] for o in splitminmax],
    #     '--', color="limegreen"
    # )
    ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), o2_data_calc, fmt="none",
        # xerr=zavrel.error("light calc"),
        yerr=zavrel.error("O2 fba"),
        elinewidth=.5, capsize=2, ecolor="limegreen"
    )
    ax.set_ylim(0)
    ax.set_xlim(0)
    p = p3 + p4 + p5
    labels = [lab.get_label() for lab in p]
    ax.legend(p, labels)
    ax.set_title("Growth rate fitted to data")

    fig.savefig("zavrel/fit_growth_o.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_growth_o_error.png", bbox_inches="tight")

    ax.clear()
    p3 = ax.plot(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), [o*fac for o in o2_points],
        'X', color="limegreen", label="O$_2$ production"
    )
    p6 = ax.plot(
        plus_blue(photons, alpha, alpha2, alpha3), [o*fac for o in o2],
        color="limegreen", label="O$_2$ production"
    )
    ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW\\cdot h}$]")
    ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    p4 = ax.plot(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), net_o2,
        "o", color='limegreen', label="zavrel data net O$_2$"
    )
    p5 = ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), net_o2, fmt="none",
        # xerr=zavrel.error("light calc"),
        yerr=zavrel.error("net O2 fba"),
        elinewidth=.5, capsize=2, ecolor="limegreen"
    )
    # p3 = ax.plot(
    #     plus_blue(photons, alpha, alpha2, alpha3), [o[0] for o in o2minmax],
    #     '--', color="limegreen", label="O$_2$ production variability"
    # )
    # ax.plot(
    #     plus_blue(photons, alpha, alpha2, alpha3), [o[1] for o in o2minmax],
    #     '--', color="limegreen"
    # )
    ax.set_ylim(0, 8)
    ax.set_xlim(0)
    p = p3 + p4
    labels = [lab.get_label() for lab in p]
    ax.legend(p, labels)
    ax.set_title("Growth rate fitted to data")

    fig.savefig("zavrel/fit_growth_net_o.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_growth_net_o_error.png", bbox_inches="tight")

    # Plot growth rate and CO2:
    ax.clear()
    p5 = ax.plot(
        plus_blue(photons, alpha, alpha2, alpha3), -1*np.array(co2),
        '-', color="red", label="CO$_2$ uptake"
    )
    p6 = ax.plot(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), co2_data_calc,
        "o", color='red', label="zavrel data CO$_2$"
    )
    p7 = ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), co2_data_calc, fmt="none",
        yerr=zavrel.error("CO2 fba"),
        elinewidth=.5, capsize=2, ecolor="red"
    )
    ax.set_ylabel("CO$_2$ flux [$\\frac{mmol}{gDW}$]")
    ax.set_ylim(0)
    ax.set_xlim(0)
    p = p5 + p6
    labels = [lab.get_label() for lab in p]
    ax.legend(p, labels)

    fig.savefig("zavrel/fit_growth_c.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_growth_c_error.png", bbox_inches="tight")

    # Plot O2 and CO2 dependent of growth rate:
    fig, ax = plt.subplots()
    # p1 = ax.plot(
    #     gro, o2,
    #     '-', color="limegreen", label="O$_2$"
    # )
    # p2 = ax.plot(
    #     growth_rate, o2_data_calc,
    #     'o', color="limegreen", label="O$_2$ data"
    # )
    # ax.set_xlim(0)
    # ax.set_ylim([0, 6])
    # ax.set_xlabel("Growth rate [$h^{-1}$]")
    # ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW\\cdot h}$]")
    # ax2 = ax.twinx()
    # p3 = ax2.plot(
    #     gro, -1*np.array(co2),
    #     '-', color="red", label="CO$_2$"
    # )
    # p2 = ax2.plot(
    #     growth_rate, co2_data_calc,
    #     'o', color="red", label="CO$_2$ data"
    # )
    # ax2.set_ylim(0, 6)
    # ax.set_xlim(0)
    # ax2.set_ylabel("CO$_2$ flux [$\\frac{mmol}{gDW}$]")
    # ax.set_title("Growth rate fitted to data")

    # fig.savefig("zavrel/fit_growth_o2out_co2.svg", bbox_inches="tight")
    # fig.savefig("zavrel/fit_growth_o2out_co2.png", bbox_inches="tight")

    # Plot dark respiration:
    # o2resp = o2resp_fac * np.array(gro)

    fig, ax = plt.subplots()
    fig.set_figheight(3)
    p1 = ax.plot(
        plus_blue(photons, alpha, alpha2, alpha3), o2resp,  # o2-h2osplit,
        '-', color="darkgreen", label="O$_2$ Respiration"
    )
    p2 = ax.plot(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), dark_resp_data,
        "o", color="darkgreen", label="zavrel data"
    )
    ax.errorbar(
        plus_blue(light+blue_light, alpha, alpha2, alpha3), dark_resp_data, fmt="none",
        # xerr=zavrel.error("light calc"),
        yerr=zavrel.error("dark resp fba"),
        elinewidth=.5, capsize=2, ecolor="darkgreen"
    )
    ax.set_xlim(0)
    ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW\\cdot h}$]")
    ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    p = p1 + p2
    labels = [lab.get_label() for lab in p]
    ax.legend(p, labels)
    ax.set_title("Growth rate fitted to data")

    fig.savefig("zavrel/fit_gro_resp.svg", bbox_inches="tight")
    fig.savefig("zavrel/fit_gro_resp_error.png", bbox_inches="tight")

# -------------------------------------------------------------------
# Fit Growth Rate And O2:
# -------------------------------------------------------------------

    # # print([o2_data[i]/growth_rate[i] for i in range(len(o2_data))])
    # factor = 1#np.mean([o2_data_calc[i]/growth_rate[i] for i in range(len(o2_data_calc))]) # scale data to match, for fitting 
    # o2_data_calc /= factor
    # (pmax, alpha, kd, kd2), _ = op.curve_fit(
    #     pd_double_fit(
    #         new_syn, photon_id, [bof_id, o2_id]
    #     ),
    #     np.hstack((light, light1)), np.hstack((growth_rate*100, o2_data_calc)),
    #     p0=(90, 0.5, 1, 1), bounds=((20, 0.2, 0.01, 0), (300, 10, 10, 10))
    # )

    # print(pmax, alpha, kd, kd2)

    # photons = np.arange(0, 1200, 0.1)*alpha
    # with new_syn as m3:
    #     atp = m3.reactions.get_by_id("R_ATP_maint")
    #     atp.upper_bound = 1000
    #     with mp.Pool(mp.cpu_count()-2) as pool:
    #         kl = get_kl(m3, pmax)
    #         gro = pool.starmap(photodamage_helper,
    #             ((m3, photon_id, "BM0009", p*-1, kl, kd, kd2) for p in photons)
    #         )
    #         o2 = pool.starmap(
    #             photodamage_helper,
    #             ((m3, photon_id, o2_id, p*-1, kl, kd, kd2) for p in photons)
    #         )
    #         co2 = pool.starmap(
    #             photodamage_helper, 
    #             ((m3, photon_id, co2_id, p*-1, kl, kd, kd2) for p in photons)
    #         )

    # # Plot Growth rate and O2:
    # fig, ax = plt.subplots()
    # p1 = ax.plot(photons/alpha, gro, 'k-', label="fitted growth curve")
    # p2 = ax.plot(light, growth_rate, "ko", label="zavrel data")
    # ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    # ax.set_ylabel("Growth rate [$h^{-1}$]")
    # ax.set_xlim(0)
    # ax.set_ylim(0)
    # ax.text(500, 0.05, (
    #     f'K$_L$ = {round(kl, 2)} $\\frac{{mmol\\ \\mathrm{{photons}}}}{{gDW\\cdot h}}$\n'
    #     f'$\\alpha$ = {round(alpha, 2)} $\\frac{{m^2}}{{gDW}}$\n'
    #     f'K$_{{d,ATP}}$ = {round(kd, 2)} %\n'
    #     f'K$_{{d,NADPH}}$ = {round(kd2, 2)} %'
    # ))
    # ax2 = ax.twinx()
    # p3 = ax2.plot(photons/alpha, [o*factor for o in o2], '-', color="limegreen", label="O$_2$ evolution")
    # ax2.set_ylim(0)
    # ax2.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW}$]")
    # p4 = ax2.plot(light1, [o*factor for o in o2_data_calc], "o", color="limegreen", label="zavrel data O$_2$")
    # p = p1 + p2 + p3 + p4
    # labels = [l.get_label() for l in p]
    # ax.legend(p, labels)
    # ax.set_title("Growth rate and O$_2$ fitted to data")

    # fig.savefig("zavrel/fit_groO2_o.svg", bbox_inches="tight")
    # fig.savefig("zavrel/fit_gro02_o.png", bbox_inches="tight")

    # # Plot growth rate and CO2:
    # ax2.clear()
    # p5 = ax2.plot(photons/alpha, -1*np.array(co2), '-', color="red", label="CO$_2$ uptake")
    # p6 = ax2.plot(light, co2_data_calc, "o", color='red', label="zavrel data CO$_2$")
    # ax2.set_ylabel("CO$_2$ flux [$\\frac{mmol}{gDW}$]")
    # ax2.set_ylim(0)
    # p = p1 + p2 + p5 + p6
    # labels = [l.get_label() for l in p]
    # ax.legend(p, labels)

    # fig.savefig("zavrel/fit_groO2_c.svg", bbox_inches="tight")
    # fig.savefig("zavrel/fit_groO2_c.png", bbox_inches="tight")

    # # Plot O2 and CO2 dependent of growth rate:
    # fig, ax = plt.subplots()
    # p1 = ax.plot(gro, o2, '-', color="limegreen", label="O$_2$")
    # p2 = ax.plot(growth_rate[:-1], o2_data_calc, 'o', color="limegreen", label="O$_2$ data")
    # ax.set_xlim(0)
    # ax.set_ylim([0,6])
    # ax.set_xlabel("Growth rate [$h^{-1}$]")
    # ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW}$]")
    # ax2 = ax.twinx()
    # p3 = ax2.plot(gro, -1*np.array(co2), '-', color="red", label="CO$_2$")
    # p2 = ax2.plot(growth_rate, co2_data_calc, 'o', color="red", label="CO$_2$ data")
    # ax2.set_ylim(0,6)
    # ax2.set_ylabel("CO$_2$ flux [$\\frac{mmol}{gDW}$]")
    # ax.set_title("Growth rate and O$_2$ fitted to data")

    # fig.savefig("zavrel/fit_groO2_o2out_co2.svg", bbox_inches="tight")
    # fig.savefig("zavrel/fit_groO2_o2out_co2.png", bbox_inches="tight")

    # # Plot dark respiration:
    # o2resp = o2resp_fac * np.array(gro)
    
    # fig, ax = plt.subplots()
    # p1 = ax.plot(photons/alpha, o2resp, '-', color="darkgreen", label="O$_2$ in Biomass")
    # p2 = ax.plot(light, dark_resp_data, "o", color="darkgreen", label="zavrel data")
    # ax.set_ylabel("O$_2$ flux [$\\frac{mmol}{gDW}$]")
    # ax.set_xlabel("Photon flux [$mmol\\cdot gDW^{-1} h^{-1}$]")
    # p = p1 + p2
    # labels = [l.get_label() for l in p]
    # ax.legend(p, labels)
    # ax.set_title("Growth rate and O$_2$ fitted to data")

    # fig.savefig("zavrel/fit_groO2_resp.svg", bbox_inches="tight")
    # fig.savefig("zavrel/fit_groO2_resp.png", bbox_inches="tight")

