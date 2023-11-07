import numpy as np
import cobra
import matplotlib.pyplot as plt
import multiprocessing as mp


def growth(model: cobra.Model, rxn_id: str, influx: float) -> float:
    '''calculate objective value in dependency of influx
    model: a cobra model to simulate
    rxn_id: the id of the reaction to vary
    influx: influx value to use

    returns: optimized objective value'''

    with model:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = influx
        return model.slim_optimize()


def get_flux(
    model: cobra.Model,
    rxn_in_id: str, rxn_out_id: str,
    influx: float
) -> float:
    '''calculate a flux in dependency of influx value
    model: a cobra model to simulate
    rxn_in_id: id of the reaction to vary
    rxn_out_id: id of the reaction which flux to calculate
    influx: influx value

    returns: flux value for rxn_out_id'''

    with model:
        rxn_in = model.reactions.get_by_id(rxn_in_id)
        # rxn_out = model.reactions.get_by_id(rxn_out_id)
        rxn_in.lower_bound = influx
        try:
            psol = cobra.flux_analysis.pfba(model, reactions=[rxn_out_id])
            # psol = model.optimize()
            return psol.fluxes[rxn_out_id]
        except cobra.exceptions.Infeasible:
            return 0


def get_fluxes(
    model: cobra.Model,
    rxn_in_id: str, rxn_out_ids: str,
    bounds: tuple
) -> tuple:
    '''calculate a set of fluxes in dependency of bounds for a reaction
    model: a cobra model to simulate
    rxn_in_id: id of the reaction to constrain
    rxn_out_ids: a list of ids of the reactions which fluxes to calculate
    bounds: tuple with lower_bound, upper_bound constrains

    returns: a tuple with flux values for rxn_out_ids'''

    with model:
        rxn_in = model.reactions.get_by_id(rxn_in_id)
        rxn_in.bounds = bounds
        try:
            psol = cobra.flux_analysis.pfba(
                model,
                reactions=rxn_out_ids,
                fraction_of_optimum=0.999
            )
            # psol = model.optimize()
            return tuple(psol.fluxes[rxn_out_id] for rxn_out_id in rxn_out_ids)
        except cobra.exceptions.Infeasible:
            return tuple(0 for rxn_out_id in rxn_out_ids)


def get_fva(
    model: cobra.Model,
    rxn_in_id: str, rxn_out_id: str,
    influx: float, p: float = None
) -> tuple:
    '''build a pair of flux minimum value and flux maximum value
    in dependency of influx value
    model: a cobra model to simulate
    rxn_in_id: id of the reaction to vary
    rxn_out_id: id of the reaction which flux variability to calculate
    influx: influx value
    p (optional): the pfba factor for the variability

    returns: tuple with (minimum, maximum) of flux for rxn_out_id'''

    with model:
        rxn_in = model.reactions.get_by_id(rxn_in_id)
        rxn_out = model.reactions.get_by_id(rxn_out_id)
        rxn_in.lower_bound = influx
        try:
            vsol = cobra.flux_analysis.flux_variability_analysis(
                model, rxn_out, pfba_factor=p, loopless=True
            )
        except cobra.exceptions.Infeasible:
            return (0, 0)
        return (vsol.minimum.iloc[0], vsol.maximum.iloc[0])


if __name__ == "__main__":
    from import_model import import_model
    model = import_model('../models_sbml/Synechocystis6803_1.sbml')

    photons = np.arange(0., 15, 0.001)
    photon_id = "EX_E1_ext_b"
    o2_id = "EX_C00007_ext_b"
    pool = mp.Pool(mp.cpu_count()) # initialize processor pool for parallel programming

    # run the simulations for plotting BOF and o2 production in dependency of photon influx:
    with model as m1:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        m1.reactions.get_by_id("EX_C00288_ext_b").lower_bound = -1 # influx of HCO3-
        with m1 as m1_1:
            biomass = pool.starmap(growth, ((m1_1, photon_id, p*-1) for p in photons))
        with m1 as m1_2:
            o2 = pool.starmap(get_flux, ((m1_2, photon_id, o2_id, p*-1) for p in photons))
        ax1.plot(photons, biomass, 'r-', label = "growth, 1 mmol/(gDW*h) HCO3")
        ax2.plot(photons, o2, 'r--', label = "O2, 1 mmol/(gDW*h) HCO3")
        ax1.set_xlabel("Photons [mmol/(gDW*h)]")
        ax1.set_ylabel("growth rate [1/h]")
        ax2.set_ylabel("O2 export [mmol/(gDW*h)]")
        m1.reactions.EX_C00031_ext_b.lower_bound = -0.1 # influx of glucose
        with m1 as m1_3:
            biomass = pool.starmap(growth, ((m1_3, photon_id, p*-1) for p in photons))
        with m1 as m1_4:
            o2 = pool.starmap(get_flux, ((m1_4, photon_id, o2_id, p*-1) for p in photons))
        ax1.plot(photons, biomass, 'g-', label = "+ 0.1 mmol Glucose")
        ax2.plot(photons, o2, 'g--', label = "O2, + 0.1 mmol Glucose")
        ax1.set_xlim(0)
        ax1.set_ylim([0,0.04])
        ax2.set_ylim([-.4,1.2])
        ax1.legend()
        ax2.legend(loc = 'lower right')
        fig.savefig('growth_mixo_carbon_lim.png')

    # run the simulations for plotting BOF and o2 production in dependency of carbon influx with different N sources:
    with model as m2:
        # zero all N influxes:
        for rxn in m2.boundary:
           if 'N' in rxn.check_mass_balance():
               rxn.lower_bound = 0.

        carbon = np.arange(0, 3., 0.0001)
        fig, ax1 = plt.subplots()
        m2.reactions.get_by_id(photon_id).lower_bound = -1000
        with m2 as m2_1:
            m2_1.reactions.get_by_id("EX_C01342_ext_b").lower_bound = -.3 # influx of NH4
            biomass = pool.starmap(growth, ((m2_1, "EX_C00288_ext_b", c*-1) for c in carbon))
        with m2 as m2_2:
            m2_2.reactions.get_by_id("EX_C00086_ext_b").lower_bound = -.3 # influx of urea
            biomass_urea = pool.starmap(growth, ((m2_1, "EX_C00288_ext_b", c*-1) for c in carbon))
        ax1.plot(carbon, biomass, 'r-', label = "growth, 0.3 mmol/(gDW*h) NH4")
        ax1.plot(carbon, biomass_urea, 'b-', label = "growth, 0.3 mmol/(gDW*h) Urea")
        ax1.set_xlabel("CO2 influx [mmol/(gDW*h)]")
        ax1.set_ylabel("growth rate [1/h]")
        ax1.set_xlim(0)
        ax1.set_ylim(0)
        ax1.legend()
        fig.savefig('growth_nitro_lim.png')

    # run the simulations for plotting min/max rubisco flux dependency of carbon influx:
    with model as m3:
        # zero all N influxes:
        for rxn in m3.boundary:
            if 'N' in rxn.check_mass_balance():
                rxn.lower_bound = 0.

        m3.reactions.get_by_id(photon_id).lower_bound = -20
        carbon = np.arange(0, 1.5, 0.001)
        fig, ax1 = plt.subplots()
        with m3 as m3_1:
            m3_1.reactions.get_by_id("EX_C01342_ext_b").lower_bound = -.3 # influx of NH4
            rubisco = pool.starmap(get_fva, ((m3_1, "EX_C00288_ext_b", "R431", c*-1) for c in carbon))
        ax1.plot(carbon, [t[0] for t in rubisco], 'g-', label = "rubisco, 0.3 mmol/(gDW*h) N (NH4)")
        ax1.plot(carbon, [t[1] for t in rubisco], 'g-')
        ax1.set_xlabel("CO2 influx [mmol/(gDW*h)]")
        ax1.set_ylabel("rubsico flux")
        ax1.set_xlim(0)
        ax1.set_ylim(0)
        ax1.legend()
        fig.savefig('growth_rubisco_fva_nitro_lim.png')
