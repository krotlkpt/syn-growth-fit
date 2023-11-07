import cobra


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
