import cobra


def import_model(filename):
    '''
    import a sbml model and add Biomass objective (R90),
    ATP maintanace reaction,
    Cytochrome oxidase respiration and set all Carbon influx to 0
    Parameters: filename - the location of the model to import
    Returns: The model as a cobrapy model
    '''
    model = cobra.io.read_sbml_model(filename)
    if len(model.reactions.query("R90")) > 0:
        old, new, latest = (True, False, False)
    elif len(model.metabolites.query("C00244_cyt")) > 0:
        old, new, latest = (False, True, False)
    else:
        old, new, latest = (False, False, True)

    # define Biomass and cytochrome respiratroy reactions:
    if old:
        bm_rxn = model.reactions.R90
        # resp_rxns = [mdodel.reactions.get_by_id(X) for X in ["R458", "R459"]]
        O2_rxns = [model.reactions.get_by_id(X) for X in ["R491"]]
    else:
        bm_rxn = model.reactions.BM0009
        # resp_rxns = [
        #     model.reactions.get_by_id(X) for X in [
        #         "PR0010", "PR0010_2", "PR0011", "PR0011_2"
        #     ]
        # ]
        O2_rxns = [model.reactions.get_by_id(X) for X in ["PR0043"]]
        if new:
            model.reactions.TR0050.add_metabolites({
                "C00244_cyt": -1,
                "C00080_cyt": 1
            })

    bm_rxn.objective_coefficient = 1.0

    if new:
        # add atp maintenance reaction:
        atp = cobra.Reaction('R_ATP_maint')
        atp.name = 'ATP maintenance reaction'
        atp.subsystem = 'cell maintenance'
        # atp.lower_bound = 1.3
        atp.lower_bound = 0
        atp.upper_bound = 1000

        atp.add_metabolites({
            model.metabolites.get_by_id('C00002_cyt'): -1.0,
            model.metabolites.get_by_id('C00001_cyt'): -1.0,
            model.metabolites.get_by_id('C00008_cyt'): 1.0,
            model.metabolites.get_by_id('C00009_cyt'): 1.0
        })

        model.add_reactions([atp])

    # define constrain 20% of evolved O2 in respiration
    prod_sum = sum([x.flux_expression for x in O2_rxns])
    # use_sum = sum([x.flux_expression for x in resp_rxns])
    # some_flux = model.problem.Constraint(
    #     use_sum - 0.05 * prod_sum,
    #     lb=0,
    #     ub=1000,
    #     name="respiration"
    # )
    # model.add_cons_vars(some_flux)

    if new:
        mehler_like = model.problem.Constraint(
            model.reactions.PR0033.flux_expression - 0.1 * prod_sum,
            lb=0,
            ub=1000,
            name="mehlerlike"
        )
    elif latest:
        mehler_like = model.problem.Constraint(
            model.reactions.MEHLER_1.flux_expression - 0.1 * prod_sum,
            lb=0,
            ub=1000,
            name="mehlerlike"
        )
    model.add_cons_vars(mehler_like)

    ps1_superox = model.problem.Constraint(
        model.reactions.PR0032.flux_expression - 0.005 * prod_sum,
        lb=0,
        ub=1000,
        name="PSI_mehler"
    )
    ps2_superox = model.problem.Constraint(
        model.reactions.PR0034.flux_expression - 0.005 * prod_sum,
        lb=0,
        ub=1000,
        name="PSII_mehler"
    )
    model.add_cons_vars([ps1_superox, ps2_superox])

    # make all carbon influx 0, open CO2
    for rxn in model.boundary:
        if 'C' in rxn.check_mass_balance():
            rxn.lower_bound = 0.
    if new:
        model.reactions.EX_C00011_ext_b.lower_bound = -1000
    elif latest:
        model.reactions.EX_co2_e.lower_bound = -1000
    # make all nitrogen influx 0, open NO3
    for rxn in model.boundary:
        if 'N' in rxn.check_mass_balance():
            rxn.lower_bound = 0.
    if new:
        model.reactions.EX_C00244_ext_b.lower_bound = -1000
    elif latest:
        model.reactions.EX_no3_e.lower_bound = -1000

    return model


def split_rubisco(model: cobra.Model, ox_percent: float) -> cobra.Model:
    '''
    Split RuBisCO reactions and apply an oxygenation constraint.

    This function splits the RuBisCO (Ribulose-1,5-bisphosphate carboxylase/
    oxygenase) reaction in a metabolic model into two separate reactions: the
    carboxylase and the oxygenase. It then applies a constraint to maintain a
    specified ratio of oxygenation to carboxylation.

    Parameters:
    model (cobra.Model): A COBRApy metabolic model containing RuBisCO reaction.
    ox_percent (float): The desired percentage of oxygenation in RuBisCO
        activity (0-100).

    Returns:
    cobra.Model: The modified metabolic model with split RuBisCO reactions
        and the applied oxygenation constraint.

    Notes:
    - The function creates two new reactions, 'rubisCo' for carboxylation and
      'rubiscO' for oxygenation, and adds them to the model.
    - The 'ox_percent' parameter specifies the desired percentage of
      oxygenation in RuBisCO activity. It is used to create a constraint on the
      model to ensure the specified ratio is maintained.
    - The RuBisCO reactions are split into two based on the presence of
      specific metabolite IDs in the model.
    '''

    mid = [met.id for met in model.metabolites]
    if "C00011_cax" in mid:
        co2_id, o2_id, h2o_id, rb15bp_id, id_3pg, id_2pglyc = (
            "C00011_cax",
            "C00007_cax",
            "C00001_cax",
            "C01182_cax",
            "C00197_cax",
            "C00988_cax"
        )
    elif "co2_cx" in mid:
        co2_id, o2_id, h2o_id, rb15bp_id, id_3pg, id_2pglyc = (
            "co2_cx",
            "o2_cx",
            "h2o_cx",
            "rb15bp_cx",
            "3pg_cx",
            "2pglyc_cx"
        )
    m = model.metabolites
    rubisCo = cobra.Reaction("RBPC")
    rubisCo.name = "RuBisCO Carboxylase"
    rubisCo.lower_bound = 0
    rubisCo.add_metabolites({
        m.get_by_id(co2_id): -1,
        m.get_by_id(h2o_id): -1,
        m.get_by_id(rb15bp_id): -1,
        m.get_by_id(id_3pg): 2
    })
    rubiscO = cobra.Reaction("RBCh")
    rubiscO.name = "RuBisCO Oxygenase"
    rubiscO.lower_bound = 0
    rubiscO.add_metabolites({
        m.get_by_id(o2_id): -1,
        m.get_by_id(rb15bp_id): -1,
        m.get_by_id(id_3pg): 1,
        m.get_by_id(id_2pglyc): 1
    })
    if len(model.reactions.query("RBPC")) == 0:
        model.add_reactions([rubisCo, rubiscO])
        model.reactions.PP0011.upper_bound = 0
    else:
        rubisCo = model.reactions.RBPC
        rubiscO = model.reactions.RBCh

    rubisco_const = model.problem.Constraint(
        rubiscO.flux_expression - ox_percent/100 * rubisCo.flux_expression,
        lb=0,
        ub=1000,
        name="rubisco"
    )
    model.add_cons_vars([rubisco_const])

    return model
