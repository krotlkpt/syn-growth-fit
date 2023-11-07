from modules.import_model import import_model
import cobra


def add_BM10(model: cobra.Model) -> cobra.Model:
    '''
    Modify biomass reaction to include glycogen as a component.

    This function modifies the biomass reaction in a metabolic model to include
    glycogen as a component of the biomass. Glycogen is scaled as a mass
    percentage to ensure that the total biomass components sum to 1 g.

    Parameters:
    model (cobra.Model): A COBRApy metabolic model with an existing biomass
        reaction.

    Returns:
    cobra.Model: The modified metabolic model with the updated biomass
        reaction.

    Notes:
    - The function creates a new metabolite 'B10_cyt' (or 'bm_glycogen_c' if
      using BiGG IDs) representing the glycogen component in the cytoplasm.
    - It also adds a new reaction 'BM0010' representing the glycogen component
      of the biomass reaction.
    - The glycogen mass is scaled as a mass percentage, ensuring that all
      biomass components sum to 1 g.
    - The function modifies the existing biomass reaction ('BM0009') to include
      the 'B10_cyt' metabolite as a component and adjusts the coefficients
      to balance the masses.

    Example:
    - If you have a biomass reaction with other components that sum to 1 g,
      this function allows you to include glycogen as a component,
      and it scales its mass based on the specified mass percentage.

    '''

    try:
        glycogen = model.metabolites.C00182_cyt
        bigg = False
    except AttributeError:
        glycogen = model.metabolites.glycogen_c
        bigg = True
    if not bigg:
        bm_glyc = cobra.Metabolite("B10_cyt", "R")
    else:
        bm_glyc = cobra.Metabolite("bm_glycogen_c", "R")
    bm_glyc.name = 'Glycogen component of biomass_cyt'
    bm_glyc.compartment = "cyt"
    model.add_metabolites([bm_glyc])

    bm10 = cobra.Reaction('BM0010')
    bm10.name = 'Glycogen component of biomass reaction'
    bm10.subsystem = 'biomass'
    bm10.lower_bound = 0
    bm10.upper_bound = 1000

    if not bigg:
        bm10.add_metabolites({
            model.metabolites.get_by_id('C00182_cyt'): -0.210311/.0341,
            model.metabolites.get_by_id('B10_cyt'): 1.0,
        })
    else:
        bm10.add_metabolites({
            model.metabolites.get_by_id('glycogen_c'): -0.210311/.0341,
            model.metabolites.get_by_id('bm_glycogen_c'): 1.0,
        })

    model.add_reactions([bm10])

    # glycogen = model.metabolites.C00182_cyt
    bof = model.reactions.BM0009
    bof.subtract_metabolites({
        glycogen: -0.21031
    })
    bof.add_metabolites({
        bm_glyc: -0.0341
    })

    return model


def update_glycogen(model: cobra.Model, percent: float = 3.41) -> cobra.Model:
    '''
    Update biomass reaction to adjust glycogen percentage.

    This function updates the biomass reaction in a metabolic model to adjust
    the percentage of glycogen in the biomass. All other biomass components
    are scaled accordingly to maintain a total of 1 g biomass.

    Parameters:
    model (cobra.Model): A COBRApy metabolic model with an existing biomass
        reaction.
    percent (float, optional): The desired percentage of glycogen in the
        biomass (default is 3.41%).

    Returns:
    cobra.Model: The modified metabolic model with updated biomass reaction.

    Notes:
    - The function checks if the 'BM0010' reaction representing glycogen is
      already in the model. If not, it calls the 'add_BM10' function to add
      it.
    - It adjusts the coefficients of all other biomass components to ensure
      that the total biomass components still sum to 1 g.
    - The 'percent' parameter specifies the desired percentage of glycogen
      in the biomass, and the function calculates the necessary adjustments
      to other components.

    Example:
    - If you want to change the glycogen percentage in the biomass, you can use
      this function with the desired percentage, and it will adjust the biomass
      composition accordingly.

    '''

    try:
        model.reactions.get_by_id("BM0010")
    except KeyError:
        model = add_BM10(model)
    bof = model.reactions.BM0009
    glyc_component = model.metabolites.get_by_id("B10_cyt")
    dif = (percent/100 + bof.metabolites[glyc_component]) * -1
    bof.subtract_metabolites({
        met: dif*bof.metabolites[met]/sum(
            [
                bof.get_coefficient(met)
                for met in bof.reactants
                if met.id not in [
                    "C00002_cyt", "C00001_cyt", "B10_cyt"
                ]
            ]
        )
        for met in bof.reactants
        if met not in [
            model.metabolites.get_by_id(metid)
            for metid in [
                "C00002_cyt", "C00001_cyt", "B10_cyt"
            ]
        ]
    })
    bof.add_metabolites({
        glyc_component: dif
    })

    return model


def update_prot_glyc(
    model: cobra.Model,
    prot_percent: float = 51, glyc_percent: float = 3.41
) -> cobra.Model:
    '''
    Update biomass reaction to adjust protein and glycogen percentages.

    This function updates the biomass reaction in a metabolic model to adjust
    both the percentage of protein and glycogen in the biomass. All other
    biomass components are scaled accordingly to maintain a total of
    1 g biomass.

    Parameters:
    model (cobra.Model): A COBRApy metabolic model with an existing biomass
        reaction.
    prot_percent (float, optional): The desired percentage of protein in the
        biomass (default is 51%).
    glyc_percent (float, optional): The desired percentage of glycogen in the
        biomass (default is 3.41%).

    Returns:
    cobra.Model: The modified metabolic model with updated biomass reaction.

    Notes:
    - The function checks if the 'BM0010' reaction representing glycogen is
      already in the model. If not, it calls the 'add_BM10' function to add it.
    - It adjusts the coefficients of protein and glycogen as well as other
      biomass components to ensure that the total biomass components still sum
      to 1 g.
    - The 'prot_percent' and 'glyc_percent' parameters specify the desired
      percentages of protein and glycogen in the biomass, and the function
      calculates the necessary adjustments to other components.
    - The function also considers whether BiGG IDs are used for metabolites and
      adjusts accordingly.

    Example:
    - If you want to change both the protein and glycogen percentages in the
      biomass, you can use this function with the desired percentages, and it
      will adjust the biomass composition accordingly.

    '''

    try:
        model.reactions.get_by_id("BM0010")
    except KeyError:
        model = add_BM10(model)
    try:
        _ = model.metabolites.C00182_cyt
        bigg = False
    except AttributeError:
        _ = model.metabolites.glycogen_c
        bigg = True
    bof = model.reactions.BM0009
    if not bigg:
        glyc_component = model.metabolites.get_by_id("B10_cyt")
        prot_component = model.metabolites.get_by_id("B1_cyt")
        extra_mets = [
            "C00002_cyt", "C00001_cyt", "B10_cyt", "B1_cyt"
        ]
    else:
        glyc_component = model.metabolites.get_by_id("bm_glycogen_c")
        prot_component = model.metabolites.get_by_id("bm_pro_c")
        extra_mets = [
            "atp_c", "h2o_c", "bm_glycogen_c", "bm_pro_c"
        ]
    dif = sum([
        glyc_dif := (glyc_percent/100 + bof.metabolites[glyc_component]) * -1,
        prot_dif := (prot_percent/100 + bof.metabolites[prot_component]) * -1
    ])
    bof.subtract_metabolites({
        met: dif*bof.metabolites[met]/sum(
            [
                bof.get_coefficient(met)
                for met in bof.reactants
                if met.id not in extra_mets
            ]
        )
        for met in bof.reactants
        if met not in [
            model.metabolites.get_by_id(metid)
            for metid in extra_mets
        ]
    })
    bof.add_metabolites({
        glyc_component: glyc_dif,
        prot_component: prot_dif
    })

    return model


if __name__ == "__main__":
    new_syn = import_model(
            '../models_sbml/Synechocystis6803_SBML_COBRA_24_4_17.sbml'
        )
    with new_syn:
        new_syn = add_BM10(new_syn)
        print(new_syn.reactions.BM0009.build_reaction_string(
            use_metabolite_names=True
        ))
        print(sum(new_syn.reactions.BM0009.metabolites.values()))
    with new_syn:
        new_syn = update_prot_glyc(new_syn, 20, 20)
        print(new_syn.reactions.BM0009.build_reaction_string(
            use_metabolite_names=True
        ))
        print(sum(new_syn.reactions.BM0009.metabolites.values()))
    with new_syn:
        new_syn = update_glycogen(new_syn, 98)
        print(new_syn.reactions.BM0009.build_reaction_string(
            use_metabolite_names=True
        ))
        print(sum(new_syn.reactions.BM0009.metabolites.values()))
