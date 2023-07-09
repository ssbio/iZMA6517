import pandas as pd
import cobra
from cobra.io import read_sbml_model
import numpy as np
import openpyxl
from cobra.flux_analysis import flux_variability_analysis, pfba


if __name__ == '__main__':
    # Load the model
    model = read_sbml_model('MODEL.xml')
    model.solver = 'gurobi'

    #turn off certain reactions if required
    TurnOffRxn = ['R0XXXX[Z,l]']
    
    for rxn in TurnOffRxn:
        model.reactions.get_by_id(rxn).upper_bound = 0
        model.reactions.get_by_id(rxn).lower_bound = 0

    #fix bounds of certain certain reactions if required
    AARxn = ['EXE_XXX']

    for rxn in AARxn:
        model.reactions.get_by_id(rxn).upper_bound = 0
        model.reactions.get_by_id(rxn).lower_bound = -2.75

    # Perform FBA to find the maximum reaction flux of biomass equation
    model.objective = 'Plant_Biomass'
    max_biomass_flux = model.optimize().objective_value

    # Fix the biomass production rate at its maximum value
    model.reactions.get_by_id('Plant_Biomass').lower_bound = max_biomass_flux
    model.reactions.get_by_id('Plant_Biomass').upper_bound = max_biomass_flux

    # Load v_omics data from Excel file
    v_omics = pd.read_excel('v_omics.xlsx', index_col=0, squeeze=True)

    # Add dummy variable and reformulation constraints
    for reaction in model.reactions:
        if reaction.id in v_omics.index:
            d = model.problem.Variable('D_' + reaction.id, lb=0, type='continuous')
            model.add_cons_vars(d)

            reformulation1 = model.problem.Constraint(
                reaction.flux_expression - d - v_omics[reaction.id],
                ub=0,
                name='reformulation1_' + reaction.id
            )
            model.add_cons_vars(reformulation1)

            reformulation2 = model.problem.Constraint(
                - reaction.flux_expression - d + v_omics[reaction.id],
                ub=0,
                name='reformulation2_' + reaction.id
            )
            model.add_cons_vars(reformulation2)

    # Add objective based on omics data
    objective = model.problem.Objective(
        sum(d for d in model.variables if d.name.startswith('D_')),
        direction='min'
    )
    model.objective = objective

    # Solve the model
    solution = model.optimize()

    # Write the results to a file
    with open('extream_result.txt', 'w') as f:
        f.write("reaction      FLUX\n")
        for reaction in model.reactions:
            f.write(f"{reaction.id}    {solution.fluxes[reaction.id]:.8f}\n")
