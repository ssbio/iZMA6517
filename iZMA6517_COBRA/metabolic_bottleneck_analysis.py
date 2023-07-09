import cobra
from cobra.io import read_sbml_model
import pandas as pd
import openpyxl

#Create an excel file at the same directory as this python file
file_name = 'mba.xlsx'

#Load model
model = read_sbml_model('iZMA6517.xml')

#Formating output file
excel_file = openpyxl.load_workbook(file_name)
page = excel_file.active
result_df = pd.DataFrame(columns=['ID', 'Biomass', 'LB', 'UB'])

#Load optimization solver
model.solver = 'glpk'

#Print wild type biomass
print(f'Initial Biomass: {round(model.optimize().objective_value, 5)}')

#Checking each reactions to eliminate possible errors
print("Creating Initial Dataframe")
for reaction in model.reactions:
    print("**********************************************************")
    print(reaction)
    print(f'ID: {reaction.id}')
    print(f'Bound: LB: {model.reactions.get_by_id(reaction.id).lower_bound}, UB: {model.reactions.get_by_id(reaction.id).upper_bound}')
    #Save the initial bounds on a proxy varaible
    proxy_v_min = model.reactions.get_by_id(reaction.id).lower_bound
    proxy_v_max = model.reactions.get_by_id(reaction.id).upper_bound
    print("Modifying bounds of the reactions.....")
    #Expansion of reaction flux, toher bounds can be used as well, depending on the nature of metabolic netwrok
    model.reactions.get_by_id(reaction.id).lower_bound = -1000
    model.reactions.get_by_id(reaction.id).upper_bound = 1000
    print(f'Modified Bound: LB: {model.reactions.get_by_id(reaction.id).lower_bound}, UB: {model.reactions.get_by_id(reaction.id).upper_bound}')
    print(model.reactions.get_by_id(reaction.id))
    #Solve the flux balance analysis problem to calculate the modified biomass growth rate after debottlenecking
    modified_fba = model.optimize()
    print(f'Modified Biomass: {round(modified_fba.objective_value, 5)}')
    print("Reverting bounds of the reactions.....")
    #Reverting the reaction flux to the original state
    model.reactions.get_by_id(reaction.id).lower_bound = proxy_v_min
    model.reactions.get_by_id(reaction.id).upper_bound = proxy_v_max
    print(f'Reverted Bound: LB: {model.reactions.get_by_id(reaction.id).lower_bound}, UB: {model.reactions.get_by_id(reaction.id).upper_bound}')
    page.append([reaction.id, modified_fba.objective_value])
    excel_file.save(filename=file_name)
    # df.loc[len(df.index)] = [reaction.id, modified_fba.objective_value, model.reactions.get_by_id(reaction.id).lower_bound, model.reactions.get_by_id(reaction.id).upper_bound]
    print("**********************************************************")