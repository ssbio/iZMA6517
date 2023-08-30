**********************EXpression disTributed REAction flux Measurement (EXTREAM)********************
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
*************************Corresponding author: Rajib Saha (rsaha2@unl.edu)**************************

Following are the steps to run EXTREAM algorithm. This instruction is written considering only two 
conditions. Control and Examined.If there are more than two conditions, the principles will remain 
the same.

1. Get the number of unique genes and number of each genes from the GPR and divide each expression by the
gene number.
2. Arrange the control condition data as arranged in the expression.xlsx file.
3. Run the expression_to_flux.py.
4. It will map transcriptomics to data reaction flux data. Discard all the reaction flux data that has
values 1000.
5. Repeat step 1-3 for the examined condition.
6. Get the highest reaction flux data between control and examined conditions and divide each reaction
flux value by that and then multiply that by a normalizing factor, which is usually 1000. For simplicity,
let's assume the value is "a".
7. If a reaction is forward then the bounds of reaction will be 0 to a. If the reaction is reversible
then the bounds will be -a to a.
8. Reactions, for which no reaction flux data found from the expression_to_flux.py, if the reaction is 
forward, set the bounds to 0 to 1000. If the reaction is reversible, set the bounds to -1000 to 1000.
9. Update the bounds of all reactions using COBRApy.
10. Reference flux (v_moics) vector will contain upper bounds of all reactions. The file can be formated as 
v_omics.xlsx.
11. Once step 1-9 are completed, run the flux balance analysis first to get the maximum biomass gorwth 
rate.
12. Run the EXTREAM.py by fixing the biomass growth rate!
****************************************************************************************************

