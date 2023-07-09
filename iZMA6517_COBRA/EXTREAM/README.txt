**********************EXpression disTributed REAction flux Measurement (EXTREAM)********************
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
*************************Corresponding author: Rajib Saha (rsaha2@unl.edu)**************************

Following are the steps to run EXTREAM algorithm. This instruction is written considering only two 
conditions. Control and Examined.If there are more than two conditions, the principles will remain 
the same.

1. Arrange the control condition data as arranged in the data.xlsx file.
2. Run the expression_to_flux.py.
3. It will map transcriptomics to data reaction flux data. Discard all the reaction flux data that has
values 1000.
4. Repeat step 1-3 for the examined condition.
5. Get the highest reaction flux data between control and examined conditions and divide each reaction
flux value by that and then multiply that by a normalizing factor, which is usually 1000. For simplicity,
let's assume the value is "a".
6. If a reaction is forward then the bounds of reaction will be 0 to a. If the reaction is reversible
then the bounds will be -a to a.
7. Reactions, for which no reaction flux data found from the expression_to_flux.py, if the reaction is 
forward, set the bounds to 0 to 1000. If the reaction is reversible, set the bounds to -1000 to 1000.
8. Update the bounds of all reactions using COBRApy.
9. Reference flux (v_moics) vector will contain upper bounds of all reactions. The file can be formated as 
v_omics.xlsx.
10. Once step 1-9 are completed, run the flux balance analysis first to get the maximum biomass gorwth 
rate.
11. Run the EXTREAM.py by fixing the biomass growth rate!
****************************************************************************************************

