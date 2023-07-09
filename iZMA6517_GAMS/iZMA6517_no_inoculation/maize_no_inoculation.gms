*************************************************************
***********Maize genome-scale metabolic model****************
*************************************************************
*****************Niaz Bahar Chowdhury************************
*************************************************************
$INLINECOM /*  */

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"
;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"

	v_omics(j)				referecen state of v(j) calculated from the GPR relation

$include "v_omics.txt"
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance(i)				steady state mass balance
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
	biomass					fixing biomass production rate
	reformulation1(j)			reformulation constraint 1
	reformulation2(j)			reformulation constraint 2
	photorespiration			setting maximum rate of photorespiration
	leafmalateshunt				shutting down leaf malate shunt
	shoot_CO2_out				shoot can survive without photosynthesis
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	Z					objective value
	D(j)					dummy variable for refomulation
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= sum(j, (D(j) - v_omics(j)));

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

biomass..			v('Root_Biomass[R]') + v('Shoot_Biomass[S]') + v('Seed_Biomass[K]') + v('Leaf_Biomass[L]') =e= 0.754;

reformulation1(j)..		v(j) =l= D(j) + v_omics(j);

reformulation2(j)..	      - v(j) =l= D(j) - v_omics(j);

photorespiration..		v('R03140[B,p]') =e= 0.0003 * v('R00024[B,p]');

leafmalateshunt..		v('R00342[M,c]') =e= 0;

shoot_CO2_out..			v('OUT_C00011[S,c]') =g= 0;

Model maize /all/;
******************************************************************

**********Solving Model*********************

maize.holdfixed = 1;

solve maize using lp maximizing Z;
********************************************
****************Output File*****************
FILE RESULTS /maize_result.txt/;

PUT RESULTS;

PUT "reaction      FLUX"/;

LOOP(j,
	
	PUT j.tl:0:100,"    ", v.l(j):20:8/;
		
);

PUTCLOSE;
**********************************************