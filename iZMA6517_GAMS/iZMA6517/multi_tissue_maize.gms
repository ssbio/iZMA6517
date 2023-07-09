*************************************************************
*****************pFBA for whole maize plant******************
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
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective					objective function
	mass_balance(i)				steady state mass balance
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
	biomass						fixing biomass
	pFBA1_reformulation(j)		pFBA reformulation constraint 1
	pFBA2_reformulation(j)		pFBA reformulation constraint 2
	photorespiration			setting up the ratio of photorespiration
	leafmalateshunt				shutting down leaf malate shunt
	shoot_CO2_out				shoot can survive without photosynthesis
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	Z					objective value
	dummy(j)				dummy variable for pFBA
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= sum(j, dummy(j));

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

biomass..			v('Root_Biomass[R]') + v('Shoot_Biomass[S]') + v('Seed_Biomass[K]') + v('Leaf_Biomass[L]') =e= 52.76;

pFBA1_reformulation(j)..	v(j) =l= dummy(j);

pFBA2_reformulation(j)..	-v(j) =l= dummy(j);

photorespiration..		v('R03140[B,p]') =e= 0.0003 * v('R00024[B,p]');

leafmalateshunt..		v('R00342[M,c]') =e= 0;

shoot_CO2_out..			v('OUT_C00011[S,c]') =g= 0;

Model multi_tissue_maize /all/;

******************************************************************

**********Solving Model*********************

solve multi_tissue_maize using lp maximizing Z;

********************************************
****************Output File*****************
FILE RESULTS /multi_tissue_maize_flux.txt/;

PUT RESULTS;

PUT "reaction      FLUX"/;

LOOP(j,
	
	PUT j.tl:0:100,"    ", v.l(j):20:7/;
		
);

PUTCLOSE;
**********************************************
*$(S('C00238[P]',j) < 0 and v.l(j) < 0)