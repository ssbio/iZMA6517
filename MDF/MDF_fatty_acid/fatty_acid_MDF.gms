*************************************************************
********MAX/MIN DRIVING FORCE ANALYSIS FOR CHLAMYDIA*********
*********************FATTY ACID PATHWAY**********************
********************Niaz Bahar Chowdhury*********************
*************************************************************

$INLINECOM /*  */
$ONEMPTY

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"

	SMP					set for concnetration sampling	/1*100/
;

*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	delta_G_o(j)				maximum flux of v(j)

$include "deltaGo.txt"

	Cmax(i)					highest concentration

$include "cmax.txt"

	Cmin(i)					Lowest concentration

$include "cmin.txt"

	R					Gas constant
/0.008314/

	T					Temperature
/323.15/
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	constraint_1(j)				constranit 1 of the formulation
	constraint_2(j)				constranit 2 of the formulation
	upper_bound(i)				upper bound
	lower_bound(i)				lower bound
*	fix_objective				forcing objective function to be positive
;
**************************************************************

*********Defining Variables***********************************
POSITIVE VARIABLES

	x(i)					concentration of a metabolite

FREE VARIABLES

	v(j)					reaction flux
	B					Maximum driving force
	deltaG(j)				driving force in the given condition
	Z					objective value
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= B;

constraint_1(j)..		- deltaG(j) =g= B;

constraint_2(j)..		deltaG(j) =e= delta_G_o(j) + R * T * sum(i, S(i,j) * x(i));

upper_bound(i)..		x(i) =l= log(Cmax(i));

lower_bound(i)..		x(i) =g= log(Cmin(i));

*fix_objective..			deltaG('R08163') =g= deltaG('R01626');

Model fatty_acid_MDF /all/;
*****************************************************************
****************Output File*****************

x.up('C00003') = 0.0;
x.lo('C00003') = 0.0;

FILE RESULTS1 /C00003_DF.txt/;

PUT RESULTS1;

RESULTS1.pw=1000;


PUT "C00003	MDF	R07758	R07759	R07760	R07761	R02222	R00100"/;

LOOP(SMP,

	x.up('C00003') = x.up('C00003') + 0.1;

	x.lo('C00003') = x.lo('C00003') + 0.1;

	PUT x.up('C00003'):0:1;

	solve fatty_acid_MDF using lp maximizing Z;

	IF((fatty_acid_MDF.ModelStat eq 4),

		ABORT "ABORTED DUE TO INFEASIBLE MODEL STATUS";
	);
	
	PUT Z.l:20:1;

	PUT deltaG.l('R07758'):20:1;

	PUT deltaG.l('R07759'):20:1;

	PUT deltaG.l('R07760'):20:1;

	PUT deltaG.l('R07761'):20:1;

	PUT deltaG.l('R02222'):20:1;

	PUT deltaG.l('R00100'):20:1/;
			
);

PUTCLOSE;
**********************************************