*************************************************************
********MAX/MIN DRIVING FORCE ANALYSIS FOR CHLAMYDIA*********
*****************C4 PHOTOSYNTHESIS PATHWAY*******************
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
/313.15/
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	constraint_1(j)				constranit 1 of the formulation
	constraint_2(j)				constranit 2 of the formulation
	upper_bound(i)				upper bound
	lower_bound(i)				lower bound
	ATPtoADP				ATP to ADP ratio fixed according to Noor et al. 2014
	NADHtoNAD				NADH to NAD ratio fixed according to Noor et al. 2014
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

ATPtoADP..			x('C00002') =e= 10 * x('C00008');

NADHtoNAD..			x('C00004') =e= 0.1 * x('C00003');

Model photosynthesis_MDF /all/;
*****************************************************************
****************Output File*****************

x.up('C00003') = 0.0;
x.lo('C00003') = 0.0;

FILE RESULTS1 /NAD_DF.txt/;

PUT RESULTS1;

RESULTS1.pw=1000;


PUT "NAD+	MDF	R00345	R00342	R00216	R00206	R00024	R01512	R01061	R01068	R00762	R01067	R01829	R01845	R01641	R01056	R01523"/;

LOOP(SMP,

	x.up('C00003') = x.up('C00003') + 0.1;

	x.lo('C00003') = x.lo('C00003') + 0.1;

	PUT x.up('C00003'):0:1;

	solve photosynthesis_MDF using lp maximizing Z;

	IF((photosynthesis_MDF.ModelStat eq 4),

		ABORT "ABORTED DUE TO INFEASIBLE MODEL STATUS";
	);
	
	PUT Z.l:20:1;

	PUT deltaG.l('R00345'):20:1;

	PUT deltaG.l('R00342'):20:1;

	PUT deltaG.l('R00216'):20:1;

	PUT deltaG.l('R00206'):20:1;

	PUT deltaG.l('R00024'):20:1;

	PUT deltaG.l('R01512'):20:1;

	PUT deltaG.l('R01061'):20:1;

	PUT deltaG.l('R01068'):20:1;

	PUT deltaG.l('R00762'):20:1;

	PUT deltaG.l('R01067'):20:1;

	PUT deltaG.l('R01829'):20:1;

	PUT deltaG.l('R01845'):20:1;

	PUT deltaG.l('R01641'):20:1;

	PUT deltaG.l('R01056'):20:1;

	PUT deltaG.l('R01523'):20:1/;
			
);

PUTCLOSE;
**********************************************