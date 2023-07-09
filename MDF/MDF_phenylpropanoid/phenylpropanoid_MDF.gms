*************************************************************
********MAX/MIN DRIVING FORCE ANALYSIS FOR CHLAMYDIA*********
*******************PHENYLPROPANOID PATHWAY*******************
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
/278.15/
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	constraint_1(j)				constranit 1 of the formulation
	constraint_2(j)				constranit 2 of the formulation
	upper_bound(i)				upper bound
	lower_bound(i)				lower bound
	NADHtoNAD				Fixing the ratio of NADPH to NADP+ according to Noor et al. 2014
	NADPHtoNADP				Fixing the ratio of NADPH to NADP+ according to Noor et al. 2014
	Other1
	Other2
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

NADHtoNAD..			x('C00004') =e= 0.1 * x('C00003');

NADPHtoNADP..			x('C00005') =e= 10 * x('C00006');

Other1..			deltaG('R01945') =l= deltaG('R05700');

Other2..			deltaG('R07432') =l= deltaG('R05700');

Model phenylpropanoid_MDF /all/;
*****************************************************************
****************Output File*****************

x.up('C00003') = 0.0;
x.lo('C00003') = 0.0;

FILE RESULTS1 /NAD_DF.txt/;

PUT RESULTS1;

RESULTS1.pw=1000;


PUT "NAD+	MDF	R01616	R07432	R04342	R01945	R01942	R02193	R05700"/;

LOOP(SMP,

	x.up('C00003') = x.up('C00003') + 0.1;

	x.lo('C00003') = x.lo('C00003') + 0.1;

	PUT x.up('C00003'):0:1;

	solve phenylpropanoid_MDF using lp maximizing Z;

	IF((phenylpropanoid_MDF.ModelStat eq 4),

		ABORT "ABORTED DUE TO INFEASIBLE MODEL STATUS";
	);
	
	PUT Z.l:20:1;

	PUT deltaG.l('R01616'):20:1;

	PUT deltaG.l('R07432'):20:1;

	PUT deltaG.l('R04342'):20:1;

	PUT deltaG.l('R01945'):20:1;

	PUT deltaG.l('R01942'):20:1

	PUT deltaG.l('R02193'):20:1

	PUT deltaG.l('R05700'):20:1/;
			
);

PUTCLOSE;
**********************************************