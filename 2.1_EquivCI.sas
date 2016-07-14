/* SAS Macro 2.1: Equivalence Trial with Normal Endpoint */
%macro equivci(nsims=1000, npergrp=200, ux=0, uy=1, delta=1.2,
				sigmax=1, sigmay=1.2, alpha=0.05);
data twogvars;
	keep xmean ymean powerci powertest;
powerci=0; powertest=0;
do isim=1 to &nsims;
	xmean=0; ymean=0; s2x=0; s2y=0;
	do iobs=1 to &npergrp;
		xnor=rannor(7362); xmean=xmean+xnor; s2x=s2x+xnor**2;
		ynor=rannor(2637); ymean=ymean+ynor; s2y=s2y+ynor**2;
	end;
	xmean=xmean*&sigmax/&npergrp+&ux;
	ymean=ymean*&sigmay/&npergrp+&uy;
	sp=((s2x*&sigmax**2+s2y*&sigmay**2)/(2*&npergrp-2))**0.5;
	se=sp/(&npergrp/2)**0.5;
	* ci method;
	icw=probit(1-&alpha)*se;
	if abs(ymean-xmean)+icw < &delta then
	powerci=powerci+1/&nsims;
	*two one-sided test method;
	t1=(xmean-ymean-&delta)/se;
	t2=(xmean-ymean+&delta)/se;
	if t1=probit(1-&alpha) & t2>probit(1-&alpha) then
	powertest=powertest+1/&nsims;
end;
	output;
run;
proc print data=twogvars(obs=1); run;
%mend equivci;

/* invoke the macro */
title "equivalence test with normal response: alpha under ho";
%equivci(nsims=10000, npergrp=1000, ux=0.2, uy=0, delta=0.2, sigmax=1,
sigmay=1, alpha=0.05);
title "equivalence test with normal response: power under ha";
%equivci(nsims=10000, npergrp=198, ux=0, uy=1, delta=1.2, sigmax=0.8,
sigmay=0.8, alpha=0.05);
