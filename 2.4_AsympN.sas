/* SAS macro 2.4: sample size for dose-response trial */
%macro asympn(endpoint="normal", alpha=0.025, power=0.8, narms=5,
				delta=0, tstd=12, tacr=4);
data contrastn; set dinput;
	keep endpoint narms alpha power totalsamplesize;
array u{&narms}; array s{&narms}; array f{&narms}; array c{&narms}; 
endpoint=&endpoint; delta=&delta; alpha=&alpha; power=&power; narms=&narms;
epi = 0; s0 = 0;
do i =1 to narms; 
	epi = epi + cfig*ufig- &delta/narms; 
end;
if &endpoint = "normal" then do;
	do i =1 to narms; 
		s0 = s0 + s{i}/narms; 
	end;
end;
if &endpoint = "binary" then do;
	do i = 1 to narms;
		sfig = (u{i}*(1-u{i}))**0.5;
		s0=s0 + s{i}/narms;
	end;
end;
if &endpoint = "survival" then do;
	do i = 1 to narms;
		sfig = ufig*(1+exp(-u{i}*&tstd)*(1-exp(u{i}*&tacr))/(&tacr*u{i}))**(-0.5);
		s0 = s0 + s{i}/narms;
	end;
end;
sumscf0 = 0; sumscf = 0;
do i = 1 to narms; sumscf0 = sumscf0 + s0**2*c{i}*c{i}/f{i}; end;
do i = 1 to narms; sumscf = sumscf + s{i}**2*c{i}*c{i}/f{i}; end;
n = ((probit(1-&alpha)*sumscf0**0.5 + probit(&power)*sumscf**0.5)/
epi)**2;
totalsamplesize = round(n);
run;
proc print;
run;
%mend asympn;

/* invoke the macro */
title " = s of how to use the sas macros";
data dinput;
array u{4}(.46, .35, .32, .3); ** responses;
array s{4}(2, 2, 2, 2); ** standard deviation for normal endpoint;
array c{4}(-4, 1, 1, 2); ** contrasts;
array f{4}(.25, .25 ,.25, .25); ** sample size fractions;
%asympn(endpoint="normal", alpha=0.025, power=0.8, narms=4);
%asympn(endpoint="binary", alpha=0.025, power=0.8, narms=4);
%asympn(endpoint="survival", alpha=0.025, power=0.8, narms=4, delta=0,
tstd=2, tacr=.5);
run;
