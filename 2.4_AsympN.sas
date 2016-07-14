/* SAS Macro 2.4: Sample Size for Dose-Response Trial */
%Macro AsympN(endpoint="normal", alpha=0.025, power=0.8, nArms=5,
delta=0, tStd=12, tAcr=4);
Data contrastN; Set dInput;
Keep Endpoint nArms alpha power TotalSampleSize;
Array u{&nArms}; Array s{&nArms}; Array f{&nArms};
Array c{&nArms}; endpoint=&endpoint; delta=&delta;
alpha=&alpha; power=&power; nArms=&nArms;
epi = 0; s0 = 0;
Do i =1 To nArms; epi = epi + cfig*ufig- &delta/nArms; End;
If &endpoint = "normal" Then Do;
Do i =1 To nArms; s0 = s0 + s{i}/nArms; End;
End;
If &endpoint = "binary" Then Do;
Do i = 1 To nArms;
sfig = (u{i}*(1-u{i}))**0.5;
s0=s0 + s{i}/nArms;
End;
End;
If &endpoint = "survival" Then Do;
Do i = 1 To nArms;
sfig = ufig*(1+exp(-u{i}*&tStd)*(1-exp(u{i}*&tAcr))/(&tAcr*u{i}))**(-0.5);
s0 = s0 + s{i}/nArms;
End;
End;
sumscf0 = 0; sumscf = 0;
Do i = 1 To nArms; sumscf0 = sumscf0 + s0**2*c{i}*c{i}/f{i}; End;
Do i = 1 To nArms; sumscf = sumscf + s{i}**2*c{i}*c{i}/f{i}; End;
n = ((PROBit(1-&alpha)*sumscf0**0.5 + Probit(&power)*sumscf**0.5)/
epi)**2;
TotalSampleSize = round(n);
Run;
Proc print;
Run;
%Mend AsympN;

/* Invoke the macro */
Title " = s of How to Use the SAS Macros";
Data dInput;
Array u{4}(.46, .35, .32, .3); ** Responses;
Array s{4}(2, 2, 2, 2); ** Standard deviation for normal endpoint;
Array c{4}(-4, 1, 1, 2); ** Contrasts;
Array f{4}(.25, .25 ,.25, .25); ** Sample size fractions;
%AsympN(endpoint="normal", alpha=0.025, power=0.8, nArms=4);
%AsympN(endpoint="binary", alpha=0.025, power=0.8, nArms=4);
%AsympN(endpoint="survival", alpha=0.025, power=0.8, nArms=4, delta=0,
tStd=2, tAcr=.5);
Run;
