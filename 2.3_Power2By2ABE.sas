/* SAS Macro 2.3: Crossover Bioequivalence Trial */
%macro power2by2abe(totaln=24, swithin=0.355, uratio=1);
data abe; 
	keep swithin uratio n power;
n=&totaln; swithin=&swithin; uratio=&uratio;
* Err df for AB/BA crossover design;
n2=n-2;
t1=tinv(1-0.05,n-2); t2=-t1;
nc1=sqrt(n)*log(uratio/0.8)/sqrt(2)/swithin;
nc2=sqrt(n)*log(uratio/1.25)/sqrt(2)/swithin;
df=sqrt(n-2)*(nc1-nc2)/(2*t1);
power=probt(t2,df,nc2)-probt(t1,df,nc1);
run;
proc print; run;
%mend power2by2abe;

/* Invoke the macro */
%power2by2abe(totaln=58, swithin=0.355, uratio=1)

