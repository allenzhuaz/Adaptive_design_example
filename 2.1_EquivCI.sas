/* SAS Macro 2.1: Equivalence Trial with Normal Endpoint */
%Macro EquivCI(nSims=1000, nPerGrp=200, ux=0, uy=1, delta=1.2,
sigmax=1, sigmay=1.2, alpha=0.05);
Data TwoGVars;
Keep xMean yMean powerCI powerTest;
powerCI=0; powerTest=0;
Do iSim=1 To &nSims;
xMean=0; yMean=0; s2x=0; s2y=0;
Do iObs=1 To &nPerGrp;
xNOR=Rannor(7362); xMean=xMean+xNor; s2x=s2x+xNor**2;
yNOR=Rannor(2637); yMean=yMean+yNor; s2y=s2y+yNor**2;
End;
xMean=xMean*&sigmax/&nPerGrp+&ux;
yMean=yMean*&sigmay/&nPerGrp+&uy;
sp=((s2x*&sigmax**2+s2y*&sigmay**2)/(2*&nPerGrp-2))**0.5;
se=sp/(&nPerGrp/2)**0.5;
* CI method;
ICW=Probit(1-&alpha)*se;
If Abs(yMean-xMean)+ICW < &delta Then
powerCI=powerCI+1/&nSims;
*Two one-sided test method;
T1=(xMean-yMean-&delta)/se;
T2=(xMean-yMean+&delta)/se;
If T1=Probit(1-&alpha) & T2>Probit(1-&alpha) Then
powerTest=powerTest+1/&nSims;
End;
Output;
Run;
Proc Print Data=TwoGVars(obs=1); Run;
%Mend EquivCI;

/* Invoke the macro */
Title "Equivalence test with normal response: Alpha under Ho";
%EquivCI(nSims=10000, nPerGrp=1000, ux=0.2, uy=0, delta=0.2, sigmax=1,
sigmay=1, alpha=0.05);
Title "Equivalence test with normal response: Power under Ha";
%EquivCI(nSims=10000, nPerGrp=198, ux=0, uy=1, delta=1.2, sigmax=0.8,
sigmay=0.8, alpha=0.05);
