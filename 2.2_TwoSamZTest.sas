/* SAS Macro 2.2: Equivalence Trial with Binary Endpoint */
%Macro TwoSamZTest(nSims=100000, nPerGrp=100, px=0.3, py=0.4,
delta=0.3, alpha=0.05);
Data TwoGVars;
KEEP powerCI powerTest;
powerCI=0; powerTest=0;
Do iSim=1 To &nSims;
PxObs=Ranbin(733,&nPerGrp,&px)/&nPerGrp;
PyObs=Ranbin(236,&nPerGrp,&py)/&nPerGrp;
se=((PxObs*(1-PxObs)+PyObs*(1-PyObs))/&nPerGrp)**0.5;
*CI method;
ICW=Probit(1-&alpha)*se;
IF Abs(PxObs-PyObs)+ICW < &delta Then
powerCI=powerCI+1/&nSims;
*Two one-sided tests method;
T1=(PyObs-PxObs-&delta)/se;
T2=(PyObs-PxObs+&delta)/se;
IF T1=Probit(1-&alpha) & T2>Probit(1-&alpha) Then
powerTest=powerTest+1/&nSims;
End;
Output;
Run;
Proc Print; Run;
%Mend TwoSamZTest;

Title "Equivalence test with binary response: Alpha under Ho";
%TwoSamZTest(nPerGrp=100, px=0.1, py=0.2, delta=0.1, alpha=0.05);
Title "Equivalence test with binary response: Power under Ha";
%TwoSamZTest(nPerGrp=100, px=0.3, py=0.3, delta=0.2, alpha=0.05);
