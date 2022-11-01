(* ::Package:: *)

BeginPackage["TeukolskyT`"]

FieldEnergyMomentum::usage = "FieldEnergyMomentum[A_\mu,Procamass,spin] :> Output Tdowndown"
Tdowndown::usage = "Tdowndown[t,r,\[Theta],\[Phi]] :> Gives energy momentum tensor of A_\mu that was put into FieldEnergyMomentum."
TeukolskyTintegrandExpl
TeukolskyTlmw
TeukolskyTintegrand
TeukolskyTintegrandExplUno

TnnExpl
tnnexpltemp
TmnExpl
tmnexpltemp
TmmExpl
tmmexpltemp

Tddtemp
AdownInT
Aexpl
dAexpl

tnntempR
tnntempI

TnnInter
TmmInter
tnntempR
TmnInter
Tlm\[Omega]
Tlm\[Omega]temp
tnnDue
tnnUno
tnntempRin
tnntempRout
Tlm\[Omega]temp
TijtoExpl
dTijtoExpl
ddTijtoExpl
TnnInterNumRe

Sbar
sinfunc
rrangeout

plotsbar

CubicTerm
CubicTermOutput

Begin["PrivateTeukolsky`"]
Needs["xAct`xCoba`"]
(*Needs["SpinWeightedSpheroidalHarmonics`","./SWSHs.m"]*)
Needs["SpinWeightedSpheroidalHarmonics`","./SpinWeightedSpheroidalHarmonics.m"]

(*xAct initializations*)
$DefInfoQ=False;
$CVVerbose=False;
$PrePrint=ScreenDollarIndices;
$ShowTimeThreshold=0.5
DefManifold[M,4,IndexRange[b,l]];
DefMetric[-1,met[-b,-c],CD,PrintAs->"g"];
DefChart[ch,M,{0,1,2,3},{t[],r[],\[Theta][],\[Phi][]}];
ch/:CIndexForm[0,ch]:="t";
ch/:CIndexForm[1,ch]:="r";
ch/:CIndexForm[2,ch]:="\[Theta]";
ch/:CIndexForm[3,ch]:="\[Phi]";
DefConstantSymbol[\[Chi]]

(*Kerr metric*)
\[CapitalSigma]=r[]^2+\[Chi]^2 Cos[\[Theta][]]^2;
\[CapitalDelta]=r[]^2-2 r[]+\[Chi]^2;
KerrBL={{-(1-(2r[])/\[CapitalSigma]),0,0,(-2 r[] \[Chi] Sin[\[Theta][]]^2)/\[CapitalSigma]},{0,\[CapitalSigma]/\[CapitalDelta],0,0},{0,0,\[CapitalSigma],0},{(-2 r[] \[Chi] Sin[\[Theta][]]^2)/\[CapitalSigma],0,0,(r[]^2+\[Chi]^2+2 r[]/\[CapitalSigma] \[Chi]^2 Sin[\[Theta][]]^2)Sin[\[Theta][]]^2}};
ComponentValue[met[-b,-c]//ToBasis[ch]//ComponentArray,KerrBL]
MetricCompute[met,ch,All]
MetricCompute[met,ch,"Christoffel"[1,-1,-1],CVSimplify->Simplify]
MetricCompute[met,ch,"Ricci"[-1,-1],CVSimplify->Simplify]
eval[A_]:=A//ToBasis[ch]//ComponentArray//ToBasis[ch]//TraceBasisDummy//ToValues//Simplify
cd=CovDOfMetric[met];

(*The field's energy momentum tensor*)
coordinates={t[]->t,r[]->r,\[Theta][]->\[Theta],\[Phi][]->\[Phi]};
DArule=Table[PDch[{i,-ch}]@Aplace[{k,-ch}]->ToExpression["D"<>ToString[i]<>ToString["A"]<>ToString[k]],{i,0,3},{k,0,3}]//Flatten;
AtoAi=Table[{Aplace[{i,-ch}]->ToExpression["A"<>ToString[i]]},{i,0,3}]//Flatten;
DefConstantSymbol[\[Mu]]
DefTensor[Aplace[-b],M]
DefTensor[FS[-b,-c],M]
DefTensor[T[-b,-c],M]
(*The 2 in front of everything makes sure that G_mn=T_mn (Einstein equation) is satisfied (without additional proportionality constant in our units), with the usual definition of the Einstein tensor G_mn*)
FS/:FS[-b_,-c_]:=CD[-b]@Aplace[-c]-CD[-c]@Aplace[-b] 
(*Tdd=(ChangeCovD[(\[Mu]^2Aplace[-h]Aplace[-c]+met[d,e]FS[-h,-d]FS[-c,-e]-1/4met[-h,-c](FS[-d,-e]FS[-f,-g]met[d,f]met[e,g]+2\[Mu]^2met[e,d]Aplace[-d]Aplace[-e])),CD,PDch]//eval)//.DArule//.AtoAi;
Export["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/Tdd_PrivateTeukolsky.mx",Tdd];*)
Tdd=Import["./Tdd_PrivateTeukolsky.mx"]//.coordinates//.{Global`\[Chi]->PrivateTeukolsky`\[Chi]}/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};

(*Input field into energy momentum*)
FieldEnergyMomentum[AdownNorm1_,Procamass_,spin_]:=Module[{Adown=AdownNorm1,\[Mu]n=Procamass,\[Chi]m=spin},
prec=SetPrecision[#,20]&;
adowntemp=AdownNorm1//.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};
AdownInT[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=adowntemp//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
Aexpl=Table[ToExpression["PrivateTeukolsky`A"<>ToString[j]]->AdownInT[t,r,\[Theta],\[Phi]][[j+1]],{j,0,3}]//Flatten;
dAexpl=Table[{ToExpression["PrivateTeukolsky`D"<>ToString[i]<>ToString["A"]<>ToString[k]]->D[AdownInT[t,r,\[Theta],\[Phi]][[k+1]],{t,r,\[Theta],\[Phi]}[[i+1]]]},{i,0,3},{k,0,3}]//Flatten//Simplify;
Tddtemp=(Tdd//.coordinates/.Aexpl/.dAexpl//.{\[Chi]->spin,\[Mu]->Procamass});
Tdowndown[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Tddtemp/.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
];

(*The cubic interaction of the Proca*)
cubic=(ChangeCovD[met[b,d]Aplace[-b]CD[-c]@Aplace[-d],CD,PDch]//eval)//.DArule//.AtoAi;
CubicTerm[AdownNorm1_,Procamass_,spin_,\[Beta]sta_,\[Beta]sto_]:=Module[{Adown=AdownNorm1,\[Mu]n=Procamass,\[Chi]m=spin,bsta=\[Beta]sta,bsto=\[Beta]sto},
prec=SetPrecision[#,20]&;
AdownInT[t_,r_,\[Theta]_,\[Phi]_]:=Adown//.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};
Aexpl=Table[{ToExpression["PrivateTeukolsky`A"<>ToString[i]]->AdownInT[t,r,\[Theta],\[Phi]][[i+1]]},{i,0,3}]//Flatten;
dAexpl=Table[{ToExpression["PrivateTeukolsky`D"<>ToString[i]<>ToString["A"]<>ToString[k]]->D[AdownInT[t,r,\[Theta],\[Phi]][[k+1]],{t,r,\[Theta],\[Phi]}[[i+1]]]},{i,0,3},{k,0,3}]//Flatten//Simplify;
cubictemp=prec@(cubic//.coordinates/.Aexpl/.dAexpl//.{\[Chi]->spin,\[Mu]->Procamass});
CubicTermOutput[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=cubictemp/.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
Dcubic[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Procamass^2 AdownInT[t,r,\[Theta],\[Phi]]+beta CubicTermOutput[t,r,\[Theta],\[Phi]]//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
];

Print["---------------------------------------------------------------------"];
Print["Generating analytical Teukolsky T: start."];
(*Kinnersley tetrad:*)
DefTensor[nK[b],M]
DefTensor[mK[b],M]
DefTensor[mbarK[b],M]
nKtemp=1/(2\[CapitalSigma]) {r[]^2+\[Chi]^2,-\[CapitalDelta],0,\[Chi]}//Simplify;
ComponentValue[nK[b]//ToBasis[ch]//ComponentArray,nKtemp];
mKtemp={I \[Chi] Sin[\[Theta][]],0,1,I/Sin[\[Theta][]]}/(Sqrt[2](r[]+I \[Chi] Cos[\[Theta][]]))//Simplify;
ComponentValue[mK[b]//ToBasis[ch]//ComponentArray,mKtemp];
mbarKtemp={-I \[Chi] Sin[\[Theta][]],0,1,-I/Sin[\[Theta][]]}/(Sqrt[2](r[]-I \[Chi] Cos[\[Theta][]]))//Simplify;
ComponentValue[mbarK[b]//ToBasis[ch]//ComponentArray,mbarKtemp];


(*Integrad for the T_lm\[Omega] for general T_nn etc.*)
DefConstantSymbol[\[Omega]T]
DefConstantSymbol[mT]
\[Rho]=1/(r[]-I \[Chi] Cos[\[Theta][]]);
\[Rho]bar=1/(r[]+I \[Chi] Cos[\[Theta][]]);
L0=(PD[{2,-ch}]@#+(mT #)/Sin[\[Theta][]]-\[Omega]T \[Chi] Sin[\[Theta][]]#)&;
Lm1op=(PD[{2,-ch}]@#+(mT #)/Sin[\[Theta][]]-\[Omega]T \[Chi] Sin[\[Theta][]]#-Cot[\[Theta][]]#)&;
Jp=(PD[{1,-ch}]@#+(I/\[CapitalDelta])((r[]^2+\[Chi]^2)\[Omega]T-\[Chi] mT)#)&;
B21=-1/2 \[Rho]^8 \[Rho]bar Lm1op@(\[Rho]^-4 L0@(\[Rho]^-2 \[Rho]bar^-1 Tnn[t[],r[],\[Theta][],\[Phi][]]));
B22=-(1/(2Sqrt[2]))(\[Rho]^8 \[Rho]bar \[CapitalDelta]^2 Lm1op@(\[Rho]^-4 \[Rho]bar^2 Jp@(\[Rho]^-2 \[Rho]bar^-2 \[CapitalDelta]^-1 Tmn[t[],r[],\[Theta][],\[Phi][]])));
B2=B21+B22;
B2s1=-1/4 \[Rho]^8 \[Rho]bar \[CapitalDelta]^2 Jp@(\[Rho]^-4 Jp@(\[Rho]^-2 \[Rho]bar Tmm[t[],r[],\[Theta][],\[Phi][]]));
B2s2=-(1/(2Sqrt[2])) \[Rho]^8 \[Rho]bar \[CapitalDelta]^2 Jp@(\[Rho]^-4 \[Rho]bar^2 \[CapitalDelta]^-1 Lm1op@(\[Rho]^-2 \[Rho]bar^-2 Tmn[t[],r[],\[Theta][],\[Phi][]]));
B2s=B2s1+B2s2;
dTij=Table[{PD[{i,-ch}]@(ToExpression["PrivateTeukolsky`T"<>ToString[j]][t[],r[],\[Theta][],\[Phi][]])->ToExpression["d"<>ToString[i]<>ToString["T"]<>ToString[j]]},{i,1,2},{j,{"nn","mn","mm"}}]//Flatten;
ddTij=Table[{PD[{k,-ch}]@(PD[{i,-ch}]@(ToExpression["PrivateTeukolsky`T"<>ToString[j]][t[],r[],\[Theta][],\[Phi][]]))->ToExpression["d"<>ToString[k]<>ToString["d"]<>ToString[i]<>ToString["T"]<>ToString[j]]},{i,1,2},{k,1,2},{j,{"nn","mn","mm"}}]//Flatten;
TeukolskyTintegrand1=4/(\[Rho]bar \[Rho]^5) (B2+B2s)//.ddTij//.dTij//.coordinates; (*The 1/Sqrt[2\[Pi]] is not part of the FT below, so we include it here!*)
Export["./TeuTintegrand.mx",TeukolskyTintegrand1];
TeukolskyTintegrand=Import["./TeuTintegrand.mx"];
Print["Generating analytical Teukolsky T: done."];
Print["---------------------------------------------------------------------"];

(*Plugging in explicit solutions for the fields*)
TeukolskyTlmw[Tfielddowndown1_,Procamass_,spin_,mTeu_,\[Omega]Teu_,rstop1_,rstart1_,\[Theta]stop1_,\[Theta]start1_]:=Module[{Tfielddowndown=Tfielddowndown1,\[Mu]1=Procamass,\[Chi]1=spin},
If[mTeu<0, Print["Error: m<0, use symmetry relations to get these source modes."];Abort[];,{}];
If[\[Omega]Teu<0, Print["Error: Re(\[Omega])<0, use symmetry relations to get these source modes."];Abort[];,{}];
prec=SetPrecision[#,15]&;
coordinates={t[]->t,r[]->r,\[Theta][]->\[Theta],\[Phi][]->\[Phi]};

para=prec@{\[Chi]->spin,\[Mu]->Procamass,mT->mTeu,\[Omega]T->\[Omega]Teu};
rstart=prec@rstart1;
rstop=prec@rstop1;
\[Theta]start=prec@\[Theta]start1;
\[Theta]stop=prec@\[Theta]stop1;
rmesh=99;
thmesh=49; (*This HAS to be an odd number to hit the last \[Theta]stop exactly*)
rdatapoint[x_]:=(rstop-rstart)/((rmesh+1)^4-1) x^4+rstart-(rstop-rstart)/((rmesh+1)^4-1);
rrange=prec@rdatapoint[Range[1,rmesh+1]];
\[Theta]range=prec@Range[\[Theta]start,\[Theta]stop,(\[Theta]stop-\[Theta]start)/thmesh];

(*Generating T_nn etc. for the given pair (\[Omega]T,mT)*)
tnnexpltemp=(nKtemp//.coordinates) . ((Tfielddowndown/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]}) . (nKtemp//.coordinates));
TnnExpl[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=(tnnexpltemp)/.{\[Chi]->spin,\[Mu]->Procamass}//.coordinates/.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
tmnexpltemp=(nKtemp//.coordinates) . ((Tfielddowndown/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]}) . (mbarKtemp//.coordinates));
TmnExpl[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=(tmnexpltemp)/.{\[Chi]->spin,\[Mu]->Procamass}//.coordinates/.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
tmmexpltemp=(mbarKtemp//.coordinates) . ((Tfielddowndown/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]}) . (mbarKtemp//.coordinates));
TmmExpl[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=(tmmexpltemp)/.{\[Chi]->spin,\[Mu]->Procamass}//.coordinates/.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};

Print["---------------------------------------------------------------------"];
Print["Generating: Tnn..."];
tnntemp=Monitor[Table[(-1)^mTeu 2\[Pi] FourierCoefficient[FourierTransform[TnnExpl[t,rrange[[j]],\[Theta]range[[k]],\[Phi]+\[Pi]]//TrigReduce,t,\[Omega]T]/.{DiracDelta[a_]:>1},\[Phi],mTeu],{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],{j,k}];
tnntempR=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tnntemp[[j,k]]//Re},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TnnInterNumRe[rr_,\[Theta]\[Theta]_]:=Interpolation[tnntempR,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
tnntempI=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tnntemp[[j,k]]//Im},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TnnInterNumIm[rr_,\[Theta]\[Theta]_]:=Interpolation[tnntempI,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
TnnInter[rr_,\[Theta]\[Theta]_]:=TnnInterNumRe[r,\[Theta]]+I TnnInterNumIm[r,\[Theta]]/.{r->rr,\[Theta]->\[Theta]\[Theta]};

Print["Generating: Tmn..."];
tmntemp=Table[(-1)^mTeu 2\[Pi] FourierCoefficient[FourierTransform[TmnExpl[t,rrange[[j]],\[Theta]range[[k]],\[Phi]+\[Pi]]//TrigReduce,t,\[Omega]T]/.{DiracDelta[a_]:>1},\[Phi],mTeu],{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}];
tmntempR=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tmntemp[[j,k]]//Re},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TmnInterNumRe[rr_,\[Theta]\[Theta]_]:=Interpolation[tmntempR,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
tmntempI=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tmntemp[[j,k]]//Im},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TmnInterNumIm[rr_,\[Theta]\[Theta]_]:=Interpolation[tmntempI,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
TmnInter[rr_,\[Theta]\[Theta]_]:=TmnInterNumRe[r,\[Theta]]+I TmnInterNumIm[r,\[Theta]]/.{r->rr,\[Theta]->\[Theta]\[Theta]};

Print["Generating: Tmm..."];
tmmtemp=Table[(-1)^mTeu 2\[Pi] FourierCoefficient[FourierTransform[TmmExpl[t,rrange[[j]],\[Theta]range[[k]],\[Phi]+\[Pi]]//TrigReduce,t,\[Omega]T]/.{DiracDelta[a_]:>1},\[Phi],mTeu],{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}];
tmmtempR=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tmmtemp[[j,k]]//Re},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TmmInterNumRe[rr_,\[Theta]\[Theta]_]:=Interpolation[tmmtempR,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
tmmtempI=Flatten[Table[{{rrange[[j]],\[Theta]range[[k]]},tmmtemp[[j,k]]//Im},{j,1,Length[rrange]},{k,1,Length[\[Theta]range]}],1];
TmmInterNumIm[rr_,\[Theta]\[Theta]_]:=Interpolation[tmmtempI,InterpolationOrder->3,Method->"Spline"][rr,\[Theta]\[Theta]];
TmmInter[rr_,\[Theta]\[Theta]_]:=TmmInterNumRe[r,\[Theta]]+I TmmInterNumIm[r,\[Theta]]/.{r->rr,\[Theta]->\[Theta]\[Theta]};
Print["Generating: Done."];
Print["---------------------------------------------------------------------"];

(*Computing integrand of Tlm\[Omega]*)
TijtoExpl=Table[{ToExpression["PrivateTeukolsky`T"<>ToString[j]][t,r,\[Theta],\[Phi]]->ToExpression["T"<>ToString[j]<>ToString["Inter"]][r,\[Theta]]},{j,{"nn","mn","mm"}}]//Flatten;
dTijtoExpl=Table[{ToExpression["PrivateTeukolsky`d"<>ToString[i]<>ToString["T"]<>ToString[j]]->D[ToExpression["T"<>ToString[j]<>ToString["Inter"]][r,\[Theta]],{r,\[Theta]}[[i]]]},{i,{1,2}},{j,{"nn","mn","mm"}}]//Flatten;
ddTijtoExpl=Table[{ToExpression["PrivateTeukolsky`d"<>ToString[i]<>ToString["d"]<>ToString[k]<>ToString["T"]<>ToString[j]]->D[D[ToExpression["T"<>ToString[j]<>ToString["Inter"]][r,\[Theta]],{r,\[Theta]}[[k]]],{r,\[Theta]}[[i]]]},{i,1,2},{k,1,2},{j,{"nn","mn","mm"}}]//Flatten;
TTtemp=(TeukolskyTintegrand/.{\[Chi]->spin,\[Mu]->Procamass,\[Omega]T->\[Omega]Teu,mT->mTeu}//.ddTijtoExpl//.dTijtoExpl//.TijtoExpl);
TeukolskyTintegrandExpl[rr_,\[Theta]\[Theta]_]:=TTtemp/.{r->rr,\[Theta]->\[Theta]\[Theta]};

(*Generating spheroidal harmonics and Tlm\[Omega] for l=|m|,|m|+1,|m|+2*)
Monitor[Do[
	Print["Generating: Spheroidal harmonic for (s=-2,\[Chi],Re(\[Omega]_Teu))."];
	thetarange=prec@Range[0+10^-4,\[Pi]-10^-4,(\[Pi]-2 10^-4)/400];
	SWSHdata=Table[{x,Sqrt[2\[Pi]]SpinWeightedSpheroidalHarmonicS[-2,lTeu,mTeu,spin \[Omega]Teu,x,0]},{x,thetarange}];
	Sbar[\[Theta]\[Theta]_]:=Interpolation[SWSHdata,InterpolationOrder->2,Method->"Hermite"][\[Theta]\[Theta]];
	alm=SpinWeightedSpheroidalEigenvalue[-2,lTeu,mTeu,spin \[Omega]Teu]+2mTeu spin \[Omega]Teu-(spin \[Omega]Teu)^2;
	Print[ToString["Eigenvalues of SWSH: "<>ToString[alm]]];
	sinfunc[x_]:=Sin[x];
	plotsbar[lTeu]=Plot[Sbar[x]sinfunc[x],{x,10^-4,\[Pi]-10^-4}];
	rrangeout=rrange;
	Print["Generating: Interpolation data..."];
	Tlm\[Omega]temp=Table[{r,NIntegrate[TeukolskyTintegrandExpl[r,\[Theta]]Sbar[\[Theta]]sinfunc[\[Theta]],{\[Theta],\[Theta]start,\[Theta]stop},MaxRecursion->50,WorkingPrecision->5]},{r,rrange}];
	Print["Generating: Tlm\[Omega]."];
	Tlm\[Omega][lTeu]=Interpolation[Tlm\[Omega]temp,InterpolationOrder->2,Method->"Hermite"];
,{lTeu,Abs[mTeu],Abs[mTeu]+1}],lTeu];
Print["Generating: Done."];
Print["---------------------------------------------------------------------"];
];

Print["TeukolskyT Functions:"];
Print["TeukolskyTlmw[Tdowndown[t,r,\[Theta],\[Phi]],Procamass,spin,Teu-m-mode,Teu-Re(\[Omega]),rstop,rstart,\[Theta]stop,\[Theta]start]"];
Print["Tlm\[Omega][r]"];
Print["---------------------------------------------------------------------"];


End[]

EndPackage[]
