(* ::Package:: *)

BeginPackage["TeukolskyTsingleMode`"]

Tlm\[Omega]data
Tnn3Inter
TeukolskyTlmwSingleMode
TeukolskyTintegrand
TnnExplSimp
TTtemp1
TTtemp2

Begin["PrivateTeukolskysingleMode`"]
Needs["xAct`xCoba`"]
Needs["SpinWeightedSpheroidalHarmonics`","./SWSHs.m"]

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
DefConstantSymbol[m]
DefConstantSymbol[\[Omega]r]

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

(*Field's Tmunu*)
coordinates={t[]->t,r[]->r,\[Theta][]->\[Theta],\[Phi][]->\[Phi]};
DefTensor[Aplace[-b],M]
DArule=Table[PDch[{i,-ch}]@Aplace[{k,-ch}]->ToExpression["D"<>ToString[i]<>ToString[A]<>ToString[k]],{i,0,3},{k,0,3}]//Flatten;
AtoAi=Table[{Aplace[{i,-ch}]->ToExpression["A"<>ToString[i]]},{i,0,3}]//Flatten;
DefConstantSymbol[\[Mu]]
DefTensor[FS[-b,-c],M]
DefTensor[T[-b,-c],M]
(*The 2 in front of everything makes sure that G_mn=T_mn (Einstein equation) is satisfied (without additional proportionality constant in our units), with the usual definition of the Einstein tensor G_mn*)
FS/:FS[-b_,-c_]:=CD[-b]@Aplace[-c]-CD[-c]@Aplace[-b];
(*Tdd=(ChangeCovD[(\[Mu]^2Aplace[-b]Aplace[-c]+met[d,e]FS[-b,-d]FS[-c,-e]-1/4met[-b,-c](FS[-d,-e]FS[-f,-g]met[d,f]met[e,g]+2\[Mu]^2met[e,d]Aplace[-d]Aplace[-e])),CD,PDch]//eval)//.DArule//.AtoAi;
Export["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/Tdd_singlemode.mx",Tdd];*)
Tdd=Import["./Tdd_singlemode.mx"]//.coordinates//.{Global`\[Chi]->PrivateTeukolskysingleMode`\[Chi]}/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};

(*Kinnersley tetrad:*)
DefTensor[nK[b],M]
DefTensor[mK[b],M]
DefTensor[mbarK[b],M]
nn=1/(2\[CapitalSigma]) {r[]^2+\[Chi]^2,-\[CapitalDelta],0,\[Chi]}//Simplify;
ComponentValue[nK[b]//ToBasis[ch]//ComponentArray,nn];
mm={I \[Chi] Sin[\[Theta][]],0,1,I/Sin[\[Theta][]]}/(Sqrt[2](r[]+I \[Chi] Cos[\[Theta][]]))//Simplify;
ComponentValue[mK[b]//ToBasis[ch]//ComponentArray,mm];
mbarmbar={-I \[Chi] Sin[\[Theta][]],0,1,-I/Sin[\[Theta][]]}/(Sqrt[2](r[]-I \[Chi] Cos[\[Theta][]]))//Simplify;
ComponentValue[mbarK[b]//ToBasis[ch]//ComponentArray,mbarmbar];

(*The analytic ground work for the Teukolksy Tlmw*)
Asincos=Table[{ToExpression["PrivateTeukolskysingleMode`A"<>ToString[i]]->ToExpression["PrivateTeukolskysingleMode`Ac"<>ToString[i]][r,\[Theta]]Cos[m \[Phi]-\[Omega]r t]+ToExpression["PrivateTeukolskysingleMode`As"<>ToString[i]][r,\[Theta]]Sin[m \[Phi]-\[Omega]r t]},{i,0,3}]//Flatten;
dAsincos=Table[{ToExpression["PrivateTeukolskysingleMode`D"<>ToString[i]<>ToString["PrivateTeukolskysingleMode`A"]<>ToString[k]]->D[ToExpression["PrivateTeukolskysingleMode`A"<>ToString[k]]//.Asincos,{t,r,\[Theta],\[Phi]}[[i+1]]]},{i,0,3},{k,0,3}]//Flatten//Simplify;
trigtoexp={Sin[m \[Phi]-t \[Omega]r]->-(1/2) I (-E^(-I (m \[Phi]-t \[Omega]r))+E^(I (m \[Phi]-t \[Omega]r))),Cos[m \[Phi]-t \[Omega]r]->1/2 (E^(-I (m \[Phi]-t \[Omega]r))+E^(I (m \[Phi]-t \[Omega]r)))};

nntemp1=(nK[b]//eval).(Tdd.(nK[c]//eval))//.Asincos//.dAsincos//.coordinates;
nntemp2=(nntemp1//.trigtoexp//Expand)//.{Exp[-2I m \[Phi]+2I t \[Omega]r]->expminus,Exp[2I m \[Phi]-2I t \[Omega]r]->expplus};
TmmExplSimp=CoefficientList[nntemp2,{expminus,expplus}];

mntemp1=(nK[b]//eval).(Tdd.(mbarK[c]//eval))//.Asincos//.dAsincos//.coordinates;
mntemp2=(mntemp1//.trigtoexp//Expand)//.{Exp[-2I m \[Phi]+2I t \[Omega]r]->expminus,Exp[2I m \[Phi]-2I t \[Omega]r]->expplus};
TmnExplSimp=CoefficientList[mntemp2,{expminus,expplus}];

mmtemp1=(mbarK[b]//eval).(Tdd.(mbarK[c]//eval))//.Asincos//.dAsincos//.coordinates;
mmtemp2=(mmtemp1//.trigtoexp//Expand)//.{Exp[-2I m \[Phi]+2I t \[Omega]r]->expminus,Exp[2I m \[Phi]-2I t \[Omega]r]->expplus};
TnnExplSimp=CoefficientList[mmtemp2,{expminus,expplus}];

(*Towards the Teukolsky T*)
DefScalarFunction[Tnn];
DefScalarFunction[Tmn];
DefScalarFunction[Tmm];
Tnnin=Tnn[r[],\[Theta][]]Exp[-2I m \[Phi][]+2I \[Omega]r t[]];
Tmnin=Tmn[r[],\[Theta][]]Exp[-2I m \[Phi][]+2I \[Omega]r t[]];
Tmmin=Tmm[r[],\[Theta][]]Exp[-2I m \[Phi][]+2I \[Omega]r t[]];
\[Rho]=1/(r[]-I \[Chi] Cos[\[Theta][]]);
\[Rho]bar=1/(r[]+I \[Chi] Cos[\[Theta][]]);
L0=(PD[{2,-ch}]@#-I/Sin[\[Theta][]]PD[{3,-ch}]@#-I \[Chi] Sin[\[Theta][]]PD[{0,-ch}]@#)&;
Lm1op=(PD[{2,-ch}]@#-I/Sin[\[Theta][]]PD[{3,-ch}]@#-I \[Chi] Sin[\[Theta][]] PD[{0,-ch}]@#-Cot[\[Theta][]]#)&;
Jp=(PD[{1,-ch}]@#-1/\[CapitalDelta]((r[]^2+\[Chi]^2)PD[{0,-ch}]@#+\[Chi] PD[{3,-ch}]@#))&;
B21=-1/2\[Rho]^8\[Rho]bar Lm1op@(\[Rho]^-4L0@(\[Rho]^-2\[Rho]bar^-1Tnnin));
B22=-(1/(2Sqrt[2]))(\[Rho]^8\[Rho]bar \[CapitalDelta]^2Lm1op@(\[Rho]^-4\[Rho]bar^2Jp@(\[Rho]^-2\[Rho]bar^-2\[CapitalDelta]^-1Tmnin)));
B2=B21+B22;
B2s1=-1/4\[Rho]^8\[Rho]bar \[CapitalDelta]^2Jp@(\[Rho]^-4Jp@(\[Rho]^-2\[Rho]bar Tmmin));
B2s2=-(1/(2Sqrt[2]))\[Rho]^8\[Rho]bar \[CapitalDelta]^2Jp@(\[Rho]^-4\[Rho]bar^2\[CapitalDelta]^-1Lm1op@(\[Rho]^-2\[Rho]bar^-2Tmnin));
B2s=B2s1+B2s2;
TeukolskyTintegrand=(4 \[Rho]^-5)/(\[Rho]bar Sqrt[2\[Pi]])(B2+B2s);(*Does not have compact support*)
Export["./TeukolskyTintegrandsingleMode.mx",TeukolskyTintegrand];
TeukolskyTintegrand=Import["./TeukolskyTintegrandsingleMode.mx"]//.coordinates//.{Global`\[Chi]->PrivateTeukolskysingleMode`\[Chi],Global`\[Omega]r->PrivateTeukolskysingleMode`\[Omega]r,Global`m->PrivateTeukolskysingleMode`m}/.{Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};


TeukolskyTlmwSingleMode[AdownNorm_,rtst_,rstop_,\[Theta]start_,\[Theta]stop_,mT1_,lT1_,spin_,\[Omega]real_]:=Module[{AdownNorm1=AdownNorm,rstart=rtst,rstop1=rstop,\[Theta]start1=\[Theta]start,\[Theta]stop1=\[Theta]stop,\[Chi]1=spin,mTeu=mT1,lTeu=lT1},
collect={Tmm[r[],\[Theta][]],Tnn[r[],\[Theta][]],Tmn[r[],\[Theta][]],Table[D[D[ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]][r,\[Theta]],{r,\[Theta]}[[k1]]],{r,\[Theta]}[[k2]]]/.{r->r[],\[Theta]->\[Theta][]},{i,{"mm","nn","mn"}},{k1,1,2},{k2,1,2}]}//Flatten;
TTtemp1=(TeukolskyTintegrand//Expand)//.{Exp[-2I \[Omega]r t[]+2I m \[Phi][]]->expminus2,Exp[2I \[Omega]r t[]-2I m \[Phi][]]->expplus2};
TTtemp2=CoefficientList[TTtemp1,{expminus2,expplus2}];

TeukolskyTFT\[Phi]Int=2\[Pi] TTtemp2[[2,1]]DiracDelta[\[Omega]T+2\[Omega]r]KroneckerDelta[2m,mT];

parameters={mT->mTeu,lT->lTeu,\[Chi]->spin,\[Omega]r->\[Omega]real};
Adownvector[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=AdownNorm//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};

sinexpl=Sin[\[Omega]r t-m \[Phi]]//.parameters;
cosexpl=Cos[\[Omega]r t-m \[Phi]]//.parameters;
Acs=Table[{ToExpression["PrivateTeukolskysingleMode`Ac"<>ToString[k]][r,\[Theta]]->(CoefficientList[Adownvector[t,r,\[Theta],\[Phi]][[k+1]],{sinexpl,cosexpl}][[1,2]]),ToExpression["PrivateTeukolskysingleMode`As"<>ToString[k]][r,\[Theta]]->(CoefficientList[Adownvector[t,r,\[Theta],\[Phi]][[k+1]],{sinexpl,cosexpl}][[2,1]])},{k,0,3}]//.{t->0}//Flatten;
Acssq=Table[{ToExpression["PrivateTeukolskysingleMode`Ac"<>ToString[k]][r,\[Theta]]^2->(CoefficientList[Adownvector[t,r,\[Theta],\[Phi]][[k+1]],{sinexpl,cosexpl}][[1,2]])^2,ToExpression["PrivateTeukolskysingleMode`As"<>ToString[k]][r,\[Theta]]^2->(CoefficientList[Adownvector[t,r,\[Theta],\[Phi]][[k+1]],{sinexpl,cosexpl}][[2,1]])^2},{k,0,3}]//.{t->0}//Flatten;
dAcs={Table[D[ToExpression["PrivateTeukolskysingleMode`As"<>ToString[i]][r,\[Theta]],{r,\[Theta]}[[k]]]->D[ToExpression["PrivateTeukolskysingleMode`As"<>ToString[i]][r,\[Theta]]//.Acs,{r,\[Theta]}[[k]]],{i,0,3},{k,1,2}],Table[D[ToExpression["PrivateTeukolskysingleMode`Ac"<>ToString[i]][r,\[Theta]],{r,\[Theta]}[[k]]]->D[ToExpression["Ac"<>ToString[i]][r,\[Theta]]//.Acs,{r,\[Theta]}[[k]]],{i,0,3},{k,1,2}]}//Flatten;
dAcssq={Table[(D[ToExpression["PrivateTeukolskysingleMode`As"<>ToString[i]][r,\[Theta]],{r,\[Theta]}[[k]]])^2->(D[ToExpression["PrivateTeukolskysingleMode`As"<>ToString[i]][r,\[Theta]]//.Acs,{r,\[Theta]}[[k]]])^2,{i,0,3},{k,1,2}],Table[(D[ToExpression["PrivateTeukolskysingleMode`Ac"<>ToString[i]][r,\[Theta]],{r,\[Theta]}[[k]]])^2->(D[ToExpression["Ac"<>ToString[i]][r,\[Theta]]//.Acs,{r,\[Theta]}[[k]]])^2,{i,0,3},{k,1,2}]}//Flatten;

\[Theta]mesh=39;
rmesh=39;
rdatapoint[x_]:=(rstop-rtst)/((rmesh+1)^4-1) x^4+rtst-(rstop-rtst)/((rmesh+1)^4-1);
rrange=rdatapoint[Range[1,rmesh+1]];
\[Theta]range=Range[\[Theta]start,\[Theta]stop,(\[Theta]stop-\[Theta]start)/\[Theta]mesh];

Print["Constructing stress-energy NP-scalars..."];
Print["Generating Tmn..."];
Tmm1func=TmmExplSimp[[1,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tmm2func=TmmExplSimp[[2,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tmm3func=TmmExplSimp[[1,2]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Print["Generating Tnn..."];
Tnn1func=TnnExplSimp[[1,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tnn2func=TnnExplSimp[[2,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tnn3func=TnnExplSimp[[1,2]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Print["Generating Tmn..."];
Tmn1func=TmnExplSimp[[1,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tmn2func=TmnExplSimp[[2,1]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Tmn3func=TmnExplSimp[[1,2]]//.coordinates//.parameters//.Acssq//.dAcssq//.Acs//.dAcs;
Print["--------------------------------------------------------------"];

Tmm1discr=ConstantArray[0,{rmesh+1,\[Theta]mesh+1}];
Tmm2discr=Tmm1discr;Tmm3discr=Tmm1discr;
Tmn1discr=Tmm1discr;Tmn2discr=Tmm1discr;Tmn3discr=Tmm1discr;
Tnn1discr=Tmm1discr;Tnn2discr=Tmm1discr;Tnn3discr=Tmm1discr;

Print["Producing data for interpolating function..."];
Monitor[Do[
rinput=rrange[[ir]];
\[Theta]input=\[Theta]range[[i\[Theta]]];
Tmm1discr[[ir,i\[Theta]]]=Tmm1func//.{r->rinput,\[Theta]->\[Theta]input};
Tmm2discr[[ir,i\[Theta]]]=Tmm2func//.{r->rinput,\[Theta]->\[Theta]input};
Tmm3discr[[ir,i\[Theta]]]=Tmm3func//.{r->rinput,\[Theta]->\[Theta]input};

Tmn1discr[[ir,i\[Theta]]]=Tmn1func//.{r->rinput,\[Theta]->\[Theta]input};
Tmn2discr[[ir,i\[Theta]]]=Tmn2func//.{r->rinput,\[Theta]->\[Theta]input};
Tmn3discr[[ir,i\[Theta]]]=Tmn3func//.{r->rinput,\[Theta]->\[Theta]input};

Tnn1discr[[ir,i\[Theta]]]=Tnn1func//.{r->rinput,\[Theta]->\[Theta]input};
Tnn2discr[[ir,i\[Theta]]]=Tnn2func//.{r->rinput,\[Theta]->\[Theta]input};
Tnn3discr[[ir,i\[Theta]]]=Tnn3func//.{r->rinput,\[Theta]->\[Theta]input};
,{ir,1,rmesh+1},{i\[Theta],1,\[Theta]mesh+1}],{ir,i\[Theta]}];

Print["Generating interpolating functions.."];
Tmm1Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tmm2Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tmm3Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];

Tmn1Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tmn2Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tmn3Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];

Tnn1Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tnn2Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Tnn3Inter=Interpolation[Flatten[Table[{{rrange[[i]],\[Theta]range[[j]]},Tmm1discr[[i,j]]},{i,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->1];
Print["Done."];
Print["--------------------------------------------------------------"];


Tiireplace=Table[{ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]][r[],\[Theta][]]->ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]<>ToString[Inter]][r,\[Theta]]},{i,{"mm","mn","nn"}},{j,1,3}]//Flatten;
dTiireplace=Table[{D[ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]][r[],\[Theta][]],{r[],\[Theta][]}[[k1]]]->D[ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]][r[],\[Theta][]]//.Tiireplace//.coordinates,{r,\[Theta]}[[k1]]]},{i,{"mm","mn","nn"}},{j,1,3},{k1,1,2}]//Flatten;
ddTiireplace=Table[{D[D[ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]][r[],\[Theta][]],{r[],\[Theta][]}[[k1]]],{r[],\[Theta][]}[[k2]]]->D[D[ToExpression["PrivateTeukolskysingleMode`T"<>ToString[i]<>ToString[j]][r[],\[Theta][]]//.Tiireplace//.coordinates,{r,\[Theta]}[[k1]]],{r,\[Theta]}[[k2]]]},{i,{"mm","mn","nn"}},{j,1,3},{k1,1,2},{k2,1,2}]//Flatten;

temp1=TeukolskyTFT\[Phi]Int//.{DiracDelta[2 \[Omega]r+\[Omega]T]KroneckerDelta[2m,mT]->Dd3};
Tfor\[Theta]int=DeleteCases[CoefficientList[temp1,{Dd3}]//Flatten,0]//.parameters//.Tiireplace//.dTiireplace//.ddTiireplace//.coordinates;
Tlm\[Omega]=ConstantArray[0,{1,rmesh+1}];

SWSHrun[mT,lT,-2,spin,\[Omega]real];

Do[
Tlm\[Omega]data[[i,ir]]=NIntegrate[Tfor\[Theta]int[[i]]SWSH[\[Theta]]Sin[\[Theta]]//.coordinates//.{\[Theta]->theta,r->rrange[[ir]]},{theta,\[Theta]start,\[Theta]stop}]
,{i,1,1},{ir,1,rmesh+1}];
]


End[]

EndPackage[]
