(* ::Package:: *)

BeginPackage["VectorFieldNorm`"]

AupNorm::usage = "A^\mu[t,r,\[Theta],\[Phi]] :> Normalized contravariant vector field solution for given variables."
AdownNorm::usage = "A_\mu[t,r,\[Theta],\[Phi]] :> Normalized covariant vector field solution for given variables."
VectorField::usage = "VectorField[Rsol(interpolating func.),Anglsol(explicit in terms of \[Theta]),spin,m,Re(\[Nu]),Im(\[Nu]),Re(\[Omega]),\[Mu],nh,MBH(for normalization)] :> AupNorm, AdownNorm, edensityNorm"
edensityNorm::usage = "\[Rho][t,r,\[Theta],\[Phi]] :> Normalized energy density of field mode."
rtst
rstop
\[Theta]start
\[Theta]stop
Anorm
massNorm
massNormint
dataout
functionplotsR
functionplotsS


BBLtemp1
BBLtemp2

Rad
dRad
ddRad
Sangl
dSangl
ddSangl

Begin["Private`"]
Needs["xAct`xCoba`"]

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

(*Defining the Kerr metric*)
\[CapitalSigma]=r[]^2+\[Chi]^2 Cos[\[Theta][]]^2;
\[CapitalDelta]=r[]^2-2 r[]+\[Chi]^2;
KerrBL={{-(1-(2r[])/\[CapitalSigma]),0,0,(-2 r[] \[Chi] Sin[\[Theta][]]^2)/\[CapitalSigma]},{0,\[CapitalSigma]/\[CapitalDelta],0,0},{0,0,\[CapitalSigma],0},{(-2 r[] \[Chi] Sin[\[Theta][]]^2)/\[CapitalSigma],0,0,(r[]^2+\[Chi]^2+2 r[]/\[CapitalSigma] \[Chi]^2 Sin[\[Theta][]]^2)Sin[\[Theta][]]^2}};
ComponentValue[met[-b,-c]//ToBasis[ch]//ComponentArray,KerrBL]
MetricCompute[met,ch,All]
MetricCompute[met,ch,"Christoffel"[1,-1,-1],CVSimplify->Simplify]
MetricCompute[met,ch,"Ricci"[-1,-1],CVSimplify->Simplify]
eval[A_]:=A//ToBasis[ch]//ComponentArray//ToBasis[ch]//TraceBasisDummy//ToValues//Simplify
cd=CovDOfMetric[met];

(*Emulation of coordinate transformation from Dolan's coordinates to Boyer-Lindquist coordinates*)
DefTensor[\[CapitalLambda][b,-c],M]
\[CapitalLambda]updownFKKStoBL={{1,0,0,s^2},{0,1,0,0},{0,0,-1/(s Sin[\[Theta][]]),0},{0,0,0,s}};
ComponentValue[\[CapitalLambda][b,-c]//ToBasis[ch]//ComponentArray,\[CapitalLambda]updownFKKStoBL];

(*Polarization tensor in Dolan's coordinates (differ from FKKS by a complex unit in front of r)*)
DefTensor[BDolan[b,c],M]
BupupFKKS={{-((-2 m r y^4+s^2 y^4-2 m r^3 y^4 \[Nu]^2+r^2 y^4 (1+s^2 \[Nu]^2)+r^4 (y^2+s^2 (-1+y^2 \[Nu]^2)))/((-2 m r+r^2+s^2) (s^2-y^2) (r^2+y^2) (1+r^2 \[Nu]^2) (-1+y^2 \[Nu]^2))),(I r^3 \[Nu])/((r^2+y^2) (1+r^2 \[Nu]^2)),-((I y^3 \[Nu])/((r^2+y^2) (-1+y^2 \[Nu]^2))),-((2 m r y^2-s^2 y^2+2 m r^3 y^2 \[Nu]^2-r^4 y^2 \[Nu]^2-r^2 (s^2+y^4 \[Nu]^2))/((-2 m r+r^2+s^2) (s^2-y^2) (r^2+y^2) (1+r^2 \[Nu]^2) (-1+y^2 \[Nu]^2)))},{-((I r^3 \[Nu])/((r^2+y^2) (1+r^2 \[Nu]^2))),(-2 m r+r^2+s^2)/((r^2+y^2) (1+r^2 \[Nu]^2)),0,-((I r \[Nu])/((r^2+y^2) (1+r^2 \[Nu]^2)))},{(I y^3 \[Nu])/((r^2+y^2) (-1+y^2 \[Nu]^2)),0,(-s^2+y^2)/((r^2+y^2) (-1+y^2 \[Nu]^2)),-((I y \[Nu])/((r^2+y^2) (-1+y^2 \[Nu]^2)))},{(-2 m r y^2+s^2 y^2-2 m r^3 y^2 \[Nu]^2+r^4 y^2 \[Nu]^2+r^2 (s^2+y^4 \[Nu]^2))/((-2 m r+r^2+s^2) (s^2-y^2) (r^2+y^2) (1+r^2 \[Nu]^2) (-1+y^2 \[Nu]^2)),(I r \[Nu])/((r^2+y^2) (1+r^2 \[Nu]^2)),(I y \[Nu])/((r^2+y^2) (-1+y^2 \[Nu]^2)),-((-2 m (r+r^3 \[Nu]^2)+(r^2+y^2) (1+r^2 \[Nu]^2+s^2 \[Nu]^2-y^2 \[Nu]^2))/((-2 m r+r^2+s^2) (s^2-y^2) (r^2+y^2) (1+r^2 \[Nu]^2) (-1+y^2 \[Nu]^2)))}}/.{m->1,s->\[Chi]};
ComponentValue[BDolan[b,c]//ToBasis[ch]//ComponentArray,BupupFKKS];

(*Polarization tensor in Boyer-Lindquist coordinates*)
DefTensor[BBL[b,c],M]
BBLtemp1=(\[CapitalLambda][f,-g]\[CapitalLambda][d,-e]BDolan[g,e]//eval)/.{y->\[Chi] Cos[\[Theta][]],s->\[Chi]}//Transpose;
BBLtemp2=(\[CapitalLambda]updownFKKStoBL . BupupFKKS) . Transpose[\[CapitalLambda]updownFKKStoBL]/.{y->\[Chi] Cos[\[Theta][]],s->\[Chi]};
BBLtest=(BBLtemp1-BBLtemp2)//Simplify;
Print["BBL test = "<>ToString[BBLtest]];
Print["Should be 0! If not, then BBL is not correct!"];
ComponentValue[BBL[b,c]//ToBasis[ch]//ComponentArray,BBLtemp2];

(*The Dolan eq. (8)*)
DefScalarFunction[R]
DefScalarFunction[S]
DefConstantSymbol[\[Omega]]
DefConstantSymbol[m]
DefConstantSymbol[\[Nu]r]
DefConstantSymbol[\[Nu]i]
DefConstantSymbol[\[Omega]r]
DefConstantSymbol[\[Omega]i]
Z[t[],r[],\[Theta][],\[Phi][]]=R[r[]]S[\[Theta][]]Exp[I m \[Phi][]]Exp[-I \[Omega] t[]];
DefTensor[A[-b],M]
A/:A[c_]:=Module[{b},BBL[c,b]CD[-b]@Z[t[],r[],\[Theta][],\[Phi][]]];
A/:A[-c_]:=Module[{b,e},met[-e,-c]BBL[e,b]CD[-b]@Z[t[],r[],\[Theta][],\[Phi][]]];

(*Analyitc form of the real part of the field A^\mu. This part takes about 50 min to evaluate, hence, we simply load in the copmuted result, since it is always the same, and comment the following out*)
(*sin=Sin[m \[Phi]-\[Omega]r t];
cos=Cos[m \[Phi]-\[Omega]r t];
complextoreal={\[Nu]\[Rule]r\[Nu]+I i\[Nu],\[Omega]\[Rule]\[Omega]r+I \[Omega]i,R[r[]]\[RuleDelayed] Rr[r[]]+I Ri[r[]],R'[r[]]\[RuleDelayed] Rr'[r[]]+I Ri'[r[]],S[\[Theta][]]\[Rule]Sr[\[Theta][]]+I Si[\[Theta][]],S'[\[Theta][]]\[Rule]Sr'[\[Theta][]]+I Si'[\[Theta][]]};

AupSimptemp=(A[b]//eval)//.complextoreal//.{t[]\[Rule]t,r[]\[Rule]r,\[Theta][]\[Rule]\[Theta],\[Phi][]\[Rule]\[Phi]};
AupSimptemp1=CoefficientList[AupSimptemp[[1]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["t-comp. done"];
AupSimptemp2=CoefficientList[AupSimptemp[[2]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["r-comp. done"];
AupSimptemp3=CoefficientList[AupSimptemp[[3]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["\[Theta]-comp. done"];
AupSimptemp4=CoefficientList[AupSimptemp[[4]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["\[Phi]-comp. done"];

AupSimp={sin AupSimptemp1[[2,1]]+cos AupSimptemp1[[1,2]],sin AupSimptemp2[[2,1]]+cos AupSimptemp2[[1,2]],sin AupSimptemp3[[2,1]]+cos AupSimptemp3[[1,2]],sin AupSimptemp4[[2,1]]+cos AupSimptemp4[[1,2]]};
Export["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/AupSimp.mx",AupSimp];*)
(*---------------------------------------------------------------------------------------------*)
(*Analyitc form of the real part of the field A_\mu. This part takes about 50 min to evaluate, hence, we simply load in the copmuted result, since it is always the same, and comment the following out*)
(*sin=Sin[m \[Phi]-\[Omega]r t];
cos=Cos[m \[Phi]-\[Omega]r t];
complextoreal={\[Nu]\[Rule]r\[Nu]+I i\[Nu],\[Omega]\[Rule]\[Omega]r+I \[Omega]i,R[r[]]\[RuleDelayed] Rr[r[]]+I Ri[r[]],R'[r[]]\[RuleDelayed] Rr'[r[]]+I Ri'[r[]],S[\[Theta][]]\[Rule]Sr[\[Theta][]]+I Si[\[Theta][]],S'[\[Theta][]]\[Rule]Sr'[\[Theta][]]+I Si'[\[Theta][]]};

AdownSimptemp=(A[-b]//eval)//.complextoreal//.{t[]\[Rule]t,r[]\[Rule]r,\[Theta][]\[Rule]\[Theta],\[Phi][]\[Rule]\[Phi]};
AdownSimptemp1=CoefficientList[AdownSimptemp[[1]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["t-comp. done"];
AdownSimptemp2=CoefficientList[AdownSimptemp[[2]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["r-comp. done"];
AdownSimptemp3=CoefficientList[AdownSimptemp[[3]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["\[Theta]-comp. done"];
AdownSimptemp4=CoefficientList[AdownSimptemp[[4]]//Re//ComplexExpand,{sin,cos}]//Simplify; Print["\[Phi]-comp. done"];

AdownSimp={sin AdownSimptemp1[[2,1]]+cos AdownSimptemp1[[1,2]],sin AdownSimptemp2[[2,1]]+cos AdownSimptemp2[[1,2]],sin AdownSimptemp3[[2,1]]+cos AdownSimptemp3[[1,2]],sin AdownSimptemp4[[2,1]]+cos AdownSimptemp4[[1,2]]};
Export["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/AdownSimp.mx",AdownSimp];*)

(*Importing the result of the commented code above*)
coordinates={t[]->t,r[]->r,\[Theta][]->\[Theta],\[Phi][]->\[Phi]};
coordconv={Global`t->t,Global`r->r,Global`\[Theta]->\[Theta],Global`\[Phi]->\[Phi]};
coordrev={t->Private`t,r->Private`r,\[Theta]->Private`\[Theta],\[Phi]->Private`\[Phi]};
AupSimp=Import[Global`$root<>"Expressions/AupSimp.mx"]//.coordinates//.{Global`m->Private`m,Global`\[Omega]r->Private`\[Omega]r,Global`\[Omega]i->\[Omega]i,Global`\[Chi]->Private`\[Chi],Global`r\[Nu]->Private`r\[Nu],Global`i\[Nu]->Private`i\[Nu],Global`Si->Private`Si,Global`Ri->Private`Ri,Global`Sr->Private`Sr,Global`Rr->Private`Rr}//.coordconv;
(*AdownSimp=(met[-b,-c]//eval).AupSimp//.coordinates//.coordconv/.{t->Private`t,r->Private`r,\[Theta]->Private`\[Theta],\[Phi]->Private`\[Phi]};*)
AdownSimp=Import[Global`$root<>"Expressions/AdownSimp.mx"]//.coordinates//.{Global`m->Private`m,Global`\[Omega]r->Private`\[Omega]r,Global`\[Omega]i->\[Omega]i,Global`\[Chi]->Private`\[Chi],Global`r\[Nu]->Private`r\[Nu],Global`i\[Nu]->Private`i\[Nu],Global`Si->Private`Si,Global`Ri->Private`Ri,Global`Sr->Private`Sr,Global`Rr->Private`Rr}//.coordconv;
(*AupSimp=(met[b,c]//eval).AdownSimp//.coordinates//.coordconv/.{t->Private`t,r->Private`r,\[Theta]->Private`\[Theta],\[Phi]->Private`\[Phi]};*)
Print["Imported analytic field modes."]

(*Generating the energy momentum tensor for this minimally coupled massive vector*)
DefTensor[Aplace[-b],M]
DArule=Table[PDch[{i,-ch}]@Aplace[{k,-ch}]->ToExpression["D"<>ToString[i]<>ToString["A"]<>ToString[k]],{i,0,3},{k,0,3}]//Flatten;
AtoAi=Table[{Aplace[{i,-ch}]->ToExpression["A"<>ToString[i]]},{i,0,3}]//Flatten;

DefConstantSymbol[\[Mu]]
DefTensor[FS[-b,-c],M]
DefTensor[T[-b,-c],M]
FS/:FS[-b_,-c_]:=CD[-b]@Aplace[-c]-CD[-c]@Aplace[-b]
(*Tdd=(ChangeCovD[(\[Mu]^2Aplace[-h]Aplace[-c]+met[d,e]FS[-h,-d]FS[-c,-e]-1/4met[-h,-c](FS[-d,-e]FS[-f,-g]met[d,f]met[e,g]+2\[Mu]^2met[e,d]Aplace[-d]Aplace[-e])),CD,PDch]//eval)//.DArule//.AtoAi; (*Tdd stands for T_\mu\nu*)
Export["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/Tdd.mx",Tdd];*)

(*The energy density is the contraction of the timelike Killing field into the T_\mu\nu. Since \xi^\mu=(1,0,0,0)^\mu, we can take simple the T_tt component as the energy density*)
Tdd=Import[Global`$root<>"Expressions/Tdd.mx"]//.coordinates//.{Global`m->Private`m,Global`\[Omega]r->Private`\[Omega]r,Global`\[Omega]i->\[Omega]i,Global`\[Chi]->Private`\[Chi],Global`r\[Nu]->Private`r\[Nu],Global`i\[Nu]->Private`i\[Nu],Global`Si->Private`Si,Global`Ri->Private`Ri,Global`Sr->Private`Sr,Global`Rr->Private`Rr}//.coordconv;
edenrule={Table[ToExpression["Global`D"<>ToString[i]<>ToString["A"]<>ToString[k]]->ToExpression["D"<>ToString[i]<>ToString["A"]<>ToString[k]],{i,0,3},{k,0,3}],Table[ToExpression["Global`A"<>ToString[i]]->ToExpression["A"<>ToString[i]],{i,0,3}]}//Flatten;
edenrule1={Table[ToExpression["Private`D"<>ToString[i]<>ToString["A"]<>ToString[k]]->ToExpression["D"<>ToString[i]<>ToString["A"]<>ToString[k]],{i,0,3},{k,0,3}],Table[ToExpression["Private`A"<>ToString[i]]->ToExpression["A"<>ToString[i]],{i,0,3}]}//Flatten;
Tud=(met[d,e]//eval) . Tdd;
Export[Global`$root<>"Tud.mx",Tud];
Tud=Import[Global`$root<>"Expressions/Tud.mx"]//.coordinates//.{Global`m->Private`m,Global`\[Omega]r->Private`\[Omega]r,Global`\[Omega]i->\[Omega]i,Global`\[Chi]->Private`\[Chi],Global`\[Mu]->Private`\[Mu],Global`r\[Nu]->Private`r\[Nu],Global`i\[Nu]->Private`i\[Nu],Global`Si->Private`Si,Global`Ri->Private`Ri,Global`Sr->Private`Sr,Global`Rr->Private`Rr}//.coordconv;
Tudtemp=Tud[[1,1]]//.edenrule/.edenrule1/.{t->Private`t,r->Private`r,\[Theta]->Private`\[Theta],\[Phi]->Private`\[Phi]};
Print["Imported analytic Tud, Tdd and energy density."];


(*Interpolation code*)
VectorField[Rsol0_InterpolatingFunction,Anglsol0_,\[Chi]0_,m0_,nh0_,r\[Nu]0_,i\[Nu]0_,\[Omega]r0_,\[Mu]0_,MBH1_]:=Module[{Rsol=Rsol0,Anglsol=Anglsol0,\[Chi]1=\[Chi]0,m1=m0,r\[Nu]1=r\[Nu]0,i\[Nu]1=i\[Nu]0,\[Omega]r1=\[Omega]r0,\[Mu]1=\[Mu]0,nh1=nh0,MBH=MBH1},
With[{precValue=TeuInterEnv`$precValue},
prec=SetPrecision[#,precValue]&;
];
parameters=prec@{\[Chi]->\[Chi]0,m->m0,r\[Nu]->r\[Nu]0,i\[Nu]->i\[Nu]0,\[Omega]i->0,\[Omega]r->\[Omega]r0,\[Mu]->\[Mu]0,nh->nh0};

(*The radial range depending on the input parameters of the system*)
Mtilde=1; (*1/(MBH1+1);*) (*we introduce Mtilde here ONLY to fix rstop to be same (or roughly the same) as in the ModesCode! nothing more!!*)
rstop=(3/2)If[m0<3,prec@(4(10(m0+nh0)Mtilde)/(\[Mu]0^2)-2),prec@((400(m0+nh0)Mtilde)/(7\[Mu]0^2)-1)]; (*The -1 makes sure that the range of the Rsol and Sangl is larger than the interpolation of the field modes*)
rplus=1+Sqrt[1-\[Chi]0^2];
\[Epsilon]r=10^-2;
rtst=prec@(\[Epsilon]r+rplus);
rglue=If[nh0==0,0.12*rstop,0.2*rstop];
Print["precision: "<>ToString[TeuInterEnv`$precValue]];
Print["rplus = "<>ToString[rplus]];
Print["rstart = "<>ToString[rtst]];
Print["epsilonr = "<>ToString[(rtst-rplus)/rplus]];
Print["rstop = "<>ToString[rstop]];
Print["r glue = "<>ToString[rglue]];

(*Generating the interpolations for a grid of 40x40 points over (r,\[Theta]) for each of the A^\mu and A_\mu components*)
Print["Interpolation: Start."];
\[Theta]mesh=If[m0==1,99,If[m0==2,59,59]];
rmesh=If[m0==1,199,If[m0==2,699,699]];
p=If[rstop<500,4,5];
Print["rmesh = "<>ToString[rmesh]];
Print["Thetamesh = "<>ToString[\[Theta]mesh]];
rdatapoint[x_]:=(rstop-rtst)/((rmesh+1)^p-1) x^p+rtst-(rstop-rtst)/((rmesh+1)^p-1);
(*rrange=Range[rtst,rstop,(rstop-rtst)/rmesh];*)
rrange=rdatapoint[Range[2,rmesh+2]];

slp=0.25; (*slope of the step function*)
step[x_]:=\[Pi]/(\[Theta]mesh+1) ((\[Theta]mesh+1)^2/2 Exp[slp(x-(\[Theta]mesh+1)/2)])/((\[Theta]mesh+1)+(\[Theta]mesh+1)/2(Exp[slp(x-(\[Theta]mesh+1)/2)]-1));
\[Theta]range=prec@step[Range[1,\[Theta]mesh+1]];
\[Theta]start=prec@(\[Theta]range//First);
\[Theta]stop=prec@(\[Theta]range//Last);
Print["ThetaStart = "<>ToString[\[Theta]start]];
Print["ThetaStop = "<>ToString[\[Theta]stop]];

R\[Theta]ranges=prec@{{rtst,rstop},{\[Theta]start,\[Theta]stop}};

Print["Exponential extension: Start."];
(*Exponential extension of the Radial function*)
Rslope[f1_,x1_,f2_,x2_]:=(f1-f2)/(x1-x2);
(*First matching point*)
r1=2rglue/3;
R1re=Rsol[r1]//Re//Log;
R1im=Rsol[r1]//Im//Log;
dR1re=D[Rsol[r],r]/.{r->r1}//Re//Log;
dR1im=D[Rsol[r],r]/.{r->r1}//Im//Log;
ddR1re=D[Rsol[r],{r,2}]/.{r->r1}//Re//Log;
ddR1im=D[Rsol[r],{r,2}]/.{r->r1}//Im//Log;
(*Second matching point*)
r2=rglue;
R2re=Rsol[r2]//Re//Log;
R2im=Rsol[r2]//Im//Log;
dR2re=D[Rsol[r],r]/.{r->r2}//Re//Log;
dR2im=D[Rsol[r],r]/.{r->r2}//Im//Log;
ddR2re=D[Rsol[r],{r,2}]/.{r->r2}//Re//Log;
ddR2im=D[Rsol[r],{r,2}]/.{r->r2}//Im//Log;
(*Slopes between matching points*)
RslopeRe=Rslope[R1re,r1,R2re,r2];
RslopeIm=Rslope[R1im,r1,R2im,r2];
dRslopeRe=Rslope[dR1re,r1,dR2re,r2];
dRslopeIm=Rslope[dR1im,r1,dR2im,r2];
ddRslopeRe=Rslope[ddR1re,r1,ddR2re,r2];
ddRslopeIm=Rslope[ddR1im,r1,ddR2im,r2];

(*Final solution*)
Rad[rr_]:=Piecewise[{{Rsol[y]//.{y->rr},rr<rglue},{Exp[RslopeRe rr]/Exp[RslopeRe rglue]Re[Rsol[rglue]]+I Exp[RslopeIm rr]/Exp[RslopeIm rglue]Im[Rsol[rglue]],rr>rglue}}];
dRad[rr_]:=Piecewise[{{D[Rsol[y],y]//.{y->rr},rr<rglue},{Exp[dRslopeRe rr]/Exp[dRslopeRe rglue]Re[D[Rsol[x],x]/.{x->rglue}]+I Exp[dRslopeIm rr]/Exp[dRslopeIm rglue]Im[D[Rsol[x],x]/.{x->rglue}],rr>rglue}}];
ddRad[rr_]:=Piecewise[{{D[Rsol[y],{y,2}]//.{y->rr},rr<rglue},{Exp[ddRslopeRe rr]/Exp[ddRslopeRe rglue]Re[D[Rsol[x],{x,2}]/.{x->rglue}]+I Exp[ddRslopeIm rr]/Exp[ddRslopeIm rglue]Im[D[Rsol[x],{x,2}]/.{x->rglue}],rr>rglue}}];
(*Angular solution*)
Sangl[\[Theta]\[Theta]_]:=Anglsol/.{Global`\[Theta]->Private`\[Theta],TeuInterEnv`\[Theta]->Private`\[Theta]}//.coordinates//.{\[Theta]->\[Theta]\[Theta]};
dSangltemp=D[Sangl[\[Theta]],\[Theta]];
dSangl[\[Theta]\[Theta]1_]:=dSangltemp//.{\[Theta]->\[Theta]\[Theta]1};
ddSangltemp=D[dSangl[\[Theta]],\[Theta]];
ddSangl[\[Theta]\[Theta]2_]:=ddSangltemp//.{\[Theta]->\[Theta]\[Theta]2};
With[{Rad=Rad, dRad = dRad, ddRad = ddRad},
solidentify={Ri[r_]:>Im[Rad[r]],Rr[r_]:>Re[Rad[r]],Ri'[r_]:>Im[dRad[r]],Rr'[r_]:>Re[dRad[r]],Rr''[r_]:>Re[ddRad[r]],Ri''[r_]:>Im[ddRad[r]],Si[\[Theta]_]:>Im[Sangl[\[Theta]]],Sr[\[Theta]_]:>Re[Sangl[\[Theta]]],Si'[\[Theta]_]:>Im[dSangl[\[Theta]]],Sr'[\[Theta]_]:>Re[dSangl[\[Theta]]],Sr''[\[Theta]_]:>Re[ddSangl[\[Theta]]],Si''[\[Theta]_]:>Im[ddSangl[\[Theta]]]};
];
Print["Solution at horizon: Test:"];
Print["Rad[rtst] = "<>ToString[Rad[rtst]]];
Print["Exponential extension: Done!"];

AdownSimpExpl={0,0,0,0};
AupSimpExpl={0,0,0,0};
dataout={0,0,0,0};

sin=Sin[m \[Phi]-t \[Omega]r];
cos=Cos[m \[Phi]-t \[Omega]r];
Print["Generating Interpolating functions for \!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\)"];
Do[
Clear[coefflistdown,coefflistup,coefflistInterdown,coefflistInterup];
coefflistdownExpl=ConstantArray[0,{rmesh+1,\[Theta]mesh+1,2}];
coefflistupExpl=ConstantArray[0,{rmesh+1,\[Theta]mesh+1,2}];
coefflistdown=DeleteCases[CoefficientList[AdownSimp[[i]]/.{\[Omega]i->0},{sin,cos}]/.{\[Nu]i->i\[Nu],\[Nu]r->r\[Nu]}/.{\[Chi]->\[Chi]0,m->m0,r\[Nu]->r\[Nu]0,i\[Nu]->i\[Nu]0,\[Omega]r->\[Omega]r0,\[Mu]->\[Mu]0,nh->nh0},0,2]//Flatten;
coefflistup=DeleteCases[CoefficientList[AupSimp[[i]]/.{\[Omega]i->0},{sin,cos}]/.{\[Nu]i->i\[Nu],\[Nu]r->r\[Nu]}/.{\[Chi]->\[Chi]0,m->m0,r\[Nu]->r\[Nu]0,i\[Nu]->i\[Nu]0,\[Omega]r->\[Omega]r0,\[Mu]->\[Mu]0,nh->nh0},0,2]//Flatten;
coefflistdown1OE = Global`OptimizedFunction[{r,\[Theta]}, coefflistdown[[1]]//.solidentify];
coefflistup1OE = Global`OptimizedFunction[{r,\[Theta]}, coefflistup[[1]]//.solidentify];
coefflistdown2OE = Global`OptimizedFunction[{r,\[Theta]}, coefflistdown[[2]]//.solidentify];
coefflistup2OE = Global`OptimizedFunction[{r,\[Theta]}, coefflistup[[2]]//.solidentify];
testval = coefflistdown1OE[2,1];
If[Not[NumericQ[testval]], Print["ERROR: coefficients not numbers after insertion of values"];Abort[];];
DistributeDefinitions[prec,precValue,coefflistdown1OE,coefflistup1OE, solidentify, rrange, \[Theta]range, \[Theta]mesh, rmesh, coefflistdown2OE,coefflistup2OE];
SetSharedVariable[coefflistdownExpl,coefflistupExpl,coefflistdownExpl,coefflistupExpl];
	ParallelDo[
	If[Not[NumericQ[coefflistdown1OE[2,1]]], Print["ERROR: Not a number in vector field norm."]];
	Clear[rinput,\[Theta]input];
	rinput=prec@rrange[[ir]];
	\[Theta]input=prec@\[Theta]range[[i\[Theta]]];
	coefflistdownExpl[[ir,i\[Theta],1]]=prec@(coefflistdown1OE[rinput,\[Theta]input]); (*cos part*)
	coefflistupExpl[[ir,i\[Theta],1]]=prec@(coefflistup1OE[rinput,\[Theta]input]); (*cos part*)
	coefflistdownExpl[[ir,i\[Theta],2]]=prec@(coefflistdown2OE[rinput,\[Theta]input]); (*sin part*)
	coefflistupExpl[[ir,i\[Theta],2]]=prec@(coefflistup2OE[rinput,\[Theta]input]); (*sin part*)
	,{ir,1,rmesh+1},{i\[Theta],1,\[Theta]mesh+1}];
Print["Data produced for: "<>ToString[i]];
dataout[[i]]=Table[{{rrange[[k]],\[Theta]range[[j]]},coefflistdownExpl[[k,j,1]]},{k,1,rmesh+1},{j,1,\[Theta]mesh+1}]; (*only the cosine part*)
coefflistInterdown=Table[Interpolation[Flatten[Table[{{rrange[[k]],\[Theta]range[[j]]},coefflistdownExpl[[k,j,l]]},{k,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->3,Method->"Spline"],{l,1,2}];
coefflistInterup=Table[Interpolation[Flatten[Table[{{rrange[[k]],\[Theta]range[[j]]},coefflistupExpl[[k,j,l]]},{k,1,rmesh+1},{j,1,\[Theta]mesh+1}],1],InterpolationOrder->3,Method->"Spline"],{l,1,2}];
AdownSimpExpl[[i]]=(cos  coefflistInterdown[[1]][r,\[Theta]]+sin coefflistInterdown[[2]][r,\[Theta]]);
AupSimpExpl[[i]]=(cos  coefflistInterup[[1]][r,\[Theta]]+sin coefflistInterup[[2]][r,\[Theta]]);
,{i,1,4}];
Print["Interpolation: Done."];

(*Final output of the the interpolating of the fields*)
Adownvector[tt1_,rr1_,\[Theta]\[Theta]1_,\[Phi]\[Phi]1_]:=AdownSimpExpl//.parameters//.{t->tt1,r->rr1,\[Theta]->\[Theta]\[Theta]1,\[Phi]->\[Phi]\[Phi]1}//Evaluate;
Aupvector[tt2_,rr2_,\[Theta]\[Theta]2_,\[Phi]\[Phi]2_]:=AupSimpExpl//.parameters//.{t->tt2,r->rr2,\[Theta]->\[Theta]\[Theta]2,\[Phi]->\[Phi]\[Phi]2}//Evaluate;

(*Testing the vector result*)
atest=Adownvector[0,2,2,2];
Print["Adownvector[0,2,2,2] = "<>ToString[atest]];

(*Integrating over the spatial 3-slice outside of the horizon for the complete mass*)
Aexpl=Table[{ToExpression["Private`A"<>ToString[i]]->Adownvector[t,r,\[Theta],\[Phi]][[i+1]]},{i,0,3}]//Flatten;
dAexpl=Table[{ToExpression["Private`D"<>ToString[i]<>ToString["A"]<>ToString[k]]->D[Adownvector[t,r,\[Theta],\[Phi]][[k+1]],{t,r,\[Theta],\[Phi]}[[i+1]]]},{i,0,3},{k,0,3}]//Flatten//Simplify;

Print["Massnormalization: Start."];
(*Mass integral*)
prec=SetPrecision[#,30]&;
TuptdowntTEMP=prec@(Tudtemp/.{\[Omega]i->0}//.parameters)//.Aexpl//.dAexpl//.coordinates;
Tuptdownt[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=TuptdowntTEMP//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
Print["Tuptdownt test = "<>ToString[prec@Tuptdownt[0,2,2,2]]];
sqrtmg[r_,\[Theta]_]:=prec@(Sin[\[Theta]](r^2+\[Chi]^2 Cos[\[Theta]]^2)//.parameters);
Print["sqrtmg test = "<>ToString[sqrtmg[2,2]]];
Print["rplus = "<>ToString[rplus]];
rintstart=prec@(rplus+10^-1.9);
rintstop=prec@Last[rrange];
Print["rintstart = "<>ToString[rintstart]];
Print["epsilonrint = "<>ToString[(rintstart-rplus)/rplus]];
Print["rintstop = "<>ToString[rintstop]];
Print["thetastart = "<>ToString[\[Theta]start]];
Print["thetastop = "<>ToString[\[Theta]stop]];
energydensity = {r,\[Theta],\[Phi]}|->-Tuptdownt[0,r,\[Theta],\[Phi]];
massNorm=NIntegrate[-Tuptdownt[0,r,\[Theta],\[Phi]]sqrtmg[r,\[Theta]],{r,rintstart,rintstop},{\[Theta],\[Theta]start,\[Theta]stop},{\[Phi],0,2\[Pi]},MaxRecursion->200,MinRecursion->10,WorkingPrecision->15];
massNormint[rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=-Tuptdownt[0,r,\[Theta],\[Phi]]sqrtmg[r,\[Theta]]/.{r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
Print["massNorm obtained: "<>ToString[massNorm]];

(*Normalization to the mass in terms of the inital BH mass*)
Anorm=SetPrecision[Sqrt[MBH1/massNorm],20];
Print["Anorm obtained: "<>ToString[Anorm]];
Print["Massnormalization: Done."];

(*Final output of the the interpolating of the fields and the final output of the package*)
edensityNorm[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=-Anorm^2 Tuptdownt[t,r,\[Theta],\[Phi]] Sqrt[(-met[{0,-ch},{0,-ch}]//eval)/.parameters//.{t[]->t,r[]->r,\[Theta][]->\[Theta],\[Phi][]->\[Phi]}]//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
(*Note the sign and the additional lapse factor in front of the meassure of the intergral for the nergy density. Compare with 1705.01544v2 eq. (10): Sqrt[-g]=NSqrt[\[Gamma]]*)
AdownNorm[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Anorm AdownSimpExpl//.coordinates/.parameters//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
AupNorm[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Anorm AupSimpExpl//.coordinates/.parameters//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
]

Print["VectorFieldNorm Functions:"];
Print["VectorField[Rsol,Angsol,spin,Field-m-mode,\[Nu]r,\[Nu]i,Re(\[Omega]),Procamass,overtone,Normalization-mass]"];
Print["FieldEnergyMomentum[AdownNorm[t,r,\[Theta],\[Phi]],Procamass,spin]"];
Print["---------------------------------------------------------------------"];
End[]

EndPackage[]



