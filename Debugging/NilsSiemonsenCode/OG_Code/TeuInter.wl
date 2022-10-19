(* ::Package:: *)

BeginPackage["TeuInter`"]

TeuCode::usage "something."

Begin["TeuInterEnv`"]

(*All of the input has to be fiven in M0 units*)
TeuCode[min_,nin_,spininitial_,procamass_,massinitial_]:=Module[{min1=min,nin1=nin,spin=spininitial,massinitial1=massinitial},
prec=SetPrecision[#,20]&;

m=min;
n=nin;
spinin=prec@spin;
Mini=prec@massinitial;
\[Mu]=prec@procamass;

spinstring=NumberForm[spinin*10^6//Round,6,DigitBlock->5,ExponentStep->6,NumberSeparator->""]
modedataoutput=prec@Import[ToString["./precision_mode_library/m"<>ToString[m]<>ToString["n"]<>ToString[n]<>ToString["_a"]<>ToString[spinstring]<>ToString["_Sm1_prec.mx"]]]//.{QNMcode`r\[Nu]->r\[Nu],QNMcode`i\[Nu]->i\[Nu],QNMcode`\[Theta]->\[Theta],QNMcode`R->R}//.{QNMmode`r\[Nu]->r\[Nu],QNMmode`i\[Nu]->i\[Nu],QNMmode`\[Theta]->\[Theta],QNMmode`R->R};

(*This is the solution on an M0 backgroud, so we don't need to convert any units!*)
pos=Position[modedataoutput[[3,All]],Nearest[modedataoutput[[3,All]],\[Mu]]//First][[1,1]];
Rsol=R/.modedataoutput[[6,pos]];
Anglsol=modedataoutput[[2,pos]];
Procamass=modedataoutput[[3,pos]];
wr=modedataoutput[[5,pos]]//Re;
wi=modedataoutput[[5,pos]]//Im;
\[Nu]r=r\[Nu]/.modedataoutput[[4,pos]];
\[Nu]i=i\[Nu]/.modedataoutput[[4,pos]];

(*Plotting things*)
out1=LogLogPlot[Interpolation[Table[{modedataoutput[[3,i]],modedataoutput[[5,i]]//Im},{i,1,pos}],InterpolationOrder->2][\[Mu]in],{\[Mu]in,0.05,Procamass},PlotRange->All];

(*We set the Procamass to the BH mass M0, since we cancel that out later on anyway*)
prec@KerrProcaEvolution[m,n,spinin,\[Mu],Mini];
out2=Mfspinf;
EprocainM0=1-Mfspinf[[1]];
EprocainMp=1/Mfspinf[[1]]-1;
VectorField[Rsol,Anglsol,spinin,1,0,\[Nu]r,\[Nu]i,wr,\[Mu],EprocainM0];
rstart=prec@rtst;
rstop=prec@(rstop/3);

(*Plotting things*)
out3=GraphicsRow[{LogPlot[Rsol[r]//Re//Abs,{r,rstart,rstop},PlotRange->All],LogPlot[Rsol[r]//Im//Abs,{r,rstart,rstop},PlotRange->All],Plot[AupNorm[0,5,hh,1] . AdownNorm[0,5,hh,1],{hh,\[Theta]start,\[Theta]stop}]}];

(*energy density*)
phimax=\[Phi]/.FindMaximum[edensityNorm[0,r,\[Pi]/2,\[Phi]],{{r,10,rstart,20},{\[Phi],\[Pi],0,2\[Pi]}}][[2]];
rmesh=79;
rdatapoint[x_]:=(rstop-rstart)/((rmesh+1)^4-1) x^4+rstart-(rstop-rstart)/((rmesh+1)^4-1);
rrange1=Delete[rdatapoint[Range[1,rmesh+1]],1];
edenseplotdata2=Table[{rr,edensityNorm[0,rr,\[Pi]/2,phimax]//Re},{rr,rrange1}];
edensefunc2=Interpolation[edenseplotdata2,InterpolationOrder->2,Method->"Spline"];
plo1=Show[LogPlot[ edensefunc2[1/Procamass^2]^-1 edensefunc2[r/Procamass^2],{r,(rstart+0.1) Procamass^2,rstop Procamass^2}],LogPlot[{Exp[-2r],10^-3},{r,(rstart) Procamass^2,rstop Procamass^2},PlotStyle->Dashed]];
plo2=Plot[edensefunc2[r],{r,(rstart),rstop},PlotRange->All];
plo3=LogPlot[edensefunc2[r],{r,(rstart),rstop},PlotRange->{{Automatic,Automatic},{10^-3,10^-10}}];
out4=GraphicsRow[{plo1,plo2,plo3}]

(*Mass test*)
Aexpl1=Table[{ToExpression["A"<>ToString[i]]->AdownNorm[t,r,\[Theta],\[Phi]][[i+1]]},{i,0,3}]//Flatten;
dAexpl1=Table[{ToExpression["D"<>ToString[i]<>ToString["A"]<>ToString[k]]->D[AdownNorm[t,r,\[Theta],\[Phi]][[k+1]],{t,r,\[Theta],\[Phi]}[[i+1]]]},{i,0,3},{k,0,3}]//Flatten//Simplify;
Tud=Import["/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/Tud.mx"]//.{Private`A0->A0,Private`A1->A1,Private`A2->A2,Private`A3->A3,Private`\[Mu]->\[Mu],Private`\[Chi]->\[Chi],Private`r->r,Private`\[Theta]->\[Theta],Private`D1A2->D1A2,Private`D2A1->D2A1,Private`D0A3->D0A3,Private`D0A2->D0A2,Private`D0A1->D0A1,Private`D1A3->D1A3,Private`D1A0->D1A0,Private`D2A3->D2A3,Private`D2A0->D2A0,Private`D3A0->D3A0,Private`D3A1->D3A1,Private`D3A2->D3A2}//.{r[]->r,\[Theta][]->\[Theta]};
Tudtemp=Tud[[1,1]];
Tuptdownt[tt_,rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=SetPrecision[Tudtemp/.{\[Omega]i->0}//.{\[Chi]->spinin},20]//.Aexpl1//.dAexpl1//.{t->tt,r->rr,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
sqrtmg[r_,\[Theta]_]:=Sin[\[Theta]](r^2+spinin^2 Cos[\[Theta]]^2);
rintstart=rstart;
rintstop=rstop;
out5=N[NIntegrate[-Tuptdownt[0,r,\[Theta],\[Phi]]sqrtmg[r,\[Theta]],{r,rintstart,rintstop},{\[Theta],0,\[Pi]},{\[Phi],0,2\[Pi]},MaxRecursion->100,WorkingPrecision->10],20];

(*Teukolsky T*)
FieldEnergyMomentum[AdownNorm[t,r,\[Theta],\[Phi]],\[Mu],spinin];
TeukolskyTlmw[Tdowndown[t,r,\[Theta],\[Phi]],\[Mu],spinin,2,2wr,rstop,rstart,\[Theta]stop,\[Theta]start];
out6=GraphicsRow[Flatten[{Table[Plot[Tlm\[Omega][i][r]//ReIm,{r,rstart,rstop},PlotRange->All],{i,2,4}],Table[Plot[Conjugate[Tlm\[Omega][i][r]]//ReIm,{r,rstart,rstop},PlotRange->All],{i,2,4}]}]];

(*MST solution*)
Monitor[Do[
Print["--------------------------------------"];
MSTgRoots[spinin,-2,j,l,2];
\[Nu]sol[w_]:=nurealMY[w]+I nuimaginaryMY[w];
plota=Plot[{nurealMY[w],nuimaginaryMY[w]},{w,0,2}];
plotb=Plot[AlmSWSH[w],{w,0,2}];
plotrenanglmom[l]=GraphicsRow[{plota,plotb}];
Print[\[Nu]sol[2wr]];
Clear[\[Nu]sol];
\[Nu]sol[w_]:=nufullBHPT[w];
Print[\[Nu]sol[2wr]];
Print[AlmSWSH[2wr]];
MSTSolution[spinin,\[Nu]sol[2wr],2wr,AlmSWSH[2wr],-2,j,l];
RinSet[l,j]=SetPrecision[Rnumerical[r],20];
RinSet[l,-j]=Conjugate[Rnumerical[r]];
RinSetMST[l,j]=RinMST[r];
BincSet[l,j]=Binc;
BrefSet[l,j]=Bref;
BincSet[l,-j]=Conjugate[Binc];
AinSetMST[l,j]=AinMST (2wr)^(-1-I 4 wr);
,{l,{2,3,4}},{j,{2}}],l];

out7=GraphicsRow[{Show[Table[LogPlot[{Re[(Tlm\[Omega][l][r]RinSet[l,mm])/(r^2-2r+spinin^2)^2]//Abs,Im[(Tlm\[Omega][l][r]RinSet[l,mm])/(r^2-2r+spinin^2)^2]//Abs},{r,rstart+1,rstop},PlotRange->All,PlotStyle->RGBColor[(l-2)/2,0,0]],{l,{2,3,4}},{mm,{2}}]//Flatten],Table[Plot[{Re[(Tlm\[Omega][l][r]RinSet[l,mm])/(r^2-2r+spinin^2)^2],Im[(Tlm\[Omega][l][r]RinSet[l,mm])/(r^2-2r+spinin^2)^2]},{r,rstart+1,rstop},PlotRange->All,PlotStyle->RGBColor[(l-2)/2,0,0]],{l,{2,3,4}},{mm,{2}}]//Flatten,Show[Table[LogPlot[{Re[RinSet[l,mm]]//Abs,Im[RinSet[l,mm]]//Abs},{r,rstart+1,rstop},PlotRange->All,PlotStyle->RGBColor[(l-2)/2,0,0]],{l,{2,3,4}},{mm,{2}}]//Flatten]}];
out8={Table[(RinSet[i,2]/.{r->40})/BrefSet[i,2]//Abs,{i,2,4}],Table[(RinSet[i,2]/.{r->40})/BincSet[i,2]//Abs,{i,2,4}]};

(*Gravitational waves*)
Monitor[Do[
Print["Mode: "<>ToString[{l,j}]];
GravitationalRadiation[BincSet[l,j],2wr,Tlm\[Omega][l][r],RinSet[l,j],spinin,j,l,rstart+1,rstop];
RePsi4[l,j]=Re[Psi4Inf[t,r,\[Theta],\[Phi]]]//ComplexExpand//TrigReduce//Simplify;
RePsi4[l,-j]=Re[Psi4mInf[t,r,\[Theta],\[Phi]]]//ComplexExpand//TrigReduce//Simplify;
PolWavelm[l,j]=PolWaveForm[t,rstar,\[Theta],\[Phi]];
PolWavelm[l,-j]=PolWaveFormm[t,rstar,\[Theta],\[Phi]];
Powerlm[l,j]=GWPower;
,{l,{2,3}},{j,{2}}],{l,j}];

(*The factor of 2 comes from including also the m<0-modes*)
TotalPower=2Sum[Powerlm[l,2],{l,{2,3}}]; 
Psi4[tt_,rstar_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Sum[RePsi4[l,mm],{l,{2,3,4}},{mm,{-2,2}}]//.{t->tt,r->rstar,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};
Polwave[tt_,rstar_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Sum[PolWavelm[l,mm],{l,{2,3,4}},{mm,{-2,2}}]//.{t->tt,r->rstar,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]};

(*Print["Generating GW modes..."];
Do[
thetarange=prec@Range[0,\[Pi],\[Pi]/400];SWSHdata=Table[{x,Sqrt[2\[Pi]]SpinWeightedSpheroidalHarmonicS[-2,lGW,2,0,x,0]},{x,thetarange}];Swsh[\[Theta]\[Theta]_]:=Interpolation[SWSHdata,InterpolationOrder\[Rule]2,Method->"Hermite"][\[Theta]\[Theta]];
hlm[lGW,2]=2\[Pi] NIntegrate[Sin[\[Theta]]Swsh[\[Theta]]Sum[PolWavelm[lGW,2]/.{t\[Rule]0,rstar\[Rule]0,\[Phi]\[Rule]0},{lGW,{2,3,4}}],{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]10]Exp[I 2wr Tret]
,{lGW,{2,3,4,5,6}}]*)
out9={TotalPower,TotalPower EprocainM0^-2,Psi4[t,rstar,\[Theta],\[Phi]],Polwave[t,rstar,\[Theta],\[Phi]]};

(*Generating output and saving*)
output={out1,out2,out3,out4,out5,out6,out7,out8,out9,out10};
Export[ToString["./GWoutput/GW_m"<>ToString[m]<>ToString["n"]<>ToString[n]<>ToString["_a"]<>ToString[spinstring]<>ToString["_mu"]<>ToString[Round[100 \[Mu]]]<>ToString[".mx"]],output];
];

End[]

EndPackage[]
