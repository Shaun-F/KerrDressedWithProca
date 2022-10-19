(* ::Package:: *)

BeginPackage["WeylPsi4`"]

GravitationalRadiation::usage = "Zlm\[Omega]Inf[Binc,\[Omega],Tlmw,RinIn,\[Chi],rstart,rstop] :> Gives amplitude of Rlmw at infinity."
Psi4Inf
Psi4mInf
GWmodes
Zlm\[Omega]Inf
hlm
PolWaveForm
PolWaveFormm
GWPower
plotsh
plotphiint
plotslm
Slmw

Begin["PrivateWeylPsi4`"]
(*Needs["SpinWeightedSpheroidalHarmonics`","./SWSHs.m"]*)
Needs["SpinWeightedSpheroidalHarmonics`","./SpinWeightedSpheroidalHarmonics.m"]

GravitationalRadiation[BincIn_,\[Omega]In_,TlmwIn_,RinIn_,\[Chi]_,mTeu_,lTeu_,rstart_,rstop_]:=Module[{BincIn1=BincIn,\[Omega]In1=\[Omega]In,TlmwIn1=TlmwIn,RinIn1=RinIn,\[Chi]1=\[Chi],rstart1=rstart,rstop1=rstop},
If[mTeu<0,Print["Error: m<0; use symmetry relations!"];Abort[];,{}];
prec=SetPrecision[#,20]&;
integrand[rr_]:=((TlmwIn RinIn)/((rr^2-2rr+\[Chi]^2)^2))/.{Global`r->rr};
Print[rstop];
Print[rstart];
Zlm\[Omega]Inf=1/(2 I \[Omega]In BincIn) NIntegrate[integrand[rr],{rr,rstart,rstop},MaxRecursion->50,WorkingPrecision->10];
Clear[Slmw];

thetarange=prec@Range[10^-4,\[Pi]-10^-4,(\[Pi]-2 10^-4)/400];
SWSHdata=Table[{x,Sqrt[2\[Pi]]SpinWeightedSpheroidalHarmonicS[-2,lTeu,mTeu,\[Chi] \[Omega]In,x,0]},{x,thetarange}];
Slmw[\[Theta]\[Theta]_]:=Interpolation[SWSHdata,InterpolationOrder->2,Method->"Hermite"][\[Theta]\[Theta]];

alm=SpinWeightedSpheroidalEigenvalue[-2,lTeu,mTeu,\[Chi] \[Omega]In]+2mTeu \[Chi] \[Omega]In-(\[Chi] \[Omega]In)^2;
Print[ToString["Eigenvalues of SWSH: "<>ToString[alm]]];
Psi4Inf[tt_,rstar_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Zlm\[Omega]Inf/Sqrt[2 \[Pi]] Slmw[\[Theta]] Exp[I \[Omega]In (rstar1-t)+I mTeu \[Phi]]/.{t->tt,rstar1->rstar,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]}//ExpToTrig;
Psi4mInf[tt_,rstar_,\[Theta]\[Theta]_,\[Phi]\[Phi]_]:=Conjugate[Zlm\[Omega]Inf]/Sqrt[2 \[Pi]] Slmw[\[Pi]-\[Theta]] Exp[-I \[Omega]In (rstar1-t)-I mTeu \[Phi]]/.{t->tt,rstar1->rstar,\[Theta]->\[Theta]\[Theta],\[Phi]->\[Phi]\[Phi]}//ExpToTrig;
PolWaveForm[tt_,rstar1_,\[Theta]_,\[Phi]_]:=-2 Zlm\[Omega]Inf/\[Omega]In^2 Slmw[\[Theta]]/Sqrt[2\[Pi]]Exp[I \[Omega]In (rstar1-tt)+I mTeu \[Phi]];
PolWaveFormm[tt_,rstar1_,\[Theta]_,\[Phi]_]:=-2 Conjugate[Zlm\[Omega]Inf]/\[Omega]In^2 Slmw[\[Pi]-\[Theta]]/Sqrt[2\[Pi]]Exp[-I \[Omega]In (rstar1-tt)-I mTeu \[Phi]];
Print["Psi4(r->Inf): Done."];

(*Radiated power and luminosity*)
GWPower=Abs[Zlm\[Omega]Inf]^2/(4\[Pi] \[Omega]In^2);
Print["Total GW energy flux: Done."];
]
Print["-------------------------------------------------"];
Print["WeylPsi4 Functions:"];
Print["GravitationalRadiation[Binc,\[Omega](Teu),Tlmw,Rin,\[Chi],m(Teu),l(Teu),rstart,rstop]:> Psi4Inf[t,rstar,\[Theta],\[Phi]], GWPower"];
Print["GWmodes[l(GW),m(GW)] :> hlm[t,rstar]"];
Print["-------------------------------------------------"];

End[]

EndPackage[]
