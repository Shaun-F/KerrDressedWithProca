(* ::Package:: *)

BeginPackage["SpinWeightedSpheroidalHarmonics`"]

SWSHrun::usage = "SWSHrun[m,l,s,\[Chi],\[Omega]] :> Outputs SWSH and Alm for the given parameters."
SWSH::usage = "SWSH[\[Theta]] is the spin weighted spheroidal harmonic for the provided indices"
Alm::usage = "Eigenvalue of the spin weighted spheroidal harmonic for the provided indices"
leaversave

Begin["PrivateSWSH`"]

SWSHrun[min_,lin_,sin_,spin_,win_]:=Module[{m=min,l=lin,\[Chi]=spin,wint=win,sin1=sin},
(*Reliability regime of this code*)
Clear[winit,wang,br,ar,rminus,rplus,k1,k2,Almtemp,alm,Alm,SWSH];
If[(lin-Abs[min])>=3&&spin>0,Print["Error: SWSH-code not optimized for this (l-|m|)-value!"]; Abort[],{}];
If[spin>0&&spin<0.5,Print["Error: SWSH-code not optimized for this \[Chi]-value!"];Abort[],{}];
If[lin>10&&spin==0,Print["Error: SWSH-code not optimized for this l-value!"]; Abort[],{}];
(*precision function*)
prec=SetPrecision[#,20]&;

(*The parameters of interest*)
winit=prec@2*win;
ar=prec@spin/2;(*Additional factor of 2 due to converting from units of M to units of M/2*)

(*Kerr metric quantities with M=1/2*)
br=prec@Sqrt[1-4*ar^2];
rminus=prec@(1-br)/2;
rplus=prec@(1+br)/2;
s=prec@sin;
k1=prec@1/2*Abs[m-s];
k2=prec@1/2*Abs[m+s];
Almtemp=prec@(lin*(lin+1)-s*(s+1));
Print["Alm-in value: "<>ToString[Almtemp//N]];

(*Leaver continued fraction*)
Nmax=100;
If[m>=0,wdensity=200;
wrange=prec@Delete[Range[0,winit,winit/wdensity],1];
Monitor[Do[
wang=prec@wrange[[i]];
\[Gamma]ang[n_]:=prec@(2*ar*wang*(n+k1+k2+s));
\[Beta]ang[n_]:=prec@(n*(n-1)+2*n*(k1+k2+1-2*ar*wang)-(2*ar*wang*(2*k1+s+1)-(k1+k2)*(k1+k2+1))-(ar^2*wang^2+s*(s+1)+alm));
\[Alpha]ang[n_]:=prec@(-2*(n+1)*(n+2*k1+1));
leavercontfrac=prec@(\[Beta]ang[0]+ContinuedFractionK[-\[Alpha]ang[n-1]\[Gamma]ang[n],\[Beta]ang[n],{n,1,Nmax}]);
leaversave[i]=prec@leavercontfrac;
Almsol=alm/.FindRoot[leavercontfrac==0,{alm,Almtemp},WorkingPrecision->MachinePrecision,AccuracyGoal->10,PrecisionGoal->6];
Almtemp=prec@Almsol;
Almout[i]=prec@Almtemp;
,{i,1,Length[wrange]}],i];
(*Output eigenvalue*)
Alm=Almout[Length[wrange]];
,
wang=prec@winit;
\[Gamma]ang[n_]:=prec@2*ar*wang*(n+k1+k2+s);
\[Beta]ang[n_]:=prec@n*(n-1)+2*n*(k1+k2+1-2*ar*wang)-(2*ar*wang*(2*k1+s+1)-(k1+k2)*(k1+k2+1))-(ar^2*wang^2+s*(s+1)+alm);
\[Alpha]ang[n_]:=prec@-2*(n+1)*(n+2*k1+1);
leavercontfrac=prec@(\[Beta]ang[0]+ContinuedFractionK[-\[Alpha]ang[n-1]\[Gamma]ang[n],\[Beta]ang[n],{n,1,Nmax}]);
leaversave[1]=prec@leavercontfrac;
Almsol=alm/.FindRoot[leavercontfrac==0,{alm,Almtemp+If[l==4,0.5,1]},WorkingPrecision->MachinePrecision,AccuracyGoal->10,PrecisionGoal->6];
Alm=prec@Almsol;
];

(*no-spin limit eigenvalues*)
If[ar==0,Alm=lin*(lin+1)-s*(s+1),If[win==0,Alm=lin*(lin+1)-s*(s+1),{}]];
Print["Alm-out value: "<>ToString[Alm//N]];

(*The angular solution*)
Sep=prec@Alm;
Nmax=25;
w=prec@winit;
sima=prec@(w*rplus-ar*m)/br;
ae[0]=prec@1;
ae[1]=prec@(-\[Beta]ang[0]/\[Alpha]ang[0]);
i=1;
While[i<Nmax,ae[i+1]=(-\[Beta]ang[i]*ae[i]-\[Gamma]ang[i]*ae[i-1])/(\[Alpha]ang[i]);i=i+1];
Sang=prec@(Exp[ar*w*u]*(1+u)^(1/2*Abs[m-s])*(1-u)^(1/2*Abs[m+s])*Sum[ae[i]*(1+u)^(i),{i,0,Nmax-1}]/.u->u/.{alm->Alm});
Normalization=prec@(1/Sqrt[NIntegrate[Sang*Conjugate[Sang],{u,-1,1},MaxRecursion->12,PrecisionGoal->10,AccuracyGoal->10]]);
SangNormal=prec@Normalization*Sang;
SWSH[\[Theta]_]:=prec@(SangNormal/.{u->Cos[\[Theta]]});

(*Rootcount test (for win>0)*)
If[ar>0,
If[m>0&&win>0&&Alm>lin*(lin+1)-s*(s+1), Print["SWSH Eigenvalue Error 1!"],{}];
If[m<0&&win>0&&Alm<lin*(lin+1)-s*(s+1),
Print["SWSH Eigenvalue Error 2!"];
If[lin<=8,Print["The Rootcount of the resulting SWSH is: "<>ToString[CountRoots[SWSH[x],{x,0.1,\[Pi]-0.1}]]<>"  (It should be: "<>ToString[lin-Abs[min]]<>")"];
,Print["Check the Sbar[\[Theta]]-plot to count roots! There should be: "<>ToString[lin-Abs[min]]];
];
,{}];
,{}];

(*Rootcount test (for win<0)*)
If[ar>0,
If[m>0&&win<0&&Alm<lin*(lin+1)-s*(s+1), Print["SWSH Eigenvalue Error 1!"];
If[lin<=8,Print["The Rootcount of the resulting SWSH is: "<>ToString[CountRoots[SWSH[x],{x,0.1,\[Pi]-0.1}]]<>"  (It should be: "<>ToString[lin-Abs[min]]<>")"];
,Print["Check the Sbar[\[Theta]]-plot to count roots! There should be: "<>ToString[lin-Abs[min]]];
];
,{}];
If[m<0&&win<0&&Alm>lin*(lin+1)-s*(s+1),
Print["SWSH Eigenvalue Error 2!"];
,{}];
,{}];
]
Print["----------------------------------------------"];
Print["SWSHs functions:"];
Print["SWSHrun[m,l,s,\[Chi],Re(\[Omega])] :> SWSH[\[Theta]] and the corresponding Alm."]


End[]

EndPackage[]
