(* ::Package:: *)

BeginPackage["MSTformalism`"]

MSTgRoots::usage = "MSTgRoots[spin,s,m,l] :> Output: nureal, nuimaginary as roots of g(\[Nu])"
nurealMY
nuimaginaryMY
nufullBHPT
fout
fmout
bout
bm\[Nu]m1out
foverfout
RLout
Rout
Lout
RLm\[Nu]m1out

MSTSolution::usage = "MSTgRoots[spin,g(\[Nu])-root as complex number,Re(\[Omega]),Alm,s,m,l] :> Output: Binc and RinOutput"
RinOutput
R\[Nu]COutput
RinHyper
RinMST
RinHyperMST
RC1
RC2
RC1Hyper
RC2Hyper
Binc
Bref
AlmSWSH
Almset
AinMST
AoutMST
gsol
Kout
K\[Nu]
K\[Nu]MST
Km\[Nu]m1
Km\[Nu]m1MST
Rnumerical

Begin["PrivateMST`"]
Needs["SpinWeightedSpheroidalHarmonics`","./SpinWeightedSpheroidalHarmonics.m"]
Needs["Teukolsky`RenormalizedAngularMomentum`","./RenormalizedAngularMomentum.m"]
prec=SetPrecision[#,20]&;

(*For the spheroidal harmonics*)
\[Alpha]leaver[n_,k1_]:=prec@(-2(n+1)(n+2k1+1));
\[Beta]leaver[n_,c_,k1_,k2_,s_]:=prec@(n(n-1)+2n(k1+k2+1-2c)-(2c(2k1+s+1)-(k1+k2)(k1+k2+1))-(c^2+s(s+1)+Alm));
\[Gamma]leaver[n_,c_,k1_,k2_,s_]:=prec@(2c(n+k1+k2+s))

(*For the MST g(\[Nu])*)
alpha[n_,\[Nu]_]:=(I \[Epsilon] \[Kappa](n+\[Nu]+1+s+I \[Epsilon])(n+\[Nu]+1+s-I \[Epsilon])(n+\[Nu]+1+I \[Tau]))/((n+\[Nu]+1)(2n+2\[Nu]+3))
beta[n_,\[Nu]_]:=-\[Lambda]-s(s+1)+(n+\[Nu])(n+\[Nu]+1)+\[Epsilon]^2+\[Epsilon](\[Epsilon]-m \[Chi])+(\[Epsilon](\[Epsilon]-m \[Chi])(s^2+\[Epsilon]^2))/((n+\[Nu])(n+\[Nu]+1))
gamma[n_,\[Nu]_]:=(-I \[Epsilon] \[Kappa](n+\[Nu]-s+I \[Epsilon])(n+\[Nu]-s-I \[Epsilon])(n+\[Nu]-I \[Tau]))/((n+\[Nu])(2n+2\[Nu]-1))


MSTgRoots[spin_,sin_,min_,lin_,wmax_]:=Module[{spin1=spin,s1=sin,m1=min,l1=lin},
\[Chi]=SetPrecision[spin,64];
s=sin;
m=min;
l=lin;
points=100;
prec=SetPrecision[#,20]&;
\[Omega]value=prec@Join[Range[0,wmax,wmax/points]];
Almguess=prec@(l(l+1)-s(s+1));
\[Nu]guess=prec@l;
gamma = PrivateMST`gamma;
beta=PrivateMST`beta;
alpha=PrivateMST`alpha;
Clear[Almset,L,R,gsol,\[Epsilon],c,one,zero,\[Tau],\[Kappa],k1,k2,\[Nu]];
If[lin-Abs[min]>2,Print["Error: The MST-code is not reliable for l-|m|>2!"]; Abort[];,{}];
If[m<0,Print["Error: The MST-code is not reliable for m<0; use symmetry relations for Rin!"]; Abort[];,{}];

\[Kappa]=prec@Sqrt[1-\[Chi]^2];
k1=prec@(1/2 Abs[m-s]);
k2=prec@(1/2 Abs[m+s]);
gfuncall=ConstantArray[0,{points}];
Do[
Clear[\[Epsilon],\[Lambda],c,\[Nu],\[Tau]];
\[Epsilon]=prec@2\[Omega]value[[i]];
\[Tau]=prec@(\[Epsilon]-m \[Chi])/\[Kappa];
c=prec@\[Chi] \[Omega]value[[i]];
\[Nu]=prec@(\[Nu]r+I \[Nu]i);

PrivateMST`\[Epsilon]=\[Epsilon];
PrivateMST`\[Lambda]=\[Lambda];
PrivateMST`\[Kappa]=\[Kappa];
PrivateMST`m=m;
PrivateMST`s=s;
PrivateMST`\[Chi]=\[Chi];

(*Angular eigvenvalues; conventions from Berti et. al.*)
(*Almset[i]=prec@FindRoot[\[Beta]leaver[0,c,k1,k2,s]+prec@ContinuedFractionK[-\[Alpha]leaver[n-1,k1]\[Gamma]leaver[n,c,k1,k2,s],\[Beta]leaver[n,c,k1,k2,s],{n,1,100}],{Alm,Almguess},WorkingPrecision->MachinePrecision];
Almguess=prec@If[lin-Abs[min]==2,If[i<=0.02*points,Almguess,Alm/.Almset[i]],If[i<=1,Almguess,Alm/.Almset[i]]];*)
Almguess=SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,lin,min,c]+2min c-c^2;
Almset[i]={Alm->Almguess};
\[Lambda]=prec@(Almguess-2c m+c^2);
(*the renormalized angular momentum; conventions from Sasaki et. al.*)
R[1]=prec@ContinuedFractionK[If[n==1,-gamma[1,\[Nu]],-alpha[n-1,\[Nu]]gamma[n,\[Nu]]],beta[n,\[Nu]],{n,1,20}];
R[0]=prec@ContinuedFractionK[If[n==1,-gamma[0,\[Nu]],-alpha[n-2,\[Nu]]gamma[n-1,\[Nu]]],beta[n-1,\[Nu]],{n,1,20}];
L[-1]=prec@ContinuedFractionK[If[n==1,-alpha[-1,\[Nu]],-alpha[-n,\[Nu]]gamma[-n+1,\[Nu]]],beta[-n,\[Nu]],{n,1,20}];
Clear[gfunc];
gfunc[\[Nu]tr_,\[Nu]ti_]:=beta[0,\[Nu]]+alpha[0,\[Nu]]R[1]+gamma[0,\[Nu]]L[-1]//.{\[Nu]i->\[Nu]ti,\[Nu]r->\[Nu]tr};
gfuncall[[i]]=gfunc[\[Nu]r,\[Nu]i];
gsol[i]=SetPrecision[If[Re[\[Nu]guess]<=1.01(lin-0.5),{With[{\[Nu]ri=If[lin==2,-0.5,lin-0.5]},{\[Nu]i->Abs[Re[\[Nu]i/.FindRoot[Re[gfunc[\[Nu]ri,\[Nu]i]]==0,{\[Nu]i,Im[\[Nu]guess]+1},MaxIterations->200]]]}],\[Nu]r->lin-0.5}//Flatten,
												 {With[{\[Nu]ii=0},FindRoot[Re[gfunc[\[Nu]r,\[Nu]ii]]==0,{\[Nu]r,Re[\[Nu]guess]},MaxIterations->200]],\[Nu]i->0}//Flatten
										          ],20];
\[Nu]guess=SetPrecision[(\[Nu]r+I \[Nu]i+10^-6 (1+I))/.gsol[i],20];
(*Tests of the solution*)
zero1[i]=gfunc[\[Nu]r/.gsol[i],\[Nu]i/.gsol[i]]; (*result has to be zero by construction*)
zero2[i]=Re[gfunc[\[Nu]r,\[Nu]i]/.gsol[i]]; (*result has to be zero by construction*)
zero3[i]=Im[gfunc[\[Nu]r,\[Nu]i]/.gsol[i]]; (*result has to be zero by construction*)
one[i]=R[0]L[-1]/.gsol[i]; (*result has to be one by construction*)
,{i,1,points}];

(*AlmSWSH[w_]:=prec@Interpolation[Table[prec@{\[Omega]value[[i]],Alm/.Almset[i]},{i,1,points}],Method->"Spline",InterpolationOrder->1][w];*)
AlmSWSH[w_]:=prec@(SpinWeightedSpheroidalEigenvalue[-2,lin,min,SetPrecision[spin*w,64]]+2min spin w-(spin w)^2);
nurealMY[w_]:=prec@Interpolation[Table[{\[Omega]value[[i]],\[Nu]r/.gsol[i]},{i,1,points}],Method->"Spline",InterpolationOrder->1][w];
nuimaginaryMY[w_]:=prec@Interpolation[Table[{\[Omega]value[[i]],\[Nu]i/.gsol[i]},{i,1,points}],Method->"Spline",InterpolationOrder->1][w];
nufullBHPT[w_]:=prec@RenormalizedAngularMomentum[-2,lin,min,\[Chi],2 SetPrecision[w,64],Method->"Monodromy"];
(*nufullBHPT[w_]:=prec@If[Im[RenormalizedAngularMomentum[-2,lin,min,\[Chi],2 SetPrecision[w,64],Method->"Monodromy"]]==0,RenormalizedAngularMomentum[-2,lin,min,\[Chi],2 SetPrecision[w,64],Method->"Monodromy"],lin-1/2+I Im[RenormalizedAngularMomentum[-2,lin,min,\[Chi],2 SetPrecision[w,64],Method->"Monodromy"]]];*)
]


(*Constructing the corresponding solution based on the irregular confluent hypergeometric functions*)
MSTSolution[spin1_,\[Nu]sol1_,wpick_,Almin_,sin_,min_,lin_,rstopCode_]:=Block[{\[Chi]=spin1,\[Nu]sol=\[Nu]sol1,\[Omega]pick=wpick,s1=sin,m1=min,l1=lin},
prec=SetPrecision[#,20]&;
prech=SetPrecision[#,100]&;
NmaxL=If[min==0,40,If[min<7,60,60]];
Nmax=30;
Clear[m,l,\[Nu],\[Lambda],\[Epsilon],\[Kappa],\[Tau],\[Epsilon]p,\[Chi],s];
m=min;
l=lin;
s=sin;
\[Lambda]=SetPrecision[(Almin-2 spin1 wpick min+(spin1 wpick)^2),100];
\[Epsilon]=prech@2wpick;
\[Kappa]=prech@Sqrt[1-spin1^2];
\[Tau]=prech@((\[Epsilon]-m spin1)/\[Kappa]);
gsolexpl=SetPrecision[{\[Nu]r->Re[\[Nu]sol1],\[Nu]i->Im[\[Nu]sol1]},100];
\[Chi]=prech@spin1;
Do[
R[j]=SetPrecision[ContinuedFractionK[If[n==1,-gamma[j,\[Nu]sol1],-alpha[n-2+j,\[Nu]sol1]gamma[n+j-1,\[Nu]sol1]],beta[n+j-1,\[Nu]sol1],{n,1,50}],64];
L[j]=SetPrecision[ContinuedFractionK[If[n==1,-alpha[j,\[Nu]sol1],-alpha[-n+j+1,\[Nu]sol1]gamma[-n+j+2,\[Nu]sol1]],beta[-n+j+1,\[Nu]sol1],{n,1,50}],64];
,{j,-NmaxL,NmaxL}];

Do[
f[0]=1;
f[n]=SetPrecision[Product[R[k],{k,1,n}],64];
f[-n]=SetPrecision[Product[L[k],{k,-n,-1}],64];
,{n,1,NmaxL}];
Do[
fm\[Nu]m1[i]=f[-i];
,{i,-NmaxL,NmaxL}];

fout=Table[f[n],{n,-NmaxL,NmaxL}];
Rout=Table[R[n],{n,-NmaxL,NmaxL}];
Lout=Table[L[n],{n,-NmaxL,NmaxL}];

(*Sasaki et. al. notation*)
xofr[r_]:=prec@((1+Sqrt[1-spin1^2]-r)/(2Sqrt[1-spin1^2]));
z=prec@(\[Epsilon] \[Kappa](1-x));
\[Epsilon]p=prec@((\[Epsilon]+\[Tau])/2);
sMST=s;
rMST=0;
IrrConflHyper[a_,c_,x_]:=Gamma[1-c]/Gamma[a-c+1] Hypergeometric1F1[a,c,x]+Gamma[c-1]/Gamma[a] x^(1-c) Hypergeometric1F1[a-c+1,2-c,x]; (*From eq. 6.5(7) in "Higher transcendental functions by Erdelyi et. al (CalTech) 1981 ||||| See also https://dlmf.nist.gov/33.2 *)
ConvenientHypergeometric[a_,b_,c_,x_]:=(Gamma[c]Gamma[b-a])/(Gamma[b]Gamma[c-a]) (1-x)^-a Hypergeometric2F1[a,c-b,a-b+1,1/(1-x)]+(Gamma[c]Gamma[a-b])/(Gamma[a]Gamma[c-b]) (1-x)^-b Hypergeometric2F1[c-a,b,b-a+1,1/(1-x)];

(*The Columbic expansion*)
\[Nu]=SetPrecision[\[Nu]sol1,100];
Print["\[Nu]sol = "<>ToString[\[Nu]]];
R\[Nu]plus=prec@(2^\[Nu] Exp[-\[Pi] \[Epsilon]]Exp[I \[Pi](\[Nu]+1-sMST)] Gamma[\[Nu]+1-sMST+I \[Epsilon]]/Gamma[\[Nu]+1+sMST-I \[Epsilon]] Exp[-I z]z^(\[Nu]+I \[Epsilon]p) (z-\[Epsilon] \[Kappa])^(-sMST-I \[Epsilon]p) Sum[I^n f[n](2z)^n IrrConflHyper[n+\[Nu]+1-sMST+I \[Epsilon],2n+2\[Nu]+2,2I z],{n,-Nmax,Nmax}]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
R\[Nu]minus=prec@(2^\[Nu] Exp[-\[Pi] \[Epsilon]]Exp[-I \[Pi](\[Nu]+1+sMST)]Exp[I z]z^(\[Nu]+I \[Epsilon]p) (z-\[Epsilon] \[Kappa])^(-sMST-I \[Epsilon]p) Sum[I^n Pochhammer[\[Nu]+1+sMST-I \[Epsilon],n]/Pochhammer[\[Nu]+1-sMST+I \[Epsilon],n] f[n](2z)^n IrrConflHyper[n+\[Nu]+1+sMST-I \[Epsilon],2n+2\[Nu]+2,-2I z],{n,-Nmax,Nmax}]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
R\[Nu]C=prec@(R\[Nu]plus+R\[Nu]minus);
Rm\[Nu]m1plus=prec@((2^\[Nu] Exp[-\[Pi] \[Epsilon]]Exp[I \[Pi](\[Nu]+1-sMST)] Gamma[\[Nu]+1-sMST+I \[Epsilon]]/Gamma[\[Nu]+1+sMST-I \[Epsilon]] Exp[-I z]z^(\[Nu]+I \[Epsilon]p) (z-\[Epsilon] \[Kappa])^(-sMST-I \[Epsilon]p) Sum[I^n fm\[Nu]m1[n](2z)^n IrrConflHyper[n+\[Nu]+1-sMST+I \[Epsilon],2n+2\[Nu]+2,2I z],{n,-Nmax,Nmax}])/.{\[Nu]->-\[Nu]-1}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
Rm\[Nu]m1minus=prec@((2^(-\[Nu]-1) Exp[-\[Pi] \[Epsilon]]Exp[-I \[Pi](-\[Nu]-1+1+sMST)]Exp[I z]z^(-\[Nu]-1+I \[Epsilon]p) (z-\[Epsilon] \[Kappa])^(-sMST-I \[Epsilon]p) Sum[I^n Pochhammer[-\[Nu]-1+1+sMST-I \[Epsilon],n]/Pochhammer[-\[Nu]-1+1-sMST+I \[Epsilon],n] fm\[Nu]m1[n](2z)^n IrrConflHyper[n+-\[Nu]-1+1+sMST-I \[Epsilon],2n+2(-\[Nu]-1)+2,-2I z],{n,-Nmax,Nmax}])/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
R\[Nu]Cm\[Nu]m1=prec@(Rm\[Nu]m1plus+Rm\[Nu]m1minus);

(*The hypergeometric expansion*)
R0\[Nu]=Exp[I \[Epsilon] \[Kappa] x](-x)^(-sMST-(I/2)(\[Epsilon]+\[Tau])) (1-x)^((I/2)(\[Epsilon]+\[Tau])+\[Nu]) Sum[f[n] (Gamma[1-sMST-I \[Epsilon]-I \[Tau]]Gamma[2n+2\[Nu]+1])/(Gamma[n+\[Nu]+1-I \[Tau]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]]) (1-x)^n ConvenientHypergeometric[-n-\[Nu]-I \[Tau],-n-\[Nu]-sMST-I \[Epsilon],-2n-2\[Nu],1/(1-x)],{n,-Nmax,Nmax}]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl;
R0m\[Nu]m1=(Exp[I \[Epsilon] \[Kappa] x](-x)^(-sMST-(I/2)(\[Epsilon]+\[Tau])) (1-x)^((I/2)(\[Epsilon]+\[Tau])+\[Nu]) Sum[fm\[Nu]m1[n] (Gamma[1-sMST-I \[Epsilon]-I \[Tau]]Gamma[2n+2\[Nu]+1])/(Gamma[n+\[Nu]+1-I \[Tau]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]]) (1-x)^n ConvenientHypergeometric[-n-\[Nu]-I \[Tau],-n-\[Nu]-sMST-I \[Epsilon],-2n-2\[Nu],1/(1-x)],{n,-Nmax,Nmax}])/.{\[Nu]->-\[Nu]-1}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl;
RC1Hyper[rr_]:=R0\[Nu]//.gsolexpl//.{x->xofr[rr]}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl;
RC2Hyper[rr_]:=R0m\[Nu]m1//.gsolexpl//.{x->xofr[rr]}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl;
RinHyper[rr_]:=R0\[Nu]+R0m\[Nu]m1//.gsolexpl//.{x->xofr[rr]}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl;

(*Matching relations*)
Cnj[n_,j_]:=(Gamma[1-sMST-2I \[Epsilon]p]Gamma[2n+2\[Nu]+1])/(Gamma[n+\[Nu]+1-I \[Tau]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]]) (Pochhammer[-n-\[Nu]-I \[Tau] ,j]Pochhammer[-n-\[Nu]-sMST-I \[Epsilon],j])/(Pochhammer[-2n-2\[Nu],j]Factorial[j]) (\[Epsilon] \[Kappa])^(-n+j) f[n];
Dnj[n_,j_]:=(-1)^n (2I)^(n+j) Gamma[n+\[Nu]+1-sMST+I \[Epsilon]]/Gamma[2n+2\[Nu]+2] Pochhammer[\[Nu]+1+sMST-I \[Epsilon],n]/Pochhammer[\[Nu]+1-sMST+I \[Epsilon],n] Pochhammer[n+\[Nu]+1-sMST+I \[Epsilon],j]/(Pochhammer[2n+2\[Nu]+2,j]Factorial[j]) f[n];
K\[Nu]1=prec@(Exp[I \[Epsilon] \[Kappa]](\[Epsilon] \[Kappa])^(sMST-\[Nu]) 2^-\[Nu] Sum[Dnj[n,rMST-n],{n,-Nmax,rMST}]^-1 Sum[Cnj[n,n-rMST],{n,rMST,Nmax}]);
K\[Nu]=prec@(((Exp[I \[Epsilon] \[Kappa]](2\[Epsilon] \[Kappa])^(sMST-\[Nu]-rMST) 2^-sMST I^rMST Gamma[1-sMST-2I \[Epsilon]p ]Gamma[rMST+2\[Nu]+2])/(Gamma[rMST+\[Nu]+1-sMST+I \[Epsilon]]Gamma[rMST+\[Nu]+1+I \[Tau]]Gamma[rMST+\[Nu]+1+sMST+I \[Epsilon]]))Sum[(-1)^n Gamma[n+rMST+2\[Nu]+1]/Factorial[n-rMST] Gamma[n+\[Nu]+1+sMST+I \[Epsilon]]/Gamma[n+\[Nu]+1-sMST-I \[Epsilon]] Gamma[n+\[Nu]+1+I \[Tau]]/Gamma[n+\[Nu]+1-I \[Tau]] f[n],{n,rMST,NmaxL}]Sum[(-1)^n/(Factorial[rMST-n]Pochhammer[rMST+2\[Nu]+2,n]) Pochhammer[\[Nu]+1+sMST-I \[Epsilon],n]/Pochhammer[\[Nu]+1-sMST+I \[Epsilon],n] f[n],{n,Range[-NmaxL,rMST]}]^-1/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
Km\[Nu]m1=prec@((((Exp[I \[Epsilon] \[Kappa]](2\[Epsilon] \[Kappa])^(sMST-(-\[Nu]-1)-rMST) 2^-sMST I^rMST Gamma[1-sMST-2I \[Epsilon]p ]Gamma[rMST+2(-\[Nu]-1)+2])/(Gamma[rMST+(-\[Nu]-1)+1-sMST+I \[Epsilon]]Gamma[rMST+(-\[Nu]-1)+1+I \[Tau]]Gamma[rMST+(-\[Nu]-1)+1+sMST+I \[Epsilon]]))Sum[(-1)^n Gamma[n+rMST+2(-\[Nu]-1)+1]/Factorial[n-rMST] Gamma[n+(-\[Nu]-1)+1+sMST+I \[Epsilon]]/Gamma[n+(-\[Nu]-1)+1-sMST-I \[Epsilon]] Gamma[n+(-\[Nu]-1)+1+I \[Tau]]/Gamma[n+(-\[Nu]-1)+1-I \[Tau]] fm\[Nu]m1[n],{n,rMST,NmaxL}]Sum[(-1)^n/(Factorial[rMST-n]Pochhammer[rMST+2(-\[Nu]-1)+2,n]) Pochhammer[(-\[Nu]-1)+1+sMST-I \[Epsilon],n]/Pochhammer[(-\[Nu]-1)+1-sMST+I \[Epsilon],n] fm\[Nu]m1[n],{n,Range[-NmaxL,rMST]}]^-1)/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
Kout=prec@((K\[Nu]-K\[Nu]1)//.gsolexpl);

(*Functions determined through matching*)
RinMatching=prec@(K\[Nu] R\[Nu]C+Km\[Nu]m1 R\[Nu]Cm\[Nu]m1);
RC1[rr_]:=prec@(K\[Nu] R\[Nu]C//.gsolexpl//.{x->xofr[rr]});
RC2[rr_]:=prec@(Km\[Nu]m1 R\[Nu]Cm\[Nu]m1//.gsolexpl//.{x->xofr[rr]});
RinOutput[rr_]:=prec@(RinMatching//.gsolexpl//.{x->xofr[rr]});

(*For the asymtptotic amplitude*)
Aplus=prec@(Exp[-(\[Pi]/2)\[Epsilon]]Exp[(\[Pi]/2)I(\[Nu]+1-sMST)]2^(-1+sMST-I \[Epsilon]) Gamma[\[Nu]+1-sMST+I \[Epsilon]]/Gamma[\[Nu]+1+sMST-I \[Epsilon]] Sum[f[n],{n,-NmaxL,NmaxL}]);
Binc=prec@(wpick^-1 (K\[Nu]-I Exp[-I \[Pi] \[Nu]]Km\[Nu]m1 Sin[\[Pi](\[Nu]-sMST+I \[Epsilon])]/Sin[\[Pi](\[Nu]+sMST-I \[Epsilon])])Aplus Exp[-I(\[Epsilon] Log[\[Epsilon]](*-\[Epsilon] (1-\[Kappa])/2*))]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
Aminus=prec@(2^(-1-sMST+I \[Epsilon])Exp[-(\[Pi]/2)I(\[Nu]+1+sMST)]Exp[-(\[Pi]/2)\[Epsilon]]Sum[(-1)^n Pochhammer[\[Nu]+1+sMST-I \[Epsilon],n]/Pochhammer[\[Nu]+1-sMST+I \[Epsilon],n]f[n],{n,-NmaxL,NmaxL}]);
Bref=prec@(wpick^(-1-2sMST)(K\[Nu]+I Exp[I \[Pi] \[Nu]]Km\[Nu]m1)Aminus Exp[I(\[Epsilon] Log[\[Epsilon]]-\[Epsilon] (1-\[Kappa])/2)]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);

(*MST's original notation*)
zp=prec@(-\[Epsilon] \[Kappa] x);
Do[
b[n]=prec@(I^n Gamma[\[Nu]+1-sMST+I \[Epsilon]]Gamma[\[Nu]+1-sMST-I \[Epsilon]]Gamma[n+\[Nu]+1+sMST+I \[Epsilon]]Gamma[n+\[Nu]+1+sMST-I \[Epsilon]]/(Gamma[\[Nu]+1+sMST+I \[Epsilon]]Gamma[\[Nu]+1+sMST-I \[Epsilon]]Gamma[n+\[Nu]+1-sMST+I \[Epsilon]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]])f[n]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
bm\[Nu]m1[n]=prec@(I^n Gamma[\[Nu]+1-sMST+I \[Epsilon]]Gamma[\[Nu]+1-sMST-I \[Epsilon]]Gamma[n+\[Nu]+1+sMST+I \[Epsilon]]Gamma[n+\[Nu]+1+sMST-I \[Epsilon]]/(Gamma[\[Nu]+1+sMST+I \[Epsilon]]Gamma[\[Nu]+1+sMST-I \[Epsilon]]Gamma[n+\[Nu]+1-sMST+I \[Epsilon]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]])fm\[Nu]m1[n]/.{\[Nu]->-\[Nu]-1}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
,{n,-NmaxL,NmaxL}];
RCin=prec@(Exp[-I zp]zp^(\[Nu]-sMST) (1+(\[Epsilon] \[Kappa])/zp)^((I/2)(\[Epsilon]-\[Tau])) 2^\[Nu] Exp[I \[Pi] (\[Nu]+1-sMST+I \[Epsilon])]Sum[b[n](-2zp)^n Gamma[n+\[Nu]+1-sMST+I \[Epsilon]]/Gamma[n+\[Nu]+1+sMST-I \[Epsilon]]IrrConflHyper[n+\[Nu]+1-sMST+I \[Epsilon],2n+2\[Nu]+2,2I zp],{n,-Nmax,Nmax}]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
RCout=prec@(Exp[I zp]zp^(\[Nu]-sMST) (1+(\[Epsilon] \[Kappa])/zp)^((I/2)(\[Epsilon]-\[Tau])) 2^\[Nu] Exp[-I \[Pi] (\[Nu]+1+sMST-I \[Epsilon])]Sum[b[n](-2zp)^n IrrConflHyper[n+\[Nu]+1+sMST-I \[Epsilon],2n+2\[Nu]+2,-2I zp],{n,-Nmax,Nmax}]/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
RCm\[Nu]m1MST=prec@(-I Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-sMST+I \[Epsilon])]/Sin[\[Pi](\[Nu]+sMST-I \[Epsilon])]RCin+I Exp[I \[Pi] \[Nu]]RCout/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
K\[Nu]MST=prec@((\[Epsilon] \[Kappa])^(-\[Nu]-rMST+sMST) 2^(-\[Nu]-rMST) (-I)^rMST Gamma[1-sMST-I \[Epsilon]-I \[Tau]]/(Gamma[1+rMST+\[Nu]+I \[Tau]]Gamma[1+rMST+\[Nu]-sMST-I \[Epsilon]]Gamma[1+rMST+\[Nu]-sMST+I \[Epsilon]])Sum[Gamma[n+\[Nu]+1+I \[Tau]]Gamma[n+rMST+2\[Nu]+1]/(Factorial[n-rMST]Gamma[n+\[Nu]+1-I \[Tau]])f[n],{n,rMST,NmaxL}]Sum[((-I)^n b[n])/(Factorial[rMST-n]Gamma[n+rMST+2\[Nu]+2]),{n,-NmaxL,rMST}]^-1/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
Km\[Nu]m1MST=prec@((\[Epsilon] \[Kappa])^(-(-\[Nu]-1)-rMST+sMST) 2^(-(-\[Nu]-1)-rMST) (-I)^rMST Gamma[1-sMST-I \[Epsilon]-I \[Tau]]/(Gamma[1+rMST+(-\[Nu]-1)+I \[Tau]]Gamma[1+rMST+(-\[Nu]-1)-sMST-I \[Epsilon]]Gamma[1+rMST+(-\[Nu]-1)-sMST+I \[Epsilon]])Sum[Gamma[n+(-\[Nu]-1)+1+I \[Tau]]Gamma[n+rMST+2(-\[Nu]-1)+1]/(Factorial[n-rMST]Gamma[n+(-\[Nu]-1)+1-I \[Tau]])fm\[Nu]m1[n],{n,rMST,NmaxL}]Sum[((-I)^n bm\[Nu]m1[n])/(Factorial[rMST-n]Gamma[n+rMST+2(-\[Nu]-1)+2]),{n,-NmaxL,rMST}]^-1/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
RinMSTp=prec@((K\[Nu]MST-I Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-sMST+I \[Epsilon])]/Sin[\[Pi](\[Nu]+sMST-I \[Epsilon])]Km\[Nu]m1MST)RCin+(K\[Nu]MST+I Exp[I \[Pi] \[Nu]]Km\[Nu]m1MST)RCout);
AinMSTp=prec@(Exp[-(I/2)\[Pi](-\[Nu]-1+sMST-I \[Epsilon])]2^(-1+sMST-I \[Epsilon]) (K\[Nu]MST-I Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-sMST+I \[Epsilon])]/Sin[\[Pi](\[Nu]+sMST-I \[Epsilon])]Km\[Nu]m1MST)Sum[b[n]I^n Gamma[n+\[Nu]+1-sMST+I \[Epsilon]]/Gamma[n+\[Nu]+1+sMST-I \[Epsilon]],{n,-NmaxL,NmaxL}]);
AoutMSTp=prec@(Exp[(I/2)\[Pi](\[Nu]+1+sMST-I \[Epsilon])]2^(-1-sMST+I \[Epsilon])(K\[Nu]MST+I Exp[I \[Pi] \[Nu]]Km\[Nu]m1MST)Sum[b[n](-I)^n,{n,-NmaxL,NmaxL}]);
RinMST[rr_]:=prec@(RinMSTp//.{x->xofr[rr]}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
AinMST=prec@(AinMSTp/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);
AoutMST=prec@(AoutMSTp/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);

(*hypegeoemtric expansion in original MST formulation*)
R0\[Nu]MST=prec@(Exp[I \[Epsilon] \[Kappa] x](-x)^(\[Nu]-sMST-(I/2)(\[Epsilon]-\[Tau]))(1-x)^((I/2)(\[Epsilon]-\[Tau]))Sum[f[n]Gamma[1-sMST-I \[Epsilon]-I \[Tau]]Gamma[2n+2\[Nu]+1]/(Gamma[n+\[Nu]+1-I \[Tau]]Gamma[n+\[Nu]+1-sMST-I \[Epsilon]])(-x)^n Hypergeometric2F1[-n-\[Nu]-I \[Tau],-n-\[Nu]+sMST+I \[Epsilon],-2n-2\[Nu],1/x],{n,-Nmax,Nmax}]);
R0m\[Nu]m1MST=prec@(Exp[I \[Epsilon] \[Kappa] x](-x)^((-\[Nu]-1)-sMST-(I/2)(\[Epsilon]-\[Tau]))(1-x)^((I/2)(\[Epsilon]-\[Tau]))Sum[fm\[Nu]m1[n]Gamma[1-sMST-I \[Epsilon]-I \[Tau]]Gamma[2n+2(-\[Nu]-1)+1]/(Gamma[n+(-\[Nu]-1)+1-I \[Tau]]Gamma[n+(-\[Nu]-1)+1-sMST-I \[Epsilon]])(-x)^n Hypergeometric2F1[-n-(-\[Nu]-1)-I \[Tau],-n-(-\[Nu]-1)+sMST+I \[Epsilon],-2n-2(-\[Nu]-1),1/x],{n,-Nmax,Nmax}]);
RinHyperMST[rr_]:=prec@((R0\[Nu]MST+R0m\[Nu]m1MST)//.{x->xofr[rr]}/.{\[Nu]->\[Nu]r+I \[Nu]i}/.gsolexpl);

(*Numerical solution*)
Clear[\[Kappa]];
setprec=MachinePrecision;
maxstep=10^9;
\[CapitalDelta][r_]:=SetPrecision[r^2-2r+spin1^2,setprec];
Kf[r_,\[Chi]t_,\[Omega]t_,mt_]:=SetPrecision[(r^2+\[Chi]t^2)\[Omega]t-\[Chi]t mt,setprec];
Lambda[Alm_,\[Chi]t_,\[Omega]t_,mt_]:=SetPrecision[Alm+\[Chi]t^2 \[Omega]t^2-2\[Chi]t mt \[Omega]t,setprec];
Diffr=\[CapitalDelta][r]^-s D[\[CapitalDelta][r]^(s+1) D[#,r],r]+((Kf[r,spin1,wpick,min]^2-2I s(r-1)Kf[r,spin1,wpick,min])/\[CapitalDelta][r]+4I s wpick r-Lambda[Almin,spin1,wpick,min])#&;
(*Collect[Diffr@f[r]/.{f\[Rule](f[(#-rp[\[Chi]])/(rp[\[Chi]]-rm[\[Chi]])]&)}/.{r\[Rule]x(rp[\[Chi]]-rm[\[Chi]])+rp[\[Chi]]}//Simplify,{f[x],f'[x],f''[x]}]*)
rp[\[Chi]t_]:=SetPrecision[1+Sqrt[1-\[Chi]t^2],20];
rm[\[Chi]t_]:=SetPrecision[1-Sqrt[1-\[Chi]t^2],20];
rstart[a_]:=SetPrecision[rp[a]+10^-5,20];
rstop[w_]:=SetPrecision[20(10(min))/w^2,20];
rstophere=(3/2)*rstopCode;
Diffx= (-Almin+2 min spin1 wpick-spin1^2 wpick^2+4 I s wpick (-x rm[spin1]+(1+x) rp[spin1])+((-min spin1+wpick (spin1^2+(x rm[spin1]-(1+x) rp[spin1])^2)) (-min spin1-2 I s (-1-x rm[spin1]+(1+x) rp[spin1])+wpick (spin1^2+(x rm[spin1]-(1+x) rp[spin1])^2)))/(spin1^2+2 x rm[spin1]-2 (1+x) rp[spin1]+(x rm[spin1]-(1+x) rp[spin1])^2))#+(2 (1+s) (1+x rm[spin1]-(1+x) rp[spin1]) )/(rm[spin1]-rp[spin1]) D[#,x]+(spin1^2+2 x rm[spin1]-2 (1+x) rp[spin1]+(x rm[spin1]-(1+x) rp[spin1])^2) /(rm[spin1]-rp[spin1])^2 D[#,{x,2}]&;
Frob[x_]:=x^(-I \[Kappa]) (1+Sum[Symbol["PrivateMST`c"<>ToString[n]]x^n,{n,1,2}]);
exptemp=Series[SetPrecision[x^(I \[Kappa]+2) (Diffx@Frob[x])//Simplify,20],{x,0,3}]==0;
exptemp1=SetPrecision[exptemp//LogicalExpand,20];
Coeff=NSolve[exptemp1,{\[Kappa],c1,c2},WorkingPrecision->setprec][[1]];
FrobOfr[r_]:=SetPrecision[Frob[x]//.Coeff//.{x->(r-rp[spin1])/(rp[spin1]-rm[spin1])}//Simplify,20];
adaptivestep=If[wpick<0.45,0.01,If[wpick<1,0.001,0.0005]];
R0=SetPrecision[Limit[FrobOfr[r],r->rstart[spin1]],20];
dR0=SetPrecision[Limit[D[FrobOfr[r],r],r->rstart[spin1]],20];
solR=NDSolve[{Diffr@R[r]==0,R[rstart[spin1]]==R0,R'[rstart[spin1]]==dR0},R,{r,rstart[spin1],rstophere},Method->{"TimeIntegration"->{"ExplicitRungeKutta","DifferenceOrder"->8}}(*,Method->{"StiffnessSwitching"}*),WorkingPrecision->setprec,MaxSteps->maxstep,PrecisionGoal->15,AccuracyGoal->15(*,MaxStepSize->adaptivestep*)]//First;
norm=SetPrecision[RinOutput[(2/3)*rstophere]/(R[(2/3)*rstophere]/.solR),100];
Rnumerical[rr_]:=SetPrecision[norm R[r]/.solR//.{r->rr},setprec];
];


Print["------------------------------------------------------"];
Print["MSTformalism Functions:"];
Print["MSTgRoots[spin,s,m,l,\[Omega]_max(max Re(\[Omega]) to which code should be evolved)]:> nureal[w] and nuimaginary[w]"];
Print["MSTSolution[spin,\[Nu]sol,Re(\[Omega]),Alm,s,m,l] :> RinOutpur[r] (based on irregular confluent hypergeo. functions)"];
Print["------------------------------------------------------"];

End[]

EndPackage[]
