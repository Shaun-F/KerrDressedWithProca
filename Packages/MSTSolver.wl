(* ::Package:: *)

(* ::Chapter:: *)
(*Setup*)


MSTSolver;
<<Teukolsky`
<<SpinWeightedSpheroidalHarmonics`

$MSTMaxN = 30; (*Max index of series expansions to use throughout*)
$RAMMaxN = 500; (*Max N to use for calculating the renormalized angular momentum*)
$MSTDumbyR = 0; (*summation halting index for Equation 165 in Sasaki et al.*)
$MSTGlobalPrecision = 50;


(*helper functions*)
ToMSTPrecision[expr_]:=SetPrecision[expr, $MSTGlobalPrecision]

IrregularConfluentHGF::usage="Irregular confluent hypergeometric function (Tricomi's hypergeometric function).";
(*See https://en.wikipedia.org/wiki/Confluent_hypergeometric_function   or   https://authors.library.caltech.edu/43491/1/Volume%201.pdf    for details*)
IrregularConfluentHGF[a_,b_,z_]:=Gamma[1-b]/Gamma[a+1-b]*Hypergeometric1F1[a,b,z] + Gamma[b-1]/Gamma[a]*z^(1-b)*Hypergeometric1F1[a+1-b,2-b,z];


(*Definitions taken from Sasaki et al*)
ClearAll[\[Epsilon],z,\[Kappa],rplus,rminus,x,\[Tau],\[Epsilon]plus,\[Epsilon]minus,\[Alpha]\[Nu]n,\[Beta]\[Nu]n,\[Gamma]\[Nu]n,Rcf, Lcf, GRoot, SeriesCoeff, ConstructHGFCoefficients, MatchingCcoeff, MatchingDcoeff, MatchingKcoeff];
\[Kappa][\[Chi]_]:=Sqrt[1-\[Chi]^2];
rplus[\[Chi]_]:=1+\[Kappa][\[Chi]];
rminus[\[Chi]_]:=1-\[Kappa][\[Chi]];
\[Epsilon][\[Omega]_]:=2*\[Omega]
z[r_,\[Omega]_]:=\[Omega]*r
ztilde[r_,\[Omega]_,\[Chi]_]:=\[Omega]*(r-rminus[\[Chi]]);
x[\[Omega]_,r_,\[Chi]_]:=(z[rplus[\[Chi]],\[Omega]]-z[r,\[Omega]])/(\[Epsilon][\[Omega]]*\[Kappa][\[Chi]]);
\[Tau][\[Omega]_,\[Chi]_,m_]:=(\[Epsilon][\[Omega]]-m*\[Chi])/\[Kappa][\[Chi]];
\[Epsilon]plus[\[Omega]_,\[Chi]_,m_]:=(\[Epsilon][\[Omega]]+\[Tau][\[Omega],\[Chi],m])/2;
\[Epsilon]minus[\[Omega]_,\[Chi]_,m_]:=(\[Epsilon][\[Omega]]-\[Tau][\[Omega],\[Chi],m])/2;
\[Alpha]\[Nu]n[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_][n_]:=Block[{res},res=(I*\[Epsilon][\[Omega]]*\[Kappa][\[Chi]](n+\[Nu]+1+s+I*\[Epsilon][\[Omega]])(n+\[Nu]+1+s-I*\[Epsilon][\[Omega]])(n+\[Nu]+1+I*\[Tau][\[Omega],\[Chi],m]))/((n+\[Nu]+1)(2*n+2*\[Nu]+3));Return[res//ToMSTPrecision]];
\[Beta]\[Nu]n[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_]:=Block[{res},res=-\[Lambda]-s(s+1)+(n+\[Nu])(n+\[Nu]+1)+\[Epsilon][\[Omega]]^2+\[Epsilon][\[Omega]](\[Epsilon][\[Omega]]-m*\[Chi])+(\[Epsilon][\[Omega]](\[Epsilon][\[Omega]]-m*\[Chi])(s^2+\[Epsilon][\[Omega]]^2))/((n+\[Nu])(n+\[Nu]+1));Return[res//ToMSTPrecision]];
\[Gamma]\[Nu]n[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_][n_]:=Block[{res},res=(-I*\[Epsilon][\[Omega]]*\[Kappa][\[Chi]](n+\[Nu]-s+I*\[Epsilon][\[Omega]])(n+\[Nu]-s-I*\[Epsilon][\[Omega]])(n+\[Nu]-I*\[Tau][\[Omega],\[Chi],m]))/((n+\[Nu])(2*n+2*\[Nu]-1));Return[res//ToMSTPrecision]];
Rcf::usage="Continued fraction representation of series coefficient ratios \!\(\*FractionBox[SubscriptBox[\(f\), \(n\)], SubscriptBox[\(f\), \(n - 1\)]]\)";
Rcf[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_,MaxN_]:=
Block[{Rcfres},Rcfres=
ContinuedFractionK[
If[iter==0, 
-\[Gamma]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n], 
-\[Alpha]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n+iter-1]\[Gamma]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n+iter]
],
\[Beta]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n+iter],
 {iter,0,MaxN}];
Return[Rcfres//ToMSTPrecision]
];
Lcf::usage="Continued fraction representation of series coefficient ratios \!\(\*FractionBox[SubscriptBox[\(f\), \(n\)], SubscriptBox[\(f\), \(n + 1\)]]\)";
Lcf[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_,MaxN_]:=
Block[{res},res=
ContinuedFractionK[
If[iter==0, 
-\[Alpha]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n], 
-\[Alpha]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n-iter]\[Gamma]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][n-iter+1]
],
\[Beta]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n-iter],
 {iter,0,MaxN}];
Return[res//ToMSTPrecision]
];

GRoot[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{expr,nsample=0},
expr = \[Beta]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][nsample] + \[Alpha]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][nsample]*Rcf[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][nsample+1, $MSTMaxN] + \[Gamma]\[Nu]n[\[Nu],\[Omega],\[Chi],m,s][nsample]*Lcf[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][nsample-1, $MSTMaxN];
Return[expr//ToMSTPrecision]
];

SeriesCoeff::usage="Construct single coefficient of the hypergeometric function expansion of the radial function using continued fraction representation of coefficient ratios";
SeriesCoeff[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_]:=
With[{f0=1},
If[n>=1,

Return[f0*Product[
Rcf[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][i,$MSTMaxN],
{i,1,n}
]];,

Return[f0*Product[Lcf[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][-i, $MSTMaxN],
{i,1,Abs[n]}
]];(*Product*)

];(*If*)
];(*With*)

ConstructHGFCoefficients[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][MaxN_]:=Association@@Table[i->SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][i], {i,-MaxN,MaxN}];

MatchingCcoeff[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_,j_]:=
Block[{MatchingCcoefftmp},
With[{\[Epsilon]val = \[Epsilon][\[Omega]], \[Epsilon]plusval = \[Epsilon]plus[\[Omega],\[Chi],m], \[Tau]val = \[Tau][\[Omega],\[Chi],m],\[Kappa]val = \[Kappa][\[Chi]],coeff = SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]},

MatchingCcoefftmp = (Gamma[1-s-2*I*\[Epsilon]plusval]*Gamma[2*n+2*\[Nu]+1])/(Gamma[n+\[Nu]+1-I*\[Tau]val]*Gamma[n+\[Nu]+1-s-I*\[Epsilon]val])*(Pochhammer[(-n-\[Nu]-I*\[Tau]val),j]Pochhammer[-n-\[Nu]-s-I*\[Epsilon]val,j])/(Pochhammer[-2*n-2*\[Nu],j]Factorial[j]) (\[Epsilon]val*\[Kappa]val)^(-n + j)*coeff;
Return[MatchingCcoefftmp//ToMSTPrecision]
];
];
MatchingDcoeff[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_][n_,j_]:=
Block[{MatchingDcoefftmp},
With[{\[Epsilon]val = \[Epsilon][\[Omega]],\[Epsilon]plusval = \[Epsilon]plus[\[Omega],\[Chi],m], \[Tau]val = \[Tau][\[Omega],\[Chi],m],\[Kappa]val = \[Kappa][\[Chi]],coeff = SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]},

MatchingDcoefftmp=(-1)^n*(2*I)^(n+j)*(Gamma[n+\[Nu]+1-s+I*\[Epsilon]val]*Pochhammer[\[Nu]+1+s-I*\[Epsilon]val,n])/(Gamma[2*n+2*\[Nu]+2]*Pochhammer[\[Nu]+1-s+I*\[Epsilon]val,n])*Pochhammer[n+\[Nu]+1-s+I*\[Epsilon]val,j]/(Pochhammer[2*n+2*\[Nu]+2,j]*Factorial[j])*coeff;
Return[MatchingDcoefftmp//ToMSTPrecision]
];
];
MatchingKcoeff[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{MatchingKcoefftmp},
With[{\[Epsilon]val = \[Epsilon][\[Omega]], \[Kappa]val = \[Kappa][\[Chi]], MaxN = $MSTMaxN, rval = $MSTDumbyR},

MatchingKcoefftmp = Exp[I*\[Epsilon]val* \[Kappa]val]*(\[Epsilon]val*\[Kappa]val)^(s-\[Nu])*2^-\[Nu]*(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-MaxN\)\), \(rval\)]\(\(MatchingDcoeff[\[Nu], \[Omega], \[Chi], m, s, \[Lambda]]\)[n, rval - n]\)\))^-1*(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = rval\), \(MaxN\)]\(\(MatchingCcoeff[\[Nu], \[Omega], \[Chi], m, s, \[Lambda]]\)[n, n - rval]\)\));
Return[MatchingKcoefftmp//ToMSTPrecision]
];
];


(* ::Chapter:: *)
(*Functions*)


RenormedAngularMomentum::usage="Renormalized angular momentum to be used in MST formalism. Uses BlackHolePerturbationsToolkit's Teukolsky package";
SetAttributes[RenormedAngularMomentum, Listable];
Options[RenormedAngularMomentum]={Method->{"Monodromy", "nmax"->$RAMMaxN}};
RenormedAngularMomentum[\[Omega]_,\[Chi]_,m_,l_,s_, OptionsPattern[]]:=
Block[{testval, RenormedAngularMomentumtmp,nsample=0},
(*Used Teukolsky package's more sophisticated algorithm to compute renormalized angular momentum*)
(*!!!!! Check convergence for range of parameter space we're interested in.*)

(*Teukolsky uses SpinWeightedSpheroidalEigenvalue to compute the SWSH eigenvalue*)
testval = RenormalizedAngularMomentum[s,l,m,\[Chi],\[Omega], Method->OptionValue[Method]];

Return[testval//ToMSTPrecision];


(*I forgot why I built the following code*)
(*
If[Im[testval]==0,
Return[testval//ToMSTPrecision]
];

If[TrueQ[\[Omega]>0.6&&Im[testval]!=0],
RenormedAngularMomentumtmp = l-testval;
Return[RenormedAngularMomentumtmp//ToMSTPrecision]
];

If[TrueQ[\[Omega]<=0.6&&Im[testval]!=0],
RenormedAngularMomentumtmp = l-1+testval;
Return[RenormedAngularMomentumtmp//ToMSTPrecision];
];
*)
];

A\[Nu]plus[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{A\[Nu]plustmp},
With[{ \[Nu]val = ToMSTPrecision[\[Nu]], coeff= Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]], \[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]]},

A\[Nu]plustmp = Exp[-\[Pi]/2*\[Epsilon]val]*Exp[\[Pi]/2 I(\[Nu]val+1-s)]2^(-1+s-I*\[Epsilon]val) Gamma[\[Nu]val+1-s+I*\[Epsilon]val]/Gamma[\[Nu]val+1+s-I*\[Epsilon]val] \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-$MSTMaxN\)\), \($MSTMaxN\)]\(coeff[n]\)\);
Return[A\[Nu]plustmp//ToMSTPrecision]
];
];



A\[Nu]minus[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_,MaxIndex_:$MSTMaxN]:=
Block[{A\[Nu]minustmp},
With[{ \[Nu]val = ToMSTPrecision[\[Nu]], coeff= Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]], \[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]]},

A\[Nu]minustmp = Exp[-\[Pi]/2*\[Epsilon]val]*Exp[-(\[Pi]/2)I(\[Nu]val+1+s)]2^(-1-s+I*\[Epsilon]val) \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-MaxIndex\)\), \(MaxIndex\)]\((
\*SuperscriptBox[\((\(-1\))\), \(n\)] coeff[n]*
\*FractionBox[\(Pochhammer[\[Nu] + 1 + s - I*\[Epsilon]val, n]\), \(Pochhammer[\[Nu] + 1 - s + I*\[Epsilon]val, n]\)])\)\);
Return[A\[Nu]minustmp//ToMSTPrecision]
];
];

IncomingAmplitudeB::usage="Calculated Amplitude \!\(\*SuperscriptBox[\(B\), \(inc\)]\) of asymptotic radial function";
IncomingAmplitudeB[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{Kcoeff1, Kcoeff2, IncomingAmplitudeBtmp, Aplus},
With[{\[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]], \[Kappa]val = ToMSTPrecision[\[Kappa][\[Chi]]]},

Kcoeff1 = MatchingKcoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]];
Kcoeff2 = MatchingKcoeff[-\[Nu]-1,\[Omega],\[Chi],m,s,\[Lambda]];
Aplus = A\[Nu]plus[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]];
IncomingAmplitudeBtmp = 1/\[Omega]*(Kcoeff1 - I*Exp[-I*\[Pi]*\[Nu]]*Sin[\[Pi](\[Nu]-s+I*\[Epsilon]val)]/Sin[\[Pi](\[Nu]+s-I*\[Epsilon]val)]*Kcoeff2)*Aplus*Exp[-I*(\[Epsilon]val*Log[\[Epsilon]val](*-(1-\[Kappa]val)/2*\[Epsilon]val*))];      (*__________________TEST___________________*)
Return[IncomingAmplitudeBtmp//ToMSTPrecision]
];
];

TransmittedAmplitudeB::usage="Calculated Amplitude \!\(\*SuperscriptBox[\(B\), \(trans\)]\) of asymptotic radial function";
TransmittedAmplitudeB[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{Kcoeff1, Kcoeff2, TransmittedAmplitudeBtmp, Aplus},
With[{\[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]],\[Epsilon]plus = ToMSTPrecision[\[Epsilon]plus[\[Omega],\[Chi],m]], \[Kappa]val = ToMSTPrecision[\[Kappa][\[Chi]]],coeff= Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]]},

TransmittedAmplitudeBtmp = ((\[Epsilon]val*\[Kappa]val)/\[Omega])^(2*s)*Exp[I*\[Kappa]val*\[Epsilon]plus(1+2*Log[\[Kappa]val]/(1+\[Kappa]val))]*\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-$MSTMaxN\)\), \($MSTMaxN\)]\(coeff[n]\)\);
Return[TransmittedAmplitudeBtmp//ToMSTPrecision]
];
];

ReflectedAmplitudeB::usage="Calculated Amplitude \!\(\*SuperscriptBox[\(B\), \(ref\)]\) of asymptotic radial function";
ReflectedAmplitudeB[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{Kcoeff1, Kcoeff2, ReflectedAmplitudeBtmp, Aminus},With[{\[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]],\[Epsilon]plus = ToMSTPrecision[\[Epsilon]plus[\[Omega],\[Chi],m]], \[Kappa]val = ToMSTPrecision[\[Kappa][\[Chi]]],coeff= Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]]},
Kcoeff1 = MatchingKcoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]];
Kcoeff2 = MatchingKcoeff[-\[Nu]-1,\[Omega],\[Chi],m,s,\[Lambda]];
Aminus = A\[Nu]minus[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]];
ReflectedAmplitudeBtmp = \[Omega]^(-1-2*s)*(Kcoeff1 + I*Exp[I*\[Pi]*\[Nu]]*Kcoeff2)*Aminus*Exp[I(\[Epsilon]val*Log[\[Epsilon]val] - (1-\[Kappa]val)/2*\[Epsilon]val)];
Return[ReflectedAmplitudeBtmp//ToMSTPrecision]
];
];

TransmittedAmplitudeC::usage="Calculated Amplitude \!\(\*SuperscriptBox[\(C\), \(trans\)]\) of asymptotic radial function";
TransmittedAmplitudeC[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_]:=
Block[{Kcoeff1, Kcoeff2, TransmittedAmplitudeCtmp, Aminus},
With[{\[Epsilon]val = ToMSTPrecision[\[Epsilon][\[Omega]]],\[Epsilon]plus = ToMSTPrecision[\[Epsilon]plus[\[Omega],\[Chi],m]], \[Kappa]val = ToMSTPrecision[\[Kappa][\[Chi]]],coeff= Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]]},
Aminus = A\[Nu]minus[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]];
TransmittedAmplitudeCtmp=\[Omega]^(-1-2*s)*Aminus*Exp[I(\[Epsilon]val*Log[\[Epsilon]val]-(1-\[Kappa]val)/2*\[Epsilon]val)];
Return[TransmittedAmplitudeCtmp//ToMSTPrecision];
];
];


Options[OuterRPlus] = {SeriesSign->1};
OuterRPlus[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_, OptionsPattern[]][r_]:=
Block[{OuterRPlustmp},
With[{\[Epsilon]val = \[Epsilon][\[Omega]],\[Epsilon]plusval = \[Epsilon]plus[\[Omega],\[Chi],m],\[Kappa]val = \[Kappa][\[Chi]], coeff = Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]]},
OuterRPlustmp = 2^\[Nu]*Exp[-\[Pi]*\[Epsilon]val]*Exp[I*\[Pi]*(\[Nu]+1-s)]*Gamma[\[Nu]+1-s+I*\[Epsilon]val]/Gamma[\[Nu]+1+s-I*\[Epsilon]val]*Exp[-I*ztilde[r,\[Omega],\[Chi]]]*ztilde[r,\[Omega],\[Chi]]^(\[Nu]+I*\[Epsilon]plusval) (ztilde[r,\[Omega],\[Chi]]-\[Epsilon]val*\[Kappa]val)^(-s-I*\[Epsilon]plusval)*
\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-$MSTMaxN\)\), \($MSTMaxN\)]\((
\*SuperscriptBox[\(I\), \(n\)]*coeff[OptionValue[SeriesSign]*n]*
\*SuperscriptBox[\((2*ztilde[r, \[Omega], \[Chi]])\), \(n\)]*IrregularConfluentHGF[n + \[Nu] + 1 - s + I*\[Epsilon]val, \ 2*n + 2*\[Nu] + 2, \ 2*I*ztilde[r, \[Omega], \[Chi]]])\)\);
Return[OuterRPlustmp//ToMSTPrecision]
]
];

Options[OuterRMinus] = {SeriesSign->1};
OuterRMinus[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_, OptionsPattern[]][r_]:=
Block[{OuterRMinustmp},
With[{\[Epsilon]val = \[Epsilon][\[Omega]],\[Epsilon]plusval = \[Epsilon]plus[\[Omega],\[Chi],m],\[Kappa]val = \[Kappa][\[Chi]], coeff = Function[{n},SeriesCoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][n]]},
OuterRMinustmp =  2^\[Nu]*Exp[-\[Pi]*\[Epsilon]val]*Exp[-I*\[Pi]*(\[Nu]+1+s)]Exp[I*ztilde[r,\[Omega],\[Chi]]]*ztilde[r,\[Omega],\[Chi]]^(\[Nu]+I*\[Epsilon]plusval) (ztilde[r,\[Omega],\[Chi]]-\[Epsilon]val*\[Kappa]val)^(-s-I*\[Epsilon]plusval)*\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-$MSTMaxN\)\), \($MSTMaxN\)]\((
\*SuperscriptBox[\(I\), \(n\)]*
\*FractionBox[\(Pochhammer[\[Nu] + 1 + s - I*\[Epsilon]val, n]\), \(Pochhammer[\[Nu] + 1 - s + I*\[Epsilon]val, n]\)] coeff[OptionValue[SeriesSign]*n]*
\*SuperscriptBox[\((2  ztilde[r, \[Omega], \[Chi]])\), \(n\)]*IrregularConfluentHGF[n + \[Nu] + 1 + s - I*\[Epsilon]val, \ 2*n + 2*\[Nu] + 2, \ \(-2\)*I*ztilde[r, \[Omega], \[Chi]]])\)\);
Return[OuterRMinustmp//ToMSTPrecision]
];
];

OuterRC::usage="Asymptotic expansion of radial function in terms of Coulombic wavefunctions satisfying the ingoing boundary condition. 
See Eq. (152) of Sasaki et. al";
Options[OuterRC] = Options[OuterRPlus];
OuterRC[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_, OptionsPattern[]][r_]:=OuterRPlus[\[Nu],\[Omega],\[Chi],m,s,\[Lambda], SeriesSign->OptionValue[SeriesSign]][r]+OuterRMinus[\[Nu],\[Omega],\[Chi],m,s,\[Lambda],SeriesSign->OptionValue[SeriesSign]][r]

OuterRIncoming::usage="Asymptotic expansion of radial function in terms of Coulombic wavefunctions satisfying the ingoing boundary condition and using matching of solution to Horizon expansion. 
See Eq. (166) of Sasaki et. al";
Options[OuterRIncoming] = Options[OuterRPlus];
OuterRIncoming[\[Nu]_,\[Omega]_,\[Chi]_,m_,s_,\[Lambda]_, OptionsPattern[]][r_]:=
MatchingKcoeff[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]]*OuterRC[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][r] + MatchingKcoeff[-\[Nu]-1,\[Omega],\[Chi],m,s,\[Lambda]]*OuterRC[-\[Nu]-1,\[Omega],\[Chi],m,s,\[Lambda]][r];


RNumerical::usage="Numerical solution of Homogeneous radial Teukolsky equation
RNumerical[\[Nu]_,\[Omega]_,\[Chi]_,l_,m_,s_,\[Lambda]_, RSTOP_]
RSTOP parameter is the radius at which to halt the numerical integration.
";
RNumerical[\[Nu]_,\[Omega]_,\[Chi]_,l_,m_,s_,\[Lambda]_, RSTOP_]:=
Block[{
rinx, 
Radial,
\[CapitalDelta],
K,
xinr,
DiffEQX,
DiffEQr,
rstart,
rstop,
Frob,
FrobSeriesEquation,
FrobSeriesEquationSet,
FrobSolution,
FrobinR,
InitialConditions,
RSolution,
SampleRadius,
Normalization,
Solution,
NormalizedRSolution

},
With[{
\[Lambda]val = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Chi]*\[Omega]],
\[CapitalDelta] =  r^2-2*r + \[Chi]^2,
K = (r^2+\[Chi]^2)*\[Omega] -m*\[Chi]},

DiffEQr=((\[CapitalDelta])^-s*D[\[CapitalDelta]^(s+1)*D[#,r],r] + ((K^2-2*I*s*(r-1)*K)/\[CapitalDelta]+4*I*s*\[Omega]*r -\[Lambda])#)& ;

(*Frobenius expansion near horizon to generate initial conditions*)
rinx = X*(rplus[\[Chi]]-rminus[\[Chi]])+rplus[\[Chi]];
xinr =( r-rplus[\[Chi]])/(rplus[\[Chi]]-rminus[\[Chi]]);
With[{\[CapitalDelta]x = \[CapitalDelta]/.r->rinx,Kx = K/.r->rinx, jacobian = D[xinr,r]},
DiffEQX = (\[CapitalDelta]x^-s*D[\[CapitalDelta]x^(s+1),X]*D[#,X]*jacobian^2 + \[CapitalDelta]x*jacobian^2*D[#,X,X] + ((Kx^2-2*I*s*(rinx-1)Kx)/\[CapitalDelta]x+4*I*s*\[Omega]*rinx - \[Lambda])#)&;
];

rstart = (rplus[\[Chi]]+10^-5)//ToMSTPrecision;
rstop = RSTOP//ToMSTPrecision;
Frob[X_]:=X^(-I*kk) (1 + C1*(X)+C2*(X)^2);
FrobSeriesEquation= Series[ExpandAll[X^(I kk+2) *DiffEQX@Frob[X]], {X,0,3}]==0;
FrobSeriesEquationSet = FrobSeriesEquation//LogicalExpand;
FrobSolution = NSolve[FrobSeriesEquationSet, {kk, C1, C2}][[1]];


FrobinR[r_]:=Frob[(r-rplus[\[Chi]])/(rplus[\[Chi]]-rminus[\[Chi]])]/.FrobSolution//ToMSTPrecision;
InitialConditions = {Limit[FrobinR[r], r->rstart], Limit[D[FrobinR[r],r], r->rstart]}//ToMSTPrecision;

RSolution = NDSolve[
{ToMSTPrecision[DiffEQr@Radial[r]==0], Radial[rstart]==InitialConditions[[1]], Derivative[1][Radial][rstart]==InitialConditions[[2]]}, 
Radial, 
{r,rstart, rstop}, 
WorkingPrecision->$MachinePrecision,
Method->{"TimeIntegration"->{"ExplicitRungeKutta","DifferenceOrder"->8}},
 MaxSteps->10^9]//First//First//Last;(*Pick out the InterpolatingFunction from the return value of NDSolve*)

SampleRadius = If[rstop<10, rstop,10];
(*Numerical integration doesnt determine normalization of radial function. Use analytic result for that*)

Block[{$MSTMaxN=10}, (*Using lower value of MaxN results in much faster computation time with little decrease in required precision*)
Normalization = OuterRIncoming[\[Nu],\[Omega],\[Chi],m,s,\[Lambda]][SampleRadius]/RSolution[SampleRadius]//SetPrecision[#, 100]&;
];
NormalizedRSolution = OperateOnInterpolationValues[RSolution, (#*Normalization&)];
MSTDebug = <|"FrobSeriesEquationSet"->FrobSeriesEquationSet,"UnnormalizedRSolution"->RSolution, "Normalization"->Normalization, "IC"-> InitialConditions, "rstart"->rstart, "rstop"->rstop, "FrobeniusSolution"->FrobSolution, "DiffEQx"->DiffEQX, "DiffEQr"->DiffEQr|>;
Return[NormalizedRSolution];
]
]
