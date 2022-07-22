(* ::Package:: *)

(* ::Chapter:: *)
(*Setup*)


TeukolskySolver;

SolutionPath = $FKKSRoot <> "Solutions/";

If["xActSetup" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]
];

If["EnergyMomentumwl" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "EnergyMomentum.wl"}]]
]

If["MSTSolver" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "MSTSolver.wl"}]]
]
Get["KerrWithProca`"]
Get["SpinWeightedSpheroidalHarmonics`"]
Get["Teukolsky`"]


(*Setup basic functions and parameters*)
\[Theta]\[Epsilon]=10^-5;

Teukrho = 1/(r[]-I*a*Cos[\[Theta][]]);
Teukrhob = 1/(r[]+I*a*Cos[\[Theta][]]);
Kerr\[CapitalDelta] = r[]^2-2*M*r[]+a^2;
Kerr\[CapitalSigma] = r[]^2+a^2*Cos[\[Theta]]^2;

L[s_][expr_]:=
	Block[{tmp},
		With[{EXPRESSION = FromxActVariables[expr]},
			tmp = (D[#,\[Theta]]- I/Sin[\[Theta]] D[#,\[Phi]]-I*a*Sin[\[Theta]]D[#,t]+s*Cot[\[Theta]]*#)&@EXPRESSION
			];
		ToxActVariables[tmp]
	];

Jp[expr_]:=
	Block[{tmp},
		With[{EXPRESSION = FromxActVariables[expr], \[CapitalDelta] = FromxActVariables[Kerr\[CapitalDelta]]},
			tmp=(D[#,r]-1/\[CapitalDelta] ((r^2+a^2)D[#,t]+a*D[#,\[Phi]]))&@EXPRESSION
		];
	ToxActVariables[tmp]
	];


(* ::Title:: *)
(*Functions*)


(*Define the Kinnersley Tetrad*)
Block[{nn,mm,mmbar},
	DefTensor[nn[\[Zeta]], \[ScriptCapitalM], PrintAs->"n"];
	DefTensor[mm[\[Zeta]], \[ScriptCapitalM], PrintAs->"m"];
	DefTensor[mmbar[\[Zeta]], \[ScriptCapitalM], PrintAs->"\!\(\*OverscriptBox[\(m\), \(_\)]\)"];
		Block[{Print},
			ComponentValue[ComponentArray[nn[{\[Zeta],sphericalchart}]], {r[]^2+a^2, -(r[]^2-2*M*r[] + a^2),0,a}/(2*(r[]^2+a^2*Cos[\[Theta][]]^2))];
			ComponentValue[ComponentArray[mm[{\[Zeta],sphericalchart}]], {I*a*Sin[\[Theta][]],0,1,I/Sin[\[Theta][]]}/(Sqrt[2]*(r[]+I*a*Cos[\[Theta][]]))];
			ComponentValue[ComponentArray[mmbar[{\[Zeta],sphericalchart}]], {-I*a*Sin[\[Theta][]],0,1,-I/Sin[\[Theta][]]}/(Sqrt[2]*(r[]-I*a*Cos[\[Theta][]]))];
		];
	KinnersleyN = ToCTensor[nn, {sphericalchart}]/.TensorValues[nn];
	KinnersleyM = ToCTensor[mm, {sphericalchart}]/.TensorValues[mm];
	KinnersleyMBar = ToCTensor[mmbar, {sphericalchart}]/.TensorValues[mmbar];
]


(*Calculate various projections of EM-tensor onto kinnersley tetrad*)
Options[Tnn]={AsSymbolic->False, AsInterpolatingFunction->False, AsOptimized->False, ToCompiled->False, WithProperties->False};
Tnn[solution_, OptionsPattern[]]:=
Block[{TLL,tmp, tmp$, tmp$OE,tmp$$, tmp$IP, tmp$r, ww$, mm$, res, rmin$, rmax$,Messengernn},
With[{TnnFilePath = FileNameJoin[{$FKKSRoot, "Expressions", "FKKSEnergyMomentumNN.mx"}]},
If[FileExistsQ[TnnFilePath],
tmp = Import[TnnFilePath];,

$Messenger="Forming analytic expression for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(nn\)]\)";
TLL = Import[$FKKSRoot<>"Expressions/FKKSEnergyMomentumRealCTensor.mx"];
tmp = FromxActVariables[TLL[-\[Zeta],-\[Xi]]KinnersleyN[\[Zeta]]KinnersleyN[\[Xi]]];
Export[TnnFilePath, tmp];
];
];
$Messenger="Mapping solution onto \!\(\*SubscriptBox[\(\[CapitalTau]\), \(nn\)]\)";
res = Evaluate[tmp/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution]];

If[OptionValue[AsSymbolic], 
Return[res];
];

$Messenger="Generating Optimized function for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(nn\)]\) with properties: Compiled->"<>ToString[OptionValue[ToCompiled]]<>" WithProperties->"<>ToString[OptionValue[WithProperties]];
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->OptionValue[ToCompiled],
WithProperties->OptionValue[WithProperties]];

If[OptionValue[AsOptimized],
Return[tmp$]
];

If[OptionValue[AsInterpolatingFunction], 
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, nn] = Subscript[(T^1), nn][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), nn][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), nn][r,\[Theta]]*)
Block[{rdom, Dmat, CoefficientSolution, reprrule, coeff1exp, coeff2exp, coeff3exp,coeff1, coeff2, coeff3, TnnDecomposed},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&;
Dmat[phi1_, phi2_,phi3_] := {{Cos[phi1],Sin[phi1],1},{Cos[phi2],Sin[phi2],1},{Cos[phi3], Sin[phi3], 1}};

With[{rpoints = solution["Solution", "R"]["Domain"]//First//Last, \[Theta]points = 200,\[Theta]sampling = {0,\[Pi]/2, \[Pi]}, \[Omega]value = solution["Solution", "\[Omega]"]//Re//Evaluate,
 mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1]
T[0,r,\[Theta],p2]
T[0,r,\[Theta],p3]

\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1

\[NoBreak])(\[NoBreak]A
B
C

\[NoBreak])
*)
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Theta]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Theta]sampling[[1]]/2], Bb->tmp$[0,r,\[Theta],\[Theta]sampling[[2]]/2], Cc->tmp$[0,r,\[Theta],\[Theta]sampling[[3]]/2]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule], ToCompiled->True];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule], ToCompiled->True];
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule], ToCompiled->True];

coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff3 = GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points},DensityFunctions->{(#^4&), Identity}];
On[CompiledFunction::cfsa];

TnnDecomposed = {t,r,\[Theta],\[Phi]}|->Evaluate[coeff1[r,\[Theta]]*Cos[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff2[r,\[Theta]]*Sin[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff3[r,\[Theta]]]
];(*end With statement*)

$Messenger = "Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function ... Done.";
Return[TnnDecomposed]
];(*end Block statement*)
]; (*end If statement*)

(*Default Return value*)
Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[{t,r,\[Theta],\[Phi]}|->FromxActVariables[res]]];
]
];



Options[Tmbmb]=Options[Tnn];
Tmbmb[solution_,  OptionsPattern[]]:=
Block[{TLL,tmp, tmp$, tmp$OE, tmp$IP, tmp$r, ww$, mm$, res,rmin$,rmax$,Messengermbmb},
With[{TmbmbFilePath = FileNameJoin[{$FKKSRoot, "Expressions", "FKKSEnergyMomentumMbarMbar.mx"}]},
If[FileExistsQ[TmbmbFilePath],
tmp = Import[TmbmbFilePath];,

$Messenger="Forming analytic expression for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\)";
TLL = Import[$FKKSRoot<>"Expressions/FKKSEnergyMomentumRealCTensor.mx"];
tmp = FromxActVariables[TLL[-\[Zeta],-\[Xi]]KinnersleyMBar[\[Zeta]]KinnersleyMBar[\[Xi]]];
Export[TmbmbFilePath, tmp];
];
];
$Messenger="Mapping solution onto \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\)";
res = Evaluate[tmp/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution]];

If[OptionValue[AsSymbolic], 
Return[res];
];

$Messenger="Generating Optimized function for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\) with properties: Compiled->"<>ToString[OptionValue[ToCompiled]]<>" WithProperties->"<>ToString[OptionValue[WithProperties]];
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->OptionValue[ToCompiled],
WithProperties->OptionValue[WithProperties]];

If[OptionValue[AsOptimized],
Return[tmp$]
];

If[OptionValue[AsInterpolatingFunction], 
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, nn] = Subscript[(T^1), nn][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), nn][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), nn][r,\[Theta]]*)
Block[{rdom, Dmat, CoefficientSolution, reprrule, coeff1exp, coeff2exp, coeff3exp,coeff1, coeff2, coeff3, TnnDecomposed},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&;
Dmat[phi1_, phi2_,phi3_] := {{Cos[phi1],Sin[phi1],1},{Cos[phi2],Sin[phi2],1},{Cos[phi3], Sin[phi3], 1}};
With[{rpoints = solution["Solution", "R"]["Domain"]//First//Last, \[Theta]points = 200,\[Theta]sampling = {0,\[Pi]/2, \[Pi]}, \[Omega]value = solution["Solution", "\[Omega]"]//Re, mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1]
T[0,r,\[Theta],p2]
T[0,r,\[Theta],p3]

\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1

\[NoBreak])(\[NoBreak]A
B
C

\[NoBreak])
*)
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Theta]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Theta]sampling[[1]]/2], Bb->tmp$[0,r,\[Theta],\[Theta]sampling[[2]]/2], Cc->tmp$[0,r,\[Theta],\[Theta]sampling[[3]]/2]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule], ToCompiled->True];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule], ToCompiled->True];
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule], ToCompiled->True];

coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff3 = GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points},DensityFunctions->{(#^4&), Identity}];
On[CompiledFunction::cfsa];

TnnDecomposed = {t,r,\[Theta],\[Phi]}|->Evaluate[coeff1[r,\[Theta]]*Cos[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff2[r,\[Theta]]*Sin[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff3[r,\[Theta]]]
];(*end With statement*)

$Messenger = "Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\) Interpolating Function ... Done.";
Return[TnnDecomposed]
];(*end Block statement*)
]; (*end If statement*)

Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[FromxActVariables[res]]];
]
]

Options[Tmbn]=Options[Tnn];
Tmbn[solution_, OptionsPattern[]]:=
Block[{TLL,tmp, tmp$, tmp$OE, tmp$IP, tmp$r, ww$, mm$, res,rmin$, rmax$,Messengermbn},
With[{TmbnFilePath = FileNameJoin[{$FKKSRoot, "Expressions", "FKKSEnergyMomentumMbarN.mx"}]},
If[FileExistsQ[TmbnFilePath],
tmp = Import[TmbnFilePath];,

$Messenger="Forming analytic expression for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] n\)]\)";
TLL = Import[$FKKSRoot<>"Expressions/FKKSEnergyMomentumRealCTensor.mx"];
tmp = FromxActVariables[TLL[-\[Zeta],-\[Xi]]KinnersleyMBar[\[Zeta]]KinnersleyN[\[Xi]]];
Export[TmbnFilePath, tmp];
];
];
$Messenger="Mapping solution onto \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] n\)]\)";
res = Evaluate[tmp/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution]];

If[OptionValue[AsSymbolic], 
Return[res];
];

$Messenger="Generating Optimized function for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] n\)]\) with properties: Compiled->"<>ToString[OptionValue[ToCompiled]]<>" WithProperties->"<>ToString[OptionValue[WithProperties]];
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->OptionValue[ToCompiled],
WithProperties->OptionValue[WithProperties]];

If[OptionValue[AsOptimized],
Return[tmp$]
];

If[OptionValue[AsInterpolatingFunction], 
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(\*OverscriptBox[\(m\), \(_\)] n\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, nn] = Subscript[(T^1), nn][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), nn][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), nn][r,\[Theta]]*)
Block[{rdom, Dmat, CoefficientSolution, reprrule, coeff1exp, coeff2exp, coeff3exp,coeff1, coeff2, coeff3, TnnDecomposed},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&;
Dmat[phi1_, phi2_,phi3_] := {{Cos[phi1],Sin[phi1],1},{Cos[phi2],Sin[phi2],1},{Cos[phi3], Sin[phi3], 1}};
With[{rpoints = solution["Solution", "R"]["Domain"]//First//Last, \[Theta]points = 200,\[Theta]sampling = {0,\[Pi]/2, \[Pi]}, \[Omega]value = solution["Solution", "\[Omega]"]//Re, mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1]
T[0,r,\[Theta],p2]
T[0,r,\[Theta],p3]

\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1

\[NoBreak])(\[NoBreak]A
B
C

\[NoBreak])
*)
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Theta]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Theta]sampling[[1]]/2], Bb->tmp$[0,r,\[Theta],\[Theta]sampling[[2]]/2], Cc->tmp$[0,r,\[Theta],\[Theta]sampling[[3]]/2]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule], ToCompiled->True];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule], ToCompiled->True];
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule], ToCompiled->True];

coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}];
coeff3 = GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points},DensityFunctions->{(#^4&), Identity}];
On[CompiledFunction::cfsa];

TnnDecomposed = {t,r,\[Theta],\[Phi]}|->Evaluate[coeff1[r,\[Theta]]*Cos[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff2[r,\[Theta]]*Sin[-2*\[Omega]value*t + 2*mvalue*\[Phi]]+coeff3[r,\[Theta]]]
];(*end With statement*)

$Messenger = "Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(\*OverscriptBox[\(m\), \(_\)] n\)]\) Interpolating Function ... Done.";
Return[TnnDecomposed]
];(*end Block statement*)
]; (*end If statement*)

Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[FromxActVariables[res]]];
]
]

(*
\[ScriptCapitalT] = 2 (-1/2\[Rho]^8*Overscript[\[Rho], _]*Subscript[L, -1](\[Rho]^-4*Subscript[L, 0]((\[Rho]^-2*Overscript[\[Rho], _]^-1*Subscript[\[ScriptCapitalT], nn]))) - 1/(2*Sqrt[2])\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2*Subscript[L, -1](\[Rho]^-4*Overscript[\[Rho], _]^2Subscript[J, +]((\[Rho]^-2*Overscript[\[Rho], _]^-2*\[CapitalDelta]^-1*Subscript[\[ScriptCapitalT], Overscript[m, _]n]))) + -1/4\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2*Subscript[J, +] (\[Rho]^-4*Subscript[J, +] ((\[Rho]^-2*Overscript[\[Rho], _]*Subscript[\[ScriptCapitalT], Overscript[m, _]Overscript[m, _]])))  -  1/(2*Sqrt[2])\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2Subscript[J, +](\[Rho]^-4*Overscript[\[Rho], _]^2*\[CapitalDelta]^-1*Subscript[L, -1]((\[Rho]^-2*Overscript[\[Rho], _]^-2*Subscript[\[ScriptCapitalT], Overscript[m, _]n]))))
*)
TeukolskySource[solution_]:=Block[{\[ScriptCapitalT], B, BPrime,expr},
Monitor[
$Messenger="Beginning Teukolsky Source Calculation...";
With[{
tmbmb = Tmbmb[solution, AsInterpolatingFunction->True],
tmbn = Tmbn[solution, AsInterpolatingFunction->True],
tnn = Tnn[solution, AsInterpolatingFunction->True], 
\[Rho] = Teukrho,
\[Rho]b = Teukrhob, 
\[CapitalDelta] = Kerr\[CapitalDelta]
},
B = -1/2 \[Rho]^8*\[Rho]b*L[-1][\[Rho]^-4*L[0][(\[Rho]^-2*\[Rho]b^-1*tnn[t,r,\[Theta],\[Phi]])]]-1/(2*Sqrt[2]) \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*L[-1][\[Rho]^-4*\[Rho]b^2 Jp[(\[Rho]^-2*\[Rho]b^-2*\[CapitalDelta]^-1*tmbn[t,r,\[Theta],\[Phi]])]];
BPrime = -1/4 \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*Jp[\[Rho]^-4*Jp[(\[Rho]^-2*\[Rho]b*tmbmb[t,r,\[Theta],\[Phi]])]]-1/(2*Sqrt[2]) \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*Jp[\[Rho]^-4*\[Rho]b^2*\[CapitalDelta]^-1*L[-1][(\[Rho]^-2*\[Rho]b^-2*tmbn[t,r,\[Theta],\[Phi]])]];
expr = (2(B+BPrime))//FromxActVariables//ApplySolutionSet[solution];
\[ScriptCapitalT] = OptimizedFunction[{t,r,\[Theta],\[Phi]},expr];
Return[\[ScriptCapitalT]];
],
$Messenger
]
]

(*
Subscript[\[ScriptCapitalT], lm\[Omega], integrand] = 2 2/Sqrt[2*\[Pi]]*\[Rho]^-5*Overscript[\[Rho], _]^-1 \[ScriptCapitalT]Subscript[(r,\[Theta]) , -2]Subscript[S^a\[Omega], lm](\[Theta])
*)
TeukolskySourceModalIntegrand[solution_,l_,m_]:=
With[{Gamma = solution[["Parameters","\[Chi]"]]*2*solution[["Solution","\[Omega]"]]},
OptimizedFunction[{r,\[Theta]},
Evaluate[ApplySolutionSet[solution][FromxActVariables[2/Sqrt[2*\[Pi]]*Teukrho^-5*Teukrhob^-1*TeukolskySource[solution][0,r,\[Theta],0]*SpinWeightedSpheroidalHarmonicS[-2,l,m,Gamma,\[Theta],0]]]]
]
]

(*
Subscript[\[ScriptCapitalT], lm\[Omega]] = \[Integral]d\[Theta]d\[Phi] sin(\[Theta]) Subscript[\[ScriptCapitalT], lm\[Omega], integrand]
*)
TeukolskySourceModal[solution_,l_,m_]:=
Block[{SourceIntegrand,SourceIntegrandFunction,IntegrationFunction,TlmwOnMesh,TlmwData, TlmwInterpolation,radialMesh,rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#, solution["Parameters", "\[Chi]"]]&},
radialMesh = Table[r, {r,rdom[[1]], rdom[[2]], 0.2}];
$Messenger = "Generating Teukoslky Integrand";
SourceIntegrand = TeukolskySourceModalIntegrand[solution,l,m];
SourceIntegrandFunction[x_?NumericQ,\[Theta]_?NumericQ]:=SourceIntegrand[x,\[Theta]];
IntegrationFunction = {x}|->NIntegrate[SourceIntegrandFunction[x,\[Theta]]*Sin[\[Theta]], {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon]}, WorkingPrecision->10, MaxRecursion->100];
$Messenger = "Evaluating integration at each radial mesh node";
DistributeDefinitions[SourceIntegrand, SourceIntegrandFunction,SpinWeightedSpheroidalHarmonicS];
TlmwOnMesh = ParallelMap[IntegrationFunction, radialMesh];
TlmwData = Thread[{radialMesh,TlmwOnMesh}];
$Messenger = "Generating Interpolation function";
TlmwInterpolation = Interpolation[TlmwData];
Return[TlmwInterpolation];
$Messenger= "Done."
]

(*
Subscript[Z^\[Infinity], lm\[Omega]] = 1/(2i\[Omega] Subscript[B^inc, lm\[Omega]])\!\(
\*SubscriptBox[
SuperscriptBox[\(\[Integral]\), \(\[Infinity]\)], 
SubscriptBox[\(r\), \(+\)]]\(
\*FractionBox[\(
\*SubscriptBox[\(\[ScriptCapitalT]\), \(lm\[Omega]\)]\((r)\)
\*SubscriptBox[
SuperscriptBox[\(R\), \(in\)], \(lm\[Omega]\)]\((r)\)\), \(
\*SuperscriptBox[\(\[CapitalDelta]\), \(2\)]\((r)\)\)]dr\)\)
*)

(*!!!!! There has to be a better way of obtaining the Incident amplitude*)

TeukolskyZInfinity[ProcaSolution_, TeukolskyLM\[Omega]_,l_,m_]:=
Block[{
DeltaFunction,
TeukRin, 
integration,
Amplitude,
Zcoefficient,
\[CapitalDelta]=Evaluate[Kerr\[CapitalDelta]//FromxActVariables//ApplySolutionSet[ProcaSolution]], 
rdom = ProcaSolution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#, ProcaSolution["Parameters", "\[Chi]"]]&,
\[Chi]v = SetPrecision[ProcaSolution["Parameters", "\[Chi]"],60],
\[Omega]v = SetPrecision[2*ProcaSolution["Solution", "\[Omega]"],60]
},
DeltaFunction = {r}|->Evaluate[\[CapitalDelta]];
TeukRin = Evaluate@TeukolskyRadial[-2,l,m,\[Chi]v,\[Omega]v , BoundaryConditions->"In", Method->{"NumericalIntegration", "Domain"->rdom}];
Amplitude = TeukolskyRadial[-2,l,m,\[Chi]v,\[Omega]v , Method->"MST"]["In"]["Amplitudes"]["Incidence"];
integration = NIntegrate[TeukolskyLM\[Omega][r]*TeukRin[r]/DeltaFunction[r]^2, {r,rdom[[1]], rdom[[2]]}];
Zcoefficient = integration/(2*I*2*ProcaSolution["Solution", "\[Omega]"]*Amplitude);
Return[Zcoefficient]
]

Options[EnergyFlux]={TeukolskyTlmw->Automatic};
EnergyFlux[ProcaSolution_,OptionsPattern[]]:=
Block[{Zcoeff,
TeukModal,
EnFlux,
teuk\[Omega] = 2*Re[ProcaSolution["Solution", "\[Omega]"]]},
If[TrueQ[OptionValue[TeukolskyTlmw]==Automatic],
TeukModal = TeukolskySourceModal[ProcaSolution,2,2],
TeukModal = OptionValue[TeukolskyTlmw]
];
Zcoeff = TeukolskyZInfinity[ProcaSolution,TeukModal, 2,2];
EnFlux = Abs[Zcoeff]^2/(4*\[Pi]*teuk\[Omega]^2)
]
