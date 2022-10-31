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
];

If["MSTSolver" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "MSTSolver.wl"}]]
];

If["KerrWithProca" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]
];

If["HelperFunctions" \[NotElement] Names["Global`*"],
  Get[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]
];

Get["SpinWeightedSpheroidalHarmonics`"]
Get["Teukolsky`"]


(*Setup basic functions and parameters*)
\[Theta]\[Epsilon]=10^-3;
$TSCoefficent\[Theta]points=50; (*Number of points to use for \[Theta] argument in calculating Interpolating function for coefficients in Tnn, Tmn, and tmm*);
$TSNIntegrateMethod = {"GlobalAdaptive", "MaxErrorIncreases"->2000, Method->{"ClenshawCurtisRule", "SymbolicProcessing"->0, "Points"->20}};
$TSNIntegratePrecision = 5;
$TSNIntegrateMaxRecursion = 50;
$TSNIntegrateMinRecursion=0;
$TSMaxRadialPoints=300;
$TSPrecision = 25;ToTSPrecision[expr_]:=SetPrecision[expr, $TSPrecision];
$TSRuntimeOptions = {"CatchMachineOverflow"->False, "CatchMachineIntegerOverflow"->False, "CompareWithTolerance"->False, "EvaluateSymbolically"->False, "RuntimeErrorHandler"->Null, "WarningMessages"->True};

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
	
Lexpl[s_][expr_]:=
	Block[{tmp},
		With[{EXPRESSION = FromxActVariables[expr]},
			tmp = (D[#,\[Theta]] + ((2*m)/Sin[\[Theta]]-2*a*\[Omega]*Sin[\[Theta]]+s*Cot[\[Theta]])*#)&@EXPRESSION
			];
		ToxActVariables[tmp]
	];

Jpexpl[expr_]:=
	Block[{tmp},
		With[{EXPRESSION = FromxActVariables[expr], \[CapitalDelta] = FromxActVariables[Kerr\[CapitalDelta]]},
			tmp=(D[#,r]+(I*#)/\[CapitalDelta] ((r^2+a^2)*2*\[Omega] - a*2*m))&@EXPRESSION
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
			ComponentValue[ComponentArray[nn[{\[Zeta],sphericalchart}]], {r[]^2+a^2, -Kerr\[CapitalDelta],0,a}/(2*Kerr\[CapitalSigma])//Evaluate];
			ComponentValue[ComponentArray[mm[{\[Zeta],sphericalchart}]], {I*a*Sin[\[Theta][]],0,1,I/Sin[\[Theta][]]}/(Sqrt[2]*(r[]+I*a*Cos[\[Theta][]]))];
			ComponentValue[ComponentArray[mmbar[{\[Zeta],sphericalchart}]], {-I*a*Sin[\[Theta][]],0,1,-I/Sin[\[Theta][]]}/(Sqrt[2]*(r[]-I*a*Cos[\[Theta][]]))];
		];
	KinnersleyN = ToCTensor[nn, {sphericalchart}]/.TensorValues[nn];
	KinnersleyM = ToCTensor[mm, {sphericalchart}]/.TensorValues[mm];
	KinnersleyMBar = ToCTensor[mmbar, {sphericalchart}]/.TensorValues[mmbar];
]


(*Calculate various projections of EM-tensor onto kinnersley tetrad*)
Options[Tnn]={AsSymbolic->False, AsInterpolatingFunction->False, AsOptimized->False, ToCompiled->False, WithProperties->False, ExtractTimeDependence->False};
Tnn[solution_?AssociationQ, OptionsPattern[]]:=
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

If[OptionValue[AsOptimized],
	$Messenger="Generating Optimized function for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(nn\)]\) with properties: Compiled->"<>ToString[OptionValue[ToCompiled]]<>" WithProperties->"<>ToString[OptionValue[WithProperties]];
	tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->OptionValue[ToCompiled], WithProperties->OptionValue[WithProperties]];
	Return[tmp$]
];

If[OptionValue[AsInterpolatingFunction], 
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, nn] = Subscript[(T^1), nn][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), nn][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), nn][r,\[Theta]]*)
Block[{rdom, coefftimestart,coefftimestop,Dmat, radialdomain, CoefficientSolution, reprrule,  coeff1exp, coeff2exp, coeff3exp, coeff1,coeff2, coeff3, Decomposed,Prefactor},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&//ToTSPrecision;
Dmat[phi1_, phi2_, phi3_] := {{Cos[2*solution["Parameters", "m"]*phi1],Sin[2*solution["Parameters", "m"]*phi1],1},{Cos[2*solution["Parameters", "m"]*phi2],Sin[2*solution["Parameters", "m"]*phi2],1}, {Cos[2*solution["Parameters", "m"]*phi3],Sin[2*solution["Parameters", "m"]*phi3],1}};
radialdomain = solution["Solution", "R"]["Domain"]//First//ToTSPrecision;
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
With[{rpoints = If[radialdomain[[-1]]<$TSMaxRadialPoints, radialdomain[[-1]], $TSMaxRadialPoints+10*Log[radialdomain[[-1]]]], \[Theta]points = $TSCoefficent\[Theta]points, \[Phi]sampling = {0, \[Pi]/(4*solution["Parameters", "m"]), \[Pi]/(8*solution["Parameters", "m"])}, \[Omega]value = solution["Solution", "\[Omega]"]//Re,
 mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1],T[0,r,\[Theta],p2],T[0,r,\[Theta],p3]\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1
\[NoBreak])(\[NoBreak]A,B,C\[NoBreak])
*)
$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}], "Generating coefficient interpolating functions"}];
coefftimestart = AbsoluteTime[];
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Phi]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Phi]sampling[[1]]], Bb->tmp$[0,r,\[Theta],\[Phi]sampling[[2]]], Cc->tmp$[0,r,\[Theta],\[Phi]sampling[[3]]]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
If[Not[OptionValue[ExtractTimeDependence]],
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff3= GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
];
On[CompiledFunction::cfsa];
coefftimestop = AbsoluteTime[];
If[OptionValue[ExtractTimeDependence],
(*Numerical prefactor in agreement with nils code*)
(*expression comes from fourier transforming TNN, which has the functional form CC[r,\[Theta]]+AA[r,\[Theta]] Cos[2  (\[Phi]+\[Pi])-2 t \[Omega]]+BB[r,\[Theta]] Sin[2 (\[Phi]+\[Pi])-2 t \[Omega]]*)
Prefactor = Sqrt[2]*\[Pi]^(3/2);
Decomposed = {r$,\[Theta]$}|->Evaluate[Prefactor*(coeff1[r$,\[Theta]$]-I*coeff2[r$,\[Theta]$])];,

Decomposed = {t$,r$,\[Theta]$,\[Phi]$}|->Evaluate[coeff1[r$,\[Theta]$]*Cos[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff2[r$,\[Theta]$]*Sin[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff3[r$,\[Theta]$]];
];
];(*end With statement*)

$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function ... Done."}], "Generating coefficient interpolating functions... Done", "time to generate coefficient interpolating functions: "<>ToString[coefftimestop-coefftimestart, InputForm]}];
Return[Decomposed]
];(*end Block statement*)
]; (*end If statement*)

(*Default Return value*)
Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[{t,r,\[Theta],\[Phi]}|->FromxActVariables[res]]];
]
];


Options[Tmbmb]=Options[Tnn];
Tmbmb[solution_?AssociationQ,  OptionsPattern[]]:=
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

If[OptionValue[AsOptimized],
	$Messenger="Generating Optimized function for \!\(\*SubscriptBox[\(\[CapitalTau]\), \(\*OverscriptBox[\(m\), \(_\)] \*OverscriptBox[\(m\), \(_\)]\)]\) with properties: Compiled->"<>ToString[OptionValue[ToCompiled]]<>" WithProperties->"<>ToString[OptionValue[WithProperties]];
	tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->OptionValue[ToCompiled], WithProperties->OptionValue[WithProperties]];
	Return[tmp$]
];

If[OptionValue[AsInterpolatingFunction], 
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(mbmb\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, mbmb] = Subscript[(T^1), mbmb][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), mbmb][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), mbmb][r,\[Theta]]*)
Block[{rdom, coefftimestart,coefftimestop,Dmat, radialdomain, CoefficientSolution, reprrule, coeff1exp, coeff2exp, coeff3exp,coeff1, coeff2, coeff3, Decomposed},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&//ToTSPrecision;
Dmat[phi1_, phi2_, phi3_] := {{Cos[2*solution["Parameters", "m"]*phi1],Sin[2*solution["Parameters", "m"]*phi1],1},{Cos[2*solution["Parameters", "m"]*phi2],Sin[2*solution["Parameters", "m"]*phi2],1}, {Cos[2*solution["Parameters", "m"]*phi3],Sin[2*solution["Parameters", "m"]*phi3],1}};
radialdomain = solution["Solution", "R"]["Domain"]//First//ToTSPrecision;
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
With[{rpoints = If[radialdomain[[-1]]<$TSMaxRadialPoints, radialdomain[[-1]], $TSMaxRadialPoints+10*Log[radialdomain[[-1]]]], \[Theta]points = $TSCoefficent\[Theta]points, \[Phi]sampling = {0, \[Pi]/(4*solution["Parameters", "m"]), \[Pi]/(8*solution["Parameters", "m"])}, \[Omega]value = solution["Solution", "\[Omega]"]//Re,
 mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1],T[0,r,\[Theta],p2],T[0,r,\[Theta],p3]\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1
\[NoBreak])(\[NoBreak]A,B,C\[NoBreak])
*)
$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}], "Generating coefficient interpolating functions"}];
coefftimestart = AbsoluteTime[];
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Phi]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Phi]sampling[[1]]], Bb->tmp$[0,r,\[Theta],\[Phi]sampling[[2]]], Cc->tmp$[0,r,\[Theta],\[Phi]sampling[[3]]]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];

coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
If[Not[OptionValue[ExtractTimeDependence]],
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff3= GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
];
On[CompiledFunction::cfsa];
coefftimestop = AbsoluteTime[];
If[OptionValue[ExtractTimeDependence],
(*Numerical prefactor in agreement with nils code*)
(*expression comes from fourier transforming TNN, which has the functional form CC[r,\[Theta]]+AA[r,\[Theta]] Cos[2  (\[Phi]+\[Pi])-2 t \[Omega]]+BB[r,\[Theta]] Sin[2 (\[Phi]+\[Pi])-2 t \[Omega]]*)
Decomposed = {r$,\[Theta]$}|->Evaluate[Sqrt[2]*\[Pi]^(3/2)*(coeff1[r$,\[Theta]$]-I*coeff2[r$,\[Theta]$])];,

Decomposed = {t$,r$,\[Theta]$,\[Phi]$}|->Evaluate[coeff1[r$,\[Theta]$]*Cos[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff2[r$,\[Theta]$]*Sin[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff3[r$,\[Theta]$]];
];
];(*end With statement*)

$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(mbmb\)]\) Interpolating Function ... Done."}], "Generating coefficient interpolating functions... Done", "time to generate coefficient interpolating functions: "<>ToString[coefftimestop-coefftimestart, InputForm]}];
Return[Decomposed]
];(*end Block statement*)
]; (*end If statement*)

Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[FromxActVariables[res]]];
]
];


Options[Tmbn]=Options[Tnn];
Tmbn[solution_?AssociationQ, OptionsPattern[]]:=
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
$Messenger = Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(mbn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}];
(*Assume Expression of the formSubscript[T, mbn] = Subscript[(T^1), mbn][r,\[Theta]] Cos[-\[Omega] t + m \[Phi]] + Subscript[(T^2), mbn][r,\[Theta]] Sin[-\[Omega] t + m \[Phi]] + Subscript[(T^3), mbn][r,\[Theta]]*)
Block[{rdom, coefftimestart,coefftimestop,Dmat, radialdomain, CoefficientSolution, reprrule, coeff1exp, coeff2exp, coeff3exp,coeff1, coeff2, coeff3, Decomposed},
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&//ToTSPrecision;
Dmat[phi1_, phi2_, phi3_] := {{Cos[2*solution["Parameters", "m"]*phi1],Sin[2*solution["Parameters", "m"]*phi1],1},{Cos[2*solution["Parameters", "m"]*phi2],Sin[2*solution["Parameters", "m"]*phi2],1}, {Cos[2*solution["Parameters", "m"]*phi3],Sin[2*solution["Parameters", "m"]*phi3],1}};
radialdomain = solution["Solution", "R"]["Domain"]//First//ToTSPrecision;
tmp$ = OptimizedFunction[{t,r,\[Theta],\[Phi]}, res, ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
With[{rpoints = If[radialdomain[[-1]]<$TSMaxRadialPoints, radialdomain[[-1]], $TSMaxRadialPoints+10*Log[radialdomain[[-1]]]], \[Theta]points = $TSCoefficent\[Theta]points, \[Phi]sampling = {0, \[Pi]/(4*solution["Parameters", "m"]), \[Pi]/(8*solution["Parameters", "m"])}, \[Omega]value = solution["Solution", "\[Omega]"]//Re,
 mvalue = solution["Parameters", "m"]},
Off[CompiledFunction::cfsa];
(*We solve the expression 
Subscript[T, nn][0,r,\[Theta],Subscript[\[Phi], sample]]=A*Cos[2*m*Subscript[\[Phi], sample]]+B*Sin[2*m*Subscript[\[Phi], sample]]+C 
for the undetermined coefficients A, B, and C by choosing three values of Subscript[\[Phi], sample] and solving the resulting matrix equation 
(\[NoBreak]T[0,r,\[Theta],p1],T[0,r,\[Theta],p2],T[0,r,\[Theta],p3]\[NoBreak]) = (\[NoBreak]Cos[2m p1]	Sin[2m p1]	1
Cos[2m p2]	Sin[2m p2]	1
Cos[2m p3]	Sin[2m p3]	1
\[NoBreak])(\[NoBreak]A,B,C\[NoBreak])
*)
$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(nn\)]\) Interpolating Function", ProgressIndicator[Appearance->"Percolate"]}], "Generating coefficient interpolating functions"}];
coefftimestart = AbsoluteTime[];
CoefficientSolution = LinearSolve[Dmat[Sequence@@\[Phi]sampling], {Aa,Bb,Cc}];
reprrule = {Aa->tmp$[0,r,\[Theta],\[Phi]sampling[[1]]], Bb->tmp$[0,r,\[Theta],\[Phi]sampling[[2]]], Cc->tmp$[0,r,\[Theta],\[Phi]sampling[[3]]]};
coeff1exp = OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[1]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff2exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[2]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];

coeff1 = GenerateInterpolation[coeff1exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
coeff2= GenerateInterpolation[coeff2exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
If[Not[OptionValue[ExtractTimeDependence]],
coeff3exp =  OptimizedFunction[{r,\[Theta]},Evaluate[CoefficientSolution[[3]]/.reprrule//ToTSPrecision], ToCompiled->True, CompilationTarget->"WVM", RuntimeOptions->$TSRuntimeOptions];
coeff3= GenerateInterpolation[coeff3exp[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
];
On[CompiledFunction::cfsa];
coefftimestop = AbsoluteTime[];
If[OptionValue[ExtractTimeDependence],
(*Numerical prefactor in agreement with nils code*)
(*expression comes from fourier transforming TNN, which has the functional form CC[r,\[Theta]]+AA[r,\[Theta]] Cos[2  (\[Phi]+\[Pi])-2 t \[Omega]]+BB[r,\[Theta]] Sin[2 (\[Phi]+\[Pi])-2 t \[Omega]]*)
Decomposed = {r$,\[Theta]$}|->Evaluate[Sqrt[2]*\[Pi]^(3/2)*(coeff1[r$,\[Theta]$]-I*coeff2[r$,\[Theta]$])];,

Decomposed = {t$,r$,\[Theta]$,\[Phi]$}|->Evaluate[coeff1[r$,\[Theta]$]*Cos[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff2[r$,\[Theta]$]*Sin[-2*\[Omega]value*t$ + 2*mvalue*\[Phi]$]+coeff3[r$,\[Theta]$]];
];
];(*end With statement*)

$Messenger = Column[{Row[{"Generating \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(mbn\)]\) Interpolating Function ... Done."}], "Generating coefficient interpolating functions... Done", "time to generate coefficient interpolating functions: "<>ToString[coefftimestop-coefftimestart, InputForm]}];
Return[Decomposed]
];(*end Block statement*)
]; (*end If statement*)

Return[OptimizedFunction[{t,r,\[Theta],\[Phi]},  Evaluate[FromxActVariables[res]]];
]
];


(*
\[ScriptCapitalT] = 2 (-1/2\[Rho]^8*Overscript[\[Rho], _]*Subscript[L, -1](\[Rho]^-4*Subscript[L, 0]((\[Rho]^-2*Overscript[\[Rho], _]^-1*Subscript[\[ScriptCapitalT], nn]))) - 1/(2*Sqrt[2])\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2*Subscript[L, -1](\[Rho]^-4*Overscript[\[Rho], _]^2Subscript[J, +]((\[Rho]^-2*Overscript[\[Rho], _]^-2*\[CapitalDelta]^-1*Subscript[\[ScriptCapitalT], Overscript[m, _]n]))) + -1/4\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2*Subscript[J, +] (\[Rho]^-4*Subscript[J, +] ((\[Rho]^-2*Overscript[\[Rho], _]*Subscript[\[ScriptCapitalT], Overscript[m, _]Overscript[m, _]])))  -  1/(2*Sqrt[2])\[Rho]^8*Overscript[\[Rho], _]*\[CapitalDelta]^2Subscript[J, +](\[Rho]^-4*Overscript[\[Rho], _]^2*\[CapitalDelta]^-1*Subscript[L, -1]((\[Rho]^-2*Overscript[\[Rho], _]^-2*Subscript[\[ScriptCapitalT], Overscript[m, _]n]))))
*)
TeukolskySource[solution_?AssociationQ]:=Block[{\[ScriptCapitalT], B, BPrime,expr},
$Messenger="Beginning Teukolsky Source Calculation...";
With[{
tmbmb = Tmbmb[solution, AsInterpolatingFunction->True, ExtractTimeDependence->True],
tmbn = Tmbn[solution, AsInterpolatingFunction->True, ExtractTimeDependence->True],
tnn = Tnn[solution, AsInterpolatingFunction->True, ExtractTimeDependence->True], 
\[Rho] = Teukrho,
\[Rho]b = Teukrhob, 
\[CapitalDelta] = Kerr\[CapitalDelta]
},
B = -1/2 \[Rho]^8*\[Rho]b*Lexpl[-1][\[Rho]^-4*Lexpl[0][(\[Rho]^-2*\[Rho]b^-1*tnn[r,\[Theta]])]]-1/(2*Sqrt[2]) \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*Lexpl[-1][\[Rho]^-4*\[Rho]b^2 Jpexpl[(\[Rho]^-2*\[Rho]b^-2*\[CapitalDelta]^-1*tmbn[r,\[Theta]])]];
BPrime = -1/4 \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*Jpexpl[\[Rho]^-4*Jpexpl[(\[Rho]^-2*\[Rho]b*tmbmb[r,\[Theta]])]]-1/(2*Sqrt[2]) \[Rho]^8*\[Rho]b*\[CapitalDelta]^2*Jpexpl[\[Rho]^-4*\[Rho]b^2*\[CapitalDelta]^-1*Lexpl[-1][(\[Rho]^-2*\[Rho]b^-2*tmbn[r,\[Theta]])]];
expr = (2(B+BPrime))//FromxActVariables//ApplySolutionSet[solution];
\[ScriptCapitalT] = OptimizedFunction[{r,\[Theta]},Evaluate@expr];
Return[\[ScriptCapitalT]];
]
];


(*
Subscript[\[ScriptCapitalT], lm\[Omega], integrand] = 2 2/Sqrt[2*\[Pi]]*\[Rho]^-5*Overscript[\[Rho], _]^-1 \[ScriptCapitalT]Subscript[(r,\[Theta]) , -2]Subscript[S^a\[Omega], lm](\[Theta])
*)
Options[TeukolskySourceModalIntegrand]={Source->Null};
TeukolskySourceModalIntegrand[solution_?AssociationQ,l_?IntegerQ,m_?IntegerQ, OptionsPattern[]]:=
With[{Gamma = solution[["Parameters","\[Chi]"]]*2*solution[["Solution","\[Omega]"]]//Re},
(*Additional factor of Sqrt[2*\[Pi]] comes from normalization of SWSH package compared to normalization defined in http://link.springer.com/10.12942/lrr-2003-6*)
(*Note, integrand contains jacobian already*)
If[TrueQ[OptionValue[Source]==Null],
retval = OptimizedFunction[{r,\[Theta]},Evaluate[ApplySolutionSet[solution][FromxActVariables[2/(Teukrhob*Teukrho^5)*TeukolskySource[solution][r,\[Theta]]*Sqrt[2*\[Pi]]*SpinWeightedSpheroidalHarmonicS[-2,l,m,Gamma,\[Theta],0]*Sin[\[Theta]]]]]],
retval = OptimizedFunction[{r,\[Theta]},Evaluate[ApplySolutionSet[solution][FromxActVariables[2/(Teukrhob*Teukrho^5)*OptionValue[Source][r,\[Theta]]*Sqrt[2*\[Pi]]*SpinWeightedSpheroidalHarmonicS[-2,l,m,Gamma,\[Theta],0]*Sin[\[Theta]]]]]]
];
Return[retval];
];


(*
Subscript[\[ScriptCapitalT], lm\[Omega]] = \[Integral]d\[Theta]d\[Phi] sin(\[Theta]) Subscript[\[ScriptCapitalT], lm\[Omega], integrand]
*)
Options[TeukolskySourceModal]={Recalculate->False, SourceIntegrand->Null};
TeukolskySourceModal[solution_?AssociationQ,l_?IntegerQ,m_?IntegerQ, OptionsPattern[]]:=
Block[{SourceIntegrand,SourceIntegrandFunction,IntegrationFunction,TlmwOnMesh,TlmwData, TlmwInterpolation,radialMesh,rpoints,
rdom = solution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#, solution["Parameters", "\[Chi]"]]&
},
rpoints = If[rdom[[-1]]<$TSMaxRadialPoints, rdom[[-1]], $TSMaxRadialPoints+10*Log[rdom[[-1]]]];
radialMesh = Table[r, {r,rdom[[1]], rdom[[2]], (rdom[[2]]-rdom[[1]])/rpoints}];
If[OptionValue[Recalculate]||!KeyExistsQ[solution, "Derived"],
$Messenger = "Generating Teukolsky Integrand";
If[TrueQ[OptionValue[SourceIntegrand]==Null],
SourceIntegrand = TeukolskySourceModalIntegrand[solution,l,m];,
SourceIntegrand = OptionValue[SourceIntegrand];
];
SourceIntegrandFunction[x_?NumericQ,\[Theta]_?NumericQ]:=SourceIntegrand[x,\[Theta]];
IntegrationFunction[x_?NumericQ]:=NIntegrate[SourceIntegrandFunction[x,\[Theta]], {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon]}, WorkingPrecision->10, MaxRecursion->100];
$Messenger = "Evaluating integration at each radial mesh node";
DistributeDefinitions[SourceIntegrand, SourceIntegrandFunction,SpinWeightedSpheroidalHarmonicS];
TlmwData = ParallelMap[({#,IntegrationFunction[#]}&), radialMesh];
(*TlmwData = Thread[{radialMesh,TlmwOnMesh}];*)
$Messenger = "Generating Interpolation function";
TlmwInterpolation = Interpolation[TlmwData];
Return[TlmwInterpolation];
$Messenger= "Done.",
Return[solution["Derived"]["Tlmw"]];
];
];


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
TeukolskyZInfinity[ProcaSolution_?AssociationQ, TeukolskyLM\[Omega]_InterpolatingFunction,l_?IntegerQ,m_?IntegerQ]:=
Block[{
DeltaFunction,
RenormedAngMom,
SWSHEigenvalue,
TeukRin,
integrand,
integration,
Amplitude,
Zcoefficient,
\[CapitalDelta]=Evaluate[Kerr\[CapitalDelta]//FromxActVariables//ApplySolutionSet[ProcaSolution]], 
rdom = ProcaSolution["Solution", "R"]["Domain"]//First//HorizonCoordToRadial[#, ProcaSolution["Parameters", "\[Chi]"]]&,
\[Chi]v = SetPrecision[ProcaSolution["Parameters", "\[Chi]"],40],
\[Omega]v = SetPrecision[2*ProcaSolution["Solution", "\[Omega]"]//Re,40],
mv = 2*ProcaSolution["Parameters", "m"]
},
DeltaFunction = {r}|->Evaluate[\[CapitalDelta]];
RenormedAngMom = RenormedAngularMomentum[\[Omega]v,\[Chi]v,l,m,-2];
SWSHEigenvalue = SpinWeightedSpheroidalEigenvalue[-2,l,m,\[Chi]v*\[Omega]v];
$Messenger="Generating Homogeneous solution";
TeukRin = RNumerical[RenormedAngMom,\[Omega]v,\[Chi]v,l,m,-2,SWSHEigenvalue,rdom[[2]]];
(*TeukRin = Evaluate@TeukolskyRadial[-2,l,m,\[Chi]v,\[Omega]v , BoundaryConditions->"In", Method->{"NumericalIntegration", "Domain"->rdom}];
Amplitude = TeukolskyRadial[-2,l,m,\[Chi]v,\[Omega]v , Method->"MST"]["In"]["Amplitudes"]["Incidence"];*)
$Messenger="Calculating Asymptotic Amplitude \!\(\*SubscriptBox[\(B\), \(inc\)]\)";
If[(\[Not]NumericQ[RenormedAngMom])||(\[Not]NumericQ[SWSHEigenvalue]),
Print["Error: Renormalized angular momentum = "<>ToString[RenormedAngMom, InputForm]<>" \n SWSHEigenvalue = "<>ToString[SWSHEigenvalue,InputForm]];
];
Amplitude = IncomingAmplitudeB[RenormedAngMom,\[Omega]v,\[Chi]v,m,-2,SWSHEigenvalue];
integrand[r_?NumericQ]:= TeukolskyLM\[Omega][r]*TeukRin[r]/DeltaFunction[r]^2;
$Messenger="Performing Integration";
integration = NIntegrate[integrand[r], 
{r,rdom+1.5}//Flatten, 
PrecisionGoal->$TSNIntegratePrecision,
MaxRecursion->$TSNIntegrateMaxRecursion,
Method->$TSNIntegrateMethod,
WorkingPrecision->$TSPrecision];
Zcoefficient = integration/(2*I*\[Omega]v*Amplitude);



Print[integration];


$Messenger="";
Return[Zcoefficient]
];


Options[EnergyFlux]={TeukolskyTlmw->Automatic, ZCoefficient->Automatic};
EnergyFlux[ProcaSolution_?AssociationQ,l_?IntegerQ,m_?IntegerQ,OptionsPattern[]]:=
Block[{Zcoeff,
TeukModal,
EnFlux,
teuk\[Omega] = 2*Re[ProcaSolution["Solution", "\[Omega]"]]},

If[TrueQ[OptionValue[TeukolskyTlmw]==Automatic],
(*True*)
TeukModal = TeukolskySourceModal[ProcaSolution,l,m],
(*False*)
TeukModal = OptionValue[TeukolskyTlmw]
];(*end if*)

If[TrueQ[OptionValue[ZCoefficient]==Automatic],
Zcoeff = TeukolskyZInfinity[ProcaSolution,TeukModal, l,m];,
Zcoeff = OptionValue[ZCoefficient];
];

EnFlux = Abs[Zcoeff]^2/(4*\[Pi]*teuk\[Omega]^2);
Return[EnFlux]
];


NPPsi4::usage = "Compute the Newman-Penrose \!\(\*SubscriptBox[\(\[Psi]\), \(4\)]\) scalar at asymptotic null infinity";
Options[NPPsi4] = {ZCoefficient->Automatic, l->2, m->2, TeukModal->Automatic};
NPPsi4[ProcaSolution_?AssociationQ, OptionsPattern[]]:=
Block[{ZCoeff, SWSH, Psi4,teuk\[Omega] = 2*Re[ProcaSolution["Solution", "\[Omega]"]], teuk\[Chi] = ProcaSolution["Parameters", "\[Chi]"]},
If[TrueQ[OptionValue[ZCoefficient]==Automatic],
ZCoeff = TeukolskyZInfinity[ProcaSolution, OptionValue[TeukModal],OptionValue[l], OptionValue[m]];,
ZCoeff = OptionValue[ZCoefficient]
];
SWSH = {\[Theta]}|->SpinWeightedSpheroidalHarmonicS[-2, OptionValue[l], OptionValue[m], teuk\[Chi]*teuk\[Omega],\[Theta],0];
Psi4 = {t,r,\[Theta],\[Phi]}|->ZCoeff/Sqrt[2*\[Pi]]*SWSH[\[Theta]]*Exp[I*teuk\[Omega]*(TortoiseR[r,teuk\[Chi],1]-t) + I*OptionValue[m]*\[Phi]];
Return[Psi4];
];
