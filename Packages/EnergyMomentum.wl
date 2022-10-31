(* ::Package:: *)

EnergyMomentumwl;
(*Import Proca Mode Solver*)
If["KerrWithProca"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]]
(*Import xAct Setup*)
If["xActSetup"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]]
(*Import Helper functions*)
If["HelperFunctions"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]]

$EMMaxRadialPoints=200;
$EMCoefficient\[Theta]points=50;


(*$Assumptions={r>0,M>0,J>0,1>a\[GreaterEqual]0,\[Pi]\[GreaterEqual]\[Theta]\[GreaterEqual]0, \[Lambda]\[Element]Complexes,\[Omega]\[Element]Complexes, \[Nu]\[Element]Complexes};*)


M=1;
\[Delta]\[Theta]=10^-4;
SolutionPath =$FKKSRoot<>"Solutions/";

$EMMinRecursion = 0;
$EMMaxRecursion = Automatic;
$EMNIntegratePrecision = 5; (*Produces sufficient accuracy*)
$EMNIntegrateMethod = {"GlobalAdaptive"};

AnalyticSolution = <|"Parameters"-><|"\[Epsilon]"->\[Epsilon],"\[Mu]Nv"->\[Mu]Nv,"m"->m,"n"->n,"\[Chi]"->\[Chi],"M"->M|>,"Solution"-><|"\[Omega]"->\[Omega],"\[Nu]"->\[Nu],"R"->R,"S"->S|>|>;



complexToRealReprRule = {R[\[Epsilon]_]:>Rr[\[Epsilon]]+I*Ri[\[Epsilon]],Derivative[1][R][\[Epsilon]_]:>Derivative[1][Rr][\[Epsilon]]+I*Derivative[1][Ri][\[Epsilon]],S[t_]:>Sr[t]+I*Si[t],Derivative[1][S][t_]:>Derivative[1][Sr][t]+I*Derivative[1][Si][t], \[Omega]:>\[Omega]r+I*\[Omega]i, \[Nu]:>\[Nu]r+I*\[Nu]i};
killImagTerms = {
Im[Rr[\[Epsilon]_]]:>0,Im[Ri[\[Epsilon]_]]:>0,
 Re[Rr[\[Epsilon]_]]:>Rr[\[Epsilon]], Re[Ri[\[Epsilon]_]]:>Ri[\[Epsilon]],
Im[Derivative[1][Rr][\[Epsilon]_]]:>0,Im[Derivative[1][Ri][\[Epsilon]_]]:>0,
 Re[Derivative[1][Rr][\[Epsilon]_]]:>Derivative[1][Rr][\[Epsilon]], Re[Derivative[1][Ri][\[Epsilon]_]]:>Derivative[1][Ri][\[Epsilon]],
Im[(Rr^\[Prime]\[Prime])[\[Epsilon]_]]:>0,Im[(Ri^\[Prime]\[Prime])[\[Epsilon]_]]->0,
 Re[(Rr^\[Prime]\[Prime])[\[Epsilon]_]]:>(Rr^\[Prime]\[Prime])[\[Epsilon]], Re[(Ri^\[Prime]\[Prime])[\[Epsilon]_]]:>(Ri^\[Prime]\[Prime])[\[Epsilon]],
Im[Sr[\[Epsilon]_]]:>0,Im[Si[\[Epsilon]_]]:>0,
 Re[Sr[\[Epsilon]_]]:>Sr[\[Epsilon]], Re[Si[\[Epsilon]_]]:>Si[\[Epsilon]],
Im[Derivative[1][Sr][\[Epsilon]_]]:>0,Im[Derivative[1][Si][\[Epsilon]_]]:>0,
 Re[Derivative[1][Sr][\[Epsilon]_]]:>Derivative[1][Sr][\[Epsilon]], Re[Derivative[1][Si][\[Epsilon]_]]:>Derivative[1][Si][\[Epsilon]],
Im[(Sr^\[Prime]\[Prime])[\[Epsilon]_]]:>0,Im[(Si^\[Prime]\[Prime])[\[Epsilon]_]]->0,
 Re[(Sr^\[Prime]\[Prime])[\[Epsilon]_]]:>(Sr^\[Prime]\[Prime])[\[Epsilon]], Re[(Si^\[Prime]\[Prime])[\[Epsilon]_]]:>(Si^\[Prime]\[Prime])[\[Epsilon]]};


(*A^\[Mu]*)
Options[FKKSProca]={SymbolicExpression->False, Optimized->False, RealPart->True, QuasiboundState->False, IndexStructure->up, AsInterpolatingFunction->False};
FKKSProca[solution_:analytic, OptionsPattern[]]:=
Block[{res,tmp,gradZ,Z,A,Aupreal, AuprealComponents,Ar, Filename,Filenamecomps,radialdomain,mv,\[Omega]v},
If[OptionValue[RealPart],
Filename = "FKKSProcaRealCTensor.mx";
Filenamecomps = "AuprealComponents.mx";,

Filename = "FKKSProcaCTensor.mx";
];

If[OptionValue[RealPart],
(*Take real part*)
With[{
ProcaFilePath =FileNameJoin[{ $FKKSRoot,"Expressions",Filename}],ProcaComponentsFilePath = FileNameJoin[{ $FKKSRoot,"Expressions",Filenamecomps}]},

(*If procafile exists, import it, else build it*)
If[
FileExistsQ[ProcaFilePath],
res = Import[ProcaFilePath];,

(*If components were already calculated, import it, else calculate them (WARNING: will take probably 20 mins)*)
If[
FileExistsQ[ProcaComponentsFilePath],
AuprealComponents = Import[ProcaComponentsFilePath];,

Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
A = (Polarization[\[Xi],\[Zeta]]Cd[-\[Zeta]]@Z)//ToBasisExpand//FromBasisExpand//Head;
CloseKernels[];LaunchKernels[4];
PrintTemporary@Dynamic[ProcaMessenger//MatrixForm];
ProcaMessenger=ConstantArray["", $KernelCount];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
KeyRule = ParallelTable[$KernelID->i, {i,0,3}];
AuprealComponents = {0,0,0,0};
SetSharedVariable[ProcaMessenger,  AuprealComponents];
DistributeDefinitions[KeyRule,complexToRealReprRule, killImagTerms,A,ProcaComponentsFilePath];
ParallelDo[
ProcaMessenger[[1+($KernelID/.KeyRule)]]=Row[{"Kernel "<>ToString[$KernelID]<>" working on iteration "<>ToString[i], Spacer[20], ProgressIndicator[Appearance->"Percolate"]}];
AuprealComponents[[i+1]] = (A[[1]][[i+1]]//.complexToRealReprRule//Re//ComplexExpand)/.killImagTerms//Simplify[#,TimeConstraint->Infinity]&;
ProcaMessenger[[1+($KernelID/.KeyRule)]]="Kernel "<>ToString[$KernelID]<>" done.";,
{i,0,3}];
Export[ProcaComponentsFilePath,AuprealComponents];
ProcaMessenger[[1+($KernelID/.KeyRule)]]="Kernel "<>ToString[$KernelID]<>" Done.";
];
DefTensor[Ar[\[Zeta]],\[ScriptCapitalM]];
Block[{Print},ComponentValue[ComponentArray[Ar[{\[Zeta],sphericalchart}]], AuprealComponents]];
res = ToCTensor[Ar, {sphericalchart}]/.TensorValues[Ar];
Export[ProcaFilePath,res];
]
];,

(*Not real part*)
With[{ProcaFilePath = $FKKSRoot<>"Expressions/FKKSProcaCTensor.mx"},
If[
FileExistsQ[ProcaFilePath],
res = Import[ProcaFilePath];,
Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
res = Head[Polarization[\[Zeta],\[Xi]]Cd[-\[Xi]]@Z//FromBasisExpand//Simplify];
Export[ProcaFilePath,res];
]
];
];
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];

If[OptionValue[RealPart],
tmp = res/.ToParamSymbols//FromxActVariables//ApplyRealSolutionSet[solution, QuasiboundState->OptionValue[QuasiboundState]];,
tmp = res/.ToParamSymbols//FromxActVariables//PrimedToSymbolic//ApplySolutionSet[solution];
];

If[OptionValue[SymbolicExpression], 
Return[tmp]
];

If[OptionValue[Optimized],
Return[OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[tmp]]]
];

If[OptionValue[AsInterpolatingFunction],
radialdomain = solution["Solution", "R"]["Domain"]//First;
mv = solution["Parameters", "m"];
\[Omega]v= solution["Solution", "\[Omega]"];
Block[{ dmat, Idmat, tmpoe, VectorSamplingValues, coeff1, coeff2,coeff1oe,coeff2oe,coeff1Interp, coeff2Interp,coeff1interpfunction,coeff2interpfunction},
	With[{\[Phi]sampling = {0, \[Pi]/(2*mv)},
			rdom = radialdomain//HorizonCoordToRadial[#,solution["Parameters", "\[Chi]"]]&,
			rpoints = If[radialdomain[[-1]]<$EMMaxRadialPoints, radialdomain[[-1]], $EMMaxRadialPoints+10*Log[radialdomain[[-1]]]],
			\[Theta]points = $EMCoefficient\[Theta]points
			},
			dmat[\[Phi]1_, \[Phi]2_]:= {{Cos[\[Phi]1], Sin[\[Phi]1]}, {Cos[\[Phi]2], Sin[\[Phi]2]}};
			Idmat[\[Phi]1_, \[Phi]2_]:={{-Csc[\[Phi]1-\[Phi]2] Sin[\[Phi]2],Csc[\[Phi]1-\[Phi]2] Sin[\[Phi]1]},{Cos[\[Phi]2] Csc[\[Phi]1-\[Phi]2],-Cos[\[Phi]1] Csc[\[Phi]1-\[Phi]2]}};
			tmpoe = OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[tmp[{\[Zeta],sphericalchart}]//ComponentArray]];
			VectorSamplingValues = {tmpoe[0,r,\[Theta],\[Phi]sampling[[1]]], tmpoe[0,r,\[Theta],\[Phi]sampling[[2]]]};
			{coeff1, coeff2} = Idmat[\[Phi]sampling[[1]], \[Phi]sampling[[2]]] . VectorSamplingValues;
			coeff1oe = OptimizedFunction[{r,\[Theta]}, Evaluate[coeff1]];
			coeff2oe = OptimizedFunction[{r,\[Theta]}, Evaluate[coeff2]];
			coeff1Interp = GenerateInterpolation[coeff1oe[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points},DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];
			coeff2Interp = GenerateInterpolation[coeff2oe[r,\[Theta]], {r,rdom[[1]], rdom[[2]],(rdom[[2]]-rdom[[1]])/rpoints}, {\[Theta],\[Theta]\[Epsilon],\[Pi]-\[Theta]\[Epsilon],\[Pi]/\[Theta]points}, DensityFunctions->{(#^4&), Identity}, InterpolationOrder->3, Method->"Spline"];	
		];
	coeff1interpfunction = {r$,\[Theta]$}|->Evaluate@Table[coeff1Interp[[i]][r$,\[Theta]$], {i,1,4}];
	coeff2interpfunction = {r$,\[Theta]$}|->Evaluate@Table[coeff2Interp[[i]][r$,\[Theta]$], {i,1,4}];
	With[{mvv = mv, \[Omega]vv=\[Omega]v, C1 = coeff1interpfunction, C2 = coeff2interpfunction},
	If[OptionValue[QuasiboundState],
	decomposed = {t$,r$,\[Theta]$,\[Phi]$}|->Evaluate[C1[r$,\[Theta]$]*Cos[-Re[\[Omega]vv]*t$ + mvv*\[Phi]$]+C2[r$,\[Theta]$]*Sin[-Re[\[Omega]vv]*t$ + mvv*\[Phi]$]],
	decomposed = {t$,r$,\[Theta]$,\[Phi]$}|->Evaluate[C1[r$, \[Theta]$]*Cos[-\[Omega]vv*t$ + mvv*\[Phi]$]+C2[r$,\[Theta]$]*Sin[-\[Omega]vv*t$ + mvv*\[Phi]$]]
	];
	];
	Return[{t,r,\[Theta],\[Phi]}|->Evaluate[CTensor[decomposed[t,r,\[Theta],\[Phi]], {sphericalchart}]]];
	];
];


Return[{t,r,\[Theta],\[Phi]}|->Evaluate[tmp]]
]


(*
Subscript[F, \[Mu]\[Nu]]=\!\(
\*SubscriptBox[\(\[Del]\), \(\[Mu]\)]
\*SubscriptBox[\(A\), \(\[Nu]\)]\)-\!\(
\*SubscriptBox[\(\[Del]\), \(\[Nu]\)]
\*SubscriptBox[\(A\), \(\[Mu]\)]\)
*)
Options[FKKSFieldStrength]={SymbolicExpression->False, Optimized->False, RealPart->True, QuasiboundState->False, AsInterpolatingFunction->False};
FKKSFieldStrength[solution_, OptionsPattern[]]:=
Block[{res,tmp,A,Filename},
If[OptionValue[RealPart],
Filename = "FKKSFieldStrengthRealCTensor.mx";,
Filename = "FKKSFieldStrengthCTensor.mx";
];
With[{FSFilePath = $FKKSRoot<>"Expressions/"<>Filename},
If[
FileExistsQ[FSFilePath],
res = Import[FSFilePath];,
A=FKKSProca[analytic, RealPart->OptionValue[RealPart], AsInterpolatingFunction->OptionValue[AsInterpolatingFunction]];
res = Head[Cd[-\[Zeta]]@A[-\[Xi]] - Cd[-\[Xi]]@A[-\[Zeta]]];
Export[FSFilePath, res];
]
];
(*Format output based on input arguments*)
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];
If[OptionValue[RealPart],
tmp = res/.ToParamSymbols//FromxActVariables//ApplyRealSolutionSet[solution, QuasiboundState->OptionValue[QuasiboundState]];,
tmp = res/.ToParamSymbols//FromxActVariables//PrimedToSymbolic//ApplySolutionSet[solution];
];

If[OptionValue[SymbolicExpression], 
Return[tmp]
];

If[OptionValue[Optimized],
Return[OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[tmp]]]
];

Return[{t,r,\[Theta],\[Phi]}|->Evaluate[tmp]]
]


(*
Subscript[\[ScriptCapitalT], \[Mu]\[Nu]]
*)
Options[FKKSEnergyMomentum]={SymbolicExpression->False, Debug->False, Optimized->False, RealPart->True, QuasiboundState->True, AsInterpolatingFunction->False};
FKKSEnergyMomentum[solution_Association, OptionsPattern[]]:=
Block[{res,tmp,tmp$, mass = \[Mu]Nv, tmpOE,F,A, Filename, chartcoords={t[],r[],\[Theta][],\[Phi][]}},
If[OptionValue[AsInterpolatingFunction],
InterpedAup = FKKSProca[solution, AsInterpolatingFunction->True, RealPart->OptionValue[RealPart], QuasiboundState->OptionValue[QuasiboundState]][Sequence@@chartcoords];
InterpedAdown = (InterpedAup[{-\[Xi],-sphericalchart}]//ToBasis[sphericalchart]//ComponentArray//ApplyRealSolutionSet[MySolution, QuasiboundState->OptionValue[QuasiboundState]]);
With[{EMFilePath=$FKKSRoot<>"Expressions/FKKSEnergyMomentumInVector.mx"},
If[FileExistsQ[EMFilePath],
EMinA = Import[EMFilePath],
A=CTensor[{A0[t[],r[],\[Theta][],\[Phi][]],A1[t[],r[],\[Theta][],\[Phi][]],A2[t[],r[],\[Theta][],\[Phi][]],A3[t[],r[],\[Theta][],\[Phi][]]},{sphericalchart}];
F = Head[Cd[-\[Zeta]]@A[-\[Xi]] - Cd[-\[Xi]]@A[-\[Zeta]]];
EMinA = F[-\[Zeta],-\[Gamma]]F[-\[Xi],\[Gamma]] + mass^2*A[-\[Zeta]]A[-\[Xi]] -1/4 met[-\[Zeta],-\[Xi]](F[-\[Iota],-\[Gamma]]F[\[Iota],\[Gamma]] + 2*mass^2*A[-\[Gamma]]A[\[Gamma]])//ToBasis[sphericalchart]//ComponentArray//Simplify;
Export[EMFilePath,EMinA];
];
];
Areprule = Table[ToExpression["A"<>ToString[i]][Sequence@@chartcoords]->InterpedAdown[[i+1]], {i,0,3}];
PDreprule = Table[PDsphericalchart[{i,-sphericalchart}][A[{j,-sphericalchart}]]->D[InterpedAdown[[j+1]], chartcoords[[i+1]]], {i,0,3},{j,0,3}];
tmp$ = EMinA/.PDreprule/.Areprule//FromxActVariables//ApplyRealSolutionSet[solution, QuasiboundState->OptionValue[QuasiboundState]];
res = OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[CTensor[tmp$, {-sphericalchart, -sphericalchart}]]];
Return[res];
,
If[OptionValue[RealPart],
Filename = "FKKSEnergyMomentumRealCTensor.mx";,
Filename = "FKKSEnergyMomentumCTensor.mx";
];
With[{EMFilePath = $FKKSRoot<>"Expressions/"<>Filename},
If[
FileExistsQ[EMFilePath],
res = Import[EMFilePath];,

F = FKKSFieldStrength[analytic, RealPart->OptionValue[RealPart]];
A=FKKSProca[analytic, RealPart->OptionValue[RealPart]];
tmp$ = F[-\[Zeta],-\[Gamma]]F[-\[Xi],\[Gamma]] + mass^2*A[-\[Zeta]]A[-\[Xi]] -1/4 met[-\[Zeta],-\[Xi]](F[-\[Iota],-\[Gamma]]F[\[Iota],\[Gamma]] + 2*mass^2*A[-\[Gamma]]A[\[Gamma]]);
res=Evaluate@Head[tmp$];
If[TrueQ[Head[res]!=CTensor],Print["ERROR: Form of EM tensor incorrect"],Export[EMFilePath,res]];
]
];
(*Format output based on input arguments*)
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];
If[OptionValue[RealPart],
tmp = res/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution, QuasiboundState->OptionValue[QuasiboundState]];,
tmp = res/.ToParamSymbols//FromxActVariables//PrimedToSymbolic//ApplySolutionSet[solution];
];
If[OptionValue[SymbolicExpression], 
Return[tmp]
];

If[OptionValue[Optimized],
Return[OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[tmp]]]
];

Return[{t,r,\[Theta],\[Phi]}|->Evaluate[tmp]]

]
]


(*
\[Rho] == Subscript[\[ScriptCapitalT]^t, t]
*)
Options[FKKSEnergyDensity]={SymbolicExpression->False,Optimized->True,ToCompiled->False, RealPart->True , WithProperties->True, QuasiboundState->True}\[Union]Options[Experimental`OptimizeExpression]\[Union]Options[Compile];
FKKSEnergyDensity[solution_Association, OptionsPattern[]]:=
Block[{res, resOE,T,OptimizedResult, tmp,UnOptimizedResult, Filename},
If[OptionValue[RealPart],
Filename = "FKKSEnergyDensityReal.mx";,
Filename = "FKKSEnergyDensity.mx";
];
With[{ENFilePath =$FKKSRoot<>"Expressions/"<>Filename},
If[FileExistsQ[ENFilePath],
res = Import[ENFilePath];,

T = FKKSEnergyMomentum[analytic, RealPart->OptionValue[RealPart], QuasiboundState->OptionValue[QuasiboundState]];
res = -1*T[{0,sphericalchart},{0,-sphericalchart}]//FromxActVariables;
Export[ENFilePath, res];
];
];
(*Format output based on input arguments*)
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];
If[OptionValue[RealPart],
tmp = res/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution, QuasiboundState->OptionValue[QuasiboundState]];,
tmp = res/.ToParamSymbols//FromxActVariables//PrimedToSymbolic//ApplySolutionSet[solution];
];
If[OptionValue[SymbolicExpression], 
Return[tmp]
];

If[OptionValue[Optimized],
Return[
Evaluate@OptimizedFunction[{t,r,\[Theta],\[Phi]}, Evaluate[tmp], 
ToCompiled->OptionValue[ToCompiled],
WithProperties->OptionValue[WithProperties],
ExcludedForms->OptionValue[ExcludedForms],
ExternalForms->OptionValue[ExternalForms],
InertForms->OptionValue[InertForms],
OptimizationLevel->OptionValue[OptimizationLevel],
OptimizationSymbol->OptionValue[OptimizationSymbol],
ScopingSymbol->OptionValue[ScopingSymbol],
CompilationOptions->OptionValue[CompilationOptions],
CompilationTarget->OptionValue[CompilationTarget],
Parallelization->OptionValue[Parallelization],
RuntimeAttributes->OptionValue[RuntimeAttributes],
RuntimeOptions->OptionValue[RuntimeOptions]
]
]
];

Return[{t,r,\[Theta],\[Phi]}|->Evaluate[tmp]]
]


(*
E = \[Integral]\[Rho]\[Sqrt](-g)dV
*)
Options[FKKSTotalEnergy]={QuasiboundState->True};
FKKSTotalEnergy[solution_Association, OptionsPattern[]]:=
Block[{TotalEnergyMessenger,
energydensity,
res, 
weight,
integrand,
result,
rmin,
IntegrationErrors, 
energydensityInterpolation,
SpacialEnergyDensity,
nradialpoints, 
n\[Theta]points, 
n\[Phi]points,
\[Chi] = solution["Parameters", "\[Chi]"],
Radialdomain = HorizonCoordToRadial[solution["Solution", "R"]["Domain"][[1]], solution["Parameters", "\[Chi]"]]
},

energydensity = FKKSEnergyDensity[solution, ToCompiled->True,QuasiboundState->OptionValue[QuasiboundState]];

SpacialEnergyDensity = {r,\[Theta],\[Phi]}|->energydensity[0,r,\[Theta],\[Phi]];

weight = {r,\[Theta]}|->Evaluate[Kerrmetweight[r,\[Theta],\[Chi]]];

integrand[r_?NumericQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ]:=SpacialEnergyDensity[r,\[Theta],\[Phi]]*weight[r,\[Theta]];

(*Safety check on integrand*)
If[!NumericQ[integrand[2,2,2]],
	Print["Error! Integrand not a number! integrand[2,2,2] = "<>ToString[integrand[2,2,2], InputForm]];
	Return[Null];
];
rmin = rplusN[solution["Parameters", "\[Chi]"]]+10^-1.9;
res = NIntegrate[
integrand[r,\[Theta],\[Phi]], 
{r,rmin,rmin,Radialdomain[[2]]}, (*Exclude r=Subscript[r, initial]*)
{\[Theta],\[Delta]\[Theta],\[Delta]\[Theta], \[Pi]-\[Delta]\[Theta], \[Pi]-\[Delta]\[Theta]},(*Exclude \[Theta]=0,\[Pi]*)
{\[Phi],0,2\[Pi]},
MaxRecursion->$EMMaxRecursion,
Method->$EMNIntegrateMethod,
MinRecursion->$EMMinRecursion,
PrecisionGoal->$EMNIntegratePrecision,
IntegrationMonitor:>((IntegrationErrors=Through[#1@"Error"])&)
];
If[Total@Errors>0.01,
Print["ERROR: integration errors greater than 1%"];
];
TotalEnergyMessenger = Row[{"Integration complete. Returning result ", ProgressIndicator[Appearance->"Percolate"]}];

result = <|"Total"->res, "Errors"->Total@IntegrationErrors|>;
Return[result]
]


FinalDimlessSpin::usage="
Final dimensionless spin of the black hole in terms of the initial mass (Relative to an energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\)!), the real part of the proca frequency, the mode number, the initial dimensionless frequency, and the ratio between the final mass and initial mass.

FinalDimlessSpin[InitialDimlessSpin_,freq_,mode_,MInitial_,MRatio_]

InitialDimlessSpin: Initial dimensionless spin of kerr black hole.
freq: Real part of the proca frequency (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))
mode: mode number of proca cloud.
MInitial: Initial mass of kerr black hole (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\)).
MRatio: Ratio between final mass \!\(\*SubscriptBox[\(M\), \(f\)]\) and inital mass \!\(\*SubscriptBox[\(M\), \(i\)]\) of kerr black hole, \!\(\*FractionBox[SubscriptBox[\(M\), \(f\)], SubscriptBox[\(M\), \(i\)]]\).
";
FinalDimlessSpin[InitialDimlessSpin_,freq_,mode_,MInitial_,MRatio_]:=mode/freq*1/MInitial*1/MRatio (1-1/MRatio)+InitialDimlessSpin/MRatio^2;


SaturationCondition::usage="
Condition that must be satisfied for superradiant amplification to cease. Condition is satisfied when this function equals 0.
\!\(\*FractionBox[\(\[Omega]\), \(m\)]\)-\!\(\*FractionBox[SubscriptBox[\(\[Chi]\), \(f\)], \(2 \*SubscriptBox[\(M\), \(f\)] \((1 - \[Sqrt]\((1 - \*SuperscriptBox[SubscriptBox[\(\[Chi]\), \(f\)], \(2\)])\))\)\)]\)=0

SaturationCondition[freq_,mode_,InitialDimlessSpin_,MInitial_,MRatio_]

freq: real part of proca frequency (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))
mode: mode number of proca cloud
InitialDimlessSpin: Initial dimensionless spin of kerr black hole
MInitial: Initial mass of kerr black hole (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))
MRatio: Ratio between final and initial mass of kerr black hole, \!\(\*FractionBox[SubscriptBox[\(M\), \(f\)], SubscriptBox[\(M\), \(i\)]]\)
";
SaturationCondition[freq_,mode_,InitialDimlessSpin_,MInitial_,MRatio_]:=
Block[{chiFinal = FinalDimlessSpin[InitialDimlessSpin,freq,mode,MInitial,MRatio]},
freq/mode-chiFinal/(2*MRatio*MInitial*(1+Sqrt[1-chiFinal^2]))
]


FKKSFinalMass::usage="
compute final mass of Kerr black hole.

FinalMass[freq_,mode_,InitialDimlessSpin_,MInitial_]

freq: real part of proca frequency (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))
mode: mode number of proca cloud
InitialDimlessSpin: Initial dimensionless spin of kerr black hole
MInitial: Initial mass of kerr black hole (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))

Returns final mass of kerr black hole (relative to energy scale \!\(\*SubscriptBox[\(M\), \(0\)]\))
";
FKKSFinalMass[freq_,mode_,InitialDimlessSpin_,MInitial_]:=
Block[{mm, MRatio, allsolutions},
allsolutions = mm/.NSolve[SaturationCondition[freq,mode,InitialDimlessSpin,MInitial,mm],mm,Reals];
MRatio = Select[allsolutions,  0<=#<=1&]//First; (*Assuming energy flows from BH to proca cloud*)
MRatio*MInitial
]


Options[FKKSNormalization]={Recalculate->False};
FKKSNormalization[solution_?AssociationQ, InitialMass_, OptionsPattern[]]:=
Block[{totalenergy, totalenergyresults,finalmass, normalization,wv = solution["Solution", "\[Omega]"]//Re, mv = solution["Parameters", "m"], \[Chi]v = solution["Parameters", "\[Chi]"]},
If[!OptionValue[Recalculate],
If[KeyExistsQ[solution, "Derived"],
If[KeyExistsQ[solution["Derived"], "Normalization"],
Return[solution["Derived", "Normalization"]];
]
];
];
totalenergyresults = FKKSTotalEnergy[solution,QuasiboundState->True];
totalenergy = totalenergyresults["Total"];
finalmass = FKKSFinalMass[wv, mv, \[Chi]v, InitialMass];
normalization = Sqrt[(InitialMass-finalmass)/totalenergy];

(*Error Analysis*)
NormalizationError = normalization*totalenergyresults["Errors"]/totalenergy//Total;
If[NormalizationError/normalization>10^-2,
	Print["Error! normalization uncertainty greater than 0.01 %"];
];

(*Renormalize energy*)
RenormalizedEnergy = totalenergy*normalization^2;

<|
	"Normalization"->normalization, 
	"NormalizationFractionalError"->NormalizationError/normalization,
	"TotalEnergy"->RenormalizedEnergy, 
	"TotalEnergyErrors"->totalenergyresults["Errors"], 
	"FinalMass"->finalmass, 
	"FinalSpin"->FinalDimlessSpin[\[Chi]v,wv,mv,InitialMass, finalmass/InitialMass]
|>
]

