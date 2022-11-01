(* ::Package:: *)

EnergyMomentumwl;
(*Import Proca Mode Solver*)
If[\[Not]MemberQ[Names["Global`*"], "KerrWithProca"], Get[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]]
(*Import xAct Setup*)
If[\[Not]MemberQ[Names["Global`*"], "xActSetup"], Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]]
(*Import Helper functions*)
If[\[Not]MemberQ[Names["Global`*"], "HelperFunctions"], Get[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]]


(*$Assumptions={r>0,M>0,J>0,1>a\[GreaterEqual]0,\[Pi]\[GreaterEqual]\[Theta]\[GreaterEqual]0, \[Lambda]\[Element]Complexes,\[Omega]\[Element]Complexes, \[Nu]\[Element]Complexes};*)


M=1;
\[Delta]\[Theta]=10^-5;
SolutionPath =$FKKSRoot<>"Solutions/";

$EMMinRecursion = 0;
$EMMaxRecursion = Automatic;
$EMNIntegratePrecision = 5; (*Produces sufficient accuracy*)
$EMNIntegrateMethod = {"GlobalAdaptive"};

AnalyticSolution = <|"Parameters"-><|"\[Epsilon]"->\[Epsilon],"\[Mu]Nv"->\[Mu]Nv,"m"->m,"n"->n,"\[Chi]"->\[Chi],"M"->M|>,"Solution"-><|"\[Omega]"->\[Omega],"\[Nu]"->\[Nu],"R"->R,"S"->S|>|>;



complexToRealReprRule = {R[\[Epsilon]_]:>Rr[\[Epsilon]]+I*Ri[\[Epsilon]],Derivative[1][R][\[Epsilon]_]:>Derivative[1][Rr][\[Epsilon]]+I*Derivative[1][Ri][\[Epsilon]],S[t_]:>Sr[t]+I*Si[t],Derivative[1][S][t_]:>Derivative[1][Sr][t]+I*Derivative[1][Si][t], \[Omega]:>\[Omega]r+I*\[Omega]i, \[Nu]:>\[Nu]r+I*\[Nu]i};


(*A^\[Mu]*)
Options[FKKSProca]={SymbolicExpression->False, Optimized->False, RealPart->True, QuasiboundState->False, IndexStructure->up};
FKKSProca[solution_:analytic, OptionsPattern[]]:=
Block[{res,tmp,gradZ,Z,A,Aupreal, AuprealComponents,Ar,assumps, Filename,Filenamecomps},

If[OptionValue[RealPart],
(*Real filenames*)
If[TrueQ[OptionValue[IndexStructure]==up],
Filename = "FKKSProcaUpRealCTensor.mx";
Filenamecomps = "AuprealComponents.mx";
Filename="FKKSProcaDownRealCTensor.mx";
Filenamecomps = "AdownrealComponents.mx"
],

(*Not real filenames*)
If[TrueQ[OptionValue[IndexStructure]==up],
Filename = "FKKSProcaUpCTensor.mx";,
Filename = "FKKSProcaDownCTensor.mx"
];
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

If[OptionValue[IndexStructure]==up,
(*up*)
Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
A = (Polarization[\[Xi],\[Zeta]]Cd[-\[Zeta]]@Z)//ToBasisExpand//FromBasisExpand//Head;
CloseKernels[];LaunchKernels[4];
PrintTemporary@Dynamic[ProcaMessenger//MatrixForm];
ProcaMessenger=ConstantArray["", $KernelCount];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
KeyRule = ParallelTable[$KernelID->i, {i,0,3}];
AuprealComponents = {0,0,0,0};
SetSharedVariable[ProcaMessenger,  AuprealComponents];
DistributeDefinitions[KeyRule,complexToRealReprRule, A,ProcaComponentsFilePath];
ParallelDo[
ProcaMessenger[[1+($KernelID/.KeyRule)]]=Row[{"Kernel "<>ToString[$KernelID]<>" working on iteration "<>ToString[i], Spacer[20], ProgressIndicator[Appearance->"Percolate"]}];
assumps = {Derivative[1][Si][\[Theta]]\[Element]Reals, Derivative[1][Sr][\[Theta]]\[Element]Reals,Derivative[1][Ri][\[Theta]]\[Element]Reals,Derivative[1][Rr][\[Theta]]\[Element]Reals};
AuprealComponents[[i+1]] = (A[[1]][[i+1]]//.complexToRealReprRule//Re//ComplexExpand)//Simplify[#,Assumptions->assumps,TimeConstraint->Infinity]&;
ProcaMessenger[[1+($KernelID/.KeyRule)]]="Kernel "<>ToString[$KernelID]<>" done.";,
{i,0,3}];
Export[ProcaComponentsFilePath,AuprealComponents];
ProcaMessenger[[1+($KernelID/.KeyRule)]]="Kernel "<>ToString[$KernelID]<>" Done.";
];,
(*down*)
Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
A = (Polarization[-\[Xi],\[Zeta]]Cd[-\[Zeta]]@Z)//ToBasisExpand//FromBasisExpand//Head;
CloseKernels[];LaunchKernels[4];
PrintTemporary@Dynamic[ProcaMessenger//MatrixForm];
ProcaMessenger=ConstantArray["", $KernelCount];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
KeyRule = ParallelTable[$KernelID->i, {i,0,3}];
AuprealComponents = {0,0,0,0};
SetSharedVariable[ProcaMessenger,  AuprealComponents];
DistributeDefinitions[KeyRule,complexToRealReprRule, A,ProcaComponentsFilePath];
ParallelDo[
ProcaMessenger[[1+($KernelID/.KeyRule)]]=Row[{"Kernel "<>ToString[$KernelID]<>" working on iteration "<>ToString[i], Spacer[20], ProgressIndicator[Appearance->"Percolate"]}];
assumps = {Derivative[1][Si][\[Theta]]\[Element]Reals, Derivative[1][Sr][\[Theta]]\[Element]Reals,Derivative[1][Ri][\[Theta]]\[Element]Reals,Derivative[1][Rr][\[Theta]]\[Element]Reals};
AuprealComponents[[i+1]] = (A[[1]][[i+1]]//.complexToRealReprRule//Re//ComplexExpand)//Simplify[#,Assumptions->assumps,TimeConstraint->Infinity]&;
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
With[{ProcaFilePath = FileNameJoin[{ $FKKSRoot,"Expressions",Filename}]},
If[
FileExistsQ[ProcaFilePath],
res = Import[ProcaFilePath];,

If[TrueQ[OptionValue[IndexStructure]==up],
(*up*)
Block[{ch=sphericalchart},
Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
res = Head[Polarization[\[Zeta],\[Xi]]Cd[-\[Xi]]@Z//FromBasisExpand//Simplify];
Export[ProcaFilePath,res];
](*End of Block*),

(*down*)
Block[{ch=sphericalchart},
Z = R[r[]]*S[\[Theta][]]*Exp[-I*\[Omega]*t[]]*Exp[I*m*\[Phi][]];
res = Head[Polarization[-\[Zeta],\[Xi]]Cd[-\[Xi]]@Z//FromBasisExpand//Simplify];
Export[ProcaFilePath,res];
](*End of Block*)
];
](*End of outer If*)
] (*End of Not Real Part*)
]; (*End of If Optionvalue[RealPart]*)

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
Subscript[F, \[Mu]\[Nu]]=\!\(
\*SubscriptBox[\(\[Del]\), \(\[Mu]\)]
\*SubscriptBox[\(A\), \(\[Nu]\)]\)-\!\(
\*SubscriptBox[\(\[Del]\), \(\[Nu]\)]
\*SubscriptBox[\(A\), \(\[Mu]\)]\)
*)
Options[FKKSFieldStrength]={SymbolicExpression->False, Optimized->False, RealPart->True, QuasiboundState->False};
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
A=FKKSProca[analytic, RealPart->OptionValue[RealPart]];
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
Options[FKKSEnergyMomentum]={SymbolicExpression->False, Debug->False, Optimized->False, RealPart->True, QuasiboundState->False};
FKKSEnergyMomentum[solution_, OptionsPattern[]]:=
Block[{res,tmp,tmp$, mass = \[Mu]Nv, tmpOE,F,A, Filename},
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

T = FKKSEnergyMomentum[analytic, RealPart->OptionValue[RealPart]];
res = -1*T[{0,sphericalchart},{0,-sphericalchart}]//FromxActVariables;
Export[ENFilePath, res];
];
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
IntegrationErrors, 
energydensityInterpolation,
SpacialEnergyDensity,
nradialpoints, 
n\[Theta]points, 
n\[Phi]points,
\[Chi] = solution["Parameters", "\[Chi]"],
Radialdomain = HorizonCoordToRadial[solution["Solution", "R"]["Domain"][[1]], solution["Parameters", "\[Chi]"]]
},

energydensity = FKKSEnergyDensity[solution, ToCompiled->True, QuasiboundState->OptionValue[QuasiboundState]];

SpacialEnergyDensity = {r,\[Theta],\[Phi]}|->energydensity[0,r,\[Theta],\[Phi]];

weight = {r,\[Theta]}|->Evaluate[Kerrmetweight[r,\[Theta],\[Chi]]];

integrand[r_?NumericQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ]:=SpacialEnergyDensity[r,\[Theta],\[Phi]]*weight[r,\[Theta]];

(*Safety check on integrand*)
If[!NumericQ[integrand[2,2,2]],
	Print["Error! Integrand not a number! integrand[2,2,2] = "<>ToString[integrand[2,2,2], InputForm]];
	Return[Null];
];
res = NIntegrate[
integrand[r,\[Theta],\[Phi]], 
{r,Radialdomain[[1]],Radialdomain[[1]],Radialdomain[[2]]}, (*Exclude r=Subscript[r, initial]*)
{\[Theta],\[Theta]\[Epsilon],\[Theta]\[Epsilon], \[Pi]-\[Theta]\[Epsilon], \[Pi]-\[Theta]\[Epsilon]},(*Exclude \[Theta]=0,\[Pi]*)
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
totalenergyresults = FKKSTotalEnergy[solution, QuasiboundState->True];
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
RenormalizedEnergyErrors = totalenergyresults["Errors"]*normalization^2;

<|
	"Normalization"->normalization, 
	"NormalizationFractionalError"->NormalizationError/normalization,
	"TotalEnergy"->RenormalizedEnergy, 
	"TotalEnergyErrors"->RenormalizedEnergyErrors, 
	"FinalMass"->finalmass, 
	"FinalSpin"->FinalDimlessSpin[\[Chi]v,wv,mv,InitialMass, finalmass/InitialMass]
|>
]
