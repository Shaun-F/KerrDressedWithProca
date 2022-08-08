(* ::Package:: *)

EnergyMomentumwl
(*Import Proca Mode Solver*)
<<KerrWithProca`
(*Import xAct Setup*)
If["xActSetup"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]]

(*$Assumptions={r>0,M>0,J>0,1>a\[GreaterEqual]0,\[Pi]\[GreaterEqual]\[Theta]\[GreaterEqual]0, \[Lambda]\[Element]Complexes,\[Omega]\[Element]Complexes, \[Nu]\[Element]Complexes};*)


M=1;
SolutionPath =$FKKSRoot<>"Solutions/";
(*
AllSolutions = getResults[SolutionPath];
\[Mu]OrderingIndex = Ordering@AllSolutions[[All, "Parameters", "\[Mu]Nv"]];
AllSolutions = AllSolutions[[\[Mu]OrderingIndex]];
*)
AnalyticSolution = <|"Parameters"-><|"\[Epsilon]"->\[Epsilon],"\[Mu]Nv"->\[Mu]Nv,"m"->m,"n"->n,"\[Chi]"->\[Chi],"M"->M|>,"Solution"-><|"\[Omega]"->\[Omega],"\[Nu]"->\[Nu],"R"->R,"S"->S|>|>;



OvertoneModeGroup[nvalue_?NumberQ, mvalue_?NumberQ]:=Select[Select[AllSolutions, #["Parameters"]["n"]==nvalue&], #["Parameters"]["m"]==mvalue&]


Kerrmetweight[r_,\[Theta]_,\[Chi]_]:=(r^2+\[Chi]^2*Cos[\[Theta]]^2)*Sin[\[Theta]]

OvertoneModeGroup[nvalue_?NumberQ, mvalue_?NumberQ]:=Select[Select[AllSolutions, #["Parameters"]["n"]==nvalue&], #["Parameters"]["m"]==mvalue&]

HorizonCoordToRadial[xN_, \[Chi]_]:=xN * (rplusN[\[Chi]] - rminusN[\[Chi]])  +  rplusN[\[Chi]];
RadialToHorizonCoord[rN_,\[Chi]_]:=(rN - rplusN[\[Chi]])/(rplusN[\[Chi]] - rminusN[\[Chi]]);

(*https://mathematica.stackexchange.com/questions/59944/extracting-the-function-from-interpolatingfunction-object/59963#59963*)
InterpolationToPiecewise[if_,x_]:=Module[{main,default,grid},grid=if["Grid"];
Piecewise[{if@"GetPolynomial"[#,x-#],x<First@#}&/@grid[[2;;-2]],if@"GetPolynomial"[#,x-#]&@grid[[-1]]]]/;if["InterpolationMethod"]=="Hermite"

(*Redundant?*)
CompileInterpolating[interp_InterpolatingFunction]:=Compile[{\[Psi]}, Evaluate[InterpolationToPiecewise[interp, \[Psi]]]]

ParamsToReprRule[solution_]:=AssociationThread[(Symbol/@(solution[ "Parameters"]//Keys))->(solution[ "Parameters"]//Values)]//Normal

ToRational[expr_]:=Rationalize[expr,0]
PrimedToSymbolic[expr_]:=expr/.{R^\[Prime]\[Prime]->ddR, Derivative[1][R]->dR};

SolToReprRule[solution_,OptionsPattern[{real->False}]]:=
With[{\[Omega]v = solution["Solution","\[Omega]"],\[Nu]v = solution["Solution","\[Nu]"], Rv=solution["Solution","R"], Sv = solution["Solution","S"], \[Chi]v = solution["Parameters","\[Chi]"]},
If[OptionValue[real],
AssociationThread[{\[Omega],\[Nu],S,R,dR,ddR}->Flatten@{Re[\[Omega]v],
											\[Nu]v,
											Sv,
											(Rv[RadialToHorizonCoord[#,\[Chi]v]]&), 
											((Derivative[1][Rv][RadialToHorizonCoord[#,\[Chi]v]]*1/(rplusN[\[Chi]v]-rminusN[\[Chi]v]))&),
											((Derivative[1][Derivative[1][Rv]][RadialToHorizonCoord[#,\[Chi]v]]*1/(rplusN[\[Chi]v]-rminusN[\[Chi]v])^2)&)
}
],

AssociationThread[{\[Omega],\[Nu],S,R,dR,ddR}->Flatten@{\[Omega]v,
											\[Nu]v,
											Sv,
											(Rv[RadialToHorizonCoord[#,\[Chi]v]]&), 
											((Derivative[1][Rv][RadialToHorizonCoord[#,\[Chi]v]]*1/(rplusN[\[Chi]v]-rminusN[\[Chi]v]))&),
											((Derivative[1][Derivative[1][Rv]][RadialToHorizonCoord[#,\[Chi]v]]*1/(rplusN[\[Chi]v]-rminusN[\[Chi]v])^2)&)
}
](*Dont forgot to convert radial function from horizon co-ordinates to Boyer-Lindquist radial coordinates! Polarization tensor is in terms of radial coordinates, so we must match the arguments for the Z fkks function.*)
]
];

ThreadThrough::argn="The arguments `1` and `2` must have equal lengths";
ThreadThrough::usage="Thread a list of functions over a list of arguments.
ThreadThrough[{f,g,h}, {x,y,z}] -> {f[x], g[y], h[z]}
";
ThreadThrough[functionlist_?ListQ, Operands_?ListQ]:=
(If[Length@functionlist!=Length@Operands,
Return@Message[MessageName[ThreadThrough,"argn"],functionlist, Operands]
];
Table[functionlist[[i]][Operands[[i]]], {i,1,Length@functionlist}]
)


ToParamSymbols:={a->\[Chi]};

ApplySolutionSet[solution_,OptionsPattern[{real->False}]][expr_]:= expr/.ToParamSymbols/.ParamsToReprRule[solution]/.SolToReprRule[solution,real->OptionValue[real]];

OptimizedFunction[vars_, expr_, OptionsPattern[{WithProperties->True,ToCompiled->False,OptimizationSymbol->aa, Experimental`OptimizeExpression, Compile}]]:=Module[{i,exproe, res,ret},
exproe = Experimental`OptimizeExpression[expr,
ExcludedForms->OptionValue[ExcludedForms],
ExternalForms->OptionValue[ExternalForms],
InertForms->OptionValue[InertForms],
OptimizationLevel->OptionValue[OptimizationLevel],
OptimizationSymbol->OptionValue[OptimizationSymbol],
ScopingSymbol->OptionValue[ScopingSymbol]
];
res = exproe/.Experimental`OptimizedExpression[f_]:>Function@@Hold[vars, f];

If[OptionValue[ToCompiled],
With[{resc = res, 
thevars = vars,
opts = {
CompilationOptions->OptionValue[CompilationOptions],
CompilationTarget->OptionValue[CompilationTarget],
Parallelization->OptionValue[Parallelization],
RuntimeAttributes->OptionValue[RuntimeAttributes],
RuntimeOptions->OptionValue[RuntimeOptions]
}
},
rett = Compile@@Hold[thevars, Evaluate[resc]@@thevars,opts];
ret =rett
];
Return[ret]
];

If[OptionValue[WithProperties],
(*some DownValue/OwnValue/ModuleNumber fuckery lies below*)

(****Assign the OwnValue of res into the DownValue of the unique GLOBAL variable defined by Unique[OptimizedFunction]. Create a further downvalue of this unique variable*)
(****ret[Properties] downvalue tells us the variable assigned to the output of OptimizedFunction (this function) is an optimized expression. Used later...*)
ret = Unique[OptimizedFunction];
With[{rett = ret},
DownValues[Evaluate@rett]=Flatten@{OwnValues[res]/.HoldPattern[res]->rett[Sequence@@Evaluate@Table[Pattern[Evaluate@vars[[i]],_]?NumericQ, {i,1,Length[vars]}]]/.HoldPattern[Function[\[Zeta]_,\[Xi]_]]:>\[Xi], Evaluate@HoldPattern[rett[Properties]]:>Optimized};
Return[Evaluate[rett]]
]
];

Return[res]
]

DetMet[r_,\[Theta]_,\[Chi]_]:= -(r^2+\[Chi]^2*Cos[\[Theta]]^2)^2*Sin[\[Theta]]^2

ErgoRadius[\[Theta]_,\[Chi]_, OptionsPattern[{coords->BL}]]:=
If[TrueQ[OptionValue[coords]==BL],1+Sqrt[1-\[Chi]^2*Cos[\[Theta]]^2],
If[TrueQ[OptionValue[coords]==Horizon],RadialToHorizonCoord[1+Sqrt[1-\[Chi]^2*Cos[\[Theta]]^2],\[Chi]],
Print["Coords must be either BL (boyer-lindquist) or Horizon"]
]
];

SetAttributes[GenerateInterpolation, HoldFirst];
GenerateInterpolation[funcform_, args___, OptionsPattern[{ DensityFunctions->Automatic, Metadata->False}]]:=
Block[{
func,
variables = Hold/@List[args][[All,1]],
StartPoints = List[args][[All,2]],
EndPoints = List[args][[All,3]],
NArgs = Length@List@args,
funcOnPoints,
densityfuncs,
RecastDensityToBoundary,
recastfunctions,
CoordinateRanges,
retval
},
func = Hold[funcform][[1,0]];(*Extract head*)
RecastDensityToBoundary[Identity,_,_][x_]:=x;
RecastDensityToBoundary[denfunc_,rstart_,rstop_][x_]:=denfunc[x]/((denfunc[rstop]-denfunc[rstart])/(rstop-rstart))+(rstart*denfunc[rstop]-rstop*denfunc[rstart])/(denfunc[rstop] -denfunc[rstart]);
SetAttributes[RecastDensityToBoundary, Listable];
SetAttributes[func, Listable];

If[TrueQ[OptionValue[DensityFunctions]==Automatic],
densityfuncs = ConstantArray[Identity,Length[variables]];,
densityfuncs = OptionValue[DensityFunctions]
];

SetAttributes[RecastDensityToBoundary, Listable];
recastfunctions = RecastDensityToBoundary[densityfuncs, StartPoints, EndPoints];
Messenger="Generating coordinate mesh";
CoordinateRanges =Table[
With[{arg = List[args][[coordIter]], rcfunc =recastfunctions[[coordIter]],var = variables[[coordIter]] },
Table[rcfunc[ReleaseHold[var]], arg]
],
{coordIter, 1, Length[variables]}];
Messenger="Mapping function over mesh. \n\t Function Sample with timing: "<>ToString[Part[func@@StartPoints//Timing,1], InputForm];
funcOnPoints = Outer[func,Sequence@@CoordinateRanges];

Messenger="Generating Interpolation function";
retval = ListInterpolation[funcOnPoints, CoordinateRanges];
If[OptionValue[Metadata],
Return[<|"Interpolation"->retval, "Mesh"->CoordinateRanges, "DensityFunctions"->densityfuncs|>],
Return[retval]
]
]

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

ApplyRealSolutionSet[solution_][expr_]:=
Block[{repr, solR=solution["Solution","R"], solS=solution["Solution", "S"]},
repr =Flatten[{
Thread[{\[Omega]r,\[Omega]i}->ReIm[solution["Solution","\[Omega]"]]],
Thread[{\[Nu]r,\[Nu]i}->ReIm[solution["Solution","\[Nu]"]]],
HoldPattern@Derivative[d_][Rr][x_]->Re[Derivative[d][solR][x]],
HoldPattern@Derivative[d_][Ri][x_]->Im[Derivative[d][solR][x]],
HoldPattern@Derivative[d_][Sr][x_]->Re[Derivative[d][solS][x]],
HoldPattern@Derivative[d_][Si][x_]->Im[Derivative[d][solS][x]],
Rr->( Re[solR[#]]&),
Ri->( Im[solR[#]]&),
Sr->( Re[solS[#]]&),
Si->( Im[solS[#]]&),
\[Chi]->solution["Parameters", "\[Chi]"],
m-> solution["Parameters", "m"],
\[Mu]Nv->solution["Parameters", "\[Mu]Nv"]
}];
expr//.repr
]


(*A^\[Mu]*)
FKKSProca[solution_:analytic, OptionsPattern[{SymbolicExpression->False, Optimized->False, RealPart->True}]]:=
Block[{res,tmp,gradZ,Z,A,Aupreal, AuprealComponents,Ar, Filename,Filenamecomps},
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
tmp = res/.ToParamSymbols//FromxActVariables//ApplyRealSolutionSet[solution];,
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


(*A^\[Mu]Subscript[A, \[Mu]]*)
FKKSProcaNorm[solution_:analytic,OptionsPattern[{SymbolicExpression->False, Optimized->False, RealPart->True}]]:=
Block[{tmp,res,A, Filename},
If[OptionValue[RealPart],
Filename = "FKKSProcaNormReal.mx";,
Filename = "FKKSProcaNorm.mx";
];
With[{FSFilePath = $FKKSRoot<>"Expressions/"<>Filename},
If[
FileExistsQ[FSFilePath],
res = Import[FSFilePath];,

A=FKKSProca[analytic, RealPart->OptionValue[RealPart]];
res = A[-\[Zeta]]A[\[Zeta]];
Export[FSFilePath, res];
]
];
(*Format output based on input arguments*)
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];
If[OptionValue[RealPart],
tmp = res/.ToParamSymbols//FromxActVariables//ApplyRealSolutionSet[solution];,
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
FKKSFieldStrength[solution_, OptionsPattern[{SymbolicExpression->False, Optimized->False, RealPart->True}]]:=
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
tmp = res/.ToParamSymbols//FromxActVariables//ApplyRealSolutionSet[solution];,
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
FKKSEnergyMomentum[solution_, OptionsPattern[{SymbolicExpression->False, Debug->False, Optimized->False, RealPart->True}]]:=
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
tmp = res/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution];,
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
FKKSEnergyDensity[solution_, OptionsPattern[{SymbolicExpression->False,Optimized->True,ToCompiled->False, RealPart->True , WithProperties->True,Experimental`OptimizeExpression, Compile}]]:=
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
tmp = res/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution];,
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
FKKSTotalEnergy[solution_]:=
Block[{TotalEnergyMessenger,energydensity,res, weight,integrationdomain, integrand,result,IntegrationErrors, energydensityInterpolation,tointerp,
\[Chi] = solution["Parameters", "\[Chi]"],
Radialdomain = HorizonCoordToRadial[solution["Solution", "R"]["Domain"][[1]], solution["Parameters", "\[Chi]"]]
},
PrintTemporary@Dynamic[TotalEnergyMessenger];

TotalEnergyMessenger = Row[{"Computing energy density", ProgressIndicator[Appearance->"Percolate"]}];
energydensity = FKKSEnergyDensity[solution];

TotalEnergyMessenger = Row[{"Generating interpolating function for energy density", ProgressIndicator[Appearance->"Percolate"]}];
tointerp = {r,\[Theta],\[Phi]}|->energydensity[0,r,\[Theta],\[Phi]];
energydensityInterpolation=Generate3DInterpolation[tointerp,Radialdomain[[1]], Radialdomain[[2]]/2, \[Theta]\[Epsilon], \[Pi]-\[Theta]\[Epsilon],0,2*\[Pi], Radialdomain[[2]], \[Theta]points, \[Phi]points];

weight = {r,\[Theta]}|->Kerrmetweight[r,\[Theta],\[Chi]];

TotalEnergyMessenger = Row[{"Integrating energy density", ProgressIndicator[Appearance->"Percolate"]}];
integrationdomain = Radialdomain;
integrand[r_?NumericQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ]:=energydensityInterpolation[r,\[Theta],\[Phi]]*weight[r,\[Theta]];
res = NIntegrate[
integrand[r,\[Theta],\[Phi]], 
{r,integrationdomain[[1]],integrationdomain[[1]],integrationdomain[[2]]}, (*Exclude r=Subscript[r, initial]*)
{\[Theta],\[Theta]\[Epsilon],\[Theta]\[Epsilon], \[Pi]-\[Theta]\[Epsilon], \[Pi]-\[Theta]\[Epsilon]},(*Exclude \[Theta]=0,\[Pi]*)
{\[Phi],0,2\[Pi]},
MaxRecursion->$MaxRecursion,
Method->$NIntegrateMethod,
MinRecursion->$MinRecursion,
PrecisionGoal->$NIntegratePrecision,
IntegrationMonitor:>((IntegrationErrors=Through[#1@"Error"])&)
];
If[Total@Errors>0.01,
Print["ERROR: integration errors greater than 1%"];
];
result = <|"Total"->res, "Errors"->Total@IntegrationErrors|>;
Return[result]
]


(*
Subscript[\[ScriptCapitalT]^\[Sigma], \[Sigma]]
*)
FKKSEnergyMomentumTrace[solution_, OptionsPattern[{SymbolicExpression->False,Optimized->True, RealPart->True}]]:=
Block[{res, resOE,fkksEM,OptimizedResult, tmp, Filename},
If[OptionValue[RealPart],
Filename = "FKKSEnergyMomentumTraceReal.mx";,
Filename = "FKKSEnergyMomentumTrace.mx";
];
With[{ENFilePath =$FKKSRoot<>"Expressions/"<>Filename},
If[FileExistsQ[ENFilePath],
res = Import[ENFilePath];,
fkksEM = FKKSEnergyMomentum[analytic, RealPart->OptionValue[RealPart]];
res = fkksEM[{\[Zeta],sphericalchart},{\[Xi],sphericalchart}]met[{-\[Zeta],-sphericalchart},{-\[Xi],-sphericalchart}]//TraceBasisDummy;
Export[ENFilePath, res];
];
];
(*Format output based on input arguments*)
If[TrueQ[solution==analytic],
Return[res/.ToParamSymbols]
];
If[OptionValue[RealPart],
tmp = res/.ToParamSymbols/.\[Omega]i->0//FromxActVariables//ApplyRealSolutionSet[solution];,
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
Block[{mm, MRatio},
MRatio=First[mm/.NSolve[SaturationCondition[freq,mode,InitialDimlessSpin,MInitial,mm],mm,Reals] ];
MRatio*MInitial
]


FKKSNormalization[solution_, InitialMass_]:=
Block[{totalenergy, finalmass, normalization,wv = solution["Solution", "\[Omega]"]//Re, mv = solution["Parameters", "m"], \[Chi]v = solution["Parameters", "\[Chi]"]},
totalenergy = FKKSTotalEnergy[solution]["Total"];
finalmass = FKKSFinalMass[wv, mv, \[Chi]v, InitialMass];
normalization = Sqrt[(InitialMass-finalmass)/totalenergy];
<|"Normalization"->normalization, "TotalEnergy"->totalenergy,"FinalMass"->finalmass|>
]
