(* ::Package:: *)

HelperFunctions;


FixProcaSolution[solution_]:=
Block[{OutputSolution},
OutputSolution = solution;
OutputSolution[["Solution", "R"]]=solution[["Solution", "R(r)"]];
OutputSolution[["Solution", "S"]]={\[Theta]}|->Evaluate[Sum[C[tempiter]*SphericalHarmonicY[Abs[solution["Parameters", "m"]]+2*tempiter+solution["Parameters", "\[Eta]"],solution["Parameters", "m"],\[Theta],\[Phi]]*Exp[-I*solution["Parameters", "m"]*\[Phi]], {tempiter,0,solution["Parameters", "KMax"]}]/.solution["Solution", "AngularKernelVector"]];
Return[OutputSolution]
];


RenormalizeProcaSolution::usage="rescale a given proca solution by a given normalization factor";
RenormalizeProcaSolution[solution_, normalizationfactor_]:=
Block[{temporarySolution=solution},
(*Rescale radial function by normalization factor*)
temporaryRadialInterpolatingFunction = temporarySolution["Solution", "R"];
temporaryRadialInterpolatingFunctionList = List@@temporaryRadialInterpolatingFunction; (*Convert interpolating function to list*)
temporaryRadialInterpolatingFunctionList[[4]] = temporaryRadialInterpolatingFunctionList[[4]]*normalizationfactor; (*Rescale the 'ValuesOnGrid' Annotation of the Radial interpolating function*)
temporarySolution["Solution", "R"] = InterpolatingFunction@@temporaryRadialInterpolatingFunctionList; (*convert back to interpolating function*)
Return[temporarySolution];
];


AppendToSolution[solution_][name_, expr_]:=Block[{TemporarySolution,AssociationKey,printData},
TemporarySolution = solution;
If[TrueQ[Head[name]==String],
AssociationKey = name,
AssociationKey=ToString[name]
];
If[\[Not]KeyExistsQ[TemporarySolution, "Derived"],
TemporarySolution["Derived"]=<||>
];
TemporarySolution["Derived",AssociationKey]=expr;
printData = solution["Parameters"][[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s", "\[Chi]", "KMax", "branch"}]];
Export[$SolutionPath<>"RunData_"<>assocToString[printData]<>".mx", TemporarySolution];
]


styleMarkdown::usage = "pretty printing";
styleMarkdown[str_String] :=
    Module[{style, applyRules},
        style[s_String, {weight_:Plain, slant_:Plain}] := StringJoin[
            "\!\(\*StyleBox[\"", s, "\", FontWeight -> ", ToString @ weight, ", FontSlant -> ",
             ToString @ slant, "]\)"];
        applyRules[s_String, rules_] := StringReplace[s, # ~~ Shortest
             @ x__ ~~ # :> style[x, #2]& @@@ rules];
        applyRules[str, {"***" -> {Bold, Italic}, "**" -> {Bold}, "*"
             -> {Plain, Italic}}]
    ];


getResults::usage = "Retrieve results from disk. second argument is optional and specifies which parameters we wish to retrieve 
					example: getResults[NotebookDirectory[]<>'Solutions/', {{'m',1},{'n',1}}]
";
getResults[absolDir_?StringQ, parameterSet : (_?ListQ | Null) : Null] :=
    Module[{localfiles = FileNames[All, absolDir]},
        indexer = StringContainsQ[#, "RunData"]& /@ localfiles;
        filenames = Pick[localfiles, indexer];
        If[parameterSet == Null,
            Return[Import /@ filenames, Module]
        ];
        If[ListQ @ parameterSet,
            If[Length @ Dimensions @ parameterSet == 1,
                ParameterSetIndicies = StringContainsQ[ToString[parameterSet
                    [[1]]] <> ToString[parameterSet[[2]]]] /@ filenames;
                Return[Import /@ Pick[filenames, ParameterSetIndicies
                    ], Module];
            ];
            If[Length @ Dimensions @ parameterSet > 1,
                ParameterSetIndicies = Table[StringContainsQ[ToString[
                    parameterSet[[i]][[1]]] <> ToString[parameterSet[[i]][[2]]]] /@ filenames,
                     {i, 1, Length @ parameterSet}];
                ParameterSetIndiciesAll = Table[AllTrue[ParameterSetIndicies
                    [[All, i]], TrueQ], {i, 1, Length @ First @ ParameterSetIndicies}];
                Return[Import /@ Pick[filenames, ParameterSetIndiciesAll
                    ], Module];
            ];
        ];
    ];


assocToString::usage = "Convert association to string in format 'keyvalue'";
assocToString[assoc_] :=
    StringReplace[StringJoin @@ Table[ToString[(assoc // Keys)[[i]], 
        InputForm] <> ToString[(assoc // Values)[[i]], InputForm], {i, 1, Length[
        assoc // Keys]}], {"/" -> "_", "\"" -> ""}];


CacheResults::usage = "Save data in .mx format to disk in local directory";
CacheResults[assoc_, filename_] :=
    (Export[NotebookDirectory[] <> ToString[filename] <> ".mx", assoc
        ]);


R2toC[vec_] :=
    vec[[1]] + I * vec[[2]];


rplusN::usage = "outer horizon in terms of numerical variables";
rplusN[\[Chi]_]:=1+Sqrt[1-\[Chi]^2];


rminusN::usage = "inner horizon in terms of numerical variables";
rminusN[\[Chi]_]:=1-Sqrt[1-\[Chi]^2];


IterabletoList[iter_] :=
    Table[iter[[i]], {i, 1, iter // Length}];


DisjunctiontoList[disj_] :=
    Map[#[[2]]&, IterabletoList[disj]];


SymbolQ[sym_] :=
    MatchQ[Head @ sym, _Symbol];


Max\[Mu]::usage = "Compute the maximum gravitational coupling for superradiance to be efficient";
Max\[Mu][parameters_] :=
    parameters["m"] * parameters["\[Chi]"] / (2 * rplusN[parameters["\[Chi]"]]);


EmptyQ::usage = "Pattern test to check if list is empty. Used in plotting sector";
EmptyQ[expr_?ListQ] :=
    If[Length[expr] == 0,
        True
        ,
        False
    ];


ToxActVariables[expr_]:=expr/.{t->t[], r->r[], \[Theta]->\[Theta][], \[Phi]->\[Phi][]};
FromxActVariables[expr_]:=expr/.{t[]->t, r[]->r, \[Theta][]->\[Theta], \[Phi][]->\[Phi]};


OvertoneModeGroup[nvalue_?NumberQ, mvalue_?NumberQ]:=Select[Select[AllSolutions, #["Parameters"]["n"]==nvalue&], #["Parameters"]["m"]==mvalue&];


DetMet[r_,\[Theta]_,\[Chi]_]:= (r^2+\[Chi]^2*Cos[\[Theta]]^2)^2*Sin[\[Theta]]^2;


Kerrmetweight[r_,\[Theta]_,\[Chi]_]:=(r^2+\[Chi]^2*Cos[\[Theta]]^2)*Sin[\[Theta]];


TortoiseR::usage="Tortoise coordinates in terms of Boyer Lindquist coordinates";
TortoiseR[r_,\[Chi]_,M_]:=Block[{rp = M*rplusN[\[Chi]], rm = M*rminusN[\[Chi]]}, r + 2*M*rp/(rp-rm)*Log[(r-rp)/(2*M)] - 2*M*rm/(rp-rm)*Log[(r-rm)/(2*M)]];


OvertoneModeGroup[nvalue_?NumberQ, mvalue_?NumberQ]:=Select[Select[AllSolutions, #["Parameters"]["n"]==nvalue&], #["Parameters"]["m"]==mvalue&];


HorizonCoordToRadial[xN_, \[Chi]_]:=xN * (rplusN[\[Chi]] - rminusN[\[Chi]])  +  rplusN[\[Chi]];


RadialToHorizonCoord[rN_,\[Chi]_]:=(rN - rplusN[\[Chi]])/(rplusN[\[Chi]] - rminusN[\[Chi]]);


(*https://mathematica.stackexchange.com/questions/59944/extracting-the-function-from-interpolatingfunction-object/59963#59963*)
InterpolationToPiecewise[if_,x_]:=Module[{main,default,grid},grid=if["Grid"];
Piecewise[{if@"GetPolynomial"[#,x-#],x<First@#}&/@grid[[2;;-2]],if@"GetPolynomial"[#,x-#]&@grid[[-1]]]]/;if["InterpolationMethod"]=="Hermite";


(*Redundant?*)
CompileInterpolating[interp_InterpolatingFunction]:=Compile[{\[Psi]}, Evaluate[InterpolationToPiecewise[interp, \[Psi]]]];


ParamsToReprRule[solution_]:=AssociationThread[(Symbol/@(solution[ "Parameters"]//Keys))->(solution[ "Parameters"]//Values)]//Normal;


ToRational[expr_]:=Rationalize[expr,0];


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
);


ApplySolutionSet[solution_,OptionsPattern[{real->False}]][expr_]:= expr/.ToParamSymbols/.ParamsToReprRule[solution]/.SolToReprRule[solution,real->OptionValue[real]];


ToParamSymbols:={a->\[Chi]};


 Options[OptimizedFunction]={WithProperties->True,ToCompiled->False,OptimizationSymbol->aa, Options@Experimental`OptimizeExpression, Options@Compile}//Flatten;
OptimizedFunction[vars_, expr_, OptionsPattern[]]:=Module[{i,exproe, res,ret},
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
];


Options[ErgoRadius] = {coords->BL};
ErgoRadius[\[Theta]_,\[Chi]_, OptionsPattern[]]:=
If[TrueQ[OptionValue[coords]==BL],1+Sqrt[1-\[Chi]^2*Cos[\[Theta]]^2],
If[TrueQ[OptionValue[coords]==Horizon],RadialToHorizonCoord[1+Sqrt[1-\[Chi]^2*Cos[\[Theta]]^2],\[Chi]],
Print["Coords must be either BL (boyer-lindquist) or Horizon"]
]
];


SetAttributes[GenerateInterpolation, HoldFirst];
Options[GenerateInterpolation] = { DensityFunctions->Automatic, Metadata->False};
GenerateInterpolation[funcform_, args___, OptionsPattern[]]:=
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
DistributeDefinitions[func, CoordinateRanges];
funcOnPoints = With[{operand = {func,Sequence@@CoordinateRanges}},
						Parallelize[Outer@@operand]
					];

Messenger="Generating Interpolation function";
retval = ListInterpolation[funcOnPoints, CoordinateRanges];
If[OptionValue[Metadata],
Return[<|"Interpolation"->retval, "Mesh"->CoordinateRanges, "DensityFunctions"->densityfuncs|>],
Return[retval]
]
];


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
];
