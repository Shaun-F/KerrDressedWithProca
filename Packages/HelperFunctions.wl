(* ::Package:: *)

(*File title to know if package has already been imported*)
HelperFunctions;


Options[SolutionToFilename]={Prelabel->Null, KeyList->{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s", "\[Chi]", "KMax", "branch"}, ExtensionType->".mx"};
SolutionToFilename[Assoc_Association, parentdirectory_, OptionsPattern[]]:=
Block[{printData, PrelabelString, parametersection,fileName,absoluteFileName},
If[KeyExistsQ[Assoc, "Parameters"],
printData = Assoc["Parameters"][[OptionValue[KeyList]]];,
printData = Assoc[[OptionValue[KeyList]]];
];
If[TrueQ[OptionValue[Prelabel]==Null],
PrelabelString = "",
PrelabelString = OptionValue[Prelabel]
];
parametersection = assocToString[printData];
fileName = PrelabelString<>"RunData_"<>parametersection<>OptionValue[ExtensionType];
absoluteFileName = FileNameJoin[{parentdirectory,fileName}];
Return[absoluteFileName]
]


RecastInterpolationFunction::usage="Recast interpolating function over a new domain obtained by operating on the old domain with DomainFunction"
RecastInterpolationFunction[int_InterpolatingFunction, DomainFunction_]:=
Block[{GridPoints = int["Grid"]//Flatten, NewGridPoints, NewData},
NewGridPoints = DomainFunction/@GridPoints;
NewData = Thread[{NewGridPoints, int/@GridPoints}];
Interpolation[NewData, InterpolationOrder->int["InterpolationOrder"], Method->int["InterpolationMethod"]]
]


TruncateInterpolatingFunction[int_InterpolatingFunction, TruncatingFunction_, OptionsPattern[{OnValues->False, OnDomain->True}]]:=
Block[{interp = int, interpdomain,NewInterp, NewValuesAndGrid = {}},
interpdomain = interp["Grid"]//Flatten;
If[OptionValue[OnValues],
For[i=1, i<=Length[interpdomain], i++, 
If[int[interpdomain[[i]]]//TruncatingFunction,
AppendTo[NewValuesAndGrid, {interpdomain[[i]],interp[interpdomain[[i]]]}];
];(*If statmement*)
];(*For loop*)
];(*If statement*)

If[OptionValue[OnDomain],
For[i=1, i<=Length[interpdomain], i++,
If[interpdomain[[i]]//TruncatingFunction,
AppendTo[NewValuesAndGrid, {interpdomain[[i]],interp[interpdomain[[i]]]}];
];(*If statement*)
](*For loops*)
];(*If statement*)

NewInterp = Interpolation[NewValuesAndGrid, Method->interp["InterpolationMethod"], InterpolationOrder->interp["InterpolationOrder"]]
]


FixProcaSolution[solution_]:=
Block[{OutputSolution},
OutputSolution = solution;
OutputSolution[["Solution", "R"]]=solution[["Solution", "R(r)"]];
OutputSolution[["Solution", "S"]]={\[Theta]}|->Evaluate[Sum[C[tempiter]*SphericalHarmonicY[Abs[solution["Parameters", "m"]]+2*tempiter+solution["Parameters", "\[Eta]"],solution["Parameters", "m"],\[Theta],\[Phi]]*Exp[-I*solution["Parameters", "m"]*\[Phi]], {tempiter,0,solution["Parameters", "KMax"]}]/.solution["Solution", "AngularKernelVector"]];
Return[OutputSolution]
];


OperateOnInterpolationValues::usage="Operate on the values inputted into the interpolation function with a listable function";
OperateOnInterpolationValues[interp_, function_]:=
Block[{Out},
Out = List@@interp;
Out[[4]]=Out[[4]]//function;
Return[InterpolatingFunction@@Out]
]


RenormalizeProcaSolution::usage="rescale a given proca solution by a given normalization factor";
RenormalizeProcaSolution[solution_, normalizationfactor_]:=
Block[{temporaryRadialInterpolatingFunction,temporaryRadialInterpolatingFunctionList, temporarySolution=solution},
(*Rescale radial function by normalization factor*)
temporaryRadialInterpolatingFunction = temporarySolution["Solution", "R"];
temporaryRadialInterpolatingFunctionList = List@@temporaryRadialInterpolatingFunction; (*Convert interpolating function to list*)
temporaryRadialInterpolatingFunctionList[[4]] = temporaryRadialInterpolatingFunctionList[[4]]*normalizationfactor; (*Rescale the 'ValuesOnGrid' Annotation of the Radial interpolating function*)
temporarySolution["Solution", "R"] = InterpolatingFunction@@temporaryRadialInterpolatingFunctionList; (*convert back to interpolating function*)
Return[temporarySolution];
];


AppendToSolution[solution_, OptionsPattern[{Prelabel->Null}]][NameValuePair_List]:=Block[{TemporarySolution,AssociationKey,AbsFilename},
TemporarySolution = solution;
If[\[Not]KeyExistsQ[TemporarySolution, "Derived"],
TemporarySolution["Derived"]=<||>
];
If[Length@Dimensions@NameValuePair>1,
(*True*)
Do[
TemporarySolution["Derived",ToString[NameValuePair[[i]][[1]]]]=NameValuePair[[i]][[2]];,
{i,1,Length@NameValuePair}
];,

(*Else*)
TemporarySolution["Derived",ToString[NameValuePair[[1]]]]=NameValuePair[[2]];
];

AbsFilename = SolutionToFilename[solution, $SolutionPath, Prelabel->OptionValue[Prelabel]];
Export[AbsFilename, TemporarySolution];
];


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
					example: getResults[NotebookDirectory[]<>'Solutions/', {{'m',1},{'n',1}, {'\[Mu]Nv', '1_20'}}]
";
getResults[absolDir_?StringQ, parameterSet : (_?ListQ | Null) : Null] :=
	Block[{m,n,\[Chi],\[Mu]Nv,\[Delta]\[Mu],\[Delta]n,\[Delta]m,\[Epsilon],eta,mode,overtone,polarization,s,orbitalnumber,l,\[Chi]v,KMax,branch, indexer, filenames,DataList,ParameterSetIndicies, ParameterSetIndiciesAll, CleansedDataList,MassOrderingIndicies},   
     Module[{localfiles = FileNames[All, absolDir]},
        indexer = StringContainsQ[#, "RunData"]& /@ localfiles;
        filenames = Pick[localfiles, indexer];
        If[parameterSet == Null,
            DataList = Import /@ filenames;
        ];
        If[ListQ @ parameterSet,
            If[Length @ Dimensions @ parameterSet == 1,
                ParameterSetIndicies = Thread[StringContainsQ[ToString[parameterSet[[1]]] <> ToString[parameterSet[[2]]]] /@ filenames || StringContainsQ[ToString[FullForm[ToExpression@parameterSet[[1]]]] <> ToString[parameterSet[[2]]]]/@filenames];
                DataList = Import /@ Pick[filenames, ParameterSetIndicies];
            ];
            If[Length @ Dimensions @ parameterSet > 1,
                ParameterSetIndicies = Table[
                                             Thread[StringContainsQ[ToString[parameterSet[[i]][[1]]] <> ToString[parameterSet[[i]][[2]]]] /@ filenames || StringContainsQ[ToString[FullForm[ToExpression@parameterSet[[i]][[1]]]] <> ToString[parameterSet[[i]][[2]]]]/@filenames],
                                             {i, 1, Length @ parameterSet}
                                             ];
                ParameterSetIndiciesAll = Table[AllTrue[ParameterSetIndicies[[All, i]], TrueQ], {i, 1, Length @ First @ ParameterSetIndicies}];
                DataList = Import /@ Pick[filenames, ParameterSetIndiciesAll];
            ];
        ];
        CleansedDataList = Apply[Intersection, 
                                 DeleteCases[DataList, #]&/@ {
                                                             <|___, "Solution"-><|___, "\[Omega]"->s_/;Not[NumericQ[s]], ___|>,___|>, 
                                                             <|___, "Solution"->s_/;StringQ[s], ___|>
                                                             }
                                 ];
        MassOrderingIndicies = Ordering@CleansedDataList[[All, "Parameters", "\[Mu]Nv"]];
        CleansedDataList[[MassOrderingIndicies]]
    ]
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


OvertoneModeGroup[nvalue_?NumberQ, mvalue_?NumberQ, SolutionSet_]:=Select[Select[SolutionSet, #["Parameters"]["n"]==nvalue&], #["Parameters"]["m"]==mvalue&];


DetMet[r_,\[Theta]_,\[Chi]_]:= (r^2+\[Chi]^2*Cos[\[Theta]]^2)^2*Sin[\[Theta]]^2;


Kerrmetweight[r_,\[Theta]_,\[Chi]_]:=(r^2+\[Chi]^2*Cos[\[Theta]]^2)*Sin[\[Theta]];


TortoiseR::usage="Tortoise coordinates in terms of Boyer Lindquist coordinates";
TortoiseR[r_,\[Chi]_,M_]:=Block[{rp = M*rplusN[\[Chi]], rm = M*rminusN[\[Chi]]}, r + 2*M*rp/(rp-rm)*Log[(r-rp)/(2*M)] - 2*M*rm/(rp-rm)*Log[(r-rm)/(2*M)]];


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
											((Derivative[2][Rv][RadialToHorizonCoord[#,\[Chi]v]]*1/(rplusN[\[Chi]v]-rminusN[\[Chi]v])^2)&)
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


ApplySolutionSet[solution_,OptionsPattern[{real->False}]][expr_]:= expr/.ToParamSymbols/.ParamsToReprRule[solution]/.SolToReprRule[solution,real->OptionValue[real]];


Options[ApplyRealSolutionSet]={QuasiboundState->False, ProperDerivative->True};
ApplyRealSolutionSet[solution_, OptionsPattern[]][expr_]:=
	Block[{repr, res},
	With[{solR=solution["Solution","R"], solS=solution["Solution", "S"], \[Chi]v=solution["Parameters", "\[Chi]"]},
	(*postfactor = 1/(rplusN[\[Chi]v]-rminusN[\[Chi]v]);*)
		With[{MapToRadial =  (RadialToHorizonCoord[#,\[Chi]v]&),postfactor=1/(rplusN[\[Chi]v]-rminusN[\[Chi]v])},
		repr = {
			HoldPattern@Derivative[d_][Rr][r_]:>Re[Derivative[d][solR][MapToRadial[r]]]*postfactor^d,
			HoldPattern@Derivative[d_][Ri][r_]:>Im[Derivative[d][solR][MapToRadial[r]]]*postfactor^d,
			HoldPattern@Derivative[d_][Sr][x_]:>Re[Derivative[d][solS][x]],
			HoldPattern@Derivative[d_][Si][x_]:>Im[Derivative[d][solS][x]],
			HoldPattern@Derivative[d_][Rr]:>(Re[Derivative[d][solR][MapToRadial[#]]]*postfactor^d&),
			HoldPattern@Derivative[d_][Ri]:>(Im[Derivative[d][solR][MapToRadial[#]]]*postfactor^d&),
			HoldPattern@Derivative[d_][Sr]:>(Re[Derivative[d][solS][#]]&),
			HoldPattern@Derivative[d_][Si]:>(Im[Derivative[d][solS][#]]&),
			\[Omega]r -> Re[solution["Solution","\[Omega]"]],
			\[Omega]i->Evaluate@If[OptionValue[QuasiboundState], 0, Im[solution["Solution","\[Omega]"]]],
			\[Nu]r -> Re[solution["Solution","\[Nu]"]], 
			\[Nu]i->Im[solution["Solution","\[Nu]"]],
			\[Chi]->solution["Parameters", "\[Chi]"], 
			a->solution["Parameters", "\[Chi]"],
			\[Chi]v->solution["Parameters", "\[Chi]"],
			m->solution["Parameters", "m"], 
			\[Mu]Nv->solution["Parameters", "\[Mu]Nv"],
			\[Mu]->solution["Parameters", "\[Mu]Nv"],
			n->solution["Parameters", "n"],
			Rr->( Re[solR[MapToRadial[#]]]&), 
			Ri->( Im[solR[MapToRadial[#]]]&), 
			Sr->( Re[solS[#]]&),
			Si->( Im[solS[#]]&)
			};
		res = expr/.repr;
		Return[res]
		];
		];
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


ToParamSymbols={a->\[Chi]};


Options[OptimizedFunction]={WithProperties->False,ToCompiled->False,OptimizationSymbol->aa, Options@Experimental`OptimizeExpression, Options@Compile}//Flatten;
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
ret = Compile@@Hold[thevars, Evaluate[resc]@@thevars,opts];
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
Options[GenerateInterpolation] = { DensityFunctions->Automatic, Metadata->False}\[Union]Options[ListInterpolation];
GenerateInterpolation[funcform_, args___, OptionsPattern[]]:=
Block[{
func,
variables = Hold/@List[args][[All,1]],
StartPoints = List[args][[All,2]],
EndPoints = List[args][[All,3]],
NArgs = Length@List@args,
funcOnPoints,
densityfuncs,
FunctionDimension, Switcher,
InterpList,
scalar, vectorial,
RecastDensityToBoundary,
recastfunctions,
CoordinateRanges,
Generator,
funcOnPointsReal,
funcOnPointsImag,
interpRe,
interpIm,
retval
},

func = Hold[funcform][[1,0]];(*Extract head*)
RecastDensityToBoundary[Identity,_,_][x_]:=x;
RecastDensityToBoundary[denfunc_,rstart_,rstop_][x_]:=denfunc[x]/((denfunc[rstop]-denfunc[rstart])/(rstop-rstart))+(rstart*denfunc[rstop]-rstop*denfunc[rstart])/(denfunc[rstop] -denfunc[rstart]);
SetAttributes[RecastDensityToBoundary, Listable];
SetAttributes[func, Listable];

(*if input function is vectorial or scalar*)
FunctionDimension = Length[func@@StartPoints];
If[TrueQ[FunctionDimension>0],
	Switcher = vectorial,
	Switcher = scalar
];

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
funcOnPoints = N@With[{operand = {func,Sequence@@CoordinateRanges}},
						Parallelize[Outer@@operand]
					];			
Switch[Switcher,

scalar,			
Messenger="Generating Interpolation function for scalar function";
Generator = ListInterpolation[#1, CoordinateRanges, Method->OptionValue[Method], InterpolationOrder->OptionValue[InterpolationOrder], PeriodicInterpolation->OptionValue[PeriodicInterpolation]]&;
If[AnyTrue[funcOnPoints, (TrueQ[Head[#]==Complex] &), Depth[funcOnPoints]-1] && TrueQ[OptionValue[Method]=="Spline"],
funcOnPointsReal = Re[funcOnPoints];
funcOnPointsImag = Im[funcOnPoints];
interpRe = Generator@funcOnPointsReal;
interpIm = Generator@funcOnPointsImag;
With[{unheldvars = ReleaseHold/@variables},
retval = Evaluate[unheldvars] |->Evaluate[interpRe[Sequence@@unheldvars] + I*interpIm[Sequence@@unheldvars]];
];,
retval = Generator@funcOnPoints;
];,


vectorial,
Messenger="Generating Interpolation function for Vectorial function";
InterpList = {};
Do[
With[{SubspaceFunctionPoints = Map[Part[#,i]&, funcOnPoints, {-(2)}]},
Generator = ListInterpolation[#1, CoordinateRanges, Method->OptionValue[Method], InterpolationOrder->OptionValue[InterpolationOrder], PeriodicInterpolation->OptionValue[PeriodicInterpolation]]&;
If[AnyTrue[SubspaceFunctionPoints, (TrueQ[Head[#]==Complex] &), Depth[SubspaceFunctionPoints]-1] && TrueQ[OptionValue[Method]=="Spline"],
funcOnPointsReal = Re[SubspaceFunctionPoints];
funcOnPointsImag = Im[SubspaceFunctionPoints];
interpRe = Generator@SubspaceFunctionPoints;
interpIm = Generator@SubspaceFunctionPoints;
With[{unheldvars = ReleaseHold/@variables},
AppendTo[InterpList, Evaluate[unheldvars] |->Evaluate[interpRe[Sequence@@unheldvars] + I*interpIm[Sequence@@unheldvars]]];
];,
AppendTo[InterpList,Generator@SubspaceFunctionPoints];
];
];(*End of With Statement*), {i, 1,FunctionDimension}
];(*End of Do*)
retval = InterpList;
];(*End of Switch*)

If[OptionValue[Metadata],
Return[<|"Interpolation"->retval, "Mesh"->CoordinateRanges, "DensityFunctions"->densityfuncs|>
],
Return[retval]
]

];
