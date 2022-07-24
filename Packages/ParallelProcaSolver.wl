(* ::Package:: *)

Needs["KerrWithProca`"]
$WorkingDirectory = "/home/shaunf/Documents/Computer/Code/projects/Massive_Vector_Field_Dynamical_Friction/ProcaAroundKerr/FKKSsolver/Mathematica/"
SolutionPath = $WorkingDirectory<>"Solutions/";
PlotPath = $WorkingDirectory<>"Plots/";
LogFile = $WorkingDirectory<>"RunLogs.txt";
PlotFileExtension = ".png";
ParameterDirectory = $WorkingDirectory<>"parameters.ini";
If[\[Not]DirectoryQ[SolutionPath],CreateDirectory[SolutionPath]]
If[\[Not]DirectoryQ[PlotPath],CreateDirectory[PlotPath]]


parameters = ToExpression/@KeyMap[ToString[ToExpression[#, InputForm, HoldForm]]&,Import[ParameterDirectory]]


(*Specify if iterations should produces logs during evaluation that get printed to output*)
WatchEvaluation = True;
TrackIteration=True;
LogRun = True;

(* ----- Prepare logging ----- *)
DeleteFile[LogFile];
CreateFile[LogFile];
If[LogRun,
LogStream = OpenWrite[LogFile];
WriteString[LogStream, DateString[]<>"\n\n"];
WriteString[LogStream, ToString[parameters, InputForm]<>"\n\n"];
];
Messenger = "\nExecuting Parameter space iteration"
Print[Messenger]


(* ----- Define Worker function ----- *)
WorkerFunction[mval_,nval_]:=
Block[{$HistoryLength=0},
With[{k=mval, j=nval},
 \[Mu]IterationCounter=1;
Do[

parameters["m"]=k;
parameters["l"]=parameters["m"]+parameters["s"];
parameters["n"]=j;

(*Setup next iteration*)
parameters["\[Mu]Nv"]=i*parameters["m"]; (*account for m-dependence of superradiant threshold*)
parameters["EndingX"]=calcEndingX[parameters];

Print["Kernel "<>ToString[$KernelID//N, InputForm]<>" executing search with \[Mu]="<>ToString[parameters["\[Mu]Nv"]//N, InputForm]<>" m="<>ToString[parameters["m"]]<>" n="<>ToString[parameters["n"]]];


\[Mu]IterationMessenger = " Current \[Mu] Iteration: "<>ToString[\[Mu]IterationCounter, InputForm]<>":"<>ToString[\[Mu]Iterations, InputForm];
\[Mu]Messenger= " Current \[Mu] value: "<>ToString[parameters["\[Mu]Nv"], InputForm];

Logger = "m Iteration Count: "<>ToString[mIterationCounter, InputForm]<>"\n"; 
Logger = Logger<> "m Value: "<>ToString[parameters["m"], InputForm]<>"\n"; 
Logger = Logger<>"n Value: "<>ToString[parameters["n"], InputForm]<>"\n";
Logger=Logger<>"\[Mu] Value: "<>ToString[parameters["\[Mu]Nv"], InputForm]<>"\n";

(* --------- Minimization using native Nelder Mead algorithm --------- *)
minimizetimestart = AbsoluteTime[];

(*Determine next guess for \[Omega] and \[Nu]*)
\[Nu]FitFunction[\[Mu]IterationCounter] = Null;
If[\[Mu]IterationCounter>=6,
(*If we're on the fifth iteration, use previous Q-values to approximate \[Omega]-value, instead of non-rel. limit*)
\[Omega]ValHistory = Table[\[Omega]Holder[iter], {iter,1,\[Mu]IterationCounter-1}];
\[Nu]ValHistory = Table[\[Nu]Holder[iter], {iter,1,\[Mu]IterationCounter-1}];
\[Mu]ValHistory = Table[\[Mu]Holder[iter], {iter,1,\[Mu]IterationCounter-1}];
\[Omega]Fit[\[Mu]IterationCounter] = Fit[Thread[{\[Mu]ValHistory, \[Omega]ValHistory}], {1,\[Zeta],\[Zeta]^2, \[Zeta]^3, \[Zeta]^4}, \[Zeta], WorkingPrecision->precision] (*InterpolatingPolynomial[Thread[{\[Mu]ValHistory, \[Omega]ValHistory}],\[Zeta]]*)(*Interp. Poly seems to fail for \[Omega]-values near \[Mu]=0.4*);
\[Nu]Fit[\[Mu]IterationCounter]=Fit[Thread[{\[Mu]ValHistory, \[Nu]ValHistory}],  {1,\[Zeta],\[Zeta]^2, \[Zeta]^3, \[Zeta]^4}, \[Zeta], WorkingPrecision->precision];
Current\[Omega]Guess = \[Omega]Fit[\[Mu]IterationCounter]/.\[Zeta]->parameters["\[Mu]Nv"];
Current\[Nu]Guess = \[Nu]Fit[\[Mu]IterationCounter]/.\[Zeta]->parameters["\[Mu]Nv"];
omegaGuess = Current\[Omega]Guess;
nuGuess = getNuValue[omegaGuess, parameters, Current\[Nu]Guess]//N;
\[Nu]FitFunction[\[Mu]IterationCounter] = Function[{\[Zeta]}, Evaluate@\[Nu]Fit[\[Mu]IterationCounter]]
,
(*use non-rel limit for initial guess of \[Omega] in first 4 iterations*)
omegaGuess = N[omegaNRNonRel[parameters] + I*omegaNINonRel[parameters]];
nuGuess = getNuValue[omegaGuess, parameters,nuNNonRel[omegaGuess,parameters]]//N;
];
omegaBoundarySize = 1/2*(omegaNRNonRel[ReplacePart[parameters,"n"->parameters["n"]+1]] - omegaGuess//Re);
omegaBoundary = {omegaGuess - omegaBoundarySize, omegaGuess+omegaBoundarySize};
(*Perform minimization*)
MinimizationResults=RadialMinimize[omegaGuess,  nuGuess, omegaBoundary,  parameters, WatchEvaluation,\[Nu]FitFunction[\[Mu]IterationCounter]];
Logger = Logger<> "Minimization Results: "<>ToString[MinimizationResults//N, InputForm];
minimizetimestop = AbsoluteTime[];
Logger = Logger<> "\nMinimization Time: "<>ToString[N[minimizetimestop-minimizetimestart], InputForm];



(* --------- Solve for kernel of angular matrix --------- *)
angulartimestart = AbsoluteTime[];
Sfunctemp = Sum[C[tempiter]*SphericalHarmonicY[Abs[parameters["m"]]+2*tempiter+parameters["\[Eta]"],parameters["m"],\[Theta],\[Phi]]*Exp[-I*parameters["m"]*\[Phi]], {tempiter,0,parameters["KMax"]}];
coeffs = SolveAngularSystem[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"],parameters];
Sfunc=FunctionInterpolation[Evaluate[Sfunctemp/.coeffs], {\[Theta],0,\[Pi]}];
angulartimestop = AbsoluteTime[];
Logger = Logger<>"\nAngular Solver Time: "<>ToString[N[angulartimestop-angulartimestart], InputForm];


(* --------- Safety checks on solution --------- *)
If[parameters["\[Mu]Nv"]^2<=Abs[MinimizationResults["\[Omega]"]]^2,
Print["Kernel "<>ToString[$KernelID, InputForm]<>" says: Error! Frequency larger than field mass! Breaking \[Mu] iteration ... (Parameter Point: m="<>ToString[parameters["m"], InputForm]<>" n="<>ToString[parameters["n"]]<>" \[Mu]="<>ToString[parameters["\[Mu]Nv"], InputForm]<>")"];
Break[];
];
If[TheSolution["R"][parameters["EndingX"]]//Log10>1,
Print["Kernel "<>ToString[$KernelID]<>" says: Error! Minimization failed! Breaking \[Mu] Iteration ... (Parameter Point: m="<>ToString[parameters["m"], InputForm]<>" n="<>ToString[parameters["n"]]<>" \[Mu]="<>ToString[parameters["\[Mu]Nv"], InputForm]<>")"];
Break[];
];

(* --------- Output solution and save to disk --------- *)
If[\[Mu]IterationCounter<6, 
\[Omega]fitOutVariable="Not enough iterations to construct fit";
\[Nu]fitOutVariable="Not enough iterations to construct fit";, 

\[Nu]fitOutVariable=\[Nu]Fit//DownValues;
 \[Omega]fitOutVariable=\[Omega]Fit//DownValues;
];
TheSolution = <|"\[Omega]" -> MinimizationResults["\[Omega]"],
			    "\[Nu]"->  MinimizationResults["\[Nu]"],
		    	"R"-> getRadialSolution[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"], parameters],
		             "S"-> Sfunc,
			    "AngularKernelVector"-> coeffs,
			    "Q"-> QValue[parameters["\[Mu]Nv"], MinimizationResults["\[Omega]"]],
			    "\[Omega]fit"-> \[Omega]fitOutVariable,
			    "\[Nu]fit"-> \[Nu]fitOutVariable
			|>;
RadialPlot=LogLinearPlot[
TheSolution["R"][r]//Abs, 
{r,parameters["StartingX"], parameters["EndingX"]}, 
PlotLabel->"|R(\!\(\*SubscriptBox[\(x\), \(N\)]\))|",
 PlotRange->All,
 AxesLabel->{"\!\(\*SubscriptBox[\(x\), \(N\)]\)", "R(\!\(\*SubscriptBox[\(x\), \(N\)]\))"},
 ImageSize->Large, 
Epilog-> Inset[Style[Grid[Partition[Table[(parameters//Keys)[[tempiter]] -> (parameters//Values)[[tempiter]], {tempiter,1, Length@parameters}],2], ItemSize->18], 10]]
];
AngularPlot =Plot[
TheSolution["S"][\[Theta]]//Abs, 
{\[Theta],0,\[Pi]}, 
PlotLabel->"|S(\[Theta])|", 
AxesLabel->{"\[Theta]", ""},
ImageSize->Large,
Epilog-> Inset[Style[Grid[Partition[Table[(parameters//Keys)[[tempiter]] -> (parameters//Values)[[tempiter]], {tempiter,1, Length@parameters}],2], ItemSize->18], 10]]
];

AppendTo[TheSolution, "RadialPlot"-> RadialPlot];
AppendTo[TheSolution, "AngularPlot"-> AngularPlot];

(*Prepare solution data and write to disk*)
printData=parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s", "\[Chi]", "KMax", "branch"}]];
CacheData = <| "Parameters"-> parameters, "Solution"-> TheSolution|>;
Logger = Logger<> "\nCaching results to: "<>assocToString[printData]<>".mx\n\n";
Export[SolutionPath<>"RunData_"<>assocToString[printData]<>".mx", CacheData];
If[LogRun,Logger>>>LogFile];

(* --------- Prepare next loop --------- *);
\[Mu]Holder[\[Mu]IterationCounter]=parameters["\[Mu]Nv"]//SetPrecision[#,parameters["precision"]]&;
\[Omega]Holder[\[Mu]IterationCounter]=MinimizationResults["\[Omega]"]//SetPrecision[#,parameters["precision"]]&;
\[Nu]Holder[\[Mu]IterationCounter]=MinimizationResults["\[Nu]"]//SetPrecision[#,parameters["precision"]]&;

\[Mu]IterationCounter++;
,
{i,parameters["\[Mu]range"]//First, parameters["\[Mu]range"]//Last,parameters["\[Delta]\[Mu]"]}];
ClearAll[\[Mu]Holder, \[Omega]Holder, \[Nu]Holder, TheSolution, AngularPlot, RadialPlot,CacheData, Logger, \[Omega]Fit, \[Nu]Fit,\[Nu]FitFunction];
ClearSystemCache[];
]
]


(* ----- Setup Parallelization ----- *)
CloseKernels[];
If[WatchEvaluation, Print["Launching slave kernels..."]];
LaunchKernels[Min[$ProcessorCount,(Last@parameters["mrange"] - First@parameters["mrange"] + 1)*(Last@parameters["nrange"] - First@parameters["nrange"] + 1)]];
If[WatchEvaluation, Print["Done. Launched "<>ToString[$KernelCount]<>" slave kernels."]]
ParallelEvaluate[Needs["KerrWithProca`"]];
DistributeDefinitions[parameters, WorkerFunction,SolutionPath];

(*Turn off irrelevant error messages*)
Off[ClebschGordan::phy];ParallelEvaluate[Off[ClebschGordan::phy]];
Off[NMinimize::nnum];ParallelEvaluate[Off[NMinimize::nnum]];


(* ----- Execute Parallelized Iterations ----- *)
runStartTime = AbsoluteTime[];
(Evaluators = Flatten@Table[ParallelSubmit[{mv, nv},WorkerFunction[mv,nv]], 
			{nv, parameters["nrange"]//First, parameters["nrange"]//Last}, 
			{mv, parameters["mrange"]//First, parameters["mrange"]//Last}
			])//Partition[#, $KernelCount/2//Floor]&//Grid
WaitAll[Evaluators]

runEndTime=AbsoluteTime[];
Print[Style["Finished! Total Time: "<>ToString[N@(Round[runEndTime-runStartTime]/60), InputForm]<>" Minutes",Blue, Italic, 24]]
(*Clean up*)
CloseKernels[];
