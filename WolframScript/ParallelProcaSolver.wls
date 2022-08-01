#!/usr/bin/env wolframscript
(* ::Package:: *)

$FKKSRoot = "/home/fell/Documents/Github/KerrDressedWithProca/";
(*$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";*)
LogFile = $FKKSRoot<>"Logs/LogFile.log";

ConfigurationParameters = Import[$FKKSRoot<>"WolframScript/Config.ini"];
LogRun = ToExpression[ConfigurationParameters["LogRun"]];
SaveToDisk = ToExpression[ConfigurationParameters["SaveToDisk"]];
NumberOfKernels = ConfigurationParameters["NumberofSlaveKernels"]//ToExpression;

If[!TrueQ[(ConfigurationParameters["SolutionFilePath"]==Automatic)],
	$SolutionPath = ConfigurationParameters["SolutionFilePath"];,
	$SolutionPath = $FKKSRoot<>"Solutions/";
];

Print["FKKSRoot: "<>$FKKSRoot];
Print["Solution Path: "];Print[$SolutionPath];
Print["Log File Path: "];Print[LogFile];
Print["/nLog Run: "<>ToString[TrueQ[LogRun], InputForm]];
Print["Save to disk: "<>ToString[TrueQ[SaveToDisk], InputForm]];


Import[FileNameJoin[{$FKKSRoot,"Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];

parameterlist = Import[$FKKSRoot<>"WolframScript/parameters.ini"]//ToExpression;
parameters = KeyMap[ToString[ToExpression[#]]&, parameterlist];
parameters["StartingX"] = parameters["\[Epsilon]"];


If[FileExistsQ[LogFile], DeleteFile[LogFile]];
CreateFile[LogFile];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];
Messenger = "\n\n-----------------------\nBeginning Parallelized Proca Solver\n";
Print[Messenger];

Options[WorkerFunction]={OverwritePrevious->False};
WorkerFunction[mval_,nval_, OptionsPattern[]]:=
	Block[{HistoryLength=0},
		With[{k=mval, j=nval},
			muIterationCounter=1;
			Logger = {"Beginning search routine:   time:" <> DateString[] <>" Parameter values: (m,n,\[Chi]) = (" <> ToString[k, InputForm] <>"," <> ToString[j, InputForm] <> "," <> ToString[parameters["\[Chi]"], InputForm] <> ")"}
			Do[
				(*Prepare parameter values*)
				parameters["m"]=k;
				parameters["l"]=parameters["m"]+parameters["s"];
				parameters["n"]=j;
				
				(*Setup iteration*)	
				parameters["\[Mu]Nv"]=i*parameters["m"]; (*account for m-dependence of superradiant threshold*)
				parameters["EndingX"]=calcEndingX[parameters]; (*Determine cutoff radius for radial equation integrator*)
				
				(*Construct name for file to save data in*)
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s", "\[Chi]", "KMax", "branch"}]];
				SolutionFileName = $SolutionPath<>"RunData_"<>assocToString[printData]<>".mx";

				(*If OverwritePrevious option is false and if the solution file already exists, skip current loop iteration*)
				If[! OptionValue[OverwritePrevious],
					If[FileExistsQ[SolutionFileName],
						Continue[]
					];
				];
				
				(*Output current iteration to terminal*)
				Print["Kernel "<>ToString[$KernelID]<>" working on parameter point (\[Mu], m, n) = ("<>ToString[parameters["\[Mu]Nv"], InputForm]<>", "<>ToString[parameters["m"], InputForm]<>" ,"<>ToString[parameters["n"], InputForm]<>")"];
				
				AppendTo[Logger,"m Iteration Count: "<>ToString[mIterationCounter, InputForm]];
				AppendTo[Logger,"m value: "<>ToString[parameters["m"], InputForm]];
				AppendTo[Logger,"n value: "<>ToString[parameters["n"], InputForm]];
				AppendTo[Logger,"mu value: "<>ToString[parameters["\[Mu]Nv"], InputForm]];
				
				(*Minimization using native Nelder Mead algorithm*)
				minimizetimestart = AbsoluteTime[];
				
				(*Determine next guess for frequency and angular eigenvalue*)
				nuFitFunction[muIterationCounter] = Null;
				If[TrueQ[muIterationCounter >= 6],
					
					(*If we're on the 6th iteration, use previous values to approximate current frequency value using fit functions*)
					omegaValHistory = Table[omegaHolder[iter], {iter, 1, muIterationCounter-1}];
					nuValHistory = Table[nuHolder[iter], {iter, 1, muIterationCounter-1}];
					muValHistory = Table[muHolder[iter], {iter, 1, muIterationCounter-1}];
					omegaFit[muIterationCounter] = Interpolation[Thread[{muValHistory, omegaValHistory}], InterpolationOrder->2];
					nuFit[muIterationCounter] = Interpolation[Thread[{muValHistory, nuValHistory}], InterpolationOrder->2];
					CurrentomegaGuess = omegaFit[muIterationCounter][parameters["\[Mu]Nv"]];
					CurrentnuGuess = nuFit[muIterationCounter][parameters["\[Mu]Nv"]];
					omegaGuess = CurrentomegaGuess;
					nuGuess = getNuValue[omegaGuess, parameters, CurrentnuGuess];
					nuFitFunction[muIterationCounter] = nuFit[muIterationCounter]; (*To be used later for the radial minimization*)
					,
					(*Use non-relativistic limit for initial guess of omega in first 5 iterations*)
					omegaGuess = N[omegaNRNonRel[parameters] + I*omegaNINonRel[parameters]];
					nuGuess = getNuValue[omegaGuess, parameters, nuNNonRel[omegaGuess, parameters]];
				];
				
				(*If guess for imaginary part is negative, bound state is not superradiantly growing -> skip iteration*)
				If[Im[omegaGuess] < 0, 
					Print["Superradiant condition fails for parameter (\[Mu],m,n,\[Chi],s)=("<>ToString[parameters["\[Mu]Nv"], InputForm] <> "," <>ToString[parameters["m"], InputForm] <> "," <>ToString[parameters["n"], InputForm] <> "," <>ToString[parameters["\[Chi]"], InputForm] <> "," <>ToString[parameters["s"], InputForm] <> ")"];
					Goto[SuperradiantConditionFailure];
				];
				
				omegaBoundarySize = (1/2)*(omegaNRNonRel[ReplacePart[parameters, "n"->parameters["n"]+1]] - Re[omegaGuess]); (*Use next overtone frequency to construct boundary for omega search*)
				omegaBoundary = {omegaGuess - omegaBoundarySize, omegaGuess + omegaBoundarySize};
				
				(*Perform Minimization*)
				MinimizationResults = RadialMinimize[omegaGuess, nuGuess, omegaBoundary, parameters, False, nuFitFunction[muIterationCounter]];
				AppendTo[Logger,"Minimization Results: "<>ToString[MinimizationResults//N, InputForm]];
				minimizetimestop = AbsoluteTime[];
				AppendTo[Logger,"Minimization Time: "<>ToString[N[minimizetimestop-minimizetimestart], InputForm]];
				
				
				(*Solve for Kernel of angular matrix*)
				angulartimestart = AbsoluteTime[];
				Sfunctemp = Sum[C[tempiter]*SphericalHarmonicY[Abs[parameters["m"]] + 2*tempiter + parameters["\[Eta]"],parameters["m"],\[Theta], \[Phi]]*Exp[-I*parameters["m"]*\[Phi]], {tempiter,0, parameters["KMax"]}];
				coeffs = SolveAngularSystem[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"],parameters];
				Sfunc = Function[{\[Theta]},Evaluate[Sfunctemp/.coeffs]];
				angulartimestop = AbsoluteTime[];
				AppendTo[Logger,"Angular Solver Time: "<>ToString[N[angulartimestop-angulartimestart], InputForm]];
				
				
				(*Format solution*)
				If[muIterationCounter < 6,
					omegafitOutVariable = "Not enough iterations to construct fit";
					nufitOutVariable = "Not enough iterations to construct fit";,
					
					nufitOutVariable = nuFit//DownValues;
					omegafitOutVariable = omegaFit//DownValues;
				];
				RadialSolution = getRadialSolution[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"], parameters];
				
				RadialPlot = LogLinearPlot[
								RadialSolution[r]//Abs,
								{r, parameters["StartingX"], parameters["EndingX"]},
								PlotLabel->"|R( \!\(\*SubscriptBox[\(x\), \(n\)]\) )|",
								PlotRange->All,
								AxesLabel->{"\!\(\*SubscriptBox[\(x\), \(n\)]\)", "R( \!\(\*SubscriptBox[\(x\), \(n\)]\) )"},
								ImageSize->Large,
								Epilog->Inset[Style[Grid[Partition[Table[(parameters//Keys)[[tempiter]]->(parameters//Values)[[tempiter]], {tempiter, 1, Length@parameters}],2], ItemSize->18], 10]]
							];
				AngularPlot = LogLinearPlot[
								Sfunc[theta]//Abs,
								{theta, 0, \[Pi]},
								PlotLabel->"|S( \[Theta] )|",
								PlotRange->All,
								AxesLabel->{"\[Theta]",""},
								ImageSize->Large,
								Epilog->Inset[Style[Grid[Partition[Table[(parameters//Keys)[[tempiter]]->(parameters//Values)[[tempiter]], {tempiter, 1, Length@parameters}],2], ItemSize->18], 10]]
							];
				
				TheSolution = <|"\[Omega]" -> MinimizationResults["\[Omega]"],
						"\[Nu]"-> MinimizationResults["\[Nu]"],
						"R"->RadialSolution,
						"S"->Sfunc,
						"AngularKernelVector"->coeffs,
						"Q"->QValue[parameters["\[Mu]Nv"], MinimizationResults["\[Omega]"]],
						"\[Omega]fit"->omegafitOutVariable,
						"\[Nu]fit"->nufitOutVariables,
						"RadialPlot"->RadialPlot,
						"AngularPlot"->AngularPlot
						|>;
						
				(*Safety checks on solution*)
				If[TrueQ[parameters["\[Mu]Nv"]^2 <= Abs[MinimizationResults["\[Omega]"]]^2],
					Print["Kernel "<>ToString[$KernelID,InputForm]<>" says: Error! Frequency larger than field mass! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];                
					Break[];
				];
				If[TrueQ[TheSolution["R"][parameters["EndingX"]]//Log10 > 1],
					Print["Kernel "<>ToString[$KernelID]<>" says: Error! Minimization failed! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];
					Break[];
				];
						
				(*Print solution to terminal*)
				Print["Result for Kernel "<>ToString[$KernelID, InputForm]<>": \n\t \[Omega] = "<>ToString[Evaluate@MinimizationResults["\[Omega]"], InputForm]];
				
				(*Prepare solution data and write to disk*)
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s","\[Chi]", "KMax", "branch"}]];
				CacheData = <| "Parameters" -> parameters, "Solution" -> TheSolution|>;
				AppendTo[Logger,"Caching results to: "<>assocToString[printData]<>".mx"];
				AppendTo[Logger, ""];AppendTo[Logger, ""];
				If[SaveToDisk,
					Export[$SolutionPath<>"RunData_"<>assocToString[printData]<>".mx",CacheData];
				];
				If[LogRun, 
					PutAppend[Sequence@@Logger, LogFile];
				];
				
				(*Prepare next loop*)
				muHolder[muIterationCounter] = parameters["\[Mu]Nv"]//SetPrecision[#, parameters["precision"]]&;
				omegaHolder[muIterationCounter] = MinimizationResults["\[Omega]"]//SetPrecision[#, parameters["precision"]]&;
				nuHolder[muIterationCounter] =  MinimizationResults["\[Nu]"]//SetPrecision[#, parameters["precision"]]&;
				
				muIterationCounter++;
				Continue[];
				
				Label[SuperradiantConditionFailure];
				AppendTo[Logger,"Superradiant condition failed for parameters (\[Mu],m,n,\[Chi],s)=(" <>ToString[parameters["\[Mu]Nv"], InputForm] <> "," <>ToString[parameters["m"], InputForm] <> "," <>ToString[parameters["n"], InputForm] <> "," <>ToString[parameters["\[Chi]"], InputForm] <> "," <>ToString[parameters["s"], InputForm] <> ")"];
				If[LogRun, PutAppend[Sequence@@Logger,LogFile]];	
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s","\[Chi]", "KMax", "branch"}]];
				CacheData = <|"Parameters" -> parameters, "Solution" -> "Superradiant condition failed"|>;
				If[SaveToDisk,Export[$SolutionPath <> "RunData_" <> assocToString[printData]<>".mx", CacheData];];
			
					
				,
				{i, parameters["\[Mu]range"]//First, parameters["\[Mu]range"]//Last, parameters["\[Delta]\[Mu]"]}
			]; (*mu iteration loop end*)
			
			(*Clear symbols used in mu iteration to prevent memory leaks*)
			ClearAll[muHolder, omegaHolder, nuHolder, TheSolution, AngularPlot, RadialPlot, CacheData, Logger, omegaFit, nuFit, nuFitFunction];
			ClearSystemCache[];
		]; (*With statement end*)
	];(*Block statement end*)
	(*End of worker function definition*)



(*Setup parallelization*)
CloseKernels[];
NumberKernels = Min[NumberOfKernels, (Last[parameters["mrange"]] - First[parameters["mrange"]]+1)*(Last@parameters["nrange"] - First@parameters["nrange"] + 1)];
Print["Launching "<>ToString[NumberKernels,InputForm]<>" Kernels."];
LaunchKernels[NumberKernels];
DistributeDefinitions[parameters, WorkerFunction, $SolutionPath, $FKKSRoot, LogFile, LogRun];
ParallelEvaluate[Import[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]];

(*turn off irrelevant error messages*)
Off[ClebschGordan::phy]; ParallelEvaluate[Off[ClebschGordan::phy]];
Off[NMinimize::nnum]; ParallelEvaluate[Off[NMinimize::nnum]];


(*Execute parallelized kerr solver*)
runStartTime = AbsoluteTime[];
Evaluators = Flatten@Table[
				ParallelSubmit[{mv, nv}, WorkerFunction[mv,nv]],
				{nv, parameters["nrange"] // First, parameters["nrange"] // Last}, 
				{mv, parameters["mrange"] // First, parameters["mrange"] // Last}
			];
WaitAll[Evaluators];



(*Final Message*)
runEndTime = AbsoluteTime[];
Print["Finished! Total Time: "<>ToString[N[Round[runEndTime-runStartTime]/60], InputForm]<>" Minutes"];	
			
(*Cleanup*)
CloseKernels[];	
