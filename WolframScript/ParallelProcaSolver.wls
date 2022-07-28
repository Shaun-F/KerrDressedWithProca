#!/usr/bin/env wolframscript
(* ::Package:: *)

LogRun=True;
SaveToDisk=True;


$FKKSRoot = "/home/shaunf/Documents/Computer/Code/projects/Massive_Vector_Field_Dynamical_Friction/ProcaAroundKerr/FKKSsolver/Mathematica/";
$SolutionPath = $FKKSRoot<>"Solutions/";
LogFile = $FKKSRoot<>"Logs/LogFile";

Import[FileNameJoin[{$FKKSRoot,"Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];

parameterlist = Import[$FKKSRoot<>"WolframScript/parameters.ini"]//ToExpression;
parameters = KeyMap[ToString[ToExpression[#]]&, parameterlist];



If[FileExistsQ[LogFile], DeleteFile[LogFile]];
CreateFile[LogFile];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];
Messenger = "Beginning Parallelized Proca Solver\n";
Print[Messenger];

WorkerFunction[mval_,nval_]:=
	Block[{HistoryLength=0},
		With[{k=mval, j=nval},
			muIterationCounter=1;
			Do[
				(*Prepare parameter values*)
				parameters["m"]=k;
				parameters["l"]=parameters["m"]+parameters["s"];
				parameters["n"]=j;
				
				(*Setup iteration*)
				parameters["\[Mu]Nv"]=i*parameters["m"]; (*account for m-dependence of superradiant threshold*)
				parameters["EndingX"]=calcEndingX[parameters]; (*Determine cutoff radius for radial equation integrator*)
				
				(*Output current iteration to terminal*)
				Print["Kernel "<>ToString[$KernelID]<>" working on parameter point (\[Mu], m, n) = ("<>ToString[parameters["\[Mu]Nv"], InputForm]<>", "<>ToString[parameters["m"], InputForm]<>" ,"<>ToString[parameters["n"], InputForm]<>")"];
				
				Logger = "m Iteration Count: "<>ToString[mIterationCounter, InputForm]<>"\n";
				Logger = Logger<>"m value: "<>ToString[parameters["m"], InputForm]<>"\n";
				Logger = Logger<>"n value: "<>ToString[parameters["n"], InputForm]<>"\n";
				Logger = Logger<>"mu value: "<>ToString[parameters["\[Mu]Nv"], InputForm]<>"\n";
				
				(*Minimization using native Nelder Mead algorithm*)
				minimizetimestart = AbsoluteTime[];
				
				(*Determine next guess for frequency and angular eigenvalue*)
				nuFitFunction[muIterationCounter] = Null;
				If[TrueQ[muIterationCounter >= 6],
					(*If we're on the 6th iteration, use previous values to approximate current frequency value using fit functions*)
					omegaValHistory = Table[omegaHolder[iter], {iter, 1, muIterationCounter-1}];
					nuValHistory = Table[nuHolder[iter], {iter, 1, muIterationCounter-1}];
					muValHistory = Table[muHolder[iter], {iter, 1, muIterationCounter-1}];
					omegaFit[muIterationCounter] = Fit[Thread[{muValHistory,omegaValHistory}], {1, zeta, zeta^2, zeta^3, zeta^4}, zeta, WorkingPrecision->parameters["precision"]];
					nuFit[muIterationCounter] = Fit[Thread[{muValHistory,nuValHistory}], {1, zeta, zeta^2, zeta^3, zeta^4}, zeta, WorkingPrecision->parameters["precision"]];
					CurrentomegaGuess = omegaFit[muIterationCounter]/.zeta->parameters["\[Mu]Nv"];
					CurrentnuGuess = nuFit[muIterationCounter]/.zeta->parameters["\[Mu]Nv"];
					omegaGuess = CurrentomegaGuess;
					nuGuess = getNuValue[omegaGuess,parameters, CurrentnuGuess];
					nuFitFunction[muIterationCounter] = Function[{zeta}, Evaluate@nuFit[muIterationCounter]]; (*To be used later for the radial minimization*),
					
					(*Use non-relativistic limit for initial guess of omega in first 5 iterations*)
					omegaGuess = N[omegaNRNonRel[parameters] + I*omegaNINonRel[parameters]];
					nuGuess = getNuValue[omegaGuess, parameters, nuNNonRel[omegaGuess, parameters]];
				];
				
				omegaBoundarySize = (1/2)*(omegaNRNonRel[ReplacePart[parameters, "n"->parameters["n"]+1]] - Re[omegaGuess]); (*Use next overtone frequency to construct boundary for omega search*)
				omegaBoundary = {omegaGuess - omegaBoundarySize, omegaGuess + omegaBoundarySize};
				
				(*Perform Minimization*)
				MinimizationResults = RadialMinimize[omegaGuess, nuGuess, omegaBoundary, parameters, False, nuFitFunction[muIterationCounter]];
				Logger = Logger<>"Minimization Results: "<>ToString[MinimizationResults//N, InputForm];
				minimizetimestop = AbsoluteTime[];
				Logger = Logger<>"\nMinimization Time: "<>ToString[N[minimizetimestop-minimizetimestart], InputForm];
				
				
				(*Solve for Kernel of angular matrix*)
				angulartimestart = AbsoluteTime[];
				Sfunctemp = Sum[C[tempiter]*SphericalHarmonicY[Abs[parameters["m"]] + 2*tempiter + parameters["\[Eta]"],parameters["m"],\[Theta], \[Phi]]*Exp[-I*parameters["m"]*\[Phi]], {tempiter,0, parameters["KMax"]}];
				coeffs = SolveAngularSystem[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"],parameters];
				Sfunc = Evaluate[Sfunctemp/.coeffs];
				angulartimestop = AbsoluteTime[];
				Logger = Logger<>"\nAngular Solver Time: "<>ToString[N[angulartimestop-angulartimestart], InputForm];
				
				
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
				If[parameters["\[Mu]Nv"]^2 <= Abs[MinimizationResults["\[Omega]"]]^2,
					Print["Kernel "<>ToString[$KernelID,InputForm]<>" says: Error! Frequency larger than field mass! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];                
					Break[];
				];
				If[TheSolution["R"][parameters["EndingX"]]//Log10 > 1,
					Print["Kernel "<>ToString[$KernelID]<>" says: Error! Minimization failed! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];
					Break[];
				];
						
				(*Print solution to terminal*)
				Print["Result for Kernel "<>ToString[$KernelID, InputForm]<>": \n\t \[Omega] = "<>ToString[Evaluate@MinimizationResults["\[Omega]"], InputForm]];
				
				(*Prepare solution data and write to disk*)
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s","\[Chi]", "KMax", "branch"}]];
				CacheData = <| "Parameters" -> parameters, "Solution" -> TheSolution|>;
				Logger = Logger<>"\nCaching results to: "<>assocToString[printData]<>".mx\n\n";
				If[SaveToDisk,
					Export[$SolutionPath<>"RunData_"<>assocToString[printData]<>".mx",CacheData];
				];
				If[LogRun, Logger >>> LogFile];
				
				(*Prepare next loop*)
				muHolder[muIterationCounter] = parameters["\[Mu]Nv"]//SetPrecision[#, parameters["precision"]]&;
				omegaHolder[muIterationCounter] = MinimizationResults["\[Omega]"]//SetPrecision[#, parameters["precision"]]&;
				nuHolder[muIterationCounter] =  MinimizationResults["\[Nu]"]//SetPrecision[#, parameters["precision"]]&;
				
				muIterationCounter++;
				
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
NumberKernels = Min[$ProcessorCount, (Last[parameters["mrange"]] - First[parameters["mrange"]]+1)*(Last@parameters["nrange"] - First@parameters["nrange"] + 1)];
LaunchKernels[NumberKernels];
DistributeDefinitions[parameters, WorkerFunction, $SolutionPath, $FKKSRoot];
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
				
				
				
				
					
					


































































				
