(* ::Package:: *)

ModeOvertoneWorker;

(*Import Helper functions*)
If["HelperFunctions"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]]
(*Import xAct Setup*)
If["xActSetup"\[NotElement]Names["Global`*"], Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]]


Options[WorkerFunction]={OverwritePrevious->False};
WorkerFunction[modeval_,overtoneval_, InitialParameters_, LogFile_, OptionsPattern[]]:=
	Block[{HistoryLength=0, 
			parameters=InitialParameters,
			ResolutionReductionPoint, muIterationCounter,
			Logger, printData,SolutionFileName,
			minimizetimestart,minimizetimestop,
			nuFitFunction,omegaValHistory,nuValHistory,muValHistory,omegaHolder,nuHolder,muHolder,omegaFit,nuFit,CurrentomegaGuess,CurrentnuGuess,omegaGuess,nuGuess,
			omegaBoundarySize,omegaBoundary,MinimizationResults,
			angulartimestart, angulartimestop, angularcoeffs, Sfunc,Sfunctemp,coeffs,
			omegafitOutVariable,nufitOutVariable,
			RadialSolution,RadialPlot,AngularPlot, TheSolution,CacheData,CacheResults
			},
		With[{k=modeval, j=overtoneval},
			(*Turn off irrelevant error messages*)
			Off[ClebschGordan::phy];
			Off[NMinimize::nnum];
			Off[InterpolatingFunction::dmval];
			muIterationCounter=1;
			ResolutionReductionPoint=Null;
			Logger = {"Beginning search routine:   time:" <> DateString[] <>" Parameter values: (m,n,\[Chi]) = (" <> ToString[k, InputForm] <>"," <> ToString[j, InputForm] <> "," <> ToString[parameters["\[Chi]"], InputForm] <> ")"};
			While[True,
				(*-----------Prepare parameter values-----------*)
				parameters["m"]=k;
				parameters["l"]=parameters["m"]+parameters["s"];
				parameters["n"]=j;
				
				(*-----------Determine new mass value-----------*)	
				parameters["\[Mu]Nv"]=(parameters["\[Mu]range"][[1]] + muIterationCounter*parameters["\[Delta]\[Mu]"])*parameters["m"]; (*account for m-dependence of superradiant threshold*)
				(*If the Imaginary part of the frequency is starting to decrease, we half the step size and remove the mode postfactor to more accurately probe the fall off*)
				If[muIterationCounter>3,
					If[Im[omegaHolder[muIterationCounter-2]]>Im[omegaHolder[muIterationCounter-1]],
						Print["Im(omega) decreasing. Reducing mass resolution..."];
						(*Determine the muIterationCounter at which the resolution should be reduced*)
						If[!NumberQ[ResolutionReductionPoint],
							ResolutionReductionPoint=muIterationCounter-1;
						];
						(*Reduce resolution and continue from last value*)
						parameters["\[Mu]Nv"]=First[parameters["\[Mu]range"]] + (muIterationCounter+ResolutionReductionPoint)*parameters["\[Delta]\[Mu]"]/2;
					];
					
					If[Im[omegaHolder[muIterationCounter-3]]>Im[omegaHolder[muIterationCounter-2]]<Im[omegaHolder[muIterationCounter]-1],
						(*If the instability rate decreases then increases, we probably jumped across the superradiant threshold gap \[Omega]=m*\[CapitalOmega]. Break the loop.*)
						Print["Im(omega) decreased then increased! Possibly jumped superradiant threshold gap. Breaking..."];
						AppendTo[Logger, "Im(omega) decreased then increased! Possibly jumped superradiant threshold gap. Breaking..."];
						PutAppend[Sequence@@Logger, LogFile];
						Break[];	
					];
				];
				(*-----------Safety check on mass value-----------*)
				If[!NumberQ[parameters["\[Mu]Nv"]],
					Print["Error! Mass value is not a number! \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]];
					AppendTo[Logger, "Error! Mass value is not a number! \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]];
					PutAppend[Sequence@@Logger, LogFile];
					Break[];
				];
				parameters["EndingX"]=calcEndingX[parameters]; (*Determine cutoff radius for radial equation integrator*)
				
				(*-----------Construct name for file to save data in-----------*)
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s", "\[Chi]", "KMax", "branch"}]];
				SolutionFileName = $SolutionPath<>"RunData_"<>assocToString[printData]<>".mx";
				
				(*-----------If OverwritePrevious option is false and if the solution file already exists, skip current loop iteration-----------*)
				If[! OptionValue[OverwritePrevious],
					If[FileExistsQ[SolutionFileName],
						Continue[]
					];
				];
				
				(*-----------Output current iteration to terminal-----------*)
				Print["Kernel "<>ToString[$KernelID]<>" working on parameter point (\[Mu], m, n) = ("<>ToString[parameters["\[Mu]Nv"], InputForm]<>", "<>ToString[parameters["m"], InputForm]<>" ,"<>ToString[parameters["n"], InputForm]<>")"];
				
				AppendTo[Logger,"m value: "<>ToString[parameters["m"], InputForm]];
				AppendTo[Logger,"n value: "<>ToString[parameters["n"], InputForm]];
				AppendTo[Logger,"mu value: "<>ToString[parameters["\[Mu]Nv"], InputForm]];
				
				(*-----------Minimization using native Nelder Mead algorithm-----------*)
				minimizetimestart = AbsoluteTime[];
				
				(*Determine next guess for frequency and angular eigenvalue*)
				nuFitFunction[muIterationCounter] = Null;
				If[TrueQ[muIterationCounter >= 6],
					(*If we're on the 6th iteration, use previous values to approximate current frequency value using fit functions*)
					omegaValHistory = Table[omegaHolder[iter], {iter, 1, muIterationCounter-1}];
					nuValHistory = Table[nuHolder[iter], {iter, 1, muIterationCounter-1}];
					muValHistory = Table[muHolder[iter], {iter, 1, muIterationCounter-1}];
					omegaFit[muIterationCounter] = Interpolation[Thread[{muValHistory, omegaValHistory}], InterpolationOrder->parameters["FrequencyInterpolatingOrder"]];
					nuFit[muIterationCounter] = Interpolation[Thread[{muValHistory, nuValHistory}], InterpolationOrder->parameters["FrequencyInterpolatingOrder"]];
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
				
				(*Perform Minimization-*)
				MinimizationResults = RadialMinimize[omegaGuess, nuGuess, omegaBoundary, parameters, False, nuFitFunction[muIterationCounter]];
				AppendTo[Logger,"Minimization Results: "<>ToString[MinimizationResults//N, InputForm]];
				minimizetimestop = AbsoluteTime[];
				AppendTo[Logger,"Minimization Time: "<>ToString[N[minimizetimestop-minimizetimestart], InputForm]];
				
				
				(*-----------Solve for Kernel of angular matrix-----------*)
				angulartimestart = AbsoluteTime[];
				Sfunctemp = Sum[C[tempiter]*SphericalHarmonicY[Abs[parameters["m"]] + 2*tempiter + parameters["\[Eta]"],parameters["m"],\[Theta], \[Phi]]*Exp[-I*parameters["m"]*\[Phi]], {tempiter,0, parameters["KMax"]}];
				coeffs = SolveAngularSystem[MinimizationResults["\[Omega]"], MinimizationResults["\[Nu]"],parameters];
				Sfunc = Function[{\[Theta]},Evaluate[Sfunctemp/.coeffs]];
				angulartimestop = AbsoluteTime[];
				AppendTo[Logger,"Angular Solver Time: "<>ToString[N[angulartimestop-angulartimestart], InputForm]];
				
				
				(*-----------Format solution-----------*)
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
						
				(*-----------Safety checks on solution-----------*)
				If[TrueQ[parameters["\[Mu]Nv"]^2 <= Abs[MinimizationResults["\[Omega]"]]^2],
					Print["Kernel "<>ToString[$KernelID,InputForm]<>" says: Error! Frequency larger than field mass! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];                
					Break[];
				];
				If[TrueQ[TheSolution["R"][parameters["EndingX"]]//Log10 > 1],
					Print["Kernel "<>ToString[$KernelID]<>" says: Error! Minimization failed! Breaking mu iteration ... (Parameter Point : m = "<>ToString[parameters["m"], InputForm]<>" n = "<>ToString[parameters["n"]]<>" \[Mu] = "<>ToString[parameters["\[Mu]Nv"], InputForm]<>") "];
					Break[];
				];
						
				(*-----------Print solution to terminal-----------*)
				Print["Result for Kernel "<>ToString[$KernelID, InputForm]<>": \n\t \[Omega] = "<>ToString[Evaluate@MinimizationResults["\[Omega]"], InputForm]];
				
				(*-----------Prepare solution data and write to disk-----------*)
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
				
				(*-----------Prepare next loop-----------*)
				muHolder[muIterationCounter] = parameters["\[Mu]Nv"]//SetPrecision[#, parameters["precision"]]&;
				omegaHolder[muIterationCounter] = MinimizationResults["\[Omega]"]//SetPrecision[#, parameters["precision"]]&;
				nuHolder[muIterationCounter] =  MinimizationResults["\[Nu]"]//SetPrecision[#, parameters["precision"]]&;
				
				muIterationCounter++;
				Continue[];
				
				(*-----------Fallback if superradiatn condition fails-----------*);
				Label[SuperradiantConditionFailure];
				AppendTo[Logger,"Superradiant condition failed for parameters (\[Mu],m,n,\[Chi],s)=(" <>ToString[parameters["\[Mu]Nv"], InputForm] <> "," <>ToString[parameters["m"], InputForm] <> "," <>ToString[parameters["n"], InputForm] <> "," <>ToString[parameters["\[Chi]"], InputForm] <> "," <>ToString[parameters["s"], InputForm] <> "). Breaking mass loop..."];
				If[LogRun, PutAppend[Sequence@@Logger,LogFile]];	
				printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s","\[Chi]", "KMax", "branch"}]];
				CacheData = <|"Parameters" -> parameters, "Solution" -> "Superradiant condition failed"|>;
				If[SaveToDisk,Export[$SolutionPath <> "RunData_" <> assocToString[printData]<>".mx", CacheData];];
				Break[];
			]; (*mu iteration loop end*)
			
			(*Clear symbols used in mu iteration to prevent memory leaks*)
			ClearAll[muHolder, omegaHolder, nuHolder, TheSolution, AngularPlot, RadialPlot, CacheData, Logger, omegaFit, nuFit, nuFitFunction];
			ClearSystemCache[];
		]; (*With statement end*)
	];(*Block statement end*)
	(*End of worker function definition*)

