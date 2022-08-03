#!/usr/bin/env wolframscript
(* ::Package:: *)

$FKKSRoot = "/home/fell/Documents/Github/KerrDressedWithProca/"; (*directory on ITP cluster to git repo*)
(*$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";*)
LogFile = $FKKSRoot<>"Logs/LogFile.log";

ConfigurationParameters = Import[$FKKSRoot<>"WolframScript/Config.ini"];
LogRun = ToExpression[ConfigurationParameters["LogRun"]];
SaveToDisk = ToExpression[ConfigurationParameters["SaveToDisk"]];
NumberOfKernels = ConfigurationParameters["NumberofSlaveKernels"]//ToExpression;
OverwriteExisting = ConfigurationParameters["OverwriteExisting"]//ToExpression;

If[!TrueQ[(ConfigurationParameters["SolutionFilePath"]==Automatic)],
	$SolutionPath = ConfigurationParameters["SolutionFilePath"];,
	$SolutionPath = $FKKSRoot<>"Solutions/";
];

Print["----------------------------------"];
Print["FKKSRoot: "<>$FKKSRoot];
Print["Solution Path: "];Print["     "<>ToString@$SolutionPath];
Print["Log File Path: "];Print["     "<>ToString@LogFile];
Print["Log Run: "<>ToString[TrueQ[LogRun], InputForm]];
Print["Save to disk: "<>ToString[TrueQ[SaveToDisk], InputForm]];
Print["Overwrite Previous Solutions: "<>ToString[TrueQ[OverwriteExisting], InputForm]];
Print["----------------------------------"];

Import[FileNameJoin[{$FKKSRoot,"Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "ModeOvertoneWorker.wl"}]];

parameterlist = Import[$FKKSRoot<>"WolframScript/parameters.ini"]//ToExpression;
parameters = KeyMap[ToString[ToExpression[#]]&, parameterlist];
parameters["StartingX"] = parameters["\[Epsilon]"];


CreateFile[LogFile, OverwriteTarget->True];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];


Messenger = "\n\n-----------------------\nBeginning Parallelized Proca Solver\n";
Print[Messenger];

(*Setup parallelization*)
CloseKernels[];
NumberKernels = Min[NumberOfKernels, (Last[parameters["mrange"]] - First[parameters["mrange"]]+1)*(Last@parameters["nrange"] - First@parameters["nrange"] + 1)];
Print["Launching "<>ToString[NumberKernels,InputForm]<>" Kernels."];
LaunchKernels[NumberKernels];
Print["Launched "<>ToString[NumberKernels, InputForm]<>" Kernels. Distributing definitions and import packages."];
DistributeDefinitions[parameters, WorkerFunction, $SolutionPath, $FKKSRoot, LogFile, LogRun];
ParallelEvaluate[Import[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]];

(*Execute parallelized kerr solver*)
Print["Executing parallel routine..."];
runStartTime = AbsoluteTime[];
Evaluators = Flatten@Table[
				ParallelSubmit[{mv, nv}, WorkerFunction[mv,nv, parameters, LogFile, OverwritePrevious->OverwriteExisting]],
				{nv, parameters["nrange"] // First, parameters["nrange"] // Last}, 
				{mv, parameters["mrange"] // First, parameters["mrange"] // Last}
			];
WaitAll[Evaluators];



(*Final Message*)
runEndTime = AbsoluteTime[];
Print["Finished! Total Time: "<>ToString[N[Round[runEndTime-runStartTime]/60], InputForm]<>" Minutes"];	
			
(*Cleanup*)
CloseKernels[];	
