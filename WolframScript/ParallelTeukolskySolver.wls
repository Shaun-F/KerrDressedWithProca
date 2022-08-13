#!/usr/bin/env wolframscript
(* ::Package:: *)

Unprotect[Print];
Print[expr_] /; ! FreeQ[expr, _KernelObject] = Null;
Print[expr_] /; ! FreeQ[expr, s_String /; StringContainsQ[s, "KernelObject"]] = Null;
Print[expr_] /; ! FreeQ[expr, s_String /; StringContainsQ[s, "remote kernel"]] = Null;
Protect[Print];


(*$FKKSRoot = "/home/fell/Documents/Github/KerrDressedWithProca/";*) (*directory on ITP cluster to git repo*)
$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";
LogFile = $FKKSRoot<>"Logs/TeukolskyLogFile.log";


ConfigurationParameters = Import[$FKKSRoot<>"WolframScript/TeukolskyConfig.ini"];
LogRun = ToExpression[ConfigurationParameters["LogRun"]];
SaveToDisk = ToExpression[ConfigurationParameters["SaveToDisk"]];
NumberOfKernels = ConfigurationParameters["NumberofSlaveKernels"]//ToExpression;
OverwriteExisting = ConfigurationParameters["OverwriteExisting"]//ToExpression;
OrbitalNumber = ConfigurationParameters["OrbitalNumber"]//ToExpression;
ModeNumber = ConfigurationParameters["ModeNumber"]//ToExpression;

If[!TrueQ[(ConfigurationParameters["SolutionFilePath"]==Automatic)],
	$SolutionPath = ConfigurationParameters["SolutionFilePath"];,
	$SolutionPath = $FKKSRoot<>"Solutions/";
];


Print["----------------------------------"];
Print["FKKSRoot: "];Print["    "<>$FKKSRoot];
Print["Solution Path: "];Print["     "<>ToString@$SolutionPath];
Print["Log File Path: "];Print["     "<>ToString@LogFile];
Print["Log Run: "<>ToString[TrueQ[LogRun], InputForm]];
Print["Save to Disk: "<>ToString[TrueQ[SaveToDisk], InputForm]];
Print["Number of Kernels: "<>ToString[NumberOfKernels]];
Print["Orbital Number: "<>ToString[OrbitalNumber]];
Print["Mode Number: "<>ToString[ModeNumber]];
Print["----------------------------------"];


Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "EnergyMomentum.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "MSTSolver.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "TeukolskySolver.wl"}]]];


CreateFile[LogFile, OverwriteTarget->True];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];


Messenger = "\n\n-----------------------\nBeginning Parallelized Teukolsky Solver\n";
Print[Messenger];


Print["Importing data"];
AllFileNames = FileNames[All, $SolutionPath];

AllData = Import/@AllFileNames;
Print["Number of solution sets imported: "<>ToString[Length@AllData, InputForm]];


Print["Launching "<>ToString[NumberOfKernels, InputForm]<>" kernels."];
CloseKernels[];
Off[PrintAsCharacter::argx];
Off[NIntegrate::slwcon];
LaunchKernels[NumberOfKernels];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
ParallelEvaluate[Off[NIntegrate::slwcon]];
If[$KernelCount==NumberOfKernels,
Print["Launched "<>ToString[NumberOfKernels, InputForm]<>" kernels."];,
Print["Failed to launch all requested kernels. Number of launched kernels: "<>ToString[$KernelCount, InputForm]];
];


Print["Data Imported. Running normalization procedure on "<>ToString[$KernelCount,InputForm]<>" kernels."];
Do[
	Block[{Normalization,NormedData, Tlmw, AsymptoticAmplitude,AsymptoticFlux, Logger},
	Logger = {};
	With[{SingleSolution = AllData[[i]]},
		IterationStartTime = AbsoluteTime[];
		(*Skip solution if already calculated*)
		If[KeyExistsQ[SingleSolution, "Derived"]&&KeyExistsQ[SingleSolution["Derived"], "Einf"],
			Print["Solution already exists, skipping iteration."];
			AppendTo[Logger, ""]; AppendTo[Logger, ""];
			AppendTo[Logger, "Solution already exists, skipping iteration."];
			AppendTo[Logger, ""]; AppendTo[Logger, ""];
			PutAppend[Logger, LogFile];
			Continue[]
		];
		Print["\n\nProgress: "<>ToString[i]<>" out of "<>ToString[Length[AllData]]];
		AppendTo[Logger, "\n\nProgress: "<>ToString[i]<>" out of "<>ToString[Length[AllData]]];
		(*----------Black hole evolution----------*)
		Print["Beginning normalization routine of solution"];
		AppendTo[Logger, "Beginning normalization routine of solution"];
		Normalization=FKKSNormalization[SingleSolution, 1, Recalculate->OverwriteExisting];

		NormedData = RenormalizeProcaSolution[SingleSolution, Normalization["Normalization"]];
		Print["Black hole evolution code complete. Proceeding to Teukolsky solver."];
		AppendTo[Logger, "Black hole evolution code complete. Proceeding to Teukolsky solver."];

		(*----------mode decomposition of Teukolsky source function----------*)
		Print["Calculating mode decomposition of Teukolsky source term"];
		AppendTo[Logger, "Calculating mode decomposition of Teukolsky source term"];
		Tlmw = TeukolskySourceModal[NormedData, OrbitalNumber, ModeNumber];


		(*----------Asymptotic amplitude----------*)
		Print["Mode decomposition calculation complete. Calculating asymptotic amplitudes"];
		AppendTo[Logger, "Mode decomposition calculation complete. Calculating asymptotic amplitudes"];
		Block[{$MSTMaxN=15},
			AsymptoticAmplitude = TeukolskyZInfinity[NormedData, Tlmw, OrbitalNumber, ModeNumber];
		];


		(*----------Asymptotic energy flux----------*)
		Print["Calculating asymptotic energy flux"];
		AppendTo[Logger, "Calculating asymptotic energy flux"];
		AsymptoticFlux = EnergyFlux[NormedData,OrbitalNumber,ModeNumber, TeukolskyTlmw->Tlmw, ZCoefficient->AsymptoticAmplitude];
		If[SaveToDisk,
			Print["Asymptotic flux calculation complete. Appending data to solution"];
			AppendToSolution[SingleSolution, Prelabel->"FullSolutionSet"][{{"Einf", AsymptoticFlux}, {"Tlmw",Tlmw}, {"Normalization", Normalization}, {"Zinf", AsymptoticAmplitude}}];
		];
		IterationStopTime=AbsoluteTime[];
		Print["Memory in use by Wolfram Kernel: "<>ToString[MemoryInUse[], InputForm]];
		Print["Time to complete iteration: "<>ToString[IterationStopTime-IterationStartTime, InputForm]];
		AppendTo[Logger, "Energy Flux: "<>ToString[AsymptoticFlux, InputForm]];
		AppendTo[Logger, "Memory in use by Wolfram Kernel: "<>ToString[MemoryInUse[], InputForm]];
		AppendTo[Logger, "Time to complete iteration: "<>ToString[IterationStopTime-IterationStartTime, InputForm]];
		PutAppend[Sequence@@Logger, LogFile];
	];(*With*)
	];(*Block*)
	ClearSystemCache[],
	{i,1,Length@AllData}
];(*do*)


Print["Teukolsky solver complete."];
