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
Print["FKKSRoot: "<>$FKKSRoot];
Print["Solution Path: "];Print["     "<>ToString@$SolutionPath];
Print["Log File Path: "];Print["     "<>ToString@LogFile];
Print["Solution Path: "];Print["     "<>ToString@$SolutionPath];
Print["Log Run: "<>ToString[TrueQ[LogRun], InputForm]];
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
	Block[{Normalization,NormedData, Tlmw, AsymptoticAmplitude,AsymptoticFlux},
	With[{SingleSolution = AllData[[i]]},
		Print["\n\nProgress: "<>ToString[i]<>" out of "<>ToString[Length[AllData]]];
		
		(*----------Black hole evolution----------*)
		Normalization=FKKSNormalization[SingleSolution, 1, Recalculate->OverwriteExisting];
		If[SaveToDisk,
			Print["Normalization complete. Appending data to solution set and normalizing radial functions."];
			AppendToSolution[SingleSolution]["Normalization", Normalization];
		];
		NormedData = RenormalizeProcaSolution[SingleSolution, Normalization["Normalization"]];
		Print["Black hole evolution code complete. Proceeding to Teukolsky solver."];


		(*----------mode decomposition of Teukolsky source function----------*)
		Print["Calculating mode decomposition of Teukolsky source term"];
		Tlmw = TeukolskySourceModal[NormedData, OrbitalNumber, ModeNumber];
		If[SaveToDisk,
			AppendToSolution[SingleSolution]["Tlmw", Tlmw];
		];


		(*----------Asymptotic amplitude----------*)
		Print["Mode decomposition calculation complete. Calculating asymptotic amplitudes"];
		Block[{$MSTMaxN=15},
			AsymptoticAmplitude = TeukolskyZInfinity[NormedData, Tlmw, OrbitalNumber, ModeNumber];
		];
		If[SaveToDisk,
			Print["Asymptotic amplitude calculation complete. Appending data to solution"];
			AppendToSolution[SingleSolution]["Zinf", AsymptoticAmplitude];
		];


		(*----------Asymptotic energy flux----------*)
		Print["Calculating asymptotic energy flux"];
		AsymptoticFlux = EnergyFlux[NormedData,OrbitalNumber,ModeNumber, TeukolskyTlmw->Tlmw, ZCoefficient->AsymptoticAmplitude];
		If[SaveToDisk,
			Print["Asymptotic flux calculation complete. Appending data to solution"];
			AppendToSolution[SingleSolution]["Einf", AsymptoticFlux];
		];
	];
	];,
	{i,1,Length@AllData}
];


Print["Teukolsky solver complete."];
