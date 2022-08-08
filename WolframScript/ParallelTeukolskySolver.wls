#!/usr/bin/env wolframscript
(* ::Package:: *)

Unprotect[CellPrint];
CellPrint[Cell[s_, "Print", label_, ShowCellLabel->True]]:=Block[{},
CellPrint[Cell[s,"Print",ShowCellLabel->False]]];
Protect[CellPrint];


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


Import[FileNameJoin[{$FKKSRoot, "Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];
Get[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]
Import[FileNameJoin[{$FKKSRoot, "Packages", "EnergyMomentum.wl"}]]
Import[FileNameJoin[{$FKKSRoot, "Packages", "TeukolskySolver.wl"}]];


CreateFile[LogFile, OverwriteTarget->True];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];


Messenger = "\n\n-----------------------\nBeginning Parallelized Teukolsky Solver\n";
Print[Messenger];


Print["Importing data"];
AllFileNames = FileNames[All, $FKKSRoot<>"Solutions/"];

AllFileNames = AllFileNames[[1;;6]];


AllData = Import/@AllFileNames;
Print["Number of solution sets imported: "<>ToString[Length@AllData, InputForm]];


Print["Launching "<>ToString[NumberOfKernels, InputForm]<>" kernels."];
CloseKernels[];
LaunchKernels[NumberOfKernels];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
Off[PrintAsCharacter::argx];
If[$KernelCount==NumberOfKernels,
Print["Launched "<>ToString[NumberOfKernels, InputForm]<>" kernels."];,
Print["Failed to launch all requested kernels. Number of launched kernels: "<>ToString[$KernelCount, InputForm]];
];


Print["Data Imported. Running normalization procedure on "<>ToString[$KernelCount,InputForm]<>" kernels."];
DistributeDefinitions[FKKSNormalization, AllData, FKKSTotalEnergy, FKKSEnergyDensity,TeukolskySourceModal];
TaskList = Table[ParallelSubmit[{j}, FKKSNormalization[AllData[[j]], 1, Recalculate->OverwriteExisting]], {j,1,Length@AllData}]
Normalizations = WaitAll[TaskList];

If[SaveToDisk,
Print["Normalization complete. Appending data to solution sets and normalizing radial functions."];
Do[AppendToSolution[AllData[[j]]]["Normalization", Normalizations[[j]]], {i,1,Length@AllData}];
];
AllNormedData = Table[RenormalizeProcaSolution[AllData[[j]], Normalizations[[j]]["Normalization"]], {j,1,Length@AllData}];
Print["Black hole evolution code complete. Proceeding to Teukolsky solver."];


Print["Calculating mode decomposition of Teukolsky source term"];
TaskListModeDecomp = Table[ParallelSubmit[{OrbitalNumber, ModeNumber,j}, TeukolskySourceModal[AllNormedData[[j]], OrbitalNumber, ModeNumber]], {j,1,Length@AllNormedData}];
AllTlmw = WaitAll[TaskListModeDecomp];
If[SaveToDisk,
Do[AppendToSolution[AllData[[j]]]["Tlmw", AllTlmw[[j]]], {j,1,Length@AllData}];
];


Print["Mode decomposition calculation complete. Calculating asymptotic amplitudes"];
Block[{$MSTMaxN=15},
TaskListAsympAmps = Table[ParallelSubmit[{AllNormedData, AllTlmw, OrbitalNumber, ModeNumber,j},TeukolskyZInfinity[AllNormedData[[j]], AllTlmw[[j]],OrbitalNumber, ModeNumber]], {j,1,Length@AllNormedData}];
AllAsymptoticAmplitudes = WaitAll[TaskListAsympAmps];
];
If[SaveToDisk,
Print["Asymptotic amplitude calculation complete. Appending data to solution"];
Do[AppendToSolution[AllData[[j]]]["Zinf", AllAsymptoticAmplitudes[[j]]], {j,1,Length@AllData}];
];


Print["Calculating asymptotic energy flux"]
TaskListAsympFlux = Table[ParallelSubmit[{AllNormedData,AllTlmw,AllAsymptoticAmplitudes,i},EnergyFlux[AllNormedData[[i]],2,2, TeukolskyTlmw->AllTlmw[[i]], ZCoefficient->AllAsymptoticAmplitudes[[i]]]], {i,1,Length@AllNormedData}];
AllAsymptoticFluxes = WaitAll[TaskListAsympFlux];
If[SaveToDisk,
Print["Asymptotic flux calculation complete. Appending data to solution"];
Do[AppendToSolution[AllData[[j]]]["Einf", AllAsymptoticFluxes[[j]]], {j,1,Length@AllData}];
];


Print["Teukolsky solver complete."];
