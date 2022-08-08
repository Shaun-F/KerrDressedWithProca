#!/usr/bin/env wolframscript
(* ::Package:: *)

Unprotect[CellPrint];
CellPrint[Cell[s_, "Print", label_, ShowCellLabel->True]]:=Block[{},
CellPrint[Cell[s,"Print",ShowCellLabel->False]]];
Protect[CellPrint];


$FKKSRoot = "/home/fell/Documents/Github/KerrDressedWithProca/"; (*directory on ITP cluster to git repo*)
(*$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";*)
LogFile = $FKKSRoot<>"Logs/TeukolskyLogFile.log";


ConfigurationParameters = Import[$FKKSRoot<>"WolframScript/TeukolskyConfig.ini"];
LogRun = ToExpression[ConfigurationParameters["LogRun"]];
SaveToDisk = ToExpression[ConfigurationParameters["SaveToDisk"]];
NumberOfKernels = ConfigurationParameters["NumberofSlaveKernels"]//ToExpression;
OverwriteExisting = ConfigurationParameters["OverwriteExisting"]//ToExpression;


Print["----------------------------------"];
Print["FKKSRoot: "<>$FKKSRoot];
Print["Solution Path: "];Print["     "<>ToString@$SolutionPath];
Print["Log File Path: "];Print["     "<>ToString@LogFile];
Print["Log Run: "<>ToString[TrueQ[LogRun], InputForm]];
Print["----------------------------------"];


Import[FileNameJoin[{$FKKSRoot,"Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "TeukolskySolver.wl"}]];


CreateFile[LogFile, OverwriteTarget->True];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];


Messenger = "\n\n-----------------------\nBeginning Parallelized Teukolsky Solver\n";
Print[Messenger];


Print["Importing data"];
AllFileNames = FileNames[All, $FKKSRoot<>"Solutions/";]
AllData = Import/@AllFileNames;
Print["Data Imported. Running normalization procedure on "<>ToString[NumberOfKernels,InputForm]<>" kernels."];

DistributeDefinitions[FKKSNormalization, AllData, FKKSTotalEnergy, FKKSEnergyDensity];
TaskList = Table[ParallelSubmit[{j}, FKKSNormalization[AllData[[j]], 1, Recalculate->True]], {j,1,Length@AllData}]
Normalizations = WaitAll[TaskList];

Print["Normalization complete. Appending data to solution sets and normalizing radial functions."];
Do[AppendToSolution[AllData[[j]]]["Normalization", Normalizations[[j]]], {i,1,Length@AllData}];
AllNormedData = Table[RenormalizeProcaSolution[AllData[[j]], Normalizations[[j]]["Normalization"]], {j,1,Length@AllData}];

Print["Black hole evolution code complete. Proceeding to Teukolsky solver."];









