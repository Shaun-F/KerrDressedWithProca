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
