#!/usr/bin/env wolframscript
(* ::Package:: *)

Unprotect[Print];
Print[expr_] /; ! FreeQ[expr, _KernelObject] = Null;
Print[expr_] /; ! FreeQ[expr, s_String /; StringContainsQ[s, "KernelObject"]] = Null;
Print[expr_] /; ! FreeQ[expr, s_String /; StringContainsQ[s, "remote kernel"]] = Null;
Protect[Print];



$FKKSRoot = "/home/fell/Documents/Github/KerrDressedWithProca/"; (*directory on ITP cluster to git repo*)
(*$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";*)
LogFile = $FKKSRoot<>"Logs/EnergyConverterLogFile.log";


ConfigurationParameters = Import[$FKKSRoot<>"WolframScript/EnergyConverterConfig.ini"];
LogRun = ToExpression[ConfigurationParameters["LogRun"]];
SaveToDisk = ToExpression[ConfigurationParameters["SaveToDisk"]];
NumberOfKernels = ConfigurationParameters["NumberofSlaveKernels"]//ToExpression;


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
Print["----------------------------------"];


Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "xActSetup.wl"}]]];
Block[{Print, CellPrint, Return},Import[FileNameJoin[{$FKKSRoot, "Packages", "EnergyMomentum.wl"}]]];


CreateFile[LogFile, OverwriteTarget->True];
If[\[Not]DirectoryQ[$SolutionPath], CreateDirectory[$SolutionPath]];


Messenger = "\n\n-----------------------\nBeginning Parallelized Energy Converter\n";
Print[Messenger];


Print["Importing data"];
AllFileNames = FileNames[All, $SolutionPath];


AllData = Import/@AllFileNames;
Print["Number of solution sets imported: "<>ToString[Length@AllData, InputForm]];


Print["Launching "<>ToString[NumberOfKernels, InputForm]<>" kernels."];
CloseKernels[];
Off[PrintAsCharacter::argx];
LaunchKernels[NumberOfKernels];
ParallelEvaluate[Off[PrintAsCharacter::argx]];
If[$KernelCount==NumberOfKernels,
Print["Launched "<>ToString[NumberOfKernels, InputForm]<>" kernels."];,
Print["Failed to launch all requested kernels. Number of launched kernels: "<>ToString[$KernelCount, InputForm]];
];

Print["Performing Energy Density Conversion"];

ParallelDo[
str="Solution "<>ToString[i]<>" of "<>ToString[Length[AllFileNames]];
Print[str];
str>>>LogFile;

GenerateEnergyDensityCSVData[AllData[[i]]],
{i,1,Length[AllFileNames]}
]


Print["Energy Density Converter Complete."]
