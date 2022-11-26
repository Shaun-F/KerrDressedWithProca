#!/usr/bin/env wolframscript
(* ::Package:: *)

$root = "/home/fell/Documents/Github/KerrDressedWithProca/Debugging/NilsSiemonsenCode/OG_Code_MyExecution/";
SetDirectory[$root];
$FKKSRoot = $root<>"../../../";
Import[$FKKSRoot<>"Packages/HelperFunctions.wl"];
$SolutionPath = $FKKSRoot<>"Solutions/";
Block[{Print},
<<HelperFunctions.wl;
<<BHEvol.m;
<<VectorFieldNorm.m;
<<SpinWeightedSpheroidalHarmonics.m;
<<RenormalizedAngularMomentum.m;
<<MSTformalism.m;
<<WeylPsi4.m;
<<TeukolskyT.m;
<<TeuInter.wl;
]


Do[TeuCode[1,0,9/10, i,1,MySolution->False], 
{i,6/100,33/100,1/100}]
