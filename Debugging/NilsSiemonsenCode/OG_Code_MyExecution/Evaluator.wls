#!/usr/bin/env wolframscript
(* ::Package:: *)

$root = "/home/shaunf/Documents/Computer/Code/projects/Massive_Vector_Field_Dynamical_Friction/ProcaAroundKerr/FKKSsolver/KerrDressedWithProca/Debugging/NilsSiemonsenCode/OG_Code_MyExecution/";
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


Do[TeuCode[1,0,9/10, i,1,MySolution->True], 
{i,6/100,33/100,1/100}]
