(* ::Package:: *)

LaunchKernels[2];
Needs["PathLoader`","./PathLoader.wl"];
ParallelNeeds["PathLoader`"];

ParallelNeeds["VectorFieldNorm`"];(*,"/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/VectorFieldNorm.m"];*)
ParallelNeeds["BHEvol`"];(*,"/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/BHEvol.m"];*)
ParallelNeeds["TeukolskyT`"];(*,"/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/TeukolskyT.m"];*)
ParallelNeeds["MSTformalism`"];(*,"/home/nils/uni/projects/superrad/mathematica_scripts/Teukolsky_code/MSTformalism.m"];*)
ParallelNeeds["WeylPsi4`"];(*,"./WeylPsi4.m"];*)
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];(*,"./SpinWeightedSpheroidalHarmonics.m"];*)

DistributeDefinitions[TeuCode];

m=1;
n=0;
spin=0.99;
Mini=1;
masses={0.1,0.15};

Table[
TeuCode[m,n,spin,masses[[i]],Mini];
,{i,1,2}];
