(* ::Package:: *)

CPU=1;
LaunchKernels[CPU];
Needs["QNMovertoneSearch`","/home/nsiemonsen/proca_superrad/precision_mode_library/overT_search_prec.wl"];

DistributeDefinitions[QNMoverSearchPrec];

spin={Range[0.6,0.995,(0.995-0.6)/(40-1)][[39]]};
(*spin={0.4};*)
minput=2;
ninput=0;

OH[a_] := a/(2 (1 + Sqrt[1 - a^2]))
musat = ((mu /. NSolve[0.91 mu (1 - mu^2/(2(minput)^2)) == minput Oh, mu][[3]]) /. {Oh->OH[aa]} // Simplify);
mustop = Re[Table[musat/.aa->spin[[i]],{i,1,CPU}]];
Print[mustop];
Print[spin];

(*spin={0.995,0.995,0.995,0.995,0.55,0.99,0.99,0.99,0.95,0.95,0.95,0.95,0.9,0.9,0.9,0.9,0.85,0.85,0.85,0.85,0.8,0.8,0.8,0.8,0.75,0.75,0.75,0.75,0.7,0.7,0.7,0.7,0.65,0.65,0.65,0.65,0.6,0.6,0.6,0.6};
ninput={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};*)

aa1=spin//Length;
aa2=min//Length;
aa3=nin//Length;
aa4=\[Mu]sta//Length;
aa5=\[Mu]sto//Length;
aa6=winit//Length;
Print[aa1];
Print[aa2];
Print[aa3];
Print[aa4];
Print[aa5];
Print[aa6];
Print["---------------------------"];


ParallelTable[
Pause[10*i];
spinstring=NumberForm[spin[[i]]*10^6//Round,6,DigitBlock->5,ExponentStep->6,NumberSeparator->""];
strmMode=OpenWrite[ToString["/home/nsiemonsen/proca_superrad/m2_modes/m"<>ToString[minput]<>ToString["n"]<>ToString[ninput]<>ToString["_a"]<>ToString[spinstring]<>ToString["_S"]<>ToString["m"]<>ToString[1]<>ToString["_AccGoal5.log"]]];
AppendTo[$Output,strmMode];
AppendTo[$Messages,strmMode];
QNMoverSearchPrec[minput,ninput,spin[[i]],0.25,mustop[[i]],2,10^-10,-1,1,0.005];
Close[strmMode];
ClearSystemCache[];
,{i,1,CPU}]

Print["Parallel Session: Done!"];
