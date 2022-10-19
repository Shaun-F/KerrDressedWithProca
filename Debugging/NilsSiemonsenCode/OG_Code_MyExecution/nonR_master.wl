(* ::Package:: *)

CPU=1;
LaunchKernels[CPU];
Needs["QNMnonRelsearch`","./nonR_search.wl"];
destinationPATH=ToString[$root<>"surrogate"];
Print["Save Directory: "<>ToString[destinationPATH]];

ninput=0;
winit=10^-9;

spin=0.99;

OH[a_]:=a/(2(1+Sqrt[1-a^2]))
musat[aa_]:=(mu/.Solve[0.91mu(1-mu^2/2)==Oh,mu][[3]]/.{Oh -> OH[a]}//Simplify)/.{a -> aa}
Y[mu_,aa_]:=(mu-0.05)/(musat[a]-0.05)//.{\[Mu] -> mu,a -> aa}

Do[
sol[a]= x/.Solve[Y[x,a]==1,x][[1]];
,{a,{spin}}];
mustop=Table[sol[a], {a, {spin}}]//Re;
mustop=0.51;(*{0.2,0.4,0.47,0.5};*)

Print["\[Mu] stop: "<>ToString@mustop];
Print["spin: "<>ToString@spin];

Table[
spinstring=NumberForm[spin*10^6//Round,6,DigitBlock->5,ExponentStep->6,NumberSeparator->""];
strmMode=OpenWrite[ToString[ToString[destinationPATH]<>"/m"<>ToString[1]<>ToString["n"]<>ToString[ninput]<>ToString["_a"]<>ToString[spinstring]<>ToString["_S"]<>ToString["m"]<>ToString[1]<>ToString["_prec"]<>ToString["_HPe.log"]]];
AppendTo[$Output,strmMode];
AppendTo[$Messages,strmMode];
br=If[spin<0.84,1,2];
QNMnonRel[1,ninput,spin,0.05,mustop,br,winit,-1,1, ModePrecision->i];
Close[strmMode];
ClearSystemCache[];
,{i,17,25}];

Print["Parallel Session: Done!"];
