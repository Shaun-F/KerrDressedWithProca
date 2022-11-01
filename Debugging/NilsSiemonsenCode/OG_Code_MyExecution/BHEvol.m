(* ::Package:: *)

BeginPackage["BHEvol`"]

KerrProcaEvolution::usage = "KerrProcaEvolution[m-mode,n-overtone,spin,procamass,Min] in terms of old BH mass units."
Mfspinf::usage = "Mfspinf[[1]] = Final BH mass in M0=1 units; Mfspinf[[2]] = Final dimensionless BH spin parameter."

Begin["PrivateBHEvol`"]

KerrProcaEvolution[min_,nin_,spinin_,procamass_,Mini_]:=Module[{mode1=min,overtone1=nin,spin1=spin,Procamass1=Procamass},
prec=SetPrecision[#,TeuInterEnv`$precValue]&;

(*Import the chosen Proca field mode*)
variableconversion={QNMcode`r\[Nu]->r\[Nu],QNMcode`i\[Nu]->i\[Nu],QNMcode`R->R,QNMcode`\[Theta]->\[Theta],QNMmode`r\[Nu]->r\[Nu],QNMmode`i\[Nu]->i\[Nu],QNMmode`R->R,QNMmode`\[Theta]->\[Theta]};
spinstring=NumberForm[spinin*10^6//Round,6,DigitBlock->5,ExponentStep->6,NumberSeparator->""];
ProcaSolution=DeleteCases[Import[Global`$root<>"precision_mode_library/m"<>ToString[min]<>ToString["n"]<>ToString[nin]<>ToString["_a"]<>ToString[spinstring]<>If[min==1,ToString["_Sm1_prec.mx"],ToString["_Sm1_prec_overT.mx"]]]//.variableconversion,0,3];
ProcaSolutionLength=Length[ProcaSolution[[3,All]]];
\[Omega]rof\[Mu]=Interpolation[Table[{ProcaSolution[[3,i]],ProcaSolution[[5,i]]//Re},{i,1,ProcaSolutionLength}],InterpolationOrder->1,Method->"Spline"];
(*The 'procamass' variable is in M0 units and must be converted to M' units by multiplying Mini. The output wr is then in M' units and must be converted back, so we mut divide by Mini!*)
wRinM0units=\[Omega]rof\[Mu][procamass Mini]/Mini;

wr=prec@wRinM0units;
dimlesspin[Mbar_]:=min/(wr Mini Mbar) (1-1/Mbar)+spinin/Mbar^2;
omegahorizon[Mbar_]:=dimlesspin[Mbar]/(2Mbar Mini(1+Sqrt[1-dimlesspin[Mbar]^2]));
massoutputinM0units=prec@(Mbar/.NSolve[omegahorizon[Mbar]==wr/min,Mbar][[1]])Mini;
dimlessspinoutput=prec@dimlesspin[massoutputinM0units/Mini]; (*We need to divide by Mini here, because the chi takes in \bar{M}, rather than stuff in M0 units! But massoutput is written in M0 units!*)

Mfspinf=prec@{massoutputinM0units,dimlessspinoutput};
];
Print["---------------------------------------------------------------------"];
Print["BHEvol Functions:"];
Print["KerrProcaEvolution[m-mode,n-overtone,\[Chi] (dim.-less), Procamass (initial BH mass units)]"];
Print["---------------------------------------------------------------------"];


End[]

EndPackage[]
