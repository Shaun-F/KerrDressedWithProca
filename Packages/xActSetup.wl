(* ::Package:: *)

xActSetup (*Do Not Remove, used for recognition of package already being imported*)
M=1;

Block[{Print},Needs["xAct`xCoba`"]]

Off[ValidateSymbol::used]
Off[ValidateSymbol::invalid]
$DefInfoQ=False;
$UndefInfoQ=False;
$CVSimplify=Identity;
(*$Assumptions= {r>0,M>0,J>0,1>a\[GreaterEqual]0,\[Pi]\[GreaterEqual]\[Theta]\[GreaterEqual]0,\[Omega]\[Element]Complexes, \[Nu]\[Element]Complexes, R[r]\[NotEqualTilde]0, S[\[Theta]]\[NotEqualTilde]0, \[Sigma]\[NotEqualTilde]0, \[CapitalSigma]\[NotEqualTilde]0, \[CapitalDelta]\[NotEqualTilde]0};*)
(*Since were only using one covariant derivative, redefine printed symbols*)
PrintAs[RiemannCD]^="R";
PrintAs[RicciCD]^="R";
PrintAs[RicciScalarCD]^="R";
PrintAs[ChristoffelCD]^="\[CapitalGamma]";
PrintAs[EinsteinCD]^="G";

concise[expr_]:=Tooltip[\[GrayCircle], expr]
ExtractBiPair[ass_, pair_]:=<|pair[[1]]->First[ass][pair[[1]]], pair[[2]]-> Last[ass][pair[[2]]]|>

(*Define basic objects of xAct*)
DefManifold[\[ScriptCapitalM],4, {\[Alpha],\[Beta],\[Rho],\[Gamma],\[Delta],\[Iota],\[Zeta],\[Xi]}]; (*Manifold of our spacetime*)
DefChart[sphericalchart, \[ScriptCapitalM],{0,1,2,3},{t[],r[],\[Theta][],\[Phi][]}, FormatBasis->{"Partials", "Differentials"}, ChartColor->Blue]; (*Spherical co-ordinate chart*)
(*Cast numerical indices to latin/greek representation*)
ch/:CIndexForm[0,ch]:="t";
ch/:CIndexForm[1,ch]:="r";
ch/:CIndexForm[2,ch]:="\[Theta]";
ch/:CIndexForm[3,ch]:="\[Phi]";

(*Define the constants of the system*)
DefConstantSymbol[M] (*Mass of BH*)
DefConstantSymbol[J](*Angular momentum of BH*)
DefConstantSymbol[a](*Dimensionless spin*)
DefConstantSymbol[G](*Newtons constant*)
DefConstantSymbol[\[Omega]](*angular frequency*)
DefConstantSymbol[\[Omega]r];
DefConstantSymbol[\[Omega]i];
DefConstantSymbol[m](*projected angular momentum quantum number*)
DefConstantSymbol[\[CapitalMu]](*Mass of proca field*);
DefConstantSymbol[q];
DefConstantSymbol[\[Pi]];
DefConstantSymbol[\[Chi]];
DefConstantSymbol[\[Nu]];
DefConstantSymbol[\[Nu]r];
DefConstantSymbol[\[Nu]i];
DefConstantSymbol[\[Mu]];

(*Define some scalar fields*)
DefScalarFunction[Z];
DefScalarFunction[S];
DefScalarFunction[R];
DefScalarFunction[Sr];
DefScalarFunction[Rr];
DefScalarFunction[Si];
DefScalarFunction[Ri];


(*Define metric object*)
DefMetric[-1,g[-\[Zeta],-\[Xi]], CD, {";", "\[Del]"},PrintAs->"g"]

(*G=c=1*)
KerrComponents = {{-1+(2 M r)/((J^2 Cos[\[Theta]]^2)/M^2+r^2),0,0,-((2 J r*Sin[\[Theta]]^2)/((J^2 Cos[\[Theta]]^2)/M^2+r^2))},{0,(J^2 Cos[\[Theta]]^2+M^2 r^2)/(J^2+M^2 r (-2 M+r)),0,0},{0,0,(J^2 Cos[\[Theta]]^2)/M^2+r^2,0},{-((2 J r Sin[\[Theta]]^2)/((J^2 Cos[\[Theta]]^2)/M^2+r^2)),0,0,Sin[\[Theta]]^2 (J^2/M^2+r^2+(2 J^2 M r Sin[\[Theta]]^2)/(J^2 Cos[\[Theta]]^2+M^2 r^2))}}/.J->a*M//FullSimplify;
IKerrComponents = Inverse[KerrComponents]//FullSimplify;
met = CTensor[ToxActVariables@KerrComponents, {-sphericalchart,-sphericalchart}];
SetCMetric[met,-sphericalchart,SignatureOfMetric->{3,1,0}];
g=met;

(*Precompute all metric-derived quantities, e.g. christoffel, riemann, etc.*)
(*MetricCompute[g,sphericalchart,All]; *)

(*Define the covariant derivative of this CTensor metric*)
Cd = CovDOfMetric[g];

(*Define killing vectors*)
DefTensor[\[CapitalTau][\[Zeta]],\[ScriptCapitalM],PrintAs->"\[ScriptCapitalT]"]
DefTensor[\[CapitalPhi][\[Zeta]], \[ScriptCapitalM],PrintAs->"\[CapitalPhi]"]
Block[{Print},ComponentValue[ComponentArray[\[CapitalTau][{\[Zeta],sphericalchart}]], {1,0,0,0}]];
Block[{Print},ComponentValue[ComponentArray[\[CapitalPhi][{\[Zeta],sphericalchart}]], {0,0,0,1}]];
TimelikeKilling = ToCTensor[\[CapitalTau],{sphericalchart}]/.TensorValues[\[CapitalTau]];
SpacelikeKilling = ToCTensor[\[CapitalPhi],{sphericalchart}]/.TensorValues[\[CapitalPhi]];


PolarizationComponents={{(-((a^2+r^2)^2/((a^2+(-2+r) r) (r^2+1/\[Nu]^2)))+(a^2 Sin[\[Theta]]^2)/(1/\[Nu]^2-a^2 Cos[\[Theta]]^2))/(\[Nu]^2 (r^2+a^2 Cos[\[Theta]]^2)),(I r (a^2+r^2) \[Nu])/((1+r^2 \[Nu]^2) (r^2+a^2 Cos[\[Theta]]^2)),-((I a^2 \[Nu] Cos[\[Theta]] Sin[\[Theta]])/((r^2+a^2 Cos[\[Theta]]^2) (-1+a^2 \[Nu]^2 Cos[\[Theta]]^2))),(a (-((a^2+r^2)/((a^2+(-2+r) r) (r^2+1/\[Nu]^2)))+1/(1/\[Nu]^2-a^2 Cos[\[Theta]]^2)))/(\[Nu]^2 (r^2+a^2 Cos[\[Theta]]^2))},{-((I r (a^2+r^2) \[Nu])/((1+r^2 \[Nu]^2) (r^2+a^2 Cos[\[Theta]]^2))),(a^2+(-2+r) r)/((1+r^2 \[Nu]^2) (r^2+a^2 Cos[\[Theta]]^2)),0,-((I a r \[Nu])/((1+r^2 \[Nu]^2) (r^2+a^2 Cos[\[Theta]]^2)))},{(I a^2 \[Nu] Cos[\[Theta]] Sin[\[Theta]])/((r^2+a^2 Cos[\[Theta]]^2) (-1+a^2 \[Nu]^2 Cos[\[Theta]]^2)),0,-(1/((r^2+a^2 Cos[\[Theta]]^2) (-1+a^2 \[Nu]^2 Cos[\[Theta]]^2))),(I a \[Nu] Cot[\[Theta]])/((r^2+a^2 Cos[\[Theta]]^2) (-1+a^2 \[Nu]^2 Cos[\[Theta]]^2))},{(a (-((a^2+r^2)/((a^2+(-2+r) r) (r^2+1/\[Nu]^2)))+1/(1/\[Nu]^2-a^2 Cos[\[Theta]]^2)))/(\[Nu]^2 (r^2+a^2 Cos[\[Theta]]^2)),(I a r \[Nu])/((1+r^2 \[Nu]^2) (r^2+a^2 Cos[\[Theta]]^2)),-((I a \[Nu] Cot[\[Theta]])/((r^2+a^2 Cos[\[Theta]]^2) (-1+a^2 \[Nu]^2 Cos[\[Theta]]^2))),(-(a^2/((a^2+(-2+r) r) (r^2+1/\[Nu]^2)))+Csc[\[Theta]]^2/(1/\[Nu]^2-a^2 Cos[\[Theta]]^2))/(\[Nu]^2 (r^2+a^2 Cos[\[Theta]]^2))}};
(*h = -(r*dr + a^2*Sin[\[Theta]]*Cos[\[Theta]]*d\[Theta])\[TensorWedge]dt + a*Sin[\[Theta]]*(r*Sin[\[Theta]]*dr + (r^2+a^2)*Cos[\[Theta]]*d\[Theta])\[TensorWedge]d\[Phi];
Bsymb = Table[Symbol["B"<>ToString[i]<>ToString[j]], {i,1,4},{j,1,4}];
sys = Table[Sum[Bsymb[[\[Alpha],\[Gamma]]]*(KerrComponents[[\[Gamma],\[Beta]]] + I*h[[\[Gamma],\[Beta]]]*\[Nu]), {\[Gamma],1,4}], {\[Beta],1,4},{\[Alpha],1,4}]//Simplify;
equ =Flatten[ Thread/@Thread[sys==DiagonalMatrix[{1,1,1,1}]]];
Coeffs = Flatten@Bsymb;
sol = Solve[equ, Coeffs, Assumptions\[Rule]True]//Simplify;
PolarizationComponents = Bsymb/.sol//FullSimplify//First;
*)
(*Define the polarization tensor object, then assign values to it*)
DefTensor[B[\[Zeta],\[Xi]], \[ScriptCapitalM], PrintAs-> "B"]; (*Define the tensor object*)
Block[{Print},ComponentValue[ComponentArray[B[{\[Zeta],sphericalchart},{\[Xi],sphericalchart}]], PolarizationComponents//ToxActVariables]];

Polarization = (ToCTensor[B,{sphericalchart,sphericalchart}]/.TensorValues[B]);

(*See arXiv:2004.09536*) 
dr = {0,1,0,0};
dt = {1,0,0,0};
d\[Theta] = {0,0,1,0};
d\[Phi] = {0,0,0,1};
DefTensor[H[\[Zeta],\[Xi]], \[ScriptCapitalM], PrintAs->"h"]
hcomps = -(r*dr + a^2*Sin[\[Theta]]*Cos[\[Theta]]*d\[Theta])\[TensorWedge]dt + a*Sin[\[Theta]]*(r*Sin[\[Theta]]*dr + (r^2+a^2)*Cos[\[Theta]]*d\[Theta])\[TensorWedge]d\[Phi];
Block[{Print},ComponentValue[ComponentArray[H[{-\[Zeta],-sphericalchart},{-\[Xi],-sphericalchart}]], hcomps//ToxActVariables]];
CKT = ToCTensor[H, {-sphericalchart, -sphericalchart}]/.TensorValues[H];
Block[{res = Head[Polarization[\[Zeta],\[Xi]](g[-\[Xi],-\[Iota]]+I*\[Nu] CKT[-\[Xi],-\[Iota]])//Simplify][[1]]},
If[TrueQ[res!=DiagonalMatrix[{1,1,1,1}]], Print["Error in polarization tensor definition"]]]
