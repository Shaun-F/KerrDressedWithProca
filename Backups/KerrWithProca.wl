(* ::Package:: *)

BeginPackage["KerrWithProca`"];


notebookPath = "/home/shaunf/Documents/Computer/Code/projects/Massive_Vector_Field_Dynamical_Friction/ProcaAroundKerr/FKKSsolver/Mathematica/";


RadialMinimize::usage = "Use Nelder Mead to find the minimum value of Log(|R(\!\(\*SubscriptBox[\(r\), \(max\)]\))|) over complex \[Omega]-space"

omegaNRNonRel::usage = "non-relativistic approximation to real part of proca frequency"

omegaNINonRel::usage = "non-relativistic approximation to imaginary part of proca frequency"

nuNNonRel::usage = "non-relativistic approximation to angular eigenvalue"

SolveAngularSystem::usage = "Solve angular FKKS equation in matrix form for kernel vectors"

getRadialSolution::usage = "Solve radial equation numerically for input \[Omega], \[Nu], and parameter association"

FrobeniusSystem::usage = "System of equations resulting from expanding radial function near horizon and solving radial equation order by order"

NelderMead2D::usage = "
NelderMead2D[objFunc, initialSimplex, Precisiongoal, Accuracygoal, Watcher]

objFunc is the HEAD of the function to be minimized and whose argument structure is such that it can be evaluated on a 2D set.\[IndentingNewLine]initialSimplex is the initial simplex on which to begin the algorithm. Must contain 3 points.\[IndentingNewLine]Precision goal sets absolute difference in values of objective function evaluated on nodes of simplex. This is one of the Termination constraints.\[IndentingNewLine]Accuracy goal sets absolute size of simplex. This is the other termination constraints.\[IndentingNewLine]BOTH termination constraints must be satisfied for algorithm to \[IndentingNewLine]terminate. Note: This unusual requirement stems from requirements of the code this algorithm is being embedded in.
Watcher is a boolean that switches wether an iteration counter should be printed during evaluation. DEFAULT: False.
"

functomin::usage = "input is frequency. Output is value of radial function at upper boundary"

Mmatrix::usage = "Angular equation in matrix form"

Mmatrixdet::usage = "Determinant of angular equation"

getNuValue::usage = "Get non-relativistic approximation to angular eigenvalue"

rplusN::usage = "outer horizon in terms of numerical variables"

rminusN::usage = "inner horizon in terms of numerical variables"

CacheResults::usage = "Save data in .mx format to disk in local directory"

assocToString::usage = "Convert association to string in format 'keyvalue'"

styleMarkdown::usage = "pretty printing"

QValue::usage = "Compute Q value (numerical variable) for given \[Mu] and \[Omega]"

QApprox::usage = "linear approximation to Q-value over \[Mu]-space"

calcEndingX::usage = "calculate the ending value of the radial coordinate given input parameter association"

calcOrbital::usage = "calculate the orbital angular momentum number given input parameter association"

getResults::usage = "Retrieve results from disk. second argument is optional and specifies which parameters we wish to retrieve 
					example: getResults[NotebookDirectory[]<>'Solutions/', {{'m',1},{'n',1}}]
"

ToPrecision::usage = "Set precision of expression according to parameter association. Usage: expr//ToPrecision[parameters]"

EmptyQ::usage = "Pattern test to check if list is empty. Used in plotting sector"

Max\[Mu]::usage = "Compute the maximum gravitational coupling for superradiance to be efficient"


Begin["KWP`"];


(*Helper functions*)

(*Some of these are deprecated*)

CacheResults[assoc_, filename_] :=
    (Export[NotebookDirectory[] <> ToString[filename] <> ".mx", assoc
        ]);

assocToString[assoc_] :=
    StringReplace[StringJoin @@ Table[ToString[(assoc // Keys)[[i]], 
        InputForm] <> ToString[(assoc // Values)[[i]], InputForm], {i, 1, Length[
        assoc // Keys]}], {"/" -> "_", "\"" -> ""}];

R2toC[vec_] :=
    vec[[1]] + I * vec[[2]];

IterabletoList[iter_] :=
    Table[iter[[i]], {i, 1, iter // Length}];

DisjunctiontoList[disj_] :=
    Map[#[[2]]&, IterabletoList[disj]];

SymbolQ[sym_] :=
    MatchQ[Head @ sym, _Symbol]

EmptyQ[expr_?ListQ] :=
    If[Length[expr] == 0,
        True
        ,
        False
    ];

Max\[Mu][parameters_] :=
    parameters["m"] * parameters["\[Chi]"] / (2 * rplusN[parameters["\[Chi]"]])

ToPrecision[parameters_][expr_] :=
    SetPrecision[expr, parameters["precision"]];

TotalQuantum[orbital_, overtone_] :=
    (overtone + orbital + 1);

InsertValues[parameters_][expr_] :=
    expr /. {\[Mu]N -> parameters["\[Mu]Nv"], m -> parameters["m"], \[Eta] -> parameters[
        "\[Eta]"], n -> parameters["n"], l -> parameters["l"], s -> parameters["s"
        ], \[Chi] -> parameters["\[Chi]"]};

GetIngoingSolution[frobsol_] :=
    Block[{
        inx =
            MapIndexed[
                    If[Re[#1[[1]]] <= 0,
                        #2
                        ,
                        Unevaluated[Sequence[]]
                    ]&
                    ,
                    frobsol[[All, All, 2]]
                ] // First
    },
        frobsol[[inx]]
    ];

JfromLS[l_, s_] :=
    Abs[l + s] // Rationalize[#, 0]&;

QValue[\[Mu]_, \[Omega]_] :=
    Sqrt[\[Mu]^2 - \[Omega]^2]

styleMarkdown[str_String] :=
    Module[{style, applyRules},
        style[s_String, {weight_:Plain, slant_:Plain}] := StringJoin[
            "\!\(\*StyleBox[\"", s, "\", FontWeight -> ", ToString @ weight, ", FontSlant -> ",
             ToString @ slant, "]\)"];
        applyRules[s_String, rules_] := StringReplace[s, # ~~ Shortest
             @ x__ ~~ # :> style[x, #2]& @@@ rules];
        applyRules[str, {"***" -> {Bold, Italic}, "**" -> {Bold}, "*"
             -> {Plain, Italic}}]
    ]

QApprox[\[Mu]_, Q1_, Q2_, \[Mu]1_, \[Mu]2_] :=
    (Q2 - Q1) / (\[Mu]2 - \[Mu]1) * \[Mu] + (Q1 * \[Mu]2 - Q2 * \[Mu]1) / (\[Mu]2 - \[Mu]1)

calcEndingX[parameters_] :=
    2 * 10 * (parameters["m"] + parameters["n"]) / parameters["\[Mu]Nv"] ^
         2

calcOrbital[parameters_] :=
    parameters["m"] + parameters["s"]

getResults[absolDir_?StringQ, parameterSet : (_?ListQ | Null) : Null] :=
    Module[{localfiles = FileNames[All, absolDir]},
        indexer = StringContainsQ[#, "RunData"]& /@ localfiles;
        filenames = Pick[localfiles, indexer];
        If[parameterSet == Null,
            Return[Import /@ filenames, Module]
        ];
        If[ListQ @ parameterSet,
            If[Length @ Dimensions @ parameterSet == 1,
                ParameterSetIndicies = StringContainsQ[ToString[parameterSet
                    [[1]]] <> ToString[parameterSet[[2]]]] /@ filenames;
                Return[Import /@ Pick[filenames, ParameterSetIndicies
                    ], Module];
            ];
            If[Length @ Dimensions @ parameterSet > 1,
                ParameterSetIndicies = Table[StringContainsQ[ToString[
                    parameterSet[[i]][[1]]] <> ToString[parameterSet[[i]][[2]]]] /@ filenames,
                     {i, 1, Length @ parameterSet}];
                ParameterSetIndiciesAll = Table[AllTrue[ParameterSetIndicies
                    [[All, i]], TrueQ], {i, 1, Length @ First @ ParameterSetIndicies}];
                Return[Import /@ Pick[filenames, ParameterSetIndiciesAll
                    ], Module];
            ];
        ];
    ];


(*Angular functions*)

qth[a_, \[Nu]_, \[Theta]_] :=
    1 - a^2 * \[Nu]^2 * Cos[\[Theta]] ^ 2;

Kth[m_, a_, \[Omega]_, \[Theta]_] :=
    m - a * \[Omega] * Sin[\[Theta]] ^ 2;

\[Sigma][a_, m_, \[Omega]_, \[Nu]_] :=
    a * \[Nu]^2 (m - a * \[Omega]) + \[Omega];

\[CapitalLambda][a_, m_, \[Omega]_, \[Nu]_, \[Mu]_] :=
    \[Mu]^2 / \[Nu]^2 - \[Sigma][a, m, \[Omega], \[Nu]] / \[Nu] + 2 * a * \[Omega] * m - a^2 * \[Omega]^2;

\[Gamma][\[Omega]_, \[Mu]_] :=
    \[Omega]^2 - \[Mu]^2;

diffeqang[s_, m_, a_, \[Omega]_, \[Nu]_] :=
    Collect[qth[a, \[Nu], \[Theta]] / Sin[\[Theta]] * D[Sin[\[Theta]] / qth[a, \[Nu], \[Theta]] * D[s, \[Theta]],
         \[Theta]] - qth (Kth[m, a, \[Omega], \[Theta]] ^ 2 / (qth[a, \[Nu], \[Theta]] * Sin[\[Theta]] ^ 2) + (2 - qth[
        a, \[Nu], \[Theta]]) / qth[a, \[Nu], \[Theta]] ^ 2 * \[Sigma][a, m, \[Omega], \[Nu]] / \[Nu] - \[Mu]^2 / \[Nu]^2) s, {s, 
        Derivative[1][s], s^\[Prime]\[Prime]}];

spherharmons = Flatten @ Table[Row[{{l, m}, Spacer[20], SphericalHarmonicY[
    l, m, \[Theta], 0]}], {l, 0, 4}, {m, -l, l}] // MatrixForm;

Off[ClebschGordan::tri];

Off[ClebschGordan::phy];

dl[l_] :=
    2 * l + 1;

sym[l1_, l2_, l3_, m_] :=
    (-1) ^ m * Sqrt[dl[l1] * dl[l2] * dl[l3] / (4 * \[Pi])] * ThreeJSymbol[
        {l1, 0}, {l2, 0}, {l3, 0}] * ThreeJSymbol[{l1, -m}, {l2, 0}, {l3, m}];

all[l_, lp_, m_] :=
    4 / 3 Sqrt[\[Pi] / 5] sym[l, 2, lp, m] + Sqrt[4 * \[Pi]] / 3 * sym[l, 0, 
        lp, m];

bll[l_, lp_, m_] :=
    2 * Sqrt[\[Pi]] / 5 * sym[l, 0, lp, m] + 8 * Sqrt[\[Pi] / 5] / 7 sym[l, 2,
         lp, m] + 16 * Sqrt[\[Pi]] / 105 * sym[l, 4, lp, m];

dll[l_, lp_, m_] :=
    Sqrt[4 \[Pi] / 3] (lp * Sqrt[(lp + 1) ^ 2 - m^2] / Sqrt[(2 lp + 1) (2
         lp + 3)] * sym[l, 1, lp + 1, m] - (lp + 1) Sqrt[lp^2 - m^2] / Sqrt[(
        2 lp + 1) (2 lp - 1)] * sym[l, 1, lp - 1, m]);

Mmatrixdef =
    Function[{l, lp, m, a, \[Omega], \[Nu], \[Mu]},
        (\[CapitalLambda][a, m, \[Omega], \[Nu], \[Mu]] - lp * (lp + 1)) * KroneckerDelta[l, lp] + 
            (-\[Nu]^2 * \[CapitalLambda][a, m, \[Omega], \[Nu], \[Mu]] + \[Nu]^2 * lp * (lp + 1) - 2 * \[Sigma][a, m, \[Omega], \[Nu]] * 
            \[Nu] + \[Gamma][\[Omega], \[Mu]]) a^2 * all[l, lp, m] - 2 * a^2 * \[Nu]^2 * dll[l, lp, m] - \[Gamma][
            \[Omega], \[Mu]] * \[Nu]^2 * a^4 * bll[l, lp, m]
    ];

Mmatrix =
    Function[{KMax, eta, m, a, \[Omega], \[Nu], \[Mu]},
        Table[Mmatrixdef[m + 2 * k + eta, m + 2 * kp + eta, m, a, \[Omega], 
            \[Nu], \[Mu]], {k, 0, KMax}, {kp, 0, KMax}]
    ];

Mmatrixdet[\[Omega]N_, \[Nu]N_, parameters_] :=
    Block[{
        Mmatrixdettemp =
            Evaluate[
                \[Nu]N ^ parameters["KMax"] * Mmatrix[parameters["KMax"],
                     parameters["\[Eta]"], parameters["m"], parameters["\[Chi]"] * M, \[Omega]N / M, \[Nu]N / 
                    M, parameters["\[Mu]Nv"] / M] //
                Det //
                Refine //
                Collect[#, M]& //
                Collect[#, \[Nu]N]& //
                Refine
            ]
    },
        Evaluate @
            Block[{eta = parameters["\[Eta]"], m = parameters["m"], \[Chi] = parameters[
                "\[Chi]"], \[Mu]N = parameters["\[Mu]Nv"]},
                Evaluate @ Mmatrixdettemp
            ]
    ];

ansatz[b1_] :=
    {1, b1, 0};

SolveAngularSystem[\[Omega]value_?NumberQ, \[Nu]value_?NumberQ, parameters_] :=
    Block[{bcoeffs = Table[C[i], {i, 0, parameters["KMax"]}], \[Eta] = parameters[
        "\[Eta]"], m = parameters["m"]},
        NumMatrix = \[Nu]N ^ -7 * Mmatrix[parameters["KMax"], \[Eta], m, \[Chi] * M,
             \[Omega]N / M, \[Nu]N / M, \[Mu]N / M] // N;
        equalto = Table[0, {i, 0, parameters["KMax"]}];
        System = (NumMatrix . bcoeffs == equalto) /. {\[Omega]N -> \[Omega]value, \[Nu]N
             -> \[Nu]value, \[Mu]N -> parameters["\[Mu]Nv"], m -> parameters["m"], \[Eta] -> parameters[
            "\[Eta]"], \[Chi] -> parameters["\[Chi]"]};
        Return[Flatten @ {C[0] -> 1, Solve[System /. C[0] -> 1, bcoeffs
            [[2 ;; parameters["KMax"] + 1]]]}];
    ];


(*Radial functions*)

qr[\[Nu]_, r_] :=
    1 + \[Nu]^2 * r^2;

\[Sigma][a_, m_, \[Omega]_, \[Nu]_] :=
    a * \[Nu]^2 (m - a * \[Omega]) + \[Omega];

Kr[a_, m_, r_] :=
    -a * m + (a^2 + r^2) \[Omega];

rplus = (rs + Sqrt[rs^2 - 4 * a^2]) / 2;

rminus = (rs - Sqrt[rs^2 - 4 * a^2]) / 2;

\[CapitalDelta] = (r - rplus) (r - rminus);

diffeqrad = D[\[CapitalDelta] / qr * D[R[r], r], r] + (Kr[a, m, r] ^ 2 / (qr * \[CapitalDelta]) +
     (2 - qr) / qr^2 * \[Sigma][a, m, \[Omega], \[Nu]] / \[Nu] - \[Mu]^2 / \[Nu]^2) R[r] == 0;

diffeqradOp = (D[\[CapitalDelta] / qr * D[#, r], r] + (Kr[a, m, r] ^ 2 / (qr * \[CapitalDelta]) +
     (2 - qr) / qr^2 * \[Sigma][a, m, \[Omega], \[Nu]] / \[Nu] - \[Mu]^2 / \[Nu]^2) #)&;

\[CapitalKappa][\[Omega]_, \[Chi]_, m_] :=
    (\[Omega] * rplusN[\[Chi]] * (rplusN[\[Chi]] + rminusN[\[Chi]]) - \[Chi] * m) / (rplusN[\[Chi]] -
         rminusN[\[Chi]]);

rplusN[\[Chi]_] :=
    1 + (1 - \[Chi]^2) ^ (1 / 2);

rminusN[\[Chi]_] :=
    1 - (1 - \[Chi]^2) ^ (1 / 2);

KerrSurfaceRotationN[\[Chi]_] :=
    (1 / 2) (\[Chi] / rplusN[\[Chi]]);

qrN = qr[\[Nu], rN * M];

rpN = rplusN[\[Chi]] // Rationalize;

rmN = rminusN[\[Chi]] // Rationalize;

rsN = rpN + rmN // Rationalize;

\[CapitalDelta]N = (rN - rpN) * (rN - rmN) // Rationalize;

\[Sigma]N =
    \[Sigma][a, m, \[Omega], \[Nu]] * M //
    Expand //
    Rationalize;

KrN =
    Kr[a, m, rN * M] / M //
    Rationalize //
    Simplify;

qrNX =
    qrN /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]] //
    Rationalize //
    FullSimplify;

KrNX =
    KrN /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]] //
    Rationalize //
    FullSimplify;

\[Sigma]NX =
    \[Sigma]N /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]]//
    Rationalize //
    FullSimplify;

\[CapitalDelta]NX =
    \[CapitalDelta]N /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]] //
    Rationalize //
    FullSimplify;

P1 = (2 * rN - rsN) / (rN - rminusN[\[Chi]]) // Rationalize;

P2 = -(2 * rN * \[Nu]N^2) / qrN;

Q1 = KrN^2 / (rN - rminusN[\[Chi]]) ^ 2;

Q2 = qrN / (rN - rminusN[\[Chi]]) * ((2 - qrN) / qrN^2 * \[Sigma]N / \[Nu]N - \[Mu]N^2 / 
    \[Nu]N^2) // Rationalize;

P1x = P1 /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]];

P2x = P2 /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]];

Q1x = Q1 /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]];

Q2x = Q2 /. rN -> xN * (rplusN[\[Chi]] - rminusN[\[Chi]]) + rplusN[\[Chi]];


DiffRadOp = (\!\(
\*SubscriptBox[\(\[PartialD]\), \(rN, rN\)]#\) + (qrN/\[CapitalDelta]N)*((2rN-rsN)/qrN-(2rN*\[Nu]N^2*\[CapitalDelta]N)/qrN^2)\!\(
\*SubscriptBox[\(\[PartialD]\), \(rN\)]#\) + (qrN/\[CapitalDelta]N)(-(\[Mu]N/\[Nu]N)^2+KrN^2/(\[CapitalDelta]N*qrN)+((2-qrN)/qrN^2)*(\[Sigma]N/\[Nu]N))#)&;
DiffRadOpx = (\!\(
\*SubscriptBox[\(\[PartialD]\), \(xN, xN\)]#\) + (P1x/xN+P2x*(rplusN[\[Chi]]-rminusN[\[Chi]]))\!\(
\*SubscriptBox[\(\[PartialD]\), \(xN\)]#\) + (Q1x/xN^2+Q2x (rplusN[\[Chi]]-rminusN[\[Chi]])/xN)#)&;

FrobeniusSeries[x_,\[Kappa]_,c1_,c2_]:=x^(-I*\[Kappa])*(1 + c1*x + c2*x^2);
(*FrobeniusSeries[xN,\[Kappa],C[1],C[2]]; (xN^(I*\[Kappa]+2)DiffRadOpx@%)//Simplify; Series[%, {xN,0,2}]\[Equal]0;
FrobeniusSystem = Simplify[LogicalExpand[%], TimeConstraint\[Rule]Infinity]
FrobeniusSystem[\[Kappa]_,c1_,c2_,\[Nu]N_,\[Omega]N_,parameters_:parameters]:=Evaluate[
Block[{\[Chi]=parameters["\[Chi]"],m=parameters["m"], \[Mu]N=parameters["\[Mu]Nv"]},
Block[{fs = FrobeniusSeries[xN,\[Kappa],c1,c2]},
Block[{expr = xN^(I*\[Kappa]+2)DiffRadOpx@fs},
Return@LogicalExpand[Series[expr, {xN,0,2}]\[Equal]0]]]]]*)
FrobeniusSystem[\[Kappa]_,c1_,c2_,\[Nu]Nt_,\[Omega]Nt_,parameters_]:=Evaluate[
Evaluate[Import[notebookPath<>"Expressions/FrobeniusSystem.mx"][\[Kappa],c1,c2,\[Omega]Nt, \[Nu]Nt, parameters["m"], parameters["\[Chi]"], parameters["\[Mu]Nv"]]]
];

getRadialSolution[w_?NumberQ, v_?NumberQ, parameters_]:=
	Module[{\[Omega]vv = w, \[Nu]vv=v,\[Chi]=parameters["\[Chi]"],m=parameters["m"],\[Mu]N=parameters["\[Mu]Nv"]},
		StartingRadius =  parameters["StartingX"];
		EndingRadius =  parameters["EndingX"];
		FrobeniusSolution =(NSolve[FrobeniusSystem[\[Kappa],c1,c2,\[Nu]vv,\[Omega]vv,parameters]//Rationalize[#,0]&, {\[Kappa],c1,c2}, VerifySolutions->True, WorkingPrecision->parameters["precision"]][[parameters["branch"]]])[[All,2]];
		FrobR0 = (Limit[FrobeniusSeries[x,Sequence@@#], x->parameters["\[Epsilon]"]]&)@ FrobeniusSolution;
		FrobRPrime0 = (Limit[D[FrobeniusSeries[x,Sequence@@#],x], x-> parameters["\[Epsilon]"]]&)@FrobeniusSolution;
		BC = Thread@{FrobR0, FrobRPrime0}//ToPrecision[parameters];
		eq = {((DiffRadOpx@R[xN]/.{a->\[Chi]*M,\[Omega]->\[Omega]N/M,\[Nu]->\[Nu]N/M,\[Mu]->\[Mu]N/M}//Collect[#,M]&//InsertValues[parameters]//ToPrecision[parameters])/.{ \[Nu]N->\[Nu]vv, \[Omega]N->\[Omega]vv})==0, R[StartingRadius]==FrobR0// ToPrecision[parameters], Derivative[1][R][StartingRadius]==FrobRPrime0// ToPrecision[parameters]};
		RadialSolutions = NDSolve[eq// ToPrecision[parameters],R, {xN, StartingRadius, EndingRadius},WorkingPrecision->parameters["precision"], MaxSteps->parameters["integratorMaxStep"], PrecisionGoal->parameters["integratorPrecisionGoal"]]//First;
		Return[R/.RadialSolutions, Module]
];


(*Non-relativistic formulas*)

omegaNRNonRel[parameters_] :=
    Block[{\[Mu]N = parameters["\[Mu]Nv"], m = parameters["m"], l = parameters[
        "l"], n = parameters["n"], s = parameters["s"]},
        \[Mu]N * (1 - \[Mu]N^2 / (2 * TotalQuantum[l, n] ^ 2)) // Rationalize[#, 0]&
    ];

gfuncN[j_, parameters_] :=
    Block[{\[Chi] = parameters["\[Chi]"], \[Mu]N = parameters["\[Mu]Nv"], m = parameters[
        "m"], l = parameters["l"], n = parameters["n"], s = parameters["s"]},
        Product[k^2 * (1 - \[Chi]^2) + (\[Chi] * m - 2 * rplusN[\[Chi]] * omegaNRNonRel[parameters]) ^ 2, {k, 1, j}] // Rationalize[#, 0]&
    ];

Ccoeff[n_, l_, j_] :=
    With[{NN = n + l + 1},
        If[NN >= 0,
            ((2 ^ (2 l + 2 j + 1) * (NN + l)!) / (NN ^ (2 l + 4) * (NN- l - 1)!)) (l! / ((l + j)! (l + j + 1)!)) ^ 2 * (1 + (2 (1 + l - j)(1 - l + j)) / (l + j)) ^ 2 // Rationalize
            ,
            Print["Error: Principal quantum number less than zero!"];
            Abort[]
        ]
    ];

omegaNINonRel[parameters_] :=
    Block[{\[Chi] = parameters["\[Chi]"], \[Mu]N = parameters["\[Mu]Nv"], m = parameters[
        "m"], l = parameters["l"], n = parameters["n"], s = parameters["s"]},
        Block[{j = Abs[m]},
            2 * rplusN[\[Chi]] * Ccoeff[n, l, j] * gfuncN[j, parameters] *(m * KerrSurfaceRotationN[\[Chi]] - omegaNRNonRel[parameters]) * \[Mu]N ^ (2 l + 2 j + 5)
        ]
    ];

nuNNonRel[parameters_] :=
    Block[{\[Chi] = parameters["\[Chi]"], \[Mu]N = parameters["\[Mu]Nv"], m = parameters[
        "m"], l = parameters["l"], n = parameters["n"], PolarizationSigned = 
        parameters["s"]},
        Block[{\[Omega]NNonRel = omegaNRNonRel[parameters] + omegaNINonRel[parameters
            ] * I},
            Evaluate @
                Switch[PolarizationSigned,
                    -1,(-\[Omega]NNonRel) / (m - \[Chi] * \[Omega]NNonRel),
                    0,(1 / (2 * \[Chi])) (m + 1 - \[Chi] * \[Omega]NNonRel + \[Sqrt]((\[Chi] * \[Omega]NNonRel - m - 1) ^ 2 + 4 * \[Chi] * \[Omega]NNonRel)),
                    1,Evaluate @Block[{X = \[Chi], \[Nu]},
                                Roots[X * \[Nu]^3 (m - X * \[Omega]NNonRel) - \[Nu]^
                                    2 ((m + 1) (m + 2) - X * \[Omega]NNonRel (2 m - X * \[Omega]NNonRel)) + \[Omega]NNonRel * 
                                    \[Nu] + \[Omega]NNonRel^2 == 0, \[Nu]][[2]]
                            ]
                ]
        ]
    ];(*Block[{\[Chi]=X},Roots[\[Chi]*\[Nu]^3(m - \[Chi]*\[Omega]NNonRel)-\[Nu]^2((m+1)(m+2) - \[Chi]*\[Omega]NNonRel(2m-\[Chi]*\[Omega]NNonRel)) + \[Omega]NNonRel*\[Nu] + \[Omega]NNonRel^2\[Equal]0,\[Nu]][[2]]//Simplify]*)

nuNNonRel[\[Omega]NNonRel_, parameters_] :=
    Block[{\[Chi] = parameters["\[Chi]"], m = parameters["m"], PolarizationSigned
         = parameters["s"]},
        Evaluate @
            Switch[PolarizationSigned,
                -1,(-\[Omega]NNonRel) / (m - \[Chi] * \[Omega]NNonRel),
                0,(1 / (2 * \[Chi])) (m + 1 - \[Chi] * \[Omega]NNonRel + \[Sqrt]((\[Chi] * \[Omega]NNonRel- m - 1) ^ 2 + 4 * \[Chi] * \[Omega]NNonRel)),
                1,Evaluate @
                        Block[{X = \[Chi], \[Nu]},
                            Roots[X * \[Nu]^3 (m - X * \[Omega]NNonRel) - \[Nu]^2 ((
                                m + 1) (m + 2) - X * \[Omega]NNonRel (2 m - X * \[Omega]NNonRel)) + \[Omega]NNonRel * \[Nu] + 
                                \[Omega]NNonRel^2 == 0, \[Nu]][[2]]
                        ]
            ]
    ];

getNuValue[\[Omega]N_?NumberQ, parameters_, ClosestGuess_?NumberQ] :=
    Block[{eta = parameters["\[Eta]"], m = parameters["m"], \[Chi] = parameters[
        "\[Chi]"], \[Mu]N = parameters["\[Mu]Nv"]},
        Module[{
            thedet =
                Function[{w, \[Nu]N},
                    Evaluate @ Mmatrixdet[w, \[Nu]N, parameters]
                ]
        },
            FindRoot[thedet[\[Omega]N, \[Nu]N] // ToPrecision[parameters], {\[Nu]N, 
                ClosestGuess}, WorkingPrecision -> parameters["precision"], MaxIterations
                 -> 1000][[1]][[2]] //
            Rationalize[#, 0]& //
            Return
        ]
    ]; (*Min@Cases[%,Power[\[Nu]N,x_?NumberQ]\[RuleDelayed]x,-1]*) 


functomin[parameters_?AssociationQ, \[Nu]fitted_:Null][\[Omega]realimaglist_?ListQ] :=
    Block[{\[Omega]value = \[Omega]realimaglist[[1]] + I * \[Omega]realimaglist[[2]]},
        Block[{
            nuroot =Evaluate @If[TrueQ[\[Nu]fitted == Null],
                                   Evaluate @ ToPrecision[parameters]@getNuValue[\[Omega]value, parameters, ToPrecision[parameters]@nuNNonRel[\[Omega]value,parameters]]
                                   ,
                                   Evaluate @ ToPrecision[parameters]@getNuValue[\[Omega]value, parameters, ToPrecision[parameters]@\[Nu]fitted]
                                 ]
              },
              retval = getRadialSolution[\[Omega]value, nuroot, parameters][parameters["EndingX"]] // Abs;
            (*Return large negative value if retval is 0. Cures errors
                 raised by taking naive logarithm*)
              If[retval == 0,
                Return[-100]
                ,
                Return[Log10[retval]]
            ];
        ]
    ];


RadialMinimize[\[Omega]Guess_?NumberQ, \[Nu]Guess_?NumberQ, OmegaBoundary_?ListQ,
     parameters_?AssociationQ, Watcher_?BooleanQ, \[Nu]fit_:Null] :=
    Module[{currval = {0, 0}, CurrentIteration = 0},
        RMMessengerGenerator =
            Function[{X, counter},
                "Current Value: " <> ToString[X, InputForm] <> "\n Current Step: "
                     <> ToString[counter, InputForm]
            ];
        RMMessenger = RMMessengerGenerator[currval, CurrentIteration]
            ;
        (*If[Watcher, PrintTemporary@Dynamic[RMMessenger]];*)
        SimplexBase = {\[Omega]Guess // Re, \[Omega]Guess // Im} // ToPrecision[parameters
            ];
        SimplexNodeOne = SimplexBase - {(OmegaBoundary[[2]] - \[Omega]Guess)
             / parameters["omegaRes"], 0} // Re;
        SimplexNodeTwo = SimplexBase - {0, 10 ^ (Log10[Abs @ Im[\[Omega]Guess
            ]] - 1)};
        InitialSimplex = {SimplexBase, SimplexNodeOne, SimplexNodeTwo
            };
        If[TrueQ[\[Nu]fit == Null],\[Nu]fitted=Null,\[Nu]fitted = \[Nu]fit[parameters["\[Mu]Nv"]]];
        result = 
                NMinimize[
                    {Hold @ functomin[parameters, \[Nu]fitted][{\[Xi], \[Psi]}], \[Psi] > 0},{\[Xi], \[Psi]},
                    Method -> {"NelderMead", "PostProcess" -> False, 
                        "InitialPoints" -> InitialSimplex},
                    WorkingPrecision -> parameters["OptimizationAccuracy"],
                    EvaluationMonitor :>(RMMessenger = RMMessengerGenerator[\[Xi] + I * \[Psi], CurrentIteration];
                                        CurrentIteration++;
                                        ),
                    MaxIterations -> parameters["MaxInterationsMinimization"]
                ][[2]];
        \[Omega]Result = (\[Xi] + I * \[Psi]) /. result;
        If[TrueQ[\[Nu]fit == Null],
            \[Nu]Result = getNuValue[\[Omega]Result, parameters, nuNNonRel[\[Omega]Result,
                 parameters]]
            ,
            \[Nu]Result = getNuValue[\[Omega]Result, parameters, \[Nu]fitted];
        ];
        Return[<|"\[Omega]" -> \[Omega]Result, "\[Nu]" -> \[Nu]Result|>, Module];
    ];


End[];


EndPackage[];
