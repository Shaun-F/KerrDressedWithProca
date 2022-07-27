(* ::Package:: *)

If["HelperFunctions"\[NotElement]Names["Global`"], Import[$FKKSRoot<>"Packages/HelperFunctions.wl"]];


RadialMinimize::usage = "Use Nelder Mead to find the minimum value of Log(|R(\!\(\*SubscriptBox[\(r\), \(max\)]\))|) over complex \[Omega]-space
RadialMinimize[\[Omega]Guess_?NumberQ, \[Nu]Guess_?NumberQ, OmegaBoundary_?ListQ,
     parameters_?AssociationQ, Watcher_?BooleanQ, \[Nu]fit_:Null] 
"

ProcaDiffOpRad::usage = "Radial differential operator from FKKS Decomposition of Proca equations of motions";

omegaNRNonRel::usage = "non-relativistic approximation to real part of proca frequency"

omegaNINonRel::usage = "non-relativistic approximation to imaginary part of proca frequency"

nuNNonRel::usage = "non-relativistic approximation to angular eigenvalue"

SolveAngularSystem::usage = "Solve angular FKKS equation in matrix form for kernel vectors"

getRadialSolution::usage = "Solve radial equation numerically for input \[Omega], \[Nu], and parameter association"

FrobeniusSystem::usage = "System of equations resulting from expanding radial function near horizon and solving radial equation order by order"

functomin::usage = "input is frequency. Output is value of radial function at upper boundary"

Mmatrix::usage = "Angular equation in matrix form"

Mmatrixdet::usage = "Determinant of angular equation"

getNuValue::usage = "Get non-relativistic approximation to angular eigenvalue"



QValue::usage = "Compute Q value (numerical variable) for given \[Mu] and \[Omega]"

QApprox::usage = "linear approximation to Q-value over \[Mu]-space"

calcEndingX::usage = "calculate the ending value of the radial coordinate given input parameter association"

calcOrbital::usage = "calculate the orbital angular momentum number given input parameter association"

ToPrecision::usage = "Set precision of expression according to parameter association. Usage: expr//ToPrecision[parameters]"


(*
Begin["KWP`"];
*)


ToPrecision[parameterAssoc_?AssociationQ][expr_] := SetPrecision[expr, parameterAssoc["precision"]];


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
    Sqrt[\[Mu]^2 - \[Omega]^2];


calcEndingX[parameters_] :=
    2 * 10 * (parameters["m"] + parameters["n"]) / parameters["\[Mu]Nv"] ^
         2;


calcOrbital[parameters_] :=
    parameters["m"] + parameters["s"];


KerrSurfaceRotationN[\[Chi]_] := (1 / 2) (\[Chi] / rplusN[\[Chi]]);


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
        Derivative[1][s], Derivative[2][s]}];

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
    Block[{NumMatrix,equalto,System,res,bcoeffs = Table[C[i], {i, 0, parameters["KMax"]}], \[Eta] = parameters[
        "\[Eta]"], m = parameters["m"]},
        NumMatrix = \[Nu]N ^ -7 * Mmatrix[parameters["KMax"], \[Eta], m, \[Chi] * M,
             \[Omega]N / M, \[Nu]N / M, \[Mu]N / M] // N;
        equalto = Table[0, {i, 0, parameters["KMax"]}];
        System = (NumMatrix . bcoeffs == equalto) /. {\[Omega]N -> \[Omega]value, \[Nu]N
             -> \[Nu]value, \[Mu]N -> parameters["\[Mu]Nv"], m -> parameters["m"], \[Eta] -> parameters[
            "\[Eta]"], \[Chi] -> parameters["\[Chi]"]};
        res = Flatten @ {C[0] -> 1, Solve[System /. C[0] -> 1, bcoeffs
            [[2 ;; parameters["KMax"] + 1]]]};
        Return[res];
    ];


(*Radial functions*)

qr[\[Nu]_, r_] :=
    1 + \[Nu]^2 * r^2;

\[Sigma][a_, m_, \[Omega]_, \[Nu]_] :=
    a * \[Nu]^2 (m - a * \[Omega]) + \[Omega];

Kr[a_, m_, r_] :=
    -a * m + (a^2 + r^2) \[Omega];


ProcaDiffOpRad = 
	Block[{Diffeq,
		   rplus = (1 + Sqrt[1 -  \[Chi]^2]),
		   rminus = (1 - Sqrt[1 - \[Chi]^2])
		   },
		   With[{\[CapitalDelta] = (r - rplus) (r - rminus), 
		         qr = 1 + \[Nu]^2 * r^2, 
		         \[Sigma]=\[Chi] * \[Nu]^2 (m -\[Chi] * \[Omega]) + \[Omega], 
		         Kr = -\[Chi]* m + (\[Chi]^2 + r^2) \[Omega]
		         },
		         Diffeq = (D[\[CapitalDelta]/qr*D[#,r],r] + (Kr^2/(qr*\[CapitalDelta])+(2-qr)/qr^2*\[Sigma]/\[Nu]-\[Mu]^2/\[Nu]^2)*#)&;
		         Diffeq/.{\[Chi]->Global`\[Chi],m->Global`m,r->Global`r,\[Mu]->\[Mu], \[Nu]->Global`\[Nu], \[Omega]->Global`\[Omega]}
		         ]
		   ];


DiffRadOpx = 
	With[{
	rplusNval = rplusN[\[Chi]] // Rationalize, 
	rminusNval = rminusN[\[Chi]] // Rationalize,
	qrN = qr[\[Nu], rN * M], (*numerical Subscript[q, r]*)
	rsN = rplusN[\[Chi]] + rminusN[\[Chi]]  // Rationalize, (*numerical Subscript[r, s]*)
	\[CapitalDelta]N = (rN - rplusN[\[Chi]]) * (rN - rminusN[\[Chi]]) // Rationalize, (*numerical \[CapitalDelta]*)
	\[Sigma]N =\[Sigma][a, m, \[Omega], \[Nu]] * M //Expand //Rationalize, (*numerical \[Sigma]*)
	KrN = Kr[a, m, rN * M] / M //Rationalize //Simplify (*numerical Subscript[K, r]*)
	},
	Block[{
	P1 = (2 * rN - rsN) / (rN - rminusNval) // Rationalize,
	P2 = -(2 * rN * \[Nu]N^2) / qrN,
	Q1 = KrN^2 / (rN - rminusNval) ^ 2,
	Q2 = qrN / (rN - rminusNval) * ((2 - qrN) / qrN^2 * \[Sigma]N / \[Nu]N - \[Mu]N^2 / \[Nu]N^2) // Rationalize
	},
	
	With[{
	P1x = P1 /. rN -> xN * (rplusNval - rminusNval) + rplusNval,
	P2x = P2 /. rN -> xN * (rplusNval - rminusNval) + rplusNval,
	Q1x = Q1 /. rN -> xN * (rplusNval - rminusNval) + rplusNval,
	Q2x = Q2 /. rN -> xN * (rplusNval - rminusNval) + rplusNval
	},
	(\!\(
\*SubscriptBox[\(\[PartialD]\), \(xN, xN\)]#\) + (P1x/xN+P2x*(rplusNval-rminusNval))\!\(
\*SubscriptBox[\(\[PartialD]\), \(xN\)]#\) + (Q1x/xN^2+Q2x (rplusNval-rminusNval)/xN)#)&
	]
	]
	];
	
DiffRadOpxsplit = 
	DiffRadOpx/.\[Omega]N->\[Omega]Nr + I*\[Omega]Ni;


FrobeniusSeries[x_,\[Kappa]_,c1_,c2_]:=x^(-I*\[Kappa])*(1 + c1*x + c2*x^2);

Options[FrobeniusSystem] = {Recalculate->False};
FrobeniusSystem[\[Kappa]_,c1_,c2_,\[Nu]Nt_,\[Omega]Nt_,parameters_, OptionsPattern[]]:=
	Block[{},
		If[OptionValue[Recalculate],
		
			Block[{expr, exprSeries, exprSystem, fs = FrobeniusSeries[xN,\[Kappa],c1,c2]},
				expr = xN^(I*\[Kappa]+2)*DiffRadOpx@fs;
				exprSeries = Series[expr,{xN,0,2}]==0;
				exprSystem = LogicalExpand[exprSeries];
				Export[$FKKSRoot<>"Expression/FrobeniusSystem.mx", exprSystem];
				Return[exprSystem]
				];,
				
			Evaluate[Import[$FKKSRoot<>"Expressions/FrobeniusSystem.mx"][\[Kappa],c1,c2,\[Omega]Nt, \[Nu]Nt, parameters["m"], parameters["\[Chi]"], parameters["\[Mu]Nv"]]]
			]
		];

getRadialSolution[w_?NumberQ, v_?NumberQ, parameters_]:=
	Block[{
		(*StartingRadius, 
		EndingRadius, 
		AnalyticFrobeniusSystem,
		FullFrobeniusSolution, 
		FrobeniusSolutionResults,
		FrobeniusSolutionSeries,
		FrobeniusSolution, 
		FrobR0, 
		FrobRPrime0, 
		BC, 
		eq, 
		RadialSolutions,*)
		\[Omega]vv = w, 
		\[Nu]vv=v,
		\[Chi]=parameters["\[Chi]"],
		m=parameters["m"],
		\[Mu]N=parameters["\[Mu]Nv"]},
		
		StartingRadius =  parameters["StartingX"];
		EndingRadius =  parameters["EndingX"];
		AnalyticFrobeniusSystem = FrobeniusSystem[\[Kappa],c1,c2,\[Nu]vv,\[Omega]vv,parameters]//Rationalize[#,0]&;
		FullFrobeniusSolution = NSolve[AnalyticFrobeniusSystem, {\[Kappa],c1,c2}, VerifySolutions->True, WorkingPrecision->parameters["precision"]];
		FrobeniusSolutionResults =(FullFrobeniusSolution[[parameters["branch"]]])[[All,2]];
		FrobeniusSolutionSeries = {xN}|->Evaluate@FrobeniusSeries[xN,Sequence@@FrobeniusSolutionResults];
		FrobR0 = Limit[FrobeniusSolutionSeries[xN], xN->parameters["\[Epsilon]"]];
		FrobRPrime0 = Limit[Evaluate@D[FrobeniusSolutionSeries[xN],xN], xN-> parameters["\[Epsilon]"]];
		BC = Thread@{FrobR0, FrobRPrime0}//ToPrecision[parameters];
		eq = {
			((DiffRadOpx@R[xN]/.{a->\[Chi]*M,\[Omega]->\[Omega]N/M,\[Nu]->\[Nu]N/M,\[Mu]->\[Mu]N/M}//Collect[#,M]&//InsertValues[parameters])/.{ \[Nu]N->\[Nu]vv, \[Omega]N->\[Omega]vv})==0, 
			R[StartingRadius]==FrobR0, 
			Derivative[1][R][StartingRadius]==FrobRPrime0
			}// ToPrecision[parameters];
		RadialSolutions = NDSolve[eq,
								 R, 
								 {xN, StartingRadius, EndingRadius},
								 WorkingPrecision->parameters["precision"], 
								 MaxSteps->parameters["RadialIntegratorMaxStep"], 
								 PrecisionGoal->parameters["RadialIntegratorPrecisionGoal"]
								 ]//First;
		Return[R/.RadialSolutions]
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
        thedet = Function[{w, \[Nu]N}, Evaluate @ Mmatrixdet[w, \[Nu]N, parameters]];
        FindRoot[thedet[\[Omega]N, \[Nu]N] // ToPrecision[parameters], {\[Nu]N, 
                ClosestGuess}, WorkingPrecision -> parameters["precision"], MaxIterations
                 -> 1000][[1]][[2]] //Rationalize[#, 0]& //Return
    ]; 


functomin[parameters_?AssociationQ, \[Nu]fitted_:Null][{\[Omega]real_?NumericQ, \[Omega]imag_?NumericQ}] :=
    Block[{\[Omega]value = ToPrecision[parameters][\[Omega]real + I * \[Omega]imag]},
        Block[{
            nuroot =Evaluate @If[TrueQ[\[Nu]fitted == Null],
                                   Evaluate @ getNuValue[\[Omega]value, parameters, ToPrecision[parameters]@nuNNonRel[\[Omega]value,parameters]]
                                   ,
                                   Evaluate @ getNuValue[\[Omega]value, parameters, ToPrecision[parameters]@\[Nu]fitted]
                                 ]
              },
              radsol = getRadialSolution[\[Omega]value, nuroot, parameters];
              retval = radsol[parameters["EndingX"]] // Abs;
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
     Block[{MinimizationConstraints,MinimizationMethod,RMMessengerGenerator,SimplexBase,SimplexNodeOne,SimplexNodeTwo,InitialSimplex,\[Nu]fitted,result,\[Omega]Result,\[Nu]Result},
    Module[{currval = {0, 0}, CurrentIteration = 0},
        RMMessengerGenerator =
            Function[{X, counter},
                "Current Value: " <> ToString[X, InputForm] <> "\n Current Step: "
                     <> ToString[counter, InputForm]
            ];
        RMMessenger = RMMessengerGenerator[currval, CurrentIteration]
            ;
        (*If[Watcher, PrintTemporary@Dynamic[RMMessenger]];*)
        SimplexBase = {\[Omega]Guess // Re, \[Omega]Guess // Im} // ToPrecision[parameters];
        SimplexNodeOne = SimplexBase - {(OmegaBoundary[[2]] - \[Omega]Guess)/parameters["omegaRes"], 0} // Re;
        SimplexNodeTwo = SimplexBase - {0, 10 ^ (Log10[Abs @ Im[\[Omega]Guess]] - 1)};
        InitialSimplex = {SimplexBase, SimplexNodeOne, SimplexNodeTwo};
        If[TrueQ[\[Nu]fit == Null],\[Nu]fitted=Null,\[Nu]fitted = \[Nu]fit[parameters["\[Mu]Nv"]]];
        
        MinimizationConstraints = ToPrecision[parameters][{100*Im[\[Omega]Guess]> \[Psi] > 10^-13}];
        MinimizationMethod = {"NelderMead", 
                              "PostProcess" -> False, 
                              "InitialPoints" -> InitialSimplex,
		                     "Tolerance"->0
                              };
        
        (*Constraint on imaginary part comes from ad hoc lower bound together with rough scaling of imaginary part with mode number*)                      
		With[{Contraints = {100*Im[\[Omega]Guess]>\[Psi]>10^((-12)*(parameters["m"]+1)/2)&&parameters["\[Mu]Nv"]^2>\[Xi]^2+\[Psi]^2&&Re[omegaBoundary[[2]]]>\[Xi]>Re[omegaBoundary[[1]]]}//ToPrecision[parameters],
				method = MinimizationMethod
		},
		
        result = 
                NMinimize[
                    {Hold @ functomin[parameters, \[Nu]fitted][{\[Xi], \[Psi]}], 
                    Sequence@@Contraints
                    },
                    {\[Xi], \[Psi]},
                    Method -> method,
                    WorkingPrecision -> parameters["OptimizationAccuracy"],
                    EvaluationMonitor :>(RMMessenger = RMMessengerGenerator[\[Xi] + I * \[Psi], CurrentIteration];
                                        CurrentIteration++;
                                        ),
                    MaxIterations -> parameters["MaxIterationsMinimization"]
                ][[2]];
        ];
        \[Omega]Result = (\[Xi] + I * \[Psi]) /. result;
        If[TrueQ[\[Nu]fit == Null],
            \[Nu]Result = getNuValue[\[Omega]Result, parameters, nuNNonRel[\[Omega]Result,
                 parameters]]
            ,
            \[Nu]Result = getNuValue[\[Omega]Result, parameters, \[Nu]fitted];
        ];
    ];
    
    Return[<|"\[Omega]" -> \[Omega]Result, "\[Nu]" -> \[Nu]Result|>];
    ];


(*
End[];
EndPackage[];
*)
