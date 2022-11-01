(* ::Package:: *)

BeginPackage["QNMnonRelsearch`"]

QNMnonRel

Begin["QNMcode`"]


Options[QNMnonRel]={ModePrecision->20};
QNMnonRel[min_,nin_,spin_,\[Mu]sta_,\[Mu]sto_,branch1_,\[Omega]ImIn_,Sin_,direc_, OptionsPattern[]]:=Module[{min1=min,nin1=nin,spin1=spin,\[Mu]sta1=\[Mu]sta,\[Mu]sto1=\[Mu]sto,branch2=branch1,wijk=\[Omega]ImIn,Sin1=Sin},

Print["Save directory: "<>ToString[Global`destinationPATH]];

$Assumptions=M>0&&x!=0;
kmax=5;
xmax=2;
\[Epsilon]=10^-6;
xstart=(\[Epsilon] M)/(rp[a]-rm[a])//Simplify//Rationalize[#,0]&;
\[Omega]nmax=(m \[Chi])/(2(1+Sqrt[1-\[Chi]^2]))//Simplify;
setprec=24;
maxstep=10^8;
modeprec=OptionValue[ModePrecision];

m=min;
\[Eta]=If[Sin==0,1,0];
nh=nin;
S=Sin;
n=Abs[m]+nh+S+1;
\[Chi]=spin//Rationalize[#,0]&;

Print["We consider the following parameters:"];
Print["m = "<>ToString[m]];
Print["n = "<>ToString[nh]];
Print["\[Chi] = "<>ToString[\[Chi]//N]];
Print["S = "<>ToString[Sin//N]];
Print["branch = "<>ToString[branch1]];
Print["Log of precision cutoff = 10^-"<>ToString[modeprec]];
Print["-------------------------------------"];

a=\[Chi] M;\[Mu]=\[Mu]n/M ;\[Omega]=(\[Omega]nr+I \[Omega]ni)/M;\[Nu]=\[Nu]n/M;
rp[a_]:=M+Sqrt[M^2-a^2];
rm[a_]:=M-Sqrt[M^2-a^2];
rstop[\[Mu]_]:=2 (10(m+nh))/\[Mu]^2//Rationalize[#,0]&;
rtstart[a_]:=(xstart(rp[a]-rm[a])+rp[a])/M//Simplify//Rationalize[#,0]&;
\[CapitalOmega]H[a_]:=a/(2M rp[a]);
\[CapitalDelta][r_,a_]:=r^2-2M r+a^2;
Kr[r_,a_,\[Omega]_,m_]:=(a^2+r^2)\[Omega]-a m;
qr[r_,\[Nu]_]:=1+\[Nu]^2 r^2;
\[CapitalLambda][\[Nu]_,\[Mu]_,a_,m_,\[Omega]_]:=\[Mu]^2/\[Nu]^2-(\[Omega]+a \[Nu]^2 (m-a \[Omega]))/\[Nu]+2a \[Omega] m-a^2 \[Omega]^2;
\[Gamma][\[Omega]_,\[Mu]_]:=Sqrt[\[Omega]^2-\[Mu]^2];
\[Sigma][\[Omega]_,a_,\[Nu]_,m_]:=\[Omega]+a \[Nu]^2 (m-a \[Omega]);
Three[l1_,l2_,l3_,m1_]:=If[l3>l1+l2||Abs[l1-l2]>l3 ||Abs[m1]>l1||Abs[m1]>l3,0,ThreeJSymbol[{l1,-m1},{l2,0},{l3,m1}]];
Bracket[l1_,l2_,l3_,m_]:=(-1)^m Sqrt[((2l1+1)(2l2+1)(2l3+1))/(4\[Pi])]Three[l1,l2,l3,0]Three[l1,l2,l3,m];
c2[l_,m_,lp_]:=(2Sqrt[\[Pi]])/3 Bracket[l,0,lp,m]+4/3 Sqrt[\[Pi]/5]Bracket[l,2,lp,m];
c4[l_,m_,lp_]:=(2Sqrt[\[Pi]])/5 Bracket[l,0,lp,m]+8/7 Sqrt[\[Pi]/5]Bracket[l,2,lp,m]+(16Sqrt[\[Pi]])/105 Bracket[l,4,lp,m];
d2[l_,m_,lp_]:=Sqrt[(4\[Pi])/3](lp Sqrt[((lp+1)^2-m^2)/((2lp+1)(2lp+3))]Bracket[l,1,lp+1,m]-(lp+1)Sqrt[(lp^2-m^2)/((2lp+1)(2lp-1))]Bracket[l,1,lp-1,m]);
Mtemp[l_,m_,lp_,\[Nu]_,\[Mu]_,a_,\[Omega]_]:=(\[CapitalLambda][\[Nu],\[Mu],a,m,\[Omega]]-lp(lp+1))KroneckerDelta[l,lp]+(-\[Nu]^2\[CapitalLambda][\[Nu],\[Mu],a,m,\[Omega]]+\[Nu]^2 lp(lp+1)-2\[Sigma][\[Omega],a,\[Nu],m]\[Nu]+\[Gamma][\[Omega],\[Mu]]^2)a^2 c2[l,m,lp]-2a^2 \[Nu]^2 d2[l,m,lp]-\[Gamma][\[Omega],\[Mu]]^2 \[Nu]^2 a^4 c4[l,m,lp]//Refine;
Mat[m_,\[Eta]_,\[Nu]_,\[Mu]_,a_,\[Omega]_]:=Table[Mtemp[Abs[m]+2k+\[Eta],m,Abs[m]+2kp+\[Eta],\[Nu],\[Mu],a,\[Omega]],{k,0,kmax},{kp,0,kmax}]//Simplify;
\[Kappa][a_,\[Omega]_,m_]:=(2 M rp[a](\[Omega]-m \[CapitalOmega]H[a]))/(rp[a]-rm[a])//Simplify;
(*Clear[\[Chi],\[Omega]ni,\[Omega]nr,\[Nu]n,\[Mu]n,m,\[Eta],a,\[Omega],\[Mu],\[Nu]]
DiffRadr:=D[\[CapitalDelta][r,a]D[#,r],r]+(Kr[r,a,\[Omega],m]^2/\[CapitalDelta][r,a]-\[CapitalLambda][\[Nu],\[Mu],a,m,\[Omega]]+2 a \[Omega] m-a^2\[Omega]^2-\[Mu]^2r^2)#-(2r \[Nu]^2)/qr[r,\[Nu]](\[CapitalDelta][r,a]D[#,r]+r \[Sigma][\[Omega],a,\[Nu],m]/\[Nu]#)&;
Collect[DiffRadr@Y[r]/.{Y\[Rule](Y[#/M]&)}/.{r\[Rule]rt M}//Simplify,{Y[rt],Y'[rt],Y''[rt]}]
Collect[(DiffRadr@f[r]/.{f\[Rule](f[(#-rp[a])/(rp[a]-rm[a])]&)}/.{r\[Rule]x(rp[a]-rm[a])+rp[a]})/(x(1+x))//Simplify,{f[x],f'[x],f''[x]}]*)
DiffRadrt=(-((2 rt^2 \[Nu]n (m \[Nu]n^2 \[Chi]-I (-1+\[Nu]n^2 \[Chi]^2) (\[Omega]ni-I \[Omega]nr)))/((1+rt^2 \[Nu]n^2) (-2 rt+rt^2+\[Chi]^2)))+(-rt^2 \[Mu]n^2-\[Mu]n^2/\[Nu]n^2+(m \[Chi]-I (rt^2+\[Chi]^2) (\[Omega]ni-I \[Omega]nr))^2/(-2 rt+rt^2+\[Chi]^2)+(I \[Omega]ni+\[Omega]nr)/\[Nu]n+\[Nu]n \[Chi] (m-\[Chi] (I \[Omega]ni+\[Omega]nr)))/(-2 rt+rt^2+\[Chi]^2)) #+(-((2 rt \[Nu]n^2)/(1+rt^2 \[Nu]n^2))+(2 (-1+rt))/(-2 rt+rt^2+\[Chi]^2)) D[#,rt]+D[#,{rt,2}]&;
DiffRadx=1/(x (1+x)) (-(\[Mu]n^2/\[Nu]n^2)-\[Mu]n^2 (1+Sqrt[1-\[Chi]^2]+2 x Sqrt[1-\[Chi]^2])^2-(2 M \[Nu]n (1+Sqrt[1-\[Chi]^2]+2 x Sqrt[1-\[Chi]^2]) (1-\[Chi]^2+Sqrt[1-\[Chi]^2]-2 x (-1+\[Chi]^2)) (m \[Nu]n^2 \[Chi]-I (-1+\[Nu]n^2 \[Chi]^2) (\[Omega]ni-I \[Omega]nr)))/(Sqrt[1-\[Chi]^2] (M+M \[Nu]n^2 (2-\[Chi]^2+2 Sqrt[1-\[Chi]^2]-4 x^2 (-1+\[Chi]^2)+4 x (1-\[Chi]^2+Sqrt[1-\[Chi]^2]))))+(I m \[Chi]+2 (1+Sqrt[1-\[Chi]^2]-2 x^2 (-1+\[Chi]^2)+2 x (1-\[Chi]^2+Sqrt[1-\[Chi]^2])) (\[Omega]ni-I \[Omega]nr))^2/(4 x (1+x) (-1+\[Chi]^2))+(I \[Omega]ni+\[Omega]nr)/\[Nu]n+\[Nu]n \[Chi] (m-\[Chi] (I \[Omega]ni+\[Omega]nr))) #+(1+2 x+(4 M x (1+x) \[Nu]n^2 (-1+\[Chi]^2) (1+Sqrt[1-\[Chi]^2]+2 x Sqrt[1-\[Chi]^2]))/(Sqrt[1-\[Chi]^2] (M+M \[Nu]n^2 (2-\[Chi]^2+2 Sqrt[1-\[Chi]^2]-4 x^2 (-1+\[Chi]^2)+4 x (1-\[Chi]^2+Sqrt[1-\[Chi]^2]))))) /(x (1+x)) D[#,x]+D[#,{x,2}]&;

prec=SetPrecision[#,setprec]&;
(*The angular equation in matrix form*)
mat=prec@(Mat[m,\[Eta],\[Nu],\[Mu],a,\[Omega]]\[Nu]n^2);
matdettemp=prec@Det[mat]; (*With and without the \[Nu]n^12 in front, which is analytically equivalent to not having it there, changes the resulting root starting from the sixth decimal*)
matdet[r\[Nu]_,i\[Nu]_]:=prec@(matdettemp//.{\[Nu]n->r\[Nu]+I i\[Nu]});
absdetmat[r\[Nu]_,i\[Nu]_]:=prec@Log[Abs[matdet[r\[Nu],i\[Nu]]]];

Rmaxfunc[wIn_,\[Nu]In_,\[Mu]In_,branch3_]:=Module[{wIn1=wIn,\[Nu]In1=\[Nu]In,\[Mu]In1=\[Mu]In},
	prec=SetPrecision[#,setprec]&;
	\[Mu]i=prec@\[Mu]In;

	wrpoint=Re[wIn]//Rationalize[#,0]&;
	wipoint=Im[wIn]//Rationalize[#,0]&;

	\[Nu]rpoint=Re[\[Nu]In]//Rationalize[#,0]&;
	\[Nu]ipoint=Im[\[Nu]In]//Rationalize[#,0]&;

	\[Nu]rootSollist=prec@NSolve[SetPrecision[(matdettemp/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint}),setprec]==0,\[Nu]n,WorkingPrecision->setprec];
	\[Nu]rootVallist=prec@Table[{\[Nu]n/.\[Nu]rootSollist[[l]]//Re,\[Nu]n/.\[Nu]rootSollist[[l]]//Im},{l,1,Length[\[Nu]rootSollist]}];
	\[Nu]roottemp=prec@Nearest[\[Nu]rootVallist,{\[Nu]rpoint,\[Nu]ipoint},WorkingPrecision->setprec]//First;
	\[Nu]root=prec@{r\[Nu]->\[Nu]roottemp[[1]],i\[Nu]->\[Nu]roottemp[[2]]};
			
	Frob[x_]:=x^(-I \[Lambda]) (1+Sum[Symbol["QNMcode`c"<>ToString[n]]x^n,{n,1,2}]);
	exptemp=Series[SetPrecision[x^(I \[Lambda] +2) (DiffRadx@Frob[x])//Simplify,setprec],{x,0,2}]==0;
	\[CapitalKappa]=(-m \[Chi]+2 (1+Sqrt[1-\[Chi]^2]) (I \[Omega]ni+\[Omega]nr))/(2 Sqrt[1-\[Chi]^2])/.{\[Omega]nr->wrpoint,\[Omega]ni->wipoint};
	exptemp1=SetPrecision[exptemp/.{\[Mu]n->\[Mu]i}/.{\[Omega]nr->wrpoint,\[Omega]ni->wipoint}/.{\[Nu]n->r\[Nu]+I i\[Nu]}/.\[Nu]root//LogicalExpand,setprec];
	(*setprecNSolve[ju_]:=Piecewise[{{setprec-1,ju\[LessEqual] 22},{setprec-2,ju>22}}];*)
	Coeff=NSolve[exptemp1,{\[Lambda],c1,c2},WorkingPrecision->setprec-2(*setprecNSolve[i]*)][[branch3]];
	\[Lambda]test=Abs[(\[Lambda]/.Coeff)+\[CapitalKappa]];
	(*If[Log[\[Lambda]test]<-15,{},Print[\[Lambda]test]Print["Error: Possible problem with value for \[Lambda] or \[Kappa]!"]Print[Precision[\[CapitalKappa]]]];*)
	FrobOfr[rt_]:=SetPrecision[Frob[x]//.Coeff//.{x->(rt M-rp[a])/(rp[a]-rm[a])}(*/.{\[Lambda]\[Rule]\[CapitalKappa]}*)/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint}//Simplify,setprec];
	R0=SetPrecision[Limit[FrobOfr[rt],rt->rtstart[a]],setprec];
	dR0=SetPrecision[Limit[D[FrobOfr[rt],rt],rt->rtstart[a]],setprec];
		
	solR=NDSolve[{(DiffRadrt@R[rt]/.{\[Nu]n->r\[Nu]+I i\[Nu]}/.\[Nu]root/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint})==0,R[rtstart[a]]==R0,R'[rtstart[a]]==dR0},R,{rt,rtstart[a],rstop[\[Mu]i]},Method->{"StiffnessSwitching"},WorkingPrecision->setprec-2,MaxSteps->maxstep]//First;

	Rmax=SetPrecision[Log[Abs[R[rt]/.solR/.{rt->rstop[\[Mu]i]}]],setprec];
	\[Nu]rootout=prec@(\[Nu]roottemp[[1]]+I \[Nu]roottemp[[2]]);
	\[Omega]out=prec@wIn;
];

(*Regime of consideration*)
\[Mu]start=\[Mu]sta;
\[Mu]stop=\[Mu]sto;
\[Omega]mesh=20;
\[Omega]precstep=2;
\[Omega]rstepsize=400;
\[Mu]mesh=126;
branch=branch1;

(*Non-relativistic limit and its values as first guesses*)
\[Omega]nonRel[\[Mu]n_,na_]:=(\[Mu]n(1-\[Mu]n^2/(2na^2)))//Rationalize[#,0]&;(*Provides only a guess for the real part.*)
If[Sin==-1,
(*This is for S=-1*)
Print["Picked S=-1 \[Nu]nonrel"];
\[Nu]nonRel[\[Omega]a_]:=-\[Omega]a/(1-\[Chi] \[Omega]a)//Rationalize[#,0]&; 
,
If[Sin==0,
(*This is for S=0*)
Print["Picked S=0 \[Nu]nonrel"];
\[Nu]nonRel[\[Omega]a_]:=1/(2\[Chi])(2-\[Chi] \[Omega]a+Sqrt[(-2+\[Chi] \[Omega]a)^2+4\[Chi] \[Omega]a]);
,
If[Sin==1, (*This is for S=+1*)
Print["Picked S=+1 \[Nu]nonrel"];
\[Nu]nonRel[\[Omega]b_]:=v//.(NSolve[a v^3(1-a \[Omega]a)-(6-a \[Omega]a(2-a \[Omega]a))v^2+\[Omega]a v+\[Omega]a^2,v][[2]]//Simplify)/.{M->1}//.{\[Omega]a->\[Omega]b};
,
Print["Error: S not in (-1,0,+1)!"];
Abort[];
]]];

(*For a given complex (!) \[Omega],this provides an initial guess for real and imaginary parts*)
\[Mu]range={0,Range[\[Mu]start,\[Mu]stop,(\[Mu]stop-\[Mu]start)/\[Mu]mesh]}//Flatten//Rationalize[#,0]&;
\[Mu]all=Take[\[Mu]range,-Length[\[Mu]range]+1];
\[Mu]rangenumber=Length[\[Mu]range]-1;

(*Initial \[Omega] guesses*)
\[Omega]intR=\[Omega]nonRel[\[Mu]start,n];
\[Omega]intI=\[Omega]ImIn;
\[Nu]int=SetPrecision[{r\[Nu]->Re[\[Nu]nonRel[\[Omega]intR+I \[Omega]intI]],i\[Nu]->Im[\[Nu]nonRel[\[Omega]intR+I \[Omega]intI]]},setprec];

(*The function to be minimized*)
funcmin[win_,\[Nu]in_,\[Mu]in_,branch_]:=Module[{win1=win,\[Nu]in1=\[Nu]in,\[Mu]in1=\[Mu]in},
Rmaxfunc[win,\[Nu]in,\[Mu]in,branch];
Rmax
];

(*Isomorphism from R^2 to C (and inverse)*)
R2toC[vec_]:=vec[[1]]+I vec[[2]];
CtoR2[num_]:={Re[num],Im[num]};

\[Nu]rootminall=ConstantArray[0,Length[\[Mu]range]];
solRminall=ConstantArray[0,Length[\[Mu]range]];
\[Omega]minall=ConstantArray[0,Length[\[Mu]range]];
\[Mu]minall=ConstantArray[0,Length[\[Mu]range]];
plotlist=ConstantArray[0,{Length[\[Mu]range],\[Omega]precstep}];

\[Nu]rootminout=ConstantArray[0,Length[\[Mu]range]];
solRminout=ConstantArray[0,Length[\[Mu]range]];
\[Omega]minout=ConstantArray[0,Length[\[Mu]range]];
\[Mu]minout=ConstantArray[0,Length[\[Mu]range]];
lminout=ConstantArray[0,Length[\[Mu]range]];

(*The only reliable value for \[Omega] and \[Nu] is starting at [[2]], since in the first step, we use the non-rel. result without (!) optimization and minimization to our equations!!! Always remove the first value from the list*)
\[Nu]rootminall[[1]]=\[Nu]int;
\[Omega]minall[[1]]=SetPrecision[\[Omega]intR+I \[Omega]intI,setprec];

(*The minimum search mesh boundaries for the real part of \[Omega]*)
rrinitial=1/2 (\[Omega]nonRel[\[Mu]start,n+1]-\[Omega]nonRel[\[Mu]start,n]);
rrdelta[u_]:=If[u==1,rrinitial,1/2 (\[Omega]nonRel[\[Mu]range[[u+1]],n+1]-\[Omega]nonRel[\[Mu]range[[u+1]],n])];

(*The minimum search mesh boundaries for the imaginary part of \[Omega]*)
iiThreshhold=10^-1;
iiinitial=1;
iidelta[u_]:=SetPrecision[If[u<4,iiinitial,If[Abs[Log10[Im[\[Omega]minall[[u]]]]-Log10[Im[\[Omega]minall[[u-1]]]]]>iiThreshhold,Abs[Log10[Im[\[Omega]minall[[u]]]]-Log10[Im[\[Omega]minall[[u-1]]]]],0.5]],setprec];

(*The code:*)
offset=0;
Clear[y];
For[i=1,i<=\[Mu]rangenumber,i++,
branchtemp=branch;
If[i<= 2,{},
If[Abs[Log10[Im[\[Omega]minall[[i]]]]-Log10[Im[\[Omega]minall[[i-1]]]]]>2,
branch=Select[{1,2},#!=branchtemp&]//First; 
\[Omega]minall[[i]]=\[Omega]minalltemp;
i=i-1;
Print["flip!"];
,
branch=branchtemp]];

Print[ToString[i]<>ToString["/"]<>ToString[\[Mu]mesh]];
\[Mu]i=\[Mu]range[[i+1]];
Print[\[Mu]i//N];

(*We always take the previous \[Omega] as the new guess, since these are much closer to the actual value than the non-rel limit. In the very first step, we use the non-rel limit as guess*)
wintr=SetPrecision[\[Omega]nonRel[\[Mu]i,n],setprec];
winti=SetPrecision[Im[\[Omega]minall[[i]]],setprec];
\[Nu]root=\[Nu]rootminall[[i]];

(*The initial iteration boaundries*)
rrboundary=rrdelta[i];
iiboundary=iidelta[i];     

(*For the very first iteration, we choose this asymmetric search interval, since the real part of the frequency is monotonically increasing with \[Mu] (if not rescaled by \[Mu]^-1)*)
(*\[Omega]rrange=SetPrecision[If[wintr+rrboundary<\[Omega]nmax,Range[wintr-rrboundary,wintr+rrboundary,(2rrboundary)/\[Omega]mesh],Range[wintr-rrboundary,\[Omega]nmax- 10^-6,Abs[wintr-rrboundary-(\[Omega]nmax- 10^-6)]/\[Omega]mesh]],setprec];*)
\[Omega]irange=SetPrecision[10^Range[Log10[winti]-iiboundary,Log10[winti]+iiboundary,(2iiboundary)/\[Omega]mesh],setprec];
(*\[Omega]rrangeplot[i]=Show[ListPlot[Table[{i,\[Omega]rrange[[k]]},{k,1,Length[\[Omega]rrange]}]],ListPlot[{{i,\[Omega]minall[[i]]//Re}},PlotStyle\[Rule]Red]];*)
\[Omega]irangeplot[i]=Show[ListLogPlot[Table[{i,\[Omega]irange[[k]]},{k,1,Length[\[Omega]irange]}]],ListLogPlot[{{i,\[Omega]minall[[i]]//Im}},PlotStyle->Red]];

	\[Delta]\[Omega]r=0;
	aa=1;
	aatest=-1;
	wrpoint=wintr;
		While[aatest<0,
			AAtest=maxR[2]-maxR[1];
			wrpoint=SetPrecision[If[aa<=2,wintr-aa \[Delta]\[Omega]r+offset,If[AAtest<0,wintr-aa \[Delta]\[Omega]r+offset,wintr+aa \[Delta]\[Omega]r+offset]],setprec];
			wipoint=winti;

			(*To pick the correct \[Nu] solution form the set of roots of the EVP, we compare it to either the non-rel. limit using the current \[Omega], or the previous (previous mass) minimal result for \[Nu]*)
			\[Nu]rpoint=prec@(r\[Nu]/.\[Nu]root);
			\[Nu]ipoint=prec@(i\[Nu]/.\[Nu]root);

			\[Nu]rootSollist=NSolve[SetPrecision[(matdettemp/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint}),setprec]==0,\[Nu]n,WorkingPrecision->setprec];
			\[Nu]rootVallist=SetPrecision[Table[{\[Nu]n/.\[Nu]rootSollist[[l]]//Re,\[Nu]n/.\[Nu]rootSollist[[l]]//Im},{l,1,Length[\[Nu]rootSollist]}],setprec];
			\[Nu]roottemp=Nearest[\[Nu]rootVallist,{\[Nu]rpoint,\[Nu]ipoint},WorkingPrecision->setprec]//First;
			\[Nu]root=SetPrecision[{r\[Nu]->\[Nu]roottemp[[1]],i\[Nu]->\[Nu]roottemp[[2]]},setprec];
			
			Frob[x_]:=x^(-I \[Lambda]) (1+Sum[Symbol["QNMcode`c"<>ToString[n]]x^n,{n,1,2}]);
			exptemp=Series[SetPrecision[x^(I \[Lambda] +2) (DiffRadx@Frob[x])//Simplify,setprec],{x,0,2}]==0;
			\[CapitalKappa]=(-m \[Chi]+2 (1+Sqrt[1-\[Chi]^2]) (I \[Omega]ni+\[Omega]nr))/(2 Sqrt[1-\[Chi]^2])/.{\[Omega]nr->wrpoint,\[Omega]ni->wipoint};
			exptemp1=SetPrecision[exptemp/.{\[Mu]n->\[Mu]i}/.{\[Omega]nr->wrpoint,\[Omega]ni->wipoint}/.{\[Nu]n->r\[Nu]+I i\[Nu]}/.\[Nu]root//LogicalExpand,setprec];
			(*setprecNSolve[ju_]:=Piecewise[{{setprec-1,ju\[LessEqual] 22},{setprec-2,ju>22}}];*)
			Coeff=NSolve[exptemp1,{\[Lambda],c1,c2},WorkingPrecision->setprec-2(*setprecNSolve[i]*)][[branch]];
			\[Lambda]test=Abs[(\[Lambda]/.Coeff)+\[CapitalKappa]];
			(*If[Log[\[Lambda]test]<-15,{},Print[\[Lambda]test]Print["Error: Possible problem with value for \[Lambda] or \[Kappa]!"]Print[Precision[\[CapitalKappa]]]];*)
			FrobOfr[rt_]:=SetPrecision[Frob[x]//.Coeff//.{x->(rt M-rp[a])/(rp[a]-rm[a])}(*/.{\[Lambda]\[Rule]\[CapitalKappa]}*)/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint}//Simplify,setprec];
			R0=SetPrecision[Limit[FrobOfr[rt],rt->rtstart[a]],setprec];
			dR0=SetPrecision[Limit[D[FrobOfr[rt],rt],rt->rtstart[a]],setprec];
		
			solR=NDSolve[{(DiffRadrt@R[rt]/.{\[Nu]n->r\[Nu]+I i\[Nu]}/.\[Nu]root/.{\[Mu]n->\[Mu]i,\[Omega]nr->wrpoint,\[Omega]ni->wipoint})==0,R[rtstart[a]]==R0,R'[rtstart[a]]==dR0},R,{rt,rtstart[a],rstop[\[Mu]i]},Method->{"StiffnessSwitching"},WorkingPrecision->setprec-2,MaxSteps->maxstep]//First;

		    solRIterationR[aa]=solR;
			\[Nu]rootIterationR[aa]=\[Nu]root;
			\[Omega]IterationR[aa]=SetPrecision[wrpoint,20];

			maxR[aa]=SetPrecision[Log[Abs[R[rt]/.solRIterationR[aa]/.{rt->rstop[\[Mu]i]}]],setprec];
			aatest=If[aa<=2,-1,maxR[aa]-maxR[aa-1]];
			\[Delta]\[Omega]r=SetPrecision[rrdelta[i]/\[Omega]rstepsize,20];
			If[aa>500,Print["aa>500: Terminated!"];Break[],{}];
		;aa++];
		If[AAtest<0,Print["-1"],Print["+1"]];
		Print[aa-1];
		\[Omega]minR=\[Omega]IterationR[aa-1];
		If[aa>50&&AAtest<0,offset=1(\[Omega]minR-\[Omega]nonRel[\[Mu]i,n]),{}];
		If[aa>100&&AAtest<0,offset=2(\[Omega]minR-\[Omega]nonRel[\[Mu]i,n]),{}];
		\[Nu]transfer=\[Nu]roottemp[[1]]+I \[Nu]roottemp[[2]];
		\[Omega]transfer=\[Omega]minR+I winti;
(*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*)
(*Initial simplex dimensions*)
prec=SetPrecision[#,modeprec]&;

Print["Simplex code..."];
rrinitial2=\[Delta]\[Omega]r/10;
iiinitial2=10^(Log10[Im[\[Omega]minall[[i]]]]-1);
\[Nu]int=prec@\[Nu]transfer;
wint=prec@\[Omega]transfer;
wvec0=prec@{Re[wint],Im[wint]};
wvec1=prec@{Re[wint]-rrinitial2,Im[wint]};
wvec2=prec@{Re[wint],Im[wint]-iiinitial2};
newsimplex=prec@{wvec0,wvec1,wvec2};

Print["Initial \[Omega] guess: "<>ToString[wint, InputForm]];
testnorm=1;
count=1;
While[
(*termination test:*)
testnorm>10^(Log10[Im[wint]]-modeprec),

(*Simplex*)
funcofsimplex=Table[funcmin[R2toC[newsimplex[[j]]],\[Nu]int,\[Mu]i,branch],{j,1,3}];
fh=prec@Max[funcofsimplex];
posfh=Position[funcofsimplex,fh][[1]];
xh=prec@newsimplex[[posfh]]//First;
fl=Min[funcofsimplex];
posfl=Position[funcofsimplex,fl][[1]];
xl=prec@newsimplex[[posfl]]//First;
fs=prec@Delete[funcofsimplex,{posfh,posfl}]//First;
posfs=Position[funcofsimplex,fs][[1]];
xs=prec@newsimplex[[posfs]]//First;
centroid=prec@1/2(xl+xs);

(*Transformation of simplex*)
\[Alpha]Simplex=1;
\[Beta]Simplex=1/2;
\[Gamma]Simplex=2;
\[Delta]Simplex=1/2;
(*Reflection point*)
xr=prec@(centroid+\[Alpha]Simplex(centroid-xh));
fr=prec@funcmin[R2toC[xr],\[Nu]int,\[Mu]i,branch];

(*Body of the method*)
If[fl<= fr&&fr<fs,
(*Print["Reflected 1"];*)
xnew=xr;
x1=xs;
x2=xl;
,
If[fr<fl,
(*expansion*)
xe=prec@(centroid+\[Gamma]Simplex(xr-centroid));
fe=prec@funcmin[R2toC[xe],\[Nu]int,\[Mu]i,branch];
If[fe<fr,
(*Print["Expanded"];*)
xnew=xe;
x1=xs;
x2=xl;
,
(*Print["Reflected 2"];*)
xnew=xr;
x1=xs;
x2=xl;
];
,
(*contraction*)
If[fr>= fs,
(*outside*)
If[fr<fh,
xc=prec@(centroid+\[Beta]Simplex(xr-centroid));
fc=prec@funcmin[R2toC[xc],\[Nu]int,\[Mu]i,branch];
If[fc<= fr,
(*Print["Contracted 1"];*)
xnew=xc;
x1=xs;
x2=xl;
,
(*Print["Shrunk 1"];*)
xnew=prec@(xl+\[Delta]Simplex(xs-xl));
x1=prec@(xl+\[Delta]Simplex(xh-xl));
x2=xl;
];
,
(*inside*)
xc=prec@(centroid+\[Beta]Simplex(xh-centroid));
fc=prec@funcmin[R2toC[xc],\[Nu]int,\[Mu]i,branch];
If[fc<fh,
(*Print["Contracted 2")*)
xnew=xc;
x1=xs;
x2=xl;
,
(*Print["Shrunk 2"];*)
xnew=prec@(xl+\[Delta]Simplex(xs-xl));
x1=prec@(xl+\[Delta]Simplex(xh-xl));
x2=xl;
];
,Print["Error1!"]];

,Print["Error2!"]];
];
];

(*new simplex*)
newsimplex=prec@({xnew,x1,x2}//Abs);
xnewprint=R2toC[xnew];
simplexset[i,count]=newsimplex;
\[Nu]int=prec@\[Nu]rootout;
testnorm=prec@Max[{Norm[xnew-x1],Norm[xnew-x2],Norm[x1-x2]}];
testnormset[i,count]=testnorm;
count=count+1;
If[count>5000,Print["Terminated loop 3"];Break[],{}];
];
	\[Omega]minalltemp=\[Omega]minall[[i]];
	funcmin[prec@R2toC[1/3(xnew+x1+x2)],\[Nu]int,\[Mu]i,branch];
	\[Mu]minall[[i+1]]=\[Mu]i;
	solRminall[[i+1]]=solR;
	\[Nu]rootminall[[i+1]]=prec@{r\[Nu]->Re[\[Nu]rootout],i\[Nu]->Im[\[Nu]rootout]};
	\[Omega]minall[[i+1]]=\[Omega]out;
	lminout[[i]]=\[Lambda]/.Coeff;
	printw=prec@\[Omega]minall[[i+1]];
	Print["\[Omega] value: "<>ToString[printw, InputForm]];
	Print["---------------------------------------------"];
(*If[count>5000,Print["Terminated 3!"];Break[],{}];*)
If[aa>500,Print["Terminated 2!"];Break[],{}];
Print["Iteration count: "<>ToString@count];

solRminout[[i]]=solRminall[[i+1]];
\[Omega]minout[[i]]=\[Omega]minall[[i+1]];
\[Nu]rootminout[[i]]=\[Nu]rootminall[[i+1]];
\[Mu]minout[[i]]=\[Mu]i;

Print["Solving angular equation..."];
b=Table[Symbol["QNMcode`b"<>ToString[i]],{i,0,kmax}];
AnglCoeffList=ConstantArray[0,{Length[\[Mu]minout]}];
AnglFuncList=ConstantArray[0,{Length[\[Mu]minout]}];

(*The b0 will parameterize all the solutions for \vec{b} in the kernel of mat. With pick a b0, such that the resulting angular solution has a global maximum of O(1); for numerical convenience.*)
Y[l_,m_,\[Theta]_]:=SphericalHarmonicY[l,m,\[Theta],\[Phi]]Exp[-I m \[Phi]]//Simplify;
Sfunc[m_,\[Theta]_,\[Eta]_]:=Sum[Symbol["QNMcode`b"<>ToString[kp]]Y[Abs[m]+2kp+\[Eta],m,\[Theta]],{kp,0,kmax}];

Do[(*Note the multiplication by some power of \[Nu]n in order to get the components of the matrix to be of order 1, rather than e-10*)
b0norm=10;
matplug=prec@(mat//Simplify)//.{\[Nu]n->r\[Nu]+I i\[Nu]}/.\[Nu]rootminout[[n]]/.{\[Mu]n->\[Mu]minout[[n]],\[Omega]nr->Re[\[Omega]minout[[n]]],\[Omega]ni->Im[\[Omega]minout[[n]]]}//Simplify;
mattemp=prec@matplug . b;
linsys=Table[mattemp[[l]]==0,{l,1,kmax+1}]//.{b0->b0norm};
solvar=Table[Symbol["QNMcode`b"<>ToString[i]],{i,0,kmax}];
AnglCoeff=Solve[linsys,solvar]//.{b0->b0norm}//Flatten;
AnglCoeffList[[n]]=AnglCoeff;
Splot[\[Theta]_]:=Sfunc[m,\[Theta],\[Eta]]//.AnglCoeff//.{b0->b0norm}//Simplify;
AnglFuncList[[n]]=Splot[\[Theta]];
,{n,1,i}];
Print["Done"];

modedataoutput={AnglCoeffList,AnglFuncList,\[Mu]minout,\[Nu]rootminout,\[Omega]minout,solRminout,lminout};
spinstring=NumberForm[spin*10^6//Round,6,DigitBlock->5,ExponentStep->6,NumberSeparator->""];

OutputFile =ToString[ToString[Global`destinationPATH]<>"/m"<>ToString[m]<>ToString["n"]<>ToString[nh]<>ToString["_a"]<>ToString[spinstring]<>ToString["_S"]<>ToString[If[Sin<0,"m","p"]]<>ToString[Abs[Sin]]<>ToString["_prec_"]<>ToString[OptionValue[ModePrecision]]<>ToString["_HPee.mx"]];
Print["Exporting data to "<>OutputFile];
Export[OutputFile,modedataoutput];

ClearSystemCache["Numerical"];

];
Print["Minimization: Done!"];
];

End[]

EndPackage[]



