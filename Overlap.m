(* ::Package:: *)

(*
Calculating the an overlap of a pair of 
waveforms (http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.084012)

Maxim Zalutskiy, 2015
*)
BeginPackage["Overlap`"]
ReadKerrQNM::usage="ReadKerrQNM[l, m, n] loads quasi-normal modes for given l, m and n.";

Coefficients::usage="Coefficients[data,lm,modes,negmodes,time_)points,"<>
"FindMaximum_options].\n Produces the output {Maximum of \[Rho]2, {\[Delta],a,\[Theta]}, {{A,\[Phi]},...},"<>
"{{standard_error_in_A, standard_error_in_\[Phi]},...}}";

\[Rho]2::usage="\[Rho]2[data,lm,modes,pmodes,\[Delta],a,\[Theta],it0].\nCalculates the overlap of two waveforms.";
OverlapPlot::usage = "OverlapPlot[data,coeffs,lm].\n Plots the overlap."
MCBootStrap::usage = "MCBootStrap[data,lm,modes,pmodes,Nt,Nstat].";
MCBootStrapCError::usage = "";
MaxOverlap::usage = "";
GetTime::usage = "";
Protect[Mass,Spin,Theta,CoeffsDataRand,\[Rho]2DataRand,GTDataRand];

Begin["Private`"]


(* Quasi-normal modes loader*)
ReadKerrQNM[l_,m_,n_]:=Module[{nname,mnamea,mname,lname,rawdat,Nelems,NL,
	Lmin,Lmax,i,a,\[Omega],Alm,Allmn},
	nname=If[n<10,"0"<>ToString[n],ToString[n]];
	mnamea=If[Abs[m]<10,"0"<>ToString[Abs[m]],ToString[Abs[m]]];
	mname=If[m<0,"-"<>mnamea,"+"<>mnamea];
	lname = If[l<10,"0"<>ToString[l],ToString[l]];
	Import[path<>"KerrQNM_"<>nname<>".h5",{"HDF5","Datasets",{"/n"<>
		nname<>"/m"<>mname<>"/L"<>lname<>".dat"}}]
]


DefineInterpolations[lm_,modes_,neg_]:=
	Module[{modesize,lmsize,i,ll,mm,nn,QNMdata,QNMsize,lmin,Need\[Omega],
	Re\[Omega],Im\[Omega],spin,k,l,Need\[ScriptCapitalA]},
	
	modesize = Length[modes];
	lmsize = Length[lm];
	For[i=1,i<=modesize,++i,
		ll=modes[[i,1]];
		mm=neg*modes[[i,2]];
		nn=modes[[i,3]];
		QNMdata = ReadKerrQNM[ll,mm,nn];
		QNMsize = Length[QNMdata];
		lmin = Max[2,Abs[mm]];
		Need\[Omega] = If[Head[\[Omega]bar[ll,mm,nn]]==InterpolatingFunction,False,True,True];
		If[Need\[Omega],
			Re\[Omega][i_]:=QNMdata[[i,2]];
			Im\[Omega][i_]:=QNMdata[[i,3]];
			spin[i_]:=QNMdata[[i,1]];
			\[Omega]bar[ll,mm,nn] = 
				Interpolation[Table[{spin[j],Re\[Omega][j]+I*Im\[Omega][j]},{j,1,QNMsize}]];
		]
		For[k=1,k<=lmsize,++k,
			l=lm[[k,1]];
			Need\[ScriptCapitalA]=If[Head[\[ScriptCapitalA][l,ll,mm,nn]]==InterpolatingFunction,False,True,
											True];				
			If[Need\[ScriptCapitalA],
				\[ScriptCapitalA][l,ll,mm,nn]=Interpolation[
					Table[{QNMdata[[j,1]],QNMdata[[j,2(l-lmin)+6]]
							+QNMdata[[j,2(l-lmin)+7]]I},{j,1,QNMsize}]];
			];
		];
	];
]


\[Psi]pos[lm_, mode_, \[Delta]_, a_, \[Theta]_, t_] := WignerD[{mode[[1]],-lm[[2]],-mode[[2]]},\[Theta]]*
				\[ScriptCapitalA][lm[[1]], mode[[1]], mode[[2]], mode[[3]]][a]*
				Exp[-I*\[Omega]bar[mode[[1]], mode[[2]], mode[[3]]][a]*t/\[Delta]]


\[Psi]neg[lm_, pmode_, \[Delta]_, a_, \[Theta]_, t_] := (-1)^(lm[[1]]+pmode[[1]])*
			   WignerD[{pmode[[1]],-lm[[2]],-pmode[[2]]},\[Theta]]*
			   Conjugate[\[ScriptCapitalA][lm[[1]], pmode[[1]], -pmode[[2]], pmode[[3]]][a]*
			   Exp[-I*\[Omega]bar[pmode[[1]], -pmode[[2]], pmode[[3]]][a]*t/\[Delta]]]


\[Psi]k[lm_, modes_, pmodes_, \[Delta]_, a_, \[Theta]_, t_] := Module[{pos,neg={}},
(*Print["\[Psi]k:3"];*)
	pos = Table[\[Psi]pos[lm,modes[[i]],\[Delta],a,\[Theta],t], {i, 1, Length[modes]}];
	If[Length[pmodes]>0,
		neg = Table[\[Psi]neg[lm,pmodes[[i]],\[Delta],a,\[Theta],t], {i, 1, Length[pmodes]}];
	];
	pos~Join~neg
]


Clm[lm_, modes_, pmodes_, \[Delta]_, a_, \[Theta]_, t_] := Module[{pos,neg=0,msize},
	msize = Length[modes];
	(*the first msize amplitutes in \[Alpha] are for the positve modes*)
	pos = Sum[modes[[i,4]]*Exp[I*modes[[i,5]]]*\[Psi]pos[lm,modes[[i]],\[Delta],a,\[Theta],t], {i, 1, msize}];
	If[Length[pmodes]>0,
		neg = Sum[modes[[i+msize,4]]*Exp[I*modes[[i+msize,5]]]*\[Psi]neg[lm,pmodes[[i]],\[Delta],a,\[Theta],t], 
			  {i, 1, Length[pmodes]}];
	];
	pos+neg
]


SVDInverse[mat_]:=Module[{u,w,v,invw,small=10^-12,\[Sigma]2},
	{u,w,v}=SingularValueDecomposition[mat];
	invw=DiagonalMatrix[Table[If[w[[i, i]] > small, 1/w[[i, i]], 0], {i, 1, Length[w]}]];
	\[Sigma]2 = Table[Sum[invw[[j,j]]*Abs[v[[i,j]]]^2,{j,1,Length@v}],{i,1,Length@v}];
	{v.invw.ConjugateTranspose[u], \[Sigma]2}
]


SVD\[Sigma]2[mat_]:=Module[{u,w,v,invw,small=10^-12,\[Sigma]2},
	{u,w,v}=SingularValueDecomposition[mat];
	\[Sigma]2 = Table[Sum[(Abs[v[[i,j]]]^2/w[[j,j]])^2,{j,1,Length@v}],{i,1,Length@v}];
	\[Sigma]2
]


Options[GetTime]={GTDataRand->False}
GetTime[data_,lm_,it0_,OptionsPattern[]]:=Module[{i,dati,l,m,time},
	If[OptionValue[GTDataRand],
		(*random indexing*)
		dati=RandomInteger[{-it0,-1},it0],
		dati=Range[-it0,-1],dati=Range[-it0,-1]
	];
	For[i=1,i<=Length[lm],i++,
		l=lm[[i,1]];
		m=lm[[i,2]];
		time[l,m]=Table[data[[i,1,dati[[j]],1]],{j,1,it0}];
	];
	{time, dati}	
]



Options[MaxOverlap]=Union[{Mass->0.9,Spin->0.6,Theta->0}, Options@FindMaximum]
MaxOverlap[data_,lm_,modes_,pmodes_,t0_,opts:OptionsPattern[]]:=Module[{time},
	time = GetTime[data,lm,t0];
	FindMaximum[\[Rho]2[data, lm, modes, pmodes, \[Delta], a, \[Theta], t0, time], 
	{{\[Delta], OptionValue[Mass]}, {a, OptionValue[Spin]}, {\[Theta], OptionValue[Theta]}}, 
	Evaluate@FilterRules[{opts},Options@FindMaximum]]
]


DataIndex[data_,lm_]:=Table[Position[data,lm[[i]]][[1,1]],{i,1,Length[lm]}]


Options[Coefficients] = 
	Union[{Mass->0.9,Spin->0.6,Theta->0,CoeffsDataRand->False}, Options@FindMaximum]
Coefficients[data_,lm_,modes_,pmodes_,it0_,opts:OptionsPattern[]]:=
	Module[{time,\[Rho]2max,alpha,alphaExp,\[Psi]kmat,\[Psi],A,B, \[Delta]max,amax,\[Theta]max,invB,
			inv\[Psi]k,\[Sigma]2,modSize,chi2,Ndata,Mvar,dind},
	
	time = GetTime[data,lm,it0,GTDataRand->OptionValue[CoeffsDataRand]];
	dind = DataIndex[data,lm];

	\[Rho]2max = FindMaximum[\[Rho]2[data,lm,modes,pmodes,\[Delta],a,\[Theta],it0,time], 
			{{\[Delta],OptionValue[Mass]},{a,OptionValue[Spin]},{\[Theta],OptionValue[Theta]}},
			Evaluate@FilterRules[{opts},Options@FindMaximum]];

	\[Delta]max = \[Rho]2max[[2, 1, 2]];	
	amax = \[Rho]2max[[2, 2, 2]];
	\[Theta]max = \[Rho]2max[[2, 3, 2]];

	\[Psi]kmat = Flatten[Table[\[Psi]k[lm[[i]], modes, pmodes, \[Delta]max, amax, \[Theta]max, #]&/@ 
				time[[1]][lm[[i, 1]], lm[[i, 2]]], {i,1,Length[lm]}], 1];
	\[Psi] = Flatten[Table[data[[dind[[j]], 1, time[[2,i]], 2]] + I*data[[dind[[j]], 
				1, time[[2,i]], 3]], {j, 1, Length[dind]}, {i,1,it0}], 1];
(*	\[Psi]\[Psi] = Re[Conjugate[\[Psi]].\[Psi]];
	\[Sigma]2\[Psi]kmat = SVD\[Sigma]2[\[Psi]kmat];*)
	A = ConjugateTranspose[\[Psi]kmat].\[Psi];
	B = ConjugateTranspose[\[Psi]kmat].\[Psi]kmat;
	
	(*out = Eigensystem[{Outer[Times,A,Conjugate[A]],B}];*)
	{invB,\[Sigma]2} = SVDInverse[B];
	alpha = invB.A;

	chi2 = Norm[\[Psi]kmat.alpha-\[Psi]]^2;

	modSize = Length[modes];
	alphaExp = Table[{modes[[i,1]],modes[[i,2]],modes[[i,3]],Abs[alpha[[i]]],
					   Arg[alpha[[i]]]}, {i,1,modSize}];

	(*For negative modes pmodes the phase is negative because the whole negative modes
	expression is conjugated in \[Psi]k*)
	alphaExp = alphaExp~Join~Table[{pmodes[[i,1]],pmodes[[i,2]],pmodes[[i,3]],
									Abs[alpha[[modSize+i]]],-Arg[alpha[[modSize+i]]]}, 
									{i,1,Length[pmodes]}];
	Ndata = 2*Length[\[Psi]];
	Mvar = 2*Length[alpha];
	{\[Rho]2max, alphaExp, ErrorConvert[alphaExp, \[Sigma]2*chi2/(Ndata-Mvar)]}
]


ErrorConvert[\[Alpha]_, \[Sigma]2_]:=
	Table[{Sqrt@Abs@\[Sigma]2[[i]], Sqrt[ArcSin[Abs@\[Sigma]2[[i]]/\[Alpha][[i,4]]]]}, 
	{i,1,Length[\[Sigma]2]}]


MCBootStrap[data_,lm_,modes_,pmodes_,Nt_,Nstat_]:=
Module[{i,coeffs,coeffsentry,rcoeffs={},d\[Delta]=0,da=0,d\[Theta]=0,err=0},

	coeffs = Coefficients[data,lm,modes,pmodes,Nt,AccuracyGoal->6,
						Gradient->{"FiniteDifference"}];
	For[i=1,i<=Nstat,i++,
		coeffsentry = Quiet[Check[Coefficients[data,lm,modes,pmodes,Nt,AccuracyGoal->6,
							(*WorkingPrecision\[Rule]48*)Gradient->{"FiniteDifference"},
							CoeffsDataRand->True],False]];

		If[Head[coeffsentry]==List, 
			Print["Gamble: ", i];
			AppendTo[rcoeffs, coeffsentry];
			d\[Delta] += (coeffs[[1,2,1,2]]-coeffsentry[[1,2,1,2]])^2;
			da += (coeffs[[1,2,2,2]]-coeffsentry[[1,2,2,2]])^2;
			d\[Theta] += (coeffs[[1,2,3,2]]-coeffsentry[[1,2,3,2]])^2;
			,Null(*empty*)
			,i--(*Coefficients function failed, repeat the iteration*)
		];	
	];
	{rcoeffs,{Sqrt[d\[Delta]/Nstat],Sqrt[da/Nstat],Sqrt[d\[Theta]/Nstat]}}
]


MCBootStrapCError[coeffs_,mc_]:=Module[{i,cffs,cSize,dcffs,Nstat},
	cffs = coeffs[[2]];
	cSize = Length@cffs;
	dcffs = Table[0,{i,1,cSize}];
	Nstat = Length[mc[[1]]];
	For[i=1,i<=Nstat,i++,
		dcffs += Table[{(mc[[1,i,2,j,4]]-cffs[[j,4]])^2,
						(mc[[1,i,2,j,5]]-cffs[[j,5]])^2},{j,1,cSize}];
	];
	Sqrt[dcffs/Nstat]
]


Options[\[Rho]2] = {\[Rho]2DataRand->False}
\[Rho]2[data_,lm_,modes_,pmodes_,\[Delta]_?NumberQ,a_?NumberQ,\[Theta]_?NumberQ,it0_,time_,OptionsPattern[]]:=
	Module[{A,B,out,\[Psi]kmat,\[Psi],\[Psi]\[Psi],dind},
	
(*	time2 = GetTime[data,lm,it0,GTDataRand\[Rule]OptionValue[\[Rho]2DataRand]];*)
	DefineInterpolations[lm,modes,1];
	DefineInterpolations[lm,pmodes,-1];

	\[Psi]kmat = Flatten[Table[\[Psi]k[lm[[i]], modes, pmodes, \[Delta], a, \[Theta], #]&/@ 
				time[[1]][lm[[i, 1]], lm[[i, 2]]], {i,1,Length[lm]}], 1];

	dind = DataIndex[data,lm];
	\[Psi] = Flatten[Table[data[[dind[[j]], 1, time[[2,i]], 2]] + I*data[[dind[[j]], 1, time[[2,i]], 3]], 
				{j, 1, Length[dind]}, {i,1, it0}], 1];

	\[Psi]\[Psi] = Re[Conjugate[\[Psi]].\[Psi]];
	A = ConjugateTranspose[\[Psi]kmat].\[Psi];
	B = ConjugateTranspose[\[Psi]kmat].\[Psi]kmat;

	out = Re[Conjugate[A].SVDInverse[B][[1]].A/\[Psi]\[Psi]];

	(*This speeds up the calculations, but one must be careful with it*)
	(*\[Rho]2max[data,lm,modes,pmodes,\[Delta],a,\[Theta],it0] = out;*)
	out
]


OverlapPlot::nodata="cannot find a plot for mode `1`";
Options[OverlapPlot]=Union[{}, Options@Plot]
OverlapPlot[data_,coeffs_,lm_,opts:OptionsPattern[]]:=
Module[{rdata,idata,size,datasize,i,dataCnum,rlp,pr,\[Delta],a,\[Theta]},
	
	datasize=Length[data];
	For[i=1,i<=datasize,i++,
		If[data[[i,2,1]]==lm[[1]] && data[[i,2,2]]==lm[[2]],
			dataCnum=i;Break[];
		];
	(*Got to the end of the list, but the requested mode is not there*)
	If[i==datasize,Message[OverlapPlot::nodata,lm];Return[]];
	];

	size=Length[data[[dataCnum,1]]];
	rdata=Table[{data[[dataCnum,1,i,1]],data[[dataCnum,1,i,2]]},{i,1,size}];
	idata=Table[{data[[dataCnum,1,i,1]],data[[dataCnum,1,i,3]]},{i,1,size}];
	
	\[Delta] = coeffs[[1,2,1,2]];
	a = coeffs[[1,2,2,2]];
	rlp = ListPlot[rdata];

	pr = Plot[Re[Clm[lm, coeffs[[2]], {}, \[Delta], a, 0, t]], {t,0,80}, PlotStyle->Red,
			   Axes->False,Frame->True,Evaluate@FilterRules[{opts},Options@Plot]];
	Show[pr,rlp]
]


End[] (*End Private*)
EndPackage[]
