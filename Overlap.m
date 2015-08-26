(* ::Package:: *)

(*
Calculating the an overlap of a pair of 
waveforms (http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.084012)

Maxim Zalutskiy, 2015
*)
BeginPackage["Overlap`"]
ReadKerrQNM::usage="ReadKerrQNM[l, m, n] loads quasi-normal modes "<>
"for given l, m and n.";

\[Rho]2max::usage="";

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
			\[Omega]bar[ll,mm,nn]=
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


\[Psi]k[lm_, modes_, pmodes_,\[Delta]_, a_, \[Theta]_, t_] := Module[{pos,neg={}},
	pos = Table[WignerD[{modes[[i, 1]],-lm[[2]],-modes[[i, 2]]},\[Theta]]*\[ScriptCapitalA][lm[[1]], modes[[i, 1]], modes[[i, 2]], modes[[i, 3]]][a]*
		Exp[-I*\[Omega]bar[modes[[i, 1]], modes[[i, 2]], modes[[i, 3]]][a]*t/\[Delta]], {i, 1, Length[modes]}];
	If[Length[pmodes]>0,
		neg = Table[(-1)^(pmodes[[i, 1]]+lm[[1]])*WignerD[{pmodes[[i, 1]],-lm[[2]],-pmodes[[i, 2]]},\[Theta]]*Conjugate[\[ScriptCapitalA][lm[[1]], pmodes[[i, 1]], -pmodes[[i, 2]], pmodes[[i, 3]]][a]*
		Exp[-I*\[Omega]bar[pmodes[[i, 1]], -pmodes[[i, 2]], pmodes[[i, 3]]][a]*t/\[Delta]]], {i, 1, Length[pmodes]}];
	];
	pos~Join~neg
]


\[Rho]2max[data_,lm_,modes_,pmodes_,\[Delta]_?NumberQ,a_?NumberQ,\[Theta]_?NumberQ,it0_]:=
	Module[{A,B,out,\[Psi]kmat,\[Psi],\[Psi]\[Psi],time,i,l,m},
	For[i=1,i<=Length[lm],i++,
		l=lm[[i,1]];
		m=lm[[i,2]];
		time[l,m]=Table[data[[i,1,j,1]],{j,-it0,-1}];
	];
	DefineInterpolations[lm,modes,1];
	DefineInterpolations[lm,pmodes,-1];

	\[Psi]kmat = Flatten[Table[\[Psi]k[lm[[i]], modes, pmodes, \[Delta], a, \[Theta], #]&/@ 
				time[lm[[i, 1]], lm[[i, 2]]], {i,1,Length[lm]}], 1];

	\[Psi] = Flatten[Table[data[[j, 1, i, 2]] + I*data[[j, 1, i, 3]], 
				{j, 1, Length[data]}, {i,-it0, -1}], 1];
	\[Psi]\[Psi] = Re[Conjugate[\[Psi]].\[Psi]];
	A = ConjugateTranspose[\[Psi]kmat].\[Psi];
	B = ConjugateTranspose[\[Psi]kmat].\[Psi]kmat;
	
	out = Re[Conjugate[A].Inverse[B].A/\[Psi]\[Psi]];
	out
	
	(*This speeds up the calculations, but one must be careful with it*)
	(*\[Rho]2max[data,lm,modes,pmodes,\[Delta],a,\[Theta],it0] = out;
	out*)
]


End[] (*End Private*)
EndPackage[]
