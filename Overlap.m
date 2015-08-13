(* ::Package:: *)

(*
Calculating the an overlap of a pair of waveforms
Maxim Zalutskiy
2015
*)
BeginPackage["Overlap`"]
ReadKerrQNM::usage="ReadKerrQNM[l, m, n] loads quasi-normal modes "<>
"for given l, m and n.";

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


DefineInterpolations[lm_,]=


\[Rho]2max[data_,lm_,modes_,pmodes_,\[Delta]_,a_,\[Theta]_,t0_]:=
Conjugate[A[data,lm,modes,pmodes,\[Delta],a,\[Theta],t0]]*
Inverse[B[modes,pmodes,\[Delta],a,\[Theta]]]*A[data,lm,modes,pmodes,\[Delta],a,\[Theta],t0]


End[] (*End Private*)
EndPackage[]
