(* ::Package:: *)

(* ::Title:: *)
(* Open XXZ chain: symmetry sectors, eigenstate entanglement, edge dynamics *)


(* Evaluate this file with Get[".../xxz_open_edges_MaTeX.m"], or evaluate its
   sections in order. All code is active; example alternatives are commented.
   EDIT sections 3 and 6. Re-evaluate section 4 after changing the model, then
   section 5 for the sweep and section 7 for dynamics. No kernel restart needed.
   MaTeX must be installed/configured as in the reference file. *)



(* ::Section::Closed:: *)
(* 1. Physical definitions and numerical implementations *)


(* L is the TOTAL number of spins, ordered {1,...,L}; L>=3; open boundaries.
   S_j^a = sigma_j^a/2, hbar=1, |0>=up and |1>=down.

   H(Delta) = J Sum[Sx_j Sx_(j+1)+Sy_j Sy_(j+1)+Delta Sz_j Sz_(j+1)]
              + h Sum[Sz_j] + bLeft Sz_1 + bRight Sz_L.
   The bond sum is j=1,...,L-1. J is SIGNED, Jxy=J and Jz=J Delta.
   There are no additional probe spins, independent edge couplings, or defects.
   The field-free homogeneous model is obtained with h=bLeft=bRight=0.
   At Delta=1, J>0 is AFM and J<0 FM. For other Delta use the actual signed
   pair (Jxy,Jz), rather than labeling the whole model by J alone.

   System AB={1,L}, dimension 4; environment E={2,...,L-1}, dimension 2^(L-2).
   H_AB=(h+bLeft)Sz_1+(h+bRight)Sz_L.
   H_E consists of the interior bonds and fields; H_AB,E is the two end bonds.
   All bonds have the same J,Delta. Edge-state order is {00,01,10,11}.

   Full ED is performed independently in all surviving standard symmetry
   sectors: magnetization m, reflection R, spin inversion F at m=0, RF at
   m=0 when applicable, and S^2 at Delta=1 for a uniform effective field.
   At Delta=-1, staggered S^2 is resolved. Nonzero uniform h keeps S,m good
   but splits SU(2) multiplets; full SU(2) requires zero fields.
   F and RF map m to -m, so they cannot be extra simultaneous labels at m!=0.
   Redundant RF is omitted when both R and F are already specified.
   A charge is used only if its parameter conditions AND commutator check pass.
   Integrable higher conserved charges and accidental degeneracies are not
   asserted to be additional resolved symmetry sectors.

   Pure-state entanglement, in NATS:
     S_half = S(Tr_{floor(L/2)+1,...,L} |psi><psi|).
     S_edges = S(Tr_{2,...,L-1} |psi><psi|).
   S_edges is entanglement between the TWO EDGES TOGETHER and the middle,
   not entanglement between spin 1 and spin L. Their entanglement is measured
   separately by concurrence and negativity. All entropy logs are natural.
   Eigenvectors inside a degenerate eigenspace are not unique. The code fixes
   a deterministic orthonormal basis inside each resolved degeneracy; these
   entropies describe that basis, not a basis-independent degenerate subspace.

   Dynamics is the exact finite-system unitary evolution, with no master
   equation, bath limit, weak-coupling, Markov, or smoothing approximation.
   For mixed full states S(rho_AB) is a reduced-state entropy, not a pure-state
   entanglement entropy. No dynamical map/Choi object is inferred from a single
   trajectory or from an initially correlated full-chain state.
*)

ClearAll["XXZE`*"];

XXZE`Require[test_, message_] := If[!TrueQ[test],
  Throw[Failure["XXZRun", <|"Message" -> message|>], "XXZFailure"]];
XXZE`RealQ[x_] := NumericQ[x] && NumberQ[N[x]] && TrueQ[Im[N[x]] == 0];
XXZE`Id[n_] := SparseArray[Band[{1,1}] -> 1, {n,n}];
XXZE`MaxAbs[a_] := If[Head[a] === SparseArray,
  Max[0, Sequence @@ Abs[a["NonzeroValues"]]], Max[0, Sequence @@ Abs[Flatten[a]]]];
XXZE`Hermitian[a_] := (a + ConjugateTranspose[a])/2;
XXZE`KetDensity[v_] := Outer[Times, v, Conjugate[v]];
XXZE`SiteOp[a_, k_, l_] := KroneckerProduct[
  XXZE`Id[2^(k-1)], SparseArray[a], XXZE`Id[2^(l-k)]];
XXZE`BondOp[a_, b_, k_, l_] := KroneckerProduct[
  XXZE`Id[2^(k-1)], SparseArray[a], SparseArray[b], XXZE`Id[2^(l-k-1)]];

XXZE`Hamiltonian[delta_?NumericQ, p_Association] := Module[
  {l=p["L"], j=p["J"], h, x, y, z},
  XXZE`Require[IntegerQ[l] && l>=3, "L must be an integer >=3."];
  XXZE`Require[And@@(XXZE`RealQ /@ {delta,j,p["h"],p["bLeft"],p["bRight"]}),
    "Delta, J and fields must be finite real numbers."];
  {x,y,z}=Table[PauliMatrix[k]/2,{k,3}];
  h=SparseArray[{}, {2^l,2^l}];
  Do[h += j (XXZE`BondOp[x,x,k,l]+XXZE`BondOp[y,y,k,l]+
      delta XXZE`BondOp[z,z,k,l]), {k,l-1}];
  Do[h += p["h"] XXZE`SiteOp[z,k,l], {k,l}];
  h += p["bLeft"] XXZE`SiteOp[z,1,l]+p["bRight"] XXZE`SiteOp[z,l,l];
  SparseArray[h]
];

(* General partial trace for a factor R satisfying rho=R.R^dagger.
   This avoids constructing the full density matrix during time evolution. *)
XXZE`ReducedFactor[fac_, sites_List, l_Integer] := Module[
  {r=Dimensions[fac][[2]], rest, a},
  rest=Complement[Range[l],sites];
  a=ArrayReshape[
    Transpose[ArrayReshape[fac,Join[ConstantArray[2,l],{r}]],
      Ordering[Join[sites,rest,{l+1}]]],
    {2^Length[sites],2^Length[rest] r}];
  a . ConjugateTranspose[a]
];
XXZE`ReducedKet[psi_, sites_List, l_Integer] :=
  XXZE`ReducedFactor[Transpose[{psi}],sites,l];
XXZE`EdgesDensity[psi_, l_Integer] := XXZE`ReducedKet[psi,{1,l},l];
XXZE`Entropy[rho_] := Module[{e},
  e=Select[Re[Eigenvalues[N[XXZE`Hermitian[rho]]]], #>0 &];
  If[e==={},0.,-Total[e Log[e]]]
];
XXZE`HalfEntropy[psi_, l_Integer] := Module[{sv,p},
  sv=SingularValueList[ArrayReshape[psi,{2^Floor[l/2],2^Ceiling[l/2]}]];
  p=Select[sv^2,#>0&]; If[p==={},0.,-Total[p Log[p]]]
];
XXZE`EdgesEntropy[psi_, l_Integer] := XXZE`Entropy[XXZE`EdgesDensity[psi,l]];
XXZE`TraceB[r_] := {{r[[1,1]]+r[[2,2]],r[[1,3]]+r[[2,4]]},
  {r[[3,1]]+r[[4,2]],r[[3,3]]+r[[4,4]]}};
XXZE`TraceA[r_] := {{r[[1,1]]+r[[3,3]],r[[1,2]]+r[[3,4]]},
  {r[[2,1]]+r[[4,3]],r[[2,2]]+r[[4,4]]}};
XXZE`PartialTransposeB[r_] := ArrayReshape[
  Transpose[ArrayReshape[r,{2,2,2,2}],{1,4,3,2}],{4,4}];
XXZE`Negativity[r_] := Total[Clip[-Re[Eigenvalues[
  XXZE`Hermitian[XXZE`PartialTransposeB[r]]]],{0,Infinity}]];
XXZE`Concurrence[r_] := Module[{ev,rows,root,flip,s},
  {ev,rows}=Eigensystem[N[XXZE`Hermitian[r]]];
  root=Transpose[rows] . (Sqrt[Clip[Re[ev],{0,Infinity}]] Conjugate[rows]);
  flip=KroneckerProduct[PauliMatrix[2],PauliMatrix[2]];
  s=Reverse[Sort[Sqrt[Clip[Re[Eigenvalues[XXZE`Hermitian[
    root . flip . Conjugate[r] . flip . root]]],{0,Infinity}]]]];
  Max[0.,First[s]-Total[Rest[s]]]
];

(* State constructors: pure vectors or arbitrary positive density matrices
   are accepted below. ProductKet normalizes each local two-component vector. *)
XXZE`ProductKet[locals_List] := Module[{v},
  XXZE`Require[Length[locals]>0 && AllTrue[locals,
    VectorQ[#,NumericQ] && Length[#]==2 && Norm[#]>0&],
    "ProductKet requires nonzero two-component local vectors."];
  v=Fold[Flatten[KroneckerProduct[#1,N[#2/Norm[#2]]]]&, {1.}, locals]; v/Norm[v]
];
XXZE`NeelKet[n_Integer] := XXZE`ProductKet[
  Table[If[OddQ[k],{1,0},{0,1}],{k,n}]];
XXZE`HelixKet[n_Integer,q_?NumericQ,theta_:Pi/2,phi_:0] := XXZE`ProductKet[
  Table[{Cos[theta/2],Exp[I(phi+(k-1)q)] Sin[theta/2]},{k,n}]];
(* HelixKet is a preparation, not a claim that it is an eigenstate of this
   open Hermitian Hamiltonian. No tuned transverse boundary fields are added. *)
XXZE`StateFactor[input_, dim_Integer, tol_, label_] := Module[
  {x=N[Normal[input]],e,v,keep,fac},
  If[VectorQ[x,NumericQ],
    XXZE`Require[Length[x]==dim && Norm[x]>0,label<>" has an invalid vector."];
    Return[Transpose[{x/Norm[x]}]]];
  XXZE`Require[MatrixQ[x,NumericQ] && Dimensions[x]=={dim,dim},
    label<>" must be a vector or density matrix of the required dimension."];
  XXZE`Require[Norm[x-ConjugateTranspose[x],"Frobenius"]<=tol && Abs[Tr[x]-1]<=tol,
    label<>" density matrix must be Hermitian with trace one."];
  {e,v}=Eigensystem[XXZE`Hermitian[x]];
  XXZE`Require[Min[Re[e]]>=-tol,label<>" density matrix has a negative eigenvalue."];
  e=Clip[Re[e],{0,Infinity}];e=e/Total[e];
  keep=Select[Range[dim],e[[#]]>0&];
  fac=Transpose[Sqrt[e[[keep]]] v[[keep]]];fac
];
(* Tensor order AB,E is rearranged into the physical order 1,2,...,L. *)
XXZE`JoinEdgeMiddle[edgeFactor_, middleFactor_, l_Integer] := Module[{a},
  a=KroneckerProduct[edgeFactor,middleFactor];
  ArrayReshape[Transpose[ArrayReshape[a,{2,2,2^(l-2),Dimensions[a][[2]]}],
    {1,3,2,4}],{2^l,Dimensions[a][[2]]}]
];

(* Symmetries and deterministic bases in degenerate eigenspaces. *)
XXZE`CanonicalRows[rows_,tol_] := Module[{p,chosen={},v},
  If[Length[rows]==1,Return[rows/Norm[First[rows]]]];
  p=Transpose[rows] . Conjugate[rows];
  Do[v=p[[All,k]];
    Do[v-=Conjugate[w] . v w,{w,chosen}];
    If[Norm[v]>tol,AppendTo[chosen,v/Norm[v]]];
    If[Length[chosen]==Length[rows],Break[]],{k,Length[p]}];
  XXZE`Require[Length[chosen]==Length[rows],"Degenerate eigenbasis construction failed."];
  chosen
];
XXZE`FixPhase[v_,tol_] := Module[{k},
  k=SelectFirst[Range[Length[v]],Abs[v[[#]]]>tol&,1];v Exp[-I Arg[v[[k]]]]
];
XXZE`Charges[l_Integer,delta_,p_] := Module[{d=2^l,r,f,fields,out=<||>,st,sx,sy,sz},
  fields=ConstantArray[p["h"],l];
  fields[[1]]+=p["bLeft"];fields[[-1]]+=p["bRight"];
  r=SparseArray[Table[{1+FromDigits[Reverse[IntegerDigits[k,2,l]],2],k+1}->1,
    {k,0,d-1}],{d,d}];
  f=SparseArray[Table[{d-k,k+1}->1,{k,0,d-1}],{d,d}];
  If[And@@MapThread[TrueQ[#1==#2]&,{fields,Reverse[fields]}],AssociateTo[out,"R"->r]];
  If[And@@(TrueQ[#==0]& /@ fields),AssociateTo[out,"F"->f]];
  If[And@@MapThread[TrueQ[#1 == -#2]&,{fields,Reverse[fields]}] &&
    !(KeyExistsQ[out,"R"] && KeyExistsQ[out,"F"]),AssociateTo[out,"RF"->r . f]];
  If[(TrueQ[delta==1] || TrueQ[delta == -1]) &&
    And@@(TrueQ[#==First[fields]]& /@ fields),
    st=If[TrueQ[delta==1],ConstantArray[1,l],Table[(-1)^(k-1),{k,l}]];
    sx=Total[Table[st[[k]] XXZE`SiteOp[PauliMatrix[1]/2,k,l],{k,l}]];
    sy=Total[Table[st[[k]] XXZE`SiteOp[PauliMatrix[2]/2,k,l],{k,l}]];
    sz=Total[Table[XXZE`SiteOp[PauliMatrix[3]/2,k,l],{k,l}]];
    AssociateTo[out,If[TrueQ[delta==1],"S2","StaggeredS2"]->(sx . sx+sy . sy+sz . sz)]];
  out
];
XXZE`SplitCharge[s_,q_,name_,tol_] := Module[{b=s["Basis"],e,v,groups,label,rows},
  {e,v}=Eigensystem[N[XXZE`Hermitian[ConjugateTranspose[b] . q . b]]];
  groups=Gather[Range[Length[e]],Abs[e[[#1]]-e[[#2]]]<tol&];
  Table[
    rows=XXZE`CanonicalRows[Orthogonalize[v[[g]]],tol];
    label=If[MemberQ[{"S2","StaggeredS2"},name],
      Round[(-1+Sqrt[1+4 Re[Mean[e[[g]]]]])/2,1/2],Round[Re[Mean[e[[g]]]]]];
    <|"Basis"->b . Transpose[rows],"Labels"->Append[s["Labels"],
      Switch[name,"S2","S","StaggeredS2","StaggeredS",_,name]->label]|>,{g,groups}]
];
XXZE`Diagonalize[delta_?NumericQ,p_Association,num_Association] := Catch[Module[
  {l=p["L"],h,charges,active,tol=num["Tolerance"],scale,ids,hm,m,names,
   spaces,out={},b,hs,e,v,ord,groups,rows,residual,leak,records,orth},
  h=XXZE`Hamiltonian[delta,p];scale=Max[1.,Norm[h,Infinity]];
  charges=XXZE`Charges[l,delta,p];
  active=Select[Keys[charges],XXZE`MaxAbs[h . charges[#]-charges[#] . h]<tol scale&];
  Do[
    m=l/2-nDown;
    ids=Select[Range[2^l],Total[IntegerDigits[#-1,2,l]]==nDown&];
    hm=N[Normal[h[[ids,ids]]]];
    names=Select[active,!MemberQ[{"F","RF"},#] || m==0&];
    names=SortBy[names,Switch[#,"S2"|"StaggeredS2",0,"R",1,_,2]&];
    spaces={<|"Basis"->IdentityMatrix[Length[ids]],"Labels"-><|"m"->m|>|>};
    Do[spaces=Flatten[XXZE`SplitCharge[#,N[Normal[charges[name][[ids,ids]]]],
        name,tol]& /@ spaces,1],{name,names}];
    Do[
      b=s["Basis"];hs=XXZE`Hermitian[ConjugateTranspose[b] . hm . b];
      {e,v}=Eigensystem[N[hs]];ord=Ordering[Re[e]];e=Re[e[[ord]]];v=v[[ord]];
      groups=Gather[Range[Length[e]],Abs[e[[#1]]-e[[#2]]]<tol scale&];
      Do[v[[g]]=XXZE`CanonicalRows[Orthogonalize[v[[g]]],tol],{g,groups}];
      rows=Transpose[b . Transpose[v]];
      rows=XXZE`FixPhase[#,tol]& /@ rows;
      residual=Max[MapThread[Norm[hm . #2-#1 #2]&,{e,rows}]];
      orth=XXZE`MaxAbs[rows . ConjugateTranspose[rows]-IdentityMatrix[Length[rows]]];
      leak=XXZE`MaxAbs[hm . b-b . hs];
      XXZE`Require[residual<100 tol scale && leak<100 tol scale && orth<100 tol,
        "Sector eigenpairs failed validation: "<>ToString[s["Labels"],InputForm]];
      AppendTo[out,<|"Labels"->s["Labels"],"Dimension"->Length[e],
        "Indices"->ids,"Eigenvalues"->e,"EigenvectorsInMagnetizationBasis"->rows,
        "Residual"->residual,"OrthogonalityError"->orth,"Leakage"->leak|>],
      {s,spaces}],{nDown,0,l}];
  XXZE`Require[Total[Lookup[out,"Dimension"]]==2^l,"The sector dimensions are incomplete."];
  <|"Delta"->delta,"Parameters"->p,"Hamiltonian"->h,"ActiveCharges"->active,
    "Sectors"->out,"DimensionSum"->Total[Lookup[out,"Dimension"]],
    "MaxResidual"->Max[Lookup[out,"Residual"]]|>
],"XXZFailure"];
XXZE`Eigenvector[ed_,sector_Integer,level_Integer] := Module[{s=ed["Sectors"][[sector]]},
  Normal[SparseArray[Thread[s["Indices"]->
    s["EigenvectorsInMagnetizationBasis"][[level]]],{2^ed["Parameters"]["L"]}]]
];
XXZE`GroundKet[ed_] := Module[{pairs,choice},
  pairs=Flatten[Table[{ed["Sectors"][[k]]["Eigenvalues"][[r]],k,r},
    {k,Length[ed["Sectors"]]},{r,ed["Sectors"][[k]]["Dimension"]}],1];
  choice=First[SortBy[pairs,First]]; XXZE`Eigenvector[ed,choice[[2]],choice[[3]]]
];
XXZE`EigenstateEntropies[ed_] := Module[{out=ed,l=ed["Parameters"]["L"],s,v},
  out["Sectors"]=Table[
    s=ed["Sectors"][[k]];
    v=Table[XXZE`Eigenvector[ed,k,r],{r,s["Dimension"]}];
    Join[s,<|"HalfEntropy"->(XXZE`HalfEntropy[#,l]& /@ v),
      "EdgesEntropy"->(XXZE`EdgesEntropy[#,l]& /@ v)|>],{k,Length[ed["Sectors"]]}];out
];
XXZE`SweepLabels[labels_] := KeySort[KeyDrop[labels,{"S","StaggeredS"}]];
XXZE`PackEigenstates[ed_] := Module[{groups},
  groups=GatherBy[ed["Sectors"],ToString[Normal[XXZE`SweepLabels[# ["Labels"]]],InputForm]&];
  Association[Table[With[{ordered=SortBy[Flatten[
    Table[Transpose[{s["Eigenvalues"],s["HalfEntropy"],s["EdgesEntropy"]}],{s,g}],1],First]},
    ToString[Normal[XXZE`SweepLabels[First[g]["Labels"]]],InputForm]->
    <|"Labels"->XXZE`SweepLabels[First[g]["Labels"]],"Energies"->ordered[[All,1]],
      "HalfEntropy"->ordered[[All,2]],"EdgesEntropy"->ordered[[All,3]]|>],{g,groups}]]
];
XXZE`SweepPoint[delta_,p_,num_] := Module[{ed=XXZE`Diagonalize[delta,p,num]},
  If[FailureQ[ed],Return[ed]];
  XXZE`PackEigenstates[XXZE`EigenstateEntropies[ed]]
];
XXZE`IsospectralQ[a_,b_,samples_,tol_] := Module[{ea,eb},
  ea=(#[a]["Energies"]& /@ samples);eb=(#[b]["Energies"]& /@ samples);
  Dimensions[ea]===Dimensions[eb] && XXZE`MaxAbs[ea-eb]<tol Max[1.,XXZE`MaxAbs[{ea,eb}]]
];

(* Evolution uses all the separately computed sector eigenvectors. Only the
   compact spectral coefficients and occupied initial sector projections are
   sent to time workers, once per chunk. The Hamiltonian is never re-diagonalized
   during evolution. Cross-sector coherences are preserved when assembling rho. *)
XXZE`EvolutionEngine[ed_,fac_] := Module[{sectors},
  sectors=Table[With[{s=s},<|"Indices"->s["Indices"],"E"->s["Eigenvalues"],
    "V"->Transpose[s["EigenvectorsInMagnetizationBasis"]],
    "C"->Conjugate[s["EigenvectorsInMagnetizationBasis"]] . fac[[s["Indices"],All]]|>],
    {s,ed["Sectors"]}];
  <|"Dimension"->Length[fac],"Columns"->Dimensions[fac][[2]],
    "L"->ed["Parameters"]["L"],"Sectors"->sectors|>
];
XXZE`Evolve[t_,engine_] := Module[{a},
  a=ConstantArray[0.+0.I,{engine["Dimension"],engine["Columns"]}];
  Do[a[[s["Indices"],All]] += s["V"] . (Exp[-I s["E"] t] s["C"]),
    {s,engine["Sectors"]}];a
];
XXZE`Observables[t_,rho_,rho0_] := Module[{a,b,sa,sb,sab,z,x,ma,mb,e},
  a=XXZE`TraceB[rho];b=XXZE`TraceA[rho];
  {sa,sb,sab}=XXZE`Entropy /@ {a,b,rho};{x,z}=PauliMatrix /@ {1,3};
  ma=Re[Tr[z . a]];mb=Re[Tr[z . b]];e=Reverse[Sort[Re[Eigenvalues[XXZE`Hermitian[rho]]]]];
  <|"Time"->t,"RhoAB"->rho,"RhoA"->a,"RhoB"->b,"EigenvaluesAB"->e,
    "Purity"->Re[Tr[rho . rho]],"EntropyA"->sa,"EntropyB"->sb,"EntropyAB"->sab,
    "MutualInfo"->sa+sb-sab,"Concurrence"->XXZE`Concurrence[rho],
    "Negativity"->XXZE`Negativity[rho],"InitialOverlap"->Re[Tr[rho0 . rho]],
    "MagnetizationA"->ma,"MagnetizationB"->mb,
    "ConnectedZZ"->Re[Tr[KroneckerProduct[z,z] . rho]]-ma mb,
    "ConnectedXX"->Re[Tr[KroneckerProduct[x,x] . rho]]-Re[Tr[x . a]] Re[Tr[x . b]],
    "TraceError"->Abs[Tr[rho]-1.],
    "HermiticityError"->Norm[rho-ConjugateTranspose[rho],"Frobenius"],
    "MinEigenvalue"->Last[e]|>
];
XXZE`TimeChunk[times_,engine_,rho0_,tol_] := Catch[Module[{a,r,records},
  records=Table[
    a=XXZE`Evolve[t,engine];r=XXZE`ReducedFactor[a,{1,engine["L"]},engine["L"]];
    XXZE`Observables[t,r,rho0],{t,times}];
  XXZE`Require[Max[Lookup[records,"TraceError"]]<tol &&
    Max[Lookup[records,"HermiticityError"]]<tol &&
    Min[Lookup[records,"MinEigenvalue"]]>=-tol,
    "A reduced state failed trace, Hermiticity, or positivity validation."];
  <|"KernelID"->$KernelID,"Results"->records|>
],"XXZFailure"];
XXZE`TimeGrid[cfg_] := Module[{t},
  XXZE`Require[XXZE`RealQ[cfg["tMax"]] && cfg["tMax"]>0,"tMax must be positive."];
  Switch[cfg["Grid"],
    "Linear",
      XXZE`Require[XXZE`RealQ[cfg["dt"]] && cfg["dt"]>0,"dt must be positive."];
      t=N[cfg["dt"] Range[0,Floor[cfg["tMax"]/cfg["dt"]]]];
      If[Abs[Last[t]-cfg["tMax"]]<10^-12 Max[1.,cfg["tMax"]],
        t[[-1]]=N[cfg["tMax"]],AppendTo[t,N[cfg["tMax"]]]];t,
    "Log",
      XXZE`Require[XXZE`RealQ[cfg["tMinPositive"]] &&
        0<cfg["tMinPositive"]<cfg["tMax"] && IntegerQ[cfg["LogPoints"]] && cfg["LogPoints"]>=2,
        "Log grid requires 0<tMinPositive<tMax and LogPoints>=2."];
      Join[{0.},Exp[Subdivide[Log[N[cfg["tMinPositive"]]],Log[N[cfg["tMax"]]],cfg["LogPoints"]-1]]],
    _,XXZE`Require[False,"Grid must be Linear or Log."]]
];



(* ::Section:: *)
(* 2. Presentation, MaTeX labels, output files, and execution helpers *)


(* This section contains no Hamiltonian or observable definitions.
   Chop is applied only to plotted coordinates, never to saved raw data.
   Plot captions carry the Hamiltonian parameters, not initial-state names.
   PNG only; WXF contains numerical data. A numbered rerun suffix prevents
   overwriting, with no random strings or times in the folder name. *)

XXZE`$TeXReady=Automatic;
XXZE`TeX[s_String,size_:16] := Module[{g},
  If[XXZE`$TeXReady===Automatic,
    XXZE`$TeXReady=Quiet[Check[Needs["MaTeX`"];True,False]];
    If[!TrueQ[XXZE`$TeXReady],Print["MaTeX unavailable; using native labels. Configure MaTeX as in the reference file."]]];
  If[TrueQ[XXZE`$TeXReady],
    g=Quiet[Check[MaTeX`MaTeX[s,FontSize->size],$Failed]];
    If[MatchQ[g,_Graphics],Return[g]];
    XXZE`$TeXReady=False;Print["MaTeX rendering failed; check TeX and Ghostscript configuration."]];
  Style[s,Black,FontFamily->"Times",FontSize->size]
];
XXZE`TeXNumber[x_] := Module[{s,parts},
  If[TrueQ[x==0],Return["0"]];
  If[TrueQ[x==Round[x]] && Abs[x]<10^6,Return[ToString[Round[x]]]];
  s=StringReplace[ToString[N[x],InputForm],RegularExpression["`[0-9.]*"]->""];
  parts=StringSplit[s,"*^"];parts[[1]]=StringReplace[parts[[1]],RegularExpression["\\.$"]->""];
  If[Length[parts]==2,parts[[1]]<>"\\times 10^{"<>parts[[2]]<>"}",First[parts]]
];
XXZE`Ticks[lo_?NumericQ,hi_?NumericQ] := Module[{a=lo,b=hi,pad},
  If[a==b,pad=Max[1.,.05 Abs[a]];a-=pad;b+=pad];
  ({#,XXZE`TeX[XXZE`TeXNumber[#],14],{.012,0}}& /@ FindDivisions[{a,b},5])
];
XXZE`LogTicks[lo_?NumericQ,hi_?NumericQ] := Module[{exponents},
  exponents=Range[Floor[lo/Log[10]],Ceiling[hi/Log[10]]];
  ({N[# Log[10]],XXZE`TeX["10^{"<>ToString[#]<>"}",14],{.012,0}}& /@ exponents)
];
XXZE`Caption[p_,deltaInfo_String,extra_String:""] := Column[
  XXZE`TeX[#,13]& /@ DeleteCases[{
    "\\mathrm{XXZ~OBC}:\\quad L="<>XXZE`TeXNumber[p["L"]]<>
      ",\\quad J="<>XXZE`TeXNumber[p["J"]]<> ",\\quad J_z=J\\Delta,\\quad S^a=\\sigma^a/2",
    "h="<>XXZE`TeXNumber[p["h"]]<> ",\\quad b_L="<>XXZE`TeXNumber[p["bLeft"]]<>
      ",\\quad b_R="<>XXZE`TeXNumber[p["bRight"]]<> ",\\quad AB=\\{1,L\\}",
    deltaInfo,extra},""],Alignment->Center,Spacings->.3];
XXZE`SectorLabel[a_] := StringRiffle[KeyValueMap[
  Switch[#1,"R","r","F","f","RF","rf","StaggeredS","S_{\\mathrm{stag}}",_,#1]<>
    "="<>XXZE`TeXNumber[#2]&,a],",\\quad "];
XXZE`Styles[n_] := Table[Switch[k,
  1,Directive[Blue,AbsoluteThickness[2]],
  2,Directive[Orange,AbsoluteThickness[2],Dashed],
  _,Directive[ColorData[97][k],AbsoluteThickness[2],Dashing[{.012,.004,.003,.004}]]],{k,n}];
XXZE`Save[file_,data_,format_] := Module[{r=Quiet[Check[Export[file,data,format],$Failed]]},
  XXZE`Require[StringQ[r],"Export failed: "<>file];r];
XXZE`BaseDirectory[] := Module[{dir},
  dir=If[StringQ[$InputFileName] && StringLength[$InputFileName]>0,
    DirectoryName[ExpandFileName[$InputFileName]],Quiet[Check[NotebookDirectory[],Directory[]]]];
  If[StringQ[dir],dir,Directory[]]];
XXZE`Token[x_] := Module[{s},
  s=StringReplace[ToString[N[x],InputForm],RegularExpression["`[0-9.]*"]->""];
  s=StringReplace[s,RegularExpression["\\.$"]->""];
  StringReplace[s,{"*^"->"e","."->"p","-"->"m","+"->""}]];
XXZE`NewDirectory[root_,tag_] := Module[{dir,k=1},
  XXZE`Require[StringQ[root] && StringQ[tag] &&
    StringMatchQ[tag,(LetterCharacter|DigitCharacter|"_"|"-")..],"Invalid OutputRoot or RunLabel."];
  dir=FileNameJoin[{root,tag}];
  While[DirectoryQ[dir],k++;dir=FileNameJoin[{root,tag<>"_run"<>ToString[k]}]];
  CreateDirectory[dir,CreateIntermediateDirectories->True];dir];

XXZE`SweepPlots[samples_,grid_,groups_,p_,cfg_,dir_] := Module[
  {curves,styles,labels,branches,data,plots,panel,fields,ylabels,titles,energies},
  fields={"Energies","HalfEntropy","EdgesEntropy"};
  ylabels={"E","S_{\\mathrm{half}}\\;(\\mathrm{nats})","S_{AB|E}\\;(\\mathrm{nats})"};
  titles={"\\text{Energy spectrum}","\\text{Half-chain eigenstate entanglement}",
    "\\text{Edges versus middle eigenstate entanglement}"};
  Do[
    styles=XXZE`Styles[Length[group]];
    labels=XXZE`TeX[XXZE`SectorLabel[First[samples][#]["Labels"]],14]& /@ group;
    plots=Table[
      curves={};data={};
      Do[energies=(#[group[[k]]][fields[[j]]]& /@ samples);
        branches=(Transpose[{N[grid],#}]& /@ Transpose[energies]);
        curves=Join[curves,branches];
        data=Join[data,ConstantArray[styles[[k]],Length[branches]]],{k,Length[group]}];
      ListLinePlot[Chop[curves,cfg["ChopTolerance"]],Joined->(j==1),
        InterpolationOrder->1,PlotStyle->data,PlotRange->All,Frame->True,Axes->False,
        FrameStyle->Directive[Black,AbsoluteThickness[1.5]],
        FrameTicks->{{XXZE`Ticks,None},{XXZE`Ticks,None}},
        FrameLabel->{XXZE`TeX["\\Delta",18],XXZE`TeX[ylabels[[j]],17]},
        PlotLabel->XXZE`TeX[titles[[j]],17],
        PlotMarkers->If[j==1,None,{Automatic,4}],ImageSize->600,
        AspectRatio->.75,Background->White],{j,3}];
    panel=Column[{
      XXZE`Caption[p,"\\Delta\\in["<>XXZE`TeXNumber[First[grid]]<>","<>
        XXZE`TeXNumber[Last[grid]]<>"];\\quad N_\\Delta="<>ToString[Length[grid]]],
      LineLegend[styles,labels],GraphicsRow[plots,ImageSize->1900,Spacings->15]},
      Alignment->Center,Spacings->.6];
    XXZE`Save[FileNameJoin[{dir,"sector_group_"<>ToString[gi]<>".png"}],panel,"PNG"],
    {gi,Length[groups]},{group,{groups[[gi]]}}]
];
XXZE`FixedSectorPlots[ed_,groups_,p_,cfg_,dir_] := Module[
  {styles,labels,series,plots,panel,s,metrics,ylabels},
  metrics={"Eigenvalues","HalfEntropy","EdgesEntropy"};
  ylabels={"E_n","S_{\\mathrm{half}}\\;(\\mathrm{nats})","S_{AB|E}\\;(\\mathrm{nats})"};
  Do[
    styles=XXZE`Styles[Length[group]];
    labels=XXZE`TeX[XXZE`SectorLabel[ed["Sectors"][[#]]["Labels"]],14]& /@ group;
    plots=Table[
      series=Table[s=ed["Sectors"][[k]];Transpose[{Range[s["Dimension"]],s[metrics[[j]]]}],{k,group}];
      ListLinePlot[Chop[series,cfg["ChopTolerance"]],Joined->(j==1),PlotStyle->styles,
        PlotMarkers->If[j==1,None,{Automatic,5}],Frame->True,Axes->False,PlotRange->All,
        FrameTicks->{{XXZE`Ticks,None},{XXZE`Ticks,None}},
        FrameStyle->Directive[Black,AbsoluteThickness[1.5]],
        FrameLabel->{XXZE`TeX["n\\quad(\\text{energy order within sector})",15],
          XXZE`TeX[ylabels[[j]],17]},ImageSize->600,AspectRatio->.75,Background->White],{j,3}];
    panel=Column[{XXZE`Caption[p,"\\Delta="<>XXZE`TeXNumber[ed["Delta"]]],
      LineLegend[styles,labels],GraphicsRow[plots,ImageSize->1900]},Alignment->Center];
    XXZE`Save[FileNameJoin[{dir,"fixed_sector_group_"<>ToString[gi]<>".png"}],panel,"PNG"],
    {gi,Length[groups]},{group,{groups[[gi]]}}]
];
XXZE`Series[records_,key_] := ({#["Time"],#[key]}& /@ records);
XXZE`MatrixSeries[records_,key_,i_,j_,part_] := ({#["Time"],part[#[key][[i,j]]]}& /@ records);
XXZE`TimePlot[records_,series_,legends_,title_,ylabel_,name_,p_,delta_,time_,cfg_,dir_] :=
 Module[{data,styles,plot,panel,xt,log,note},
  styles=Table[Directive[ColorData[97][k],AbsoluteThickness[2]],{k,Length[series]}];
  Do[
    data=If[log,Select[#,First[#]>0&]& /@ series,series];
    data=Chop[data,cfg["ChopTolerance"]];
    (* Avoid Chop sending a positive logarithmic time to zero. *)
    If[log,data=Select[#,First[#]>0&]& /@ data];
    If[AllTrue[data,Length[#]>0&],
      xt=If[log,XXZE`LogTicks,XXZE`Ticks];
      plot=ListLinePlot[data,Frame->True,Axes->False,PlotRange->All,
        ScalingFunctions->If[log,{"Log",None},None],
        FrameStyle->Directive[Black,AbsoluteThickness[1.6]],
        FrameTicks->{{XXZE`Ticks,None},{xt,None}},
        FrameLabel->{XXZE`TeX["t\\quad(\\hbar=1)",20],XXZE`TeX[ylabel,18]},
        PlotStyle->styles,InterpolationOrder->1,
        PlotLegends->Placed[LineLegend[styles,XXZE`TeX[#,15]& /@ legends,
          LegendLayout->"Column"],Right],ImageSize->1100,AspectRatio->.48,Background->White];
      note="t\\in[0,"<>XXZE`TeXNumber[time["tMax"]]<>"];\\quad N_t="<>
        ToString[Length[records]]<>";\\quad \\text{grid: "<>time["Grid"]<>"}"<>
        If[time["Grid"]=="Linear",";\\quad \\delta t="<>XXZE`TeXNumber[time["dt"]],
          ";\\quad t_{\\min,+}="<>XXZE`TeXNumber[time["tMinPositive"]]];
      panel=Column[{XXZE`TeX[title,22],XXZE`Caption[p,"\\Delta="<>XXZE`TeXNumber[delta],note],plot},
        Alignment->Center,Spacings->.6];
      XXZE`Save[FileNameJoin[{dir,name<>If[log,"_logtime",""]<>".png"}],panel,"PNG"]],
    {log,If[TrueQ[cfg["ExportLogTime"]],{False,True},{False}]}]
];
XXZE`DynamicsPlots[records_,p_,delta_,time_,cfg_,dir_] := Module[
  {basis={"00","01","10","11"},pairs=Subsets[Range[4],{2}],emit,legends,series},
  emit[name_,title_,ylabel_,data_,labels_] :=
    XXZE`TimePlot[records,data,labels,title,ylabel,name,p,delta,time,cfg,dir];
  emit["01_populations","\\text{Edge populations}","\\rho_{AB;ss}(t)",
    Table[XXZE`MatrixSeries[records,"RhoAB",k,k,Re],{k,4}],
    ("\\rho_{"<>#<> ","<>#<>"}"& /@ basis)];
  Do[
    legends=("\\operatorname{"<>part<>"}\\rho_{"<>basis[[#[[1]]]]<>","<>basis[[#[[2]]]]<>"}"& /@ pairs);
    If[part=="Abs",legends=("|\\rho_{"<>basis[[#[[1]]]]<>","<>basis[[#[[2]]]]<>"}|"& /@ pairs)];
    series=(XXZE`MatrixSeries[records,"RhoAB",#[[1]],#[[2]],
      Switch[part,"Re",Re,"Im",Im,"Abs",Abs]]& /@ pairs);
    emit["02_coherences_"<>part,"\\text{All six edge coherences: "<>part<>"}",
      "\\rho_{AB;ij}(t)",series,legends],{part,{"Re","Im","Abs"}}];
  emit["03_entropies_mutual_information","\\text{Reduced-state entropies and mutual information}",
    "S,I\\;(\\mathrm{nats})",XXZE`Series[records,#]& /@ {"EntropyA","EntropyB","EntropyAB","MutualInfo"},
    {"S(\\rho_A)","S(\\rho_B)","S(\\rho_{AB})","I(A:B)=S_A+S_B-S_{AB}"}];
  emit["04_edge_entanglement","\\text{Entanglement between the two edge spins}","C,\\mathcal N",
    XXZE`Series[records,#]& /@ {"Concurrence","Negativity"},{"C(\\rho_{AB})","\\mathcal N(\\rho_{AB})"}];
  emit["05_purity_initial_overlap","\\text{Purity and initial reduced-state overlap}","\\text{Purity / overlap}",
    XXZE`Series[records,#]& /@ {"Purity","InitialOverlap"},
    {"\\operatorname{Tr}\\rho_{AB}(t)^2","\\operatorname{Tr}[\\rho_{AB}(0)\\rho_{AB}(t)]"}];
  emit["06_magnetizations","\\text{Local edge magnetizations}","\\langle\\sigma^z\\rangle",
    XXZE`Series[records,#]& /@ {"MagnetizationA","MagnetizationB"},
    {"\\langle\\sigma_1^z\\rangle","\\langle\\sigma_L^z\\rangle"}];
  emit["07_connected_correlations","\\text{Connected edge correlations}","\\langle\\sigma_1^k\\sigma_L^k\\rangle_c",
    XXZE`Series[records,#]& /@ {"ConnectedZZ","ConnectedXX"},
    {"\\langle Z_1Z_L\\rangle-\\langle Z_1\\rangle\\langle Z_L\\rangle",
     "\\langle X_1X_L\\rangle-\\langle X_1\\rangle\\langle X_L\\rangle"}];
  emit["08_density_validation","\\text{Reduced-state numerical checks}","\\text{Error / minimum eigenvalue}",
    XXZE`Series[records,#]& /@ {"TraceError","HermiticityError","MinEigenvalue"},
    {"|\\operatorname{Tr}\\rho_{AB}-1|","\\|\\rho_{AB}-\\rho_{AB}^\\dagger\\|_F","\\lambda_{\\min}(\\rho_{AB})"}];
  emit["09_density_eigenvalues","\\text{All four edge-state eigenvalues}","\\lambda_k(\\rho_{AB})",
    Table[({#["Time"],#["EigenvaluesAB"][[k]]}& /@ records),{k,4}],
    Table["\\lambda_"<>ToString[k],{k,4}]];
  Do[
    legends={"\\operatorname{"<>part<>"}\\rho_{A;01}","\\operatorname{"<>part<>"}\\rho_{B;01}"};
    If[part=="Abs",legends={"|\\rho_{A;01}|","|\\rho_{B;01}|"}];
    emit["10_local_coherences_"<>part,"\\text{Single-edge coherences: "<>part<>"}",
      "\\rho_{Q;01}(t)",Table[XXZE`MatrixSeries[records,key,1,2,
        Switch[part,"Re",Re,"Im",Im,"Abs",Abs]],{key,{"RhoA","RhoB"}}],legends],
    {part,{"Re","Im","Abs"}}]
];

(* Execution helpers: no nested parallelism. Sweep points are parallel tasks;
   dynamics uses a fixed ED with at most one copied engine per time chunk. *)
XXZE`StartKernels[cfg_] := Module[{available},
  If[!TrueQ[cfg["Enabled"]],Return[0]];
  XXZE`Require[IntegerQ[cfg["Kernels"]] && cfg["Kernels"]>=1,"Kernels must be a positive integer."];
  available=Length[Kernels[]];
  If[available<cfg["Kernels"],Quiet[LaunchKernels[cfg["Kernels"]-available]]];
  available=Min[cfg["Kernels"],Length[Kernels[]]];
  If[available<cfg["Kernels"],Print["Available parallel kernels: ",available,
    " (requested ",cfg["Kernels"],"). Check the local Mathematica kernel/license configuration."]];
  available
];
XXZE`ParallelSweep[grid_,p_,num_,cfg_] := Module[{n,chunks,jobs,res},
  n=Min[Length[grid],XXZE`StartKernels[cfg]];
  If[n<2,Print["Eigenstate sweep: serial."];Return[XXZE`SweepPoint[#,p,num]& /@ grid]];
  DistributeDefinitions[XXZE`SweepPoint];
  chunks=Table[grid[[1+Floor[(k-1) Length[grid]/n];;Floor[k Length[grid]/n]]],{k,n}];
  Print["Eigenstate sweep: ",n," parallel Delta chunks."];
  jobs=Table[With[{gg=g,pp=p,nn=num},ParallelSubmit[
    XXZE`SweepPoint[#,pp,nn]& /@ gg]],{g,chunks}];
  res=WaitAll[jobs];Flatten[res,1]
];
XXZE`ParallelTimes[times_,engine_,rho0_,num_,cfg_] := Module[
  {n,reserve,limit,chunks,jobs,packets,records,reference,error},
  reserve=4. (ByteCount[engine]+ByteCount[rho0])+256. 1024.^2;
  limit=Max[1,Floor[cfg["MemoryBudgetGB"] 1024.^3/reserve]];
  n=Min[Length[times],limit,XXZE`StartKernels[cfg]];
  If[n<2,Print["Time evolution: serial."];
    packets=XXZE`TimeChunk[times,engine,rho0,num["ValidationTolerance"]];
    XXZE`Require[AssociationQ[packets],"Serial time evolution failed."];
    Return[packets["Results"]]];
  DistributeDefinitions[XXZE`TimeChunk];
  chunks=Table[times[[1+Floor[(k-1) Length[times]/n];;Floor[k Length[times]/n]]],{k,n}];
  Print["Time evolution: ",n," parallel chunks; ",Length[times]," samples."];
  jobs=Table[With[{tt=t,ee=engine,rr=rho0,tol=num["ValidationTolerance"]},
    ParallelSubmit[XXZE`TimeChunk[tt,ee,rr,tol]]],{t,chunks}];
  packets=WaitAll[jobs];
  XXZE`Require[AllTrue[packets,AssociationQ],"A parallel time chunk failed."];
  records=Flatten[Lookup[packets,"Results"],1];
  XXZE`Require[Lookup[records,"Time"]===times,"Parallel time ordering failed."];
  reference=XXZE`TimeChunk[{First[times],Last[times]},engine,rho0,num["ValidationTolerance"]];
  XXZE`Require[AssociationQ[reference],"Serial endpoint validation failed."];
  error=Max[MapThread[Norm[#1["RhoAB"]-#2["RhoAB"],"Frobenius"]&,
    {{First[records],Last[records]},reference["Results"]}]];
  XXZE`Require[error<num["ValidationTolerance"],"Parallel/serial endpoint mismatch."];
  records
];



(* ::Section:: *)
(* 3. Choose the model and numerical/plot settings *)


model = <|
  "L" -> 6,                   (* Total spins, including both edges *)
  "J" -> 1,                   (* Signed exchange; -1 is allowed *)
  "h" -> 0,                   (* Uniform longitudinal field *)
  "bLeft" -> 0,               (* Additional field on site 1 *)
  "bRight" -> 0               (* Additional field on site L *)
|>;

deltaFixed = 1/2;              (* Full ED kept for dynamics and eigenstate plots *)
deltaSweep = Range[-5, 5, 1/50];(* Exact grid includes Delta=-1 and +1 *)

numerics = <|
  "Tolerance" -> 10^-10,
  "ValidationTolerance" -> 10^-8,
  "SpectrumMatchingTolerance" -> 10^-8
|>;
parallelSettings = <|"Enabled" -> True,"Kernels" -> 8,"MemoryBudgetGB" -> 20.|>;
plotSettings = <|"ChopTolerance" -> 10^-10,"ExportLogTime" -> False|>;
outputRoot = FileNameJoin[{XXZE`BaseDirectory[],"xxz_open_edges_results"}];
runLabel = "run01";
runSweep = True;
saveData = True;



(* ::Section:: *)
(* 4. Fixed-Delta ED: each sector independently, all eigenpairs retained *)


fixedRun = Catch[Module[{tag,groups,energySets},
  tag=runLabel<>"_L"<>ToString[model["L"]]<>"_J"<>XXZE`Token[model["J"]]<>
    "_h"<>XXZE`Token[model["h"]]<>"_bL"<>XXZE`Token[model["bLeft"]]<>
    "_bR"<>XXZE`Token[model["bRight"]];
  runDirectory=XXZE`NewDirectory[outputRoot,tag];
  {fixedEDSeconds,fixedED}=AbsoluteTiming[XXZE`Diagonalize[deltaFixed,model,numerics]];
  XXZE`Require[AssociationQ[fixedED],"Fixed-Delta diagonalization failed: "<>ToString[fixedED,InputForm]];
  fixedED=XXZE`EigenstateEntropies[fixedED];
  sectorTable=MapIndexed[<|"Sector"->First[#2],"Labels"->#1["Labels"],
    "Dimension"->#1["Dimension"],"Emin"->First[#1["Eigenvalues"]],
    "Emax"->Last[#1["Eigenvalues"]],"Residual"->#1["Residual"]|>&,fixedED["Sectors"]];
  Print["Fixed Delta=",deltaFixed,"; dimension=",fixedED["DimensionSum"],
    "; ED seconds=",fixedEDSeconds];Print[Dataset[sectorTable]];
  energySets=Lookup[fixedED["Sectors"],"Eigenvalues"];
  groups=Gather[Range[Length[energySets]],
    Length[energySets[[#1]]]==Length[energySets[[#2]]] &&
      XXZE`MaxAbs[energySets[[#1]]-energySets[[#2]]]<
        numerics["SpectrumMatchingTolerance"] Max[1.,XXZE`MaxAbs[energySets]]&];
  XXZE`FixedSectorPlots[fixedED,groups,model,plotSettings,runDirectory];
  XXZE`Save[FileNameJoin[{runDirectory,"parameters.txt"}],ToString[
    <|"Model"->model,"DeltaFixed"->deltaFixed,"DeltaSweep"->deltaSweep,
      "Numerics"->numerics,"ParallelSettings"->parallelSettings,"PlotSettings"->plotSettings,
      "Conventions"->"OBC; total L; S=Pauli/2; AB={1,L}; natural logarithms; deterministic degenerate basis",
      "FixedEDSeconds"->fixedEDSeconds,"SectorTable"->sectorTable|>,InputForm],"Text"];
  If[TrueQ[saveData],XXZE`Save[FileNameJoin[{runDirectory,"fixed_sector_eigenpairs.wxf"}],fixedED,"WXF"]];
  <|"Directory"->runDirectory,"SectorGroups"->groups|>
],"XXZFailure"];
If[FailureQ[fixedRun],Print[fixedRun]];



(* ::Section:: *)
(* 5. Delta sweep: spectra and both eigenstate entropies, equivalent sectors together *)


(* Every sweep point computes all eigenpairs in every sector. Eigenvectors are
   used to compute the two entropies and then released; the compact sweep saves
   energies and entropies. Full vectors at deltaFixed are in fixedED/the WXF.
   Get vectors at another Delta with XXZE`Diagonalize[delta,model,numerics].
   At +/-1, total-spin subblocks are merged only for the Delta-axis plots.
   Entropy panels use points because energy-order labels need not track the
   same eigenstate through crossings. Isospectral grouping is checked on the
   chosen grid; it is not a proof of equivalence at all real Delta. *)

sweepRun = If[TrueQ[runSweep] && AssociationQ[fixedRun],Catch[Module[{keys,valid,groups,dir},
  XXZE`Require[Length[deltaSweep]>=2 && VectorQ[deltaSweep,XXZE`RealQ] &&
    And@@Thread[Differences[N[deltaSweep]]>0],"deltaSweep must be a strictly increasing real list."];
  {sweepSeconds,sweepResults}=AbsoluteTiming[
    XXZE`ParallelSweep[deltaSweep,model,numerics,parallelSettings]];
  XXZE`Require[AllTrue[sweepResults,AssociationQ],"A Delta point failed; inspect sweepResults."];
  keys=Keys[First[sweepResults]];
  valid=AllTrue[sweepResults,Sort[Keys[#]]==Sort[keys] &&
    Total[Length[# ["Energies"]]& /@ Values[#]]==2^model["L"]&];
  XXZE`Require[valid,"The sweep has inconsistent sector labels or incomplete dimensions."];
  groups=Gather[keys,XXZE`IsospectralQ[#1,#2,sweepResults,numerics["SpectrumMatchingTolerance"]]&];
  dir=XXZE`NewDirectory[runDirectory,"delta_sweep"];
  XXZE`SweepPlots[sweepResults,deltaSweep,groups,model,plotSettings,dir];
  If[TrueQ[saveData],XXZE`Save[FileNameJoin[{dir,"eigenstate_sweep.wxf"}],
    <|"Model"->model,"DeltaGrid"->deltaSweep,"Seconds"->sweepSeconds,
      "SectorGroups"->groups,"Results"->sweepResults|>,"WXF"]];
  Print["Sweep seconds=",sweepSeconds,"; sectors=",Length[keys],"; panels=",Length[groups]];
  <|"Directory"->dir,"SectorGroups"->groups|>
],"XXZFailure"],Missing["SweepSkipped"]];
If[FailureQ[sweepRun],Print[sweepRun]];



(* ::Section:: *)
(* 6. Choose the initial state and time grid *)


(* ONE active definition per state; alternatives below are commented.
   Mode="Factorized": edge input has dimension 4; middle input 2^(L-2).
   Each input may be a pure vector or a positive trace-one density matrix.
   AB may be internally entangled. Mode="Full": input has dimension 2^L,
   or is a 2^L x 2^L density matrix, in the physical site order {1,...,L}.
   Full mode accepts initially correlated/entangled edges and middle.
   delta and ed are supplied to all generators; no hidden external state cache.
*)

initialStateSettings = <|"Mode" -> "Factorized"|>;
edgeInitial[delta_,p_,ed_] := {1,1,1,1}/2;       (* |+>_1 tensor |+>_L *)
middleInitial[delta_,p_,ed_] := XXZE`NeelKet[p["L"]-2];
fullInitial[delta_,p_,ed_] := XXZE`GroundKet[ed]; (* Used only in Full mode *)

(* Edge alternatives: replace the ONE active edgeInitial definition.
   edgeInitial[delta_,p_,ed_] := {0,1,1,0}/Sqrt[2];   (* Bell pair *)
   edgeInitial[delta_,p_,ed_] := {0,0,1,0};           (* |10> *)
   edgeInitial[delta_,p_,ed_] := IdentityMatrix[4]/4; (* Mixed edges *)

   Middle alternatives:
   middleInitial[delta_,p_,ed_] := XXZE`ProductKet[ConstantArray[{1,0},p["L"]-2]];
   middleInitial[delta_,p_,ed_] := XXZE`HelixKet[p["L"]-2,Pi/3];
   middleInitial[delta_,p_,ed_] := myMiddleState;

   Full-state examples: set initialStateSettings["Mode"]="Full" first.
   fullInitial[delta_,p_,ed_] := XXZE`Eigenvector[ed,3,1];
     (* sector 3, first eigenvector; inspect sectorTable for its labels *)
   fullInitial[delta_,p_,ed_] := XXZE`NeelKet[p["L"]];
   fullInitial[delta_,p_,ed_] := XXZE`HelixKet[p["L"],Pi/3];
   fullInitial[delta_,p_,ed_] := myFullState;
   fullInitial[delta_,p_,ed_] := IdentityMatrix[2^p["L"]]/2^p["L"];
   A normalized eigenvector evolves only by a phase, so all its reduced-state
   observables are stationary. Degenerate ground-state selection is the first
   vector in the first lowest-energy sector; choose another explicitly if needed.
*)

timeSettings = <|
  "Grid" -> "Linear",        (* Linear or Log *)
  "tMax" -> 50.,
  "dt" -> 0.1,
  "tMinPositive" -> 0.01,    (* Used only for Log *)
  "LogPoints" -> 501         (* Positive log points; t=0 is also included *)
|>;



(* ::Section:: *)
(* 7. Evolve at deltaFixed, calculate edge observables, export MaTeX PNGs *)


dynamicsRun = If[AssociationQ[fixedRun],Catch[Module[
  {fac,ef,mf,input,engine,rho0,times,checks,dir,a0,aEnd,initialRecord,purityFull},
  XXZE`Require[fixedED["Parameters"]===model && TrueQ[fixedED["Delta"]==deltaFixed],
    "Model or Delta changed: rerun section 4 before dynamics."];
  times=XXZE`TimeGrid[timeSettings];
  Switch[initialStateSettings["Mode"],
    "Factorized",
      input=<|"Edge"->edgeInitial[deltaFixed,model,fixedED],
        "Middle"->middleInitial[deltaFixed,model,fixedED]|>;
      ef=XXZE`StateFactor[input["Edge"],4,numerics["ValidationTolerance"],"Edges"];
      mf=XXZE`StateFactor[input["Middle"],2^(model["L"]-2),numerics["ValidationTolerance"],"Middle"];
      fac=XXZE`JoinEdgeMiddle[ef,mf,model["L"]],
    "Full",
      input=<|"Full"->fullInitial[deltaFixed,model,fixedED]|>;
      fac=XXZE`StateFactor[input["Full"],2^model["L"],numerics["ValidationTolerance"],"Full system"],
    _,XXZE`Require[False,"Initial mode must be Factorized or Full."]];
  rho0=XXZE`ReducedFactor[fac,{1,model["L"]},model["L"]];
  engine=XXZE`EvolutionEngine[fixedED,fac];
  {dynamicsSeconds,dynamicsResults}=AbsoluteTiming[
    XXZE`ParallelTimes[times,engine,rho0,numerics,parallelSettings]];
  XXZE`Require[ListQ[dynamicsResults] && AllTrue[dynamicsResults,AssociationQ],
    "Time evolution returned an invalid result."];
  a0=XXZE`Evolve[0.,engine];aEnd=XXZE`Evolve[Last[times],engine];
  checks=<|"InitialFactorError"->Norm[a0-fac,"Frobenius"],
    "InitialEdgeError"->Norm[First[dynamicsResults]["RhoAB"]-rho0,"Frobenius"],
    "MaxTraceError"->Max[Lookup[dynamicsResults,"TraceError"]],
    "MaxHermiticityError"->Max[Lookup[dynamicsResults,"HermiticityError"]],
    "MinEigenvalue"->Min[Lookup[dynamicsResults,"MinEigenvalue"]],
    "EndpointNormError"->Abs[Tr[ConjugateTranspose[aEnd] . aEnd]-1.]|>;
  XXZE`Require[Max[Values[KeyDrop[checks,"MinEigenvalue"]]]<numerics["ValidationTolerance"] &&
    checks["MinEigenvalue"]>=-numerics["ValidationTolerance"],"Dynamics validation failed."];
  purityFull=Re[Tr[(ConjugateTranspose[fac] . fac) . (ConjugateTranspose[fac] . fac)]];
  dir=XXZE`NewDirectory[runDirectory,"dynamics_Delta"<>XXZE`Token[deltaFixed]];
  XXZE`Save[FileNameJoin[{dir,"initial_states.wxf"}],
    <|"Mode"->initialStateSettings["Mode"],"Inputs"->input,"FullStateFactor"->fac,"RhoAB0"->rho0|>,"WXF"];
  XXZE`Save[FileNameJoin[{dir,"parameters.txt"}],ToString[
    <|"Model"->model,"Delta"->deltaFixed,"TimeSettings"->timeSettings,"TimeGrid"->times,
      "StateMode"->initialStateSettings["Mode"],"GlobalPurity"->purityFull,
      "StateDefinitions"->{DownValues[edgeInitial],DownValues[middleInitial],DownValues[fullInitial]},
      "Numerics"->numerics,"ParallelSettings"->parallelSettings,"PlotSettings"->plotSettings,
      "Seconds"->dynamicsSeconds,"Checks"->checks,
      "Notes"->"Pauli magnetizations; natural logs; reduced overlap is survival probability only for pure rhoAB0. Reduced entropy is bipartite entanglement only for a pure full state."|>,InputForm],"Text"];
  If[TrueQ[saveData],XXZE`Save[FileNameJoin[{dir,"edge_dynamics.wxf"}],
    <|"Model"->model,"Delta"->deltaFixed,"Checks"->checks,"Results"->dynamicsResults|>,"WXF"]];
  XXZE`DynamicsPlots[dynamicsResults,model,deltaFixed,timeSettings,plotSettings,dir];
  Print["Dynamics seconds=",dynamicsSeconds,"; checks=",checks];Print["Saved: ",dir];
  <|"Directory"->dir,"Checks"->checks,"GlobalPurity"->purityFull|>
],"XXZFailure"],Missing["FixedEDFailed"]];
If[FailureQ[dynamicsRun],Print[dynamicsRun]];



(* ::Section:: *)
(* 8. Access results and reuse definitions *)


(* fixedED["Sectors"][[k]]["Labels"]
   fixedED["Sectors"][[k]]["Eigenvalues"]
   XXZE`Eigenvector[fixedED,k,r]                 (* Full 2^L vector *)
   fixedED["Sectors"][[k]]["HalfEntropy"]
   fixedED["Sectors"][[k]]["EdgesEntropy"]
   dynamicsResults[[q]]["RhoAB"] // MatrixForm
   dynamicsResults[[q]]["RhoA"] // MatrixForm
   dynamicsResults[[q]]["RhoB"] // MatrixForm
   dynamicsResults[[q]]["EigenvaluesAB"]
   XXZE`EdgesDensity[XXZE`Eigenvector[fixedED,k,r],model["L"]]

   Fresh full ED at another anisotropy, without changing plotting/state settings:
   otherED=XXZE`Diagonalize[1,model,numerics];
   otherED=XXZE`EigenstateEntropies[otherED];
   All functions are globally accessible under XXZE`; no Begin/End required.
   For another initial state at the SAME model, edit section 6 and rerun 7.
   For another model or deltaFixed, rerun 4 and then 5/7 as needed.
*)
