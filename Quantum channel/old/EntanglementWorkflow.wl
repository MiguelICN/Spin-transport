(* ::Package:: *)

BeginPackage["EntanglementWorkflow`"]

(*Definitions to compute entanglement entropy fast MACHINE PRECISION ONLY*)
ExpandBlockVector::usage = "ExpandBlockVector[v, map, type, l] expands a projected state back to the full Hilbert space."
revIndex::usage = "revIndex[i, l] returns the reverse-binary index of i with L bits."
ProjectBlockVector::usage = "ProjectBlockVector[state, map, type, l] projects a full state into parity block."
ExpandFullState::usage = "ExpandFullState[evenState, oddState, evenMap, oddMap, l] reconstructs full state from even/odd sectors."
makeExpansionMatrix::usage = "makeExpansionMatrix[map, l, type] returns the sparse matrix that maps block basis to full basis."
expandFullStateFastTemplate::usage = "Template for fast compiled expander."
VNentropy::usage = "von Neumann entanglement entropy."
new::usage = "new[colVec, l] returns the entanglement entropy of a bipartite system."
compiledEvolverSector::usage = "compiledEvolverSector[eigenvec, eigenval] returns a compiled time evolution function."
buildMaps::usage = "buildMaps[l] generates the even and odd parity index maps."
buildExpandFullStateFast::usage = "Builds a compiled fast expander for use in entropy evaluation."


VNentropy::usage "Von Neumann entropy."
sdiag::usage "Diagonal entropy of Eigenvectors."
pageEntropy::usage "Page value for tensor product Hilbert spaces."
haar::usage "Haar state of dimension 2^L."
pur::usage "Purity"
decoherence::usage "Total decoherence of a purity curve"


Begin["`Private`"]

Needs["Developer`"]

ClearAll[ExpandBlockVector]
ExpandBlockVector[v_, map_, type_, l_] := 
 Module[{full = ConstantArray[0, 2^l], k = 1, a, b, coeff},
  Do[
   If[Length[pair] == 1,
    full[[pair[[1]] + 1]] = v[[k]],
    {a, b} = pair;
    coeff = v[[k]]/Sqrt[2];
    If[type === "even",
     full[[a + 1]] = coeff; full[[b + 1]] = coeff,
     full[[a + 1]] = coeff; full[[b + 1]] = -coeff]
   ];
   k++,
   {pair, map}];
  full
 ]

ClearAll[revIndex]
revIndex[i_, l_] := FromDigits[Reverse[IntegerDigits[i, 2, l]], 2]

ClearAll[ProjectBlockVector]
ProjectBlockVector[state_, map_, type_, l_] := 
 Module[{res, k = 1, a, b},
  res = ConstantArray[0, Length[map]];
  Do[
   If[Length[pair] == 1,
    res[[k]] = state[[pair[[1]] + 1]],
    {a, b} = pair;
    res[[k]] = If[type === "even",
      (state[[a + 1]] + state[[b + 1]])/Sqrt[2],
      (state[[a + 1]] - state[[b + 1]])/Sqrt[2]]
   ];
   k++,
   {pair, map}];
  res
 ]

ClearAll[ExpandFullState]
ExpandFullState[evenState_, oddState_, evenMap_, oddMap_, l_] := 
 ExpandBlockVector[evenState, evenMap, "even", l] + 
  ExpandBlockVector[oddState, oddMap, "odd", l]

ClearAll[makeExpansionMatrix]
makeExpansionMatrix[map_, l_, type_] := 
 SparseArray[
  Flatten[
   Table[
    If[Length[pair] == 1,
     {{pair[[1]] + 1, k} -> 1},
     Which[
      type === "even",
      {{pair[[1]] + 1, k} -> 1/Sqrt[2], {pair[[2]] + 1, k} -> 1/Sqrt[2]},
      type === "odd",
      {{pair[[1]] + 1, k} -> 1/Sqrt[2], {pair[[2]] + 1, k} -> -1/Sqrt[2]}
     ]
    ],
    {k, Length[map]}, {pair, {map[[k]]}}], 2],
  {2^l, Length[map]}]

ClearAll[expandFullStateFastTemplate]
expandFullStateFastTemplate[evenMat_, oddMat_][evenState_, oddState_] := 
 evenMat . evenState + oddMat . oddState

ClearAll[buildExpandFullStateFast]
buildExpandFullStateFast[evenMap_, oddMap_, l_] := 
 Module[{evenMat, oddMat},
  evenMat = makeExpansionMatrix[evenMap, l, "even"];
  oddMat = makeExpansionMatrix[oddMap, l, "odd"];
  Compile[{{evenState, _Complex, 1}, {oddState, _Complex, 1}},
   expandFullStateFastTemplate[evenMat, oddMat][evenState, oddState],
   CompilationTarget -> "WVM",
   RuntimeOptions -> {"EvaluateSymbolically" -> False}]
 ]

ClearAll[VNentropy]
VNentropy[0 | 0.] = 0;
VNentropy[x_] := -x Log[x]

ClearAll[new]
new[colVec_, l_] := 
 Module[{m},
  m = ArrayReshape[colVec, {2^(l/2), 2^(l/2)}];
  Total[Map[VNentropy, Chop[SingularValueList[m]]^2]]
 ]

ClearAll[compiledEvolverSector]
compiledEvolverSector[eigenvec_, eigenval_] := 
 With[{evT = Transpose[eigenvec], evC = Conjugate[eigenvec]},
  Compile[{{state, _Complex, 1}, {tlist, _Real, 1}},
   Module[{ck},
    ck = evC . state;
    Table[evT . (ck Exp[-I eigenval t]), {t, tlist}]
   ],
   CompilationTarget -> "WVM",
   RuntimeOptions -> {"EvaluateSymbolically" -> False}]
 ]

ClearAll[buildMaps]
buildMaps[l_] := 
 Module[{visited, evenMap = {}, oddMap = {}, j},
  visited = ConstantArray[False, 2^l, SparseArray];
  Do[
   If[!visited[[i + 1]],
    j = revIndex[i, l];
    visited[[i + 1]] = True;
    visited[[j + 1]] = True;
    If[i == j,
     AppendTo[evenMap, {i}],
     AppendTo[evenMap, {i, j}]; AppendTo[oddMap, {i, j}]
    ]
   ],
   {i, 0, (2^l) - 1}];
  {evenMap, oddMap}
 ]


sdiag[L_,La_,Lb_,ini_]:=Module[{state=Dyad[ini],dima=2^(La),dimb=2^(Lb),rho},
rho=MatrixPartialTrace[state,2,{dima,dimb}];
Chop[Total[Map[VNentropy,Chop[Eigenvalues[rho]]]]]
];

pageEntropy[La_,Lb_]:=PolyGamma[0,(2^La)(2^Lb)+1]-PolyGamma[0,(2^Lb)+1]-((2^La)-1)/(2*(2^Lb));

haar[L_Integer]:=Normalize[Table[RandomVariate[NormalDistribution[0,1]]+I RandomVariate[NormalDistribution[0,1]],{2^L}]];

pur[L_,La_,Lb_,ini_,eigenval_,eigenvec_,t_]:=Module[{state=StateEvolution[t,ini,eigenval,eigenvec],rho,\[Lambda],purity,dima=2^(La),dimb=2^(Lb),vn},
rho=MatrixPartialTrace[Dyad[state],2,{dima,dimb}];
Tr[Chop[rho . rho]]
];

(*Compute total decoherence via trapezoidal rule*)
decoherence[purity_,tlist_]:=Module[{pur=Chop[purity]},(1/tlist[[-1]])Total[Table[(pur[[i]]+pur[[i+1]])/2*(tlist[[i+1]]-tlist[[i]]),{i,1,Length[purity]-1}]]]


End[]
EndPackage[]
