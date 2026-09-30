(* ::Package:: *)
(* ::Title:: *)
(*XXZ boundary probes - complete dynamical map*)

(* Edit this file, then evaluate Get["/your/path/xxz_edge_dynamical_map.m"].
   To rerun individual definition/input/execution cells after loading, first use
   Begin["EdgeChannel`Private`"]; and finish with End[].
   The public result functions below remain accessible after Get finishes. *)

BeginPackage["EdgeChannel`"];
ClearAll["EdgeChannel`*","EdgeChannel`Private`*"];
ChoiZ::usage = "ChoiZ: see editable sections and examples in this file.";
SuperoperatorZ::usage = "SuperoperatorZ: see editable sections and examples in this file.";
ApplyChannelZ::usage = "ApplyChannelZ: see editable sections and examples in this file.";
ChoiExchange::usage = "ChoiExchange: see editable sections and examples in this file.";
SuperoperatorExchange::usage = "SuperoperatorExchange: see editable sections and examples in this file.";
ApplyChannelExchange::usage = "ApplyChannelExchange: see editable sections and examples in this file.";
channelResultsZ::usage = "channelResultsZ: see editable sections and examples in this file.";
channelResultsExchange::usage = "channelResultsExchange: see editable sections and examples in this file.";
Begin["`Private`"];

(* ::Section:: *)
(*Conventions and workflow*)
(* L is the number of CHAIN spins. Total size is L+2, ordered A,1,...,L,B.
   Probe basis: {00,01,10,11}; Z|0> = +|0>, Z|1> = -|1>.
   hbar=1; S^alpha = sigma^alpha/2 on the chain.
   HC[Delta] = Jxy/4 Sum[XX+YY+Delta ZZ] + w/2 Sum[Z] + ed/2 Z_d.
   Thus Jz=Jxy Delta; Delta itself is never divided by 10.
   HA=omegaA ZA/2, HB=omegaB ZB/2.
   Z interaction: gA ZA Z1/4 + gB ZL ZB/4
                = (gA/2) ZA S1z + (gB/2) ZB SLz.
   Exchange: gA(sigmaA+ S1- + sigmaA- S1+), likewise at B,
             with sigma+ = (X+iY)/2 and S+ = (X+iY)/2.
   OLD FILE CONVERSION: gExchange=4 lambdaOld; omega=2 DeltaProbeOld.
   Initial state ALWAYS factorizes between AB and C: rhoAB0 (x) rhoC0.
   AB may be internally entangled. Chain vectors may be arbitrary entangled
   pure states; nonthermal density matrices are also accepted.
   No master equation, temperature, weak-coupling or Markov approximation.

   Each file is independent: no QMB library or C compiler is needed.
   Definitions -> editable parameters/state -> execution -> PNG exports.
   Every named function takes delta explicitly; numerical results retain it.
   Z and Exchange inputs are separate; set either Run...=False to skip it.
   Input functions are evaluated anew at EACH Delta, including ground states.
   Results retain only small reduced matrices, not full evolving wavefunctions.
   Diagonalization uses total-Z sectors; no large matrices copied to subkernels.
   Exact diagonalization still grows exponentially: start at L=6 or L=8.
   Figures are raw sampled trajectories, joined by lines; no smoothing.
   All logarithms in entropies are natural (nats).
*)

(* ::Section:: *)
(*1 - Numerical utilities; validation errors stop the current run*)

RunAbort[delta_, message_] := Throw[Failure["InvalidRun", <|"Delta" -> delta,
  "Message" -> message|>], "RunFailure"];
Require[delta_, test_, message_] := If[!TrueQ[test], RunAbort[delta, message]];
RealScalarQ[delta_, x_] := NumberQ[N[x]] && TrueQ[Im[N[x]] == 0];
SparseId[delta_, n_] := SparseArray[Band[{1, 1}] -> 1., {n, n}];
Pauli[delta_, k_] := N[PauliMatrix[k]];
SiteOp[delta_, op_, j_, n_] := KroneckerProduct[
  SparseId[delta, 2^(j-1)], SparseArray[op], SparseId[delta, 2^(n-j)]];
BondOp[delta_, a_, b_, j_, n_] := KroneckerProduct[
  SparseId[delta, 2^(j-1)], SparseArray[a], SparseArray[b], SparseId[delta, 2^(n-j-1)]];
KetDensity[delta_, v_] := Outer[Times, v, Conjugate[v]];
HermitianPart[delta_, m_] := (m + ConjugateTranspose[m])/2;
VNEntropy[delta_, m_] := Module[{ev},
  ev = Select[Re[Eigenvalues[HermitianPart[delta, m]]], # > 0. &];
  If[ev === {}, 0., -Total[ev Log[ev]]]
];
TraceOutB[delta_, r_] := {{r[[1,1]]+r[[2,2]], r[[1,3]]+r[[2,4]]},
  {r[[3,1]]+r[[4,2]], r[[3,3]]+r[[4,4]]}};
TraceOutA[delta_, r_] := {{r[[1,1]]+r[[3,3]], r[[1,2]]+r[[3,4]]},
  {r[[2,1]]+r[[4,3]], r[[2,2]]+r[[4,4]]}};
PartialTransposeB[delta_, r_] := ArrayReshape[
  Transpose[ArrayReshape[r, {2,2,2,2}], {1,4,3,2}], {4,4}];
Negativity[delta_, r_] := Total[Clip[-Re[Eigenvalues[
  HermitianPart[delta, PartialTransposeB[delta, r]]]], {0., Infinity}]];
Concurrence[delta_, r_] := Module[{v, rows, root, flip, ev},
  {v, rows} = Eigensystem[HermitianPart[delta, r]];
  root = Transpose[rows].(Sqrt[Clip[Re[v], {0., Infinity}]] Conjugate[rows]);
  flip = KroneckerProduct[Pauli[delta, 2], Pauli[delta, 2]];
  ev = Reverse[Sort[Sqrt[Clip[Re[Eigenvalues[HermitianPart[delta,
    root.flip.Conjugate[r].flip.root]]], {0., Infinity}]]]];
  Max[0., ev[[1]]-Total[Rest[ev]]]
];

(* Density factor R obeys rho=R.R^dagger. Vectors are normalized explicitly;
   matrices must already have unit trace and be positive within ValidationTolerance.
   Only negative roundoff eigenvalues are clipped; all positive weights retained.
   For a matrix, numerical clipping/renormalization is reported in the run metadata. *)
StateFactor[delta_, input_, dim_, tol_, label_] := Module[{x=N[Normal[input]], ev, rows, keep, fac, rho},
  If[VectorQ[x, NumericQ],
    Require[delta, Length[x] == dim && Norm[x] > 0., label <> ": invalid vector dimension or zero vector"];
    fac = Transpose[{x/Norm[x]}];
    Return[<|"Factor"->fac, "InputTraceOrNorm"->Norm[x], "CorrectionNorm"->0.|>]
  ];
  Require[delta, MatrixQ[x, NumericQ] && Dimensions[x] == {dim,dim}, label <> ": expected vector or square density matrix"];
  Require[delta, Norm[x-ConjugateTranspose[x], "Frobenius"] <= tol && Abs[Tr[x]-1.] <= tol,
    label <> ": matrix must be Hermitian with trace 1"];
  {ev, rows} = Eigensystem[HermitianPart[delta, x]];
  Require[delta, Min[Re[ev]] >= -tol, label <> ": density matrix is not positive"];
  ev = Clip[Re[ev], {0., Infinity}]; ev = ev/Total[ev];
  keep = Flatten[Position[ev, _?(# > 0. &)]];
  fac = Transpose[Sqrt[ev[[keep]]] rows[[keep]]];
  rho = fac.ConjugateTranspose[fac];
  <|"Factor"->fac, "InputTraceOrNorm"->Tr[x], "CorrectionNorm"->Norm[rho-x,"Frobenius"]|>
];

(* ::Section:: *)
(*2 - Hamiltonians and initial chain states*)

ChainHamiltonian[delta_?NumericQ, p_Association] := Module[{n=p["L"], h, x,y,z},
  {x,y,z}=Table[Pauli[delta,k],{k,3}];
  h=SparseArray[{}, {2^n,2^n}];
  Do[h += p["Jxy"]/4 (BondOp[delta,x,x,j,n]+BondOp[delta,y,y,j,n]+
    delta BondOp[delta,z,z,j,n]),{j,n-1}];
  Do[h += p["w"]/2 SiteOp[delta,z,j,n],{j,n}];
  h += p["ed"]/2 SiteOp[delta,z,p["DefectSite"],n];
  Re[h]
];
FullHamiltonian[delta_?NumericQ, model_, p_Association] := Module[{n=p["L"]+2,h,x,y,z},
  {x,y,z}=Table[Pauli[delta,k],{k,3}];
  h=KroneckerProduct[SparseId[delta,2],ChainHamiltonian[delta,p],SparseId[delta,2]]+
    p["omegaA"]/2 SiteOp[delta,z,1,n]+p["omegaB"]/2 SiteOp[delta,z,n,n];
  If[model=="Z",
    h += p["gA"]/4 BondOp[delta,z,z,1,n]+p["gB"]/4 BondOp[delta,z,z,n-1,n],
    h += p["gA"]/2 (BondOp[delta,x,x,1,n]+BondOp[delta,y,y,1,n])+
         p["gB"]/2 (BondOp[delta,x,x,n-1,n]+BondOp[delta,y,y,n-1,n])
  ]; Re[h]
];
ProductKet[delta_, localKets_List] := Module[{v},
  Require[delta, Length[localKets]>0 && And@@(VectorQ[#,NumericQ] && Length[#]==2 && Norm[#]>0 & /@ localKets),
    "Each product-state factor must be a nonzero two-component vector"];
  v=Fold[Flatten[KroneckerProduct[#1,#2/Norm[#2]]] &, {1.}, N[localKets]];
  v/Norm[v]
];
NeelKet[delta_, n_Integer] := ProductKet[delta, Table[If[OddQ[j],{1.,0.},{0.,1.}],{j,n}]];

(* Number of |1> spins labels conserved magnetization sectors. The ground-state
   routine searches ALL sectors, including fully polarized sectors. For degenerate
   ground spaces it projects a fixed computational basis vector into that space:
   choose smallest number of |1> spins first, then smallest nonzero projection.
   GroundTolerance defines numerical degeneracy. Override ChainState... to choose
   another vector/symmetry sector or a coherent superposition in that manifold. *)
SectorIndices[delta_, n_] := GatherBy[Range[2^n], Total[IntegerDigits[#-1,2,n]] &];
SortedEigensystem[delta_, h_] := Module[{e,v,order},
  {e,v}=Eigensystem[Normal[N[h]]]; order=Ordering[Re[e]];
  {Re[e[[order]]],v[[order]]}
];
GroundKet[delta_?NumericQ, p_Association] := Module[
  {h=ChainHamiltonian[delta,p],groups,records,emin,candidates,chosen,rows,weights,k,v,result,tol},
  groups=SortBy[SectorIndices[delta,p["L"]], Total[IntegerDigits[First[#]-1,2,p["L"]]] &];
  records=Table[With[{es=SortedEigensystem[delta,h[[ix,ix]]]},
    <|"Indices"->ix,"E"->es[[1]],"Rows"->es[[2]]|>],{ix,groups}];
  emin=Min[Flatten[(#["E"] & /@ records)]];
  tol=p["GroundTolerance"] Max[1.,Abs[emin]];
  candidates=Select[records,First[#["E"]]<=emin+tol &]; chosen=First[candidates];
  rows=Pick[chosen["Rows"], (#<=emin+tol & /@ chosen["E"])];
  weights=Total[Abs[rows]^2]; k=First[FirstPosition[weights,_?(#>p["GroundTolerance"] &)]];
  v=Transpose[rows].Conjugate[rows[[All,k]]]; v=v/Norm[v];
  result=ConstantArray[0.+0.I,2^p["L"]]; result[[chosen["Indices"]]]=v; result
];

(* ::Section:: *)
(*3 - State-vector engines, diagonalized once per Delta*)

(* Preparing four input columns per chain-state component yields an isometry.
   An arbitrary two-probe input is applied afterwards, without another evolution.
   For dephasing only four conditional CHAIN Hamiltonians are needed. *)
BuildInitialVector[delta_, vAB_, vC_, n_] := Flatten[Transpose[
  TensorProduct[ArrayReshape[vAB,{2,2}],vC],{1,3,2}]];
BuildEngine[delta_?NumericQ, model_, p_, chainFactor_] := Module[
  {n=p["L"],r=Dimensions[chainFactor][[2]],signs,hc,blocks,h,initial,groups},
  signs={{1,1},{1,-1},{-1,1},{-1,-1}};
  If[model=="Z",
    hc=ChainHamiltonian[delta,p]; groups=SectorIndices[delta,n];
    blocks=Table[
      h=hc+signs[[s,1]] p["gA"]/4 SiteOp[delta,Pauli[delta,3],1,n]+
        signs[[s,2]] p["gB"]/4 SiteOp[delta,Pauli[delta,3],n,n];
      PrepareBlocks[delta,h,chainFactor,groups],{s,4}];
    <|"Delta"->delta,"Model"->model,"L"->n,"RankC"->r,"Blocks"->blocks,
      "LocalEnergies"->(signs.{p["omegaA"],p["omegaB"]}/2)|>,
    h=FullHamiltonian[delta,model,p];
    initial=Transpose[Flatten[Table[BuildInitialVector[delta,UnitVector[4,s],
      chainFactor[[All,k]],n],{s,4},{k,r}],1]];
    blocks=PrepareBlocks[delta,h,initial,SectorIndices[delta,n+2]];
    <|"Delta"->delta,"Model"->model,"L"->n,"RankC"->r,"Blocks"->blocks|>
  ]
];
PrepareBlocks[delta_, h_, initial_, groups_] := Module[{data,es,c},
  data=Reap[Do[
    (* Skip only exactly zero initial blocks; no small-amplitude truncation. *)
    If[Max[Abs[Flatten[initial[[ix,All]]]]]>0.,
      es=SortedEigensystem[delta,h[[ix,ix]]];
      c=Conjugate[es[[2]]].initial[[ix,All]];
      Sow[<|"Indices"->ix,"E"->es[[1]],"V"->Transpose[es[[2]]],"C"->c|>]
    ],{ix,groups}]][[2]];
  <|"Dimension"->Length[h],"Columns"->Dimensions[initial][[2]],
    "Sectors"->If[data==={},{},First[data]]|>
];
EvolveBlocks[delta_, t_?NumericQ, data_] := Module[{out},
  out=ConstantArray[0.+0.I,{data["Dimension"],data["Columns"]}];
  Do[out[[b["Indices"],All]]=b["V"].(Exp[-I b["E"] t] b["C"]),{b,data["Sectors"]}]; out
];
DephasingFactors[delta_?NumericQ,t_?NumericQ,engine_] := Module[{vectors},
  Require[delta, engine["Delta"]==delta && engine["Model"]=="Z", "Engine/Delta/model mismatch"];
  vectors=Table[Exp[-I engine["LocalEnergies"][[s]] t]
    Flatten[EvolveBlocks[delta,t,engine["Blocks"][[s]]]],{s,4}];
  vectors.ConjugateTranspose[vectors]
];
(* Return four 4 x (dC rankC) matrices W_i, one for each AB input |i>.
   E(|i><j|)=W_i.W_j^dagger. Chain columns include sqrt(probability). *)
OutputBlocks[delta_?NumericQ,t_?NumericQ,engine_] := Module[{out,n=engine["L"],r=engine["RankC"]},
  Require[delta, engine["Delta"]==delta, "Engine/Delta mismatch"];
  out=EvolveBlocks[delta,t,engine["Blocks"]];
  Table[ArrayReshape[Transpose[ArrayReshape[out[[All,(s-1) r+1;;s r]],
    {2,2^n,2,r}],{1,3,2,4}],{4,2^n r}],{s,4}]
];

(* ::Section:: *)
(*4 - Time grids, reproducible run records and PNG labels*)

ValidateParameters[delta_, p_] := Module[{},
  Require[delta, RealScalarQ[delta,delta], "Delta must be a real number"];
  Require[delta, IntegerQ[p["L"]] && p["L"]>=1, "L must be a positive integer"];
  Require[delta, IntegerQ[p["DefectSite"]] && 1<=p["DefectSite"]<=p["L"], "DefectSite must lie in 1..L"];
  Require[delta, And@@(RealScalarQ[delta,p[#]] & /@ {"Jxy","w","ed","omegaA","omegaB","gA","gB","tMax","dt"}),
    "Hamiltonian parameters and times must be real"];
  Require[delta,p["tMax"]>0 && p["dt"]>0,"tMax and dt must be positive"];
  Require[delta,p["ValidationTolerance"]>0 && p["GroundTolerance"]>0 && p["InverseTolerance"]>0,
    "All tolerances must be positive"];
];
TimeGrid[delta_,p_] := Module[{t=N[Range[0,Floor[p["tMax"]/p["dt"]],1] p["dt"]]},
  If[Abs[Last[t]-p["tMax"]]<=10^-12 Max[1.,p["tMax"]],
    t[[-1]]=N[p["tMax"]], t=Append[t,N[p["tMax"]]]]; t
];
BaseDirectory[delta_] := Module[{dir},
  dir=If[StringQ[$InputFileName] && StringLength[$InputFileName]>0,
    DirectoryName[ExpandFileName[$InputFileName]], Quiet[Check[NotebookDirectory[],Directory[]]]];
  If[StringQ[dir],dir,Directory[]]
];
SafeName[delta_, x_] := StringReplace[ToString[x,InputForm],
  Except[LetterCharacter|DigitCharacter|"-"|"_"] -> "p"];
CreateRunDirectory[delta_,p_,model_,kind_] := Module[{dir},
  dir=FileNameJoin[{p["OutputRoot"],kind<>"_"<>model,
    "L"<>ToString[p["L"]]<>"_D"<>SafeName[delta,delta]<>"_"<>StringTake[CreateUUID[],8]}];
  Require[delta,StringQ[Quiet[Check[CreateDirectory[dir,CreateIntermediateDirectories->True],$Failed]]],
    "Could not create output directory: "<>dir]; dir
];
SaveFile[delta_,file_,data_,format_] := Module[{result},
  result=Quiet[Check[Export[file,data,format],$Failed]];
  Require[delta,StringQ[result],"Export failed: "<>file]; result
];
(* The exact state arrays and state-function definitions go in initial_states.wxf
   and parameters.txt. Every PNG identifies them by SHA256 plus a human label.
   Large arbitrary states cannot reasonably be printed component by component
   on a figure. Their hash names the exact saved input, including all phases. *)
StateHash[delta_, x_] := IntegerString[Hash[N[x],"SHA256"],16,64];
PlotCaption[delta_,p_,model_,stateInfo_] := Column[{
  Row[{"XXZ OBC; A-1-...-L-B; L=",p["L"],"; Ntotal=",p["L"]+2,
    "; Delta=",delta,"; Jxy=",p["Jxy"],"; Jz=",p["Jxy"] delta}],
  Row[{"w=",p["w"],"; ed=",p["ed"],"; defect site=",p["DefectSite"],
    "; omegaA=",p["omegaA"],"; omegaB=",p["omegaB"],"; gA=",p["gA"],"; gB=",p["gB"]}],
  If[model=="Z","Hint=(gA ZA Z1 + gB ZL ZB)/4",
    "Hint=gA(sigmaA+ S1- + sigmaA- S1+) + gB(SL+ sigmaB- + SL- sigmaB+)"],
  "HC=Jxy Sum(XX+YY+Delta ZZ)/4 + w Sum(Z)/2 + ed Z_d/2; HQ=omegaQ ZQ/2; hbar=1",
  Row[{"t=0..",p["tMax"],"; dt=",p["dt"]," (final step shortened if needed); raw samples; nats"}],
  Row[{"validation tol=",p["ValidationTolerance"],"; ground tol=",p["GroundTolerance"],
    "; inverse tol=",p["InverseTolerance"],"; Sz-sector diagonalization"}],
  Row[{"Chain: ",p["ChainStateLabel"],"; ",stateInfo["ChainSummary"]}],
  "chain SHA256: "<>stateInfo["ChainHash"],
  stateInfo["ProbeSummary"],stateInfo["ProbeHashLine"]
},Alignment->Center,Spacings->.25];
ExportPlot[delta_,p_,model_,info_,dir_,name_,series_,legends_,quantity_] := Module[{plot,panel},
  plot=ListLinePlot[series,Frame->True,Axes->False,PlotRange->All,
    FrameLabel->{"Time t (hbar=1)",quantity},LabelStyle->Directive[Black,14],
    PlotLegends->Placed[LineLegend[legends],Right],ImageSize->1100,AspectRatio->.48];
  panel=Column[{Style[quantity,16,Bold],Style[PlotCaption[delta,p,model,info],11],plot},
    Alignment->Center,Spacings->1];
  SaveFile[delta,FileNameJoin[{dir,name<>".png"}],panel,"PNG"]
];
ExtractSeries[delta_,records_,key_] := ({#["Time"],#[key]} & /@ records);
MatrixSeries[delta_,records_,key_,i_,j_,part_] := ({#["Time"],part[#[key][[i,j]]]} & /@ records);
DataAt[delta_?NumericQ,t_?NumericQ,runs_,key_] := Module[{run,record},
  run=Select[runs,AssociationQ[#] && TrueQ[#["Delta"]==delta] &];
  If[Length[run]!=1,Return[Missing["DeltaNotFound",delta]]];
  record=Select[First[run]["Results"],Abs[#["Time"]-t]<=10^-12 Max[1.,Abs[t]] &];
  If[Length[record]!=1,Missing["TimeNotOnStoredGrid",t],First[record][key]]
];

(* ::Section:: *)
(*5 - Complete dynamical map: normalized Choi state and row-vectorized superoperator*)
(* chi = (1/4) Sum_ij |i><j|_in (x) E(|i><j|)_out, ordered INPUT then OUTPUT.
   Tr chi=1; Tr_out chi=I_4/4. E(rho)=4 Tr_in[(rho^T (x) I) chi].
   Row-vectorization: Flatten[E(rho)] = Superoperator.Flatten[rho].
   Choi purity/entropy characterize the map, not entanglement between probes.
   Process MI uses the bipartition (Ain,Aout)|(Bin,Bout).
   A nonzero product distance detects a nonproduct channel; it does NOT alone
   establish entangling ability, non-Markovianity, or a transport coefficient.
*)
ChoiFromFactors[delta_,f_] := Module[{chi=ConstantArray[0.+0.I,{16,16}],ix},
  ix={1,6,11,16}; chi[[ix,ix]]=f/4; chi
];
ChoiMatrix[delta_?NumericQ,t_?NumericQ,engine_] := Module[{w},
  If[engine["Model"]=="Z",ChoiFromFactors[delta,DephasingFactors[delta,t,engine]],
    w=OutputBlocks[delta,t,engine];
    ArrayFlatten[Table[w[[i]].ConjugateTranspose[w[[j]]],{i,4},{j,4}]]/4]
];
ChoiToSuper[delta_,chi_] := ArrayReshape[Table[
  4 chi[[4(i-1)+a,4(j-1)+b]],{a,4},{b,4},{i,4},{j,4}],{16,16}];
SuperToChoi[delta_,s_] := ArrayReshape[Table[
  s[[4(a-1)+b,4(i-1)+j]]/4,{i,4},{a,4},{j,4},{b,4}],{16,16}];
ApplySuper[delta_,s_,rho_] := ArrayReshape[s.Flatten[rho],{4,4}];
DynamicalMap[delta_?NumericQ,t_?NumericQ,engine_,rho_] :=
  ApplySuper[delta,ChoiToSuper[delta,ChoiMatrix[delta,t,engine]],rho];
TraceChoiOutput[delta_,chi_] := Table[Sum[chi[[4(i-1)+a,4(j-1)+a]],{a,4}],{i,4},{j,4}];
ProcessMarginals[delta_,chi_] := Module[{q,a,b},
  q=ArrayReshape[Transpose[ArrayReshape[chi,{2,2,2,2,2,2,2,2}],
    {1,3,2,4,5,7,6,8}],{4,4,4,4}];
  a=Table[Sum[q[[i,k,j,k]],{k,4}],{i,4},{j,4}];
  b=Table[Sum[q[[k,i,k,j]],{k,4}],{i,4},{j,4}];
  {a,b,ArrayReshape[q,{16,16}]}
];
ChannelObservables[delta_,t_,chi_,s_,p_] := Module[{ev,a,b,q,phi,transfer},
  ev=Reverse[Sort[Re[Eigenvalues[HermitianPart[delta,chi]]]]];
  {a,b,q}=ProcessMarginals[delta,chi]; phi=Flatten[IdentityMatrix[4]]/2;
  transfer=Table[Re[Tr[KroneckerProduct[IdentityMatrix[2],Pauli[delta,k]].
    ApplySuper[delta,s,KroneckerProduct[Pauli[delta,k],IdentityMatrix[2]]]]]/4,{k,3}];
  <|"Delta"->delta,"Time"->t,"Choi"->chi,"Superoperator"->s,"Eigenvalues"->ev,
    "ChoiPurity"->Re[Tr[chi.chi]],"ChoiEntropy"->VNEntropy[delta,chi],
    "ChoiRank"->Count[ev,_?(#>p["ValidationTolerance"] &)],
    "ProcessMutualInfo"->VNEntropy[delta,a]+VNEntropy[delta,b]-VNEntropy[delta,chi],
    "ProductDistance"->Norm[q-KroneckerProduct[a,b],"Frobenius"],
    "IdentityProcessFidelity"->Re[Conjugate[phi].chi.phi],
    "TransferX"->transfer[[1]],"TransferY"->transfer[[2]],"TransferZ"->transfer[[3]],
    "TraceError"->Abs[Tr[chi]-1.],"MinEigenvalue"->Last[ev],
    "HermiticityError"->Norm[chi-ConjugateTranspose[chi],"Frobenius"],
    "TPError"->Norm[TraceChoiOutput[delta,chi]-IdentityMatrix[4]/4,"Frobenius"],
    "UnitalError"->Norm[ApplySuper[delta,s,IdentityMatrix[4]]-IdentityMatrix[4],"Frobenius"]|>
];
(* Adjacent sampled intermediate maps V(t,s)=E(t) E(s)^(-1).
   A negative minimum Choi eigenvalue, beyond numerical error, witnesses failure
   of CP divisibility. Positive values on this grid do not prove divisibility at
   all times. If s is singular/ill-conditioned the test is UNRESOLVED, not zero.
   No pseudoinverse is substituted for a nonunique intermediate map. *)
IntermediateChecks[delta_,res_,p_] := Module[{prev,sv,ratio,v,ev,check},
  MapIndexed[Function[{rec,index},
    If[First[index]==1 || !TrueQ[p["CheckDivisibility"]],
      Join[rec,<|"IntermediateMinEigenvalue"->Missing["NotTested"],
        "IntermediateConditionRatio"->Missing["NotTested"]|>],
      prev=res[[First[index]-1]]["Superoperator"];
      sv=SingularValueList[prev]; ratio=Last[sv]/First[sv];
      If[ratio<=p["InverseTolerance"],
        ev=Missing["SingularOrIllConditioned"],
        v=Transpose[LinearSolve[Transpose[prev],Transpose[rec["Superoperator"]]]];
        ev=Min[Re[Eigenvalues[HermitianPart[delta,SuperToChoi[delta,v]]]]]
      ];
      Join[rec,<|"IntermediateMinEigenvalue"->ev,"IntermediateConditionRatio"->ratio|>]
    ]],res]
];
ChannelPlots[delta_,p_,model_,info_,dir_,res_] := Module[{pairs,basis,valid},
  ExportPlot[delta,p,model,info,dir,"01_choi_purity_and_identity_fidelity",
    (ExtractSeries[delta,res,#] & /@ {"ChoiPurity","IdentityProcessFidelity"}),
    {"Tr[chi(t)^2]","<Phi4|chi(t)|Phi4>, |Phi4>=Sum_i |ii>/2"},"Normalized Choi purity and process fidelity to identity"];
  ExportPlot[delta,p,model,info,dir,"02_choi_entropy",
    {ExtractSeries[delta,res,"ChoiEntropy"]},{"S(chi)=-Tr[chi ln(chi)]"},"Normalized Choi entropy (nats)"];
  ExportPlot[delta,p,model,info,dir,"03_all_16_choi_eigenvalues",
    Table[({#["Time"],#["Eigenvalues"][[k]]} & /@ res),{k,16}],
    Table["lambda_"<>ToString[k],{k,16}],"All 16 eigenvalues of chi, sorted descending at each time"];
  ExportPlot[delta,p,model,info,dir,"04_process_mutual_information",
    {ExtractSeries[delta,res,"ProcessMutualInfo"]},
    {"S(chi_A)+S(chi_B)-S(chi_AB)"},"Process mutual information: (Ain,Aout) | (Bin,Bout), nats"];
  ExportPlot[delta,p,model,info,dir,"05_nonproduct_channel",
    {ExtractSeries[delta,res,"ProductDistance"]},
    {"||chi_(Ain,Aout,Bin,Bout) - chi_A tensor chi_B||_F"},"Distance from the product of Choi marginals (Frobenius)"];
  ExportPlot[delta,p,model,info,dir,"06_pauli_response_A_to_B",
    (ExtractSeries[delta,res,#] & /@ {"TransferX","TransferY","TransferZ"}),
    {"R_Bx,Ax","R_By,Ay","R_Bz,Az"},
    "R_Bk,Ak=Tr[(I tensor sigma_k) E(sigma_k tensor I)]/4"];
  ExportPlot[delta,p,model,info,dir,"07_channel_validation",
    (ExtractSeries[delta,res,#] & /@ {"TraceError","TPError","HermiticityError","MinEigenvalue"}),
    {"|Tr(chi)-1|","||Tr_out(chi)-I/4||_F","||chi-chi^dagger||_F","lambda_min(chi)"},"Complete positivity and trace-preservation checks"];
  ExportPlot[delta,p,model,info,dir,"08_unitality",
    {ExtractSeries[delta,res,"UnitalError"]},{"||E(I_4)-I_4||_F"},"Unitality deviation (nonzero is allowed for exchange)"];
  If[TrueQ[p["CheckDivisibility"]],
    valid=Select[ExtractSeries[delta,res,"IntermediateMinEigenvalue"],NumericQ[Last[#]] &];
    If[valid=!={},
      (* Points only: do not join across unresolved time intervals. *)
      SaveFile[delta,FileNameJoin[{dir,"09_intermediate_CP_test.png"}],Column[{
        Style["Adjacent-time CP test; omitted points are unresolved, not zero",16,Bold],
        Style[PlotCaption[delta,p,model,info],11],
        ListPlot[valid,Frame->True,Axes->False,PlotRange->All,
          FrameLabel->{"Later time t","min eigenvalue of normalized Choi[V(t,s)]"},
          PlotLabel->"s = previous sampled time; negative beyond numerical error witnesses non-CP divisibility",
          ImageSize->1100,LabelStyle->Directive[Black,14]]},Alignment->Center],"PNG"],
      Print["No well-conditioned adjacent intermediate maps; see stored Missing values."]]
  ];
  If[model=="Z",
    pairs=Subsets[Range[4],{2}]; basis={"00","01","10","11"};
    Do[ExportPlot[delta,p,model,info,dir,"10_dephasing_factors_"<>part,
      (MatrixSeries[delta,res,"DephasingFactors",#[[1]],#[[2]],Switch[part,"Re",Re,"Im",Im,"Abs",Abs]] & /@ pairs),
      (part<>" F["<>basis[[#[[1]]]]<>","<>basis[[#[[2]]]]<>"]" & /@ pairs),
      part<>" F_ij(t), with rho_ij(t)=F_ij(t) rho_ij(0)"],{part,{"Re","Im","Abs"}}];
    ExportPlot[delta,p,model,info,dir,"11_collective_coherence_difference",
      {({#["Time"],Abs[#["DephasingFactors"][[1,4]]]-Abs[#["DephasingFactors"][[2,3]]]} & /@ res)},
      {"|F_00,11|-|F_01,10|"},"Parallel minus antiparallel coherence magnitude (no weak-coupling interpretation assumed)"]
  ];
];
RunChannel[delta_?NumericQ,model_,p_Association,chainFunction_] := Catch[Module[
  {cs,cf,engine,times,res,chi,s,f,info,dir,checks,elapsed,metadata},
  ValidateParameters[delta,p]; times=TimeGrid[delta,p];
  Print["Channel: ",model,"; Delta=",delta,"; L=",p["L"],"; t=0..",p["tMax"],"; points=",Length[times]];
  elapsed=First[AbsoluteTiming[
    cs=chainFunction[delta,p]; cf=StateFactor[delta,cs,2^p["L"],p["ValidationTolerance"],"Chain"];
    engine=BuildEngine[delta,model,p,cf["Factor"]];
    res=Table[
      If[model=="Z",f=DephasingFactors[delta,t,engine];chi=ChoiFromFactors[delta,f],
        chi=ChoiMatrix[delta,t,engine]];
      s=ChoiToSuper[delta,chi];
      Join[ChannelObservables[delta,t,chi,s,p],If[model=="Z",<|"DephasingFactors"->f|>,<||>]],{t,times}];
    res=IntermediateChecks[delta,res,p];
  ]];
  checks=<|"MaxTraceError"->Max[Lookup[res,"TraceError"]],"MaxTPError"->Max[Lookup[res,"TPError"]],
    "MaxHermiticityError"->Max[Lookup[res,"HermiticityError"]],
    "MinEigenvalue"->Min[Lookup[res,"MinEigenvalue"]],
    "IdentityAtZeroError"->Norm[First[res]["Superoperator"]-IdentityMatrix[16],"Frobenius"]|>;
  Require[delta,Max[checks["MaxTraceError"],checks["MaxTPError"],checks["MaxHermiticityError"],
    checks["IdentityAtZeroError"]]<=p["ValidationTolerance"] && checks["MinEigenvalue"]>=-p["ValidationTolerance"],
    "Channel validation failed: "<>ToString[checks,InputForm]];
  info=<|"ChainSummary"->("factor columns="<>ToString[Dimensions[cf["Factor"]][[2]]]),
    "ChainHash"->StateHash[delta,cf["Factor"]],
    "ProbeSummary"->"Full map for arbitrary rho_AB(0), independent of the fixed chain input",
    "ProbeHashLine"->"Choi convention: input tensor output; |Phi4>=Sum_i |ii>/2; Tr chi=1"|>;
  dir=CreateRunDirectory[delta,p,model,"channel"];
  metadata=<|"Delta"->delta,"Model"->model,"Parameters"->p,"TimeGrid"->times,
    "Checks"->checks,"EvolutionSeconds"->elapsed,"StateInfo"->info,
    "ChainInputNormOrTrace"->cf["InputTraceOrNorm"],"ChainCorrectionNorm"->cf["CorrectionNorm"],
    "ChainStateDefinition"->DownValues[chainFunction],
    "ChoiConvention"->"input tensor output, normalized to 1", "Vectorization"->"row-major Flatten"|>;
  SaveFile[delta,FileNameJoin[{dir,"initial_states.wxf"}],
    <|"Delta"->delta,"ChainFactor"->cf["Factor"],"ChainInput"->cs|>,"WXF"];
  SaveFile[delta,FileNameJoin[{dir,"parameters.txt"}],ToString[metadata,InputForm],"Text"];
  If[TrueQ[p["SaveData"]],SaveFile[delta,FileNameJoin[{dir,"channel_data.wxf"}],
    <|"Metadata"->metadata,"Results"->res|>,"WXF"]];
  ChannelPlots[delta,p,model,info,dir,res];
  Print["Completed in ",Round[elapsed,.01]," s (evolution); PNGs: ",dir];
  <|"Delta"->delta,"Model"->model,"Parameters"->p,"Directory"->dir,"Checks"->checks,"Results"->res|>
],"RunFailure"];

(* ::Chapter:: *)
(*I - SIGMA-Z / PURE DEPHASING*)

(* ::Section:: *)
(*A - FREE INPUT PARAMETERS: edit this entire block BEFORE execution*)

(* L counts only chain spins. Total spins = L+2.
   deltaListZ: any real anisotropies, e.g. {0.5,1.,1.5}, no integer encoding.
   Couplings, fields, chain size and time interval below belong ONLY to this section.
   tMax is the FULL physical evolution time; dt is the time sampling step.
   Both sections default to the SAME chain input and Hamiltonian parameters.
   Change ChainStateZ below for a product state, ground state, or arbitrary input.
   Ground states are of the ISOLATED chain HC[Delta], recalculated for every Delta.
   State labels describe your preparation; update them when changing the functions.
   initial_states.wxf always records the actual arrays, even when SaveData=False.
   Setting SaveData=False suppresses only the trajectory WXF, not PNG exports.
*)
RunZ = True;
deltaListZ = {0.5, 1., 1.5};
parametersZ = <|
  "L" -> 6,
  "Jxy" -> 1.,             (* Jz is ALWAYS Jxy*Delta *)
  "w" -> 0.,              (* uniform chain field: w Sum Z/2 *)
  "ed" -> 0.,             (* defect field ed Z_d/2; zero switches it off *)
  "DefectSite" -> 3,       (* integer 1..L, even if ed=0 *)
  "omegaA" -> 0., "omegaB" -> 0., (* HQ=omegaQ sigmaQz/2 *)
  "gA" -> 0.1, "gB" -> 0.1,      (* convention stated above and on every PNG *)
  "tMax" -> 50.,          (* FULL EVOLUTION TIME: t=0 through t=50 *)
  "dt" -> 0.1,            (* 501 points here; final tMax always included *)
  "ChainStateLabel" -> "Neel |0101...>, site 1 = |0>",
  "GroundTolerance" -> 10.^-10,
  "ValidationTolerance" -> 10.^-8,
  "InverseTolerance" -> 10.^-8,
  "CheckDivisibility" -> True, (* adjacent-time channel test *)
  "SaveData" -> True,
  "OutputRoot" -> FileNameJoin[{BaseDirectory[0.], "xxz_boundary_results"}]
|>;

(* ::Section:: *)
(*B - FREE INITIAL STATES: functions of Delta; edit BEFORE execution*)

Clear[ChainStateZ];
ChainStateZ[delta_?NumericQ, p_Association] := NeelKet[delta, p["L"]];

(* REPLACE the active definition above by ONE of these examples:

   ChainStateZ[delta_?NumericQ,p_Association] := GroundKet[delta,p];
   (* Also set parametersZ["ChainStateLabel"]="Isolated-chain ground state";
      selection within degenerate ground spaces is documented at GroundKet. *)

   ChainStateZ[delta_?NumericQ,p_Association] :=
     ProductKet[delta, ConstantArray[{1.,0.},p["L"]]];  (* all |0> *)

   ChainStateZ[delta_?NumericQ,p_Association] :=
     ProductKet[delta, Table[{Cos[Pi/6],Exp[I Pi j/7] Sin[Pi/6]},{j,p["L"]}]];
   (* Example product state with local polar angle Pi/3 and phase Pi*j/7. *)

   ChainStateZ[delta_?NumericQ,p_Association] := myChainVector[delta];
   (* Supply exactly 2^L amplitudes, or a 2^L by 2^L density matrix.
      Define myChainVector[delta_] yourself; vectors are normalized automatically.
      To use a fixed state across a sweep, return the same vector for every delta.
      Random input: use BlockRandom[SeedRandom[12345]; ...] for reproducibility
      and use the SAME generator/seed in both files when comparing their outputs. *)
*)

(* No probe input is required to COMPUTE a channel: all 16 input operators
   are represented. An arbitrary probe state can be applied afterwards below. *)

(* ::Section:: *)
(*C - EXECUTION: all free parameters and the FULL TIME are specified above*)

Clear[channelResultsZ,failedZ];
channelResultsZ = {};
If[TrueQ[RunZ],
  Print["Starting Z: L=",parametersZ["L"],"; total spins=",parametersZ["L"]+2,
    "; full time=",parametersZ["tMax"],"; dt=",parametersZ["dt"],"; Delta list=",deltaListZ];
  AbsoluteTiming[
    channelResultsZ=Table[RunChannel[delta,"Z",parametersZ,ChainStateZ],
      {delta,deltaListZ}];
  ] // Print;
  failedZ=Select[channelResultsZ,FailureQ];
  If[failedZ=!={},Print["FAILED runs: ",failedZ]];
];

(* ::Section:: *)
(*D - ACCESS RESULTS: functions of Delta and t; use a stored grid time*)

ChoiZ[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,channelResultsZ,"Choi"];
SuperoperatorZ[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,channelResultsZ,"Superoperator"];
ApplyChannelZ[delta_?NumericQ,t_?NumericQ,rho0_] := Module[{s},
  s=SuperoperatorZ[delta,t]; If[MissingQ[s],s,ApplySuper[delta,s,rho0]]
];
(* Examples after execution:
   ChoiZ[1.,10.] // MatrixForm
   ApplyChannelZ[1.,10.,KetDensity[1.,{1.,1.,1.,1.}/2]] // MatrixForm
   ApplyChannelZ[1.,10.,KetDensity[1.,{0.,0.,1.,0.}]] // MatrixForm
   Compare with RhoAB in the other file only at identical chain input, model,
   Delta, couplings, fields, L and time. The channel does not depend on rho0.
   Exact quantities and parameters appear on EVERY PNG.
   Rerun parameters + states + execution after editing inputs. *)

(* ::Chapter:: *)
(*II - EXCHANGE-COUPLED PROBES*)

(* ::Section:: *)
(*A - FREE INPUT PARAMETERS: edit this entire block BEFORE execution*)

(* L counts only chain spins. Total spins = L+2.
   deltaListExchange: any real anisotropies, e.g. {0.5,1.,1.5}, no integer encoding.
   Couplings, fields, chain size and time interval below belong ONLY to this section.
   tMax is the FULL physical evolution time; dt is the time sampling step.
   Both sections default to the SAME chain input and Hamiltonian parameters.
   Change ChainStateExchange below for a product state, ground state, or arbitrary input.
   Ground states are of the ISOLATED chain HC[Delta], recalculated for every Delta.
   State labels describe your preparation; update them when changing the functions.
   initial_states.wxf always records the actual arrays, even when SaveData=False.
   Setting SaveData=False suppresses only the trajectory WXF, not PNG exports.
*)
RunExchange = True;
deltaListExchange = {0.5, 1., 1.5};
parametersExchange = <|
  "L" -> 6,
  "Jxy" -> 1.,             (* Jz is ALWAYS Jxy*Delta *)
  "w" -> 0.,              (* uniform chain field: w Sum Z/2 *)
  "ed" -> 0.,             (* defect field ed Z_d/2; zero switches it off *)
  "DefectSite" -> 3,       (* integer 1..L, even if ed=0 *)
  "omegaA" -> 0., "omegaB" -> 0., (* HQ=omegaQ sigmaQz/2 *)
  "gA" -> 0.1, "gB" -> 0.1,      (* convention stated above and on every PNG *)
  "tMax" -> 50.,          (* FULL EVOLUTION TIME: t=0 through t=50 *)
  "dt" -> 0.1,            (* 501 points here; final tMax always included *)
  "ChainStateLabel" -> "Neel |0101...>, site 1 = |0>",
  "GroundTolerance" -> 10.^-10,
  "ValidationTolerance" -> 10.^-8,
  "InverseTolerance" -> 10.^-8,
  "CheckDivisibility" -> True, (* adjacent-time channel test *)
  "SaveData" -> True,
  "OutputRoot" -> FileNameJoin[{BaseDirectory[0.], "xxz_boundary_results"}]
|>;

(* ::Section:: *)
(*B - FREE INITIAL STATES: functions of Delta; edit BEFORE execution*)

Clear[ChainStateExchange];
ChainStateExchange[delta_?NumericQ, p_Association] := NeelKet[delta, p["L"]];

(* REPLACE the active definition above by ONE of these examples:

   ChainStateExchange[delta_?NumericQ,p_Association] := GroundKet[delta,p];
   (* Also set parametersExchange["ChainStateLabel"]="Isolated-chain ground state";
      selection within degenerate ground spaces is documented at GroundKet. *)

   ChainStateExchange[delta_?NumericQ,p_Association] :=
     ProductKet[delta, ConstantArray[{1.,0.},p["L"]]];  (* all |0> *)

   ChainStateExchange[delta_?NumericQ,p_Association] :=
     ProductKet[delta, Table[{Cos[Pi/6],Exp[I Pi j/7] Sin[Pi/6]},{j,p["L"]}]];
   (* Example product state with local polar angle Pi/3 and phase Pi*j/7. *)

   ChainStateExchange[delta_?NumericQ,p_Association] := myChainVector[delta];
   (* Supply exactly 2^L amplitudes, or a 2^L by 2^L density matrix.
      Define myChainVector[delta_] yourself; vectors are normalized automatically.
      To use a fixed state across a sweep, return the same vector for every delta.
      Random input: use BlockRandom[SeedRandom[12345]; ...] for reproducibility
      and use the SAME generator/seed in both files when comparing their outputs. *)
*)

(* No probe input is required to COMPUTE a channel: all 16 input operators
   are represented. An arbitrary probe state can be applied afterwards below. *)

(* ::Section:: *)
(*C - EXECUTION: all free parameters and the FULL TIME are specified above*)

Clear[channelResultsExchange,failedExchange];
channelResultsExchange = {};
If[TrueQ[RunExchange],
  Print["Starting Exchange: L=",parametersExchange["L"],"; total spins=",parametersExchange["L"]+2,
    "; full time=",parametersExchange["tMax"],"; dt=",parametersExchange["dt"],"; Delta list=",deltaListExchange];
  AbsoluteTiming[
    channelResultsExchange=Table[RunChannel[delta,"Exchange",parametersExchange,ChainStateExchange],
      {delta,deltaListExchange}];
  ] // Print;
  failedExchange=Select[channelResultsExchange,FailureQ];
  If[failedExchange=!={},Print["FAILED runs: ",failedExchange]];
];

(* ::Section:: *)
(*D - ACCESS RESULTS: functions of Delta and t; use a stored grid time*)

ChoiExchange[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,channelResultsExchange,"Choi"];
SuperoperatorExchange[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,channelResultsExchange,"Superoperator"];
ApplyChannelExchange[delta_?NumericQ,t_?NumericQ,rho0_] := Module[{s},
  s=SuperoperatorExchange[delta,t]; If[MissingQ[s],s,ApplySuper[delta,s,rho0]]
];
(* Examples after execution:
   ChoiExchange[1.,10.] // MatrixForm
   ApplyChannelExchange[1.,10.,KetDensity[1.,{1.,1.,1.,1.}/2]] // MatrixForm
   ApplyChannelExchange[1.,10.,KetDensity[1.,{0.,0.,1.,0.}]] // MatrixForm
   Compare with RhoAB in the other file only at identical chain input, model,
   Delta, couplings, fields, L and time. The channel does not depend on rho0.
   Exact quantities and parameters appear on EVERY PNG.
   Rerun parameters + states + execution after editing inputs. *)

End[];
EndPackage[];
