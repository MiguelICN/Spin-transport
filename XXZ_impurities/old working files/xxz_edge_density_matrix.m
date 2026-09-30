(* ::Package:: *)

(* ::Title:: *)
(*XXZ boundary probes - reduced density matrix*)


(* Edit this file, then evaluate Get["/your/path/xxz_edge_density_matrix.m"].
   To rerun individual definition/input/execution cells after loading, first use
   Begin["EdgeDensity`Private`"]; and finish with End[].
   The public result functions below remain accessible after Get finishes. *)

BeginPackage["EdgeDensity`"];
ClearAll["EdgeDensity`*","EdgeDensity`Private`*"];
RhoABZ::usage = "RhoABZ: see editable sections and examples in this file.";
RhoAZ::usage = "RhoAZ: see editable sections and examples in this file.";
RhoBZ::usage = "RhoBZ: see editable sections and examples in this file.";
RhoABExchange::usage = "RhoABExchange: see editable sections and examples in this file.";
RhoAExchange::usage = "RhoAExchange: see editable sections and examples in this file.";
RhoBExchange::usage = "RhoBExchange: see editable sections and examples in this file.";
densityResultsZ::usage = "densityResultsZ: see editable sections and examples in this file.";
densityResultsExchange::usage = "densityResultsExchange: see editable sections and examples in this file.";
Begin["`Private`"];


SetDirectory[NotebookDirectory[]];


(* ::Section::Closed:: *)
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



(* ::Section::Closed:: *)
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
  root = Transpose[rows] . (Sqrt[Clip[Re[v], {0., Infinity}]] Conjugate[rows]);
  flip = KroneckerProduct[Pauli[delta, 2], Pauli[delta, 2]];
  ev = Reverse[Sort[Sqrt[Clip[Re[Eigenvalues[HermitianPart[delta,
    root . flip . Conjugate[r] . flip . root]]], {0., Infinity}]]]];
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
  rho = fac . ConjugateTranspose[fac];
  <|"Factor"->fac, "InputTraceOrNorm"->Tr[x], "CorrectionNorm"->Norm[rho-x,"Frobenius"]|>
];



(* ::Section::Closed:: *)
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
  v=Transpose[rows] . Conjugate[rows[[All,k]]]; v=v/Norm[v];
  result=ConstantArray[0.+0.I,2^p["L"]]; result[[chosen["Indices"]]]=v; result
];



(* ::Section::Closed:: *)
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
      "LocalEnergies"->(signs . {p["omegaA"],p["omegaB"]}/2)|>,
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
      c=Conjugate[es[[2]]] . initial[[ix,All]];
      Sow[<|"Indices"->ix,"E"->es[[1]],"V"->Transpose[es[[2]]],"C"->c|>]
    ],{ix,groups}]][[2]];
  <|"Dimension"->Length[h],"Columns"->Dimensions[initial][[2]],
    "Sectors"->If[data==={},{},First[data]]|>
];
EvolveBlocks[delta_, t_?NumericQ, data_] := Module[{out},
  out=ConstantArray[0.+0.I,{data["Dimension"],data["Columns"]}];
  Do[out[[b["Indices"],All]]=b["V"] . (Exp[-I b["E"] t] b["C"]),{b,data["Sectors"]}]; out
];
DephasingFactors[delta_?NumericQ,t_?NumericQ,engine_] := Module[{vectors},
  Require[delta, engine["Delta"]==delta && engine["Model"]=="Z", "Engine/Delta/model mismatch"];
  vectors=Table[Exp[-I engine["LocalEnergies"][[s]] t]
    Flatten[EvolveBlocks[delta,t,engine["Blocks"][[s]]]],{s,4}];
  vectors . ConjugateTranspose[vectors]
];
(* Return four 4 x (dC rankC) matrices W_i, one for each AB input |i>.
   E(|i><j|)=W_i.W_j^dagger. Chain columns include sqrt(probability). *)
OutputBlocks[delta_?NumericQ,t_?NumericQ,engine_] := Module[{out,n=engine["L"],r=engine["RankC"]},
  Require[delta, engine["Delta"]==delta, "Engine/Delta mismatch"];
  out=EvolveBlocks[delta,t,engine["Blocks"]];
  Table[ArrayReshape[Transpose[ArrayReshape[out[[All,(s-1) r+1;;s r]],
    {2,2^n,2,r}],{1,3,2,4}],{4,2^n r}],{s,4}]
];



(* ::Section::Closed:: *)
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



(* ::Section::Closed:: *)
(*5 - Reduced density matrix and observables*)


RhoAB[delta_?NumericQ,t_?NumericQ,engine_,rho0_] := Module[{f,w},
  If[engine["Model"]=="Z", f=DephasingFactors[delta,t,engine]; rho0 f,
    w=OutputBlocks[delta,t,engine];
    Sum[rho0[[i,j]] w[[i]] . ConjugateTranspose[w[[j]]],{i,4},{j,4}]]
];
DensityObservables[delta_,t_,r_,rho0_] := Module[{a,b,sa,sb,sab,za,zb,xx,zz,ma,mb},
  a=TraceOutB[delta,r]; b=TraceOutA[delta,r];
  {sa,sb,sab}=VNEntropy[delta,#] & /@ {a,b,r};
  za=KroneckerProduct[Pauli[delta,3],IdentityMatrix[2]];
  zb=KroneckerProduct[IdentityMatrix[2],Pauli[delta,3]];
  ma=Re[Tr[za . r]]; mb=Re[Tr[zb . r]];
  zz=Re[Tr[za . zb . r]]-ma mb;
  xx=Re[Tr[KroneckerProduct[Pauli[delta,1],Pauli[delta,1]] . r]]-
    Re[Tr[Pauli[delta,1] . a]] Re[Tr[Pauli[delta,1] . b]];
  <|"Delta"->delta,"Time"->t,"RhoAB"->r,"RhoA"->a,"RhoB"->b,
    "Purity"->Re[Tr[r . r]],"EntropyA"->sa,"EntropyB"->sb,"EntropyAB"->sab,
    "MutualInfo"->sa+sb-sab,"Concurrence"->Concurrence[delta,r],
    "Negativity"->Negativity[delta,r],"InitialOverlap"->Re[Tr[rho0 . r]],
    "MagnetizationA"->ma,"MagnetizationB"->mb,"ConnectedZZ"->zz,"ConnectedXX"->xx,
    "TraceError"->Abs[Tr[r]-1.],"HermiticityError"->Norm[r-ConjugateTranspose[r],"Frobenius"],
    "MinEigenvalue"->Min[Re[Eigenvalues[HermitianPart[delta,r]]]]|>
];
DensityPlots[delta_,p_,model_,info_,dir_,res_] := Module[{pairs,basis,labels},
  basis={"00","01","10","11"}; pairs=Subsets[Range[4],{2}];
  ExportPlot[delta,p,model,info,dir,"01_populations",
    Table[MatrixSeries[delta,res,"RhoAB",i,i,Re],{i,4}],
    ("rho["<>#<> ","<>#<>"]" & /@ basis),"Probe populations rho_AB[s,s](t)"];
  Do[
    labels=(part<>" rho["<>basis[[#[[1]]]]<>","<>basis[[#[[2]]]]<>"]" & /@ pairs);
    ExportPlot[delta,p,model,info,dir,"02_coherences_"<>part,
      (MatrixSeries[delta,res,"RhoAB",#[[1]],#[[2]],Switch[part,"Re",Re,"Im",Im,"Abs",Abs]] & /@ pairs),
      labels,part<>" of all six independent probe coherences"],{part,{"Re","Im","Abs"}}];
  ExportPlot[delta,p,model,info,dir,"03_entropies_and_mutual_information",
    (ExtractSeries[delta,res,#] & /@ {"EntropyA","EntropyB","EntropyAB","MutualInfo"}),
    {"S(rho_A)","S(rho_B)","S(rho_AB)","I(A:B)=S_A+S_B-S_AB"},"Von Neumann entropies and mutual information (nats)"];
  ExportPlot[delta,p,model,info,dir,"04_entanglement",
    (ExtractSeries[delta,res,#] & /@ {"Concurrence","Negativity"}),
    {"Wootters concurrence C(rho_AB)","N=(||rho_AB^(T_B)||_1-1)/2"},"Two-probe entanglement"];
  ExportPlot[delta,p,model,info,dir,"05_purity_and_initial_overlap",
    (ExtractSeries[delta,res,#] & /@ {"Purity","InitialOverlap"}),
    {"Tr[rho_AB(t)^2]","Tr[rho_AB(0) rho_AB(t)]"},
    "Purity and initial overlap (survival probability only for pure rho_AB(0))"];
  ExportPlot[delta,p,model,info,dir,"06_local_magnetizations",
    (ExtractSeries[delta,res,#] & /@ {"MagnetizationA","MagnetizationB"}),
    {"<sigma_A^z>","<sigma_B^z>"},"Local probe magnetizations (Pauli convention)"];
  ExportPlot[delta,p,model,info,dir,"07_connected_correlations",
    (ExtractSeries[delta,res,#] & /@ {"ConnectedZZ","ConnectedXX"}),
    {"<ZA ZB>-<ZA><ZB>","<XA XB>-<XA><XB>"},"Equal-time connected PROBE correlations"];
  ExportPlot[delta,p,model,info,dir,"08_density_validation",
    (ExtractSeries[delta,res,#] & /@ {"TraceError","HermiticityError","MinEigenvalue"}),
    {"|Tr(rho_AB)-1|","||rho_AB-rho_AB^dagger||_F","lambda_min(rho_AB)"},"Reduced-state numerical checks"];
];
RunDensity[delta_?NumericQ,model_,p_Association,chainFunction_,probeFunction_] := Catch[Module[
  {cs,ps,cf,pf,rho0,engine,times,res,r,info,dir,checks,elapsed,metadata},
  ValidateParameters[delta,p]; times=TimeGrid[delta,p];
  Print["Density: ",model,"; Delta=",delta,"; L=",p["L"],"; t=0..",p["tMax"],"; points=",Length[times]];
  elapsed=First[AbsoluteTiming[
    cs=chainFunction[delta,p]; ps=probeFunction[delta,p];
    cf=StateFactor[delta,cs,2^p["L"],p["ValidationTolerance"],"Chain"];
    pf=StateFactor[delta,ps,4,p["ValidationTolerance"],"AB"];
    rho0=pf["Factor"] . ConjugateTranspose[pf["Factor"]];
    engine=BuildEngine[delta,model,p,cf["Factor"]];
    res=Table[r=RhoAB[delta,t,engine,rho0]; DensityObservables[delta,t,r,rho0],{t,times}];
  ]];
  checks=<|"MaxTraceError"->Max[Lookup[res,"TraceError"]],
    "MaxHermiticityError"->Max[Lookup[res,"HermiticityError"]],
    "MinEigenvalue"->Min[Lookup[res,"MinEigenvalue"]],
    "InitialStateError"->Norm[First[res]["RhoAB"]-rho0,"Frobenius"],
    "DephasingPopulationDrift"->If[model=="Z",Max[Abs[Flatten[
      (Diagonal[#["RhoAB"]]-Diagonal[rho0] & /@ res)]]],0.]|>;
  Require[delta,Max[checks["MaxTraceError"],checks["MaxHermiticityError"],
    checks["InitialStateError"],checks["DephasingPopulationDrift"]]<=p["ValidationTolerance"] &&
    checks["MinEigenvalue"]>=-p["ValidationTolerance"],"Reduced-state validation failed: "<>ToString[checks,InputForm]];
  info=<|"ChainSummary"->("factor columns="<>ToString[Dimensions[cf["Factor"]][[2]]]),
    "ChainHash"->StateHash[delta,cf["Factor"]],
    "ProbeSummary"->Row[{"AB input: ",p["ProbeStateLabel"],"; rho_AB(0)=",MatrixForm[rho0]}],
    "ProbeHashLine"->("AB SHA256: "<>StateHash[delta,rho0])|>;
  dir=CreateRunDirectory[delta,p,model,"density"];
  metadata=<|"Delta"->delta,"Model"->model,"Parameters"->p,"TimeGrid"->times,
    "Checks"->checks,"EvolutionSeconds"->elapsed,"StateInfo"->info,
    "ChainInputNormOrTrace"->cf["InputTraceOrNorm"],"ChainCorrectionNorm"->cf["CorrectionNorm"],
    "ProbeInputNormOrTrace"->pf["InputTraceOrNorm"],"ProbeCorrectionNorm"->pf["CorrectionNorm"],
    "ChainStateDefinition"->DownValues[chainFunction],"ProbeStateDefinition"->DownValues[probeFunction]|>;
  SaveFile[delta,FileNameJoin[{dir,"initial_states.wxf"}],<|"Delta"->delta,
    "ChainFactor"->cf["Factor"],"RhoAB0"->rho0,"ChainInput"->cs,"ProbeInput"->ps|>,"WXF"];
  SaveFile[delta,FileNameJoin[{dir,"parameters.txt"}],ToString[metadata,InputForm],"Text"];
  If[TrueQ[p["SaveData"]],SaveFile[delta,FileNameJoin[{dir,"density_data.wxf"}],
    <|"Metadata"->metadata,"Results"->res|>,"WXF"]];
  DensityPlots[delta,p,model,info,dir,res];
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
deltaListZ = {0.01, 0.5, 1., 1.5,10.};
parametersZ = <|
  "L" -> 6,
  "Jxy" -> 1.,             (* Jz is ALWAYS Jxy*Delta *)
  "w" -> 0.,              (* uniform chain field: w Sum Z/2 *)
  "ed" -> 0.,             (* defect field ed Z_d/2; zero switches it off *)
  "DefectSite" -> 3,       (* integer 1..L, even if ed=0 *)
  "omegaA" -> 1., "omegaB" -> 1., (* HQ=omegaQ sigmaQz/2 *)
  "gA" -> 0.1, "gB" -> 0.1,      (* convention stated above and on every PNG *)
  "tMax" -> 200.,          (* FULL EVOLUTION TIME: t=0 through t=50 *)
  "dt" -> 0.1,            (* 501 points here; final tMax always included *)
  "ChainStateLabel" -> "Neel |0101...>, site 1 = |0>",
  "ProbeStateLabel" -> "|+>_A tensor |+>_B",
  "GroundTolerance" -> 10.^-10,
  "ValidationTolerance" -> 10.^-8,
  "InverseTolerance" -> 10.^-8,
  "CheckDivisibility" -> False, (* unused in density file *)
  "SaveData" -> True,
  "OutputRoot" -> FileNameJoin[{BaseDirectory[0.], "xxz_boundary_results"}]
|>;



(* ::Section::Closed:: *)
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

Clear[ProbeStateZ];
ProbeStateZ[delta_?NumericQ,p_Association] := {1.,1.,1.,1.}/2;
(* Free two-probe input: any four-vector or 4x4 density matrix.
   Basis is {00,01,10,11}. Examples:
   ProbeStateZ[delta_?NumericQ,p_Association] := {0.,0.,1.,0.}; (* |10> *)
   ProbeStateZ[delta_?NumericQ,p_Association] := {0.,1.,1.,0.}/Sqrt[2.]; (* Bell *)
   ProbeStateZ[delta_?NumericQ,p_Association] := IdentityMatrix[4]/4.;
   For exchange, |10> is a useful local magnetization-transfer preparation.
   Keep |++> to compare dephasing and exchange with the same input state.
   Update ProbeStateLabel above. This state is independent of the chain initially. *)



(* ::Section:: *)
(*C - EXECUTION: all free parameters and the FULL TIME are specified above*)


Clear[densityResultsZ,failedZ];
densityResultsZ = {};
If[TrueQ[RunZ],
  Print["Starting Z: L=",parametersZ["L"],"; total spins=",parametersZ["L"]+2,
    "; full time=",parametersZ["tMax"],"; dt=",parametersZ["dt"],"; Delta list=",deltaListZ];
  AbsoluteTiming[
    densityResultsZ=Table[RunDensity[delta,"Z",parametersZ,ChainStateZ,ProbeStateZ],
      {delta,deltaListZ}];
  ] // Print;
  failedZ=Select[densityResultsZ,FailureQ];
  If[failedZ=!={},Print["FAILED runs: ",failedZ]];
];



(* ::Section::Closed:: *)
(*D - ACCESS RESULTS: functions of Delta and t; use a stored grid time*)


RhoABZ[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsZ,"RhoAB"];
RhoAZ[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsZ,"RhoA"];
RhoBZ[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsZ,"RhoB"];
(* Examples after execution: RhoABZ[1.,10.] // MatrixForm
   densityResultsZ[[1]]["Results"] is the first Delta's full list of observables.
   Exact definitions and all parameters appear on EVERY exported PNG.
   Rerun parameters + states + execution to change L, fields, couplings or Delta.
   No stale diagonalization is reused after a parameter edit. *)



(* ::Chapter:: *)
(*II - EXCHANGE-COUPLED PROBES*)


(* ::Section::Closed:: *)
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
deltaListExchange = {0.01,.5, 1., 1.5,10.};
parametersExchange = <|
  "L" -> 6,
  "Jxy" -> 1.,             (* Jz is ALWAYS Jxy*Delta *)
  "w" -> 0.,              (* uniform chain field: w Sum Z/2 *)
  "ed" -> 0.,             (* defect field ed Z_d/2; zero switches it off *)
  "DefectSite" -> 3,       (* integer 1..L, even if ed=0 *)
  "omegaA" -> 1., "omegaB" -> 1., (* HQ=omegaQ sigmaQz/2 *)
  "gA" -> 0.1, "gB" -> 0.1,      (* convention stated above and on every PNG *)
  "tMax" -> 200.,          (* FULL EVOLUTION TIME: t=0 through t=50 *)
  "dt" -> 0.1,            (* 501 points here; final tMax always included *)
  "ChainStateLabel" -> "Neel |0101...>, site 1 = |0>",
  "ProbeStateLabel" -> "|+>_A tensor |+>_B",
  "GroundTolerance" -> 10.^-10,
  "ValidationTolerance" -> 10.^-8,
  "InverseTolerance" -> 10.^-8,
  "CheckDivisibility" -> False, (* unused in density file *)
  "SaveData" -> True,
  "OutputRoot" -> FileNameJoin[{BaseDirectory[0.], "xxz_boundary_results"}]
|>;



(* ::Section::Closed:: *)
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

Clear[ProbeStateExchange];
ProbeStateExchange[delta_?NumericQ,p_Association] := {1.,1.,1.,1.}/2;
(* Free two-probe input: any four-vector or 4x4 density matrix.
   Basis is {00,01,10,11}. Examples:
   ProbeStateExchange[delta_?NumericQ,p_Association] := {0.,0.,1.,0.}; (* |10> *)
   ProbeStateExchange[delta_?NumericQ,p_Association] := {0.,1.,1.,0.}/Sqrt[2.]; (* Bell *)
   ProbeStateExchange[delta_?NumericQ,p_Association] := IdentityMatrix[4]/4.;
   For exchange, |10> is a useful local magnetization-transfer preparation.
   Keep |++> to compare dephasing and exchange with the same input state.
   Update ProbeStateLabel above. This state is independent of the chain initially. *)



(* ::Section:: *)
(*C - EXECUTION: all free parameters and the FULL TIME are specified above*)


Clear[densityResultsExchange,failedExchange];
densityResultsExchange = {};
If[TrueQ[RunExchange],
  Print["Starting Exchange: L=",parametersExchange["L"],"; total spins=",parametersExchange["L"]+2,
    "; full time=",parametersExchange["tMax"],"; dt=",parametersExchange["dt"],"; Delta list=",deltaListExchange];
  AbsoluteTiming[
    densityResultsExchange=Table[RunDensity[delta,"Exchange",parametersExchange,ChainStateExchange,ProbeStateExchange],
      {delta,deltaListExchange}];
  ] // Print;
  failedExchange=Select[densityResultsExchange,FailureQ];
  If[failedExchange=!={},Print["FAILED runs: ",failedExchange]];
];



(* ::Section::Closed:: *)
(*D - ACCESS RESULTS: functions of Delta and t; use a stored grid time*)


RhoABExchange[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsExchange,"RhoAB"];
RhoAExchange[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsExchange,"RhoA"];
RhoBExchange[delta_?NumericQ,t_?NumericQ] := DataAt[delta,t,densityResultsExchange,"RhoB"];
(* Examples after execution: RhoABExchange[1.,10.] // MatrixForm
   densityResultsExchange[[1]]["Results"] is the first Delta's full list of observables.
   Exact definitions and all parameters appear on EVERY exported PNG.
   Rerun parameters + states + execution to change L, fields, couplings or Delta.
   No stale diagonalization is reused after a parameter edit. *)

End[];
EndPackage[];
