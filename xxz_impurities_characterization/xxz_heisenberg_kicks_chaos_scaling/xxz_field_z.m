(* ::Title:: *)
(* Impurity-field kicks along z *)

(* ::Section:: *)
(* 1. MODEL, DRIVING, AND ALL USER PARAMETERS *)

(* ::Subsection:: *)
(* Definitions: physical model and the fixed driving class *)
ClearAll["Global`*"];

(* Order: {A,1,...,Lchain,B}; spins = sigma/2; hbar=1.
   OPEN XXZ chain with TWO additional endpoint spins:
   H0 = J Sum_chain (Sjx S(j+1)x + Sjy S(j+1)y + Delta Sjz S(j+1)z)
        + OmegaA sAz + OmegaB sBz
        + gXYA (sAx S1x+sAy S1y) + gZA sAz S1z
        + gXYB (SLx sBx+SLy sBy) + gZB SLz sBz.

   THIS FILE: LOCAL FIELD pulses on the two extra impurities, along z.
   VA = s_A^z
   VB = s_B^z
   H(t) = H0 + Sum_n delta(t-n T) (alphaA VA+alphaB VB),
   alphaA=KickScale KickA; alphaB=KickScale KickB; T=Period.
   Pulses are simultaneous and have the same period.
   A pulse area is a rotation angle for a FIELD pulse. For a BOND pulse,
   it is an integrated coupling strength. Static and pulse couplings are separate.

   Parameter-independent sector guarantee: total magnetization Mz.
   Reflection/spin flips/SU(2)/local charges can depend on parameter choices.
   They are tested BEFORE Floquet block exponentiation and diagonalization.
   x/y-field files have no unconditional unitary-sector guarantee with generic
   nonzero probe z-fields. Zero Omegas add the corresponding global flip.

   Ideal delta kicks and all eigenvalues are retained at model level.
   Arithmetic has finite precision. No time steps, Trotter/Magnus expansion,
   bath approximation, unfolding, or spectral truncation. *)

driveTag="field_z";
driveLabel="Impurity-field kicks along z";
drivingLocation="field";
fieldAxis="z";              (* fixed by this file *)
bondWeights={0,0,0};           (* fixed {XX,YY,ZZ} pulse components *)
expectedSymmetries={"Mz"};

Lchain=6;                         (* chain ONLY: 6+2 = 8 total spins *)
requestedKernels=8;
workingPrecision=MachinePrecision; (* optional 30; then use exact inputs *)
saveEigenvectors=False;          (* True: every sector eigenvector *)
keepFloquetBlocks=False;         (* True: retain sector matrices *)
exportPlots=True;
numberOfBins=24;
minLevelsForRPlot=20;            (* plot only; every block remains in the table *)

p=<|"J"->1,"Delta"->1/2,"OmegaA"->3/5,"OmegaB"->9/10,
   "gXYA"->1/5,"gXYB"->1/5,"gZA"->1/10,"gZB"->1/10,
   "Period"->1,"KickA"->Pi/2,"KickB"->Pi/3,"KickScale"->1|>;
scanParameter="KickScale";       (* any key of p *)
scanValues=Range[0,2,1/10];       (* overrides the baseline scanParameter *)
extraSweepRules[x_]:=<||>;       (* optional simultaneous changes *)

(* EXAMPLES: edit above; rerun the entire file.
   Delta: scanParameter="Delta"; scanValues={0,1/2,1,2,5}.
   Period: scanParameter="Period"; scanValues=Range[1/5,2,1/10].
   Only pulse A: set KickB->0 in p; scan "KickA".
   Reflection: match OmegaB=OmegaA, gXYB=gXYA, gZB=gZA, KickB=KickA.
   x/y-field flip symmetry: set OmegaA=OmegaB=0.
   Pure ZZ STATIC attachments: set gXYA=gXYB=0.
   Sweep BOTH XY attachments: scanParameter="gXYA";
      extraSweepRules[x_]:=<|"gXYB"->x|>.
   Homogeneous STATIC enlarged XXZ:
      gXYA=gXYB=J, gZA=gZB=J Delta, OmegaA=OmegaB=0.
   When scanning Delta, keep matching through
      extraSweepRules[x_]:=<|"gZA"->p["J"] x,"gZB"->p["J"] x|>.
   Axes/weights identify THIS file. Use another file for another driving class.
   Convention: g(s+S-+s-S+)=2g(sxSx+sySy), so gXY=2g in that convention. *)

(* ::Subsection:: *)
(* Computations: validate parameters *)
tolerance=If[workingPrecision===MachinePrecision,10^-10,
   10^-Floor[workingPrecision/2]];
gapTolerance=10 tolerance;
nspin=Lchain+2;dim=2^nspin;
If[!IntegerQ[Lchain] || Lchain<2 || !MemberQ[Keys[p],scanParameter] ||
   !ListQ[scanValues] || Length[scanValues]==0,
   Print["Use Lchain>=2, a key of p, and nonempty scanValues."];Abort[]];
pSweep=(Join[p,<|scanParameter->#|>,extraSweepRules[#]]&)/@scanValues;
If[!TrueQ[And@@((And@@(NumericQ/@Values[#]) &&
      And@@(TrueQ[Im[#]==0]&/@Values[#]) && #["Period"]>0)&/@pSweep)],
   Print["Parameters must be real numeric; Period must be positive."];Abort[]];
If[!IntegerQ[numberOfBins] || numberOfBins<1 ||
   !IntegerQ[minLevelsForRPlot] || minLevelsForRPlot<3,Abort[]];
Print[driveLabel,"; chain=",Lchain,"; total spins=",nspin,
   "; dimension=",dim,"; expected symmetries=",expectedSymmetries];

(* ::Section:: *)
(* 2. EIGHT KERNELS AND EXACT SPARSE OPERATORS *)

(* ::Subsection:: *)
(* Definitions *)
parallelMap[f_,items_]:=ParallelMap[f,items,
   Method->"CoarsestGrained",DistributedContexts->None];
identity[d_]:=IdentityMatrix[d,SparseArray];
spin=<|"x"->PauliMatrix[1]/2,"y"->PauliMatrix[2]/2,"z"->PauliMatrix[3]/2|>;
single[a_,i_]:=KroneckerProduct[identity[2^(i-1)],SparseArray[a],identity[2^(nspin-i)]];
adjacent[a_,i_]:=KroneckerProduct[identity[2^(i-1)],SparseArray[a],identity[2^(nspin-i-1)]];
maxEntry[m_SparseArray]:=Max[0,Abs[m["NonzeroValues"]]];
maxEntry[m_]:=Max[Abs[Flatten[m]]];
relativeCommutator[s_,a_]:=maxEntry[s.a-a.s]/
   (Max[1,maxEntry[s]] Max[1,maxEntry[a]]);
H0[q_]:=q["J"](bulkXY+q["Delta"] bulkZZ)+q["OmegaA"] zA+q["OmegaB"] zB+
   q["gXYA"] edgeXYA+q["gXYB"] edgeXYB+q["gZA"] edgeZZA+q["gZB"] edgeZZB;
pulseMatrices[q_]:=Module[{a,b},
   a=MatrixExp[N[-I q["KickScale"] q["KickA"] localVA,workingPrecision]];
   b=MatrixExp[N[-I q["KickScale"] q["KickB"] localVB,workingPrecision]];
   If[drivingLocation=="field",{single[a,1],single[b,nspin]},
      {adjacent[a,1],adjacent[b,nspin-1]}]];
kickMatrix[q_]:=Module[{a,b},{a,b}=pulseMatrices[q];b.a];

(* ::Subsection:: *)
(* Computations *)
(* BEGIN PARALLEL STARTUP *)
If[$KernelCount<requestedKernels,LaunchKernels[requestedKernels-$KernelCount]];
If[$KernelCount!=requestedKernels,
   Print["Requested ",requestedKernels," kernels; running ",$KernelCount,
      ". Check the license/configuration or close surplus kernels."];Abort[]];
ParallelEvaluate[$HistoryLength=0;];
Print["Parallel kernel IDs: ",ParallelEvaluate[$KernelID]];
(* END PARALLEL STARTUP *)
localXY=KroneckerProduct[spin["x"],spin["x"]]+KroneckerProduct[spin["y"],spin["y"]];
localZZ=KroneckerProduct[spin["z"],spin["z"]];
localBond=bondWeights[[1]] KroneckerProduct[spin["x"],spin["x"]]+
   bondWeights[[2]] KroneckerProduct[spin["y"],spin["y"]]+bondWeights[[3]] localZZ;
bulkXY=Total[Table[adjacent[localXY,i],{i,2,nspin-2}]];
bulkZZ=Total[Table[adjacent[localZZ,i],{i,2,nspin-2}]];
zA=single[spin["z"],1];zB=single[spin["z"],nspin];
edgeXYA=adjacent[localXY,1];edgeXYB=adjacent[localXY,nspin-1];
edgeZZA=adjacent[localZZ,1];edgeZZB=adjacent[localZZ,nspin-1];
localVA=If[drivingLocation=="field",spin[fieldAxis],localBond];localVB=localVA;

(* ::Section:: *)
(* 3. EXACT FLOQUET FACTORIZATION *)

(* ::Subsection:: *)
(* Definitions *)
(* Exact pulse jump: Exp[-I(alphaA VA+alphaB VB)].
   Disjoint supports imply [VA,VB]=0 exactly, hence
   UF=KB.KA.Exp[-I H0 T].
   The 4x4/16x16 pulse-support check is the same identity at any chain length.
   Reflection may preserve KB.KA without preserving KA and KB individually.
   We therefore restrict the COMBINED kick, not the individual kicks. *)
verifyFactorization[q_]:=Module[{d,va,vb,a,b,ka,kb,combined},
   d=Length[localVA];va=KroneckerProduct[localVA,IdentityMatrix[d]];
   vb=KroneckerProduct[IdentityMatrix[d],localVB];
   a=q["KickScale"] q["KickA"];b=q["KickScale"] q["KickB"];
   ka=MatrixExp[N[-I a va,workingPrecision]];
   kb=MatrixExp[N[-I b vb,workingPrecision]];
   combined=MatrixExp[N[-I(a va+b vb),workingPrecision]];
   <|"ExactCommutation"->(va.vb==vb.va),
     "FactorizationError"->maxEntry[combined-kb.ka]|>];
(* ::Subsection:: *)
(* Computations *)
factorizationChecks=verifyFactorization/@pSweep;
factorizationCheck=<|"ExactCommutation"->And@@Lookup[factorizationChecks,"ExactCommutation"],
   "FactorizationError"->Max[Lookup[factorizationChecks,"FactorizationError"]],
   "ParameterPointsChecked"->Length[pSweep]|>;
Print["Exact pulse factorization: ",factorizationCheck];
If[!TrueQ[factorizationCheck["ExactCommutation"]] ||
   factorizationCheck["FactorizationError"]>tolerance,Abort[]];

(* ::Section:: *)
(* 4. DECLARE AND VERIFY SECTORS BEFORE FLOQUET DIAGONALIZATION *)

(* ::Subsection:: *)
(* Definitions *)
(* Mz=total magnetization; BulkMz excludes A,B; sAz/sBz=probe magnetizations;
   Pz=product sigma_z; S2=total spin squared; Fx/Fy=products sigma_x/sigma_y;
   Reflection reverses the full chain including interchange of A and B.
   We require [S,H0]=[S,KB.KA]=0. These common symmetries permit exact
   sector-by-sector EXPONENTIATION as well as diagonalization.
   Floquet-only symmetries that do not preserve H0 are not exhaustively
   searched. This finite candidate list is not exhaustive. *)
candidateNames={"Mz","BulkMz","sAz","sBz","Pz","S2","Fx","Fy","Reflection"};
globalFlip[a_]:=KroneckerProduct@@ConstantArray[SparseArray[PauliMatrix[a]],nspin];
reflectionMatrix[]:=SparseArray[Table[
   {FromDigits[Reverse[IntegerDigits[b,2,nspin]],2]+1,b+1}->1,
   {b,0,dim-1}],{dim,dim}];
symmetryCheck[h_,k_]:=Association@Table[name-><|
   "H0Error"->relativeCommutator[symmetryOperators[name],h],
   "KickError"->relativeCommutator[symmetryOperators[name],k]|>,{name,candidateNames}];
refineSector[block_,name_,s_]:=Module[
   {b=block["Basis"],sb,values,vectors,labels,groups},
   sb=ConjugateTranspose[b].(s.b);
   If[maxEntry[s.b-b.sb]/Max[1,maxEntry[s]]>tolerance,Return[{block}]];
   {values,vectors}=Eigensystem[N[sb,workingPrecision]];
   labels=Round[4 Re[values]]/4;
   (* Labels are known quarter-integer symmetry charges. Floquet phases
      and eigenvalues are NEVER rounded. *)
   groups=GatherBy[Range[Length[values]],labels[[#]]&];
   If[Length[groups]==1,Return[{block}]];
   Table[<|"Labels"->Join[block["Labels"],<|name->labels[[First[g]]]|>],
      "Basis"->b.Transpose[Orthogonalize[vectors[[g]]]]|>,{g,groups}]];
buildSectorBasis[names_]:=Module[{blocks,b},
   blocks={<|"Labels"-><||>,"Basis"->IdentityMatrix[dim]|>};
   Do[blocks=Flatten[refineSector[#,name,symmetryOperators[name]]&/@blocks,1],
      {name,names}];
   b=Join[Sequence@@Lookup[blocks,"Basis"],2];
   If[Dimensions[b]!={dim,dim} ||
      maxEntry[ConjugateTranspose[b].b-IdentityMatrix[dim]]>10 tolerance,
      Return[Failure["InvalidBasis",<||>]]];
   blocks];
sectorBasisCache[names_List]:=sectorBasisCache[names]=buildSectorBasis[names];

(* ::Subsection:: *)
(* Computations: fixed candidate matrices *)
totalX=Total[Table[single[spin["x"],i],{i,nspin}]];
totalY=Total[Table[single[spin["y"],i],{i,nspin}]];
totalZ=Total[Table[single[spin["z"],i],{i,nspin}]];
symmetryOperators=<|"Mz"->totalZ,"BulkMz"->totalZ-zA-zB,"sAz"->zA,"sBz"->zB,
   "Pz"->globalFlip[3],"S2"->totalX.totalX+totalY.totalY+totalZ.totalZ,
   "Fx"->globalFlip[1],"Fy"->globalFlip[2],"Reflection"->reflectionMatrix[]|>;
(* Each point will verify its sectors BEFORE exponentiating/diagonalizing UF.
   Only the tiny local pulse exponentials have been computed so far. *)

(* ::Section:: *)
(* 5. COMPUTATIONS INSIDE VERIFIED SECTORS, IN PARALLEL *)

(* ::Subsection:: *)
(* Definitions *)
freeBlockCache=<||>;
staticKeys={"J","Delta","OmegaA","OmegaB","gXYA","gXYB","gZA","gZB","Period"};
freeSignature[q_]:=Lookup[q,staticKeys];
sectorFree[q_,names_,block_,hblock_]:=Module[{key},
   If[!reuseSectorFree,Return[MatrixExp[N[-I q["Period"] hblock,workingPrecision]]]];
   key=ToString[{names,Normal[block["Labels"]]},InputForm];
   If[!KeyExistsQ[freeBlockCache,key],
      AssociateTo[freeBlockCache,key->MatrixExp[N[-I q["Period"] hblock,workingPrecision]]]];
   freeBlockCache[key]];
circularRatios[eigenvalues_]:=Module[{phases,gaps,next,d,ratios},
   d=Length[eigenvalues];phases=Sort[Mod[-Arg[eigenvalues]+Pi,2 Pi]-Pi];
   If[d<3,Return[<|"Phases"->phases,"Gaps"->{},"Ratios"->{},
      "MeanR"->Missing["TooFewLevels"],"Status"->"Fewer than 3 levels"|>]];
   gaps=Differences[Append[phases,First[phases]+2 Pi]];
   If[Min[gaps]<=gapTolerance,
      Return[<|"Phases"->phases,"Gaps"->gaps,"Ratios"->{},
         "MeanR"->Missing["DegenerateOrNearDegenerate"],
         "Status"->"Check degeneracy/additional symmetry/precision"|>]];
   next=RotateLeft[gaps];
   ratios=MapThread[Min[#1,#2]/Max[#1,#2]&,{gaps,next}];
   <|"Phases"->phases,"Gaps"->gaps,"Ratios"->ratios,"MeanR"->Mean[ratios],
     "Status"->"OK"|>];
solveBlock[q_,names_,block_,h_,k_]:=Module[
   {b=block["Basis"],hb,kb,leakH,leakK,u0,uf,values,vectors,answer},
   If[block["Labels"]==<||>,
      hb=Normal[N[h,workingPrecision]];kb=k;leakH=0;leakK=0,
      hb=ConjugateTranspose[b].(h.b);kb=ConjugateTranspose[b].(k.b);
      leakH=maxEntry[h.b-b.hb]/Max[1,maxEntry[h]];leakK=maxEntry[k.b-b.kb]];
   If[Max[leakH,leakK]>10 tolerance,
      Return[Failure["NonInvariantSector",<|"Sector"->block["Labels"]|>]]];
   (* Both H0 and combined K preserve this sector: no projection approximation. *)
   u0=sectorFree[q,names,block,hb];uf=Normal[kb.u0];
   If[saveEigenvectors,{values,vectors}=Eigensystem[uf],values=Eigenvalues[uf]];
   answer=Join[<|"Sector"->block["Labels"],"Dimension"->Length[values],
      "Eigenvalues"->values,"H0Leakage"->leakH,"KickLeakage"->leakK,
      "UnitarityError"->maxEntry[ConjugateTranspose[uf].uf-IdentityMatrix[Length[uf]]],
      "ModulusError"->Max[Abs[Abs[values]-1]]|>,circularRatios[values]];
   answer=Join[answer,<|"Quasienergies"->answer["Phases"]/q["Period"]|>];
   If[saveEigenvectors,answer=Join[answer,<|"Eigenvectors"->vectors|>]];
   If[keepFloquetBlocks,answer=Join[answer,<|"FloquetBlock"->uf|>]];
   answer];
solvePoint[q_]:=Module[{h,k,checks,names,blocks,results},
   h=H0[q];k=kickMatrix[q];checks=symmetryCheck[h,k];
   names=Select[candidateNames,Max[Values[checks[#]]]<=tolerance&];
   If[!SubsetQ[names,expectedSymmetries],
      Return[Failure["ExpectedSymmetryFailed",<|"Checks"->checks|>]]];
   blocks=sectorBasisCache[names];If[FailureQ[blocks],Return[blocks]];
   results=solveBlock[q,names,#,h,k]&/@blocks;
   If[AnyTrue[results,FailureQ],Return[First[Select[results,FailureQ]]]];
   <|"Parameters"->q,"ParameterValue"->q[scanParameter],"KernelID"->$KernelID,
     "SymmetryChecks"->checks,"ConservedCandidates"->names,"Blocks"->results|>];

(* ::Subsection:: *)
(* Computations *)
reuseSectorFree=TrueQ[And@@(freeSignature[#]==freeSignature[pSweep[[1]]]&/@pSweep)];
Print["Free sector propagators reused: ",reuseSectorFree,
   "; additional symmetries checked at EVERY parameter point."];
DistributeDefinitions[solvePoint,solveBlock,sectorFree,circularRatios,symmetryCheck,
   sectorBasisCache,buildSectorBasis,refineSector,kickMatrix,pulseMatrices,H0,
   single,adjacent,identity,maxEntry,relativeCommutator,workingPrecision,
   tolerance,gapTolerance,reuseSectorFree,saveEigenvectors,keepFloquetBlocks,
   scanParameter,expectedSymmetries,candidateNames,symmetryOperators,nspin,dim,
   drivingLocation,localVA,localVB,bulkXY,bulkZZ,zA,zB,edgeXYA,edgeXYB,edgeZZA,edgeZZB];
(* BEGIN CACHE RESET *)
ParallelEvaluate[freeBlockCache=<||>;Clear[sectorBasisCache];
   sectorBasisCache[names_List]:=sectorBasisCache[names]=buildSectorBasis[names];];
(* END CACHE RESET *)
{sweepSeconds,pointResults}=AbsoluteTiming[parallelMap[solvePoint,pSweep]];
If[AnyTrue[pointResults,FailureQ],Print[Select[pointResults,FailureQ]];Abort[]];
allBlocks=Flatten[Lookup[pointResults,"Blocks"],1];
If[Max[Lookup[allBlocks,"UnitarityError"]]>10 tolerance ||
   Max[Lookup[allBlocks,"ModulusError"]]>10 tolerance ||
   !TrueQ[And@@(Total[Lookup[#["Blocks"],"Dimension"]]==dim&/@pointResults)],
   Print["Sector spectral checks failed."];Abort[]];
Print["Sweep seconds: ",sweepSeconds,"; kernels used: ",
   Sort[DeleteDuplicates[Lookup[pointResults,"KernelID"]]]];
sectorSummary=Flatten[Table[
   {i,pointResults[[i]]["ParameterValue"],pointResults[[i]]["ConservedCandidates"],
      #["Sector"],#["Dimension"],#["MeanR"],#["Status"]}&/@pointResults[[i]]["Blocks"],
   {i,Length[pointResults]}],1];
Print[Grid[Prepend[sectorSummary,
   {"Point","Scan value","Conserved candidates","Sector","Dimension","<r>","Status"}],
   Frame->All,Alignment->Left]];
Print["Candidate list is not exhaustive. Sectors are NEVER pooled for <r>."];

(* ::Section:: *)
(* 6. PLOTS: EVERY FIGURE DECLARES ALL PARAMETERS USED *)

(* ::Subsection:: *)
(* Definitions *)
formatValue[x_]:=If[NumericQ[x],ToString[x,InputForm],x];
parameterValues[qs_,key_]:=Module[{v=DeleteDuplicates[Lookup[qs,key]]},
   If[Length[v]==1,formatValue[First[v]],
      Grid[Partition[PadRight[formatValue/@v,6 Ceiling[Length[v]/6],""],6],
         Alignment->Left,Spacings->{0.4,0.2}]]];
parameterHeader[qs_,note_]:=Column[{
   Style[driveLabel<>" | open XXZ | S=sigma/2",13,Bold],
   Row[{"Lchain=",Lchain,", total spins=",nspin,", dimension=",dim,
      ", pulses at A and B simultaneously"}],
   Grid[Partition[(Row[{#," = ",parameterValues[qs,#]}]&/@Keys[p]),3],
      Alignment->Left,Spacings->{1,0.4}],
   Style[note,11]},Alignment->Left];
allPhases[i_]:=Sort[Flatten[Lookup[pointResults[[i]]["Blocks"],"Phases"]]];
allQuasienergies[i_]:=Sort[Flatten[Lookup[pointResults[[i]]["Blocks"],"Quasienergies"]]];
densityPoints[i_]:=Module[{density},
   density=BinCounts[N[allQuasienergies[i]],{-energyBound,energyBound,binWidth}]/
      (dim binWidth);
   Flatten[Table[{{-energyBound+(b-1) binWidth,density[[b]]},
      {-energyBound+b binWidth,density[[b]]}},{b,numberOfBins}],1]];
makeDiscreteDOS[i_]:=Total[DiracDelta/@(energy-allQuasienergies[i])]/dim;
makeDiscreteSectorDOS[i_,j_]:=With[{b=pointResults[[i]]["Blocks"][[j]]},
   Total[DiracDelta/@(energy-b["Quasienergies"])]/b["Dimension"]];
sectorDOSPlot[i_,j_]:=Module[{b,density,points},
   b=pointResults[[i]]["Blocks"][[j]];
   density=BinCounts[N[b["Quasienergies"]],{-energyBound,energyBound,binWidth}]/
      (b["Dimension"] binWidth);
   points=Flatten[Table[{{-energyBound+(k-1) binWidth,density[[k]]},
      {-energyBound+k binWidth,density[[k]]}},{k,numberOfBins}],1];
   ListLinePlot[points,Frame->True,Axes->False,
      FrameLabel->{"quasienergy epsilon","normalized sector DOS"},
      PlotLabel->parameterHeader[{pSweep[[i]]},
         "Sector "<>ToString[Normal[b["Sector"]],InputForm]<>"; dimension="<>
         ToString[b["Dimension"]]<>"; bins="<>ToString[numberOfBins]],
      PlotRange->{{-energyBound,energyBound},{0,All}},ImageSize->1050]];
(* Optional after the sweep: sectorDOSPlot[pointIndex,blockIndex].
   makeDiscreteDOS / makeDiscreteSectorDOS return the exact discrete DOS
   from the computed eigenvalues. Bins affect the display only. *)
sectorKey[labels_]:=ToString[Normal[labels],InputForm];

(* ::Subsection:: *)
(* Computations *)
spectrumPoints=Flatten[Table[({pointResults[[i]]["ParameterValue"],#}&/@allPhases[i]),
   {i,Length[pointResults]}],1];
spectrumPlot=ListPlot[N[spectrumPoints],Frame->True,Axes->False,
   FrameLabel->{scanParameter,"quasienergy phase = epsilon T"},
   PlotLabel->parameterHeader[pSweep,"Full spectrum: all verified sectors"],
   PlotRange->{All,{-Pi,Pi}},PlotStyle->Directive[Blue,PointSize[0.004]],ImageSize->1050];
Print[spectrumPlot];

dosIndices=DeleteDuplicates[Round[Subdivide[1,Length[pointResults],3]]];
energyBound=N[Max[Pi/#["Period"]&/@pSweep]](1+10^-12);
binWidth=2 energyBound/numberOfBins;
dosPlots=Table[ListLinePlot[densityPoints[i],Frame->True,Axes->False,
   FrameLabel->{"quasienergy epsilon","normalized DOS"},
   PlotLabel->parameterHeader[{pSweep[[i]]},
      "DOS of all sectors; point "<>ToString[i]<>"; bins="<>ToString[numberOfBins]],
   PlotRange->{{-energyBound,energyBound},{0,All}},PlotStyle->Directive[Blue,Thick],
   ImageSize->1050],{i,dosIndices}];
Scan[Print,dosPlots];

rRows=Flatten[Table[
   ({sectorKey[#["Sector"]],N[pointResults[[i]]["ParameterValue"]],N[#["MeanR"]]}&/@
      Select[pointResults[[i]]["Blocks"],
         NumericQ[#["MeanR"]] && #["Dimension"]>=minLevelsForRPlot&]),
   {i,Length[pointResults]}],1];
rGroups=GatherBy[rRows,First];
ratioPlot=If[rGroups=={},
   Graphics[Text["No valid sectors with the requested minimum level count"],
      PlotLabel->parameterHeader[pSweep,"No <r> data; inspect sectorSummary"],ImageSize->1050],
   ListPlot[(#[[All,{2,3}]]&/@rGroups),Frame->True,Axes->False,
      PlotLegends->(First[#][[1]]&/@rGroups),PlotRange->{All,{0,1}},
      FrameLabel->{scanParameter,"mean circular gap ratio <r>"},
      PlotLabel->parameterHeader[pSweep,
         "Each series is ONE sector; at least "<>ToString[minLevelsForRPlot]<>
         " levels. Reference: Poisson~0.3863, COE~0.53."],ImageSize->1050]];
Print[ratioPlot];

(* ::Section:: *)
(* 7. EXPORTS AND ACCESS TO VERIFIED BLOCKS *)

(* ::Subsection:: *)
(* Definitions *)
baseDirectory[]:=Module[{base},
   base=If[StringQ[$InputFileName] && $InputFileName!="",DirectoryName[$InputFileName],
      Quiet[Check[NotebookDirectory[],Directory[]]]];
   If[StringQ[base],base,Directory[]]];
(* pointResults[[i]]["Blocks"][[j]] contains one sector result.
   keepFloquetBlocks=True retains its complete ["FloquetBlock"].
   Full-space basis:
      sectorBasisCache[pointResults[[i]]["ConservedCandidates"]][[j]]["Basis"].
   Eigenvectors (if requested) are rows in that sector basis.
   No tiny/degenerate sector is deleted. Its status remains in the table.
   minLevelsForRPlot controls the plot only. *)

(* ::Subsection:: *)
(* Computations: keep previous runs and export parameter-labeled PNGs *)
If[exportPlots,
   rootDirectory=FileNameJoin[{baseDirectory[],"xxz_driving_results",driveTag,
      "Lchain"<>ToString[Lchain],scanParameter}];
   If[!DirectoryQ[rootDirectory],
      CreateDirectory[rootDirectory,CreateIntermediateDirectories->True]];
   runNumber=1;
   While[DirectoryQ[FileNameJoin[{rootDirectory,"run"<>IntegerString[runNumber,10,3]}]],
      runNumber++];
   outputDirectory=FileNameJoin[{rootDirectory,"run"<>IntegerString[runNumber,10,3]}];
   CreateDirectory[outputDirectory];
   Export[FileNameJoin[{outputDirectory,"spectrum.png"}],spectrumPlot,ImageResolution->160];
   Export[FileNameJoin[{outputDirectory,"mean_r_by_sector.png"}],ratioPlot,ImageResolution->160];
   Do[Export[FileNameJoin[{outputDirectory,"DOS_point"<>ToString[dosIndices[[j]]]<>".png"}],
      dosPlots[[j]],ImageResolution->160],{j,Length[dosPlots]}];
   Export[FileNameJoin[{outputDirectory,"parameters.csv"}],
      Prepend[(Lookup[#,Keys[p]]&/@pSweep),Keys[p]]];
   Export[FileNameJoin[{outputDirectory,"sector_statistics.csv"}],
      Prepend[({#[[1]],#[[2]],ToString[#[[3]],InputForm],
         ToString[#[[4]],InputForm],#[[5]],#[[6]],#[[7]]}&/@sectorSummary),
         {"Point","ScanValue","ConservedCandidates","Sector","Dimension","MeanR","Status"}]];
   Export[FileNameJoin[{outputDirectory,"results.wxf"}],<|
      "Driving"->driveLabel,"Location"->drivingLocation,"FieldAxis"->fieldAxis,
      "BondWeights"->bondWeights,"Lchain"->Lchain,"Parameters"->pSweep,
      "ScanParameter"->scanParameter,"WorkingPrecision"->workingPrecision,
      "FactorizationCheck"->factorizationCheck,"PointResults"->pointResults|>];
   Print["Saved: ",outputDirectory]];

