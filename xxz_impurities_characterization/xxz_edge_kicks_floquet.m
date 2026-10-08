(* ::Title:: *)
(* XXZ endpoint kicks: exact Floquet spectra on eight parallel kernels *)

(* ::Text:: *)
(* Open this .m as a notebook, or run Get["full/path/xxz_edge_kicks_floquet.m"].
   Evaluate top to bottom after changing Section 1. Each section separates
   DEFINITIONS from COMPUTATIONS. Timing is printed for the expensive stages.
   Spin order: {A, 1, ..., Lchain, B}; S = sigma/2; hbar = 1.
   Both protocols have simultaneous instantaneous kicks with the SAME period.
   No time discretization, Magnus/Trotter expansion, or spectral truncation.
   Numerical roundoff remains. Binning changes only the plotted DOS.
   Parallelism is over independent parameter points, not individual eigenvalues. *)

(* ::Section:: *)
(* 1. PARAMETERS: all user choices are here *)

(* ::Subsection:: *)
(* Definitions: edit these values *)

ClearAll["Global`*"];

Lchain = 6;                         (* 6 chain + 2 probes = 8 TOTAL spins.
                                      Lchain=8 means 10 total spins, dim=1024. *)
requestedKernels = 8;
workingPrecision = MachinePrecision; (* fast double precision; use 30 or 40
                                        only if extra precision is needed *)
saveEigenvectors = False;          (* False: ALL eigenvalues, enough for DOS.
                                      True: also compute/store ALL eigenvectors *)
exportPlots = True;
numberOfBins = 24;

(* Physical model:
   H0 = J Sum_chain [Sx.Sx + Sy.Sy + Delta Sz.Sz]
      + OmegaA sAz + OmegaB sBz
      + gXYA (sAx S1x + sAy S1y) + gZA sAz S1z
      + gXYB (SLx sBx + SLy sBy) + gZB SLz sBz.
   The chain itself has no applied magnetic field in this script. *)
p = <|
   "J" -> 1,                       (* bulk XY exchange; bulk ZZ = J Delta *)
   "Delta" -> 1/2,                 (* bulk anisotropy *)
   "OmegaA" -> 3/5, "OmegaB" -> 9/10, (* static probe z-field strengths *)
   "gXYA" -> 1/5, "gXYB" -> 1/5,  (* STATIC endpoint XY exchange strengths *)
   "gZA" -> 1/10, "gZB" -> 1/10,  (* STATIC endpoint ZZ strengths *)
   "Period" -> 1,                  (* T>0; frequency = 2 Pi/T *)
   "KickA" -> Pi/2, "KickB" -> Pi/3, (* pulse areas at the two ends *)
   "KickScale" -> 1                (* multiplies BOTH pulse areas *)
   |>;

(* Field protocol: QA = KickScale KickA sA^(fieldAxisA), and likewise for B.
   Effective rotation angles are KickScale KickA and KickScale KickB. *)
fieldAxisA = "x";                  (* independently choose "x", "y", or "z" *)
fieldAxisB = "x";

(* Bond protocol: QA = KickScale KickA
                  (wx sAx S1x + wy sAy S1y + wz sAz S1z), likewise for B.
   Pulse areas are dimensionless integrated coupling strengths.
   These pulse couplings are additional to the STATIC couplings in p. *)
bondWeights = {1, 1, 0};           (* {wx,wy,wz}; XY exchange pulses *)
(* {1,0,0}: XX pulse; {0,1,0}: YY pulse; {0,0,1}: ZZ pulse;
   {1,1,1}: Heisenberg pulse. One choice applies to both endpoint bonds. *)

(* Change ONE parameter in p while keeping the others fixed.
   The sweep overrides that parameter's baseline value in p.
   Reuse of exp(-I H0 T) is automatic for pulse-only sweeps. *)
scanParameter = "KickScale";
scanValues = Range[0, 2, 1/10];
(* Optional: couple several changes to the SAME sweep variable x.
   Leave <||> for an ordinary one-parameter sweep.
   Example below for equal attachment strengths; see recipes. *)
extraSweepRules[x_] := <||>;

(* READY-TO-USE RECIPES. Replace settings above; rerun the entire file.
   A. Bulk anisotropy: scanParameter="Delta"; scanValues={0,1/2,1,2,5}.
   B. Period: scanParameter="Period"; scanValues=Range[1/5,2,1/10].
   C. Only kick A: set p["KickB"]=0; sweep "KickA".
      Both protocols retain the same H0.
   D. Transverse axes: fieldAxisA="x"; fieldAxisB="y".
      Field axes affect ONLY Field. bondWeights affects ONLY Bond.
   E. Pure ZZ attachment: p["gXYA"]=0; p["gXYB"]=0.
      With z-field kicks, [K,H0]=0; x/y-field kicks generally do not commute.
   F. Sweep both XY attachments together:
      scanParameter="gXYA"; scanValues=Range[0,1,1/20];
      extraSweepRules[x_]:=<|"gXYB"->x|>.
   G. Homogeneous static enlarged XXZ chain:
      p["gXYA"]=p["J"]; p["gXYB"]=p["J"];
      p["gZA"]=p["J"] p["Delta"]; p["gZB"]=p["J"] p["Delta"];
      p["OmegaA"]=0; p["OmegaB"]=0.
      If scanning Delta, keep endpoint ZZ terms matched with:
      extraSweepRules[x_]:=<|"gZA"->p["J"] x,"gZB"->p["J"] x|>.
   H. High precision: workingPrecision=30; keep exact inputs (1/5, not 0.2).
   Normalization: g (s+ S- + s- S+) = 2 g (sx Sx + sy Sy).
   To reproduce that convention in the attached notebook, set gXY=2 g.
   Lchain, axes, weights, bins, and numerical controls are NOT keys of p:
   edit them directly and rerun. They cannot be scanned with scanParameter. *)

(* ::Subsection:: *)
(* Computations: validate choices and construct sweep inputs *)

checkTolerance = If[workingPrecision === MachinePrecision, 10^-10,
   10^-Floor[workingPrecision/2]];
If[!IntegerQ[Lchain] || Lchain < 2,
   Print["Lchain must be at least 2: the kicked endpoint bonds must be disjoint."];
   Abort[]];
If[!MemberQ[Keys[p], scanParameter] || !ListQ[scanValues] ||
   Length[scanValues] == 0,
   Print["Choose scanParameter from Keys[p], and a nonempty scanValues list."];
   Abort[]];
If[!MemberQ[{"x","y","z"},fieldAxisA] || !MemberQ[{"x","y","z"},fieldAxisB] ||
   Length[bondWeights] != 3 || !IntegerQ[numberOfBins] || numberOfBins < 1,
   Print["Check field axes, three bondWeights, and positive numberOfBins."];
   Abort[]];
pSweep = (Join[p, <|scanParameter -> #|>, extraSweepRules[#]] &) /@ scanValues;
If[!TrueQ[And @@ ((And @@ (NumericQ /@ Values[#]) &&
         And @@ (TrueQ[Im[#] == 0] & /@ Values[#]) && #["Period"] > 0) & /@
      Join[{p}, pSweep])] ||
   !TrueQ[And @@ (NumericQ[#] && Im[#] == 0 & /@ bondWeights)],
   Print["All physical parameters must be real numeric, with Period > 0."];
   Abort[]];
nspin = Lchain + 2; dim = 2^nspin;
Print["Chain spins = ", Lchain, "; total spins = ", nspin,
   "; full Floquet dimension = ", dim, "; precision = ", workingPrecision];

(* ::Section:: *)
(* 2. EIGHT-KERNEL PARALLEL SETUP *)

(* ::Subsection:: *)
(* Definitions *)
parallelMap[f_, items_] := ParallelMap[f, items,
   Method -> "CoarsestGrained", DistributedContexts -> None];

(* ::Subsection:: *)
(* Computations *)
(* BEGIN PARALLEL STARTUP *)
If[$KernelCount < requestedKernels, LaunchKernels[requestedKernels - $KernelCount]];
If[$KernelCount != requestedKernels,
   Print["Requested ", requestedKernels, " kernels, but ", $KernelCount,
      " are running. Check your parallel-kernel license/configuration.",
      " If extra kernels are running, close the surplus and rerun."];
   Abort[]];
ParallelEvaluate[$HistoryLength = 0;];
Print["Using ", $KernelCount, " parallel kernels. IDs: ", ParallelEvaluate[$KernelID]];
(* END PARALLEL STARTUP *)

(* ::Section:: *)
(* 3. SPARSE STATIC OPERATORS AND SMALL PULSE MATRICES *)

(* ::Subsection:: *)
(* Definitions *)
identity[d_] := IdentityMatrix[d, SparseArray];
localSpin = <|"x" -> PauliMatrix[1]/2, "y" -> PauliMatrix[2]/2,
   "z" -> PauliMatrix[3]/2|>;
single[a_, site_] := KroneckerProduct[
   identity[2^(site-1)], SparseArray[a], identity[2^(nspin-site)]];
adjacent[a_, site_] := KroneckerProduct[
   identity[2^(site-1)], SparseArray[a], identity[2^(nspin-site-1)]];

H0[q_Association] := q["J"] (bulkXY + q["Delta"] bulkZZ) +
   q["OmegaA"] zA + q["OmegaB"] zB +
   q["gXYA"] edgeXYA + q["gXYB"] edgeXYB +
   q["gZA"] edgeZZA + q["gZB"] edgeZZB;
staticKeys = {"J","Delta","OmegaA","OmegaB","gXYA","gXYB","gZA","gZB"};
freeSignature[q_] := Lookup[q, Append[staticKeys, "Period"]];
freeMatrix[q_] := MatrixExp[Normal[N[-I q["Period"] H0[q], workingPrecision]]];

localKicks[q_, kind_] := Module[{v = pulseGenerators[kind]},
   {MatrixExp[N[-I q["KickScale"] q["KickA"] v[[1]], workingPrecision]],
    MatrixExp[N[-I q["KickScale"] q["KickB"] v[[2]], workingPrecision]]}];
kickMatrices[q_, kind_] := Module[{ka, kb},
   {ka,kb} = localKicks[q,kind];
   If[kind == "Field", {single[ka,1], single[kb,nspin]},
      {adjacent[ka,1], adjacent[kb,nspin-1]}]];
floquetFromFree[q_, kind_, u0_] := Module[{ka,kb},
   {ka,kb} = kickMatrices[q,kind];
   Normal[kb.(ka.u0)]];
maxEntry[m_] := Max[Abs[Flatten[Normal[m]]]];

(* ::Subsection:: *)
(* Computations: build once *)
setupSeconds = First[AbsoluteTiming[
   localXY = KroneckerProduct[localSpin["x"],localSpin["x"]] +
      KroneckerProduct[localSpin["y"],localSpin["y"]];
   localZZ = KroneckerProduct[localSpin["z"],localSpin["z"]];
   localBond = bondWeights[[1]] KroneckerProduct[localSpin["x"],localSpin["x"]] +
      bondWeights[[2]] KroneckerProduct[localSpin["y"],localSpin["y"]] +
      bondWeights[[3]] localZZ;
   bulkXY = Total[Table[adjacent[localXY,i],{i,2,nspin-2}]];
   bulkZZ = Total[Table[adjacent[localZZ,i],{i,2,nspin-2}]];
   zA = single[localSpin["z"],1]; zB = single[localSpin["z"],nspin];
   edgeXYA = adjacent[localXY,1]; edgeXYB = adjacent[localXY,nspin-1];
   edgeZZA = adjacent[localZZ,1]; edgeZZB = adjacent[localZZ,nspin-1];
   pulseGenerators = <|"Field" -> {localSpin[fieldAxisA],localSpin[fieldAxisB]},
      "Bond" -> {localBond,localBond}|>;
   protocols = {"Field","Bond"};
   ]];
Print["Static-operator setup seconds: ", setupSeconds];

(* ::Section:: *)
(* 4. EXACT FLOQUET FACTORIZATION *)

(* ::Subsection:: *)
(* Definitions: proof and checks *)
(* The exact delta-pulse jump is Exp[-I(QA+QB)].
   UF = Exp[-I(QA+QB)].Exp[-I H0 T].
   The field pulses act on distinct spins; bond pulses act on disjoint pairs.
   Thus [QA,QB]=0 and UF = KB.KA.Exp[-I H0 T] exactly.
   The checks below use the complete pulse support: 2+2 or 4+4 tensor factors,
   hence 4x4 or 16x16 matrices. Remaining spins are identity spectators.
   This is the SAME operator identity at every chain length, and avoids
   exponentiating a large matrix just to verify a local pulse identity.
   It does NOT assert [K,H0]=0, or UF=Exp[-I(T H0+QA+QB)]. *)

verifyFactorization[kind_] := Module[
   {v = pulseGenerators[kind], d, va, vb, alpha, beta, ka, kb, together},
   d = Length[v[[1]]];
   va = KroneckerProduct[v[[1]],IdentityMatrix[d]];
   vb = KroneckerProduct[IdentityMatrix[d],v[[2]]];
   alpha = p["KickScale"] p["KickA"]; beta = p["KickScale"] p["KickB"];
   ka = MatrixExp[N[-I alpha va,workingPrecision]];
   kb = MatrixExp[N[-I beta vb,workingPrecision]];
   together = MatrixExp[N[-I(alpha va+beta vb),workingPrecision]];
   <|"SupportDimension" -> d^2,
     "ExactCommutation" -> (va.vb == vb.va),
     "KickFactorizationError" -> maxEntry[together-kb.ka]|>];

(* ::Subsection:: *)
(* Computations: distribute definitions, then verify in parallel *)
DistributeDefinitions[verifyFactorization,pulseGenerators,p,workingPrecision,maxEntry];
factorizationChecks = AssociationThread[protocols,
   parallelMap[verifyFactorization,protocols]];
Print["Exact factorization checks: ", factorizationChecks];
If[!TrueQ[And @@ Lookup[Values[factorizationChecks],"ExactCommutation"]] ||
   Max[Lookup[Values[factorizationChecks],"KickFactorizationError"]] > checkTolerance,
   Print["Factorization checks failed."]; Abort[]];

(* ::Section:: *)
(* 5. FULL FLOQUET SPECTRA: PARAMETER POINTS RUN IN PARALLEL *)

(* ::Subsection:: *)
(* Definitions *)
(* Each worker receives one q, computes/reuses Ufree, and solves BOTH protocols.
   Changing static couplings or T requires one Ufree per parameter point,
   reused between protocols. Pulse-only sweeps need Ufree just ONCE.
   Never cache every large Ufree on every kernel.
   All dim eigenvalues are computed even when saveEigenvectors=False. *)
solveProtocol[q_, kind_, u0_] := Module[
   {uf, eigenvalues, eigenvectors, phases, order, answer},
   uf = floquetFromFree[q,kind,u0];
   If[saveEigenvectors,
      {eigenvalues,eigenvectors} = Eigensystem[uf],
      eigenvalues = Eigenvalues[uf]];
   phases = Mod[-Arg[eigenvalues]+Pi,2 Pi]-Pi; order = Ordering[phases];
   answer = <|"ParameterValue" -> q[scanParameter], "Period" -> q["Period"],
      "Eigenvalues" -> eigenvalues[[order]], "Phases" -> phases[[order]],
      "Quasienergies" -> phases[[order]]/q["Period"],
      "EigenvalueModulusError" -> Max[Abs[Abs[eigenvalues]-1]]|>;
   If[saveEigenvectors, answer = Join[answer,
      <|"Eigenvectors" -> eigenvectors[[order]],
        "EigenvectorResidual" -> maxEntry[
           uf.Transpose[eigenvectors] -
           Transpose[eigenvectors].DiagonalMatrix[eigenvalues]]|>]];
   answer];
solvePoint[q_] := Module[{u0},
   u0 = If[reuseFree, sharedFree, freeMatrix[q]];
   <|"KernelID" -> $KernelID, "Parameters" -> q,
     "Field" -> solveProtocol[q,"Field",u0],
     "Bond" -> solveProtocol[q,"Bond",u0]|>];

(* ::Subsection:: *)
(* Computations *)
reuseFree = TrueQ[And @@ (freeSignature[#] == freeSignature[p] & /@ pSweep)];
{freeSeconds, sharedFree} = AbsoluteTiming[If[reuseFree,freeMatrix[p],None]];
Print[If[reuseFree, "One shared free propagator computed in " <> ToString[freeSeconds] <>
   " seconds.", "H0 or T changes: one free propagator per parameter point."]];

DistributeDefinitions[solvePoint,solveProtocol,floquetFromFree,kickMatrices,localKicks,
   single,adjacent,identity,freeMatrix,H0,bulkXY,bulkZZ,zA,zB,edgeXYA,edgeXYB,
   edgeZZA,edgeZZB,pulseGenerators,nspin,dim,workingPrecision,saveEigenvectors,
   scanParameter,reuseFree,sharedFree,maxEntry];
Print["Diagonalizing ", Length[pSweep], " parameter points for both protocols."];
{sweepSeconds, pointResults} = AbsoluteTiming[parallelMap[solvePoint,pSweep]];
spectra = Association@Table[kind -> Lookup[pointResults,kind],{kind,protocols}];
usedKernelIDs = Sort[DeleteDuplicates[Lookup[pointResults,"KernelID"]]];
spectrumChecks = Association@Table[kind -> <|
   "AllLevelsRetained" -> (And @@
      (Length[#["Eigenvalues"]] == dim & /@ spectra[kind])),
   "MaxModulusError" -> Max[Lookup[spectra[kind],"EigenvalueModulusError"]],
   "MaxEigenvectorResidual" -> If[saveEigenvectors,
      Max[Lookup[spectra[kind],"EigenvectorResidual"]],Missing["NotRequested"]]|>,
   {kind,protocols}];
Print["Sweep seconds: ", sweepSeconds, "; kernel IDs used: ", usedKernelIDs];
Print["Spectrum checks: ", spectrumChecks];
If[!TrueQ[And @@ Lookup[Values[spectrumChecks],"AllLevelsRetained"]] ||
   Max[Lookup[Values[spectrumChecks],"MaxModulusError"]] > checkTolerance ||
   (saveEigenvectors &&
      Max[Lookup[Values[spectrumChecks],"MaxEigenvectorResidual"]] > checkTolerance),
   Print["Spectrum checks failed. Increase precision or check inputs."]; Abort[]];

(* ::Section:: *)
(* 6. INSPECT THE BASELINE FULL MATRICES AND CHECK UNITARITY *)

(* ::Subsection:: *)
(* Definitions *)
(* These are full dim x dim matrices at p, not necessarily a sweep point.
   Run MatrixForm[UFField] or MatrixForm[UFBond] to inspect their entries.
   Full unitarity is checked once per protocol, not at every sweep point. *)
checkFullUnitarity[uf_] :=
   maxEntry[ConjugateTranspose[uf].uf-IdentityMatrix[Length[uf]]];

(* ::Subsection:: *)
(* Computations *)
baselineFree = If[reuseFree,sharedFree,freeMatrix[p]];
UFField = floquetFromFree[p,"Field",baselineFree];
UFBond = floquetFromFree[p,"Bond",baselineFree];
Hstatic = H0[p];
DistributeDefinitions[checkFullUnitarity,maxEntry];
unitarityChecks = AssociationThread[protocols,
   parallelMap[checkFullUnitarity,{UFField,UFBond}]];
Print["Baseline full-matrix unitarity errors: ", unitarityChecks];
If[Max[Values[unitarityChecks]] > checkTolerance,
   Print["Unitarity check failed."]; Abort[]];

(* ::Section:: *)
(* 7. QUASIENERGY SPECTRUM AND DOS *)

(* ::Subsection:: *)
(* Definitions *)
(* lambda=Exp[-I epsilon T]; phase=epsilon T in [-Pi,Pi);
   epsilon in [-Pi/T,Pi/T). Full multiplicities are retained.
   Normalized discrete DOS: rho(e)=(1/dim) Sum_a delta(e-epsilon_a).
   makeDiscreteDOS["Field",i] returns that expression on demand.
   Plot binning is display only, with no smoothing or pooling across q. *)
phasePoints[kind_] := Flatten[Table[
   ({scanValues[[i]],#} & /@ spectra[kind][[i]]["Phases"]),
   {i,Length[scanValues]}],1];
makeDiscreteDOS[kind_,i_] :=
   Total[DiracDelta /@ (energy-spectra[kind][[i]]["Quasienergies"])]/dim;
binDensities[kind_,i_] := BinCounts[N[spectra[kind][[i]]["Quasienergies"]],
   {-energyBound,energyBound,binWidth}]/(dim binWidth);
densityPoints[kind_,i_] := Module[{density=binDensities[kind,i]},
   Flatten[Table[{{-energyBound+(b-1) binWidth,density[[b]]},
      {-energyBound+b binWidth,density[[b]]}},{b,numberOfBins}],1]];

(* ::Subsection:: *)
(* Computations: inexpensive plotting stays on the main kernel *)
spectrumPlots = Table[ListPlot[N[phasePoints[kind]],Frame->True,Axes->False,
   FrameLabel->{scanParameter,"quasienergy phase = epsilon T"},
   PlotLabel->kind<>" kicks",PlotRange->{All,{-Pi,Pi}},
   PlotStyle->Directive[If[kind=="Field",Blue,Red],PointSize[0.004]],
   ImageSize->500],{kind,protocols}];
spectrumFigure = GraphicsRow[spectrumPlots,ImageSize->1050];
Print[spectrumFigure];
dosIndices = DeleteDuplicates[Round[Subdivide[1,Length[scanValues],3]]];
energyBound = N[Max[Pi/#["Period"] & /@ pSweep]] (1+10^-12);
(* Tiny outward bin padding avoids losing a roundoff-shifted zone endpoint;
   it does not move or discard any eigenvalue. *)
binWidth = 2 energyBound/numberOfBins;
dosPlots = Table[ListLinePlot[
   {densityPoints["Field",i],densityPoints["Bond",i]},
   PlotStyle->{Directive[Blue,Thick],Directive[Red,Thick,Dashed]},
   PlotLegends->{"Field kicks","Bond kicks"},Frame->True,Axes->False,
   FrameLabel->{"quasienergy epsilon","normalized density rho(epsilon)"},
   PlotLabel->Row[{scanParameter," = ",N[scanValues[[i]]],
      ", T = ",N[spectra["Field"][[i]]["Period"]]}],
   PlotRange->{{-energyBound,energyBound},{0,All}},ImageSize->600],
   {i,dosIndices}];
Scan[Print,dosPlots];

(* ::Section:: *)
(* 8. PNG EXPORT AND ACCESS TO RESULTS *)

(* ::Subsection:: *)
(* Definitions *)
outputBaseDirectory[] := Module[{base},
   base = If[StringQ[$InputFileName] && $InputFileName != "",
      DirectoryName[$InputFileName],Quiet[Check[NotebookDirectory[],Directory[]]]];
   If[StringQ[base],base,Directory[]]];

(* ::Subsection:: *)
(* Computations *)
If[exportPlots,
   outputDirectory = FileNameJoin[{outputBaseDirectory[],
      "xxz_floquet_results","Lchain"<>ToString[Lchain]<>"_scan_"<>scanParameter}];
   If[!DirectoryQ[outputDirectory],
      CreateDirectory[outputDirectory,CreateIntermediateDirectories->True]];
   Export[FileNameJoin[{outputDirectory,"spectrum.png"}],spectrumFigure,
      ImageResolution->180];
   Do[Export[FileNameJoin[{outputDirectory,
      "DOS_panel"<>ToString[i]<>".png"}],dosPlots[[i]],ImageResolution->180],
      {i,Length[dosPlots]}];
   Print["PNG plots saved in: ",outputDirectory]];

(* Results:
   factorizationChecks; spectrumChecks; unitarityChecks; sweepSeconds;
   spectra["Field"][[i]], spectra["Bond"][[i]], pointResults[[i]]["Parameters"];
   Hstatic, UFField, UFBond, spectrumFigure, dosPlots.
   Eigenvectors are present only if saveEigenvectors=True.
   makeDiscreteDOS["Bond",i] constructs the full discrete DOS on demand.
   The no-kick protocols coincide exactly. DOS alone does not determine chaos.
   Existing eight kernels are reused; kernels remain open for the next run. *)

(* ::Section:: *)
(* 9. SYMMETRY SECTORS AND CIRCULAR MEAN GAP RATIOS AT THE BASELINE p *)

(* ::Subsection:: *)
(* Definitions: controls and symmetry operators *)
(* This section analyzes UFField and UFBond at p, not the complete sweep.
   It can be rerun by itself once the earlier baseline matrices exist.
   Known clear cases:
     XY/Heisenberg/ZZ bond pulses: total Mz conserved (wx=wy).
     XX or YY bond pulses: total z parity conserved, Mz generally broken.
     z-field pulses: total Mz conserved.
     x-field pulses with OmegaA=OmegaB=0: global x flip conserved.
     Equal endpoint parameters and equal pulses: reflection conserved.
   To preserve reflection, set OmegaB=OmegaA, gXYB=gXYA, gZB=gZA,
   KickB=KickA and, for Field, fieldAxisB=fieldAxisA; rerun the file.
   Special pulse angles can have extra Floquet symmetries; checks use UF itself.
   The candidate list is not a proof that no other symmetry exists.
   SU(2) cases, pure-ZZ attachments, and zero kicks have extra tested charges.
   Noncommuting symmetries are refined only within blocks they preserve. *)

symmetryTolerance = checkTolerance;
gapTolerance = 10 checkTolerance;
candidateNames = {"Mz","BulkMz","sAz","sBz","Pz","S2","Fx","Fy","Reflection"};

globalFlip[axis_] := KroneckerProduct @@
   ConstantArray[SparseArray[PauliMatrix[axis]],nspin];
makeReflection[] := SparseArray[
   Table[{FromDigits[Reverse[IntegerDigits[b,2,nspin]],2]+1,b+1}->1,
      {b,0,dim-1}],{dim,dim}];
relativeCommutatorError[s_,uf_] :=
   maxEntry[s.uf-uf.s]/Max[1,maxEntry[s]];

(* A sector basis B is a dim x d matrix with orthonormal columns.
   Only the symmetry matrices are diagonalized to obtain B.
   Ublock=B^dagger.UF.B is a restriction to an invariant subspace;
   no Floquet eigenvalues are discarded.
   All candidate charge eigenvalues are multiples of 1/4. Round below merely
   labels their known quantized values; Floquet phases are NEVER rounded. *)
refineSector[block_,name_,s_] := Module[
   {b=block["Basis"],sb,values,vectors,labels,groups},
   sb=ConjugateTranspose[b].(s.b);
   If[maxEntry[s.b-b.sb]/Max[1,maxEntry[s]]>symmetryTolerance,
      Return[{block}]];
   {values,vectors}=Eigensystem[N[sb,workingPrecision]];
   labels=Round[4 Re[values]]/4;
   groups=GatherBy[Range[Length[values]],labels[[#]]&];
   If[Length[groups]==1,Return[{block}]];
   Table[<|"Labels"->Join[block["Labels"],<|name->labels[[First[g]]]|>],
      "Basis"->b.Transpose[Orthogonalize[vectors[[g]]]]|>,{g,groups}]];

buildSectors[uf_,names_] := Module[
   {blocks={<|"Labels"-><||>,"Basis"->IdentityMatrix[dim]|>}},
   Do[blocks=Flatten[
      refineSector[#,name,symmetryOperators[name]]& /@ blocks,1],{name,names}];
   blocks];

circularRatios[eigenvalues_] := Module[{phases,gaps,next,ratios,d},
   d=Length[eigenvalues];
   phases=Sort[Mod[-Arg[eigenvalues]+Pi,2 Pi]-Pi];
   If[d<3,Return[<|"Phases"->phases,"Gaps"->{},"Ratios"->{},
      "MeanR"->Missing["TooFewLevels"],"MinimumGap"->Missing["TooFewLevels"],
      "Status"->"Fewer than 3 levels"|>]];
   gaps=Differences[Append[phases,First[phases]+2 Pi]];
   If[Min[gaps]<=gapTolerance,
      Return[<|"Phases"->phases,"Gaps"->gaps,"Ratios"->{},
         "MeanR"->Missing["DegenerateOrNearDegenerate"],
         "MinimumGap"->Min[gaps],
         "Status"->"Check degeneracy, unresolved symmetry, or precision"|>]];
   next=RotateLeft[gaps];
   ratios=MapThread[Min[#1,#2]/Max[#1,#2]&,{gaps,next}];
   <|"Phases"->phases,"Gaps"->gaps,"Ratios"->ratios,"MeanR"->Mean[ratios],
     "MinimumGap"->Min[gaps],"Status"->"OK"|>];

analyzeSector[job_] := Module[{uf,b,ub,eigenvalues,leak},
   uf=If[job["Protocol"]=="Field",UFField,UFBond];
   b=job["Basis"]; ub=ConjugateTranspose[b].(uf.b);
   leak=maxEntry[uf.b-b.ub];
   eigenvalues=Eigenvalues[N[ub,workingPrecision]];
   Join[<|"Protocol"->job["Protocol"],"Sector"->job["Labels"],
      "Dimension"->Length[ub],"Eigenvalues"->eigenvalues,
      "LeakageError"->leak,
      "UnitarityError"->maxEntry[ConjugateTranspose[ub].ub-
         IdentityMatrix[Length[ub]]]|>,circularRatios[eigenvalues]]];

(* ::Subsection:: *)
(* Computations: candidate checks and simultaneous block construction *)
totalX=Total[Table[single[localSpin["x"],i],{i,nspin}]];
totalY=Total[Table[single[localSpin["y"],i],{i,nspin}]];
totalZ=Total[Table[single[localSpin["z"],i],{i,nspin}]];
symmetryOperators=<|
   "Mz"->totalZ, "BulkMz"->totalZ-zA-zB,
   "sAz"->zA,"sBz"->zB, "Pz"->globalFlip[3],
   "S2"->totalX.totalX+totalY.totalY+totalZ.totalZ,
   "Fx"->globalFlip[1],"Fy"->globalFlip[2],"Reflection"->makeReflection[]|>;
symmetryErrors=Association@Table[kind->Association@Table[
   name->relativeCommutatorError[symmetryOperators[name],
      If[kind=="Field",UFField,UFBond]],{name,candidateNames}],{kind,protocols}];
conservedCandidates=Association@Table[kind->
   Select[candidateNames,symmetryErrors[kind][#]<=symmetryTolerance&],
   {kind,protocols}];
Print["Candidate [S,UF] errors: ",symmetryErrors];
Print["Conserved candidates: ",conservedCandidates];

{sectorBasisSeconds,sectorBases}=AbsoluteTiming[Association@Table[
   kind->buildSectors[If[kind=="Field",UFField,UFBond],conservedCandidates[kind]],
   {kind,protocols}]];
If[!TrueQ[And @@ Table[
      Total[Length[#["Basis"][[1]]]& /@ sectorBases[kind]]==dim,{kind,protocols}]],
   Print["Sector dimensions do not add to the full dimension."];Abort[]];

(* Computations: diagonalize independent symmetry blocks on the same 8 kernels *)
basisChecks=Association@Table[kind->Module[
   {q=Join[Sequence@@Lookup[sectorBases[kind],"Basis"],2]},
   maxEntry[ConjugateTranspose[q].q-IdentityMatrix[dim]]],{kind,protocols}];
Print["Complete sector-basis orthogonality errors: ",basisChecks];
If[Max[Values[basisChecks]]>10 symmetryTolerance,
   Print["Sector bases do not form a complete orthonormal basis."];Abort[]];
sectorJobs=Flatten[Table[
   Join[<|"Protocol"->kind|>,#]& /@ sectorBases[kind],{kind,protocols}],1];
DistributeDefinitions[analyzeSector,circularRatios,maxEntry,UFField,UFBond,
   workingPrecision,gapTolerance];
{sectorSeconds,sectorStatistics}=AbsoluteTiming[parallelMap[analyzeSector,sectorJobs]];
If[Max[Lookup[sectorStatistics,"LeakageError"]] > symmetryTolerance ||
   Max[Lookup[sectorStatistics,"UnitarityError"]] > 10 symmetryTolerance,
   Print["A proposed block is not invariant/unitary. Do not use its ratios."];Abort[]];

sectorTable=Prepend[
   ({#["Protocol"],#["Sector"],#["Dimension"],#["MeanR"],#["Status"]}& /@
      sectorStatistics),{"Protocol","Sector quantum numbers","Dimension","<r>","Status"}];
Print["Basis seconds: ",sectorBasisSeconds,"; sector diagonalization seconds: ",
   sectorSeconds];
Print[Grid[sectorTable,Frame->All,Alignment->Left]];
Print["Keep different-sector spectra separate. Small blocks have large fluctuations."];

(* Access: sectorStatistics[[i]]["Eigenvalues"], ["Phases"], ["Gaps"],
   ["Ratios"], ["MeanR"], ["LeakageError"]; sectorBases[kind][[i]]["Basis"].
   Include the gap across the 2 Pi phase boundary and pair gaps cyclically.
   Degeneracies are flagged, never removed or merged.
   Poisson reference ~0.386; COE ~0.53; CUE ~0.60 (large sectors).
   For the real x/z-field or real bond generators, the symmetric time frame
   exp(-I H0 T/2).exp(-I QA-I QB).exp(-I H0 T/2) is complex symmetric.
   A chaotic block preserved by that antiunitary symmetry is compared to COE.
   A conserved magnetization or reflection alone does not establish integrability.
   At zero kick strength, additional integrability/conservation can be relevant. *)

(* Optional CSV beside the PNGs. Raw floating-point means, no fitted/smoothed data. *)
If[exportPlots,
   Export[FileNameJoin[{outputDirectory,"baseline_sector_gap_ratios.csv"}],
      Prepend[({#["Protocol"],ToString[#["Sector"],InputForm],#["Dimension"],
         #["MeanR"],#["Status"]}& /@ sectorStatistics),
         {"Protocol","Sector","Dimension","MeanR","Status"}]]];


