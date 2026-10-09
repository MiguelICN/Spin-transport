(* ::Package:: *)

(* ::Title:: *)
(* Edge Heisenberg-kicked XXZ chain: spectral statistics in the largest symmetry sector *)


(* ::Section::Closed:: *)
(* 0. MODEL, DRIVE, SECTORS AND DIAGNOSTICS (definitions only) *)


(* Spins ordered {A, 1, ..., L, B};  N = L + 2;  S = sigma/2;  hbar = 1.
   Computational basis: site 1 (= impurity A) is the most significant bit; bit 0 = spin up;
   M = Sum_i (1/2 - bit_i).

   STATIC HAMILTONIAN
     H0 = J Sum_{j=1}^{L-1} [ S_j^x S_{j+1}^x + S_j^y S_{j+1}^y + Delta S_j^z S_{j+1}^z ]
        + OmegaA S_A^z + OmegaB S_B^z
        + gXYA (S_A^x S_1^x + S_A^y S_1^y) + gZA S_A^z S_1^z
        + gXYB (S_L^x S_B^x + S_L^y S_B^y) + gZB S_L^z S_B^z .

   DRIVE: simultaneous isotropic Heisenberg delta kicks on the two impurity bonds,
     H(t) = H0 + alpha Sum_n delta(t - n T) V,     V = S_A.S_1 + S_L.S_B .

   FLOQUET OPERATOR (exact; the two bond kicks act on disjoint sites and commute)
     UF(alpha, T) = K(alpha) . Exp[-I H0 T],
     K(alpha) = (c0 + c1 P_A)(c0 + c1 P_B),  P_X = singlet projector of bond X = 1/4 - S.S,
     c0 = Exp[-I alpha/4],  c1 = Exp[3 I alpha/4] - Exp[-I alpha/4].
   alpha -> alpha + 2 Pi gives UF -> -UF (quasienergies shift by Pi/T): alpha in [0, 2 Pi).

   FAST DIAGONALIZATION (same idea as the fixed-kick-basis trick of the QKT workflow).
     P_A and P_B commute, so in their joint eigenbasis W (labels nA, nB in {0, 1}) the kick is
     diagonal: K = diag[(c0 + c1 nA)(c0 + c1 nB)].  For each size (and each T) the fixed core
       C = W^dagger . Exp[-I H0 T] . W
     is built ONCE; for every alpha, UF in this basis is the row-scaled matrix
       diag[k(alpha)] . C      (no matrix exponential, no matrix product per alpha).
     The core is cleared after each size.

   SYMMETRY SECTORS (verified numerically in Section 4)
     "default" (Delta = 1/2, OmegaA = OmegaB != 0, gXYA != gXYB): only U(1) (Mz).
        Largest sector M = 0 (N even) or M = 1/2 (N odd):  D = Binomial[N, Floor[N/2]].
     "su2" (Delta = 1, Omega = 0, gXY = gZ at each end, gA != gB): full SU(2).
        Largest sector (Sstar, M = Sstar), Sstar maximizing d(N, S) = Binomial[N, N/2 - S] - Binomial[N, N/2 - S - 1]
        (ties -> larger S, whose M = S subspace is smaller and cheaper to resolve).

     L  :   5    6    7    8    9    10    11    12    13     14
     N  :   7    8    9   10   11    12    13    14    15     16
     default D:  35   70  126  252  462   924  1716  3432  6435  12870
     su2 D    :  14   28   48   90  165   297   572  1001  2002   3640
     D_default ~ 2^N Sqrt[2/(Pi N)];  D_su2 / D_default ~ 0.3 at N = 12-16 (slowly decreasing).
     D >= 500:  default L >= 10 (L = 9 gives 462);  su2 L >= 11.

   DIAGNOSTICS (largest sector, dimension D)
     Ensemble: for nominal alpha0, R realizations alpha_i = alpha0 + Subdivide[-w, w, R - 1];
       adaptive w from the level velocities v_k = d theta_k/d alpha = <k|V|k> at alpha0, such that
       levels move ~ targetDisplacement mean spacings relative to each other across the window.
     DOS, P(r), <r>, SFF: pooled over the R realizations (eigenvalues only).
     Entanglement: ONLY the central realization alpha0 (Schur eigenvectors); half-chain von
       Neumann entropy (natural log), cut after Floor[N/2] spins (impurity A + first bulk spins);
       reference: random states of the same sector (mean +- sd). No fitted or averaged lines.
     SFF (as in the QKT spectral workflow): per realization
         K_i(n) = |Sum_k Exp[-I n theta_k]|^2 / D,   n = 1 ... tauMax D,  tau = n / D,
       raw eigenphases (no unfolding), Gaussian-smoothed in tau (width sffSmoothingTau), then
       mean and standard error over realizations.  References: COE analytic
         K = 2 tau - tau Log[1 + 2 tau] (tau <= 1), 2 - tau Log[(2 tau + 1)/(2 tau - 1)] (tau > 1),
       and a numerical COE ensemble of the same dimension D processed identically.
       Caveat: the SFF uses raw phases, as for the QKT whose quasienergy DOS is flat. If the DOS
       plot is not flat (small W T), the early-time SFF contains the DOS contribution.
     <r>: circular gap ratios; references Poisson 0.3863, COE 0.5307, CUE 0.5996.

   SECTION 11 (T scan): SU(2) point, single L, fixed alpha (or fixed alpha/T), all diagnostics vs T,
     plus the overlap of the (symmetric-gauge) Floquet eigenvectors with the eigenvectors of the
     time-averaged Hamiltonian H_avg = H0 + (alpha/T) V. *)


(* ::Section:: *)
(* 1. USER PARAMETERS *)


ClearAll["Global`*"];

parameterPoints = <|
   "default" -> <|"Label" -> "anisotropic point", "J" -> 1, "Delta" -> 1/2,
      "OmegaA" -> 3/5, "OmegaB" -> 3/5, "gXYA" -> 1/5, "gXYB" -> 2/5, "gZA" -> 1/10, "gZB" -> 1/10|>,
   "su2" -> <|"Label" -> "SU(2) point", "J" -> 1, "Delta" -> 1,
      "OmegaA" -> 0, "OmegaB" -> 0, "gXYA" -> 4/5, "gXYB" -> 13/10, "gZA" -> 4/5, "gZB" -> 13/10|>|>;
pointsToRun = {"default", "su2"};

(* bulk sizes per point; chosen so that the largest-sector dimension is ~ 500 - 3500 *)
lValues = <|"default" -> {9, 10, 11, 12}, "su2" -> {10, 11, 12, 13}|>;

period = 4;                            (* T for the alpha scans: W T >> 2 Pi (genuinely Floquet) *)
nAlpha = 41;
alphaValues = N[2 Pi Range[0, nAlpha - 1]/nAlpha];

(* ensemble *)
ensembleSize = 21;                     (* R *)
windowMode = "Adaptive";               (* "Adaptive" or "Fixed" *)
targetDisplacement = 4;                (* adaptive: relative level motion across the window, in mean spacings *)
maxHalfWidth = 0.25;                   (* cap on the adaptive half-width w *)
fixedHalfWidth = 0.01;                 (* w used when windowMode = "Fixed" *)

(* spectral form factor *)
sffTauMax = 3.0; sffSmoothingTau = 0.01;
includeNumericalCOE = True; coeEnsembleSize = 20;

(* T scan (Section 11): SU(2) point, one size *)
runTScan = True;
tScanPoint = "su2"; tScanL = 12;                      (* D = 1001 *)
tScanMode = "FixedAlpha";                             (* "FixedAlpha" or "FixedAlphaOverT" *)
tScanAlpha = N[Pi/2];                                 (* used in "FixedAlpha" *)
tScanAlphaOverT = 0.4;                                (* used in "FixedAlphaOverT": alpha = a T (valid for |alpha| <= Pi) *)
tScanValues = N[Exp[Subdivide[Log[1/4], Log[8], 30]]];   (* 31 log-spaced periods, 0.25 ... 8 *)

(* other numerics *)
nRandomStates = 30; dosBins = 48; rBins = 25;
testAlphas = {1.234, 4.321}; tolerance = 10^-9; fullCheckMaxSpins = 9;
useParallel = True; requestedKernels = 8;   (* ~ 8 x (16 D^2) bytes for the shared core; reduce for D > 3000 *)
exportPNG = True; imageResolution = 120;
outputRoot = FileNameJoin[{If[StringQ[$InputFileName] && $InputFileName != "",
      DirectoryName[$InputFileName], Quiet[Check[NotebookDirectory[], Directory[]]]],
    "xxz_heisenberg_kicks_chaos_scaling"}];


(* ::Section::Closed:: *)
(* 2. COLORS AND TEXT STYLES *)


inkPrimary = RGBColor["#0b0b0b"]; inkSecondary = RGBColor["#52514e"]; inkMuted = RGBColor["#8a8986"];
dataBlue = RGBColor["#2a78d6"]; dataBlueLight = RGBColor["#9ec5f4"]; dataBlueDark = RGBColor["#184f95"];
bandGray = RGBColor["#d9d8d4"];
sizeRamp6 = RGBColor /@ {"#86b6ef", "#6d9bd3", "#5480b8", "#3c679e", "#254e84", "#0d366b"};
(* one hue, light (smallest L) -> dark (largest L); interpolated in LAB for other counts *)
sizeColor[L_, ls_List] := Module[{k = First[FirstPosition[ls, L]], n = Length[ls]},
   If[n == 6, sizeRamp6[[k]],
    ColorConvert[Blend[ColorConvert[{sizeRamp6[[1]], sizeRamp6[[-1]]}, "LAB"],
      If[n == 1, 1, (k - 1)/(n - 1)]], "RGB"]]];
refStyles = <|"Poisson" -> Directive[inkSecondary, AbsoluteThickness[1.6], AbsoluteDashing[{6, 4}]],
   "COE" -> Directive[inkPrimary, AbsoluteThickness[1.6]],
   "CUE" -> Directive[inkSecondary, AbsoluteThickness[1.6], AbsoluteDashing[{1.5, 3}]],
   "COEanalytic" -> Directive[inkPrimary, AbsoluteThickness[1.6], AbsoluteDashing[{6, 4}]],
   "COEnumeric" -> Directive[inkMuted, AbsoluteThickness[2]]|>;
refR = <|"Poisson" -> 2 Log[2.] - 1, "COE" -> 0.5307, "CUE" -> 0.5996|>;
baseStyle = {FontFamily -> "Helvetica", FontSize -> 12, FontColor -> inkPrimary};
frameStyle = Directive[inkSecondary, AbsoluteThickness[0.8]];
(* NOTE: with duplicated graphics options the FIRST occurrence wins, so ImageSize and GridLines
   are not part of commonOptions; they are given explicitly in every plot. *)
commonOptions = {Frame -> True, Axes -> False, FrameStyle -> frameStyle,
   FrameTicksStyle -> Directive[11, inkSecondary], LabelStyle -> Directive[12, inkPrimary],
   AspectRatio -> Full, BaseStyle -> baseStyle};
gridOptions = {GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.92]]};


(* ::Section::Closed:: *)
(* 3. OPERATORS *)


sigma = SparseArray /@ {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
spinHalf = sigma/2;
localXY = KroneckerProduct[spinHalf[[1]], spinHalf[[1]]] + KroneckerProduct[spinHalf[[2]], spinHalf[[2]]];
localZZ = KroneckerProduct[spinHalf[[3]], spinHalf[[3]]];
localSS = localXY + localZZ;
op1[m_, i_, n_] := KroneckerProduct[IdentityMatrix[2^(i - 1), SparseArray], m,
   IdentityMatrix[2^(n - i), SparseArray]];
op2[m_, i_, n_] := KroneckerProduct[IdentityMatrix[2^(i - 1), SparseArray], m,
   IdentityMatrix[2^(n - i - 1), SparseArray]];
buildH0[q_, n_] := q["J"] Sum[op2[localXY + q["Delta"] localZZ, i, n], {i, 2, n - 2}] +
   q["OmegaA"] op1[spinHalf[[3]], 1, n] + q["OmegaB"] op1[spinHalf[[3]], n, n] +
   q["gXYA"] op2[localXY, 1, n] + q["gZA"] op2[localZZ, 1, n] +
   q["gXYB"] op2[localXY, n - 1, n] + q["gZB"] op2[localZZ, n - 1, n];
singletA[n_] := IdentityMatrix[2^n, SparseArray]/4 - op2[localSS, 1, n];
singletB[n_] := IdentityMatrix[2^n, SparseArray]/4 - op2[localSS, n - 1, n];
kickCoefficients[alpha_] := {Exp[-I alpha/4], Exp[3 I alpha/4] - Exp[-I alpha/4]};
kickFull[alpha_, n_] := Module[{c = kickCoefficients[alpha], id = IdentityMatrix[2^n, SparseArray]},
   (c[[1]] id + c[[2]] singletB[n]) . (c[[1]] id + c[[2]] singletA[n])];
totalSpin[n_] := Table[Sum[op1[spinHalf[[a]], i, n], {i, n}], {a, 3}];
sPlusTotal[n_] := Sum[op1[SparseArray[{{0, 1}, {0, 0}}], i, n], {i, n}];   (* S^+ = Sx + I Sy, real *)
(* Cholesky-QR, applied twice: orthonormalizes the columns of a nearly orthonormal x to machine
   precision without any eigensolver (no degenerate-cluster issues). *)
cholOrthonormalize[x_] := Module[{y = x, g},
   Do[g = ConjugateTranspose[y] . y; g = (g + ConjugateTranspose[g])/2;
    y = Transpose[LinearSolve[Transpose[CholeskyDecomposition[g]], Transpose[y]]], {2}]; y];
(* Hermitian eigendecomposition via the Schur form: {eigenvalues, Q with eigenvectors as COLUMNS}.
   Q is unitary to machine precision even inside (near-)degenerate subspaces. Eigensystem is NOT
   used for Hermitian matrices: for the 1001 x 1001 matrix PA + 2 PB (eigenvalues 0..3, multiplicities
   526/200/200/75) it returned eigenvectors with max|Q^T Q - 1| = 1.0 (Mathematica 15.0.1). *)
hermitianEigensystem[h_] := Module[{q, t}, {q, t} = SchurDecomposition[(h + ConjugateTranspose[h])/2];
   {Re[Diagonal[t]], q}];
globalProduct[a_, n_] := KroneckerProduct @@ ConstantArray[sigma[[a]], n];
reflection[n_] := SparseArray[Table[
    {FromDigits[Reverse[IntegerDigits[b, 2, n]], 2] + 1, b + 1} -> 1, {b, 0, 2^n - 1}], {2^n, 2^n}];
candidateOperators[n_] := Module[{s = totalSpin[n], fx = globalProduct[1, n],
    fy = globalProduct[2, n], r = reflection[n]},
   <|"Mz" -> s[[3]], "Sx_total" -> s[[1]], "S2" -> s[[1]] . s[[1]] + s[[2]] . s[[2]] + s[[3]] . s[[3]],
     "Pz" -> globalProduct[3, n], "Fx" -> fx, "Fy" -> fy, "R" -> r, "R.Fx" -> r . fx, "R.Fy" -> r . fy|>];
maxAbs[m_SparseArray] := Max[0, Abs[m["NonzeroValues"]]];
maxAbs[m_] := Max[Abs[Flatten[m]]];
relativeCommutator[a_, b_] := maxAbs[a . b - b . a]/(Max[1, maxAbs[a]] Max[1, maxAbs[b]]);
magnetizationIndices[n_, m_] := Select[Range[2^n], Total[1/2 - IntegerDigits[# - 1, 2, n]] == m &];
multiplicity[n_, s_] := Binomial[n, n/2 - s] - If[n/2 - s - 1 >= 0, Binomial[n, n/2 - s - 1], 0];
largestSpin[n_] := Max[MaximalBy[Range[Mod[n, 2]/2, n/2], multiplicity[n, #] &]];
sectorDimension["default", L_] := Binomial[L + 2, Floor[(L + 2)/2]];
sectorDimension["su2", L_] := multiplicity[L + 2, largestSpin[L + 2]];


(* ::Section::Closed:: *)
(* 4. SYMMETRY STATEMENT AND VERIFICATION *)


(* "default": conserved Mz and Pz (a function of Mz); broken Sx_total, S2 (Delta != 1, Omega != 0),
     Fx, Fy (Omega != 0), R, R.Fx, R.Fy (gXYA != gXYB).
   "su2": conserved SU(2) (Mz, Sx_total, S2; Fx, Fy are SU(2) elements); broken R, R.Fx, R.Fy.
   Test: relative commutators [Q, H0] (exact rational arithmetic, exact 0 if conserved) and
   [Q, K(alpha)] (machine precision, ~1e-17 if conserved). Threshold: tolerance.
   For large N the full-space test is expensive; it is run up to verifyMaxSpins and the symmetry
   content does not depend on N. *)
verifyMaxSpins = 12;
declaredConserved = <|"default" -> {"Mz", "Pz"},
   "su2" -> {"Mz", "Sx_total", "S2", "Pz", "Fx", "Fy"}|>;
declaredBroken = <|"default" -> {"Sx_total", "S2", "Fx", "Fy", "R", "R.Fx", "R.Fy"},
   "su2" -> {"R", "R.Fx", "R.Fy"}|>;

verifySymmetries[pointName_, n_] := Module[{q = parameterPoints[pointName], h, ks, cands, table, ok},
   h = buildH0[q, n]; ks = N[kickFull[#, n]] & /@ testAlphas; cands = candidateOperators[n];
   table = AssociationMap[Function[name, <|
        "H0" -> relativeCommutator[cands[name], h],
        "K" -> Max[relativeCommutator[cands[name], #] & /@ ks]|>], Keys[cands]];
   ok = And @@ (Max[Values[table[#]]] < tolerance & /@ declaredConserved[pointName]) &&
      And @@ (Max[Values[table[#]]] > 10^-3 & /@ declaredBroken[pointName]);
   <|"Table" -> table, "Consistent" -> ok|>];

symmetryGrid[sym_, pointName_] := Grid[Prepend[KeyValueMap[{#1,
       If[MemberQ[declaredConserved[pointName], #1], Style["conserved", Darker[Green]], Style["broken", Red]],
       N[#2["H0"]], N[#2["K"]]} &, sym["Table"]],
    Style[#, Bold] & /@ {"operator", "declared", "rel. [Q, H0]", "rel. [Q, K]"}],
   Frame -> All, Alignment -> Left, BaseStyle -> baseStyle];

Print[Grid[Prepend[Table[{L, L + 2, sectorDimension["default", L], sectorDimension["su2", L],
      "S* = " <> ToString[largestSpin[L + 2], InputForm]}, {L, 5, 14}],
    Style[#, Bold] & /@ {"L", "N", "D default (M = 0 or 1/2)", "D su2 (S = M = S*)", ""}],
   Frame -> All, BaseStyle -> baseStyle]];


(* ::Section::Closed:: *)
(* 5. LARGEST SECTOR, KICK EIGENBASIS, FIXED CORE *)


(* Sector basis (columns, full-space coefficients restricted to the M-subspace 'Indices'),
   static H0, and the joint eigenbasis W of the two bond singlet projectors. *)
largestSector[pointName_, n_] := Module[
   {q = parameterPoints[pointName], type, label, m, idx, b, s2, sPlus, idxUp, vals, vecs,
    sel, sStar, h, hm, hs, pa, pb, pas, pbs, leak, checks, cands, jm, jv, jw, lab, groups, w, nA, nB},
   type = If[pointName === "su2", "SU2", "Mz"];
   h = N[buildH0[q, n]];
   If[type === "Mz",
    m = If[EvenQ[n], 0, 1/2]; idx = magnetizationIndices[n, m];
    b = IdentityMatrix[Length[idx]]; label = "M = " <> ToString[m, InputForm];
    checks = <|"BasisOrthonormality" -> 0.|>,
    (* S = M = Sstar states are exactly the highest-weight states of the M = Sstar subspace:
       basis = null space of S^+ : [M = Sstar] -> [M = Sstar + 1].  Obtained from the SVD (no degenerate
       eigensolver), then purified, b -> b - S^- [S^+ S^-]^-1 S^+ b, where S^+ S^- >= 2 Sstar + 2 on
       M = Sstar + 1 (well conditioned), and re-orthonormalized by Cholesky-QR. *)
    sStar = largestSpin[n]; m = sStar; idx = magnetizationIndices[n, m];
    idxUp = magnetizationIndices[n, m + 1];
    label = "S = M = " <> ToString[sStar, InputForm];
    If[idxUp === {},
     b = IdentityMatrix[Length[idx]] // N; sPlus = ConstantArray[0., {1, Length[idx]}],
     sPlus = Normal[N[sPlusTotal[n][[idxUp, idx]]]];
     {vals, sel, vecs} = SingularValueDecomposition[sPlus];
     sel = Count[Diagonal[sel], x_ /; x > 10^-6];                 (* rank of S^+ *)
     b = vecs[[All, sel + 1 ;;]]; Clear[vals, vecs];
     s2 = LinearSolve[sPlus . Transpose[sPlus]];                    (* factorized once *)
     Do[b = cholOrthonormalize[b - Transpose[sPlus] . s2[sPlus . b]], {2}]];
    s2 = totalSpin[n]; s2 = Normal[N[Total[# . # & /@ s2][[idx, idx]]]];
    checks = <|"BasisOrthonormality" -> maxAbs[ConjugateTranspose[b] . b - IdentityMatrix[Length[b[[1]]]]],
      "S2Eigen" -> maxAbs[s2 . b - sStar (sStar + 1) b],
      "HighestWeight" -> maxAbs[sPlus . b],
      "DimensionVsMultiplicity" -> Length[b[[1]]] - multiplicity[n, sStar]|>;
    Clear[s2, sPlus]];
   leak = <|"H0OutsideM" -> maxAbs[h[[Complement[Range[2^n], idx], idx]]]|>;
   hm = Normal[h[[idx, idx]]]; pa = Normal[N[singletA[n][[idx, idx]]]];
   pb = Normal[N[singletB[n][[idx, idx]]]];
   hs = ConjugateTranspose[b] . hm . b; hs = (hs + ConjugateTranspose[hs])/2;
   pas = ConjugateTranspose[b] . pa . b; pbs = ConjugateTranspose[b] . pb . b;
   pas = (pas + ConjugateTranspose[pas])/2; pbs = (pbs + ConjugateTranspose[pbs])/2;
   leak = Join[leak, <|"H0Sector" -> maxAbs[hm . b - b . hs], "PASector" -> maxAbs[pa . b - b . pas],
      "PBSector" -> maxAbs[pb . b - b . pbs]|>];
   (* joint eigenbasis of the commuting projectors: eigenvalue of PA + 2 PB = nA + 2 nB *)
   (* Schur vectors of PA + 2 PB, then each block is projected with its exact projector
      Q(nA,nB) = [PA or 1-PA].[PB or 1-PB] and re-orthonormalized (Cholesky-QR), twice. *)
   jm = pas + 2 pbs; {jv, jw} = hermitianEigensystem[jm]; jw = Transpose[jw];   (* rows = vectors *)
   lab = Round[jv]; groups = GatherBy[Range[Length[lab]], lab[[#]] &];
   w = Transpose[Join @@ (Module[{l = lab[[First[#]]], qa, qb, qq, x = Transpose[jw[[#]]]},
         qa = If[Mod[l, 2] == 1, pas, IdentityMatrix[Length[pas]] - pas];
         qb = If[Quotient[l, 2] == 1, pbs, IdentityMatrix[Length[pbs]] - pbs];
         qq = qa . qb; qq = (qq + ConjugateTranspose[qq])/2;
         Do[x = cholOrthonormalize[qq . x], {2}]; Transpose[x]] & /@ groups)];
   lab = Join @@ (lab[[#]] & /@ groups); Clear[jw];
   nA = Mod[lab, 2]; nB = Quotient[lab, 2];
   checks = Join[checks, <|"KickBasisUnitarity" -> maxAbs[ConjugateTranspose[w] . w - IdentityMatrix[Length[w]]],
      "PADiagonal" -> maxAbs[ConjugateTranspose[w] . pas . w - DiagonalMatrix[N[nA]]],
      "PBDiagonal" -> maxAbs[ConjugateTranspose[w] . pbs . w - DiagonalMatrix[N[nB]]]|>];
   <|"Point" -> pointName, "N" -> n, "L" -> n - 2, "Type" -> type, "Label" -> label,
     "Indices" -> idx, "BW" -> b . w, "D" -> Length[hs],
     "HW" -> ConjugateTranspose[w] . hs . w,          (* H0 in the kick eigenbasis *)
     "nA" -> nA, "nB" -> nB, "VDiag" -> N[1/2 - nA - nB],
     "Checks" -> Join[checks, leak]|>];

kickPhases[alpha_, sd_] := With[{c = kickCoefficients[alpha]}, (c[[1]] + c[[2]] sd["nA"]) (c[[1]] + c[[2]] sd["nB"])];
coreFromH[hw_, T_] := MatrixExp[-I T hw];                       (* built once per size / per T *)
floquetMatrix[alpha_, core_, sd_] := kickPhases[alpha, sd] core;  (* diag(k) . core, row scaling *)

crossChecks[sd_, core_, T_] := Module[{n = sd["N"], q = parameterPoints[sd["Point"]], a = First[testAlphas],
    out, full, sec, bw = sd["BW"], gen},
   gen = N[(op2[localSS, 1, n] + op2[localSS, n - 1, n])[[sd["Indices"], sd["Indices"]]]];
   out = <|"VDiagonal" -> maxAbs[ConjugateTranspose[bw] . Normal[gen] . bw - DiagonalMatrix[sd["VDiag"]]],
     "KickClosedForm" -> maxAbs[MatrixExp[-I a ConjugateTranspose[bw] . Normal[gen] . bw] -
        DiagonalMatrix[kickPhases[a, sd]]]|>;
   If[n <= fullCheckMaxSpins,
    full = Eigenvalues[Normal[N[kickFull[a, n]]] . MatrixExp[-I T Normal[N[buildH0[q, n]]]]];
    sec = Eigenvalues[floquetMatrix[a, core, sd]];
    out = Join[out, <|"SectorInFullSpectrum" -> Max[Min[Abs[full - #]] & /@ sec]|>]];
   out];


(* ::Section::Closed:: *)
(* 6. SPECTRAL TOOLS: GAP RATIOS, ENTANGLEMENT, SFF, COE REFERENCE *)


gapRatios[theta_] := Module[{th = Sort[Mod[theta, 2 Pi]], g, pairs},
   g = Differences[Append[th, First[th] + 2 Pi]];
   pairs = Select[Transpose[{g, RotateLeft[g]}], Max[#] > 10^-13 &];
   Min[#]/Max[#] & /@ pairs];

entropyOf[vec_, sd_] := Module[{n = sd["N"], nA, psi, sv, p},
   nA = Floor[n/2];
   psi = Normal[SparseArray[Thread[sd["Indices"] -> vec], 2^n]];
   sv = SingularValueList[ArrayReshape[psi, {2^nA, 2^(n - nA)}]];
   p = Select[sv^2, # > 10^-15 &];
   -Total[p Log[p]]];

(* discrete Floquet SFF of one realization: |Tr U^t|^2 / D, t = 1 ... tMax *)
sffDiscrete = Compile[{{phases, _Real, 1}, {tMax, _Integer}},
   Module[{n = Length[phases], z, power, out, t},
    z = Exp[-I phases]; power = ConstantArray[1.0 + 0.0 I, n]; out = ConstantArray[0.0, tMax];
    For[t = 1, t <= tMax, t++, power = power z; out[[t]] = Abs[Total[power]]^2/n];
    out], CompilationTarget -> "WVM", RuntimeOptions -> "Speed"];

(* mean and standard error over realizations of the Gaussian-smoothed SFF *)
processSFF[phaseLists_, d_] := Module[{tMax = Round[sffTauMax d], all, sm},
   all = sffDiscrete[#, tMax] & /@ phaseLists;
   sm = GaussianFilter[#, sffSmoothingTau d] & /@ all;
   <|"Tau" -> N[Range[tMax]/d], "SFFMean" -> Mean[sm],
     "SFFSE" -> If[Length[sm] > 1, StandardDeviation[sm]/Sqrt[Length[sm]], 0 sm[[1]]]|>];

(* numerical COE ensemble of dimension d: U = W^T W, W Haar (CUE) *)
coeMatrix[d_] := Module[{z, q, r, ph, w},
   z = (RandomVariate[NormalDistribution[], {d, d}] + I RandomVariate[NormalDistribution[], {d, d}])/Sqrt[2.];
   {q, r} = QRDecomposition[z]; ph = Diagonal[r]/Abs[Diagonal[r]];
   w = ConjugateTranspose[q] . DiagonalMatrix[ph];
   Transpose[w] . w];
coeReference[d_] := Module[{phs},
   phs = Table[Mod[-Arg[Eigenvalues[coeMatrix[d]]], 2 Pi], {coeEnsembleSize}];
   Join[processSFF[phs, d], <|"MeanR" -> Mean[Flatten[gapRatios /@ phs]], "Size" -> coeEnsembleSize|>]];

randomReference[sd_] := Module[{vals},
   vals = Table[With[{c = RandomVariate[NormalDistribution[], sd["D"]] +
          I RandomVariate[NormalDistribution[], sd["D"]]},
       entropyOf[sd["BW"] . Normalize[c], sd]], {nRandomStates}];
   <|"Mean" -> Mean[vals], "SD" -> StandardDeviation[vals]|>];


(* ::Section::Closed:: *)
(* 7. ONE PARAMETER POINT: ENSEMBLE AROUND alpha0 FOR A GIVEN CORE (= GIVEN T) *)


solveEnsemble[alpha0_, T_, core_, sd_, extra_ : <||>] := Module[{d = sd["D"], f, q, t, lam, theta, eps,
    order, vel, sigmaV, w, alphas, phaseLists, rLists, rPooled, rMeans, ent, sff, havg, phi, psiSym, ovl},
   (* central realization: Schur (orthonormal eigenvectors) *)
   f = floquetMatrix[alpha0, core, sd];
   {q, t} = SchurDecomposition[f];
   lam = Diagonal[t]; theta = Mod[-Arg[lam], 2 Pi]; eps = -Arg[lam]/T; order = Ordering[eps];
   vel = sd["VDiag"] . (Abs[q]^2);                       (* d theta_k / d alpha = <k|V|k> *)
   sigmaV = StandardDeviation[vel];
   w = If[windowMode === "Adaptive",
     Min[maxHalfWidth, targetDisplacement (2 Pi/d)/(2 Max[sigmaV, 10^-12])], fixedHalfWidth];
   alphas = If[ensembleSize == 1, {alpha0}, alpha0 + Subdivide[-w, w, ensembleSize - 1]];
   (* eigenvalues only for all realizations *)
   phaseLists = Table[If[Abs[a - alpha0] < 10^-14, Sort[theta],
       Sort[Mod[-Arg[Eigenvalues[floquetMatrix[a, core, sd]]], 2 Pi]]], {a, alphas}];
   rLists = gapRatios /@ phaseLists; rPooled = Flatten[rLists]; rMeans = Mean /@ rLists;
   (* entanglement: central realization only *)
   ent = Table[entropyOf[sd["BW"] . q[[All, k]], sd], {k, order}];
   sff = processSFF[phaseLists, d];
   (* optional: overlap with the eigenvectors of H_avg = H0 + (alpha0/T) V (symmetric gauge) *)
   ovl = If[TrueQ[extra["Overlap"]],
     havg = sd["HW"] + (alpha0/T) DiagonalMatrix[sd["VDiag"]]; havg = (havg + ConjugateTranspose[havg])/2;
     phi = Last[hermitianEigensystem[havg]];               (* columns *)
     psiSym = Conjugate[kickPhases[alpha0/2, sd]] q;     (* K^(-1/2) psi *)
     Mean[Max /@ Transpose[Abs[ConjugateTranspose[phi] . psiSym]^2]], Missing[]];
   <|"alpha0" -> alpha0, "T" -> T, "HalfWidth" -> w, "R" -> Length[alphas], "Alphas" -> alphas,
     "Displacement" -> 2 w sigmaV/(2 Pi/d), "LevelVelocitySD" -> sigmaV,
     "EpsPooled" -> Flatten[(Mod[# + Pi, 2 Pi] - Pi)/T & /@ phaseLists],
     "rPooled" -> rPooled, "MeanR" -> Mean[rPooled],
     "MeanRError" -> If[Length[rMeans] > 1, StandardDeviation[rMeans]/Sqrt[Length[rMeans]],
       StandardDeviation[rPooled]/Sqrt[Length[rPooled]]],
     "RealizationMeanR" -> rMeans,
     "CentralEps" -> eps[[order]], "CentralEntropy" -> ent,
     "Tau" -> sff["Tau"], "SFFMean" -> sff["SFFMean"], "SFFSE" -> sff["SFFSE"],
     "OverlapHavg" -> ovl,
     "SchurOffDiagonal" -> maxAbs[t - DiagonalMatrix[lam]],
     "UnitarityError" -> maxAbs[ConjugateTranspose[f] . f - IdentityMatrix[d]]|>];


(* ::Section::Closed:: *)
(* 8. PLOTS FOR ONE PARAMETER POINT *)


poissonR[r_] := 2/(1 + r)^2;
coeR[r_] := 27/4 (r + r^2)/(1 + r + r^2)^(5/2);
cueR[r_] := 81 Sqrt[3]/(2 Pi) (r + r^2)^2/(1 + r + r^2)^4;
coeSFF[tau_?NumericQ] := If[tau <= 1, 2 tau - tau Log[1 + 2 tau], 2 - tau Log[(2 tau + 1)/(2 tau - 1)]];
fmt[x_] := ToString[x, InputForm];
num[x_, d_] := ToString[NumberForm[N[x], {Max[1, Ceiling[Log10[Abs[N[x]] + 1]]] + d, d}]];

modelLine[pointName_, T_] := Module[{q = parameterPoints[pointName]},
   Row[{"J = ", fmt[q["J"]], ",   \[CapitalDelta] = ", fmt[q["Delta"]],
     ",   \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(A\)]\) = ", fmt[q["OmegaA"]],
     ",   \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(B\)]\) = ", fmt[q["OmegaB"]],
     ",   \!\(\*SubscriptBox[\(g\), \(XY\)]\) (A, B) = (", fmt[q["gXYA"]], ", ", fmt[q["gXYB"]], ")",
     ",   \!\(\*SubscriptBox[\(g\), \(Z\)]\) (A, B) = (", fmt[q["gZA"]], ", ", fmt[q["gZB"]], ")",
     ",   T = ", num[T, 3]}]];
systemLine[sd_] := Row[{parameterPoints[sd["Point"], "Label"], "   \[CenterDot]   L = ", sd["L"],
    " bulk spins (N = ", sd["N"], " incl. impurities A, B)   \[CenterDot]   largest sector ", sd["Label"],
    " (D = ", sd["D"], " levels per realization)"}];
ensembleLine[res_] := Row[{"R = ", res["R"], " kick strengths \[Alpha] \[Element] [",
    num[res["alpha0"] - res["HalfWidth"], 3], ", ", num[res["alpha0"] + res["HalfWidth"], 3],
    "]  (central \[Alpha] = ", num[res["alpha0"], 3], ", half-width w = ", num[res["HalfWidth"], 3],
    ");  levels move \[TildeTilde] ", num[res["Displacement"], 1], " mean spacings across the window"}];
titleBlock[what_, extra_, sd_, res_] := Column[Join[{
     Style[what, 15, Bold, inkPrimary],
     Style[systemLine[sd], 12, inkSecondary],
     Style[modelLine[sd["Point"], res["T"]], 12, inkSecondary],
     Style[ensembleLine[res], 12, inkSecondary]},
    Style[#, 12, inkSecondary] & /@ extra], Spacings -> 0.25, Alignment -> Left];

plotsFor[res_, sd_, ref_, coe_] := Module[{T = res["T"], emax, nA = Floor[sd["N"]/2], dos, ent, pr, sff,
    tau = res["Tau"], m = res["SFFMean"], se = res["SFFSE"], coeLines, legItems},
   emax = Pi/T;
   dos = Labeled[Legended[Show[
       Histogram[res["EpsPooled"], {-emax, emax, 2 emax/dosBins}, "PDF",
        ChartStyle -> Directive[EdgeForm[Directive[White, AbsoluteThickness[0.6]]], dataBlueLight],
        PlotRange -> {{-emax, emax}, {0, All}}, ImageSize -> {720, 430}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
        FrameLabel -> {"quasienergy \[Epsilon]  (units of 1/T)", "density of states \[Rho](\[Epsilon])"}],
       Graphics[{refStyles["Poisson"], Line[{{-emax, T/(2 Pi)}, {emax, T/(2 Pi)}}]}]],
      {SwatchLegend[{dataBlueLight}, {"all R realizations pooled"}, LegendMarkerSize -> 14],
       LineLegend[{refStyles["Poisson"]}, {"uniform T/(2\[Pi])"}]}],
     titleBlock["Quasienergy density of states", {}, sd, res], Top];
   ent = Labeled[Legended[Show[
       Graphics[{bandGray, Rectangle[{-emax, ref["Mean"] - ref["SD"]}, {emax, ref["Mean"] + ref["SD"]}]}],
       ListPlot[Transpose[{res["CentralEps"], res["CentralEntropy"]}],
        PlotStyle -> Directive[dataBlue, Opacity[0.6], PointSize[0.006]]],
       Graphics[{{refStyles["Poisson"], Line[{{-emax, ref["Mean"]}, {emax, ref["Mean"]}}]},
         {refStyles["CUE"], Line[{{-emax, nA Log[2.]}, {emax, nA Log[2.]}}]}}],
       PlotRange -> {{-emax, emax}, {0, 1.05 nA Log[2.]}}, ImageSize -> {720, 430}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
       FrameLabel -> {"quasienergy \[Epsilon]  (units of 1/T)",
         "entanglement entropy  \!\(\*SubscriptBox[\(S\), \(A\)]\)"}],
      {PointLegend[{dataBlue}, {Row[{"Floquet eigenstates at \[Alpha] = ", num[res["alpha0"], 3], " (single diagonalization)"}]},
        LegendMarkerSize -> 8],
       LineLegend[{refStyles["Poisson"], refStyles["CUE"]},
        {"random states of the same sector (\[PlusMinus] sd: gray band)", Row[{"maximum  ", nA, " ln 2"}]}]}],
     Column[{Style["Half-chain entanglement entropy of the Floquet eigenstates", 15, Bold, inkPrimary],
       Style[systemLine[sd], 12, inkSecondary], Style[modelLine[sd["Point"], T], 12, inkSecondary],
       Style[Row[{"single realization at \[Alpha] = ", num[res["alpha0"], 4], ";  bipartition after ", nA,
          " spins (impurity A + first ", nA - 1, " bulk spins);  <\!\(\*SubscriptBox[\(S\), \(A\)]\)> / <\!\(\*SubscriptBox[\(S\), \(random\)]\)> = ",
          num[Mean[res["CentralEntropy"]]/ref["Mean"], 3]}], 12, inkSecondary]}, Spacings -> 0.25], Top];
   pr = Labeled[Legended[Show[
       Histogram[res["rPooled"], {0, 1, 1/rBins}, "PDF",
        ChartStyle -> Directive[EdgeForm[Directive[White, AbsoluteThickness[0.6]]], dataBlueLight],
        PlotRange -> {{0, 1}, {0, 2.1}}, ImageSize -> {720, 430}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
        FrameLabel -> {"gap ratio  r = min(\!\(\*SubscriptBox[\(s\), \(n\)]\), \!\(\*SubscriptBox[\(s\), \(n + 1\)]\)) / max(\!\(\*SubscriptBox[\(s\), \(n\)]\), \!\(\*SubscriptBox[\(s\), \(n + 1\)]\))", "probability density  P(r)"}],
       Plot[poissonR[x], {x, 0, 1}, PlotStyle -> refStyles["Poisson"]],
       Plot[coeR[x], {x, 0, 1}, PlotStyle -> refStyles["COE"]],
       Plot[cueR[x], {x, 0, 1}, PlotStyle -> refStyles["CUE"]]],
      {SwatchLegend[{dataBlueLight}, {Row[{"R realizations pooled (", Length[res["rPooled"]], " ratios)"}]},
        LegendMarkerSize -> 14],
       LineLegend[{refStyles["Poisson"], refStyles["COE"], refStyles["CUE"]},
        {"Poisson  (<r> = 0.386)", "COE  (<r> = 0.531)", "CUE  (<r> = 0.600)"}]}],
     titleBlock["Distribution of consecutive level-spacing ratios",
      {Row[{"pooled  <r> = ", num[res["MeanR"], 3], " \[PlusMinus] ", num[res["MeanRError"], 3],
         "   (error: sd of the ", res["R"], " realization means / \[Sqrt]R)"}]}, sd, res], Top];
   coeLines = If[AssociationQ[coe], {ListLinePlot[Transpose[{coe["Tau"], coe["SFFMean"]}],
       PlotStyle -> refStyles["COEnumeric"]]}, {}];
   legItems = If[AssociationQ[coe],
     {{Directive[dataBlue, AbsoluteThickness[2]], Directive[dataBlue, Opacity[0.25], AbsoluteThickness[8]],
       refStyles["COEnumeric"], refStyles["COEanalytic"]},
      {"mean over R realizations", "\[PlusMinus] standard error",
       Row[{"numerical COE, D = ", sd["D"], " (", coe["Size"], " matrices)"}], "COE analytic"}},
     {{Directive[dataBlue, AbsoluteThickness[2]], Directive[dataBlue, Opacity[0.25], AbsoluteThickness[8]],
       refStyles["COEanalytic"]}, {"mean over R realizations", "\[PlusMinus] standard error", "COE analytic"}}];
   sff = Labeled[Legended[Show[
       ListLinePlot[{Transpose[{tau, m - se}], Transpose[{tau, m + se}]},
        PlotStyle -> {Directive[dataBlue, Opacity[0]], Directive[dataBlue, Opacity[0]]},
        Filling -> {1 -> {2}}, FillingStyle -> Directive[dataBlue, Opacity[0.18]],
        PlotRange -> {{0, sffTauMax}, {0, 1.6}}, ImageSize -> {720, 430}, Evaluate[Sequence @@ commonOptions],
        GridLines -> {{1}, None}, GridLinesStyle -> Directive[GrayLevel[0.75], Dotted],
        FrameLabel -> {"\[Tau] = n / D   (n = number of periods; Heisenberg time at \[Tau] = 1)",
          "spectral form factor  K(\[Tau])"}],
       Sequence @@ coeLines,
       ListLinePlot[Transpose[{tau, m}], PlotStyle -> Directive[dataBlue, AbsoluteThickness[2]]],
       Plot[coeSFF[x], {x, 0.001, sffTauMax}, PlotStyle -> refStyles["COEanalytic"]]],
      LineLegend @@ legItems],
     titleBlock["Spectral form factor",
      {Row[{"K(\[Tau]) = |Tr \!\(\*SuperscriptBox[\(U\), \(n\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\) / D per realization (raw eigenphases), Gaussian-smoothed in \[Tau] (width ",
         sffSmoothingTau, "), then averaged"}]}, sd, res], Top];
   <|"DOS" -> dos, "Entanglement" -> ent, "Pr" -> pr, "SFF" -> sff|>];

exportPlots[res_, tag_, sd_, ref_, coe_, dir_] := Module[{pl = plotsFor[res, sd, ref, coe]},
   KeyValueMap[Export[FileNameJoin[{dir, tag <> "_" <> #1 <> ".png"}], #2,
      ImageResolution -> imageResolution] &, pl]];

slim[res_, ref_] := <|KeyTake[res, {"alpha0", "T", "HalfWidth", "R", "Displacement", "LevelVelocitySD",
      "MeanR", "MeanRError", "RealizationMeanR", "rPooled", "OverlapHavg", "SchurOffDiagonal", "UnitarityError"}],
   "EntropyRatio" -> Mean[res["CentralEntropy"]]/ref["Mean"]|>;


(* ::Section::Closed:: *)
(* 8b. RUN SETUP (needed by Sections 9 AND 11; evaluate before either) *)


If[useParallel,
  If[$KernelCount < requestedKernels, LaunchKernels[requestedKernels - $KernelCount]];
  ParallelEvaluate[$HistoryLength = 0];
  Print["Parallel kernels: ", $KernelCount]];

workerFunctions := DistributeDefinitions[solveEnsemble, hermitianEigensystem, floquetMatrix, kickPhases, kickCoefficients,
   gapRatios, entropyOf, sffDiscrete, processSFF, maxAbs, plotsFor, exportPlots, slim, titleBlock,
   systemLine, modelLine, ensembleLine, fmt, num, poissonR, coeR, cueR, coeSFF,
   parameterPoints, ensembleSize, windowMode, targetDisplacement, maxHalfWidth, fixedHalfWidth,
   sffTauMax, sffSmoothingTau, dosBins, rBins, exportPNG, imageResolution, gridOptions, inkPrimary, inkSecondary,
   inkMuted, dataBlue, dataBlueLight, dataBlueDark, bandGray, refStyles, baseStyle, frameStyle, commonOptions];

releaseShared[] := (Clear[coreShared, sectorShared, refShared, coeShared, dirShared];
   If[useParallel, ParallelEvaluate[Clear[coreShared, sectorShared, refShared, coeShared, dirShared];
     ClearSystemCache[]]]; ClearSystemCache[]);

If[!AssociationQ[coeCache], coeCache = <||>];       (* numerical COE references, keyed by D *)
If[!AssociationQ[summary], summary = <||>];
If[!DirectoryQ[outputRoot], CreateDirectory[outputRoot, CreateIntermediateDirectories -> True]];


(* ::Section:: *)
(* 9. RUN THE alpha SCANS (both points, all sizes, period T = period). Set pointsToRun = {} to skip. *)


Do[
  Do[
   n = L + 2;
   Print[Style[parameterPoints[pointName, "Label"] <> ",  L = " <> ToString[L] <> " (N = " <> ToString[n] <>
      "),  T = " <> num[period, 3], 14, Bold]];
   If[n <= verifyMaxSpins,
    sym = verifySymmetries[pointName, n]; Print[symmetryGrid[sym, pointName]];
    If[!sym["Consistent"], Print["Symmetry statement NOT confirmed. Stop."]; Abort[]],
    sym = Missing["NotRun"]];
   sd = largestSector[pointName, n];
   core = coreFromH[sd["HW"], period];                      (* built once for this size *)
   cc = crossChecks[sd, core, period];
   Print["Largest sector ", sd["Label"], ", D = ", sd["D"], ";  checks: ", sd["Checks"], ";  cross-checks: ", cc];
   If[Max[Values[sd["Checks"]] /. x_Integer :> Abs[x]] > 10^-8 || Max[Values[cc]] > 10^-8,
    Print["Sector checks failed. Stop."]; Abort[]];
   ref = randomReference[sd];
   If[includeNumericalCOE && !KeyExistsQ[coeCache, sd["D"]],
    Print["  numerical COE reference, D = ", sd["D"], " ..."];
    coeCache[sd["D"]] = coeReference[sd["D"]];
    Print["  COE check: <r> = ", num[coeCache[sd["D"]]["MeanR"], 4], " (expected 0.5307)"]];
   dir = FileNameJoin[{outputRoot, pointName, "L" <> ToString[L]}];
   If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
   coreShared = core; sectorShared = KeyDrop[sd, {"HW"}]; refShared = ref;
   coeShared = If[includeNumericalCOE, coeCache[sd["D"]], None]; dirShared = dir;
   If[useParallel, workerFunctions; DistributeDefinitions[coreShared, sectorShared, refShared, coeShared,
      dirShared, alphaValues, period]];
   job = Function[k, Module[{res = solveEnsemble[alphaValues[[k]], period, coreShared, sectorShared]},
      If[exportPNG, exportPlots[res, sectorShared["Point"] <> "_L" <> ToString[sectorShared["L"]] <>
         "_a" <> IntegerString[k - 1, 10, 2] <> "_alpha" <> num[res["alpha0"], 4],
         sectorShared, refShared, coeShared, dirShared]];
      slim[res, refShared]]];
   {sec, results} = AbsoluteTiming[If[useParallel,
      ParallelMap[job, Range[Length[alphaValues]], Method -> "FinestGrained", DistributedContexts -> None],
      job /@ Range[Length[alphaValues]]]];
   Print["  ", Length[results], " alphas x R = ", ensembleSize, " in ", Round[sec, 0.1],
    " s;  max Schur off-diagonal = ", Max[Lookup[results, "SchurOffDiagonal"]],
    ";  max unitarity error = ", Max[Lookup[results, "UnitarityError"]],
    ";  w in [", num[Min[Lookup[results, "HalfWidth"]], 3], ", ", num[Max[Lookup[results, "HalfWidth"]], 3],
    "];  level displacement in [", num[Min[Lookup[results, "Displacement"]], 1], ", ",
    num[Max[Lookup[results, "Displacement"]], 1], "] spacings"];
   Export[FileNameJoin[{dir, "results.wxf"}], <|"Sector" -> KeyDrop[sd, {"HW", "BW"}],
      "Symmetry" -> sym, "CrossChecks" -> cc, "RandomReference" -> ref, "Results" -> results|>];
   summary[{pointName, L}] = <|"D" -> sd["D"], "Sector" -> sd["Label"], "Results" -> results|>;
   Clear[core, sd]; releaseShared[],               (* delete the core of this size *)
   {L, lValues[pointName]}],
  {pointName, pointsToRun}];


(* ::Section::Closed:: *)
(* 10. SCALING PLOTS (one hue; darker = larger L) *)


refLines[xmin_, xmax_, labelX_] := Table[{{refStyles[k], Line[{{xmin, refR[k]}, {xmax, refR[k]}}]},
    Text[Style[k, 11, inkSecondary], {labelX, refR[k]}, {-1, 0}]}, {k, {"Poisson", "COE", "CUE"}}];
scalingHeader[what_, pointName_, extra_] := Column[Join[{Style[what, 15, Bold, inkPrimary],
     Style[Row[{parameterPoints[pointName, "Label"], "   \[CenterDot]   largest symmetry sector   \[CenterDot]   L = ",
        Row[lValues[pointName], ", "], " bulk spins"}], 12, inkSecondary],
     Style[modelLine[pointName, period], 12, inkSecondary]}, Style[#, 12, inkSecondary] & /@ extra],
   Spacings -> 0.25, Alignment -> Left];
sizeLegend[pointName_] := With[{ls = lValues[pointName]},
   LineLegend[Table[Directive[sizeColor[L, ls], AbsoluteThickness[2]], {L, ls}],
    Table[Row[{"L = ", L, "   (D = ", summary[{pointName, L}]["D"], ")"}], {L, ls}],
    LegendMarkers -> Table[{Graphics[{sizeColor[L, ls], Disk[]}], 9}, {L, ls}],
    LegendLabel -> Style["bulk size", 12, inkSecondary]]];
windowNote = If[windowMode === "Adaptive",
   Row[{"each point: R = ", ensembleSize, " realizations pooled in an adaptive \[Alpha]-window (target \[TildeTilde] ",
     targetDisplacement, " mean spacings of relative level motion, w \[LessEqual] ", maxHalfWidth, ")"}],
   Row[{"each point: R = ", ensembleSize, " realizations pooled in \[Alpha] \[PlusMinus] ", fixedHalfWidth}]];

Do[
  ls = lValues[pointName];
  rPlot = Labeled[Legended[Show[
      Graphics[refLines[0, 2 Pi, 2 Pi + 0.05]],
      ListPlot[Table[With[{rs = summary[{pointName, L}]["Results"]},
         Transpose[{Lookup[rs, "alpha0"], MapThread[Around, {Lookup[rs, "MeanR"], Lookup[rs, "MeanRError"]}]}]],
        {L, ls}], Joined -> True,
       PlotStyle -> Table[Directive[sizeColor[L, ls], AbsoluteThickness[1.6]], {L, ls}],
       PlotMarkers -> Table[{Graphics[{sizeColor[L, ls], Disk[]}], 7}, {L, ls}],
       IntervalMarkersStyle -> Table[Directive[sizeColor[L, ls], AbsoluteThickness[1]], {L, ls}]],
      PlotRange -> {{0, 2 Pi}, {0.25, 0.65}}, PlotRangePadding -> {{0.05, 0.6}, 0.01}, PlotRangeClipping -> False,
      Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions], ImageSize -> {900, 480},
      FrameTicks -> {{Automatic, None}, {Table[{k Pi/4, If[k == 0, "0", k Pi/4]}, {k, 0, 8}], None}},
      FrameLabel -> {"kick strength  \[Alpha]", "mean gap ratio  <r>"}], sizeLegend[pointName]],
    scalingHeader["Mean gap ratio vs kick strength", pointName,
     {windowNote, "error bars: sd of the realization means / \[Sqrt]R;  dashed / solid / dotted: Poisson / COE / CUE"}], Top];
  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_meanR_vs_alpha.png"}], rPlot, ImageResolution -> 150];
  scal = Table[With[{rs = summary[{pointName, L}]["Results"]},
     <|"L" -> L, "rAvg" -> Mean[Rest[Lookup[rs, "MeanR"]]],
       "rAvgErr" -> StandardDeviation[Rest[Lookup[rs, "MeanR"]]]/Sqrt[Length[rs] - 1],
       "rStatic" -> First[Lookup[rs, "MeanR"]], "S" -> Mean[Rest[Lookup[rs, "EntropyRatio"]]],
       "rAll" -> Flatten[Rest[Lookup[rs, "rPooled"]]]|>], {L, ls}];
  panelR = Labeled[Legended[Show[
      Graphics[refLines[First[ls] - 0.3, Last[ls] + 0.3, Last[ls] + 0.35]],
      Graphics[Table[{sizeColor[s["L"], ls], AbsolutePointSize[11], Point[{s["L"], s["rAvg"]}],
         AbsoluteThickness[1.2], Line[{{s["L"], s["rAvg"] - s["rAvgErr"]}, {s["L"], s["rAvg"] + s["rAvgErr"]}}],
         AbsolutePointSize[11], Point[{s["L"] + 0.15, s["rStatic"]}],
         White, AbsolutePointSize[6], Point[{s["L"] + 0.15, s["rStatic"]}]}, {s, scal}]],
      PlotRange -> {{First[ls] - 0.5, Last[ls] + 1.2}, {0.3, 0.65}}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
      ImageSize -> {520, 400}, FrameLabel -> {"bulk size  L", "mean gap ratio  <r>"},
      FrameTicks -> {{Automatic, None}, {ls, None}}],
     PointLegend[{inkSecondary, inkSecondary}, {Row[{"average over the ", nAlpha - 1, " values \[Alpha] \[NotEqual] 0"}],
       "\[Alpha] = 0 (open): U = Exp[-I H0 T]"},
      LegendMarkers -> {{Graphics[Disk[]], 10}, {Graphics[{EdgeForm[inkSecondary], White, Disk[]}], 10}}]],
    Column[{Style["(a)  <r> vs system size", 13, Bold],
      Style[Row[{"error: sd over \[Alpha] / \[Sqrt]", nAlpha - 1, " (filled)"}], 11, inkSecondary]}], Top];
  panelS = Labeled[Show[
      Graphics[{refStyles["COE"], Line[{{First[ls] - 0.3, 1}, {Last[ls] + 0.3, 1}}],
        Text[Style["random states", 11, inkSecondary], {Last[ls] + 0.35, 1}, {-1, 0}]}],
      Graphics[Table[{sizeColor[s["L"], ls], AbsolutePointSize[11], Point[{s["L"], s["S"]}]}, {s, scal}]],
      PlotRange -> {{First[ls] - 0.5, Last[ls] + 1.6}, {0, 1.05}}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
      ImageSize -> {520, 400}, FrameLabel -> {"bulk size  L", "<\!\(\*SubscriptBox[\(S\), \(A\)]\)> / <\!\(\*SubscriptBox[\(S\), \(random\)]\)>"},
      FrameTicks -> {{Automatic, None}, {ls, None}}],
    Column[{Style["(b)  eigenstate entanglement / random-state value", 13, Bold],
      Style["central realization at each \[Alpha] \[NotEqual] 0, averaged over \[Alpha]", 11, inkSecondary]}], Top];
  panelP = Labeled[Legended[Show[
      Plot[poissonR[x], {x, 0, 1}, PlotStyle -> refStyles["Poisson"]],
      Plot[coeR[x], {x, 0, 1}, PlotStyle -> refStyles["COE"]],
      Plot[cueR[x], {x, 0, 1}, PlotStyle -> refStyles["CUE"]],
      ListLinePlot[Table[With[{h = HistogramList[s["rAll"], {0, 1, 1/rBins}, "PDF"]},
         Transpose[{MovingAverage[h[[1]], 2], h[[2]]}]], {s, scal}],
       PlotStyle -> Table[Directive[sizeColor[s["L"], ls], AbsoluteThickness[2]], {s, scal}]],
      PlotRange -> {{0, 1}, {0, 2.1}}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions], ImageSize -> {520, 400},
      FrameLabel -> {"gap ratio  r", "P(r)"}], sizeLegend[pointName]],
    Column[{Style["(c)  P(r) pooled over all \[Alpha] \[NotEqual] 0", 13, Bold],
      Style["dashed / solid / dotted: Poisson / COE / CUE surmise", 11, inkSecondary]}], Top];
  scalPlot = Labeled[Grid[{{panelR, panelS, panelP}}, Spacings -> 3, Alignment -> Top],
    scalingHeader["Finite-size scaling of the spectral statistics", pointName, {windowNote}], Top];
  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_scaling_summary.png"}], scalPlot, ImageResolution -> 150];
  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_meanR_table.csv"}],
   Prepend[Flatten[Table[With[{rs = summary[{pointName, L}]["Results"]},
       {L, summary[{pointName, L}]["D"], #["alpha0"], #["HalfWidth"], #["R"], #["Displacement"], #["MeanR"],
          #["MeanRError"], #["EntropyRatio"]} & /@ rs], {L, ls}], 1],
    {"L", "D", "alpha0", "half_width_w", "R", "level_displacement_spacings", "meanR_pooled",
     "meanR_error", "entropy_ratio_central"}]];
  Print[rPlot]; Print[scalPlot],
  {pointName, pointsToRun}];


(* ::Section::Closed:: *)
(* 11. T SCAN: SU(2) POINT, ONE SIZE, ALL DIAGNOSTICS AS A FUNCTION OF THE PERIOD T *)


(* Fixed: model parameters of tScanPoint, L = tScanL, and either alpha (tScanMode = "FixedAlpha")
   or alpha/T (tScanMode = "FixedAlphaOverT", i.e. fixed time-averaged Hamiltonian; only meaningful
   while |alpha| <= Pi because alpha is defined modulo 2 Pi).
   H0 is diagonalized ONCE in the kick eigenbasis; for each T the core is Q.Exp[-I E T].Q^dagger. *)
If[runTScan,
  n = tScanL + 2;
  Print[Style["T scan: " <> parameterPoints[tScanPoint, "Label"] <> ",  L = " <> ToString[tScanL] <>
     ",  mode " <> tScanMode, 14, Bold]];
  sdT = largestSector[tScanPoint, n];
  {eT, qT} = hermitianEigensystem[sdT["HW"]];                       (* columns = eigenvectors of H0 *)
  ccT = Join[crossChecks[sdT, qT . (Exp[-I eT Last[tScanValues]] ConjugateTranspose[qT]), Last[tScanValues]],
    <|"H0EigenUnitarity" -> maxAbs[ConjugateTranspose[qT] . qT - IdentityMatrix[Length[qT]]],
      "CoreVsMatrixExp" -> maxAbs[qT . (Exp[-I eT Last[tScanValues]] ConjugateTranspose[qT]) -
         coreFromH[sdT["HW"], Last[tScanValues]]]|>];
  Print["Largest sector ", sdT["Label"], ", D = ", sdT["D"], ";  checks: ", sdT["Checks"], ";  cross-checks: ", ccT];
  If[Max[Values[sdT["Checks"]] /. x_Integer :> Abs[x]] > 10^-8 || Max[Values[ccT]] > 10^-8,
   Print["Sector checks failed. Stop."]; Abort[]];
  widthT = Max[eT] - Min[eT];
  refT = randomReference[sdT];
  If[includeNumericalCOE && !KeyExistsQ[coeCache, sdT["D"]], coeCache[sdT["D"]] = coeReference[sdT["D"]]];
  dirT = FileNameJoin[{outputRoot, "Tscan_" <> tScanPoint <> "_L" <> ToString[tScanL]}];
  If[!DirectoryQ[dirT], CreateDirectory[dirT, CreateIntermediateDirectories -> True]];
  alphaOfT[T_] := If[tScanMode === "FixedAlphaOverT", tScanAlphaOverT T, tScanAlpha];
  If[tScanMode === "FixedAlphaOverT" && Max[Abs[alphaOfT /@ tScanValues]] > Pi,
   Print["Warning: alpha = a T exceeds Pi for some T; H_avg is then not uniquely defined."]];
  sectorShared = sdT;                                   (* includes HW, needed for the overlap *)
  refShared = refT; coeShared = If[includeNumericalCOE, coeCache[sdT["D"]], None]; dirShared = dirT;
  qShared = qT; eShared = eT;
  If[useParallel, workerFunctions; DistributeDefinitions[sectorShared, refShared, coeShared, dirShared,
     qShared, eShared, tScanValues, alphaOfT, tScanMode, tScanAlpha, tScanAlphaOverT]];
  jobT = Function[k, Module[{T = tScanValues[[k]], core, res},
     core = qShared . (Exp[-I eShared T] ConjugateTranspose[qShared]);
     res = solveEnsemble[alphaOfT[T], T, core, sectorShared, <|"Overlap" -> True|>];
     If[exportPNG, exportPlots[res, "Tscan_" <> sectorShared["Point"] <> "_L" <> ToString[sectorShared["L"]] <>
        "_t" <> IntegerString[k - 1, 10, 2] <> "_T" <> num[T, 3], sectorShared, refShared, coeShared, dirShared]];
     slim[res, refShared]]];
  {secT, resultsT} = AbsoluteTiming[If[useParallel,
     ParallelMap[jobT, Range[Length[tScanValues]], Method -> "FinestGrained", DistributedContexts -> None],
     jobT /@ Range[Length[tScanValues]]]];
  Print["  ", Length[resultsT], " periods in ", Round[secT, 0.1], " s;  max Schur off-diagonal = ",
   Max[Lookup[resultsT, "SchurOffDiagonal"]], ";  bandwidth W of H0 in the sector = ", num[widthT, 3]];

  tHeader[what_, extra_] := Column[Join[{Style[what, 15, Bold, inkPrimary],
      Style[Row[{parameterPoints[tScanPoint, "Label"], "   \[CenterDot]   L = ", tScanL, " (N = ", n,
         ")   \[CenterDot]   largest sector ", sdT["Label"], " (D = ", sdT["D"], ")   \[CenterDot]   ",
         If[tScanMode === "FixedAlpha", Row[{"fixed \[Alpha] = ", num[tScanAlpha, 3]}],
          Row[{"fixed \[Alpha]/T = ", num[tScanAlphaOverT, 3]}]]}], 12, inkSecondary],
      Style[Row[{"J = 1,  \[CapitalDelta] = ", fmt[parameterPoints[tScanPoint, "Delta"]], ",  \!\(\*SubscriptBox[\(g\), \(XY\)]\) = \!\(\*SubscriptBox[\(g\), \(Z\)]\) (A, B) = (",
         fmt[parameterPoints[tScanPoint, "gXYA"]], ", ", fmt[parameterPoints[tScanPoint, "gXYB"]], "),  \[CapitalOmega] = 0;  bandwidth of H0 in the sector W = ",
         num[widthT, 2]}], 12, inkSecondary]}, Style[#, 12, inkSecondary] & /@ extra], Spacings -> 0.25, Alignment -> Left];
  wLine[ymin_, ymax_] := Graphics[{inkMuted, AbsoluteDashing[{2, 3}], Line[{{Log[2 Pi/widthT], ymin}, {Log[2 Pi/widthT], ymax}}],
     Text[Style["W T = 2\[Pi]", 11, inkSecondary], {Log[2 Pi/widthT], ymax}, {-1.1, 1}]}];
  (* the log-x plot must come first in Show so that its T tick labels are used *)
  panelTR = Labeled[Show[
     ListLogLinearPlot[Transpose[{Lookup[resultsT, "T"], MapThread[Around, {Lookup[resultsT, "MeanR"], Lookup[resultsT, "MeanRError"]}]}],
      Joined -> True, PlotStyle -> Directive[dataBlue, AbsoluteThickness[1.6]], PlotMarkers -> {Graphics[{dataBlue, Disk[]}], 7}],
     Graphics[refLines[Log[Min[tScanValues]], Log[Max[tScanValues]], Log[Max[tScanValues]] + 0.05]],
     wLine[0.3, 0.65],
     PlotRange -> {{Log[Min[tScanValues]], Log[Max[tScanValues]]}, {0.3, 0.65}}, PlotRangePadding -> {{0.05, 0.5}, 0.01},
     PlotRangeClipping -> False, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions], ImageSize -> {560, 400},
     FrameLabel -> {"period  T", "mean gap ratio  <r>"}],
    Column[{Style["(a)  <r> vs T", 13, Bold], Style[Row[{"R = ", ensembleSize, " realizations pooled; error: sd of realization means / \[Sqrt]R"}], 11, inkSecondary]}], Top];
  panelTO = Labeled[Show[
     ListLogLinearPlot[Transpose[{Lookup[resultsT, "T"], Lookup[resultsT, "OverlapHavg"]}],
      Joined -> True, PlotStyle -> Directive[dataBlue, AbsoluteThickness[1.6]], PlotMarkers -> {Graphics[{dataBlue, Disk[]}], 7}],
     wLine[0, 1.02],
     PlotRange -> {{Log[Min[tScanValues]], Log[Max[tScanValues]]}, {0, 1.02}}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
     ImageSize -> {560, 400}, FrameLabel -> {"period  T", "overlap with \!\(\*SubscriptBox[\(H\), \(avg\)]\) eigenstates"}],
    Column[{Style["(b)  how static is the Floquet operator?", 13, Bold],
      Style["mean over Floquet eigenstates (symmetric gauge) of max |<\[Phi]|\[Psi]>\!\(\*SuperscriptBox[\(|\), \(2\)]\), \!\(\*SubscriptBox[\(H\), \(avg\)]\) = \!\(\*SubscriptBox[\(H\), \(0\)]\) + (\[Alpha]/T) V", 11, inkSecondary]}], Top];
  panelTS = Labeled[Show[
     ListLogLinearPlot[Transpose[{Lookup[resultsT, "T"], Lookup[resultsT, "EntropyRatio"]}],
      Joined -> True, PlotStyle -> Directive[dataBlue, AbsoluteThickness[1.6]], PlotMarkers -> {Graphics[{dataBlue, Disk[]}], 7}],
     Graphics[{refStyles["COE"], Line[{{Log[Min[tScanValues]], 1}, {Log[Max[tScanValues]], 1}}]}],
     wLine[0, 1.05],
     PlotRange -> {{Log[Min[tScanValues]], Log[Max[tScanValues]]}, {0, 1.05}}, Evaluate[Sequence @@ commonOptions], Evaluate[Sequence @@ gridOptions],
     ImageSize -> {560, 400}, FrameLabel -> {"period  T", "<\!\(\*SubscriptBox[\(S\), \(A\)]\)> / <\!\(\*SubscriptBox[\(S\), \(random\)]\)>"}],
    Column[{Style["(c)  eigenstate entanglement / random-state value", 13, Bold],
      Style["central realization at each T; solid line: random states", 11, inkSecondary]}], Top];
  tPlot = Labeled[Grid[{{panelTR, panelTO, panelTS}}, Spacings -> 3, Alignment -> Top],
    tHeader["Spectral statistics vs the period T", {windowNote}], Top];
  Export[FileNameJoin[{dirT, "Tscan_summary.png"}], tPlot, ImageResolution -> 150];
  Export[FileNameJoin[{dirT, "Tscan_table.csv"}],
   Prepend[{#["T"], #["alpha0"], #["HalfWidth"], #["R"], #["Displacement"], #["MeanR"], #["MeanRError"],
       #["OverlapHavg"], #["EntropyRatio"], widthT #["T"]/(2 Pi)} & /@ resultsT,
    {"T", "alpha", "half_width_w", "R", "level_displacement_spacings", "meanR_pooled", "meanR_error",
     "overlap_Havg", "entropy_ratio_central", "W_T_over_2pi"}]];
  Export[FileNameJoin[{dirT, "results.wxf"}], <|"Sector" -> KeyDrop[sdT, {"HW", "BW"}], "Bandwidth" -> widthT,
     "RandomReference" -> refT, "Results" -> resultsT|>];
  Print[tPlot];
  Clear[sdT, eT, qT]; Clear[qShared, eShared]; releaseShared[]];

Print["Output written to ", outputRoot];
