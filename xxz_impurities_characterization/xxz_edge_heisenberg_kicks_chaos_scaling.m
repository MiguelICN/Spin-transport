(* ::Package:: *)

(* ::Title:: *)
(* Edge Heisenberg-kicked XXZ chain: ensemble-pooled spectral statistics in the largest symmetry sector *)


(* ::Section:: *)
(* 0. MODEL, DRIVE AND DIAGNOSTICS (definitions only) *)


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
     UF(alpha) = K(alpha) . Exp[-I H0 T],   K(alpha) = Exp[-I alpha V].
   One bond: Exp[-I alpha S.S] = c0 1 + c1 P_singlet, P_singlet = 1/4 - S.S,
     c0 = Exp[-I alpha/4], c1 = Exp[3 I alpha/4] - Exp[-I alpha/4].
   alpha -> alpha + 2 Pi gives UF -> -UF (quasienergies shift by Pi/T): alpha in [0, 2 Pi).

   ENSEMBLE (used for EVERY spectral statistic reported in the plots).
     For each nominal alpha0, R realizations alpha_i = alpha0 + Subdivide[-w, w, R - 1].
     All quasienergies, gap ratios and eigenvector entropies of the R realizations are pooled.
     w is either fixed or adaptive. Adaptive w uses the level velocities
       v_k = d theta_k / d alpha = <k| V |k>   (theta = -Arg[eigenvalue], exact for unitary UF)
     at alpha0 and sets the relative level displacement across the window,
       Delta = 2 w sd(v) / (2 Pi / D)   [in units of the mean level spacing],
     to targetDisplacement, capped by maxHalfWidth. Delta is printed on every plot: if it
     is < 1 the realizations are strongly correlated and the pool is effectively ~1 sample.

   DIAGNOSTICS (largest symmetry sector, dimension D per realization)
     DOS: pooled histogram of quasienergies eps = -Arg[lambda]/T in (-Pi/T, Pi/T].
     Entanglement: half-chain von Neumann entropy (natural log) of every Floquet eigenvector,
       cut after N_A = Floor[N/2] spins (impurity A + first N_A - 1 bulk spins), all realizations
       pooled; reference: random states of the same sector (mean +- sd).
     P(r): pooled circular gap ratios r_n = min(s_n, s_{n+1})/max(s_n, s_{n+1}) (no unfolding
       needed); surmises: Poisson 2/(1+r)^2, COE/CUE from Atas et al., PRL 110, 084101 (2013).
       <r> = pooled mean; error bar = sd of the R realization means / Sqrt[R].
       Large-D references: Poisson 0.3863, COE 0.5307, CUE 0.5996.
     SFF: K(n) = <|Tr U^n|^2>_ens / D, plotted vs tau = n / D (Heisenberg time n_H = D periods).
       Optional unfolding (sffUnfold): theta -> 2 Pi F(theta), F = cumulative of the
       ensemble-averaged DOS smoothed with a periodic Gaussian of width
       unfoldKernelSpacings mean spacings. The disconnected part |<Tr U^n>|^2 / D is drawn too.
       COE: K = 2 tau - tau Log[1 + 2 tau] (tau <= 1), 2 - tau Log[(2 tau + 1)/(2 tau - 1)]
       (tau > 1); Poisson: K = 1.
   Expected class if chaotic: COE (H0 and V are real; with the half-kick time origin UF is
   complex symmetric). *)


(* ::Section:: *)
(* 1. USER PARAMETERS *)


ClearAll["Global`*"];

Lvalues = Range[5, 10];                (* bulk chain lengths; N = L + 2 spins *)
nAlpha = 41;                           (* nominal alpha0 values in [0, 2 Pi) *)
alphaValues = N[2 Pi Range[0, nAlpha - 1]/nAlpha];
period = 1;                            (* T *)

parameterPoints = <|
   "default" -> <|"Label" -> "anisotropic point", "J" -> 1, "Delta" -> 1/2,
      "OmegaA" -> 3/5, "OmegaB" -> 3/5, "gXYA" -> 1/5, "gXYB" -> 2/5, "gZA" -> 1/10, "gZB" -> 1/10|>,
   "su2" -> <|"Label" -> "SU(2) point", "J" -> 1, "Delta" -> 1,
      "OmegaA" -> 0, "OmegaB" -> 0, "gXYA" -> 4/5, "gXYB" -> 13/10, "gZA" -> 4/5, "gZB" -> 13/10|>|>;
pointsToRun = {"default", "su2"};

(* ensemble *)
ensembleSize = 21;                     (* R; odd => alpha0 itself is one realization *)
windowMode = "Adaptive";               (* "Adaptive" or "Fixed" *)
targetDisplacement = 4;                (* adaptive: relative level motion across the window, in spacings *)
maxHalfWidth = 0.25;                   (* cap on the adaptive half-width w *)
fixedHalfWidth = 0.05;                 (* w used when windowMode = "Fixed" *)
pooledEntanglement = True;             (* False: eigenvectors only for the realization closest to alpha0 *)

(* spectral form factor *)
sffMaxTau = 3; sffLogBins = 40;
sffUnfold = True; unfoldKernelSpacings = 25; unfoldGrid = 4096;

(* other numerics *)
nRandomStates = 30; dosBins = 48; rBins = 25; entropyBins = 48;
fullCheckMaxSpins = 9; testAlphas = {1.234, 4.321}; tolerance = 10^-9;
useParallel = True; requestedKernels = 8;
exportPNG = True; imageResolution = 120;
outputRoot = FileNameJoin[{If[StringQ[$InputFileName] && $InputFileName != "",
      DirectoryName[$InputFileName], Quiet[Check[NotebookDirectory[], Directory[]]]],
    "xxz_heisenberg_kicks_chaos_scaling"}];


(* ::Section:: *)
(* 2. COLORS AND TEXT STYLES *)


(* Text inks; data color (one blue); a second accent only for the disconnected SFF part.
   Reference curves (Poisson/COE/CUE) are neutral ink distinguished by dash pattern.
   System-size ramp: one hue, light (smallest L) -> dark (largest L); the 6-step ramp is
   evenly spaced in OKLab between #86b6ef and #0d366b (validated: monotone lightness,
   adjacent Delta L >= 0.06, light end >= 2:1 contrast on white). *)
inkPrimary = RGBColor["#0b0b0b"]; inkSecondary = RGBColor["#52514e"]; inkMuted = RGBColor["#8a8986"];
dataBlue = RGBColor["#2a78d6"]; dataBlueLight = RGBColor["#9ec5f4"]; dataBlueDark = RGBColor["#184f95"];
accentOrange = RGBColor["#eb6834"]; bandGray = RGBColor["#d9d8d4"];
sizeRamp6 = RGBColor /@ {"#86b6ef", "#6d9bd3", "#5480b8", "#3c679e", "#254e84", "#0d366b"};
sizeColor[L_] := Module[{k = First[FirstPosition[Lvalues, L]], n = Length[Lvalues]},
   If[n == 6, sizeRamp6[[k]],
    ColorConvert[Blend[ColorConvert[{sizeRamp6[[1]], sizeRamp6[[-1]]}, "LAB"],
      If[n == 1, 1, (k - 1)/(n - 1)]], "RGB"]]];
refStyles = <|"Poisson" -> Directive[inkSecondary, AbsoluteThickness[1.6], AbsoluteDashing[{6, 4}]],
   "COE" -> Directive[inkPrimary, AbsoluteThickness[1.6]],
   "CUE" -> Directive[inkSecondary, AbsoluteThickness[1.6], AbsoluteDashing[{1.5, 3}]]|>;
refR = <|"Poisson" -> 2 Log[2.] - 1, "COE" -> 0.5307, "CUE" -> 0.5996|>;
baseStyle = {FontFamily -> "Helvetica", FontSize -> 12, FontColor -> inkPrimary};
frameStyle = Directive[inkSecondary, AbsoluteThickness[0.8]];


(* ::Section:: *)
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
globalProduct[a_, n_] := KroneckerProduct @@ ConstantArray[sigma[[a]], n];
reflection[n_] := SparseArray[Table[
    {FromDigits[Reverse[IntegerDigits[b, 2, n]], 2] + 1, b + 1} -> 1, {b, 0, 2^n - 1}], {2^n, 2^n}];
candidateOperators[n_] := Module[{s = totalSpin[n], fx = globalProduct[1, n],
    fy = globalProduct[2, n], r = reflection[n]},
   <|"Mz" -> s[[3]], "Sx_total" -> s[[1]], "S2" -> s[[1]].s[[1]] + s[[2]].s[[2]] + s[[3]].s[[3]],
     "Pz" -> globalProduct[3, n], "Fx" -> fx, "Fy" -> fy, "R" -> r, "R.Fx" -> r.fx, "R.Fy" -> r.fy|>];
maxAbs[m_SparseArray] := Max[0, Abs[m["NonzeroValues"]]];
maxAbs[m_] := Max[Abs[Flatten[m]]];
relativeCommutator[a_, b_] := maxAbs[a.b - b.a]/(Max[1, maxAbs[a]] Max[1, maxAbs[b]]);
magnetizationIndices[n_, m_] := Select[Range[2^n], Total[1/2 - IntegerDigits[# - 1, 2, n]] == m &];
multiplicity[n_, s_] := Binomial[n, n/2 - s] - If[n/2 - s - 1 >= 0, Binomial[n, n/2 - s - 1], 0];


(* ::Section:: *)
(* 4. SYMMETRY STATEMENT AND VERIFICATION *)


(* "default": conserved Mz and Pz (a function of Mz); broken Sx_total, S2 (Delta != 1, Omega != 0),
     Fx, Fy (Omega != 0), R, R.Fx, R.Fy (gXYA != gXYB). Largest sector M = 0 (N even), +1/2 (N odd).
   "su2": conserved SU(2) (Mz, Sx_total, S2; Fx, Fy are SU(2) elements); broken R, R.Fx, R.Fy
     (gA != gB). Largest sector (Sstar, M = Sstar), Sstar maximizing the multiplicity.
   The test: relative commutator [Q, H0] (exact rational arithmetic) and [Q, K(alpha)]
   (machine precision, at testAlphas). Zero (H0) or <= tolerance (K) = conserved. *)
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


(* ::Section:: *)
(* 5. LARGEST SECTOR *)


largestSector[pointName_, n_] := Module[
   {q = parameterPoints[pointName], type, label, m, idx, b, s2, sPlus, idxUp, vals, vecs,
    sel, sStar, h, hm, hs, pa, pb, pas, pbs, leak, checks, cands},
   type = If[pointName === "su2", "SU2", "Mz"];
   h = N[buildH0[q, n]];
   If[type === "Mz",
    m = If[EvenQ[n], 0, 1/2]; idx = magnetizationIndices[n, m];
    b = IdentityMatrix[Length[idx]]; label = "M = " <> ToString[m, InputForm];
    checks = <|"BasisOrthonormality" -> 0.|>,
    sStar = First[MaximalBy[Range[Mod[n, 2]/2, n/2], multiplicity[n, #] &]];
    m = sStar; idx = magnetizationIndices[n, m];
    cands = candidateOperators[n];
    s2 = Normal[N[cands["S2"][[idx, idx]]]];
    {vals, vecs} = Eigensystem[(s2 + Transpose[s2])/2];
    sel = Flatten[Position[Abs[vals - sStar (sStar + 1)], x_ /; x < 10^-8, {1}]];
    b = Transpose[Orthogonalize[vecs[[sel]]]];
    idxUp = magnetizationIndices[n, m + 1];
    sPlus = N[(cands["Sx_total"] + I totalSpin[n][[2]])[[idxUp, idx]]];
    label = "S = M = " <> ToString[sStar, InputForm];
    checks = <|"BasisOrthonormality" -> maxAbs[ConjugateTranspose[b].b - IdentityMatrix[Length[sel]]],
      "S2Eigen" -> maxAbs[s2.b - sStar (sStar + 1) b],
      "HighestWeight" -> If[Length[idxUp] > 0, maxAbs[sPlus.b], 0.],
      "DimensionVsMultiplicity" -> Length[sel] - multiplicity[n, sStar]|>];
   leak = <|"H0OutsideM" -> maxAbs[h[[Complement[Range[2^n], idx], idx]]]|>;
   hm = Normal[h[[idx, idx]]]; pa = Normal[N[singletA[n][[idx, idx]]]];
   pb = Normal[N[singletB[n][[idx, idx]]]];
   hs = ConjugateTranspose[b].hm.b; hs = (hs + ConjugateTranspose[hs])/2;
   pas = ConjugateTranspose[b].pa.b; pbs = ConjugateTranspose[b].pb.b;
   leak = Join[leak, <|"H0Sector" -> maxAbs[hm.b - b.hs], "PASector" -> maxAbs[pa.b - b.pas],
      "PBSector" -> maxAbs[pb.b - b.pbs]|>];
   <|"Point" -> pointName, "N" -> n, "L" -> n - 2, "Type" -> type, "Label" -> label,
     "Indices" -> idx, "Basis" -> b, "D" -> Length[hs], "U0" -> MatrixExp[-I period hs],
     "PA" -> pas, "PB" -> pbs, "PAB" -> pas.pbs,
     "V" -> IdentityMatrix[Length[hs]]/2 - pas - pbs,     (* S_A.S_1 + S_L.S_B in the sector *)
     "Checks" -> Join[checks, leak]|>];

kickSector[alpha_, sd_] := Module[{c = kickCoefficients[alpha]},
   c[[1]]^2 IdentityMatrix[sd["D"]] + c[[1]] c[[2]] (sd["PA"] + sd["PB"]) + c[[2]]^2 sd["PAB"]];

crossChecks[sd_] := Module[{n = sd["N"], q = parameterPoints[sd["Point"]], a = First[testAlphas],
    gen, kExp, out, full, sec},
   gen = N[(op2[localSS, 1, n] + op2[localSS, n - 1, n])[[sd["Indices"], sd["Indices"]]]];
   kExp = MatrixExp[-I a ConjugateTranspose[sd["Basis"]].Normal[gen].sd["Basis"]];
   out = <|"KickClosedForm" -> maxAbs[kExp - kickSector[a, sd]],
     "VSector" -> maxAbs[ConjugateTranspose[sd["Basis"]].Normal[gen].sd["Basis"] - sd["V"]]|>;
   If[n <= fullCheckMaxSpins,
    full = Eigenvalues[Normal[N[kickFull[a, n]]].MatrixExp[-I period Normal[N[buildH0[q, n]]]]];
    sec = Eigenvalues[kickSector[a, sd].sd["U0"]];
    out = Join[out, <|"SectorInFullSpectrum" -> Max[Min[Abs[full - #]] & /@ sec]|>]];
   out];


(* ::Section:: *)
(* 6. ONE NOMINAL ALPHA: ENSEMBLE, POOLED STATISTICS *)


gapRatios[theta_] := Module[{th = Sort[theta], g, pairs},
   g = Differences[Append[th, First[th] + 2 Pi]];
   pairs = Select[Transpose[{g, RotateLeft[g]}], Max[#] > 0 &];
   Min[#]/Max[#] & /@ pairs];

entropyOf[vec_, sd_] := Module[{n = sd["N"], nA, psi, sv, p},
   nA = Floor[n/2];
   psi = Normal[SparseArray[Thread[sd["Indices"] -> vec], 2^n]];
   sv = SingularValueList[ArrayReshape[psi, {2^nA, 2^(n - nA)}]];
   p = Select[sv^2, # > 10^-15 &];
   -Total[p Log[p]]];

(* one realization: Schur decomposition (orthonormal Floquet eigenvectors) *)
oneRealization[a_, sd_, withVectors_] := Module[{u, q, t, lam, eps, order, ent},
   u = kickSector[a, sd].sd["U0"];
   {q, t} = SchurDecomposition[u];
   lam = Diagonal[t]; eps = -Arg[lam]/period; order = Ordering[eps];
   ent = If[withVectors, Table[entropyOf[sd["Basis"].q[[All, j]], sd], {j, order}], {}];
   <|"alpha" -> a, "Eps" -> eps[[order]], "Theta" -> Mod[eps[[order]] period, 2 Pi],
     "Entropy" -> ent, "Velocity" -> Re[Diagonal[ConjugateTranspose[q].sd["V"].q]][[order]],
     "SchurOffDiagonal" -> maxAbs[t - DiagonalMatrix[lam]],
     "UnitarityError" -> maxAbs[ConjugateTranspose[u].u - IdentityMatrix[sd["D"]]]|>];

(* ensemble-averaged smoothed DOS on the circle and the unfolding map theta -> 2 Pi F(theta) *)
unfoldData[thetaPooled_, d_] := Module[{g = unfoldGrid, h, s, kernel, dens, cdf, centers, f},
   h = N[BinCounts[thetaPooled, {0, 2 Pi, 2 Pi/g}]];
   s = Max[1., unfoldKernelSpacings (2 Pi/d)/(2 Pi/g)];
   kernel = Table[Exp[-Min[k, g - k]^2/(2 s^2)], {k, 0, g - 1}]; kernel = kernel/Total[kernel];
   dens = Re[InverseFourier[Fourier[h] Fourier[kernel]]] Sqrt[g];   (* cyclic convolution *)
   dens = Clip[dens, {0, Infinity}];
   cdf = Accumulate[dens]/Total[dens];
   centers = (Range[g] - 1/2) 2 Pi/g;
   f = Interpolation[Transpose[{Prepend[Range[g] 2 Pi/g, 0.], Prepend[cdf, 0.]}], InterpolationOrder -> 1];
   <|"Map" -> f, "DensityTheta" -> Transpose[{centers, dens g/(2 Pi Total[dens])}]|>];

binnedMean[pairs_, nb_, lo_, hi_] := Module[{w = (hi - lo)/nb, groups},
   groups = GatherBy[pairs, Min[nb, Floor[(#[[1]] - lo)/w] + 1] &];
   SortBy[{Mean[#[[All, 1]]], Mean[#[[All, 2]]]} & /@ Select[groups, Length[#] >= 3 &], First]];

solveAlpha[alpha0_, sd_] := Module[{d = sd["D"], centre, sigmaV, w, alphas, reals, iCentre,
    thetaPooled, unf, thetaUnf, ns, traces, meanTr, rLists, rPooled, rMeans, disp, densEps},
   centre = oneRealization[alpha0, sd, True];
   sigmaV = StandardDeviation[centre["Velocity"]];
   w = If[windowMode === "Adaptive",
     Min[maxHalfWidth, targetDisplacement (2 Pi/d)/(2 Max[sigmaV, 10^-12])], fixedHalfWidth];
   alphas = If[ensembleSize == 1, {alpha0}, alpha0 + Subdivide[-w, w, ensembleSize - 1]];
   iCentre = First[Ordering[Abs[alphas - alpha0], 1]];
   reals = Table[If[k == iCentre && Abs[alphas[[k]] - alpha0] < 10^-14, centre,
      oneRealization[alphas[[k]], sd, pooledEntanglement || k == iCentre]], {k, Length[alphas]}];
   disp = 2 w sigmaV/(2 Pi/d);
   (* pooled gap ratios *)
   rLists = gapRatios /@ Lookup[reals, "Theta"];
   rPooled = Flatten[rLists]; rMeans = Mean /@ rLists;
   (* unfolding and SFF *)
   thetaPooled = Flatten[Lookup[reals, "Theta"]];
   unf = unfoldData[thetaPooled, d];
   thetaUnf = If[sffUnfold, 2 Pi unf["Map"] /@ # & /@ Lookup[reals, "Theta"], Lookup[reals, "Theta"]];
   ns = Range[Ceiling[sffMaxTau d]];
   traces = Total[Exp[-I Outer[Times, #, ns]]] & /@ thetaUnf;
   meanTr = Mean[traces];
   densEps = SortBy[{Mod[#[[1]] + Pi, 2 Pi] - Pi, #[[2]]} & /@ unf["DensityTheta"], First];
   densEps[[All, 1]] = densEps[[All, 1]]/period; densEps[[All, 2]] = densEps[[All, 2]] period;
   <|"alpha0" -> alpha0, "HalfWidth" -> w, "Alphas" -> alphas, "R" -> Length[alphas],
     "Displacement" -> disp, "LevelVelocitySD" -> sigmaV,
     "EpsPooled" -> Flatten[Lookup[reals, "Eps"]], "DensityCurve" -> densEps,
     "rPooled" -> rPooled, "MeanR" -> Mean[rPooled],
     "MeanRError" -> If[Length[rMeans] > 1, StandardDeviation[rMeans]/Sqrt[Length[rMeans]],
       StandardDeviation[rPooled]/Sqrt[Length[rPooled]]],
     "RealizationMeanR" -> rMeans,
     "EntropyPairs" -> Flatten[Transpose[{#["Eps"], #["Entropy"]}] & /@ Select[reals, #["Entropy"] =!= {} &], 1],
     "Tau" -> N[ns/d], "SFF" -> Mean[Abs[traces]^2]/d, "SFFDisconnected" -> Abs[meanTr]^2/d,
     "MaxSchurOffDiagonal" -> Max[Lookup[reals, "SchurOffDiagonal"]],
     "MaxUnitarityError" -> Max[Lookup[reals, "UnitarityError"]]|>];

randomReference[sd_] := Module[{vals},
   vals = Table[With[{c = RandomVariate[NormalDistribution[], sd["D"]] +
          I RandomVariate[NormalDistribution[], sd["D"]]},
       entropyOf[sd["Basis"].Normalize[c], sd]], {nRandomStates}];
   <|"Mean" -> Mean[vals], "SD" -> StandardDeviation[vals]|>];


(* ::Section:: *)
(* 7. PLOTS FOR ONE NOMINAL ALPHA *)


poissonR[r_] := 2/(1 + r)^2;
coeR[r_] := 27/4 (r + r^2)/(1 + r + r^2)^(5/2);
cueR[r_] := 81 Sqrt[3]/(2 Pi) (r + r^2)^2/(1 + r + r^2)^4;
coeSFF[tau_?NumericQ] := If[tau <= 1, 2 tau - tau Log[1 + 2 tau], 2 - tau Log[(2 tau + 1)/(2 tau - 1)]];
fmt[x_] := ToString[x, InputForm];
num[x_, d_] := ToString[NumberForm[N[x], {Max[1, Ceiling[Log10[Abs[N[x]] + 1]]] + d, d}]];

modelLine[pointName_] := Module[{q = parameterPoints[pointName]},
   Row[{"J = ", fmt[q["J"]], ",   \[CapitalDelta] = ", fmt[q["Delta"]],
     ",   \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(A\)]\) = ", fmt[q["OmegaA"]],
     ",   \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(B\)]\) = ", fmt[q["OmegaB"]],
     ",   \!\(\*SubscriptBox[\(g\), \(XY\)]\) (A, B) = (", fmt[q["gXYA"]], ", ", fmt[q["gXYB"]], ")",
     ",   \!\(\*SubscriptBox[\(g\), \(Z\)]\) (A, B) = (", fmt[q["gZA"]], ", ", fmt[q["gZB"]], ")",
     ",   T = ", fmt[period]}]];
systemLine[sd_] := Row[{parameterPoints[sd["Point"], "Label"], "   \[CenterDot]   L = ", sd["L"],
    " bulk spins (N = ", sd["N"], " incl. impurities A, B)   \[CenterDot]   largest sector ", sd["Label"],
    " (D = ", sd["D"], " levels per realization)"}];
ensembleLine[res_] := Row[{"pooled over R = ", res["R"], " kick strengths \[Alpha] \[Element] [",
    num[res["alpha0"] - res["HalfWidth"], 3], ", ", num[res["alpha0"] + res["HalfWidth"], 3],
    "]  (nominal \[Alpha] = ", num[res["alpha0"], 3], ", half-width w = ", num[res["HalfWidth"], 3],
    ");  levels move \[TildeTilde] ", num[res["Displacement"], 1], " mean spacings across the window"}];

titleBlock[what_, extra_, sd_, res_] := Column[Join[{
     Style[what, 15, Bold, inkPrimary],
     Style[systemLine[sd], 12, inkSecondary],
     Style[modelLine[sd["Point"]], 12, inkSecondary],
     Style[ensembleLine[res], 12, inkSecondary]},
    Style[#, 12, inkSecondary] & /@ extra], Spacings -> 0.25, Alignment -> Left];

commonOptions = {Frame -> True, Axes -> False, FrameStyle -> frameStyle, FrameTicksStyle -> Directive[11, inkSecondary],
   LabelStyle -> Directive[12, inkPrimary], ImageSize -> {720, 430}, AspectRatio -> Full,
   GridLines -> Automatic, GridLinesStyle -> Directive[GrayLevel[0.92]], BaseStyle -> baseStyle};

plotsFor[res_, sd_, ref_] := Module[{emax = Pi/period, nA = Floor[sd["N"]/2], dos, ent, pr, sff,
    kpts, kdis, entMean},
   (* DOS *)
   dos = Labeled[Legended[Show[
       Histogram[res["EpsPooled"], {-emax, emax, 2 emax/dosBins}, "PDF",
        ChartStyle -> Directive[EdgeForm[Directive[White, AbsoluteThickness[0.6]]], dataBlueLight],
        PlotRange -> {{-emax, emax}, {0, All}}, Evaluate[Sequence @@ commonOptions],
        FrameLabel -> {"quasienergy \[Epsilon]  (units of 1/T)", "density of states \[Rho](\[Epsilon])"}],
       ListLinePlot[res["DensityCurve"], PlotStyle -> Directive[dataBlueDark, AbsoluteThickness[2]]],
       Graphics[{refStyles["Poisson"], Line[{{-emax, period/(2 Pi)}, {emax, period/(2 Pi)}}]}]],
      {SwatchLegend[{dataBlueLight}, {"pooled histogram"}, LegendMarkerSize -> 14],
       LineLegend[{Directive[dataBlueDark, AbsoluteThickness[2]], refStyles["Poisson"]},
        {Row[{"smoothed DOS (Gaussian, ", unfoldKernelSpacings, " spacings)"}], "uniform T/(2\[Pi])"}]}],
     titleBlock["Quasienergy density of states", {}, sd, res], Top];
   (* entanglement *)
   entMean = binnedMean[res["EntropyPairs"], entropyBins, -emax, emax];
   ent = Labeled[Legended[Show[
       Graphics[{bandGray, Rectangle[{-emax, ref["Mean"] - ref["SD"]}, {emax, ref["Mean"] + ref["SD"]}]}],
       ListPlot[res["EntropyPairs"], PlotStyle -> Directive[dataBlue, Opacity[0.25], PointSize[0.005]]],
       ListLinePlot[entMean, PlotStyle -> Directive[dataBlueDark, AbsoluteThickness[2]]],
       Graphics[{{refStyles["Poisson"], Line[{{-emax, ref["Mean"]}, {emax, ref["Mean"]}}]},
         {refStyles["CUE"], Line[{{-emax, nA Log[2.]}, {emax, nA Log[2.]}}]}}],
       PlotRange -> {{-emax, emax}, {0, 1.05 nA Log[2.]}}, Evaluate[Sequence @@ commonOptions],
       FrameLabel -> {"quasienergy \[Epsilon]  (units of 1/T)",
         Row[{"entanglement entropy  \!\(\*SubscriptBox[\(S\), \(A\)]\)"}]}],
      {PointLegend[{dataBlue}, {"Floquet eigenstates (all realizations)"}, LegendMarkerSize -> 8],
       LineLegend[{Directive[dataBlueDark, AbsoluteThickness[2]], refStyles["Poisson"], refStyles["CUE"]},
        {Row[{"mean in ", entropyBins, " \[Epsilon]-bins"}], "random states of the same sector (\[PlusMinus] sd: gray band)",
         Row[{"maximum  ", nA, " ln 2"}]}]}],
     titleBlock["Half-chain entanglement entropy of the Floquet eigenstates",
      {Row[{"bipartition after ", nA, " spins: impurity A + first ", nA - 1, " bulk spins | rest;  ",
         "<\!\(\*SubscriptBox[\(S\), \(A\)]\)> / <\!\(\*SubscriptBox[\(S\), \(random\)]\)> = ",
         num[Mean[res["EntropyPairs"][[All, 2]]]/ref["Mean"], 3]}]}, sd, res], Top];
   (* P(r) *)
   pr = Labeled[Legended[Show[
       Histogram[res["rPooled"], {0, 1, 1/rBins}, "PDF",
        ChartStyle -> Directive[EdgeForm[Directive[White, AbsoluteThickness[0.6]]], dataBlueLight],
        PlotRange -> {{0, 1}, {0, 2.1}}, Evaluate[Sequence @@ commonOptions],
        FrameLabel -> {"gap ratio  r = min(\!\(\*SubscriptBox[\(s\), \(n\)]\), \!\(\*SubscriptBox[\(s\), \(n + 1\)]\)) / max(\!\(\*SubscriptBox[\(s\), \(n\)]\), \!\(\*SubscriptBox[\(s\), \(n + 1\)]\))", "probability density  P(r)"}],
       Plot[poissonR[x], {x, 0, 1}, PlotStyle -> refStyles["Poisson"]],
       Plot[coeR[x], {x, 0, 1}, PlotStyle -> refStyles["COE"]],
       Plot[cueR[x], {x, 0, 1}, PlotStyle -> refStyles["CUE"]]],
      {SwatchLegend[{dataBlueLight}, {Row[{"pooled data (", Length[res["rPooled"]], " ratios)"}]},
        LegendMarkerSize -> 14],
       LineLegend[{refStyles["Poisson"], refStyles["COE"], refStyles["CUE"]},
        {"Poisson  (<r> = 0.386)", "COE  (<r> = 0.531)", "CUE  (<r> = 0.600)"}]}],
     titleBlock["Distribution of consecutive level-spacing ratios",
      {Row[{"pooled  <r> = ", num[res["MeanR"], 3], " \[PlusMinus] ", num[res["MeanRError"], 3],
         "   (error: sd of the ", res["R"], " realization means / \[Sqrt]R)"}]}, sd, res], Top];
   (* SFF *)
   kpts = Select[Transpose[{res["Tau"], res["SFF"]}], #[[2]] > 0 &];
   kdis = Select[Transpose[{res["Tau"], res["SFFDisconnected"]}], #[[2]] > 0 &];
   sff = Labeled[Legended[Show[
       ListLogLogPlot[kpts, PlotStyle -> Directive[inkMuted, Opacity[0.35], PointSize[0.004]],
        PlotRange -> {{1/sd["D"], sffMaxTau}, All}, Evaluate[Sequence @@ commonOptions],
        FrameLabel -> {"\[Tau] = n / D   (n = number of periods; Heisenberg time n = D)",
          "spectral form factor  K(\[Tau]) / D"}],
       LogLogPlot[{coeSFF[x], 1}, {x, 1/sd["D"], sffMaxTau}, PlotStyle -> {refStyles["COE"], refStyles["Poisson"]}],
       ListLogLogPlot[Exp[binnedMean[Log[kdis], sffLogBins, Log[1./sd["D"]], Log[sffMaxTau]]],
        Joined -> True, PlotStyle -> Directive[accentOrange, AbsoluteThickness[1.6]]],
       ListLogLogPlot[Exp[binnedMean[Log[kpts], sffLogBins, Log[1./sd["D"]], Log[sffMaxTau]]],
        Joined -> True, PlotStyle -> Directive[dataBlue, AbsoluteThickness[2.2]]]],
      LineLegend[{Directive[dataBlue, AbsoluteThickness[2.2]], Directive[inkMuted, AbsoluteThickness[3]],
        Directive[accentOrange, AbsoluteThickness[1.6]], refStyles["COE"], refStyles["Poisson"]},
       {"<|Tr \!\(\*SuperscriptBox[\(U\), \(n\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)> / D  (log-binned)",
        "same, every n", "|<Tr \!\(\*SuperscriptBox[\(U\), \(n\)]\)>\!\(\*SuperscriptBox[\(|\), \(2\)]\) / D  (disconnected part)",
        "COE", "Poisson"}]],
     titleBlock["Spectral form factor in Heisenberg units",
      {If[sffUnfold, Row[{"phases unfolded with the ensemble-averaged smoothed DOS (Gaussian width ", unfoldKernelSpacings,
          " spacings); average over the R realizations"}], "no unfolding; average over the R realizations"]},
      sd, res], Top];
   <|"DOS" -> dos, "Entanglement" -> ent, "Pr" -> pr, "SFF" -> sff|>];

exportAlpha[res_, k_, sd_, ref_, dir_] := Module[{pl = plotsFor[res, sd, ref], tag},
   tag = sd["Point"] <> "_L" <> ToString[sd["L"]] <> "_a" <> IntegerString[k - 1, 10, 2] <>
     "_alpha" <> num[res["alpha0"], 4];
   KeyValueMap[Export[FileNameJoin[{dir, tag <> "_" <> #1 <> ".png"}], #2,
      ImageResolution -> imageResolution] &, pl]];

slim[res_, ref_] := <|KeyTake[res, {"alpha0", "HalfWidth", "R", "Displacement", "LevelVelocitySD", "MeanR",
      "MeanRError", "RealizationMeanR", "rPooled", "MaxSchurOffDiagonal", "MaxUnitarityError"}],
   "EntropyRatio" -> If[res["EntropyPairs"] === {}, Missing[], Mean[res["EntropyPairs"][[All, 2]]]/ref["Mean"]]|>;

processAlpha[k_] := Module[{res = solveAlpha[alphaValues[[k]], sectorShared]},
   If[exportPNG, exportAlpha[res, k, sectorShared, refShared, dirShared]];
   slim[res, refShared]];


(* ::Section:: *)
(* 8. RUN: SYMMETRIES, SECTORS, ALL ALPHAS (plots are produced on the worker kernels) *)


If[useParallel,
  If[$KernelCount < requestedKernels, LaunchKernels[requestedKernels - $KernelCount]];
  Print["Parallel kernels: ", $KernelCount]];

summary = <||>;
If[!DirectoryQ[outputRoot], CreateDirectory[outputRoot, CreateIntermediateDirectories -> True]];

Do[
  Do[
   n = L + 2;
   sym = verifySymmetries[pointName, n];
   Print[Style[parameterPoints[pointName, "Label"] <> ",  L = " <> ToString[L] <> " (N = " <> ToString[n] <> ")", 14, Bold]];
   Print[symmetryGrid[sym, pointName]];
   If[!sym["Consistent"], Print["Symmetry statement NOT confirmed. Stop."]; Abort[]];
   sd = largestSector[pointName, n]; cc = crossChecks[sd];
   Print["Largest sector ", sd["Label"], ", D = ", sd["D"], ";  sector checks: ", sd["Checks"],
    ";  cross-checks: ", cc];
   If[Max[Values[sd["Checks"]] /. x_Integer :> Abs[x]] > 10^-8 || Max[Values[cc]] > 10^-8,
    Print["Sector checks failed. Stop."]; Abort[]];
   ref = randomReference[sd];
   dir = FileNameJoin[{outputRoot, pointName, "L" <> ToString[L]}];
   If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
   sectorShared = sd; refShared = ref; dirShared = dir;
   If[useParallel,
    DistributeDefinitions[sectorShared, refShared, dirShared, alphaValues, period, parameterPoints,
     ensembleSize, windowMode, targetDisplacement, maxHalfWidth, fixedHalfWidth, pooledEntanglement,
     sffMaxTau, sffLogBins, sffUnfold, unfoldKernelSpacings, unfoldGrid, dosBins, rBins, entropyBins,
     exportPNG, imageResolution, inkPrimary, inkSecondary, inkMuted, dataBlue, dataBlueLight, dataBlueDark,
     accentOrange, bandGray, refStyles, baseStyle, frameStyle, commonOptions,
     processAlpha, solveAlpha, oneRealization, unfoldData, binnedMean, gapRatios, entropyOf, kickSector,
     kickCoefficients, maxAbs, plotsFor, exportAlpha, slim, titleBlock, systemLine, modelLine, ensembleLine,
     fmt, num, poissonR, coeR, cueR, coeSFF]];
   {sec, results} = AbsoluteTiming[If[useParallel,
      ParallelMap[processAlpha, Range[Length[alphaValues]], Method -> "FinestGrained", DistributedContexts -> None],
      processAlpha /@ Range[Length[alphaValues]]]];
   Print["  ", Length[results], " nominal alphas x R = ", ensembleSize, " realizations in ", Round[sec, 0.1],
    " s;  max Schur off-diagonal = ", Max[Lookup[results, "MaxSchurOffDiagonal"]],
    ";  max unitarity error = ", Max[Lookup[results, "MaxUnitarityError"]],
    ";  window half-width w in [", num[Min[Lookup[results, "HalfWidth"]], 3], ", ", num[Max[Lookup[results, "HalfWidth"]], 3],
    "];  level displacement across window in [", num[Min[Lookup[results, "Displacement"]], 1], ", ",
    num[Max[Lookup[results, "Displacement"]], 1], "] spacings"];
   Export[FileNameJoin[{dir, "results.wxf"}], <|"Sector" -> KeyDrop[sd, {"U0", "PA", "PB", "PAB", "Basis", "V"}],
      "Symmetry" -> sym, "CrossChecks" -> cc, "RandomReference" -> ref, "Results" -> results|>];
   summary[{pointName, L}] = <|"D" -> sd["D"], "Sector" -> sd["Label"], "Results" -> results|>,
   {L, Lvalues}],
  {pointName, pointsToRun}];


(* ::Section:: *)
(* 9. SCALING PLOTS (one hue; darker = larger L) *)


refLines[xmin_, xmax_, labelX_] := Table[{{refStyles[k], Line[{{xmin, refR[k]}, {xmax, refR[k]}}]},
     Text[Style[k, 11, inkSecondary], {labelX, refR[k]}, {-1, 0}]}, {k, {"Poisson", "COE", "CUE"}}];

scalingHeader[what_, pointName_, extra_] := Column[Join[{Style[what, 15, Bold, inkPrimary],
     Style[Row[{parameterPoints[pointName, "Label"], "   \[CenterDot]   largest symmetry sector   \[CenterDot]   L = ",
        First[Lvalues], " \[Ellipsis] ", Last[Lvalues], " bulk spins"}], 12, inkSecondary],
     Style[modelLine[pointName], 12, inkSecondary]}, Style[#, 12, inkSecondary] & /@ extra],
   Spacings -> 0.25, Alignment -> Left];

sizeLegend[pointName_] := LineLegend[Table[Directive[sizeColor[L], AbsoluteThickness[2]], {L, Lvalues}],
   Table[Row[{"L = ", L, "   (D = ", summary[{pointName, L}]["D"], ")"}], {L, Lvalues}],
   LegendMarkers -> Table[{Graphics[{sizeColor[L], Disk[]}], 9}, {L, Lvalues}],
   LegendLabel -> Style["bulk size", 12, inkSecondary]];

windowNote = If[windowMode === "Adaptive",
   Row[{"each point: R = ", ensembleSize, " realizations pooled in an adaptive \[Alpha]-window (target \[TildeTilde] ",
     targetDisplacement, " mean spacings of relative level motion, w \[LessEqual] ", maxHalfWidth, ")"}],
   Row[{"each point: R = ", ensembleSize, " realizations pooled in \[Alpha] \[PlusMinus] ", fixedHalfWidth}]];

Do[
  (* <r> vs alpha *)
  rPlot = Labeled[Legended[Show[
      Graphics[refLines[0, 2 Pi, 2 Pi + 0.05]],
      ListPlot[Table[With[{rs = summary[{pointName, L}]["Results"]},
         Transpose[{Lookup[rs, "alpha0"], MapThread[Around, {Lookup[rs, "MeanR"], Lookup[rs, "MeanRError"]}]}]],
        {L, Lvalues}], Joined -> True,
       PlotStyle -> Table[Directive[sizeColor[L], AbsoluteThickness[1.6]], {L, Lvalues}],
       PlotMarkers -> Table[{Graphics[{sizeColor[L], Disk[]}], 7}, {L, Lvalues}],
       IntervalMarkersStyle -> Table[Directive[sizeColor[L], AbsoluteThickness[1]], {L, Lvalues}]],
      PlotRange -> {{0, 2 Pi}, {0.25, 0.65}}, PlotRangePadding -> {{0.05, 0.6}, 0.01}, PlotRangeClipping -> False,
      Evaluate[Sequence @@ commonOptions], ImageSize -> {900, 480},
      FrameTicks -> {{Automatic, None}, {Table[{k Pi/4, If[k == 0, "0", k Pi/4]}, {k, 0, 8}], None}},
      FrameLabel -> {"kick strength  \[Alpha]", "mean gap ratio  <r>"}], sizeLegend[pointName]],
    scalingHeader["Mean gap ratio vs kick strength", pointName,
     {windowNote, "error bars: sd of the realization means / \[Sqrt]R;  dashed / solid / dotted: Poisson / COE / CUE"}], Top];
  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_meanR_vs_alpha.png"}], rPlot, ImageResolution -> 150];

  (* alpha-averaged <r> and entanglement vs L; pooled P(r) over all alpha != 0 per L *)
  scal = Table[With[{rs = summary[{pointName, L}]["Results"]},
     <|"L" -> L, "rAvg" -> Mean[Rest[Lookup[rs, "MeanR"]]],
       "rAvgErr" -> StandardDeviation[Rest[Lookup[rs, "MeanR"]]]/Sqrt[Length[rs] - 1],
       "rStatic" -> First[Lookup[rs, "MeanR"]], "rStaticErr" -> First[Lookup[rs, "MeanRError"]],
       "S" -> Mean[DeleteMissing[Rest[Lookup[rs, "EntropyRatio"]]]],
       "rAll" -> Flatten[Rest[Lookup[rs, "rPooled"]]]|>], {L, Lvalues}];
  panelR = Labeled[Legended[Show[
      Graphics[refLines[First[Lvalues] - 0.3, Last[Lvalues] + 0.3, Last[Lvalues] + 0.35]],
      Graphics[Table[{sizeColor[s["L"]], AbsolutePointSize[11], Point[{s["L"], s["rAvg"]}],
         AbsoluteThickness[1.2], Line[{{s["L"], s["rAvg"] - s["rAvgErr"]}, {s["L"], s["rAvg"] + s["rAvgErr"]}}],
         AbsolutePointSize[11], Point[{s["L"] + 0.15, s["rStatic"]}],
         White, AbsolutePointSize[6], Point[{s["L"] + 0.15, s["rStatic"]}]}, {s, scal}]],
      PlotRange -> {{First[Lvalues] - 0.5, Last[Lvalues] + 1.2}, {0.3, 0.65}}, Evaluate[Sequence @@ commonOptions],
      ImageSize -> {520, 400}, FrameLabel -> {"bulk size  L", "mean gap ratio  <r>"},
      FrameTicks -> {{Automatic, None}, {Lvalues, None}}],
     PointLegend[{inkSecondary, inkSecondary}, {Row[{"average over the ", nAlpha - 1, " values \[Alpha] \[NotEqual] 0"}], "static, \[Alpha] = 0 (open)"},
      LegendMarkers -> {{Graphics[Disk[]], 10}, {Graphics[{EdgeForm[inkSecondary], White, Disk[]}], 10}}]],
    Column[{Style["(a)  <r> vs system size", 13, Bold], Style[Row[{"error: sd over \[Alpha] / \[Sqrt]", nAlpha - 1, " (filled)"}], 11, inkSecondary]}], Top];
  panelS = Labeled[Show[
      Graphics[{refStyles["COE"], Line[{{First[Lvalues] - 0.3, 1}, {Last[Lvalues] + 0.3, 1}}],
        Text[Style["random states", 11, inkSecondary], {Last[Lvalues] + 0.35, 1}, {-1, 0}]}],
      Graphics[Table[{sizeColor[s["L"]], AbsolutePointSize[11], Point[{s["L"], s["S"]}]}, {s, scal}]],
      PlotRange -> {{First[Lvalues] - 0.5, Last[Lvalues] + 1.6}, {0, 1.05}}, Evaluate[Sequence @@ commonOptions],
      ImageSize -> {520, 400}, FrameLabel -> {"bulk size  L", "<\!\(\*SubscriptBox[\(S\), \(A\)]\)> / <\!\(\*SubscriptBox[\(S\), \(random\)]\)>"},
      FrameTicks -> {{Automatic, None}, {Lvalues, None}}],
    Column[{Style["(b)  eigenstate entanglement / random-state value", 13, Bold],
      Style["averaged over all eigenstates and all \[Alpha] \[NotEqual] 0", 11, inkSecondary]}], Top];
  panelP = Labeled[Legended[Show[
      Plot[poissonR[x], {x, 0, 1}, PlotStyle -> refStyles["Poisson"]],
      Plot[coeR[x], {x, 0, 1}, PlotStyle -> refStyles["COE"]],
      Plot[cueR[x], {x, 0, 1}, PlotStyle -> refStyles["CUE"]],
      ListLinePlot[Table[With[{h = HistogramList[s["rAll"], {0, 1, 1/rBins}, "PDF"]},
         Transpose[{MovingAverage[h[[1]], 2], h[[2]]}]], {s, scal}],
       PlotStyle -> Table[Directive[sizeColor[s["L"]], AbsoluteThickness[2]], {s, scal}]],
      PlotRange -> {{0, 1}, {0, 2.1}}, Evaluate[Sequence @@ commonOptions], ImageSize -> {520, 400},
      FrameLabel -> {"gap ratio  r", "P(r)"}], sizeLegend[pointName]],
    Column[{Style["(c)  P(r) pooled over all \[Alpha] \[NotEqual] 0", 13, Bold],
      Style["dashed / solid / dotted: Poisson / COE / CUE surmise", 11, inkSecondary]}], Top];
  scalPlot = Labeled[Grid[{{panelR, panelS, panelP}}, Spacings -> 3, Alignment -> Top],
    scalingHeader["Finite-size scaling of the spectral statistics", pointName, {windowNote}], Top];
  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_scaling_summary.png"}], scalPlot, ImageResolution -> 150];

  Export[FileNameJoin[{outputRoot, pointName, pointName <> "_meanR_table.csv"}],
   Prepend[Flatten[Table[With[{rs = summary[{pointName, L}]["Results"]},
       {L, summary[{pointName, L}]["D"], #["alpha0"], #["HalfWidth"], #["R"], #["Displacement"], #["MeanR"],
          #["MeanRError"], #["EntropyRatio"]} & /@ rs], {L, Lvalues}], 1],
    {"L", "D", "alpha0", "half_width_w", "R", "level_displacement_spacings", "meanR_pooled",
     "meanR_error", "entropy_ratio"}]];
  Print[rPlot]; Print[scalPlot],
  {pointName, pointsToRun}];
Print["Output written to ", outputRoot];
