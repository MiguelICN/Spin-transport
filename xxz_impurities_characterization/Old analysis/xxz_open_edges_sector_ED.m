(* ::Title:: *)
(* Open XXZ chain: end spins as system, interior spins as environment *)

(* ::Section:: *)
(* 1. Model, conventions, and free parameters *)
(* Sites are {1,...,L}, L>=3; S^a = PauliMatrix[a]/2, hbar=1.
   H = J Sum[Sx[j] Sx[j+1]+Sy[j] Sy[j+1]+Delta Sz[j] Sz[j+1],j=1..L-1]
       + h Sum[Sz[j],j=1..L] + bLeft Sz[1] + bRight Sz[L].
   No bond connects L to 1. No additional probe spins are attached.
   The homogeneous, field-free XXZ model has only {L,J,Delta}; set all fields
   to zero. The optional fields are longitudinal, real, and static.
   J is SIGNED: at Delta=1, J>0 is antiferromagnetic, J<0 ferromagnetic.
   At general Delta, the signed couplings are Jxy=J and Jz=J Delta.
   Changing J alone reverses only the exchange part if fields are fixed.
   Reversing J,h,bLeft,bRight together gives -H and reverses energy order.
   Do not infer the ground state at arbitrary Delta from the sign of J alone.

   System S={1,L}, dimension 4. Environment E={2,...,L-1}, dimension 2^(L-2).
   H_S=(h+bLeft) Sz[1]+(h+bRight) Sz[L].
   H_E has interior bonds j=2..L-2 and interior uniform fields.
   H_SE consists of bonds (1,2) and (L-1,L), with the SAME J and Delta.
   For L=3 the environment has one spin and no interior bond.
   This partition alone specifies neither an initial environment state nor
   a dynamical map. A channel additionally requires a fixed rho_E and an
   initially factorized rho_S tensor rho_E (in the reordered S,E basis).

   ED below means full finite-matrix diagonalization, numerical to machine
   precision, with every eigenpair retained; it is not symbolic radical algebra.
   Charges resolved: Sz_total, reflection R when allowed, spin inversion F in
   m=0 when allowed, and RF in m=0 when allowed. F maps m to -m, so F and Sz
   cannot be simultaneous labels outside m=0. At Delta=1, additionally resolve
   S_total^2 if the field is uniform. At Delta=-1 use staggered S_total^2.
   Accidental degeneracies and extra integrable conserved quantities are not
   claimed to be a complete classification of all parameter-dependent charges.
*)

ClearAll[XXZSpin, XXZHamiltonian, XXZCharges, XXZSplitCharge,
  XXZSectorED, XXZEdgeDensity, XXZRun, XXZValidate, XXZZeroQ];

Options[XXZHamiltonian] = {"UniformField" -> 0, "LeftField" -> 0,
  "RightField" -> 0};
XXZHamiltonian::params = "Use integer L>=3 and finite real J, Delta, h, bLeft, bRight.";
XXZSpin[l_Integer, a_Integer, j_Integer] :=
  KroneckerProduct @@ ReplacePart[
    ConstantArray[SparseArray[IdentityMatrix[2]], l],
    j -> SparseArray[PauliMatrix[a]/2]];

XXZHamiltonian[l_Integer, j_?NumericQ, delta_?NumericQ, OptionsPattern[]] :=
 Module[{h = OptionValue["UniformField"], bl = OptionValue["LeftField"],
   br = OptionValue["RightField"], sx, sy, sz},
  If[l < 3 || !And @@ (TrueQ[Element[#, Reals]] & /@ {j, delta, h, bl, br}),
    Message[XXZHamiltonian::params]; Return[$Failed]];
  sx = Table[XXZSpin[l, 1, k], {k, l}];
  sy = Table[XXZSpin[l, 2, k], {k, l}];
  sz = Table[XXZSpin[l, 3, k], {k, l}];
  SparseArray[j Sum[sx[[k]].sx[[k+1]] + sy[[k]].sy[[k+1]] +
      delta sz[[k]].sz[[k+1]], {k, l-1}] + h Total[sz] + bl First[sz] + br Last[sz]]
 ];

(* ::Section:: *)
(* 2. Symmetry operators and sector construction *)
(* Reflection reverses site order. F=Product[PauliX[j]] flips every spin.
   Uniform nonzero fields preserve S^2 at Delta=1 but break SU(2) rotations
   of H: the good labels are (S,m), with Zeeman-split multiplets.
   Full SU(2) requires zero effective fields.
   At Delta=-1, U=Product[PauliZ[j],j even] gives
      U H_exchange(J,-1) U^dagger = H_exchange(-J,1).
   Thus U S_total^2 U^dagger, not ordinary S_total^2, is the extra charge.
   Reflection and F commute with this squared staggered spin as well.
*)
XXZZeroQ[a_, tol_] := Max[Abs[Flatten[Normal[N[a]]]]] < tol;
XXZCharges[l_Integer, delta_] := Module[{d = 2^l, r, f, sx, sy, sz, s2,
   stagger, result},
  r = SparseArray[Table[{1+FromDigits[Reverse[IntegerDigits[k,2,l]],2], k+1}->1,
      {k,0,d-1}], {d,d}];
  f = SparseArray[Table[{d-k,k+1}->1,{k,0,d-1}],{d,d}];
  result = <|"R"->r, "F"->f, "RF"->r.f|>;
  If[TrueQ[delta == 1] || TrueQ[delta == -1],
    stagger = If[TrueQ[delta == -1], Table[(-1)^(k-1),{k,l}], ConstantArray[1,l]];
    sx = Total[Table[stagger[[k]] XXZSpin[l,1,k],{k,l}]];
    sy = Total[Table[stagger[[k]] XXZSpin[l,2,k],{k,l}]];
    sz = Total[Table[XXZSpin[l,3,k],{k,l}]];
    s2 = SparseArray[sx.sx + sy.sy + sz.sz];
    AssociateTo[result, If[TrueQ[delta == 1],"S2","StaggeredS2"] -> s2]];
  result
 ];

(* A sector basis is stored as columns in its magnetization subspace.
   Diagonalize each commuting charge within the already-resolved subspaces.
   Orthonormalize within equal-charge groups before further splitting. *)
XXZSplitCharge[sector_, q_, name_, tol_] := Module[{b = sector["Basis"],
    vals, vecs, groups, sub, label},
  {vals,vecs} = Eigensystem[N[ConjugateTranspose[b].q.b]];
  groups = Gather[Range[Length[vals]], Abs[vals[[#1]]-vals[[#2]]] < tol &];
  Table[
    sub = Transpose[Orthogonalize[vecs[[g]]]];
    label = If[MemberQ[{"S2","StaggeredS2"},name],
      Round[(-1+Sqrt[1+4 Re[Mean[vals[[g]]]]])/2,1/2],
      Round[Re[Mean[vals[[g]]]]]];
    <|"Basis"->b.sub, "Labels"->Append[sector["Labels"],
        If[name=="S2","S",If[name=="StaggeredS2","StaggeredS",name]]->label]|>,
    {g,groups}]
 ];

Options[XXZSectorED] = Join[Options[XXZHamiltonian], {"Tolerance"->10^-9}];
XXZSectorED::check = "A sector residual exceeded tolerance: `1`. No result returned.";
XXZSectorED[l_Integer, j_?NumericQ, delta_?NumericQ, OptionsPattern[]] :=
 Module[{h, charges, tol=OptionValue["Tolerance"], d=2^l, out={}, ids, hm,
    m, candidates, sectors, q, es, vals, rows, ord, b, fullRows, residual,
    leakage, active, pars, scale, eigTol, fields, eligible},
  pars = <|"L"->l,"J"->j,"Delta"->delta,
    "UniformField"->OptionValue["UniformField"],
    "LeftField"->OptionValue["LeftField"],"RightField"->OptionValue["RightField"]|>;
  h = XXZHamiltonian[l,j,delta,"UniformField"->pars["UniformField"],
    "LeftField"->pars["LeftField"],"RightField"->pars["RightField"]];
  If[h===$Failed,Return[$Failed]];
  scale = Max[1, Norm[h,Infinity]]; eigTol = 100 tol scale;
  charges = XXZCharges[l,delta];
  (* Enforce parameter conditions before numerical verification: a tiny
     nonzero field is not silently interpreted as an exact symmetry. *)
  fields = ConstantArray[pars["UniformField"],l];
  fields[[1]] += pars["LeftField"]; fields[[-1]] += pars["RightField"];
  eligible = Select[Keys[charges], Switch[#,
    "R", And@@MapThread[TrueQ[#1==#2]&,{fields,Reverse[fields]}],
    "F", And@@(TrueQ[#==0]& /@ fields),
    "RF", And@@MapThread[TrueQ[#1 == -#2]&,{fields,Reverse[fields]}],
    "S2"|"StaggeredS2", And@@(TrueQ[#==First[fields]]& /@ fields),
    _,False]&];
  active = Select[eligible, XXZZeroQ[h.charges[#]-charges[#].h,tol scale]&];
  Do[
    m = l/2-nDown;
    ids = Select[Range[d],Total[IntegerDigits[#-1,2,l]]==nDown&];
    hm = N[Normal[h[[ids,ids]]]];
    candidates = Select[active, !MemberQ[{"F","RF"},#] || m==0&];
    (* S^2 first, then spatial/discrete parities. *)
    candidates = SortBy[candidates, Switch[#,"S2"|"StaggeredS2",0,"R",1,"F",2,_,3]&];
    sectors = {<|"Basis"->IdentityMatrix[Length[ids]],"Labels"-><|"m"->m|>|>};
    Do[
      q = N[Normal[charges[name][[ids,ids]]]];
      sectors = Flatten[XXZSplitCharge[#,q,name,tol]& /@ sectors,1],
      {name,candidates}];
    Do[
      b = s["Basis"];
      leakage = Max[Abs[Flatten[hm.b-b.(ConjugateTranspose[b].hm.b)]]];
      {vals,rows} = Eigensystem[N[ConjugateTranspose[b].hm.b]];
      ord = Ordering[Re[vals]]; vals = Re[vals[[ord]]];
      (* Normalize/orthogonalize only within degenerate ENERGY eigenspaces. *)
      rows = rows[[ord]];
      Do[rows[[g]] = Orthogonalize[rows[[g]]],
        {g,Gather[Range[Length[vals]],Abs[vals[[#1]]-vals[[#2]]] < tol scale&]}];
      rows = Transpose[b.Transpose[rows]];
      residual = Max[MapThread[Norm[hm.#2-#1 #2]&,{vals,rows}]];
      If[Max[residual,leakage]>eigTol,
        Message[XXZSectorED::check,{s["Labels"],residual,leakage}];Return[$Failed]];
      fullRows = Table[Normal[SparseArray[Thread[ids->v],{d}]],{v,rows}];
      AppendTo[out,<|"Labels"->s["Labels"],"Dimension"->Length[vals],
        "BasisIndices"->ids,"SectorBasis"->b,"Eigenvalues"->Chop[vals],
        "Eigenvectors"->fullRows,"Residual"->residual,"Leakage"->leakage|>],
      {s,sectors}],
    {nDown,0,l}];
  <|"Parameters"->pars,"Hamiltonian"->h,"ActiveCharges"->active,
    "Sectors"->out,"Eigenvalues"->Sort[Flatten[Lookup[out,"Eigenvalues"]]],
    "DimensionSum"->Total[Lookup[out,"Dimension"]]|>
 ];

(* ::Section:: *)
(* 3. Edge reduced state: exact trace of the interior *)
(* Input is a normalized pure state in the usual site order {1,...,L}.
   Output order is {|up up>,|up down>,|down up>,|down down>} on {1,L}.
   The two edges need not be independent: their environment is shared. *)
XXZEdgeDensity[psi_List,l_Integer] := Module[{a},
  If[Length[psi]!=2^l,Return[$Failed]];
  a = ArrayReshape[Transpose[ArrayReshape[psi,{2,2^(l-2),2}],{1,3,2}],
    {4,2^(l-2)}];
  Chop[a.ConjugateTranspose[a]]
 ];

(* ::Section:: *)
(* 4. Independent verification of the complete ED *)
XXZValidate[result_Association,tol_:10^-8] := Module[{h=result["Hamiltonian"],
    v, evals, d, checks},
  d = Length[h];
  v = Flatten[Lookup[result["Sectors"],"Eigenvectors"],1];
  evals = Flatten[Lookup[result["Sectors"],"Eigenvalues"]];
  checks = <|"DimensionSum"->result["DimensionSum"],
    "ExpectedDimension"->d,
    "OrthonormalityError"->Max[Abs[Flatten[v.ConjugateTranspose[v]-IdentityMatrix[d]]]],
    "EigenpairError"->Max[MapThread[Norm[h.#2-#1 #2]&,{evals,v}]],
    "TraceError"->Abs[Total[evals]-Tr[h]],
    "TraceH2Error"->Abs[Total[evals^2]-Tr[h.h]]|>;
  (* Independent unsplit diagonalization only for small chains. *)
  If[d<=256,AssociateTo[checks,"FullSpectrumError"->
    Max[Abs[Sort[Re[Eigenvalues[N[Normal[h]]]]]-Sort[evals]]]]];
  AssociateTo[checks,"Passed"->(result["DimensionSum"]==d &&
    checks["OrthonormalityError"]<tol &&
    checks["EigenpairError"]<tol Max[1,Norm[h,Infinity]] &&
    checks["TraceError"]<tol Max[1,Norm[h,Infinity]] d &&
    checks["TraceH2Error"]<tol Max[1,Norm[h,Infinity]]^2 d &&
    Lookup[checks,"FullSpectrumError",0]<tol Max[1,Norm[h,Infinity]])];
  checks
 ];

(* ::Section:: *)
(* 5. Run/export function: all parameters appear on every plot *)
(* Call XXZRun independently for each point in any parameter sweep.
   No moving average, unfolding, or processing of the spectrum.
   A fresh numbered folder prevents overwriting without random names.
   WXF retains every sector basis, eigenvalue, eigenvector, and edge state.
   Full ED storage and validation scale exponentially; begin with L=5.
*)
Options[XXZRun] = Join[Options[XXZSectorED],{"OutputDirectory"->Automatic}];
XXZRun[l_Integer,j_?NumericQ,delta_?NumericQ,OptionsPattern[]] := Module[
  {res,seconds,check,root,base,dir,n=1,table,plot,data,label,edges},
  {seconds,res}=AbsoluteTiming[XXZSectorED[l,j,delta,
    "UniformField"->OptionValue["UniformField"],
    "LeftField"->OptionValue["LeftField"],"RightField"->OptionValue["RightField"],
    "Tolerance"->OptionValue["Tolerance"]]];
  If[res===$Failed,Return[$Failed]];
  check=XXZValidate[res]; Print["Validation: ",check];
  If[!TrueQ[check["Passed"]],Return[$Failed]];
  table=MapIndexed[<|"Sector"->First[#2],"Labels"->#1["Labels"],
    "Dimension"->#1["Dimension"],"Emin"->Min[#1["Eigenvalues"]],
    "Emax"->Max[#1["Eigenvalues"]]|>&,res["Sectors"]];
  Print[Dataset[table]];
  root=OptionValue["OutputDirectory"];
  If[root===Automatic,
    root=If[$FrontEnd===Null,Directory[],Quiet[Check[NotebookDirectory[],Directory[]]]];
    If[!StringQ[root],root=Directory[]];
    root=FileNameJoin[{root,"xxz_open_edges_results"}]];
  If[!StringQ[root],root=FileNameJoin[{Directory[],"xxz_open_edges_results"}]];
  If[!DirectoryQ[root],CreateDirectory[root,CreateIntermediateDirectories->True]];
  base="L"<>ToString[l]<>If[TrueQ[delta==1],"_isotropic","_anisotropic"];
  dir=FileNameJoin[{root,base<>"_run"<>ToString[n]}];
  While[DirectoryQ[dir],n++;dir=FileNameJoin[{root,base<>"_run"<>ToString[n]}]];
  CreateDirectory[dir];
  label="Open XXZ; S={1,L}; "<>StringRiffle[
    KeyValueMap[ToString[#1]<>"="<>ToString[#2,InputForm]&,res["Parameters"]],"; "];
  data=Table[({k,#}& /@ res["Sectors"][[k]]["Eigenvalues"]),
    {k,Length[res["Sectors"]]}];
  plot=ListPlot[data,Joined->False,PlotRange->All,Frame->True,
    FrameLabel->{"Symmetry sector index","Energy E"},PlotLabel->label,
    PlotStyle->Black,ImageSize->1000,
    FrameTicks->{{Automatic,None},{Range[Length[data]],None}}];
  Export[FileNameJoin[{dir,"spectrum_all_sectors.png"}],plot];
  edges=Map[XXZEdgeDensity[#,l]& /@ #["Eigenvectors"]&,res["Sectors"]];
  AssociateTo[res,<|"EdgeReducedStates"->edges,"Validation"->check,
    "SectorTable"->table,"DiagonalizationSeconds"->seconds|>];
  Export[FileNameJoin[{dir,"complete_sector_ED.wxf"}],res];
  Export[FileNameJoin[{dir,"parameters_and_sectors.txt"}],
    label<>"\nSpin convention S=Pauli/2; boundary=open; no attached probes.\n"<>
    "Computational site order={1,...,L}; edge order={1,L}.\n"<>
    "DiagonalizationSeconds="<>ToString[seconds,InputForm]<>"\n"<>
    "ActiveCharges="<>ToString[res["ActiveCharges"],InputForm]<>"\n"<>
    "Validation="<>ToString[check,InputForm]<>"\n"<>
    StringRiffle[ToString[#,InputForm]& /@ table,"\n"],"Text"];
  Print["Saved: ",dir]; res
 ];

(* ::Section:: *)
(* 6. Delta != 1: full ED in each surviving symmetry sector *)
(* Edit these free inputs. J is not replaced by Abs[J].
   deltaAnisotropic can be negative, zero, or positive, but cannot be 1.
   Delta=-1 is treated with its extra staggered total-spin charge.
   At zero fields generic sectors are (m,R) and additionally F in m=0.
*)
chainLength=5;
exchangeJ=1;
deltaAnisotropic=1/2;
uniformField=0;
leftBoundaryField=0;
rightBoundaryField=0;

If[TrueQ[deltaAnisotropic==1],
  Print["Section 6 requires Delta != 1."],
  anisotropicResult=XXZRun[chainLength,exchangeJ,deltaAnisotropic,
    "UniformField"->uniformField,"LeftField"->leftBoundaryField,
    "RightField"->rightBoundaryField]];

(* ::Section:: *)
(* 7. Delta = 1: full ED in each surviving symmetry sector *)
(* At zero fields: (S,m,R); F at m=0 is redundant with S for ordinary SU(2)
   but is explicitly labeled when present. Energies repeat for m=-S..S.
   A uniform field keeps (S,m,R), shifts E by h*m, and splits multiplets.
   Boundary field differences can destroy S^2 and/or reflection: then retain
   only the charges that actually commute. Delta=1 alone does not ensure SU(2).
   Set isotropicJ=-1 to examine the ferromagnetic isotropic chain.
*)
isotropicLength=5;
isotropicJ=1;
isotropicUniformField=0;
isotropicLeftField=0;
isotropicRightField=0;

isotropicResult=XXZRun[isotropicLength,isotropicJ,1,
  "UniformField"->isotropicUniformField,"LeftField"->isotropicLeftField,
  "RightField"->isotropicRightField];

(* ::Section:: *)
(* 8. Access results without re-diagonalizing *)
(* anisotropicResult["Sectors"][[k]]["Labels"]
   anisotropicResult["Sectors"][[k]]["Eigenvalues"]
   anisotropicResult["Sectors"][[k]]["Eigenvectors"][[r]]
   anisotropicResult["EdgeReducedStates"][[k,r]]
   isotropicResult["SectorTable"] // Dataset
   loaded=Import[".../complete_sector_ED.wxf"];

   For a reproducible sign comparison, run XXZSectorED[L,+Abs[J],Delta,...]
   and XXZSectorED[L,-Abs[J],Delta,...] with the same stated fields.
   At zero fields they have identical eigenspaces with reversed energies;
   their ground states need not coincide.
*)
