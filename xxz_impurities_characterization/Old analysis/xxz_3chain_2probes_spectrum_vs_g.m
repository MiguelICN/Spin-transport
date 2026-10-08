(* Five-spin order: {A,1,2,3,B}.  S = sigma/2. *)
ClearAll["Global`*"];
Nspin = 5; dim = 2^Nspin; J = 1.; delta = 0.5; omegaA = 0.6; omegaB = 0.9;
gList = Range[0., 0.5, 0.005];

id = IdentityMatrix[2]; sx0 = PauliMatrix[1]/2; sy0 = PauliMatrix[2]/2; sz0 = PauliMatrix[3]/2;
op[a_, j_] := KroneckerProduct @@ ReplacePart[ConstantArray[id, Nspin], j -> a];
sx[j_] := op[sx0, j]; sy[j_] := op[sy0, j]; sz[j_] := op[sz0, j];
sp[j_] := sx[j] + I sy[j]; sm[j_] := sx[j] - I sy[j];

Hchain = J (sx[2].sx[3] + sy[2].sy[3] + delta sz[2].sz[3] +
             sx[3].sx[4] + sy[3].sy[4] + delta sz[3].sz[4]);
HZ = omegaA sz[1] + omegaB sz[5];
Vzz = sz[1].sz[2] + sz[4].sz[5];
Vex = sp[1].sm[2] + sm[1].sp[2] + sp[4].sm[5] + sm[4].sp[5];

(* 1. ZZ coupling, no Zeeman splitting *)
ezz0 = Table[Sort[Chop[Eigenvalues[N[Hchain + g Vzz]]]], {g, gList}];
p1 = ListPlot[Flatten[Table[{gList[[i]], ezz0[[i, j]]}, {i, Length[gList]}, {j, dim}], 1],
   PlotRange -> All, AxesLabel -> {"g", "E"}, PlotLabel -> "ZZ, no Zeeman"];

(* 2. ZZ coupling, with probe Zeeman splitting *)
ezzZ = Table[Sort[Chop[Eigenvalues[N[Hchain + HZ + g Vzz]]]], {g, gList}];
p2 = ListPlot[Flatten[Table[{gList[[i]], ezzZ[[i, j]]}, {i, Length[gList]}, {j, dim}], 1],
   PlotRange -> All, AxesLabel -> {"g", "E"}, PlotLabel -> "ZZ, with Zeeman"];

(* 3. Exchange coupling, no Zeeman splitting *)
eex0 = Table[Sort[Chop[Eigenvalues[N[Hchain + g Vex]]]], {g, gList}];
p3 = ListPlot[Flatten[Table[{gList[[i]], eex0[[i, j]]}, {i, Length[gList]}, {j, dim}], 1],
   PlotRange -> All, AxesLabel -> {"g", "E"}, PlotLabel -> "Exchange, no Zeeman"];

(* 4. Exchange coupling, with probe Zeeman splitting *)
eexZ = Table[Sort[Chop[Eigenvalues[N[Hchain + HZ + g Vex]]]], {g, gList}];
p4 = ListPlot[Flatten[Table[{gList[[i]], eexZ[[i, j]]}, {i, Length[gList]}, {j, dim}], 1],
   PlotRange -> All, AxesLabel -> {"g", "E"}, PlotLabel -> "Exchange, with Zeeman"];

GraphicsGrid[{{p1, p2}, {p3, p4}}, ImageSize -> Large]
