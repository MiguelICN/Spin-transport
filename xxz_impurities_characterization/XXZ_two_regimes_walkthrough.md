# Walkthrough: independent ZZ and exchange edge implementations

Updated 26 October 2026: added ZZ section 1B (isolated-bulk susceptibility), exact eigenvector input and optional controlled dynamics comparison. The exchange file is unchanged by this revision. Previous update, 16 October 2026: incorporates the independent assessment and preserves the subsequent accelerated chapter-3 channel implementation. The audit examined the earlier QST revision; its statements about one-time-at-a-time channel evolution do not describe the current chapter 3.

Formatting correction, 26 October 2026: the ZZ `.m` is assembled directly from the verified uploaded package. Standalone explanatory comments use `Text` cell markers rather than `Input` markers. Section 1B is ordered as heading, definitions (1B.1), computation (1B.2), selection (1B.3), and plots/export (1B.4). The optional comparison 2.8 follows the chapter-2 export and precedes chapter 3. Evaluate the opening context cell, shared definitions, and 1B.1 before 1B.2; a missing-definitions guard now stops with that instruction rather than producing cascading `Lookup` errors. A real `Get`-based definition load in the `XXZZZ` context was tested, followed by a six-state susceptibility calculation. Front-end opening/rendering on Windows has not been reproduced here.

This document describes `XXZ_ZZ_independent.m` and `XXZ_Exchange_independent.m`. The numerical conventions and reusable routines come from the supplied `xxz_open_edges_revised_benchmark.nb`. Chapter 3 now follows the supplied revision-2 `edge_map_implementation_guide(1).md` and accompanying channel properties, adapted to the existing Hamiltonian and length conventions. The previous channel helpers are retained for the independent QST chapter.

Each self-contained Wolfram `.m` file contains executable input cells, section markers, its own definitions and four chapters:

1. Static eigenpairs, eigenstate entropies, computational-basis IPR, and optional sector DOS/gap diagnostics. ZZ additionally has section 1B: isolated-bulk boundary response and eigenvector selection.
2. Exact unitary evolution of a chosen global initial state and the reduced edge density matrix.
3. The edge dynamical map for a fixed bulk preparation, including its Choi state and channel diagnostics.
4. Quantum-state texture (QST): independent edge-state dynamics and channel texturizing-power calculations.

The `.m` files replace the previously delivered notebook format for this revision. They contain ordinary executable Wolfram code, not commented-out input cells. The two static scans are independent. Chapters 2 and 3 require the notebook's initialization definitions, but no static computation. Chapter 3 also requires no evaluation of chapter 2. The files begin in separate named Wolfram contexts (`XXZZZ` and `XXZExchange`). Evaluate the Begin cell before evaluating their definitions and selected run cells. Evaluating the entire file with `Get` executes all run cells; it is not a definitions-only loader. Neither notebook loads the other or loads the benchmark at runtime.

The opening cell now saves the caller's `$ContextPath` in an explicitly qualified symbol and sets the search path to `{"System`"}` after entering the file context. Existing `Global` definitions therefore cannot capture unqualified helper names. The final cell restores the saved path and calls `End[]`. Evaluate the opening cell separately before definitions/run cells and the closing cell when leaving the file, including after an aborted interactive run. No kernel reset or clearing of the user's `Global` symbols is required. Loading through `Get` reads the separate source expressions in sequence; passing the whole source through a single pre-parsed `ToExpression` is not the supported loading method.


## 1. Physical models and conventions

There are **L physical spin-1/2 sites**, with edges $A=1$, $B=L$, and bulk $M=\{2,\ldots,L-1\}$. The bulk contains $L-2$ spins. Site 1 is the leftmost tensor factor. Bit 0 means up and bit 1 means down. The edge basis is $00,01,10,11$.

The Hamiltonian convention is $S^a=\sigma^a/2$, with $\hbar=1$. All entropies use natural logarithms. $J$ and $g$ have energy units; $g$ is not $g/J$. Weak coupling means $|g|\ll |J|$.

Both notebooks use the same field-free bulk:

$$
H_M=J\sum_{j=2}^{L-2}
\left(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y+\Delta S_j^zS_{j+1}^z\right).
$$

The two end bonds are:

| Notebook | Edge interaction |
|---|---|
| ZZ | $g(S_1^zS_2^z+S_{L-1}^zS_L^z)$ |
| Exchange | $g(S_1^xS_2^x+S_1^yS_2^y+S_{L-1}^xS_L^x+S_{L-1}^yS_L^y)$ |

Thus $H=H_M+H_{\rm edge}$. There are no additional magnetic fields on any site, in any direction. There is no ZZ term on an exchange end bond, and no XY term on a ZZ end bond. The exchange matrix element for $|\uparrow\downarrow\rangle\leftrightarrow|\downarrow\uparrow\rangle$ is $g/2$.

Each run has its own association containing only `L`, `J`, `Delta`, and `g`. The validator accepts integer $L\ge4$, finite real parameters, and nonzero $J,g$. Positive and negative $g$ are allowed. An exactly zero $\Delta$ is allowed. An extra model key, including a field parameter, is rejected rather than silently ignored.

### 1.1 ZZ symmetry sectors

For every allowed $g,\Delta$, both edge z components and the bulk magnetization are conserved:

$$
s_A=S_1^z,\qquad s_B=S_L^z,\qquad
m_M=\sum_{j=2}^{L-1}S_j^z.
$$

The code retains $m=s_A+s_B+m_M$ as a redundant, useful label. A fixed-edge computational block has dimension

$$
\binom{L-2}{(L-2)/2-m_M}.
$$

Reflection $R$ and spin flip $F=\prod_j\sigma_j^x$ are symmetries of the full Hamiltonian, but they need not preserve an individual fixed-edge block. The compatible labels are:

| Fixed-edge block | Spatial label resolved inside it |
|---|---|
| $s_A=s_B$ | Reflection $R=\pm1$ |
| $s_A=-s_B$, $m_M=0$ | Combined $RF=\pm1$ |
| $s_A=-s_B$, $m_M\ne0$ | Neither of those operators preserves the block individually |

Reflection and spin flip can still relate distinct partner blocks. Those blocks remain separately represented and diagonalized; matching spectra do not cause their eigenstates or statistics to be pooled.

Conditioning on the edge values gives the exact bulk operator

$$
H_{s_As_B}=H_M+g s_A S_2^z+g s_B S_{L-1}^z.
$$

These conditional boundary terms arise from the two-spin interactions. They are not externally added fields in the original $L$-spin Hamiltonian.

**Every fixed-edge ZZ block is an integrable open XXZ chain with diagonal boundary fields, for every allowed $g$ and $\Delta$.** The full model is a direct sum of these integrable blocks. Its gap ratios are an integrable control, not a Poisson-to-GOE crossover probe. Integrability does not force a finite, partly symmetry-resolved sample to have the exact Poisson mean: additional charges, exact degeneracies and finite-size effects remain relevant. The exchange model is a different boundary interaction and can display a crossover; its finite-size statistics still require checks.

The integrable diagonal-boundary construction follows E. K. Sklyanin, *Boundary conditions for integrable quantum systems*, J. Phys. A **21**, 2375 (1988), [doi:10.1088/0305-4470/21/10/015](https://doi.org/10.1088/0305-4470/21/10/015). General open XXZ integrable boundaries are also treated in [arXiv:1401.4901](https://arxiv.org/abs/1401.4901). These references establish the integrable construction; the audit's finite-size gap-ratio values are separate numerical observations, not a universal theorem about each finite block.

The implementation checks the actual conditional-bulk charges and refines a block only when the charge preserves that block and commutes with its Hamiltonian. This retains the benchmark's small-system exception: for a two-site bulk, equal conditional boundary fields are uniform, so the ordinary or staggered bulk Casimir can supply an extra label at $\Delta=+1$ or $-1$, respectively. Such child labels remain distinct in the results. For `StaggeredSBulk`, the Casimir is $(\widetilde S_M^x)^2+(\widetilde S_M^y)^2+(S_M^z)^2$, with $\widetilde S_M^{x,y}=\sum_{j=2}^{L-1}(-1)^j S_j^{x,y}$ and **uniform** $S_M^z=\sum_{j=2}^{L-1}S_j^z$. Reversing the overall staggering sign leaves the Casimir unchanged.

### 1.2 Exchange symmetry sectors

The exchange interaction conserves total $m$, but it can change either edge spin and the bulk magnetization separately. The compatible labels are:

- $m\ne0$: $m,R$.
- $m=0$: $m,R,F$.

Spin flip maps $m$ to $-m$; it is not an independent label within a nonzero-$m$ block. For odd $L$, magnetizations are half-integers, so there is no $m=0$ block. Every permitted sector, including both signs of $m$, is included in a static eigensystem.

The XY-only end bonds prevent full-chain SU(2) symmetry from being restored simply by choosing $\Delta=\pm1$. Homogeneous-chain total-spin labels are therefore not imposed on this model.

### 1.3 What changes at special anisotropies

At $\Delta=0$, the exchange Hamiltonian is free fermionic; each conditional ZZ bulk Hamiltonian is also free fermionic. The listed spin and spatial sectors remain valid and complete as a decomposition of the Hilbert space. Additional mode occupations are conserved. **The notebook resolves the benchmark's compatible spin/spatial charges; it does not claim to diagonalize every integrability charge or assign every fermion occupation as a separate sector label.** Consequently a listed sector can still contain exact degeneracies.

No energy-degenerate vectors are mixed across the resolved sectors. Inside a remaining degenerate eigenspace, the eigensolver chooses a basis. Individual eigenstate entropies and IPR can depend on that basis, especially at special parameter values. The spectrum, complete unitary evolution, and correctly clustered time-averaged map do not depend on this choice.

At $\Delta=\pm1/2$, the notebook identifies the previously observed special bulk-degeneracy control. A sub-Poisson weak-coupling gap ratio at this point is not interpreted as a generic Poisson-to-GOE crossover. The default comparison anisotropy is $9/20=0.45$.

## 2. Chapter 1: static calculations

### 2.1 Independent scans in the two directions

Sections 1.2–1.5 compute a $g$ scan at fixed $\Delta$, with `statGModel`, `statGGrid`, and `statGRun`.

Sections 1.6–1.9 independently compute a $\Delta$ scan at fixed $g$, with `statDeltaModel`, `statDeltaGrid`, and `statDeltaRun`.

Each point is a separate static Hamiltonian. No state is dynamically transported between points. An eigenvalue's index means its energy rank within its sector at that point; it is not an overlap-tracked eigenstate identity through a crossing.

The initial static defaults are $L=8,J=1$. The $g$ scan uses 31 samples from 0.01 to 1 at $\Delta=0.45$; the anisotropy scan uses 31 samples from $-1$ to 2 at $g=0.1$, including zero anisotropy. These are editable demonstration sizes and ranges, not imposed physical restrictions.

### 2.2 Numerical diagonalization

`regimeHamiltonian` assembles the model from sparse embedded spin operators. `regimeSetup` constructs the compatible subspaces, reusing cached parent symmetry bases for a given $L$ and tolerance. Parameter-dependent ZZ refinements are checked for the current model.

The setup audits computational-index coverage, symmetry-basis orthonormality, and invariance under the bulk and interaction parts. Thus an invalid fixed-edge partition cannot silently drop a U(1)-conserving edge-flip term. The real-sector guard also checks Hermiticity, reality, and total-$S^z$ conservation.

For each sector basis $B_s$, the calculation forms

$$
H_s=B_s^\dagger H B_s
$$

within its computational-index block, then performs dense `Eigensystem`. The benchmark's `diagonalizeBlock` sorts eigenvalues, reconstructs eigenvectors in the computational coordinates of that block, fixes their real overall sign, and checks eigenvector residuals, orthogonality, and subspace leakage. All sector dimensions must sum to $2^L$.

These are machine-precision eigensystems, as in the benchmark. Small nonzero energy differences are retained: the diagonalization does not `Chop`, deduplicate, or unfold levels. Increasing the precision of a parameter alone does not turn this implementation into an arbitrary-precision solver.

### 2.3 Quantities attached to every eigenvector

For each normalized eigenvector $|n\rangle$, the retained quantities are:

| Quantity | Definition and meaning |
|---|---|
| `Eigenvalues` | $E_n$, in the Hamiltonian's energy units |
| `Vectors` | The eigenvectors, stored with the sector's full-chain computational `Indices` |
| `HalfEntropy` | Entanglement of sites $1,\ldots,\lfloor L/2\rfloor$ with the remaining sites |
| `EdgesEntropy` | Entanglement of the joint edge pair $AB$ with the bulk $M$ |
| `IPR` | $\sum_z|\langle z|n\rangle|^4$ in the physical computational basis |

The half-chain entropy uses squared singular values of the reshaped eigenket. The edge/bulk entropy uses the eigenvalues of the $4\times4$ reduced edge state. Both calculate $-\sum p\log p$.

The IPR is bounded by $2^{-L}\le\mathrm{IPR}\le1$. It measures concentration in the computational basis, not in the symmetry-adapted basis. A large IPR by itself is not a localization or chaos diagnosis.

For the fixed-edge eigenbasis in the ZZ notebook,

$$
|n\rangle=|s_A\rangle\otimes|\phi_n\rangle_M\otimes|s_B\rangle,
$$

so `EdgesEntropy` vanishes. The code checks this. The half-chain entropy and IPR can nevertheless vary. Superpositions of degenerate eigenstates from different edge blocks can have nonzero edge/bulk entanglement, but they are not the basis selected for the static ZZ eigensystem.

Both scans retain all eigenvectors. `eigenket[run["Points"][[k]],s,n]` reconstructs the full ket for point $k$, sector $s$, and level $n$. `sectorTable` reports the actual labels, dimensions, energy range, smallest positive gap, and residual.

Every static panel includes the resolved sector labels and its actual dimension. If a label occurs with different dimensions across the scan, the displayed dimension set records that fact. The static plot settings also allow selecting specific `SectorKeys`; all sectors and eigenpairs remain computed and retained.

The current ZZ static plots overlay symmetry-related sectors using symmetry-paired eigenvectors. Each member retains its own labels and dimension; odd-chain antiparallel half-chain entropies retain their two spin-flip families. Export counts depend on groups and parameter-dependent refinements. Compact quantity-first filenames are mapped back to full plot keys in `ExportIndex_all_sectors.csv`. Raw spectra remain sector-resolved. Caption formatting uses rounded sampling spacings; it does not alter sampling times.

### 2.4 DOS before adjacent-gap statistics

Sections 1.10–1.11 provide optional postprocessing of either completed static scan. `diagSector` must identify exactly one resolved sector at every selected point; the code refuses to merge multiple refinements.

The DOS uses the entire raw sector spectrum. Its histogram is normalized to unit integral, and the empirical CDF counts every level, including degeneracies. A rank window $\{a,b\}$ selects sorted levels

$$
n_{\min}=\max(1,\lfloor aD\rfloor+1),\qquad
n_{\max}=\min(D,\lceil bD\rceil).
$$

The DOS record retains the ranks and their actual energy boundaries. A central rank window selects a fraction of the levels, not a fixed fraction of the energy bandwidth. Changing the rank window requires refreshing the DOS diagnostic before computing the corresponding means; no ED rerun is needed.

For consecutive spacings $\delta_n=E_{n+1}-E_n$, the ratio is

$$
\widetilde r_n=\frac{\min(\delta_n,\delta_{n+1})}{\max(\delta_n,\delta_{n+1})}.
$$

A gap is unresolved when it is at or below `GapResolution` times $\max(1,\max_n|E_n|)$. If any adjacent-gap pair in the selected window is unresolved, the mean is `Missing`; those levels are never deleted to manufacture a resolved sample. The record includes individual ratios, unresolved counts, mean spacing, the mean ratio, the fraction below 0.05, and two uncertainty diagnostics.

`IIDSEM` treats ratios as independent. `BlockSEM` takes the standard error of means of complete, nonoverlapping consecutive blocks of `BlockLength` ratios. The final incomplete block is omitted only from that error estimate; the reported mean still uses every resolved ratio. These are descriptive finite-spectrum estimates, not calibrated confidence intervals. No automatic crossover is inferred.

The $L=8$ default selected sectors have dimensions 10 (ZZ) and 23 (exchange), yielding only 4 and 11 central-window ratios. The revised diagnostics therefore withhold `MeanR` and report `InsufficientStatistics` when the sector has fewer than `diagMinimumDimension=500` states or the selected window fewer than `diagMinimumRatios=100` **numeric, resolved** ratios. `ResolvedRatioCount` and `RatioEntryCount` distinguish usable ratios from all adjacent-gap-pair entries. Unresolved entries remain `Missing`; they no longer count toward the sample threshold. `FractionRBelow005` is present in the insufficient-statistics branch as `Missing["InsufficientStatistics"]`. This does not override the existing rule that unresolved gaps make a reported mean unreliable: the underlying gap diagnostic still withholds that mean when its resolution check fails. DOS is still available. Lowering the thresholds is an explicit demonstration choice; meeting them is a minimum sample-size check, not a guarantee of asymptotic statistics or a calibrated confidence interval. The static size remains small by default to avoid launching a large all-sector scan. For quantitative diagnostics, increase the size/select a larger sector and inspect its actual dimension before interpreting the result. In practice the corresponding default symmetry families need about $L=16$ for ZZ or $L=14$ for exchange to pass the 500-state requirement; compatible refinements and the chosen magnetization can change these counts.

## 2B. ZZ section 1B: isolated-bulk boundary susceptibility

This section sits between static characterization and reduced-state dynamics. It uses initialization definitions only: neither static scan has to run. It was added only to the attached ZZ implementation; the exchange implementation has not been changed by this revision.

The proposed target is **edge–edge entanglement generation**, measured by concurrence and negativity. Susceptibility measures eigenvector sensitivity and is a candidate predictor to test, not an entanglement witness or a demonstrated monotonic predictor. Edge–bulk entanglement is a different quantity: for pure global preparations it equals the reduced edge entropy, already available as `EntropyAB`.

### Operators and normalization

The unperturbed Hamiltonian is the field-free isolated bulk $H_0$, with $N=L-2$ spins. Define $Z_L=S_2^z$ and $Z_R=S_{L-1}^z$, each with eigenvalues $\pm1/2$. The conditional Hamiltonians are

$$H_a=H_0+g(s_A Z_L+s_B Z_R),\qquad s_A,s_B=\pm\tfrac12.$$

For an isolated eigenstate $|n\rangle$, the code computes

$$\chi_{ij}^{n}=\operatorname{Re}\sum_{m\ne n}
\frac{\langle n|Z_i|m\rangle\langle m|Z_j|n\rangle}{(E_m-E_n)^2},\qquad
\chi_n^{(\pm)}=\sum_{m\ne n}\frac{|\langle m|(Z_L\pm Z_R)|n\rangle|^2}{(E_m-E_n)^2}.$$

Thus $\chi^{(\pm)}=\chi_{LL}+\chi_{RR}\pm2\chi_{LR}$, with units of inverse energy squared. For parallel/antiparallel conditional fields, $1-F_n^{(\pm)}=g^2\chi_n^{(\pm)}/4+O(g^3)$. This factor of four follows from the spin-1/2 convention; the code does not replace $S^z$ by $\sigma^z$.

### Complete transition space and eigenstate identifiers

`bulkEigenCatalog[model,mu,tolerance]` diagonalizes the entire chosen **bulk magnetization** space. It first resolves compatible reflection, spin-flip and, at special anisotropies, spin-Casimir labels using the existing charge routines. It records the actual labels, sector dimension, energy-order `Level`, catalog `ID`, and computational-basis vector for every eigenstate. Dense diagonalization is confined to resolved sectors; bulk operators are initially sparse. This remains exact diagonalization with exponential Hilbert-space cost, not a large-system tensor-network method.

Response sums include every target eigenstate in that magnetization space, **including different parity sectors**. $Z_L+Z_R$ preserves reflection; $Z_L-Z_R$ reverses it. At zero magnetization both reverse spin-flip parity. Restricting the sum to the source's complete parity label would incorrectly remove allowed transitions. Other bulk magnetizations do not contribute because the boundary operators conserve magnetization.

`boundaryResponse` retains the matrix elements of both boundary operators, the full real $2\times2$ fidelity metric, both scalar susceptibilities, and the variances of $Z_L\pm Z_R$. The variances are zeroth spectral moments; they distinguish instantaneous quench sensitivity from the inverse-frequency-squared weighting of susceptibility.

`boundaryPeaks` returns the discrete distribution as a table of frequencies $E_m-E_n$, target IDs/labels and weights $|\langle m|V_\pm|n\rangle|^2$. Both positive and negative frequencies are retained. The plotted weights are delta-peak weights, not a binned DOS or broadened density. No additional normalization or missing factor of $2\pi$ is implicit.

### Degeneracies and numerical resolution

An eigenstate with any other energy within `susGapTolerance Max[1,Max[Abs[energies]]]` in the same bulk magnetization space receives `Missing["DegenerateOrUnresolved"]` for the tensor and both susceptibilities. Its partner IDs remain visible. Such states are excluded from susceptibility ranking and the nondegenerate fidelity check. No small denominator is replaced by a finite cutoff, and no degenerate level is silently deleted from a reported susceptibility.

This is deliberately conservative: some degeneracies cannot mix under a particular perturbation, but this revision does not implement perturbation-dependent degenerate-subspace susceptibility. Cross-magnetization degeneracies do not invalidate the calculation because the perturbations cannot connect them. A degenerate catalog eigenvector can still be used as a precisely saved initial state for exact dynamics; its eigenvector statistics remain basis-dependent.

### Running the sensitivity experiment

1. Evaluate the shared definitions and section 1B.1 definitions.
2. Set `susModel` (`L`, `J`, `Delta`), `susMu`, ED and gap tolerances. There is no edge coupling in this catalog.
3. Compute the catalog and response. Choose `susRankKey="ChiPlus"` or `"ChiMinus"`; optionally restrict `susEnergyWindow={Emin,Emax}` to compare eigenstates at similar energies. The ranked ID list runs from least to most sensitive.
4. Set `susSelectedID` to any displayed ID. `susSelectedBulk` stores the actual full bulk ket, source model, sector labels and level. IDs are local to that catalog, not labels that track an eigenstate across parameter changes.
5. Inspect its transition peaks and `susFidelity`. The fidelity check diagonalizes $H_0+(g/2)V_\pm$ in the complete fixed-magnetization eigenbasis and takes the eigenvector with largest overlap with the original state. `MatchedOverlap` exposes branch ambiguity; test decreasing `susCheckG` and do not interpret a low-overlap match as controlled perturbation theory.

Figures show each scalar susceptibility against energy, the selected state's transition weights, and actual versus perturbative infidelity. Raw vectors, labels, matrix elements, rankings, selection and finite-field checks are retained in the numerical export.

### Choosing the exact bulk eigenvector in chapter 2

Chapter 2 remains independent of section 1B. Its `dynBulkChoice` offers:

- `"Neel"`: preserves the supplied file's original default preparation.
- `"Eigenstate"`: set `dynBulkMu`, `dynBulkSector` and `dynBulkLevel`. This constructs the isolated-bulk catalog for the chapter-2 model and requires exactly one match. Include compatible extra spin labels at $\Delta=\pm1$ if needed.
- `"SelectedEigenvector"`: uses `susSelectedBulk`, including its stored vector. Set chapter-2 `L,J,Delta` to match the selected catalog; `g` is independently adjustable. For example, `dynModel=Join[susSelectedBulk["Model"],<|"g"->0.1|>]` before assembling the initial factor.

These correspond to new `bulkFactor` types `"Eigenstate"` and `"Eigenvector"`. The latter checks source-model agreement, vector normalization and the eigenvalue residual against the current $H_0$. The channel and QST bulk constructors also inherit these types through their existing fallback to `bulkFactor`; no channel algorithm has been changed.

Use `dynEdgeState=ConstantArray[1/2,4]` for the separable $|+\rangle_A|+\rangle_B$ input when testing entanglement **generation**. The supplied default Bell edge state is retained for continuity, but its concurrence already starts at one and measures preservation/decay instead. Keep edge input, model, time grid and energy window controlled when comparing bulk eigenstates.

The assembled full initial factor is projected into the **coupled model's** symmetry sectors before any Hamiltonian block is diagonalized. The isolated-bulk parity labels alone do not determine these occupied sectors. Chapter 2 now defaults to `ProjectionCutoff=10^-28` (a squared projection-weight threshold), removing numerical zero projections. The engine still reports every sector weight, skipped total weight and reconstruction error, and rejects discarded weight above `EDTolerance^2`. This is tolerance-controlled pruning, not an assertion that floating-point projections vanish exactly. Coherences between occupied sectors remain in the global unitary evolution.

### Optional controlled dynamics comparison

Section 2.8 is disabled by default (`testEigenstates=False`). Set its independent `testModel`, `testMu`, `testIDs`, `testTimes` and `testEdge`, then enable it. It builds its own catalog and response, evolves each chosen eigenvector only in occupied coupled-model sectors, and retains full records and projection diagnostics. It compares susceptibility with the **sampled-window maxima** of concurrence, negativity and edge entropy and produces concurrence trajectories labeled by state ID. IDs copied from section 1B require matching $L,J,\Delta,\mu$.

Choose several low/high-susceptibility states in a narrow energy window, then vary weak $g$ and observation time. Maxima over a finite grid are not all-time maxima. High susceptibility may indicate dephasing/edge–bulk entanglement without generating edge–edge entanglement; the comparison is intended to reveal precisely this distinction.

### Validation of this addition

The new functions were executed in Wolfram. For a four-spin bulk at $\Delta=0.45$, the six-state zero-magnetization catalog is complete; a weak-field example gives infidelity $4.80692818\times10^{-7}$ versus the prediction $4.80693067\times10^{-7}$ at $g=10^{-3}$. The degeneracy guard identifies unresolved states at $\Delta=0$.

A total-$L=6$ saved-eigenvector preparation with product $|++\rangle$ edges agrees with independent full-Hamiltonian matrix-exponential evolution at $t=0,0.1,2$ to $4.4\times10^{-16}$. Only 4 of 28 available coupled-model sectors were diagonalized; discarded weight was $9.2\times10^{-34}$ or smaller. These are small-system checks, not evidence that susceptibility predicts entanglement. The optional two-eigenstate comparison was also executed at total L=4. Plot-data assembly generated all six expected figure entries with rendering stubbed; full-size scans, parallel paths and MaTeX rendering were not validated in this revision.


## 3. Chapter 2: dynamics of a chosen global state

### 3.1 Independent preparation

`dynModel`, `dynTime`, `dynSettings`, and `dynInitial` belong only to this chapter. The initial default is $L=8,J=1,\Delta=0.45,g=0.1$, with a Bell edge pair and a Néel bulk.

There is one active `dynInitial` assignment. Two modes are available:

- `Product`: an arbitrary normalized edge vector or density matrix, together with an independently chosen bulk preparation.
- `Full`: an arbitrary normalized full-chain vector or density matrix in physical site order, including correlated and entangled edge/bulk states.

`stateFactor` represents a density matrix as $\rho=AA^\dagger$. A pure state uses one column. A mixed state is checked for normalization, Hermiticity, and positivity, then factored using its eigensystem. `joinStates` restores physical site order when combining an edge factor and a bulk factor.

Bulk choices are Néel, polarized, helix, infinite temperature, Gibbs, thermal typicality, and a user-supplied state. Gibbs and typicality preparations use the field-free $L-2$-site bulk Hamiltonian $H_M$, without either end interaction.

For Gibbs preparations the weights are $e^{-\beta E}/Z$, with an energy shift for numerical stability. Typicality uses a specified number of complex Gaussian columns filtered by $e^{-\beta H_M/2}$, followed by a common normalization. It is a finite-sample approximation to the requested thermal state. Its subsequent evolution and channel are exact for the sampled, normalized bulk density matrix. The seed and sample count are explicit.

### 3.2 Helix inputs

The reused state constructor is

$$
|\mathrm{helix}\rangle=
\bigotimes_{j=1}^{n}\left[
\cos(\theta/2)|\uparrow\rangle+
 e^{i[\phi+(j-1)q]}\sin(\theta/2)|\downarrow\rangle\right].
$$

`q`, `theta`, and `phi` are independently editable. For a bulk helix, $\phi$ is the azimuth at physical site 2; for a full-chain helix, it is the azimuth at site 1. The notebook contains explicit commented replacements for both cases.

These are initial states of the field-free edge-coupled model. No matched transverse field or rigid-rotation protocol is imported from the separate helix file. The helix is not presumed to be an energy eigenstate.

### 3.3 Project first, diagonalize occupied sectors

`occupiedEngine` computes

$$
A_s=B_s^\dagger A_0,\qquad
w_s=\frac{\|A_s\|_F^2}{\|A_0\|_F^2}
$$

before diagonalizing a sector Hamiltonian. Every sector's label, dimension, weight, and diagonalization decision are recorded. Only sectors above `ProjectionCutoff` are diagonalized. The default cutoff is zero, so only exactly zero projections are skipped; numerical roundoff can therefore retain very small weights. If a positive cutoff is chosen, the discarded total weight is reported and must not exceed `EDTolerance^2`.

The setup still constructs/audits the sector bases needed for the projections; it does not solve the Hamiltonian eigensystem in unpopulated sectors.

For an occupied sector, let $V_s$ contain its energy eigenvectors and $C_s=V_s^\dagger A_0$ its initial coefficients. Evolution is

$$
A(t)=\sum_s V_s\,\mathrm{diag}(e^{-iE_{s,n}t})\,C_s,
$$

embedded into the full physical basis. The sectors are summed at the amplitude/factor level before taking a partial trace. Thus coherences between different magnetization or parity sectors are preserved. Projection does not replace the initial state by an incoherent sector mixture.

The engine checks reconstruction at $t=0$. `stateInputDiagnostics` also reports purity, mean energy, energy variance, and the Frobenius norm of $[H,\rho(0)]$. The commutator norm is evaluated in a QR-compressed span of the initial factor and its Hamiltonian image. A stationary initial density matrix produces a visible diagnostic rather than unexplained constant observables.

### 3.4 ZZ simplification

The exact propagator has the controlled form

$$
U(t)=\sum_s|s\rangle\langle s|_{AB}\otimes e^{-iH_s t},
\qquad s\in\{00,01,10,11\}.
$$

The ZZ engine works inside these conditional-bulk blocks and their compatible refinements. It never diagonalizes an exchange-like block that combines different fixed edge spins.

For any initial state, including a correlated one, the four joint edge populations are constants, both edge z magnetizations are constants, and the boundary spin currents vanish. Edge coherences can evolve through overlaps between conditional bulk trajectories. Fixed populations therefore do not imply a stationary reduced density matrix.

The ZZ boundary-current entries are assigned exactly zero because $[H,S_1^z]=[H,S_L^z]=0$. This is an analytic consequence of the model, not an omitted current calculation or a numerical cutoff. Exchange currents are evaluated from the two-site current operators.

### 3.5 Retained state observables

At every time the notebook stores $\rho_{AB}(t)$, $\rho_A(t)$, $\rho_B(t)$, and the four eigenvalues of $\rho_{AB}$. This includes all populations and complex coherences.

The scalar observables retain the benchmark's meanings:

| Family | Quantities |
|---|---|
| Entropies | $S_A,S_B,S_{AB}$ and $I(A:B)=S_A+S_B-S_{AB}$ |
| Edge entanglement | Concurrence and negativity of the two-qubit edge state |
| Purities | $\mathrm{Tr}\rho_A^2$, $\mathrm{Tr}\rho_B^2$, $\mathrm{Tr}\rho_{AB}^2$ |
| Local observables | All six $\langle\sigma_A^a\rangle,\langle\sigma_B^a\rangle$ |
| Correlations | All nine $\langle\sigma_A^a\sigma_B^b\rangle$ and their connected counterparts |
| Reduced overlap | $\mathrm{Tr}[\rho_{AB}(0)\rho_{AB}(t)]$ |
| Full overlap | $\mathrm{Tr}[\rho(0)\rho(t)]/\mathrm{Tr}[\rho(0)^2]$ |
| Currents | Left and right boundary spin currents |
| Validation | Trace error, Hermiticity error, and minimum edge-state eigenvalue |

For a pure global state, `FullOverlap` is the survival probability $|\langle\psi(0)|\psi(t)\rangle|^2$. For a mixed state it is a normalized Hilbert–Schmidt overlap, not Uhlmann fidelity. Likewise the reduced overlap is not generally a fidelity.

For a pure global state, $S_{AB}$ is the edge/bulk entanglement entropy. For a mixed global state, it is the entropy of the edge marginal and is not, by itself, an entanglement measure.

Local expectations and correlations use **Pauli matrices**. The Hamiltonian and currents use $S=\sigma/2$. In the exchange notebook,

$$
j_{j\to j+1}=g(S_j^xS_{j+1}^y-S_j^yS_{j+1}^x)
$$

on either boundary bond. The sign convention gives $d\langle S_1^z\rangle/dt=-\langle j_{1\to2}\rangle$ and $d\langle S_L^z\rangle/dt=\langle j_{L-1\to L}\rangle$. Currents come from operator expectations, not finite differences.

Concurrence uses the benchmark's factor-based singular-value algorithm. Negativity sums the magnitudes of the negative eigenvalues of the partial transpose. Reduced-state traces and positivity are checked during the trajectory.

## 4. Chapter 3: dynamics of the dynamical map

### 4.1 Preparation, conventions, and independence

This chapter has its own `mapModel`, `mapBulk`, `mapTime`, `mapSettings`, and `mapProbes`. Its defaults remain $L=6,J=1,\Delta=0.45,g=0.1$, with a Néel bulk. Evaluate the definitions and this chapter's inputs/setup; no static scan, state trajectory, or QST result is required.

The acceleration guide uses $L$ for the **bulk** length. These files continue to use $L$ for the **total** length, so the bulk has $L-2$ spins. The Hamiltonians are unchanged: $S^\alpha=\sigma^\alpha/2$, the exchange end-bond flip amplitude is $g/2$, and no magnetic fields are introduced. Thus the guide's exchange parameter maps to $g/2$ here; its ZZ parameter maps to $g$. The guide's field-retuning examples are excluded. In particular, $\Delta=1$ in the bulk does not give full-chain SU(2) symmetry with the present XY-only exchange end bonds.

For a fixed bulk preparation $\rho_M=BB^\dagger$, the channel is

$$
\Phi_t(X)=\operatorname{Tr}_M[U(t)(X\otimes\rho_M)U(t)^\dagger].
$$

`chBuild` takes the model, bulk factor, and settings. **There is no initial edge state argument.** The channel acts on all edge inputs. `mapProbes` supplies two edge density matrices for the BLP witness and a one-spin B state for signalling/retention; those probes are applied after construction and cannot affect the Choi matrix. Correlated edge/bulk initial states remain a chapter-2 task: this product-preparation channel is not a universal map for their edge marginals.

The existing normalization and ordering are retained:

$$
C_t=\sum_{i,j=1}^4|i\rangle\langle j|\otimes\Phi_t(|i\rangle\langle j|),
\qquad \operatorname{Tr}C_t=4,\qquad \chi_t=C_t/4.
$$

`Choi` stores $C_t$; `ChoiState` stores $\chi_t$, with input-reference index first and output-edge index second. `Transfer` uses row-major flattening, $\operatorname{Flatten}[\Phi_t(X)]=S_t\operatorname{Flatten}[X]$. The engine checks $S_0=I_{16}$, Choi positivity, and output partial trace $\operatorname{Tr}_{\rm out}C_t=I_4$.

The new `ch`-prefixed functions belong to chapter 3. The legacy channel helpers remain available for the independent QST chapter. The preceding channel-only revision left chapters 1, 2, and 4 byte-for-byte unchanged. The present audit patch changes the context wrapper, a gap-diagnostic count, the legacy map embedding used by QST, and chapter 4. The accelerated `ch`-prefixed channel implementation and its run cells are retained.

### 4.2 Bulk factors, symmetries, and reusable setup

All previous preparation choices remain available, including helices, arbitrary supplied pure/mixed bulk states, Gibbs states, infinite temperature, and the existing `"Typical"` ensemble. Gibbs factors are now assembled by bulk magnetization sectors, so each column has definite magnetization and exact structural zeros. This represents the same Gibbs density matrix; positive and negative finite inverse temperature are supported with an energy shift that prevents exponential overflow.

A separate optional preparation is

```mathematica
mapBulk=<|"Type"->"SectorTypical","Beta"->1.,
  "SamplesPerSector"->2,"Seed"->1234|>;
```

It draws independent complex Gaussian vectors in every bulk magnetization sector, applies $e^{-\beta H_M/2}$, and normalizes the entire factor together. It is a stochastic thermal estimator with a sector-sparse factor. It has a different finite-sample ensemble from the existing global-vector `"Typical"` preparation; the latter has not been silently redefined. Increase samples and check seed dependence when using either approximation.

With `"Symmetries"->Automatic`, real eigensystems are obtained in magnetization blocks, using reflection parity when it acts inside the relevant block and reusing spin-flip-related eigenvectors. `None` disables these extra shortcuts while retaining magnetization blocks. The complete occupied subspace is retained, including all parity components. These channel blocks need not reproduce the fine labels used for static eigenstate classification: channel propagation retains the whole necessary space, and is invariant under a basis change within a degenerate energy space.

The setup checks realness, Hermiticity, magnetization conservation, reflection invariance where used, eigenpair residuals, and normalization using `EDTolerance`. Displayed sector weights sum to one; the underlying all-edge-input construction has total squared norm four, so its block weights are divided by four.

Only exactly zero factor columns are pruned inside a block. This chapter requires `ProjectionCutoff->0`; it does not approximate the channel by discarding weakly populated sectors. This refers to the supplied factor after its preparation/validation, not to an assertion that a stochastic thermal factor is an exact Gibbs state.

The eigensystem and projected coefficients are constructed once. Changing only the time grid requires rerunning the time-evaluation cell, not diagonalizing again. Changing the model, bulk factor, or setup settings requires rerunning setup.

### 4.3 Exchange: active-column propagation in bounded batches

For a bulk factor with $r$ columns, there are $4r$ edge-input/bulk-factor combinations. `chExEngine` identifies which combinations have support in each total-magnetization block and constructs only those block columns. It does not first allocate a full $2^L\times4r$ initial array. A bulk factor supported in one magnetization sector populates three total-magnetization sectors through the four possible edge inputs.

Within each occupied block, propagation uses

$$
A_q(t)=V_q e^{-iE_qt}C_q.
$$

Several times are packed into one matrix product. Since $V_q$ is real, the product is evaluated as two real matrix multiplications, $V_q\Re Z+iV_q\Im Z$. Complex initial bulk phases remain fully represented in $C_q$; no bulk-sector dephasing is performed for exchange.

The block amplitudes are scattered directly into

$$
W_{(i,\beta),(b,k)}(t)=\langle\beta,b|U(t)|i,b_k\rangle,
\qquad C_t=W_tW_t^\dagger.
$$

This Gram construction preserves positivity without sixteen separate matrix-unit evolutions. The analytic derivative uses the same batches with phases multiplied by $-iE$, followed by $\dot C=\dot W W^\dagger+W\dot W^\dagger$.

`BatchTimes` bounds the outer batch, and `MaxEntries` additionally limits the number of entries in one batched W array. One time slice requires $16\,2^{L-2}r$ complex entries. If even one slice exceeds the limit, the code stops with an explicit instruction to raise the limit or reduce factor rank; it does not silently exceed it. Derivative arrays, eigensystems, coefficients, and temporary products require additional memory.

### 4.4 ZZ: conditional-bulk propagation or overlap contraction

The fixed edge configurations $s$ give four conditional bulk Hamiltonians $H_s$. The exact channel is

$$
X_s(t)=e^{-iH_st}B,\qquad
K_{ss'}(t)=\operatorname{Tr}[X_s(t)X_{s'}(t)^\dagger],\qquad
\Phi_t(X)=K(t)\circ X.
$$

`DephasingFactors` stores $K$ (denoted $F_t$ in the new plot labels). Its diagonal is one. This is a pure-dephasing, unital channel for every bulk preparation, with Choi rank at most four. The Choi matrix is supported on $|s\rangle_{\rm in}|s\rangle_{\rm out}$ and equals $K$ on that support.

Only occupied **bulk** magnetization blocks are required. Coherences of $\rho_M$ between distinct bulk magnetizations cannot contribute after the bulk trace, because both conditional Hamiltonians preserve that magnetization. This exact simplification is specific to this ZZ construction.

Two interchangeable methods compute the same $K$ and analytic $\dot K$:

- `"ZZMethod"->"Direct"`: propagate active factor columns in real-arithmetic batches, then contract the four conditional states. `MaxColumns` limits time-count times active factor columns in each block; a single block slice must fit. This is usually favorable for pure or low-rank preparations.
- `"ZZMethod"->"Overlap"`: precompute six conditional-eigenbasis overlap matrices. For block $q$, with $C_{s,q}=V_{s,q}^\dagger B_q$,
  $$
  G^{ss'}_{nm}=(V_{s,q}^{T}V_{s',q}^{*})_{nm}
                 (C_{s,q}C_{s',q}^{\dagger})_{nm},
  \qquad
  K_{ss'}(t)=\sum_{q,n,m}G^{ss'}_{nm}e^{-i(E_{s,n}-E_{s',m})t}.
  $$
  The time loop no longer depends on the bulk-factor rank. Derivatives multiply each term by $-i(E_{s,n}-E_{s',m})$.

`Automatic` compares the guide's arithmetic-cost estimates for the requested grid and selects overlap only when it is cheaper and its cache fits `MaxOverlapEntries`. The six caches contain $6\sum_qd_q^2$ entries. Selection is an estimate, not a measured guarantee of the faster method on every machine. The chosen method and estimates are retained in `mapChannelRun`. An explicitly requested oversized overlap cache stops with a clear message.

The overlap cache is built once per trajectory call, not at every time. `BatchTimes` also bounds overlap evaluation and optional time-reversal evaluation. The direct and overlap methods can be forced for numerical comparisons.

### 4.5 Observables and the complex-bulk time-reversal check

The complete Choi and transfer matrices and analytic transfer derivative are retained at every time. Previously available observables retain their meanings and normalization.

| Diagnostic | Meaning |
|---|---|
| Choi eigenvalues, CP/TP errors | Numerical validity of the channel construction |
| `KrausRank` | Number of Choi eigenvalues above the relative cutoff |
| `MapEntropy` | $S(C_t/4)$, entropy exchange for a maximally entangled reference/edge input |
| Choi PT minimum and negativity | Entanglement between the input reference and output edge pair |
| Minimum transfer singular value | Conditioning/invertibility of the channel |
| Generator and canonical rates | Local divisibility diagnostics where the map is invertible |
| RHP diagnostics | Instantaneous and finite-step tests of CP divisibility |
| BLP distance and accumulation | Distinguishability revival for the selected pair of edge inputs |
| Unitality and mixed-input entropy | Response to the maximally mixed edge state |
| U(1), reflection, flip, RF covariance errors | Compatibility of the map with those transformations |
| `TimeReversalError` | Correct conjugate-bulk relation below, at sampled check times |
| Double-quantum response | Complex $00,11$ output matrix element from input $|00\rangle\langle11|$ |
| Signalling and retention | Largest Bloch response singular value at B and A for the specified B preparation |
| ZZ `MinDephasingEigenvalue` | $\lambda_{\min}(K_t)/4$, the Choi minimum restricted to its four-dimensional support |
| ZZ `IntermediateSupportMinimum` | Minimum support eigenvalue of the normalized intermediate Choi matrix, when the earlier map is invertible |

ZZ entropy and Choi support eigenvalues use a $4\times4$ problem, with twelve exact zero Choi eigenvalues appended where the full spectrum is required. The partial-transpose spectrum is $1/4$ four times and $\pm|K_{ss'}|/4$ for the six unordered pairs. The minimum transfer singular value is $\min_{ss'}|K_{ss'}|$. These shortcuts avoid repeated $16\times16$ spectral decompositions.

Map entropy vanishes for a unitary channel and is $\log16$ for complete depolarization. Choi negativity concerns the $4\otimes4$ reference/output split, not entanglement of A with B. Negative partial transpose proves non-entanglement-breaking behavior; PPT alone does not prove separability for a general $4\otimes4$ Choi state.

Covariance generally requires both the Hamiltonian and bulk preparation to respect the transformation. A helix can break U(1) covariance in exchange. ZZ always preserves computational-basis coherences independently, and is U(1)-covariant even for a non-invariant bulk preparation. Without suitable covariance in exchange, the double-quantum matrix element remains defined but need not describe a closed one-dimensional coherence mode.

Because H is real in the computational basis, the valid time-reversal identity is

$$
S_{-t}^{[\rho_M]}=\left(S_t^{[\rho_M^*]}\right)^*.
$$

The same-bulk identity is justified only when $\rho_M=\rho_M^*$. The previous comparison could therefore report a nonzero residual for a valid complex helix preparation. `chTimeReversal` conjugates the projected initial coefficients, reuses the same eigensystems, and tests the correct relation. `"Endpoints"` checks first, middle, and last times by default; `"All"` checks every time; `"None"` skips the diagnostic. Unsampled rows carry `Missing["NotSampled"]`. This diagnostic is outside the main trajectory loop.

### 4.6 Generator, canonical rates, and memory witnesses

Where $\sigma_{\min}(S_t)$ passes `MinimumSingularValue`,

$$
\mathcal L_t=\dot S_tS_t^{-1}.
$$

Exchange uses a right-side `LinearSolve`. A cached $225\times256$ contraction matrix projects this generator onto the traceless two-qubit Pauli basis $G_{\mu\nu}=\sigma_\mu\otimes\sigma_\nu/2$. Thus constructing the $15\times15$ GKS matrix is a single matrix-vector product, instead of rebuilding 225 Kronecker products at every time.

ZZ uses entrywise division $\ell=\dot K\oslash K$. If Q has orthonormal rows spanning the subspace perpendicular to $(1,1,1,1)$, its nontrivial canonical rates are the three eigenvalues of $Q\ell Q^T$; twelve zero rates complete the same fifteen-rate convention. The full generator remains available as $\operatorname{diag}(\operatorname{Flatten}\ell)$.

No pseudoinverse is substituted near singularity. The generator, rates, and affected intermediate-map diagnostics become `Missing`. With the stated basis normalization,

$$
\mathrm{RHPInstantaneous}(t)=\frac12\sum_k\max(0,-\gamma_k(t)).
$$

Finite-step RHP is $\max(0,\|C[S_{t_{k+1}}S_{t_k}^{-1}]\|_1/4-1)$; the code separately stores this quantity divided by the time step. For ZZ, the intermediate map is the Schur multiplier $K_{t_{k+1}}\oslash K_{t_k}$, so its support minimum and trace norm require only a $4\times4$ matrix. Finite-step values depend on the grid and are not exact derivatives. Missing intervals are not replaced with zero.

BLP uses $D(t)=\tfrac12\|\Phi_t(\rho_1-\rho_2)\|_1$ and accumulates positive sampled increments in one pass. It is a grid-dependent lower bound for the selected pair, not an optimization over all pairs. ZZ computational-basis probes give a constant witness; the input cell retains coherence-sensitive alternatives.

Signalling and retention maximize the corresponding $3\times3$ Bloch response over the A-input direction at fixed $\rho_B$. `SignalOnset` remains the interpolated first crossing of the chosen threshold, not a fitted propagation velocity or a strict causal boundary.

### 4.7 Spectral infinite-time average

`chAverage` computes the spectral infinite-time average, not a late-time sample. The absolute energy threshold is `AverageEnergyTolerance` times $\max(1,\max|E|)$ over the relevant energies. Clusters have bounded total width; adjacent near-equalities cannot form arbitrarily wide chains. Cluster membership, widths, and tolerances are retained. Check tolerance dependence when small physical splittings approach this threshold.

**Exchange:** energies are clustered globally across all occupied total-magnetization blocks. Degenerate states in different blocks must remain coherent within a common energy projector. Singleton contributions are accumulated in vectorized form: writing the eigenvector in edge/bulk form as $V_n$ and its input coefficients as $B_n$, each contributes $(B_nB_n^\dagger)\otimes(V_nV_n^\dagger)$ to the unnormalized Choi matrix. A non-singleton cluster builds one coherent W factor from all its members and adds $WW^\dagger$. This avoids one full evolving array per isolated level while preserving cross-sector degeneracies. Highly degenerate spectra still require substantial work and transient storage. The one-slice Choi-factor size is checked against `MaxEntries` before averaging.

**ZZ:** the bulk trace eliminates cross-bulk-magnetization terms exactly. Each occupied bulk sector is therefore clustered across its four conditional spectra, preserving degeneracies between different edge configurations. Matching-energy entries of the six overlap matrices give the averaged off-diagonal $K_{ss'}$; the diagonal remains one. The cache is checked against `MaxOverlapEntries`. No full global factor is allocated for every energy level.

Both results retain `Choi`, `ChoiState`, `Transfer`, TP error, map entropy, and Choi partial-transpose minimum. Exchange additionally reports singleton/degenerate cluster counts; ZZ reports clusters by bulk sector and matched off-diagonal energy-pair counts. Normalizations remain identical to the instantaneous channel.

### 4.8 Numerical checks and measured small-system timing

The revised code was executed in Wolfram. For each regime, 72 cases combined $L=4,5$, $\Delta=0,-1,0.45$, $g=-0.2,0.3$, complex random rank-two/Gibbs/infinite-temperature bulk factors, and symmetry shortcuts on/off. Reference Choi matrices used independently assembled full Hamiltonians, matrix-unit embeddings, matrix exponentials, and explicit bulk traces.

| Check | ZZ maximum error | Exchange maximum error |
|---|---:|---:|
| 72 direct-reference Choi comparisons | $7.78\times10^{-16}$ | $1.34\times10^{-14}$ |
| Analytic derivative versus previous implementation | $9.01\times10^{-16}$ | $1.50\times10^{-15}$ |
| Spectral average versus previous implementation in the $L=4$ subset | $8.89\times10^{-16}$ | $8.65\times10^{-15}$ |
| Negative-temperature Gibbs density from the new factor versus the old factor | $1.67\times10^{-16}$ | $1.67\times10^{-16}$ |
| $L=5$, $\Delta=1$ complex-helix average versus explicit full-space energy projectors | $1.25\times10^{-15}$ | $1.41\times10^{-14}$ |
| Retained numerical diagnostics at $L=6$, Gibbs bulk, 21 times | $4.52\times10^{-15}$ | $7.52\times10^{-15}$ |

Forced ZZ overlap and direct evolution agreed within $2.24\times10^{-16}$ for Choi matrices and $1.39\times10^{-16}$ for derivatives in the complex-helix check. At $L=8$ with an infinite-temperature bulk, the forced ZZ methods also agreed within $2.23\times10^{-16}$ for Choi matrices and $7.32\times10^{-17}$ for derivatives; automatic selection chose overlap on a 128-point grid. The conjugate-bulk time-reversal residual was zero to machine arithmetic in these tests. The sector-typical constructor reproduced the same seeded factor and normalized it within $4.45\times10^{-16}$.

A single small timing run at $L=6$, $\Delta=0.45$, $g=0.3$, Gibbs $\beta=0.8$, and 21 times gave:

| Computation | ZZ previous / revised | Exchange previous / revised |
|---|---:|---:|
| Full trajectory diagnostics | 1.044 s / 0.143 s | 1.036 s / 0.116 s |
| Spectral time average | 0.0294 s / 0.00310 s | 0.0154 s / 0.00473 s |

These are illustrative timings in the test kernel, not a large-system benchmark or guaranteed speedup. ED is outside the trajectory timings; the revised GKS cache was warmed. The previous trajectory also recomputed a same-bulk negative-time map at every point, whereas the revised time-reversal check is separate and sampled. Thus the timing improvement includes both faster algebra and removal of that repeated diagnostic work.

The actual numerical input/setup/trajectory/time-reversal/average cells were also executed through temporary `.m` files in fresh kernels at $L=5$. Both completed with TP errors below $7.78\times10^{-16}$ and left an unrelated state-dynamics model unchanged. Scope comparisons for that preceding channel-only revision confirmed that every byte before chapter 3 and from chapter 4 onward was unchanged. That historical statement does not apply to the subsequent QST/context patches described below. Full `.m` syntax checks, near-singular-map handling, and memory-limit rejection passed in both regimes. New channel plot wrappers were exercised with native text labels; MaTeX integration remains untested in the available kernel.

## 5. Chapter 4: quantum-state texture of the edge spins

QST here means **quantum-state texture**, not state transfer or tomography. Its computational basis is the same edge basis $00,01,10,11$ used throughout the files. The chapter follows the supplied `QST_chapter_edges.md`. It contains two independently initialized calculations: A, a particular reduced-state trajectory, and B, a channel calculation. Neither uses evaluated results of chapters 1–3, and neither requires the other QST calculation.

### 5.1 Texture quantities and basis conventions

The fixed textureless edge ket and density matrix are

$$
|f_1\rangle=\tfrac12(1,1,1,1)^T=|+\rangle_A|+\rangle_B,
\qquad f_1=|f_1\rangle\langle f_1|.
$$

Only the edges are prepared in this state. The bulk has its own specification; the QST preparation helper rejects a bulk indistinguishable from its all-plus pure state within the validation tolerance. It does not silently substitute another bulk preparation.

For a density matrix of dimension $D$,

$$
\Sigma(\rho)=\operatorname{Re}\sum_{ij}\rho_{ij},\qquad
T_F(\rho)=1-\frac{\Sigma(\rho)}D,\qquad
R(\rho)=-\ln\frac{\Sigma(\rho)}D.
$$

The pair uses $D=4$ and either individual edge uses $D=2$. These quantities depend on the chosen basis. For the pair, $\Sigma(\rho_{AB})/4=\langle f_1|\rho_{AB}|f_1\rangle$: the probability of the all-zero outcome after a Hadamard on each edge. The code computes the matrix-element probability separately as a consistency check.

Grand sums and textures are retained without clipping. An exactly zero grand sum produces `Infinity` for rugosity. A negative grand sum, including a small negative value caused by roundoff, produces `Missing["NegativeGrandSum",...]` for rugosity while retaining the raw grand sum. Positivity and physical-range checks use the chosen numerical tolerance; they do not modify the stored matrices.

### 5.2 A: open dynamics starting from textureless edges

The global preparation is

$$
\rho(0)=f_1\otimes\rho_E,
\qquad \rho_{AB}(t)=\operatorname{Tr}_E[U(t)\rho(0)U(t)^\dagger].
$$

The editable bulk options are Néel, polarized, planar or conical helix, polar helix, a supplied ket/density matrix, infinite temperature, Gibbs, and Gibbs typicality. A polar helix is constructed from explicit `ThetaList` and `PhiList` arrays for the $L-2$ bulk sites, with local kets $(\cos(\theta_j/2),e^{i\phi_j}\sin(\theta_j/2))$. This is a state constructor, not an assertion of stationarity or an exact phantom-state condition. The other preparations reuse chapter 2's factorization routines unchanged.

`qstAFactor` combines the fixed edge ket with the bulk factor in physical site order. Initial edge density matrix and global purity are retained. `occupiedEngine` projects this factor onto every available symmetry basis before diagonalization, records each weight, and diagonalizes only occupied blocks. `SectorCut` is the same projection threshold as earlier chapters: the total discarded weight must still satisfy the engine's ED-tolerance bound. Setting a cutoff does not authorize a large uncontrolled loss of norm.

In ZZ, the sectors are conditional bulk spaces with fixed edge configurations. Four bulk evolutions generate the edge coherences; populations remain fixed. In exchange, the full-chain magnetization/parity sectors are used, preserving their mutual coherent amplitudes. In both cases the reduced dynamics is obtained from exact global unitary phases, without a Markov or master-equation approximation.

The new QST batching wrapper calls `evolveBatch` and the existing `stateRecord` unchanged. It bounds the number of time slices in a temporary full-factor array and optionally maps batches across parallel kernels. All usual state observables survive in each record, including entropies, purity, concurrence, negativity, overlaps, local Pauli expectations, correlations and currents. The wrapper adds:

| Record family | Meaning |
|---|---|
| `SigmaAB`, `SigmaA`, `SigmaB` | Pair and single-edge grand sums |
| `TextureFAB`, `TextureFA`, `TextureFB` | Corresponding fidelity textures |
| `RugosityAB`, `RugosityA`, `RugosityB` | Corresponding logarithmic textures |
| `Nonmultiplicativity` | $\Sigma_{AB}-\Sigma_A\Sigma_B$ |
| `L1Coherence` | $\sum_{i\ne j}|\rho_{AB,ij}|$, not a texture monotone substituted for QST |
| `F1Probability` | Directly evaluated $\langle f_1|\rho_{AB}|f_1\rangle$ |
| Comparison values | $C_{\ell_1}/3$, $1-\operatorname{Tr}\rho_{AB}^2$, and $S(\rho_{AB})/\ln4$ |

The independently editable fine short-time grid uses the same initialized engine and does not repeat ED. Main and fine trajectories have separate time snapshots and validations.

The optional infinite-time state average uses global energy clusters across all occupied sectors, including different fixed-edge configurations in ZZ. For an eigenvector $|n\rangle$, let $V_n$ be its $4\times2^{L-2}$ edge/bulk coefficient matrix and $c_n$ its row of initial-factor coefficients. The revised `qstAverageState` evaluates

$$
\overline\rho_{AB}=\sum_G\sum_{n,m\in G}
(c_n\cdot c_m^*)V_nV_m^\dagger.
$$

For a singleton cluster, the contribution is $\|c_n\|^2V_nV_n^\dagger$. For a degenerate cluster, the routine adds all diagonal contributions and each off-diagonal pair together with its conjugate transpose. This is exactly the reduced density matrix of the coherent energy-projected factor, without allocating a full $2^L\times r$ array for every cluster. Cross-sector degeneracies remain coherent. Embedding each eigenvector costs $O(2^L)$ storage/work, and coefficient contractions still depend on factor rank $r$; large degenerate clusters retain a quadratic number of pairs. The routine reports the averaged edge matrix, $\overline\Sigma_{AB}$, singleton/degenerate cluster counts, widths, and tolerances.

Both this embedding and `mapLevelFactors` explicitly assign complex machine entries (`+0.I`) into their complex zero arrays. This prevents Mathematica from unpacking the arrays when an eigenvector is real. The external audit tested these optimizations: its exchange $L=12$ Néel map average changed from 14.0 s to 0.90 s, and its exchange $L=10$ infinite-temperature state average from 3.98 s to 0.13 s, with unchanged results to its reported precision. These are the auditor's timings for its test engine, not measurements of the complete revised files or guarantees for a user's machine. The average is still evaluated once per requested cluster tolerance (three by default).

The average rugosity is a **separate finite-window estimate**: trapezoidal integration of the sampled $R(t)$ with respect to physical time, divided by the window length. A logarithmic sampling grid is not averaged with equal sample weights. A sampled exact zero is reported as `Missing["SingularSamplesRequireQuadrature"]` for this estimate, since an isolated logarithmic divergence need not make its time integral divergent. Invalid values are also reported, never dropped. Rugosity of the averaged state is not used as a surrogate.

### 5.3 B: channel texturizing power and response operator

Section B independently fixes the model and bulk state and constructs the channel through `mapBuild`. It does not use A's factor, eigensystem, time grid, or edge trajectory. There is no edge-state parameter in channel construction. The unique free input $f_1$ is applied afterward, defining

$$
\Omega_{T_F}(t)=1-\frac{\Sigma[\Phi_t(f_1)]}{4},
\qquad \Omega_R(t)=-\ln\frac{\Sigma[\Phi_t(f_1)]}{4}.
$$

Three independently contracted formulas are retained and compared at every time:

$$
\Sigma[\Phi_t(f_1)]
=\frac14\mathbf1_{16}^T S_t\mathbf1_{16}
=\frac14\sum_{ij}(C_t)_{ij}.
$$

Here $C_t$ is the unnormalized Choi matrix with trace four, not the normalized Choi state. The factors of four follow from `Flatten[f1]=ones16/4`. The first route actually forms `applyTransfer[S,f1]`; the other two use the transfer and Choi arrays.

The response operator is

$$
\mathcal R_t=4\Phi_t^*(f_1),\qquad
\Sigma[\Phi_t(\rho)]=\operatorname{Tr}(\mathcal R_t\rho).
$$

With the row-major flattening convention it is `Transpose[ArrayReshape[Total[S],{4,4}]]`. The transpose is required by the trace contraction; no conjugation is inserted in that formula. Hermiticity and positivity are checked independently.

The ordered eigenvalues of $\mathcal R_t$ give the extreme achievable output grand sums. The code retains the full eigenvectors and their $f_1$ overlaps. Inside a degenerate eigenspace, an individual eigenvector and its overlap are basis dependent; the extreme eigenvalue remains well defined. These are bounds over input edge states, not over bulk preparations or Hamiltonians.

Two channel properties are tested separately:

| Property | Exact condition | Recorded distances |
|---|---|---|
| Free channel | $\Phi_t(f_1)=f_1$ | Max-entry and trace norms of $\Phi_t(f_1)-f_1$ |
| Texture preservation | $\mathcal R_t=4f_1$ | Max-entry and trace norms of $\mathcal R_t-4f_1$ |

The trace norms here have **no extra factor of one half**. Also recorded is $\operatorname{Tr}\mathcal R_t-4$; its vanishing means that the output grand sum averaged over any orthonormal edge basis equals one. It does not imply texture preservation by itself.

B retains $S_t,C_t,C_t/4,\Phi_t(f_1)$, the response eigendata, all three route values, and the comparison diagnostics: map entropy, unitality, mixed-input entropy, covariance violations, minimum Choi eigenvalue, TP error, and numerical Kraus rank. In ZZ it also retains the $4\times4$ conditional Gram matrix $\Gamma$, with $S_t=\operatorname{diag}(\operatorname{Flatten}\Gamma)$.

The three routes are algebraically equivalent rearrangements of the same transfer matrix. Their agreement validates index/reshuffling conventions; it is not an independent physical channel calculation.

For ZZ, the exact reductions are

$$
\mathcal R_t=\Gamma_t^T,\qquad
\Phi_t(f_1)=\Gamma_t/4,\qquad
\Omega_{T_F}(t)=\frac34-\frac18\sum_{s<s'}\Re\Gamma_{ss'}(t).
$$

The transpose here is an ordinary transpose, not a conjugate transpose. For a common bulk preparation, section A's edge state is therefore $\Gamma_t/4$, and the response spectrum equals the spectrum of $\Gamma_t$. Chapter 4 re-expresses the same six complex dephasing overlaps used in chapter 3; it does not add independent dynamical information in ZZ. It supplies their interpretation in terms of texture and its operational input/output quantities. `texturizingRecord` now stores the three residuals `ZZResponseGammaError`, `ZZOutputGammaError`, and `ZZOmegaGammaError`, and `qstMapChecks` requires them to pass.

The optional spectral average reuses the optimized `timeAveragedMap`. It evaluates $\overline\Omega_{T_F}=1-\mathbf1^T\overline S\mathbf1/16$ and the spectrum of $\overline{\mathcal R}$. This spectrum describes the averaged response, not the time average of instantaneous extreme eigenvalues. Sampled-window $\Omega_R$ is integrated separately using the same nonlinear-average precautions as A.

### 5.4 QST validation and independence

A checks the initial $f_1$, physical grand-sum bounds, reduced-state trace/Hermiticity/positivity, and the identities

$$
\Sigma_{AB}=1+\langle\sigma_1^x\rangle+\langle\sigma_L^x\rangle+
\langle\sigma_1^x\sigma_L^x\rangle,
\quad \Sigma_A=1+\langle\sigma_1^x\rangle,
\quad \Sigma_B=1+\langle\sigma_L^x\rangle.
$$

B checks $S_0=I_{16}$, zero initial texturizing power, $\mathcal R_0=4f_1$, route agreement, the response expectation on $f_1$, CP/TP, and response Hermiticity/positivity. A separate reproducible random-density-matrix test checks the response identity at the first, middle and last times. Its sample count and seed can change without repeating ED or evolution.

Both sections compare batched evolution with single-time evaluation at endpoints. Their built-in small-system reference uses full-chain `MatrixExp`, explicit matrix-unit embeddings and an index-sum partial trace at first/middle/last times when $L\le6$; larger runs report it as skipped. Because it uses `regimeHamiltonian`, this reference tests the projection, evolution, and partial-trace pipeline, not independent Hamiltonian assembly. The random-input response test likewise checks the response convention against the stored transfer matrix. Independent Hamiltonian/channel checks are documented separately in §6.

Both A and B now report `SectorReconstructionError` in both regimes. ZZ additionally tests fixed edge populations, diagonal transfer structure, the unit diagonal of $\Gamma$, and the three exact texture identities above. The old `ZZCurrentError` was tautological because the ZZ currents are analytically set to zero. It has been replaced with `ZZEdgeConservationError`, computed from the actual operators as $\max_{j\in\{1,L\}}\|[H,S_j^z]\|_{\max}$. This check runs once per validation call, rather than evaluating prescribed zero current records at every time.

The optional final A/B comparison now checks the complete $4\times4$ matrices $\rho_{AB}^{(A)}(t)$ and `OutputF1` from B, as well as their grand sums. It compares **used** model, bulk specification/factor, actual time grid, and numerical settings, and verifies record times against the snapshots. Parallel/batch settings are excluded because they do not change the intended result. Editing input variables after a run cannot mislabel old records as a new computation.

After a successful comparison, section 4.C automatically writes a separate `ZZ_QST_AB_comparison` or `Exchange_QST_AB_comparison` bundle under B's output parent, with both used snapshots and the comparison result. Earlier A/B bundles may still say `NotCompared`; the separate comparison is the saved, authoritative result for those snapshots. Re-exporting all plots is unnecessary. A and B remain independently executable; the final comparison runs only when both results are available.

Changing the model or bulk requires rebuilding the respective engine. Changing only a time grid starts at that section's evolution cell; the engine is reused. Cluster tolerances affect only the average cells. Random-test count/seed affect only the response validation. Section A's initial factor and section B's map factor never depend on the other section's evaluated variables.


## 6. Computational scope and validation

The main static cost is dense ED inside each symmetry block. A block of dimension $d_s$ requires roughly quadratic storage and cubic diagonalization time. All requested eigenvectors are retained at all scan points, so a long static scan can consume substantially more memory than one eigensystem.

Dynamics diagonalizes only occupied sector Hamiltonians, then evaluates exact spectral phases. Chapter 2 processes time samples one at a time; the revised chapter 3 and the existing QST implementation process bounded batches. Reduced records are retained, avoiding storage of every full global trajectory factor.

For a chapter-3 exchange channel with bulk-factor rank $r$, one Choi factor W has $16\,2^{L-2}r$ complex entries, the same count as $2^L\times4r$. At $L=12$ and full bulk rank, that is about 268 MB (256 MiB) per array. The default `MaxEntries->2^22` would reject that one-slice size; raise it deliberately if resources allow. Derivative arrays, intermediate products, eigensystems, and stored records add to peak memory. Batch limits cap particular arrays, **not total process memory**.

The new ZZ direct engine works block by block, with conditional arrays of size $d_q\times(r_q\times\text{batch times})$. Its overlap engine instead caches $6\sum_qd_q^2$ entries and uses batched energy phases. These costs differ from the old unblocked ZZ full-factor estimate of about 67 MB at $L=12$ and full rank; that old figure was per array, not a peak estimate. The optional average also has its own transient storage costs.

Chapter 3 retains small channel/diagnostic records for every time, so output storage still grows with the number of samples even when transient propagation is batched. The time-reversal check runs separately, using the same eigensystems; by default only three times are checked. QST retains its independent batched evolution engine and memory controls, with the revised spectral-average contractions described in §5.2. Its parallel-worker limits do not cap eigensystem duplication, total records, or process-wide memory. Typicality reduces factor rank by approximating the bulk state; it does not remove the ED cost.

The QST transient-factor batch rule remains

$$
n_{\rm batch}=\min\!\left(n_{\rm requested},\left\lfloor
\frac{\mathrm{MemoryMB}\,2^{20}}{4 n_k\,16\,2^L\,\mathrm{Columns}}
\right\rfloor\right),
$$

where $n_k$ is the active worker count (one in serial mode). In exchange B at $L=12$ with full-rank bulk, `Columns=4096`; one factor is 256 MiB, and the four-array allowance requires `MemoryMB>=1024` for one serial slice. The floor scales with active workers. The default remains 512 for the small default system, with this larger-size requirement documented beside it. A too-small budget now reports the calculated minimum explicitly. It does not promise that the whole run will fit: worker copies, eigensystems, and retained records are additional.

QST now catches a worker/batch abort and rejects malformed, incomplete, or failed batch returns before assembling a trajectory. This prevents a partial result from being treated as a complete run. It is a defensive code change; successful execution on real parallel subkernels remains unverified. QST plot captions consistently say `Exchange`; logarithmic tick generation supplies within-decade ticks, and QST labels reuse the shared cached MaTeX loader rather than retrying package loading for every label.

### 6.1 Earlier implementation checks

The reorganized definitions were executed in a Wolfram kernel and compared with independently assembled full-space references:

| Check | ZZ | Exchange |
|---|---:|---:|
| Complete-spectrum comparisons | 30 | 30 |
| Maximum Hamiltonian discrepancy | 0 | 0 |
| Maximum full-spectrum error | $7.99\times10^{-15}$ | $6.44\times10^{-15}$ |
| Maximum evolved-factor error | $1.67\times10^{-16}$ | $2.36\times10^{-16}$ |
| Maximum reduced-state error | $2.31\times10^{-16}$ | $3.05\times10^{-16}$ |
| Maximum Choi error | $6.55\times10^{-16}$ | $6.13\times10^{-16}$ |
| Maximum infinite-time-average error | $3.34\times10^{-16}$ | $7.70\times10^{-15}$ |
| Maximum trace-preservation error | $4.45\times10^{-16}$ | $6.67\times10^{-16}$ |

The spectrum checks used $L=4,5,6$, $\Delta=0,0.45,0.5,1,-1$, and both signs of $g$. They checked dimension completeness and the small-bulk ZZ Casimir refinement. State tests included complex product states, full correlated pure states, and a complex mixed global state. Channel tests included helix, Gibbs, and infinite-temperature bulks, with explicit matrix-unit embeddings and explicit bulk traces as references.

Additional time-average checks at $L=4$ included $\Delta=0,1,-1,1/2$, preserving cross-sector degeneracies. The QR-compressed initial-state commutator norm also agreed with an explicit full-density-matrix commutator to machine precision.

The analytic transfer derivative agreed with a central finite difference to about $2.3\times10^{-11}$. The GKS implementation recovered the known canonical rates 0.6 and 0.8 for a two-qubit amplitude-damping/dephasing test with the stated basis normalization.

A polarized state diagonalized exactly one occupied sector in each notebook. The actual state- and map-chapter run cells were also tested in fresh kernels with unrelated static variables left untouched; the map chapter left an unrelated state-dynamics model untouched. Initial-factor reconstruction errors in these checks were below $4\times10^{-16}$.

The earlier notebook input boxes were parsed in Wolfram and compared with the tested source expressions. DOS normalization, resolved gap diagnostics, graphics construction, and WXF round-trips were checked. The deliverables are executable `.m` sources. The preceding revisions received the syntax and small-run checks described here; the validation limits of the current audit patch are stated in §6.3.

The independent audit supplied on 15 October reports 40 spectrum comparisons per regime at $L=4$–7, direct-state checks at $L=8$, and channel checks at $L=6$ for five bulk preparations. It found no numerical errors in the prior implementation. Its numerical findings are attributed to that audit; they are not claimed as newly rerun tests here. The revisions address its interpretation, sample-size and performance findings without changing the Hamiltonians, state definitions or channel conventions.

For the preceding QST update, QST was run in a Wolfram kernel for **32 cases per regime**: $L=4,5$, $\Delta=0,0.45,1,-1$, both signs of $g$, and conical-helix or infinite-temperature bulks. The worst combined state/channel/reference/identity error was $1.34\times10^{-15}$ (ZZ) and $4.45\times10^{-15}$ (exchange). The new compact time-average calculation agreed with the prior full-cluster-factor construction within $2.24\times10^{-16}$ in both regimes, including degenerate clusters spanning sectors. Exact-zero rugosity and rejection of the excluded all-plus bulk were checked.

The test kernel has no MaTeX installation or available parallel subkernels. Native-fallback graphics and PNG export were exercised, including singular-rugosity markers; mathematical label mappings and sector-dimension labels were checked. The parallel path reports the lack of kernels explicitly. All 131 mathematical observable/sector label strings compiled with standard LaTeX (`amsmath`, `amssymb`). The actual QST A run cells completed a small ZZ case, including PNG/WXF/parameters export and WXF readback. The actual QST B run cells completed an independent exchange/infinite-temperature case with no A results present. Successful parallel execution, MaTeX integration, full default scans, and large-$L$ performance remain untested. The current chapter-3 tests and limited timing measurements are reported separately in §4.8; they do not extend those earlier validation claims.


### 6.2 Independent assessment supplied on 16 October

The supplied `XXZ_two_regimes_assessment_2026-10-16.md` reports no numerical error in the averaging/QST revision it examined. It used the functions verbatim but supplied independently built fixed-magnetization engines, rather than rerunning these files' sector setup. It examined the version before the latest accelerated chapter-3 revision, so its “unchanged channel code” and one-time-at-a-time performance observations are historical.

The auditor ran 30 cases per regime at $L=6$, with $\Delta=0,0.45,0.5,\pm1$, both signs of coupling, and Néel, complex-helix, and rank-three complex mixed bulks. Every case had energy degeneracies spanning sectors. Full-space references used bit-wise Hamiltonians, energy projectors, explicit embeddings, and index-sum partial traces. Its principal results were:

| Quantity | ZZ | Exchange |
|---|---:|---:|
| Averaged Choi versus full-space projectors | $8.0\times10^{-14}$ | $\le6.8\times10^{-14}$ against corrected/30-digit references |
| Averaged edge state | $1.9\times10^{-14}$ | $\le3.2\times10^{-15}$ for the rechecked free-fermion cases |
| Three grand-sum routes, spread | $8.9\times10^{-16}$ | $1.3\times10^{-15}$ |
| Random-input response identity | $4.5\times10^{-16}$ | $4.5\times10^{-16}$ |

The initially reported exchange discrepancies near $10^{-10}$ at $\Delta=0$ came from the auditor's machine-precision full-space eigenvectors not being orthonormal within degenerate spaces. Rechecking with 30-digit or cluster-orthonormalized references removed that discrepancy; it was not a demonstrated code error. The exact ZZ identities above were separately checked by the auditor to about $10^{-16}$.

Its principal findings are addressed as follows:

| Finding | Current response |
|---|---|
| M1: names captured from `Global` | Explicit context-path isolation and restoration |
| M2: packed complex vectors unpacked during map averaging | Complex assignment in `mapLevelFactors` |
| M3: full initial-factor allocation per state-average cluster | Singleton/pair contractions, preserving cross-sector degeneracies |
| Tautological ZZ current check | Actual Hamiltonian/edge-magnetization commutators |
| Weak, unsaved A/B check | Full-matrix plus grand-sum comparison; separate automatic comparison export |
| Different parallel settings prevented comparison | Parallel/batch settings excluded from compatibility keys |
| Missing ZZ B reconstruction check | Added, consistent across both files |
| Caption, tick, and MaTeX-loader issues | `Exchange` label, within-decade log ticks, shared cached loader |
| Mixed-bulk memory floor | Documented beside the setting; calculated minimum included in the failure message |
| Unresolved ratios counted toward sample threshold | Only numeric ratios count; count metadata and consistent missing fraction supplied |
| Parallel execution not tested | Failed/aborted result handling strengthened; execution remains unverified |

The audit's speed measurements in §5.2 are attributed to it. Neither its reported tests nor the earlier implementation checks should be read as execution of the new files in their entirety.

### 6.3 Checks performed for this audit patch

The Wolfram evaluator returned connection/internal errors for the requested regression suite and even for a minimal `$Version`/arithmetic test. Consequently, **this patch has not received a fresh Mathematica execution check**. No new Mathematica timing, packed-array measurement, context-isolation run, parallel run, or export round-trip is claimed. The context and complex-assignment changes follow the forms the auditor explicitly tested.

Independent NumPy/SciPy calculations checked the state-average contraction against full-space energy-projector partial traces for **90 cases per regime**: $L=4,5,6$, $\Delta=0,0.45,0.5,\pm1$, $g=-0.3,0.1$, and complex pure, complex rank-three mixed, or infinite-temperature bulk factors. All 90 cases in each regime contained cross-sector degeneracies. Maximum differences were $7.19\times10^{-14}$ for ZZ and $6.03\times10^{-13}$ for exchange. The ZZ $\Gamma/4$ state and texture formula agreed with direct unitary partial traces within $2.23\times10^{-16}$. These tests verify the mathematical contraction and conventions in a separate implementation; they do not execute the Wolfram source.

The generated `.m` files passed a local balanced-delimiter/string/comment check, and no template placeholders remain. A source comparison confirms that chapter 3 differs from the latest channel revision only by the `+0.I` correction in the legacy averaging helper used by QST: the accelerated channel definitions and run cells are unchanged. The underlying Hamiltonian assembly, symmetry-sector construction, and state-evolution definitions were not modified. These source checks are narrower than a Wolfram parser or runtime test. MaTeX, real parallel subkernels, full default scans, and production-size runtimes remain outside the verified scope.
