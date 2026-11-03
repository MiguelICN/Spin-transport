# Numerical methodology: phantom XXZ bulk with two edge probes

Companion to `XXZ_Phantom_edge_dynamics.m`, revised following the independent audit of 2 November 2026. This document describes the implemented methods, distinguishes tested results from performance estimates, and records the audit response in Section 13.

The workflow computes the reduced dynamics of two additional edge spins and their complete quantum channel, with an initially pure, field-dressed bulk eigenstate. It treats ZZ and exchange boundary interactions in independent chapters. It does not perform static parameter scans, QST calculations, bulk-state optimization, or infinite-time averaging.

Its central numerical choice is **exact diagonalization followed by direct spectral evaluation at each requested time**. No differential equation is integrated. Consequently, time samples can be evaluated independently in parallel, apart from diagnostics that compare neighbouring samples or accumulate changes over a trajectory.

## 1. Model, indexing and normalization

`N` is the number of **bulk** spins. `L = N + 2` includes the two additional probes. Thus, `N -> 4` means six spins in the full calculation; `N -> 10` means twelve.

The physical tensor order is

\[
A,\ M_1,\ldots,M_N,\ B.
\]

In the formulas below, \(S_j^\alpha\) acts on bulk spin \(M_j\), while \(s_A^\alpha,s_B^\alpha\) act on the probes. Bulk sites \(j=1,N\) are physical sites \(2,L-1\), respectively.

The conventions are

\[
S^\alpha=s^\alpha=\frac{\sigma^\alpha}{2},\qquad \hbar=1,
\qquad |0\rangle=|\uparrow\rangle,\quad |1\rangle=|\downarrow\rangle.
\]

The edge basis is \(\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}\). Entropies use natural logarithms and are measured in nats.

### 1.1 Field-dressed bulk

The isolated bulk Hamiltonian used both for preparation verification and during subsequent evolution is

\[
\begin{aligned}
H_M={}&J\sum_{j=1}^{N-1}
\left(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y
+\Delta S_j^zS_{j+1}^z\right)\\
&+\frac J2S_1^x
+\frac J2\left(\cos\varphi\,S_N^x+\sin\varphi\,S_N^y\right),
\qquad \Delta=\cos\gamma.
\end{aligned}
\]

`Gamma` is the helix pitch \(\gamma\). The code derives `Delta` from it; these are not independent inputs in this workflow. `Phi` is derived from the selected helix family and boundary branch.

The two transverse fields act **only on the ends of the bulk**. They remain present during the dynamics. Their values do not depend on the additional probe states or on \(g\). There are no applied fields on the probes and no direct probe–probe interaction.

### 1.2 Probe interactions

For ZZ coupling,

\[
V_{ZZ}=g\left(s_A^zS_1^z+S_N^zs_B^z\right).
\]

For exchange coupling,

\[
V_{\mathrm{ex}}=g\left(s_A^xS_1^x+s_A^yS_1^y
+S_N^xs_B^x+S_N^ys_B^y\right).
\]

The latter gives a flip-flop matrix element \(g/2\), since

\[
s^xS^x+s^yS^y=\frac12(s^+S^-+s^-S^+).
\]

Both \(J\) and \(g\) have energy units. The interaction coefficient is \(g\), not implicitly \(gJ\). The full Hamiltonian is \(H=H_M+V\), with the bulk term embedded in the physical tensor order. The implementation requires nonzero real \(J\) and \(g\).

## 2. Initial bulk state and its independent verification

### 2.1 Analytic helices

For handedness \(\eta=\pm1\), `phantomKet` constructs

\[
|\psi_1^\eta\rangle=
\bigotimes_{j=1}^{N}
\frac{|0\rangle-e^{i\eta(j-2)\gamma}|1\rangle}{\sqrt2},
\]

or

\[
|\psi_2^\eta\rangle=
\bigotimes_{j=1}^{N}
\frac{|0\rangle+e^{i\eta j\gamma}|1\rangle}{\sqrt2}.
\]

The selected boundary angle is reduced modulo \(2\pi\):

| Family | `BoundaryBranch` | Right-field angle \(\varphi\) | Total bulk eigenenergy |
|---|---|---|---|
| 1 | `PMC` | \(\eta(N-3)\gamma\) | \((N-3)J\cos\gamma/4\) |
| 2 | `PMC` | \(\eta(N+1)\gamma\) | \((N+1)J\cos\gamma/4\) |
| 1 or 2 | `Alternative` | \(\eta(N-1)\gamma+\pi\) | \((N-1)J\cos\gamma/4\) |

These are energies of the isolated, field-dressed bulk. They do not include the interaction with the additional probes. A bulk eigenstate need not remain an eigenstate after those probes are attached at nonzero \(g\).

### 2.2 Supplied numerical states

`chooseBulkKet` accepts three sources:

| Choice | What is used |
|---|---|
| `"UploadedPhantom"` | The supplied numerical phantom vector, embedded in the `.m` file. |
| `"UploadedSpinCurMax"` | The separate supplied approximate current-maximizing vector, also embedded. |
| `"Analytic"` | The product helix generated from the chosen family, handedness and pitch. |

The supplied vectors exist for \(N=4,\ldots,10\) and require family 1, positive handedness, `PMC`, and \(\gamma=\pi/3\). They are normalized before use. The analytic constructor supports other allowed sizes and parameter choices without external state files.

`UploadedSpinCurMax` is generally **not the same product state** as `UploadedPhantom`. It occupies the same degenerate bulk energy level in the supplied examples. Its individual bond and field expectation values need not equal those of the product helix, even when their total energies coincide.

The independent audit diagonalized the projected bond-current operator within this degenerate level. It found that `UploadedSpinCurMax` approximates the positive-current maximizer for **J = +1**, with reported infidelity about 3.5×10⁻⁶–1.1×10⁻⁴ over the seven supplied sizes. It is an accurate energy eigenstate, but not the exact current-maximizing vector. These are different accuracy questions. Because the physical current includes J, changing J's sign reverses the current and interchanges its extrema. The workflow preserves the uploaded vector; it does not silently replace it with a newly optimized state.

### 2.3 Energy conventions

The code uses the spin-operator Hamiltonian written in Section 1, without an added constant energy shift. In this convention, the dimensional longitudinal coefficient is

\[
\Delta_z=J\cos\gamma,
\]

and the product-helix energy per bulk bond is \(\Delta_z/4\).

For the supplied family-1 states at \(J=1\), \(\gamma=\pi/3\), the screenshot table gives

\[
E_M=\frac{N-3}{8}.
\]

The energy numbers in the supplied `.txt` files instead equal \(-4E_M\):

| Bulk \(N\) | Full \(L\) | Energy in the implemented convention | Supplied `.txt` energy |
|---:|---:|---:|---:|
| 4 | 6 | 0.125 | −0.5 |
| 5 | 7 | 0.250 | −1.0 |
| 6 | 8 | 0.375 | −1.5 |
| 7 | 9 | 0.500 | −2.0 |
| 8 | 10 | 0.625 | −2.5 |
| 9 | 11 | 0.750 | −3.0 |
| 10 | 12 | 0.875 | −3.5 |

The verification reports both conventions and the converted diagonalized energy. It does not rescale the evolution Hamiltonian to force agreement with a differently normalized energy list. The default remains **J = 1, g = 0.1 in spin-½ units**, following the screenshot Hamiltonian. Plot captions explicitly display \(S^\alpha=\sigma^\alpha/2\), \(\hbar=1\), J, g and the dimensionless \(\Delta=\cos\gamma\).

There is a precise conditional conversion. If the generator Hamiltonian is assumed to be \(\lambda H_M(J=1)+cI\), with the same field-to-bond ratio and a size-independent shift c, the listed energies fix \(\lambda=-4\) and \(c=0\). Setting `"J" -> -4` therefore reproduces these bulk energies and the uniformly rescaled bulk Hamiltonian. The energy data alone do not establish the generator's probe interaction, current definition, or a unique Hamiltonian outside that assumed family.

For the **entire coupled system**, uniform rescaling by −4 requires both

\[
J_{\rm new}=-4J_{\rm old},\qquad g_{\rm new}=-4g_{\rm old}.
\]

Then \(U_{\rm new}(t)=U_{\rm old}(-4t)\). Changing J alone changes the relative coupling g/J and is not merely a change of time unit. If an independently specified signed probe coefficient multiplies Pauli products, \(g_P\sigma^\alpha\sigma^\alpha\), its corresponding code coefficient is \(g=4g_P\). No sign or rescaling is inferred automatically from the uploaded energy list.

### 2.4 Verification by diagonalization

`bulkEigenAudit` diagonalizes \(H_M\) and checks the original chosen vector, rather than replacing it with an arbitrary eigenvector returned by diagonalization. It records:

1. The expected table energy \(E_*\).
2. The expectation value \(\langle\psi_M|H_M|\psi_M\rangle\).
3. The residual \(\|H_M|\psi_M\rangle-E_*|\psi_M\rangle\|_2\).
4. The closest diagonalized energy.
5. The dimension of the numerically identified eigenspace around \(E_*\).
6. The projection weight in that entire eigenspace.
7. Bond-sum and boundary-field energy contributions, and overlap with the analytic product helix.

If several diagonalized levels satisfy

\[
|E_n-E_*|\le\epsilon_{\mathrm{ED}}\max(1,\|H_M\|_\infty),
\]

all their eigenvectors enter the projector. The tested projection weight is

\[
\sum_{n\text{ in this eigenspace}}
|\langle E_n|\psi_M\rangle|^2.
\]

This avoids an incorrect failure caused by arbitrary basis rotations inside degenerate eigenspaces. The residual and projection checks use scaled tolerances, with a factor of 100 as implemented in `bulkEigenAudit`.

The preliminary verification section can check both supplied vectors at every size \(N=4,\ldots,10\), and all eight combinations of family, handedness and boundary branch at the selected verification size. These checks are separate from each chapter's own check of its actual chosen bulk state. Running all-size verification is not part of the reported trajectory timing.

## 3. What is—and is not—used as an evolution sector

### 3.1 Magnetization is measured initially, but not assumed conserved

The transverse fields generally make

\[
[H_M,\sum_jS_j^z]\ne0,
\qquad
[H,s_A^z+\sum_jS_j^z+s_B^z]\ne0.
\]

`magnetizationAudit` examines Hamiltonian matrix elements connecting different computational-basis magnetizations and reports the largest such element. It also reports the initial weights

\[
w_\mu=\operatorname{Tr}(P_\mu\rho_0).
\]

**Those initial magnetization slices are not propagated independently.** An initially empty magnetization slice cannot be discarded merely because its initial weight is zero; the transverse fields may populate it later.

For the equatorial product helix, the bulk weights are binomial:

\[
w_\mu=2^{-N}\binom{N}{N/2-\mu}.
\]

Every allowed bulk magnetization has nonzero weight. The separate supplied current-maximizing vector need not have the same distribution.

### 3.2 Exact computational-basis components

`exactComponents` builds a graph whose vertices are computational-basis states. A nonzero off-diagonal Hamiltonian element connects two vertices. Connected components identify invariant blocks in this basis.

`occupiedBlocks` then:

1. Computes the input-factor weight in each component.
2. Skips components with exactly zero computed weight.
3. Retains all nonzero input columns within each occupied component.
4. Diagonalizes the retained block and projects the input factor into its eigenbasis.

There is no user-defined cutoff that drops small but nonzero occupation weights. Coherences between retained components are preserved when the evolved amplitudes are assembled.

This is a safe computational-basis connectivity reduction, **not a complete classification of every possible symmetry or conserved charge**. It does not construct a minimal Krylov subspace or automatically diagonalize additional non-diagonal symmetry operators. The old field-free magnetization/reflection/spin-flip classification is not reused without justification.

### 3.3 Exact ZZ reduction

With ZZ probe interactions, both probe \(z\) components remain conserved even though bulk magnetization does not. The full Hamiltonian separates into four conditional bulk Hamiltonians:

\[
H_{s_A,s_B}=H_M+g(s_AS_1^z+s_BS_N^z),
\qquad s_A,s_B\in\{+\tfrac12,-\tfrac12\}.
\]

The typical numerical task is therefore four complex Hermitian diagonalizations of dimension \(2^N\), rather than one of dimension \(2^{N+2}\). Any further disconnected computational components are handled by the same occupation test.

All four probe inputs are required to build an input-independent channel. A computational-basis probe preparation may use only one block for its particular state evolution, but that does not justify omitting the other channel inputs. The combined workflow prepares the channel and state together.

`FixedEdgeOccupation` distinguishes the physical input's edge weights from the inputs required for the Choi calculation. The `InputWeight` in a ZZ component report is normalized within the assigned edge-input block; four entries equal to one are not four probabilities of one in a single normalized physical state.

### 3.4 Exchange reduction

Exchange interactions do not conserve the two probe \(z\) components individually. With the transverse bulk fields, the computational-basis graph is normally connected, and the retained block has dimension

\[
D=2^{N+2}.
\]

The code uses complex Hermitian full-space ED in this case. Parallelizing time samples does not change that Hilbert-space dimension or remove the exponential cost of diagonalization.

## 4. Diagonalization and spectral propagation

### 4.1 Complex Hermitian ED

Hamiltonians are first built as sparse matrices. `complexED` checks Hermiticity relative to \(\max(1,\|H\|_\infty)\), then converts the retained block to a dense complex matrix for `Eigensystem`. Energies are sorted in ascending order; the corresponding eigenvectors become columns of \(V\).

The implementation preserves complex eigenvectors. It does not discard imaginary components or apply the legacy real-vector sign-fixing routine. Complex amplitudes are essential for the transverse \(y\) field and helical states.

The block audit checks residuals and mutual orthogonality for the first, middle and last returned eigenvectors. This is a sampled ED check, not an exhaustive residual check of every eigenvector. The separately chosen bulk ket receives the full residual and eigenspace-projection test described earlier.

### 4.2 Four input columns, one reusable engine

For a fixed normalized bulk ket, define the \(D\times4\) input matrix

\[
F_0(:,i)=|i\rangle_{AB}\otimes|\psi_M\rangle,
\]

with tensor factors reordered into physical site order. Its four columns represent four separate edge basis inputs. It is not itself a trace-one physical state factor: \(\|F_0\|_F^2=4\).

For each retained block \(\nu\), the engine stores its basis indices, active input columns, energies \(E_\nu\), eigenvectors \(V_\nu\), sparse Hamiltonian \(H_\nu\), and projected coefficients

\[
C_\nu=V_\nu^\dagger F_{0,\nu}.
\]

At time \(t\),

\[
F_\nu(t)=V_\nu\operatorname{diag}(e^{-iE_\nu t})C_\nu.
\]

The retained pieces are scattered into their physical rows and input columns to produce \(F(t)\).

`evolveInputs` evaluates a batch of times together. It concatenates the phase-weighted coefficients into a matrix and uses dense matrix multiplication for the expensive multiplication by \(V_\nu\).

The derivative is evaluated analytically as

\[
\dot F_\nu(t)=-iH_\nu F_\nu(t),
\]

using the stored sparse Hamiltonian. This removes a second dense eigenvector multiplication. A spectral-derivative fallback remains for an older engine without the stored `H` entry.

There is no ODE solver, Trotter approximation, or finite-difference time derivative. The sampling spacing controls which times are observed, not an integration step. A full propagator \(U(t)\) is not constructed in production propagation, although dense eigenvector matrices are stored.

### 4.3 Linear and logarithmic sampling

For a linear grid, `timeGrid` starts at zero and generates multiples of `dt` up to `tMax`. If the endpoint is not a multiple of the spacing, it appends the endpoint; the final interval is then shorter. Near-equal endpoints are reconciled numerically.

For a logarithmic grid, it includes zero separately and generates `LogPoints` positive, geometrically spaced times between `tMinPositive` and `tMax`. The positive samples are uniform in log time. The complete list therefore has `LogPoints + 1` samples.

Every listed time is evaluated spectrally. Finite-step channel witnesses use the actual difference between adjacent times, including the first positive logarithmic time and any shorter final linear interval.

## 5. Reduced edge-state calculation

The initial preparation supported here is

\[
\rho_0=\rho_{AB}(0)\otimes|\psi_M\rangle\langle\psi_M|.
\]

The bulk is pure and fixed. The edges may be supplied as a pure vector or a valid mixed density matrix. The current workflow does not expose the earlier general mixed-bulk or correlated full-chain input interfaces.

`stateFactor` writes the edge input as

\[
\rho_{AB}(0)=EE^\dagger,
\]

where \(E\) has four rows. For a density-matrix input it checks dimensions, Hermiticity, unit trace and positivity within tolerance, clips tolerated negative eigenvalues, and constructs a normalized eigenfactor from the positive eigenvalues.

The evolving full-chain factor is

\[
A(t)=F(t)E.
\]

To trace out the bulk, `partitionFactor` reorders and reshapes amplitudes into a matrix \(B(t)\) with four edge rows and environment/factor columns. Then

\[
\rho_{AB}(t)=B(t)B(t)^\dagger.
\]

The implementation does not propagate a full \(D\times D\) density matrix. Single-edge states are obtained by tracing one edge out of the resulting \(4\times4\) matrix.

| Recorded quantity | Numerical definition |
|---|---|
| Entropies | \(S(\rho)=-\sum_{\lambda>0}\lambda\ln\lambda\). |
| Purities | \(\operatorname{Tr}\rho^2\). |
| Mutual information | \(S(\rho_A)+S(\rho_B)-S(\rho_{AB})\). |
| Concurrence | A factor-based spin-flip calculation: QR reduction of the edge factor followed by singular values, avoiding the less stable eigenvalues of \(\rho\widetilde\rho\). |
| Negativity | Sum of negative eigenvalue magnitudes of the partial transpose. |
| Local components | \(\langle\sigma_A^\alpha\rangle,\langle\sigma_B^\alpha\rangle\); these use Pauli operators, not \(S^\alpha=\sigma^\alpha/2\). |
| Correlations | \(\langle\sigma_A^\alpha\sigma_B^\beta\rangle\), and the connected version after subtracting the product of means. |
| `InitialOverlap` | \(\operatorname{Tr}[\rho_{AB}(0)\rho_{AB}(t)]\); this is not general Uhlmann fidelity. |
| `FullOverlap` | \(\operatorname{Tr}[\rho_0\rho(t)]/\operatorname{Tr}\rho_0^2\), evaluated from factors. For a pure initial global state this is its survival probability. |
| Boundary currents | Exchange-bond \(S^z\) currents from two-site reductions. For ZZ they are set to their exact value zero because each probe's \(S^z\) commutes with the Hamiltonian. |

The default edge vector is `ConstantArray[1/2,4]`, namely the separable state \(|+\rangle_A|+\rangle_B\). Choosing an initially entangled edge state instead studies its subsequent evolution, rather than entanglement generation from separable edges.

## 6. Input-independent channel and Choi matrix

### 6.1 Definition and normalization

For the same fixed bulk state,

\[
\Lambda_t(X)=\operatorname{Tr}_M\left[U(t)
\bigl(X\otimes|\psi_M\rangle\langle\psi_M|\bigr)U(t)^\dagger\right].
\]

It is independent of the edge input selected for the reduced-state experiment. That input enters \(A(t)=F(t)E\), but not the channel construction from \(F(t)\).

The code uses the unnormalized Choi matrix

\[
C_t=\sum_{i,j=0}^{3}|i\rangle\langle j|_R
\otimes\Lambda_t(|i\rangle\langle j|),
\qquad \operatorname{Tr}C_t=4.
\]

`ChoiState` is \(\omega_t=C_t/4\), with trace one. The reference \(R\) is a mathematical device for representing the channel; no extra physical spin is introduced into the Hamiltonian or preparation.

### 6.2 Exchange channel from an amplitude factor

For environment basis state \(m\), input edge state \(i\), and output edge state \(o\), define

\[
W_{(i,o),m}(t)=\langle o,m|U(t)|i,\psi_M\rangle.
\]

`choiFactor` obtains this \(16\times2^N\) matrix by reshaping the propagated input columns. The channel is then built as

\[
C_t=W(t)W(t)^\dagger,
\qquad
\dot C_t=\dot W W^\dagger+W\dot W^\dagger.
\]

This constructs the entire channel using four propagated input kets. It avoids sixteen separately prepared matrix-unit evolutions and avoids forming full-system density matrices.

`choiToTransfer` reshuffles indices to form the \(16\times16\) matrix \(\mathsf S_t\) satisfying

\[
\operatorname{vec}_{\rm row}[\Lambda_t(X)]
=\mathsf S_t\operatorname{vec}_{\rm row}(X).
\]

The same reshuffling gives \(\dot{\mathsf S}_t\). `applyTransfer` uses this row-vectorization convention consistently.

### 6.3 ZZ channel from four conditional trajectories

Define

\[
|\phi_a(t)\rangle=e^{-iH_a t}|\psi_M\rangle,
\qquad \Gamma_{ab}(t)=\langle\phi_b(t)|\phi_a(t)\rangle.
\]

The exact channel is the Schur map

\[
\Lambda_t(X)=\Gamma(t)\circ X.
\]

`channelValue` constructs \(\Gamma\) from the four conditional kets and evaluates its derivative from their analytic derivatives. It embeds \(\Gamma\) into rows and columns `1,6,11,16` of the unnormalized Choi matrix. The other Choi entries vanish. The transfer matrix is simply

\[
\mathsf S_t=\operatorname{diag}(\operatorname{Flatten}\Gamma(t)).
\]

This structure permits several exact reductions:

- The nonzero Choi eigenvalues are the eigenvalues of \(\Gamma\); twelve structural zeros complete the Choi spectrum.
- Choi rank is at most four, and can be smaller. It is one at \(t=0\).
- Choi-state entropy is \(S(\Gamma/4)\).
- Choi partial-transpose eigenvalues are \(\Gamma_{aa}/4\) and \(\pm|\Gamma_{ab}|/4\), for \(a<b\).
- The smallest transfer singular value is \(\min_{a,b}|\Gamma_{ab}|\).

These identities replace unnecessary generic \(16\times16\) operations in the ZZ diagnostics. They remain valid with the transverse fields on the bulk, because the probes still couple diagonally in their own \(z\) basis.

**Choi rank and transfer invertibility are different properties.** The identity channel at \(t=0\) has Choi rank one but an invertible, rank-16 transfer matrix.

## 7. Channel diagnostics and numerical conditioning

### 7.1 Basic checks and channel characterizations

`channelRecord` computes Choi positivity, trace preservation, Choi normalization, numerical Kraus rank, Choi entropy and reference–output negativity. The rank threshold is

\[
\lambda_k(C_t)>
\epsilon_{\mathrm{rank}}\max(1,\lambda_{\max}(C_t)).
\]

It also records unitality, the entropy of the maximally mixed input's output, covariance deviations, the factor acting on \(|00\rangle\langle11|\), and selected signalling/retention responses.

Probe-dependent quantities use `Probes`, not the physical input chosen in `EdgeState`. `prepareChannelProbes` precomputes the input difference and Bloch contraction matrices once. Changing these probes changes the witnesses being evaluated; it does not change the channel.

For ZZ, unitality and edge \(z\)-rotation covariance are exact consequences of the Schur structure. Reflection and spin-flip covariance are measured, rather than assumed for the field-dressed preparation. The default computational-basis B probe gives zero A-to-B signalling response in ZZ, and the largest A-retention singular value remains one because the \(z\) component is conserved; neither diagnostic exhausts the channel's action on coherent inputs.

Choi negativity concerns entanglement across **reference | output**, not physical edge A | edge B.

### 7.2 Generator and canonical rates

When \(\sigma_{\min}(\mathsf S_t)\) is above `MinimumSingularValue`, the generator matrix is

\[
\mathsf L_t=\dot{\mathsf S}_t\mathsf S_t^{-1}.
\]

The generic implementation uses a right linear solve, rather than explicitly forming a matrix inverse. For ZZ, it evaluates the equivalent diagonal expression directly:

\[
\mathsf L_t=\operatorname{diag}
\left(\operatorname{Flatten}\left[\dot\Gamma(t)/\Gamma(t)\right]\right),
\]

where division is entrywise.

The canonical decoherence matrix is extracted in the orthonormal traceless two-qubit Pauli basis

\[
G_{\mu\nu}=\frac{\sigma^\mu\otimes\sigma^\nu}{2},
\qquad (\mu,\nu)\ne(0,0).
\]

The general extraction is a fixed linear transformation of `Flatten[Generator]`. `gksBasis` caches that transformation once, replacing the repeated construction and tracing of 225 operator products at every time sample. After Hermitization to remove roundoff asymmetry, its eigenvalues are the canonical rates.

For ZZ, only the traceless diagonal subspace is required. `gksDiagonalBasis` restricts the calculation to

\[
\left\{\frac{I\otimes\sigma^z}{2},
\frac{\sigma^z\otimes I}{2},
\frac{\sigma^z\otimes\sigma^z}{2}\right\}.
\]

The code diagonalizes a \(3\times3\) matrix, pads with twelve zero rates, and sorts all fifteen rates. Sorting is performed separately at each time; a plotted rate index does not track a fixed eigenoperator through crossings.

If the channel is below the invertibility threshold, the generator, canonical rates and instantaneous RHP witness are stored as `Missing`. The channel and Choi state remain defined. Small coherence factors can make logarithmic derivatives large, so sharp rate peaks must be interpreted together with transfer conditioning.

### 7.3 Instantaneous and finite-step RHP witnesses

With the operator normalization above and system dimension four, the instantaneous witness is

\[
g_{\mathrm{RHP}}(t)=\frac12\sum_k\max(0,-\gamma_k(t)).
\]

For adjacent samples \(t_{k-1},t_k\), `finishChannel` also constructs

\[
V_{k,k-1}=\Lambda_{t_k}\Lambda_{t_{k-1}}^{-1},
\]

and evaluates

\[
w_k=\max\left(0,\frac{\|C[V_{k,k-1}]\|_1}{4}-1\right),
\qquad r_k=\frac{w_k}{t_k-t_{k-1}}.
\]

The previous sample's singular value determines whether this calculation is allowed. The first sample has no preceding interval and is stored as `Missing` for these quantities.

In ZZ, the intermediate map is Schur multiplication by

\[
Q_k=\Gamma(t_k)/\Gamma(t_{k-1}).
\]

Its trace norm can be obtained from the \(4\times4\) matrix \(Q_k\), without constructing its full \(16\times16\) Choi matrix. At a fixed uniform spacing, the finite-step witness and its rate differ only by that spacing. Their agreement with the instantaneous witness is a small-spacing statement at regular times, not an identity for a finite interval near a pole.

Positive RHP witnesses diagnose failure of CP divisibility. They do not mean that \(\Lambda_t\), obtained from the original global unitary, has ceased to be a physical channel.

### 7.4 Selected-pair BLP witness

The trace-distance probe is

\[
D_k=\frac12\|\Lambda_{t_k}(\rho_1-\rho_2)\|_1.
\]

The accumulated sampled backflow is

\[
\mathcal N_k=\sum_{j=2}^{k}\max(0,D_j-D_{j-1}).
\]

This is evaluated with a single cumulative sum. It does not repeatedly recompute every preceding difference at each time, and does not optimize over input pairs. Finite sampling can miss intervening increases.

### 7.5 Channel entanglement-generation diagnostic for ZZ

The additional ZZ quantity is

\[
\mathcal C_{\mathrm{gen}}(t)=
\max_{\rho\in\mathrm{Sep}}C[\Lambda_t(\rho)]
=C[\Gamma(t)/4].
\]

For a pure product edge input with amplitudes \((a,b)\otimes(c,d)\), diagonal filtering gives output concurrence \(4|abcd|C[\Gamma/4]\). Equatorial product inputs attain the maximum, and mixtures of product states cannot improve it by convexity. Thus no numerical optimization over separable edge inputs is performed for this ZZ quantity.

This is a property of the channel with the **chosen fixed bulk state**. It is not an optimization over bulk preparations, and it should not be identified with RHP, BLP, or Choi negativity. A channel can have strong memory without generating appreciable edge entanglement.

## 8. Parallel methodology

### 8.1 Work assigned to the master and workers

The master performs bulk-state verification, Hamiltonian construction, block ED, coefficient preparation, worker initialization, final assembly, global accumulation and acceptance checks. Plotting and export also remain on the master.

Workers perform time propagation, analytic derivatives, reduced-state observables, Choi construction and per-sample channel diagnostics. Each participating kernel receives the required engine and preparation data once, rather than with every time sample. Only reduced records are returned, not the full-chain amplitude trajectory.

`DistributeDefinitions` distributes the calculation functions and their dependencies. The entire interactive context, existing output records and all embedded input-state tables are not deliberately broadcast as a single context payload.

### 8.2 Coarse groups and batch boundaries

The requested times are partitioned into batches of `BatchTimes` samples. Except for the first batch, each batch also receives the immediately preceding time sample. Workers can then evaluate finite-step RHP across the boundary using the correct adjacent pair.

After calculation, the extra leading sample is discarded from that batch's returned records. The master sorts batches by their original index, concatenates the results, and checks that the number and exact ordering of all times match the requested list. It recomputes the global BLP accumulation from the assembled distances.

Batches are assigned round-robin to coarse worker groups. `ParallelEvaluate` executes those groups concurrently. This is not a claim that the source uses `ParallelMap` with a particular scheduling option.

Each worker calls `configureWorkerThreads[]`, which requests

```mathematica
SetSystemOptions[
 "ParallelOptions" -> {"MKLThreadNumber" -> 1, "ParallelThreadNumber" -> 1}
]
```

and reads back **both** settings using `SystemOptions["ParallelOptions" -> optionName]`. `Performance["WorkerReport"]` is a list of associations containing `DefinitionsReady`, `KernelID`, and `ThreadLimit`; the latter records `Confirmed`, `MKLThreadSetting`, and `ParallelThreadSetting`. All workers must confirm both values equal 1 before trajectory evaluation. A failed confirmation prints the actual settings, clears task payloads and stops the parallel run. It does not continue under an unverified threading policy. Explicit serial execution remains available with `ParallelEnabled -> False`.

This corrects the previous top-level `"MKLThreadNumber"` option path, which was rejected in Mathematica 15.0.1. Setting and reading back the nested options was tested in the available Wolfram kernel. Actual multi-worker execution remains untested here. Readback establishes the configured options, not an operating-system measurement of every library thread. The master's threading options are not modified by production worker setup.

### 8.3 Pool size, memory budget and fallback

Both chapters currently default to:

```mathematica
"BatchTimes" -> 32,
"ParallelEnabled" -> True,
"ParallelKernels" -> 12,
"ParallelMemoryMB" -> 16384,
"ParallelVerify" -> True
```

Choose a worker count no larger than the machine's physical cores, preferably leaving one core for the master. Eight kernels can be requested instead of twelve if appropriate for the hardware; twelve is a request, not a claim about available physical cores. The code does not infer physical cores from a possibly logical processor count. The actual worker count is limited by the request, number of batches, estimated working memory and available/licensed subkernels. If fewer than two workers are usable, the calculation prints a reason and uses optimized serial execution. It never reports twelve workers merely because twelve were requested.

For task size \(B\) MiB, full dimension \(D\), and maximum batch length \(b=\texttt{BatchTimes}+1\), the pool selector estimates

\[
M_{\rm master}=4B,
\qquad
M_{\rm worker}=256+3B+
\frac{10\times16\times D\times4\times b}{2^{20}}
\quad\text{MiB}.
\]

These estimates determine a permitted worker count under `ParallelMemoryMB`. They are allocation heuristics, not an operating-system memory cap. Resident earlier runs, Mathematica overhead, result buffers, assembly copies and other applications can increase actual memory usage.

Worker task data are cleared after completion or an interrupted parallel evaluation. Existing kernels are left open for reuse.

### 8.4 Startup comparison and performance report

With `ParallelVerify -> True`, the master and the first worker independently calculate the first up to three samples. The full returned state and channel records are compared before the main parallel trajectory is accepted. Both missing values are treated as matching; numerical entries are compared by their maximum absolute difference.

This is a startup consistency check on **one worker and a few samples**, not a proof of every worker result over the entire trajectory. All returned batches subsequently undergo completeness, ordering and physical acceptance checks.

Inspect:

```mathematica
zzRun["Performance"]
exRun["Performance"]
```

The report includes actual execution mode, requested and used kernel counts, kernel IDs, memory estimates, worker setup and batch-computation times, assembly time, batch count and startup comparison error. `WorkerSetupSeconds` covers definition and data distribution; kernel launch occurs earlier and is not separately timed in that field. The startup comparison is also outside `BatchComputeSeconds`.

The enclosing `zzSeconds` or `exSeconds` measures the complete `runDynamics` call, including these overheads. Engine-setup timing, preliminary verification, reference checks, plotting and export are separate.

### 8.5 Measuring actual serial-versus-parallel performance

The optional `benchmarkDynamics` helper compares the same engine, physical edge input, full time list, tolerances and witness probes in optimized serial and parallel modes. Disabled `If[False, ...]` input cells appear after each chapter's dynamics cell; change False to True to run the corresponding benchmark:

```mathematica
zzBenchmark = benchmarkDynamics[
 zzEngine, zzEdgeStateUsed, zzTimes, zzSettingsUsed, zzProbesUsed];
Print[Dataset[{zzBenchmark}]];
```

Use the corresponding `ex...` symbols for exchange. Select an appropriate `ParallelKernels` value in the saved configuration passed to this helper. It launches the requested kernel count before either timed run; participating workers are still limited by the normal memory/batch rules. The helper can leave additional idle kernels open. It is disabled by default because it runs the trajectory twice and temporarily holds both returned trajectories in memory. It does not change production records.

The report contains `KernelLaunchSeconds`, `SerialSeconds`, `ParallelRequestedSeconds`, `ParallelExecuted`, `Speedup`, `MaximumRecordDifference`, `RecordsAgree`, `SampleCount`, and the actual `ParallelPerformance`. `RecordsAgree` compares every returned state/channel record entry, including every batch boundary, at `ValidationTolerance`. This compares the cleaned public records; the ordinary acceptance checks still use raw records. If no actual parallel execution occurred, `Speedup` is `Missing["ParallelNotExecuted"]`, never a ratio advertised as parallel acceleration.

Both run timings include normal preparation, definition/payload distribution where applicable, startup checks, record transfer, assembly and cleaning. They exclude kernel launch, prebuilt-engine ED, plotting and exports. Launching beforehand removes cold process startup from the comparison; it does not remove all first-use/cache effects. For reliable machine-specific measurements, repeat and compare medians, and inspect the reported worker count and record agreement. A separate plot/export timing is needed for whole-chapter performance.

The audit measured per-sample diagnostics of roughly 3–6 ms and serial record cleaning around 0.35–0.50 ms in its two-core environment; these are externally reported profiling results, not predicted desktop runtimes. Many small-system runs are dominated by diagnostics, communication or plotting. Dense exchange propagation becomes more important as N increases. No fixed parallel speedup is promised, and the audit's estimates of 2–5× are not measured multi-kernel speedups.

## 9. Tolerances and cleaning numerical zeros

| Setting | Current default | Role |
|---|---:|---|
| `EDTolerance` | \(10^{-10}\) | Scaled Hermiticity, ED and bulk-eigenspace checks; individual checks include the scaling factors described above. |
| `ValidationTolerance` | \(10^{-8}\) | State/channel acceptance and master–worker startup comparison. |
| `RankTolerance` | \(10^{-10}\) | Relative threshold for numerical Choi rank. |
| `MinimumSingularValue` | \(10^{-8}\) | Threshold for allowing inverse-based channel diagnostics. |
| `RecordChop` | \(10^{-12}\) | Cleans small numerical components in completed records. |
| `DenseWorkingMemoryMB` | 4096 | Guard on the estimated working memory of a single dense block diagonalization. |
| `DirectReferenceMaximumL` | 6 | Largest total length used for the automatic full-propagator reference check. |

`cleanRecords` applies `Chop` after the diagnostics and raw acceptance checks. It leaves `Time` unchanged. It does not modify sampling times, the Hamiltonian, ED energies, input occupation tests or spectral propagation. The saved `Checks` describe the unchopped calculations.

The number of displayed digits in captions is independent of numerical precision. The production matrices use machine-precision arithmetic; six-digit labels do not mean that six-digit arithmetic is used.

## 10. Scaling and practical limits

For one connected block, dense ED has quadratic eigenvector storage and approximately cubic diagonalization cost. Parallel time evaluation accelerates repeated propagation and diagnostics; it does not distribute the dense ED across the twelve subkernels.

The following storage figures assume dense machine-complex eigenvectors, at 16 bytes per entry:

| Bulk \(N\) | Full \(L\) | Typical ZZ block dimension | Typical exchange dimension | Four ZZ eigenvector matrices | Exchange eigenvector matrix |
|---:|---:|---:|---:|---:|---:|
| 4 | 6 | 16 | 64 | 0.0156 MiB | 0.0625 MiB |
| 8 | 10 | 256 | 1024 | 4 MiB | 16 MiB |
| 10 | 12 | 1024 | 4096 | 64 MiB | 256 MiB |

These are **eigenvector-array sizes, not peak memory**. Dense ED work arrays, worker copies and retained results add to them. `complexED` uses a separate conservative estimate of eight dense complex block-sized arrays before allowing a diagonalization.

The complete state and channel records are retained for all samples. Their memory grows linearly with the sample count, and temporary copies can coexist during parallel return and assembly. Increasing `tMax` while keeping the number of samples fixed is different from increasing the number of samples at fixed spacing. Very large time grids are not streamed to disk by the present implementation.

For a small system, kernel launch, data distribution and diagnostics can dominate the actual propagation. Increasing the kernel count therefore does not imply a proportional speedup. For long times, spectral phase accuracy also depends on errors in energies multiplied by time; sparse analytic propagation does not remove floating-point limitations.

## 11. Validation performed and limits of the evidence

### 11.1 Automatic checks in the workflow

Each accepted run verifies:

- The chosen isolated bulk state's table energy and eigenspace projection.
- Identity of the channel at \(t=0\).
- Reduced-state trace and positivity within tolerance.
- Choi positivity, trace and trace-preservation constraints.
- Agreement between two contractions of the same propagated input data: applying the computed channel to the chosen edge input and directly tracing the evolved state factor.
- Completeness and ordering of the assembled time samples.
- The master–worker startup comparison when parallel execution is available and that check is enabled.

These internal consistency tests share the propagation engine. Passing CP/TP or reconstruction checks alone does not establish that every Hamiltonian term or long-time phase is correct.

For total \(L\le6\) by default, `smallReferenceCheck` additionally constructs `MatrixExp[-I t H]` at \(t=0,0.137,1\). It compares evolved state factors and an explicit environment-index sum for the Choi entries with the production implementation. This is the small-system exception to the production rule of not forming \(U(t)\).

### 11.2 Development checks already completed

The supplied vectors at \(N=4,\ldots,10\) were checked against independently constructed bulk Hamiltonians and diagonalized eigenspaces. The largest supplied-state residual was approximately \(5\times10^{-15}\). Wolfram tests also checked the table branches and small-system state/channel evolution.

For the optimized update, all records were compared with the user-tested uploaded implementation at selected times, including samples around sharp channel features. In an 81-sample ZZ benchmark with \(N=4\), the previous serial calculation took approximately 4.104 seconds and the optimized serial calculation 0.386 seconds, a 10.6-fold speedup. State-record differences were below \(2\times10^{-15}\), and channel-record differences below \(2\times10^{-14}\) in that benchmark.

The optimized exchange calculation was also exercised at total \(L=10\), with a 1024-dimensional retained block. Its sampled state/channel consistency and CP/TP checks passed at approximately machine precision. This test was not a full independent propagator comparison at \(L=10\), nor a full-length production trajectory benchmark.

The development Wolfram environment could not launch parallel subkernels. The optimized serial path, batch-boundary handling and explicit unavailable-kernel fallback were exercised, but a genuine 8–12-kernel execution and speedup were **not measured there**. The startup comparison is included to detect worker inconsistencies on the user's actual installation. Plotting and export definitions were preserved from the user-tested file; the optimized update did not establish a new independent MaTeX-rendering test.

The benchmark is evidence for removing serial overhead, not a promised runtime for another machine, system size, time grid or worker count.

## 12. Reading and reusing a computed run

The two chapters have independent `zz...` and `ex...` settings. Each stores snapshots of the model, bulk choice, edge input, time settings and numerical settings used for the calculation. Editing an input cell later does not alter a stored run.

Typical data access is:

```mathematica
zzRecords[[k]]["RhoAB"]
zzMapRecords[[k]]["Choi"]
zzMapRecords[[k]]["ChoiState"]
zzMapRecords[[k]]["DephasingFactors"]
zzRun["Checks"]
zzRun["Performance"]
```

For exchange, use the corresponding `ex...` names. `DephasingFactors` is specific to the ZZ construction.

To apply an existing channel to a different valid edge input while retaining the same bulk preparation:

```mathematica
applyTransfer[zzMapRecords[[k]]["Transfer"], rhoEdge]
```

Changing the bulk state, model parameters or time samples requires a corresponding new calculation. Changing only the physical edge input does not change the already computed channel, although the file's combined run interface normally calculates state and channel records together.

Definitions and results live in the private Wolfram context **XXZPhantomDynamics**. Evaluate the opening context cell and all shared definitions before interacting with the chapter cells. The final `End[]` cell restores the caller's context and search path; it does not erase results. After leaving that context, access a result explicitly, for example:

```mathematica
XXZPhantomDynamics`zzRun["Performance"]
```

Load the updated shared definitions before using the parallel run cells. Re-evaluating only the input settings in a session containing the old definitions does not install the optimized methods.


## 13. Response to the independent audit (2 November 2026)

The audit found no physics or numerical correctness defects in the tested cases. Its independent Hamiltonian/MatrixExp comparisons support both couplings and both supplied N = 4 preparations. Its larger-size profiles and parameter checks have the narrower scope stated in the audit; they are not full large-system trajectory validation.

Changes made in this revision:

- Corrected the worker thread-option path; both configured thread limits are read back on every worker, saved in the performance report, and required before parallel computation.
- Added an optional full-trajectory serial/parallel benchmark, with explicit handling of serial fallback and no automatic duplicate production calculation.
- Made the existing spin normalization explicit in captions and clarified the relation between the supplied energies, J, g and time units. Default physical parameters and uploaded vectors are unchanged.
- Identified `UploadedSpinCurMax` as an approximate positive-current maximizer under the audit's J = +1 convention, distinct from its much more accurate energy-eigenstate property.
- Kept all existing observables, full record fields, sampling times, marker plots and exports. The audit's compact-record proposal would change the public data interface; moving cleaning before acceptance would change validation semantics unless separately refactored. Neither is required to correct the thread defect.
- Did not add a new exact-current-maximizing preparation or remove existing helper functions. Those are optional extensions, not corrections to the tested supplied-state experiment. Plot downsampling, joined curves and reduced exports were not imposed.

The audit also verified the alternative boundary-angle branch at generic pitches, both families and both handedness signs. This supports its existing implementation; it does not require changing the supplied screenshots or claiming an analytic proof from numerical tests.

Validation of this revision is recorded below. Genuine multi-kernel speedup and MaTeX rendering must still be checked on the desktop installation; neither the original developer environment nor the audit environment could launch subkernels.

- In the available Wolfram kernel, the new thread-setting helper returned `Confirmed -> True`, with both nested option readbacks equal to 1.
- Compared all returned state and channel record fields before/after this patch at nine times (including t = 35.1, 35.2, 35.3 and 333), in three-sample batches, for both regimes at N = 4. Both maximum record differences were exactly zero in that test.
- Small-system direct-reference checks passed: maximum state-factor / Choi errors were 4.14×10⁻¹⁵ / 4.33×10⁻¹⁴ for ZZ and 6.90×10⁻¹⁵ / 1.33×10⁻¹⁴ for exchange.
- Exercised the new benchmark's unavailable-subkernel path: it returned `ParallelExecuted -> False`, `Speedup -> Missing["ParallelNotExecuted"]`, and full-record agreement. This is fallback validation, not execution of the parallel branch.
- Confirmed that embedded input data are unchanged and explanatory Text cells contain only Mathematica comments. The added caption normalization uses ordinary MaTeX expressions; actual MaTeX rendering was not rerun here.
