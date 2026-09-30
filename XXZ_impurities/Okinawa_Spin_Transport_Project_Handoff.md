# Okinawa — Spin Transport: project handoff

**Status date:** 22 September 2026.  **Scope:** This note records the discussions, documents, code, and results identified in this project and its earlier boundary-spin discussions. It distinguishes a result reported by the user from an implementation present in a file and from a proposed experiment. It does not report new simulations or verify unpublished novelty.

## 1. Research question and its development

The broad question is how much can be learned about an interacting spin chain by preparing and measuring one or two controllable spins at its ends. The usual geometry has an open spin-1/2 XXZ chain of **L chain sites** with separate system spins **A** and **B**, ordered `A,1,...,L,B`. The chief numerical objects are the two-spin reduced state `rho_AB(t)` and, when the chain preparation is fixed and the two-spin input is varied, the complete dynamical map `E_t` acting on A and B.

The user's questions have included: which properties of the chain and its transport can be inferred from endpoints alone; what `rho_AB(t)` and the full channel add beyond a single spin's coherence; whether the end spins can act as source and detector; how integrability-preserving versus generic impurity couplings affect an observable experiment; and whether open-chain phantom/helix states can be identified from edge-only time traces. These are distinct research questions. No unique transport coefficient or phantom-state witness has yet been established by the numerical outputs reviewed here.

An earlier [three-month proposal](#documents-and-code-to-open-first) narrowed one possible collaboration to two comparisons: (i) full single-qubit process tomography versus one optimized relaxation/noise observable for inferring a defined chain property, and (ii) two probes used as a source and detector of nonlocal propagation. It proposed an XX/free-fermion large-chain benchmark followed by interacting XXZ with MPS, with boundary integrability as an optional backaction comparison. **This was a proposal, not a completed numerical program.**

## 2. Physical model and conventions in the current Mathematica scripts

The current `EdgeDensity` and `EdgeChannel` packages state the following conventions (with `hbar=1`, chain spin `S=σ/2`, and the probes' Pauli matrices denoted `σ_A,σ_B`):

```text
H_C = (Jxy/4) Σ_{j=1}^{L-1} (X_j X_{j+1} + Y_j Y_{j+1} + Δ Z_j Z_{j+1})
      + (w/2) Σ_j Z_j + (ed/2) Z_DefectSite
H_Q = (ω_A/2) Z_A + (ω_B/2) Z_B
H_int,Z = (g_A/4) Z_A Z_1 + (g_B/4) Z_L Z_B
H_int,exchange = g_A(σ_A^+ S_1^- + σ_A^- S_1^+)
               + g_B(S_L^+ σ_B^- + S_L^- σ_B^+).
```

Here `Δ=Jz/Jxy`; it is entered as a real number, not as a scaled integer. `ed` is a **Hamiltonian** field defect, separate from a defect inserted into an initial helix. The scripts document an old-to-new coupling conversion `gExchange = 4 lambdaOld` and `omega = 2 DeltaProbeOld`; check conventions before comparing early reports with recent runs.

For a fixed chain state `rho_C(0)` and an initially factorized system–environment state, the exact reduced map is

```text
E_t(X) = Tr_C[ exp(-i H t) (X ⊗ rho_C(0)) exp(+i H t) ],
rho_AB(t) = E_t(rho_AB(0)).
```

The probes A and B may be initially entangled with **each other**, while the scripts assume `rho_AB(0) ⊗ rho_C(0)`. The chain input can be a pure vector or a density matrix. The numerical construction does not impose a thermal bath, weak coupling, Markov approximation, or master equation. Conversely, a separate chain preparation (for example, another helix pitch) defines a **different** `E_t`.

The `Z` option is longitudinal pure dephasing: probe computational-basis populations are conserved. It is implemented through four conditional chain Hamiltonians and their probe-coherence factors. The exchange option transfers excitation between each probe and the adjacent chain spin; its map is built from full-system evolution. Neither option, by default, adds the tuned boundary fields required to make a finite open-chain phantom helix an exact eigenstate. The ordinary `Z` and exchange Hamiltonians are also not the same as the integrability-preserving boundary family discussed in the impurity paper.

## 3. Completed earlier calculations and reported observations

The user-supplied **`spin project.pdf`**, dated 24 August 2026, documents a prior exact-spectral-decomposition study of two flip-flop-coupled end spins with an open XXZ chain. Its stated preparation is a Bell pair on A,B and a thermal chain state, with inverse temperature `β`; the zero-temperature and infinite-temperature limits are compared. It defines and displays reduced-state purity, entropies, mutual information, concurrence, pair overlap/survival, global initial-state overlap, and a Choi-state description of the complete two-qubit map. It describes sweeps over `L=3,...,7`, `Δ` in approximately `[0.5,1.5]`, several `β` values, and a separate weaker-coupling set. The report contains example plots at `L=7` for `β=∞` and `β=0`, plus Choi-purity and Choi-eigenvalue figures.

The report says that the example anisotropy sweeps did not substantially change the displayed qualitative behavior. It discusses a loss of A–B concurrence while correlations with the chain develop. These are **observations and interpretations in that report**; they are not a newly reanalyzed dataset here. In particular, for a mixed thermal global input, reduced-state purity is not by itself a pure-state entanglement measure, and the report's global `Tr[rho(0)rho(t)]` is an initial-state overlap rather than the usual pure-state survival probability. The report makes the latter distinction explicitly.

There are more recent exported example images in the project folder, including `01_populations.png`, `01_populations(1).png`, `04_entanglement.png`, `08_density_validation.png`, `01_choi_purity_and_identity_fidelity.png`, and `07_channel_validation.png`. Their existence establishes that results were plotted for reduced states and channels; these images alone do not establish a transport regime or a successful edge-only helix classifier. For example, the readable labels on `01_populations(1).png` specify the exchange model, `L=6`, `Δ=0.01`, `g_A=g_B=0.1`, a Néel chain, `|++>` probes, and `t=0,...,200` at `dt=0.1`.

The user also reported a later **actual helix-related Z-model run** at `L=8` (10 total spins), `t=0,...,1000`, `dt=0.1` (10,001 samples), using 12 parallel time chunks/kernels. Runs for `Δ=0.01, 0.5, 1` completed and reported a zero parallel-versus-serial endpoint matrix error; `Δ=1.5` and `10` failed the `PhantomHelixKet` range check `-1 <= Δ <= 1`. The reported time-evolution timings for the three completed anisotropies were approximately 15.69 s, 4.26 s, and 4.22 s, respectively. This verifies a completed computation and the expected invalid-input failure, **not** a measured helix lifetime or a finding about transport. The result directory reported in the conversation used `density_Z/run02/L8_Delta...`; current source defaults and run labels should be checked independently before reproducing that run.

## 4. Implemented numerical infrastructure

Two independent Wolfram Language packages in the `Okinawa - Spin transport` Library folder cover the currently active calculations:

| File | Implemented purpose | Points to check when using it |
| --- | --- | --- |
| `xxz_edge_density_matrix_MATEX_parallel_v3(1).m` | Fixed-input `rho_AB(t)`, `rho_A(t)`, `rho_B(t)` for `Z` or exchange coupling; exported trajectories of populations, six coherences, entropy, mutual information, concurrence, negativity, purity, initial overlap, local magnetizations, connected **probe** correlations, and numerical checks. Includes helix constructors in its `Z` section. | The saved `Z` section has five successive `Clear[ChainStateZ]` definitions; **only the last, central-phase-defect definition, is active**. Its saved `deltaListZ` still includes `1.5` and `10`, which are invalid for that active matched-helix constructor. Its saved exchange section remains Néel by default. |
| `xxz_edge_dynamical_map_MATEX_parallel_v3.m` | Complete two-probe channel for a fixed chain preparation, as a normalized `16×16` Choi state and a superoperator; includes `ChoiZ`, `ChoiExchange`, `ApplyChannelZ`, and `ApplyChannelExchange`. Exports Choi purity/entropy/eigenvalues, process mutual information, distance from a product channel, cross-probe Pauli response, unitality, and CPTP checks. | This saved file has default Néel chain states; it has **not** been shown to contain the later helix subsection. Its saved defaults differ from the density file's, including chain length and some probe frequencies. Match **every** physical and initial-state parameter before comparing channels with density trajectories. |

Earlier versions `xxz_edge_density_matrix_MATEX_v2.m` and `xxz_edge_dynamical_map_MATEX_v2.m` remain available. The v3 scripts diagonalize by conserved total-Z sectors, reuse each eigensystem, divide time samples among parallel workers, and check selected parallel results against serial endpoints. Their metadata and selected trajectories can be written to WXF; PNG plots are raw sampled trajectories, not smoothing-based fits. A separate earlier `spin_channels_statevector.m` and several notebooks (`spin channels testing.nb`, `spin channels clean.nb`, `hamiltonian.nb`, `hamiltonian _xxz.nb`, `choi_thermal_sweep.nb`) are also stored; this handoff does not claim that every notebook has been brought into the current v3 workflow.

The initial chain `NeelKet`, `GroundKet` (isolated-chain ground state), arbitrary product inputs, and pure or mixed chain preparations are implemented as alternatives. The density package exposes an independently editable `ProbeStateZ` and `ProbeStateExchange`; its saved defaults use `|++>_A,B`. The channel package does not require choosing a probe input to compute the entire map. In the density package, the currently saved default parameters include `L=8`, `Jxy=1`, `w=ed=0`, `ω_A=ω_B=1`, `g_A=g_B=0.1`, `tMax=1000`, and `dt=0.1` for both coupling sections. **Defaults in a source file are not proof that each combination was executed.**

The channel script has an adjacent-time CP-divisibility diagnostic when intermediate maps are well-conditioned; a singular or ill-conditioned intermediate map is recorded as unresolved. Choi purity is a map diagnostic; it is neither A–B concurrence nor, on its own, a proof of non-Markovianity. Connected `ZZ` and `XX` values from the density script are correlations **between the probes**, not direct measurements of `⟨S_1^z(t) S_L^z(0)⟩` in the chain.

## 5. Phantom Bethe states and spin helices: discussion and code

The project discussed phantom Bethe roots, the difference between fixed-magnetization Bethe states and factorized spin-helix superpositions, chiral currents, periodic-chain commensurability, open-chain boundary conditions, and whether these states should be described as scars in different settings. The working easy-plane product input is

```text
|helix(q,θ,φ)⟩ = ⊗_{j=1}^L [ cos(θ/2)|0⟩_j
                  + exp(i[φ+(j-1)q]) sin(θ/2)|1⟩_j ],
q = ± arccos(Δ) for a matched transverse phantom pitch, |Δ| <= 1.
```

`SpinHelixKet` implements an arbitrary real pitch and polar angle. `PhantomHelixKet` selects matched `q` and either chirality. `HelixCentralDefectKet` starts from a matched transverse helix, adds a specified phase to one interior spin, and leaves the first and last **chain** spins initially identical to those of the unmodified helix. The discussed variants are `q=±arccos Δ`, `θ=π/3`, fixed `q=π/3` while `Δ` changes, and a central phase shift `π/2`. A fixed-pitch real helix can be initialized at arbitrary real `Δ`, but the matched constructor intentionally rejects `|Δ|>1`.

The user's **specific edge-only question** is whether dynamics visible only on the boundary spins can distinguish a phantom/helix preparation from ordinary or deliberately defected chain states. Static endpoint preparations alone cannot certify the entire bulk state: the central-defect construction is an explicit control with initially unchanged endpoint chain spins. A proposed comparison therefore uses time-dependent single-edge and joint-edge signals, pitch detuning, both chiralities, a tilted helix, and defects. Suggested quantities include endpoint coherences, `rho_AB(t)`, joint-versus-marginal trace distance, connected correlations, and onset times. **No validated edge-only certification threshold, exponent, or length scaling has yet been reported.**

The important physics constraint is that `Δ=cos q` is the bulk matching condition. An ordinary finite open XXZ chain without suitable end fields has boundary-induced evolution even for a matched nontrivial helix; appropriate open-chain boundary terms can make special helices exact eigenstates. The separately attached probes and the script's `g_A,g_B` terms further change the full Hamiltonian. Thus neither the saved density script nor the proposed exchange subsection makes the entire `A+chain+B` preparation an exact phantom eigenstate automatically. A uniform `w` field, a local Hamiltonian `ed` defect, and an initial phase defect are different controls.

In the immediately preceding discussion, a **paste-in** phantom/helix subsection was supplied for `ChainStateExchange` with one active central-defect choice and commented alternatives, plus a suggested valid list `deltaListExchange={-0.5,0,0.5}` and separate run labels. That subsection was **provided in chat**; it has not been shown to have been inserted into, run in, or saved as a new version of the existing exchange file. Likewise the suggested first test `L=12, Δ=0.5, q≈π/3, θ=π/2` with pitch/chirality/defect comparisons remains a proposed test.

## 6. Boundary impurity, integrability, and transport discussions

The paper **“Edge modes and boundary impurities in the anisotropic Heisenberg spin chain,”** *Phys. Rev. B* **111**, 174430 (2025), motivated questions about the gapped antiferromagnetic XXZ chain, fractional `±1/4` edge magnetization, spin-1/2 impurities, integrable and generic boundary attachments, Kondo screening, bound modes, and a nonintegrable midgap regime. The discussion covered the paper's use of Bethe ansatz, exact diagonalization, and DMRG, and its integrable boundary parameter `d_q`. The current purely transverse flip-flop probe coupling should **not** be assumed to lie in that interacting model's integrable impurity family. The user's question about an integrability-breaking edge perturbation producing prethermal behavior was discussed as a possible research direction; no corresponding time-scale law has been measured in this project.

The user asked whether Bethe integrability makes `rho_A(t)` or `rho_AB(t)` at `L=50` as directly tractable as a Gaussian transverse-field Ising evolution at much larger `L`. The discussion concluded that Bethe-solvable spectra and thermodynamics do not, by themselves, yield an inexpensive generic real-time, two-probe reduced density matrix: full interacting nonequilibrium dynamics may still require MPS/TDVP or other numerical methods. No `L=50` two-probe calculation was reported.

The observable-level questions include edge magnetizations, endpoint relaxation/coherence, response after perturbing A, signal arrival at B, frequency-dependent transmission, connected probe correlations, and whether finite-size trends support ballistic, diffusive, or anomalous propagation. Source–detector response and direct calculation of relevant chain correlators were proposed to make any inferred property checkable. In the exchange model, probe polarization dynamics can be related to local exchanged magnetization under appropriate continuity conditions; in the pure-`Z` coupling, longitudinal probe populations are conserved, so this is a different sensing protocol. No extraction of a diffusion constant or transport exponent from the present small-system outputs is documented.

Two single-qubit XXZ sensing works were compared with this approach: **“Harnessing spin-qubit decoherence to probe strongly interacting quantum systems,”** arXiv:2410.22003 and *Phys. Rev. B* **111**, L161115 (2025), and **“Thermal robustness of sensing quantum phases of matter via qubit probes,”** arXiv:2609.12964. These papers concern information encoded in a **single** qubit coupled to an XXZ chain, including anisotropy/phase information and thermal robustness. The discussions distinguished exact conditional-unitary reduced dynamics from weak-coupling analytical or master-equation approximations. The project's two **edge** probes were raised as a way to ask about spatially separated information and the full joint channel; that extra information has not yet been quantified against a single-probe baseline in reported project results.

## 7. Other papers discussed, with their status in this project

| Source | Role in the discussion | Project-specific result status |
| --- | --- | --- |
| `2102.03295v2.pdf`, `2102.03299v2.pdf` (and later copies) | Phantom Bethe roots, periodic and open XXZ chiral states, special open-boundary criteria and semi-phantom states. | Studied as theory background; only the product-state initializers and limited runs above are documented project computations. |
| `s41567-022-01651-7.pdf`, Nature Physics (2022) | Experimental long-lived helices and anisotropy sensitivity; finite open-end effects; discussion of higher-dimensional/nonintegrable analogues. | No project experiment or scar-spectrum analysis reported. |
| `2206.06771v3.pdf` | Proposed boundary-dissipative preparation of phantom Bethe states. | The current scripts implement **unitary** dynamics, not this Lindblad preparation. |
| `PhysRevX.14.041070.pdf`, “Engineering Hierarchical Symmetries” | Recently discussed Floquet/prethermal symmetry ladder `SU(2) → U(1) → Z2 → E` and whether reduced edge channels could reveal successive regimes. | Paper analysis and possible application only; no driven edge-probe implementation or result reported. |
| `Two Qubits at the Ends of an XXZ Chain Mediated Entanglement, Non-Markovianity, and Phantom-Helix Probes in Exact Diagonalization.pdf` | An externally generated research-ideas review collected in the project folder. | Treat its feasibility, novelty, and proposed directions as **suggestions**, not independently verified claims or completed work. |

## 8. Documents and code to open first

For a new agent, start with these exact names; most recent scripts and papers are in the Library folder **`Okinawa - Spin transport`**. The earlier report and impurity paper are stored as **`spin project.pdf`** and **`PhysRevB.111.174430 (1).pdf`** outside that project folder.

1. `spin project.pdf` — earlier physical setup, thermal Bell-pair computations, example results, and Choi formulation. Pay attention to which quantities were actually calculated there.
2. `project.pdf` / `project(4).tex` — the September three-month sensing/source–detector **proposal**, including constraints on measurable targets and proposed XX and MPS stages.
3. `xxz_edge_density_matrix_MATEX_parallel_v3(1).m` — most recent identified reduced-state script with helix definitions and current `Z`/exchange defaults.
4. `xxz_edge_dynamical_map_MATEX_parallel_v3.m` — separate full-channel implementation; check its defaults separately.
5. The named exported plots and their accompanying `parameters.txt`, `initial_states.wxf`, and saved trajectory WXF **when available in the user's local results directories**; the project folder contains some PNG copies, not a complete verified inventory of the underlying WXF data.
6. The XXZ impurity, central-probe, and phantom/helix papers identified above when the associated physical claim matters.

## 9. Unresolved points and evidential limits

- Decide which physical property is to be inferred, and compare endpoint measurements with the **same property's independently computed chain reference**. The full channel is available numerically, but its advantage over an optimized single observable has not been established.
- Analyze actual stored helix and non-helix trajectories with matching `L`, Hamiltonian, probe preparation, coupling, time grid, and output metadata before claiming sensitivity to chirality, pitch, or a remote defect.
- Separate ordinary-open-chain boundary decay from any exact tuned-boundary helix claim; such tuned boundaries are not in the default Hamiltonian. A matched product chain input is not automatically an eigenstate of the probe-coupled system.
- Treat the draft proposals for source–detector transmission, integrable-boundary comparisons, prethermalization, non-Markovianity, and hierarchical Floquet symmetries as **unperformed** unless later files or run records show otherwise.
- Avoid interpreting a finite `L=6` or `L=8` recurrence or a Choi-purity trend as a thermodynamic transport law, a diffusion coefficient, or a standalone proof of information backflow.

**Practical status:** There is a defined exact two-probe model; an earlier report includes thermal Bell-pair and Choi results for short chains; newer Mathematica reduced-state and complete-channel programs plus example plots exist; an `L=8` Z-coupled matched-helix run has been reported at three allowed anisotropies; and multiple physically distinct research tests have been proposed. There is no documented completed finite-size transport classification or validated edge-only phantom-state witness in the material reviewed for this handoff.
