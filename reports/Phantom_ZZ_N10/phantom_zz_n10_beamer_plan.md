# Plan: Beamer presentation of the `Phantom_ZZ_N10` run

This plan is for an agent who will write a new Beamer deck explaining **every dynamical quantity** in this folder. Read it all before starting. Section 1 summarises the analysis of the run that the plan is based on. Sections 2–6 are the instructions.

---

## 0. Deliverable and ground rules

- **File:** `reports/phantom_zz_n10_beamer.tex`. Do not modify `friday_beamer.tex`.
- **Preamble:** copy the preamble of `reports/friday_beamer.tex`: document class, colours, macros `\ket`, `\bra`, `\Tr`, `\Heff`, `\dmat`, `\diag`, `\spec`. Use the same notation as that deck.
- **Figures:** use relative paths. Add `\graphicspath{{./}{Phantom_ZZ_N10/}}` and include files as e.g. `\includegraphics{Phantom_ZZ_N10/Channel_Gamma_Abs.png}`.
- **Build** from `reports/`:
  1. `pdflatex -synctex=1 -interaction=nonstopmode -halt-on-error phantom_zz_n10_beamer.tex`
  2. `bibtex phantom_zz_n10_beamer`
  3. Run `pdflatex` twice more.
  4. Finish with **zero `Overfull \vbox`** warnings.
- **Slide-density rules.** These are hard limits. The user explicitly asked not to overload slides.
  - At most **one main figure per slide**. Two figures are allowed only if both are small and directly compared.
  - At most **4 bullets** per slide, each ≤ 2 lines.
  - Define each quantity **on the slide where it is first plotted**, with exactly one displayed formula. Never collect many definitions on one slide.
  - Notation slides have no figures. Result slides have a figure, a one-line definition, and ≤ 3 observations.
  - Do not use `shrink` to make content fit. Split the slide instead.
- **Do not invent numbers.** Every number on a slide must come from Section 1 of this plan, from `settings.txt`, or from data you extract yourself (Section 5). Mark eyeballed values with "≈".
- **References:** reuse `references.bib` (`Hall2014`, `Wootters1998`). Add a new entry only if you have verified it.

---

## 1. Analysis of the run (facts the deck relies on)

### 1.1 Parameters (from `settings.txt` and the plot headers)

| Quantity | Value |
|---|---|
| Bulk size | N = 10 bulk spins, L = 12 with the two edges |
| Couplings | J = 1, γ = π/3, Δ = cos γ = 1/2 |
| Coupling type and strength | ZZ, g = 0.1 |
| Initial bulk state | phantom helix Ψ₁⁺: family 1, handedness +1, branch `PMC`, φ = π/3, state `UploadedPhantom` |
| Bulk energy | E = 7/8 = bonds 9/8 + left field −1/8 + right field −1/8 |
| Degeneracy of the bulk level | 8 (verified: eigenstate residual 3·10⁻¹⁶, weight in eigenspace 1) |
| Initial edge state | \|++⟩, ρ_AB(0) = all entries 1/4 |
| Time grid | t ∈ [0, 500], δt = 0.01, 50001 samples |
| BLP pair | ρ₁ = \|++⟩⟨++\| and ρ₂ = \|+−⟩⟨+−\| (`Probes`) |
| Signalling probe on B | ρ_B = \|1⟩⟨1\| |
| Invertibility threshold | `MinimumSingularValue` = 10⁻⁸, the cutoff for computing the generator |

Numerical checks from `settings.txt`:

| Check | Value |
|---|---|
| State trace error | 7·10⁻¹⁵ |
| Min state eigenvalue | −8·10⁻¹⁶ |
| Choi TP error | 5·10⁻¹⁴ |
| Choi trace error | 3·10⁻¹⁴ |
| Min Choi eigenvalue | −3·10⁻¹⁵ |

### 1.2 The model, in the code's conventions

The source of truth is `xxz_helix_states/XXZ_Phantom_edge_dynamics_methodology.md`, §1–7.

- **Spins and tensor order.** Spin operators are S^α = s^α = σ^α/2, with ħ = 1, |0⟩ = |↑⟩ and |1⟩ = |↓⟩. The tensor order is A, M₁, …, M_N, B.
- **Bulk Hamiltonian, with transverse boundary fields:**
  H_M = J Σ_j (S^x_j S^x_{j+1} + S^y_j S^y_{j+1} + Δ S^z_j S^z_{j+1}) + (J/2) S^x_1 + (J/2)(cos φ S^x_N + sin φ S^y_N).
  **These fields break bulk M_z conservation.** The setting `MagnetizationConservedExactly → False`. The helix has binomial weight over all M_z sectors. The deck must state this. It differs from the Friday deck, which wrote the fields as "optional".
- **ZZ coupling:** V = g(s^z_A S^z_1 + S^z_N s^z_B). This gives the conditional Hamiltonians H_a = H_M + g(s_A S^z_1 + s_B S^z_N), with ‖V_a‖ = g/2 = 0.05.
- **Conditional bulk states and the channel:** |φ_a(t)⟩ = e^{−iH_a t}|ψ_M⟩, Γ_ab(t) = ⟨φ_b|φ_a⟩, and Λ_t(X) = Γ(t)∘X.

### 1.3 Key structural fact: with the |++⟩ input, the state *is* the channel

Because ρ_AB(0) has every entry equal to 1/4:

  **ρ_AB(t) = Γ(t)/4.**

The nonzero part of the Choi state is the same matrix: the Choi eigenvalues are the eigenvalues of Γ, plus 12 zeros. So almost every plotted quantity, State_* and Channel_* alike, is a function of the **six complex numbers Γ_ab(t)**, a < b. Exact identities (from the code, `channelRecord`, and from Γ/4 = ρ):

| Plot | Exact expression in terms of Γ |
|---|---|
| `Channel_MapEntropy` | S(Γ/4), i.e. **identical to** `State_EntropyAB` |
| `Channel_EntanglingConcurrence` | C[Γ/4], i.e. **identical to** `State_Concurrence` |
| `Channel_MinDephasingEigenvalue` | λ_min(Γ) = 4·`State_MinEigenvalue` |
| `Channel_ChoiEigenvalues` | eig(Γ) (4 curves) plus 12 zeros, i.e. 4·`State_EdgeDensityEigenvalues` |
| `Channel_MinSingularValue` | min_{a<b} \|Γ_ab\| |
| `Channel_ChoiPTMinimum` | −max_{a<b} \|Γ_ab\|/4 |
| `Channel_ChoiNegativity` | ¼ Σ_{a<b} \|Γ_ab\|, an ℓ₁-coherence (1.5 at t = 0) |
| `Channel_DoubleQuantum*` | Γ_{00,11} |
| `State_EdgeCoherences_*` | Γ_ab/4 |
| `State_LocalAX`, `State_LocalAY` | ⟨σ^x_A⟩ = ½ Re(Γ_{00,10} + Γ_{01,11}) and ⟨σ^y_A⟩ = −½ Im(Γ_{00,10} + Γ_{01,11}); likewise for B with (Γ_{00,01}, Γ_{10,11}). **Verify the sign of the y-component numerically.** |
| `State_Purity` | Tr ρ² = (4 + 2Σ_{a<b}\|Γ_ab\|²)/16 |

**Make one slide of this.** It is the key to reading all the plots and removes about 40 % of them as duplicates.

### 1.4 Quantities that are trivial because both edge s^z are conserved

The following are exactly constant. The deck should state each value and its one-line reason, not show the plot.

| Plot | Exact value | Reason |
|---|---|---|
| `State_EdgePopulations` | 1/4 each | populations frozen |
| `State_MagnetizationA/B`, `State_LocalAZ/BZ` | 0 | populations frozen |
| `State_CorrelationZZ`, `State_ConnectedZZ` | 0 | populations frozen |
| `State_CurrentLeft/Right` | 0 | [H, s^z_A] = [H, s^z_B] = 0 |
| `Channel_SignalAtoB` | 0 | the probe on B is \|1⟩⟨1\| |
| `Channel_RetentionA` | 1 | z-component of A is conserved |
| `Channel_UnitalityError` | 0 | Schur map |
| `Channel_MixedInputEntropy` | ln 4 ≈ 1.386 | Schur map |
| `Channel_U1CovarianceError` | 0 | Schur map |
| `Channel_DephasingPopulationError` | 0 | Γ_aa = 1 |
| `Channel_KrausRank` | 4 for t > 0 (1 at t = 0) | rank Γ ≤ 4 |

### 1.5 Pure numerical-consistency checks (flat at 0)

`State_TraceError`, `State_HermiticityError`, `Channel_TPError`, `Channel_ChoiTraceError`, `Channel_MinChoiEigenvalue`.

**Problem:** these are plotted on Mathematica's default ±1 axis, which only shows "≪ 1". Do not include these PNGs. Put their maxima from `settings.txt` (§1.1) in a table instead. If you can read `results.wxf`, put max over t of each error there too (§5).

### 1.6 Symmetries visible in the data

Observed in the plots:
- Im Γ_{00,11} ≡ 0, so the double-quantum factor is real.
- \|ρ_A\| = \|ρ_B\|, ⟨σ^x_A⟩ = ⟨σ^x_B⟩ and ⟨σ^y_A⟩ = −⟨σ^y_B⟩.
- In Re Γ and \|Γ\|, the curves of Γ_{00,01} and Γ_{00,10} lie on top of Γ_{01,11} and Γ_{10,11}. In Im Γ the same pairs are sign mirrors.
- Correlations: ZX = −XZ, YZ = ZY, and ConnectedYX ≈ −ConnectedXY.
- **RF covariance is exact** (`RFCovarianceError` = 0). Reflection alone and spin-flip alone are **broken**: both reach ≈ 0.79, and their two curves look identical.

Interpretation, to verify before presenting:
- The Friday deck appendix A.5 has two antiunitary symmetries of the helix preparation, Θ = (∏σ^x)K (spin flip) and T = e^{−iφM_z} ℛ K (reflection).
- Their product ΘT is **unitary**. It maps V_a ↦ V_{σ(ā)}, where σ swaps A and B and the bar flips both edge spins. That gives Γ_ab(t) = Γ_{σ(ā)σ(b̄)}(t), which is exactly ℛℱ-covariance.
- Reality of Γ_{00,11} is Friday A.5, Proposition 3.
- Reflection and flip alone each come with complex conjugation and time reversal, so they are not channel covariances.

Present this on **one** slide titled "Symmetries of Γ(t)". Before writing it, confirm the identity Γ_ab = Γ_{σ(ā)σ(b̄)} numerically from the data, or from the plots if the data is unavailable. If it fails, present only the observations.

**Consequence for the deck:** only **3 independent complex functions** remain:
- Γ_{01,10}, the "flip-flop" pair;
- Γ_{00,11}, real, the "double-quantum" pair;
- one of the four "single-flip" factors, Γ_{00,01} for example (the other three follow by symmetry).

**Verify** this count from the data before stating it.

### 1.7 The sharp feature at t ≈ 222.6, and what does *not* show it

This is the central physics message. It is the same as the Friday slide "Why only the rates spike".

**Spikes at t ≈ 222.6:**
- `Channel_CanonicalRates`: γ₁₅ ≈ +9.5, γ₁ ≈ −9.7.
- `Channel_MinimumRate`: ≈ −9.7.
- `Channel_RHPInstantaneous`: ≈ 4.8.
- `Channel_RHPFiniteStepRate`: ≈ 4.9.
- `Channel_RHPFiniteStep`: ≈ 0.049, which is the rate × δt.
- Smaller bumps of ≈ ±0.3–0.4 at t ≈ 77, 361 and 374.

**Dips to ≈ 0, smooth, at the same times:**
- `Channel_MinSingularValue`.
- \|Γ_{01,10}\| (`Channel_Gamma_Abs`, `State_EdgeCoherences_Abs`). min\|Γ_{01,10}\| = 1.07·10⁻⁴ at t = 222.60, from `friday_beamer_math_audit.md`. **Re-check this from the data.** The value is above the 10⁻⁸ threshold, so the rates were computed there.
- Near-zeros at t ≈ 75, 360 and 373 match the small rate bumps.

**Smooth through t ≈ 222.6:**
- Entropies, purities and mutual information.
- Concurrence and negativity (zero there).
- All correlations and local Bloch components.
- Choi spectrum, Choi negativity and map entropy.
- BLP distance and accumulated BLP. The accumulated BLP has no step at 222; it sits inside a slow rise over 207–297.
- Initial and full overlaps.
- Covariance errors.

Message: **only quantities built from Λ_t⁻¹ (the generator, its rates, the RHP witnesses) inherit the spike.** Every quantity that is a function of the state or the map is smooth.

### 1.8 Other notable features

- **Entanglement window.**
  - C(ρ_AB) is exactly 0 until t ≈ 279. It peaks at ≈ 0.117 at t ≈ 298 and returns to 0 at t ≈ 345, with a small bump ≈ 0.004 at t ≈ 478–492.
  - Negativity has the same shape, peaking at ≈ 0.058.
  - The peak coincides with Re Γ_{01,10} ≈ −0.82 at t ≈ 300, the purity bump (0.46), the dip in S_AB (0.88), and the peak in I(A:B) (0.49).
  - Link to Part I of the Friday deck: C = C[Γ/4] is the separable-input entanglement capacity, and no Bell state is reached.
- **Recurring times.** t ≈ 297–300 and t ≈ 435 are extrema in many plots.
- **Kinks are artefacts of sorting or minimisation, not physics.** Examples: `State_MinEigenvalue` at 76, 252, 311, 413 and 444; `Channel_MinSingularValue`; `Channel_MinDephasingEigenvalue`; the onset and offset of the concurrence. They come from taking a minimum or a max(0, ·) over crossing branches. Say so once.
- **Unsorted curves.** `State_EdgeDensityEigenvalues` shows λ₃ and λ₄ crossing, so these curves are not sorted. Mention it, or re-plot them sorted.
- **Overlaps.**
  - `State_InitialOverlap` = Tr[ρ_AB(0)ρ_AB(t)] falls 1 → ≈ 0.41 by t ≈ 100, then to ≈ 0.26 at t ≈ 420, then rises to ≈ 0.44.
  - `State_FullOverlap` is the survival probability of the full initial state |++⟩⊗|ψ_M⟩. It falls to ≈ 0.06 at t ≈ 440.
  - The full overlap equals (1/16)\|Σ_a ⟨ψ_M\|φ_a(t)⟩\|², a **Loschmidt-type echo** of the bulk. This connects to the Friday Loschmidt-echo slide. **Verify this formula** before using it.
- **Non-Markovianity summary.** BLP for the chosen pair shows backflow, with accumulated 𝒩 ≈ 1.12 by t ≈ 447. The steps lie at t ≈ 92–125, 207–297 and 375–447. These are not at the rate spike: BLP sees broad recoherence, while RHP sees the near-singularity. Both witness non-Markovianity, and they measure different things.

### 1.9 Label problems in the existing PNGs (decide: fix or footnote)

1. Edge operators are labelled σ₁ and σ_L. Readers may take these as chain sites. Add a footnote: "σ₁ ≡ edge A, σ_L ≡ edge B".
2. `State_EdgePopulations` uses an unexplained subscript "ss".
3. Hidden curves: 2 of the 6 curves are invisible in `Gamma_Re`, `Gamma_Abs` and `EdgeCoherences_Re/Abs`. State this explicitly, and cite the symmetry slide.
4. `CurrentLeft` and `CurrentRight` define the current with opposite signs. These plots are dropped anyway.
5. The y-axis of `State_InitialOverlap` starts at 0.2.
6. The plots never define S_t, ℛ, ℱ, M_{A→B} or ρ₁,₂. Define them on the slides.

---

## 2. Deck structure (target 26–30 main slides plus an appendix)

### Part 0 — Opening (2 slides)

| # | Title | Content |
|---|---|---|
| 0.1 | Title | "Phantom-helix bulk, ZZ edge coupling: full dynamical characterization (N = 10)" |
| 0.2 | Outline and one-sentence message | The four parts below. Message: *everything is a function of Γ(t); only inverse-based quantities spike.* |

### Part A — Notation and general framework (9 slides, no data plots)

One concept per slide. No figures.

| # | Title | Content |
|---|---|---|
| A1 | System and conventions | Tensor order A, M₁…M_N, B. S = σ/2, ħ = 1, \|0⟩ = \|↑⟩. Edge basis (00, 01, 10, 11). Entropies in nats. Local plots use Pauli σ, not S. |
| A2 | Bulk Hamiltonian with boundary fields | The H_M formula from §1.2. Say the transverse end fields make the helix an exact eigenstate **and** break M_z conservation. |
| A3 | Initial bulk state: the phantom helix | The product formula ψ₁^η = ⊗_j (\|0⟩ − e^{iη(j−2)γ}\|1⟩)/√2. Table: N = 10, γ = π/3, φ = π/3, E = 7/8 (with its contributions), degeneracy 8, verification residual. One bullet: binomial M_z weights. |
| A4 | ZZ coupling and conditional Hamiltonians | V, the conserved edge s^z, U_t = Σ_a \|a⟩⟨a\| ⊗ U_a(t), the table of H_a, and ‖V_a‖ = 0.05. |
| A5 | Conditional bulk states and the dephasing matrix | \|φ_a(t)⟩, Γ_ab = ⟨φ_b\|φ_a⟩. Properties: Gram matrix, Γ ⪰ 0, Γ_aa = 1, \|Γ_ab\| ≤ 1. |
| A6 | The exact channel and its representations | Λ_t(X) = Γ∘X. Choi matrix C_t (Tr = 4) and Choi state C_t/4. Transfer matrix S_t = diag(vec Γ) (define it here, since the plots use it). Kraus rank = rank Γ. |
| A7 | With the \|++⟩ input, state = channel | ρ_AB(t) = Γ/4. Point forward to the dictionary slide. |
| A8 | Time-local generator | Needs invertibility: σ_min(S_t) = min\|Γ_ab\|. L_t = Ṡ S⁻¹ = diag(Γ̇/Γ), η_ab = d/dt ln Γ_ab = −κ + iω. Note the code's 10⁻⁸ threshold. |
| A9 | Canonical master equation and rates | Basis {σ^z⊗I, I⊗σ^z, σ^z⊗σ^z}/2. d = Gᵀ η G, rates = spec d plus 12 zeros, sorted at every t. **Caveat:** sorting means a rate index does not follow a fixed operator through crossings. Bound: \|γ_k\| ≤ ‖η‖. |
| A10 | Non-Markovianity witnesses | RHP instantaneous g = ½ Σ max(0, −γ_k). Finite step: w_k = max(0, ‖C[Λ_{t_k}Λ_{t_{k−1}}⁻¹]‖₁/4 − 1), with rate r_k = w_k/δt. BLP: D(t) = ½‖Λ_t(ρ₁ − ρ₂)‖₁, accumulated 𝒩 = Σ max(0, ΔD), with ρ₁ = \|++⟩ and ρ₂ = \|+−⟩. |

A10 may be split into two slides (RHP | BLP) if it exceeds the density rules.

### Part B — Reading guide (3 slides)

| # | Title | Content |
|---|---|---|
| B1 | Dictionary: every quantity is a function of Γ | The table of §1.3, reduced to ≤ 8 rows. Put the rest in the appendix. |
| B2 | What is exactly constant, and why | The table of §1.4. No plots. Footnote: these PNGs exist but are omitted. |
| B3 | Numerical checks | The table of §1.5, with numbers. No plots. |

### Part C — Results (one quantity or theme per slide, 13–15 slides)

Each slide: figure on the left (~58 % width); on the right, a one-line definition plus ≤ 3 observations taken from §1.6–1.8.

| # | Title | Figure | Key observations |
|---|---|---|---|
| C1 | Conditional bulk overlaps \|Γ_ab(t)\| | `Channel_Gamma_Abs.png` | Γ_{01,10} explores the most, with near-zeros at ≈ 75, 222, 360 and 373. Γ_{00,11} stays in 0.48–0.95. Two curves are hidden by symmetry. |
| C2 | Phases of Γ: Re and Im | `Channel_Gamma_Re.png` and `Channel_Gamma_Im.png`, side by side, small | Im Γ_{00,11} ≡ 0. Re Γ_{01,10} reaches −0.82 at t ≈ 300 (a sign change). |
| C3 | Symmetries of Γ(t) | `Channel_RFCovarianceError.png` or `Channel_ReflectionCovarianceError.png` | ℛℱ exact; ℛ and ℱ alone broken. Explanation per §1.6, then "3 independent factors". |
| C4 | Edge coherences | `State_EdgeCoherences_Abs.png` | Equal to Γ/4; populations frozen at 1/4. |
| C5 | Single-edge dynamics | `State_LocalAX.png` (optionally `State_LocalAY.png`) | Formulas from §1.3. A = B for x, A = −B for y. Each edge dephases towards its z-axis, so ρ_A → nearly I/2. |
| C6 | Entropies and purity | `State_EntropyAB.png` (with `State_EntropyA.png` small) | S_A → ln 2 by t ≈ 125. S_AB ≈ 1.1–1.29, dipping to 0.88 at t ≈ 297. Equals the map entropy (dictionary). |
| C7 | Edge mutual information | `State_MutualInfo.png` | Peaks at 172, 297 and 462. Correlations build up through the bulk with no direct A–B term. |
| C8 | Two-point correlations | `State_CorrelationXX.png` and `State_CorrelationYY.png` | YY reaches −0.83 at t ≈ 300. Relations ZX = −XZ, YZ = ZY. Mixed z-correlations are nonzero even though ⟨σ^z⟩ = 0; connected = raw for z. |
| C9 | Entanglement between the edges | `State_Concurrence.png` | Zero until 279, max 0.117 at 298, zero after 345. It equals C[Γ/4], the maximum over separable inputs. No Bell state (link to the Friday Part I). Negativity has the same shape. |
| C10 | Return probabilities | `State_FullOverlap.png` (with `State_InitialOverlap.png` small) | Definitions. Full overlap = Loschmidt-type echo (after verifying the formula). |
| C11 | Choi spectrum and map entropy | `Channel_ChoiEigenvalues.png` | 4 nonzero eigenvalues = eig Γ; Kraus rank 4. Map entropy = S_AB. |
| C12 | Reference–output negativity | `Channel_ChoiNegativity.png` | Equals ¼Σ\|Γ_ab\|. It measures how far the map is from fully dephasing; it is not A–B entanglement. |
| C13 | Invertibility of Λ_t | `Channel_MinSingularValue.png` | σ_min = min\|Γ_ab\|. Near-zeros at 75, 222, 360 and 373. **Add a log-scale zoom** (§5) showing 1.07·10⁻⁴ at 222.60. |
| C14 | Canonical rates | `Channel_CanonicalRates.png` | Only γ₁ and γ₁₅ move. Spike ±9.6 at 222.6; small bumps at 77, 361 and 374. Mirror shape. |
| C15 | RHP witnesses | `Channel_RHPInstantaneous.png` (with `Channel_RHPFiniteStep.png` small) | Instantaneous ≈ 4.8 and finite-step rate ≈ 4.9 agree at δt = 0.01. CP-divisibility is violated mainly at the spike. |
| C16 | BLP witness | `Channel_BLPDistance.png` and `Channel_BLPAccumulated.png` | Broad backflow episodes; no step at 222.6. BLP and RHP detect different things. |

### Part D — Synthesis (3 slides)

| # | Title | Content |
|---|---|---|
| D1 | What spikes and what does not at t ≈ 222.6 | A two-column table from §1.7: spikes vs smooth. One sentence on why (ratio vs linear quantity; Friday slide "Why only the rates spike"). |
| D2 | Zoom figure | **New figure** (§5): t ∈ [200, 245], three stacked panels sharing the t-axis: (i) log₁₀\|Γ_{01,10}\|; (ii) γ₁ and γ₁₅; (iii) S_AB, C and BLP D(t). The visual proof of D1. |
| D3 | Summary and open questions | 4 bullets: Γ is the whole story; conservation kills the z-sector; only inverse-based quantities spike; degeneracy hypothesis (Friday deck) as the origin. Open questions: exact zeros vs near-zeros; comparison run with a non-degenerate bulk state. |

### Appendix

| # | Title | Content |
|---|---|---|
| App-1 | Full plot inventory | A table of all 83 PNGs: filename → slide used, *or* reason omitted ("duplicate of …", "exactly constant", "numerical check"). Use §1.3–1.5 and §3. |
| App-2 | Remaining duplicates | e.g. `State_EntropyB`, `State_PurityB`, `State_LocalBX/BY`, `State_Negativity`, `State_ConnectedXZ`, `Channel_FlipCovarianceError`, `Channel_DoubleQuantum*`, `Channel_MinimumRate`, `Channel_RHPFiniteStepRate`: a small grid of thumbnails, at most 4 per slide. |
| App-3 | Derivation of d = Gᵀ η G | Optional: copy from `friday_beamer.tex` (frames "Appendix: Derivation of the canonical form (I)/(II)"). |

---

## 3. Plot-to-slide map (all 83 PNGs)

**Used in Part C:**
- `Channel_`: `Gamma_Abs`, `Gamma_Re`, `Gamma_Im`, `RFCovarianceError` or `ReflectionCovarianceError`, `ChoiEigenvalues`, `ChoiNegativity`, `MinSingularValue`, `CanonicalRates`, `RHPInstantaneous`, `RHPFiniteStep`, `BLPDistance`, `BLPAccumulated`.
- `State_`: `EdgeCoherences_Abs`, `LocalAX`, `LocalAY`, `EntropyAB`, `EntropyA`, `MutualInfo`, `CorrelationXX`, `CorrelationYY`, `Concurrence`, `FullOverlap`, `InitialOverlap`.

**Duplicates** (appendix, with a reason):
- `State_EntropyB`, `State_PurityB`, `State_LocalBX`, `State_LocalBY`, `State_Negativity`.
- `State_LocalCoherences_Re/Im/Abs`.
- All `State_Connected*` involving z, plus `State_ConnectedYY` and `State_ConnectedYX`.
- `State_CorrelationZX`, `State_CorrelationZY`.
- `Channel_MapEntropy`, `Channel_EntanglingConcurrence`, `Channel_MinDephasingEigenvalue`, `Channel_ChoiPTMinimum`, `Channel_DoubleQuantumReal/Imag/Abs`.
- `Channel_FlipCovarianceError`, `Channel_MinimumRate`, `Channel_RHPFiniteStepRate`.
- `State_EdgeCoherences_Re/Im`, `State_EdgeDensityEigenvalues`, `State_MinEigenvalue`.

**Optional in Part C, if space allows:** `State_Purity`, `State_PurityA`, `State_CorrelationXY`, `State_CorrelationXZ`, `State_CorrelationYZ`, `State_ConnectedXY`, `State_ConnectedXX`.

**Exactly constant** (B2 table only): see §1.4.

**Numerical checks** (B3 table only): see §1.5.

---

## 4. Wording rules (consistency with the Friday deck)

- Say **"conditional bulk states"** for |φ_a(t)⟩, never "histories" or "trajectories".
- Say **"near-zero"** and **"tall finite spike"** for t ≈ 222.6. Do not call it a "divergence": min\|Γ\| = 10⁻⁴ ≠ 0.
- The canonical-rate claims must not contradict `friday_beamer.tex` slides "Canonical master equation: the result" and "Why only the rates spike". Reuse their formulas.
- Write "Λ_t is CPTP at all times" wherever non-Markovianity is discussed.

---

## 5. Additional work: data extraction and new figures

These need Mathematica (`wolframscript`) to read `results.wxf`, which is 400 MB. If Mathematica is unavailable, skip §5, keep the eyeballed "≈" values, and say so in the deck footnote and in your final report.

1. **Load the data.** `r = Import["reports/Phantom_ZZ_N10/results.wxf"]`. Inspect `Keys[r]` and the record structure: records carry keys like "Time", "DephasingFactors", "CanonicalRates", "MinSingularValue", …; see `channelRecord` in `xxz_helix_states/XXZ_Phantom_edge_dynamics.m`.
2. **Extract these numbers:**
   - max over t of every §1.5 check.
   - The location and value of min\|Γ_{01,10}\| over [200, 245], expected 222.60 and 1.07·10⁻⁴.
   - Peak γ₁ and γ₁₅.
   - Concurrence onset, peak and offset.
3. **Verify these identities numerically** to 10⁻¹⁰:
   - ρ_AB = Γ/4.
   - The §1.3 table.
   - The §1.6 symmetry Γ_ab = Γ_{σ(ā)σ(b̄)}, and Im Γ_{00,11} = 0.
   - The full-overlap echo formula in §1.8.

   Report any failure to the user instead of presenting the claim.
4. **New figures**, saved as PNG in `reports/Phantom_ZZ_N10/new/`, with the plotting script saved alongside:
   - `zoom_222.png`: the three stacked panels of slide D2.
   - `minsv_log.png`: log₁₀ σ_min(S_t) over [0, 500], for slide C13.
   - Optional: `independent_gammas.png`, showing only the 3 independent Γ's (\|·\| and arg).

---

## 6. Acceptance checklist (the agent must confirm each)

- [ ] Compiles with no errors and no Overfull vbox. No frame uses `shrink`.
- [ ] Every result slide has ≤ 1 main figure and ≤ 4 bullets, and defines its quantity in one formula.
- [ ] Part A has no figures, and every symbol used later is defined there or on its slide.
- [ ] S_t, ℛ, ℱ, M_{A→B}, ρ₁, ρ₂ and ρ_B are all defined.
- [ ] Every one of the 83 PNGs appears in the App-1 inventory, either with the slide that uses it or with a reason it is omitted.
- [ ] The transverse boundary fields and the non-conservation of M_z are stated.
- [ ] Synthesis slide D1 is present and agrees with §1.7.
- [ ] All numbers either come from §1 or are re-extracted, and eyeballed values carry "≈".
- [ ] Unverified interpretations (the §1.6 symmetry explanation, the §1.8 echo formula) are labelled as such if they were not checked in §5.
- [ ] `friday_beamer.tex` is untouched.
