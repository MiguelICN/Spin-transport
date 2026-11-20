# Improvement plan for `friday_beamer.tex` (today's meeting)

## Goal of today's talk

**Part I (keep as is):** Starting from product states of the edge spins, the Ising-coupled channel cannot create a Bell state, for any bulk state.

**Part II (refocus):** The reduced map $\Lambda_t$ is always a well-defined CPTP map, because the global evolution is unitary. Asking for a time-local generator, or a canonical master equation, is an *extra* requirement. It can fail, depending on the initial bulk state.

**The observation that motivates Part II.** With the phantom-helix bulk state, the canonical rates show large, isolated spikes. With the single pure eigenstates used in earlier runs (the isolated-chain ground state), no such spikes appeared.

**The hypothesis to present.** The relevant difference is not the energy or the magnetization of the bulk state. It is whether the state lies in a **degenerate (or nearly degenerate) level** of $H_E$. A two-fold degenerate level is already enough to produce the effect.

**Register.** Physical, motivated by the observation, and justified by short exact statements. There are no iff conditions, no proofs and no perturbation theory in the main slides. The detailed analysis stays in the appendix and is labelled as ongoing.

---

## 1. Diagnosis of the current deck

| # | Issue | Where | Fix |
|---|---|---|---|
| P1 | Part II opens with three slides of derivation (canonical form I, II, rates). That is heavier than today needs. | Slides 12–14 | Compress to one "result" slide; move the derivation to the appendix (§3.2) |
| P2 | The **observation** that motivates the question (helix gives spikes, earlier eigenstates did not) is never stated. The numerical slide shows only the helix and ends with "why?". | Numerical slide | New observation slide (§3.3) |
| P3 | There is no short statement of **where large rates can come from at all**. The audience needs one line that turns "large rates" into "the bulk state". | Before the hypothesis | New slide (§3.4) |
| P4 | The hypothesis slides ("Core physical argument", "Minimal illustration") use first-order degenerate perturbation theory ($W_a=P_\mathcal{D}V_aP_\mathcal{D}$, "dressing $\sim\lVert V\rVert/\Delta E$"). The audit showed it is quantitatively off for the helix runs. | Part III, slides 1–2 | Replace with an exact two-level picture (§3.5–3.6) |
| P5 | Nothing shows that the helix and the earlier eigenstates actually **differ in degeneracy**. That is the one piece of evidence the hypothesis needs, and it is cheap (§3.7). | — | New evidence slide (§3.7) |
| P6 | "Finite peaks vs genuine divergences" goes into the exact-zero question (Lemma 4, symmetry). That is beyond today's scope. | Part III, slide 3 | Move to the appendix; keep one bullet (§3.6) |
| P7 | The summary table has rows unrelated to today's message ("phase gate / $\sigma^z\sigma^z$ phase"). It also has no Part I line and no "open questions". | Summary | Rewrite (§3.8) |
| P8 | Missing caveats that someone will ask about: (i) degeneracy is necessary, not sufficient; (ii) only degeneracies the boundary field can couple count ($V_a$ conserves $M_z$); (iii) mixtures over many levels can dephase by a different mechanism. | — | One line each, on §3.5 and §3.8 |

---

## 2. New structure

| # | Slide | Action |
|---|---|---|
| 0 | Title | Small edit to the subtitle (§3.1) |
| 1–11 | **Part I: setup, conditional Hamiltonians, Case I, Bell obstruction, Case II, Schur channel, $\lvert{++}\rangle$ no-go, general no-go, caveat** | **Keep unchanged** (minor edits in §4) |
| 12 | Exact reduced map vs time-local generator | Keep, add one sentence (§3.2) |
| 13 | Canonical master equation: the result | **Merge** of the three derivation slides (§3.2) |
| 14 | Observation: helix vs earlier eigenstates | **New**, reusing the existing figure (§3.3) |
| 15 | Where can large rates come from? | **New** (§3.4) |
| 16 | Hypothesis: degeneracy of the initial bulk level | **Rewrite** of "Core physical argument" (§3.5) |
| 17 | Simplest example: a two-level bulk | **Rewrite** of "Minimal illustration" (§3.6) |
| 18 | Is the hypothesis consistent with the runs? | **New** (§3.7) |
| 19 | Summary and open questions | **Rewrite** (§3.8) |
| App. | Derivation (I), (II), rate definition, finite-vs-divergent, A.1–A.9 | Move or keep after `\appendix` |

That gives 11 + 8 main slides. Part II stays light: every slide has at most one displayed equation.

---

## 3. Slide-by-slide content

### 3.1 Title
Keep the title. Change the subtitle to:
*Partial results: Bell-state no-go and the canonical master equation — why the initial bulk state matters.*

### 3.2 Slides 12–13: the map always exists; the generator is an extra assumption
**Slide 12** (existing "Parametric dynamical map vs time-local master equation"): keep it. Add one sentence, boxed:
> Unitarity of the global evolution guarantees $\Lambda_t$. It does **not** guarantee $\Lambda_t^{-1}$, and the time-local generator needs $\Lambda_t^{-1}$.

Also add the one-line reason that $\Lambda_t$ is always CPTP: $\Gamma(t)$ is the Gram matrix of the vectors $U_a(t)\sqrt{\omega_E}$, so $\Gamma\succeq0$ with unit diagonal, for every $\omega_E$ and every $t$.

**Slide 13** (new, replaces the three derivation slides): give only the result.
- Coherences obey $\dot\rho_{ab}=\eta_{ab}(t)\rho_{ab}$, with $\eta_{ab}=\dot\Gamma_{ab}/\Gamma_{ab}$.
- The canonical form has three dephasing channels, $L_k\in\operatorname{span}\{\sigma^z_A,\sigma^z_B,\sigma^z_A\sigma^z_B\}$. Its rates are $\{\tilde\gamma_k\}=\spec\bigl(G^{\mathsf T}\eta(t)G\bigr)$, and the other 12 rates are exactly zero.
- One sentence: *the rates are built from the logarithmic derivatives $\dot\Gamma_{ab}/\Gamma_{ab}$*.
- Footnote: "Derivation: Appendix".

### 3.3 Slide 14: the observation
Left: the existing `Channel_CanonicalRates.png` (helix, $N=10$, 8-fold degenerate level).
Right, bullets:
- **Helix bulk state:** large isolated spikes ($|\gamma|\approx9.7$ at $t\approx223$ for $N=10$; up to $\approx67$ for $N=4$).
- **Earlier runs, isolated-chain ground state as bulk state:** rates stayed small and smooth.
  - If you have a saved rates or $|\Gamma_{ab}|$ plot from a ground-state run, put it side by side. The `10_dephasing_factors_Abs.png` files under `XXZ_impurities/xxz_boundary_results_a/channel_Z/` may serve, if one of them is a ground-state run.
  - Otherwise state the observation in words.
- Both bulk states are **pure** and **eigenstates** of the isolated chain.
- **Question:** what property of the initial bulk state separates the two behaviours?

### 3.4 Slide 15: where can large rates come from?
A single exact line, no derivation beyond it. For a pure bulk state, $\Gamma_{ab}(t)=\langle\phi_b(t)|\phi_a(t)\rangle$ with $|\phi_a(t)\rangle=e^{-iH_at}|\psi_0\rangle$, and
$$\frac{\dot\Gamma_{ab}}{\Gamma_{ab}}=\frac{i\,\langle\phi_b|(V_b-V_a)|\phi_a\rangle}{\langle\phi_b|\phi_a\rangle},\qquad |\text{numerator}|\le\lVert V_a-V_b\rVert\sim g .$$
Reading, in a box:
- The **numerator** is always small: at most the edge coupling, for every bulk state.
- So the rates can become large **only if the denominator becomes small**, i.e. if two conditional bulk histories become (nearly) orthogonal.
- The question becomes: **which initial bulk states can a weak boundary field rotate into an orthogonal state?**

This one identity is the mathematical justification for even expecting irregular rates, and for tying them to the initial bulk state. It holds for mixed $\omega_E$ too, with the trace in place of the overlap.

### 3.5 Slide 16: the hypothesis (rewrite "Core physical argument")
Keep the two-column layout. Remove all perturbative language ($W_a$, "zeroth order in $1/\Delta E$", "dressing by $\lVert V\rVert/\Delta E$"). Use plain physics instead:

| Isolated eigenstate | Eigenstate in a degenerate level |
|---|---|
| The nearest other level (that the field can couple to) is far away compared with the coupling $g$. | Other states sit at the *same* energy. |
| Mixing them in costs energy, so the weak field can only slightly deform $\psi_0$. Each branch mostly picks up a phase. | The field can mix them in at **no energy cost** and rotate $\psi_0$ by a large angle inside the level. |
| Histories stay nearly parallel, $|\Gamma_{ab}|\approx1$, and the rates stay small. | Opposite edge configurations push in opposite directions ($V_{\bar a}=-V_a$), so histories can become nearly orthogonal and the rates spike. |

Bottom line, boxed: *The energy and the magnetization of $\psi_0$ do not decide this. What decides it is whether $\psi_0$ has degenerate partners that the boundary field connects it to.*

Two short caveats, in a footnote:
- **Necessary, not sufficient.** If $\psi_0$ is already a state that the field does not rotate (an eigenvector of the field within the level), nothing happens.
- **What "degenerate" means here.** $V_a$ conserves $M_z$, so only degeneracies within the same $M_z$ sector count. "Degenerate" in practice means gap ≲ coupling.

### 3.6 Slide 17: simplest example, a two-level bulk (rewrite "Minimal illustration")
This answers "is a 2-dimensional degenerate sector enough?": yes. The model is exactly solvable, so it is not perturbation theory.
- **Model.** Keep only two bulk levels, $H_E=\tfrac{\Delta}{2}\sigma^z$, with boundary field $V_a=+\delta\sigma^x$ and $V_{\bar a}=-\delta\sigma^x$. Start in an eigenstate.
- **Exact result:**
$$\Gamma_{a\bar a}(t)=1-2\,\frac{\delta^2}{\Omega^2}\sin^2(\Omega t),\qquad \Omega=\sqrt{\Delta^2/4+\delta^2}.$$
- **Gapped**, $\Delta\gg\delta$: $\Gamma$ stays $\approx1$, and the rates stay small.
- **Degenerate**, $\Delta=0$: $\Gamma=\cos(2\delta t)$ reaches zero periodically, and the rates spike.
- The level energy $E_0$ never appears.
- **Figure:** $|\Gamma(t)|$ and $\kappa=-\tfrac{d}{dt}\ln|\Gamma|$ for $\Delta/\delta=20,\,4,\,0$. Write a ~20-line Python script and save it with the figure.
- **One bullet on spike height.** If $\psi_0$ sits unevenly inside the degenerate level, $|\Gamma|$ dips to a small nonzero value instead of zero, and the spike is tall but finite. That is what the helix data show ($\min|\Gamma_{01,10}|\approx10^{-4}$). Whether exact zeros (true divergences) occur is an open question: see the appendix.
- **One bullet on $d>2$.** Several frequencies interfere, and the regular comb becomes irregular, with occasional deep dips. The mechanism is unchanged.

(I checked numerically: the formula above agrees with direct evaluation to $10^{-15}$. In the degenerate case with weight imbalance $c=\langle\sigma^x\rangle$, the spike height is $\delta(1-c^2)/|c|$. Keep that formula off the main slide, or as a footnote at most.)

### 3.7 Slide 18: is the hypothesis consistent with the runs?
This is the evidence slide. It is a comparison of numbers, not a derivation. Show a two-row table of gap$/\lVert V\rVert$, with $\lVert V\rVert=(g_A+g_B)/4=0.05$:

| Bulk state | Degeneracy of its level | Gap to nearest coupled level | gap$/\lVert V\rVert$ | Spikes? |
|---|---|---|---|---|
| Isolated-chain ground state, $L=6$–$10$, $\Delta=0.5,1$ | 1 (within its $M_z$ sector) | $0.33$–$0.56$ | **7–11** | No |
| Phantom helix, $N=4$–$10$ | $d=2$–$8$ | $0.002$–$0.13$ | **0.04–2.7** | Yes |

- The helix row is from appendix A.1.
- I computed the ground-state row today for an open XXZ chain with $J=1$, $S=\sigma/2$, no bulk boundary fields, taking the gap within the ground state's $M_z$ sector. **Check it against the conventions of the actual earlier runs before presenting it.**
- Note for odd $L$: the ground state is 2-fold degenerate ($M_z=\pm\tfrac12$). The field cannot connect the two partners because it conserves $M_z$, so the state behaves as isolated. This illustrates "only coupled degeneracies count".

Optional extra bullet: for $N=4,5$ ($d=2$), the observed spike spacing ($\approx70$–$76$) matches the rotation frequency set by the coupling inside the level (appendix A.3). Label it as an estimate.

Bottom line: *The two families differ precisely in degeneracy, and the spikes follow the degeneracy.*

### 3.8 Slide 19: summary and open questions
**Part I:** No edge state with $C<1$, in particular no product state, reaches $|\Phi^\pm\rangle$ or $|\Psi^\pm\rangle$, for any bulk state.

**Part II, as a logical chain:**
1. Global unitarity ⇒ $\Lambda_t$ is CPTP and smooth for every bulk state and every $t$.
2. A canonical master equation needs $\Lambda_t^{-1}$, i.e. $\Gamma_{ab}(t)\neq0$.
3. The rates are $\dot\Gamma/\Gamma$, and the numerator is always $\lesssim g$. Large rates therefore require near-orthogonal conditional bulk histories.
4. Hypothesis: that can happen only if the initial bulk state has degenerate partners the boundary field couples to. Isolated eigenstates are protected.
5. Consistent with the data: ground states (gap/‖V‖ ≈ 7–11) show no spikes; the helix (gap/‖V‖ ≲ 3) does.

**Open questions (ongoing work):**
- When are the spikes true divergences (exact zeros) rather than finite near-zeros?
- Quantitative theory for $d>2$.
- Mixed bulk states over many levels, which can dephase by phase averaging (a different mechanism).
- A controlled test (below).

Drop the current table, or cut it to the rows "conditional histories / $|\Gamma|$ / rates / depends on $E_0$ or $M_z$?".

---

## 4. Minor edits to Part I (content unchanged)
- Setup slide, outline bullet: rephrase item (2) as "canonical master equation, and why the initial bulk state decides whether its rates stay regular".
- "Case II: $|{++}\rangle$": in the parenthetical "Exceptions…", keep it as is; it is already corrected.
- "Caveat: the no-go is basis-specific": its last sub-bullet ("approximately CZ for an isolated non-degenerate bulk eigenstate") ties in with Part II. Add "(cf. Part II: isolated eigenstates keep $|\Gamma|\approx1$)". Otherwise unchanged.

## 5. Moved to the appendix (unchanged)
- "Derivation of the canonical form (I)" and "(II)"
- "Canonical master equation and definition of the canonical rates"
- "When do the canonical rates blow up?" (its key point now lives in §3.4)
- "Finite peaks vs genuine divergences"
- A.1–A.9 as they are

## 6. Optional, if time allows: a controlled test
Run the existing notebook with the **same** $N$, $g$ and $H_E$, and two bulk eigenstates:
- (a) one that is isolated within its $M_z$ sector;
- (b) one in a degenerate level, with similar energy or $\langle M_z\rangle$ if possible.

Show $|\Gamma_{ab}(t)|$ and the rates side by side. This would turn "consistent with" into a direct test that varies only the degeneracy. It is not needed for today: §3.7 already carries the argument.

## 7. Order of work (if time is short)
1. §3.3, §3.4 and §3.5 (LaTeX only, ~45 min). These carry the message.
2. §3.7 table: verify the ground-state gaps against the actual run conventions (~15 min).
3. §3.6 two-level figure (~15 min).
4. Reshuffle into the appendix, and write §3.8 (~15 min).
