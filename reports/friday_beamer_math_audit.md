# Mathematical audit of `friday_beamer.tex`

I checked every displayed equation and every mathematical claim in [friday_beamer.tex](friday_beamer.tex), slide by slide. I also tested the claims about the phantom-helix runs with an independent NumPy recomputation that uses the conventions of `XXZ_Phantom_edge_dynamics.m`. That recomputation reproduces the plotted canonical rates exactly: $\gamma_{\min}=-9.05$ and $\gamma_{\max}=9.11$ at $t=222.6$ for $N=10$; $-65.8$ and $+58.5$ at $t=1582.7$ for $N=4$; $\pm2.06$ for $N=5$.

Notation follows the slides. $a\in\{00,01,10,11\}$, $\bar a$ is the opposite configuration ($\overline{00}=11$, $\overline{01}=10$), $V_a=g_As^a_AS^z_1+g_Bs^a_BS^z_N$, $H_a=H_E+V_a$, $U_a(t)=e^{-iH_at}$, and
$\Gamma_{ab}(t)=\operatorname{Tr}_E[\omega_E U_b^\dagger(t)U_a(t)]$.

---

## 0. Summary

| # | Slide (frame title) | Statement | Verdict |
|---|---|---|---|
| E1 | *Case 1: Pure initial bulk state* | "non-degenerate … $g\ll\Delta E$ (here $g=0.1\ll J=1$)" ⇒ $\lvert\Gamma_{ab}\rvert=1-O(g/\Delta E)$ | **Error.** The hypothesis is checked against the wrong scale. The relevant gap is to the nearest other level of $H_E$, and in the actual runs it is $0.0020$–$0.133$, so it is *not* $\gg\lVert V_a\rVert=g/2$. Corrected rigorous bound in §1. |
| E2 | *Case 2: degenerate sector* | "**Exact** reduction … $\Gamma^{\mathcal D}_{a\bar a}=\sum_m p_m e^{-2i\delta_m t}$" | **Error.** This is first-order degenerate perturbation theory only. The second-order term is even in $V$ and spoils the reduction. Numerically, the deviation from exact $\Gamma$ is $0.27$–$0.52$. See §2. |
| E3 | *Minimal mechanism ($d=2$)* | $\eta^{\mathcal D}_{a\bar a}=-\Delta\delta\tan(\Delta\delta t)=\frac1{t-t_n}+O(t-t_n)$ | **Error.** The term $-i(\delta_1+\delta_2)$ is missing, so the remainder is $O(1)$, not $O(t-t_n)$. See §3. |
| E4 | *Numerical rates*; *$d=2$*; *Comparison* | The rates "diverge"; $d=2$ is "sufficient" for poles; $N=4,5$ "exhibit periodic rate poles" | **Error.** Beyond first order, the exact spin-flip symmetry gives only $\Gamma_{a\bar a}(t)=\Gamma_{a\bar a}(-t)$, not reality, so exact real-time zeros of $\Gamma_{01,10}$ are non-generic (codimension 2). In every run ($N=4,\dots,10$), $\min_t\lvert\Gamma_{ab}\rvert\ge1.07\times10^{-4}>0$, and the peaks are **finite**: resolved maxima are $66.7$, $2.0$ and $9.7$. See §4. |
| E5 | *Comparison: $d=2$ vs $d>2$* | Zero condition for $d\ge3$ is "polygon inequality $\max p_m\le\frac12$ (open set)"; "$d>2$ removes the fine-tuning" | **Error.** The polygon inequality is necessary but not sufficient. For generic frequencies an exact zero has measure zero. Corrected statement in §5. |
| E6 | *Derivation of the canonical form (I)* | "Because $\mathcal L_t(E_{ab})=\eta_{ab}E_{ab}$ never mixes entries, $H_{\rm eff}$ is diagonal and $\mathsf d$ is supported on $F_{1,2,3}$" | **Logical gap.** The conclusion is true, but it does not follow from the stated reason alone. A rigorous proof is in §6. |
| m1 | *Case II: Starting in $\lvert{++}\rangle$* | "For $a\ne b$, $U_b^\dagger U_a\ne I_E$" | Minor. False if $g_A=0$ or $g_B=0$, and $\ne I$ does not imply damping. |
| m2 | *Comparison* | "$t\approx224$" | Minor. The deepest near-zero is at $t=222.60$ ($\lvert\Gamma_{01,10}\rvert=1.07\times10^{-4}$). |
| m3 | *Parametric map vs. time-local* | "depends **entirely** on $\omega_E$" | Overstated. It also depends on $H_E$, $g_A$ and $g_B$. |
| m4 | *Canonical rates* | "$\gamma_1\ge0$ iff CP-divisible" | Correct but degenerate: with 12 padded zeros, $\gamma_1=\min(0,\tilde\gamma_{\min})\le0$ always. |
| m5 | *General Schur no-go* | local-filter optimality "\cite{Wootters1998}" | Attribution. Wootters defines $C$; the optimality result is the report's Lemma (Sec. 3.9). |
| R | *Why $\lvert{++}\rangle$ cannot evolve into a Bell state* | Bell fidelity $\le\frac12$ | **Correct**, but the bound is basis-specific. It does **not** exclude maximal entanglement. Remark in §7. |

**Checked and correct.** The conditional-Hamiltonian table, $U_t=\sum_aP_a\otimes U_a$, both Case I slides, the Case II global and reduced states, the Schur form, $\Gamma\succeq0$, $\Gamma_{aa}=1$ and $\lvert\Gamma_{ab}\rvert\le1$, the Bell fidelities $(1\pm\operatorname{Re}\Gamma)/4$, the general no-go ($r=\lvert\Phi_\theta\rangle\langle\Phi_\theta\rvert$ forced), $\det\Lambda_t=\prod_{a<b}\lvert\Gamma_{ab}\rvert^2$, Step 1, the matrix $G$, the elementwise and matrix forms of $\eta$, the inversion $\mathsf d=G^{\mathsf T}\eta G$ and $h_a=-\frac14\sum_b\operatorname{Im}\eta_{ab}$, the canonical master equation, the rate definition (it matches `diagonalCanonicalRates`), the Rayleigh bound, the local pole $\kappa\simeq-r/(t-t_*)$, the $d=2$ formula for $\Gamma$ and its zero condition, and the $\Theta$ pairing within first order.

---

## 1. E1: the weak-coupling hypothesis for a non-degenerate eigenstate

**What the slide says.** If $\lvert\psi_0\rangle=\lvert E_0\rangle$ is non-degenerate with gap $\Delta E\gg g$, then $\lvert\Gamma_{ab}(t)\rvert=1-O(g/\Delta E)$ for all $t$. The slide adds "here $g=0.1\ll J=1$".

**The mistake.** The bound is controlled by $v:=\lVert V_a\rVert$ divided by the **distance from $E_0$ to the rest of $\operatorname{spec}(H_E)$**, not by $g/J$. For spin-$\tfrac12$,
$$\lVert V_a\rVert=\tfrac{g_A}{4}+\tfrac{g_B}{4}=\tfrac g2=0.05 .$$
In the phantom runs, the distance from the helix level to the nearest other level of $H_E$ is:

| $N$ | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| $d=\dim\mathcal D$ | 2 | 2 | 4 | 4 | 4 | 8 | 8 |
| gap $\Delta E$ | 0.133 | 0.098 | 0.018 | 0.048 | 0.0031 | 0.0067 | 0.0020 |

So $v/\Delta E$ ranges from $0.38$ to $25$. The regime "$g\ll\Delta E$" is **never** satisfied in these runs. The same failure also invalidates the first-order treatment of Case 2 (§2).

**Correct, rigorous statement.**

> **Lemma 1.** Let $H_E\lvert E_0\rangle=E_0\lvert E_0\rangle$ be non-degenerate with $\Delta E=\operatorname{dist}(E_0,\operatorname{spec}H_E\setminus\{E_0\})$. Let $v_a=\lVert V_a\rVert<\Delta E/2$ and $\varepsilon_a:=v_a/(\Delta E-v_a)$. Then for all $t\in\mathbb R$,
> $$\lvert\Gamma_{ab}(t)\rvert\;\ge\;\cos\theta_a\cos\theta_b\cos(\theta_a+\theta_b)-\varepsilon_a-\varepsilon_b-\varepsilon_a\varepsilon_b,\qquad \sin\theta_a\le\varepsilon_a .$$

*Proof.*

1. **One eigenvalue near $E_0$.** By Weyl's inequality, $H_a$ has exactly one eigenvalue $E_a$ with $\lvert E_a-E_0\rvert\le v_a$, with eigenvector $\lvert E_0^a\rangle$. Every other eigenvalue $\lambda$ satisfies $\lvert\lambda-E_0\rvert\ge\Delta E-v_a$.
2. **Overlap with the perturbed eigenvector.** Write $\lvert E_0\rangle=c_a\lvert E_0^a\rangle+\lvert r_a\rangle$ with $r_a\perp E_0^a$. Then
$$(H_a-E_0)\lvert E_0\rangle=c_a(E_a-E_0)\lvert E_0^a\rangle+(H_a-E_0)\lvert r_a\rangle .$$
The two terms are orthogonal, the left side equals $V_a\lvert E_0\rangle$ with norm $\le v_a$, and $\lVert(H_a-E_0)r_a\rVert\ge(\Delta E-v_a)\lVert r_a\rVert$. Hence $\lVert r_a\rVert\le\varepsilon_a$.
3. **Evolved state.** It follows that $\lvert\phi_a(t)\rangle=c_ae^{-iE_at}\lvert E_0^a\rangle+e^{-iH_at}\lvert r_a\rangle$, with $\lvert c_a\rvert=\cos\theta_a$ and $\sin\theta_a=\lVert r_a\rVert$.
4. **Bound on the overlap.** Expand $\langle\phi_b\vert\phi_a\rangle$ and bound the three cross terms by Cauchy–Schwarz ($\le\varepsilon_a$, $\le\varepsilon_b$, $\le\varepsilon_a\varepsilon_b$). For the main term, use $\lvert\langle E_0^b\vert E_0^a\rangle\rvert\ge\cos(\theta_a+\theta_b)$, the triangle inequality for angles between rays. $\square$

For $\varepsilon_a=\varepsilon_b=\varepsilon$, the right-hand side is positive for $\varepsilon\lesssim0.3$, i.e. for $\Delta E\gtrsim4.3\,v=0.22$. **None** of the computed systems satisfies this. The qualitative conclusion, that a well-isolated non-degenerate eigenstate gives no zeros, remains true, but the slide must state the correct condition and must not claim that the phantom runs are in that regime.

---

## 2. E2: the degenerate-sector reduction is first order, not exact

**Claim on the slide.** With $W_a=P_{\mathcal D}V_aP_{\mathcal D}$, "Exact reduction":
$\Gamma^{\mathcal D}_{a\bar a}(t)=\operatorname{Tr}[\omega_E e^{iH^{\rm eff}_{\bar a}t}e^{-iH^{\rm eff}_at}]=\sum_mp_me^{-2i\delta_mt}$.

**The algebra is correct for $H_a^{\rm eff}=E_0+W_a$.** Because $W_{\bar a}=-W_a$, the product $e^{-iW_{\bar a}t}{}^\dagger e^{-iW_at}=e^{-2iW_at}$. What is wrong is the word *exact*. A Schrieffer–Wolff expansion with $Q=1-P_{\mathcal D}$ gives
$$H^{\rm eff}_a=E_0+W_a+K_a+O\!\left(\tfrac{v^3}{\Delta E^2}\right),\qquad K_a=P V_aQ\,(E_0-H_E)^{-1}QV_aP .$$
$K_a$ is **quadratic** in $V_a$, so $K_{\bar a}=K_a$. Therefore:

- $H^{\rm eff}_a=E_0+W_a+K_a$ and $H^{\rm eff}_{\bar a}=E_0-W_a+K_a$;
- $[H^{\rm eff}_a,H^{\rm eff}_{\bar a}]=2[W_a,K_a]\neq0$ in general;
- the product of exponentials no longer collapses to $e^{-2iW_at}$.

Together with the $O(v/\Delta E)$ leakage of $\lvert\psi_0\rangle$ out of the dressed sector, the first-order formula is controlled only for $v\ll\Delta E$ and $\lvert t\rvert\ll\Delta E/v^2$.

**Numerical confirmation.** The table below gives $\max_t\lvert\Gamma^{\rm exact}-\Gamma^{\rm 1st}\rvert$ for $t\le100$ (and $t\le300$):

| $N$ | pair $(00,11)$ | pair $(01,10)$ |
|---|---|---|
| 4 | 0.27 (0.27) | 0.025 (0.026) |
| 5 | 0.31 (0.32) | 0.047 (0.048) |
| 6 | 0.36 (0.36) | 0.041 (0.22) |
| 8 | 0.45 (0.45) | 0.16 (0.24) |
| 10 | 0.52 (0.52) | 0.079 (0.23) |

A striking case: for $(00,11)$ the computation gives $W_{00}\equiv0$ on $\mathcal D$ for every $N$, so first order predicts $\Gamma_{00,11}\equiv1$. The exact value drops to $0.46$–$0.73$. The plotted window $t\le2000$ is far beyond the first-order validity time. **Corrected wording:** "To first order in $V$ (valid for $\lVert V\rVert\ll\Delta E$, $t\ll\Delta E/\lVert V\rVert^2$)…".

---

## 3. E3: the logarithmic derivative for $d=2$

Let $\sigma=\delta_1+\delta_2$ and $\Delta\delta=\delta_1-\delta_2$. With $p_1=p_2=\tfrac12$,
$$\Gamma^{\mathcal D}_{a\bar a}(t)=e^{-i\sigma t}\cos(\Delta\delta\,t)
\;\Rightarrow\;
\eta^{\mathcal D}_{a\bar a}=\frac{\dot\Gamma}{\Gamma}=-i\sigma-\Delta\delta\tan(\Delta\delta\,t).$$
Near $t_n=\frac{(2n+1)\pi}{2\lvert\Delta\delta\rvert}$, put $x=\Delta\delta(t-t_n)$. Then $\tan(\Delta\delta t)=-\cot x$, and so
$$\boxed{\eta^{\mathcal D}_{a\bar a}(t)=\frac1{t-t_n}-i(\delta_1+\delta_2)-\frac{\Delta\delta^2}{3}(t-t_n)+O\big((t-t_n)^3\big).}$$
The slide's formula holds only when $\delta_1+\delta_2=0$, which is the $\Theta$-paired case. The damping rate $\kappa=-\operatorname{Re}\eta$ is unaffected; the phase rate is $\omega=-\sigma$.

---

## 4. E4: exact zeros versus near-zeros, and whether the rates actually diverge

### 4.1 What the exact symmetries really imply

Let $K$ denote complex conjugation in the computational basis, and let $\mathcal R$ be the site reflection $j\mapsto N+1-j$. Define
$$\Theta=\Big(\textstyle\prod_j\sigma^x_j\Big)K,\qquad T=e^{-i\varphi M_z}\,\mathcal R\,K,\qquad M_z=\textstyle\sum_jS^z_j .$$

**Facts about $\Theta$.** $\Theta$ maps $(S^x,S^y,S^z)\mapsto(S^x,S^y,-S^z)$. Consequently:

- $\Theta H_E\Theta^{-1}=H_E$: the XXZ bonds and both in-plane boundary fields are invariant;
- $\Theta V_a\Theta^{-1}=V_{\bar a}$;
- $\Theta\lvert\Psi_1^+\rangle=e^{i\chi}\lvert\Psi_1^+\rangle$: each equatorial site spinor maps to itself up to a phase.

**Facts about $T$.** $K$ maps the in-plane angle $\theta\mapsto-\theta$, $\mathcal R$ swaps sites $1\leftrightarrow N$, and $e^{-i\varphi M_z}$ adds $\varphi$. Consequently:

- the field directions are exchanged correctly, $0\mapsto\varphi$ and $\varphi\mapsto0$, so $TH_ET^{-1}=H_E$;
- $TS^z_1T^{-1}=S^z_N$;
- the helix angles map as $\theta_j=\pi+(j-2)\gamma\mapsto\varphi-\theta_{N+1-j}=\theta_j-2\pi$, so $T\lvert\Psi_1^+\rangle\propto\lvert\Psi_1^+\rangle$.

All of these identities were verified numerically to $\le2\times10^{-15}$ for $N=4,5,7,10$.

> **Proposition 2 ($\Theta$ gives evenness only).** For all $t$, $\Gamma_{a\bar a}(t)=\Gamma_{a\bar a}(-t)$.
>
> *Proof.* Since $\Theta$ is antiunitary, $\Theta U_a(t)\Theta^{-1}=e^{+iH_{\bar a}t}=U_{\bar a}(-t)$, and $\langle\Theta x\vert\Theta y\rangle=\langle y\vert x\rangle$. With $x=U_a(t)\psi$ and $y=U_{\bar a}(t)\psi$:
> $$\Gamma_{a\bar a}(t)=\langle U_{\bar a}(t)\psi\vert U_a(t)\psi\rangle=\langle\Theta U_a(t)\psi\vert\Theta U_{\bar a}(t)\psi\rangle=\langle U_{\bar a}(-t)\psi\vert U_a(-t)\psi\rangle=\Gamma_{a\bar a}(-t).$$
> The phases $e^{i\chi}$ cancel. $\square$
>
> $\Theta$ does **not** make $\Gamma_{a\bar a}$ real: its complex conjugate is $\Gamma_{\bar aa}(t)$, which is a different function.

> **Proposition 3 ($T$ makes $\Gamma_{00,11}$ real; requires $g_A=g_B$).** Under $T$, $V_a\mapsto V_{\sigma(a)}$, where $\sigma$ swaps $01\leftrightarrow10$ and fixes $00,11$. The same steps as above give $\Gamma_{ab}(t)=\Gamma_{\sigma(b)\sigma(a)}(-t)$. For $(00,11)$, combining this with Proposition 2 gives $\Gamma_{00,11}(t)=\Gamma_{00,11}(t)^*\in\mathbb R$. For $(01,10)$ it only reproduces Proposition 2. $\square$

The numerics agree: $\max_t\lvert\operatorname{Im}\Gamma_{00,11}\rvert\le3\times10^{-13}$, while $\Gamma_{01,10}$ is genuinely complex. The time-reversal identity $\lvert\Gamma(t)-\Gamma(-t)\rvert$ holds to $\le3\times10^{-14}$.

### 4.2 Why exact zeros of $\Gamma_{01,10}$ are non-generic

A real-time zero of a complex analytic $\Gamma$ requires two real equations, $\operatorname{Re}\Gamma(t_*)=\operatorname{Im}\Gamma(t_*)=0$, in one unknown $t_*$.

**Lemma 4.** Suppose $\Gamma(t_*)=0$ and $\Gamma'(t_*)=v\ne0$. Perturb $\Gamma$ to $\Gamma+\delta\Gamma$, with $\delta\Gamma$ smooth and of size $\epsilon$. Then
$$\min_{t\approx t_*}\lvert\Gamma+\delta\Gamma\rvert=\lvert\operatorname{Im}(\bar v\,\delta\Gamma(t_*))\rvert/\lvert v\rvert+O(\epsilon^2).$$
This is nonzero unless $\delta\Gamma(t_*)\parallel v$, which is a codimension-1 condition.

*Proof.* Near $t_*$, $\Gamma+\delta\Gamma\approx v(t-t_*)+\delta\Gamma(t_*)$, a straight line in $\mathbb C$. Its distance from $0$ is the component of $\delta\Gamma(t_*)$ orthogonal to $v$. $\square$

If instead a symmetry forces $\Gamma=e^{i\alpha(t)}R(t)$ with $R$ real, a zero is a sign change of $R$, which is stable (codimension 1). This is the case for $\Gamma_{00,11}$ (Proposition 3). It also holds within first order for $\Gamma_{01,10}$, where $\Gamma^{\mathcal D}=\cos(2\delta t)$. Beyond first order, only Proposition 2 survives for $(01,10)$, so the first-order zeros are lifted to near-zeros of depth set by the corrections of §2.

### 4.3 Numerical facts, independent recomputation for $t\in[0,2000]$

Method: a $\delta t=0.1$ grid, followed by golden-section refinement of each minimum.

| $N$ | $\min_t\lvert\Gamma_{01,10}\rvert$ (at $t$) | $\Gamma$ at the min | $\min_t\lvert\Gamma_{00,11}\rvert$ | resolved $\max\lvert\gamma_k\rvert$ |
|---|---|---|---|---|
| 4 | $5.1\times10^{-4}$ (1582.70) | $10^{-5}-5.1\cdot10^{-4}i$ | 0.734 | **66.7** |
| 5 | $1.5\times10^{-2}$ (1020.44) | $3\cdot10^{-5}-1.5\cdot10^{-2}i$ | 0.681 | **2.0** |
| 6 | $2.9\times10^{-3}$ | — | 0.634 | — |
| 7 | $8.9\times10^{-4}$ | — | 0.589 | — |
| 8 | $6.3\times10^{-4}$ | — | 0.539 | — |
| 9 | $4.6\times10^{-4}$ | — | 0.504 | — |
| 10 | $1.07\times10^{-4}$ (222.60) | $-2\cdot10^{-5}-1.05\cdot10^{-4}i$ | 0.455 | **9.7** |

"Resolved" means evaluated on a $5\times10^{-5}$ grid around the peak. These maxima equal the plotted values, so **the plotted peaks are smooth and finite, not sampled poles**.

Note the geometry at each minimum: $\operatorname{Re}\Gamma\approx0$ while $\operatorname{Im}\Gamma\ne0$. This is exactly what Lemma 4 predicts for a first-order real zero crossing that has been lifted by higher-order complex corrections. $\Gamma_{00,11}$ is real but never changes sign in the window, so it produces no zero either. All other pairs also stay strictly away from zero. The smallest value is $9.7\times10^{-4}$, for $(00,10)$ at $N=7$; the overall minimum over all pairs and all $N$ is the $1.07\times10^{-4}$ at $N=10$.

### 4.4 Corrected statements for the slides

- "$\Gamma_{ab}$ is a finite sum of exponentials, so any zeros are isolated, and rates diverge **only** at exact zeros." This part is correct.
- "For complex $\Gamma_{ab}$, exact real-time zeros are non-generic. They occur robustly only when a symmetry makes $\Gamma$ real up to a known phase."
- "In the phantom-helix runs ($N=4$–$10$, $t\le2000$) there are **no exact zeros**. The tall spikes are near-zeros, with height $\sim\lvert\dot\Gamma\rvert/\min\lvert\Gamma\rvert$: $66.7$ ($N=4$), $2.0$ ($N=5$), $9.7$ ($N=10$, $t=222.6$)."
- The $d=2$ slide should say that the **period** of the spikes is $\pi/\lvert\Delta\delta\rvert$, which follows from first order: $\delta=0.0224$ gives $70.3$ for $N=4$, and $\delta=0.0208$ gives $75.6$ for $N=5$. Their **depth**, and hence their finite height, is set by beyond-first-order corrections. Replace "periodic rate poles" with "periodic rate peaks".

---

## 5. E5: the $d\ge3$ zero condition

**Lemma 5 (closed polygon).** There exist phases $\theta_m$ with $\sum_mp_me^{i\theta_m}=0$ iff $\max_mp_m\le\tfrac12$.

This is only **necessary** for $\Gamma^{\mathcal D}(t)=\sum_mp_me^{-2i\delta_mt}$ to vanish, because the phases are restricted to the line $\theta_m=-2\delta_mt$.

**Proposition 6.** Let $d\ge3$, with the differences $\{\delta_m-\delta_1\}_{m\ge2}$ rationally independent.

1. The phase flow is dense in $\mathbb T^{d-1}$ (Kronecker), so $\inf_t\lvert\Gamma^{\mathcal D}(t)\rvert=0$ iff $\max p_m\le\tfrac12$.
2. The zero set $Z_p\subset\mathbb T^{d-1}$ has codimension $2$. For fixed $p$, the frequency vectors whose line meets $Z_p$ form a countable union of codimension-1 sets, which has Lebesgue measure zero.
3. Hence, for almost every frequency vector, $\Gamma^{\mathcal D}$ **never vanishes**. It only approaches $0$ arbitrarily closely over long times.

**Symmetric case.** With $\Theta$ pairing $\pm\delta_\alpha$ and equal weights $p_\alpha$, plus a weight $p_0$ on $\delta=0$, we get $\Gamma^{\mathcal D}=p_0+2\sum_\alpha p_\alpha\cos(2\delta_\alpha t)\in\mathbb R$. With $\{\delta_\alpha\}$ rationally independent, $\inf_t\Gamma^{\mathcal D}=2p_0-1$. Sign changes, and hence exact first-order zeros, exist iff $p_0<\tfrac12$. For $d=2$, $p_0=0$.

**Corrected takeaways.**

- *Takeaway 1* stays: within first order, $d=2$ already produces zeros, through two-branch interference with $p_1=p_2$ enforced by $\Theta$.
- *Takeaway 2* must change: $d>2$ does **not** remove fine-tuning for generic states. Without symmetry, exact zeros remain measure-zero for every $d$; what $d>2$ changes is the **pattern**, from a single-frequency comb to quasi-periodic, multi-frequency near-zeros.
- The sparseness of tall spikes at $N=10$ arises because only accidental deep near-zeros, such as $1.07\times10^{-4}$ at $t=222.6$, give large finite peaks.

---

## 6. E6: why $H_{\rm eff}$ is diagonal and $\mathsf d$ lives on $\{F_1,F_2,F_3\}$

The slide's reason ("$\mathcal L_t$ never mixes entries") is not by itself a proof. Rigorous argument:

> **Proposition 7.** Suppose $\mathcal L_t(E_{ab})=\eta_{ab}E_{ab}$ with $\eta_{aa}=0$, and $\mathcal L_t$ is written in Hall's canonical form, with $H$ traceless and $\mathsf d$ Hermitian in an orthonormal traceless basis. Then $H$ is diagonal and $\mathsf d$ is supported on $\operatorname{span}\{F_1,F_2,F_3\}$.

*Proof.*

1. **Covariance.** For every diagonal unitary $D_\theta=\operatorname{diag}(e^{i\theta_a})$, $\mathcal L_t(D_\theta XD_\theta^\dagger)=D_\theta\mathcal L_t(X)D_\theta^\dagger$, since $D_\theta E_{ab}D_\theta^\dagger=e^{i(\theta_a-\theta_b)}E_{ab}$ and $\mathcal L_t$ multiplies each $E_{ab}$ by a scalar.
2. **Transport by uniqueness.** Conjugating the canonical form shows that $\operatorname{Ad}_{D_\theta^\dagger}\circ\mathcal L_t\circ\operatorname{Ad}_{D_\theta}$ has canonical data $(D_\theta^\dagger HD_\theta,\ \mathcal C_\theta)$. Here $\mathcal C=\sum_{jk}\mathsf d_{jk}\lvert O_j\rangle\rangle\langle\langle O_k\rvert$ is the basis-free Kossakowski operator on traceless operators, and $\mathcal C_\theta=\operatorname{Ad}_{D_\theta^\dagger}\mathcal C\operatorname{Ad}_{D_\theta}$. By uniqueness of the canonical form, $D_\theta^\dagger HD_\theta=H$ and $\mathcal C_\theta=\mathcal C$ for all $\theta$.
3. **$H$ is diagonal.** The commutant of all $D_\theta$ is the diagonal algebra.
4. **Block structure of $\mathcal C$.** $\operatorname{Ad}_{D_\theta}$ acts on traceless operators with characters $e^{i(\theta_a-\theta_b)}$ on $E_{ab}$ ($a\ne b$) and trivially on the diagonal traceless subspace $\operatorname{span}\{F_\mu\}$. These characters are pairwise distinct, so $\mathcal C$ commuting with all of them must be block-diagonal:
$$\mathcal C=\mathcal C_{\rm diag}\oplus\bigoplus_{a\ne b}c_{ab}\lvert E_{ab}\rangle\rangle\langle\langle E_{ab}\rvert,\qquad c_{ab}\in\mathbb R .$$
5. **Populations are fixed.** Since $\eta_{bb}=0$, we have $\mathcal L_t(E_{bb})=0$. On diagonal inputs, the Hamiltonian part and the $\mathcal C_{\rm diag}$ part vanish, because $F\rho F-\frac12\{F^2,\rho\}=0$ when $[F,\rho]=0$. Each jump term gives $D[E_{ab}](E_{bb})=E_{aa}-E_{bb}$. Hence
$$0=\mathcal L_t(E_{bb})=\sum_{a\ne b}c_{ab}(E_{aa}-E_{bb}),$$
and linear independence of the $E_{aa}$ forces $c_{ab}=0$.

Therefore $\mathcal C=\mathcal C_{\rm diag}$. $\square$

The rest of the derivation (Steps 3–5) is then correct as written. In particular, the elementwise formula for $\eta_{ab}$, $\eta=G\mathsf dG^{\mathsf T}+(-ih-q)\mathbf1^{\mathsf T}+\mathbf1(ih-q)^{\mathsf T}$, the diagonal check $\eta_{aa}=\mathbf g^{a\mathsf T}S\mathbf g^a-2q_a=0$, and both inversion formulas were all re-derived and hold.

---

## 7. Remark (not an error): the Bell no-go is basis-specific

The slides correctly prove that the four **computational-basis** Bell states are unreachable from any $C<1$ input, and that $\langle\Phi^\pm\vert\rho_f\vert\Phi^\pm\rangle\le\frac12$ and $\langle\Psi^\pm\vert\rho_f\vert\Psi^\pm\rangle\le\frac12$. However, the Schur structure does **not** forbid maximal entanglement.

Example: $\Gamma=uu^\dagger$ with $u=(1,1,1,-1)$ is a valid correlation matrix (rank 1, unit diagonal). Its channel is $\Lambda(X)=\mathrm{CZ}\,X\,\mathrm{CZ}$, and
$$\rho_f=\mathrm{CZ}\lvert{++}\rangle\langle{++}\rvert\mathrm{CZ},\qquad C=1 .$$
Its Bell fidelities are exactly $\frac12$ (e.g. $\langle\Phi^-\vert\rho_f\vert\Phi^-\rangle=\frac12$), which saturates the slide's bound, and it is one local Hadamard away from $\lvert\Phi^+\rangle$.

So the defensible conclusion for Friday is: *the Ising channel cannot reshape the edge populations into those of $\lvert\Phi^\pm\rangle$ or $\lvert\Psi^\pm\rangle$, so it cannot output these states from any $C<1$ input, and its maximal separable-input entanglement is $C[\Gamma/4]$.* Whether a given bulk $\omega_E$ realizes a near-rank-1 $\Gamma$ with non-additive phase $\pi$ is a dynamical question, not one excluded by the channel structure. Section 9 answers it for the present setup.

---

## 9. Can the time-dependent Schur channel be exactly a CZ gate in our setup?

### 9.1 Unitary Schur channels are diagonal phase gates

$\Lambda_t$ is unitary iff $\operatorname{rank}\Gamma(t)=1$ (report, Sec. 3.7), i.e. iff $\Gamma_{ab}(t)=u_a(t)u_b(t)^*$ with $\lvert u_a\rvert=1$. Then $\Lambda_t(X)=D_uXD_u^\dagger$ with $D_u=\operatorname{diag}(u)$, and every such $D_u$ factorizes as
$$D_u=e^{i\alpha}\,\big(e^{-i\beta_A\sigma^z_A/2}\otimes e^{-i\beta_B\sigma^z_B/2}\big)\,\operatorname{diag}(1,1,1,e^{i\Phi}),\qquad
\Phi=\arg\!\big(u_{00}u_{11}u_{01}^*u_{10}^*\big).$$
Acting on $\lvert{++}\rangle$, the output is the pure state $\tfrac12(u_{00},u_{01},u_{10},u_{11})$, whose concurrence is
$$C=2\lvert\psi_{00}\psi_{11}-\psi_{01}\psi_{10}\rvert=\tfrac12\lvert u_{00}u_{11}-u_{01}u_{10}\rvert=\lvert\sin(\Phi/2)\rvert .$$
It is CZ-equivalent ($C=1$) iff $\Phi=\pi$, and then $\rho_f=\mathrm{CZ}\lvert{++}\rangle\langle{++}\rvert\mathrm{CZ}=\lvert\chi\rangle\langle\chi\rvert$ with $\lvert\chi\rangle=\tfrac1{\sqrt2}(\lvert0+\rangle+\lvert1-\rangle)$. This is maximally entangled but has populations $\tfrac14$, so it is **not** one of $\lvert\Phi^\pm\rangle,\lvert\Psi^\pm\rangle$.

### 9.2 Result A: a channel that is unitary on any time interval is never entangling

> **Proposition 8.** Assume $g_A,g_B\neq0$. If $\Lambda_t$ is unitary for all $t$ in some open interval $I$, then it is unitary for all $t\in\mathbb R$ and $\Phi(t)\equiv0$. So $\Lambda_t$ is a product of local $z$ rotations and never creates $A|B$ entanglement.

*Proof.*

1. **Extension to all times.** $\lvert\Gamma_{ab}(t)\rvert^2$ is real-analytic in $t$, being a finite sum of exponentials times its conjugate. If it equals $1$ on $I$, it equals $1$ on $\mathbb R$. So $\Lambda_t$ is unitary for all $t$.
2. **Collinearity of the Gram vectors.** Let $A_a(t)=U_a(t)\sqrt{\omega_E}$, so that $\Gamma_{ab}=\operatorname{Tr}(A_b^\dagger A_a)$ is a Gram matrix in Hilbert–Schmidt space with $\lVert A_a\rVert_2=1$. Rank one means $A_a(t)=u_a(t)M(t)$ for a common $M(t)$ with $\lVert M\rVert_2=1$. Choose $u_{00}\equiv1$, $M=A_{00}$; then $u_a=\operatorname{Tr}(A_{00}^\dagger A_a)$ is smooth.
3. **Eigenvector condition.** Differentiate $A_a=u_aM$ using $i\dot A_a=H_aA_a$ and $i\dot M=H_{00}M$. This gives $i\dot u_aM=u_a(H_a-H_{00})M$, i.e.
$$(V_a-V_{00})M(t)=\lambda_a(t)M(t),\qquad \lambda_a=i\dot u_a/u_a .$$
4. **Fixed eigenvalues.** For $a=10$, $V_{10}-V_{00}=-g_AS^z_1$; for $a=01$, $V_{01}-V_{00}=-g_BS^z_N$. Hence the range of $M(t)$ lies in an eigenspace of $S^z_1$ and of $S^z_N$, with eigenvalues $m_1(t),m_N(t)\in\{\pm\frac12\}$. Since each $\lambda_a$ is continuous and takes values in a discrete set, $m_1$ and $m_N$ are constant.
5. **The phase cancels.** Therefore $u_a(t)=\exp\{-it[g_A(s^a_A-\tfrac12)m_1+g_B(s^a_B-\tfrac12)m_N]\}$, which is affine in $s_A$ and $s_B$. In $\Phi$ the signs are $(+,-,-,+)$ for $(00,01,10,11)$, and $\sum_a\pm s^a_A=\sum_a\pm s^a_B=\sum_a\pm1=0$. Hence $\Phi(t)\equiv\Phi(0)=0$. $\square$

**Consequence.** Entanglement between the edges can only be generated while $\Lambda_t$ is **non-unitary**, i.e. while the bulk carries away which-path information. A CZ channel valid on any time window is impossible in this model.

### 9.3 Exact CZ at an isolated time is non-generic

$\operatorname{rank}\Gamma(t_*)=1$ at a single $t_*$ requires all four Gram vectors $A_a(t_*)$ to realign simultaneously. For a pure bulk state this means $\lvert\langle\phi_b(t_*)\vert\phi_a(t_*)\rangle\rvert=1$ for every pair, where each $\lvert\Gamma_{ab}\rvert$ reaches its absolute maximum, together with $\Phi(t_*)=\pi$: many real conditions on one real unknown. Quasi-periodicity of $\Gamma(t)$ guarantees approximate returns to $\Gamma(0)=\mathbf1\mathbf1^{\mathsf T}$, which is the identity channel with $\Phi=0$. It does **not** guarantee returns to a rank-one $\Gamma$ with $\Phi=\pi$. An exact CZ at an isolated time therefore requires fine-tuned spectral relations and is not expected.

### 9.4 An approximate CZ is possible: weak coupling to an isolated bulk eigenstate

Take $\omega_E=\lvert E_0\rangle\langle E_0\rvert$ non-degenerate, in the setting of Lemma 1, with $\varepsilon=v/(\Delta E-v)\ll1$ and $v=\lVert V_a\rVert$.

*Step 1 (the channel is nearly unitary, uniformly in time).* From Lemma 1, $\lvert\phi_a(t)\rangle=c_ae^{-iE_at}\lvert E_0^a\rangle+e^{-iH_at}\lvert r_a\rangle$, where $\lVert r_a\rVert\le\varepsilon$ and $e^{-iH_at}r_a\perp E_0^a$. Because $\lvert\langle E_0^b\vert x\rangle\rvert\le\sin\angle(E_0^a,E_0^b)=O(\varepsilon)$ for every unit $x\perp E_0^a$, all three cross terms in $\langle\phi_b\vert\phi_a\rangle$ are $O(\varepsilon^2)$. Hence
$$\Gamma_{ab}(t)=u_a(t)u_b(t)^*+O(\varepsilon^2)\quad\text{uniformly in }t,\qquad u_a(t)=\tfrac{c_a}{\lvert c_a\rvert}e^{-iE_at}.$$

*Step 2 (the entangling phase).* $\Phi(t)=-t\,J_{\rm eff}+\Phi_0$, where $J_{\rm eff}:=E_{00}+E_{11}-E_{01}-E_{10}$ is built from the **exact** dressed energies. Every function of $s_A,s_B\in\{\pm\frac12\}$ reduces to $\alpha+\beta s_A+\gamma s_B+\delta s_As_B$, and only $\delta$ survives the alternating sum. Rayleigh–Schrödinger perturbation theory with $\langle n\vert V_a\vert0\rangle=g_As_AA_n+g_Bs_BB_n$, where $A_n=\langle n\vert S^z_1\vert0\rangle$ and $B_n=\langle n\vert S^z_N\vert0\rangle$, then gives
$$J_{\rm eff}=g_Ag_B\,\chi_{1N}+O\!\left(\frac{g^4}{\Delta E^3}\right),\qquad
\chi_{1N}=2\operatorname{Re}\sum_{n\ne0}\frac{\langle0\vert S^z_1\vert n\rangle\langle n\vert S^z_N\vert0\rangle}{E_0-E_n}.$$
The first-order terms, and the third-order terms (which carry $s_A^2=s_B^2=\frac14$), are affine and drop out.

*Step 3 (fidelity and concurrence).* At $t_*$ with $\Phi(t_*)=\pi$, i.e. $t_*\simeq\pi/\lvert g_Ag_B\chi_{1N}\rvert$, the fidelity with $\lvert\chi_u\rangle=\frac12 D_u\lvert{++}\rangle$ is $F=\frac1{16}\sum_{ab}u_a^*\Gamma_{ab}u_b=1-O(\varepsilon^2)$. The two-qubit inequality "fully entangled fraction $\le(1+C)/2$" (Verstraete–Verschelde) then gives
$$C[\rho_f(t_*)]\ge2F-1=1-O(\varepsilon^2).$$

**Answer.**

- *Exactly CZ:* no. It is impossible on any interval (Proposition 8) and non-generic at isolated times (§9.3).
- *Approximately CZ:* yes, in principle, if the bulk starts in an **isolated non-degenerate eigenstate** with $\chi_{1N}\neq0$, on the long time scale $t_*\simeq\pi/(g_Ag_B\lvert\chi_{1N}\rvert)$, which scales as $\Delta E/g^2$ because $\chi_{1N}=O(1/\Delta E)$. This regime is the opposite of the degenerate helix runs. A thermal (mixed) bulk washes it out, because $J_{\rm eff}$ depends on the eigenstate.
- *Even then:* the output is $\mathrm{CZ}\lvert{++}\rangle$, not a computational-basis Bell state. Obtaining $\lvert\Phi^+\rangle$ needs an additional local Hadamard on $B$. The basis-specific no-go of the slides stands.

---

## 10. Finite peaks versus genuine divergences: the analytic argument

### 10.1 What is exact

For all $t$, the canonical rates are finite iff every $\Gamma_{ab}(t)\neq0$, because $\eta_{ab}=\dot\Gamma_{ab}/\Gamma_{ab}$ and the rates are eigenvalues of $G^{\mathsf T}\eta G$. Conversely, if $\Gamma_{ab}(t_*)=0$ with order $r$, then $\kappa_{ab}\simeq-r/(t-t_*)$ and the bound $\gamma_{\min}\le\kappa_{ab}\le\gamma_{\max}$ forces divergence.

Each $\Gamma_{ab}$ extends to an **entire** function of complex $t$ (a finite exponential sum, not identically zero since $\Gamma_{ab}(0)=1$), so its complex zeros are isolated. The whole question is therefore whether any zero lies **exactly on the real axis**.

### 10.2 Real-on-the-axis versus complex: Schwarz reflection

- **$\Gamma_{00,11}$ is real for real $t$** (Proposition 3). By the Schwarz reflection principle, $\Gamma_{00,11}(\bar t)=\overline{\Gamma_{00,11}(t)}$, so its complex zeros come in conjugate pairs. As a parameter $\lambda$ (e.g. $g$) varies, a pair can collide on the real axis and split into two real zeros. Real zeros, i.e. sign changes, therefore occupy **open** regions of parameter space and are robust: if they occur, they are genuine divergences.
- **$\Gamma_{01,10}$ is only even**, $\Gamma(-t)=\Gamma(t)$ (Proposition 2). Its zero set is symmetric under $t\mapsto-t$ but **not** under $t\mapsto\bar t$. By the implicit function theorem, a simple zero $t_k(\lambda)$ moves analytically with $\lambda$. A real zero requires the single real condition $\operatorname{Im}t_k(\lambda)=0$, which is met only on a codimension-1 set of parameter values (isolated points along a one-parameter family). At **generic fixed parameters**, $\Gamma_{01,10}$ has no real zero.

### 10.3 Shape of a near-zero: why the peaks are tall but finite

Let $t_0$ be the point of closest approach to $0$, where $\Gamma(t_0)=e^{i\beta}i\epsilon'$ with $\epsilon'\in\mathbb R$, $v=\dot\Gamma(t_0)=\lvert v\rvert e^{i\beta}$, and $\tau=t-t_0$. Then $\Gamma(t)\simeq e^{i\beta}(\lvert v\rvert\tau+i\epsilon')$, and
$$\kappa(t)=-\frac{\lvert v\rvert^2\tau}{\lvert v\rvert^2\tau^2+\epsilon'^2},\qquad
\omega(t)=-\frac{\lvert v\rvert\epsilon'}{\lvert v\rvert^2\tau^2+\epsilon'^2}.$$

- $\kappa$ is a bounded, antisymmetric "dispersive" peak with extrema $\mp\lvert v\rvert/(2\lvert\epsilon'\rvert)$ at $\tau=\pm\lvert\epsilon'\rvert/\lvert v\rvert$.
- $\omega$ is a Lorentzian of height $\lvert v\rvert/\lvert\epsilon'\rvert$ and area $-\pi\operatorname{sgn}\epsilon'$: the phase of $\Gamma$ winds smoothly by $\pi$.
- By the Rayleigh bound, the canonical rates inherit peaks of height $O(\lvert v\rvert/\lvert\epsilon'\rvert)$.
- A **true divergence** is the singular limit $\epsilon'\to0$: $\kappa\to-1/\tau$, and $\omega\to-\pi\operatorname{sgn}(\epsilon')\delta(\tau)$, a discontinuous phase jump.

### 10.4 Application to the plotted peaks

- Within first-order degenerate theory, $\Theta$ makes $\Gamma^{\mathcal D}_{01,10}$ real ($=\cos 2\delta t$ for $d=2$), so the slides' poles are genuine **in that model**.
- In the exact dynamics only evenness survives for $(01,10)$, so each first-order zero is displaced to a near-zero with $\epsilon'$ equal to the component of the beyond-first-order correction orthogonal to $v$ (Lemma 4). The generic expectation is therefore a **tall but finite** peak, of height $\sim\lvert\dot\Gamma\rvert/\epsilon'$.
- A genuine divergence at the actual parameters would be an accidental coincidence. It cannot be excluded by a theorem about a specific transcendental function, but it is non-generic.
- The notebooks' `Missing["NearSingularMap"]` flag ($\lvert\Gamma\rvert<10^{-8}$) never fires, which is consistent with this.

### 10.5 The dependence on the bulk preparation

If $H_E$ is real in the computational basis (no $S^y$ fields) and $\psi_0$ is real, complex conjugation $K$ is a symmetry with $KV_aK=V_a$, so $\Gamma_{ab}(t)^*=\Gamma_{ab}(-t)$ for every pair. If in addition $\Theta\psi_0\propto\psi_0$, Proposition 2 makes **every** $\Gamma_{a\bar a}$ real. Its zero crossings, and hence the divergences of the canonical rates, are then robust.

The complex equatorial helix (phases $e^{i(j-2)\gamma}$, plus the $S^y$ component of the right field) breaks $K$. It leaves only $\Theta$ and the reflection $T$, so robust reality survives for $(00,11)$ but not for $(01,10)$. Whether the canonical generator has genuine poles, or only tall finite peaks, is therefore a property of the **initial bulk state** and its symmetries. This is consistent with the point that the time-local description is an extra, preparation-dependent assumption.

---

## 11. Changes applied to `friday_beamer.tex`

All the corrections E1–E6 and m1–m5 were applied. The full perturbative material (Lemma 1, first-order degenerate theory with its validity conditions, the $d=2$ and $d\ge3$ analyses, the symmetry propositions, the near-zero shape, and §9) was moved to an appendix. The main part keeps only the qualitative argument contrasting an isolated (non-degenerate) bulk eigenstate with a degenerate one.

---

## 8. Reproducibility

The scripts are in the conversation scratch folder:

- [check_gamma_zeros.py](file:///C:/Users/Miguel/.gemini/antigravity-cli/brain/934a3a58-6980-409b-9820-55ee19336bee/scratch/check_gamma_zeros.py): spectrum and gap, $\Theta$ check, exact $\Gamma_{ab}(t)$ and $\dot\Gamma_{ab}(t)$ on the notebook grid, refined minima, canonical rates, and the first-order comparison;
- [check_T_symmetry.py](file:///C:/Users/Miguel/.gemini/antigravity-cli/brain/934a3a58-6980-409b-9820-55ee19336bee/scratch/check_T_symmetry.py): the reflection antiunitary $T$;
- [check_peak_height.py](file:///C:/Users/Miguel/.gemini/antigravity-cli/brain/934a3a58-6980-409b-9820-55ee19336bee/scratch/check_peak_height.py): resolved peak heights.

All use $J=1$, $g_A=g_B=0.1$, $\gamma=\pi/3$, family 1, $+$, PMC, and bulk $\lvert\Psi_1^+\rangle$, matching `settings.txt` of `Phantom_ZZ_N10_run2`.
