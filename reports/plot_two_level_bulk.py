"""Two-level bulk model: exact |Gamma(t)| and damping rate kappa(t) = -d/dt ln|Gamma|."""
import pathlib
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "legend.fontsize": 8.6,
})

tau = np.linspace(0.0, np.pi, 4000)  # dimensionless time delta * t
cases = [
    (20.0, "#5f6e7d", "--", r"$\Delta/\delta = 20$ (gapped)"),
    (4.0,  "#006e7d", "-.", r"$\Delta/\delta = 4$"),
    (0.0,  "#b83227", "-",  r"$\Delta/\delta = 0$ (degenerate)"),
]

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(4.9, 3.85), sharex=True,
    gridspec_kw={"height_ratios": [1.0, 1.05], "hspace": 0.14},
)

for r, color, ls, label in cases:
    omega = np.sqrt(0.25 * r**2 + 1.0)
    gamma = 1.0 - (2.0 / omega**2) * np.sin(omega * tau)**2
    dgamma = -(2.0 / omega) * np.sin(2.0 * omega * tau)
    kappa = np.where(np.abs(gamma) > 2e-2, -dgamma / gamma, np.nan)
    ax1.plot(tau, np.abs(gamma), ls=ls, color=color, lw=1.7, label=label)
    ax2.plot(tau, kappa, ls=ls, color=color, lw=1.7)

for x0 in (0.25 * np.pi, 0.75 * np.pi):
    ax1.axvline(x0, color="#b83227", lw=0.7, ls=":", alpha=0.55)
    ax2.axvline(x0, color="#b83227", lw=0.7, ls=":", alpha=0.55)

ax1.set_ylabel(r"$|\Gamma_{a\bar a}(t)|$")
ax1.set_ylim(-0.04, 1.08)
ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=2, frameon=True, framealpha=0.95)
ax2.axhline(0.0, color="#888888", lw=0.6)
ax2.set_ylabel(r"$\kappa(t)/\delta$")
ax2.set_xlabel(r"$\delta\,t$")
ax2.set_ylim(-12.0, 12.0)
ax2.set_xlim(0.0, np.pi)
ax2.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi],
               [r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])

out = pathlib.Path(__file__).resolve().parent
fig.savefig(out / "two_level_bulk.pdf", bbox_inches="tight")
fig.savefig(out / "two_level_bulk.png", dpi=220, bbox_inches="tight")
plt.close(fig)
