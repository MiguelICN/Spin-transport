"""Generate auxiliary figures for phantom_zz_n10_beamer.tex from results.wxf (or cached .npz).

Outputs in reports/Phantom_ZZ_N10/new/:
  - zoom_222.png          : 3-panel zoom on t in [200, 245] around the near-zero at t = 222.60
  - minsv_log.png         : log10 sigma_min(S_t) over [0, 500] with inset around t = 222.60
  - independent_gammas.png: |Gamma_ab(t)| and arg(Gamma_ab(t))/pi for the 4 independent pairs
"""
import argparse
import pathlib
import struct
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

DTYPE_MAP = {
    0x00: np.dtype('<i1'), 0x01: np.dtype('<i2'), 0x02: np.dtype('<i4'), 0x03: np.dtype('<i8'),
    0x10: np.dtype('<u1'), 0x11: np.dtype('<u2'), 0x12: np.dtype('<u4'), 0x13: np.dtype('<u8'),
    0x22: np.dtype('<f4'), 0x23: np.dtype('<f8'), 0x33: np.dtype('<c8'), 0x34: np.dtype('<c16'),
}
SKIP_KEYS = {'Choi', 'ChoiState', 'Transfer', 'TransferDerivative', 'Generator', 'DephasingDerivative'}


class MiniWXFReader:
    def __init__(self, buf: memoryview):
        self.buf = buf
        self.pos = 0

    def read_varint(self) -> int:
        shift = 0
        res = 0
        while True:
            b = self.buf[self.pos]
            self.pos += 1
            res |= (b & 0x7F) << shift
            if not (b & 0x80):
                return res
            shift += 7

    def skip_expr(self):
        tok = self.buf[self.pos]
        self.pos += 1
        if tok == 0x66:
            n = self.read_varint()
            self.skip_expr()
            for _ in range(n):
                self.skip_expr()
        elif tok in (0x73, 0x53, 0x42, 0x49, 0x52):
            self.pos += self.read_varint()
        elif tok == 0x43:
            self.pos += 1
        elif tok == 0x6A:
            self.pos += 2
        elif tok == 0x69:
            self.pos += 4
        elif tok in (0x4C, 0x72):
            self.pos += 8
        elif tok == 0x41:
            for _ in range(self.read_varint()):
                self.pos += 1
                self.skip_expr()
                self.skip_expr()
        elif tok in (0xC1, 0xC2):
            vt = self.buf[self.pos]
            self.pos += 1
            cnt = 1
            for _ in range(self.read_varint()):
                cnt *= self.read_varint()
            self.pos += cnt * DTYPE_MAP[vt].itemsize

    def read_expr(self):
        tok = self.buf[self.pos]
        self.pos += 1
        if tok == 0x66:
            n = self.read_varint()
            h = self.read_expr()
            args = [self.read_expr() for _ in range(n)]
            if h == ('sym', 'List'):
                return args
            if h == ('sym', 'Complex'):
                return complex(args[0], args[1])
            if h == ('sym', 'Rational'):
                return args[0] / args[1]
            if h == ('sym', 'Missing'):
                return np.nan
            return ('func', h, args)
        elif tok == 0x73:
            n = self.read_varint()
            s = bytes(self.buf[self.pos:self.pos + n]).decode('utf-8')
            self.pos += n
            return {'True': True, 'False': False, 'Null': None, 'Pi': np.pi}.get(s, ('sym', s))
        elif tok == 0x53:
            n = self.read_varint()
            s = bytes(self.buf[self.pos:self.pos + n]).decode('utf-8')
            self.pos += n
            return s
        elif tok == 0x43:
            v = struct.unpack_from('<b', self.buf, self.pos)[0]; self.pos += 1; return v
        elif tok == 0x6A:
            v = struct.unpack_from('<h', self.buf, self.pos)[0]; self.pos += 2; return v
        elif tok == 0x69:
            v = struct.unpack_from('<i', self.buf, self.pos)[0]; self.pos += 4; return v
        elif tok == 0x4C:
            v = struct.unpack_from('<q', self.buf, self.pos)[0]; self.pos += 8; return v
        elif tok == 0x72:
            v = struct.unpack_from('<d', self.buf, self.pos)[0]; self.pos += 8; return v
        elif tok in (0x49, 0x52):
            n = self.read_varint()
            s = bytes(self.buf[self.pos:self.pos + n]).decode('ascii'); self.pos += n
            return int(s) if tok == 0x49 else float(s)
        elif tok == 0x41:
            d = {}
            for _ in range(self.read_varint()):
                self.pos += 1
                k = self.read_expr()
                if k in SKIP_KEYS:
                    self.skip_expr()
                else:
                    d[k] = self.read_expr()
            return d
        elif tok in (0xC1, 0xC2):
            vt = self.buf[self.pos]; self.pos += 1
            dims = tuple(self.read_varint() for _ in range(self.read_varint()))
            dt = DTYPE_MAP[vt]
            cnt = int(np.prod(dims))
            arr = np.frombuffer(self.buf, dtype=dt, count=cnt, offset=self.pos).reshape(dims).copy()
            self.pos += cnt * dt.itemsize
            return arr
        raise ValueError(f"Unknown token {tok:#x}")


def load_series(here: pathlib.Path, npz_path: str | None = None):
    if npz_path and pathlib.Path(npz_path).exists():
        z = np.load(npz_path)
        return (
            z['times'], z['gamma'], z['rates'], z['C_MinSingularValue'],
            z['S_EntropyAB'], z['S_Concurrence'], z['C_BLPDistance'],
        )
    wxf_path = here.parent / "results.wxf"
    raw = wxf_path.read_bytes()
    d = MiniWXFReader(memoryview(raw[2:])).read_expr()
    times = d['Times']
    recs, mrecs = d['Records'], d['MapRecords']
    gamma = np.array([r['DephasingFactors'] for r in mrecs])
    rates = np.array([r['CanonicalRates'] for r in mrecs])
    minsv = np.array([r['MinSingularValue'] for r in mrecs])
    s_ab = np.array([r['EntropyAB'] for r in recs])
    conc = np.array([r['Concurrence'] for r in recs])
    blp_d = np.array([r['BLPDistance'] for r in mrecs])
    return times, gamma, rates, minsv, s_ab, conc, blp_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=None, help="Optional pre-extracted .npz path")
    ap.add_argument("--outdir", default=None, help="Optional output directory")
    args = ap.parse_args()

    here = pathlib.Path(__file__).resolve().parent
    outdir = pathlib.Path(args.outdir) if args.outdir else here
    outdir.mkdir(parents=True, exist_ok=True)

    times, gamma, rates, minsv, s_ab, conc, blp_d = load_series(here, args.npz)

    plt.rcParams.update({
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "legend.fontsize": 8.3,
        "axes.edgecolor": "#192b3c",
        "xtick.color": "#192b3c",
        "ytick.color": "#192b3c",
    })

    # Locate exact near-zero around 222.6
    sl_zoom = (times >= 200.0) & (times <= 245.0)
    tz = times[sl_zoom]
    g0110_z = np.abs(gamma[sl_zoom, 1, 2])
    imin = int(np.argmin(g0110_z))
    t0 = float(tz[imin])
    gmin = float(g0110_z[imin])

    # 1. zoom_222.png
    exp10 = int(np.floor(np.log10(gmin)))
    mant = gmin / (10.0 ** exp10)
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.3, 4.45), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.15, 1.0], "hspace": 0.12},
    )
    ax1.plot(tz, np.log10(g0110_z), color="#b83227", lw=1.7, label=r"$\log_{10}|\Gamma_{01,10}(t)|$")
    ax1.scatter([t0], [np.log10(gmin)], color="#b83227", s=18, zorder=5)
    ax1.annotate(
        rf"$|\Gamma_{{01,10}}({t0:.2f})| = {mant:.2f}\times 10^{{{exp10}}}$" + "\n" + r"(above $10^{-8}$ cutoff)",
        xy=(t0, np.log10(gmin)), xytext=(t0 + 3.8, -3.65),
        fontsize=8.2, color="#192b3c",
        arrowprops={"arrowstyle": "->", "color": "#b83227", "lw": 0.9},
    )
    ax1.set_ylabel(r"$\log_{10}|\Gamma_{01,10}|$")
    ax1.set_ylim(-4.35, -0.5)
    ax1.legend(loc="lower left", frameon=True, framealpha=0.92)
    ax1.grid(True, alpha=0.22, lw=0.5)

    ax2.axhline(0.0, color="#888888", lw=0.6)
    ax2.plot(tz, rates[sl_zoom, 14], color="#006e7d", lw=1.5, label=r"$\gamma_{15}(t)$ (max rate)")
    ax2.plot(tz, rates[sl_zoom, 0], color="#b83227", lw=1.5, label=r"$\gamma_{1}(t)$ (min rate)")
    ax2.set_ylabel(r"$\gamma_k(t)$")
    ax2.set_ylim(-10.8, 10.8)
    ax2.legend(loc="upper right", ncol=1, frameon=True, framealpha=0.92)
    ax2.grid(True, alpha=0.22, lw=0.5)

    ax3.plot(tz, s_ab[sl_zoom], color="#192b3c", lw=1.7, label=r"$S(\rho_{AB})$ (nats)")
    ax3.plot(tz, blp_d[sl_zoom], color="#006e7d", lw=1.7, ls="--", label=r"BLP $D(t)$")
    ax3.plot(tz, conc[sl_zoom], color="#d97706", lw=1.6, ls="-.", label=r"$C(\rho_{AB}) \equiv 0$")
    ax3.set_ylabel("Smooth obs.")
    ax3.set_xlabel(r"$t$")
    ax3.set_ylim(-0.06, 1.36)
    ax3.set_xlim(200.0, 245.0)
    ax3.legend(loc="center right", ncol=3, frameon=True, framealpha=0.92)
    ax3.grid(True, alpha=0.22, lw=0.5)

    for ax in (ax1, ax2, ax3):
        ax.axvline(t0, color="#5f6e7d", lw=0.85, ls=":")

    fig.savefig(outdir / "zoom_222.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # 2. minsv_log.png
    fig, ax = plt.subplots(figsize=(5.8, 3.1))
    ax.plot(times, np.log10(minsv), color="#006e7d", lw=1.35,
            label=r"$\log_{10}\sigma_{\min}(\mathsf{S}_t) = \log_{10}\min_{a<b}|\Gamma_{ab}(t)|$")
    ax.axhline(-8.0, color="#5f6e7d", lw=0.95, ls="--", label=r"Code cutoff $10^{-8}$")
    ax.scatter([t0], [np.log10(gmin)], color="#b83227", s=22, zorder=5)
    ax.annotate(
        rf"$t={t0:.2f}$: ${mant:.2f}\times 10^{{{exp10}}}$",
        xy=(t0, np.log10(gmin)), xytext=(72, -4.85),
        fontsize=8.5, color="#b83227",
        arrowprops={"arrowstyle": "->", "color": "#b83227", "lw": 0.9},
    )
    for tx in (75.37, 361.09, 373.80):
        idx_near = int(np.argmin(np.abs(times - tx)))
        w = slice(max(0, idx_near - 300), min(len(times), idx_near + 300))
        im_loc = w.start + int(np.argmin(minsv[w]))
        ax.scatter([times[im_loc]], [np.log10(minsv[im_loc])], color="#192b3c", s=13, zorder=4)

    ax.set_xlim(0.0, 500.0)
    ax.set_ylim(-8.6, 0.35)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\log_{10}\sigma_{\min}(\mathsf{S}_t)$")
    ax.legend(loc="lower left", frameon=True, framealpha=0.93, fontsize=7.9)
    ax.grid(True, alpha=0.22, lw=0.5)

    # Inset resolving the smooth bottom around t = 222.60
    ins = inset_axes(ax, width="33%", height="38%", loc="lower right", borderpad=1.65)
    sl_in = (times >= 222.35) & (times <= 222.85)
    ins.plot(times[sl_in], minsv[sl_in] * 1e4, color="#b83227", lw=1.4)
    ins.scatter([t0], [gmin * 1e4], color="#b83227", s=14)
    ins.set_title(r"Zoom ($\times 10^{-4}$)", fontsize=7.6, pad=2.5)
    ins.set_xticks([222.4, 222.6, 222.8])
    ins.tick_params(labelsize=7)
    ins.set_xlim(222.35, 222.85)
    ins.grid(True, alpha=0.25, lw=0.4)

    fig.savefig(outdir / "minsv_log.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # 3. independent_gammas.png (4 independent functions under RF symmetry)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.0, 3.85), sharex=True,
        gridspec_kw={"height_ratios": [1.05, 1.0], "hspace": 0.14},
    )
    specs = [
        ((0, 3), "#4c9a2a", "-",  r"$\Gamma_{00,11}\in\mathbb{R}$ (double-quantum)"),
        ((1, 2), "#b83227", "-",  r"$\Gamma_{01,10}$ (flip-flop)"),
        ((0, 1), "#5b4bc4", "-",  r"$\Gamma_{00,01}=\Gamma_{01,11}^*$"),
        ((0, 2), "#b45309", "--", r"$\Gamma_{00,10}=\Gamma_{10,11}^*$"),
    ]
    for (i, j), col, ls, lbl in specs:
        gij = gamma[:, i, j]
        ax1.plot(times, np.abs(gij), color=col, ls=ls, lw=1.45, label=lbl)
        phase_unwrapped = np.unwrap(np.angle(gij)) / np.pi
        ax2.plot(times, phase_unwrapped, color=col, ls=ls, lw=1.35)

    ax1.set_ylabel(r"$|\Gamma_{ab}(t)|$")
    ax1.set_ylim(-0.03, 1.05)
    ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=2, frameon=True, framealpha=0.94, fontsize=8.0)
    ax1.grid(True, alpha=0.22, lw=0.5)

    ax2.set_ylabel(r"Unwrapped $\arg\Gamma_{ab}/\pi$")
    ax2.set_xlabel(r"$t$")
    ax2.set_xlim(0.0, 500.0)
    ax2.grid(True, alpha=0.22, lw=0.5)

    fig.savefig(outdir / "independent_gammas.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("Generated zoom_222.png, minsv_log.png, independent_gammas.png in", outdir)


if __name__ == "__main__":
    main()
