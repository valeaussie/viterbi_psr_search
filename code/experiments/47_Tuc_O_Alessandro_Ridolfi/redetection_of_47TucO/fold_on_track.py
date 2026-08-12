#!/usr/bin/env python3
"""
fold_on_track.py
Fold a barycentred dedispersed PRESTO time series (.dat) with the identical
in-house folder, choosing the phase model with --model:

  track  : per-segment Viterbi frequencies (track.dat), held piecewise-constant
  poly   : polynomial fit F0..Fn from fit_summary.dat (smooth f(t))
  kepler : Keplerian fit f_spin,K,Pb,T0 from fit_summary.dat (smooth f(t))

Same data, same binning, same S/N definition, NO search: the only thing that
changes between runs is the phase model. The .dat is barycentred, so no frame
conversion is applied.

Usage:
  python fold_on_track.py --model track \
      --dat FILE.dat --inf FILE.inf --track track.dat
  python fold_on_track.py --model poly \
      --dat FILE.dat --inf FILE.inf --track track.dat --summary fit_summary.dat
  python fold_on_track.py --model kepler \
      --dat FILE.dat --inf FILE.inf --track track.dat --summary fit_summary.dat
"""
import argparse, math, sys
import numpy as np


def read_inf(path):
    tsamp = mjd = None
    with open(path) as fh:
        for line in fh:
            low = line.lower()
            if "width of each time series bin" in low:
                tsamp = float(line.split("=")[-1])
            elif "epoch of observation" in low:
                mjd = float(line.split("=")[-1])
    if tsamp is None:
        sys.exit("ERROR: could not read tsamp from .inf")
    return tsamp, mjd


def read_track(path, freq_col, time_col):
    data = np.loadtxt(path, comments="#", ndmin=2)
    ncol = data.shape[1]
    print(f"track.dat: {data.shape[0]} rows x {ncol} cols")
    for c in range(ncol):
        col = data[:, c]
        print(f"  col {c}: min={col.min():.6g} max={col.max():.6g} median={np.median(col):.6g}")
    if time_col is None:
        cand = [c for c in range(ncol) if np.all(np.diff(data[:, c]) > 0)]
        time_col = cand[0] if cand else 0
    if freq_col is None:
        cand = [c for c in range(ncol) if c != time_col
                and 1.0 < np.median(data[:, c]) < 2000.0]
        freq_col = cand[0] if cand else 1
    print(f"Using time_col={time_col}, freq_col={freq_col} "
          f"(override with --time-col / --freq-col if wrong)")
    return data[:, time_col].astype(float), data[:, freq_col].astype(float)


def read_summary(path):
    """Parse 'key value [err]' lines from fit_summary.dat into {key: float}."""
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    d[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return d


def load_dat_trimmed(path):
    x = np.fromfile(path, dtype=np.float32)
    nz = np.nonzero(x)[0]
    if nz.size:
        x = x[: nz[-1] + 1]
    return x


# ---------------- phase models (all return phase in turns) ----------------
def make_track_phase(tmid, ftrack):
    Tcoh = tmid[1] - tmid[0]
    edges = np.concatenate(([tmid[0] - Tcoh / 2.0], tmid + Tcoh / 2.0))
    cum = np.concatenate(([0.0], np.cumsum(ftrack * Tcoh)))
    nseg = ftrack.size

    def phase(t):
        seg = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, nseg - 1)
        return cum[seg] + ftrack[seg] * (t - edges[seg])
    return phase


def make_poly_phase(coeffs, tref):
    # coeffs = [F0, F1, ...] as derivatives; f(t)=sum Fk/k! (t-tref)^k
    # phase(t)=sum Fk/(k+1)! (t-tref)^(k+1)
    def phase(t):
        dt = t - tref
        out = np.zeros_like(t, dtype=np.float64)
        for k, fk in enumerate(coeffs):
            out += fk / math.factorial(k + 1) * dt ** (k + 1)
        return out
    return phase


def make_kepler_phase(fspin, K, Pb, T0):
    w = 2.0 * np.pi / Pb
    def phase(t):
        return fspin * t + (K / w) * np.cos(w * (t - T0))
    return phase


def model_freq(name, t, tmid, ftrack, coeffs, tref, kep):
    """f(t) for the model, used only for the sanity check against the track."""
    if name == "track":
        Tcoh = tmid[1] - tmid[0]
        idx = np.clip(((t - (tmid[0] - Tcoh / 2.0)) / Tcoh).astype(int), 0, ftrack.size - 1)
        return ftrack[idx]
    if name == "poly":
        dt = t - tref
        return sum(coeffs[k] / math.factorial(k) * dt ** k for k in range(len(coeffs)))
    if name == "kepler":
        return kep["f_spin"] - kep["K"] * np.sin(2 * np.pi * (t - kep["T0"]) / kep["Pb"])


def fold(x, tsamp, phase_fn, nbin, nsub, chunk=5_000_000):
    N = x.size
    Tobs = N * tsamp
    prof_sum = np.zeros(nbin); prof_cnt = np.zeros(nbin)
    sub = np.zeros((nsub, nbin)); subcnt = np.zeros((nsub, nbin))
    for s0 in range(0, N, chunk):
        s1 = min(s0 + chunk, N)
        j = np.arange(s0, s1)
        t = (j + 0.5) * tsamp
        phi = phase_fn(t)
        b = np.minimum(((phi % 1.0) * nbin).astype(np.int64), nbin - 1)
        w = x[s0:s1].astype(np.float64)
        prof_sum += np.bincount(b, weights=w, minlength=nbin)
        prof_cnt += np.bincount(b, minlength=nbin)
        si = np.minimum((t / Tobs * nsub).astype(np.int64), nsub - 1)
        flat = si * nbin + b
        sub += np.bincount(flat, weights=w, minlength=nsub * nbin).reshape(nsub, nbin)
        subcnt += np.bincount(flat, minlength=nsub * nbin).reshape(nsub, nbin)
    return prof_sum / np.maximum(prof_cnt, 1), sub / np.maximum(subcnt, 1)


def profile_snr(prof):
    s = np.sort(prof)
    k = max(1, int(0.6 * prof.size))
    base = s[:k]
    off = base.std()
    return 0.0 if off == 0 else (prof.max() - np.median(base)) / off

def ridge_phase(subprof):
    """Peak phase (turns) per sub-integration from the first harmonic, and amplitude."""
    nsub, nbin = subprof.shape
    k = np.arange(nbin)
    c = (subprof * np.exp(-2j * np.pi * k / nbin)[None, :]).sum(axis=1)
    return -np.angle(c) / (2 * np.pi), np.abs(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["track", "poly", "kepler"], default="track")
    ap.add_argument("--dat", required=True)
    ap.add_argument("--inf", required=True)
    ap.add_argument("--track", required=True)
    ap.add_argument("--summary", default=None, help="fit_summary.dat (needed for poly/kepler)")
    ap.add_argument("--nbin", type=int, default=64)
    ap.add_argument("--nsub", type=int, default=32)
    ap.add_argument("--freq-col", type=int, default=None)
    ap.add_argument("--time-col", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--diagnose", action="store_true")
    ap.add_argument("--alignf0", action="store_true",
                    help="fit and remove the best constant frequency offset (F0 search)")
    ap.add_argument("--nsub-align", type=int, default=256,
                    help="fine sub-integrations used for the F0 alignment")
    a = ap.parse_args()

    tsamp, mjd = read_inf(a.inf)
    tmid, ftrack = read_track(a.track, a.freq_col, a.time_col)
    x = load_dat_trimmed(a.dat)
    Tobs = x.size * tsamp

    coeffs = []; tref = 0.0; kep = {}
    if a.model in ("poly", "kepler"):
        if a.summary is None:
            sys.exit("ERROR: --summary fit_summary.dat is required for poly/kepler")
        d = read_summary(a.summary)
        if a.model == "poly":
            order = int(d.get("poly_order_used", 5))
            coeffs = [d[f"poly_f{k}"] for k in range(order + 1)]
            tref = d["poly_tref_s"]
            phase_fn = make_poly_phase(coeffs, tref)
        else:
            kep = dict(f_spin=d["kepler_f_spin_hz"], K=d["kepler_K_hz"],
                       Pb=d["kepler_Pb_s"], T0=d["kepler_T0_s"])
            phase_fn = make_kepler_phase(kep["f_spin"], kep["K"], kep["Pb"], kep["T0"])
    else:
        phase_fn = make_track_phase(tmid, ftrack)

    # sanity check: does the model sit on the track?
    fm = model_freq(a.model, tmid, tmid, ftrack, coeffs, tref, kep)
    dfW = 1.0 / (tmid[1] - tmid[0])
    resid = (fm - ftrack).std()
    print(f"tsamp={tsamp:.6g}s  Tobs={Tobs:.1f}s  NT={ftrack.size}  Tcoh={(tmid[1]-tmid[0]):.2f}s")
    print(f"model={a.model}  model-vs-track RMS={resid:.3e} Hz ({resid/dfW:.2f} bins)")

    if a.alignf0:
        nsa = a.nsub_align - (a.nsub_align % a.nsub)   # multiple of nsub
        prof0, sub_fine = fold(x, tsamp, phase_fn, a.nbin, nsa)
        tk = (np.arange(nsa) + 0.5) * Tobs / nsa
        tref_a = Tobs / 2.0
        m = np.arange(a.nbin // 2 + 1)
        P = np.fft.rfft(sub_fine, axis=1)
        c1 = P[:, 1]
        grid = np.linspace(-8.0, 8.0, 16001) / Tobs   # +-8 turns
        score = np.abs((c1[None, :] * np.exp(2j * np.pi * grid[:, None] * (tk - tref_a)[None, :])).sum(1))
        df_best = grid[int(np.argmax(score))]
        Pshift = P * np.exp(2j * np.pi * m[None, :] * df_best * (tk - tref_a)[:, None])
        sub_al = np.fft.irfft(Pshift, n=a.nbin, axis=1)
        prof = sub_al.sum(0)
        fac = nsa // a.nsub
        subprof = sub_al[: a.nsub * fac].reshape(a.nsub, fac, a.nbin).sum(1)
        snr = profile_snr(prof)
        print(f"\n>>> {a.model.upper()} FOLD (F0-aligned) peak S/N = {snr:.2f}")
        print(f"    removed df = {df_best:.3e} Hz ({df_best*Tobs:.2f} turns); "
              f"unaligned S/N was {profile_snr(prof0):.2f}")
    else:
        prof, subprof = fold(x, tsamp, phase_fn, a.nbin, a.nsub)
        snr = profile_snr(prof)
        print(f"\n>>> {a.model.upper()} FOLD peak S/N = {snr:.2f}")

    diag = None
    if a.diagnose:
        f0 = ftrack[ftrack.size // 2] if a.model == "track" else \
             (coeffs[0] if a.model == "poly" else kep["f_spin"])
        tk = (np.arange(a.nsub) + 0.5) * Tobs / a.nsub
        ph, amp = ridge_phase(subprof)
        phu = np.unwrap(ph * 2 * np.pi) / (2 * np.pi)
        A = np.vstack([np.ones_like(tk), tk - Tobs / 2, (tk - Tobs / 2) ** 2]).T
        coef = np.linalg.lstsq(A * amp[:, None], phu * amp, rcond=None)[0]
        curv = phu - (coef[0] + coef[1] * (tk - Tobs / 2))
        print("\n--- drift diagnostic ---")
        print(f"linear residual freq  df = {coef[1]:.3e} Hz "
              f"({coef[1]*Tobs:.2f} turns over Tobs)")
        print(f"fractional offset     df/f0 = {coef[1]/f0:.3e} "
              f"(bary/topo v/c ~ 2.72e-5)")
        print(f"curvature RMS (orbital residual) = {curv.std():.3f} turns")
        diag = (tk, phu, coef)

    out = a.out or f"fold_on_track_{a.model}.pdf"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 6),
                                   gridspec_kw={"height_ratios": [1, 2]})
    ph = (np.arange(a.nbin) + 0.5) / a.nbin
    ax1.plot(np.concatenate([ph, ph + 1]), np.concatenate([prof, prof]), lw=1.0, color="k")
    ax1.set_xlim(0, 2); ax1.set_ylabel("flux"); ax1.set_xticklabels([])
    ax1.margins(y=0.15)
    ax1.text(0.02, 0.95, f"{a.model}  S/N = {snr:.1f}", transform=ax1.transAxes,
             va="top", ha="left")
    ax2.imshow(np.hstack([subprof, subprof]), aspect="auto", origin="lower",
               extent=[0, 2, 0, Tobs], cmap="viridis")
    if diag is not None:
        tk, phu, coef = diag
        fit = (coef[0] + coef[1] * (tk - Tobs / 2) + coef[2] * (tk - Tobs / 2) ** 2) % 1.0
        for off in (0.0, 1.0):
            ax2.plot((phu % 1.0) + off, tk, "r.", ms=4)
            ax2.plot(fit + off, tk, "w-", lw=1.0)
    ax2.set_xlabel("pulse phase"); ax2.set_ylabel("time (s)")
    fig.tight_layout(); fig.savefig(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()