#!/usr/bin/env python3
"""Detection threshold L_th via the exponential-tail method
(O'Leary et al. 2026, Eqs 13-17), for the num_harm = 1 search.

Noise model: each spectrogram bin is exponential(1) (normalised Fourier
power under Gaussian noise). For each of Nreal realisations of an
(Nbin x NT) spectrogram we run the same +/-1-bin Viterbi and keep the
final accumulated log-likelihood of every terminating bin, delta[:, -1].
Pooling all Nbin x Nreal values gives the noise-only p(L); we fit its
tail above the 99.99th percentile and invert for L_th at per-subband
false-alarm probability alpha'.
"""
import argparse
import numpy as np
import csv
import os




def viterbi_delta_final(spec):
    """Return delta[:, -1] for the +/-1-bin kernel (vectorised, delta only).
    Verified identical to the repo's viterbi() on random spectrograms."""
    Nf, NT = spec.shape
    d = spec[:, 0].copy()
    for t in range(1, NT):
        best = d.copy()                           # stay (f)
        best[1:] = np.maximum(best[1:], d[:-1])   # predecessor f-1
        best[:-1] = np.maximum(best[:-1], d[1:])  # predecessor f+1
        d = best + spec[:, t]
    return d


def noise_loglikes(NT, Nbin, Nreal, rng):
    out = np.empty(Nbin * Nreal)
    for i in range(Nreal):
        spec = rng.exponential(size=(Nbin, NT))
        out[i * Nbin:(i + 1) * Nbin] = viterbi_delta_final(spec)
    return out


def fit_Lth(L, Nbin, Nreal, alpha_prime, tail_pct):
    """Eqs 13-17. Returns dict with L_th, lambda, L_tail, k, N_tail."""
    L_tail = np.percentile(L, tail_pct)
    tail = L[L > L_tail]
    N_tail = tail.size
    lam = N_tail / np.sum(tail - L_tail)                      # Eq 14 (MLE)
    NQ = Nbin
    inner = Nreal * NQ * (1.0 - (1.0 - alpha_prime) ** (1.0 / NQ)) / N_tail
    L_th = L_tail - (1.0 / lam) * np.log(inner)               # Eq 17
    k = N_tail / (Nreal * NQ)
    return dict(L_th=L_th, lam=lam, L_tail=L_tail, k=k, N_tail=N_tail)


def bootstrap_sigma(L, Nbin, Nreal, alpha_prime, tail_pct, nboot, rng):
    Lr = L.reshape(Nreal, Nbin)
    vals = np.empty(nboot)
    for b in range(nboot):
        idx = rng.integers(0, Nreal, size=Nreal)
        vals[b] = fit_Lth(Lr[idx].ravel(), Nbin, Nreal, alpha_prime, tail_pct)["L_th"]
    return np.std(vals)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--NT", type=int, required=True, help="segments (time steps)")
    ap.add_argument("--Nbin", type=int, required=True,
                    help="frequency bins N_Q (= rows of that run's loglikes.dat)")
    ap.add_argument("--Nreal", type=int, default=10000)
    ap.add_argument("--alpha-prime", type=float, default=0.1)
    ap.add_argument("--tail-pct", type=float, default=99.99)
    ap.add_argument("--nboot", type=int, default=0, help="bootstrap resamples (0 = skip)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None, help="append one CSV row of results")
    ap.add_argument("--Tcoh", type=float, default=float("nan"), help="T_coh (s), for the output row")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    L = noise_loglikes(args.NT, args.Nbin, args.Nreal, rng)
    r = fit_Lth(L, args.Nbin, args.Nreal, args.alpha_prime, args.tail_pct)
    sig = bootstrap_sigma(L, args.Nbin, args.Nreal, args.alpha_prime,
                          args.tail_pct, args.nboot, rng) if args.nboot else float("nan")

    print(f"NT={args.NT}  Nbin(N_Q)={args.Nbin}  Nreal={args.Nreal}  alpha'={args.alpha_prime}")
    print(f"  L_tail (p{args.tail_pct}) = {r['L_tail']:.3f}   N_tail = {r['N_tail']}   k = {r['k']:.3e}")
    print(f"  lambda_MLE = {r['lam']:.4f}")
    print(f"  L_th = {r['L_th']:.3f}" + (f" +/- {sig:.3f}" if args.nboot else ""))
    if args.out:
        new = not os.path.exists(args.out)
        with open(args.out, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["N_T","T_coh_s","Nbin","alpha_prime","L_tail",
                            "lambda_MLE","k","N_tail","L_th","sigma_Lth"])
            w.writerow([args.NT, f"{args.Tcoh:.6f}", args.Nbin, args.alpha_prime,
                        f"{r['L_tail']:.4f}", f"{r['lam']:.4f}", f"{r['k']:.3e}",
                        r["N_tail"], f"{r['L_th']:.4f}", f"{sig:.4f}"])