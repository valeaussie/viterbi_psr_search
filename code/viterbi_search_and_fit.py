#!/usr/bin/env python

import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.fftpack as fft
from scipy.optimize import curve_fit
import argparse
import configparser
import ast
import os
from viterbi import viterbi, backtrace


"""
process_timeseries_het.py

This script performs a Viterbi-based search for a drifting narrowband signal
(e.g. a pulsar or continuous-wave signal) in a time-series dataset, using a
heteroscedastic (time-varying noise) spectrogram statistic.

The pipeline operates as follows:

    1. Load a time-domain signal from a binary file.
    2. Divide the time series into segments of length Tsft.
    3. For each segment, compute a Short Fourier Transform (SFT) in a specified
       frequency band, using heterodyning and optional resampling/decimation.
    4. Assemble these per-segment spectra into a time-frequency spectrogram.
    5. Apply a heteroscedastic normalisation / likelihood construction so that
       each timestep is weighted by its local noise properties (i.e. the noise
       variance may change with time).
    6. Apply the Viterbi algorithm to identify the most likely frequency path
       through the resulting (noise-weighted) spectrogram.
    7. Fit the recovered path with two independent models:
         (a) a polynomial Taylor expansion (default order 5 -> f0..f5) at the
             path midpoint,
         (b) a circular-orbit (Kepler) Doppler model giving (f_spin, K, Pb,
             T0).  Local f0..f5 are then evaluated analytically at the path
             midpoint from the fitted curve.
       Both fits are written to disk on every run, together with a summary
       flag indicating whether the path is consistent with a circular binary
       orbit or with an isolated source.
"""


def _load_params_file(path: str) -> dict:
    """Load key=value params file into a dict, converting types when possible."""
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"Params file not found: {path}")

    with open(path, "r") as f:
        raw = f.read()

    # configparser requires a section header
    raw = "[params]\n" + raw

    cfg = configparser.ConfigParser()
    cfg.optionxform = str  # preserve case sensitivity of keys
    cfg.read_string(raw)

    out = {}
    for k, v in cfg["params"].items():
        v = v.strip()
        try:
            out[k] = ast.literal_eval(v)
        except Exception:
            out[k] = v
    return out


def compute_sft(data, tsamp, Tsft, fh, bw, offset=0, padding_factor=1):
    print(len(data))
    count = int(Tsft / tsamp)

    data_to_sft = np.concatenate((
        (data[offset:offset + count] - np.mean(data[offset:offset + count])),
        np.zeros(int(count * (padding_factor - 1)))
    ))
    print(len(data_to_sft))
    ts = np.linspace(0, (len(data_to_sft) - 1) * tsamp, len(data_to_sft))
    het_ts = np.exp(-1j * 2 * np.pi * fh * ts)

    het_data = het_ts * data_to_sft

    new_tsamp = 1. / bw
    decim_factor = int(new_tsamp / tsamp)
    new_tsamp = tsamp * decim_factor

    het_decim_data = sig.resample(het_data, int(count * padding_factor / decim_factor))

    if len(het_decim_data) % 2 == 1:
        het_decim_data = het_decim_data[:-1]

    het_data_fft = fft.fft(het_decim_data)

    norm = np.var(het_decim_data) * len(het_decim_data)

    fs = np.linspace(-bw / 2, bw / 2, len(het_decim_data))
    return 1 / norm * np.abs(fft.fftshift(het_data_fft)) ** 2, fs


# ---------------------------------------------------------------------------
# Fitting helpers
# ---------------------------------------------------------------------------

def fit_polynomial_at_midpoint(t, path, order, df_bin):
    """
    Fit a polynomial of given order to the full Viterbi path, with the
    independent variable shifted to the path midpoint so that the returned
    coefficients are local Taylor coefficients at tref = t_midpoint.

    If the path has fewer points than (order + 2), the order is reduced
    automatically to len(path) - 2 to keep the fit overdetermined and to
    allow a covariance estimate.

    Returns
    -------
    tref : float
        Reference time (path midpoint), in seconds from the start of the
        dedispersed time series.
    coeffs_taylor : np.ndarray, shape (order_used + 1,)
        Taylor coefficients [f0, f1, f2, ..., f_{order_used}] at tref.  These
        satisfy
            f(t) = f0 + f1*(t-tref) + (1/2!)*f2*(t-tref)^2 + ...
    errs_taylor : np.ndarray, shape (order_used + 1,)
        1-sigma uncertainties on the Taylor coefficients.
    chi2, dof, red_chi2 : floats
        Goodness of fit using a quantisation uncertainty of df_bin per point.
    path_fit : np.ndarray
        Polynomial evaluated at the input times.
    order_used : int
        The polynomial order actually used (may be less than `order` if the
        path was too short).
    """
    npoints = len(path)
    order_used = int(order)
    if npoints < order_used + 2:
        order_used = max(0, npoints - 2)
        print(f"WARNING: only {npoints} path points; reducing polynomial "
              f"order from {order} to {order_used}.")

    imid = npoints // 2
    tref = float(t[imid])
    x = t - tref

    sigma_f = df_bin * np.ones_like(path)

    # np.polyfit returns highest-order first: [a_n, a_{n-1}, ..., a_1, a_0]
    # The model is y = a_n x^n + ... + a_1 x + a_0 (no factorial).
    # Taylor convention is f(t) = sum_k (1/k!) f_k (t-tref)^k, so
    #   f_k = k! * a_k.
    poly_coeffs, cov = np.polyfit(x, path, order_used, w=1.0 / sigma_f, cov=True)
    poly_errs = np.sqrt(np.diag(cov))

    # Reorder from highest-first to lowest-first, then convert to Taylor
    a_low_to_high = poly_coeffs[::-1]
    a_errs_low_to_high = poly_errs[::-1]

    coeffs_taylor = np.array([
        math.factorial(k) * a_low_to_high[k] for k in range(order_used + 1)
    ])
    errs_taylor = np.array([
        math.factorial(k) * a_errs_low_to_high[k] for k in range(order_used + 1)
    ])

    # Pad with NaNs if order_used < order so callers always get the same length
    if order_used < order:
        pad = np.full(order - order_used, np.nan)
        coeffs_taylor = np.concatenate([coeffs_taylor, pad])
        errs_taylor = np.concatenate([errs_taylor, pad])

    path_fit = np.polyval(poly_coeffs, x)
    resid = path - path_fit
    chi2 = float(np.sum((resid / sigma_f) ** 2))
    dof = npoints - (order_used + 1)
    red_chi2 = chi2 / dof if dof > 0 else float("nan")

    return tref, coeffs_taylor, errs_taylor, chi2, dof, red_chi2, path_fit, order_used


def kepler_model(t, f_spin, K, Pb, T0):
    """Circular-orbit apparent-frequency model.

    f_app(t) = f_spin - K * sin(2 * pi * (t - T0) / Pb)

    Parameters
    ----------
    f_spin : intrinsic spin frequency (Hz)
    K      : apparent-frequency semi-amplitude (Hz),
             K = f_spin * v_orb_los / c
    Pb     : orbital period (s)
    T0     : zero-crossing reference time (s, measured from the start of the
             dedispersed time series); not the same as the standard T_asc
             unless conventions are aligned externally
    """
    return f_spin - K * np.sin(2.0 * np.pi * (t - T0) / Pb)


def kepler_derivs(tref_val, f_spin, K, Pb, T0, order):
    """Analytic derivatives of the circular-orbit model at tref_val.

    Returns Taylor coefficients [f0, f1, f2, ..., f_order] satisfying
        f(t) = sum_k (1/k!) f_k (t - tref)^k

    For the model f(t) = f_spin - K sin(omega (t - T0)) with omega = 2*pi/Pb,
    the k-th derivative is
        f^{(k)}(t) = -K * omega^k * sin(omega (t - T0) + k*pi/2).

    The Taylor coefficient is just f^{(k)}(tref) directly (no factorial).
    """
    omega = 2.0 * np.pi / Pb
    phase = omega * (tref_val - T0)
    # k=0 is the function value itself
    f0_val = f_spin - K * np.sin(phase)
    out = [f0_val]
    for k in range(1, order + 1):
        out.append(-K * (omega ** k) * np.sin(phase + 0.5 * np.pi * k))
    return np.array(out)


def write_psrfold_candfile(out_path, dm,
                           poly_taylor, poly_tref_s, mjd_start,
                           kep_popt, kep_taylor_at_mid, tref_mid_s,
                           label):
    """
    Write a minimal psrfold_fil candfile (upstream 6-column format).

    The upstream psrfold_fil candfile schema is:
        #id dm acc F0 F1 S/N
    No F2, no orbital columns. Higher-order coefficients cannot be passed
    through the candfile and must be supplied via direct CLI flags
    (e.g. --f2) instead.

    Both rows pass apparent-frequency Taylor parameters (F0, F1) at their
    respective reference times:
        Row 0: polynomial fit at poly_tref_s
        Row 1: Kepler-derived Taylor expansion at tref_mid_s

    A sidecar .info file records the metadata (pepoch values, F2, label,
    Keplerian parameters) that cannot be embedded in the candfile.
    """
    SEC_PER_DAY = 86400.0

    # Polynomial row: apparent F0, F1 at poly_tref_s
    poly_F0 = float(poly_taylor[0])
    poly_F1 = (float(poly_taylor[1])
               if len(poly_taylor) > 1 and not np.isnan(poly_taylor[1])
               else 0.0)
    poly_F2 = (float(poly_taylor[2])
               if len(poly_taylor) > 2 and not np.isnan(poly_taylor[2])
               else 0.0)
    pepoch_poly_mjd = mjd_start + poly_tref_s / SEC_PER_DAY

    # Kepler row: apparent F0, F1 at tref_mid_s (NOT rest-frame f_spin)
    f_spin_fit, K_fit, Pb_s_fit, T0_s_fit = kep_popt
    kep_F0 = float(kep_taylor_at_mid[0])
    kep_F1 = float(kep_taylor_at_mid[1])
    kep_F2 = float(kep_taylor_at_mid[2])
    pepoch_kep_mjd = mjd_start + tref_mid_s / SEC_PER_DAY

    with open(out_path, "w") as fh:
        fh.write("#id dm acc F0 F1 S/N\n")
        fh.write(
            f"0 {dm:.6f} 0.0 {poly_F0:.12f} {poly_F1:.12e} 1.0\n"
        )
        fh.write(
            f"1 {dm:.6f} 0.0 {kep_F0:.12f} {kep_F1:.12e} 1.0\n"
        )

    # Sidecar metadata file
    info_path = out_path + ".info"
    with open(info_path, "w") as fh:
        fh.write("# Sidecar metadata for the psrfold_fil candfile\n")
        fh.write(f"# Verdict from automatic classification: {label}\n")
        fh.write("#\n")
        fh.write("# psrfold_fil --candfile takes only the 6-column upstream\n")
        fh.write("# schema (#id dm acc F0 F1 S/N). F2 and orbital parameters\n")
        fh.write("# are NOT read from the candfile. To fold with F2 (or\n")
        fh.write("# higher), drop --candfile and pass the parameters as\n")
        fh.write("# direct CLI flags using the values below.\n")
        fh.write("#\n")
        fh.write("# Row 0 (polynomial fit at path midpoint):\n")
        fh.write(f"#   --pepoch {pepoch_poly_mjd:.10f}\n")
        fh.write(f"#   --f0 {poly_F0:.12f}\n")
        fh.write(f"#   --f1 {poly_F1:.12e}\n")
        fh.write(f"#   --f2 {poly_F2:.12e}\n")
        fh.write("#\n")
        fh.write("# Row 1 (Kepler-derived Taylor expansion at path midpoint):\n")
        fh.write(f"#   --pepoch {pepoch_kep_mjd:.10f}\n")
        fh.write(f"#   --f0 {kep_F0:.12f}\n")
        fh.write(f"#   --f1 {kep_F1:.12e}\n")
        fh.write(f"#   --f2 {kep_F2:.12e}\n")
        fh.write("#\n")
        fh.write("# Kepler fit parameters (informational, not used for folding):\n")
        fh.write(f"#   f_spin = {f_spin_fit:.12f} Hz (rest-frame)\n")
        fh.write(f"#   K      = {K_fit:.6e} Hz\n")
        fh.write(f"#   Pb     = {Pb_s_fit:.6f} s ({Pb_s_fit / SEC_PER_DAY:.10f} d)\n")
        fh.write(f"#   T0     = {T0_s_fit:.6f} s from start of data\n")

def fit_kepler(t, path, df_bin, Tsft):
    """Fit a circular-orbit model to the full Viterbi path.

    To mitigate local-minimum issues when the observation contains multiple
    orbital cycles (so that several values of Pb give similar residuals), the
    fit is started from a small grid of Pb_guess values spanning a few
    decades, and the best result is kept.

    Returns popt, perr, chi2, dof, red_chi2, path_fit.
    """
    f_spin_guess = float(np.mean(path))
    K_guess = 0.5 * (np.max(path) - np.min(path))
    Tobs = t[-1] - t[0] + Tsft

    # Grid of period guesses spanning short orbits (Pb << Tobs) to long ones
    # (Pb >> Tobs).  At least one of these should put the local optimiser in
    # the basin of the true minimum.
    Pb_grid = np.geomspace(max(2.0 * Tsft, 100.0),
                           max(20.0 * Tobs, 100.0 * Tsft),
                           20)

    bounds = (
        [f_spin_guess - 0.1, 0.0, Tsft, -100.0 * Tobs],
        [f_spin_guess + 0.1,
         5.0 * (np.max(path) - np.min(path) + 1e-6),
         100.0 * Tobs,
         +100.0 * Tobs]
    )

    sigma_f_const = df_bin

    best = None
    for Pb_guess in Pb_grid:
        # Try a couple of phase guesses for T0 to cover the sin/cos ambiguity
        for phase_frac in (0.0, 0.25, 0.5, 0.75):
            T0_guess = float(t[int(np.argmax(path))]) - phase_frac * Pb_guess
            # Clip T0_guess into bounds
            T0_guess = float(np.clip(T0_guess, bounds[0][3] + 1.0,
                                     bounds[1][3] - 1.0))
            p0 = [f_spin_guess, K_guess, Pb_guess, T0_guess]
            try:
                popt, pcov = curve_fit(
                    kepler_model, t, path,
                    p0=p0, bounds=bounds, maxfev=20000
                )
                resid = path - kepler_model(t, *popt)
                chi2_try = float(np.sum((resid / sigma_f_const) ** 2))
                if best is None or chi2_try < best[2]:
                    best = (popt, pcov, chi2_try)
            except Exception:
                continue

    if best is None:
        # All starts failed; return a no-orbit guess
        print("WARNING: every Kepler-fit starting point failed.")
        popt = np.array([f_spin_guess, 0.0, Tobs, 0.0])
        perr = np.array([np.nan] * 4)
        path_fit = kepler_model(t, *popt)
    else:
        popt, pcov, _ = best
        perr = np.sqrt(np.diag(pcov))
        path_fit = kepler_model(t, *popt)

    sigma_f = df_bin * np.ones_like(path)
    resid = path - path_fit
    chi2 = float(np.sum((resid / sigma_f) ** 2))
    dof = len(path) - len(popt)
    red_chi2 = chi2 / dof if dof > 0 else float("nan")

    return popt, perr, chi2, dof, red_chi2, path_fit


def classify_orbit(K_fit, K_err, df_bin, red_chi2_kepler, red_chi2_poly):
    """Decide whether the path is consistent with a circular binary orbit.

    Heuristic. Returns (label, reasons) where label is one of
        "circular_orbit" or "isolated_or_unresolved".

    Conditions for "circular_orbit":
      (a) K is significantly larger than zero (K > 5 * K_err and K_err finite).
      (b) K is at least 2 bin widths (orbit actually resolved by the spectrogram).
      (c) The Kepler fit is acceptable (red_chi2_kepler < 5).
      (d) The Kepler fit is at least as good as the polynomial fit
          (red_chi2_kepler <= 1.5 * red_chi2_poly).  If a polynomial fits the
          path with similar quality, there is no compelling evidence for an
          orbit.
    """
    reasons = []
    cond_a = (np.isfinite(K_err) and K_err > 0 and K_fit > 5.0 * K_err)
    cond_b = (K_fit > 2.0 * df_bin)
    cond_c = (np.isfinite(red_chi2_kepler) and red_chi2_kepler < 5.0)
    cond_d = (
        np.isfinite(red_chi2_kepler) and np.isfinite(red_chi2_poly)
        and red_chi2_kepler <= 1.5 * red_chi2_poly
    )

    reasons.append(f"K = {K_fit:.3e} Hz (err = {K_err:.3e})")
    reasons.append(f"2 * bin width = {2.0 * df_bin:.3e} Hz")
    reasons.append(f"red_chi2_kepler = {red_chi2_kepler:.3f}")
    reasons.append(f"red_chi2_poly   = {red_chi2_poly:.3f}")
    reasons.append(f"K significant (>5 sigma): {cond_a}")
    reasons.append(f"K resolved (>2 bins):     {cond_b}")
    reasons.append(f"Kepler fit acceptable:    {cond_c}")
    reasons.append(f"Kepler not worse than poly: {cond_d}")

    if cond_a and cond_b and cond_c and cond_d:
        return "circular_orbit", reasons
    return "isolated_or_unresolved", reasons


def write_summary(out_path, label, reasons,
                  poly_tref, poly_taylor, poly_errs,
                  poly_chi2, poly_dof, poly_red_chi2,
                  kep_popt, kep_perr,
                  kep_chi2, kep_dof, kep_red_chi2,
                  kep_taylor_at_mid, tref_mid,
                  poly_order, poly_order_used):
    """Write a single human-readable summary file with both fits and a verdict."""
    f_spin_fit, K_fit, Pb_fit, T0_fit = kep_popt
    f_spin_err, K_err, Pb_err, T0_err = kep_perr

    if label == "circular_orbit":
        verdict_text = (
            "VERDICT: This pulsar is probably in a circular orbit. The Kepler\n"
            "fit is well constrained and describes the path better than the\n"
            "polynomial fit. The recommended fold parameters are the\n"
            "Kepler-derived Taylor expansion at tref_mid."
        )
    else:
        verdict_text = (
            "VERDICT: This pulsar is probably NOT in a circular orbit (or any\n"
            "orbit is not resolved by the current Tsft / observation). The\n"
            "polynomial fit at the path midpoint is the more appropriate\n"
            "source of fold parameters."
        )

    with open(out_path, "w") as fh:
        fh.write("# Viterbi-search fit summary\n")
        fh.write("# All Taylor coefficients f_k satisfy:\n")
        fh.write("#   f(t) = sum_k (1/k!) * f_k * (t - tref)^k\n")
        fh.write("# i.e. f_k = k-th derivative of f at tref.\n")
        fh.write("\n")
        fh.write(verdict_text + "\n\n")

        fh.write("# Diagnostic conditions\n")
        for r in reasons:
            fh.write(f"#   {r}\n")
        fh.write("\n")

        # Polynomial block
        fh.write("# ---------------------------------------------------------\n")
        fh.write(f"# Polynomial fit at path midpoint\n")
        fh.write(f"# Requested order: {poly_order}; order actually used: "
                 f"{poly_order_used}\n")
        fh.write("# ---------------------------------------------------------\n")
        fh.write(f"poly_order_requested {poly_order}\n")
        fh.write(f"poly_order_used {poly_order_used}\n")
        fh.write(f"poly_tref_s {poly_tref:.12f}\n")
        for k in range(poly_order + 1):
            if np.isnan(poly_taylor[k]):
                fh.write(f"poly_f{k} nan nan\n")
            else:
                fh.write(f"poly_f{k} {poly_taylor[k]:.15e} "
                         f"{poly_errs[k]:.6e}\n")
        fh.write(f"poly_chi2 {poly_chi2:.6f}\n")
        fh.write(f"poly_dof {poly_dof}\n")
        fh.write(f"poly_red_chi2 {poly_red_chi2:.6f}\n\n")

        # Kepler block
        fh.write("# ---------------------------------------------------------\n")
        fh.write("# Kepler circular-orbit fit\n")
        fh.write("# Model: f_app(t) = f_spin - K * sin(2*pi*(t - T0)/Pb)\n")
        fh.write("# ---------------------------------------------------------\n")
        fh.write(f"kepler_f_spin_hz {f_spin_fit:.15f} {f_spin_err:.6e}\n")
        fh.write(f"kepler_K_hz {K_fit:.15e} {K_err:.6e}\n")
        fh.write(f"kepler_Pb_s {Pb_fit:.6f} {Pb_err:.6e}\n")
        fh.write(f"kepler_T0_s {T0_fit:.6f} {T0_err:.6e}\n")
        fh.write(f"kepler_chi2 {kep_chi2:.6f}\n")
        fh.write(f"kepler_dof {kep_dof}\n")
        fh.write(f"kepler_red_chi2 {kep_red_chi2:.6f}\n")
        fh.write(f"# Kepler-derived Taylor expansion at path midpoint\n")
        fh.write(f"kepler_tref_s {tref_mid:.12f}\n")
        for k in range(poly_order + 1):
            fh.write(f"kepler_f{k} {kep_taylor_at_mid[k]:.15e}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str,
                        help="Input data file (PRESTO-format float32 .dat)")
    parser.add_argument('--tsamp', type=float,
                        help="Time series sampling period (in seconds)")
    parser.add_argument('--Tsft', type=float,
                        help="Length of each Short Fourier Transform (s)")
    parser.add_argument('--Nsft', type=int,
                        help="Number of SFTs to compute", required=False, default=-1)
    parser.add_argument('--f0', type=float,
                        help="Lower edge of the search band (Hz)")
    parser.add_argument('--padding-factor', type=float,
                        help="Zero-padding factor before each FFT", default=1)
    parser.add_argument('--bw', type=float,
                        help="Width of the search band (Hz)")
    parser.add_argument('--spec-flo', type=float,
                        help="Lowest frequency in spectrogram plot", required=False)
    parser.add_argument('--spec-fhi', type=float,
                        help="Highest frequency in spectrogram plot", required=False)
    parser.add_argument('--plot-path',
                        help="Overplot best Viterbi path on spectrogram",
                        action='store_true')
    parser.add_argument('--no-plot-path',
                        help="Do not overplot best Viterbi path on spectrogram",
                        action='store_false', dest='plot_path')
    parser.add_argument('--out_prefix', type=str,
                        help="Prefix for output files", default="search")
    parser.add_argument('--top_paths', type=int,
                        help="Save the top N Viterbi paths and log likelihoods",
                        required=False, default=1)
    parser.add_argument('--save-delta', action='store_true', dest='save_delta')
    parser.add_argument('--num-harm', type=int,
                        help="Number of harmonics to sum incoherently", default=1)
    parser.add_argument('--poly-order', type=int,
                        help="Polynomial order for Taylor-expansion fit "
                             "(default 5: gives f0..f5)", default=5)
    parser.add_argument('--params', type=str, default=None,
                        help="Path to .params file. CLI args override file values.")

    pre_args, _ = parser.parse_known_args()
    params_dict = _load_params_file(pre_args.params)
    parser.set_defaults(**params_dict)
    args = parser.parse_args()

    required_fields = ['infile', 'tsamp', 'Tsft', 'f0', 'bw']
    missing = [field for field in required_fields if getattr(args, field) is None]
    if missing:
        raise ValueError(f"Missing required parameters: {missing}")

    tsamp = args.tsamp
    data = np.fromfile(args.infile, dtype=np.single)
    Tfull = len(data) * tsamp
    Tsft = args.Tsft
    if args.Nsft == -1:
        Nsft = int(Tfull // Tsft)
    elif args.Nsft > 0:
        Nsft = args.Nsft
    else:
        print("--Nsft must be > 0!")
        sys.exit(1)

    downsample_tsamp = tsamp
    downsample_data = data

    bw = args.bw
    f0 = args.f0
    poly_order = int(args.poly_order)

    # ------------------------------------------------------------------
    # Build spectrogram
    # ------------------------------------------------------------------
    _, fs = compute_sft(downsample_data, downsample_tsamp, Tsft,
                        f0 + bw / 2., bw, padding_factor=args.padding_factor)
    spec = np.zeros((len(fs), Nsft), dtype=float)
    unsummed_spec = np.zeros((len(fs), args.num_harm), dtype=float)
    for i in range(0, Nsft):
        print(f"Doing timestep {i + 1} of {Nsft}")
        for harm in range(1, args.num_harm + 1):
            harm_spec, _ = np.abs(compute_sft(
                downsample_data, downsample_tsamp, Tsft,
                (f0 + bw / 2) * harm, bw * harm,
                offset=int(i * Tsft / downsample_tsamp),
                padding_factor=args.padding_factor
            ))
            harm_spec = harm_spec[:harm * len(fs)]
            avgd_harm_spec = np.mean(harm_spec.reshape((-1, harm)), axis=1)
            unsummed_spec[:, harm - 1] = avgd_harm_spec
        spec[:, i] = np.sum(unsummed_spec, axis=1)

    # ------------------------------------------------------------------
    # Run Viterbi
    # ------------------------------------------------------------------
    delta, backptrs = viterbi(spec)

    if args.save_delta:
        np.savetxt(f"{args.out_prefix}_delta.dat", delta)

    top_n_delta_idxs = np.argsort(delta[:, -1])[-args.top_paths:]
    with open(f"{args.out_prefix}_paths.dat", "w") as f:
        for idx in top_n_delta_idxs:
            path_top = " ".join([str(x)
                                 for x in f0 + bw / 2 + fs[backtrace(backptrs, idx)]])
            print(f"{delta[idx, -1]} {path_top}", file=f)

    path = f0 + bw / 2 + fs[backtrace(backptrs, np.argmax(delta[:, -1]))]
    t = Tsft / 2. + np.arange(Nsft) * Tsft
    df_bin = float(fs[1] - fs[0])

    print("Recovered Viterbi path:")
    print(path)

    # ------------------------------------------------------------------
    # Polynomial fit at the path midpoint
    # ------------------------------------------------------------------
    (poly_tref, poly_taylor, poly_errs,
     poly_chi2, poly_dof, poly_red_chi2,
     poly_path_fit, poly_order_used) = fit_polynomial_at_midpoint(
        t, path, poly_order, df_bin
    )

    print(f"\nPolynomial fit (order {poly_order_used}) at tref = "
          f"{poly_tref:.3f} s (path midpoint):")
    for k in range(poly_order + 1):
        if np.isnan(poly_taylor[k]):
            print(f"  f{k} = (not fitted, polynomial order too low)")
        else:
            print(f"  f{k} = {poly_taylor[k]:.12e}  +/- {poly_errs[k]:.3e}")
    print(f"  chi2 = {poly_chi2:.3f}, dof = {poly_dof}, "
          f"reduced chi2 = {poly_red_chi2:.3f}")

    # ------------------------------------------------------------------
    # Kepler fit (always run)
    # ------------------------------------------------------------------
    (kep_popt, kep_perr, kep_chi2, kep_dof, kep_red_chi2,
     kep_path_fit) = fit_kepler(t, path, df_bin, Tsft)
    f_spin_fit, K_fit, Pb_fit, T0_fit = kep_popt
    f_spin_err, K_err, Pb_err, T0_err = kep_perr

    # Kepler-derived Taylor expansion at path midpoint
    tref_mid = float(t[len(path) // 2])
    kep_taylor_at_mid = kepler_derivs(
        tref_mid, f_spin_fit, K_fit, Pb_fit, T0_fit, poly_order
    )

    print("\nKepler-fit orbital parameters:")
    print(f"  f_spin  = {f_spin_fit:.12f} +/- {f_spin_err:.3e} Hz")
    print(f"  K       = {K_fit:.6e}    +/- {K_err:.3e} Hz")
    print(f"  Pb      = {Pb_fit:.3f}   +/- {Pb_err:.3e} s")
    print(f"  T0      = {T0_fit:.3f}   +/- {T0_err:.3e} s "
          f"(from start of dedispersed time series)")
    print(f"  chi2    = {kep_chi2:.3f}, dof = {kep_dof}, "
          f"reduced chi2 = {kep_red_chi2:.3f}")
    print(f"\nKepler-derived Taylor expansion at tref = {tref_mid:.3f} s "
          f"(path midpoint):")
    for k in range(poly_order + 1):
        print(f"  f{k} = {kep_taylor_at_mid[k]:.12e}")

    # ------------------------------------------------------------------
    # Decide which fit to trust and write summary
    # ------------------------------------------------------------------
    label, reasons = classify_orbit(
        K_fit, K_err, df_bin, kep_red_chi2, poly_red_chi2
    )
    print(f"\nClassification: {label}")
    for r in reasons:
        print(f"  {r}")

    write_summary(
        f"{args.out_prefix}_fit_summary.dat",
        label, reasons,
        poly_tref, poly_taylor, poly_errs,
        poly_chi2, poly_dof, poly_red_chi2,
        kep_popt, kep_perr,
        kep_chi2, kep_dof, kep_red_chi2,
        kep_taylor_at_mid, tref_mid,
        poly_order, poly_order_used,
    )

# Write a psrfold_fil candfile with both candidates
    # NOTE: requires --dm and --mjd-start as inputs.  If not provided,
    # the candfile is skipped.
    if args.dm is not None and args.mjd_start is not None:
        write_psrfold_candfile(
            f"{args.out_prefix}_psrfold.candfile",
            dm=args.dm,
            poly_taylor=poly_taylor,
            poly_tref_s=poly_tref,
            mjd_start=args.mjd_start,
            kep_popt=kep_popt,
            kep_taylor_at_mid=kep_taylor_at_mid,
            tref_mid_s=tref_mid,
            label=label,
        )
        print(f"\nWrote {args.out_prefix}_psrfold.candfile (--pepoch values "
              f"are listed as comments inside the file).")
    else:
        print("\nSkipping psrfold candfile (need --dm and --mjd-start).")


    # Save the Viterbi track plus both fitted curves
    np.savetxt(
        f"{args.out_prefix}_track.dat",
        np.column_stack([t, path, poly_path_fit, kep_path_fit]),
        header="time_s viterbi_frequency_hz polynomial_fit_hz kepler_fit_hz"
    )

    # ------------------------------------------------------------------
    # Plot spectrogram with Viterbi path + both fits overlaid
    # ------------------------------------------------------------------
    fig, ax = plt.subplots()
    ax.imshow(spec, aspect='auto', extent=[0, Tfull, f0 + bw, f0], cmap='gray_r')
    ax.autoscale(False)

    if args.plot_path:
        ax.plot(t, path, 'r', linewidth=2, label='Viterbi path')
        ax.plot(t, poly_path_fit, color='tab:blue', linestyle='--',
                linewidth=1.5, label=f'Polynomial (order {poly_order})')
        ax.plot(t, kep_path_fit, color='tab:orange', linestyle=':',
                linewidth=1.5, label='Kepler fit')
        ax.legend(loc='best', fontsize=8)

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Frequency (Hz)')
    if args.spec_flo is not None and args.spec_fhi is not None:
        ax.set_ylim(args.spec_flo, args.spec_fhi)
    plt.savefig(f'{args.out_prefix}_spectrogram.png', dpi=600,
                bbox_inches='tight')
    plt.clf()

    # Log-likelihood vs frequency at the final time step
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Log likelihood')
    plt.plot(f0 + bw / 2 + fs, delta[:, -1])
    plt.savefig(f'{args.out_prefix}_loglikes.png', dpi=600,
                bbox_inches='tight')
    plt.clf()

    print(f"\nLoglike: {np.max(delta[:, -1])}")
