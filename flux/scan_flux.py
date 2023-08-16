#!/usr/bin/env python

import subprocess
import numpy as np
from generate_data import generate_data, clean_data

process_script = "/fred/oz002/ldunn/viterbi_psr_search/simulatesearch_test/python_test/process_timeseries_het.py"
noise_script = "/fred/oz002/ldunn/viterbi_psr_search/simulatesearch_test/python_test/do_noise.py"

noise_params = {"tobs": 7200}
sig_params = {"width": 0.05, "flux": None}
parfile_params = {"freq": 173, "fdot": 0, "asini": 0.01, "orbitP": 4500, "ecc": 0}

noise_template = "./noise_template.params"
sig_template = "./sig_template.params"
par_template = "./template.par"

fluxes = np.logspace(np.log10(0.0002), np.log10(0.002), 10)

Tsft = 440

Nsft = int(7200//Tsft)

lls = []

for flux in fluxes:
    
    out_prefix = f"flux_{flux:.3E}"
    sig_params["flux"] = flux

    generate_data(sig_template, noise_template, par_template, sig_params, noise_params, parfile_params, out_prefix)

    process_cmd = f"python {process_script} --in {out_prefix}_data.dat --f0 150 --bw 50 --Tsft {Tsft} --Nsft {Nsft} --tsamp 6.4e-5"
    print(process_cmd)
    ll = float((subprocess.run(process_cmd, capture_output=True, text=True, shell=True).stdout).split('\n')[-2].split(":")[-1])
    lls.append(ll)

    clean_data(out_prefix)

print(lls)

np.savetxt("flux_results.dat", lls)
