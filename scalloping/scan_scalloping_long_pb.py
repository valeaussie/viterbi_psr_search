#!/usr/bin/env python

import subprocess
import numpy as np
from generate_data import generate_data, clean_data

process_script = "/fred/oz002/ldunn/viterbi_psr_search/paper/code/process_timeseries_het.py"
noise_script = "/fred/oz002/ldunn/viterbi_psr_search/simulatesearch_test/python_test/do_noise.py"

noise_params = {"tobs": 7200}
sig_params = {"width": 0.05, "flux": 0.001}
parfile_params = {"freq": 173, "fdot": 0, "asini": 0.01, "orbitP": 4500, "ecc": 0, "tasc": 54501.467101261028311}

noise_template = "./noise_template.params"
sig_template = "./sig_template.params"
par_template = "./template.par"


Tsft = 440

freq_offsets = np.linspace(0, 1.0/Tsft, 11)

Nsft = int(7200//Tsft)

lls = []

for freq_offset in freq_offsets:
    
    out_prefix = f"freq_offset_{freq_offset:.3E}"
    parfile_params['freq'] = 173 + freq_offset
    parfile_params['tasc'] += np.random.uniform(-1,1)*parfile_params['orbitP']/86400

    generate_data(sig_template, noise_template, par_template, sig_params, noise_params, parfile_params, out_prefix)

    process_cmd = f"python {process_script} --in {out_prefix}_data.dat --f0 150 --bw 50 --Tsft {Tsft} --Nsft {Nsft} --tsamp 6.4e-5 --padding-factor 2"
    print(process_cmd)
    ll = float((subprocess.run(process_cmd, capture_output=True, text=True, shell=True).stdout).split('\n')[-2].split(":")[-1])
    lls.append(ll)

    clean_data(out_prefix)

print(lls)

np.savetxt("scalloping_results.dat", lls)
