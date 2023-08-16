#!/usr/bin/env python

import subprocess
import numpy as np
from generate_data import generate_data, clean_data
import itertools
import pathlib
import sys

noise_params = {"tobs": 7200}
sig_params = {"width": 0.05, "flux": 0.001}
parfile_params = {"freq": 173, "fdot": 0, "asini": 0.01, "orbitP": 4500, "ecc": 0, "tasc": 54501.467101261028311}

noise_template = "./noise_template.params"
sig_template = "./sig_template.params"
par_template = "./template.par"


tobss = [1800, 3600, 7200, 14400]
porb_fracs = [10, 5, 2, 1, 0.5]
fluxes = np.logspace(np.log10(0.0002), np.log10(0.002), 5)
eccentricities = [0, 0.01, 0.1, 0.5]

paramsets = list(itertools.product(tobss, porb_fracs, fluxes, eccentricities))
print(len(paramsets))
tobs, porb_frac, flux, ecc = paramsets[int(sys.argv[1])]

noise_params['tobs'] = tobs
sig_params['flux'] = flux
parfile_params['orbitP'] = tobs*porb_frac/86400
parfile_params['ecc'] = ecc

out_id = f"{tobs}_{porb_frac}_{ecc}_{flux:.3E}"
out_dir = f"data_{out_id}"
pathlib.Path(f"data_{out_id}").mkdir(exist_ok=True)
out_prefix = f"{out_dir}/{out_id}"
print(f"Doing {out_prefix}")
generate_data(sig_template, noise_template, par_template, sig_params, noise_params, parfile_params, out_prefix)

clean_data(out_prefix, leave_ts=True)
