#!/usr/bin/env python

import subprocess
import os
import argparse

def populate_sig_params(template, width, flux):
    params = template
    params = params.replace("PULSEWIDTH", str(width))
    params = params.replace("FLUX", str(flux))
    return params.strip()

def populate_noise_params(template, tobs):
    params = template
    params = params.replace("TOBS", str(tobs))
    return params.strip()

def populate_parfile(template, freq, fdot, orbitP, asini, ecc, tasc):
    params = template
    params = params.replace("FREQ", str(freq))
    params = params.replace("FDOT", str(fdot))
    params = params.replace("ORBITP", str(orbitP))
    params = params.replace("ASINI", str(asini))
    params = params.replace("ECCENTRICITY", str(ecc))
    params = params.replace("TASC", str(tasc))
    return params.strip()

def generate_data(sig_template, noise_template, parfile_template,
        sig_params, noise_params, parfile_params, out_prefix):
    
    parfile = f"{out_prefix}.par"
    print(populate_parfile(open(parfile_template, "r").read(), **parfile_params), file=open(parfile, "w"))

    sig_paramfile = f"{out_prefix}_signal.params"
    print(populate_sig_params(open(sig_template, "r").read(), **sig_params), file=open(sig_paramfile, "w"))

    noise_paramfile = f"{out_prefix}_noise.params"
    print(populate_noise_params(open(noise_template, "r").read(), **noise_params), file=open(noise_paramfile, "w"))

    tempo2_cmd = f"tempo2 -f {parfile} -pred \"PKS 58456.3 58456.9 600 2400 16 16 900\""
    sys_noise_cmd = f"simulateSystemNoise -p {noise_paramfile} -o {out_prefix}_noise.dat"
    sim_psr_cmd = f"simulateComplexPsr -noDMsmear -p {noise_paramfile} -p {sig_paramfile} -o {out_prefix}_psr.dat"
    make_sf_cmd = f"createSearchFile -f {out_prefix}_psr.dat -f {out_prefix}_noise.dat -p {noise_paramfile} -o {out_prefix}_data.sf"

    prepdata_cmd = f"prepdata -o {out_prefix}_data -dm 2.64 {out_prefix}_data.sf"

    print(tempo2_cmd)
    subprocess.run(tempo2_cmd, shell=True)

    print(sys_noise_cmd)
    subprocess.run(sys_noise_cmd, shell=True)
    
    print(sim_psr_cmd)
    subprocess.run(sim_psr_cmd, shell=True)

    print(make_sf_cmd)
    subprocess.run(make_sf_cmd, shell=True)

    print(prepdata_cmd)
    subprocess.run(prepdata_cmd, shell=True)

def clean_data(out_prefix):
    os.remove(f"{out_prefix}_noise.dat")
    os.remove(f"{out_prefix}_psr.dat")
    os.remove(f"{out_prefix}_data.dat")
    os.remove(f"{out_prefix}_data.sf")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tobs', type=float, required=True, help="Total length of observation (s)")
    parser.add_argument('--pulse_width', required=True, type=float, help="Width of gaussian pulse (turns)")
    parser.add_argument('--flux', required=True, type=float, help='Flux density of pulsar (Jy)')
    parser.add_argument('--freq', required=True, type=float, help="Pulse frequency (Hz")
    parser.add_argument('--fdot', required=True, type=float, help="Pulsar spindown rate (Hz/s")
    parser.add_argument('--orbitP', required=True, type=float, help="Orbital period (s)")
    parser.add_argument('--asini', required=True, type=float, help="Projected semimajor axis (lt-s)")
    parser.add_argument('--ecc', required=True, type=float, default=0, help="Eccentricity parameter")

    parser.add_argument('--noise_template', required=True, type=str, help="Location of template noise parameter file")
    parser.add_argument('--sig_template', required=True, type=str, help="Location of template signal parameter file")
    parser.add_argument('--parfile_template', required=True, type=str, help="Location of template .par file")

    parser.add_argument('--out_prefix', required=True, type=str, help="Output prefix")

    args = parser.parse_args()

    sig_params = {"flux": args.flux, "width": args.pulse_width}
    noise_params = {"tobs": args.tobs}
    psr_params = {"freq": args.freq, "fdot": args.fdot, "orbitP": args.orbitP/86400, "asini": args.asini, "ecc": args.ecc}
    
    generate_data(args.sig_template, args.noise_template, args.parfile_template,
            sig_params, noise_params, psr_params, args.out_prefix)
