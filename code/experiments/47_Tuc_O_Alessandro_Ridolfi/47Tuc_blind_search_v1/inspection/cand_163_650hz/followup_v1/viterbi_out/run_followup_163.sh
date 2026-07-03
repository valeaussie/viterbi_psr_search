#!/bin/bash

CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"

source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh

TSAMP=7.65607476635514e-05
MJD_START=59391.128954238112783
TOBS=14395.4
F0=647.0
BW=6.0

DMS="25.10 25.20 25.30"
NTS="32 64 128 256"

for DM in ${DMS}; do
    DMTAG=$(printf "%.2f" ${DM})
    DAT="${EXP}/stage2_dedisp/DM${DMTAG}_bary/47Tuc_DM${DMTAG}_bary.dat"

    for NT in ${NTS}; do
        TSFT=$(python3 -c "print(${TOBS}/${NT})")
        OUTDIR="${EXP}/inspection/cand_163_650hz/followup_v1/viterbi_out/DM${DMTAG}_Nt${NT}"
        PREFIX="${OUTDIR}/followup_DM${DMTAG}_Nt${NT}"
        mkdir -p ${OUTDIR}

        cat > ${OUTDIR}/followup_DM${DMTAG}_Nt${NT}.params << EOF
infile         = "${DAT}"
tsamp          = ${TSAMP}
Tsft           = ${TSFT}
Nsft           = -1
f0             = ${F0}
bw             = ${BW}
padding_factor = 1
num_harm       = 1
top_paths      = 1
save_delta     = False
out_prefix     = "${PREFIX}"
plot_path      = True
dm             = ${DM}
mjd_start      = ${MJD_START}
EOF

        echo "Running DM=${DM} Nt=${NT}..."
        cd ${OUTDIR}
        python3 ${VITERBI} --params ${OUTDIR}/followup_DM${DMTAG}_Nt${NT}.params
        echo "Done DM=${DM} Nt=${NT}"
    done
done

echo "All done."