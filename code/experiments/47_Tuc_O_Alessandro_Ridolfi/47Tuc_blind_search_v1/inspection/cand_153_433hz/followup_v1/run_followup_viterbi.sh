#!/bin/bash

CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"
VITERBI_DIR="${EXP}/inspection/cand_153_433hz/followup_v1/viterbi"
OUT_DIR="${EXP}/inspection/cand_153_433hz/followup_v1/viterbi_out"
mkdir -p ${OUT_DIR}

TSAMP=7.656074766355e-05
MJD_START=59391.125996527684038
F0=430.32
BW=6.0

DMS="25.13 25.18 25.23 25.28 25.33"
NTS="32 64 128 256"
TOBS=14395.4

for DM in ${DMS}; do
    DMTAG=$(printf "%.2f" ${DM})
    DAT="${VITERBI_DIR}/J0000-00_cfbf00000_Plan1_1_DM${DMTAG}.dat"

    for NT in ${NTS}; do
        TSFT=$(python3 -c "print(${TOBS}/${NT})")
        PREFIX="${OUT_DIR}/DM${DMTAG}_Nt${NT}/followup_DM${DMTAG}_Nt${NT}"
        PARAMS="${OUT_DIR}/DM${DMTAG}_Nt${NT}/followup_DM${DMTAG}_Nt${NT}.params"
        mkdir -p "${OUT_DIR}/DM${DMTAG}_Nt${NT}"

        cat > ${PARAMS} << EOF
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

        echo "Running DM=${DM} Nt=${NT} Tsft=${TSFT}..."
        python3 ${VITERBI} --params ${PARAMS}
        echo "Done DM=${DM} Nt=${NT}"
    done
done

echo "All done."