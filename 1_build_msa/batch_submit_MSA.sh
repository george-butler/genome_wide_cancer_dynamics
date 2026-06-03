#!/usr/bin/env bash
set -euo pipefail

OUTDIR="$1"

find "${OUTDIR}" -maxdepth 1 -type d -name "MSA_*" | sort -V | while read -r i
do
    genes="${i}/split_gene_list.txt"

    if [[ ! -s "${genes}" ]]; then
        echo "Missing or empty: ${genes}"
        continue
    fi

    N=$(wc -l < "${genes}")

    echo "Submitting ${i} (${N} genes)"

    sbatch --array=1-${N} scripts/submit_MSA.sh "${genes}" "${i}"
done
