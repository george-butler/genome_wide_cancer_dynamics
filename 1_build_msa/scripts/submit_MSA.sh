#!/bin/bash
#SBATCH --job-name=msa_array
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G

set -euo pipefail


GENE_LIST="$1"
INPUT_DIR="$2"

#SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/scripts"
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source activate MSAbuilder
# Where your per-gene pipeline writes outputs (adjust if your wrapper uses something else)
OUT_BASE="msa_out"

mkdir -p logs "${OUT_BASE}" bench

gene=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENE_LIST")
if [[ -z "${gene:-}" ]]; then
  echo "No gene for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

# Per-gene logs
OUT_LOG="logs/${gene}.out"
ERR_LOG="logs/${gene}.err"
exec >"$OUT_LOG" 2>"$ERR_LOG"

echo "[INFO] $(date -Iseconds) START gene=${gene}"
echo "[INFO] job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "[INFO] cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}"

export TMPDIR="${SLURM_TMPDIR:-/tmp}"

# Run your wrapper (NO nohup)
# If your wrapper accepts threads, pass $SLURM_CPUS_PER_TASK
# Otherwise it will still run fine.
set +e
srun --cpu-bind=cores \
  ${SCRIPTS_DIR}/submit_msa_and_tree.sh \
  "$gene" "$INPUT_DIR"
rc=$?
set -e

echo "[INFO] $(date -Iseconds) END gene=${gene} exit_code=${rc}"

# ---- Collect MACSE benchmarks (best-effort) ----
# This assumes MACSE was run with out_dir somewhere under the gene output dir
# and out_file_prefix == gene, so files look like:
#   ${gene}_benchmark.tsv, ${gene}_pipeline.log, ${gene}_bench_summary.txt
#
# If your wrapper writes them in a known dir, set MACSE_OUT_DIR explicitly.
#
# Heuristic search (cheap): look within OUT_BASE/gene for matching files.
gene_dir="${OUT_BASE}/${gene}"
if [[ -d "$gene_dir" ]]; then
  find "$gene_dir" -maxdepth 4 -type f \
    \( -name "${gene}_benchmark.tsv" -o -name "${gene}_pipeline.log" -o -name "${gene}_bench_summary.txt" \) \
    -exec cp -f {} bench/ \; || true
fi

# Track failures for easy retry
if [[ "$rc" -ne 0 ]]; then
  echo -e "${gene}\t${rc}" >> logs/failed_genes.tsv
fi

exit "$rc"

