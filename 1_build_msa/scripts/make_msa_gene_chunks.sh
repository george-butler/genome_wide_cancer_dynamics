#!/usr/bin/env bash
set -euo pipefail

HUMAN_PROTEINS="$1"              # input fasta
OUTDIR="$2"                      # output directory
CHUNK_SIZE="${3:-1000}"          # genes per split dir

mkdir -p "${OUTDIR}"

# Output files
GENE_SIZES="${OUTDIR}/human_proteins.gene_sizes.tsv"
SORTED_GENES="${OUTDIR}/human_proteins.sorted_by_size.gene_list.txt"

###############################################################################
# STEP 1: Get protein lengths and sort by size
###############################################################################

awk '
BEGIN{OFS="\t"}

/^>/{
    if(id!="")
        print id, len

    id=$1
    sub(/^>/,"",id)
    len=0
    next
}

{
    gsub(/[ \t\r\n]/,"")
    len += length($0)
}

END{
    if(id!="")
        print id, len
}
' "${HUMAN_PROTEINS}" \
| sort -k2,2nr > "${GENE_SIZES}"

###############################################################################
# STEP 2: Extract ordered gene list
###############################################################################

cut -f1 "${GENE_SIZES}" > "${SORTED_GENES}"

###############################################################################
# STEP 3: Split into MSA_1 ... MSA_N directories
###############################################################################

awk -v chunk="${CHUNK_SIZE}" -v outdir="${OUTDIR}" '
{
    batch = int((NR-1)/chunk) + 1

    dir = outdir "/MSA_" batch
    file = dir "/split_gene_list.txt"

    if (!(dir in made)) {
        system("mkdir -p " dir)
        made[dir]=1
    }

    print $0 >> file
}
' "${SORTED_GENES}"

###############################################################################
# DONE
###############################################################################

echo "Done."
echo "Gene sizes: ${GENE_SIZES}"
echo "Sorted gene list: ${SORTED_GENES}"
echo "MSA directories created under: ${OUTDIR}"
