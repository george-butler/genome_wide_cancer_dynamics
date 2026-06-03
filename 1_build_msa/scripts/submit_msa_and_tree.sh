#!/bin/bash
#SBATCH --job-name=msa_job        # Job name
#SBATCH --output=%x_%j.out             # Standard output and error log
#SBATCH --error=%x_%j.err              # Standard error log
#SBATCH --nodes=1                     # Number of tasks (e.g., cores)
#SBATCH --ntasks-per-node=10              # Number of CPU cores per task
#SBATCH --time=03:00:00                # Maximum runtime (e.g., 24 hours)
#SBATCH --mem=5G

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_gene> <input_dir> <ref_proteins_path>"
    exit 1
fi

# Infer launch directory
LAUNCH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set MACSE path
MACSE_PATH="${LAUNCH_DIR}/../MACSE_V2_PIPELINES"

# Check if MACSE exists
if [[ ! -d "${MACSE_PATH}" ]]; then
    echo "ERROR: MACSE_V2_PIPELINES not found at:"
    echo "  ${MACSE_PATH}"
    echo
    echo "Please install/download MACSE_V2_PIPELINES from:"
    echo "  https://github.com/ranwez/MACSE_V2_PIPELINES"
    exit 1
fi

# Input parameters
GENE=$1
INPUT_DIR=$2

REF_DIR=$3
DIRPATH="miniprot_all"

ALL_CDS=${INPUT_DIR}/all_combined_CDS.fasta
ALL_PROT=${INPUT_DIR}/all_combined_proteins.fasta
ALL_GENEIDS=${INPUT_DIR}/all_combined_gene_ids.txt
#REPORT=${INPUT_DIR}/all_gene_report.txt

OUTGROUPS="all_amphibians.txt" 

#output_dir=${INPUT_DIR}/MSA/${GENE}
output_dir=${INPUT_DIR}/${GENE}

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

input_fasta=$output_dir/input_proteins.fasta
input_cds_fasta=$output_dir/input_cds.fasta

if [[ ! -s ${input_fasta} ]] ; then
# New script
grep ":${GENE}$" ${ALL_GENEIDS} | grep -f ${OUTGROUPS} | sort | head -n 1 > $output_dir/outgroup.txt
grep ":${GENE}$" ${ALL_GENEIDS} | grep -v -f ${OUTGROUPS} | sort > $output_dir/geneids.txt

formatted_date=$(date "+%Y-%m-%d::%H:%M:%S")
echo "Step 1 : Making input files ... - $formatted_date"
# unordered
seqtk subseq ${ALL_CDS} $output_dir/outgroup.txt > ${input_cds_fasta}
seqtk subseq ${ALL_CDS} $output_dir/geneids.txt >> ${input_cds_fasta}

seqtk subseq ${ALL_PROT} $output_dir/outgroup.txt > ${input_fasta}
seqtk subseq ${ALL_PROT} $output_dir/geneids.txt >> ${input_fasta}

sed -i '/^>/ s/:/ /g' ${input_fasta}
sed -i '/^>/ s/:/ /g' ${input_cds_fasta}
fi

# Check if outgroup.txt is empty
if [[ ! -s $output_dir/outgroup.txt ]]; then
    #echo "Error: $output_dir/outgroup.txt is empty. Exiting."
    echo -e "::${GENE}\tNO_OUTGROUP"
    exit 1
fi

# Count the number of sequences in input_cds_fasta
num_sequences=$(grep -c "^>" "$input_fasta")

# Check if the number of sequences is less than 50
if [[ "$num_sequences" -lt 50 ]]; then
    #echo "Error: Number of sequences in $input_fasta is less than 50. Exiting."
    echo -e "::${GENE}\tLT_50"
    exit 1
fi


# Define output file names
msa_output="${output_dir}/ALIGN_${GENE}/${GENE}_final_unmask_align_AA.aln"
msa_output_fixed="${output_dir}/ALIGN_${GENE}/${GENE}_final_unmask_align_AA_fixed.aln"
iqtree_output="${output_dir}/iqtree"


rm -rf ${output_dir}/ALIGN_${GENE}/

formatted_date=$(date "+%Y-%m-%d::%H:%M:%S")

if [ ! -s "${msa_output}" ] ; then
# Step 1: Perform multiple sequence alignment using MAFFT
echo "Step 2 : Started Running  MSA... - $formatted_date"
${MACSE_PATH}/OMM_MACSE/S_OMM_MACSE_V12.01.sh --in_seq_file ${input_cds_fasta} --out_dir ${output_dir}/ALIGN_${GENE} --out_file_prefix ${GENE} --genetic_code_number 1 --min_percent_NT_at_ends 0.2
formatted_date=$(date "+%Y-%m-%d::%H:%M:%S")

echo "Complete! $formatted_date"
echo -e "::$GENE\tFINISHED\t$formatted_date"
fi

sed 's/[!,*]/-/g' ${output_dir}/ALIGN_${GENE}/${GENE}_final_homol_AA.aln > ${output_dir}/ALIGN_${GENE}/${GENE}_final_homol_AA_fixed.aln

ls ${output_dir}/ALIGN_${GENE}/* | grep -v "_export.aln\|_hmmCleaner_mask2_detail.aln\|_maskHomolog_stat.csv\|readme_output.txt\|*_fixed.aln"  | xargs -I % rm -f %

