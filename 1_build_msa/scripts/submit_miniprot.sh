#!/bin/bash
#SBATCH --job-name=miniprot_job        # Job name
#SBATCH --output=%x_%j.out             # Standard output and error log
#SBATCH --error=%x_%j.err              # Standard error log
#SBATCH --nodes=1                     # Number of tasks (e.g., cores)
#SBATCH --ntasks-per-node=4              # Number of CPU cores per task
#SBATCH --time=05:00:00                # Maximum runtime (e.g., 24 hours)

# Parse command-line arguments
query_fa=$1    # FASTA file with query sequences
db=$2          # Database file
outdir=$3      # Output directory

#SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/scripts"
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source activate MSAbuilder

# Ensure output directory exists
mkdir -p $outdir
awk '$1 ~ /^>/ { print $1}' ${db} |  sed 's/>//' > $outdir/gene_list.txt

GENE=$( basename ${db} ".fa" | awk -F"_" '{ print $2 }' )

if [[ ${query_fa} != *VGP/* ]]; then
    SP=$( dirname ${query_fa} | awk -F"/" '{ split($(NF-3),p,"_") ; g=toupper(substr(p[1],1,1)); s=substr(p[2],1,3) ; print $(NF-3)"_"$NF}' )
else
    SP=$( dirname ${query_fa} | awk -F"/" '{ split($(NF),p,"_") ; g=toupper(substr(p[1],1,1)); s=substr(p[2],1,3) ; print $NF"_"$(NF-1)}' )
fi

OUT=aln_miniprot_genes.gff
OUT_FIXED=aln_miniprot_genes_fixed.gff
OUT_STATS=aln_miniprot_genes_fixed_stats.txt

GFF=$outdir/all_genes.gff
FA=${query_fa}

# Run Miniprot
if [ ! -s $outdir/${OUT} ]; then
echo "Step 1 : Running miniprot"
miniprot -t 8 --gff $query_fa ${db} > $outdir/${OUT}
fi

# Print completion message
echo "INFO : Miniprot job completed. Output saved to $outdir/${OUT}"

if [ ! -s $outdir/all_genes.CDS.fasta ]; then
grep -v "#PAF" $outdir/${OUT} > $outdir/${OUT_FIXED} 

awk -F"\t" '$3 == "mRNA" && !($9 ~ /frameshift/) { split($9,k,";"); match($9,"Identity=([^;]+)"); idty=substr($9,RSTART+9,RLENGTH-9); match($9,"ID=([^;]+)"); mid=substr($9,RSTART+3,RLENGTH-3); match($9,"Rank=([^;]+)"); rank=substr($9,RSTART+5,RLENGTH-5); match($9,"Positive=([^;]+)"); pstv=substr($9,RSTART+9,RLENGTH-9); match($9,"Target=([^$]+)"); target=substr($9,RSTART+7,RLENGTH-7); split(target,gene," "); f1=(2 * pstv * idty ) / ( pstv + idty ) ; print gene[1]"\t"idty"\t"pstv"\t"rank"\t"f1"\t"mid }' $outdir/${OUT_FIXED} | awk -F"\t" '$2 >= 0.75 && $3 > 0.75 { print }' | sort -k1,1 -V -k4,4n > $outdir/${OUT_STATS}

awk -F'\t' '{if ($1 in min) {if ($4 < min[$1]) {min[$1] = $4; last[$1] = $NF}} else {min[$1] = $4; last[$1] = $NF}} END {for (key in min) print key, min[key], last[key]}' $outdir/${OUT_STATS} | awk '{ print $1"\t"$NF }' > $outdir/geneid_mapping.txt

awk -F"\t" '{ print $NF }' $outdir/geneid_mapping.txt > $outdir/all_geneids.txt

grep -f $outdir/all_geneids.txt $outdir/${OUT_FIXED} > $outdir/all_genes.gff

awk 'NR==FNR {genes[$1]; next} $3 == "mRNA" {split($9, a, ";"); for (i in a) if (a[i] ~ /Target=/) {found[substr(a[i], 8)] = 1}} END {found_count=0; not_found_count=0; for (gene in genes) {if (gene in found) {found_count++} else {not_found[not_found_count++] = gene}}; print "Total genes found:"found_count; print "Total genes not found:", not_found_count; if (not_found_count > 0) {print "Genes not found:"; for (j = 0; j < not_found_count; j++) print not_found[j]}}' $outdir/gene_list.txt $outdir/all_genes.gff > $outdir/miniprot_report.txt 

sed -i 's/ /_/g' $outdir/all_genes.gff

agat_sp_extract_sequences.pl -gff ${GFF} -f ${FA} --cis --cfs -p -o $outdir/all_genes.proteins_tmp.fasta
agat_sp_extract_sequences.pl -gff ${GFF} -f ${FA} -p -o $outdir/all_genes.proteins_tmp_with_start_stop.fasta
agat_sp_extract_sequences.pl -gff ${GFF} -f ${FA} -t cds -o $outdir/all_genes.CDS_tmp.fasta

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {map[$2]=$1; next} {if ($1 ~ /^>/) {split($1, a, " "); header=a[1]; header=substr(header, 2); if (header in map) {print ">" map[header]} else {print $1}} else {print}}' $outdir/geneid_mapping.txt $outdir/all_genes.proteins_tmp.fasta > $outdir/all_genes.proteins.fasta

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {map[$2]=$1; next} {if ($1 ~ /^>/) {split($1, a, " "); header=a[1]; header=substr(header, 2); if (header in map) {print ">" map[header]} else {print $1}} else {print}}' $outdir/geneid_mapping.txt $outdir/all_genes.CDS_tmp.fasta > $outdir/all_genes.CDS.fasta


sed -i "s/^>/>${SP} /g" $outdir/all_genes.proteins.fasta
sed -i "s/^>/>${SP} /g" $outdir/all_genes.CDS.fasta

rm -f $outdir/all_genes.CDS_tmp.fasta
rm -f $outdir/all_genes.proteins_tmp.fasta

fi

if [ ! -s $outdir/all_genes.proteins_tmp_with_start_stop.fasta ] ; then
   sed -i 's/ /_/g' $outdir/all_genes.gff
   agat_sp_extract_sequences.pl -gff ${GFF} -f ${FA} -p -o $outdir/all_genes.proteins_tmp_with_start_stop.fasta
fi

if [ ! -s ${outdir}/psuedogene_report_ori.txt ] ; then
   ${SCRIPTS_DIR}/find_psuedogenes.sh ${outdir}/all_genes.proteins_tmp_with_start_stop.fasta ${SP}
   mv ${outdir}/psuedogene_report.txt ${outdir}/psuedogene_report_ori.txt
fi

# Added updates to make sure a single query gene loci is attached to only one target ( non overlapping filter ) 
REMOVED_TSV=$outdir/genes_removed.tsv
SUMMARY_TSV=$outdir/genes_summary.tsv
OUT_LATEST_GFF=$outdir/aln_miniprot_genes_final.gff

if [ ! -s $outdir/all_genes.proteins_latest.fasta ]; then
${SCRIPTS_DIR}/filter_miniprot_best_per_target.py --min_identity 0.75 --min_positive 0.75 --removed_tsv ${REMOVED_TSV} --summary_tsv ${SUMMARY_TSV} $outdir/${OUT_FIXED} ${OUT_LATEST_GFF}

sed -i 's/Target=\([^;]*\)/Target=\1/; s/ /_/g' ${OUT_LATEST_GFF}

awk -F"\t" '$11 ~ /kept/ { print $1"\t"$2 }' ${SUMMARY_TSV} | sort | uniq > $outdir/geneid_mapping_latest.txt

agat_sp_extract_sequences.pl -gff ${OUT_LATEST_GFF} -f ${FA} --cis --cfs -p -o $outdir/all_genes.proteins_latest_tmp.fasta
agat_sp_extract_sequences.pl -gff ${OUT_LATEST_GFF} -f ${FA} -p -o $outdir/all_genes.proteins_with_start_stop_latest_tmp.fasta
agat_sp_extract_sequences.pl -gff ${OUT_LATEST_GFF} -f ${FA} -t cds -o $outdir/all_genes.CDS_latest_tmp.fasta

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {map[$2]=$1; next} {if ($1 ~ /^>/) {split($1, a, " "); header=a[1]; header=substr(header, 2); if (header in map) {print ">" map[header]} else {print $1}} else {print}}' $outdir/geneid_mapping_latest.txt $outdir/all_genes.proteins_latest_tmp.fasta > $outdir/all_genes.proteins_latest.fasta

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {map[$2]=$1; next} {if ($1 ~ /^>/) {split($1, a, " "); header=a[1]; header=substr(header, 2); if (header in map) {print ">" map[header]} else {print $1}} else {print}}' $outdir/geneid_mapping_latest.txt $outdir/all_genes.CDS_latest_tmp.fasta > $outdir/all_genes.CDS_latest.fasta

sed -i "s/^>/>${SP}:/g" $outdir/all_genes.proteins_latest.fasta
sed -i "s/^>/>${SP}:/g" $outdir/all_genes.CDS_latest.fasta

rm -f $outdir/all_genes.CDS_latest_tmp.fasta
rm -f $outdir/all_genes.proteins_latest_tmp.fasta
fi 

[ -f $outdir/psuedogene_report.txt ] && mv $outdir/psuedogene_report.txt $outdir/psuedogene_report.txt.bkp_$(date +%Y%m%d_%H%M%S)

if [ ! -s ${outdir}/psuedogene_report_latest.txt ] ;then
# Find psuedogenes
${SCRIPTS_DIR}/find_psuedogenes.sh ${outdir}/all_genes.proteins_with_start_stop_latest_tmp.fasta ${SP}

cp ${outdir}/psuedogene_report.txt ${outdir}/psuedogene_report_latest.txt

mkdir -p ${outdir}/changed_mapping
awk -F"\t" 'NR==FNR { a[$1] = $2 ;next } ($1 in a ) &&  ( $2 != a[$1] ) { print }' $outdir/geneid_mapping.txt $outdir/geneid_mapping_latest.txt > $outdir/changed_mapping/changed_mapping_latest.txt

awk -F"\t" 'NR==FNR { a[$1] = $2 ;next } ($1 in a ) &&  ( $2 != a[$1] ) { print }' $outdir/geneid_mapping_latest.txt $outdir/geneid_mapping.txt > $outdir/changed_mapping/changed_mapping_ori.txt

grep -f <( awk -F"\t" '{ print $2}' $outdir/changed_mapping/changed_mapping_latest.txt ) ${OUT_LATEST_GFF} > ${outdir}/changed_mapping/changed_mapping.gff

seqtk subseq $outdir/all_genes.proteins_tmp_with_start_stop.fasta <( awk -F"\t" '{ print $2}' $outdir/changed_mapping/changed_mapping_ori.txt ) > ${outdir}/changed_mapping/changed_mapping_ori.proteins.fasta
seqtk subseq $outdir/all_genes.proteins_with_start_stop_latest_tmp.fasta <( awk -F"\t" '{ print $2}' $outdir/changed_mapping/changed_mapping_latest.txt ) > ${outdir}/changed_mapping/changed_mapping_latest.proteins.fasta

${SCRIPTS_DIR}/find_psuedogenes.sh ${outdir}/changed_mapping/changed_mapping_latest.proteins.fasta ${SP}_changed_latest

mv ${outdir}/changed_mapping/psuedogene_report.txt ${outdir}/changed_mapping/changed_psuedogene_report_latest.txt

${SCRIPTS_DIR}/find_psuedogenes.sh ${outdir}/changed_mapping/changed_mapping_ori.proteins.fasta ${SP}_changed_ori

mv ${outdir}/changed_mapping/psuedogene_report.txt ${outdir}/changed_mapping/changed_psuedogene_report_ori.txt

{
  head -n 1 ${outdir}/changed_mapping/changed_psuedogene_report_latest.txt
  {
    tail -n +2 ${outdir}/changed_mapping/changed_psuedogene_report_latest.txt
    tail -n +2 ${outdir}/changed_mapping/changed_psuedogene_report_ori.txt
  } | sort -u
} > ${outdir}/changed_mapping/psuedogene_report.txt

fi 
