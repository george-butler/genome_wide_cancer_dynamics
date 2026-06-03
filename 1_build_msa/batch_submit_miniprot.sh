#!/usr/bin/env bash

for i in `cat data/genomes.txt`; 
do 
d=$( dirname ${i} ); 
sbatch scripts/submit_miniprot.sh ${i} data/test_prot.fasta ${d}/miniprot_all ; 
done
