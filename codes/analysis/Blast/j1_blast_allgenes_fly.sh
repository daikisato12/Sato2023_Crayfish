#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORKDIR=${data}/blastp
QUERY_PEPFILE=${data}/ref/GCF_020424385.1_ASM2042438v2_translated_cds_longest.aa
DROS_PEPFILE=${data}/ref/Drosophila_melanogaster.BDGP6.32.pep.all.fa
OUTFILE="GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros"

cd $WORKDIR
makeblastdb -in $DROS_PEPFILE -dbtype prot -hash_index
blastp -evalue 1e-5 -outfmt 7 -query $QUERY_PEPFILE -db $DROS_PEPFILE -out $OUTFILE.blastp
#cat $OUTFILE.blastp | awk '/hits found/{getline;print}' | grep -v "#" > $OUTFILE.best_hits.txt

perl /codes/analysis/Blast/1_mklist_blastp_res.pl

echo ending at
date
