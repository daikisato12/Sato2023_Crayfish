#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORK_DIR=/home/daikisato177/Data/crayfish/Ref_Xu2021/v2/blastp
QUERY_PEPFILE=/home/daikisato177/Data/crayfish/Ref_Xu2021/v2/GCF_020424385.1_ASM2042438v2_translated_cds_longest.aa
DROS_PEPFILE=/home/daikisato177/Data/crayfish/Ref_Drosophila/Drosophila_melanogaster.BDGP6.32.pep.all.fa
OUTFILE=GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros

cd $WORK_DIR
#makeblastdb -in $NEW_PEPFILE -dbtype prot -hash_index
blastp -evalue 1e-5 -outfmt 7 -query $QUERY_PEPFILE -db $DROS_PEPFILE -out $OUTFILE.blastp
#cat $OUTFILE.blastp | awk '/hits found/{getline;print}' | grep -v "#" > $OUTFILE.best_hits.txt

echo ending at
date
