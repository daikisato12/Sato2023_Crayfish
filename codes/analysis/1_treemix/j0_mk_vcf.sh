#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 10

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads
cd $WORK_DIR

LIST=(`cat /home/daikisato177/Projects/Crayfish/Data/Pool/samplelist.txt`)
LIST2=""
for line in "${LIST[@]}"; do
        LIST2="$LIST2 $line.sort.filtered.bam";
done

samtools mpileup -C50 -q20 -Q20 -uf ${REF} -t AD,INFO/AD,DP,DV,DPR,INFO/DPR,DP,DP4,SP -v ${LIST2} | bcftools call -vmO z -o ./mpileup/pooled_7pops.vcf.gz

echo ending at
date
