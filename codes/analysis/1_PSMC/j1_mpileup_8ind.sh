#!/bin/sh
#$ -S /bin/bash
#$ -t 1-8
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=/home/daikisato177/Projects/Crayfish/Data/One_ind
cd $FASTQDIR

FILENAME=`ls -1 *_R1.fastq.gz | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 13- | rev | uniq)
REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads_1ind
cd $WORK_DIR
samtools mpileup -C 50 -q 20 -Q 20 -d 200 -u -f $REF ${SAMPLE}.sort.filtered.bam |\
bcftools call -m -v - | bcftools filter -g 10 -G 10 -O v |\
bcftools filter -s FAIL -e '%QUAL<20 || INFO/DP<=3 || INFO/DP>=30 || INFO/MQ<=30 || INFO/MQB<1e-20 || INFO/RPB<0.0001 || INFO/BQB<0.0001 || (INFO/DP4[1]+INFO/DP4[2]<=2) || (INFO/DP4[3]+INFO/DP4[4]<=2)' -O v -> ${SAMPLE}.sort.filtered.vcf

echo ending at
date
