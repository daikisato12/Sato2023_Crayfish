#!/bin/sh
#$ -S /bin/bash
#$ -t 1-8
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=/home/daikisato177/Projects/Crayfish/Data/One_ind/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa
FASTQ1=${FASTQDIR}/${SAMPLE}_1.fq.gz
FASTQ2=${FASTQDIR}/${SAMPLE}_2.fq.gz

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads_1ind
cd $WORK_DIR
bwa-mem2 mem -t 8 ${REF} ${FASTQ1} ${FASTQ2} > ${SAMPLE}.sam
samtools sort -@ 8 -O bam -o ${SAMPLE}.sort.bam ${SAMPLE}.sam
rm ${SAMPLE}.sam

elprep filter ${SAMPLE}.sort.bam ${SAMPLE}.sort.filtered.bam \
	--output-type bam \
	--mark-duplicates --remove-duplicates --filter-mapping-quality 20 --clean-sam \
	--nr-of-threads 8 --sorting-order coordinate
samtools index -@ 8 ${SAMPLE}.sort.filtered.bam
rm ${SAMPLE}.sort.bam

echo ending at
date
