#!/bin/sh
#$ -S /bin/bash
#$ -t 1-7
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/poolseq/fastq/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

REF=${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORKDIR=${data}/poolseq/mapped_reads/mpileup
cd $WORKDIR

samtools mpileup -C 50 -q 20 -Q 20 -f ${REF} ../${SAMPLE}.sort.filtered.bam > ${SAMPLE}.pileup
perl /codes/software/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ${SAMPLE}.pileup --indel-window 5 --min-count 6  --output ${SAMPLE}.indels.gtf
perl /codes/software/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ${SAMPLE}.indels.gtf --input ${SAMPLE}.pileup --output ${SAMPLE}.indelfiltered.pileup
# rm ${SAMPLE}.pileup

echo ending at
date
