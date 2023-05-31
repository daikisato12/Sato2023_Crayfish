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

if [ $SAMPLE = "triunfo2" ]; then
  POOLSIZE=7
elif [ $SAMPLE = "NewOrleans" ]; then
  POOLSIZE=6
else
  POOLSIZE=8
fi

WINDOWSIZE="5kb"

perl /codes/software/popoolation_1.2.2/Variance-sliding.pl --input ${SAMPLE}.indelfiltered.subsampled.pileup --output ./pi/${SAMPLE}.indelfiltered.subsampled.pi.${WINDOWSIZE} --fastq-type sanger --measure pi --window-size 5000 --step-size 5000 --min-count 2 --min-qual 20 --pool-size $POOLSIZE
echo "Done calculating pi."

echo ending at
date
