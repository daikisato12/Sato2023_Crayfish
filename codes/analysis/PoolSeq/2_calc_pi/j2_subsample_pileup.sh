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
 
# already done when making bam files
# bedtools genomecov -ibam ../${SAMPLE}.sort.filtered.bam -g $REF > ../${SAMPLE}.sort.filtered.depth
# grep -n 'LG' ../${SAMPLE}.sort.filtered.depth | perl -lane 'if(defined $hash{$F[1]}){$hash{$F[1]} = $hash{$F[1]}+$F[4] if($F[1] != 0)}else{$hash{$F[1]} = $F[4] if($F[1] != 0)}; END{foreach my $key (sort {$hash{$b} <=> $hash{$a}} keys %hash){print $key; last}}' > ../${SAMPLE}.sort.filtered.depth.mode

TARGET_COV=`cat ../${SAMPLE}.sort.filtered.depth.mode | awk '{print $1}'`
MAX_COV=$((2 * TARGET_COV))

echo $SAMPLE
echo "target coverage: $TARGET_COV"
echo "max coverage: $MAX_COV"

if [ $SAMPLE = "triunfo2" ]; then
  POOLSIZE=7
elif [ $SAMPLE = "NewOrleans" ]; then
  POOLSIZE=6
else
  POOLSIZE=8
fi

perl /codes/software/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --method withoutreplace --max-coverage $MAX_COV --target-coverage $TARGET_COV --input ${SAMPLE}.indelfiltered.pileup --output ${SAMPLE}.indelfiltered.subsampled.pileup
echo "Done making subsampled pileup."

echo ending at
date
