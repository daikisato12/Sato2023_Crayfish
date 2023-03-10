#!/bin/sh
#$ -S /bin/bash
#$ -t 1-7
#$ -cwd
#qsub j1_calc_pi_211205_window500bp.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=/home/daikisato177/Projects/Crayfish/Data/Pool/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup
cd $WORK_DIR

#bedtools genomecov -ibam ../${SAMPLE}.sort.filtered.bam -g $REF > ../${SAMPLE}.sort.filtered.depth
#grep -n 'genome' ../${SAMPLE}.sort.filtered.depth | perl -lane '$hash{$F[4]} = $F[1] if($F[1] != 0); END{foreach my $key (sort {$b <=> $a} keys %hash){print $hash{$key}; last}}' >../${SAMPLE}.sort.filtered.depth.mode

#TARGET_COV=`cat ../${SAMPLE}.sort.filtered.depth.mode | awk '{print $1}'`
#MAX_COV=$((2 * TARGET_COV))

#echo $SAMPLE
#echo "target coverage: $TARGET_COV"
#echo "max coverage: $MAX_COV"

if [ $SAMPLE = "triunfo2" ]; then
  POOLSIZE=7
elif [ $SAMPLE = "NewOrleans" ]; then
  POOLSIZE=6
else
  POOLSIZE=8
fi

#perl /home/daikisato177/Projects/Crayfish/Build/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --method withoutreplace --max-coverage $MAX_COV --target-coverage $TARGET_COV --input ${SAMPLE}.indelfiltered.pileup --output ${SAMPLE}.indelfiltered.subsampled.pileup
perl /home/daikisato177/Projects/Crayfish/Build/popoolation_1.2.2/Variance-sliding.pl --input ${SAMPLE}.indelfiltered.subsampled.pileup --output ${SAMPLE}.indelfiltered.subsampled.pi.1kb --fastq-type sanger --measure pi --window-size 1000 --step-size 1000 --min-count 2 --min-qual 20 --pool-size $POOLSIZE
perl /home/daikisato177/Projects/Crayfish/Build/popoolation_1.2.2/Variance-sliding.pl --input ${SAMPLE}.indelfiltered.subsampled.pileup --output ${SAMPLE}.indelfiltered.subsampled.theta.1kb --fastq-type sanger --measure theta --window-size 1000 --step-size 1000 --min-count 2 --min-qual 20 --pool-size $POOLSIZE
perl /home/daikisato177/Projects/Crayfish/Build/popoolation_1.2.2/Variance-sliding.pl --input ${SAMPLE}.indelfiltered.subsampled.pileup --output ${SAMPLE}.indelfiltered.subsampled.d.1kb --fastq-type sanger --measure D --window-size 1000 --step-size 1000 --min-count 2 --min-qual 20 --pool-size $POOLSIZE

echo ending at
date
