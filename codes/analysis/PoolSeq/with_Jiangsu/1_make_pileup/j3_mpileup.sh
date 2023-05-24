#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

BAMDIR=${data}/poolseq/mapped_reads
cd $BAMDIR

FILES=`ls -1 $BAMDIR/*.sort.filtered.bam`
FILELIST=$(eval echo '$FILES')
POOLNAME="pooled_8pops"
REF=${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

# make pileup
sambamba mpileup -t 10 -o ./mpileup/${POOLNAME}.pileup ${FILELIST} --samtools -C50 -q20 -Q20 -f ${REF}
# make vcf based on pileup file
sambamba mpileup ${FILELIST} -o ./mpileup/${POOLNAME}.vcf.gz --samtools -C50 -q20 -Q20 -f ${REF} --bcftools call -vmO z

perl /codes/software/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ./mpileup/${POOLNAME}.pileup --indel-window 5 --min-count 6  --output ./mpileup/${POOLNAME}.indels.gtf
perl /codes/software/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ./mpileup/${POOLNAME}.indels.gtf --input ./mpileup/${POOLNAME}.pileup --output ./mpileup/${POOLNAME}.indelfiltered.pileup
#rm ./mpileup/pooled_7pops.pileup
java -ea -Xmx7g -jar /codes/software/popoolation2_1201/mpileup2sync.jar --input ./mpileup/${POOLNAME}.indelfiltered.pileup --output ./mpileup/${POOLNAME}.indelfiltered.sync --fastq-type sanger --min-qual 20 --threads 10
#rm ./mpileup/pooled_7pops.indelfiltered.pileup
perl /codes/software/popoolation2_1201/subsample-synchronized.pl --method withoutreplace --max-coverage 50 --target-coverage 10 --input ./mpileup/${POOLNAME}.indelfiltered.sync --output ./mpileup/${POOLNAME}.indelfiltered.subsampled.sync
#rm ./mpileup/pooled_7pops.indelfiltered.sync

echo ending at
date
echo ending at
date
