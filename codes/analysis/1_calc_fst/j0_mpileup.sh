#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


REF=${data}/${ref}/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORK_DIR=${data}/${mapped_reads_pool}
cd $WORK_DIR
mkdir -p ./mpileup

LIST=(`cat /codes/analysis/0_mapping/samplelist.txt`)
LIST2=""
for line in "${LIST[@]}"; do
	LIST2="$LIST2 $line.sort.filtered.bam";
done
echo $LIST2

sambamba mpileup -t 10 -o ./mpileup/pooled_7pops.pileup ${LIST2} --samtools -C50 -q20 -Q20 -f ${REF} 
perl ${build}/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ./mpileup/pooled_7pops.pileup --indel-window 5 --min-count 6  --output ./mpileup/pooled_7pops.indels.gtf
perl ${build}/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ./mpileup/pooled_7pops.indels.gtf --input ./mpileup/pooled_7pops.pileup --output ./mpileup/pooled_7pops.indelfiltered.pileup
# rm ./mpileup/pooled_7pops.pileup
${build}/jre1.8.0_241/bin/java -ea -Xmx7g -jar ${build}/popoolation2_1201/mpileup2sync.jar --input ./mpileup/pooled_7pops.indelfiltered.pileup --output ./mpileup/pooled_7pops.indelfiltered.sync --fastq-type sanger --min-qual 20 --threads 10
# rm ./mpileup/pooled_7pops.indelfiltered.pileup

echo ending at
date
