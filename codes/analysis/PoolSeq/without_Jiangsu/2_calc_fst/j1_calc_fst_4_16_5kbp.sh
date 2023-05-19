#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORKDIR=${data}/poolseq/mapped_reads/mpileup
cd $WORKDIR

POOLNAME="pooled_7pops"
OUTFILE=${POOLNAME}.indelfiltered.mincov4_maxcov16_5kbp

perl /codes/software/popoolation2_1201/fst-sliding.pl --input ${POOLNAME}.indelfiltered.sync --output ./fst/$OUTFILE.fst --suppress-noninformative --min-coverage 4 --max-coverage 16 --min-covered-fraction 0.5 --window-size 5000 --step-size 5000 --pool-size 7
echo "Done calculating fst."

echo ending at
date
