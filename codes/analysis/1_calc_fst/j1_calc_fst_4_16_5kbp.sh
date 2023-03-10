#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=${data}/${mapped_reads_pool}
cd $WORK_DIR
mkdir -p ./fst

OUTFILE=pooled_7pops_indelfiltered_mincov4_maxcov16_1kbp
perl ${build}/popoolation2_1201/fst-sliding.pl --input ./mpileup/pooled_7pops.indelfiltered.sync --output ./fst/$OUTFILE.fst --suppress-noninformative --min-coverage 4 --max-coverage 16 --min-covered-fraction 0.5 --window-size 5000 --step-size 5000 --pool-size 7
## Fst results are stored in /data/analyzed_data/Fst/pooled_7pops_indelfiltered_mincov4_maxcov16_5kbp.fst ##

echo ending at
date
