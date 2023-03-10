#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads
cd $WORK_DIR

OUTFILE=pooled_7pops_indelfiltered_mincov4_maxcov16_1kbp
mkdir -p ./fst/
perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/fst-sliding.pl --input ./mpileup/pooled_7pops.indelfiltered.sync --output ./fst/$OUTFILE.fst --suppress-noninformative --min-coverage 4 --max-coverage 16 --min-covered-fraction 0.5 --window-size 1000 --step-size 1000 --pool-size 7

WORK_DIR2=/home/daikisato177/Projects/Crayfish/Analysis/Perl_script
cd $WORK_DIR2
mkdir -p ../Result/220323
perl mk_LG_scaf_2022.pl ../../Data_2022/Mapped_reads/fst/$OUTFILE.fst 0 ../../Data_2022/Mapped_reads/fst/$OUTFILE\_annotatedCHR.txt #fst=0
perl mkformat_PBS_2022.pl fst/$OUTFILE\_annotatedCHR.txt Aomori kamakura2 tonden2 220323/$OUTFILE\_annotatedCHR_Aomori_kamakura2_tonden2_PBS
perl mkformat_PBS_2022.pl fst/$OUTFILE\_annotatedCHR.txt Aomori okinawa2 tonden2 220323/$OUTFILE\_annotatedCHR_Aomori_okinawa2_tonden2_PBS
perl mkformat_PBS_2022.pl fst/$OUTFILE\_annotatedCHR.txt kamakura2 okinawa2 tonden2 220323/$OUTFILE\_annotatedCHR_kamakura2_okinawa2_tonden2_PBS

echo ending at
date
