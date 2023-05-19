qsub -pe def_slot 10 j1_fasterq_dump.sh
qsub -pe def_slot 10 j2_sample_fastq.sh
cd ${data}/poolseq/fastq/subsample
cat *_subsample20M_1.fastq.gz > ../rawdata/jiangsu_pseudopool_1.fastq.gz
cat *_subsample20M_2.fastq.gz > ../rawdata/jiangsu_pseudopool_2.fastq.gz