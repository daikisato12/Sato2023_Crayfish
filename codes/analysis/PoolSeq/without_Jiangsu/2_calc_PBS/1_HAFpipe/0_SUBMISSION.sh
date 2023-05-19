qsub -l s_vmem=5G -l mem_req=5G j1_prepare_dataset.sh
qsub -l s_vmem=5G -l mem_req=5G j2_HAFpipe_t1-2.sh
qsub -l s_vmem=5G -l mem_req=5G j3_HAFpipe_t3-4.sh
qsub -l s_vmem=20G -l mem_req=20G j4_merge_res.sh