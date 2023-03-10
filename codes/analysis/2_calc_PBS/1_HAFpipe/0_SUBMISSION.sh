qsub -l s_vmem=5G -l mem_req=5G j1_HAFpipe_t1-2.sh
qsub -l s_vmem=5G -l mem_req=5G -hold_jid 16886699 j2_HAFpipe_t3-4.sh
qsub -l medium -l s_vmem=20G -l mem_req=20G j3_merge_res.sh
