qsub -l s_vmem=40G -l mem_req=40G j1_mpileup_eachpool_pi.sh
qsub -l s_vmem=20G -l mem_req=20G j2_subsample_pileup.sh
qsub -l s_vmem=10G -l mem_req=10G j3_calc_pi_window5kbp.sh
