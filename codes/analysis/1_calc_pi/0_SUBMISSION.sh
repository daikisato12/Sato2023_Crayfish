qsub -l epyc -l s_vmem=40G -l mem_req=40G j0_mpileup_eachpool_pi.sh
qsub -l medium -l s_vmem=20G -l mem_req=20G j1_calc_depth.sh
qsub -l medium -l s_vmem=10G -l mem_req=10G -hold_jid 14929679 j2_calc_pi_window1kbp.sh
