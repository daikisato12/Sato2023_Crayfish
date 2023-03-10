qsub -l medium -l s_vmem=20G -l mem_req=20G -pe def_slot 10 j0_mpileup.sh
qsub -l medium -l s_vmem=20G -l mem_req=20G -hold_jid 15003976 j1_calc_fst_4_16_1kbp.sh