#!/bin/bash

statsf=SimWT_W1_V2_paramST_CSVs_v45.csv

str="file_name,nb_dups,nb_runs,dups_per_run,threshold,threshold_clear,threshold_hybrid,maximum_sum_row,nb_subfunc_over_nb_dups,nb_neofunc_over_nb_dups,nb_cons_over_nb_dups,nb_pseudo_over_nb_dups,nb_spec_over_nb_dups,nb_dblneo__dups,clear_fates %,hybrid_fates %,superhybrid_fates %,subfunc_portion,neofunc_portion,cons_portion,pseudo_portion,spec_portion,dblneo_portion"
echo $str > $statsf


for d in SimWT_W1_V2_paramST_CSVs_v45/ ; do

#for d in *WT_W0.5_V1_paramST_CSVs_v45/ ; do
    cd $d
    echo "$d"
    #rm concat*
    #rm trian*
    nb_runs=$(ls -d *dups_fates* | wc -l); 
    bash ../script_sim_2.sh concat* ../$statsf $d $nb_runs > allstats_t1_c95_h2.txt
    cd ..
done