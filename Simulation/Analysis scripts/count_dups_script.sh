#!/bin/bash

#to once
#f=concat_SimWT_W0.5_V2_param3_CSVs.csv
f="$1"
fstats="$2"
d="$3"
nb_runs="$4"
fstats2="$5"

tot=`wc -l $f | awk '{ print $1 }'`
((tot=tot-1))

echo "hello"

nb_subfunc=0
nb_neofunc=0
nb_cons=0
nb_pseudo=0
nb_spec=0 
nb_dblneo=0

sum_subfunc=0
sum_neofunc=0
sum_cons=0
sum_pseudo=0
sum_spec=0
sum_dblneo=0 

nb_lines=0
nb_empty_lines=0

clear=0
hybrid=0
superhybrid=0

threshold=0.98
threshold_clear=0.95
threshold_hybrid=0.2


echo "threshold  " $threshold
echo "threshold_clear " $threshold_clear
echo "threshold_hybrid " $threshold_hybrid
echo


re='^[0-9]+([.][0-9]+)?$'
max_row=1


write_csv(){
    echo \"$1\",\"$2\",\"$3\" >> $fstats
}

write_csv2(){
    echo \"$1\",\"$2\",\"$3\",\"$4\" >> $fstats2
}


arr=()
run_number=0
IFS=","
tot_nb_dups=0

while read -r original_node g_m g_h g_w first_descendant_id a_m a_h a_w second_descendant_id b_m b_h b_w P_subfunctionlization P_conservation P_newfunctionlizatoin P_pseudogenization P_specialization P_dblneo orig_file
do

if [ -z "$original_node" ]
then
      ((nb_empty_lines=nb_empty_lines+1))
      (IFS=$'\n'; sort <<< "${arr[*]}") | uniq -c 
      nb_dups_run=$((IFS=$'\n'; sort <<< "${arr[*]}") | uniq -c | wc -l);
      echo "nb_dups_per_run: " $nb_dups_run
      echo
      ((run_number=run_number+1))
      ((tot_nb_dups=tot_nb_dups+nb_dups_run))
      write_csv $d $run_number $nb_dups_run 
      arr=()
else


((nb_lines=nb_lines+1))

if ((nb_lines > 1 )); then


arr+=($original_node)


#echo $nb_lines
fi

fi                                   
done < $f



((tot=nb_lines-1))
echo
echo $nb_lines
echo

avg_nb_dups=`echo $tot_nb_dups/${nb_runs} |bc -l`

#str="$f, $nb_lines, $nb_empty_lines, $dups_per_run, $threshold, $threshold_clear, $threshold_hybrid, $max_row,$subfunc,$neofunc,$cons,$pseudo,$spec,$dblneo,$clear_fates,$hybrid_fates,$superhybrid_fates,$subfunc_portion,$neofunc_portion,$cons_portion,$pseudo_portion,$spec_portion,$dblneo_portion"
write_csv2 $d $avg_nb_dups $nb_runs $tot_nb_dups