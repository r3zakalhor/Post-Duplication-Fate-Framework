#!/bin/bash

#to once
#f=concat_SimWT_W0.5_V2_param3_CSVs.csv
f="$1"
fstats="$2"
d="$3"
nb_runs="$4"

tot=`wc -l $f | awk '{ print $1 }'`
((tot=tot-1))

#echo $tot

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

IFS=","


while read -r original_node g_m g_h g_w first_descendant_id a_m a_h a_w second_descendant_id b_m b_h b_w P_subfunctionlization P_conservation P_newfunctionlizatoin P_pseudogenization P_specialization P_dblneo orig_file
do

if [ -z "$original_node" ]
then
      ((nb_empty_lines=nb_empty_lines+1))
else


((nb_lines=nb_lines+1))

if ((nb_lines > 1 )); then

#s=(($P_subfunctionlization +  $P_newfunctionlizatoin + $P_conservation + $P_pseudogenization + $P_specialization |bc -l))

if ! [[ $P_subfunctionlization =~ $re ]] ; then
   P_subfunctionlization=$(echo $P_subfunctionlization | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)
fi
if ! [[ $P_newfunctionlizatoin =~ $re ]] ; then
   P_newfunctionlizatoin=$(echo $P_newfunctionlizatoin | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)

fi
if ! [[ $P_conservation =~ $re ]] ; then
   P_conservation=$(echo $P_conservation | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)
fi
if ! [[ $P_pseudogenization =~ $re ]] ; then
   P_pseudogenization=$(echo $P_pseudogenization | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)
fi
if ! [[ $P_specialization =~ $re ]] ; then
   P_specialization=$(echo $P_specialization | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)
fi

if ! [[ $P_dblneo =~ $re ]] ; then
   P_dblneo=$(echo $P_dblneo | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=100;/'| bc)
fi



sum=`echo $P_subfunctionlization + $P_newfunctionlizatoin + $P_conservation + $P_pseudogenization + $P_specialization + $P_dblneo| bc -l`

sum_subfunc=`echo $P_subfunctionlization + $sum_subfunc |bc -l `
sum_neofunc=`echo $P_newfunctionlizatoin + $sum_neofunc |bc -l `
sum_cons=`echo $P_conservation + $sum_cons |bc -l `
sum_pseudo=`echo $P_pseudogenization + $sum_pseudo|bc -l `
sum_spec=`echo $P_specialization + $sum_spec|bc -l `
sum_dblneo=`echo $P_dblneo + $sum_dblneo|bc -l `

set -a fates
fates[0]=$P_subfunctionlization
fates[1]=$P_newfunctionlizatoin
fates[2]=$P_conservation
fates[3]=$P_pseudogenization
fates[4]=$P_specialization
fates[5]=$P_dblneo

bigger02=0
for i in 0 1 2 3 4 5;
do
#echo -n ${fates[$i]} " " 
if (( $(echo "${fates[$i]} >= $threshold_clear" |bc -l) )); then
   ((clear=clear+1))
fi
if (( $(echo "${fates[$i]} >= $threshold_hybrid" |bc -l) )); then
   ((bigger02=bigger02+1))
fi
done


if (( $(echo "$bigger02 == 2" |bc -l) )); then
   ((hybrid=hybrid+1))
fi

if (( $(echo "$bigger02 > 2" |bc -l) )); then
   ((superhybrid=superhybrid+1))
fi



#echo

if (( $(echo "$sum > $max_row" |bc -l) )); then
   #echo "Sum of fates is bigger than 1 for line " $nb_lines ": " $sum
   max_row=`echo $sum`
fi   

 
if (( $(echo "$P_subfunctionlization >= $threshold" |bc -l) )); then
   ((nb_subfunc=nb_subfunc+1))
fi   
     
if (( $(echo "$P_newfunctionlizatoin >= $threshold" |bc -l) )); then
   ((nb_neofunc=nb_neofunc+1))
fi   
     
if (( $(echo "$P_conservation >= $threshold" |bc -l) )); then
   ((nb_cons=nb_cons+1))
fi   
     
if (( $(echo "$P_pseudogenization >= $threshold" |bc -l) )); then
   ((nb_pseudo=nb_pseudo+1))
fi   

if (( $(echo "$P_specialization >= $threshold" |bc -l) )); then
   ((nb_spec=nb_spec+1))
fi   

if (( $(echo "$P_dblneo >= $threshold" |bc -l) )); then
   ((nb_dblneo=nb_dblneo+1))
fi   

#echo $nb_lines
fi

fi                                   
done < $f

((tot=nb_lines-1))
echo
echo $nb_lines
echo

echo "maximum sum of a row is " $max_row
echo

echo -n "% subfunc "
subfunc=`echo $nb_subfunc/${tot} |bc -l`
echo $subfunc

echo -n "% neofunc "
neofunc=`echo $nb_neofunc/${tot} |bc -l`
echo $neofunc

echo -n "% cons "
cons=`echo $nb_cons/${tot} |bc -l`
echo $cons

echo -n "% pseudo "
pseudo=`echo $nb_pseudo/${tot} |bc -l`
echo $pseudo

echo -n "% spec "
spec=`echo $nb_spec/${tot} |bc -l`
echo $spec


echo -n "% dblneo "
dblneo=`echo $nb_dblneo/${tot} |bc -l`
echo $dblneo

echo

echo -n "% clear fates "
clear_fates=`echo ${clear}/${tot} |bc -l`
echo $clear_fates

echo -n "% hybrid fates "
hybrid_fates=`echo ${hybrid}/${tot} |bc -l`
echo $hybrid_fates

echo -n "% superhybrid fates "
superhybrid_fates=`echo ${superhybrid}/${tot} |bc -l`
echo $superhybrid_fates


echo

sum_all=`echo $sum_subfunc + $sum_neofunc + $sum_cons + $sum_pseudo + $sum_spec + $sum_dblneo| bc -l`

#echo "sum_all " $sum_all
#echo "sum_subfunc " $sum_subfunc
#echo "sum_cons " $sum_cons
#echo "sum_neofunc " $sum_neofunc
#echo "sum_pseudo " $sum_pseudo
#echo "sum_spec " $sum_spec
echo 

echo -n "sum_subfunc/sum_all " 
subfunc_portion=`echo $sum_subfunc/${tot} |bc -l`
echo $subfunc_portion

echo -n "sum_neofunc/sum_all " 
neofunc_portion=`echo $sum_neofunc/${tot} |bc -l`
echo $neofunc_portion

echo -n "sum_cons/sum_all " 
cons_portion=`echo $sum_cons/${tot} |bc -l`
echo $cons_portion

echo -n "sum_pseudo/sum_all " 
pseudo_portion=`echo $sum_pseudo/${tot} |bc -l`
echo $pseudo_portion

echo -n "sum_spec/sum_all " 
spec_portion=`echo $sum_spec/${tot} |bc -l`
echo $spec_portion


echo -n "sum_dblneo/sum_all " 
dblneo_portion=`echo $sum_dblneo/${tot} |bc -l`
echo $dblneo_portion


echo -n "nb dups/nb runs " 
dups_per_run=`echo $nb_lines/${nb_empty_lines} |bc -l`
echo $dups_per_run


write_csv(){
    echo \"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\",\"$7\",\"$8\",\"$9\",\"${10}\",\"${11}\",\"${12}\",\"${13}\",\"${14}\",\"${15}\",\"${16}\",\"${17}\",\"${18}\",\"${19}\",\"${20}\",\"${21}\",\"${22}\",\"${23}\" >> $fstats
}

str="$f, $nb_lines, $nb_empty_lines, $dups_per_run, $threshold, $threshold_clear, $threshold_hybrid, $max_row,$subfunc,$neofunc,$cons,$pseudo,$spec,$dblneo,$clear_fates,$hybrid_fates,$superhybrid_fates,$subfunc_portion,$neofunc_portion,$cons_portion,$pseudo_portion,$spec_portion,$dblneo_portion"
write_csv $d $nb_lines $nb_runs $dups_per_run $threshold $threshold_clear $threshold_hybrid $max_row $subfunc $neofunc $cons $pseudo $spec $dblneo $clear_fates $hybrid_fates $superhybrid_fates $subfunc_portion $neofunc_portion $cons_portion $pseudo_portion $spec_portion $dblneo_portion
