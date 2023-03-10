#!/bin/bash

#to once
f=concat_testset.csv


tot=`wc -l $f | awk '{ print $1 }'`
((tot=tot-1))

echo $t

nb_subfunc=0
nb_neofunc=0
nb_cons=0
nb_pseudo=0
nb_spec=0 
nb_lines=0



threshold=$1

IFS=","

while read -r original_node g_m g_h g_w first_descendant_id a_m a_h a_w second_descendant_id b_m b_h b_w P_subfunctionlization P_conservation P_newfunctionlizatoin P_pseudogenization P_specialization orig_file
do

((nb_lines=nb_lines+1))

if ((nb_lines > 1 )); then



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

echo $nb_lines
fi
                                   
done < $f

echo "threshold  " $threshold
echo
echo "subfunc"

echo $nb_subfunc/${tot} |bc -l

echo "neofunc"

echo $nb_neofunc/${tot} |bc -l

echo "cons"
echo $nb_cons/${tot} |bc -l

echo "pseudo"
echo $nb_pseudo/${tot} |bc -l

echo "spec"
echo $nb_spec/${tot} |bc -l




