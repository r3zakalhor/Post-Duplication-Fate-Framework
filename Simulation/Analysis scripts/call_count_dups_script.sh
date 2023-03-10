#!/bin/bash

statsf=nb_dups_run.csv
statsf2=avg_nb_dups.csv

str="file_name, run number, nb_dups in run"
str2="file_name, avg_nb_dups ,nb_runs, total_nb_dups"

echo $str > $statsf
echo $str2 > $statsf2

for d in *W0.01*v45/ ; do
    cd $d
    echo "$d"
    nb_runs=$(ls -d *dups_fates* | wc -l); 
    bash ../count_dups.sh concat* ../$statsf $d $nb_runs ../$statsf2> nb_dups.txt
    cd ..
done

for d in *W0.1*v45/ ; do
    cd $d
    echo "$d"
    nb_runs=$(ls -d *dups_fates* | wc -l); 
    bash ../count_dups.sh concat* ../$statsf $d $nb_runs ../$statsf2> nb_dups.txt
    cd ..
done

for d in *W0.5*v45/ ; do
    cd $d
    echo "$d"
    nb_runs=$(ls -d *dups_fates* | wc -l); 
    bash ../count_dups.sh concat* ../$statsf $d $nb_runs ../$statsf2> nb_dups.txt
    cd ..
done

for d in *W1*v45/ ; do
    cd $d
    echo "$d"
    nb_runs=$(ls -d *dups_fates* | wc -l); 
    bash ../count_dups.sh concat* ../$statsf $d $nb_runs ../$statsf2> nb_dups.txt
    cd ..
done
