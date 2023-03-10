#!/bin/bash
input="temp.in"
NumOfWIldTypes=6
NumOfReplications=20
NumOfGenerations=1100000
mkdir SimWT_W0.1_V1_param1_CSVs
# $IFS removed to allow the trimming #
#####################################
for (( i=5; i<=$NumOfWIldTypes-1; i++ ))
do
	for (( k=1; k<=$NumOfReplications; k++ ))
	do
		counter=0
		echo -n "" > param.in
		while read -r line
		do
		  let counter++
		  #echo $counter
		  if [ $counter = 7 ];
		  then
			j=0
			for word in $line; do
				if [ $j == 0 ];
				then
					word1="$word"
					#echo "$word1"
				fi
				if [ $j == 1 ];
				then
					RANDOM=$(date +%s)/k;
					#echo $RANDOM;
					word2="$RANDOM"
				fi
				let j++
			done
			newline="${word1} ${word2}"
			echo "$newline" >> param.in
			#echo "yes!"

		  #elif [ $counter = 39 ] && [ $i = 5 ];
		  #then
			#j=0
			#for word in $line; do
			#	if [ $j == 0 ];
			#	then
			#		word3="$word"
					#echo "$word1"
			#	fi
			#	if [ $j == 1 ];
			#	then
			#		word4=1.0
			#	fi
			#	let j++
			#done
			#newline="${word3} ${word4}"
			#echo "$newline" >> param.in
			#echo "yes!"

		  else
		  echo "$line" >> param.in
		   #echo "$line"
		  fi
		done < "$input"

		# run aevol here
		../../build/bin/aevol_create -C WT_W0.1_V1.txt
		../../build/bin/aevol_run -n $NumOfGenerations
		../../build/bin/aevol_post_lineage
		../../build/bin/aevol_post_protein_map lin*
		./CentralizedFateClassifier -m default proteins_list_after_events.csv -e $NumOfGenerations-100000
		file="WT_W0.1_V1_${k}"
		echo "$file"
		mkdir "$file"

		cp proteins_list_after_events.csv SimWT_W0.1_V1_param1_CSVs/"${file}proteins_list_after_events.csv"
		cp dups_fates_avg_probablities.csv SimWT_W0.1_V1_param1_CSVs/"${file}dups_fates_avg_probablities.csv"
		cp dups_fates_probablities.csv SimWT_W0.1_V1_param1_CSVs/"${file}dups_fates_probablities.csv"
		cp newick.txt SimWT_W0.1_V1_param1_CSVs/"${file}newick.txt"

		mv proteins_list_after_events.csv "$file"/
		mv lin* "$file"/
		mv dups_fates_avg_probablities.csv "$file"/
		mv dups_fates_probablities.csv "$file"/
		mv newick.txt "$file"/
		mv param.in "$file"/
	done
done
