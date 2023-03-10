#!/bin/bash

#set -x

wtfile=$1
paramfile=$2
aevolpath=$3

workdir="work_${wtfile}_${paramfile}" 

rm -rf "${workdir}"



mkdir -p "${workdir}"

cp "${wtfile}" "${workdir}/${wtfile}"
cp "${paramfile}" "${workdir}/param.in"

cd "${workdir}"

#print the param file in the workdir, but with backup step = 100, save it in a tmp file, then overwrite param file with that tmp
#note that with less than 100 generations, aevol_post_lineage does not want to work for some reason
awk '{if (/BACKUP_STEP*/){print "BACKUP_STEP             100"}else{print $0}}' param.in > tmp
mv tmp param.in 


#run aevol for 100 generations
${aevolpath}/aevol_create -C "${wtfile}"
${aevolpath}/aevol_run -n 100

${aevolpath}/aevol_post_lineage


#do the fitness analysis
${aevolpath}/aevol_post_protein_map lineage-b000000000-e000000100-i3-r-1.ae -f -g 1



echo "Sanity check: checking if WT_analyzed.txt = ${wtfile}"
diff --ignore-space-change --strip-trailing-cr --ignore-blank-lines ${wtfile} WT_analyzed.txt && echo "OK, WTs are identical" || echo "ERROR, WTs differ"



echo "DONE"

cd ..



