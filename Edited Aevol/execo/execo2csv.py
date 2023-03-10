#!/usr/bin/env python
import os
from pprint import pprint
import sys
import getopt
import numpy as np
import csv
import linecache

def main(argv):
    result_dir="defaultdir"
    output_csv="default.csv"
    try:
        opts, extraparams = getopt.getopt(sys.argv[1:],"hi:o:",["input-dir=","output-csv="]) 
    except getopt.GetoptError:
        print 'test.py -i <inputdirectory> -o <outputfile.csv>'
        sys.exit(2)

    for o,p in opts:
        if o == '-h':
            print 'execoToCSV.py -i <inputdir> -o <outputcsv>'
            print 'inputdir: a input directory where the EXECO results are stored'
            print 'outputdir: an output csv file where the results will be stored'
            sys.exit(2)
        elif o in ['-i','--input-dir']:
            result_dir = p
        elif o in ['-o','--output-csv']:
            output_csv = p


    if result_dir == "defaultdir":
        print "You must define a input directory where the EXECO results are stored"
        sys.exit()

    if output_csv == "default.csv":
        print "You must define an output csv file where the results will be stored"
        sys.exit()

    print 'Input result directory:',result_dir
    print 'Output CSV file: ',output_csv

    list_comb = os.listdir(result_dir)
   
    with open(output_csv, 'wb') as csvfile:
      csvwriter = csv.writer(csvfile)
      with open('toredo_xp.csv','wb') as csvfnotfinish:
	fnotfinish = csv.writer(csvfnotfinish)
	
	for comb_dir in list_comb:
	  print comb_dir
	  k = comb_dir.replace('/',' ').split('-')
	  i = iter(k)
	  params = dict(zip(i,i))
	  
	  if os.path.isfile(result_dir+'/'+comb_dir+"/logger_csv.log"):
	    with open(result_dir+'/'+comb_dir+"/logger_csv.log") as f:
		for line in f:
		  if line.split(",")[0] == 'SELECTION' and len(line.split(",")) == 3:
		    csvwriter.writerow([params['seed'],params['experiment'],params['fuzzy'],params['compilator'],
			  params['parallel'],params['number_of_generation'],int(line.split(",")[1]),int(line.split(",")[2])])
	    
	    
if __name__ == "__main__":
   main(sys.argv[1:])