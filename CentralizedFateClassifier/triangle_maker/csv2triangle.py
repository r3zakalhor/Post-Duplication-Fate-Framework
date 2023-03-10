import os
from os import listdir
from os.path import isfile, join
import sys, getopt
import shutil

def main(argv):

    central_exe = "/home/manuel/git/edited-aevol/CentralizedFateClassifier/CentralizedFateClassifier"
    outfilename = "out.csv"
    filter_diff = False
    include_repeats = True

    try:
        opts, args = getopt.getopt(argv,"ivodr",["infile=", "version=", "outfile=", "diff", "norepeats"])
    except getopt.GetoptError:
        print ('argument error, check code for syntax')
        sys.exit(2)
    
    version = "4"
    for opt, arg in opts:
        if opt in ("--infile"):
            infile = arg
        elif opt == "--version":
            version = arg
        elif opt == "--outfile":
            outfilename = arg
        elif opt == "--diff":
            filter_diff = True
        elif opt == "--norepeats":
            include_repeats = False
            
    column_names = []
    rows = []
    firstrow_str = ""
    is_first_line = True
    
    with open(infile) as file:
        for line in file:
            line = line.replace("\n", "").replace("\r", "")
            
            pz = line.split(",")
            if is_first_line:
                is_first_line = False
                for p in pz:
                    column_names.append(p.strip())
                firstrow_str = line
            elif line != "":
                
                if len(pz) > 0:
                    row = dict()
                    for (i, p) in enumerate(pz):
                        row[ column_names[i] ] = p
                    rows.append(row)
    
    
    outfile = open(outfilename, 'w')
    outfile.write(firstrow_str + "\n")
    
    seencmd = set()
    
    for row in rows:
        #we call CentralizedFateClassifier -m triangleviewer gm gh gw am ah aw bm bh bw
        #print(row)
        cmd = central_exe + " -m triangleviewer " + row["g.m"] + " " + row["g.h"] + " " + row["g.w"] + " " + row["a.m"] + " "
        cmd += row["a.h"] + " " + row["a.w"] + " " + row["b.m"] + " " + row["b.h"] + " " + row["b.w"]
        cmd += " " + version
        cmd += " >> tmp.txt"
        
        if not include_repeats and cmd in seencmd:
            continue
            
        seencmd.add(cmd)
        #print(cmd)
        os.system(cmd)
        
        params = dict()
        with open("tmp.txt") as file:
            for line in file:
                line = line.replace("\n", "").replace("\r", "")
                pz = line.split("=")
                if len(pz) == 2:
                   params[pz[0]] = pz[1]
        
        diffsubfunc = float(row['P_subfunctionlization']) - float(params['P_subfunc'])
        diffcons = float(row['P_conservation']) - float(params['P_cons'])
        
        if not filter_diff or abs(diffsubfunc) >= 0.1 or abs(diffcons) >= 0.1:
            #new row must have the form 
            #original node, g.m, g.h, g.w, first descendant id, a.m, a.h, a.w, second descendant id, b.m, b.h, b.w, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization 
            
            
            outfile.write(f"{row['original node']},{row['g.m']},{row['g.h']},{row['g.w']},{row['first descendant id']},{row['a.m']},{row['a.h']},{row['a.w']},{row['second descendant id']},")
            outfile.write(f"{row['b.m']},{row['b.h']},{row['b.w']},")
            outfile.write(f"{params['P_subfunc']},{params['P_cons']},{params['P_neo']},{params['P_pseudo']},{params['P_spec']}")
            print(f"delta subfunc={float(row['P_subfunctionlization']) - float(params['P_subfunc'])} cons={float(row['P_conservation']) - float(params['P_cons'])}")
        
            outfile.write("\n")
        
    outfile.close()
    
    

if __name__ == "__main__":
   main(sys.argv[1:])