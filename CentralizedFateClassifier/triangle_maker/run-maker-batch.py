import os
from os import listdir
from os.path import isfile, join
import sys, getopt
import shutil


def cut(s, length = 5):
    if s[0] == "-":
        return s[0:length + 1]
    else:
        return s[0:length]


def recalc_csv(indir, outdir, version = "5"):
    
    files = [f for f in listdir(indir) if isfile(join(indir, f))]
    for f in files:
        if "dups_fates_probablities.csv" in f:
            cmd = "python csv2csv.py --infile=" + join(indir, f) + " --version=" + version + " --outfile=" + join(outdir, f)
            print(cmd)
            os.system(cmd)
    

def main(argv):

    #TODO: hardcoded indir
    indir = ""  #"C:\\Users\\Manue\\Desktop\\tmp\\csvs_1_param"
    outdir = ""

    dov5 = False 
    doTriangles = False
    doPng = False
    noTriangleOverwrite = False

    outfilename_cat = "concat.csv"
    outfilename_triangle = "triangles.csv"
    
    
    
    try:
        opts, args = getopt.getopt(argv,"hitopvt",["help", "indir=", "dotriangles", "outdir=", "dopng", "dov5", "notriangleoverwrite"])
    except getopt.GetoptError:
        print ('argument error, check code for syntax')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ['-h', '--help']:
            print('python run-maker.py -indir=[indir]\nWill create concat.csv (concatenation of a dup csvs in indir) and triangles.csv (to use with trianglemaker.py)')
            
            sys.exit()
        elif opt == "--indir":
            indir = arg
        elif opt == "--outdir":
            outdir = arg
        elif opt == "--dotriangles":
            doTriangles = True
        elif opt == "--dov5":
            dov5 = True
        elif opt == "--notriangleoverwrite":
            noTriangleOverwrite = True
        elif opt == "--dopng":
            doPng = True
    
    
    if outdir != "" and not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outfilename_cat = join(outdir, "concat_" + indir + ".csv")
    outfilename_triangle = join(outdir, "triangles_" + indir + ".csv")

    outfile_cat = open(outfilename_cat, 'w')
    outfile_triangle = open(outfilename_triangle, 'w')


    outfile_cat.write("original node, g.m, g.h, g.w, first descendant id, a.m, a.h, a.w, second descendant id, b.m, b.h, b.w, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization, orig_file\n")
    outfile_triangle.write("subfunc, neofunc, cons, pseudo, spec, mg, hg, wg, ma, ha, wa, mb, hb, wb, fitness_g, fitness_a, fitness_b\n")

    csvs_dir = indir
    if dov5:
        print("Recalculating probabilities from " + indir + " to " + outdir)
        recalc_csv(indir, outdir, "5")
        csvs_dir = outdir


    files = [f for f in listdir(csvs_dir) if isfile(join(csvs_dir, f))]
    for f in files:
        if "dups_fates_probablities.csv" in f:
            column_names = []
            rows = []
            is_first_line = True
            with open(join(csvs_dir, f)) as file:
                for line in file:
                    line = line.replace("\n", "").replace("\r", "")
                    
                    
                    
                    pz = line.split(",")
                    if is_first_line:
                        is_first_line = False
                        for p in pz:
                            column_names.append(p.replace(" ", ""))
                    elif line != "":
                        outfile_cat.write(line + ", " + f + "\n")
                        if len(pz) > 0:
                            row = dict()
                            for (i, p) in enumerate(pz):
                                row[ column_names[i] ] = p
                            rows.append(row)
                    else:
                        rows.append( dict() )
            
            if len(rows) > 0:
                outfile_cat.write("\n")
            
            for (r, row) in enumerate(rows):
                if len(row) == 0:
                    outfile_triangle.write("\n")
                else:
                    outputlen = 6
                    outfile_triangle.write(cut(row['P_subfunctionlization']) + "," + cut(row['P_newfunctionlizatoin']) + "," + cut(row['P_conservation']) + "," + cut(row['P_pseudogenization']) + ",")
                    outfile_triangle.write(cut(row['P_specialization']) + ",")
                    outfile_triangle.write(cut(row['g.m']) + "," + cut(row['g.h']) + "," + cut(row['g.w']) + ",")
                    outfile_triangle.write(cut(row['a.m']) + "," + cut(row['a.h']) + "," + cut(row['a.w']) + ",")
                    outfile_triangle.write(cut(row['b.m']) + "," + cut(row['b.h']) + "," + cut(row['b.w']) + ",")
                    outfile_triangle.write("0,0,0\n")   #fake fitnesses

    outfile_cat.close()
    outfile_triangle.close()

    print("Wrote files " + outfilename_cat + " and " + outfilename_triangle)
    
    if not doTriangles:
        print("To generate triangles, run")
        print("python triangle_maker.py --batch=" + outfilename_triangle + " --outdir=[outdir]")
    else:
        outtriangles_dir = join(outdir, "triangle_figs")
        print("Generating triangles in directory " + outtriangles_dir)
        
        if not os.path.exists(outtriangles_dir):
            os.makedirs(outtriangles_dir)
        
        cmd_triangles = "python triangle_maker_multiproc.py --batch=" + outfilename_triangle + " --outdir=" + outtriangles_dir
        if noTriangleOverwrite:
           cmd_triangles += " --nooverwrite"
        
        os.system(cmd_triangles)
    
        #shutil.copyfile(outfilename_cat, join(indir, outfilename_cat))  
        #shutil.copyfile(outfilename_triangle, join(indir, outfilename_triangle))  
    
        if doPng:
            outpng_dir = join(outtriangles_dir, "png")
            
            print("Generating pngs in directory " + outpng_dir)
        
            if not os.path.exists(outpng_dir):
                os.makedirs(outpng_dir)
        
            prevcwd = os.getcwd()
            
            os.chdir(outtriangles_dir)
            
            cmd = "python " + join(prevcwd, "batchconvert.py")
            print(cmd)
            os.system(cmd)
            
            os.chdir(prevcwd)
            
            
            triangle_files = os.listdir(outtriangles_dir)
            for file in triangle_files:
                if file.endswith('.png'):
                    shutil.move(os.path.join(outtriangles_dir,file), os.path.join(outpng_dir,file))
            
            
    

if __name__ == "__main__":
   main(sys.argv[1:])
