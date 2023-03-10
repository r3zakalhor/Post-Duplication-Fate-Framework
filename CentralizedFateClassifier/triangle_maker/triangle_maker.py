import os
import sys, getopt
import shutil
import os.path


#!/usr/bin/python
# python triangle_maker.py --mg 0.34 --hg 0.78 --wg 0.03 --ma 0.26 --ha 1 --wa 0.025 --mb 0.4 --hb 1 --wb 0.025 --pneo 0 --psubfunc 0 --pcons 0 --pspec 1 --ppseudo 0
mg, hg, wg = 0.333333, 0.779528, 0.0222222
ma, ha, wa = 0.253968, 1, 0.0247312
mb, hb, wb = 0.253968, 1, 0.0247312

pneo, psubfunc, pcons, pspec = 0.1, 0.3, 0.2, 0.295
ppseudo = 0.1
ppseudoboth = 0.05


outdir = ""
overwrite = True


def parse_csv(filename):
    global outdir
    global overwrite 
    column_names = []
    rows = []
    is_first_line = True
    with open(filename) as file:
        for line in file:
            line = line.replace("\n", "").replace("\r", "")
            
            pz = line.split(",")
            if is_first_line:
                is_first_line = False
                for p in pz:
                    column_names.append(p.replace(" ", ""))
            elif line != "":
                if len(pz) > 0:
                    row = dict()
                    for (i, p) in enumerate(pz):
                        row[ column_names[i] ] = p
                    rows.append(row)
    
    for (r, row) in enumerate(rows):
        outfilename = ""
        if outdir != '':
            outfilename += outdir + "/" 
        outfilename += str(r + 2) + ".pdf"
        
        if overwrite or not os.path.isfile(outfilename):
            
            cmd = f"python triangle_maker.py --mg {row['mg']} --hg {row['hg']} --wg {row['wg']} --ma {row['ma']} --ha {row['ha']} --wa {row['wa']} --mb {row['mb']} --hb {row['hb']} --wb {row['wb']} --pneo {row['neofunc']} --psubfunc {row['subfunc']} --pcons {row['cons']} --pspec {row['spec']} --ppseudo {row['pseudo']}"
            cmd += " --outfile=" + outfilename 
        
            print(cmd)
            os.system(cmd)


def main(argv):
    gcolor = 'green'
    acolor = 'blue'
    bcolor = 'red'
    parse_csv_filename = ''
    global outdir
    global overwrite

    template_filename = "tikztriangle_template.tex"
    template_filename_out = "tikz_test_out.tex"

    inputfile = ''
    outputfilename = ''
    
    startfile = False
    
    try:
        opts, args = getopt.getopt(argv,"hs",["help", "startfile", "mg=", "hg=", "wg=", "ma=", "ha=", "wa=", "mb=", "hb=", "wb=", "pneo=", "psubfunc=", "pcons=", "pspec=", "ppseudo=", "gcolor=", "acolor=", "bcolor=", "batch=", "outfile=", "outdir=", "nooverwrite"])
    except getopt.GetoptError:
        print ('ERROR, format is \n  triangle_maker.py --mg value1 --hg value2 --wg value3 --ma value4 --ha value5 --wa value6 --mb value7 --hb value8 --wb value9 --pneo pro1 --psubfunc pro2 --pcons pro3 --pspec pro4 --ppseudo pro5 --gcolor color1 --acolor color2 --bcolor color3')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ['-h', '--help']:
            print('To generate one triangle file:')
            print ('triangle_maker.py --mg value1 --hg value2 --wg value3 --ma value4 --ha value5 --wa value6 --mb value7 --hb value8 --wb value9 --pneo pro1 --psubfunc pro2 --pcons pro3 --pspec pro4 --ppseudo pro5 --gcolor color1 --acolor color2 --bcolor color3')
            print('To generate triangle files in batch from a csv:')
            print('triangle_maker.py --batch=csvfile --outdir=outdir')
            sys.exit()
        elif opt in ("--batch"):
            parse_csv_filename = arg
        elif opt in ["--mg"]:
            mg = float(arg)
        elif opt in ["--startfile"]:
            startfile = True
        elif opt in ["--hg"]:
            hg = float(arg)
        elif opt in ["--wg"]:
            wg = float(arg)
        elif opt in ["--ma"]:
            ma = float(arg)
        elif opt in ["--ha"]:
            ha = float(arg)
        elif opt in ["--wa"]:
            wa = float(arg)
        elif opt in ["--mb"]:
            mb = float(arg)
        elif opt in ["--hb"]:
            hb = float(arg)
        elif opt in ["--wb"]:
            wb = float(arg)
        elif opt in ["--pneo"]:
            pneo = float(arg)
        elif opt in ["--nooverwrite"]:
            overwrite = False
            print("Will not overwrite existing files")
        elif opt in ["--psubfunc"]:
            psubfunc = float(arg)
        elif opt in ["--pcons"]:
            pcons = float(arg)
        elif opt in ["--pspec"]:
            pspec = float(arg)
        elif opt in ["--ppseudo"]:
            ppseudo = float(arg)
        elif opt in ["--gcolor"]:
            gcolor = arg
        elif opt in ["--acolor"]:
            acolor = arg
        elif opt in ["--bcolor"]:
            bcolor = arg
        elif opt in ["--outfile"]:
            outputfilename = arg
        elif opt in ["--outdir"]:
            outdir = arg
   
    if parse_csv_filename != '':
        parse_csv(parse_csv_filename)
        sys.exit()
   
    
    
    def adjust(val, origmin, origmax, destmin, destmax):
        return destmin + (val - origmin) * (destmax - destmin) / (origmax - origmin)
    
    
    
    mgtmp = mg
    matmp = ma
    mbtmp = mb
    
    if matmp + wa < mbtmp - wb and matmp + wa < mgtmp - wg:
        matmp = max(matmp, min(mbtmp - wb, mgtmp - wg) - wa - 0.001)
    if mbtmp + wb < matmp - wa and mbtmp + wb < mgtmp - wg:
        mbtmp = max(mbtmp, min(matmp - wa, mgtmp - wg) - wb - 0.001)
    if mgtmp + wg < mbtmp - wb and mgtmp + wg < matmp - wa:
        mgtmp = max(mgtmp, min(mbtmp - wb, matmp - wa) - wg - 0.001)
        
    if matmp - wa > mbtmp + wb and matmp - wa > mgtmp + wg:
        matmp = min(matmp, max(mbtmp + wb, mgtmp + wg) + wa + 0.001)
    if mbtmp - wb > matmp + wa and mbtmp - wb > mgtmp + wg:
        mbtmp = min(mbtmp, max(matmp + wa, mgtmp + wg) + wb + 0.001)
    if mgtmp - wg > mbtmp + wb and mgtmp - wg > matmp + wa:
        mgtmp = min(mgtmp, max(mbtmp + wb, matmp + wa) + wg + 0.001)
    
    minx = min(mgtmp - wg, matmp - wa, mbtmp - wb)
    maxx = max(mgtmp + wg, matmp + wa, mbtmp + wb)
    
    if minx == maxx:
        maxx = minx + 0.001

    
    hscale = 8 / (maxx - minx)
    mg_adjusted = adjust(mgtmp, minx, maxx, 1, 9)
    ma_adjusted = adjust(matmp, minx, maxx, 1, 9)
    mb_adjusted = adjust(mbtmp, minx, maxx, 1, 9)
    
    if ma_adjusted == mg_adjusted or ma_adjusted == mb_adjusted:
        ma_adjusted += 0.01
        
    if mb_adjusted == mg_adjusted or mb_adjusted == ma_adjusted:
        mb_adjusted -= 0.01
    
    wg_adjusted = (9 - 1) * wg / (maxx - minx)
    wa_adjusted = (9 - 1) * wa / (maxx - minx)
    wb_adjusted = (9 - 1) * wb / (maxx - minx)
    

    
    
    miny = min(0, hg, ha, hb) 
    maxy = max(0, hg, ha, hb)
    
    if miny == maxy:
        maxy = miny + 0.001
    
    vscale = 8 / (maxy - miny)
    hg_adjusted = adjust(hg, miny, maxy, 1, 9)
    ha_adjusted = adjust(ha, miny, maxy, 1, 9)
    hb_adjusted = adjust(hb, miny, maxy, 1, 9)
    
    zero_height = 0
    if miny < 0:
        zero_height = adjust(0, miny, maxy, 1, 9)
    
   
   

    
    # Read in the file
    with open(template_filename, 'r') as file :
         filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('<!mg>', str(mg))
    filedata = filedata.replace('<!hg>', str(hg))
    filedata = filedata.replace('<!wg>', str(wg))
    filedata = filedata.replace('<!ma>', str(ma))
    filedata = filedata.replace('<!ha>', str(ha))
    filedata = filedata.replace('<!wa>', str(wa))
    filedata = filedata.replace('<!mb>', str(mb))
    filedata = filedata.replace('<!hb>', str(hb))
    filedata = filedata.replace('<!wb>', str(wb))
    
    filedata = filedata.replace('<!mg_adjusted>', str(mg_adjusted))
    filedata = filedata.replace('<!hg_adjusted>', str(hg_adjusted))
    filedata = filedata.replace('<!wg_adjusted>', str(wg_adjusted))
    filedata = filedata.replace('<!ma_adjusted>', str(ma_adjusted))
    filedata = filedata.replace('<!ha_adjusted>', str(ha_adjusted))
    filedata = filedata.replace('<!wa_adjusted>', str(wa_adjusted))
    filedata = filedata.replace('<!mb_adjusted>', str(mb_adjusted))
    filedata = filedata.replace('<!hb_adjusted>', str(hb_adjusted))
    filedata = filedata.replace('<!zero_height>', str(zero_height))
    
    filedata = filedata.replace('<!wb_adjusted>', str(wb_adjusted))

    filedata = filedata.replace('<!hscale>', str(hscale))
    filedata = filedata.replace('<!vscale>', str(vscale))

    filedata = filedata.replace('<!pneo>', str(pneo))
    filedata = filedata.replace('<!psubfunc>', str(psubfunc))
    filedata = filedata.replace('<!pcons>', str(pcons))
    filedata = filedata.replace('<!pspec>', str(pspec))
    filedata = filedata.replace('<!ppseudo>', str(ppseudo))
    filedata = filedata.replace('<!ppseudoboth>', str(ppseudoboth))

    filedata = filedata.replace('<!gcolor>', gcolor)
    filedata = filedata.replace('<!acolor>', acolor)
    filedata = filedata.replace('<!bcolor>', bcolor)

    # Write the file out again
    with open(template_filename_out, 'w') as file:
         file.write(filedata)

    command = "pdflatex -shell-escape \"" + template_filename_out + "\""
    os.system(command)

    pdffilename = template_filename_out.replace(".tex", "") + "-figure0.pdf"

    if outputfilename != '':
        shutil.copyfile(pdffilename, outputfilename) 
    else:
        outputfilename = pdffilename
    
    if startfile:
        os.startfile(outputfilename)



if __name__ == "__main__":
   main(sys.argv[1:])






