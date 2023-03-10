import os 
from multiprocessing import Pool

#converts all pdfs in current working directory to png.  Requires inkscape 
#example
# cd out 
# python ../batchconvert.py

def convertfile(f):
    cmd = "inkscape --pdf-poppler --pdf-page=1 --export-type=svg --export-text-to-path --export-area-drawing --export-filename=" + f + ".svg --export-dpi=200 " + f
    os.system(cmd)        
    print(cmd)
    cmd = "inkscape --pdf-poppler --pdf-page=1 --export-text-to-path --export-area-drawing --export-filename=" + f.replace(".pdf", ".png") + " --export-dpi=200 " + f
    os.system(cmd)
    print(cmd)
    os.remove(f + ".svg")

if __name__ == '__main__':    
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    
    with Pool(10) as p:
       p.map(convertfile, files)
#convertfile(files[0])