How I generate all triangles in a set of directories:
for file in runs_7jan2023/*; do python run-maker.py --indir=$file --outdir=${file}_v5 --dotriangles --dopng --dov5; done


How I take the Concats directory and recalculate everything
rm -r Concat/*/concat*
for file in Concat/*; do python run-maker.py --indir=$file --outdir=${file} --dov45; done

If you have a dups_fates_probablities csv or a bunch of files of the form dups_fates_probablities.csv*, and want to generate all triangles figures, do the following:
- create a new directory and put all the csvs in it (there can be other files in the directory, they will be ignored, only filenames of the form dups_fates_probablities.csv* are considered)
- suppose the new directory is called mydir
- call python run-maker.py --indir=mydir --dotriangles
and wait 1-2 hours.
The concatenated csv and a triangles csv will be copied to mydir, the triangles will be added in mydir/trianglefigs
This generates pdfs.  If you want to generate pngs, do 
cd mydir/trianglefigs
python ../../batchconvert.py



The run_maker.py script takes as input a directory.  It scans the dups_fates_probablities.csv* files and concatenates their content, without repeating the header rows.
It creates two files: concat_mydir.csv, which is just a pure concatenation, and triangles_mydir.csv, which has the same elements by with columns rearranged to be used by triangle_maker.py.

This triangle_maker.py script takes a csv as input (columns must be the same as the provided tmp.csv)
Then for each line = 2 .. nblines, it outputs a pdf of the triangle on that line.
The file name of the pdf is the [line_number].pdf.  Thus the script will create 2.pdf, 3.pdf, and so on (note that we start on line 2.
The files go in the [outdir] argument (or cwd if unspecified).


Existing files will be mercilessly overwritten.
Example call:

python triangle_maker.py --batch=tmp.csv --outdir=out
