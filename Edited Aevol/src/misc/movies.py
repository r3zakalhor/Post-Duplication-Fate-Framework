#!/usr/bin/env python

import os
import numpy
import sys
import getopt
import png
import glob

def treat_generation(gfile,namepng,xsize,ysize,minv=0.,maxv=1.,colors=(0,255,0)):
    (rscale,gscale,bscale)=colors
    table=numpy.zeros((xsize,ysize))
    f = open(gfile,'r')
    for line in f:
        if len(line.split())==3:
            [x,y,z]=line.split()
            value=float(z)
            if value<minv:
                value=minv
            if value>maxv:
                value=maxv
            table[int(x),int(y)]=value
    f.close()
    pix=[tuple(reduce(lambda x,y:x+y, [(int((j-minv)/(maxv-minv)*rscale),int((j-minv)/(maxv-minv)*gscale),int((j-minv)/(maxv-minv)*bscale)) for j in i])) for i in table]

    pngfile=open(namepng,'wb')
    wr=png.Writer(xsize,ysize)
    wr.write(pngfile,pix)
    pngfile.close()
        
def printhelp():
    name=os.path.basename(sys.argv[0])
    print name +': a script to produce a movie from an aevol simulation with dumps activated.'
    print 'You need to have ffmpeg installed and in your path.'
    print 'Usage:'
    print '\t' + name + ' -i inputname [-s firstgen] [-e lastgen] [-p patchsize] [-m minvalue] [-M maxvalue] [-c color] [-o] [-k] sim1 sim2 .. simN'
    print 'or'
    print '\t' + name + ' -h'
    print 'Example:'
    print '\t' + " ./movies.py -i fitness_metabolic -p 40 -c '255 255 255' -m 0 -M 1 -s 30 -e 100 ~/myaevolsimulation"
    print 'See source code for a ful description of parameters and default values'
    exit(0)


opts, args = getopt.getopt(sys.argv[1:], "i:s:e:p:m:M:c:kho", ["inputname","start","end","patchsize","minvalue","maxvalue","color","keep","help","overwrite"])

start=None
end=None
inputname=None
patchsize=20
minvalue=0.
maxvalue=1.
colors=(0,255,0)
keepinter=False
overwrite=False

for o, a in opts:
    if o=='-i' or o=='--inputname': # Name of the output to treat (basename of the matching dump files, for example 'fitness_metabolic'). 
        inputname=a
    elif o=='-s' or o =='--start': # (optional) First generation to treat
        start=int(a)
    elif o=='-e' or o=='--end': # (optional) Last generation to treat
        end=int(a)
    elif o=='-p' or o=='--patchsize': # (optional) Size of each location in pixels (default 20)
        patchsize=int(a)
    elif o=='-m' or o=='--minvalue': # (optional) Value that will be scaled to 0 (black)
        minvalue=float(a)
    elif o=='-M' or o=='--maxvalue': # (optional) Value that will be scaled to 1 (full color specified by -c option)
        maxvalue=float(a)
    elif o=='-c' or o=='--color': # (optional) Color to give to maxvalue
        if (a=='red'):
            colors=(255,0,0)
        elif (a=='green'):
            colors=(0,255,0)
        elif (a=='blue'):
            colors=(0,0,255)
        elif len(a.split())==3:
            colors=map(lambda x: int(x), a.split())
        else:
            print 'Error: unrecognized argument given to -c option'
            exit(-1)
    elif o=='-k' or o=='--keep': # do not erase png files
        keepinter=True
    elif o=='-o' or o=='--overwrite': # re-generate png files if already existing
        overwrite=True
    elif o=='-h' or o=='--help':
        printhelp()
    else:
        print 'Error: unrecognized option ' + o
        exit(-1)

if not inputname:
    print 'Error: you must specifiy option -i (--inputname)'
    exit(-1)

smallestgen=1000000
biggestgen=None
warned=False

for arg in args: # If several folders are given as arguments treat all of them separatly
    first=True
    for filename in os.listdir(arg+'/stats/dump/'):
        if filename.startswith(inputname):
            ngen=int(filename.split('_').pop().split('.')[0])
            if ngen>biggestgen:
                biggestgen=ngen
            if ngen<smallestgen:
                smallestgen=ngen
            if first: # Use this first file to infer size of the grid
                mx=-1
                my=-1
                first=False
                temp=open(arg+'/stats/dump/'+inputname+'_%06d.out'%ngen)
                for line in temp:
                    if len(line)<=3 or line.startswith('#'):
                        continue
                    x=int(line.split()[0])
                    y=int(line.split()[0])
                    if x>mx:
                        mx=x
                    if y>my:
                        my=y
                temp.close()
                print "Detected x size: %d, detected y size: %d"%(mx+1,my+1)
            if (start and ngen<start) or (end and ngen>end): # do not treat generations that are not in given interval if any
                continue
            namepng='_%s%06d.png'%(inputname,ngen)
            if not os.path.exists(namepng) or overwrite:
                treat_generation(arg+'/stats/dump/'+inputname+'_%06d.out'%ngen,namepng,mx+1,my+1,minvalue,maxvalue,colors)
            elif not warned:
                print ''
                print "Warning: already existing png files found, we will use them instead of re-generating. If you want to re-generate them (to use different options), you first need to delete existing png files or to use option --overwrite"
                print ''
                warned=True
    if not start:
        start=smallestgen
    if not end:
        end=biggestgen
    os.system("ffmpeg -start_number " + str(start) + " -i _" + inputname + "%06d.png -vframes " + str(end-start) + " -r 6  -vcodec png -s "+str((mx+1)*patchsize)+"x"+str((my+1)*patchsize)+" -sws_flags neighbor -sws_dither none " + arg + '/' + inputname + ".mov")
    if not keepinter:
        os.system('rm _' + inputname + '*.png')

