#!/usr/local/bin/python

import os
import sys
import stat
import glob
import random
import tarfile
import shutil

#~ if len(sys.argv) != 2:
  #~ sys.stdout.write( 'ERROR : missing \'generation number\' operand\n' )
#~ else:
  #~ num_gener = sys.argv[1]
  
  
if 1 == 1:
  ################################################################################################################
  ## Retreive the simulation name, the number of the last saved generation and the SPS directory of the experiment
  ################################################################################################################
  # Get simulation name
  sim_name_file = open( 'sim_name.txt', 'r' )
  sim_name = (sim_name_file.readline()).rstrip( '\n' )
  sim_name_file.close()
  
  # Get number of the last saved generation
  last_gener_file = open( 'last_gener.txt', 'r' )
  last_gener = (last_gener_file.readline()).rstrip( '\n' )
  last_gener_file.close()
  
  # Get SPS dir
  sps_dir_file = open( 'SPS_dir.txt', 'r' )
  sps_dir = (sps_dir_file.readline()).rstrip( '\n' )
  sps_dir_file.close()
  
  
  
  ####################################################################################
  ## Open the general archive: a tgz file to put our "big backup"
  ####################################################################################
  general_archive_name = sim_name + '_' + last_gener + '.tgz'
  general_archive = tarfile.open( general_archive_name, 'w:gz' )
  
  
  
  ##################################################################################
  ## Save a copy of all the stat files in a targzipped "stat_<last_gener>" directory
  ##################################################################################
  # mkdir stat_<last_gener>
  os.makedirs( 'stat_' + last_gener )
  
  # cp stat_*.out stat_<last_gener>
  for filename in glob.glob( os.path.join( '.', 'stat_*.out' ) ):
    shutil.copy( filename, 'stat_' + last_gener + '/' + filename )
  
  # tar zcf stat_<last_gener>.tgz stat_<last_gener>
  stat_archive = tarfile.open( 'stat_' + last_gener + '.tgz', 'w:gz' )
  stat_archive.add( 'stat_' + last_gener )
  stat_archive.close()
  
  # rm -rf stat_<last_gener>
  shutil.rmtree( 'stat_' + last_gener )
  
  
  
  ####################################################################################
  ## Put stat_<last_gener>.tgz in the general archive
  ####################################################################################
  general_archive.add( 'stat_' + last_gener + '.tgz' )
  general_archive.add( 'backup' )
  general_archive.add( 'tree' )
  
  # add *.txt to general archive
  for filename in glob.glob( os.path.join( '.', '*.txt' ) ):
    general_archive.add( filename )
  general_archive.close()
  
  # rm (stat, backup, tree)_<last_gener>.tgz
  os.remove( 'stat_' + last_gener + '.tgz' )
  
  # rm -f backup/* (keeping only the last backup)
  os.rename( 'backup', 'backup.tmp' )
  os.makedirs( 'backup' )
  os.rename( 'backup.tmp/gen_' + '{0:0>6}'.format( last_gener ) + '.ae', 'backup/gen_' + '{0:0>6}'.format( last_gener ) + '.ae' )
  shutil.rmtree( 'backup.tmp' )
  # rm -f tree/*
  shutil.rmtree( 'tree' )
  os.makedirs( 'tree' )
  
  # move general archive to SPS
  shutil.copy( general_archive_name, sps_dir )
  os.remove( general_archive_name )
  