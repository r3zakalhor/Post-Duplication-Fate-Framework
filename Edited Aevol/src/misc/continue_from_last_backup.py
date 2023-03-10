#!/usr/local/bin/python
###!/usr/bin/env python3.1

import os
import sys
import stat
import glob
import random
import tarfile
import shutil



if len(sys.argv) != 2:
  sys.exit( 'USAGE : ' + sys.argv[0] + ' <end_of_job_file.tgz>\n' )

# Keep track of the Current Working Directory
my_cwd = os.getcwd()

# Enter the experiment base dir (if not already inside)
if sys.argv[1].find('/') != -1:
  part_string           = sys.argv[1].rpartition('/')
  exp_base_path         = part_string[0]
  end_of_job_file_name  = part_string[2]
  os.chdir( exp_base_path )
else:
 end_of_job_file_name   = sys.argv[1]
 
# Untar end_of_job_file
end_of_job_file = tarfile.open( end_of_job_file_name, 'r:gz')
end_of_job_file.extractall()
end_of_job_file.close()

# Reconstruct sim_dir
part_string   = end_of_job_file_name.rpartition( '_seed' )
part_string2  = part_string[2].rpartition( '.' )
sim_dir       = part_string[0] + '/seed' + part_string2[0]

script_dir = 'scripts'


# Get rid of the ending '/' in sim_dir and script_dir (if any)
sim_dir     = sim_dir.rstrip( '/' )
script_dir  = script_dir.rstrip( '/' )

# Get exp name
exp_name_file = open( sim_dir + '/exp_name.txt', 'r' )
exp_name = (exp_name_file.readline()).rstrip( '\n' )
exp_name_file.close()

# Get simulation name
sim_name_file = open( sim_dir + '/sim_name.txt', 'r' )
sim_name = (sim_name_file.readline()).rstrip( '\n' )
sim_name_file.close()

# Get number of the last saved generation
last_gener_file = open( sim_dir + '/last_gener.txt', 'r' )
last_gener = int( (last_gener_file.readline()).rstrip( '\n' ) )
last_gener_file.close()

# Get number of generations to do 
nb_gener_file = open( 'nb_gener.txt', 'r' )
nb_gener = int( (nb_gener_file.readline()).rstrip( '\n' ) )
nb_gener_file.close()


if last_gener >= nb_gener:
  sys.exit( sys.argv[1] + ': Simulation is finished\n' )



# **********************************************************************
# Create the "out_files.tgz" archive
# **********************************************************************
os.chdir( sim_dir )
out_files = tarfile.open( 'out_files.tgz', 'w:gz')
for filename in glob.glob( os.path.join( '.', '*.out' ) ):
  out_files.add( filename )
out_files.close()
os.chdir( '../..' )


# **********************************************************************
# Create the new script file
# **********************************************************************
os.chdir( script_dir )
old_script_file_name = sim_name + '.sh'
# NB : '{0:0>6}'.format( foo ) forces foo to span 6 digits (filled with 0s on the left)
new_script_file_name = sim_name + '_' + '{0:0>6}'.format( last_gener ) + '.sh'

old_script_file = open( old_script_file_name, 'r' )
new_script_file = open( new_script_file_name, 'w' )

# Copy the content of old_script_file into new_script_file with a few modifications : 
#   * Instead of copying the param.in file, we need to retreive the last backup and all the
#     .out files (mainly the stat files)
#   * The simulation must be run for <nb_gener> - <last_gener> simulation starting from 
#     generation <last_gener>
for cur_line in old_script_file:
  if cur_line.find( '# Copy the param.in file' ) == 0:
    new_script_file.write( '# Retreive backed files necessary to continue the simulation\n' )
  elif cur_line.find( 'cp $SPS_SIM_DIR"/param.in"' ) == 0:
    new_script_file.write( 'mkdir backup tree\n' )
    new_script_file.write( 'SPS_LAST_BACKUP=$SPS_SIM_DIR"/backup/gen_' + '{0:0>6}'.format( last_gener ) + '.ae"\n' )
    new_script_file.write( 'cp $SPS_LAST_BACKUP "backup/"\n' )
    new_script_file.write( 'cp $SPS_SIM_DIR"/out_files.tgz" .\n' )
    new_script_file.write( 'tar zxf out_files.tgz\n' )
    new_script_file.write( 'rm out_files.tgz\n' )
  elif cur_line.find( '$EXEC' ) == 0:
    new_script_file.write( '$EXEC -f backup/gen_' + '{0:0>6}'.format( last_gener ) + '.ae -n ' + str(nb_gener-last_gener) + ' > /dev/null\n' )
  else:
    new_script_file.write( cur_line )

old_script_file.close()
new_script_file.close()
os.chmod(new_script_file_name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH) # chmod 744

# **********************************************************************
#                           SUBMIT NEW JOB
# **********************************************************************
os.system( 'qsub ' + new_script_file_name )


# Go back to the directory we were in at the beginning
os.chdir( my_cwd )
