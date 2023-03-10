#! /usr/bin/python

import os
import shutil
import sys
import stat
import random
import tarfile

# ****************************** SETTINGS ******************************
exp_name            = 'my_exp'
nb_seeds            = 5
rearrangement_rates = ['1e-6', '2e-6', '5e-6', '1e-5', '2e-5', '5e-5', '1e-4']
mutation_rates      = ['1e-6', '2e-6', '5e-6', '1e-5', '2e-5', '5e-5', '1e-4']
nb_gener            = 50000

aevol_version       = 'aevol-3.1.0'
aevol_dir           = '/sps/liris/aevol/' + aevol_version
aevol               = aevol_dir + '/src/aevol_IN2P3'
backup_script       = aevol_dir + '/src/misc/make_big_backup_IN2P3.py'
exp_dir             = '/sps/liris/dparsons/sim/' + exp_name
# **********************************************************************

# Anti overwrite
if os.path.exists(exp_name):
  sys.exit( 'prepare_experiment: creating experiment \'' + exp_name + '\': Directory exists\n' )
if not os.path.exists( 'param.in' ):
  sys.exit( 'prepare_experiment: creating experiment \'' + exp_name + '\': No param.in file provided\n' )
  


# Create the experiment directory and its 'scripts' directory
os.makedirs( exp_name )
os.chdir( exp_name )
os.makedirs( 'scripts' )
submit_all_sh_file_name = 'scripts/submit_all.sh'
submit_all_sh_file      = open( submit_all_sh_file_name, 'w' )
aevol_version_info_file = open( 'aevol_version_info.txt', 'w' )
aevol_version_info_file.write( aevol_version + '\n' )
aevol_version_info_file.close()
nb_gener_info_file = open( 'nb_gener.txt', 'w' )
nb_gener_info_file.write( str(nb_gener) + '\n' )
nb_gener_info_file.close()


# Initialize random generator
random.seed()


# Create the directory, param file and script for each simulation
for rear_rate in rearrangement_rates:
  for mut_rate in mutation_rates:
    for i in range(1, nb_seeds+1):
      
      # **********************************************************************
      # 1) Set a few variables
      # **********************************************************************
      sim_dir           = 'rear_' + rear_rate + 'mut_' + mut_rate + '/seed' + str(i)
      sps_sim_dir       = exp_dir + '/' + sim_dir
      sim_name          = 'rear_' + rear_rate + 'mut_' + mut_rate + '_seed' + str(i)
      job_name          = sim_name
      
      
      # **********************************************************************
      # 2) Create simulation directory and parameter file (using template)
      # **********************************************************************
      os.makedirs( sim_dir )
      
      template_param_file = open('../param.in', 'r')
      new_param_file = open( sim_dir + '/param.in', 'w' )
      
      for cur_line in template_param_file:
        if cur_line.find('SEED') == 0:
          new_param_file.write('SEED ' + str(random.randint(1, 10000000)) + '\n')
        elif cur_line.find('POINT_MUTATION_RATE') == 0:
          new_param_file.write('POINT_MUTATION_RATE   ' + mut_rate + '\n')
        elif cur_line.find('SMALL_INSERTION_RATE') == 0:
          new_param_file.write('SMALL_INSERTION_RATE  ' + mut_rate + '\n')
        elif cur_line.find('SMALL_DELETION_RATE') == 0:
          new_param_file.write('SMALL_DELETION_RATE   ' + mut_rate + '\n')
        elif cur_line.find('DUPLICATION_RATE') == 0:
          new_param_file.write('DUPLICATION_RATE      ' + rear_rate + '\n')
        elif cur_line.find('DELETION_RATE') == 0:
          new_param_file.write('DELETION_RATE         ' + rear_rate + '\n')
        elif cur_line.find('TRANSLOCATION_RATE') == 0:
          new_param_file.write('TRANSLOCATION_RATE    ' + rear_rate + '\n')
        elif cur_line.find('INVERSION_RATE') == 0:
          new_param_file.write('INVERSION_RATE        ' + rear_rate + '\n')
        else:
          new_param_file.write(cur_line)
          
      new_param_file.close()
      template_param_file.close()
      
      
      # **********************************************************************
      # 3) Build the scripts that will be launched on the cluster
      # **********************************************************************
      script_file_name  = 'scripts/' + job_name + '.sh'
      script_file       = open(script_file_name, 'w')
      script_file.write('#!/bin/sh\n')
      script_file.write('#PBS -N ' + job_name + ' # Job name\n')
      script_file.write('#PBS -q T                             # Execution Class\n')
      script_file.write('#PBS -l platform=LINUX                # Execution Platform\n')
      script_file.write('#PBS -l u_sps_liris                   # We use the Semi Permanent Storage\n')
      script_file.write('#PBS -l T=4286000                     # Duration (in normalized units)\n')
      script_file.write('#PBS -l M=4096MB                      # Virtual Memory in MB\n')
      script_file.write('#PBS -l scratch=30GB                  # Scratch size in MB (disk space on calculator)\n')
      script_file.write('#PBS -l spool=500KB                   # Spool size in KB (space for stdout and stderr)\n')
      script_file.write('#PBS -l model=Xeon                    # Processor\n')
      script_file.write('#PBS -eo                              # Redirect stderr to stdout\n')
      script_file.write('#PBS -mb                              # Send Mail on Begin\n')
      script_file.write('#PBS -me                              # Send Mail on End\n')
      script_file.write('\n\n\n')
      script_file.write('# Simulation directories -> Where you have your input files and where you want your output files to be copied\n')
      script_file.write('SPS_BASE_DIR="' + exp_dir + '"\n')
      script_file.write('SIM_DIR="' + sim_dir + '"\n')
      script_file.write('SIM_NAME="' + sim_name + '"\n')
      script_file.write('SPS_SIM_DIR="' + sps_sim_dir + '"\n')
      script_file.write('\n')
      script_file.write('# Binary program path.\n')
      script_file.write('EXEC="' + aevol + '"\n')
      script_file.write('\n')
      script_file.write('# Create the simulation dir and cd into it\n')
      script_file.write('mkdir -p $SIM_DIR             # Create the simulation subdirectory on the calculator\n')
      script_file.write('cd $SIM_DIR                   # cd on the calculator\n')
      script_file.write('\n')
      script_file.write('# Copy the param.in file\n')
      script_file.write('cp $SPS_SIM_DIR"/param.in" .  # Copy param.in from SPS to the calculator\n')
      script_file.write('\n')
      script_file.write('# Write the sim_name.txt, SPS_dir.txt and cpu_info.txt files\n')
      script_file.write('echo ' + exp_name + ' > exp_name.txt\n')
      script_file.write('echo ' + sim_name + ' > sim_name.txt\n')
      script_file.write('echo ' + sps_sim_dir + ' > SPS_dir.txt\n')
      script_file.write('cp /proc/cpuinfo cpu_info.txt\n')
      script_file.write('\n')
      script_file.write('# Copy the big_backup script\n')
      script_file.write('cp ' + backup_script + ' make_big_backup.py\n')
      script_file.write('\n')
      script_file.write('# Run the program "$EXEC" for nb_gener generations\n')
      script_file.write('$EXEC -n ' + str(nb_gener) + ' > /dev/null\n')
      script_file.write('\n')
      script_file.write('# Make a .tar.gz from the result (with the directory structure $SIM_DIR)\n')
      script_file.write('cd ../..\n')
      script_file.write('tar zcf $SIM_NAME.tgz $SIM_DIR --remove-files\n')
      script_file.write('mv  $SIM_NAME.tgz $SPS_BASE_DIR/\n')
      script_file.close()
      
      submit_all_sh_file.write('qsub ' + exp_dir + '/' + script_file_name + '\n')
      os.chmod(script_file_name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)


submit_all_sh_file.close()
os.chmod(submit_all_sh_file_name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)


# ***************************************************************
# 4) Put the freshly created files in a tarball (and delete them)
# ***************************************************************
os.chdir( '..' )
dirs_archive = tarfile.open( exp_name + '.tar.gz', 'w:gz')
dirs_archive.add( exp_name )
dirs_archive.close()


shutil.rmtree( exp_name )
