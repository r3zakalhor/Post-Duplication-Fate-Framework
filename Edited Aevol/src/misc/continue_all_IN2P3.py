#!/usr/local/bin/python

import os
import sys
import glob

for simu in glob.glob( os.path.join( '.', '*.tgz' ) ):
  simu_name = simu.rpartition( '.tgz' )[0]
  continue_script = sys.argv[0].rpartition( '/' )[0] + '/continue_from_last_backup.py'
  
  print( 'Preparing ' + simu_name + '...\n' )
  os.system( continue_script + ' ' + simu )
  os.system( 'mv '+ simu_name+'.tgz ' + simu_name+'.tgz.bak' )

