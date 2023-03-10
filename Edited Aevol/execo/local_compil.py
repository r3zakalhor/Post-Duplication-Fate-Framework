#!/usr/bin/env python
# encoding=utf8

import os, time, datetime
from tempfile import mkstemp
import copy
import hashlib
import shutil

from threading import Thread
from execo import Put, Remote, Get, sleep, default_connection_params, Host, format_date, format_duration, SshProcess, Process
from execo.log import style
from execo import logger as ex_log
from execo.time_utils import timedelta_to_seconds
from execo_g5k import get_host_attributes, get_planning, compute_slots, \
    find_first_slot, distribute_hosts, get_jobs_specs, g5k_configuration, \
    wait_oargrid_job_start, oargridsub, oargriddel, get_oargrid_job_nodes, \
    Deployment, deploy, get_oargrid_job_info, get_host_cluster, get_cluster_site, get_host_site, get_g5k_sites, get_oar_job_info, \
    get_cluster_site, OarSubmission, \
    oarsub, get_oar_job_nodes, wait_oar_job_start, oardel, get_host_attributes
from execo_engine import Engine, logger, ParamSweeper, sweep, slugify


class compil_aevol(Engine):

    def __init__(self):
        """Overloading class initialization with parent and adding options"""
        super(compil_aevol, self).__init__()

    def run(self):
        """ """
        try:
            # Creation of the main iterator which is used for the first control loop.
            self.define_parameters()
	    
            # While there are combinations to treat
            while len(self.sweeper.get_remaining()) > 0:
	      comb = self.sweeper.get_next()
              if comb:
		self.workflow(comb)
	finally:
	  logger.info("Compilation DONE")
	  
    def define_parameters(self):
        """ """
        parameters = {
          'seed' : [51456165, 33263658, 7158785, 456847894, 1223144, 878944, 121145, 3587842, 2784564, 68984554],
          'mutation' : ['5e-5','1e-5','5e-6','1e-6','5e-7'],
          'env' : ['const','lat_3','lat_13','lat_123','lat_all'],
          'selection' : [500,1000,1500,2000]

        }
        sweeps = sweep(parameters)
        self.sweeper = ParamSweeper(os.path.join(self.result_dir, "sweeps"), sweeps)
        logger.info('Number of parameters combinations %s', len(self.sweeper.get_remaining()))

    def workflow(self, comb):
        """ """
        comb_ok = False
        logger.info(slugify(comb) + \
                             ' starts to compile')
        try:
	    export = "source /opt/intel/bin/compilervars.sh intel64; "
	    
	    src_directory = "/home/arrouan/workspace/aevol/git/world/aevol/"
	    
	    bin_directory = "/home/arrouan/workspace/aevol/compiled_binary/"
    
	    configure_option = "--with-tracing --without-x"
	    
	    if comb['blas'] == 'openblas':
	      configure_option += " --with-blas"
	    elif comb['blas'] == 'mkl':
	      configure_option += " --with-mkl"
	    elif comb['blas'] == 'atlas':
	      configure_option += " --with-atlas"
	        
	    if comb['experiment'] == 'raevol':
	      configure_option += " --with-raevol"
	      
	    if comb['compilator'] == 'intel':
	      configure_option += " CXX=icc"
	      
	    full_bin_directory = bin_directory + comb['experiment']+'_'+comb['compilator']+'_'+comb['parallel']+'_'+comb['blas']
	    
	    try:
	      os.mkdir(full_bin_directory)
            except:
	      for f in os.listdir(full_bin_directory):
		os.remove(full_bin_directory + "/" + f)
		
	    p = Process(export+'cd '+src_directory+'; autoreconf; ./configure '+configure_option+'; make clean; make; cp src/aevol_run '+full_bin_directory+'/; cp src/aevol_create '+full_bin_directory+'/')
	    p.shell = True
	    #
	    p.run()
	    
	    print p.stdout
	    
            comb_ok = True
        finally:
            if comb_ok:
                self.sweeper.done(comb)
	        logger.info(slugify(comb) + \
                             ' has been done')
            else:
                self.sweeper.cancel(comb)
                logger.warning(slugify(comb) + \
                            ' has been canceled')
        logger.info(style.step('%s Remaining'),
                        len(self.sweeper.get_remaining()))

    def is_job_alive(self):
        rez=get_oar_job_info(self.oar_job_id)
        if (rez["start_date"]+rez["walltime"] > time.time()):
            return True
        else:
            return False

if __name__ == "__main__":
    engine = compil_aevol()
    engine.start()
