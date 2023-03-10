#!/usr/bin/env python
# encoding=utf8

import os, time, datetime
from tempfile import mkstemp
import copy
import hashlib
import shutil

from threading import Thread
from execo import Put, Remote, Get, sleep, default_connection_params, Host, format_date, format_duration, SshProcess
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


class raevol_matrix(Engine):

    def __init__(self):
        """Overloading class initialization with parent and adding options"""
        super(raevol_matrix, self).__init__()
        self.options_parser.set_usage("usage: %prog ")
        self.options_parser.set_description("Execo Engine that can be used to" + \
                "perform automatic virtual machines experiments")
        self.options_parser.add_option("-n", dest="n_nodes",
                    help="maximum number of nodes used",
                    type="int",
                    default=200)
        self.options_parser.add_option("-w", dest="walltime",
                    help="walltime for the reservation",
                    type="string",
                    default="02:00:00")
        self.options_parser.add_option("-j", dest="oargrid_job_id",
                    help="oargrid_job_id to relaunch an engine",
                    type=int, default=None)
        self.options_parser.add_option("-k", dest="keep_alive",
                    help="keep reservation alive ..",
                    action="store_true")
        self.options_parser.add_option("-u", dest="selected_cluster",
                    help="run on a specific cluster.",
                    type="string",
                    default="taurus")
        self.options_parser.add_option("-o", dest="outofchart",
                    help="Run the engine outside days",
                    action="store_true")
	self.options_parser.add_option("-s", dest="storage5k_job_id",
                    help="storage5k_job_id to store the data",
                    type=int)

    def run(self):
        """ """
        if self.options.oargrid_job_id is not None:
            self.oar_job_id = self.options.oargrid_job_id
        else:
            self.oar_job_id = None
            
        self.list_of_clusters = ['parasilo','paravance','parapluie','paranoia']

        try:
            # Creation of the main iterator which is used for the first control loop.
            self.define_parameters()
	    self.working_dir = '/data/jorouzaudcornabas_'+str(self.options.storage5k_job_id)
        
            job_is_dead = False
            # While there are combinations to treat
            while len(self.sweeper.get_remaining()) > 0:
                # If no job, we make a reservation and prepare the hosts for the experiments
                if self.oar_job_id is None:
                    self.submit_all_available_best_effort(self.list_of_clusters, self.options.walltime)
                    # self.make_reservation_local()
                # Wait that the job starts
                logger.info('Waiting that the job start '+str(self.oar_job_id))
                wait_oar_job_start(self.oar_job_id)
                # Retrieving the hosts and subnets parameters
                self.hosts = get_oar_job_nodes(self.oar_job_id)
                # Hosts deployment and configuration
                default_connection_params['user'] = 'jorouzaudcornabas'

                logger.info("Start hosts configuration")
                ex_log.setLevel('INFO')
                #===============================================================
                # deployment = Deployment(hosts = self.hosts, 
                #             env_file='/home/sirimie/env/mywheezy-x64-base.env')
                # self.hosts, _ = deploy(deployment)
                #===============================================================
                if len(self.hosts) == 0:
                    break

                # Initializing the resources and threads
                available_hosts = self.hosts
                        
                threads = {}
                
                # Creating the unique folder for storing the results
                comb_dir = self.result_dir + '/logs'
                if  not os.path.exists(comb_dir):
					os.mkdir(comb_dir)

                logger.info("Starting the thread "+str(self.is_job_alive())+" "+str(len(threads.keys())))
                # Checking that the job is running and not in Error
                while self.is_job_alive() or len(threads.keys()) > 0:
                    job_is_dead = False

                    while self.options.n_nodes > len(available_hosts):
                        tmp_threads = dict(threads)
                        for t in tmp_threads:
                            if not t.is_alive():
                                available_hosts.append(tmp_threads[t]['host'])
                                del threads[t]
                        sleep(5)
                        if not self.is_job_alive():
                            job_is_dead = True
                            break
                    if job_is_dead:
                        break

                    # Getting the next combination
                    comb = self.sweeper.get_next()
                    if not comb:
                        while len(threads.keys()) > 0:
                            tmp_threads = dict(threads)
                            for t in tmp_threads:
                                if not t.is_alive():
                                    del threads[t]
                            logger.info('Waiting for threads to complete')
                            sleep(20)
                        break
                    
                    host = available_hosts[0]
                    available_hosts = available_hosts[1:]
                    logger.info("Launching thread")
                    t = Thread(target=self.workflow,
                               args=(comb, host, comb_dir))
                    threads[t] = {'host': host}
                    t.daemon = True
                    t.start()

                if not self.is_job_alive():
                    job_is_dead = True

                if job_is_dead:
                    self.oar_job_id = None
            
                
        finally:
            if self.oar_job_id is not None:
                if not self.options.keep_alive:
                    logger.info('Deleting job')
                    oardel([self.oar_job_id])
                else:
                    logger.info('Keeping job alive for debugging')


    def define_parameters(self):
        """ """
        parameters = {
	  'seed' : [51456165, 33263658, 7158785, 456847894, 1223144, 878944, 121145, 3587842],
	  'mutation' : ['5e-4','1e-4','5e-5','5e-6'],
	  'env' : ['const','lat_3','lat_all'],
	  'selection' : [750,2000,4000]
        }
        sweeps = sweep(parameters)
        self.sweeper = ParamSweeper(os.path.join(self.result_dir, "sweeps"), sweeps)
        logger.info('Number of parameters combinations %s', len(self.sweeper.get_remaining()))

    def make_reservation(self):
        """ """
        logger.info('Performing reservation')
        starttime = int(time.time() + timedelta_to_seconds(datetime.timedelta(minutes=1)))
        planning = get_planning(elements=[self.options.selected_cluster],
                            starttime=starttime)
        slots = compute_slots(planning, self.options.walltime)
        wanted = {self.options.selected_cluster: 0}
        start_date, end_date, resources = find_first_slot(slots, wanted)
        wanted[self.options.selected_cluster] = resources[self.options.
                                                    selected_cluster]
        actual_resources = distribute_hosts(resources, wanted)

        job_specs = get_jobs_specs(actual_resources, name='Aevol_diff_area')
        logger.info("try to reserve " + str(actual_resources))
        self.oargrid_job_id , _= oargridsub(job_specs,
                          walltime = end_date - start_date,
                          job_type = ['besteffort"'
                                      'allow_classic_ssh'])
        logger.info("Reservation done")

    def make_reservation_local(self):
        """Perform a reservation of the required number of nodes, with 4000 IP.
        """
        logger.info('Performing reservation')
        starttime = int(time.time() + timedelta_to_seconds(datetime.timedelta(minutes=1)))
        endtime = int(starttime + timedelta_to_seconds(datetime.timedelta(days=3,
                                                                 minutes=1)))
        self.cluster = self.options.selected_cluster
        startdate, n_nodes = self._get_nodes(starttime, endtime)
        
        while not n_nodes:
            logger.info('No enough nodes found between %s and %s, ' + \
                        'increasing time window',
                        format_date(starttime), format_date(endtime))
            starttime = endtime
            endtime = int(starttime + timedelta_to_seconds(datetime.timedelta(days=3,
                                                                minutes=1)))
            startdate, n_nodes = self._get_nodes(starttime, endtime)
            if starttime > int(time.time() + timedelta_to_seconds(
                                                datetime.timedelta(weeks=6))):
                logger.error('There are not enough nodes on %s for your ' + \
                             'experiments, abort ...', self.cluster)
                exit()
        startdate = []
        jobs_specs = get_jobs_specs({self.cluster: n_nodes},
                                    name=self.__class__.__name__)
        sub = jobs_specs[0][0]
        tmp = str(sub.resources).replace('\\', '')
        sub.resources = tmp.replace('"', '')
        sub.walltime = self.options.walltime
        sub.additional_options = '-t allow_classic_ssh -t besteffort'
        (self.oar_job_id, self.frontend) = oarsub(jobs_specs)[0]
        logger.info('Startdate: besteffort, n_nodes: %s', 
                    str(n_nodes))
        
    def _get_nodes(self, starttime, endtime):
        """ """
        planning = get_planning(elements=[self.cluster],
                                starttime=starttime,
                                endtime=endtime,
                                out_of_chart=self.options.outofchart)
        slots = compute_slots(planning, self.options.walltime)
        startdate = slots[0][0]
        i_slot = 0
        n_nodes = slots[i_slot][2][self.cluster]
        logger.info("nodes %s in %s at %s", str(n_nodes),
                    str(self.cluster), format_date(startdate))
        while n_nodes < self.options.n_nodes:
            logger.debug(slots[i_slot])
            startdate = slots[i_slot][0]
            n_nodes = slots[i_slot][2][self.cluster]
            i_slot += 1
            if i_slot == len(slots) - 1:
                return False, False
        return startdate, n_nodes
      

    def workflow(self, comb, host, comb_dir):
        """ """
        comb_ok = False
        thread_name = style.Thread(str(host).split('.')[0]) + ': '
        logger.info(thread_name + 'Starting combination ' + slugify(comb))

        try:
	    self.export = "source ~/aevol_binary/intel/linux/mkl/bin/mklvars.sh intel64; "
	    
	    bucketname = self.working_dir+'/raevol_5_mut_lat/'+slugify(comb)+'/'
	      
	    if os.path.isdir(bucketname) and os.path.exists(bucketname+'/last_gener.txt'):
	      logger.info(thread_name + "Resuming AEVOL from NFS backup")
	      
	      gen_file = open(bucketname+'/last_gener.txt', 'r')                                                                                                                                                                                   
                                                                                                                                                                                                                                                   
              last_gen = gen_file.read()
	      
	      if int(last_gen) < 500000:
                logger.info(thread_name + "Resuming AEVOL Run from "+str(int(last_gen)))
                rem = Remote(self.export+'cd '+bucketname+'; /home/jorouzaudcornabas/aevol_binary/aevol/src/aevol_run -p 16'
		      + ' -e 300000 -r '+last_gen+' >> aevol_run.log',
			  [host]).run()
                if rem.ok:
                    comb_ok = True
              else:
                  comb_ok = True
	    else:
	      Remote('mkdir -p '+bucketname,[host]).run()
		
	      param_file = '/home/jorouzaudcornabas/aevol_binary/execo/mut_lat/param_tmpl.in'
	      
	      logger.info(thread_name + 'Generate config file ' + param_file)
	      
	      f_template = open(param_file)
	      fd, outfile = mkstemp(dir='/tmp/', prefix=slugify(comb) + '_param')
	      f = os.fdopen(fd, 'w')
		
	      for line in f_template:
		line = line.replace('SEED_NUMBER',str(comb['seed']))
		line = line.replace('FUZZY_VERSION',str(comb['fuzzy']))
		if comb['move']:
                    line = line.replace('FIRST_GAUSSIAN_MEDIAN','0.25') 
                    line = line.replace('THIRD_GAUSSIAN_MEDIAN','0.65') 
                else:
                    line = line.replace('FIRST_GAUSSIAN_MEDIAN','0.2') 
                    line = line.replace('THIRD_GAUSSIAN_MEDIAN','0.6')
                line = line.replace('GAUSSIAN_HEIGHT',str(comb['height']))
		f.write(line)
		
	      f_template.close()
	      f.close()
	      
	      put_file = Put([host], [outfile],remote_location=bucketname).run()
	      if not put_file.ok:
		exit()
		
	      os.remove(outfile)
	      
	      Remote('cd '+bucketname+'; cp ' + outfile.split('/')[-1] + ' param.in',
			  [host]).run()
		
	      logger.info(thread_name + "Launching AEVOL Create")
	      Remote(self.export+'cd '+bucketname+'; /home/jorouzaudcornabas/aevol_diff_area/aevol/src/aevol_create > aevol_create.log',
			  [host]).run()
	      
	      logger.info(thread_name + "Launching AEVOL Run")
	      rem = Remote(self.export+'cd '+bucketname+'; /home/jorouzaudcornabas/aevol_diff_area/aevol/src/aevol_run -p 16 -n 500000 > aevol_run.log',
			  [host]).run()
              if rem.ok:
                    comb_ok = True
            
	    logger.info(thread_name + 'Get results ' + comb_dir + "/" + slugify(comb))
  
            #try:
                #os.mkdir(comb_dir + "/" + slugify(comb))
            #except:
                #logger.warning(thread_name +
                    #'%s already exists, removing existing files', comb_dir + "/" + slugify(comb))
		#shutil.rmtree(comb_dir+ "/" + slugify(comb))
		#try:
		  #os.mkdir(comb_dir + "/" + slugify(comb))
		#except:
		  #logger.warning(thread_name +
                    #'%s already exists, recreating directory', comb_dir + "/" + slugify(comb))

	    #get_results = Get([host], [bucketname+ "/aevol_create.log", bucketname+ "/aevol_run.log", bucketname+'/stats/'],
                            #local_location=comb_dir + "/" + slugify(comb)).run()

            #for p in get_results.processes:
                #if not p.ok:
                    #logger.error(thread_name +
                        #': Unable to retrieve the files for combination %s',
                        #slugify(comb))
                    #exit()

        finally:
            if comb_ok:
                self.sweeper.done(comb)
		# shutil.rmtree(bucketname)
                logger.info(thread_name + ': ' + slugify(comb) + \
                             ' has been done')
            else:
                self.sweeper.cancel(comb)
                logger.warning(thread_name + ': ' + slugify(comb) + \
                            ' has been canceled')
        logger.info(style.step('%s Remaining'),
                        len(self.sweeper.get_remaining()))

    def is_job_alive(self):
        rez=get_oar_job_info(self.oar_job_id)
        if rez['state'] == 'Error':
          return False

        if (rez["start_date"]+rez["walltime"] > time.time()):
            return True
        else:
            return False
        
    def get_immediately_available_nodes(self,list_of_clusters, walltime, ignore_besteffort = False):
        planning = get_planning(list_of_clusters, ignore_besteffort=False)
        slots = compute_slots(planning, walltime)
        wanted = { cluster: 0 for cluster in list_of_clusters }
        start_date, end_date, resources = find_first_slot(slots, wanted)
        actual_resources = { resource: n_nodes
                            for resource, n_nodes in resources.iteritems()
                            if resource in list_of_clusters and n_nodes > 0 }
        return start_date, end_date, actual_resources

    def submit_all_available_best_effort(self,list_of_clusters, walltime):
        start_date, end_date, resources = self.get_immediately_available_nodes(list_of_clusters, walltime, ignore_besteffort = False)
        job_specs = get_jobs_specs(resources)
        for j,f in job_specs:
            j.job_type = "besteffort"
            j.additional_options = "-t allow_classic_ssh"
        #jobs = oarsub(job_specs)
        (self.oar_job_id, self.frontend) = oarsub(job_specs)[0]
        print job_specs
        return self.oar_job_id


if __name__ == "__main__":
    engine = raevol_matrix()
    engine.start()
