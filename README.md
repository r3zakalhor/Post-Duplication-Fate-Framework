# Post-Duplication Fate Framework




#### To infer post-duplication fates of an aevol simulation, we need to do the following steps in order:

## 1) Compile edited Aevol:

#### Aevol is a computational platform that allows for the study and manipulation of populations of digital organisms evolving under different conditions. To compile Aevol use:

  - $cd Edited Aevol
  - $mkdir build
  - $cd build
  - $cmake -DCMAKE_BUILD_TYPE=Release ..
  - $make



In this process, you'll need development packages. Here are the basic dependencies on a Debian system:

- build-essential
- cmake
- zlib1g-dev
- libboost-dev
- libx11-dev  
- libboost-filesystem1.67-dev


One option to install them all at once is to run the following command:

- $sudo apt install build-essential cmake zlib1g-dev libboost-dev libx11-dev libboost-filesystem1.67-dev

## 2) Compile Fate Classifier:

#### Centralized Fate Classifier is used to reconstruct the gene tree of the aevol populations and generation, also to calculate the probabilities of post-replication fates. To compile Centralized Fate Classifier:

- $cd CentralizedFateClassifier
- $g++ -std=c++11 CentralizedFateClassifier.cpp GeneTreeConstructor.h probabilitiescalculation.h -o CentralizedFateClassifier

## 3) Launch Simulation :

#### To run a simulation in aevol we need to do following steps:

 ##### A) We need to prepare a param.in file that contains aevol simulation parameters (such as mutation rates, BACKUP_STEP, Target function/fitness funcations/environment). A list of parameter files according to environments are available in "Simulation/param files/" (e.g. param1.in represents environment 1). Also, we might need a preevolved genome as initial genome. A list of initial genomes (wildtypes) are available in "Simulation/wildtypes/" (e.g. WT_W0.1_V1.txt represents a wild type with a maximum width of 0.1, and V stands for version). To create a simulation, place a param file, a wildtype and CentralizedFateClassifier file in CentralizedFateClassifier directory to a new directory "Edited Aevol/examples/new simulation". For example:
 
- $cd Edited Aevol/examples/new simulation
- $ls
- output: param1.in WT_W0.1_V1.txt CentralizedFateClassifier
 
  
 ##### B) Second, run following command to run a simulation: 
 
- $../../build/bin/aevol_create -C WT_W0.1_V1.txt
- $../../build/bin/aevol_run -n 1100000
- $../../build/bin/aevol_post_lineage
- $../../build/bin/aevol_post_protein_map lin*
- $./CentralizedFateClassifier -m default proteins_list_after_events.csv -e 1000000

An alternative option that we often use (and that proved efficient) is to generate Wild-Types on long evolutions (typically 10 million generations) and then to initialize all the experiments from these WTs. So, instead of using ../../build/bin/aevol_create we can use ../../build/bin/aevol_create -C WT0.txt (a list of wild types that evolved for 10 million generations are in CentralizedFateClassifier directory)
  
  3)Third, run this command to run the simulation by specifying the number of generation using -n: ../../build/bin/aevol_run -n 100000 (note: just make sure the number of generations here should be greater or equal to the BACKUP_STEP in the param.in) 
    After using this command, we have list of outputs that we can use to find the list of events during generations for the best individual.
    
  4)To reconstruct the lineage of a given individual we used: ../../build/bin/aevol_post_lineage (It rebuilds it for the best individual without specifying a individual)
  
  5)Finally, using this command, we can map all proteins after each event for the given lineage file generated in the previous step: ../../build/bin/aevol_post_protein_map lineage-b000000000-e000010000-i0-r0.ae
    The output is a CSV file containing a list of proteins mapped after each event.

But at the bash script, we want to run serveral simulations with different fitness function. So, we use all these command to run simulations in bash script but with using different param.in file at each simulation. 

we can use bash script located in the CentralizedFateClassifier directory to run a bunch of aevol simulations.

C) Eventually, we have a CSV file of mapped proteins, with using following command we can build gene tree and probabilities of each duplication fate for each duplication event that happend in the simulation.
   CentralizedFateClassifier -m default csvfilename (csvfilename is the CSV file generated by the previous command)
   The outputs are three files: One newick gene Tree and Two CSV file which one of them is list of duplication fates probabilities of a duplication node and all combination of left and right descendts and another CSV file is avrage of duplication fates probabilities for each duplication event.

