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

## 3) Launch a Simulation:

#### To run a simulation in aevol we need to do following steps:

 ##### A) We need to prepare a param.in file that contains aevol simulation parameters (such as mutation rates, BACKUP_STEP, Target function/fitness funcations/environment). A list of parameter files according to environments are available in "Simulation/param files/" (e.g. param1.in represents environment 1). Also, we might need a preevolved genome as initial genome. A list of initial genomes (wildtypes) are available in "Simulation/wildtypes/" (e.g. WT_W0.1_V1.txt represents a wild type with a maximum width of 0.1, and V stands for version). To create a simulation, place a param file, a wildtype and CentralizedFateClassifier file in CentralizedFateClassifier directory to a new directory "Edited Aevol/examples/new simulation". For example:
 
- $cd Edited Aevol/examples/new simulation
- $ls
- output: param1.in WT_W0.1_V1.txt CentralizedFateClassifier
 
  
 ##### B) Second, run following command to run a simulation: 
 
- $../../build/bin/aevol_create -C WT_W0.1_V1.txt
- $../../build/bin/aevol_run -n 1100000      (-n 1100000 specified the number of generations)
- $../../build/bin/aevol_post_lineage   (To reconstruct the lineage of an individual)
- $../../build/bin/aevol_post_protein_map lin*  (To map all proteins after each event for the given lineage file generated from previous command)
- $./CentralizedFateClassifier -m default proteins_list_after_events.csv -e 1000000  (To reconstruct gene tree and calculate the probabilities of post-replication fates of 1000000 first generations from proteins_list_after_events.csv that is generated from previous command)

###### Outputs: dups_fates_probablities.csv which contains the probabilities of post duplication fates of ancestors and two it's descendants. Also, newick.txt represents the reconstructed gene tree.

## 4) Launch Bunch of Simulations:

#### There is a prepared script "Simulation/running script.sh" (place it at new simulation directory) that by specify the wildtype and param file in the script, we are able to launch many simulation (by changing seed value in param file) for specfic wildtype and param file. To run script:

- $bash ./running script.sh

###### Outputs: The output is a directory "SimWT_W0.1_V1_param1_CSVs" which contains all csv files and newick files.

## 5) Concat All CSV files:

#### To concat all CSV files we need to use python code in "Simulation/Concat simulation files" as follows:

- $run-maker.py --indir=Aevol/examples/new simulation/SimWT_W0.1_V1_param1_CSVs --outdir=Aevol/examples/new simulation/SimWT_W0.1_V1_param1_CSVs_v45 --dov45

###### Outputs: The output is a CSV file "Aevol/examples/new simulation/SimWT_W0.1_V1_param1_CSVs_v45/SimWT_W0.1_V1_param1_CSVs_v45.csv" which contains all csv files.

## 6) Analysis:

#### Finally by using scripts in "Simulation/Analysis scripts" we can analysis results SimWT_W0.1_V1_param1_CSVs_v45.csv.
