###############################################################################
#                         The Workflow example
###############################################################################
# This is the workflow example. It provides an example of one of the many
# different workflows that can be used for experiments with Aevol.
#
# The main idea underlying this workflow is to parallel wet lab experiments,
# which are conducted on evolved organisms.
# To use already evolved organisms for Aevol experiments, one can either use
# an evolved genome provided by the community or evolve his own.
# This example does the latter (more complete) case.
#
#
###################
# RECOMMENDATIONS #
###################
#
# It is strongly recommended you read the corresponding section in the user
# manual while following this example since it contains an explanation for
# each step of the example.
#
# The following set of commands should work as is if you have installed aevol
# with make install, which is recommended.
#
# If you haven't installed aevol, you will need to specify where to find the
# executables, which should be in <aevol_dir>/src for the 4 main executables
# and in <aevol_dir>/src/post_treatments for the rest, <aevol_dir> being the
# main aevol directory you have downloaded.
# E.g. if your <aevol_dir> is /home/login/aevol, then the command
# aevol_run -n 5000 will become /home/login/aevol/src/aevol_run -n 5000
# and ae_misc_lineage -e 10000 will become
# /home/login/aevol/src/post_treatments/ae_misc_lineage -e 10000
#
###############################################################################


# ========== Wild-Type generation ==========
#

#mkdir wild_type 
cd wild_type
aevol_create
aevol_run -n 5000   




# ========== Experimental setup ==========
#

cd ..
# Propagate the experiment, meaning prepare directories for different
# runs starting from the wild type 
mydirnamesA="line01 line02 line03 line04 line05"
mydirnamesB="line06 line07 line08 line09 line10"
for mydir in $mydirnamesA $mydirnamesB
do
   echo $mydir
   aevol_propagate -g 5000 -i wild_type -o $mydir -S $RANDOM
done  

# For each experiment, create a file with the parameters to change
echo "# New environment 
    ENV_GAUSSIAN  0.5   0.2   0.05
    ENV_GAUSSIAN  0.5   0.4   0.05
    ENV_GAUSSIAN  0.5   0.8   0.05
    ENV_VARIATION none" > newparam-groupA.in

echo "# New environment
    ENV_GAUSSIAN  0.5   0.2   0.05
    ENV_GAUSSIAN  0.5   0.4   0.05
    ENV_GAUSSIAN  0.5   0.8   0.05
    ENV_VARIATION none
   # New rearrangement rates
    DUPLICATION_RATE          1e-5
    DELETION_RATE                1e-5
    TRANSLOCATION_RATE    1e-5
    INVERSION_RATE              1e-5" > newparam-groupB.in


# Apply these modifications
#
for mydir in $mydirnamesA 
do
   cd $mydir
   echo $mydir
   aevol_modify  --gener 0 --file ../newparam-groupA.in
   cd ..
done 

for mydir in $mydirnamesB 
do
   cd $mydir
   echo $mydir
   aevol_modify  --gener 0 --file ../newparam-groupB.in
   cd ..
done 




# ========== Run the simulations ==========
#
for mydir in $mydirnamesA $mydirnamesB
do
   cd mydir
   aevol_run -n 20000
   cd ..
done 



# ========== Analyse the outcome ==========
#
# A set of post-treatment tools is available to help analyse the outcome.
# 

# ---------- aevol_misc_view_generation ----------
#
# The simplest miscellaneous tool is view_generation. It allows one to
# visualize a generation using the exact same graphical outputs used in
# aevol_run.
# However, since it relies on graphics, it is only available when aevol is
# compiled with x enabled (which is the default).
#
aevol_misc_view_generation -g 10000


# ---------- aevol_misc_create_eps ----------
#
# Similarly, one can obtain eps outputs for a given generation with the 
# create_eps tool. Files will be outputted in eps_files_xxxxxx (with 
# xxxxxx the generation number)
#
aevol_misc_create_eps -g 10000


# ---------- aevol_misc_robustness ----------
#
# The robustness tool computes the replication statistics of all the 
# individuals of a given generation, like the proportion of neutral, beneficial, 
# deleterious offsprings. This is done by simulating nbchildren replications 
# for each individual (1000 replications by default), with its mutation, 
# rearrangement and transfer rates. Depending on those rates and genome 
# size, there can be several mutations per replication. Those global statistics 
# are written in stat/robustness_numgener.out, with one line per individual
# in the specified generation.
# The program also outputs detailed statistics for one of the individuals (the 
# best one by default). The detailed statistics for this individual are written in
# stats/replication_numgener.out, with one line per simulated child of this 
# particular individual.
#
aevol_misc_robustness -g 10000 



# ---------- aevol_misc_lineage ----------
#
# One can reconstruct the lineage of an evolved individual.
# This will generate a lineage file whose name will look like
# lineage-b000000-e010000-i999-r1000.ae containing the complete mutational
# history of a given individual of a given generation.
# This file can then be used as the input for subsequent post-treatments.
#
aevol_misc_lineage -e 10000


# ---------- aevol_misc_ancstats ----------
#
# Statistics of a lineage can be obtained with this tool.
# The generated stats are outputted in stats/ancstats/
#
aevol_misc_ancstats -f lineage-b000000-e010000-i*-r1000.ae




# ---------- aevol_misc_fixed_mutations ----------
#
# This tool outputs the list of mutational events that occurred
# on the lineage given as input.
# The generated list is outputted in stats/
#
aevol_misc_fixed_mutations -f lineage-b000000-e010000-i*-r1000.ae


# ---------- aevol_misc_gene_families ----------
#
# This tool outputs the history of each gene family in the lineage 
# given as input. 
# The generated gene trees are outputted in gene_trees/
# This analysis can be longer than ancstats or fixed_mutations,
# it can take from a few minutes to a few hours depending on
# gene number evolution in the lineage.

aevol_misc_gene_families -f lineage-b000000-e010000-i*-r1000.ae
