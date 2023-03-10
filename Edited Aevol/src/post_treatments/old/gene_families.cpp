// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************




// =================================================================
//                              Libraries
// =================================================================
#include <list>
#include <algorithm>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>  // for the permission symbols used with mkdir

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// TODO : predire la position de tous les genes apres la mutation
// puis verifier dans la liste des genes de l'indiv apres evaluation
// s'il y a bien toujours un gene a cet endroit


// TODO : transfers
// TODO : in gene tree file, print gene sequence, final position, list of proteins and M, W, H


enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

void update_pointers_in_trees(std::list<ae_gene_tree>& gene_trees, GeneticUnit * unit);
void anticipate_mutation_effect_on_genes_in_trees(std::list<ae_gene_tree>& gene_trees, const ae_mutation * mut, int32_t unitlen_before);
void register_actual_mutation_effect_on_genes_in_trees(std::list<ae_gene_tree>& gene_trees, const ae_mutation* mut, GeneticUnit * unit, int32_t gener, double impact_on_metabolic_error);
void search_protein_in_gene_trees(std::list<ae_gene_tree>& gene_trees, ae_protein * prot, ae_gene_tree ** resultTree, bool& found);
void set_end_gener_if_active_leaves_in_trees(std::list<ae_gene_tree>& gene_trees, int32_t gener);
void write_gene_trees_to_files(std::list<ae_gene_tree>& gene_trees, int32_t end_gener);
void print_gene_trees_to_screen(std::list<ae_gene_tree>& gene_trees); // For debug purposes






int main(int argc, char** argv)
{
  // The input file (lineage.ae or lineage.rae) must contain the following information:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int32_t)
  // - end gener                                                  (int32_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (ae_individual::write_to_backup)
  // - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)
  // - ...
  // - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)




  // =====================
  //  Parse command line
  // =====================

  // Default values
  check_type  check               = LIGHT_CHECK;
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  double      tolerance           = 0;


  const char * short_options = "hVvncf:t:";
  static struct option long_options[] =
  {
    {"help",      no_argument,       NULL, 'h'},
    {"version",   no_argument,       NULL, 'V' },
    {"verbose",   no_argument,       NULL, 'v'},
    {"nocheck",   no_argument,       NULL, 'n'},
    {"fullcheck", no_argument,       NULL, 'c'},
    {"file",      required_argument, NULL, 'f'},
    {"tolerance",   required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'n' : check = NO_CHECK;                  break;
      case 'c' : check = FULL_CHECK;                break;
      case 'f' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -f or --file : missing argument.\n");
          exit(EXIT_FAILURE);
        }

        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf(lineage_file_name, "%s", optarg);
        break;
      }
      case 't' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -t or --tolerance : missing argument.\n");
          exit(EXIT_FAILURE);
        }
        check = ENV_CHECK;
        tolerance = atof(optarg);
        break;
      }
      default :
      {
        fprintf(stderr, "ERROR : Unknown option, check your syntax.\n");
        print_help(argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }

  if (lineage_file_name == NULL)
  {
    fprintf(stderr, "ERROR : Option -f or --file missing. \n");
    exit(EXIT_FAILURE);
  }


  printf("\n");
  printf("WARNING : Parameter change during simulation is not managed in general.\n");
  printf("          Only changes in environmental target done with aevol_modify are handled.\n");
  printf("\n");


  // =======================
  //  Open the lineage file
  // =======================

  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL)
  {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
    exit(EXIT_FAILURE);
  }

  int64_t t0;
  int64_t t_end;
  int32_t  final_index, final_indiv_rank;
  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_index, sizeof(final_index));
  gzread(lineage_file, &final_indiv_rank,   sizeof(final_indiv_rank));

  if (verbose)
  {
    printf("\n\n");
    printf("================================================================================\n");
    printf(" Gene families of the ancestors of indiv. #%" PRId32 " (t=%" PRId64 " to %" PRId64 ")\n",
            final_index, t0, t_end);
    printf("================================================================================\n");
  }




  // =========================
  //  Open the experience manager
  // =========================


  ae_exp_manager* exp_manager = new ae_exp_manager();
  exp_manager->load(t0, true, false);
  Environment* env = new Environment(*(exp_manager->env())); // independent copy

  int64_t backup_step = exp_manager->backup_step();


  // ==============================
  //  Prepare the initial ancestor
  // ==============================

  ae_individual * indiv = new ae_individual(exp_manager, lineage_file);
  indiv->evaluate(env);
  indiv->compute_statistical_data();

  if (verbose)
  {
    printf("Initial fitness     = %e\n", indiv->fitness());
    printf("Initial genome size = %" PRId32 "\n", indiv->total_genome_size());
  }


  // ===========================================================================
  //  Prepare the initial gene trees (one for each coding RNA in the ancestor)
  // ===========================================================================

  // Each initial gene in this ancestral genome will be the root of a gene tree,
  // where the paralogs (gene copies created by duplication) will be monitored

  std::list<ae_gene_tree> gene_trees;

  for (const auto& unit: indiv->genetic_unit_list())
  {
    for (const auto& prot: unit.protein_list(LEADING))
      gene_trees.emplace_back(t0, prot);
    for (const auto& prot: unit.protein_list(LAGGING))
      gene_trees.emplace_back(t0, prot);
  }

  // ===============================================================================
  //  Replay the mutation to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir and genome are correct each
  //   time a backup is available)
  // ===============================================================================
  ae_replication_report* rep = NULL;

  std::list<GeneticUnit>::const_iterator unit;

  ae_individual* stored_indiv = NULL;
  std::list<GeneticUnit>::const_iterator stored_unit;

  int32_t index, genetic_unit_number, unitlen_before;
  double metabolic_error_before, metabolic_error_after, impact_on_metabolic_error;


  ae_gene_tree * genetree = NULL;
  bool found = false;

  ae_exp_manager* exp_manager_backup = NULL;
  Environment* backup_env = NULL;

  bool check_now = false;


  aevol::Time::plusplus();
  while (time() <= t_end)
  {
    env->build();

    rep = new ae_replication_report(lineage_file, indiv);
    index = rep->id(); // who are we building...
    indiv->set_replication_report(rep);

    // Check now?
    check_now = ((check == FULL_CHECK && ae_utils::mod(time(), backup_step) == 0) ||
                 (check == ENV_CHECK && ae_utils::mod(time(), backup_step) == 0) ||
                 (check == LIGHT_CHECK && time() == t_end));


    if (verbose)
      printf("Rebuilding ancestor at time %" PRId64 " (index %" PRId32 ")...",
          time(), index);

    env->build();
    env->apply_variation();
    indiv->reevaluate(env);

    // Check, and possibly update, the environment according to the backup files (update necessary if the env. was modified by aevol_modify at some point)
    if (ae_utils::mod(time(), backup_step) == 0)
      {
        char env_file_name[255];
        sprintf(env_file_name, "./" ENV_FNAME_FORMAT, time());
        gzFile env_file = gzopen(env_file_name, "r");
        backup_env = new Environment();
        backup_env->load(env_file);

        if (! env->is_identical_to(*backup_env, tolerance))
          {
            printf("Warning: At t=%" PRId64 ", the replayed environment is not the same\n", time());
            printf("         as the one saved at t=%" PRId64 "... \n", time());
            printf("         with tolerance of %lg\n", tolerance);
            printf("Replacing the replayed environment by the one stored in the backup.\n");
            delete env;
            env = new Environment(*backup_env);
          }
        delete backup_env;
      }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution, or if some translocations occurred between different genetic units

    genetic_unit_number = 0;
    auto unit = indiv->genetic_unit_list_nonconst().begin();


    if (check_now)
    {
      exp_manager_backup = new ae_exp_manager();
      exp_manager_backup->load(time(), true, false);
      // TODO: disabled tmp
      // stored_indiv = new ae_individual(* (ae_individual *)exp_manager_backup->indiv_by_id(index), false);
      stored_unit = stored_indiv->genetic_unit_list().begin();
    }


    for (const auto& dnarep: rep->dna_replic_reports()) {
      assert(unit != indiv->genetic_unit_list().cend());
      unit->dna()->set_replic_report(dnarep);
      update_pointers_in_trees(gene_trees, &*unit); // because of the reevaluate at each new generation (envir. variation possible)

      // ***************************************
      //             Transfer events
      // ***************************************

      // TO DO
     /*
      mnode = dnarep->HT()->first();
      while (mnode != NULL)
      {
        mut = (ae_mutation *) mnode->obj();

        metabolic_error_before = indiv->dist_to_target_by_feature(METABOLISM);
        unitlen_before = unit.dna()->length();
        unit.compute_nb_of_affected_genes(mut, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);


        unit.dna()->undergo_this_mutation(mut);
        indiv->reevaluate(env);


        metabolic_error_after = indiv->dist_to_target_by_feature(METABOLISM);
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;


        mut->generic_description_string(mut_descr_string);
        fprintf(output, "%"PRId64 " %"PRId32 " %s %"PRId32 " %.15f  %"PRId32 " %"PRId32 " %"PRId32 " \n",\
                 time(), genetic_unit_number, \
                 mut_descr_string, unitlen_before, \
                 impact_on_metabolic_error, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);



        mnode = mnode->next();
      } */


      // ***************************************
      //           Rearrangement events
      // ***************************************

      for (const auto& mut: dnarep->rearrangements()) {
        metabolic_error_before = indiv->dist_to_target_by_feature(METABOLISM);
        unitlen_before = unit->dna()->length();
        anticipate_mutation_effect_on_genes_in_trees(gene_trees, &mut, unitlen_before);

        unit->dna()->undergo_this_mutation(&mut);

        indiv->reevaluate(env);
        metabolic_error_after = indiv->dist_to_target_by_feature(METABOLISM);
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        register_actual_mutation_effect_on_genes_in_trees(gene_trees, &mut, &*unit, time(), impact_on_metabolic_error);

        /* New genes that have been created "from scratch", i.e. not by duplication => new gene tree */
        for (const auto& prot: unit->protein_list(LEADING)) {
          search_protein_in_gene_trees(gene_trees, prot, &genetree, found);
          if (not found)
            gene_trees.emplace_back(time(), prot, &mut);
        }
        for (const auto& prot: unit->protein_list(LAGGING)) {
          search_protein_in_gene_trees(gene_trees, prot, &genetree, found);
          if (not found)
            gene_trees.emplace_back(time(), prot, &mut);
        }
         // print_gene_trees_to_screen(gene_trees);// DEBUG
         // indiv->print_protein_list(); // DEBUG
      }


      // ***************************************
      // Local events (point mutations & small indels)
      // ***************************************

      for (const auto& mut: dnarep->mutations()) {
        metabolic_error_before = indiv->dist_to_target_by_feature(METABOLISM);
        unitlen_before = unit->dna()->length();
        anticipate_mutation_effect_on_genes_in_trees(gene_trees, &mut, unitlen_before);
        unit->dna()->undergo_this_mutation(&mut);

        indiv->reevaluate(env);
        metabolic_error_after = indiv->dist_to_target_by_feature(METABOLISM);
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        register_actual_mutation_effect_on_genes_in_trees(gene_trees, &mut, &*unit, time(), impact_on_metabolic_error);

        /* New genes that have been created "from scratch", i.e. not by duplication => new gene tree */
        for (const auto& prot: unit->protein_list(LEADING)) {
          search_protein_in_gene_trees(gene_trees, prot, &genetree, found);
          if (not found)
            gene_trees.emplace_back(time(), prot, &mut);
        }
        for (const auto& prot: unit->protein_list(LAGGING)) {
          search_protein_in_gene_trees(gene_trees, prot, &genetree, found);
          if (not found)
            gene_trees.emplace_back(time(), prot, &mut);
        }
      }

      if (check_now && ae_utils::mod(time(), backup_step) == 0)
      {
        if (verbose)
        {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }

        assert(stored_unit != stored_indiv->genetic_unit_list().end());

        char * str1 = new char[unit->dna()->length() + 1];
        memcpy(str1, unit->dna()->data(), \
               unit->dna()->length()*sizeof(char));
        str1[unit->dna()->length()] = '\0';

        char * str2 = new char[(stored_unit->dna())->length() + 1];
        memcpy(str2, (stored_unit->dna())->data(), (stored_unit->dna())->length()*sizeof(char));
        str2[(stored_unit->dna())->length()] = '\0';

        if(strncmp(str1,str2, (stored_unit->dna())->length())==0)
        {
          if (verbose) printf(" OK\n");
        }
        else
        {
          if (verbose) printf(" ERROR !\n");
          fprintf(stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf(stderr, "the one saved at generation %" PRId64 "... ", t0);
          fprintf(stderr, "Rebuilt unit : %zu bp\n %s\n", strlen(str1), str1);
          fprintf(stderr, "Stored unit  : %zu bp\n %s\n", strlen(str2), str2);
          delete [] str1;
          delete [] str2;
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete exp_manager_backup;
          delete exp_manager;
          exit(EXIT_FAILURE);
        }

        delete [] str1;
        delete [] str2;

        ++stored_unit;
      }
      ++unit;
      genetic_unit_number ++;
    }

    assert(unit ==  indiv->genetic_unit_list().cend());


    if (verbose) printf(" OK\n");

    delete rep;

    if (check_now && ae_utils::mod(time(), backup_step) == 0)
    {
      assert(stored_unit == stored_indiv->genetic_unit_list().end());
      delete stored_indiv;
      delete exp_manager_backup;
    }

    aevol::Time::plusplus();
  }

  set_end_gener_if_active_leaves_in_trees(gene_trees, t_end);
  write_gene_trees_to_files(gene_trees, t_end);

  gzclose(lineage_file);
  delete exp_manager;
  delete indiv;
  delete env;

  exit(EXIT_SUCCESS);

}






void search_protein_in_gene_trees(std::list<ae_gene_tree>& gene_trees, ae_protein * prot, ae_gene_tree ** resultTree, bool& found)
{
  fflush(stdout);
  auto it = find_if(gene_trees.begin(),
                    gene_trees.end(),
                    [prot](ae_gene_tree& m)
                    { return m.search_in_leaves(prot) == NULL; });
  *resultTree = &*it;
  found = (it != gene_trees.end());
}



void update_pointers_in_trees(std::list<ae_gene_tree>& gene_trees, GeneticUnit * unit)
{
  for (auto& n: gene_trees)
    n.update_pointers_in_tree_leaves(unit);
}


void anticipate_mutation_effect_on_genes_in_trees(
    std::list<ae_gene_tree>& gene_trees,
    const ae_mutation* mut,
    int32_t unitlen_before)
{
  for (auto& n: gene_trees)
    n.anticipate_mutation_effect_on_genes_in_tree_leaves(mut, unitlen_before);
}

void register_actual_mutation_effect_on_genes_in_trees(
    std::list<ae_gene_tree>& gene_trees,
    const ae_mutation* mut,
    GeneticUnit* unit,
    int32_t gener,
    double impact_on_metabolic_error)
{
  for (auto& n: gene_trees)
    n.register_actual_mutation_effect_on_genes_in_tree_leaves(mut, unit, gener, impact_on_metabolic_error);
}



void set_end_gener_if_active_leaves_in_trees(std::list<ae_gene_tree>& gene_trees, int32_t gener)
{
  for (auto& n: gene_trees)
    n.set_end_gener_if_active_leaves(gener);
}


void print_gene_trees_to_screen(std::list<ae_gene_tree>& gene_trees)
{
  int32_t tree_number = 0;
  for (auto& n: gene_trees) {
    printf("Content of tree %d :\n", tree_number);
    n.print_to_screen();
    tree_number++;
  }
}


void write_gene_trees_to_files(std::list<ae_gene_tree>& gene_trees, int32_t end_gener)
{
  // Prepare the directory for the outputs files related to the gene trees
  char directory_name[] = "gene_trees";

  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0)
    {
      //       struct stat status;
      //       stat(directory_name, &status);
      //       if (status.st_mode & S_IFDIR) cout << "The directory exists." << endl;
      //       else cout << "This path is a file." << endl;

      if (access(directory_name, X_OK | W_OK) != 0)
        {
          fprintf(stderr, "Error: cannot enter or write in directory %s.\n", directory_name);
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      // Create the directory with permissions : rwx r-x r-x
      if (mkdir(directory_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
        {
          fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
          exit(EXIT_FAILURE);
        }
    }


  FILE * tree_statistics_file = fopen("gene_trees/gene_tree_statistics.txt", "w");
  if (tree_statistics_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file gene_trees/gene_tree_statistics.txt.\n");
      exit(EXIT_FAILURE);
    }
  fprintf(tree_statistics_file, "treeID creationType beginGener endGener nbNodes nbInternalNodes nbLeaves nbActiveLeaves\n");


  FILE * nodeattr_tabular_file = fopen("gene_trees/nodeattr_tabular.txt", "w");
  if (nodeattr_tabular_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file gene_trees/nodeattr_tabular.txt.\n");
      exit(EXIT_FAILURE);
    }
  fprintf(nodeattr_tabular_file, "treeID nodeID parentID leftchildID rightchildID nodeCreationDate dnaCreationDate nodeStatus geneLossOrDupDate strand shineDalPos nbProm meanTriangle widthTriangle heightTriangle concentration nbMutTot nbLocalmutUpstreamNeutral nbLocalmutUpstreamBenef nbLocalmutUpstreamDelet nbLocalmutCdsNeutral nbLocalmutCdsBenef nbLocalmutCdsDelet nbRearUpstreamNeutral nbRearUpstreamBenef nbRearUpstreamDelet nbRearCdsNeutral nbRearCdsBenef nbRearCdsDelet\n");

  char topol_file_name[128];
  char node_attr_file_name[128];
  int32_t tree_number = 0;

  for (auto& n: gene_trees) {
      fprintf(tree_statistics_file, "%" PRId32 " ", tree_number);
      if (n.creation_type() == INITIALIZATION)  fprintf(tree_statistics_file, "INITIALIZATION ");
      else if (n.creation_type() == LOCAL_MUTATION)  fprintf(tree_statistics_file, "LOCAL_MUTATION ");
      else if (n.creation_type() == REARRANGEMENT)  fprintf(tree_statistics_file, "REARRANGEMENT ");
      else fprintf(tree_statistics_file, "TRANSFER ");
      fprintf(tree_statistics_file, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "\n", n.begin_gener(), n.end_gener(), n.total_nb_nodes(), n.nb_internal_nodes(), n.nb_leaves(), n.nb_active_leaves());
      sprintf(topol_file_name, "gene_trees/genetree" TIMESTEP_FORMAT "-topology.tre", tree_number);
      sprintf(node_attr_file_name, "gene_trees/genetree" TIMESTEP_FORMAT "-nodeattr.txt", tree_number);
      n.write_to_files(topol_file_name, node_attr_file_name, end_gener);
      n.write_nodes_in_tabular_file(tree_number, nodeattr_tabular_file);
      tree_number ++;
    }

  fclose(tree_statistics_file);
  fclose(nodeattr_tabular_file);
}




/*!
  \brief

*/
void print_help(char* prog_path)
{
  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                 Gene families post-treatment program                 * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rgene_families -h\n");
  printf("or :    rgene_families [-vn] -f lineage_file \n");
#else
  printf("Usage : gene_families -h\n");
  printf("or :    gene_families [-vn] -f lineage_file \n");
#endif
  printf("\n");
  printf("This program traces the evolution of gene families on a lineage.\n");
  printf("of lineage_file\n");
  printf("\n");
  printf("WARNING: This program should not be used for simulations run with several genetic units\n");
  printf("(eg chromosome + plasmids). It works only for simulations with only a single chromosome.\n");
  printf("\n");
  printf("\t-h or --help       : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose    : Be verbose, listing generations as they are \n");
  printf("\t                       treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome and environment checks every\n");
  printf("\t                       <BACKUP_STEP> generations. Default behaviour is\n");
  printf("\t                       lighter as it only perform sthese checks at the\n");
  printf("\t                       ending generation.\n");
  printf("\n");
  printf("\t-f lineage_file or --file lineage_file : \n");
  printf("\t                       Compute the fixed mutations of the individuals within lineage_file.\n");
  printf("\n");
  printf("\t-t tolerance or --tolerance tolerance : \n");
  printf("\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf("\n");
}
