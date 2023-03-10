// ============================================================================
//                                   Includes
// ============================================================================
#include <cerrno>
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <zlib.h>
#include <sys/stat.h>
#include <getopt.h>
#include <assert.h>

#include "aevol.h"

using namespace aevol;

int main(int argc, char** argv) {
  int64_t luca_time;
  int64_t last_gen;
  int64_t t0 = 0;

  printf("\n\n");
  printf("====================================\n");
  printf(" Loading the info of last save ... \n");

  std::ifstream light_tree_file;
  light_tree_file.open(LIGHTTREE_FILE_NAME);
  if(light_tree_file and light_tree_file.peek() != EOF) {
    light_tree_file >> luca_time;
  }
  else {
    exit(EXIT_FAILURE);
  }

  last_gen = OutputManager::last_gener();

  printf("OK\n");
  printf("====================================\n");

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // Check that the tree was recorded
  if (not exp_manager->record_light_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                          "evolution, could not reconstruct the lineage");
  }
  int64_t backup_step = exp_manager->backup_step();

  //the table will be one index bigger for sementic. the index 0 will never be used
  ReplicationReport** reports = new ReplicationReport*[last_gen+1];
  reports[0] = nullptr;

  printf("\n\n");
  printf("====================================\n");
  printf(" Loading the branch tree file ... ");
  char input_file_name[255];

  std::map<int64_t, std::map<int32_t, ReplicationReport*>> replics;

  sprintf(input_file_name, "lightTree/tree_branches.ae");
  gzFile branches_file = gzopen( input_file_name, "r" );
  if (branches_file == Z_NULL) {
    printf("ERROR : Could not read tree file %s\n", input_file_name);
    exit(EXIT_FAILURE);
  }

  int c = gzgetc(branches_file);
  while(c > -1) {
    gzungetc(c, branches_file);
    ReplicationReport* rep = new ReplicationReport(branches_file, nullptr);
    int64_t generation;
    gzread(branches_file, &generation, sizeof(generation));
    //std::cout << "Load report from id : " << rep->id() << " generation : " << generation << '\n';
    replics[generation][rep->id()] = rep;
    c = gzgetc(branches_file);
  }
  gzclose(branches_file);
  printf("OK\n");
  printf("====================================\n");

  int32_t best_id = 0;
  int32_t best_rank = 0;

  // Specify the index of indiv
  int final_indiv_index = 0;

  for(auto rep : replics[last_gen]) {
    if (rep.first == final_indiv_index) {
      // best_rank = rep.second->rank();
      best_id = rep.first;
    }
  }

  //Find the best indiv

  /*for(auto rep : replics[last_gen]) {
    if(rep.second->rank() > best_rank) {
      best_rank = rep.second->rank();
      best_id = rep.first;
    }
  }*/
  //His rank have to be the number of indivs
  std::cout << "Best indiv is : " << best_id << " with rank : " << best_rank << '\n';
  int32_t anc_id = best_id;

  //Fill the end of all the replic with the lineage of the best
  for(int64_t gen = last_gen ; gen >= luca_time and gen > 0 ; gen--) {
    reports[gen] = replics[gen][anc_id];
    anc_id = reports[gen]->parent_id();
  }

  printf("\n\n");
  printf("====================================\n");
  printf(" Loading the trunc tree files ... \n");

  //Fill the rest of the lineage
  int64_t current_gen = backup_step;

  sprintf(input_file_name, "lightTree/tree_trunc" TIMESTEP_FORMAT ".ae", current_gen);
  gzFile trunc_file = gzopen( input_file_name, "r" );
  if (trunc_file == Z_NULL) {
    printf("ERROR : Could not read tree file %s\n", input_file_name);
    exit(EXIT_FAILURE);
  }

  //printf("Luca is %d\n",luca_time);
  for(int64_t gen = 1 ; gen < luca_time ; ) {
      //printf("Searching at gen %d\n",gen);
    c = gzgetc(trunc_file);
    if(c < 0) {
      gzclose(trunc_file);
      current_gen += backup_step;
      if(current_gen > last_gen) {
          //printf("Current gen %d -- %d\n",current_gen,last_gen);
        break;
      }
      sprintf(input_file_name, "lightTree/tree_trunc" TIMESTEP_FORMAT ".ae", current_gen);
      trunc_file = gzopen( input_file_name, "r" );
      if (trunc_file == Z_NULL) {
        printf("ERROR : Could not read tree file %s\n", input_file_name);
        exit(EXIT_FAILURE);
      }
      continue;
    }
    gzungetc(c, trunc_file);
    ReplicationReport* rep = new ReplicationReport(trunc_file, nullptr);
    reports[gen] = rep;
      //printf("Getting TRUNC the replication report for the ancestor at generation %" PRId64 " => %d (%d)\n", gen,rep->id(),rep->parent_id());
    gen++;
  }
  printf("OK\n");
  printf("====================================\n");

  printf("\n\n\n");
  printf("============================================================ \n");
  printf(" Write the replication reports in the output file... \n");

  // =======================
  //  Open the output file
  // =======================
  char output_file_name[255];

  sprintf(output_file_name, "soft_lineage" TIMESTEP_FORMAT ".ae", last_gen);

  gzFile soft_lineage_file = gzopen(output_file_name, "w");
  if (soft_lineage_file == nullptr) {
    fprintf(stderr, "File %s could not be created.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }

    printf("Loading individual %lld\n",reports[1]->parent_id());

  const Individual& initial_ancestor = *(exp_manager->indiv_by_id(reports[1]->parent_id()));

  // Write file "header"
  gzwrite(soft_lineage_file, &t0, sizeof(t0));
  gzwrite(soft_lineage_file, &last_gen, sizeof(last_gen));
  gzwrite(soft_lineage_file, &best_id, sizeof(best_id));
  gzwrite(soft_lineage_file, &best_rank, sizeof(best_rank));

  initial_ancestor.grid_cell()->save(soft_lineage_file);

  std::ofstream of;
  of.open("lignée.txt");

  for(int64_t i=1 ; i <= last_gen ; i++) {
    of << i << " : " << reports[i]->parent_id() << "  " << reports[i]->id() << std::endl;
    reports[i]->write_to_tree_file(soft_lineage_file);
  }
  of.close();
  gzclose(soft_lineage_file);
  printf("OK\n");
  printf("====================================\n");

  delete [] reports;
  delete exp_manager;

}
