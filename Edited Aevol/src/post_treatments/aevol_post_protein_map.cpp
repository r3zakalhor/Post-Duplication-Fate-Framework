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

// The input file is produced by the lineage post-treatment, please refer to it
// for e.g. the file format/content

// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <zlib.h>
#include <err.h>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>

#include <iomanip>
#include <limits>

#include "aevol.h"
#include "Protein_list.h"
#include "../macros.h"
//#include "GeneTree.h"


#include "../../../CentralizedFateClassifier/probabilitiescalculation.h"


bool SWITCH_DEBUG_MODE = false;
bool FATE_PROB_MODE = false;
int SWITCH_DEBUG_GENERATION = 1;
bool SWITCH_DEBUG_SHOW_SEQUENCES = false;


using namespace aevol;
using namespace std;

class ProteinMap {
public:
  ProteinMap() {
    
  }
  ProteinMap(int nb_prot, int dna_length) {
    nb_prot_ = nb_prot;
    start_pos_.resize(nb_prot_);
    length_.resize(nb_prot_);
    basal_level_.resize(nb_prot_);
    hamming_dist_.resize(nb_prot_);
    dist_next_prot_.resize(nb_prot_);
    dna_length_ = dna_length;
  }
  
  void add_protein(int32_t start_pos, int32_t length, int32_t basal_level, int32_t hamming_dist, int32_t dist_next_prot) {
    
    start_pos_[cpt_] = start_pos;
    length_[cpt_] = length;
    basal_level_[cpt_] = basal_level;
    hamming_dist_[cpt_] = hamming_dist;
    dist_next_prot_[cpt_] = dist_next_prot;
    cpt_++;
  }
  
  std::vector<int32_t> start_pos_;
  std::vector<int32_t> length_;
  std::vector<double> basal_level_;
  std::vector<int32_t> hamming_dist_;
  std::vector<int32_t> dist_next_prot_;
  int nb_prot_;
  int cpt_ = 0;
  int dna_length_;
};

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);
ProteinMap* compute_protein_map(Individual* indiv);
Protein_List* check_SWITCH(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos);
Protein_List* check_S_INS(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos, int32_t len);
Protein_List* check_S_DEL(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos, int32_t len);
Protein_List* check_DUP(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2,  int32_t pos3);
Protein_List* check_DEL(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2);
Protein_List* check_TRANS(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2, int32_t pos3, int32_t pos4, bool invert);
Protein_List* check_INV(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2);

void try_all_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos);
void compute_fates_of_all_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos, ofstream& fate_list_file);
void try_all_duplications_and_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos, ofstream& fate_list_file, int32_t neutral_pos);

// Command-line option variables
static char* lineage_file_name = nullptr;
static bool verbose = false;

static long pt_begin = -1;
static long pt_end = -1;
long int i_prot = 0;











int find_first_neutral_insertion_location(Individual* initial_indiv)
{
	initial_indiv->Evaluate();
	//  initial_indiv->compute_statistical_data();
	//  initial_indiv->compute_non_coding();

	Individual* mutant = NULL;
	int32_t pos1;
	int32_t mut_length;
	int32_t initial_len;
	double metabolic_error_after_switch = -1.0;
        int32_t u = 0;


	//   int32_t result_tab[30000];
	//   int32_t neutral_pos[30000];

	for (const auto& gu: initial_indiv->genetic_unit_list()) {
		initial_len = gu.dna()->length();

		printf (" (genome length: %ld)\n",initial_len);


		// initialize neutral_pos table
		int32_t nb_neutral_pos = 0;



		for (pos1 = 0;pos1 < initial_len;pos1++)
		{
		  
		  // testing insertion of a terminator
		  mutant = new Individual(*initial_indiv);
		  mutant->genetic_unit(u).dna()->do_small_insertion(pos1,30,"000000000000000111111111111111");
		  // Evaluate the mutant, compute its statistics
		  mutant->ReevaluateInContext(initial_indiv->habitat());
		  metabolic_error_after_switch = mutant->dist_to_target_by_feature(METABOLISM);
		  delete mutant;

		  
		  if (metabolic_error_after_switch == initial_indiv->dist_to_target_by_feature(METABOLISM))
		  {
			return pos1;
		  }
		  
		}
		++u;
	}
	
	return -1;
}
























void update_protein_list_parameters(GeneticUnit& gen_unit, Protein_List* protein_list)
{
  
  std::list<Protein>& proteins_lead = gen_unit.protein_list(LEADING);
  for (const auto& prot: proteins_lead)
  {
    for (int i = 0; i < protein_list->nb_prot_; i++){
      if(prot.shine_dal_pos() == protein_list->pro_list_[i].start_pos_){
        protein_list->pro_list_[i].m = prot.mean();
        protein_list->pro_list_[i].w = prot.width();
        protein_list->pro_list_[i].h = prot.height();
        protein_list->pro_list_[i].c = prot.concentration();
        protein_list->pro_list_[i].functional = prot.is_functional() ? "functional" : "non functional";
        protein_list->pro_list_[i].strand = "leading";
      }
    }
  }
  //TODO : EVIL DUPLICATED CODE!
  std::list<Protein>& proteins_lag = gen_unit.protein_list(LAGGING);
  for (const auto& prot: proteins_lag)
  {
    for (int i = 0; i < protein_list->nb_prot_; i++){
      if(prot.shine_dal_pos() == protein_list->pro_list_[i].start_pos_){
        protein_list->pro_list_[i].m = prot.mean();
        protein_list->pro_list_[i].w = prot.width();
        protein_list->pro_list_[i].h = prot.height();
        protein_list->pro_list_[i].c = prot.concentration();
        protein_list->pro_list_[i].functional = prot.is_functional() ? "functional" : "non functional";
        protein_list->pro_list_[i].strand = "lagging";
      } 
    }
  }
  
}








ExpManager* exp_manager;
int32_t wanted_rank = -1;
int32_t wanted_index = -1;
int64_t num_gener = 0;
int32_t mutation_type = 0;
int32_t nb_mutants = -1;







int main(int argc, char* argv[]) {

  interpret_cmd_line_options(argc, argv);
  
  printf("\n"
           "WARNING : Parameter change during simulation is not managed in general.\n"
           "          Only changes in environmental target done with aevol_modify are handled.\n"
           "\n");
  
  
  
  if (SWITCH_DEBUG_MODE)
  {
    cout<<"******************************************"<<endl;
    cout<<"Switch debugging activated.  Will output all switch results for all genes at generation "<<SWITCH_DEBUG_GENERATION<<endl;
    cout<<"******************************************"<<endl;
    
  }

  if (FATE_PROB_MODE) {
      cout << "******************************************" << endl;
      cout << "Fates calculation mode activated.  Will output fate probabilities of all possible switchs of each gene at generation " << SWITCH_DEBUG_GENERATION << endl;
      cout << "******************************************" << endl;
  }
  
  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
    exit(EXIT_FAILURE);
  }
  
  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;
  
  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank));
  
  if (pt_begin == -1) pt_begin = t0;
  if (pt_end == -1) pt_end = t_end;
  
  if (verbose) {
    printf("\n\n""===============================================================================\n");
    printf(" Statistics of the ancestors of indiv. %" PRId32
             " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
             final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }
  
  
  
  // =============================
  //  Open the experiment manager
  // =============================
  exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);
  
  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared())
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                            "for per grid-cell phenotypic target");
  auto phenotypicTargetHandler =
    exp_manager->world()->phenotypic_target_handler();
  if (not (phenotypicTargetHandler->var_method() == NO_VAR))
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                            "for variable phenotypic targets");
  
  int64_t backup_step = exp_manager->backup_step();
  
  
  // =========================
  //  Open the output file(s)
  // =========================
  // Create missing directories
  int status;
  status = mkdir("stats/ancestor_stats/", 0755);
  if ((status == -1) && (errno != EEXIST)) {
    err(EXIT_FAILURE, "stats/ancestor_stats/");
  }
  
  // =========================
  //  Create data structure
  // =========================
  ProteinMap** list_prot_map = new ProteinMap*[t_end-t0+1];
  
  
  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  auto* indiv = grid_cell->individual();
  GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
  indiv->Evaluate();
  indiv->compute_statistical_data();
  indiv->compute_non_coding();
  
  
  // RK list of proteins in first generation
  //****************************************************
  Protein_List* protein_list = new Protein_List(indiv->protein_list().size(),indiv->amount_of_dna());
  ProteinMap* list_prot_temp = new ProteinMap(indiv->protein_list().size(),indiv->amount_of_dna());
  //****************************************************
  //int i = time();
  list_prot_map[time()] = compute_protein_map(indiv);
  list_prot_temp = compute_protein_map(indiv);
  /*for (int i_prot = 0; i_prot < list_prot_map[i]->nb_prot_; i_prot++) {
   std::cout<<0<<","<<i_prot<<","<<list_prot_map[i]->start_pos_[i_prot]
            <<","<<list_prot_map[i]->length_[i_prot]
            <<","<<list_prot_map[i]->basal_level_[i_prot]
            <<","<<list_prot_map[i]->hamming_dist_[i_prot]
            <<","<<list_prot_map[i]->dist_next_prot_[i_prot]
            <<","<<list_prot_map[i]->nb_prot_<<","<<list_prot_map[i]->dna_length_<<std::endl;
  }*/
  std::cout << "num of prot:" << list_prot_map[time()]->nb_prot_ << std::endl;
  //protein_list->hello();
  //******************************************************
  for (; i_prot < list_prot_map[time()]->nb_prot_; i_prot++) {
    //protein_list->hello(list_prot_map[time()]->start_pos_[i_prot], list_prot_map[time()]->length_[i_prot], list_prot_map[time()]->basal_level_[i_prot], list_prot_map[time()]->hamming_dist_[i_prot], list_prot_map[time()]->dist_next_prot_[i_prot], std::to_string(i_prot));
    protein_list->initialize_protein(list_prot_map[time()]->start_pos_[i_prot], list_prot_map[time()]->length_[i_prot], list_prot_map[time()]->basal_level_[i_prot],
                                     list_prot_map[time()]->hamming_dist_[i_prot], list_prot_map[time()]->dist_next_prot_[i_prot], i_prot, i_prot);
  }
  
  indiv->Reevaluate();
  /*gen_unit.clear_transcribed_proteins();
  gen_unit.reset_expression();
  gen_unit.locate_promoters();
  gen_unit.do_transcription();
  gen_unit.do_translation();
  gen_unit.compute_phenotypic_contribution();*/
  std::list<Protein>& proteins_lead = gen_unit.protein_list(LEADING);
  for (const auto& prot: proteins_lead)
  {
    for (int i = 0; i < protein_list->nb_prot_; i++){
      if(prot.shine_dal_pos() == protein_list->pro_list_[i].start_pos_){
        protein_list->pro_list_[i].m = prot.mean();
        protein_list->pro_list_[i].w = prot.width();
        protein_list->pro_list_[i].h = prot.height();
        protein_list->pro_list_[i].c = prot.concentration();
        protein_list->pro_list_[i].functional = prot.is_functional() ? "functional" : "non functional";
        protein_list->pro_list_[i].strand = "leading";
      }
    }
  }
  //TODO : EVIL DUPLICATED CODE!
  std::list<Protein>& proteins_lag = gen_unit.protein_list(LAGGING);
  for (const auto& prot: proteins_lag)
  {
    for (int i = 0; i < protein_list->nb_prot_; i++){
      if(prot.shine_dal_pos() == protein_list->pro_list_[i].start_pos_){
        protein_list->pro_list_[i].m = prot.mean();
        protein_list->pro_list_[i].w = prot.width();
        protein_list->pro_list_[i].h = prot.height();
        protein_list->pro_list_[i].c = prot.concentration();
        protein_list->pro_list_[i].functional = prot.is_functional() ? "functional" : "non functional";
        protein_list->pro_list_[i].strand = "lagging";
      } 
    }
  }
  
  std::ofstream proteins_list_file;
  proteins_list_file.open("proteins_list_after_events.csv",std::ofstream::trunc);
  proteins_list_file<<"event_type,protein_id, protein parent id,shine_dal,length,concentration,hamming_dist,dist_next_protein,nb_proteins,dna_length,m,w,h,c,f, generation"<<std::endl;
  for (int i = 0; i < protein_list->nb_prot_; i++) {
    proteins_list_file<<"initial"<<","<<protein_list->pro_list_[i].prot_id_
                      <<","<<protein_list->pro_list_[i].prot_parent_id_
                      <<","<<protein_list->pro_list_[i].start_pos_
                      <<","<<protein_list->pro_list_[i].length_
                      <<","<<protein_list->pro_list_[i].basal_level_
                      <<","<<protein_list->pro_list_[i].hamming_dist_
                      <<","<<protein_list->pro_list_[i].dist_next_prot_
                      <<","<<protein_list->nb_prot_
                      <<","<<protein_list->dna_length_
                      <<","<<protein_list->pro_list_[i].m
                      <<","<<protein_list->pro_list_[i].w
                      <<","<<protein_list->pro_list_[i].h*protein_list->pro_list_[i].c
                      <<","<<protein_list->pro_list_[i].c
                      <<","<<protein_list->pro_list_[i].functional
	    	      <<","<<time()<<std::endl;
  }
  
  std::ofstream proteins_map_file_event;
  proteins_map_file_event.open("proteins_map_event.csv",std::ofstream::trunc);
  proteins_map_file_event<<"mut_type,protein_id,shine_dal,length,concentration,hamming_dist,dist_next_protein,nb_proteins,dna_length"<<std::endl;
  for (int i = 0 ; i < list_prot_temp->nb_prot_ ; i++) {
    proteins_map_file_event<<"initial"<<","<<i<<","<<list_prot_temp->start_pos_[i]
                           <<","<<list_prot_temp->length_[i]
                           <<","<<list_prot_temp->basal_level_[i]
                           <<","<<list_prot_temp->hamming_dist_[i]
                           <<","<<list_prot_temp->dist_next_prot_[i]
                           <<","<<list_prot_temp->nb_prot_<<","<<list_prot_temp->dna_length_<<std::endl;
  }
  //std::cout << "num of prot_temp:" << list_prot_temp->nb_prot_ << std::endl;
  //******************************************************
  
  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;
  int32_t index;
  ExpManager* exp_manager_backup = nullptr;
  int32_t unitlen_before;
  double metabolic_error_before;
  double impact_on_metabolic_error;
  char mut_descr_string[255];
  
  
  aevol::AeTime::plusplus();
  while (time() <= t_end)
  {
#ifdef __REGUL
    printf("Protein map is not supported yet\n");
    exit(-1);
#else
    rep = new ReplicationReport(lineage_file, indiv);
#endif
    
    index = rep->id(); // who we are building...
    
    if (verbose)
      printf("Rebuilding ancestor at generation %" PRId64
               " (index %" PRId32 ")...", time(), index);
    
    indiv->Reevaluate();
    
    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;
    
    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();
    
    // TODO(dpa) The following 3 for loops should be factorized.
    // However, this is not as easy as it sounds :-D
    // see std::list::splice
    
    //********************************
    std::cout <<"at gen:   "<< time() <<"\n";
    
    if (SWITCH_DEBUG_MODE && time() == SWITCH_DEBUG_GENERATION)
    {
      cout<<"Reached generation "<< time() <<".  Outputting all switch effects on every gene."<<endl;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++)
      {

        try_all_switches(indiv, gen_unit, protein_list->pro_list_[i].start_pos_);
      }
      cout<<"Done, press any key to continue."<<endl;
      string x = "";
      cin>>x;
    }

    if (FATE_PROB_MODE && time() == SWITCH_DEBUG_GENERATION)
    {
        cout << "Reached generation " << time() << ".  Outputting fate probabilities of all switch effects on every gene." << endl;
        
        
        //output the genome at this moment, so we know what was analyzed
        std::ofstream wtfile;
        wtfile.open("WT_analyzed.txt");
        const char* seq = gen_unit.dna()->data();

	for (int i = 0; i < gen_unit.dna()->length(); i++)
	{
		wtfile<<seq[i];
	}
	wtfile.close();
	
        
        
        std::ofstream fate_list_file;
        fate_list_file.open("fate_probs_list_all_switches_file.csv", std::ofstream::trunc);
        fate_list_file << "protpos, insertionpos, flipoffset1, flipoffset2, dead, subfunc, neofunc, cons, pseudo, spec, cons_loss, neo_loss, mg, hg, wg, ma, ha, wa, mb, hb, wb,ig_a,ia_g,ig_b,ib_g,ig_a_plus_b,ia_plus_b_g,pa,pb,fitness_g, fitness_postdup,fitness_flip_a, fitness_flip_b, fitness_flip_ab,metabolicerr_g,metabolicerr_postdup,metabolicerr_flip_a,metabolicerr_flip_b,metabolicerr_flip_ab" << std::endl;
        cout<<"Searching for neutral insertion location"<<endl;
        
        int32_t neutral_pos = find_first_neutral_insertion_location(indiv);
        
        if (neutral_pos == -1)
        {
        	cout<<"Error: no neutral position found!"<<endl;
        	return 0;
        }
        
        for (int i = 0; i < protein_list->nb_prot_; i++)
        {
        	
		try_all_duplications_and_switches(indiv, gen_unit, protein_list->pro_list_[i].start_pos_, fate_list_file, neutral_pos);
            //compute_fates_of_all_switches(indiv, gen_unit, protein_list->pro_list_[i].start_pos_, fate_list_file);
        }
        
        return 0;
        
    }
    
    //**********************************
    dnarep.iter_muts([&](const auto& mut) {
      // Apply mutation
      const Mutation& mut1 = *mut;
      //gen_unit.dna()->undergo_this_mutation(*mut);
      std::string mut_type;
      //*********************************
      std::vector<int32_t> output_mut_values(4);
      bool output_mut_bool;
      tie(output_mut_values[0], output_mut_values[1], output_mut_values[2], output_mut_values[3], output_mut_bool) = gen_unit.dna()->undergo_this_mutation(*mut);
      indiv->Reevaluate();
      indiv->compute_statistical_data();
      indiv->compute_non_coding();
      list_prot_temp = compute_protein_map(indiv);
      switch(mut1.mut_type())
      {
      case SWITCH : {
        mut_type = "SWITCH";
        std::cout << "SWITCH\n";
        // make sure it will compute protein list after each event!
        //if(list_prot_temp->nb_prot_ != protein_list->nb_prot_)
        // change protein list based on the new list (output_mut_values[0])
        protein_list = check_SWITCH(protein_list, list_prot_temp, output_mut_values[0]);
        break;
      }
      case S_INS: {
        mut_type = "S_INS";
        std::cout << "S_INS\n"; 
        std::cout << "pos1 and pos2:" <<output_mut_values[0]<<" & " <<output_mut_values[1]<< std:: endl;
        // call small insertion function to check any protein does change or not (output_mut_values[0] & [1])
        protein_list = check_S_INS(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1]);
        break;
      }
      case S_DEL: {
        mut_type = "S_DEL";
        std::cout << "S_DEL\n";
        std::cout << "pos1 and pos2:" <<output_mut_values[0]<<" & " <<output_mut_values[1]<< std:: endl;
        // call small deletion function to check any protein does change or not (output_mut_values[0] & [1])
        protein_list = check_S_DEL(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1]);
        break;
      }
      case DUPL: {
        mut_type = "DUPL";
        std::cout << "DUPL\n";
        std::cout << "pos1 and pos2:" <<output_mut_values[0]<<"&" <<output_mut_values[1]<<"pos3:" << output_mut_values[2]<< std:: endl;
        protein_list = check_DUP(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1],  output_mut_values[2]);
        break;
      }
      case DEL:{
        mut_type = "DEL";
        std::cout << "DEL\n";
        std::cout << "pos1 and pos2:" <<output_mut_values[0]<<"&" <<output_mut_values[1] << std:: endl;
        protein_list = check_DEL(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1]);
        break;
      }
      case TRANS:{
        mut_type = "TRANS";
        std::cout << "TRANS\n";
        std::cout << "pos1 and pos2: " <<output_mut_values[0]<<" & " <<output_mut_values[1] << " pos3 and pos 4: " <<output_mut_values[2]<<" & " <<output_mut_values[3] << " bool: "<< output_mut_bool<<std:: endl;
        protein_list = check_TRANS(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1], output_mut_values[2], output_mut_values[3], output_mut_bool);
        break;
      }
      case INV:{
        mut_type = "INV";
        std::cout << "INV\n"; 
        protein_list = check_INV(protein_list, list_prot_temp, output_mut_values[0], output_mut_values[1]);
        std::cout << "pos1 and pos2:" <<output_mut_values[0]<<" & " <<output_mut_values[1] << std:: endl;
        break;
      }
      case INS_HT:{
        mut_type = "INS_HT";
        std::cout << "INS_HT\n";
        std::cout<< "ERROR!!! INS_HT"<< std::endl;
        return;
        break;
      }
      case REPL_HT:{
        mut_type = "REPL_HT";
        std::cout << "REPL_HT\n";
        std::cout<< "ERROR!!! REPL_HT"<< std::endl;
        return;
        break;
      }
      }
      indiv->Reevaluate();
      /*gen_unit.clear_transcribed_proteins();
      gen_unit.reset_expression();
      gen_unit.locate_promoters();
      gen_unit.do_transcription();
      gen_unit.do_translation();
      gen_unit.compute_phenotypic_contribution();*/
      update_protein_list_parameters(gen_unit, protein_list);
      
      for (int i = 0; i < protein_list->nb_prot_; i++) {
        if(i == 0){
          proteins_list_file<<mut_type<<","<<protein_list->pro_list_[i].prot_id_
                            <<","<<protein_list->pro_list_[i].prot_parent_id_
                            <<","<<protein_list->pro_list_[i].start_pos_
                            <<","<<protein_list->pro_list_[i].length_
                            <<","<<protein_list->pro_list_[i].basal_level_
                            <<","<<protein_list->pro_list_[i].hamming_dist_
                            <<","<<protein_list->pro_list_[i].dist_next_prot_
                            <<","<<protein_list->nb_prot_
                            <<","<<protein_list->dna_length_
                            <<","<<protein_list->pro_list_[i].m
                            <<","<<protein_list->pro_list_[i].w
                            <<","<<protein_list->pro_list_[i].h*protein_list->pro_list_[i].c
                            <<","<<protein_list->pro_list_[i].c
                            <<","<<protein_list->pro_list_[i].functional
			    <<","<<time()<<std::endl;
        }
        else{
          proteins_list_file<<" "<<","<<protein_list->pro_list_[i].prot_id_
                            <<","<<protein_list->pro_list_[i].prot_parent_id_
                            <<","<<protein_list->pro_list_[i].start_pos_
                            <<","<<protein_list->pro_list_[i].length_
                            <<","<<protein_list->pro_list_[i].basal_level_
                            <<","<<protein_list->pro_list_[i].hamming_dist_
                            <<","<<protein_list->pro_list_[i].dist_next_prot_
                            <<","<<protein_list->nb_prot_
                            <<","<<protein_list->dna_length_
                            <<","<<protein_list->pro_list_[i].m
                            <<","<<protein_list->pro_list_[i].w
                            <<","<<protein_list->pro_list_[i].h*protein_list->pro_list_[i].c
                            <<","<<protein_list->pro_list_[i].c
                            <<","<<protein_list->pro_list_[i].functional
		            <<","<<time()<<std::endl;
        }
      }
      
      
      for (int i = 0 ; i < list_prot_temp->nb_prot_ ; i++) {
        if(i == 0){
          proteins_map_file_event<<mut_type<<","<<i<<","<<list_prot_temp->start_pos_[i]
                                 <<","<<list_prot_temp->length_[i]
                                 <<","<<list_prot_temp->basal_level_[i]
                                 <<","<<list_prot_temp->hamming_dist_[i]
                                 <<","<<list_prot_temp->dist_next_prot_[i]
                                 <<","<<list_prot_temp->nb_prot_<<","<<list_prot_temp->dna_length_<<std::endl;
        }
        else{
          proteins_map_file_event<<" "<<","<<i<<","<<list_prot_temp->start_pos_[i]
                                 <<","<<list_prot_temp->length_[i]
                                 <<","<<list_prot_temp->basal_level_[i]
                                 <<","<<list_prot_temp->hamming_dist_[i]
                                 <<","<<list_prot_temp->dist_next_prot_[i]
                                 <<","<<list_prot_temp->nb_prot_<<","<<list_prot_temp->dna_length_<<std::endl;
        }
      }
      //************************************
    });
    
    // 3) All the mutations have been replayed, we can now evaluate the new individual
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();
    
    list_prot_map[time()] = compute_protein_map(indiv);
    
    if (verbose) printf(" OK\n");
    
    delete rep;
    
    aevol::AeTime::plusplus();
  }
  
  gzclose(lineage_file);
  
  std::ofstream proteins_map_file;
  proteins_map_file.open("proteins_map_gen.csv",std::ofstream::trunc);
  proteins_map_file<<"generation,protein_id,shine_dal,length,concentration,hamming_dist,dist_next_protein,nb_proteins,dna_length"<<std::endl;
  
  for (int i = t0; i <= t_end; i++) {
    for (int j = 0; j < list_prot_map[i]->nb_prot_; j++) {
      proteins_map_file<<i<<","<<j<<","<<list_prot_map[i]->start_pos_[j]
                       <<","<<list_prot_map[i]->length_[j]
                       <<","<<list_prot_map[i]->basal_level_[j]
                       <<","<<list_prot_map[i]->hamming_dist_[j]
                       <<","<<list_prot_map[i]->dist_next_prot_[j]
                       <<","<<list_prot_map[i]->nb_prot_<<","<<list_prot_map[i]->dna_length_<<std::endl;
    }
  }
  
  proteins_map_file.flush();
  proteins_map_file.close();
  
  // Additional outputs
  
  delete exp_manager;
  delete indiv;
  
  return EXIT_SUCCESS;
}







Protein* find_protein_at_startpos(Individual* indiv, GeneticUnit& gen_unit, int startpos)
{
	Protein* protein = nullptr;


	    for (auto& strand_id: {LEADING, LAGGING}) {
	      auto& strand = gen_unit.protein_list(strand_id);
	      for (auto& p: strand) {
		
		if (p.shine_dal_pos() == startpos)
		{
			protein = &p;
		}
	      }
	    }
	  
	  return protein;
}



string get_codon_type(char c1, char c2, char c3)
{
	string s = "   ";
	s[0] = c1;
	s[1] = c2;
	s[2] = c3;
	int8_t v = std::stoi(s, nullptr, 2);

	switch(v)
	{
		case CODON_M0: return "M0";
		case CODON_M1: return "M1";
		case CODON_W0: return "W0";
		case CODON_W1: return "W1";
		case CODON_H0: return "H0";
		case CODON_H1: return "H1";
		case CODON_START: return "BE";
		case CODON_STOP: return "ST";
	}
	return "XX";
}



char revcomp(char b)
{
	if (b == '0')
		return '1';
	else
		return '0';
}


void output_protein_sequence(Individual* indiv, GeneticUnit& gen_unit, int32_t shine_dal_pos, int32_t first_translated_pos, int32_t last_translated_pos, Strand strand, bool output_details)
{
	const char* seq = gen_unit.dna()->data();

	string details = "";

	int codoncpt = 0;
	//TODO: WILL CRASH ON CIRCULAR GENES
	if (strand == LEADING)
	{
		for (int32_t pos = shine_dal_pos; pos <= last_translated_pos + 3; ++pos)
		{
			if (pos == shine_dal_pos + SHINE_DAL_SIZE)
			{
		  		cout<<" ";
		  		details += " ";
		  	}
		  	
		  	if (pos == first_translated_pos)
		  	{
		  		cout<<" ";
		  		details += " ";
		  	}
		  	
		  	
		  	cout<<seq[pos];
		  	
		  	if (pos < first_translated_pos)
		  		details += " "; 
		  	
		  	if (pos >= first_translated_pos)
		  	{
		  		codoncpt++;
		  		if (codoncpt == 3)
		  		{
		  			codoncpt = 0;
		  			cout<<" ";
		  			
		  			details += get_codon_type(seq[pos - 2], seq[pos - 1], seq[pos]) + "  ";
		  		}
		  	}
		 }
	 }
	 else
	{
		//we must print the reverse complement
		for (int32_t pos = shine_dal_pos; pos >= last_translated_pos - 3; --pos)
		{

			if (pos == shine_dal_pos - SHINE_DAL_SIZE)
			{
		  		cout<<" ";
		  		details += " ";
		  	}
		  	
		  	if (pos == first_translated_pos)
		  	{
		  		cout<<" ";
		  		details += " ";
		  	}
		  	
		  	
		  	
		  	cout<<revcomp(seq[pos]);
		  		
		  	if (pos > first_translated_pos)
		  		details += " "; 
		  		
		  	if (pos <= first_translated_pos)
		  	{
		  		codoncpt++;
		  		if (codoncpt == 3)
		  		{
		  			codoncpt = 0;
		  			cout<<" ";
		  			
		  			details += get_codon_type(revcomp(seq[pos + 2]), revcomp(seq[pos + 1]), revcomp(seq[pos])) + "  ";
		  		}
		  	}
		 }
	 }
	 if (output_details)
	 	cout<<endl<<details;
}

















void print_protein_locations(Individual* indiv, GeneticUnit& gen_unit, bool details = false)
{
	vector<int> sss;
        		for (auto& strand_id: {LEADING, LAGGING}) {
			      auto& strand = gen_unit.protein_list(strand_id);
			      for (auto& p: strand) {
				
				sss.push_back(p.shine_dal_pos());
				if (details)
				{
					cout<<p.height() * p.concentration()<<endl;
				}
			      }
			    }
			    sort(sss.begin(), sss.end());
			    for (auto p : sss)
			    {
			    	cout<<p<<" ";
			    }
			    cout<<endl;
}










void print_dna_sequence(Individual* indiv, GeneticUnit& gen_unit, int start, int end)
{
	
    	const char* seq = gen_unit.dna()->data();

	for (int i = start; i <= end; i++)
	{
			cout<<seq[i];
	}
	cout<<endl;
	
}


void get_rna_positions(Individual* indiv, GeneticUnit& gen_unit, Protein* protein, int32_t &firstpos, int32_t &lastpos)
{
	bool is_leading_strand = (protein->strand() == LEADING);
	Rna* rna = *protein->rna_list().begin();
        {

            int32_t promoter_pos = rna->promoter_pos();
            int32_t transcription_start = rna->first_transcribed_pos();
            int32_t transcription_end = rna->last_transcribed_pos();
            int32_t rna_length = transcription_end - promoter_pos + 1;

            if (is_leading_strand)
            {
                firstpos = promoter_pos;
                lastpos = transcription_end;
            }
            else
            {
                
            }
        }
}







void output_ab_data(long double fitness_orig, long double fitness_dup, long double fitness_a, long double fitness_b, long double fitness_ab, long double gm, long double gh, long double gw, Protein* protein_a, Protein* protein_b, ofstream& fate_list_file,
long double orig_metabolic_error, long double metabolic_error_after_dup, long double metabolic_error_flip_a, long double metabolic_error_flip_b, long double metabolic_error_flip_ab)
{
	long double g_points[3];
	long double a_points[3];
	long double b_points[3];
	
	
	g_points[0] = gm;
	g_points[1] = gh;
	g_points[2] = gw;
	
	if (protein_a)
	{
		a_points[0] = protein_a->mean();
		a_points[1] = protein_a->height() * protein_a->concentration();
		a_points[2] = protein_a->width();
	}
	
	if (protein_b)
	{
		b_points[0] = protein_b->mean();
		b_points[1] = protein_b->height() * protein_b->concentration();
		b_points[2] = protein_b->width();
	}
	
	
	if (protein_a && protein_b)
	{
		vector<long double> params = ProbabilitiesCalculation(g_points, a_points, b_points);

		long double subfunc = params[0];
		long double neofunc = params[2];
		long double cons = params[1];
		long double pseudo = params[3];
		long double spec = params[4];
		long double cons_loss = params[17];
		long double neo_loss = params[18];

		long double ig_a = min(1, params[5]);
		long double ig_b = min(1, params[6]);
		long double ia_g = min(1, params[7]);
		long double ib_g = min(1, params[8]);

		long double ig_a_plus_b = min(1, params[9]);
		long double ia_plus_b_g = min(1, params[10]);
		long double pa = params[11];
		long double pb = params[12];



		fate_list_file << std::fixed;
		fate_list_file << std::setprecision(4);
		fate_list_file << subfunc 
			<< "," << neofunc
			<< "," << cons
			<< "," << pseudo
			<< "," << spec
				<< "," << cons_loss 
				<< "," << neo_loss
		
		<< "," << g_points[0] << "," << g_points[1] << "," << g_points[2] 
		<< "," << a_points[0] << "," << a_points[1] << "," << a_points[2]
		<< "," << b_points[0] << "," << b_points[1] << "," << b_points[2];
		
		fate_list_file << std::fixed;
		fate_list_file << std::setprecision(4);
		

		

		fate_list_file << ", "<<ig_a<<", "<<ia_g<<", "<<ig_b<<", "<<ib_g<<", "<<ig_a_plus_b<<", "<<ia_plus_b_g<<", "<<pa<<", "<<pb;
	}		
	else
	{
		fate_list_file << " , , , , , , , , , , , , , , , , , , , , , , , ";
	}	

	
	long double tenpower = 10000000000;
	
	fate_list_file << std::scientific;
	fate_list_file << std::setprecision(24);
		fate_list_file<< ", " << fitness_orig * tenpower <<", "<< fitness_dup * tenpower <<", "<< fitness_a * tenpower <<", " << fitness_b * tenpower << ", " << fitness_ab * tenpower;
	
	fate_list_file << std::fixed;
	fate_list_file << std::setprecision(8);
	
	fate_list_file << ", " << orig_metabolic_error << ", " << metabolic_error_after_dup << ", " << metabolic_error_flip_a << ", " <<
			  ", " << metabolic_error_flip_b << ", " << metabolic_error_flip_ab;
	
	fate_list_file<<std::endl;
	

}
















void try_all_duplications_and_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos, ofstream& fate_list_file, int32_t neutral_pos)
{

    

    indiv->Reevaluate();
	
    indiv->compute_fitness(indiv->phenotypic_target());
    double orig_fitness = indiv->fitness();
	
    double orig_metabolic_error = indiv->dist_to_target_by_feature(METABOLISM);	

    /*
    //this was temp code to check whether initial genes we found corresponded to genes of the WT.  I leave it here in case it becomes useful again.
    string wt6 = "101101001101010010110111011011010101010010001011001000110100000110011010100111110010010011010010110100101000111011000010010100011010100111001001101110111111001010110110110110001010010100000101111010111111111110111001110000011100011000110110000000101101101110110110110001110110010100001101010000100110101111110000100111001011101100110010111010100111001000001110111110010011011100100001110011001010100110110110110111110001000000011000110100100100101111110001010001011001110100101000101111000110111001110011101110111100000100100111010101010101100110001001101001011011010000001101010011001111111011001001111100100111101110011111110100111011111110011000101110111111111111011000000100000100110011000000011001101010011111101001001111010100011000010000011111001001001010011100100000111000110011010100111100000010011100111010100110110001001101010110110101010010101111101000001011010011011011000011101010100010011110101010001001011001110101011011011100001101010011001111111011001000101011100110000001100110101001111100100100100100111010110110011000100010011010011100111010100110011101100001100100010001101011011100100000101100110101001111000100100111110101111000011000010101011001110010110000001001101111100000110101001100111111101100100111110011011001001101000110111011010110010000010110000000011101101011111011000001101010011100110110011110111100101101001110001111010101000111000100100111000010100110001101001011011101101000101011101100001011110011100100111100110110110000011010100110011111110011110001110100011101110111110101111111001100110010001100101000111001100100101111110110100110010101001111011000100010110000001000100110101000111001100100100000010101000111110110001010011001000000011001101010011110000001001101111000011011001010111101111001000100110011110000110110001110001001111100100110111100011001110110010110101110110111011000011010100110011111110110001101101100001001001101000000111110100010011100101101111010111101101010010001101111010001010110101111110000101010100010001010110001000101010001100011010011011101110011001001010011100010010001111101011011100110010100100010101011010111000000101001101110010000100110111000110110111000101010111101111100010010011101001101111100010101011100110101101111110110010000000111010000101011010110101000001101101100110110011001001100110100010100001111010111100001110101111110011001100000011101010111110110010101010011100100000101100110101001111110100100100101111101011011110000001100101000001100101110101010101101101111111101001001110011101111100010010010000001001101011111000001001010001010110001011010001010010101000101100101111010110110000010000110010001111100100011110100100101110011011110000010011100111010111011100101110010101110110100101100010110011010101000110100110001001011010011110101010010011100010011000100001000111011011110011011001100100010100100111000011010011010100010001100110110111000111010100110111011111100111100110101010001111111001111101101101000000110100110100110101001100111111101100011101110111111101010110110000000101010111101110101010010011101011100100111000000011100010101010100000100111000110011110001010011100001111010001100000010011100111010100110110101001101010110000101001011000100101101001110010111011111110110101111101001111000111001001111110000000010000010111100101010010001010100000100101001010000110111110011110100001101000010101110111011100011000101111101000110111100110101100011000100110100001111001000110011011011001100001101010011001111111010011000111110110111000000010000100000111010101010010011100010100100010000010100010010010010010111000000010101010101111100110010111011011001100001101010011001111111011001100110111011100000000010100011111001111111100110011101110100100110010000000110011010100111100000010011001011001010111010101110010010111000001010111110111001101101100000111110111100111010100110100101100101000100000101101000011000000100010011010100011111110010010001111100100101101011001111111110000110111011010001100111001001100010010100001101101101100001110101001101110111111001011000110011110100101001011111000010110001100011011011100011101010011011101110010100101101010110110101000011000010000011111001001001001010101010010001110010010011100111010100011001011101100100100101100111111010110100011001011111011000101010000001011010001101101011110000000001011001000101111000011111101000010111010100000110000001110000110001011000100010001101010111110001110111101100110111101110000100001101100000001000011001111110101011110110001101000000101110110011000010100011111110100110101111101001100110011001101000001011000001101011001001000010000011111000010100100111111010100010010100011111000111001000110010000110101001001100101101111111000110011001001111100100001001111100111010111110010011011000101101100101100110000001000100001010000010101111111011101000011011111111110100101011001111010011000001111101011011011100001101010011001111111001011000111001001101111000000110101001100111000101110110011001000000110111001111011001000110011000101101110010011100100011010011001000111011111100110100010110011000010010010110011000100110111110101011111111100110111100000011010100110011111100100110110010001000101001011011110101000000100101100101110011110001011100011101011111111101100000110010010011101100100000001100110101001110111110100100100100101001101011001011100100000001100110001110101001111100100100111100010111011011100100100001100011001110110110101100000001100110101001111100100100110011011111111111010001110100000001100110101001111000100100101000001011000110110111100101010000010101000100000011001111111000000010001111000000010100010010100010111011011010101011010001001001111011010110111000101110011011100110000110110110000011010100011100110011111011000110001011010000110110111000101010111101111100010010011101011100010010001100001100100001001000011110011111111111100010111010101000110110111010001000010101011111110010010101000011101000111010101001001110100011110111000001101101000001110101010001001100111100111100101100001000100000000010001100010111001110101011111010100100011011010011110100011111101001100001110110011101110000000110011010100110100111111000010011001011111001010010111000001001011100011101111010100011011101110111100111110010111000001110011110110110000001111000011110110011001100110010001000000110011001001010011010100111100101000110011101001011001011111110010111110000010101101111101010100001111001100100110000010100011011100111101010100010001001001000110110010011111111001011011010010101001001111100000101111000100000101001101100000010110001011101110000011010110111100001101101110000110101001100111111001001000101100011011001000111011011110000001000011010101100110001101010111110100101010110111111101001001101100111011101011001101110110010010100001111100000001011011011001101111010001110100111100111100111001000010001010011001101010011101110011001001001100011100101111100000001000110001111111011110111010000111000110010110001111101110100010011111101111100001101101100110100010111011100110001101010100100100100111101001111001111110010010011110010110011100011111110001000110000100000111110010010011000011101101101101000001100101100100000001100110101001111100100100101010101001001001101010011000111010111110000001111110011010000011001111100000110000011011100001011100010001101100001100100100100101001101101110110110011000101001101000101111000011110011111111100110000100110000001101000110100101000001010010010010101110110110100000010110011101011011000011011100110111111010010010111000011110011000000000110010001101011101001100000011101011010011110110000011001011000011100100101101001010111100001011000011110011010101100101000000011001010001100011011010000001101010011001111111011000111000010010011000111111110000100011101001001000111100110111011110111001001000101100011100010100001110000000011011010111011111011100011101100010000110000110110101001101000001000000100001010000011011010001011110100110111010101000111100100100110100101111101000001010110101011110100111011100011011101110001111000001111100110001101000100111011010111001110010110000000110011010100111110000010010010001000110011001111011100111010101110001111101101101001011100010011000111001000000011001101010011111001001001010001110010000100010011010101011001101000001100000010001001101010001111001001001001010010001111011000110110101001110111101101100111101010010100110001000011111011100011000100110111000000001010000010011100111011101100111101001001110001101110000001011000100011111000010001001110101110010101011100001010010001101000100110010100001000001111101000100101001001010011101101110100011001110111101011110101110001101000111100101001111111001001110010010111011011100110110100000011010100110011111110110010001111000110110001100001111010001010010001111011010000001001001111110110110101010001111111000011001001111101001011001001011000111110111100010010111001000111001001110001100011011101110100100110010000000110011010100111100000010011010011110101110010110110001000101111011111110010000111011001110010110111000000001101101100000110101001100111111101100110110100011101001010101011010101000100010010000000110110100000100111001110011011001111010010100100011001111010001110111000101010111011010001110101000110100001111011000110010010100011001011110000001001101001101111101110111000100101000110011001010000111000101101101110000110101001100111111101111110100011001101101110000110101001100111110110111101100011100100110111101011001111111110011111100110101100011001010001001000001011010000111101110111011101000001011001011110110000000001010000010011100111001101100111101001010010001100111101000111011100010101010010000010000010100011011011100011101011001110110011001101001011000000001001000000110001101001101101110101001010010101010011001111111110110011010010001101101010101000011001101111000001101010001101101010110110110110000011111011110011101001011110010101111110011000111000010010100111111010010011110010100100000111001010100000100110000000101011100011010001000100101101100111001010101010001110001101111000011000001101010011001110111011001010000100000101010000110110111010111111010010011100101101101000011101010001010101010101100110110111000111010100111111010111000111010101100000110110000110101011110110111110000011010100110011011111110010110001110100100111101101100110011101000011100100111101010100010100010011100000011111101011010011110101010010011100011011001100001000110000010000000100101110100011101011100101111101010101010101001011100011000111100000001111101010110010010011010000111101011010101111111100001001010010110101001000110101111010001110111000101010100100000100101000110110001100111111010001111011110000111110110001110110111000010010000001001000000011001101010011111000001001100100010110010000000101010100100010001000110101001011110100110010101100110100010100111110101111101011010101001100000010011001110101011101101011100000101100110001111101101110000000011000111011001110010010001101111101111100010010001100110101001011101001000001010101110010101010101011101110111000010101010010101100111000011001010101111011011011000001101010011001111111011001101101000100101000110001110101111000000111100101110110010110110101100111001000110010101100000000001001101010001111001001001001000000010100011001011111010110011101000000001011011010110000001111111111000110111010111011110000110100001000101010000110100011011111010001001100111100100101000101001100100010011101111100101100110111000001011101110000000110011010100111110000010010011110101001010001001010010101100101000101101010000000011111100100010110101111000110000101111101011001111010100110101111001001011111101101101010011000011101101101101010000001101111000000101010111011101100111111000011010111110111101011101010101001101111010000011101110100011110010001110110000100010111100";
    
    gen_unit.dna()->do_insertion(0, wt6.c_str(), wt6.length());
    gen_unit.dna()->do_deletion(wt6.length(), gen_unit.dna()->length());
    
	
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();
	
    print_protein_locations(indiv, gen_unit, true);
    
    //print_dna_sequence(indiv, gen_unit, 0, gen_unit.dna()->length());
    
    int xz;
    cin >> xz;
	*/
	
    bool is_leading_strand = true;

    
    int orig_genome_length = gen_unit.dna()->length();
    char* orig_genome_seq = new char[orig_genome_length];
    for (int i = 0; i < gen_unit.dna()->length(); ++i)
    {
    	orig_genome_seq[i] = gen_unit.dna()->data()[i];
    }

    

    Protein* protein = find_protein_at_startpos(indiv, gen_unit, startpos);

    if (!protein)
    {
        cout << "PROTEIN AT POS=" << startpos << " NOT FOUND" << endl;
        
        return;
    }

    is_leading_strand = (protein->strand() == LEADING);
    long double origm = protein->mean();
    long double origw = protein->width();
    long double origh = protein->height() * protein->concentration();
    int32_t shine_dal_pos_orig = protein->shine_dal_pos();
    int32_t first_translated_pos_orig = protein->first_translated_pos();
    int32_t last_translated_pos_orig = protein->last_translated_pos();
    int32_t nbcodons = protein->length();
    
    int protein_length = last_translated_pos_orig - shine_dal_pos_orig + 1 + 3;
    
    
    //TODO TODO
    if (!is_leading_strand)
    {
	return;
    }



    
    cout << "*********************************************************************" << endl;
    cout << endl << "ANALYZING PROTEIN AT STARTPOS=" << startpos << endl;
    cout<<"Init fitness="<<orig_fitness<<"   init metabolic_err="<<orig_metabolic_error<<endl;
    cout << "Inserting at neutral pos = "<<neutral_pos<<endl;

    cout << "PROM_SEQ=" << PROM_SEQ << "    PROM_MAX_DIFF=" << (int)PROM_MAX_DIFF << " SHINE_DAL_SEQ=" << SHINE_DAL_SEQ << endl << endl;
    

    cout << "startpos=" << (is_leading_strand ? "leading" : "lagging") << endl;
    cout << "shine_dal_pos=" << shine_dal_pos_orig << "   first_translated_pos=" << first_translated_pos_orig << "   last_translated_pos=" << last_translated_pos_orig << endl;
    cout << "m=" << origm << ", ";
    cout << "w=" << origw << ", ";
    cout << "h=" << origh << endl;
    
    cout<<"PROTEIN SEQUENCE: ";
    print_dna_sequence(indiv, gen_unit, startpos, startpos + protein_length - 1);

    long double g[3] = { origm, origh, origw };
    long double a[3];
    long double b[3];
    bool a_dead, b_dead;

    cout << endl;

    
	
	
    if (true || SWITCH_DEBUG_SHOW_SEQUENCES)
    {
    	cout<<"PROTEIN LOCATIONS BEFORE DUP"<<endl;	
    	print_protein_locations(indiv, gen_unit);
    }
    
    
    int32_t posstart_rna = 0;
    int32_t posend_rna = 0;	
    int32_t pos_insertion = 0; 
    

    get_rna_positions(indiv, gen_unit, protein, posstart_rna, posend_rna);

    std::string to_insert_prefix = "000000000000000111111111111111";
    std::string to_insert = to_insert_prefix;
    for (int i = 0; i < PROM_SIZE; ++i)
    {
    	to_insert += gen_unit.dna()->data()[posstart_rna + i];
    }
    for (int i = 0; i < protein_length; ++i)
    {
    	to_insert += gen_unit.dna()->data()[startpos + i];
    }
    to_insert += "000000000000000111111111111111";
    
    
    
    int32_t startpos_afterdup1 = startpos;
    if (neutral_pos < startpos)
    	startpos_afterdup1 = startpos + to_insert.length();
    
    int32_t startpos_afterdup2 = neutral_pos + to_insert_prefix.length() + PROM_SIZE; 
    
    
    if (!is_leading_strand)
    {
        return;
    }
    
    
    
    cout<<"RNA SEQ"<<endl;
    print_dna_sequence(indiv, gen_unit, posstart_rna, posend_rna);
    
    cout<<"TO INSERT LENGTH="<<to_insert.length()<<endl<<to_insert<<endl;
    
    
    //we will do switches and dups, which will delete the current protein object.  To avoid problems, we set it to nullptr
    protein = nullptr;
    
    
    
    
    /*if (startpos == 5559)
    {
        cout<<"SEQ BEFORE"<<endl;
	print_dna_sequence(indiv, gen_unit, 0, gen_unit.dna()->length());
    }
    else
    return;*/
    

    bool output_mut_bool;
    output_mut_bool = gen_unit.dna()->do_small_insertion(neutral_pos,to_insert.length(), to_insert.data());
    if (!output_mut_bool)
    {
       cout<<"DODUP FAILED"<<endl;
       int zz;
       cin >>zz;
    }
    
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    indiv->compute_fitness(indiv->phenotypic_target());
    long double fitness_after_dup = indiv->fitness();
    long double metabolic_error_after_dup = indiv->dist_to_target_by_feature(METABOLISM);

    if (true || SWITCH_DEBUG_SHOW_SEQUENCES)
    {
    	//cout<<"SEQ AFTER DUP"<<endl;
    	//print_dna_sequence(indiv, gen_unit, 0, gen_unit.dna()->length());
    	
    	cout<<"PROTEIN LOCATIONS AFTER DUP"<<endl;	
    	print_protein_locations(indiv, gen_unit);
    }
    

    

    int cptdead = 0;
    int cptalive = 0;
    
    Protein* protein_atest = find_protein_at_startpos(indiv, gen_unit, startpos_afterdup1);
    Protein* protein_btest = find_protein_at_startpos(indiv, gen_unit, startpos_afterdup2);
    
    if (!protein_atest || !protein_btest)
    {
    	cout<<"protein not found after tandem dup: prot a : " << (protein_atest ? "found" : "not found") << " prot b : " << (protein_btest ? "found" : "not found")<<endl;
    	cout<<"this gene will be skipped"<<endl;
    	
    	int gg;
    	cin >> gg;
    	
    	if (gg == 0)
    		    	print_dna_sequence(indiv, gen_unit, 0, gen_unit.dna()->length());
    	if (gg == 1)
    	{
    		print_protein_locations(indiv, gen_unit);
    		cout<<"rnastart="<<posstart_rna<<"   startpos="<<startpos<<"   "<<"   pos_insertion="<<pos_insertion<<"   startpos2 = "<<startpos_afterdup2<<endl;
    	}
    }
    else
    {
		//we'll do switches so these will be invalidated
	    protein_atest = nullptr;
	    protein_btest = nullptr;
    	
	
	    
	    //try to mutate every position that encodes the protein, we skip the shine_dal and spacer since they don't encode anything and kill the protein.  
	    //We also skip the stop codon, hence the -3
	    for (int offset1 = SHINE_DAL_SIZE + SHINE_START_SPACER - 1; offset1 <= protein_length - 3; ++offset1)
	    {
			
	    	//offset1 refers to position in original gene copy
	    	int pos_to_flip1 = startpos_afterdup1 + offset1;
	    	
	    	if (!is_leading_strand)
	    	{
	    		pos_to_flip1 = startpos_afterdup1 - offset1;
	    	}
	    	
	    	
	    	output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip1);
		if (!output_mut_bool)
		{
		    cout << "ERR, output_mut_bool = false" << endl;
		}

		indiv->Reevaluate();
		indiv->compute_statistical_data();
		indiv->compute_non_coding();
			
		indiv->compute_fitness(indiv->phenotypic_target());
		long double fitness_flip_a = indiv->fitness();
		long double metabolic_error_flip_a = indiv->dist_to_target_by_feature(METABOLISM);
	    	

		//we only try unordered pairs of mutations, hence offset2 starts at offset1
		for (int offset2 = offset1; offset2 <= protein_length - 3; ++offset2)
		{
			int pos_to_flip2 = startpos_afterdup2 + offset2;
	    	
		    	if (!is_leading_strand)
		    	{
		    		pos_to_flip2 = startpos_afterdup2 - offset2;
		    	}
		    	
		    	
		    	output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip2);

			if (!output_mut_bool)
			{
			    cout << "ERR, output_mut_bool = false" << endl;
			}

			indiv->Reevaluate();
			indiv->compute_statistical_data();
			indiv->compute_non_coding();
				
			indiv->compute_fitness(indiv->phenotypic_target());
			long double fitness_flip_ab = indiv->fitness();
		    	long double metabolic_error_flip_ab = indiv->dist_to_target_by_feature(METABOLISM);
		    	
		    	
			//revert the a flip to get the b only fitness
			output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip1);
			indiv->Reevaluate();
			indiv->compute_statistical_data();
			indiv->compute_non_coding();
				
			indiv->compute_fitness(indiv->phenotypic_target());
			long double fitness_flip_b = indiv->fitness();
			long double metabolic_error_flip_b = indiv->dist_to_target_by_feature(METABOLISM);
			
			//replay the a flip for future iterations on offset2
			output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip1);
			indiv->Reevaluate();
			indiv->compute_statistical_data();
			indiv->compute_non_coding();
			
			
			
			//at this point, we've flipped two bits --> check the fate
		    	Protein* protein_a = find_protein_at_startpos(indiv, gen_unit, startpos_afterdup1);
		    	Protein* protein_b = find_protein_at_startpos(indiv, gen_unit, startpos_afterdup2);
		
			string deadnote = "";
			if (!protein_a)
			{
				//cout<<"a is dead"<<endl;
				deadnote += "a";
			}
			if (!protein_b)
			{
				//cout<<"b is dead"<<endl;
				deadnote += "b";
			}
			
			if (protein_a && protein_b)
			{
				cptalive++;	
				deadnote = "ab_alive";	
			}
			else 
			{
				cptdead++;
				deadnote += "_dead";
			}
			
			
			fate_list_file << startpos <<", " << neutral_pos <<", "<< offset1<<", "<<offset2<<", "<<deadnote<<", ";
			output_ab_data(orig_fitness, fitness_after_dup, fitness_flip_a, fitness_flip_b, fitness_flip_ab, origm, origh, origw, protein_a, protein_b, fate_list_file, 
					orig_metabolic_error, metabolic_error_after_dup, metabolic_error_flip_a, metabolic_error_flip_b, metabolic_error_flip_ab);
			
			
			/*if (startpos == 5559 && offset1 == 22 && offset2 == 28)
			{
				cout<<"SEQ AFTER"<<endl;
				print_dna_sequence(indiv, gen_unit, 0, gen_unit.dna()->length());
				
				double mutant_metabolic_error = indiv->dist_to_target_by_feature(METABOLISM);	
				
				cout<<"error before="<<orig_metabolic_error<<" error after="<<mutant_metabolic_error<<endl;
				int xxx; cin>>xxx;
			}*/
			
			
			//revert the b flip
			output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip2);
			indiv->Reevaluate();
			indiv->compute_statistical_data();
			indiv->compute_non_coding();
			
			//cout<<"St,st1, st2,plen,lead,o1,pos1,o2,pos2,...."<<endl;
			//cout<<startpos<<", "<<startpos_afterdup1<<", "<<startpos_afterdup2<<", "<<protein_length<<", "<<is_leading_strand<<", "<<offset1<<", "<<pos_to_flip1<<", "<<offset2<<", "<<pos_to_flip2<<", ";
			//cout<<orig_fitness<<", "<<fitness_after_dup<<", "<<fitness_flip_a<<", "<<fitness_flip_b<<", "<<fitness_flip_ab<<endl;
			

		}	//end of offset2 loop
	    	
	    	
		//revert the a flip
		output_mut_bool = gen_unit.dna()->do_switch(pos_to_flip1);
		indiv->Reevaluate();
		indiv->compute_statistical_data();
		indiv->compute_non_coding();
	    	
	    }   //end of offset1 loop */

    }	//end of if a, b found


	cout<<"cptdead="<<cptdead<<"  CPTALIVE="<<cptalive<<endl;
	//int gg;  cin >> gg;

    //REVERT THE DUP BY DELETING SECOND SEGMENT
    gen_unit.dna()->do_deletion(neutral_pos, neutral_pos + to_insert.length());
    
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    //sanity check, make sure genome has reverted to original state
    if (gen_unit.dna()->length() != orig_genome_length)
    {
	cout<<"ERROR : Genome length has changed after reverting everything, before="<<orig_genome_length<<"  after="<<gen_unit.dna()->length()<<endl;
	int xx;
	cin >> xx;
    }
    else
    {
	for (int i = 0; i < gen_unit.dna()->length(); ++i)
	{
		if (orig_genome_seq[i] != gen_unit.dna()->data()[i])
		{
			cout<<"ERROR : Genome differs at position "<<i<<" after reverting everything"<<endl;
			int xx;
			cin >> xx;
		}
    	}
    }
    
    delete [] orig_genome_seq;

        

}













































void compute_fates_of_all_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos, ofstream& fate_list_file)
{
    indiv->Reevaluate();
	
	indiv->compute_fitness(indiv->phenotypic_target());
	double orig_fitness = indiv->fitness();
	cout<<"fitness="<<orig_fitness<<endl;
    /*gen_unit.locate_promoters();
    gen_unit.do_transcription();
    gen_unit.do_translation();
    gen_unit.compute_phenotypic_contribution();*/



    bool is_leading_strand = true;

    char* subseq;

    cout << "*********************************************************************" << endl;
    cout << endl << "ANALYZING PROTEIN AT STARTPOS=" << startpos << endl;
    cout<<"Init fitness="<<orig_fitness<<endl;

    Protein* protein = find_protein_at_startpos(indiv, gen_unit, startpos);

    if (!protein)
    {
        cout << "PROTEIN AT POS=" << startpos << " NOT FOUND" << endl;
        return;
    }

    is_leading_strand = (protein->strand() == LEADING);
    long double origm = protein->mean();
    long double origw = protein->width();
    long double origh = protein->height() * protein->concentration();
    int32_t shine_dal_pos_orig = protein->shine_dal_pos();
    int32_t first_translated_pos_orig = protein->first_translated_pos();
    int32_t last_translated_pos_orig = protein->last_translated_pos();
    int32_t nbcodons = protein->length();

    int32_t shine_dal_pos_a;
    int32_t first_translated_pos_a;
    int32_t last_translated_pos_a;
    Strand strand_a;

    int32_t shine_dal_pos_b;
    int32_t first_translated_pos_b;
    int32_t last_translated_pos_b;
    Strand strand_b;

    cout << "PROM_SEQ=" << PROM_SEQ << "    PROM_MAX_DIFF=" << (int)PROM_MAX_DIFF << "SHINE_DAL_SEQ=" << SHINE_DAL_SEQ << endl << endl;
    //cout<<"TERM_SIZE="<<(int)TERM_SIZE<<"   TERM_LOOP_SIZE="<<(int)TERM_LOOP_SIZE<<"   TERM_STEM_SIZE="<<(int)TERM_STEM_SIZE<<endl;


    cout << "strand=" << (is_leading_strand ? "leading" : "lagging") << endl;

    cout << "shine_dal_pos=" << shine_dal_pos_orig << "   first_translated_pos=" << first_translated_pos_orig << "   last_translated_pos=" << last_translated_pos_orig << endl;
    cout << "m=" << origm << ", ";
    cout << "w=" << origw << ", ";
    cout << "h=" << origh << endl;

    long double g[3] = { origm, origh, origw };
    long double a[3];
    long double b[3];
    bool a_dead, b_dead;

    cout << endl;

    if (SWITCH_DEBUG_SHOW_SEQUENCES)
    {
        const char* seq_orig = gen_unit.dna()->data();
        for (Rna* rna : protein->rna_list())
        {

            int32_t promoter_pos = rna->promoter_pos();
            int32_t transcription_start = rna->first_transcribed_pos();
            int32_t transcription_end = rna->last_transcribed_pos();
            int32_t rna_length = transcription_end - promoter_pos + 1;

            cout << "RNA   strand=" << (rna->strand() == LEADING ? "leading" : "lagging") << " (val=" << rna->strand() << ")  promoter_pos=" << promoter_pos;
            cout << "   transcript_length=" << rna->transcript_length() << "   transcription_start=" << transcription_start << "   transcription_end=" << transcription_end << endl;



            if (is_leading_strand)
            {
                for (int32_t pos = promoter_pos; pos <= transcription_end; ++pos)
                {
                    if (pos == promoter_pos + PROM_SIZE)
                        cout << "|";
                    if (pos == protein->shine_dal_pos())
                        cout << " [";
                    if (pos == protein->last_translated_pos() + 1)
                        cout << "] ";
                    if (pos == transcription_end - TERM_SIZE + 1)
                        cout << "|";

                    cout << seq_orig[pos];
                }
            }
            else
            {
                for (int32_t pos = promoter_pos; pos >= transcription_end; --pos)
                {
                    if (pos == promoter_pos - PROM_SIZE)
                        cout << "|";
                    if (pos == protein->shine_dal_pos())
                        cout << " [";
                    if (pos == protein->last_translated_pos() - 1)
                        cout << "] ";
                    if (pos == transcription_end + TERM_SIZE - 1)
                        cout << "|";

                    cout << revcomp(seq_orig[pos]);
                }
            }
            cout << endl;

        }

        cout << "PROTEIN SEQUENCE";
        if (!is_leading_strand)
            cout << "  (reverse complement)";
        cout << endl;
        output_protein_sequence(indiv, gen_unit, shine_dal_pos_orig, first_translated_pos_orig, last_translated_pos_orig, (is_leading_strand ? LEADING : LAGGING), true);
        cout << endl << endl;
    }


    //we will do switches, which will delete the current protein object.  To avoid problems, we set it to nullptr
    protein = nullptr;
	
    long double fitness_a, fitness_b;

    for (int offset = 0; offset <= nbcodons * (int)(CODON_SIZE); ++offset)
    {
        //go backwards if lagging strand
        int pos_to_affect = startpos + offset * (is_leading_strand ? 1 : -1);


        bool output_mut_bool;
        output_mut_bool = gen_unit.dna()->do_switch(pos_to_affect);


        
        if (!output_mut_bool)
        {
            cout << "ERR, output_mut_bool = false" << endl;
        }

        indiv->Reevaluate();
        indiv->compute_statistical_data();
        indiv->compute_non_coding();


		
	indiv->compute_fitness(indiv->phenotypic_target());
	fitness_a = indiv->fitness();
		


        /*gen_unit.locate_promoters();
        gen_unit.do_transcription();
        gen_unit.do_translation();
        gen_unit.compute_phenotypic_contribution();*/


        //cout << "a " << pos_to_affect << ", " << endl;

        Protein* protein_after = find_protein_at_startpos(indiv, gen_unit, startpos);

        if (protein_after)
        {
            std::cout << std::fixed;
            std::cout << std::setprecision(10);

            /*cout << protein_after->mean() << ", " << protein_after->width() << ", " << protein_after->height();
            cout << ", \t\t" <<
                protein_after->mean() - origm << ", " <<
                protein_after->width() - origw << ", " <<
                protein_after->height() - origh;

            cout << ",\t\t";*/

            //long double g[3] = { origm, origh, origw };
            a[0] = protein_after->mean();
            a[1] = protein_after->height()* protein_after->concentration();
            a[2] = protein_after->width();
        
            //long double b[3] = { protein_after->mean(), protein_after->height(), protein_after->width() };
            //vector<long double> params = ProbabilitiesCalculation(g, a, b);

            //long double ig_a = min(1, params[5]);
            //long double ia_g = min(1, params[7]);
            //long double pa = min(1, params[3]);



            std::cout << std::fixed;
            std::cout << std::setprecision(4);

            //cout << ig_a << ", " << ia_g << ", " << pa;

            //cout << endl;

            if (SWITCH_DEBUG_SHOW_SEQUENCES)
            {   
                cout << "NOT DEAD" << endl;
                shine_dal_pos_a = protein_after->shine_dal_pos();
                first_translated_pos_a = protein_after->first_translated_pos();
                last_translated_pos_a = protein_after->last_translated_pos();
                strand_a = protein_after->strand();
                output_protein_sequence(indiv, gen_unit, protein_after->shine_dal_pos(), protein_after->first_translated_pos(), protein_after->last_translated_pos(), protein_after->strand(), true);
                cout << endl << endl;
            }
            a_dead = false;
        }
        else
        {
            //cout << "DEAD" << endl;
            if (SWITCH_DEBUG_SHOW_SEQUENCES)
            {
                output_protein_sequence(indiv, gen_unit, shine_dal_pos_orig, first_translated_pos_orig, last_translated_pos_orig, (is_leading_strand ? LEADING : LAGGING), true);
                cout << endl << endl;
            }
            a_dead = true;
        }



        protein_after = nullptr;



        //REVERT THE SWITCH
        output_mut_bool = gen_unit.dna()->do_switch(pos_to_affect);

        if (!output_mut_bool)
        {
            cout << "ERR, output_mut_bool = false" << endl;
        }

        indiv->Reevaluate();
        indiv->compute_statistical_data();
        indiv->compute_non_coding();

        /*gen_unit.locate_promoters();
        gen_unit.do_transcription();
        gen_unit.do_translation();
        gen_unit.compute_phenotypic_contribution();*/

        Protein* protein_after_revert = find_protein_at_startpos(indiv, gen_unit, startpos);
        if (!protein_after_revert)
            cout << "ERROR: reverted the switch but gene is dead, NOT SUPPOSED TO HAPPEN" << endl;
        protein_after_revert = nullptr;

        //////////////////////////////////////////////////////////////////////
        if (!a_dead) {
            int temp_offset = offset;
            if (is_leading_strand) {
                temp_offset++;
            }
            else {
                temp_offset--;
            }

            for (int offset2 = offset; offset2 <= nbcodons * (int)(CODON_SIZE)-1; ++offset2)
            {
                //go backwards if lagging strand
                int pos_to_affect1 = startpos + offset2 * (is_leading_strand ? 1 : -1);
                if (is_leading_strand) {
                    pos_to_affect1++;
                }
                else {
                    pos_to_affect1--;
                }

                bool output_mut_bool1;
                output_mut_bool1 = gen_unit.dna()->do_switch(pos_to_affect1);

                if (!output_mut_bool1)
                {
                    cout << "ERR, output_mut_bool = false" << endl;
                }

                indiv->Reevaluate();
                indiv->compute_statistical_data();
                indiv->compute_non_coding();

                /*gen_unit.locate_promoters();
                gen_unit.do_transcription();
                gen_unit.do_translation();
                gen_unit.compute_phenotypic_contribution();*/


                //cout << "b " << pos_to_affect1 << ", " << endl;

                protein_after = find_protein_at_startpos(indiv, gen_unit, startpos);

                if (protein_after)
                {
                    std::cout << std::fixed;
                    std::cout << std::setprecision(10);

                    /*cout << protein_after->mean() << ", " << protein_after->width() << ", " << protein_after->height();
                    cout << ", \t\t" <<
                        protein_after->mean() - origm << ", " <<
                        protein_after->width() - origw << ", " <<
                        protein_after->height() - origh;

                    cout << ",\t\t";*/

                    //long double g[3] = { origm, origh, origw };
                    //long double a[3] = { protein_after->mean(), protein_after->height(), protein_after->width() };
                    b[0] = protein_after->mean();
                    b[1] = protein_after->height()* protein_after->concentration();
                    b[2] = protein_after->width();
                    //vector<long double> params = ProbabilitiesCalculation(g, a, b);

                    //long double ig_a = min(1, params[5]);
                    //long double ia_g = min(1, params[7]);
                    //long double pa = min(1, params[3]);



                    std::cout << std::fixed;
                    std::cout << std::setprecision(4);

                    //cout << ig_a << ", " << ia_g << ", " << pa;

                    //cout << endl;

                    if (SWITCH_DEBUG_SHOW_SEQUENCES)
                    {
                        cout << "NOT DEAD" << endl;
                        shine_dal_pos_b = protein_after->shine_dal_pos();
                        first_translated_pos_b = protein_after->first_translated_pos();
                        last_translated_pos_b = protein_after->last_translated_pos();
                        strand_b = protein_after->strand();
                        output_protein_sequence(indiv, gen_unit, protein_after->shine_dal_pos(), protein_after->first_translated_pos(), protein_after->last_translated_pos(), protein_after->strand(), true);
                        cout << endl << endl;
                    }
                    b_dead = false;
                }
                else
                {
                    //cout << "DEAD" << endl;
                    if (SWITCH_DEBUG_SHOW_SEQUENCES)
                    {
                        output_protein_sequence(indiv, gen_unit, shine_dal_pos_orig, first_translated_pos_orig, last_translated_pos_orig, (is_leading_strand ? LEADING : LAGGING), true);
                        cout << endl << endl;
                    }
                    b_dead = true;
                }



                protein_after = nullptr;



                //REVERT THE SWITCH
                output_mut_bool1 = gen_unit.dna()->do_switch(pos_to_affect1);

                if (!output_mut_bool1)
                {
                    cout << "ERR, output_mut_bool = false" << endl;
                }

                indiv->Reevaluate();
                indiv->compute_statistical_data();
                indiv->compute_non_coding();
				
				
				indiv->compute_fitness(indiv->phenotypic_target());
				fitness_b = indiv->fitness();

                /*gen_unit.locate_promoters();
                gen_unit.do_transcription();
                gen_unit.do_translation();
                gen_unit.compute_phenotypic_contribution();*/

                protein_after_revert = find_protein_at_startpos(indiv, gen_unit, startpos);
                if (!protein_after_revert)
                    cout << "ERROR: reverted the switch but gene is dead, NOT SUPPOSED TO HAPPEN" << endl;


                protein_after_revert = nullptr;



                std::cout << std::fixed;
                std::cout << std::setprecision(10);

                /*cout << "mg: " << g[0] << ", hg: " << g[1] << ", wg: " << g[2] << endl;
                if (a_dead) {
                    cout << "a is DEAD!" << endl;
                }
                else {
                    cout << "ma: " << a[0] << ", ha: " << a[1] << ", wa: " << a[2] << endl;
                }
                if (b_dead) {
                    cout << "b is DEAD!" << endl;
                }
                else {
                    cout << "mb: " << b[0] << ", hb: " << b[1] << ", wb: " << b[2] << endl;
                }*/

                if (!a_dead && !b_dead) {
                    vector<long double> params = ProbabilitiesCalculation(g, a, b);

                    long double subfunc = params[0];
                    long double neofunc = params[2];
                    long double cons = params[1];
                    long double pseudo = params[3];
                    long double spec = params[4];
					long double cons_loss = params[17];
					long double neo_loss = params[18];

					long double ig_a = params[5];
					long double ig_b = params[6];
					long double ia_g = params[7];
					long double ib_g = params[8];

					long double ig_a_plus_b = params[9];
					long double ia_plus_b_g = params[10];
					long double pa = params[11];
					long double pb = params[12];
					
	
					
					fate_list_file << std::fixed;
					fate_list_file << std::setprecision(4);
                    fate_list_file << subfunc 
                        << "," << neofunc
                        << "," << cons
                        << "," << pseudo
                        << "," << spec
						<< "," << cons_loss 
						<< "," << neo_loss
                        << "," << g[0] << "," << g[1] << "," << g[2] 
                        << "," << a[0] << "," << a[1] << "," << a[2]
                        << "," << b[0] << "," << b[1] << "," << b[2];
						
						

						fate_list_file << std::fixed;
						fate_list_file << std::setprecision(10);
						
						long double tenpower = 1000000000000000;
						
						fate_list_file<< ", " << orig_fitness * tenpower<<", "<< fitness_a * tenpower <<", "<<fitness_b * tenpower;
						
						fate_list_file << ", "<<ig_a<<", "<<ia_g<<", "<<ig_b<<", "<<ib_g<<", "<<ig_a_plus_b<<", "<<ia_plus_b_g<<", "<<pa<<", "<<pb;
						
						fate_list_file<<std::endl;
						
						
						
					if (subfunc > 0.90)
					{			
						
						cout << "subfunc, neofunc, cons, pseudo, spec, mg, hg, wg, ma, ha, wa, mb, hb, wb , fitness_g, fitness_a, fitness_b       " << endl;
						std::cout << std::fixed;
						std::cout << std::setprecision(4);
						cout << subfunc << ", " << neofunc << ", " << cons << ", " << pseudo << ", " << spec << ", " << g[0] << ", " << g[1] << ", " << g[2] << 
																																							  ", " << a[0] << ", " << a[1] << ", " << a[2] <<
																																							  ", " << b[0] << ", " << b[1] << ", " << b[2];
						cout << std::fixed;
						cout << std::setprecision(10);
						
						cout<< ", " << orig_fitness * tenpower <<", "<< fitness_a * tenpower <<", "<<fitness_b * tenpower;
						cout<<endl;
					}
                    /*if (subfunc >= 0.0001)
                        cout << "subfun observed!" << endl;
                    if (subfunc > neofunc && subfunc > cons && subfunc > pseudo && subfunc > spec)
                        cout << "fate is subfun!" << endl;*/
                }
                else {
                    //cout << "a or b is dead, NO fate probability!" << endl << endl;
                }
            }
        }
        else {
            //cout << "a is dead, So b is not important, NO fate probability!" << endl << endl;
        }
        ///////////////////////////////////////////////////////////////////////////////

    }

}




void try_all_switches(Individual* indiv, GeneticUnit& gen_unit, int startpos)
{
    indiv->Reevaluate();
    /*gen_unit.locate_promoters();
    gen_unit.do_transcription();
    gen_unit.do_translation();
    gen_unit.compute_phenotypic_contribution();*/


  bool is_leading_strand = true;
  
  char* subseq;
  
  
  
  
  cout<<"*********************************************************************"<<endl;
  cout<<endl<<"ANALYZING PROTEIN AT STARTPOS="<<startpos<<endl;
  

  Protein* protein = find_protein_at_startpos(indiv, gen_unit, startpos);

  if (!protein)
  {
  	cout<<"PROTEIN AT POS="<<startpos<<" NOT FOUND"<<endl;
  	return;
  }
  
  is_leading_strand = (protein->strand() == LEADING);
  long double origm = protein->mean();
  long double origw = protein->width();
  long double origh = protein->height() * protein->concentration();
  int32_t shine_dal_pos_orig = protein->shine_dal_pos();
  int32_t first_translated_pos_orig = protein->first_translated_pos();
  int32_t last_translated_pos_orig = protein->last_translated_pos();
  int32_t nbcodons = protein->length();
  
  
  cout<<"PROM_SEQ="<<PROM_SEQ<<"    PROM_MAX_DIFF="<<(int)PROM_MAX_DIFF<<"SHINE_DAL_SEQ="<<SHINE_DAL_SEQ<<endl<<endl;
  //cout<<"TERM_SIZE="<<(int)TERM_SIZE<<"   TERM_LOOP_SIZE="<<(int)TERM_LOOP_SIZE<<"   TERM_STEM_SIZE="<<(int)TERM_STEM_SIZE<<endl;


  cout<<"startpos="<<(is_leading_strand ? "leading" : "lagging")<<endl;
  cout<<"shine_dal_pos="<<shine_dal_pos_orig<<"   first_translated_pos="<<first_translated_pos_orig<<"   last_translated_pos="<<last_translated_pos_orig<<endl;
  cout<<"m="<<origm<<", ";
  cout<<"w="<<origw<<", ";
  cout<<"h="<<origh<<endl;

  
  cout<<endl;
   
  if (SWITCH_DEBUG_SHOW_SEQUENCES)
  {
	  const char* seq_orig = gen_unit.dna()->data();
	  for (Rna* rna : protein->rna_list())
	  {

	  	int32_t promoter_pos = rna->promoter_pos();
	  	int32_t transcription_start = rna->first_transcribed_pos();
	  	int32_t transcription_end = rna->last_transcribed_pos();
	  	int32_t rna_length = transcription_end - promoter_pos + 1;

		cout<<"RNA   strand="<<(rna->strand() == LEADING ? "leading":"lagging") << " (val="<<rna->strand()<<")  promoter_pos="<<promoter_pos;
		cout<<"   transcript_length="<<rna->transcript_length()<<"   transcription_start="<<transcription_start<<"   transcription_end="<<transcription_end<<endl;

	  	
	  	
	  	if (is_leading_strand)
	  	{
	  		for (int32_t pos = promoter_pos; pos <= transcription_end; ++pos)
	  		{
	  			if (pos == promoter_pos + PROM_SIZE)
	  				cout<<"|";
	  			if (pos == protein->shine_dal_pos())
	  				cout<<" [";
	  			if (pos == protein->last_translated_pos() + 1)
	  				cout<<"] ";
	  			if (pos == transcription_end - TERM_SIZE + 1)
	  				cout<<"|";
	  				
	  			cout<<seq_orig[pos];
	  		}	
	  	}
	  	else 
	  	{
	  		for (int32_t pos = promoter_pos; pos >= transcription_end; --pos)
	  		{
	  			if (pos == promoter_pos - PROM_SIZE)
	  				cout<<"|";
	  			if (pos == protein->shine_dal_pos())
	  				cout<<" [";
	  			if (pos == protein->last_translated_pos() - 1)
	  				cout<<"] ";
	  			if (pos == transcription_end + TERM_SIZE - 1)
	  				cout<<"|";
	  				
	  			cout<<revcomp(seq_orig[pos]);
	  		}	
	  	}
	  	cout<<endl;

	  }
	  
	  cout<<"PROTEIN SEQUENCE";
	  if (!is_leading_strand)
	  	cout<<"  (reverse complement)";
	  cout<<endl;
	  output_protein_sequence(indiv, gen_unit, shine_dal_pos_orig, first_translated_pos_orig, last_translated_pos_orig, (is_leading_strand ? LEADING : LAGGING), true);
	  cout<<endl<<endl;
  }
	  
  
  cout<<"pos_altered,m_after,w_after,h_after,";
  
  cout<<"\t\t delta_m,      delta_w,      delta_h,      \t\t ig_a,ia_g,pseudo_a"<<endl;
  
  //subseq = gen_unit.dna()->subsequence(start_pos, start_pos + nbcodons * CODON_SIZE * (is_leading_strand ? 1 : -1), (is_leading_strand ? LEADING : LAGGING));
  //delete [] subseq;
  
  
  //we will do switches, which will delete the current protein object.  To avoid problems, we set it to nullptr
  protein = nullptr;

  for (int offset = 0; offset <= nbcodons * (int)(CODON_SIZE); ++offset)
  {
    //go backwards if lagging strand
    int pos_to_affect = startpos + offset * (is_leading_strand ? 1 : -1);
    

    bool output_mut_bool;
    output_mut_bool = gen_unit.dna()->do_switch(pos_to_affect);
    
    if (!output_mut_bool)
    {
      cout<<"ERR, output_mut_bool = false"<<endl;
    }
    
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();
 
    /*gen_unit.locate_promoters();
    gen_unit.do_transcription();
    gen_unit.do_translation();
    gen_unit.compute_phenotypic_contribution();*/
    
    
    cout<<pos_to_affect<<", ";
    
    Protein* protein_after = find_protein_at_startpos(indiv, gen_unit, startpos);
    
    if (protein_after)
    {
        std::cout << std::fixed;
	std::cout << std::setprecision(10);
    
    	cout<<protein_after->mean()<<", "<<protein_after->width()<<", "<<protein_after->height() * protein->concentration();
    	cout<<", \t\t"<<
    		protein_after->mean() - origm<<", "<<
    		protein_after->width() - origw<<", "<<
    		protein_after->height() * protein->concentration() - origh;
    	
    	cout<<",\t\t";
    	
	long double g[3] = {origm, origh, origw};
	long double a[3] = {protein_after->mean(), protein_after->height() * protein->concentration(), protein_after->width()};
	long double b[3] = {protein_after->mean(), protein_after->height() * protein->concentration(), protein_after->width()};
	vector<long double> params = ProbabilitiesCalculation(g, a, b);
        
        long double ig_a = min(1, params[5]);
        long double ia_g = min(1, params[7]);
        long double pa = min(1, params[3]);
        
        
  
  	std::cout << std::fixed;
	std::cout << std::setprecision(4);
	
	cout<<ig_a<<", "<<ia_g<<", "<<pa;      
    	
    	cout<<endl;
    	
    	if (SWITCH_DEBUG_SHOW_SEQUENCES)
    	{
    		output_protein_sequence(indiv, gen_unit, protein_after->shine_dal_pos(), protein_after->first_translated_pos(), protein_after->last_translated_pos(), protein_after->strand(), true);
    		cout<<endl<<endl;
    	}

    }
    else 
    {
    	//cout<<"DEAD"<<endl;
    	if (SWITCH_DEBUG_SHOW_SEQUENCES)
    	{
    		output_protein_sequence(indiv, gen_unit, shine_dal_pos_orig, first_translated_pos_orig, last_translated_pos_orig, (is_leading_strand ? LEADING : LAGGING), true);
    		cout<<endl<<endl;
    	}
    
    }
    
    
      
    protein_after = nullptr;
    
    
    
    //REVERT THE SWITCH
    output_mut_bool = gen_unit.dna()->do_switch(pos_to_affect);
    
    if (!output_mut_bool)
    {
      cout<<"ERR, output_mut_bool = false"<<endl;
    }
    
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    /*gen_unit.locate_promoters();
    gen_unit.do_transcription();
    gen_unit.do_translation();
    gen_unit.compute_phenotypic_contribution();*/
    
    Protein* protein_after_revert = find_protein_at_startpos(indiv, gen_unit, startpos);
    if (!protein_after_revert)
      cout<<"ERROR: reverted the switch but gene is dead, NOT SUPPOSED TO HAPPEN"<<endl;
    
    
    protein_after_revert = nullptr;
  }
  
  
}


















void interpret_cmd_line_options(int argc, char* argv[]) {
  // =====================
  //  Parse command line
  // =====================
  const char * short_options = "hVb:e:v:fdg:s";
  static struct option long_options[] = {
    {"help",                no_argument, NULL, 'h'},
    {"version",             no_argument, NULL, 'V'},
    {"verbose",             no_argument, NULL, 'v'},
    {"begin",     required_argument, nullptr, 'b'},
    {"end",       required_argument, nullptr, 'e'},
    {"activate fates calculation",       no_argument, NULL, 'f'},
    {"activate switch debugging",       no_argument, NULL, 'd'},
    {"switch debug generation",       required_argument, nullptr, 'g'},
    {"switch debug show sequences",       no_argument, nullptr, 's'},
    {0, 0, 0, 0}
  };
  
  int option;
  while((option = getopt_long(argc, argv, short_options,
                              long_options, nullptr)) != -1) {
    switch(option) {
    case 'h':
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V':
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'v':
      verbose = true;
      break;
    case 'b':
      pt_begin = atol(optarg);
      break;
    case 'e':
      pt_end = atol(optarg);
      break;
      //ML EDIT
    case 'f':
      FATE_PROB_MODE = true;
      break;
    case 'd':
      SWITCH_DEBUG_MODE = true;
      break;
    case 'g':
      SWITCH_DEBUG_GENERATION = atoi(optarg);
      break;
    case 's':
      SWITCH_DEBUG_SHOW_SEQUENCES = true;
      break;
      //END ML EDIT
    default:
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }
  
  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }
  
  lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);
}

void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
    prog_name = prog_path;
  }
  
  printf("******************************************************************************\n");
  printf("*                                                                            *\n");
  printf("*                        aevol - Artificial Evolution                        *\n");
  printf("*                                                                            *\n");
  printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
  printf("*                                                                            *\n");
  printf("******************************************************************************\n");
  printf("\n");
  printf("%s: create an experiment with setup as specified in PARAM_FILE.\n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s LINEAGE_FILE [-FMv]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -v, --verbose\n\tbe verbose\n");
}

ProteinMap* compute_protein_map(Individual* indiv) {
  ProteinMap* pmap = new ProteinMap(indiv->protein_list().size(),indiv->amount_of_dna());
  
  // Make a copy of each genetic unit's protein list
  for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
    // append all proteins from `gen_unit` to `protein_list_`
    for (auto& strand_id: {LEADING, LAGGING}) {
      auto& strand = gen_unit.protein_list(strand_id);
      int pos_next = std::prev(gen_unit.protein_list(strand_id).end())->shine_dal_pos();
      bool first = true;
      int pos_first = gen_unit.protein_list(strand_id).begin()->shine_dal_pos();
      int pos_prev = -1;
      for (auto& p: strand) {
        int dist = -1;
        if (first) {
          if (strand_id == LEADING)
            dist = p.shine_dal_pos() + (indiv->amount_of_dna() - pos_next);
          else
            dist = (indiv->amount_of_dna() - p.shine_dal_pos()) + pos_next;
          
          first = false;
        } else if (p.shine_dal_pos() == std::prev(gen_unit.protein_list(strand_id).end())->shine_dal_pos()) {
          if (strand_id == LEADING)
            dist = (indiv->amount_of_dna() - p.shine_dal_pos()) + pos_first;
          else
            dist = p.shine_dal_pos() + (indiv->amount_of_dna() - pos_first);
        } else {
          if (strand_id == LEADING)
            dist = p.shine_dal_pos() - pos_prev;
          else
            dist = pos_prev - p.shine_dal_pos();
        }
        pos_prev = p.shine_dal_pos();
        int8_t prom_dist;
        gen_unit.is_promoter(LEADING, (*p.rna_list().begin())->promoter_pos(),
                             prom_dist);

        pmap->add_protein(p.shine_dal_pos(),p.length(),p.concentration(),prom_dist,dist);
      }
    }
  }
  
  return pmap;
}

//******************************************************************************************
//*

Protein_List* check_SWITCH(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos){
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(protein_list->pro_list_[i].start_pos_ == list_prot_temp->start_pos_[j]){
        protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                     list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
        break;
      }
    }
  }
  // remove suddenly deleted proteins
  bool found;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    found = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        found = true;
      }
    }
    if(!found){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  // add suddenly new proteins
  bool mustadd;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    mustadd = true;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        mustadd = false;
      }
    }
    if(mustadd){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_S_INS(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos, int32_t len){
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    if(protein_list->pro_list_[i].start_pos_ > pos){
      for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
        if(protein_list->pro_list_[i].start_pos_ + len == list_prot_temp->start_pos_[j]){
          protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                       list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          break;
        }
      }
    }
    else{
      for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
        if(protein_list->pro_list_[i].start_pos_ == list_prot_temp->start_pos_[j]){
          protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                       list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          break;
        }
      }
    }
  }
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  // add suddenly new proteins
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_S_DEL(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos, int32_t len){
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    if(protein_list->pro_list_[i].start_pos_ > pos){
      for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
        if(protein_list->pro_list_[i].start_pos_ - len == list_prot_temp->start_pos_[j]){
          protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                       list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          break;
        }
      }
    }
    else if(protein_list->pro_list_[i].start_pos_ < pos){
      for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
        if(protein_list->pro_list_[i].start_pos_ == list_prot_temp->start_pos_[j]){
          protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                       list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          break;
        }
      }
    }
  }
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  // add suddenly new proteins
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_DUP(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2,  int32_t pos3){
  int32_t seg_length;
  vector<int> js;
  vector<int> jss;
  if (pos1 < pos2){
    seg_length = pos2 - pos1;
    for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
      if(protein_list->pro_list_[i].start_pos_ > pos1 && (protein_list->pro_list_[i].start_pos_ + protein_list->pro_list_[i].length_) < pos2){
        //int j = protein_list->pro_list_[i].start_pos_ + pos3;
        for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
          if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - pos1 + pos3){
            //std::cout<< list_prot_temp->start_pos_[j]<< std::endl;
            js.push_back(j);
            jss.push_back(protein_list->pro_list_[i].prot_id_);
          }
        }
        if(protein_list->pro_list_[i].start_pos_ > pos3){
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + seg_length){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          //is match
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
      }
      else{
        if(protein_list->pro_list_[i].start_pos_ > pos3){
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + seg_length){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          //is match
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
      }
    }
    for(int i = 0; i < js.size() ; i++){
      protein_list->insert_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]],
                                   list_prot_temp->hamming_dist_[js[i]], list_prot_temp->dist_next_prot_[js[i]], i_prot, jss[i]);
      i_prot++;
    }
  } 
  else{ // pos1 > pos2
    int32_t tmp1_len = protein_list->dna_length_ - pos1;
    int32_t tmp2_len = pos2;
    seg_length = tmp1_len + tmp2_len;
    //cout<< "right plaaaaaaceee"<<std::endl;
    for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
      if((protein_list->pro_list_[i].start_pos_ > 0 and (protein_list->pro_list_[i].start_pos_ + protein_list->pro_list_[i].length_) < pos2) or (protein_list->pro_list_[i].start_pos_ > pos1 and (protein_list->pro_list_[i].start_pos_ + protein_list->pro_list_[i].length_) < protein_list->dna_length_) ){
        //int j = protein_list->pro_list_[i].start_pos_ + pos3;
        //cout<< "right plaaaaaaceee"<<list_prot_temp->start_pos_[0]<<std::endl;
        for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
          if(protein_list->pro_list_[i].start_pos_ > pos1){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - pos1 + pos3){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
          else{
            if(list_prot_temp->start_pos_[j] == protein_list->dna_length_ - pos1 + pos3 + protein_list->pro_list_[i].start_pos_){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        if(protein_list->pro_list_[i].start_pos_ > pos3){
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + seg_length){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          //is match
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
      }
      else{
        if(protein_list->pro_list_[i].start_pos_ > pos3){
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + seg_length){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          //is match
          for(int j = 0; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                           list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
            }
          }
        }
      }
    }
    for(int i = 0; i < js.size() ; i++){
      protein_list->insert_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]],
                                   list_prot_temp->hamming_dist_[js[i]], list_prot_temp->dist_next_prot_[js[i]], i_prot, jss[i]);
      i_prot++;
    }
  }
  
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  // add suddenly new proteins
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_DEL(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2){
  int len;
  vector<int> js;
  vector<int> jsl;
  if(pos1 < pos2){
    len = pos2 - pos1;
    for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
      if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos2){
        js.push_back(protein_list->pro_list_[i].start_pos_);
        jsl.push_back(protein_list->pro_list_[i].length_);
        //protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      }
      else if (protein_list->pro_list_[i].start_pos_ > pos2){
        for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
          if(protein_list->pro_list_[i].start_pos_ - len == list_prot_temp->start_pos_[j])
            protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                         list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
        }
      }
      else{
        for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
          if(protein_list->pro_list_[i].start_pos_ == list_prot_temp->start_pos_[j])
            protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                         list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
        }
      }
    }
  }
  else{ // pos1 > pos2
    for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
      if(protein_list->pro_list_[i].start_pos_ < pos2){
        js.push_back(protein_list->pro_list_[i].start_pos_);
        jsl.push_back(protein_list->pro_list_[i].length_);
        //protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      }
      else if(protein_list->pro_list_[i].start_pos_ > pos1){
        js.push_back(protein_list->pro_list_[i].start_pos_);
        jsl.push_back(protein_list->pro_list_[i].length_);
        //protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      }
      else{
        for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
          if(protein_list->pro_list_[i].start_pos_ - pos2 == list_prot_temp->start_pos_[j])
            protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j],
                                         list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
        }
      }
    }
  }
  for(int i = 0; i < js.size() ; i++){
    std::cout<< "rem pro in del :"<< js[i]<< std::endl;
    protein_list->remove_protein(js[i], jsl[i]);
  }
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_TRANS(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2, int32_t pos3, int32_t pos4, bool invert){
  //    1) Note that in Case 1 (without inversion), whichever position
  //       comes first, translocating segment [pos_1->pos_2] to pos_4
  //       through pos_3 is always equivalent to rearrange the sequences from
  //       an ABCDE order to ADCBE
  //    2) In Case 2, depending on which position comes first, we may have the
  //       following rearrangements :
  //       (ABCDE => ADB'C'E) or (ABCDE => AC'D'BE)
  //       where X' stands for "inverted X"
  //
  //  Case 1 : Without inversion
  //
  //         A      B        C       D       E
  //      |----->=======[>=======>-------[>-------|
  //          pos_1   pos_3    pos_2   pos_4
  //                         |
  //                         V
  //         A      D        C       B        E
  //      |----->-------[>=======>=======[>-------|
  //
  //
  //         A      B        C       D       E
  //      |=====>-------[>------->=======[>=======|
  //          pos_2   pos_4    pos_1   pos_3
  //                         |
  //                         V
  //         A      D        C       B        E
  //      |=====>=======[>------->-------[>=======|
  //
  //
  //         A      B        C       D       E
  //      |====[>========>-------[>------->=======|
  //          pos_3    pos_2    pos_4   pos_1
  //                         |
  //                         V
  //         A       D       C        B        E
  //      |=====[>------->-------[>=======[>=======|
  //
  //
  //         A      B        C       D       E
  //      |----[>-------->=======[>=======>-------|
  //          pos_4    pos_1    pos_3   pos_2
  //                         |
  //                         V
  //         A       D       C        B        E
  //      |-----[>=======>=======[>-------[>-------|
  //
  //  Case 2 : With inversion
  //
  //    Case 2.A
  //
  //         A      B        C       D        E
  //      |----->=======[>=======>-------<]-------|
  //          pos_1   pos_3    pos_2   pos_4
  //                         |
  //                         V
  //         A      D        B'      C'       E
  //      |----->-------<]=======<=======<]-------|
  //
  //
  //         A      B        C       D       E
  //      |=====>-------[>------->=======<]=======|
  //          pos_2   pos_4    pos_1   pos_3
  //                         |
  //                         V
  //         A      D        B'      C'       E
  //      |=====>=======<]-------<-------<]=======|
  //
  //    Case 2.B
  //
  //         A      B        C       D       E
  //      |====[>========>-------<]------->=======|
  //          pos_3    pos_2    pos_4   pos_1
  //                         |
  //                         V
  //         A       C'      D'       B       E
  //      |=====[>-------<-------[>=======>=======|
  //
  //
  //         A      B        C       D       E
  //      |----<]-------->=======[>=======>-------|
  //          pos_4    pos_1    pos_3   pos_2
  //                         |
  //                         V
  //         A       C'      D'       B       E
  //      |-----<]=======>=======<]------->-------|
  // Determine which position comes first and do the corresponding rearrangement
  // TODO(dpa) use min from std
  int32_t pos_min = Utils::min(pos1,
                               Utils::min(pos2, Utils::min(pos3, pos4)));
  vector<int> js;
  vector<int> jss;
  int len_A, len_B, len_D, len_C, len_E;
  if (not invert) {
    //         A      B        C       D       E
    //      |----->=======[>=======>-------[>-------|
    //          pos_1   pos_3    pos_2   pos_4
    if (pos_min == pos1) {
      //ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
      len_B = pos3 - pos1;
      len_D = pos4 - pos2;
      len_C = pos2 - pos3;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_D - len_B){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |=====>-------[>------->=======[>=======|
    //          pos_2   pos_4    pos_1   pos_3
    else if (pos_min == pos2) {
      //ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
      len_B = pos4 - pos2;
      len_D = pos3 - pos1;
      len_C = pos1 - pos4;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_D - len_B){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |====[>========>-------[>------->=======|
    //          pos_3    pos_2    pos_4   pos_1
    else if (pos_min == pos3) {
      //ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
      len_B = pos3 - pos2;
      len_D = pos4 - pos1;
      len_C = pos4 - pos2;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_D - len_B){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |----[>-------->=======[>=======>-------|
    //          pos_4    pos_1    pos_3   pos_2
    else { // if (pos_min == pos4)
      //ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
      len_B = pos1 - pos4;
      len_D = pos2 - pos3;
      len_C = pos3 - pos1;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_D - len_B){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
  }
  
  else { // invert
    //         A      B        C       D        E
    //      |----->=======[>=======>-------<]-------|
    //          pos_1   pos_3    pos_2   pos_4
    //                         |
    //                         V
    //         A      D        B'      C'       E
    //      |----->-------<]=======<=======<]-------|
    if (pos_min == pos1) {
      //ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
      len_A = pos1;
      len_B = pos3 - pos1;
      len_C = pos2 - pos3;
      len_D = pos4 - pos2;
      len_E = protein_list->dna_length_ - pos4;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos3 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos2 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_B + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |=====>-------[>------->=======<]=======|
    //          pos_2   pos_4    pos_1   pos_3
    //                         |
    //                         V
    //         A      D        B'      C'       E
    //      |=====>=======<]-------<-------<]=======|
    else if (pos_min == pos2) {
      //ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
      len_A = pos2;
      len_B = pos4 - pos2;
      len_C = pos1 - pos4;
      len_D = pos3 - pos1;
      len_E = protein_list->dna_length_ - pos3;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos4 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos1 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_B + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ - len_B - len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |====[>========>-------<]------->=======|
    //          pos_3    pos_2    pos_4   pos_1
    //                         |
    //                         V
    //         A       C'      D'       B       E
    //      |=====[>-------<-------[>=======>=======|
    else if (pos_min == pos3) {
      //ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
      len_A = pos3;
      len_B = pos2 - pos3;
      len_C = pos4 - pos2;
      len_D = pos1 - pos4;
      len_E = protein_list->dna_length_ - pos1;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos2 and protein_list->pro_list_[i].start_pos_ < pos4){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos4 - protein_list->pro_list_[i].start_pos_ - 1 + len_A ){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos1 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
    //         A      B        C       D       E
    //      |----<]-------->=======[>=======>-------|
    //          pos_4    pos_1    pos_3   pos_2
    //                         |
    //                         V
    //         A       C'      D'       B       E
    //      |-----<]=======>=======<]------->-------|
    else { // if (pos_min == pos_4)
      //ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
      len_A = pos4;
      len_B = pos1 - pos4;
      len_C = pos3 - pos1;
      len_D = pos2 - pos3;
      len_E = protein_list->dna_length_ - pos2;
      for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
        if(protein_list->pro_list_[i].start_pos_ > pos4 and protein_list->pro_list_[i].start_pos_ < pos1){// B
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_ + len_C + len_D){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos1 and protein_list->pro_list_[i].start_pos_ < pos3){// C
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos3 - protein_list->pro_list_[i].start_pos_ - 1 + len_A ){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else if(protein_list->pro_list_[i].start_pos_ > pos3 and protein_list->pro_list_[i].start_pos_ < pos2){// D
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == pos2 - protein_list->pro_list_[i].start_pos_ - 1 + len_A + len_C){
              js.push_back(j);
              jss.push_back(protein_list->pro_list_[i].prot_id_);
            }
          }
        }
        else{
          for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
            if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_)
              protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], 
                                           list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_,  protein_list->pro_list_[i].prot_id_);
          }
        }
      }
      for(int i = 0 ; i < js.size() ; i++){
        protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]], list_prot_temp->hamming_dist_[js[i]], 
                                     list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
      }
    }
  }
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}

Protein_List* check_INV(Protein_List* protein_list, ProteinMap* list_prot_temp, int32_t pos1, int32_t pos2){
  //if(pos1 == pos2) return protein_list;
  vector<int> js;
  vector<int> jss;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    if(protein_list->pro_list_[i].start_pos_ > pos1 and (protein_list->pro_list_[i].start_pos_ + protein_list->pro_list_[i].length_) < pos2){
      for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
        if(list_prot_temp->start_pos_[j] == pos2 - protein_list->pro_list_[i].start_pos_ + pos1 - 1){
          js.push_back(j);
          jss.push_back(protein_list->pro_list_[i].prot_id_);
          // protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_);
        }
        else if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
          protein_list->modify_protein(list_prot_temp->start_pos_[j], list_prot_temp->length_[j], list_prot_temp->basal_level_[j], 
                                       list_prot_temp->hamming_dist_[j], list_prot_temp->dist_next_prot_[j], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
        }
      }
    }
    else{
      protein_list->modify_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i], 
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], protein_list->pro_list_[i].prot_id_, protein_list->pro_list_[i].prot_id_);
      
    }
  }
  //std::cout<< "num of inv:" << js.size() << std::endl;
  for(int i = 0; i < js.size() ; i++){
    protein_list->modify_protein(list_prot_temp->start_pos_[js[i]], list_prot_temp->length_[js[i]], list_prot_temp->basal_level_[js[i]],
                                 list_prot_temp->hamming_dist_[js[i]], list_prot_temp->dist_next_prot_[js[i]], jss[i], jss[i]);
    
  }
  // remove suddenly deleted proteins
  bool remove;
  for(int i = 0 ; i < protein_list->nb_prot_ ; i++){
    remove = false;
    for(int j = 0 ; j < list_prot_temp->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[j] == protein_list->pro_list_[i].start_pos_){
        remove = true;
      }
    }
    if(!remove){
      protein_list->remove_protein(protein_list->pro_list_[i].start_pos_, protein_list->pro_list_[i].length_);
      i--;
    }
  }
  bool add;
  for(int i = 0 ; i < list_prot_temp->nb_prot_ ; i++){
    add = false;
    for(int j = 0 ; j < protein_list->nb_prot_ ; j++){
      if(list_prot_temp->start_pos_[i] == protein_list->pro_list_[j].start_pos_){
        add = true;
      }
    }
    if(!add){
      protein_list->insert_protein(list_prot_temp->start_pos_[i], list_prot_temp->length_[i], list_prot_temp->basal_level_[i],
                                   list_prot_temp->hamming_dist_[i], list_prot_temp->dist_next_prot_[i], i_prot, i_prot);
      i_prot++;
    }
  }
  protein_list->sort_protein_list();
  protein_list->dna_length_ = list_prot_temp->dna_length_;
  return protein_list;
}
