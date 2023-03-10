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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================

#include <unordered_map>
#include <algorithm>
#include <fstream>
// =================================================================
//                            Project Files
// =================================================================
#include "Individual_R.h"
#include "ExpManager.h"
#include "HybridFuzzy.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_individual_R                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

/*
 * Used at initialization
*/
Individual_R::Individual_R(ExpManager* exp_m,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng,
                       std::shared_ptr<MutationParams> param_mut,
                       double w_max,
                       int32_t min_genome_length,
                       int32_t max_genome_length,
                       bool allow_plasmids,
                       int32_t id,
                       const char* strain_name,
                       int32_t age) : Individual(exp_m,mut_prng,stoch_prng,param_mut,w_max,min_genome_length,
                                                 max_genome_length,allow_plasmids,id,strain_name,age) {

  _indiv_age = 0;
  _networked = false;
  _dist_sum = 0;

  _inherited_protein_list.resize(0);

      fitness_tab_ = nullptr;


}

Individual_R::Individual_R(const Individual_R& other)
    : Individual( other )
{
  _indiv_age = 0;
  _networked = false;
  _dist_sum = 0;

  fitness_tab_ = nullptr;
  if (exp_m_->exp_s()->get_with_heredity()) {
    _inherited_protein_list.clear();

    for (const auto& prot : other.protein_list_) {
      if (prot->concentration() >
          other.exp_m_->exp_s()->get_protein_presence_limit() && !((Protein_R*)prot)->is_signal()) {
        Protein_R* inherited_prot = new Protein_R(prot->get_gen_unit(),
                                                  (Protein_R&) *prot);
        inherited_prot->set_inherited(true);
        _inherited_protein_list.push_back(inherited_prot);
      }
    }
    //printf("Size of inherited is %d\n",_inherited_protein_list.size());
  }
}

Individual_R::Individual_R( Individual_R* parent, int32_t id,
                            std::shared_ptr<JumpingMT> mut_prng,
                            std::shared_ptr<JumpingMT> stoch_prng)
        : Individual( parent, id, mut_prng, stoch_prng )
{
  //~ printf( "ae_individual_R( parent ) : I have %d inherited proteins\n", parent->get_protein_list()->get_nb_elts() );
    _indiv_age = 0;
    _networked = false;
    _dist_sum = 0;
  fitness_tab_ = nullptr;
  if (exp_m_->exp_s()->get_with_heredity()) {
    //_inherited_protein_list.resize(parent->protein_list_.size());
    //printf("%llu -- Parent protein size: %d\n",id, parent->protein_list_.size());

    for (const auto& prot : parent->protein_list_) {
      if (prot->concentration() >
          parent->exp_m_->exp_s()->get_protein_presence_limit() && ! ((Protein_R*)prot)->is_signal()) {
        Protein_R* inherited_prot = new Protein_R(prot->get_gen_unit(),
                                                  (Protein_R&) *prot, this->exp_m_);
        inherited_prot->set_inherited(true);
        _inherited_protein_list.push_back(inherited_prot);
        //printf("Add new herited protein\n");
      }
    }

    //printf("%llu -- Inherited protein size is %d\n",id,_inherited_protein_list.size());
  }
}

Individual_R::Individual_R(ExpManager* exp_m, gzFile backup_file) : Individual( exp_m, backup_file )
{
    _indiv_age = 0;
  _networked = false;
  fitness_tab_ = nullptr;
  if( exp_m_->exp_s()->get_with_heredity() )
  {
    // Retreive inherited proteins
    int32_t nb_inherited_proteins = 0;
    gzread( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );

    //printf("Nb prot %d\n",nb_inherited_proteins);

    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
	  _inherited_protein_list.push_back( new Protein_R( backup_file ) );
    }
  }  
}

// =================================================================
//                             Destructors
// =================================================================
Individual_R::~Individual_R( void ) noexcept
{

  if (exp_m_->exp_s()->get_with_heredity()) {
    for (unsigned int i = 0; i < _inherited_protein_list.size(); i++)
      delete _inherited_protein_list[i];

    _inherited_protein_list.clear();
  }

  for (unsigned int i = 0; i < _rna_list_coding.size(); i++) {
    _rna_list_coding[i] = NULL;
  }

  std::unordered_map<int,Protein_R * >::iterator iter = signal_list.begin();

  while (iter != signal_list.end())
  {
    Protein_R* clone = iter->second;
    delete clone;
    iter++;
  }

  signal_list.clear();

  _rna_list_coding.clear();

  delete [] fitness_tab_;
}

// =================================================================
//                            Public Methods
// =================================================================
Individual_R* Individual_R::CreateClone(const Individual_R* dolly, int32_t id) {
  Individual_R* indiv = new Individual_R(*dolly);
  indiv->set_id(id);
  return indiv;
}

void Individual_R::Evaluate(bool no_signal) {
		EvaluateInContext(dynamic_cast<const Habitat_R&> (habitat()), no_signal);
}



void Individual_R::EvaluateOneAfterAnother(const Habitat_R& habitat) {
    if (evaluated_ == true)
        return; // Individual has already been evaluated, nothing to do.

    for (const auto &prot : protein_list_) {
        ((Protein_R *) prot)->reset_concentration();
    }

    if (!_networked) {
        init_indiv(habitat);
    }

    _global_dist_sum = 0;

    double fitness_sum = 0;
    std::set<int> *eval = exp_m_->exp_s()->get_list_eval_step();

    fitness_tab_ = new double[habitat.number_of_phenotypic_target_models()];

    // i is thus the age of the individual
    for (int16_t env_i = 0; env_i < habitat.number_of_phenotypic_target_models(); env_i++) {

        //Set the concentration of signals for this age
        for (auto prot1 : signal_list) {
            prot1.second->set_concentration(0.0);
        }

        for (Protein_R *prot2 : habitat.phenotypic_target_model(env_i).signals()) {
            signal_list[prot2->get_id()]->set_concentration(0.9);
        }

        _dist_sum = 0;
        exp_m_->exp_s()->set_nb_degradation_step(10);
        exp_m_->exp_s()->set_nb_indiv_age(10);

        //printf("Eval env %d for %d (%d)\n",env_i,exp_m_->exp_s()->get_nb_indiv_age(), exp_m_->exp_s()->get_nb_degradation_step());
        for (int16_t i = 1; i <= exp_m_->exp_s()->get_nb_indiv_age(); i++) {

            for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
                one_step();
            }

            // If we have to evaluate the individual at this age
            if (eval->find(i) != eval->end()) {
                //        printf("Eval !\n");
                eval_step_one_after_another(habitat, env_i);
            }
        }


        final_step_one_after_another(habitat, env_i);

        fitness_sum += fitness();
        fitness_tab_[env_i] = fitness();

        //printf("%d -- %d : %e %e -- %e\n",id(),env_i,fitness(),fitness_sum,
        //       habitat.phenotypic_target_model(env_i).area_by_feature(METABOLISM));
    }


    fitness_ = fitness_sum / ((double) habitat.number_of_phenotypic_target_models());
    dist_to_target_by_feature_[METABOLISM] =
            _global_dist_sum / ((double) habitat.number_of_phenotypic_target_models());
    /*printf("%e %d : %e\n",_global_dist_sum,habitat.number_of_phenotypic_target_models(),dist_to_target_by_feature_[METABOLISM]);*/
}


void Individual_R::EvaluateInContext(const Habitat_R& habitat, bool no_signal) {
    if (habitat.phenotypic_target_handler().var_method() == ONE_AFTER_ANOTHER)
        EvaluateOneAfterAnother(dynamic_cast<const Habitat_R &> (habitat));
    else {
        //printf("dist sum %d : %lf\n",id(),_dist_sum);
        if (evaluated_ == true)
            return; // Individual has already been evaluated, nothing to do.
        for (const auto &prot : protein_list_) {
            ((Protein_R *) prot)->reset_concentration();
        }
    }

    if (!_networked) {
        init_indiv(habitat);
    }
    /*for (const auto& prot : protein_list_) {
      printf("AT INIT %d ID %d Concentration of %d is %lf\n",AeTime::time(),id(),
             ((Protein_R*)prot)->get_id(),prot->concentration());
    }*/

    /*if (habitat.phenotypic_target_handler().var_method() == ONE_AFTER_ANOTHER) {
    _dist_sum = 0;
      exp_m_->exp_s()->set_nb_degradation_step(10);
      exp_m_->exp_s()->set_nb_indiv_age(10);

    //printf("Eval env %d for %d (%d)\n",env_i,exp_m_->exp_s()->get_nb_indiv_age(), exp_m_->exp_s()->get_nb_degradation_step());
    for (int16_t i = 1; i <= exp_m_->exp_s()->get_nb_indiv_age(); i++) {

      for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
        one_step();
      }

      // If we have to evaluate the individual at this age
      if (eval->find(i) != eval->end()) {
  //        printf("Eval !\n");
        eval_step_one_after_another(habitat, env_i);
      }
    }


    final_step_one_after_another(habitat, env_i);

    fitness_sum += fitness();
    fitness_tab_[env_i] = fitness();

    //printf("%d -- %d : %e %e -- %e\n",id(),env_i,fitness(),fitness_sum,
    //       habitat.phenotypic_target_model(env_i).area_by_feature(METABOLISM));
  }


  fitness_ = fitness_sum/((double)habitat.number_of_phenotypic_target_models());
    dist_to_target_by_feature_[METABOLISM] = _global_dist_sum / ((double)habitat.number_of_phenotypic_target_models());
    *//*printf("%e %d : %e\n",_global_dist_sum,habitat.number_of_phenotypic_target_models(),dist_to_target_by_feature_[METABOLISM]);*//*

    } else {*/
    _dist_sum = 0;


//  if (id_==543 && AeTime::time() == 5895){
//    std::vector<Protein*> protein_vector;
//         for (auto prot: protein_list_)
//                 protein_vector.push_back(prot);

//   std::sort(protein_vector.begin(), protein_vector.end(),
//             [](Protein*a, Protein*b) { return a->shine_dal_pos() < b->shine_dal_pos();});

//    for (auto prot: protein_vector)
//      printf("%d -- CPU -- Protein %d : %.18e\n", 0, prot->first_translated_pos(),
//             prot->concentration());
//  }
    std::set<int> *eval = exp_m_->exp_s()->get_list_eval_step();
    // i is thus the age of the individual
           // printf("Evaluate for %d\n",exp_m_->exp_s()->get_nb_indiv_age());
    for (int16_t i = 1; i <= exp_m_->exp_s()->get_nb_indiv_age(); i++) {
        //Set the concentration of signals for this age
        for (auto prot1 : signal_list) {
            prot1.second->set_concentration(0.0);
        }
        if (no_signal)
            for (Protein_R *prot2 : habitat.phenotypic_target(i).signals()) {
                signal_list[prot2->get_id()]->set_concentration(0.9);
            }


        for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
            one_step();
        }

      //   if (id_ == 0) {

      // std::ofstream concentrationfile;
      // concentrationfile.open("stats/online_concentration_visu.csv",std::ofstream::app);

      //   // Save concentration
      //   int32_t id = 0;
      //   for (auto prot: protein_list_) {
      //     {
      //         concentrationfile<<AeTime::time()<<","<<id_<<","<<i<<","<<id<<","<<prot->concentration()<<","<<(((Protein_R*)prot)->is_signal()? 1 : 0)<<std::endl;
      //         id++;
      //     }
      //   }

      //   concentrationfile.close();
      //   }
//  if (id_==68 && AeTime::time() == 4){
//    std::vector<Protein*> protein_vector;
//         for (auto prot: protein_list_)
//                 protein_vector.push_back(prot);

//   std::sort(protein_vector.begin(), protein_vector.end(),
//             [](Protein*a, Protein*b) { return a->shine_dal_pos() < b->shine_dal_pos();});

//    for (auto prot: protein_vector)
//      printf("%d -- CPU -- Protein %d : %.18e\n", i, prot->first_translated_pos(),
//             prot->concentration());
//  }

        // if (id() == 326)
        // for (const auto& prot : protein_list_) {
        //   if (!((Protein_R*)prot)->is_signal())
        //     printf("CLASSIC -- %d -- %d -- Protein %d : %.18e\n",i,id(),prot->first_translated_pos(),prot->concentration());
        // }

        // If we have to evaluate the individual at this age
        if (eval->find(i) != eval->end()) {//|| (id_==543 && AeTime::time() == 5895)) {// ||( (id_ == 70) && (AeTime::time()>=1570))){
            //if (id_ % 1024 == 1) printf("Eval at %d\n",i);
            eval_step(habitat, i);
            // if (id_==68 && AeTime::time() == 4) {
            // if (AeTime::time() == 1 && id_ == 41) {
              // if (id_==206)
              //   printf("%d -- CPU -- Evaluate Network at %d (env id %d) :: %lf %lf -- %lf\n",id_,i,
              //     habitat.phenotypic_target( i ).get_id(),
              //      _dist_sum,dist_to_target_by_feature_[METABOLISM],
              //      habitat.phenotypic_target( i ).fuzzy()->get_geometric_area());
            //   // phenotype_->print();
            //   // habitat.phenotypic_target( i ).fuzzy()->print();
            // }
        }
    }


    final_step(habitat, exp_m_->exp_s()->get_nb_indiv_age());
}

//    protein_list_.clear();
//    protein_list_ = _initial_protein_list;


    //initialized = false;

//}

void Individual_R::EvaluateInContext(const Habitat& habitat) {
    if (habitat.phenotypic_target_handler().var_method() == ONE_AFTER_ANOTHER)
        EvaluateOneAfterAnother(dynamic_cast<const Habitat_R&> (habitat));
    else
        EvaluateInContext(dynamic_cast<const Habitat_R&> (habitat));
}

void Individual_R::init_indiv() {
  init_indiv(dynamic_cast<const Habitat_R&> (grid_cell_->habitat()));
}

void Individual_R::init_indiv(const Habitat_R& habitat)
{
  // ---------------------------------------------------------------------------
  // 1) Transcription - Translation - Folding - make_protein_list
  // ---------------------------------------------------------------------------
  transcribed_ = false;
  translated_ = false;
  folded_ = false;

  do_transcription_translation_folding();

  if (phenotype_ != NULL) {
    delete phenotype_;
    delete phenotype_activ_;
    delete phenotype_inhib_;

    phenotype_ = NULL;
    phenotype_activ_ = NULL;
    phenotype_inhib_ = NULL;
  }
  phenotype_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_activ_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_inhib_ = FuzzyFactory::fuzzyFactory->create_fuzzy();

  //----------------------------------------------------------------------------
  // 2) Make a list of all the rna present in the individual
  //    and initialise the concentrations of the proteins
  //----------------------------------------------------------------------------
  make_rna_list();

  int local_prot_id = 0;
  for (auto prot : protein_list_)
    ((Protein_R*)prot)->set_local_id(local_prot_id++);

  _initial_protein_list = protein_list_;

  //_protein_list.insert(_protein_list.end(), habitat.signals().begin(), habitat.signals().end());
  for(Protein_R* prot : habitat.signals()) {
    Protein_R* cloned = new Protein_R(prot);

    signal_list.insert({cloned->get_id(),cloned});
    protein_list_.push_back(cloned);
  }

  //----------------------------------------------------------------------------
  // 3) Create influence graph (including the signals)
  //----------------------------------------------------------------------------
  //printf("Protein %ld RNA %ld\n",protein_list_.size(),rna_list_.size());
  set_influences();

  _networked = true;
  //initialized = true;



}

void Individual_R::one_step( void )
{
  //----------------------------------------------------------------------------
  // 4) Make the individual "live its life" and compute partial phenotypes and
  //    fitnesses
  //----------------------------------------------------------------------------
  update_concentrations();

}

void Individual_R::eval_step( const Habitat_R& habitat, int16_t age ) {
  update_phenotype();
  distance_to_target_computed_ = false;
  phenotype_computed_ = true;

  for (int i=0; i<NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0;
  }

    // if (id_==68 && AeTime::time() == 4) {
    //     printf("Compute distance\n");
    // }
  compute_distance_to_target( habitat.phenotypic_target( age ) );

  // if (id_==0)
  //   printf("At %d dist is %lf on %d -- %lf\n",age,dist_to_target_by_feature_[METABOLISM],
  //         habitat.phenotypic_target( age ).get_id(), habitat.phenotypic_target( age ).area_by_feature(METABOLISM));

  _dist_sum += dist_to_target_by_feature_[METABOLISM];

  /*if (id_ % 1024 == 1)
    printf("Dist to target à l'age %d du nouveau clone : %f -- %f\n",
           age, dist_to_target_by_feature_[METABOLISM],_dist_sum);*/
}

void Individual_R::eval_step_one_after_another( const Habitat_R& habitat, int16_t env_id ) {
  update_phenotype();
  distance_to_target_computed_ = false;
  phenotype_computed_ = true;

  for (int i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0;
  }

  compute_distance_to_target(habitat.phenotypic_target_model(env_id));


  _dist_sum += dist_to_target_by_feature_[METABOLISM];
/*
  printf("Meta error %e (%e) out of %e\n",dist_to_target_by_feature_[METABOLISM],_dist_sum,habitat.phenotypic_target_model(env_id).area_by_feature(METABOLISM));
*/
}

void Individual_R::final_step( const Habitat_R& habitat, int16_t age ) {
  dist_to_target_by_feature_[METABOLISM] = _dist_sum / (double) (exp_m_->exp_s()->get_list_eval_step()->size());
  //printf("AT %d ID %d Meta error : %lf \n",AeTime::time(),id(),dist_to_target_by_feature_[METABOLISM]);


  fitness_computed_=false;
  // yoram attention il peut y avoir des soucis si on utilise des environnements segmentés ici
  compute_fitness(habitat.phenotypic_target( age ));

  /*if (id_ % 1024 == 1)
    printf("Nombre final d'évaluations : %d -- %f -- %e -- %e\n", exp_m_->exp_s()->get_list_eval_step()->size(),
           dist_to_target_by_feature_[METABOLISM], fitness_, fitness_by_feature_[METABOLISM]);*/

  phenotype_computed_ = true;
  // if (id_ == 206)
  //   printf("%d -- CPU -- Finalize Network :: %lf %lf (%ld)\n",id_,dist_to_target_by_feature_[METABOLISM], fitness_,
  //           (exp_m_->exp_s()->get_list_eval_step()->size()));
}

void Individual_R::final_step_one_after_another( const Habitat_R& habitat, int16_t env_id ) {
    dist_to_target_by_feature_[METABOLISM] = _dist_sum / (double) (exp_m_->exp_s()->get_list_eval_step()->size());
    _global_dist_sum += dist_to_target_by_feature_[METABOLISM];
    distance_to_target_computed_ = true;
/*
    printf("FINAL -> Meta error %e (%e) out of %e\n",dist_to_target_by_feature_[METABOLISM],_dist_sum,habitat.phenotypic_target_model(env_id).area_by_feature(METABOLISM));
*/

    fitness_computed_ = false;
    // yoram attention il peut y avoir des soucis si on utilise des environnements segmentés ici
    compute_fitness(habitat.phenotypic_target_model(env_id));

    phenotype_computed_ = true;
}

void Individual_R::set_influences()
// Compute the influence of each protein over each coding RNA
// As non-coding RNAs are completely inert, we don't care about their concentration
// so we don't care if proteins activate or inhibit their transcription.
{
  /*for (const auto& prot : protein_list_) {
   if (id() == 12608 && ((Protein_R*)prot)->is_signal())
      printf("ID 12608 set_influence Signal is %ld\n",((Protein_R*)prot)->get_id());
  }*/

	  for(auto& rna : _rna_list_coding) {
		  rna->set_influences( protein_list_, id() );
	  }
}

void Individual_R::update_concentrations( void )
{
	// Compute all the changes that will be applied to the concentrations
	// Concentrations must not be changed at this stage
  for (auto& prot : protein_list_) {
    if (!((Protein_R*)prot)->is_signal()) ((Protein_R*)prot)->compute_delta_concentration(exp_m_);

	}

	// Apply the changes in concentrations we have just computed
//        int j = 0;
  for (auto& prot : protein_list_) {
		if (!((Protein_R*)prot)->is_signal()) {

                  // if (id_==70 && AeTime::time() == 1595)
                  // if (id_==188 && AeTime::time() == 524)
                  //   printf("CPU -- Protein %d : %lf + %lf\n",prot->shine_dal_pos(),prot->concentration(),((Protein_R*)prot)->_delta_concentration);

                  ((Protein_R*)prot)->update_concentration();
                }

//    if (id_ == 3) printf("Protein %d :: %lf DELTA %lf\n",j,prot->concentration(),((Protein_R*)prot)->_delta_concentration);
//    j++;
	}
}

// Multiply the concentration of each protein by <factor>
void Individual_R::multiply_concentrations( double factor )
{
  for (auto& prot : protein_list_) {
	 	  ((Protein_R*)prot)->multiply_concentration( factor );
	}
}

/*int8_t Individual_R::get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  int8_t new_quadon = new_get_quadon(gen_unit,strand,pos);
  int8_t old_quadon = old_get_quadon(gen_unit,strand,pos);

  if (new_quadon != old_quadon) {
    printf("New/Old quadon are different %d %d\n",new_quadon,old_quadon);
    exit(-1);
  }
}*/

/*int8_t Individual_R::old_get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->dna()->data();
  int32_t  len    = gen_unit->dna()->length();
  int8_t quadon   = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[Utils::mod((pos+i),len)] == '1' )
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }
  else  // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {

      if ( dna[Utils::mod((pos-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }




  return quadon;
}*/

int8_t Individual_R::get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->dna()->data();
  int32_t  len    = gen_unit->dna()->length();
  int8_t quadon[4];

  if ( strand == LEADING )
  {
#ifdef __SIMD
#pragma omp simd
#endif
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      quadon[i] = (dna[((pos+i) % len < 0 ? (pos+i) % len + len : (pos+i) % len)] == '1') ? 1 << (QUADON_SIZE - i - 1) : 0;
    }
  }
  else  // ( strand == LAGGING )
  {
#ifdef __SIMD
#pragma omp simd
#endif
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      quadon[i] = (dna[((pos-i) % len < 0 ? (pos-i) % len + len : (pos-i) % len)] != '1') ? 1 << (QUADON_SIZE - i - 1) : 0;
    }
  }

  return quadon[0]+quadon[1]+quadon[2]+quadon[3];
}

void Individual_R::save( gzFile backup_file ) const
{
  Individual::save( backup_file );
  // Test if there is heredity, and if the generation is the first one (no inherited protein list).
  if (this->exp_m_->exp_s()->get_with_heredity() )
  {
    // Write inherited proteins
    int32_t nb_inherited_proteins = _inherited_protein_list.size();
    gzwrite( backup_file, &nb_inherited_proteins,  sizeof(int32_t) );

    for (auto& prot : _inherited_protein_list) {
    	prot->save( backup_file );
    }
  }
}
// =================================================================
//                           Protected Methods
// =================================================================
void Individual_R::make_protein_list( void )
{
	  Individual::make_protein_list();

    if (this->exp_m_->exp_s()->get_with_heredity()) {
      for (auto& prot : _inherited_protein_list) {
        protein_list_.push_back(prot);
      }
    }
}

void Individual_R::make_rna_list( void )
{
  Individual::make_rna_list();
  _rna_list_coding = {};

  int local_rna_id = 0;

  // Parse the newly created RNA list and copy the coding RNAs in _rna_list_coding.
  for (const auto& gen_unit: genetic_unit_list_) {
    //GeneticUnit* genu = const_cast<GeneticUnit*>(&gen_unit);
    // Create proxies
    const auto& rna_list = gen_unit.rna_list();
    //const auto& lead = rna_list[LEADING];
    //const auto& lagg = rna_list[LAGGING];

    // append pointers to rna material to local _rna_list
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: rna_list[strand]) {
        //TODO Ugly fix, change it to avoid memory usage double
        if (rna.is_coding()) {
          Rna_R* prna =  const_cast<Rna_R*>(&rna);
          //printf("COPY OR NOT : %ld == %ld",prna->get_id(),((Rna_R)rna).get_id());
          _rna_list_coding.push_back(
             prna);//new Rna_R(genu, rna));

          prna->set_local_id(local_rna_id++);
        }
    }

  }
}

void Individual_R::update_phenotype( void )
{
  // We will use two fuzzy sets :
  //   * _phenotype_activ for the proteins realising a set of functions
  //   * _phenotype_inhib for the proteins inhibitting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets

  delete phenotype_activ_;
  delete phenotype_inhib_;
  delete phenotype_;
  // printf("%d -- Delete phenotype\n",id_);

  phenotype_activ_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_inhib_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  // printf("%d -- Allocate phenotype : %u %u %u\n",id_,((Fuzzy*)phenotype_activ_)->points().size(),
  //     ((Fuzzy*)phenotype_inhib_)->points().size(),((Fuzzy*)phenotype_)->points().size());

  std::vector<Protein *> protein_vector;
  for (auto prot : protein_list_) { 
      protein_vector.emplace_back(prot);
        //protein_vector.back()->concentration_ = 0.0;
      //for (auto rna : prot.rna_list()) {
      //    protein_vector.back()->add_RNA(rna);
      //}
    // }
  }

        bool verbose = false;
        // if (AeTime::time() ==447 &&  indiv_->grid_cell()->x() * indiv_->exp_m()->grid_height() + indiv_->grid_cell()->y()==966) {
        //   verbose = true;
        // }


  sort(protein_vector.begin(), protein_vector.end(),
       [](Protein *a, Protein *b) { return *a < *b;});
  for(auto prot : protein_vector) {
    if ( ((Protein_R*)prot)->is_functional() )
    {
      // if (id_==524 && AeTime::time() == 188)
      //  printf("Add triangle %lf %lf %lf (%lf %lf)\n",((Protein_R*)prot)->mean(),
      //         ((Protein_R*)prot)->width(),
      //         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(),
      //         ((Protein_R*)prot)->height(), ((Protein_R*)prot)->concentration() );

      // if (id_==131 && AeTime::time() == 4)
      //  printf("CPU -- Add triangle %lf %lf %lf (%lf %lf)\n",
      //         ((Protein_R*)prot)->mean(),
      //         ((Protein_R*)prot)->width(),
      //         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(),
      //         ((Protein_R*)prot)->height(), ((Protein_R*)prot)->concentration() );

      if ( ((Protein_R*)prot)->height() > 0 )
      {
//    	  added=true;
        bool verbose = false;
        // if (id_==543 && AeTime::time() == 5895)
        //   verbose = true;
        phenotype_activ_->add_triangle(  ((Protein_R*)prot)->mean(),
                                         ((Protein_R*)prot)->width(),
                                         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(), verbose );

//        if (id_==120)
//          printf("Add triangle ACTIV %f %f %f (%f %f)\n",((Protein_R*)prot)->mean(),
//                 ((Protein_R*)prot)->width(),
//                 ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(),
//                 ((Protein_R*)prot)->height(), ((Protein_R*)prot)->concentration() );
      }
      else
      {
        bool verbose = false;
        // if (id_==543 && AeTime::time() == 5895)
        //   verbose = true;
        phenotype_inhib_->add_triangle(  ((Protein_R*)prot)->mean(),
                                         ((Protein_R*)prot)->width(),
                                         ((Protein_R*)prot)->height() * ((Protein_R*)prot)->concentration(), verbose );

      }

      // if (id_==41 && AeTime::time() == 1)
      //   printf("CPU -- Phenotype : %lf %lf\n",phenotype_activ_->get_geometric_area(),phenotype_inhib_->get_geometric_area());
    }
  }


    phenotype_activ_->clip(AbstractFuzzy::max,   Y_MAX);
    phenotype_inhib_->clip(AbstractFuzzy::min, - Y_MAX);

    phenotype_activ_->simplify();
    phenotype_inhib_->simplify();

    phenotype_->add(*phenotype_activ_);
    phenotype_->add(*phenotype_inhib_);
                // if (id_==543 && AeTime::time() == 5895) {printf("BEFORE CLIP\n"); phenotype_->print();}

    phenotype_->clip(AbstractFuzzy::min, Y_MIN);
    // if (id_==543 && AeTime::time() == 5895) {printf("BEFORE SIMPLIFY\n"); phenotype_->print();}
    phenotype_->simplify();
// if (id_==543 && AeTime::time() == 5895) {phenotype_->print();}
//  _phenotype->simplify();
//  if (id_==326)
//         printf("CPU -- Phenotype : %lf %lf :: %lf\n",phenotype_activ_->get_geometric_area(),phenotype_inhib_->get_geometric_area(),phenotype_->get_geometric_area());
}

void    Individual_R::clear_everything_except_dna_and_promoters() {
    _networked = false;
    _rna_list_coding.clear();
    _dist_sum = 0.0;

    Individual::clear_everything_except_dna_and_promoters();
}

void Individual_R::create_csv(char *directory_name) {

  for (const auto& prot : protein_list_) {
    ((Protein_R*) prot)->reset_concentration();
  }

  init_indiv();


  char filename[128];

  snprintf(filename, 127, "%s/best_concentration.csv", directory_name);

  FILE * drawingfile = fopen(filename, "w");

  if (drawingfile == NULL)
  {
    fprintf(stderr, "Error: could not create the file %s\n", filename);
  }


  fprintf(drawingfile, "Generation,Protein,Concentration\n");

  //int16_t life_time =  exp_m_->exp_s()->get_nb_indiv_age();
  //set the concentrations of proteins to their initial value
  double* concentrations = new double[protein_list_.size()]; // initialise le tableau de concentrations.
  //  int16_t prot_index = 0;
  int i = 0;
  for (const auto& prot : protein_list_) {
    concentrations[i] = ((Protein_R*)prot)->concentration();
    i++;
  }


  //std::set<int>* eval = exp_m_->exp_s()->get_list_eval_step();
  // i is thus the age of the individual
  for (int16_t i = 1; i <= exp_m_->exp_s()->get_nb_indiv_age(); i++) {
    //Set the concentration of signals for this age
    printf("Prot id from signal_list is ");
    for (auto prot1 : signal_list) {
      printf("%ld ",prot1.second->get_id());
      prot1.second->set_concentration(0.0);
    }
    printf("\n");

    printf("Prot id from phenotypic target is ");
    for (Protein_R* prot2 : dynamic_cast<const Habitat_R&>(habitat()).phenotypic_target(i).signals()) {
      printf("%ld ",prot2->get_id());
      signal_list[prot2->get_id()]->set_concentration(0.9);
    }
    printf("\n");

    for (int j = 0; j < exp_m()->exp_s()->get_nb_degradation_step(); j++)
      update_concentrations();

    //affichage des points n+1 dans la concentration
    //prot_index = 0;
    int proti = 0;
    //printf("Age[%d] : ",i);
    for (const auto& prot : protein_list_) {
      if (((Protein_R*)prot)->is_signal()) printf("Protein List %ld -> %f (signal: %d) (%p)\n",
                                                  ((Protein_R*)prot)->get_id(),prot->concentration(),((Protein_R*)prot)->is_signal(),prot);
    

      // morceau ajouté pour colorer les protéines en fonctions de leur paramètres
      concentrations[proti]=((Protein_R*)prot)->concentration();
      fprintf(drawingfile, "%d,%d,%lf\n",i,proti,prot->concentration());
      proti++;
    }

    update_phenotype();

    for (auto prot : signal_list) {
      printf("after %d prot %ld : %f\n",i,prot.second->get_id(),prot.second->concentration());
    }

    printf("Creating the EPS file with the phenotype of the chosen individual at step %d... ",i);
    fflush(stdout);
    draw_phenotype(dynamic_cast<const Habitat_R&>(habitat()).phenotypic_target( i ), directory_name,i);
    printf("OK\n");
    //printf("\n");
  }

    for (const auto& prot : protein_list_) {
      if (((Protein_R*)prot)->is_signal()) printf("Protein List %ld -> %f (signal: %d)\n",((Protein_R*)prot)->get_id(),prot->concentration(),((Protein_R*)prot)->is_signal() );
    }

  delete[] concentrations;
}


void Individual_R::draw_phenotype(const PhenotypicTarget& target, char* directoryName, int generation)
{
  char filename[128];
  snprintf(filename, 127, "%s/best_phenotype_%d.csv", directoryName,generation);
  FILE * drawingfile = fopen(filename, "w");

  if (drawingfile == NULL)
  {
    fprintf(stderr, "Error: could not create the file %s\n", filename);
    return;
  }


  fprintf(drawingfile, "Generation,IndividualOrPhenotype,X,Y\n");

  if (exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)phenotype_activ())->points())
      fprintf(drawingfile, "%d,0,%lf,%lf\n",generation,p.x,p.y);
  else {
    for (int i=0; i < ((HybridFuzzy*)phenotype())->get_pheno_size(); i++) {

      float xi =  (i / (float)
          ((HybridFuzzy*) phenotype())->get_pheno_size());
      fprintf(drawingfile, "%d,0,%lf,%lf\n",generation,xi,((HybridFuzzy*) phenotype())->points()[i]);
    }
  }


  // ------------------
  //  draw environment
  // ------------------
  if (exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)target.fuzzy())->points())
      fprintf(drawingfile, "%d,1,%lf,%lf\n", generation,p.x, p.y);
  else
    for (int i=0; i < ((HybridFuzzy*)target.fuzzy())->get_pheno_size(); i++) {
      float xi =  (i / (float)
          ((HybridFuzzy*)target.fuzzy())->get_pheno_size());
      fprintf(drawingfile, "%d,1,%lf,%lf\n",generation,xi, ((HybridFuzzy*) target.fuzzy())->points()[i]);
    }

  fclose(drawingfile);
}

} // namespace aevol
