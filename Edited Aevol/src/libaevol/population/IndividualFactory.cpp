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




// ============================================================================
//                                   Includes
// ============================================================================
#include "IndividualFactory.h"
#include "ExpManager.h"
#include <limits>

namespace aevol {


//##############################################################################
//                                                                             #
//                           Class IndividualFactory                           #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
/**
 * Create an individual with random sequences
 */
Individual* IndividualFactory::create_random_individual(
    ExpManager* exp_m,
    int32_t id,
    std::shared_ptr<MutationParams> param_mut,
    std::shared_ptr<JumpingMT> mut_prng,
    std::shared_ptr<JumpingMT> stoch_prng,
    const Habitat& habitat,
    double w_max,
    int32_t min_genome_length,
    int32_t max_genome_length,
    int32_t chromosome_initial_length,
    bool allow_plasmids,
    bool plasmid_initial_gene,
    int32_t plasmid_initial_length,
    char* strain_name,
    std::shared_ptr<JumpingMT> local_prng,
    bool better_than_flat)
{
  // Create a genome-less individual with the provided parameters

  #ifndef __REGUL
  Individual *indiv = new Individual(exp_m,
                         mut_prng,
                         stoch_prng,
                         param_mut,
                         w_max,
                         min_genome_length,
                         max_genome_length,
                         allow_plasmids,
                         id,
                         strain_name,
                         0);
  #else
  Individual_R* indiv = new Individual_R(exp_m,
                                     mut_prng,
                                     stoch_prng,
                                     param_mut,
                                     w_max,
                                     min_genome_length,
                                     max_genome_length,
                                     allow_plasmids,
                                     id,
                                     strain_name,
                                     0);
  #endif

  // Give it a randomly generated genome
  indiv->add_GU(indiv, chromosome_initial_length, local_prng);

  // If it was requested that the generated individual be better than a flat
  // one, test whether it is and if not re-generate until the condition is
  // satisfied
  double env_metabolic_area;
  if (better_than_flat) {
      //printf("ONE GOOD GENE\n");

#ifdef __REGUL
      dynamic_cast<Habitat_R*>(const_cast<Habitat*>(&habitat))->ApplyVariation();
    env_metabolic_area = dynamic_cast<Habitat_R*>(const_cast<Habitat*>(&habitat))->phenotypic_target_handler().
      mean_environmental_area(exp_m->exp_s());
#else
    env_metabolic_area = habitat.phenotypic_target_handler().mean_environmental_area();
#endif

    indiv->EvaluateInContext(habitat);
    //TESTING
    //exit(EXIT_FAILURE);

    #ifdef BASE_2
    double r_compare = round((indiv->dist_to_target_by_feature(METABOLISM)-env_metabolic_area) * 1E6) / 1E6;
    #elif BASE_4
    double r_compare = round((indiv->dist_to_target_by_feature(METABOLISM)-env_metabolic_area) * 1E4) / 1E4;
    #endif
      //printf("Dist to target (%lf) du nouveau clone : %e (%e) --> %lu\n", env_metabolic_area, indiv->dist_to_target_by_feature(METABOLISM),r_compare,indiv->protein_list().size());

    // indiv->dist_to_target_by_feature(METABOLISM) >= env_metabolic_area
    int32_t counter = 0;
    while (r_compare >= 0.0) {
      counter++;
      if (counter % 1000 == 0)
        printf("Counter %d\n",counter);
#ifdef __REGUL
        indiv->evaluated_ = false;
      indiv->set_networked(false);
#endif

      // Replace the former chromosome by a new random one and re-evaluate the
      // individual
      indiv->remove_GU(0);
      indiv->add_GU(indiv, chromosome_initial_length, local_prng);
      indiv->EvaluateInContext(habitat);
      //debug :

      #ifdef BASE_2
      r_compare = round((indiv->dist_to_target_by_feature(METABOLISM)-env_metabolic_area) * 1E6) / 1E6;
      #elif BASE_4 
      r_compare = round((indiv->dist_to_target_by_feature(METABOLISM)-env_metabolic_area) * 1E4) / 1E4;
      #endif
      //printf("Dist to target (%lf) du nouveau clone : %lf (%lf) --> %lu\n", env_metabolic_area, indiv->dist_to_target_by_feature(METABOLISM),r_compare,indiv->protein_list().size());
    }
  }
  if (allow_plasmids) // We create a plasmid
  {
    if (plasmid_initial_gene) {
      // The plasmid is generated independently from the chromosome
      indiv->add_GU(indiv, plasmid_initial_length, local_prng);

      if (better_than_flat) {
        indiv->EvaluateInContext(habitat);

        while (indiv->genetic_unit(1).
            dist_to_target_by_feature(METABOLISM) >= env_metabolic_area) {
          indiv->remove_GU(1);
          indiv->add_GU(indiv, plasmid_initial_length, local_prng);
          indiv->EvaluateInContext(habitat);
        }
      }
    }
    else {
      // The plasmid is a copy of the chromosome
      char* plasmid_genome = new char[chromosome_initial_length + 1];
      strncpy(plasmid_genome,
              indiv->genetic_unit_list().back().sequence(),
              chromosome_initial_length + 1);
      indiv->add_GU(plasmid_genome, chromosome_initial_length);
    }
  }

  // Insert a few IS in the sequence
  /*if (ae_common::init_params->init_method() & WITH_INS_SEQ)
  {
    // Create a random sequence
    int32_t seq_len = 50;
    char* ins_seq = new char[seq_len+1];
    int16_t nb_insert = 50;
    int16_t nb_invert = 50;

    for (int32_t i = 0 ; i < seq_len ; i++)
    {
      ins_seq[i] = '0' + ae_common::sim->prng->random(NB_BASE);
    }
    ins_seq[seq_len] = '\0';


    // Insert the sequence at random positions
    Mutation* mut1 = NULL;
    for (int16_t i = 0 ; i < nb_insert ; i++)
    {
      mut1 = indiv->genetic_unit(0)->dna()->do_insertion(ins_seq, seq_len);
      delete mut1;
    }


    // Invert the sequence and insert it at random positions
    char* inverted_seq = new char[seq_len+1];
    for (int32_t i = 0 ; i < seq_len ; i++)
    {
      inverted_seq[i] = (ins_seq[seq_len-1-i] == '1') ? '0' : '1';
    }
    inverted_seq[seq_len] = '\0';

    for (int16_t i = 0 ; i < nb_invert ; i++)
    {
      mut1 = indiv->genetic_unit(0)->dna()->do_insertion(inverted_seq, seq_len);
      delete mut1;
    }

    delete [] ins_seq;
    delete [] inverted_seq;
  }*/

  // If the individual hasn't been evaluated yet, do it
  if (not better_than_flat) {
    indiv->EvaluateInContext(habitat);
  }

  // Compute the "good" individual's statistics
  indiv->compute_statistical_data();

  return indiv;
}

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
