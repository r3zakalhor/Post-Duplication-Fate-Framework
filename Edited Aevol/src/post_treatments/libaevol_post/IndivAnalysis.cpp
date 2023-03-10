//
// Created by dparsons on 31/05/16.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "IndivAnalysis.h"

#include "aevol.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
IndivAnalysis::IndivAnalysis(const Individual& indiv) : Individual(indiv) {

};

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
/**
 * Compute reproduction theoretical proportion of neutral offsprings.
 *
 * Compute the theoretical proportion of neutral offsprings given Carole's
 * formula, based on the mutations and rearrangement rates and not on multiple
 * replications.
 *
 * \return theoretical proportion of neutral offsprings
 */
double IndivAnalysis::compute_theoritical_f_nu() {
  // We first have to collect information about genome structure.
  // Abbreviations are chosen according to Carole's formula.
  // Please notice that compared to the formula we have the beginning
  // and ends of neutral regions instead of 'functional regions'
  GeneticUnit& chromosome = genetic_unit_list_.front();
  int32_t L = chromosome.dna()->length();
  int32_t N_G = chromosome.nb_neutral_regions(); // which is not exactly Carole's original definition
  int32_t* b_i = chromosome.beginning_neutral_regions();
  int32_t* e_i = chromosome.end_neutral_regions();
  int32_t lambda = chromosome.nb_bases_in_neutral_regions();
  int32_t l = L - lambda; // nb bases in 'functional regions'

  int32_t* lambda_i = NULL;  // nb bases in ith neutral region
  if (N_G > 0) // all the chromosome may be functional
  {
    lambda_i = new int32_t[N_G];

    for (int32_t i = 0; i < N_G - 1; i++) {
      lambda_i[i] = e_i[i] - b_i[i] + 1;
    }
    if (b_i[N_G - 1] > e_i[N_G -
                           1]) // last neutral region is overlapping on the beginning of chromosome
    {
      lambda_i[N_G - 1] = (e_i[N_G - 1] + L) - b_i[N_G - 1] + 1;
    }
    else // no overlap
    {
      lambda_i[N_G - 1] = e_i[N_G - 1] - b_i[N_G - 1] + 1;
    }
  }

  // we now compute the probabilities of neutral reproduction for
  // each type of mutation and rearrangement and update Fv
  double Fv = 1;

  // mutation + insertion + deletion
  double nu_local_mutation = 1 - ((double) l) / L;
  Fv = pow(1 - point_mutation_rate() * (1 - nu_local_mutation), L);
  Fv *= pow(1 - small_insertion_rate() * (1 - nu_local_mutation), L);
  Fv *= pow(1 - small_deletion_rate() * (1 - nu_local_mutation), L);

  // inversion ~ two local mutations
  double nu_inversion = nu_local_mutation * nu_local_mutation;
  Fv *= pow(1 - inversion_rate() * (1 - nu_inversion), L);

  // translocation ~ inversion + insertion (mathematically)
  Fv *= pow(
      1 - translocation_rate() * (1 - nu_inversion * nu_local_mutation), L);

  // long deletion
  double nu_deletion = 0; // if N_G == 0, a deletion is always not neutral
  for (int32_t i = 0; i < N_G; i++) {
    nu_deletion += lambda_i[i] * (lambda_i[i] + 1);
  }
  nu_deletion /= ((double) 2 * L * L);
  Fv *= pow(1 - deletion_rate() * (1 - nu_deletion), L);

  // duplication ~ big deletion + insertion
  Fv *= pow(1 - duplication_rate() * (1 - nu_deletion * nu_local_mutation),
            L);

  if (lambda_i != NULL) delete[] lambda_i;

  return Fv;
}

/**
 *
 */
void IndivAnalysis::compute_experimental_f_nu(
    int32_t nb_indiv,
    std::shared_ptr<JumpingMT> prng,
    FILE* output_summary /* = nullptr*/,
    bool verbose /* = false*/,
    bool full_output /* = false*/) {
  double nb_pos = 0;
  double cumul_delta_err_pos = 0;
  double cumul_delta_fitness_pos = 0;
  double nb_neg = 0;
  double cumul_delta_err_neg = 0;
  double cumul_delta_fitness_neg = 0;
  double max_pos = 0;
  double max_fitness_pos = 0;
  double max_neg = 0;
  double max_fitness_neg = 0;
  double nb_neutral_genetic = 0;
  double nb_neutral_phenotypic = 0;
  int32_t nb_events = 0;

  double parent_metabolic_error = dist_to_target_by_feature(METABOLISM);
  double parent_fitness = fitness();

  fprintf(output_summary,
	  "%" PRId64 " ",AeTime::time());
  
  #pragma omp parallel for
  for (int32_t i = 0; i < nb_indiv; i++) {
    Individual mutant(this, 0, prng, prng);
    int32_t nb_events = 0;

    // Perform transfer, rearrangements and mutations
    if (not mutant.allow_plasmids()) {
      const GeneticUnit* chromosome = &(mutant.genetic_unit_list().front());
        
      nb_events = chromosome->dna()->perform_mutations(id_);
    }
    else {
      printf("WARNING: Mutational Robustness does not handle multiple "
                 "Genetic Units\n");
    }
    if (nb_events == 0)
      {
  #pragma omp atomic
	nb_neutral_phenotypic++;

	if (full_output) {
    #pragma omp critical(full_output)
    {
      fprintf(output_summary,"%.15e ",0.0);
    }
	}
	
      }
    else
      {
	mutant.EvaluateInContext(habitat());
	double new_metabolic_error = mutant.dist_to_target_by_feature(
								      METABOLISM);
	double new_fitness = mutant.fitness();
	
	if (new_metabolic_error == parent_metabolic_error) {
	  #pragma omp atomic
    nb_neutral_phenotypic++;

    #pragma omp atomic
	  nb_neutral_genetic++;
	}
	if (new_metabolic_error > parent_metabolic_error) {
	  
    #pragma omp atomic
    nb_neg++;

    #pragma omp critical(max_neg)
    {
	  if ((new_metabolic_error - parent_metabolic_error) > max_neg) {
	    max_neg =
	      new_metabolic_error -
	      parent_metabolic_error;
	  }
    }

    #pragma omp critical(max_fit_neg) 
    {
	  if ((new_fitness - parent_fitness) < max_fitness_neg) {
	    max_fitness_neg = new_fitness - parent_fitness;
	  }
    }

    #pragma omp atomic
	  cumul_delta_err_neg += new_metabolic_error - parent_metabolic_error;

    #pragma omp atomic
	  cumul_delta_fitness_neg += new_fitness - parent_fitness;
	}
	if (new_metabolic_error < parent_metabolic_error) {
    #pragma omp atomic
	  nb_pos++;
	  #pragma omp critical(max_pos)
    {
    if ((new_metabolic_error - parent_metabolic_error) < max_pos) {
      max_pos =
	      new_metabolic_error -
            parent_metabolic_error;
	  }
    }

    #pragma omp critical(max_fit_neg)
    {
	  if ((new_fitness - parent_fitness) > max_fitness_neg) {
      max_fitness_pos = new_fitness - parent_fitness;
	    
	  }
    }	
	  
    #pragma omp atomic
	  cumul_delta_err_pos += new_metabolic_error - parent_metabolic_error;

	  #pragma omp atomic
    cumul_delta_fitness_pos += new_fitness - parent_fitness;
	}
	
	if (full_output) {
    #pragma omp critical(full_output)
    {
	  fprintf(output_summary,"%.15e ",(new_fitness - parent_fitness));
    }
	}
      }
  }
  
  
  if (full_output) {
    fprintf(output_summary,"\n");
    }
 

  if (verbose) {
    printf("f+: %f   f0_Ph: %f   f0_Gen: %f   f-:%f\n", nb_pos, nb_neutral_phenotypic, nb_neutral_genetic, nb_neg);
  }

  if (!full_output) {
    fprintf(output_summary,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
	    nb_pos / nb_indiv, nb_neutral_phenotypic / nb_indiv, nb_neutral_genetic / nb_indiv, nb_neg / nb_indiv,
        cumul_delta_err_pos / nb_pos, cumul_delta_err_neg / nb_neg,
        max_pos, max_neg,
        cumul_delta_fitness_pos / nb_pos, cumul_delta_fitness_neg / nb_neg,
        max_fitness_pos, max_fitness_neg);
  }
}

// ============================================================================
//                            Non inline accessors
// ============================================================================


void IndivAnalysis::compute_experimental_mutagenesis(
    int32_t nb_indiv,
    int32_t mutation_type,
    std::shared_ptr<JumpingMT> prng,
    FILE* output_summary /* = nullptr*/,
    bool verbose /* = false*/,
    bool full_output /* = false*/			     ) {
  double nb_pos = 0;
  double cumul_delta_err_pos = 0;
  double cumul_delta_fitness_pos = 0;
  double nb_neg = 0;
  double cumul_delta_err_neg = 0;
  double cumul_delta_fitness_neg = 0;
  double max_pos = 0;
  double max_fitness_pos = 0;
  double max_neg = 0;
  double max_fitness_neg = 0;
  double nb_neutral_genetic = 0;
  double nb_neutral_phenotypic = 0;

    double parent_metabolic_error = dist_to_target_by_feature(METABOLISM);
  double parent_fitness = fitness();

    fprintf(output_summary,
	  "%" PRId64 " ",AeTime::time());
  
  for (int32_t i = 0; i < nb_indiv; i++) {
    Individual mutant(this, 0, prng, prng);
    // Perform one mutation of the specified type
    if (not mutant.allow_plasmids()) {
      const GeneticUnit* chromosome = &(mutant.genetic_unit_list().front());
        switch (mutation_type)
	{
	case SWITCH:{
	  chromosome->dna()->do_switch();
	  break;
	}
	case S_INS: {
	  chromosome->dna()->do_small_insertion();
	  break;
	}
	case S_DEL: {
            chromosome->dna()->do_small_deletion();
          break;
        }
        case DUPL: {
	  chromosome->dna()->do_duplication();
          break;
        }
        case DEL: {
              chromosome->dna()->do_deletion();
          break;
        }
        case TRANS: {
              chromosome->dna()->do_translocation();
          break;
        }
        case INV: {
             chromosome->dna()->do_inversion();
          break;
        }
        
	default: {
	  fprintf(stderr, "Error, unexpected mutation type\n");
	  break;
	}
      }

      
        mutant.EvaluateInContext(habitat());
    double new_metabolic_error = mutant.dist_to_target_by_feature(
        METABOLISM);
    double new_fitness = mutant.fitness();

    if (new_metabolic_error == parent_metabolic_error) {
      nb_neutral_phenotypic++;
      nb_neutral_genetic++;
      //    printf(" %d",taille_pre-taille_pos);
    }
    if (new_metabolic_error > parent_metabolic_error) {
      nb_neg++;
      if ((new_metabolic_error - parent_metabolic_error) > max_neg) {
        max_neg =
            new_metabolic_error -
            parent_metabolic_error;
      }
      if ((new_fitness - parent_fitness) < max_fitness_neg) {
	max_fitness_neg = new_fitness - parent_fitness;

      }	
      cumul_delta_err_neg += new_metabolic_error - parent_metabolic_error;
      cumul_delta_fitness_neg += new_fitness - parent_fitness;
    }
    if (new_metabolic_error < parent_metabolic_error) {
      nb_pos++;
      if ((new_metabolic_error - parent_metabolic_error) < max_pos) {
        max_pos =
            new_metabolic_error -
            parent_metabolic_error;
      }
      if ((new_fitness - parent_fitness) > max_fitness_neg) {
	max_fitness_pos = new_fitness - parent_fitness;

      }	
 
      cumul_delta_err_pos += new_metabolic_error - parent_metabolic_error;
      cumul_delta_fitness_pos += new_fitness - parent_fitness;
    }
      
 	if (full_output) {
	  fprintf(output_summary,"%.15e ",(new_fitness - parent_fitness));
	}
      }
  }

   if (full_output) {
    fprintf(output_summary,"\n");
    }

  if (verbose) {
    printf("f+: %f   f0_Ph: %f   f0_Gen: %f   f-:%f\n", nb_pos, nb_neutral_phenotypic, nb_neutral_genetic, nb_neg);
  }

  if (!full_output) {
    fprintf(output_summary,
        "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
	    nb_pos / nb_indiv, nb_neutral_phenotypic / nb_indiv, nb_neutral_genetic / nb_indiv, nb_neg / nb_indiv,
        cumul_delta_err_pos / nb_pos, cumul_delta_err_neg / nb_neg,
        max_pos, max_neg,
        cumul_delta_fitness_pos / nb_pos, cumul_delta_fitness_neg / nb_neg,
        max_fitness_pos, max_fitness_neg);
  }
}

  
// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
