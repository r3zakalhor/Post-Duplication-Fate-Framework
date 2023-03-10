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



// =================================================================
//                            Project Files

#include "ReplicationReport.h"

#include "7/Dna_7.h"
#include "7/Individual_7.h"
#include "7/ExpManager_7.h"
#include "AeTime.h"
#include "DnaReplicationReport.h"
#include "ExpManager.h"
#include "GridCell.h"
#include "Individual.h"
#include "Mutation.h"
#include "Observable.h"
#include "Tree.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                         Class ReplicationReport                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
#ifdef __REGUL
ReplicationReport::ReplicationReport(Individual_R* indiv,
                                     const Individual_R* parent,
                                     Individual_R* donor /*= NULL*/)
#else
ReplicationReport::ReplicationReport(Individual* indiv,
                                     const Individual* parent,
                                     Individual* donor /*= NULL*/)
#endif
{
  indiv_ = indiv;

  id_   = indiv->id();
  // rank_ = indiv->rank();

  parent_id_ = parent->id();
  // donor_id_ is set further down

  genome_size_        = 0;
  metabolic_error_    = 0.0;
  nb_genes_activ_     = 0;
  nb_genes_inhib_     = 0;
  nb_non_fun_genes_   = 0;
  nb_coding_RNAs_     = 0;
  nb_non_coding_RNAs_ = 0;

  parent_metabolic_error_ = parent->dist_to_target_by_feature(METABOLISM);
  parent_secretion_error_ = parent->dist_to_target_by_feature(SECRETION);
  parent_genome_size_     = parent->total_genome_size();
  mean_align_score_       = 0.0;

  if (donor == NULL)
  {
    donor_id_               = -1;
    donor_metabolic_error_  = 0.0;
    donor_secretion_error_	= 0.0;
    donor_genome_size_      = 0;
  }
  else
  {
    donor_id_              = donor->id();
    donor_metabolic_error_ = donor->dist_to_target_by_feature(METABOLISM);
    donor_secretion_error_ = donor->dist_to_target_by_feature(SECRETION);
    donor_genome_size_     = donor->total_genome_size();
  }
}


// Creates an independent copy of the original report
ReplicationReport::ReplicationReport(const ReplicationReport& other) :
    dna_replic_report_(other.dna_replic_report_)
{
  parent_id_  = other.parent_id_;
  donor_id_   = other.donor_id_;

  id_   = other.id_;
  // #ifdef HAVE_MPI
  // rank_ = other.rank_;
  // #endif

  genome_size_        = other.genome_size_;
  metabolic_error_    = other.metabolic_error_;
  nb_genes_activ_     = other.nb_genes_activ_;
  nb_genes_inhib_     = other.nb_genes_inhib_;
  nb_non_fun_genes_   = other.nb_non_fun_genes_;
  nb_coding_RNAs_     = other.nb_coding_RNAs_;
  nb_non_coding_RNAs_ = other.nb_non_coding_RNAs_;

  parent_metabolic_error_ = other.parent_metabolic_error_;
  parent_secretion_error_ = other.parent_secretion_error_;
  donor_metabolic_error_  = other.donor_metabolic_error_;
  donor_secretion_error_  = other.donor_secretion_error_;
  parent_genome_size_     = other.parent_genome_size_;
  donor_genome_size_      = other.donor_genome_size_;
  mean_align_score_       = other.mean_align_score_;
  // remote_                 = other.remote_;
}


#ifdef __REGUL
ReplicationReport::ReplicationReport(gzFile tree_file, Individual_R* indiv)
#else
ReplicationReport::ReplicationReport(gzFile tree_file, Individual* indiv)
#endif
{
  indiv_ = indiv;

  gzread(tree_file, &id_,        sizeof(id_));
  // printf("REPREP ID %d\n",id_);
  // #ifdef HAVE_MPI
  // gzread(tree_file, &rank_,      sizeof(rank_));
  // #endif
  gzread(tree_file, &parent_id_, sizeof(parent_id_));
  gzread(tree_file, &donor_id_,  sizeof(donor_id_));

  gzread(tree_file, &genome_size_,         sizeof(genome_size_));
  gzread(tree_file, &metabolic_error_,     sizeof(metabolic_error_));
  gzread(tree_file, &nb_genes_activ_,      sizeof(nb_genes_activ_));
  gzread(tree_file, &nb_genes_inhib_,      sizeof(nb_genes_inhib_));
  gzread(tree_file, &nb_non_fun_genes_,    sizeof(nb_non_fun_genes_));
  gzread(tree_file, &nb_coding_RNAs_,      sizeof(nb_coding_RNAs_));
  gzread(tree_file, &nb_non_coding_RNAs_,  sizeof(nb_non_coding_RNAs_));
  // #ifdef HAVE_MPI
  // gzread(tree_file, &remote_,  sizeof(remote_));
  // #endif
  

  // printf("Indiv %d (P %d) Remote ? %d (Rank %d)\n",id_,parent_id_,remote_,rank_);

   dna_replic_report_.read_from_tree_file(tree_file);

   dna_replic_report_.compute_stats();

  parent_metabolic_error_ = -1;
  parent_secretion_error_ = -1;
  donor_metabolic_error_  = -1;
  parent_genome_size_     = -1;
  donor_genome_size_      = -1;
  mean_align_score_       = 0.0;
}


// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================
/**
 * Set the individual corresponding to this replication report and
 * the characteristics of its parent
 *
 * This should be called as soon as a replication is started (just after calling
 * the offspring constructor and before doing the mutations)
 */
#ifdef __REGUL
void ReplicationReport::init(Tree* tree, Individual_R* offspring, Individual_R* parent, int indiv_id, int parent_id)
#else
void ReplicationReport::init(Tree* tree, Individual* offspring, Individual* parent, int indiv_id, int parent_id)
#endif
{

    dna_replic_report_.clear();
  indiv_ = offspring;

  id_ = indiv_id;
  parent_id_ = parent_id;

  genome_size_        = 0;
  metabolic_error_    = 0.0;
  nb_genes_activ_     = 0;
  nb_genes_inhib_     = 0;
  nb_non_fun_genes_   = 0;
  nb_coding_RNAs_     = 0;
  nb_non_coding_RNAs_ = 0;

  parent_metabolic_error_ = parent->dist_to_target_by_feature(METABOLISM);
  parent_secretion_error_ = parent->dist_to_target_by_feature(SECRETION);
  parent_genome_size_     = parent->total_genome_size();
  mean_align_score_       = 0.0;

  // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
//  indiv_->addObserver(this, MUTATION);
//  indiv_->addObserver(tree, END_REPLICATION);
}

void ReplicationReport::init(Tree* tree,
                             Individual_7* offspring,
                             Individual_7* parent, int indiv_id,
                                int parent_id, bool remote, int rank)
{
        dna_replic_report_.clear();

      simd_indiv_ = offspring;

      id_ = (unsigned long long) indiv_id;
      parent_id_ = (unsigned long long) parent_id;

      // printf("%d -- %d -- GID %d PID %d\n",AeTime::time(),tree->exp_m_->rank(),indiv_id,parent_id);

      // rank_ = rank;

      genome_size_        = 0;
      metabolic_error_    = 0.0;
      nb_genes_activ_     = 0;
      nb_genes_inhib_     = 0;
      nb_non_fun_genes_   = 0;
      nb_coding_RNAs_     = 0;
      nb_non_coding_RNAs_ = 0;

#ifndef HAVE_MPI
      parent_metabolic_error_ = parent->metaerror;
      parent_secretion_error_ = 0.0;
      parent_genome_size_     = parent->dna_->length();
      mean_align_score_       = 0.0;
#endif

      donor_id_ = -1;
      donor_metabolic_error_ = -1;
      donor_secretion_error_ = -1;
      donor_genome_size_ = -1;

      dna_replic_report_.clear();

      //remote_ = remote;

      // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
      //simd_indiv_->addObserver(this, MUTATION);
      //simd_indiv_->addObserver(tree, END_REPLICATION);
}

#ifdef __REGUL
    void ReplicationReport::init(LightTree* tree, Individual_R* offspring, Individual_R* parent, int indiv_id, int parent_id)
#else
void ReplicationReport::init(LightTree* tree, Individual* offspring, Individual* parent, int indiv_id, int parent_id)
#endif
    {

        indiv_ = offspring;

        id_ = indiv_id;
        parent_id_ = parent_id;

        genome_size_        = 0;
        metabolic_error_    = 0.0;
        nb_genes_activ_     = 0;
        nb_genes_inhib_     = 0;
        nb_non_fun_genes_   = 0;
        nb_coding_RNAs_     = 0;
        nb_non_coding_RNAs_ = 0;

        parent_metabolic_error_ = parent->dist_to_target_by_feature(METABOLISM);
        parent_secretion_error_ = parent->dist_to_target_by_feature(SECRETION);
        parent_genome_size_     = parent->total_genome_size();
        mean_align_score_       = 0.0;

        // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
        indiv_->addObserver(this, MUTATION);
        indiv_->addObserver(tree, END_REPLICATION);
    }

    void ReplicationReport::init(LightTree* tree,
                                 Individual_7* offspring,
                                 Individual_7* parent, int indiv_id,
                                 int parent_id)
    {

        simd_indiv_ = offspring;

        id_ = (unsigned long long) indiv_id;
        parent_id_ = (unsigned long long) parent_id;
        

        // rank_ = 0;

        genome_size_        = 0;
        metabolic_error_    = 0.0;
        nb_genes_activ_     = 0;
        nb_genes_inhib_     = 0;
        nb_non_fun_genes_   = 0;
        nb_coding_RNAs_     = 0;
        nb_non_coding_RNAs_ = 0;

        parent_metabolic_error_ = parent->metaerror;
        parent_secretion_error_ = 0.0;
        parent_genome_size_     = parent->dna_->length();
        mean_align_score_       = 0.0;

        // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
        simd_indiv_->addObserver(this, MUTATION);
        simd_indiv_->addObserver(tree, END_REPLICATION);
    }

/**
 * Method called at the end of the replication of an individual.
 * Actions such as finalize the calculation of average values can be done here.
 */
#ifdef __REGUL
void ReplicationReport::signal_end_of_replication(Individual_R* indiv) {
#else
  void ReplicationReport::signal_end_of_replication(Individual* indiv) {
#endif
  // TODO <david.parsons@inria.fr> tmp patch
  if (indiv_ == NULL) indiv_ = indiv;

  // Retrieve data from the individual
  genome_size_        = indiv_->total_genome_size();
  metabolic_error_    = indiv_->dist_to_target_by_feature(METABOLISM);
  nb_genes_activ_     = indiv_->nb_genes_activ();
  nb_genes_inhib_     = indiv_->nb_genes_inhib();
  nb_non_fun_genes_   = indiv_->nb_functional_genes();
  nb_coding_RNAs_     = indiv_->nb_coding_RNAs();
  nb_non_coding_RNAs_ = indiv_->nb_non_coding_RNAs();
}


void ReplicationReport::signal_end_of_replication(Individual_7* indiv) {
      // TODO <david.parsons@inria.fr> tmp patch
      if (simd_indiv_ == NULL) simd_indiv_ = indiv;

      // Retrieve data from the individual
      genome_size_        = simd_indiv_->dna_->length();
      metabolic_error_    = simd_indiv_->metaerror;
      nb_genes_activ_     = simd_indiv_->nb_genes_activ;
      nb_genes_inhib_     = simd_indiv_->nb_genes_inhib;
      nb_non_fun_genes_   = simd_indiv_->nb_func_genes;
      nb_coding_RNAs_     = simd_indiv_->nb_coding_RNAs;
      nb_non_coding_RNAs_ = simd_indiv_->nb_non_coding_RNAs;
}
/**
 * Method called at the end of a generation.
 * Actions such as update the individuals' ranks can be done here.
 */
void ReplicationReport::signal_end_of_generation() {
    // if (!ExpManager_7::standalone_simd) {
    //     rank_ = indiv_->rank();
    // }
}

void ReplicationReport::write_to_tree_file(gzFile tree_file)
{
  // Store individual identifiers and rank
      // printf("RR ID %d PID %d\n",id_,parent_id_);

  gzwrite(tree_file, &id_,         sizeof(id_));

  //   int32_t rankx = -1;
  // if (ExpManager_7::standalone_simd) {
  //     rankx = 0;
  // } else {
  //     rankx = rank_;
  // }
// #ifdef HAVE_MPI
//   gzwrite(tree_file, &rank_,       sizeof(rank_));
// #endif

  gzwrite(tree_file, &parent_id_,  sizeof(parent_id_));
  gzwrite(tree_file, &donor_id_,   sizeof(donor_id_));

  gzwrite(tree_file, &genome_size_,         sizeof(genome_size_));
  gzwrite(tree_file, &metabolic_error_,     sizeof(metabolic_error_));
  gzwrite(tree_file, &nb_genes_activ_,      sizeof(nb_genes_activ_));
  gzwrite(tree_file, &nb_genes_inhib_,      sizeof(nb_genes_inhib_));
  gzwrite(tree_file, &nb_non_fun_genes_,    sizeof(nb_non_fun_genes_));
  gzwrite(tree_file, &nb_coding_RNAs_,      sizeof(nb_coding_RNAs_));
  gzwrite(tree_file, &nb_non_coding_RNAs_,  sizeof(nb_non_coding_RNAs_));
  // #ifdef HAVE_MPI
  // gzwrite(tree_file, &remote_,  sizeof(remote_));
  // #endif

   dna_replic_report_.write_to_tree_file(tree_file);
}


// =================================================================
//                           Protected Methods
// =================================================================





// =================================================================
//                          Non inline accessors
// =================================================================
void ReplicationReport::update(Observable& o, ObservableEvent e, void* arg) {
//        printf("Receive ??? events\n");
  switch (e) {
    case MUTATION :
        //printf("Receive mutation events\n");
//#pragma omp critical
      //{
          dna_replic_report_.add_mut(reinterpret_cast<Mutation *>(arg));
      //}
      break;
    default :
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
  }
}
} // namespace aevol
