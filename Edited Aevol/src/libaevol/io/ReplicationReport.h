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


#ifndef AEVOL_REPLICATION_REPORT_H_
#define AEVOL_REPLICATION_REPORT_H_


// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>
#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "DnaReplicationReport.h"
#include "ae_enums.h"
#include "Observer.h"
#include "ObservableEvent.h"

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class Tree;
class LightTree;




class ReplicationReport : public Observer {
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ReplicationReport() = default;
#ifdef __REGUL
    ReplicationReport(Individual_R* indiv,
                      const Individual_R* parent,
                      Individual_R* donor = NULL);
#else
  ReplicationReport(Individual* indiv,
                      const Individual* parent,
                      Individual* donor = NULL);
#endif
    // Creates a completely independent copy of the original report
    ReplicationReport(const ReplicationReport& other);

#ifdef __REGUL
    ReplicationReport(gzFile tree_file, Individual_R * indiv);
#else
  ReplicationReport(gzFile tree_file, Individual * indiv);
#endif
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ReplicationReport() = default;

    // =================================================================
    //                              Accessors
    // =================================================================
#ifdef __REGUL
    inline Individual_R * indiv() const;
#else
  inline Individual * indiv() const;
#endif
    unsigned long long id() { return id_; };
    inline int32_t  genome_size() const;
    inline unsigned long long  parent_id() const;
    inline double   parent_metabolic_error() const;
    inline double   parent_secretion_error() const;
    inline int32_t  parent_genome_size() const;

    inline double   mean_align_score() const;
    inline int32_t	donor_id() const;
    inline double   donor_metabolic_error() const;
    inline int32_t  donor_genome_size() const;
    int32_t nb(MutationType t) const {
      return dna_replic_report_.nb(t);
    }

    // TODO <david.parsons@inria.fr> re-constify
    // => const DnaReplicationReport& dna_replic_report() const
    DnaReplicationReport& dna_replic_report() {
      return dna_replic_report_;
    }

#ifdef __REGUL
    void            set_indiv(Individual_R * indiv);
#else
  void            set_indiv(Individual * indiv);
#endif
    inline void     set_parent_id(unsigned long long parent_id);
    inline void     set_parent_metabolic_error(double parent_metabolic_error);
    inline void     set_parent_secretion_error(double parent_secretion_error);
    inline void     set_parent_genome_size(int32_t parent_genome_size);
    inline void     set_donor_id(int32_t donor_id);
    inline void     set_donor_metabolic_error(double donor_metabolic_error);
    inline void     set_donor_secretion_error(double donor_secretion_error);
    inline void     set_donor_genome_size(int32_t donor_genome_size);

    // =================================================================
    //                            Public Methods
    // =================================================================
#ifdef __REGUL
    void init(Tree* tree, Individual_R* offspring, Individual_R* parent, int indiv_id, int parent_id);
    void init(LightTree* tree, Individual_R* offspring, Individual_R* parent, int indiv_id, int parent_id);
    void signal_end_of_replication(Individual_R* indiv);
#else
  void init(Tree* tree, Individual* offspring, Individual* parent, int indiv_id, int parent_id);
    void init(LightTree* tree, Individual* offspring, Individual* parent, int indiv_id, int parent_id);
    void signal_end_of_replication(Individual* indiv);
#endif

    void init(Tree* tree,
              Individual_7* offspring,
              Individual_7* parent, int indiv_id, int parent_id, bool remote = false, int rank = -1);
    void init(LightTree* tree,
              Individual_7* offspring,
              Individual_7* parent, int indiv_id, int parent_id);
    void signal_end_of_replication(Individual_7* indiv);

    void signal_end_of_generation();
    void write_to_tree_file(gzFile tree_file);


  void update(Observable& o, ObservableEvent e, void* arg) override;

  Individual_7* simd_indiv_ = nullptr;

#ifdef HAVE_MPI
  int32_t rank_ = -1;
#endif

    unsigned long long id_ = 0;
    unsigned long long parent_id_ = 0;

  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
#ifdef __REGUL
    Individual_R* indiv_ = nullptr;
#else
Individual* indiv_ = nullptr;
#endif

    // int32_t rank_ = -1;

    int32_t genome_size_ = -1;
    double metabolic_error_ = -1;
    int16_t nb_genes_activ_ = -1;
    int16_t nb_genes_inhib_ = -1;
    int16_t nb_non_fun_genes_ = -1;
    int16_t nb_coding_RNAs_ = -1;
    int16_t nb_non_coding_RNAs_ = -1;

    // List of each genetic unit's replication report
    DnaReplicationReport dna_replic_report_;

    double parent_metabolic_error_ = -1;
    double parent_secretion_error_ = -1;
    int32_t parent_genome_size_ = -1;

    int32_t	donor_id_ = -1;
    double donor_metabolic_error_ = -1;
    double donor_secretion_error_ = -1;
    int32_t donor_genome_size_ = -1;

    // CK: I think that the attributes below are obsolete
    // (HT events are now stored in the ae_dna_replic_reports)

    double mean_align_score_;

    // bool remote_;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
#ifdef __REGUL
inline Individual_R *ReplicationReport::indiv() const
#else
inline Individual *ReplicationReport::indiv() const
#endif
{
  return indiv_;
}

inline int32_t ReplicationReport::genome_size() const
{
  return genome_size_;
}

unsigned long long ReplicationReport::parent_id() const
{
  return parent_id_;
}

double ReplicationReport::parent_metabolic_error() const
{
  return parent_metabolic_error_;
}

double ReplicationReport::parent_secretion_error() const
{
  return parent_secretion_error_;
}

int32_t ReplicationReport::parent_genome_size() const
{
  return parent_genome_size_;
}

inline int32_t	ReplicationReport::donor_id() const
{
  return donor_id_;
}

inline double   ReplicationReport::donor_metabolic_error() const
{
  return donor_metabolic_error_;
}

inline int32_t  ReplicationReport::donor_genome_size() const
{
  return donor_genome_size_;
}



inline double ReplicationReport::mean_align_score() const
{
  return mean_align_score_;
}

#ifdef __REGUL
inline void ReplicationReport::set_indiv(Individual_R * indiv)
#else
inline void ReplicationReport::set_indiv(Individual * indiv)
#endif
{
  indiv_ = indiv;
}

void ReplicationReport::set_parent_id(unsigned long long parent_id)
{
  parent_id_ = parent_id;
}

void ReplicationReport::set_parent_metabolic_error(double parent_metabolic_error)
{
  parent_metabolic_error_ = parent_metabolic_error;
}

void ReplicationReport::set_parent_secretion_error(double parent_secretion_error)
{
  parent_secretion_error_ = parent_secretion_error;
}

void ReplicationReport::set_parent_genome_size(int32_t parent_genome_size)
{
  parent_genome_size_ = parent_genome_size;
}

inline void ReplicationReport::set_donor_id(int32_t donor_id)
{
  donor_id_ = donor_id;
}

inline void  ReplicationReport::set_donor_metabolic_error(double donor_metabolic_error)
{
  donor_metabolic_error_ = donor_metabolic_error;
}

inline void ReplicationReport::set_donor_secretion_error(double donor_secretion_error)
{
  donor_secretion_error_ = donor_secretion_error;
}

inline void ReplicationReport::set_donor_genome_size(int32_t donor_genome_size)
{
  donor_genome_size_ = donor_genome_size;
}

} // namespace aevol

#endif // AEVOL_REPLICATION_REPORT_H_
