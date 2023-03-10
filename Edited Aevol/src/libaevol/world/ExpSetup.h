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

#ifndef AEVOL_EXP_SETUP_H_
#define AEVOL_EXP_SETUP_H_

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdlib>
#include <set>

#include "Selection.h"
#include "Stats.h"
#include "Logging.h"
#include "Tree.h"
#include "Dump.h"
#include "JumpingMT.h"

namespace aevol {

// ===========================================================================
//                          Class declarations
// ===========================================================================
class ParamLoader;

class ExpSetup {
  friend class ExpManager;

  public:
  // =======================================================================
  //                             Constructors
  // =======================================================================
  ExpSetup() = delete;
  ExpSetup(const ExpSetup&) = delete;
  ExpSetup(ExpManager* exp_m);

  // =======================================================================
  //                             Destructors
  // =======================================================================
  virtual ~ExpSetup() noexcept;

  // =======================================================================
  //                         Accessors: getters
  // =======================================================================
  inline int get_fuzzy_flavor( void ) const;
  inline int get_simd_metadata_flavor( void ) const;

  // ----------------------------------------------------- Selection context
  Selection * sel() const { return sel_; }

  // --------------------------------------------------------------- Transfer
  bool with_HT() const { return with_HT_; };
  bool repl_HT_with_close_points( void ) const { return repl_HT_with_close_points_; };
  double HT_ins_rate() const { return HT_ins_rate_; };
  double HT_repl_rate() const { return HT_repl_rate_; };
  double repl_HT_detach_rate() const { return repl_HT_detach_rate_; };

  // --------------------------------------------------------------- Plasmids
  // See comments in ExpManager.h on how plasmids are handled
  bool   with_plasmids() const { return with_plasmids_; }
  double prob_plasmid_HT() const { return prob_plasmid_HT_; }
  double tune_donor_ability() const { return tune_donor_ability_; }
  double tune_recipient_ability() const { return tune_recipient_ability_; }
  bool   donor_cost() const { return donor_cost_; }
  bool   recipient_cost() const { return recipient_cost_; }
  bool   swap_GUs() const { return swap_GUs_; }

  // -------------------------------------------------------------- Secretion
  bool   with_secretion() const { return with_secretion_; }
  double secretion_contrib_to_fitness() const { return secretion_contrib_to_fitness_; }
  double secretion_cost() const { return secretion_cost_; }

  #ifdef BASE_4
  // -------------------------------------------------------------- Terminators
  int8_t terminator_polya_sequence_length() const { return terminator_polya_sequence_length_; }

  // -------------------------------------------------------------- MWH bases configuration
  const int8_t * aa_base_m() const { return aa_base_m_; }
  const int8_t * aa_base_w() const { return aa_base_w_; }
  const int8_t * aa_base_h() const { return aa_base_h_; }

  int8_t aa_base_m_size() const { return aa_base_m_size_; }
  int8_t aa_base_w_size() const { return aa_base_w_size_; }
  int8_t aa_base_h_size() const { return aa_base_h_size_; }
  #endif
  // ------------------------------------------------------------ Parameter for SIMD
  std::shared_ptr<MutationParams> mut_params() { return mut_params_;}
  int32_t min_genome_length() const { return min_genome_length_; }
  int32_t max_genome_length() const { return max_genome_length_; }


#ifdef __REGUL
    inline bool   get_with_heredity( void ) const;
    inline double get_degradation_rate( void ) const;
    inline int    get_nb_degradation_step( void ) const;
    inline double get_protein_presence_limit( void ) const;

    inline double get_hill_shape( void ) const;
    inline double get_hill_shape_n( void ) const;
    inline double get_hill_shape_theta( void ) const;

    inline int get_nb_indiv_age( void ) const;
    std::set<int>* get_list_eval_step( void ) const {
      return _list_eval_step;
    }
#endif

#ifdef HAVE_MPI
  int nb_rank() const { return nb_rank_; }
  int32_t global_pop_size() const { return global_pop_size_; }
  int32_t global_grid_width() const { return global_grid_width_; }
  int32_t global_grid_height() const { return global_grid_height_; }

  int32_t rank_width() const { return rank_width_; }
  int32_t rank_height() const { return rank_height_; }
#endif
  // =======================================================================
  //                         Accessors: setters
  // =======================================================================
  inline void set_fuzzy_flavor( int fuzzy_flavor );
    inline void set_simd_metadata_flavor(int metadata_flavor );
  // --------------------------------------------------------------- Transfer
  void set_with_HT(bool with_HT) { with_HT_ = with_HT; }
  void set_repl_HT_with_close_points(bool repl_HT_with_close_points) { repl_HT_with_close_points_ = repl_HT_with_close_points; }
  void set_HT_ins_rate(double HT_ins_rate) { HT_ins_rate_ = HT_ins_rate; }
  void set_HT_repl_rate(double HT_repl_rate) { HT_repl_rate_ = HT_repl_rate; }
  void set_repl_HT_detach_rate(double repl_HT_detach_rate) { repl_HT_detach_rate_ = repl_HT_detach_rate; }

  // --------------------------------------------------------------- Plasmids
  void set_with_plasmids(bool with_p) { with_plasmids_ = with_p; }
  void set_prob_plasmid_HT(double prob_p_HT) { prob_plasmid_HT_ = prob_p_HT; }
  void set_tune_donor_ability(double tune_donor_ability) { tune_donor_ability_ = tune_donor_ability; }
  void set_tune_recipient_ability(double tune_recipient_ability) { tune_recipient_ability_ = tune_recipient_ability; }
  void set_donor_cost(double donor_cost) { donor_cost_ = donor_cost; }
  void set_recipient_cost(double recipient_cost) { recipient_cost_ = recipient_cost; }
  void set_swap_GUs(bool swap_GUs) { swap_GUs_ = swap_GUs; }

  // -------------------------------------------------------------- Secretion
  void set_with_secretion(bool with_secretion) { with_secretion_ = with_secretion; }
  void set_secretion_contrib_to_fitness(double secretion_contrib) { secretion_contrib_to_fitness_ = secretion_contrib; }
  void set_secretion_cost(double secretion_cost) { secretion_cost_ = secretion_cost; }

  #ifdef BASE_4
  // -------------------------------------------------------------- Terminators
  void set_terminator_polya_sequence_length(int8_t sequence_length) { terminator_polya_sequence_length_ = sequence_length; }

  // -------------------------------------------------------------- MWH bases configuration
  void set_aa_base_m(int8_t *base_m);
  void set_aa_base_w(int8_t *base_w);
  void set_aa_base_h(int8_t *base_h);

  void set_aa_base_m_size(int8_t size) { aa_base_m_size_ = size; }
  void set_aa_base_w_size(int8_t size) { aa_base_w_size_ = size; }
  void set_aa_base_h_size(int8_t size) { aa_base_h_size_ = size; }
  #endif

    // ------------------------------------------------------------ Parameter for SIMD
    void set_min_genome_length(int32_t min_genome_length) { min_genome_length_ = min_genome_length; }
    void set_max_genome_length(int32_t max_genome_length) { max_genome_length_ = max_genome_length; }

    void set_mutation_parameters(std::shared_ptr<MutationParams> mut_params) { mut_params_ = mut_params; }

#ifdef __REGUL
    inline void set_with_heredity( bool with_heredity );
    inline void set_degradation_rate( double degradation_rate );
    inline void set_nb_degradation_step( int nb_degradation_step );
    inline void set_protein_presence_limit( double protein_presence_limit );

    inline void set_hill_shape( double hill_shape );
    inline void set_hill_shape_n( double hill_shape_n );
    inline void set_hill_shape_theta( double hill_shape_theta );

    inline void set_nb_indiv_age( int nb_indiv_age );

    inline void set_list_eval_step(std::set<int> list_eval_step);
#endif

#ifdef HAVE_MPI
void set_nb_rank(int nb_rank) { nb_rank_ = nb_rank; }
void set_global_pop_size(int global_pop_size) { global_pop_size_ = global_pop_size; }
void set_global_grid_width(int global_grid_width) { global_grid_width_ = global_grid_width; }
void set_global_grid_height(int global_grid_height) { global_grid_height_ = global_grid_height; }

void set_rank_width(int rank_width) { rank_width_ = rank_width; }
void set_rank_height(int rank_height) { rank_height_ = rank_height; }
#endif

  // =======================================================================
  //                            Public Methods
  // =======================================================================
  void write_setup_file(gzFile exp_setup_file) const;
  void save(gzFile backup_file) const;
  void load(gzFile setup_file, gzFile backup_file, bool verbose);
  /// Make the individuals reproduce
  void step_to_next_generation() { sel_->step_to_next_generation(); }
#ifdef __REGUL
    // Regulation
    void     init_binding_matrix( bool random_binding_matrix, double binding_zeros_percentage,
    		       std::shared_ptr<JumpingMT> prng);

    void     read_binding_matrix_from_backup(gzFile binding_matrix_file);
    void     write_binding_matrix_to_backup(gzFile binding_matrix_file) const;

    void     write_binding_matrix_to_file(FILE* binding_matrix_file) const;
    void     print_binding_matrix( void );

    ProteinConcentration get_binding_matrix( int row, int column ) const;
#endif

  // =======================================================================
  //                           Public Attributes
  // =======================================================================

#ifdef __REGUL
    ProteinConcentration  _binding_matrix[MAX_QUADON][MAX_CODON]; //
#endif

 protected :
  // =======================================================================
  //                           Protected Methods
  // =======================================================================
  virtual void display() {};

  // =======================================================================
  //                          Protected Attributes
  // =======================================================================
  ExpManager* exp_m_;

  int fuzzy_flavor_ = 0;
  int simd_metadata_flavor_;
  // ----------------------------------------------------- Selection context
  Selection* sel_;

  // --------------------------------------------------- Transfer parameters
  bool   with_HT_ = false;
  bool   repl_HT_with_close_points_ = false;
  double HT_ins_rate_;
  double HT_repl_rate_;
  double repl_HT_detach_rate_;

  // --------------------------------------------------- Plasmids parameters
  bool   with_plasmids_ = false;
  double prob_plasmid_HT_; // Base transfer ability independent of evolvable donor and recipient ability
  double tune_donor_ability_; // How much the individuals can tune their ability to send plasmids
  double tune_recipient_ability_; // How much the individuals can tune their ability to receive plasmids
  double donor_cost_;
  double recipient_cost_;
  bool   swap_GUs_; // Whether plasmid HT is uni- or bidirectional

  // -------------------------------------------------- Secretion parameters
  bool   with_secretion_;
  double secretion_contrib_to_fitness_;
  double secretion_cost_;

  #ifdef BASE_4
  // -------------------------------------------------- Terminator parameters
  int8_t terminator_polya_sequence_length_;

  // -------------------------------------------------- MWH bases configuration
  int8_t aa_base_m_[NB_AMINO_ACIDS];
  int8_t aa_base_w_[NB_AMINO_ACIDS];
  int8_t aa_base_h_[NB_AMINO_ACIDS];

  int8_t aa_base_m_size_;
  int8_t aa_base_w_size_;
  int8_t aa_base_h_size_;
  #endif

    // Mutation rates etc...
    std::shared_ptr<MutationParams> mut_params_;

    int32_t min_genome_length_;
    int32_t max_genome_length_;

#ifdef __REGUL
    // Binding matrix


    bool    _with_heredity;
    double  _protein_presence_limit;

    double  _degradation_rate;
    int     _nb_degradation_step;

    double _hill_shape_n;
    double _hill_shape;
    double _hill_shape_theta;

    int    _nb_indiv_age;

    std::set<int>* _list_eval_step;
#endif

  #ifdef HAVE_MPI
  int nb_rank_ = 1;
  
  // Define global grid size
  int32_t global_pop_size_;
  int32_t global_grid_width_;
  int32_t global_grid_height_;

  int32_t rank_width_;
  int32_t rank_height_;
  #endif

};


// =====================================================================
//                           Getters' definitions
// =====================================================================

inline int ExpSetup::get_fuzzy_flavor( void ) const
{
  return fuzzy_flavor_;
}

inline int ExpSetup::get_simd_metadata_flavor( void ) const
{
   return simd_metadata_flavor_;
}

#ifdef __REGUL
inline bool ExpSetup::get_with_heredity( void ) const
{
  return _with_heredity;
}

inline double ExpSetup::get_degradation_rate( void ) const
{
  return _degradation_rate;
}


inline int ExpSetup::get_nb_degradation_step( void ) const
{
  return _nb_degradation_step;
}

inline double ExpSetup::get_protein_presence_limit( void ) const
{
  return _protein_presence_limit;
}

inline double ExpSetup::get_hill_shape( void ) const
{
  return _hill_shape;
}

inline double ExpSetup::get_hill_shape_n( void ) const
{
  return _hill_shape_n;
}

inline double ExpSetup::get_hill_shape_theta( void ) const
{
  return _hill_shape_theta;
}

inline int ExpSetup::get_nb_indiv_age( void ) const
{
  return _nb_indiv_age;
}

#endif

// =====================================================================
//                           Setters' definitions
// =====================================================================
// --------------------------------------------------------------- Transfer
inline void ExpSetup::set_fuzzy_flavor( int fuzzy_flavor )
{
  fuzzy_flavor_ = fuzzy_flavor;
}

inline void ExpSetup::set_simd_metadata_flavor(int metadata_flavor)
{
  simd_metadata_flavor_ = metadata_flavor;
}

#ifdef __REGUL
inline void ExpSetup::set_with_heredity( bool with_heredity )
{
  _with_heredity = with_heredity;
}

inline void ExpSetup::set_degradation_rate( double degradation_rate )
{
  _degradation_rate = degradation_rate;
}

inline void ExpSetup::set_nb_degradation_step( int degradation_step )
{
  _nb_degradation_step = degradation_step;
}

inline void ExpSetup::set_protein_presence_limit( double protein_presence_limit )
{
  _protein_presence_limit = protein_presence_limit;
}

inline void ExpSetup::set_hill_shape( double hill_shape )
{
  _hill_shape = hill_shape;
}

inline void ExpSetup::set_hill_shape_theta( double hill_shape_theta )
{
  _hill_shape_theta = hill_shape_theta;
}

inline void ExpSetup::set_hill_shape_n( double hill_shape_n )
{
  _hill_shape_n = hill_shape_n;
}


inline void ExpSetup::set_nb_indiv_age( int nb_indiv_age )
{
  _nb_indiv_age = nb_indiv_age;
}

inline void ExpSetup::set_list_eval_step( std::set<int> list_eval_step )
{
  _list_eval_step = new std::set<int>(list_eval_step);
}
#endif

// =====================================================================
//                       functions' definition
// =====================================================================

#ifdef __REGUL
inline ProteinConcentration ExpSetup::get_binding_matrix( int row, int column ) const
{
  if (row >= MAX_QUADON || column >= MAX_CODON) {
    printf("OUT OF BOUND BINDING R %d C %d\n",row,column);
    exit(-1);
  }

  return _binding_matrix[row][column];
}
#endif

} // namespace aevol

#endif // AEVOL_EXP_SETUP_H_
