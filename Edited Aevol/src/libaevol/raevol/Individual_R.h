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


#ifndef AEVOL_INDIVIDUAL_R_H_
#define  AEVOL_INDIVIDUAL_R_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <vector>
#include <set>
#include <unordered_map>
// =================================================================
//                            Project Files
// =================================================================
#include "Individual.h"
#include "Rna_R.h"
#include "Protein_R.h"
#include "Habitat_R.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class Individual_R : public virtual Individual
{
    friend class cuda_struct;


  public :

    // =================================================================
    //                             Constructors
    // =================================================================
	Individual_R() = delete;
  Individual_R(const Individual_R& other);
  Individual_R(ExpManager * exp_m,
              std::shared_ptr<JumpingMT> mut_prng,
              std::shared_ptr<JumpingMT> stoch_prng,
              std::shared_ptr<MutationParams> param_mut,
              double w_max,
              int32_t min_genome_length,
              int32_t max_genome_length,
              bool allow_plasmids,
              int32_t id,
              const char* strain_name,
              int32_t age);
	Individual_R(  Individual_R* parent, int32_t id,
                 std::shared_ptr<JumpingMT> mut_prng,
                 std::shared_ptr<JumpingMT> stoch_prng);
  Individual_R(ExpManager* exp_m, gzFile backup_file);

  static Individual_R* CreateClone(const Individual_R* dolly, int32_t id);
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Individual_R( void ) noexcept;

    // =================================================================
    //                              Accessors
    // =================================================================
    inline std::vector<Rna_R *> get_rna_list_coding( void ) const;

    // =================================================================
    //                            Public Methods
    // =================================================================
    /**
      * Main evaluation method
      */
     virtual void Evaluate(bool no_signal = true);
     /**
      * Evaluate within the provided context
      */
     virtual void EvaluateInContext(const Habitat_R& habitat, bool no_signal = true);
     virtual void EvaluateInContext(const Habitat& habitat);
     virtual void EvaluateOneAfterAnother(const Habitat_R& habitat);

     virtual void init_indiv();
     virtual void init_indiv(const Habitat_R& habitat);
     virtual void one_step( void );
     virtual void eval_step(const Habitat_R& habitat, int16_t age);
     virtual void eval_step_one_after_another( const Habitat_R& habitat, int16_t env_id );
     virtual void final_step(const Habitat_R& habitat, int16_t age);
     virtual void final_step_one_after_another( const Habitat_R& habitat, int16_t env_id );
        //virtual void reevaluate();
     //virtual void clear_everything_except_dna_and_promoters();
     //void do_transcription_translation_folding();
    // void do_transcription();
     //void do_translation();
     //void do_folding();
     //void compute_phenotype();
     //void compute_distance_to_target(const PhenotypicTarget& target);
    //void update_phenotype();
     void update_phenotype( void );

    double fitness(int env_id) { return fitness_tab_[env_id]; }
    double fitness() { return Individual::fitness(); }

    void create_csv(char* directory_name);
    void draw_phenotype(const PhenotypicTarget& target, char* directoryName, int generation);
    void clear_dist_sum() { _dist_sum = 0; };


    void    set_influences( void );
    void    update_concentrations( void );
    void    multiply_concentrations( double factor );
    int8_t  get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos );
    //int8_t  new_get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos );
    //int8_t  old_get_quadon( const GeneticUnit* gen_unit, Strand strand, int32_t pos );
    virtual void    save( gzFile backup_file ) const;
    void    clear_everything_except_dna_and_promoters();

    inline std::vector<Protein_R*> get_inherited_protein_list( void) const;
    inline void set_networked( bool networked );

    void dump_network_to_file();

    std::vector<Rna_R *> _rna_list_coding;   // Please note that these RNAs are
    // actually managed via genetic units.

    // =================================================================
    //                           Public Attributes
    // =================================================================

    std::list<Protein*> _initial_protein_list;

    std::unordered_map<int,Protein_R*> signal_list;
  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_individual_R(const ae_individual &model)
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================
    virtual void make_protein_list( void );
    virtual void make_rna_list( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    std::vector<Protein_R*> _inherited_protein_list;


    int _indiv_age;
    bool _networked;
    double _dist_sum;
    double _global_dist_sum;
    double* fitness_tab_;

};

// =====================================================================
//                          Accessors definitions
// =====================================================================
inline std::vector<Protein_R*> Individual_R::get_inherited_protein_list( void ) const
{
  return _inherited_protein_list;
}

inline std::vector<Rna_R *> Individual_R::get_rna_list_coding( void ) const
{
  return _rna_list_coding;
}

inline void Individual_R::set_networked( bool networked )
{
  _networked = networked;
}

} // namespace aevol

#endif // AEVOL_INDIVIDUAL_R_H_
