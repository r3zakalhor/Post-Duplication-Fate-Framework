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

#include <cinttypes>

#ifndef AEVOL_MACROS_H_
#define AEVOL_MACROS_H_

#ifdef BASE_2
constexpr int8_t NB_BASE = 2; // WARNING :  A lot of stuff has been optimized for binary genomes
                              //            Changing the value of NB_BASE implies verifying the existing code
                              //            and make changes where necessary
#elif BASE_4
constexpr int8_t NB_BASE = 4; // WARNING :  A lot of stuff has been optimized for binary genomes
                              //            Changing the value of NB_BASE implies verifying the existing code
                              //            and make changes where necessary

                              // Note : opposed bases only have a different LSB
constexpr char   BASE_A = '0';
constexpr char   BASE_T = '1';
constexpr char   BASE_C = '2';
constexpr char   BASE_G = '3';

inline bool is_complementary_base(char base1, char base2) {
  // opposite bases only differ in their LSB
  return ( (base1 & 0xFE) == (base2 & 0xFE) ) && ( (base1 & 1) != (base2 & 1) );
}

inline char get_complementary_base(char base) {
  // If LSB is 0, set to 1, else set to 0
  return base ^ 1;
}
#endif
// NB The following strings are not easily replaced with `constexpr
// const char*` because they are meant to be concatenated by the
// preprocessor.

// Backup directories and file name formats
#define TIMESTEP_FORMAT "%09" PRId64
// Experimental Setup
#define EXP_S_DIR                 "exp_setup"
#define EXP_S_FNAME_BASE          "exp_setup_" TIMESTEP_FORMAT
#ifdef HAVE_MPI
#define EXP_S_FNAME_FORMAT        EXP_S_DIR "/" EXP_S_FNAME_BASE "_%d.ae"
#else
#define EXP_S_FNAME_FORMAT        EXP_S_DIR "/" EXP_S_FNAME_BASE ".ae"
#endif
#define EXP_S_CONST_FNAME_BASE        "exp_setup_const"
#ifdef HAVE_MPI
#define EXP_S_CONST_FNAME_FORMAT      EXP_S_DIR "/" EXP_S_CONST_FNAME_BASE "_%d.ae"
#else
#define EXP_S_CONST_FNAME_FORMAT      EXP_S_DIR "/" EXP_S_CONST_FNAME_BASE ".ae"
#endif
// Output Profile
#define OUT_P_DIR                 "output_profile"
#define OUT_P_FNAME_BASE          "output_profile"
#ifdef HAVE_MPI
#define OUT_P_FNAME_FORMAT        OUT_P_DIR "/" OUT_P_FNAME_BASE "_%d.ae"
#else
#define OUT_P_FNAME_FORMAT        OUT_P_DIR "/" OUT_P_FNAME_BASE ".ae"
#endif
#define OUT_P_CUR_FNAME           "output_profile.ae"
// Spatial Structure
#define WORLD_DIR             "world"
#ifdef HAVE_MPI
#define WORLD_FNAME_BASE      "world_" TIMESTEP_FORMAT "_%d"
#else
#define WORLD_FNAME_BASE      "world_" TIMESTEP_FORMAT
#endif
#define WORLD_FNAME_FORMAT    WORLD_DIR "/" WORLD_FNAME_BASE".ae"
// Stats
#define STATS_DIR   "stats"
// Tree
#define TREE_DIR    "tree"
#define LIGHTTREE_DIR "lightTree"
// LightTree file in newick format
constexpr auto LIGHTTREE_FILE_NAME = "light_tree.txt";
// Last gener file
constexpr auto LAST_GENER_FNAME = "last_gener.txt";
// Best last organism file
#define BEST_LAST_ORG_FNAME "best_last_org.txt"

#define FIXED_POPULATION_SIZE // Some calculation can be spared if we know that the size of the population is fixed

constexpr int8_t PROM_SIZE = 22;

#ifdef BASE_2
constexpr auto PROM_SEQ = "0101011001110010010110";

constexpr const char* SHINE_DAL_SEQ = "011011";

constexpr int8_t PROM_MAX_DIFF  = 4;

#elif BASE_4
constexpr auto PROM_SEQ = "3322022010310333112013";


constexpr const char* SHINE_DAL_SEQ = "301212";
constexpr const char* SHINE_DAL_SEQ_LAG = "210303";

constexpr int8_t PROM_MAX_DIFF  = 8;
#endif


constexpr int8_t TERM_STEM_SIZE = 4;
constexpr int8_t TERM_LOOP_SIZE = 3;
constexpr int8_t TERM_SIZE      = 2 * TERM_STEM_SIZE + TERM_LOOP_SIZE;

constexpr int8_t SHINE_DAL_SIZE = 6;
constexpr int8_t SHINE_START_SPACER = 4;

constexpr int8_t CODON_SIZE  = 3;
constexpr int8_t CODON_START = 0b000;
constexpr int8_t CODON_STOP  = 0b001;
constexpr int8_t CODON_M0    = 0b100;
constexpr int8_t CODON_M1    = 0b101;
constexpr int8_t CODON_W0    = 0b010;
constexpr int8_t CODON_W1    = 0b011;
constexpr int8_t CODON_H0    = 0b110;
constexpr int8_t CODON_H1    = 0b111;


constexpr int8_t PROT_START_SIZE = SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE;
constexpr int32_t DO_TRANSLATION_LOOP = SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE;

#ifdef __REGUL
constexpr int8_t MAX_CODON   = 1 << CODON_SIZE;
constexpr int8_t QUADON_SIZE = 4;
constexpr int8_t MAX_QUADON  = 1 << QUADON_SIZE;
#endif

#ifdef __PROXY_POW_APPROX
constexpr int32_t LOOKUP_TABLE_SIZE = 10000000;
#endif

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 1.0;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double H_MIN = -1.0;
constexpr double H_MAX = 1.0;
constexpr double W_MIN = 0.0;
// W_MAX is defined through a parameter
constexpr double FUZZY_ROUNDING = 1000000000000.0;

constexpr int8_t SC_MATCH_BONUS   = 1;
constexpr int8_t SC_MISMATCH_COST = 2;

#ifdef BASE_4

constexpr int8_t NB_AMINO_ACIDS = 21; // not 20, because of STOP

typedef int8_t AminoAcid;
constexpr AminoAcid PHENYLALANINE = 0,
                    LEUCINE       = 1,
                    ISOLEUCINE    = 2,
                    METHIONINE    = 3,
                    VALINE        = 4,
                    SERINE        = 5,
                    PROLINE       = 6,
                    THREONINE     = 7,
                    ALANINE       = 8,
                    TYROSINE      = 9,
                    STOP          = 10,
                    HISTIDINE     = 11,
                    GLUTAMINE     = 12,
                    ASPARAGINE    = 13,
                    LYSINE        = 14,
                    ASPARTIC_ACID = 15,
                    GLUTAMIC_ACID = 16,
                    CYSTEINE      = 17,
                    TRYPTOPHAN    = 18,
                    ARGININE      = 19,
                    GLYCINE       = 20;

constexpr char* AMINO_ACIDS_NAMES[NB_AMINO_ACIDS] = {
    (char*) "phe",
    (char*) "leu",
    (char*) "iso",
    (char*) "met",
    (char*) "val",
    (char*) "ser",
    (char*) "pro",
    (char*) "thr",
    (char*) "ala",
    (char*) "tyr",
    (char*) "sto",
    (char*) "his",
    (char*) "glu",
    (char*) "asp",
    (char*) "lys",
    (char*) "asa",
    (char*) "gla",
    (char*) "cys",
    (char*) "try",
    (char*) "arg",
    (char*) "gly"
};

template<int8_t max_diff>
struct BasalLevel_ {
  double basalLevel[max_diff + 1];

  constexpr BasalLevel_() : basalLevel() {
    for(auto diff = 0; diff <= max_diff; diff++) {
      // formula :
      // if x < 2 : 1.04 ^ -[ (diff - 2)^2 ]
      //     else : 1.65 ^ (-diff+2)

      basalLevel[diff] = 1.0;
      /* GB simplification 
      if(diff - 2 >= 0) {
        int8_t exponent = -diff + 2; // always integer
        for(auto i = 0; i > exponent; i--) {
          basalLevel[diff] /= 1.65;   // because negative exponent
        }
      }
      else {
        int8_t exponent = -(diff - 2) * (diff - 2); // always integer
        for(auto i = 0; i > exponent; i--) {
          basalLevel[diff] /= 1.04;    // because negative exponent
        }
      }
      fin GB simplification*/
      basalLevel[diff] = 1 - (double)diff / (max_diff + 1); // GB new version
    }
  }

  inline const double& operator [](int i) const {
    return basalLevel[i];
  }
};


constexpr auto BASAL_LEVELS = BasalLevel_<PROM_MAX_DIFF>();
#endif
#endif // AEVOL_MACROS_H_
