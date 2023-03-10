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
//                              Includes
// =================================================================
#include <inttypes.h>
#include <cstring>

#include <list>
#include <memory>
#include <iostream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "Individual.h"
#include "macros.h"
#include "GeneticUnit.h"
#include "Rna.h"
#include "Protein.h"
#include "MutationParams.h"
#include "../libaevol/macros.h"
#include "ExpSetup.h"

using namespace aevol;

using std::list;
using std::istream;
using std::ostream;
using std::string;
using std::cout;
using std::endl;




class TranscriptionTranslationTest : public testing::Test
{
 protected:
  Individual* indiv;

  void check_genome(const string& dir, int generation);
};

using protein_line = struct {
  int32_t id;
  string c_or_p;
  string strand;
  int32_t pos;
  int32_t len;
  int32_t lpos;
  string sequence;
  double m;
  double w;
  double h;
  double c;
  int f;
  int32_t prom_pos;
  int32_t rna_len;
  double basal_level;
};

istream& operator>> (istream& is, protein_line& prot_line) {
  return is >> prot_line.id >> prot_line.c_or_p >> prot_line.strand >>
      prot_line.pos >> prot_line.len >> prot_line.lpos >>
      prot_line.sequence >> prot_line.m >> prot_line.w >> prot_line.h >>
      prot_line.c >> prot_line.f >> prot_line.prom_pos >> prot_line.rna_len >>
      prot_line.basal_level;
}
ostream& operator<< (ostream& os, protein_line& prot_line) {
  return os << prot_line.id << " " << prot_line.c_or_p << " " <<
      prot_line.strand << " " << prot_line.pos << " " <<
      prot_line.len << " " << prot_line.lpos << " " <<
      prot_line.sequence << " " << prot_line.m << " " <<
      prot_line.w << " " << prot_line.h << " " <<
      prot_line.c << " " << prot_line.f << " " <<
      prot_line.prom_pos << " " << prot_line.rna_len << " " <<
      prot_line.basal_level;
}



// Temporary helper derived from
// extract.cpp:analyse_gu(GeneticUnit*, int32_t, FILE*, Environment*)
list<protein_line> analyse_gu(const GeneticUnit& gen_unit, int32_t gen_unit_number)
{
  list<protein_line> proteins;

  Promoters2Strands llrnas = gen_unit.rna_list();
  for(auto lrnas : llrnas) {
    for (auto rna : lrnas) {
      for (auto prot : rna.transcribed_proteins()) {
        double mean = prot->mean();
        int nfeat = -1;
        protein_line prot_line;

        prot_line.id = gen_unit.indiv()->id();
        prot_line.c_or_p = gen_unit_number != 0 ? "PLASMID" : "CHROM";
        prot_line.strand = prot->strand() == LEADING ? "LEADING" : "LAGGING";
        prot_line.pos = prot->first_translated_pos();
        prot_line.len = prot->length();
        prot_line.lpos = prot->last_translated_pos();
        prot_line.sequence = prot->AA_sequence('_');
        prot_line.m = mean;
        prot_line.w = prot->width();
        prot_line.h = prot->height();
        prot_line.c = prot->concentration();
        prot_line.f = nfeat;
        prot_line.prom_pos = rna.promoter_pos();
        prot_line.rna_len = rna.transcript_length();
        prot_line.basal_level = rna.basal_level();

        proteins.push_back(prot_line);
      }
    }
  }

  return proteins;
}

void expect_equal(const list<protein_line> expected_proteins,
                  const list<protein_line> actual_proteins) {
  auto exp_prot = expected_proteins.begin();
  for (auto act_prot = actual_proteins.begin() ;
      exp_prot != expected_proteins.end() &&
      act_prot != actual_proteins.end() ;
      exp_prot++, act_prot++) {
    EXPECT_EQ(exp_prot->c_or_p, act_prot->c_or_p);
    EXPECT_EQ(exp_prot->strand, act_prot->strand);
    EXPECT_EQ(exp_prot->pos, act_prot->pos);
    EXPECT_EQ(exp_prot->len, act_prot->len);
    EXPECT_EQ(exp_prot->lpos, act_prot->lpos);
    EXPECT_EQ(exp_prot->sequence, act_prot->sequence);
    EXPECT_NEAR(exp_prot->m, act_prot->m, 0.000001);
    EXPECT_NEAR(exp_prot->w, act_prot->w, 0.000001);
    EXPECT_NEAR(exp_prot->h, act_prot->h, 0.000001);
    EXPECT_NEAR(exp_prot->c, act_prot->c, 0.000001);
    EXPECT_EQ(exp_prot->f, act_prot->f);
    EXPECT_EQ(exp_prot->prom_pos, act_prot->prom_pos);
    EXPECT_EQ(exp_prot->rna_len, act_prot->rna_len);
    EXPECT_FLOAT_EQ(exp_prot->basal_level, act_prot->basal_level);
  }
}

void TranscriptionTranslationTest::check_genome(const string& dir, int generation) {

  ExpSetup* exp_s = new ExpSetup(nullptr);
  for (int i = 0; i <= 1; i++) {
    exp_s->set_fuzzy_flavor(i);

    if (FuzzyFactory::fuzzyFactory == NULL)
      FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s);
    else {
      delete FuzzyFactory::fuzzyFactory;
      FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s);
    }
    // Read genome from input file
    std::filebuf fb;
    string genome;
    std::ostringstream ostr;
    ostr << std::setfill('0') << std::setw(6) << generation;
    string gener = ostr.str();

    fb.open(dir + "/sequence_" + gener, std::ios::in);
    ASSERT_TRUE(fb.is_open());

    istream genome_file(&fb);
    std::getline(genome_file, genome);
    fb.close();

    // Construct individual with the genome we've read
    MutationParams params_mut;
    indiv = new Individual(nullptr, nullptr, nullptr,
                           std::make_shared<MutationParams>(params_mut), 0.1,
                           10, 1000, false, 1, "anon-strain-1", 0);
    char* raw_genome = new char[genome.size() + 1];
    strcpy(raw_genome, genome.c_str());
    indiv->add_GU(raw_genome, genome.size());

    // Do transcription and translation
    indiv->do_transcription();
    indiv->do_translation();

    auto actual_proteins = analyse_gu(indiv->genetic_unit(0), 0);


    // Read the expected results
    list<protein_line> expected_proteins;
    fb.open(dir + "/proteins_" + gener, std::ios::in);
    ASSERT_TRUE(fb.is_open());

    string buffer;
    istream proteins_file(&fb);
    // Flush header
    do std::getline(proteins_file, buffer);
    while (buffer[0] != 'i');

    protein_line prot_line;
    while (proteins_file) {
      proteins_file >> prot_line;
      expected_proteins.push_back(prot_line);
    }
    // (very dirty) Get rid of last line (added twice)
    expected_proteins.pop_back();


    // cout << "*************** EXPECTED ********************" << endl;
    // for (auto prot_line : expected_proteins) {
    //   cout << prot_line << endl;
    // }
    // cout << "**************** ACTUAL *********************" << endl;
    // for (auto prot_line : actual_proteins) {
    //   cout << prot_line << endl;
    // }

    expect_equal(expected_proteins, actual_proteins);
  }
}

TEST_F(TranscriptionTranslationTest, TestIndivVirus6) {
  for (int i = 10000 ; i <= 200000 ; i += 10000)
    check_genome("TranscriptionTranslationTest_files/virus6", i);
}

TEST_F(TranscriptionTranslationTest, TestIndivVirus7) {
  for (int i = 10000 ; i <= 200000 ; i += 10000)
    check_genome("TranscriptionTranslationTest_files/virus7", i);
}

TEST_F(TranscriptionTranslationTest, TestIndivVirus8) {
  for (int i = 10000 ; i <= 200000 ; i += 10000)
    check_genome("TranscriptionTranslationTest_files/virus8", i);
}

TEST_F(TranscriptionTranslationTest, TestIndivVirus9) {
  for (int i = 10000 ; i <= 200000 ; i += 10000)
    check_genome("TranscriptionTranslationTest_files/virus9", i);
}

TEST_F(TranscriptionTranslationTest, TestIndivVirus10) {
  for (int i = 10000 ; i <= 200000 ; i += 10000)
    check_genome("TranscriptionTranslationTest_files/virus10", i);
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
