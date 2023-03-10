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
#include <stdio.h>
#include <sys/stat.h>
#include <err.h>
#include <errno.h>


// =================================================================
//                            Project Files
// =================================================================
#include "Dump.h"
#include "ExpManager.h"
#include "Individual.h"
#include "GeneticUnit.h"
#ifdef __REGUL
  #include "raevol/Protein_R.h"
#endif

namespace aevol {





//##############################################################################
//                                                                             #
//                                Class Dump                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Dump::Dump(ExpManager * exp_m)
{
  exp_m_ = exp_m;
  int status;
  status = mkdir("stats/dump/", 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, "stats/dump/");
  }
}

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================

const char* DUMP_FORMAT = "\t%d\t%d\t%f\n";

void Dump::write_current_generation_dump()
{
  //  printf("Begin dump\n");
  write_fitness_total();
  write_secretion_present();
  write_fitness_metabolic();
  write_secreted_amount();
  write_individual_probes();
  //  printf("End dump\n");
}

void Dump::write_fitness_total()
{
  sprintf(filename_buffer,
      "stats/dump/fitness_total_" TIMESTEP_FORMAT ".out",
      AeTime::time());
  current_file = fopen(filename_buffer, "w+");
  double** map = exp_m_->world()->total_fitness_grid();
  fprintf(current_file, "#\tX\tY\tfitness_total(X, Y)\n");

  for(int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    for(int16_t y = 0 ; y < exp_m_->grid_height() ; y++)
    {
      fprintf(current_file, DUMP_FORMAT, x, y, map [x][y]);
    }
    fprintf(current_file, "\n");
  }
  fflush(current_file);
  fclose(current_file);

  // Has been allocated in World::total_fitness_grid()
  for (int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    delete [] map[x];
  }
  delete [] map;
}

void Dump::write_secreted_amount()
{
  sprintf(filename_buffer,
      "stats/dump/secreted_amount_" TIMESTEP_FORMAT ".out",
      AeTime::time()) ;
  current_file = fopen(filename_buffer, "w+");

  double** map = exp_m_->world()->secreted_amount_grid();
  fprintf(current_file, "#\tX\tY\tsecreted_amount(X, Y)\n");
  for(int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    for(int16_t y = 0 ; y < exp_m_->grid_height() ; y++)
    {
      fprintf(current_file, DUMP_FORMAT, x, y, map [x][y]);
    }
    fprintf(current_file, "\n");
  }
  fflush(current_file);
  fclose(current_file);
  for (int16_t x = 0; x < exp_m_->grid_width() ; x++)
  {
    delete [] map[x];
  }
  delete [] map;
}

void Dump::write_fitness_metabolic()
{
  sprintf(filename_buffer,
      "stats/dump/fitness_metabolic_" TIMESTEP_FORMAT ".out",
      AeTime::time());
  current_file = fopen(filename_buffer, "w+");

  double** map = exp_m_->world()->metabolic_fitness_grid();
  fprintf(current_file, "#\tX\tY\tfitness_metabolic(X, Y)\n");
  for(int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    for(int16_t y = 0 ; y < exp_m_->grid_height() ; y++)
    {
      fprintf(current_file, DUMP_FORMAT, x, y, map [x][y]);
    }
    fprintf(current_file, "\n");
  }
  fflush(current_file);
  fclose(current_file);
  for (int16_t x = 0; x < exp_m_->grid_width() ; x++)
  {
    delete [] map[x];
  }
  delete [] map;
}

void Dump::write_secretion_present()
{
  sprintf(filename_buffer,
      "stats/dump/secretion_present_" TIMESTEP_FORMAT ".out",
      AeTime::time());
  current_file = fopen(filename_buffer, "w+");

  double** map = exp_m_->world()->secretion_present_grid();
  fprintf(current_file, "#\tX\tY\tsecretion_present(X, Y)\n");
  for(int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    for(int16_t y = 0 ; y < exp_m_->grid_height() ; y++)
      fprintf(current_file, DUMP_FORMAT, x, y, map [x][y]);
    fprintf(current_file, "\n");
  }
  fflush(current_file);
  fclose(current_file);
  for (int16_t x = 0; x < exp_m_->grid_width() ; x++)
  {
    delete [] map[x];
  }
  delete [] map;
}

/*!
  \brief Write the probes (5 int and 5 double) of each individual at a given generation
*/
void Dump::write_individual_probes()
{
  sprintf(filename_buffer,
      "stats/dump/individual_probes_" TIMESTEP_FORMAT ".out",
      AeTime::time());
  current_file = fopen(filename_buffer, "w");

  fprintf(current_file, "Id\tInt_Probe_1\tInt_Probe_2\tInt_Probe_3\tInt_Probe_4\tInt_Probe_5\tDouble_Probe_1\tDouble_Probe_2\tDouble_Probe_3\tDouble_Probe_4\tDouble_Probe_5\n");

  for(int16_t x = 0 ; x < exp_m_->grid_width() ; x++)
  {
    for(int16_t y = 0 ; y < exp_m_->grid_height() ; y++)
    {
      fprintf(current_file, "%llu",
          exp_m_->world()->indiv_at(x,y)->id());
      int32_t* int_probes =
          exp_m_->world()->indiv_at(x,y)->int_probes();
      double* double_probes =
          exp_m_->world()->indiv_at(x,y)->double_probes();
      for(int16_t i=0; i<5; i++)
        fprintf(current_file, "\t%" PRId32, int_probes[i]);
      for(int16_t i=0; i<5; i++)
        fprintf(current_file, "\t%f", double_probes[i]);
      fprintf(current_file, "\n");
    }
  }
  fflush(current_file);
  fclose(current_file);
}

// =================================================================
//                           Protected Methods
// =================================================================

} // namespace aevol
