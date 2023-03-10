#ifndef __CUDA_Individual
#define __CUDA_Individual

#include <cstdint>

#include "ExpManager.h"

using namespace aevol;

void cuda_init();

/**
 * Data transfer functions
 */
void transfer_in(ExpManager* exp_m, bool init_all_struct);
void transfer_out(ExpManager* exp_m, bool delete_all_struct);

/**
 * Run all kernel for one generation -- all individuals
 */
void run_a_step(int nb_indiv,float w_max, double selection_pressure, bool first_gen);

/**
 * Debug functions
 */
void print_debug_promoters_start(ExpManager* exp_m);
void print_debug_rna(ExpManager* exp_m);
void print_debug_protein(ExpManager* exp_m);
void print_debug_phenotype(ExpManager* exp_m);
void print_debug_fitness(ExpManager* exp_m);


void print_debug_promoters_start(ExpManager* exp_m, int i);
void print_debug_rna(ExpManager* exp_m, int i);
void print_debug_protein(ExpManager* exp_m, int i);
void print_debug_phenotype(ExpManager* exp_m, int i);
#endif
