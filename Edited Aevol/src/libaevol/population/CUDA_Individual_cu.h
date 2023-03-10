//
// Created by arrouan on 16/02/17.
//

#ifndef RAEVOL_CUDA_CUDA_INDIVIDUAL_CU_H
#define RAEVOL_CUDA_CUDA_INDIVIDUAL_CU_H

#include <cstdint>
#include <cstdio>

constexpr const char* PROM_SEQ_LEAD = "0101011001110010010110";
constexpr const char* PROM_SEQ_LAG  = "1010100110001101101001";


constexpr const char* SHINE_DAL_SEQ_LEAD = "011011000";
constexpr const char* SHINE_DAL_SEQ_LAG  = "100100111";

constexpr const char* PROTEIN_END_LEAD  = "001";
constexpr const char* PROTEIN_END_LAG   = "110";

constexpr int32_t PROMOTER_ARRAY_SIZE = 10000;
constexpr int32_t TERMINATOR_ARRAY_SIZE = 100000;

constexpr int32_t BUCKET_MAX_SIZE = 20;

constexpr int32_t RNA_LIST_INCR_SIZE = 200;
constexpr int32_t RNA_LIST_PROTEIN_INCR_SIZE = 250;

constexpr int32_t PROTEIN_LIST_INCR_SIZE = 400;
/**
 * Structure
 */
struct promoterStruct {
    int32_t pos = -1;
    int8_t error = -1;
    bool leading_or_lagging; // TRUE = leading / FALSE = lagging
};

typedef struct promoterStruct pStruct;


/**
 * Host structure
 *
 */
char **host_dna;
int8_t **host_dna_lead_promoter;
int8_t **host_dna_lag_promoter;
int8_t **host_dna_lead_term;
int8_t **host_dna_lag_term;
float **host_phenotype;
pStruct** host_dynPromoterList;
pStruct** host_dynTerminatorList;

/**
 * Structure to transfer from host to device memory
 */
// All DNAs
char** dna;

// All DNA size
size_t* dna_size;

// Current maximum DNA size
size_t* max_dna_size;
int host_max_dna_size;

// Hamming distance for promoter
int8_t** dna_lead_promoter;
int8_t** dna_lag_promoter;

// Terminator
int8_t** dna_lead_term;
int8_t** dna_lag_term;

// Number of promoters
int* nb_promoters;

// Maximum number of elements in RNA List
int* max_nb_elements_rna_list;

// Maximum number of elements in Protein List
int* max_nb_elements_protein_list;

// Environment (static first)
float* target;

/**
 * Internal structure (that stay on device)
 */

// Protein structures
struct pProtein {
    int32_t protein_start;
    int32_t protein_end;
    int32_t protein_length;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    float m;
    float w;
    float h;
    double e;
    bool is_functional;
};

typedef struct pProtein cProtein;

cProtein** protein_list;

int32_t* max_nb_protein;

int32_t* idx_protein;

int32_t* nb_protein;

// RNA structures
struct pRNA {
    int32_t begin;
    int32_t end;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    double e;
    uint32_t* start_prot;
    int32_t max_protein_elements;
    int32_t nb_protein;
    int32_t length;
    int32_t start_lenght;
};

typedef struct pRNA cRNA;

cRNA*** rna;

size_t* nb_rna;

int32_t* max_nb_rna;

int32_t* idx_rna;



// Fuzzy structures
float** phenotype;
float** delta;


pStruct** dynPromoterList;
pStruct** dynTerminatorList;

/**
 * Structure to transfer from device to host
 */
float* metaerror;
double* fitness;

/**
 * Kernel functions
 */

// Initialize (or reinitialize) internal struct before new generation
__global__ void init_array(int* nb_promoters);

// DNA -> RNA
__global__ void search_start_RNA(size_t* dna_size, char** dna,
                                 int8_t** dna_lead_promoter,
                                 int8_t** dna_lag_promoter, int* nb_promoters,
                                 pStruct** dynPromoterList, int block_size);
__global__ void search_stop_RNA(size_t* dna_size, char** dna,
                                int8_t** dna_lead_term, int8_t** dna_lag_term,
                                int block_size);
__global__ void search_start_stop_RNA_bucket(size_t* dna_size, char** dna, int8_t** dna_lead_promoter,
                                  int8_t** dna_lag_promoter, int* nb_promoters,
                                  pStruct** dynPromoterList, int8_t** dna_lead_term, int8_t** dna_lag_term,
                                             int bucket_size, int block_size);
__global__ void init_RNA_struct(int pop_size, cRNA*** rna, int* nb_promoters,
                                int32_t* max_nb_rna,int32_t* idx_rna, int32_t* max_nb_elements_rna_list);
__global__ void internal_init_RNA_struct(cRNA*** rna, int32_t* max_nb_elements_rna_list, int indiv_id);
__global__ void compute_RNA(int pop_size, int8_t** dna_lead_promoter,
                            int8_t** dna_lag_promoter, int8_t** dna_lead_term,
                            int8_t** dna_lag_term, char** dna, size_t* dna_size,
                            cRNA*** rna, int32_t* idx_rna, int* nb_promoters, pStruct** dynPromoterList,
                            int thread_size, int thread_dim);
__global__ void max_rna(int32_t* idx_rna, int32_t* max_nb_rna);

// RNA -> Protein
__global__ void compute_start_protein(int32_t* idx_rna, cRNA*** rna,
                                      char** dna,size_t* dna_size, int32_t* nb_protein, int threads_size, int thread_dim);
__global__ void init_protein_struct(int pop_size, int32_t* nb_protein,
                                    cProtein** protein_list, int32_t* idx_protein,
                                    cRNA*** rna, int32_t* idx_rna, int32_t* max_nb_protein, int* max_nb_elements_protein_list,
                                    int threads_size, int thread_dim);
__global__ void compute_protein(cRNA*** rna, cProtein** protein_list, int32_t* idx_protein,
                                size_t* dna_size,char** dna, int32_t* idx_rna,int threads_size, int thread_dim, int block_size);
__global__ void max_protein(int32_t* max_nb_protein, int32_t* idx_protein);
__global__ void translate_protein(float w_max, int32_t* idx_protein,
                                  cProtein** protein_list,
                                  char** dna, size_t* dna_size,int threads_size, int thread_dim);

// Protein -> Phenotype
__global__ void compute_phenotype(int32_t* idx_protein, cProtein** protein_list,
                                  float** phenotype,int threads_size, int thread_dim);

// Phenotype -> Metaerror -> Fitness
__global__ void compute_metaerror_fitness(double selection_pressure,float** phenotype,
                                          float* target,
                                          float* metaerror, double* fitness);

__global__ void free_list(cProtein** protein_list,
                          cRNA*** rna, int32_t* idx_protein,int32_t* idx_rna);

// Debug kernels
__global__ void debug_dna(size_t* dna_size, char** dna);

__global__ void debug_promoter_start(size_t* dna_size, pStruct** dynPromoterList,
                                     int* nb_promoters,int indiv_id);
__global__ void debug_promoter_stop(size_t* dna_size,
                                    int8_t** dna_lead_term,
                                    int8_t** dna_lag_term, int* nb_promoters, int indiv_id);
__global__ void debug_rna(size_t* dna_size,
                          int8_t** dna_lead_term,
                          int8_t** dna_lag_term,
                          cRNA*** rna,int32_t* idx_rna,
                          int indiv_id);
__global__ void debug_protein(int32_t* idx_protein,
                              cProtein** protein_list, char** dna,
                              int indiv_id);

__global__ void display_size_dna(size_t* dna_size);
__global__ void debug_phenotype(float** phenotype,float* target, float* metaerror, double* fitness,
                                int indiv_id);
__global__ void debug_fitness(float** phenotype,float* target,
                              float* metaerror, double* fitness,
                              int indiv_id);
#endif //RAEVOL_CUDA_CUDA_INDIVIDUAL_CU_H
