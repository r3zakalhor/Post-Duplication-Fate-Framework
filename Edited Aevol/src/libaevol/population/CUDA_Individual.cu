#include "CUDA_Individual.h"
#include "CUDA_Individual_cu.h"
#include <cstdint>
#include <stdio.h>
#include <unistd.h>
#include<cuda.h>
#include<cuda_profiler_api.h>

#include "ExpManager.h"
#include "HybridFuzzy.h"

#define DEBUG 1
// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n",
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

#define HANDLE_ERROR( err ) ( HandleError( err, __FILE__, __LINE__ ) )

static void HandleError( cudaError_t err, const char *file, int line )
{
  if (err != cudaSuccess)
  {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
            file, line );
    exit( EXIT_FAILURE );
  }
}

void cuda_init() {
  size_t limit = 9000000000;
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, limit);
  cudaDeviceSetLimit(cudaLimitStackSize, 1024*1024);
  cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 64*1024*1024);
}

void transfer_in(ExpManager* exp_m, bool init_all_struct) {
  if (init_all_struct) {
    cudaDeviceReset();
    cudaDeviceSynchronize();
  }

  size_t uCurAvailMemoryInBytes;
  size_t uTotalMemoryInBytes;
  CUresult result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }

  host_dna = (char**)malloc(exp_m->nb_indivs()*sizeof(char*));
  checkCuda(cudaMalloc((void***)&dna,exp_m->nb_indivs()*sizeof(char*)));

  //host_dna_lead_promoter = (int8_t**)malloc(exp_m->nb_indivs()*sizeof(int8_t*));
  //checkCuda(cudaMalloc((void***)&dna_lead_promoter,exp_m->nb_indivs()*sizeof(int8_t*)));
  //host_dna_lag_promoter = (int8_t**)malloc(exp_m->nb_indivs()*sizeof(int8_t*));
  //checkCuda(cudaMalloc((void***)&dna_lag_promoter,exp_m->nb_indivs()*sizeof(int8_t*)));

  host_dna_lead_term = (int8_t**)malloc(exp_m->nb_indivs()*sizeof(int8_t*));
  checkCuda(cudaMalloc((void***)&dna_lead_term,exp_m->nb_indivs()*sizeof(int8_t*)));
  host_dna_lag_term = (int8_t**)malloc(exp_m->nb_indivs()*sizeof(int8_t*));
  checkCuda(cudaMalloc((void***)&dna_lag_term,exp_m->nb_indivs()*sizeof(int8_t*)));

  host_phenotype = (float**)malloc(exp_m->nb_indivs()*sizeof(float*));
  checkCuda(cudaMalloc((void***)&phenotype,exp_m->nb_indivs()*sizeof(float*)));

  checkCuda(cudaMalloc((void**)&dna_size,
             exp_m->nb_indivs() * sizeof(size_t)));
  size_t* host_dna_size = (size_t*)malloc(exp_m->nb_indivs()*sizeof(size_t));


  checkCuda(cudaMalloc((void**)&nb_promoters,
             exp_m->nb_indivs() * sizeof(int)));
  checkCuda(cudaMemset(nb_promoters, 0, exp_m->nb_indivs() * sizeof(int)));

  if (init_all_struct) {
    printf("Init struct");
    checkCuda(cudaMalloc((void**) &max_nb_elements_rna_list,
                         exp_m->nb_indivs() * sizeof(int)));
    checkCuda(cudaMemset(max_nb_elements_rna_list, 0, exp_m->nb_indivs() * sizeof(int)));

    checkCuda(cudaMalloc((void**) &max_nb_elements_protein_list,
                         exp_m->nb_indivs() * sizeof(int)));
    checkCuda(cudaMemset(max_nb_elements_protein_list, 0, exp_m->nb_indivs() * sizeof(int)));
  }

  checkCuda(cudaMalloc((void**)&metaerror,
             exp_m->nb_indivs() * sizeof(float)));


  checkCuda(cudaMalloc((void**)&fitness,
             exp_m->nb_indivs() * sizeof(double)));


  checkCuda(cudaMalloc((void**)&nb_protein,
             exp_m->nb_indivs() * sizeof(int32_t)));
  checkCuda(cudaMemset(nb_protein, 0, exp_m->nb_indivs() * sizeof(int32_t)));


  checkCuda(cudaMalloc((void**)&max_nb_rna,
                       sizeof(int32_t)));
  checkCuda(cudaMemset(max_nb_rna, 0, sizeof(int32_t)));


  checkCuda(cudaMalloc((void**)&max_nb_protein,
                       sizeof(int32_t)));
  checkCuda(cudaMemset(max_nb_protein, 0, sizeof(int32_t)));

  checkCuda(cudaMalloc((void**)&idx_rna,
             exp_m->nb_indivs() * sizeof(int32_t)));
  checkCuda(cudaMemset(idx_rna, 0, exp_m->nb_indivs() * sizeof(int32_t)));

  if (init_all_struct) checkCuda(cudaMalloc((void***)&protein_list,exp_m->nb_indivs()*sizeof(cProtein*)));

  checkCuda(cudaMalloc((void**)&idx_protein,
                       exp_m->nb_indivs() * sizeof(int32_t)));
  checkCuda(cudaMemset(idx_protein, 0, exp_m->nb_indivs() * sizeof(int32_t)));

  if (init_all_struct) checkCuda(cudaMalloc((void***)&rna,exp_m->nb_indivs()*sizeof(cRNA**)));

  checkCuda(cudaMalloc((void***)&dynPromoterList,
                       exp_m->nb_indivs()*sizeof(pStruct*)));
  host_dynPromoterList = (pStruct**)malloc(exp_m->nb_indivs()*sizeof(pStruct*));

  checkCuda(cudaMalloc((void***)&dynTerminatorList,
                       exp_m->nb_indivs()*sizeof(pStruct*)));
  host_dynTerminatorList = (pStruct**)malloc(exp_m->nb_indivs()*sizeof(pStruct*));

  result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }


  int x,y,max_dna=-1;

  for (int i = 0; i < exp_m->nb_indivs(); i++) {

    x = i / exp_m->grid_height();
    y = i % exp_m->grid_height();

    checkCuda(cudaMalloc((void**) &host_dna[i],
               exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(char)));

    checkCuda(cudaMemcpy(host_dna[i], exp_m->world()->indiv_at(x, y)->genetic_unit(0).dna()->data(),
               exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(char),
               cudaMemcpyHostToDevice));

    host_dna_size[i] = exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length();
    max_dna = max_dna < exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() ?
              exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() : max_dna;

    //checkCuda(cudaMalloc((void**) &host_dna_lead_promoter[i],
    //           exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(int8_t)));

    //checkCuda(cudaMalloc((void**) &host_dna_lag_promoter[i],
    //           exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(int8_t)));

    checkCuda(cudaMalloc((void**) &host_dna_lead_term[i],
               exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(int8_t)));

    checkCuda(cudaMalloc((void**) &host_dna_lag_term[i],
               exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(int8_t)));


    checkCuda(cudaMalloc((void**) &host_phenotype[i], 300 * sizeof(float)));
    checkCuda(cudaMemset(host_phenotype[i], 0.0, 300 * sizeof(float)));

    checkCuda(cudaMalloc((void**) &host_dynPromoterList[i], PROMOTER_ARRAY_SIZE * sizeof(pStruct)));
    checkCuda(cudaMalloc((void**) &host_dynTerminatorList[i], exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length() * sizeof(pStruct)));
  }

  host_max_dna_size= max_dna;

  checkCuda(cudaMemcpy(dna,host_dna,exp_m->nb_indivs()*sizeof(char*),cudaMemcpyHostToDevice));

  //checkCuda(cudaMemcpy(dna_lag_promoter,host_dna_lag_promoter,exp_m->nb_indivs()*sizeof(int8_t*),cudaMemcpyHostToDevice));
  //checkCuda(cudaMemcpy(dna_lead_promoter,host_dna_lead_promoter,exp_m->nb_indivs()*sizeof(int8_t*),cudaMemcpyHostToDevice));


  checkCuda(cudaMemcpy(phenotype,host_phenotype,exp_m->nb_indivs()*sizeof(float*),cudaMemcpyHostToDevice));

  checkCuda(cudaMemcpy(dna_lag_term,host_dna_lag_term,exp_m->nb_indivs()*sizeof(int8_t*),cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(dna_lead_term,host_dna_lead_term,exp_m->nb_indivs()*sizeof(int8_t*),cudaMemcpyHostToDevice));

  checkCuda(cudaMemcpy(dynPromoterList,host_dynPromoterList,exp_m->nb_indivs()*sizeof(pStruct*),cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(dynTerminatorList,host_dynTerminatorList,exp_m->nb_indivs()*sizeof(pStruct*),cudaMemcpyHostToDevice));

  checkCuda(cudaMemcpy(dna_size,
             host_dna_size, exp_m->nb_indivs() * sizeof(size_t), cudaMemcpyHostToDevice));

  /*free(host_dna_size);
  free(host_phenotype);
  free(host_dna_lag_term);
  free(host_dna_lead_term);
  free(host_dna_lag_promoter);
  free(host_dna_lead_promoter);
  free(host_dna);*/

  checkCuda(cudaMalloc((void**) &target,
             300 * sizeof(float)));

  float target_host[300];
  for (int i = 0; i < 300; i++) {
    target_host[i] = ((HybridFuzzy*)exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[i];
  }

  checkCuda(cudaMemcpy(target,
                       target_host,
             300 * sizeof(float), cudaMemcpyHostToDevice));

}

void print_debug_promoters_start(ExpManager* exp_m) {
  print_debug_promoters_start(exp_m,737);
}

void print_debug_promoters_start(ExpManager* exp_m, int i) {
  int x,y;

  //for (int i = 0; i < exp_m->nb_indivs(); i++) {
    x = i / exp_m->grid_height();
    y = i % exp_m->grid_height();

    for (auto& strand_id: {LEADING, LAGGING}) {
      if (strand_id == LEADING) printf("Individual %d (CPU) Promoters : LEADING ",i);
      else printf("Individual %d (CPU) Promoters : LAGGING ",i);
      auto& strand = exp_m->world()->indiv_at(x, y)->genetic_unit(0).rna_list()[strand_id];
      for (auto rna = strand.begin(); rna != strand.end(); ++rna) {
        printf("%d ",rna->promoter_pos());
      }
      printf("\n");
    }

    debug_promoter_start<<<1,1>>>(dna_size,dynPromoterList,nb_promoters,
        i);
printf("Terminator !!! \n");
  debug_promoter_stop<<<1,1>>>(dna_size,dna_lead_term,dna_lag_term,nb_promoters,
      i);
  cudaDeviceSynchronize();
  //}
  printf("Terminator END !!! \n");
}

void print_debug_rna(ExpManager* exp_m) {
  print_debug_rna(exp_m,737);
}

void print_debug_rna(ExpManager* exp_m, int i) {
  int x,y;
  //for (int i = 0; i < exp_m->nb_indivs(); i++) {
    x = i / exp_m->grid_height();
    y = i % exp_m->grid_height();

  int prot_idx = 0;

int rna_idx = 0;
    for (auto& strand_id: {LEADING, LAGGING}) {
      printf("%d -- Individual %d (%d %d CPU) \n",AeTime::time(),i,x,y);
      auto& strand = exp_m->world()->indiv_at(x, y)->genetic_unit(0).rna_list()[strand_id];
      for (auto rna = strand.begin(); rna != strand.end(); ++rna) {
        printf("RNA %d : ",rna_idx);

        if (strand_id == LEADING) printf("LEADING ");
        else printf("LAGGING ");

        printf("%d (%d) %d %f (size %d)\n",rna->promoter_pos(),
               rna->first_transcribed_pos(),
            rna->last_transcribed_pos(),
            rna->basal_level(),rna->transcript_length());

/*
        for (auto prot : rna->transcribed_proteins()) {
          printf("Protein CPU (%d) %d : %d %d (%d %d) %f %f %f\n",prot_idx,
                 rna_idx,
                 prot->shine_dal_pos(),
                 prot->last_STOP_base_pos(),
                 prot->first_translated_pos(),
                 prot->last_translated_pos(),
                 prot->mean(),
                 prot->height(),
                 prot->width());
          prot_idx++;
        }*/
      rna_idx++;
      }

    }

    debug_rna<<<1,1>>>(dna_size,dna_lead_promoter,dna_lag_promoter,rna,idx_rna,
        i);
  //}

}

void print_debug_protein(ExpManager* exp_m) {
  print_debug_protein(exp_m,737);
}

void print_debug_protein(ExpManager* exp_m, int i) {
  int x,y;

  //for (int i = 0; i < exp_m->nb_indivs(); i++) {
  x = i / exp_m->grid_height();
  y = i % exp_m->grid_height();

  int prot_idx = 0;
  printf("Protein list size on CPU %d (%d %d): %d (%lu)\n",i,x,y,exp_m->world()->indiv_at(x, y)->protein_list().size(),
         exp_m->world()->indiv_at(x, y)->genetic_unit(0).seq_length());
  for (auto prot : exp_m->world()->indiv_at(x, y)->protein_list()) {
    printf("Protein CPU (%d) %d : %d %d (%d %d) %lf %lf %lf isfunctional %d\n",prot_idx,
           prot->rna_list().size(),
           prot->shine_dal_pos(),
           prot->last_STOP_base_pos(),
           prot->first_translated_pos(),
           prot->last_translated_pos(),
           prot->mean(),
           prot->height(),
           prot->width(),
           prot->is_functional());
   /* int cod_idx = 0;
    for (auto cod : prot->AA_list()) {
      printf("COD[%d] = %d\n",cod_idx,cod->value());
      cod_idx++;
    }*/
  }
  debug_protein<<<1,1>>>(idx_protein,protein_list,dna,
      i);
  //}

}

void print_debug_phenotype(ExpManager* exp_m) {
  print_debug_phenotype(exp_m,737);
}

void print_debug_phenotype(ExpManager* exp_m, int i) {
  int x,y;
  //for (int i = 0; i < exp_m->nb_indivs(); i++) {
  x = i / exp_m->grid_height();
  y = i % exp_m->grid_height();

  exp_m->world()->indiv_at(x, y)->phenotype()->print();

  debug_phenotype<<<1,1>>>(phenotype,target,metaerror,fitness,
      i);
}

void print_debug_fitness(ExpManager* exp_m) {
  int x,y;
  int i = 15;
  //for (int i = 0; i < exp_m->nb_indivs(); i++) {
  x = i / exp_m->grid_height();
  y = i % exp_m->grid_height();

  //printf("%d %d\n",x,y);

  /*debug_fitness<<<1,1>>>(phenotype, target,
      metaerror, fitness,
      i);*/
}

void transfer_out(ExpManager* exp_m, bool delete_all_struct) {
  printf("Starting transfert\n");

  float* host_metaerror = (float*)malloc(exp_m->nb_indivs()*sizeof(float));
  double* host_fitness = (double*)malloc(exp_m->nb_indivs()*sizeof(double));
  /*size_t uCurAvailMemoryInBytes;
  size_t uTotalMemoryInBytes;
  CUresult result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }*/
  cudaDeviceSynchronize();
  //printf("Transfering Metaerror\n");
  checkCuda(cudaMemcpy(host_metaerror,
             metaerror, exp_m->nb_indivs() * sizeof(float), cudaMemcpyDeviceToHost));
  //printf("Transfering Fitness\n");
  checkCuda(cudaMemcpy(host_fitness,
             fitness, exp_m->nb_indivs() * sizeof(double), cudaMemcpyDeviceToHost));
  printf("Transfer END\n");

  bool error_detected=false;

  int x,y;
  for (int i = 0; i < exp_m->nb_indivs(); i++) {
    x = i / exp_m->grid_height();
    y = i % exp_m->grid_height();

    double fit_1 = exp_m->world()->indiv_at(x, y)->dist_to_target_by_feature(METABOLISM);
    double fit_2 = host_metaerror[i];
    float i_fit_1 = roundf(fit_1*100);
    float i_fit_2 = roundf(fit_2*100);


    if (i_fit_1 != i_fit_2) {
      printf(
          "ERROR -- Individual %d : Metaerror (CPU/GPU) : %e -- %e || Fitness (CPU/GPU) : %e -- %e\n",
          i,
          exp_m->world()->indiv_at(x, y)->dist_to_target_by_feature(METABOLISM),
          host_metaerror[i],
          exp_m->world()->indiv_at(x, y)->fitness(), host_fitness[i]);


      /*if (i == 0) {
        print_debug_promoters_start(exp_m,i);
        print_debug_rna(exp_m,i);
        print_debug_protein(exp_m,i);
      }*/

      //print_debug_promoters_start(exp_m,i);
      //print_debug_rna(exp_m,i);
      //print_debug_protein(exp_m,i);
      //print_debug_phenotype(exp_m,i);

      //char c=getchar();
      //printf("Read %c\n",c);
      //if (c=='q') { error_detected = true;}
    }
  }

  //free_list<<<1024,1>>>(protein_list,rna,idx_protein,idx_rna);

 for (int i = 0; i < exp_m->nb_indivs(); i++) {
    HANDLE_ERROR(cudaFree(host_dna[i]));
   //HANDLE_ERROR(cudaFree(host_dna_lead_promoter[i]));
   //HANDLE_ERROR(cudaFree(host_dna_lag_promoter[i]));
   HANDLE_ERROR(cudaFree(host_dna_lead_term[i]));
   HANDLE_ERROR(cudaFree(host_dna_lag_term[i]));
   HANDLE_ERROR(cudaFree(host_phenotype[i]));
   HANDLE_ERROR(cudaFree(host_dynPromoterList[i]));
   HANDLE_ERROR(cudaFree(host_dynTerminatorList[i]));
  }

  HANDLE_ERROR(cudaFree(nb_promoters));
  if (delete_all_struct) HANDLE_ERROR(cudaFree(max_nb_elements_rna_list));

  HANDLE_ERROR(cudaFree(nb_protein));
  HANDLE_ERROR(cudaFree(metaerror));
  HANDLE_ERROR(cudaFree(fitness));
  HANDLE_ERROR(cudaFree(idx_protein));
  HANDLE_ERROR(cudaFree(idx_rna));
  if (delete_all_struct) HANDLE_ERROR(cudaFree(protein_list));
  if (delete_all_struct) HANDLE_ERROR(cudaFree(rna));
  HANDLE_ERROR(cudaFree(target));

  HANDLE_ERROR(cudaFree(dna));
  //HANDLE_ERROR(cudaFree(dna_lead_promoter));
  //HANDLE_ERROR(cudaFree(dna_lag_promoter));
  HANDLE_ERROR(cudaFree(dna_lead_term));
  HANDLE_ERROR(cudaFree(dna_lag_term));
  HANDLE_ERROR(cudaFree(phenotype));
  HANDLE_ERROR(cudaFree(dynPromoterList));
  HANDLE_ERROR(cudaFree(dynTerminatorList));

  /*if (error_detected)
    exit(-42);*/
}

void run_a_step(int nb_indiv,float w_max, double selection_pressure, bool first_gen) {
  cuda_init();
  nb_indiv = 1024;


  /*limit=0;
  cudaDeviceGetLimit(&limit, cudaLimitStackSize);
  printf("cudaLimitStackSize: %u\n", (unsigned)limit);
  cudaDeviceGetLimit(&limit, cudaLimitPrintfFifoSize);
  printf("cudaLimitPrintfFifoSize: %u\n", (unsigned)limit);
  cudaDeviceGetLimit(&limit, cudaLimitMallocHeapSize);
  printf("cudaLimitMallocHeapSize: %u\n", (unsigned)limit);

  size_t uCurAvailMemoryInBytes;
  size_t uTotalMemoryInBytes;
  CUresult result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }*/

  printf("Nb individual %d / max DNA size %d\n",nb_indiv,host_max_dna_size);



  int block_size = 1 + host_max_dna_size / 65000;
  int y_dim_size = host_max_dna_size / block_size;
  int x_dim_size = 1024 * block_size;

  //int bucket_size = 1;
  int tmp_block = 1 + (( y_dim_size * 52 ) / 1024);

  int bucket_size = 19;
  int thread_number = bucket_size * 52;

  y_dim_size = 1 + y_dim_size / 19;

  dim3 dimGrid(x_dim_size,y_dim_size);

  //printf("Dim grid %d %d\n",x_dim_size,y_dim_size);
  //cudaDeviceSynchronize();
 // debug_dna<<<1,1>>>(dna_size, dna);
//  cudaDeviceSynchronize();
  //init_array<<<1024,1>>>(nb_promoters);
  //search_start_RNA<<<dimGrid,44>>>(dna_size,dna,dna_lead_promoter,dna_lag_promoter,nb_promoters,dynPromoterList,block_size);

  //debug_promoter_start<<<1,1>>>(dna_size,dna_lead_promoter,dna_lag_promoter,nb_promoters,
  //    1);
  //return;
  //cudaDeviceSynchronize();
  //search_stop_RNA<<<dimGrid,8>>>(dna_size,dna,dna_lead_term,dna_lag_term,block_size);
  //cudaDeviceSynchronize();
  /*result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }*/

  //debug_promoter_stop<<<1,1>>>(dna_size,dna_lead_promoter,dna_lag_promoter,nb_promoters,
  //    1);


  printf("X Dim %d Y Dim %d Thread %d Bucket %d Block %d\n",x_dim_size,y_dim_size,thread_number,bucket_size,block_size);

  search_start_stop_RNA_bucket<<<dimGrid,thread_number>>>(dna_size,dna,dna_lead_promoter,dna_lag_promoter,
      nb_promoters,dynPromoterList,dna_lead_term,dna_lag_term,bucket_size,block_size);

  //cudaDeviceSynchronize();

  init_RNA_struct<<<nb_indiv,1>>>(nb_indiv,rna,nb_promoters,max_nb_rna,idx_rna,max_nb_elements_rna_list);
  //cudaDeviceSynchronize();

  int max_promoters_host;
  HANDLE_ERROR(cudaMemcpy(&max_promoters_host,
                          max_nb_rna, sizeof(int32_t), cudaMemcpyDeviceToHost));

  //display_size_dna<<<1024,1>>>(dna_size);

  int threads_size = 1 + max_promoters_host / 1024;
  y_dim_size = max_promoters_host / threads_size;
  x_dim_size = 1024 * threads_size;

  printf("Max RNA %d (%d %d)\n",max_promoters_host,x_dim_size,y_dim_size);

  compute_RNA<<<x_dim_size,y_dim_size>>>(nb_indiv,dna_lead_promoter,dna_lag_promoter,dna_lead_term,dna_lag_term,
                            dna,dna_size,rna,idx_rna,nb_promoters,dynPromoterList,threads_size,y_dim_size);


  //cudaDeviceSynchronize();
  //debug_rna<<<1,1>>>(dna_size,dna_lead_promoter,dna_lag_promoter,rna,idx_rna,
  //    1020);

  //printf("Max RNA is computed\n");
  max_rna<<<nb_indiv,1>>>(idx_rna,max_nb_rna);
  //cudaDeviceSynchronize();
//  printf("Max RNA is transfering\n");

  int max_rna_host;
  HANDLE_ERROR(cudaMemcpy(&max_rna_host,
             max_nb_rna, sizeof(int32_t), cudaMemcpyDeviceToHost));

  //display_size_dna<<<1024,1>>>(dna_size);

  threads_size = 1 + max_rna_host / 1024;
  y_dim_size = max_rna_host / threads_size;
  x_dim_size = 1024 * threads_size;

  printf("Max RNA %d (%d %d)\n",max_rna_host,x_dim_size,y_dim_size);

  compute_start_protein<<<x_dim_size,y_dim_size>>>(idx_rna,rna,dna,dna_size,nb_protein,
      threads_size,y_dim_size);
  //cudaDeviceSynchronize();

/*
  result = cuMemGetInfo( &uCurAvailMemoryInBytes, &uTotalMemoryInBytes );
  if( result == CUDA_SUCCESS )
  {
    printf( "Total Memory: %d MB, Free Memory: %d MB\n",
            uTotalMemoryInBytes / ( 1024 * 1024 ),
            uCurAvailMemoryInBytes / ( 1024 * 1024 ));
  }
*/
  init_protein_struct<<<x_dim_size,y_dim_size>>>(nb_indiv,nb_protein,protein_list,
      idx_protein,rna,idx_rna,max_nb_protein,max_nb_elements_protein_list,threads_size,y_dim_size);
  //cudaDeviceSynchronize();

  int max_nb_protein_host;
  HANDLE_ERROR(cudaMemcpy(&max_nb_protein_host,
             max_nb_protein, sizeof(int32_t), cudaMemcpyDeviceToHost));

  threads_size = 1 + max_nb_protein_host / 1024;
  y_dim_size = max_nb_protein_host / threads_size;
  x_dim_size = max_rna_host * threads_size;

  block_size = 1 + x_dim_size / 65000;
  x_dim_size = x_dim_size / block_size;
  int z_dim_size = 1024 * block_size;

  dim3 dimGrid2(z_dim_size,x_dim_size);

  printf("Max Protein %d (%d %d %d)\n",max_nb_protein_host,z_dim_size,x_dim_size,y_dim_size);

  //display_size_dna<<<1024,1>>>(dna_size);
  //cudaDeviceSynchronize();

  compute_protein<<<dimGrid2,y_dim_size>>>(rna,protein_list,idx_protein,dna_size,dna,idx_rna, threads_size,y_dim_size,block_size);
  //cudaDeviceSynchronize();

  max_protein<<<nb_indiv,1>>>(max_nb_protein,idx_protein);
  HANDLE_ERROR(cudaMemcpy(&max_nb_protein_host,
             max_nb_protein, sizeof(int32_t), cudaMemcpyDeviceToHost));


  threads_size = 1 + max_nb_protein_host / 1024;
  y_dim_size = max_nb_protein_host / threads_size;
  x_dim_size = 1024 * threads_size;

  printf("Max Protein Updated %d (%d %d)\n",max_nb_protein_host,x_dim_size,y_dim_size);

  translate_protein<<<x_dim_size,y_dim_size>>>(w_max,idx_protein,protein_list,dna,dna_size,threads_size,y_dim_size);
  compute_phenotype<<<x_dim_size,y_dim_size>>>(idx_protein,protein_list,phenotype,threads_size,y_dim_size);
  compute_metaerror_fitness<<<nb_indiv,300>>>(selection_pressure,phenotype,target,metaerror,fitness);
  cudaDeviceSynchronize();

  /*if (AeTime::time()==14) {
    cudaDeviceSynchronize();
    cudaProfilerStop();
    exit(-1);
  }*/

}

__global__
void init_array(int* nb_promoters) {
  int indiv_id = blockIdx.x;

  nb_promoters[indiv_id] = 0;
  //atomicAdd(nb_promoters+ indiv_id,1);
  //atomicAdd(&nb_promoters[indiv_id],1);
  //printf("%d : %d\n",indiv_id,nb_promoters[indiv_id]);
}

__global__
void search_start_RNA(size_t* dna_size, char** dna, int8_t** dna_lead_promoter,
                      int8_t** dna_lag_promoter, int* nb_promoters,
                      pStruct** dynPromoterList, int block_size) {
  int indiv_id = blockIdx.x / block_size;
  int block_id = blockIdx.x % block_size;
  int pos_block_size = blockIdx.y;

  int dna_pos = gridDim.y*block_id+pos_block_size;

    __shared__ int dist_leading[22];
    __shared__ int dist_lagging[22];

    int motif_id = threadIdx.x;
    bool leading_or_lagging = true;

    if (dna_size[indiv_id] < PROM_SIZE) {
      //printf("START -- SMALL SIZE\n");
      /*if (dna_pos < dna_size[indiv_id] && threadIdx.x == 0) {
        dna_lead_promoter[indiv_id][dna_pos] = 0;
        dna_lag_promoter[indiv_id][dna_pos] = 0;
      }*/
    } else if (dna_pos < dna_size[indiv_id]) {
      if (motif_id >= 22) {
        // LAGGING
        motif_id -= 22;

        int pos = dna_pos - motif_id < 0 ?
                  dna_size[indiv_id] + (dna_pos - motif_id) :
                  dna_pos - motif_id;


        char s_motif=PROM_SEQ_LAG[motif_id];

        /*if (pos < 0 || pos > dna_size[indiv_id])
          printf("Checking DNA at %d (dna pos %d motif %d length %lu) for indiv %d (block %d)\n",
                 pos,dna_pos,motif_id,dna_size[indiv_id],indiv_id,block_id);*/

        char s_dna = dna[indiv_id][pos];

        dist_lagging[motif_id] =
            s_motif == s_dna ? 0 : 1;

        leading_or_lagging = true;
      } else {
        // LEADING
        int pos = dna_pos + motif_id >= dna_size[indiv_id] ?
                  dna_pos + motif_id - dna_size[indiv_id] : dna_pos + motif_id;

        dist_leading[motif_id] =
            PROM_SEQ_LEAD[motif_id] == dna[indiv_id][pos] ? 0 : 1;

        leading_or_lagging = false;
      }

      __syncthreads();

      if (threadIdx.x == 0) {
        int dist_lead = dist_leading[0] + dist_leading[1] + dist_leading[2] +
                        dist_leading[3] +
                        dist_leading[4] + dist_leading[5] + dist_leading[6] +
                        dist_leading[7] + dist_leading[8] +
                        dist_leading[9] + dist_leading[10] + dist_leading[11] +
                        dist_leading[12] + dist_leading[13] +
                        dist_leading[14] + dist_leading[15] + dist_leading[16] +
                        dist_leading[17] + dist_leading[18] +
                        dist_leading[19] + dist_leading[20] + dist_leading[21];

        int dist_lag = dist_lagging[0] + dist_lagging[1] + dist_lagging[2] +
                       dist_lagging[3] +
                       dist_lagging[4] + dist_lagging[5] + dist_lagging[6] +
                       dist_lagging[7] + dist_lagging[8] +
                       dist_lagging[9] + dist_lagging[10] + dist_lagging[11] +
                       dist_lagging[12] + dist_lagging[13] +
                       dist_lagging[14] + dist_lagging[15] + dist_lagging[16] +
                       dist_lagging[17] + dist_lagging[18] +
                       dist_lagging[19] + dist_lagging[20] + dist_lagging[21];


        //dna_lead_promoter[indiv_id][dna_pos] = dist_lead > 4 ? -1 : dist_lead;
 //       int nb_pro = dist_lead <= 4 ? 1 : 0;

        if (dist_lead <= 4) {
          int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

          dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
          dynPromoterList[indiv_id][rna_idx].leading_or_lagging = true;
          dynPromoterList[indiv_id][rna_idx].error = dist_lead;
        }

        //dna_lag_promoter[indiv_id][dna_pos] = dist_lag > 4 ? -1 : dist_lag;
//        nb_pro += dist_lag <= 4 ? 1 : 0;

        //    if (indiv_id == 410 && dna_pos == 10)
        //      printf("Promoter found at %d ; %d %d \n",dna_pos,dist_lead,dist_lag);
        if (dist_lag <= 4) {
          int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

          dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
          dynPromoterList[indiv_id][rna_idx].leading_or_lagging = false;
          dynPromoterList[indiv_id][rna_idx].error = dist_lag;
        }

      }

  }
}


__global__
void search_start_RNA_bucket(size_t* dna_size, char** dna, int8_t** dna_lead_promoter,
                      int8_t** dna_lag_promoter, int* nb_promoters,
                      pStruct** dynPromoterList, int bucket_size, int block_size) {
  int indiv_id = blockIdx.x / block_size;
  int block_id = blockIdx.x % block_size;
  int pos_block_size = blockIdx.y;

  int dna_pos = gridDim.y*block_id+pos_block_size;

  __shared__ int dist_leading[BUCKET_MAX_SIZE][22];
  __shared__ int dist_lagging[BUCKET_MAX_SIZE][22];

  int motif_id = threadIdx.x % 22;
  int dna_global_offset = threadIdx.x / 22;
  dna_pos+=dna_global_offset;

  if (dna_pos < dna_size[indiv_id] && dna_size[indiv_id] >= PROM_SIZE) {
      if (motif_id >= 22) {
        // LAGGING
        int t_motif_id = motif_id - 22;
        dist_lagging[dna_global_offset][t_motif_id] =
            PROM_SEQ_LAG[t_motif_id] == dna[indiv_id][dna_pos - t_motif_id < 0 ?
                                                    dna_size[indiv_id] + (dna_pos - t_motif_id) :
                                                    dna_pos - t_motif_id] ? 0 : 1;
      } else {
        // LEADING
        dist_leading[dna_global_offset][motif_id] =
            PROM_SEQ_LEAD[motif_id] == dna[indiv_id][dna_pos + motif_id >= dna_size[indiv_id] ?
                                                     dna_pos + motif_id - dna_size[indiv_id] : dna_pos + motif_id] ? 0 : 1;
      }


    __syncthreads();

    if (motif_id == 0) {

      int dist_lead = dist_leading[dna_global_offset][0] +
                      dist_leading[dna_global_offset][1] +
                      dist_leading[dna_global_offset][2] +
                      dist_leading[dna_global_offset][3] +
                      dist_leading[dna_global_offset][4] +
                      dist_leading[dna_global_offset][5] +
                      dist_leading[dna_global_offset][6] +
                      dist_leading[dna_global_offset][7] +
                      dist_leading[dna_global_offset][8] +
                      dist_leading[dna_global_offset][9] +
                      dist_leading[dna_global_offset][10] +
                      dist_leading[dna_global_offset][11] +
                      dist_leading[dna_global_offset][12] +
                      dist_leading[dna_global_offset][13] +
                      dist_leading[dna_global_offset][14] +
                      dist_leading[dna_global_offset][15] +
                      dist_leading[dna_global_offset][16] +
                      dist_leading[dna_global_offset][17] +
                      dist_leading[dna_global_offset][18] +
                      dist_leading[dna_global_offset][19] +
                      dist_leading[dna_global_offset][20] +
                      dist_leading[dna_global_offset][21];

      if (dist_lead <= 4) {
        int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

        dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
        dynPromoterList[indiv_id][rna_idx].leading_or_lagging = true;
        dynPromoterList[indiv_id][rna_idx].error = dist_lead;
      }
    }

    if (motif_id == 22) {
      int dist_lag = dist_lagging[dna_global_offset][0] +
                     dist_lagging[dna_global_offset][1] +
                     dist_lagging[dna_global_offset][2] +
                     dist_lagging[dna_global_offset][3] +
                     dist_lagging[dna_global_offset][4] +
                     dist_lagging[dna_global_offset][5] +
                     dist_lagging[dna_global_offset][6] +
                     dist_lagging[dna_global_offset][7] +
                     dist_lagging[dna_global_offset][8] +
                     dist_lagging[dna_global_offset][9] +
                     dist_lagging[dna_global_offset][10] +
                     dist_lagging[dna_global_offset][11] +
                     dist_lagging[dna_global_offset][12] +
                     dist_lagging[dna_global_offset][13] +
                     dist_lagging[dna_global_offset][14] +
                     dist_lagging[dna_global_offset][15] +
                     dist_lagging[dna_global_offset][16] +
                     dist_lagging[dna_global_offset][17] +
                     dist_lagging[dna_global_offset][18] +
                     dist_lagging[dna_global_offset][19] +
                     dist_lagging[dna_global_offset][20] +
                     dist_lagging[dna_global_offset][21];

      if (dist_lag <= 4) {
        int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

        dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
        dynPromoterList[indiv_id][rna_idx].leading_or_lagging = false;
        dynPromoterList[indiv_id][rna_idx].error = dist_lag;
      }
    }
  }
}


__global__
void search_start_stop_RNA_bucket(size_t* dna_size, char** dna, int8_t** dna_lead_promoter,
                             int8_t** dna_lag_promoter, int* nb_promoters,
                             pStruct** dynPromoterList,
                             int8_t** dna_lead_term, int8_t** dna_lag_term,
                             int bucket_size, int block_size) {
  int indiv_id = blockIdx.x / block_size;
  int block_id = blockIdx.x % block_size;
  int pos_block_size = blockIdx.y;

  int org_dna_pos = (gridDim.y*block_id+pos_block_size)*bucket_size;

  __shared__ int prom_dist_leading[BUCKET_MAX_SIZE][26];
  __shared__ int prom_dist_lagging[BUCKET_MAX_SIZE][26];


  __shared__ int term_dist_leading[BUCKET_MAX_SIZE][4];
  __shared__ int term_dist_lagging[BUCKET_MAX_SIZE][4];

  __shared__ int cached_dna[62];

  //__shared__ int dist_lead[BUCKET_MAX_SIZE];
  //__shared__ int dist_lag[BUCKET_MAX_SIZE];

  int motif_id = threadIdx.x % 52;
  int dna_global_offset = threadIdx.x / 52;

  int cached_dna_pos = 21+dna_global_offset;
  int dna_pos = org_dna_pos+dna_global_offset;


  if (dna_pos < dna_size[indiv_id] && dna_size[indiv_id] >= PROM_SIZE) {

    if (threadIdx.x < 62) {

      int load_pos = org_dna_pos-21+threadIdx.x;
      load_pos = load_pos < 0 ? dna_size[indiv_id] + load_pos : load_pos;
      load_pos = load_pos >= dna_size[indiv_id] ?
                     load_pos - dna_size[indiv_id] :
                     load_pos;
      //printf("Loading into shared memory at %d : DNA pos %d (w/o mod %d size %lu)\n",threadIdx.x,load_pos,org_dna_pos-21+threadIdx.x,dna_size[indiv_id]);
      cached_dna[threadIdx.x] = dna[indiv_id][load_pos];
    }



    /*if (indiv_id == 0 && blockIdx.x == 0 && blockIdx.y == 0)
      printf("Thread %d Block X %d block Y %d -- Motif ID %d (offset %d) : DNA POS %d (cached %d min %d max %d)\n",
             threadIdx.x,blockIdx.x,blockIdx.y,motif_id,dna_global_offset,dna_pos,cached_dna_pos,cached_dna_pos-21,cached_dna_pos+21);*/

    __syncthreads();

    if (motif_id >= 26 && motif_id < 48) {
      // LAGGING
      int t_motif_id = motif_id - 26;
      prom_dist_lagging[dna_global_offset][t_motif_id] =
          PROM_SEQ_LAG[t_motif_id] == cached_dna[cached_dna_pos-t_motif_id] ? 0 : 1;
    } else if (motif_id < 22) {
      // LEADING
      prom_dist_leading[dna_global_offset][motif_id] =
          PROM_SEQ_LEAD[motif_id] == cached_dna[cached_dna_pos+motif_id] ? 0 : 1;
    } else if (motif_id >= 22 && motif_id < 26) {
      int t_motif_id = motif_id - 22;
      // LEADING
      term_dist_leading[dna_global_offset][t_motif_id] =
          cached_dna[cached_dna_pos+t_motif_id] != cached_dna[cached_dna_pos-t_motif_id+10] ? 1 : 0;
    } else {
      int t_motif_id = motif_id - 48;
      term_dist_lagging[dna_global_offset][t_motif_id] =
          cached_dna[cached_dna_pos-t_motif_id] != cached_dna[cached_dna_pos+t_motif_id-10] ? 1 : 0;
    }

    __syncthreads();

    if (motif_id % 52 == 0) {

      int dist_lead = prom_dist_leading[dna_global_offset][0] +
          prom_dist_leading[dna_global_offset][1] +
          prom_dist_leading[dna_global_offset][2] +
          prom_dist_leading[dna_global_offset][3] +
          prom_dist_leading[dna_global_offset][4] +
          prom_dist_leading[dna_global_offset][5] +
          prom_dist_leading[dna_global_offset][6] +
          prom_dist_leading[dna_global_offset][7] +
          prom_dist_leading[dna_global_offset][8] +
          prom_dist_leading[dna_global_offset][9] +
          prom_dist_leading[dna_global_offset][10] +
          prom_dist_leading[dna_global_offset][11] +
          prom_dist_leading[dna_global_offset][12] +
          prom_dist_leading[dna_global_offset][13] +
          prom_dist_leading[dna_global_offset][14] +
          prom_dist_leading[dna_global_offset][15] +
          prom_dist_leading[dna_global_offset][16] +
          prom_dist_leading[dna_global_offset][17] +
          prom_dist_leading[dna_global_offset][18] +
          prom_dist_leading[dna_global_offset][19] +
          prom_dist_leading[dna_global_offset][20] +
          prom_dist_leading[dna_global_offset][21];

      if (dist_lead <= 4) {
        int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

        dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
        dynPromoterList[indiv_id][rna_idx].leading_or_lagging = true;
        dynPromoterList[indiv_id][rna_idx].error = dist_lead;

//        if (indiv_id == 0) printf("New Start RNA Found ! POS %d RNA Idx %d\n",dna_pos,rna_idx);
      }
    }
    else if (motif_id % 52 == 22) {
      int dist_lead = term_dist_leading[dna_global_offset][0] +
          term_dist_leading[dna_global_offset][1] +
          term_dist_leading[dna_global_offset][2] +
          term_dist_leading[dna_global_offset][3];
      dna_lead_term[indiv_id][dna_pos] = dist_lead == 4 ? 1 : 0;
  /*    if (dna_lead_term[indiv_id][dna_pos] == 4)
        printf("New STOP RNA Found ! POS %d RNA Idx %d\n",dna_pos);*/
    }
    else if (motif_id % 52 == 26) {
      int dist_lag = prom_dist_lagging[dna_global_offset][0] +
          prom_dist_lagging[dna_global_offset][1] +
          prom_dist_lagging[dna_global_offset][2] +
          prom_dist_lagging[dna_global_offset][3] +
          prom_dist_lagging[dna_global_offset][4] +
          prom_dist_lagging[dna_global_offset][5] +
          prom_dist_lagging[dna_global_offset][6] +
          prom_dist_lagging[dna_global_offset][7] +
          prom_dist_lagging[dna_global_offset][8] +
          prom_dist_lagging[dna_global_offset][9] +
          prom_dist_lagging[dna_global_offset][10] +
          prom_dist_lagging[dna_global_offset][11] +
          prom_dist_lagging[dna_global_offset][12] +
          prom_dist_lagging[dna_global_offset][13] +
          prom_dist_lagging[dna_global_offset][14] +
          prom_dist_lagging[dna_global_offset][15] +
          prom_dist_lagging[dna_global_offset][16] +
          prom_dist_lagging[dna_global_offset][17] +
          prom_dist_lagging[dna_global_offset][18] +
          prom_dist_lagging[dna_global_offset][19] +
          prom_dist_lagging[dna_global_offset][20] +
          prom_dist_lagging[dna_global_offset][21];

      if (dist_lag <= 4) {
        int rna_idx = atomicAdd(nb_promoters + indiv_id, 1);

        dynPromoterList[indiv_id][rna_idx].pos = dna_pos;
        dynPromoterList[indiv_id][rna_idx].leading_or_lagging = false;
        dynPromoterList[indiv_id][rna_idx].error = dist_lag;
      }
    }
    else if (motif_id % 52 == 48) {
      int dist_lag = term_dist_lagging[dna_global_offset][0] +
          term_dist_lagging[dna_global_offset][1] +
          term_dist_lagging[dna_global_offset][2] +
          term_dist_lagging[dna_global_offset][3];
      dna_lag_term[indiv_id][dna_pos] = dist_lag == 4 ? 1 : 0;
    }
  }
}


__global__
void search_stop_RNA(size_t* dna_size, char** dna, int8_t** dna_lead_term, int8_t** dna_lag_term, int block_size) {

  int indiv_id = blockIdx.x / block_size;
  int block_id = blockIdx.x % block_size;
  int pos_block_size = blockIdx.y;

  int dna_pos = gridDim.y * block_id + pos_block_size;

  __shared__ int dist_leading[4];
  __shared__ int dist_lagging[4];

  int motif_id = threadIdx.x;

  if (dna_size[indiv_id] < PROM_SIZE) {
    //printf("STOP -- SMALL SIZE\n");
    /*if (dna_pos < dna_size[indiv_id] && threadIdx.x == 0) {
      dna_lead_term[indiv_id][dna_pos] = 0;
      dna_lag_term[indiv_id][dna_pos] = 0;
    }*/
  } else if (dna_pos < dna_size[indiv_id]) {

    if (motif_id >= 4) {
      // LAGGING
      motif_id -= 4;
      int pos_1 = dna_pos - motif_id < 0 ?
                  dna_size[indiv_id] + (dna_pos - motif_id) : dna_pos -
                                                              motif_id;
      int pos_2 = dna_pos + motif_id - 10 < 0 ?
                  dna_size[indiv_id] + (dna_pos + motif_id - 10) : dna_pos +
                                                                   motif_id -
                                                                   10;

      dist_lagging[motif_id] =
          dna[indiv_id][pos_1] != dna[indiv_id][pos_2] ? 1 : 0;
    } else {
      // LEADING
      int pos_1 = dna_pos + motif_id >= dna_size[indiv_id] ?
                  (dna_pos + motif_id) - dna_size[indiv_id] : dna_pos +
                                                              motif_id;
      int pos_2 = dna_pos - motif_id + 10 >= dna_size[indiv_id] ?
                  10 + dna_pos - motif_id - dna_size[indiv_id] : dna_pos -
                                                                 motif_id +
                                                                 10;
      dist_leading[motif_id] =
          dna[indiv_id][pos_1] != dna[indiv_id][pos_2] ? 1 : 0;
    }

    __syncthreads();

    if (threadIdx.x == 0) {
      int dist_lead = dist_leading[0] + dist_leading[1] + dist_leading[2] +
                      dist_leading[3];
      int dist_lag = dist_lagging[0] + dist_lagging[1] + dist_lagging[2] +
                     dist_lagging[3];

      dna_lead_term[indiv_id][dna_pos] = dist_lead == 4 ? 1 : 0;
      dna_lag_term[indiv_id][dna_pos] = dist_lag == 4 ? 1 : 0;

      //printf("%d Found STOP at %d\n",indiv_id,dna_pos);

    }
  }
}

__global__
void internal_init_RNA_struct(cRNA*** rna, int32_t* max_nb_elements_rna_list, int indiv_id) {
  int offset = blockIdx.x*1024;
  int rna_idx = offset+threadIdx.x;

  if (rna_idx < max_nb_elements_rna_list[indiv_id]) {
    rna[indiv_id][rna_idx] = (cRNA*) malloc(sizeof(cRNA));
    rna[indiv_id][rna_idx]->max_protein_elements = 200;

    rna[indiv_id][rna_idx]->start_prot = (uint32_t*) malloc(
        rna[indiv_id][rna_idx]->max_protein_elements * sizeof(uint32_t));
  }
}

__global__
void init_RNA_struct(int pop_size, cRNA*** rna, int* nb_promoters, int32_t* max_nb_rna,int32_t* idx_rna, int32_t* max_nb_elements_rna_list) {
  int indiv_id = blockIdx.x;

  if (nb_promoters[indiv_id] > 0) {
    if (nb_promoters[indiv_id] >= max_nb_elements_rna_list[indiv_id]) {
      // Increase RNA List size
      for (int i=0; i < max_nb_elements_rna_list[indiv_id];i++) {
        free(rna[indiv_id][i]->start_prot);
        free(rna[indiv_id][i]);
      }

      free(rna[indiv_id]);
      int before_cpt=max_nb_elements_rna_list[indiv_id];
      max_nb_elements_rna_list[indiv_id] = (1+((int32_t)nb_promoters[indiv_id]/RNA_LIST_INCR_SIZE))*RNA_LIST_INCR_SIZE;

      rna[indiv_id] = (cRNA**) malloc((max_nb_elements_rna_list[indiv_id]+1)*sizeof(cRNA*));


      int block_offset = max_nb_elements_rna_list[indiv_id]/1024 + 1;

      internal_init_RNA_struct<<<block_offset,1024>>>(rna,max_nb_elements_rna_list,indiv_id);
      //cudaDeviceSynchronize();
      //
      /*for (int i=0; i < max_nb_elements_rna_list[indiv_id]; i++) {
        rna[indiv_id][i] = (cRNA*) malloc(sizeof(cRNA));
        //printf("Malloc POINTER RNA %d indiv %d : %p\n",max_nb_elements_rna_list[indiv_id],indiv_id,l_rna);
        rna[indiv_id][i]->max_protein_elements = 250;//RNA_LIST_PROTEIN_INCR_SIZE;

        rna[indiv_id][i]->start_prot =  (int*) malloc(
            rna[indiv_id][i]->max_protein_elements  * sizeof(int));
      }*/
/*      printf("Malloc Increase DONE RNA %d indiv %d (before %d current %d)\n",max_nb_elements_rna_list[indiv_id],indiv_id,
             before_cpt,nb_promoters[indiv_id]);*/
    } else if (nb_promoters[indiv_id] < max_nb_elements_rna_list[indiv_id]/2 && max_nb_elements_rna_list[indiv_id] - RNA_LIST_INCR_SIZE > 0) {
      // Decrease RNA List size
      for (int i=0; i < max_nb_elements_rna_list[indiv_id];i++) {
        free(rna[indiv_id][i]->start_prot);
        free(rna[indiv_id][i]);
      }
      free(rna[indiv_id]);
      //max_nb_elements_rna_list[indiv_id] -= RNA_LIST_INCR_SIZE;
      int before_cpt=max_nb_elements_rna_list[indiv_id];
      max_nb_elements_rna_list[indiv_id]  = max_nb_elements_rna_list[indiv_id] - RNA_LIST_INCR_SIZE == 0 ?
                                            RNA_LIST_INCR_SIZE :  max_nb_elements_rna_list[indiv_id] - RNA_LIST_INCR_SIZE;

      rna[indiv_id] = (cRNA**) malloc((max_nb_elements_rna_list[indiv_id]+1)*sizeof(cRNA*));

      int block_offset = max_nb_elements_rna_list[indiv_id]/1024 + 1;

      internal_init_RNA_struct<<<block_offset,1024>>>(rna,max_nb_elements_rna_list,indiv_id);
      //cudaDeviceSynchronize();

      /*for (int i=0; i < max_nb_elements_rna_list[indiv_id];i++) {
        rna[indiv_id][i] = (cRNA*) malloc(sizeof(cRNA));
        rna[indiv_id][i]->max_protein_elements = 250;
        rna[indiv_id][i]->start_prot =  (int*) malloc(
            rna[indiv_id][i]->max_protein_elements  * sizeof(int));
        //RNA_LIST_PROTEIN_INCR_SIZE;
      }*/
      /*printf("Malloc Decrease DONE RNA %d indiv %d (before %d current %d)\n",max_nb_elements_rna_list[indiv_id],indiv_id,
              before_cpt,nb_promoters[indiv_id]);*/
    }

    atomicMax(max_nb_rna,nb_promoters[indiv_id]+1);
    idx_rna[indiv_id] = 0;
  }
}

__global__
void compute_RNA(int pop_size, int8_t** dna_lead_promoter,
                 int8_t** dna_lag_promoter, int8_t** dna_lead_term,
                 int8_t** dna_lag_term, char** dna, size_t* dna_size,
                 cRNA*** rna, int32_t* idx_rna,  int* nb_promoters,
                 pStruct** dynPromoterList, int threads_size, int thread_dim) {

  int indiv_id = blockIdx.x / threads_size;
  int block_id = blockIdx.x % threads_size;
  int pos_block_size = threadIdx.x;

  int rna_idx = thread_dim*block_id+pos_block_size;

  if (dna_size[indiv_id] >= PROM_SIZE && rna_idx < nb_promoters[indiv_id]) {
    if (dynPromoterList[indiv_id][rna_idx].leading_or_lagging) {
      // LEADING
      // Search for terminator
      int k = dynPromoterList[indiv_id][rna_idx].pos + 22;
      k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;
      int k_end = k;
      do {
        if (dna_lead_term[indiv_id][k] == 1) {
          int32_t rna_end =
              k + 10 >= dna_size[indiv_id] ? k + 10 - dna_size[indiv_id] :
              k +
              10;

          int32_t rna_length = 0;

          if (dynPromoterList[indiv_id][rna_idx].pos > rna_end)
            rna_length = dna_size[indiv_id] - dynPromoterList[indiv_id][rna_idx].pos + rna_end;
          else
            rna_length = rna_end - dynPromoterList[indiv_id][rna_idx].pos;

          if (rna_length < 19) {
            break;
          }

          int idx = atomicAdd(idx_rna + indiv_id, 1);

          rna[indiv_id][idx]->begin = dynPromoterList[indiv_id][rna_idx].pos;
          rna[indiv_id][idx]->end = rna_end;
          rna[indiv_id][idx]->length = rna_length;
          rna[indiv_id][idx]->leading_lagging = !dynPromoterList[indiv_id][rna_idx].leading_or_lagging;

          rna[indiv_id][idx]->e = 1.0 -
                     fabs(((float) dynPromoterList[indiv_id][rna_idx].error)) /
                     5.0;

          /*if (rna_length > rna[indiv_id][idx]->max_protein_elements) {
            // Increase size
            free(rna[indiv_id][idx]->start_prot);
            int before_cpt = rna[indiv_id][idx]->max_protein_elements;
            rna[indiv_id][idx]->max_protein_elements=(1+((int32_t)rna_length/RNA_LIST_PROTEIN_INCR_SIZE))*RNA_LIST_PROTEIN_INCR_SIZE;
            rna[indiv_id][idx]->start_prot = (uint32_t*) malloc(
                (rna[indiv_id][idx]->max_protein_elements + 1) * sizeof(uint32_t));
            printf("Malloc Increase DONE RNA Protein List %d indiv %d -- %d (before %d current %d)\n",rna_length,
                   indiv_id,rna_idx,
                   before_cpt,rna[indiv_id][idx]->max_protein_elements);
          } else if ((rna_length < rna[indiv_id][idx]->max_protein_elements/2) && (rna[indiv_id][idx]->max_protein_elements - RNA_LIST_PROTEIN_INCR_SIZE > 0)) {
            // Decrease size
            free(rna[indiv_id][idx]->start_prot);
            int before_cpt = rna[indiv_id][idx]->max_protein_elements;
            rna[indiv_id][idx]->max_protein_elements-=RNA_LIST_PROTEIN_INCR_SIZE;
            rna[indiv_id][idx]->start_prot = (int*) malloc(
                (rna[indiv_id][idx]->max_protein_elements + 1) * sizeof(int));
            if (indiv_id == 828)  printf("Malloc Decrease DONE RNA Protein List %d indiv %d -- %d (before %d current %d)\n",rna_length,
                   indiv_id,rna_idx,
                   before_cpt,rna[indiv_id][idx]->max_protein_elements);
          }*/

          rna[indiv_id][idx]->start_lenght = 1;
          rna[indiv_id][idx]->nb_protein = 0;

          break;
        }

        k++;
        k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;
      } while (k != k_end);

    } else {
      // LAGGING

      // Search for terminator
      int k = dynPromoterList[indiv_id][rna_idx].pos - 22;
      k = k < 0 ? dna_size[indiv_id] + k : k;
      int k_end = k;
      do {

        if (dna_lag_term[indiv_id][k] == 1) {
          int32_t rna_end = k - 10 < 0 ? dna_size[indiv_id] + (k - 10) : k - 10;

          int32_t rna_length = 0;

          if (dynPromoterList[indiv_id][rna_idx].pos < rna_end)
            rna_length = dynPromoterList[indiv_id][rna_idx].pos + dna_size[indiv_id] - rna_end;
          else
            rna_length = dynPromoterList[indiv_id][rna_idx].pos - rna_end;

          if (rna_length < 19) {
            break;
          }

          int idx = atomicAdd(idx_rna + indiv_id, 1);

          /*printf("Indiv %d -- Setting RNA %d to begin at %d (promoter idx %d out of %d)\n",
                 indiv_id,idx,dynPromoterList[indiv_id][rna_idx].pos,rna_idx,nb_promoters[indiv_id]);*/

          rna[indiv_id][idx]->begin = dynPromoterList[indiv_id][rna_idx].pos;
          rna[indiv_id][idx]->end = rna_end;
          rna[indiv_id][idx]->length = rna_length;

          rna[indiv_id][idx]->leading_lagging = !dynPromoterList[indiv_id][rna_idx].leading_or_lagging;


          rna[indiv_id][idx]->e =
              1.0 - ((float) dynPromoterList[indiv_id][rna_idx].error) / 5.0;


          /*if (rna_length > rna[indiv_id][idx]->max_protein_elements) {
            // Increase size
            free(rna[indiv_id][idx]->start_prot);
            int before_cpt = rna[indiv_id][idx]->max_protein_elements;
            rna[indiv_id][idx]->max_protein_elements=(1+((int32_t)rna_length/RNA_LIST_PROTEIN_INCR_SIZE))*RNA_LIST_PROTEIN_INCR_SIZE;
            rna[indiv_id][idx]->start_prot = (uint32_t*) malloc(
                (rna[indiv_id][idx]->max_protein_elements + 1) * sizeof(uint32_t));
            printf("Malloc Increase DONE RNA Protein List %d indiv %d -- %d (before %d current %d)\n",
                   rna_length,
                   indiv_id,rna_idx,
                   before_cpt,rna[indiv_id][idx]->max_protein_elements);
          } else if ((rna_length < rna[indiv_id][idx]->max_protein_elements/2) && (rna[indiv_id][idx]->max_protein_elements - RNA_LIST_PROTEIN_INCR_SIZE > 0)) {
            // Decrease size
            free(rna[indiv_id][idx]->start_prot);
            int before_cpt = rna[indiv_id][idx]->max_protein_elements;
            rna[indiv_id][idx]->max_protein_elements-=RNA_LIST_PROTEIN_INCR_SIZE;
            rna[indiv_id][idx]->start_prot = (int*) malloc(
                (rna[indiv_id][idx]->max_protein_elements + 1) * sizeof(int));
            if (indiv_id == 828)   printf("Malloc Decrease DONE RNA Protein List %d indiv %d -- %d (before %d current %d)\n",
                   rna_length,
                   indiv_id,rna_idx,
                   before_cpt,rna[indiv_id][idx]->max_protein_elements);
          }*/

          rna[indiv_id][idx]->start_lenght = 1;
          rna[indiv_id][idx]->nb_protein = 0;

          break;
        }

        k--;
        k = k < 0 ? dna_size[indiv_id] + k : k;
      } while (k != k_end);

    }
  }
}

__global__ void max_rna(int32_t* idx_rna, int32_t* max_nb_rna) {
  int indiv_id = blockIdx.x;

  atomicMax(max_nb_rna, idx_rna[indiv_id]);
}

__global__ void compute_start_protein(int32_t* idx_rna, cRNA*** rna,
                                      char** dna,size_t* dna_size, int32_t* nb_protein, int threads_size, int thread_dim) {
//  int indiv_id = blockIdx.x;
//  int rna_idx = threadIdx.x;

  int indiv_id = blockIdx.x / threads_size;
  int block_id = blockIdx.x % threads_size;
  int pos_block_size = threadIdx.x;

  int rna_idx = thread_dim*block_id+pos_block_size;


  /*if (indiv_id == 737 && rna_idx == 47)
    printf("Searching for protein start on %d (%d %d) with %d\n",rna_idx,indiv_id,block_id,pos_block_size);*/

  if (rna_idx < idx_rna[indiv_id]) {

    int c_pos = rna[indiv_id][rna_idx]->begin;
    if (rna[indiv_id][rna_idx]->length > 22) {
      if (rna[indiv_id][rna_idx]->leading_lagging == 0) {
        c_pos += 22;
        c_pos =
            c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
      } else {
        c_pos -= 22;
        //if (indiv_id == 915 && rna_idx < 2) printf("C_POS -22 : %d\n", c_pos);
        c_pos = c_pos < 0 ? ((int) dna_size[indiv_id]) + c_pos : c_pos;
        /*if (indiv_id == 737 && rna_idx == 47)
          printf("MOD C_POS -22 : %d SIZE - MOD : %d SIZE: %lu\n",
                 c_pos, (int) (dna_size[indiv_id] - c_pos),
                 dna_size[indiv_id]);*/

      }

      // TODO IF SIZE is smaller than 9 return

      /*if (indiv_id == 737) {
        printf("Starting to search in between %d (%d -- %lu) and %d\n",
               c_pos, rna[indiv_id][rna_idx]->begin, dna_size[indiv_id],
               rna[indiv_id][rna_idx]->end);
      }*/
      while (c_pos != rna[indiv_id][rna_idx]->end) {

        bool start = false;
        int t_pos, k_t;

        if (rna[indiv_id][rna_idx]->leading_lagging == 0) {
          // Search for Shine Dalgarro + START codon on LEADING
          for (int k = 0; k < 9; k++) {
            k_t = k >= 6 ? k + 4 : k;
            t_pos = c_pos + k_t >= dna_size[indiv_id] ? c_pos + k_t -
                                                        dna_size[indiv_id] :
                    c_pos + k_t;

            if (dna[indiv_id][t_pos] == SHINE_DAL_SEQ_LEAD[k]) {
              start = true;
            } else {
              start = false;
              break;
            }
          }

        } else {
          // Search for Shine Dalgarro + START codon on LAGGING
          for (int k = 0; k < 9; k++) {
            k_t = k >= 6 ? k + 4 : k;
            t_pos =
                c_pos - k_t < 0 ? dna_size[indiv_id] - c_pos - k_t : c_pos -
                                                                     k_t;

            /*if (indiv_id == 737 && rna_idx == 47)
              printf("Search protein start at %d : %d %d -- %c / %c\n",c_pos,t_pos,k_t,
                     dna[indiv_id][t_pos],SHINE_DAL_SEQ_LAG[k]);*/

            if (dna[indiv_id][t_pos] == SHINE_DAL_SEQ_LAG[k]) {
              start = true;
            } else {
              start = false;
              break;
            }
          }
        }

        if (start) {
          int prot_idx = atomicAdd(&(rna[indiv_id][rna_idx]->nb_protein), 1);

          rna[indiv_id][rna_idx]->start_prot[prot_idx] = c_pos;

          /*if (indiv_id == 737 && rna_idx == 47)
            printf("RNA : %d %d || Current : %d -- %lu\n",
                   rna[indiv_id][rna_idx]->begin,
                   rna[indiv_id][rna_idx]->end,
                   c_pos, dna_size[indiv_id]);*/


          atomicAdd(&nb_protein[indiv_id], 1);
        }

        if (rna[indiv_id][rna_idx]->leading_lagging == 0) {
          c_pos++;
          c_pos =
              c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
        } else {
          c_pos--;
          c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
        }

        /*if (indiv_id == 915 && rna_idx < 2) {
          printf("POS %d END %d\n", c_pos, rna[indiv_id][rna_idx]->end);
        }*/
      }
    } /*else if (indiv_id == 915 && rna_idx < 2) {
      printf("TOO SMALL !!!\n");
    }*/
  }


}

__global__ void init_protein_struct(int pop_size, int32_t* nb_protein,
                                    cProtein** protein_list, int32_t* idx_protein,
                                    cRNA*** rna, int32_t* idx_rna, int32_t* max_nb_protein, int* max_nb_elements_protein_list,
                                    int threads_size, int thread_dim) {
  //int indiv_id = blockIdx.x;
  //int rna_idx = threadIdx.x;

  int indiv_id = blockIdx.x / threads_size;
  int block_id = blockIdx.x % threads_size;
  int pos_block_size = threadIdx.x;

  int rna_idx = thread_dim*block_id+pos_block_size;

  if (rna_idx == 0 && nb_protein[indiv_id] > 0) {
    //printf("%d -- Number of Protein %d (array size %d)\n",indiv_id,nb_protein[indiv_id],max_nb_elements_protein_list[indiv_id]);
    /*if (max_nb_elements_protein_list[indiv_id] == 0) {
      max_nb_elements_protein_list[indiv_id] = (1+((int32_t)nb_protein[indiv_id]/PROTEIN_LIST_INCR_SIZE))*PROTEIN_LIST_INCR_SIZE;

      protein_list[indiv_id] = (cProtein*) malloc((max_nb_elements_protein_list[indiv_id]+1)*sizeof(cProtein));

      printf("Malloc Init DONE Protein %d indiv %d (current %d)\n",max_nb_elements_protein_list[indiv_id],indiv_id,
             nb_protein[indiv_id]);
    } else */if (nb_protein[indiv_id] >= max_nb_elements_protein_list[indiv_id]) {
      free(protein_list[indiv_id]);
      int before_cpt=max_nb_elements_protein_list[indiv_id];
      max_nb_elements_protein_list[indiv_id] = (1+((int32_t)nb_protein[indiv_id]/PROTEIN_LIST_INCR_SIZE))*PROTEIN_LIST_INCR_SIZE;

      protein_list[indiv_id] = (cProtein*) malloc((max_nb_elements_protein_list[indiv_id]+1)*sizeof(cProtein));

      /*printf("Malloc Increase DONE Protein %d indiv %d (before %d current %d)\n",max_nb_elements_protein_list[indiv_id],indiv_id,
             before_cpt,nb_protein[indiv_id]);*/
    } else if (nb_protein[indiv_id] < max_nb_elements_protein_list[indiv_id]/2 && (max_nb_elements_protein_list[indiv_id] - PROTEIN_LIST_INCR_SIZE > 0)) {
      free(protein_list[indiv_id]);
      //max_nb_elements_rna_list[indiv_id] -= RNA_LIST_INCR_SIZE;
      int before_cpt=max_nb_elements_protein_list[indiv_id];
      max_nb_elements_protein_list[indiv_id]  = max_nb_elements_protein_list[indiv_id] - PROTEIN_LIST_INCR_SIZE == 0 ?
                                            PROTEIN_LIST_INCR_SIZE :  max_nb_elements_protein_list[indiv_id] - PROTEIN_LIST_INCR_SIZE;

      protein_list[indiv_id] = (cProtein*) malloc((max_nb_elements_protein_list[indiv_id]+1)*sizeof(cProtein));

      /*printf("Malloc Decrease DONE Protein %d indiv %d (before %d current %d)\n",max_nb_elements_protein_list[indiv_id],indiv_id,
             before_cpt,nb_protein[indiv_id]);*/
    }
    //printf("%d -- END OF Number of Protein %d (array size %d)\n",indiv_id,nb_protein[indiv_id],max_nb_elements_protein_list[indiv_id]);

    //protein_list[indiv_id] = (cProtein**) malloc((nb_protein[indiv_id])*sizeof(cProtein*));
    idx_protein[indiv_id] = 0;
  }

  if (rna_idx < idx_rna[indiv_id])
    atomicMax(max_nb_protein,nb_protein[indiv_id]);
}

__global__ void display_size_dna(size_t* dna_size) {
  int indiv_id = blockIdx.x;

  //cProtein* l_protein = (cProtein*) malloc(sizeof(22));
  //if (indiv_id ==0) printf("Address : %lu\n", (unsigned long) l_protein);
}

__global__ void compute_protein(cRNA*** rna, cProtein** protein_list, int32_t* idx_protein,
                                size_t* dna_size,char** dna,int32_t* idx_rna, int threads_size, int thread_dim, int block_size) {
  //int indiv_id = blockIdx.x;

  int debug_iid = 915;
  //int rna_idx = blockIdx.y;

  //int protein_idx = threadIdx.x;

  int indiv_id = blockIdx.x / block_size;
  int block_id = blockIdx.x % block_size;
  int pos_block_size = blockIdx.y;

  int block_x = gridDim.y*block_id+pos_block_size;

  int rna_idx = block_x / threads_size;
  block_id = block_x % threads_size;
  pos_block_size = threadIdx.x;

  int protein_idx = thread_dim*block_id+pos_block_size;

  if (rna_idx < idx_rna[indiv_id])
    if (protein_idx < rna[indiv_id][rna_idx]->nb_protein) {
      int start_protein_pos = rna[indiv_id][rna_idx]->leading_lagging == 0 ? rna[indiv_id][rna_idx]->start_prot[protein_idx] + 13 : rna[indiv_id][rna_idx]->start_prot[protein_idx] - 13;
      int length = -1;

      if (rna[indiv_id][rna_idx]->leading_lagging == 0) {
        start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                            start_protein_pos - dna_size[indiv_id]
                                                                    : start_protein_pos;

        if (rna[indiv_id][rna_idx]->start_prot[protein_idx] < rna[indiv_id][rna_idx]->end) {
          length = rna[indiv_id][rna_idx]->end - rna[indiv_id][rna_idx]->start_prot[protein_idx];
        } else {
          length = dna_size[indiv_id] - rna[indiv_id][rna_idx]->start_prot[protein_idx] + rna[indiv_id][rna_idx]->end + 1;
        }

        length -= 13;

        /*if (indiv_id == 737 && rna_idx == 47)
          printf("----> START %d %d END %d LENGTH %d Indiv ID %d RNA IDX %d (block.X %d block_x %d "
                     "block.Y %d pos_block %d) Protein IDX %d (block id %d pos block size %d threads size %d)\n",rna[indiv_id][rna_idx]->start_prot[protein_idx],start_protein_pos,
                 rna[indiv_id][rna_idx]->end,length,indiv_id,rna_idx,blockIdx.x,block_x,blockIdx.y,
                 pos_block_size,protein_idx,block_id,pos_block_size,threads_size);*/

      } else {
        start_protein_pos = start_protein_pos < 0 ?
                            dna_size[indiv_id] + start_protein_pos
                                                  : start_protein_pos;

        if (rna[indiv_id][rna_idx]->start_prot[protein_idx] > rna[indiv_id][rna_idx]->end) {
          length = rna[indiv_id][rna_idx]->start_prot[protein_idx] - rna[indiv_id][rna_idx]->end;
        } else {
          length = rna[indiv_id][rna_idx]->start_prot[protein_idx] +  dna_size[indiv_id] - rna[indiv_id][rna_idx]->end;
        }


        length -= 13;
        /*if (indiv_id == 737 && rna_idx == 47)
          printf("----> START %d %d END %d LENGTH %d Indiv ID %d RNA IDX %d (block.X %d block_x %d "
                     "block.Y %d pos_block %d) Protein IDX %d (block id %d pos block size %d threads size %d)\n",rna[indiv_id][rna_idx]->start_prot[protein_idx],start_protein_pos,
                 rna[indiv_id][rna_idx]->end,length,indiv_id,rna_idx,blockIdx.x,block_x,blockIdx.y,
                 pos_block_size,protein_idx,block_id,pos_block_size,threads_size);*/
      }

      bool is_protein = false;

      /*if (indiv_id == debug_iid) {
        printf("LENGTH is %d\n",length);
      }*/
      length+=1;
      length = length - (length%3);
      /*if (indiv_id == debug_iid) {
        printf("LENGTH UPDATED is %d\n",length);
      }*/

      for (int loop_i = 0; length - loop_i >= 2; loop_i+=3) {//start_protein_pos != rna[indiv_id][rna_idx]->end) {
        int t_k;

        /*if (indiv_id == debug_iid)
          printf("Starting loop id %d\n",loop_i);*/

        if (rna[indiv_id][rna_idx]->leading_lagging == 0) {
          start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                              start_protein_pos - dna_size[indiv_id]
                                                                      : start_protein_pos;
          is_protein=false;


          /*if (indiv_id == debug_iid)
            printf("Starting search at %d\n",start_protein_pos);*/

          for (int k = 0; k < 3; k++) {
            t_k = start_protein_pos+k >= dna_size[indiv_id] ?
                  start_protein_pos - dna_size[indiv_id] + k :
                  start_protein_pos + k;

            /*printf("%d-%d-%d :: %lu : %d || %c (%d)\n", // %c (%d) /
                   indiv_id,rna_idx,protein_idx,
                   dna_size[indiv_id],
                   t_k,
            //       dna[indiv_id][t_k],t_k,
                   PROTEIN_END_LEAD[k],k);*/

            if (dna[indiv_id][t_k] == PROTEIN_END_LEAD[k]) {
              //ab=1;
              is_protein=true;
            } else {
              //ab=0;
              is_protein=false;
              break;
            }
          }

          /*if (indiv_id == 737 && rna_idx == 47)
            printf("Protein %d : %d : %d : %d : %d -> %d |||| %d %d --  E %d ERNA %d\n",indiv_id,rna_idx,protein_idx,loop_i,
                   start_protein_pos,is_protein,
                   length - loop_i,length,t_k,rna[indiv_id][rna_idx]->end);*/
          //printf("Protein %d : %d : %d : %d : %d -> %d (%d)\n",indiv_id,rna_idx,protein_idx,loop_i,start_protein_pos,is_protein,ab);

          //cProtein* test = (cProtein*) malloc(20);
          //cRNA* l_rna = (cRNA*)malloc(sizeof(cRNA));
          if (is_protein ) {
            int prot_length = -1;
            if (rna[indiv_id][rna_idx]->start_prot[protein_idx]+13 < t_k) {
              prot_length = t_k - (rna[indiv_id][rna_idx]->start_prot[protein_idx]+13);
            } else {
              prot_length = dna_size[indiv_id] - (rna[indiv_id][rna_idx]->start_prot[protein_idx]+13) + t_k;
            }

            if (prot_length >= 3) {
              //cProtein* l_protein = (cProtein*) malloc(sizeof(cProtein));
              //cProtein* l_protein = (cProtein*)malloc(sizeof(cProtein));
              //cRNA* l_rna = (cRNA*)malloc(sizeof(cRNA));
              int idx = atomicAdd(idx_protein + indiv_id, 1);

              protein_list[indiv_id][idx].protein_start = rna[indiv_id][rna_idx]->start_prot[protein_idx];
              protein_list[indiv_id][idx].protein_end = t_k;
              protein_list[indiv_id][idx].e = rna[indiv_id][rna_idx]->e;
              protein_list[indiv_id][idx].leading_lagging = rna[indiv_id][rna_idx]->leading_lagging;
              protein_list[indiv_id][idx].protein_length = prot_length;

              /*if (indiv_id == 737)
                printf("Address %d-%d-%d -- %d : %p (%d %d) -- (%d %d) -- %d\n",
                       indiv_id, rna_idx, protein_idx, idx, (void*) l_protein,
                       l_protein->protein_start, l_protein->protein_end,
                       rna[indiv_id][rna_idx]->begin,
                       rna[indiv_id][rna_idx]->end, prot_length);*/
              // = l_protein;
            }
            break;
          }



          start_protein_pos+=3;
          start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                start_protein_pos - dna_size[indiv_id]
                                                                       : start_protein_pos;

          //if (indiv_id == debug_iid) printf("New position %d %d LOOPSIZE %d\n",start_protein_pos,loop_i,length-loop_i);
        } else {

          is_protein=false;
          start_protein_pos = start_protein_pos < 0 ?
                              dna_size[indiv_id] + start_protein_pos
                                                    : start_protein_pos;



          for (int k = 0; k < 3; k++) {
            t_k = start_protein_pos-k < 0 ?
                  dna_size[indiv_id] + (start_protein_pos - k) :
                  start_protein_pos - k;

            /*if (indiv_id == 1020) printf("%d-%d-%d :: %lu : %d || %c (%d) // %c (%d)\n", //
                   indiv_id,rna_idx,protein_idx,
                   dna_size[indiv_id],
                   t_k,
                   dna[indiv_id][t_k],t_k,
                  PROTEIN_END_LAG[k],k);*/

            if (dna[indiv_id][t_k] == PROTEIN_END_LAG[k]) {
              //ab=1;
                is_protein=true;
            } else {
              //ab=0;
                is_protein=false;
                break;
            }
          }


          /*if (indiv_id == debug_iid)
           printf("Protein %d : %d : %d : %d : %d -> %d |||| %d %d\n",indiv_id,rna_idx,protein_idx,loop_i,
                  start_protein_pos,is_protein,
                  length - loop_i,length);*/

          //cProtein* test = (cProtein*) malloc(20);
          //cRNA* l_rna = (cRNA*)malloc(sizeof(cRNA));

          if (is_protein) {
            int prot_length=-1;
            if (rna[indiv_id][rna_idx]->start_prot[protein_idx]-13 > t_k) {
              prot_length = (rna[indiv_id][rna_idx]->start_prot[protein_idx]-13) - t_k;
            } else {
              prot_length = (rna[indiv_id][rna_idx]->start_prot[protein_idx]-13) +  dna_size[indiv_id] - t_k;
            }

            if (prot_length >= 3) {
              //cProtein* l_protein = (cProtein*) malloc(sizeof(cProtein));
              int idx = atomicAdd(idx_protein + indiv_id, 1);
              //cProtein* l_protein = (cProtein*) malloc(sizeof(int));
              protein_list[indiv_id][idx].protein_start = rna[indiv_id][rna_idx]->start_prot[protein_idx];
              protein_list[indiv_id][idx].protein_end = t_k;

              /*if (indiv_id == 878) {
                for (int k = 0; k < 3; k++) {
                  t_k = start_protein_pos - k < 0 ?
                        dna_size[indiv_id] - start_protein_pos - k :
                        start_protein_pos - k;
                  printf("[%d -- %c] ",t_k,dna[indiv_id][t_k]);
                }
                printf("\n");
              }*/
              protein_list[indiv_id][idx].protein_length = prot_length;
              protein_list[indiv_id][idx].e = rna[indiv_id][rna_idx]->e;
              protein_list[indiv_id][idx].leading_lagging = rna[indiv_id][rna_idx]->leading_lagging;



              /*if (indiv_id == 737)
                printf("Address %d-%d-%d -- %d : %p (%d %d) -- (%d %d) -- %d %lu\n",
                       indiv_id, rna_idx, protein_idx, idx, (void*) l_protein,
                       l_protein->protein_start, l_protein->protein_end,
                       rna[indiv_id][rna_idx]->begin,
                       rna[indiv_id][rna_idx]->end, prot_length,dna_size[indiv_id]);*/
              //protein_list[indiv_id][idx] = l_protein;
            }
            break;
          }
          //if (is_protein) cProtein* l_protein = (cProtein*) malloc(1*sizeof(cProtein));
          //int* test = (int*)malloc(sizeof(int));


          start_protein_pos = start_protein_pos-3;
          start_protein_pos = start_protein_pos < 0 ?
                                dna_size[indiv_id] + start_protein_pos
                                                      : start_protein_pos;

          //if (indiv_id == debug_iid) printf("New position %d %d\n",start_protein_pos,loop_i);
        }
      }
    }
}

__global__ void max_protein(int32_t* max_nb_protein, int32_t* idx_protein) {
  int indiv_id = blockIdx.x;

  atomicMax(max_nb_protein, idx_protein[indiv_id]-1);
}

__global__ void translate_protein(float w_max, int32_t* idx_protein,
                                  cProtein** protein_list,
                                  char** dna, size_t* dna_size, int threads_size, int thread_dim) {
  int indiv_id = blockIdx.x / threads_size;
  int block_id = blockIdx.x % threads_size;
  int pos_block_size = threadIdx.x;

  int protein_idx = thread_dim*block_id+pos_block_size;

  if (protein_idx < idx_protein[indiv_id]) {
    // Translate RNA to codon
    int c_pos = protein_list[indiv_id][protein_idx].protein_start, t_pos;
    int end_pos = protein_list[indiv_id][protein_idx].protein_end;
    if (protein_list[indiv_id][protein_idx].leading_lagging == 0) {
      c_pos += 13;
      end_pos -=3;

      c_pos = c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
      end_pos = end_pos < 0 ? dna_size[indiv_id] + end_pos : end_pos;
    } else {
      c_pos -= 13;
      end_pos +=3;

      end_pos = end_pos >= dna_size[indiv_id] ? end_pos - dna_size[indiv_id] : end_pos;
      c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
    }

    /*
    if (indiv_id == 410)
      printf("Protein %d translate from %d to %d (%d)\n",protein_idx, c_pos,end_pos,protein_list[indiv_id][protein_idx]->leading_lagging);
    */
    int8_t value = 0;
    int8_t codon_list[64] = {};
    int8_t codon_idx = 0;
    int32_t count_loop = 0;

    bool contin = true;
    if (protein_list[indiv_id][protein_idx].leading_lagging == 0) {
      // LEADING

      while (count_loop<protein_list[indiv_id][protein_idx].protein_length/3 && codon_idx < 64) {
        value = 0;
        for (int8_t i = 0; i < 3; i++) {
          t_pos = c_pos + i >= dna_size[indiv_id] ? c_pos + i - dna_size[indiv_id] : c_pos + i;
          if (dna[indiv_id][t_pos] == '1' ) value += 1 << (CODON_SIZE - i - 1);
        }
        codon_list[codon_idx] = value;
        codon_idx++;

        count_loop++;
        c_pos+=3;
        c_pos = c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
      }
    } else {
      // LAGGING
      while (count_loop<protein_list[indiv_id][protein_idx].protein_length/3 && codon_idx < 64) {
        value = 0;
        for (int8_t i = 0; i < 3; i++) {
          t_pos = c_pos - i < 0 ? dna_size[indiv_id] + (c_pos - i) : c_pos - i;
          if (dna[indiv_id][t_pos] != '1' ) value += 1 << (CODON_SIZE - i - 1);
        }
        codon_list[codon_idx] = value;
        codon_idx++;

        count_loop++;

        c_pos-=3;
        c_pos = c_pos < 0 ? c_pos + dna_size[indiv_id] : c_pos;
      }
    }

    double M = 0.0;
    double W = 0.0;
    double H = 0.0;

    int32_t nb_m = 0;
    int32_t nb_w = 0;
    int32_t nb_h = 0;

    bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
    bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
    bool bin_h = false;


    for (int i = 0; i < codon_idx; i++) {
      switch (codon_list[i])
      {
        case CODON_M0 :
        {
          // M codon found
          nb_m++;

          // Convert Gray code to "standard" binary code
          bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ M <<= 1;
          M *= 2;

          // Add this nucleotide's contribution to M
          if (bin_m) M += 1;

          break;
        }
        case CODON_M1 :
        {
          // M codon found
          nb_m++;

          // Convert Gray code to "standard" binary code
          bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
          //~ M <<= 1;
          M *= 2;

          // Add this nucleotide's contribution to M
          if (bin_m) M += 1;

          break;
        }
        case CODON_W0 :
        {
          // W codon found
          nb_w++;

          // Convert Gray code to "standard" binary code
          bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ W <<= 1;
          W *= 2;

          // Add this nucleotide's contribution to W
          if (bin_w) W += 1;

          break;
        }
        case CODON_W1 :
        {
          // W codon found
          nb_w++;

          // Convert Gray code to "standard" binary code
          bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ W <<= 1;
          W *= 2;

          // Add this nucleotide's contribution to W
          if (bin_w) W += 1;

          break;
        }
        case CODON_H0 :
        case CODON_START : // Start codon codes for the same amino-acid as H0 codon
        {
          // H codon found
          nb_h++;

          // Convert Gray code to "standard" binary code
          bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ H <<= 1;
          H *= 2;

          // Add this nucleotide's contribution to H
          if (bin_h) H += 1;

          break;
        }
        case CODON_H1 :
        {
          // H codon found
          nb_h++;

          // Convert Gray code to "standard" binary code
          bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ H <<= 1;
          H *= 2;

          // Add this nucleotide's contribution to H
          if (bin_h) H += 1;

          break;
        }
      }
    }



    //  ----------------------------------------------------------------------------------
    //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
    //  ----------------------------------------------------------------------------------
    protein_list[indiv_id][protein_idx].m = nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
    protein_list[indiv_id][protein_idx].w = nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
    protein_list[indiv_id][protein_idx].h = nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

    //  ------------------------------------------------------------------------------------
    //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
    //  ------------------------------------------------------------------------------------
    // x_min <= M <= x_max
    // w_min <= W <= w_max
    // h_min <= H <= h_max
    protein_list[indiv_id][protein_idx].m  = (X_MAX - X_MIN) * protein_list[indiv_id][protein_idx].m + X_MIN;
    protein_list[indiv_id][protein_idx].w  = (w_max - W_MIN) * protein_list[indiv_id][protein_idx].w + W_MIN;
    protein_list[indiv_id][protein_idx].h  = (H_MAX - H_MIN) * protein_list[indiv_id][protein_idx].h + H_MIN;

    if ( nb_m == 0 || nb_w == 0 || nb_h == 0 || protein_list[indiv_id][protein_idx].w == 0.0 ||
        protein_list[indiv_id][protein_idx].h == 0.0 )
    {
      protein_list[indiv_id][protein_idx].is_functional = false;
    }
    else
    {
      protein_list[indiv_id][protein_idx].is_functional = true;
    }
  }
}


__global__ void compute_phenotype(int32_t* idx_protein, cProtein** protein_list,
                                  float** phenotype, int threads_size, int thread_dim) {
  int indiv_id = blockIdx.x / threads_size;
  int block_id = blockIdx.x % threads_size;
  int pos_block_size = threadIdx.x;

  int protein_idx = thread_dim*block_id+pos_block_size;

  if (protein_idx < idx_protein[indiv_id]) {
    if ( fabs(protein_list[indiv_id][protein_idx].w) < 1e-15 ||
        fabs(protein_list[indiv_id][protein_idx].h) < 1e-15 ) return;

    if (protein_list[indiv_id][protein_idx].is_functional) {

      // Compute triangle points' coordinates
      float x0 = protein_list[indiv_id][protein_idx].m -
                 protein_list[indiv_id][protein_idx].w;
      float x1 = protein_list[indiv_id][protein_idx].m;
      float x2 = protein_list[indiv_id][protein_idx].m +
                 protein_list[indiv_id][protein_idx].w;

      /*if (indiv_id == 991)
        printf("Protein %d : %f %f %f\n",protein_idx,protein_list[indiv_id][protein_idx]->m,
               protein_list[indiv_id][protein_idx]->w,protein_list[indiv_id][protein_idx]->h);*/

      int ix0 = (int) (x0 * 300);
      int ix1 = (int) (x1 * 300);
      int ix2 = (int) (x2 * 300);

      if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
      if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
      if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

      // Compute the first equation of the triangle
      float incY = (protein_list[indiv_id][protein_idx].h *
                    protein_list[indiv_id][protein_idx].e) / (ix1 - ix0);
      int count = 1;
      // Updating value between x0 and x1

      for (int i = ix0 + 1; i < ix1; i++) {
        atomicAdd(&(phenotype[indiv_id][i]), incY * (count++));
      }

      atomicAdd(&phenotype[indiv_id][ix1],
                (protein_list[indiv_id][protein_idx].h *
                 protein_list[indiv_id][protein_idx].e));

      // Compute the second equation of the triangle
      incY = (protein_list[indiv_id][protein_idx].h *
              protein_list[indiv_id][protein_idx].e) / (ix2 - ix1);
      count = 1;

      // Updating value between x1 and x2
      for (int i = ix1 + 1; i < ix2; i++) {
        atomicAdd(&phenotype[indiv_id][i],
                  ((protein_list[indiv_id][protein_idx].h *
                    protein_list[indiv_id][protein_idx].e) -
                   (incY * (count++))));
      }
    }
  }
}

__global__ void compute_metaerror_fitness(double selection_pressure,float** phenotype,
                                          float* target,
                                          float* metaerror, double* fitness) {
  int indiv_id = blockIdx.x;

  int fuzzy_idx = threadIdx.x;

  __shared__ float delta[300];

  if (phenotype[indiv_id][fuzzy_idx] > 1) phenotype[indiv_id][fuzzy_idx] = 1;
  if (phenotype[indiv_id][fuzzy_idx] < 0) phenotype[indiv_id][fuzzy_idx] = 0;

  delta[fuzzy_idx] = phenotype[indiv_id][fuzzy_idx] - target[fuzzy_idx];

  //if (indiv_id == 15) printf("DELTA[%d] = %f\n",fuzzy_idx,delta[fuzzy_idx]);

 /* if (threadIdx.x == 0) {
    metaerror[indiv_id] = 0;
  }*/

  __syncthreads();

  /*if (threadIdx.x < 299) {
    atomicAdd(metaerror+indiv_id,
              ((fabs(delta[fuzzy_idx]) +
                fabs(delta[fuzzy_idx + 1])) / (600.0)));
  }

  __syncthreads();*/

  if (threadIdx.x == 0) {
    metaerror[indiv_id] = 0;

    for (int i = 0; i < 299; i++) {
      metaerror[indiv_id] +=
                ((fabs(delta[i]) +
                  fabs(delta[i + 1])) / (600.0));
    }

    fitness[indiv_id] = exp(
        -selection_pressure * ((double)metaerror[indiv_id]));
  }
}

__global__ void free_list(cProtein** protein_list,
                          cRNA*** rna, int32_t* idx_protein,int32_t* idx_rna) {

  int indiv_id = blockIdx.x;

}

__global__ void debug_dna(size_t* dna_size, char** dna) {
  for (int i = 0; i < 1024; i++) {
    printf("DNA SIZE %d : %lu\n",i,dna_size[i]);

    for (size_t pos = 0; pos < dna_size[i]; pos++)
      dna[i][pos];
  }
}

__global__ void debug_promoter_start(size_t* dna_size,
                                     pStruct** dynPromoterList,
                                     int* nb_promoters, int indiv_id) {
  //int indiv_id = blockIdx.x;

 // for (int i = 0; i < 1024; i++)
 //   printf("DNA SIZE %d : %lu\n",i,dna_size[i]);
  printf("RNA %d : %d\n",indiv_id,nb_promoters[indiv_id]);

  printf("Individual %d (GPU) Promoters : LEADING ",indiv_id);
  // LEADING
  for (int idx = 0; idx < nb_promoters[indiv_id]; idx++) {
    if (dynPromoterList[indiv_id][idx].leading_or_lagging)
      printf("%d ",dynPromoterList[indiv_id][idx].pos);
  }
  printf("\n");

  printf("Individual %d (GPU) Promoters : LAGGING ",indiv_id);
  // LAGGING
  for (int idx = 0; idx < nb_promoters[indiv_id]; idx++) {
    if (!dynPromoterList[indiv_id][idx].leading_or_lagging)
      printf("%d ",dynPromoterList[indiv_id][idx].pos);
  }
  printf("END\n");
}


__global__ void debug_promoter_stop(size_t* dna_size,
                                     int8_t** dna_lead_term,
                                     int8_t** dna_lag_term, int* nb_promoters, int indiv_id) {
  //int indiv_id = blockIdx.x;

  // for (int i = 0; i < 1024; i++)
  //   printf("DNA SIZE %d : %lu\n",i,dna_size[i]);
  //printf("RNA %d : %d\n",indiv_id,nb_promoters[indiv_id]);

  printf("Term from CPU\n");
  printf("Individual %d (GPU) STOPs : LEADING ",indiv_id);
  // LEADING
  for (int pos = 0; pos < dna_size[indiv_id]; pos++) {
    if (dna_lead_term[indiv_id][pos] > 0) {
      printf("%d ",pos);
    }
  }
  printf("\n");

  printf("Individual %d (GPU) STOPs : LAGGING ",indiv_id);
  // LAGGING
  for (int pos = 0; pos < dna_size[indiv_id]; pos++) {
    if (dna_lag_term[indiv_id][pos] > 0) {
      printf("%d ",pos);
    }
  }
  printf("\n");
}

__global__ void debug_rna(size_t* dna_size,
                                     int8_t** dna_lead_term,
                                     int8_t** dna_lag_term,
                          cRNA*** rna,int32_t* idx_rna,
                          int indiv_id) {
  //int indiv_id = blockIdx.x;

  printf("Individual %d (GPU) %lu : \n",indiv_id,dna_size[indiv_id]);
  // LEADING
  for (int rna_idx = 0; rna_idx < idx_rna[indiv_id]; rna_idx++) {
    printf("RNA %d : ",rna_idx);
    if (rna[indiv_id][rna_idx]->leading_lagging == 0)
      printf("LEADING ");
    else
      printf("LAGGING ");

    printf("%d %d %f\n",rna[indiv_id][rna_idx]->begin,rna[indiv_id][rna_idx]->end,rna[indiv_id][rna_idx]->e);
  }
  printf("\n");
}

__global__ void debug_protein(int32_t* idx_protein,
                          cProtein** protein_list, char** dna,
                          int indiv_id) {
  printf("Individual %d (GPU) -- %d : \n",indiv_id,idx_protein[indiv_id]);

  for (int prot_idx = 0; prot_idx < idx_protein[indiv_id]; prot_idx++) {
    printf("Protein %d : %d %d %lf %lf %lf\n",prot_idx,
           protein_list[indiv_id][prot_idx].protein_start,
           protein_list[indiv_id][prot_idx].protein_end,
           protein_list[indiv_id][prot_idx].m,
           protein_list[indiv_id][prot_idx].h,
           protein_list[indiv_id][prot_idx].w);
  }

  printf("\n");

}

__global__ void debug_phenotype(float** phenotype,float* target, float* metaerror,
                                double* fitness,
                              int indiv_id) {
  printf("Individual %d (GPU) : \n", indiv_id);

  double delta_1[300];
  double metaerror_1 = 0;
  double fitness_1 = 0;

  for (int i = 0; i < 300; i++) {
    delta_1[i] = (double)phenotype[indiv_id][i] - (double)target[i];
  }

  for (int i = 0; i < 299; i++) {
    metaerror_1 +=
        (double)((fabs((double)delta_1[i]) +
          fabs((double)delta_1[i + 1])) / (double)(600.0));
  }

    //if (phenotype[indiv_id][i] != 0) {

  for (int i = 0; i < 300; i++) {

      if (phenotype[indiv_id][i] != 0) printf("[%d : %f]\n", i, phenotype[indiv_id][i]);


  }

  fitness_1 = exp(
      -1000.0 * (double)metaerror_1);

  printf("METAERROR %f --  %f // %e -- %e\n",metaerror_1,metaerror[indiv_id],
         fitness[indiv_id],fitness_1);
}

__global__ void debug_fitness(float** phenotype,float* target,
                              float* metaerror, double* fitness,
                                int indiv_id) {

  printf("Individual %d (GPU) : \n", indiv_id);


}
