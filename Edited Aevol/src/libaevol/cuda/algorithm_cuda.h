//
// Created by arrouan on 01/09/16.
//

//#ifdef __OPENMP_GPU
#ifndef RAEVOL_CUDA_ALGORITHM_CUDA_H
#define RAEVOL_CUDA_ALGORITHM_CUDA_H

#include <vector>
#include <list>
#include  <functional>

#include "Protein.h"
#include "raevol/Protein_R.h"
#include "raevol/Rna_R.h"
#include "ae_enums.h"
#include "Point.h"

using std::vector;
using std::list;

namespace aevol {
  class algorithm_cuda {
   public:
      static std::list<Protein_R>::iterator find_if_protein(
          std::list<Protein_R>::iterator first,
          std::list<Protein_R>::iterator last,
          int32_t shine_dal_pos);

      static std::list<Rna_R>::iterator find_if_rna_1(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos, aevol::Strand);

      static std::list<Rna_R>::iterator find_if_rna_2(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);

      static std::list<Rna_R>::iterator find_if_rna_3(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);

      static std::list<Rna_R>::reverse_iterator find_if_rna_3(
          std::list<Rna_R>::reverse_iterator first,
          std::list<Rna_R>::reverse_iterator last,
          int32_t pos);

      static std::list<Rna_R>::iterator find_if_rna_4(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);

      static std::list<Point>::const_iterator find_if_point_1(
          std::list<Point>::const_iterator first,
          std::list<Point>::const_iterator last,
          ProteinConcentration x);

      static std::list<Point>::iterator find_if_point_2(
          std::list<Point>::iterator first,
          std::list<Point>::iterator last,
          list<Point>::iterator p);

      static std::list<Point>::iterator find_if_point_3(
          std::list<Point>::iterator first,
          std::list<Point>::iterator last,
          list<Point>::iterator p);

      static std::list<Point>::const_iterator find_if_point_4(
          std::list<Point>::const_iterator first,
          std::list<Point>::const_iterator last,
          ProteinConcentration x);

      static std::list<Point>::iterator find_if_point_5(
          std::list<Point>::iterator first,
          std::list<Point>::iterator last,
          ProteinConcentration x);
  };
}


#endif //RAEVOL_CUDA_ALGORITHM_CUDA_H
//#endif
