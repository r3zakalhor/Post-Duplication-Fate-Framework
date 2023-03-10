//
// Created by arrouan on 01/09/16.
//

#include <list>
#include "algorithm_cuda.h"

#include "raevol/Protein_R.h"
#include "raevol/Rna_R.h"

namespace aevol {
//#ifdef __OPENMP_GPU
std::list<Protein_R>::iterator algorithm_cuda::find_if_protein(
    std::list<Protein_R>::iterator first, std::list<Protein_R>::iterator last,
    int32_t shine_dal_pos) {
  for (; first != last; ++first) {
    if ((*first).shine_dal_pos() == shine_dal_pos) {
      return first;
    }
  }
  return last;
}

std::list<Rna_R>::iterator algorithm_cuda::find_if_rna_1(
    std::list<Rna_R>::iterator first, std::list<Rna_R>::iterator last,
    int32_t from_pos, aevol::Strand strand) {
  for (; first != last; ++first) {
    if (strand == LEADING) {
      if ((*first).promoter_pos() >= from_pos)
        return first;
    }
    else {
      if ((*first).promoter_pos() < from_pos)
        return first;
    }
  }
  return last;
}

std::list<Rna_R>::iterator algorithm_cuda::find_if_rna_2(
    std::list<Rna_R>::iterator first,
    std::list<Rna_R>::iterator last,
    int32_t pos) {

  for (; first != last; ++first) {
      if ((*first).promoter_pos() < pos)
        return first;
  }
  return last;
}

std::list<Rna_R>::iterator algorithm_cuda::find_if_rna_3(
    std::list<Rna_R>::iterator first,
    std::list<Rna_R>::iterator last,
    int32_t pos) {

  for (; first != last; ++first) {
    if ((*first).promoter_pos() >= pos)
      return first;
  }
  return last;
}

std::list<Rna_R>::reverse_iterator algorithm_cuda::find_if_rna_3(
    std::list<Rna_R>::reverse_iterator first,
    std::list<Rna_R>::reverse_iterator last,
    int32_t pos) {

  for (; first != last; ++first) {
    if ((*first).promoter_pos() >= pos)
      return first;
  }
  return last;
}

std::list<Rna_R>::iterator algorithm_cuda::find_if_rna_4(
    std::list<Rna_R>::iterator first,
    std::list<Rna_R>::iterator last,
    int32_t pos) {

  for (; first != last; ++first) {
    if ((*first).promoter_pos() <= pos)
      return first;
  }
  return last;
}

std::list<Point>::const_iterator algorithm_cuda::find_if_point_1(
    std::list<Point>::const_iterator first,
    std::list<Point>::const_iterator last,
    ProteinConcentration x) {

  for (; first != last; ++first) {
    if ((*first).x >= x)
      return first;
  }
  return last;
}


std::list<Point>::iterator algorithm_cuda::find_if_point_2(
    std::list<Point>::iterator first,
    std::list<Point>::iterator last,
    list<Point>::iterator p) {

  for (; first != last; ++first) {
    if ((*first).x != (*p).x)
      return first;
  }
  return last;
}

std::list<Point>::iterator algorithm_cuda::find_if_point_3(
    std::list<Point>::iterator first,
    std::list<Point>::iterator last,
    list<Point>::iterator p) {

  for (; first != last; ++first) {
    if ((*first).y != (*p).y)
      return first;
  }
  return last;
}

std::list<Point>::const_iterator algorithm_cuda::find_if_point_4(
    std::list<Point>::const_iterator first,
    std::list<Point>::const_iterator last,
    ProteinConcentration x) {

  for (; first != last; ++first) {
    if ((*first).x > x)
      return first;
  }
  return last;
}

std::list<Point>::iterator algorithm_cuda::find_if_point_5(
    std::list<Point>::iterator first,
    std::list<Point>::iterator last,
    ProteinConcentration x) {

  for (; first != last; ++first) {
    if ((*first).x > x)
      return first;
  }
  return last;
}
}
//#endif
