//
// Created by arrouan on 31/07/15.
//

#ifndef AEVOL_HYBRIDFUZZY_H
#define AEVOL_HYBRIDFUZZY_H

#define PHENO_SIZE 300

#include <vector>

#include "macros.h"
#include "Point.h"
#include "AbstractFuzzy.h"
#include "Fuzzy.h"

namespace aevol {

class HybridFuzzy : public AbstractFuzzy
{
 public:
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  HybridFuzzy();
  HybridFuzzy(const HybridFuzzy& f);
  HybridFuzzy(const gzFile backup);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~HybridFuzzy();

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void save(gzFile backup) const;
  void load(gzFile backup);
  void reset();
  void simplify();
  void add_triangle(ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose = false);
  void add(const AbstractFuzzy& f);
  void sub(const AbstractFuzzy& f,  bool verbose = false);
  void add_point(ProteinConcentration x, ProteinConcentration y);

  void clip(clipping_direction direction, ProteinConcentration bound);

  bool compare(Fuzzy* fuzz);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  ProteinConcentration* points() const { return _points; };

  ProteinConcentration get_geometric_area(bool verbose = false) const;

  ProteinConcentration get_geometric_area(ProteinConcentration start_segment, ProteinConcentration end_segment) const;

  ProteinConcentration get_y(ProteinConcentration x) const;
  // get_x should be moved out of fuzzy class as it really applies to pair of points

  bool is_identical_to(const AbstractFuzzy& fs, ProteinConcentration tolerance) const;

  int get_pheno_size() const { return _pheno_size; };

  void print() const;
  inline void clear()  {reset();};
  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

 protected:
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================



  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  ProteinConcentration* _points = NULL;
  int _pheno_size = PHENO_SIZE+1;
};

ProteinConcentration trapezoid_area(const Point& p1, const Point& p2);
} // namespace aevol


#endif //AEVOL_HYBRIDFUZZY_H
