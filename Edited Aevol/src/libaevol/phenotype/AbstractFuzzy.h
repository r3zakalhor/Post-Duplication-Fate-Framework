//
// Created by arrouan on 30/07/15.
//

#ifndef AEVOL_ABSTRACTFUZZY_H
#define AEVOL_ABSTRACTFUZZY_H


#include <list>

#include "macros.h"
#include "Point.h"

namespace aevol {

class AbstractFuzzy
{
 public:
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
//  AbstractFuzzy() {};
//  AbstractFuzzy(const AbstractFuzzy& f) {};
//  AbstractFuzzy(const gzFile backup) {};

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~AbstractFuzzy() {};

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void save(gzFile backup) const = 0;
  virtual void load(gzFile backup) = 0;
  virtual void reset() = 0;
  virtual void simplify() = 0;
  virtual void add_triangle(ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose = false)  = 0;
  virtual void add(const AbstractFuzzy& f)  = 0;
  virtual void sub(const AbstractFuzzy& f, bool verbose = false) = 0;
  virtual void add_point(ProteinConcentration x, ProteinConcentration y) = 0;

  /// `clipping_direction` is only used for `clip` function's keyword.
  enum clipping_direction: bool {min, max};
  virtual void clip(clipping_direction direction, ProteinConcentration bound) = 0;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  virtual ProteinConcentration get_geometric_area(bool verbose = false) const = 0;
  virtual ProteinConcentration get_geometric_area(ProteinConcentration start_segment, ProteinConcentration end_segment) const = 0;

  virtual bool is_identical_to(const AbstractFuzzy& fs, ProteinConcentration tolerance) const = 0;

  virtual void print() const = 0;

  virtual void clear() = 0;

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

};

} // namespace aevol


#endif //AEVOL_ABSTRACTFUZZY_H
