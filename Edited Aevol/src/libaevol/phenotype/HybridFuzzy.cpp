//
// Created by arrouan on 31/07/15.
//


//#include <math.h>

#ifdef __BLAS__
#include <cblas.h>
#endif

#include <iostream>
#include "HybridFuzzy.h"

namespace aevol {


HybridFuzzy::HybridFuzzy( )
{
  _pheno_size = PHENO_SIZE;
  _points = new ProteinConcentration[_pheno_size];
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;
}

HybridFuzzy::HybridFuzzy( const HybridFuzzy& model )
{
  _pheno_size = PHENO_SIZE;
  _points = new ProteinConcentration[_pheno_size];

#ifdef __BLAS__
#ifdef __FLOAT_CONCENTRATION
  cblas_scopy(_pheno_size,model._points,1,_points,1);
#else
  cblas_dcopy(_pheno_size,model._points,1,_points,1);
#endif
#else
	for (int i=0; i < _pheno_size; i++)
		_points[i] = model._points[i];
#endif
}

HybridFuzzy::HybridFuzzy( gzFile backup_file )
{
  _pheno_size = PHENO_SIZE;
  _points = new ProteinConcentration[_pheno_size];
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;

  load( backup_file );
}

// =================================================================
//                             Destructors
// =================================================================

HybridFuzzy::~HybridFuzzy( void )
{
  if (_points != NULL) delete [] _points;
  _points = NULL;
}

// =================================================================
//                            Public Methods
// =================================================================
void HybridFuzzy::simplify( )
{

}

void HybridFuzzy::reset( )
{
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;
}

void HybridFuzzy::add_triangle( ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose )
{
  if ( fabs(width) < 1e-15 || fabs(height) < 1e-15 ) return;

  // Compute triangle points' coordinates
  ProteinConcentration x0 = mean - width;
  ProteinConcentration x1 = mean;
  ProteinConcentration x2 = mean + width;

  // Compute the first equation of the triangle
  // Updating value between x0 and x1
  //printf("Update between %d (%f) and %d (%f) \n",(int) std::ceil(x0*299.0), x0, (int) std::ceil(x1*299.0),x1);

    int loop_A_start = (int) std::ceil(x0 * 299.0);
    loop_A_start = loop_A_start < 0 ? 0 : loop_A_start;
    loop_A_start = loop_A_start > 299 ? 299 : loop_A_start;

    int loop_A_end = (int) std::ceil(x1 * 299.0);
    loop_A_end = loop_A_end < 0 ? 0 : loop_A_end;
    loop_A_end = loop_A_end > 299 ? 299 : loop_A_end;

    for (int i = loop_A_start; i < loop_A_end; i++) {
                  _points[i] += (((i / 299.0) - x0) / (x1 - x0)) * height;
                  //printf("Update point %d : %f (%f %f %f)\n", i, (((i / 299.0) - x0) / (x1 - x0)) * height,
                  //       ((i / 299.0) - x0), (x1 - x0), height);
              }


  // Compute the second equation of the triangle
  // Updating value between x1 and x2
          //printf("Update between %d (%f) and %d (%f)\n",(int) std::ceil(x1*299.0),x1, (int) std::ceil(x2*299.0),x2);

    int loop_B_start = (int) std::ceil(x1 * 299.0);
    loop_B_start = loop_B_start < 0 ? 0 : loop_B_start;
    loop_B_start = loop_B_start > 299 ? 299 : loop_B_start;

    int loop_B_end = (int) std::ceil(x2 * 299.0);
    if (loop_B_end > 299) _points[299] += height * ((x2 - 1.0) / (x2 - x1));

    loop_B_end = loop_B_end < 0 ? 0 : loop_B_end;
    loop_B_end = loop_B_end > 299 ? 299 : loop_B_end;

              for (int i = loop_B_start; i < loop_B_end; i++) {
                  _points[i] += height * ((x2 - (i / 299.0)) / (x2 - x1));
                  //printf("Update point %d : %f (%f %f %f)\n", i, height * ((x2 - (i / 299.0)) / (x2 - x1)),
                  //       (x2 - (i / 299.0)), (x2 - x1), height);
              }

}

void HybridFuzzy::add( const AbstractFuzzy& f )
{
  const HybridFuzzy to_add = (HybridFuzzy&)(f);
#ifdef __BLAS__
#ifdef __FLOAT_CONCENTRATION
  cblas_saxpy(_pheno_size, 1.0, to_add.points(), 1, _points, 1);
#else
  cblas_daxpy(_pheno_size, 1.0, to_add.points(), 1, _points, 1);
#endif
#else
		for (int i = 0; i < _pheno_size; i++) {
			if (to_add._points[i] != 0) _points[i] = _points[i] + to_add._points[i];
		}
#endif
}

void HybridFuzzy::sub( const AbstractFuzzy& f, bool verbose)
{
  const HybridFuzzy to_sub = (HybridFuzzy&)(f);
#ifdef __BLAS__
#ifdef __FLOAT_CONCENTRATION
  cblas_saxpy(_pheno_size, -1.0, to_sub.points(), 1, _points, 1);
#else
  cblas_daxpy(_pheno_size, -1.0, to_sub.points(), 1, _points, 1);
#endif
#else
		for (int i = 0; i < _pheno_size; i++) {
			if (to_sub._points[i] !=0 ) _points[i] = _points[i] - to_sub._points[i];
		}
#endif
}

ProteinConcentration HybridFuzzy::get_geometric_area(bool verbose) const
{
  return get_geometric_area(X_MIN,X_MAX);
}

ProteinConcentration HybridFuzzy::get_geometric_area( ProteinConcentration start_segment, ProteinConcentration end_segment ) const
{
  ProteinConcentration area = 0;

  int istart_segment = (int) (start_segment  * _pheno_size);
  int iend_segment = (int) (end_segment  * _pheno_size);

  if (istart_segment < 0) istart_segment = 0; else if (istart_segment > (_pheno_size-1)) istart_segment = _pheno_size-1;
  if (iend_segment < 0) iend_segment = 0; else if (iend_segment > (_pheno_size-1)) iend_segment = _pheno_size-1;
  for (int i = istart_segment; i < iend_segment; i++) {

      printf("TA [ %lf %lf ] [ %lf %lf ] = %lf\n",((double)i)/300.0,((double)i+1)/300.0,_points[i],_points[i+1],
             ((fabs(_points[i]) + fabs(_points[i+1])) / (2.0*_pheno_size)));
    area+=((fabs(_points[i]) + fabs(_points[i+1])) / (2.0*_pheno_size));
  }

  return area;
}



bool HybridFuzzy::is_identical_to( const AbstractFuzzy& f, ProteinConcentration tolerance  ) const
{
  const HybridFuzzy fs = (HybridFuzzy&)(f);
  // Since list::size() has constant complexity since C++ 11, checking
  // size is an inexpensive first step.
  if (get_pheno_size() != fs.get_pheno_size())
    return false;

  for (int i = 0; i < _pheno_size; i++)
    if (fabs(_points[i] - fs.points()[i]) > tolerance * (fabs(_points[i]) + fabs(fs.points()[i])) or
        fabs(_points[i] - fs.points()[i]) > tolerance * (fabs(_points[i]) + fabs(fs.points()[i])))
      return false;
  return true;
}


void HybridFuzzy::save( gzFile backup_file ) const
{
  gzwrite(backup_file, &_pheno_size, sizeof(_pheno_size));
  std::cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << std::endl;
  std::cout << __FILE__ << ":" << __LINE__ << ":" << _pheno_size << std::endl;

  for (int i = 0; i < _pheno_size; i++)
    gzwrite(backup_file, &_points[i], sizeof(_points[i]));
}


void HybridFuzzy::load( gzFile backup_file ) {
  gzread(backup_file, &_pheno_size, sizeof(_pheno_size));
  std::cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) <<
                                                     std::endl;
  std::cout << __FILE__ << ":" << __LINE__ << ":" << _pheno_size << std::endl;
  for (int i = 0; i < _pheno_size; i++) {
    gzread(backup_file, &_points[i], sizeof(_points[i]));
  }
}


void HybridFuzzy::clip(clipping_direction direction, ProteinConcentration bound) {

  if (direction == clipping_direction::min)
    for (int i = 0; i < _pheno_size; i++)
      _points[i] = _points[i] < bound ? bound : _points[i];
  else if (direction == clipping_direction::max)
    for (int i = 0; i < _pheno_size; i++)
      _points[i] = _points[i] > bound ? bound : _points[i];
}

void HybridFuzzy::add_point(ProteinConcentration x, ProteinConcentration y) {
  int ix = (int) ( x * _pheno_size);
  if (ix < _pheno_size)
    _points[ix] = y;
}

// =================================================================
//                           Protected Methods
// =================================================================
ProteinConcentration HybridFuzzy::get_y( ProteinConcentration x ) const
{
  int ix = (int) ( x * _pheno_size);

  ProteinConcentration retValue = _points[ix];

  return retValue;
}

void HybridFuzzy::print() const
{
  for (int i = 0; i < _pheno_size; i++)
    if (_points[i]!=0) printf("[%d : %f]\n",i,_points[i]);
  printf("\n");
}

bool HybridFuzzy::compare(Fuzzy* fuzz) {
    bool is_diff = false;
  for (int i = 0; i < 300; i++) {
      double hf = roundf(_points[i] * 10000);
      double vf = roundf(fuzz->y(i/299.0) * 10000);
    if (hf != vf) {
      printf("FUZ[%d] (%f) -> HF %f VF %f\n",i,i/299.0,_points[i],fuzz->y(i/299.0));
      is_diff = true;
    }
  }
  return is_diff;
}
}
