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
#include <string.h>
#include <assert.h>



// =================================================================
//                            Project Files
// =================================================================
#include "MersenneTwister.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                               Class MersenneTwister                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
MersenneTwister::MersenneTwister(const uint32_t& oneSeed)
{
  seed(oneSeed);
}

MersenneTwister::MersenneTwister(gzFile backup_file)
{
  uint32_t loadArray[SAVE];
  gzread(backup_file, loadArray, SAVE * sizeof(loadArray[0]));
  load(loadArray);
}

MersenneTwister::MersenneTwister(const MersenneTwister & model)
{
  memcpy(state, model.state, N * sizeof(*state));
  pNext = state + (model.pNext - model.state);
  left = model.left;
  assert((pNext - state) / sizeof(*state) + left == N);
}

// =================================================================
//                             Destructors
// =================================================================
MersenneTwister::~MersenneTwister()
{
}

// =================================================================
//                            Public Methods
// =================================================================
int32_t MersenneTwister::binomial_random(int32_t nb_drawings, double prob)
{
  // Returns an integer value that is a random deviate drawn from
  // a binomial distribution of <nb> trials each of probability <prob>.

  int32_t nb_success;

  // The binomial distribution is invariant under changing
  // ProbSuccess to 1-ProbSuccess, if we also change the answer to
  // NbTrials minus itself; we ll remember to do this below. // TODO : is it done?
  double p;
  if (prob <= 0.5) p = prob;
  else p = 1.0 - prob;

  // mean of the deviate to be produced
  double mean = nb_drawings * p;


  if (nb_drawings < 25)
  {
    // Use the direct method while NbTrials is not too large.
    // This can require up to 25 calls to the uniform random.
    nb_success = 0;
    for (int32_t j = 1 ; j <= nb_drawings ; j++)
    {
      if (random() < p) nb_success++;
    }
  }
  else if (mean < 1.0)
  {
    // If fewer than one event is expected out of 25 or more trials,
    // then the distribution is quite accurately Poisson. Use direct Poisson method.

    double g = exp(-mean);
    double t = 1.0;
    int32_t j;
    for (j = 0; j <= nb_drawings ; j++)
    {
      t = t * random();
      if (t < g) break;
    }

    if (j <= nb_drawings) nb_success = j;
    else nb_success = nb_drawings;
  }

  else
  {
    // Use the rejection method.

    double en     = nb_drawings;
    double oldg   = gammln(en + 1.0);
    double pc     = 1.0 - p;
    double plog   = log(p);
    double pclog  = log(pc);

    // rejection method with a Lorentzian comparison function.
    double sq = sqrt(2.0 * mean * pc);
    double angle, y, em, t;
    do
    {
      do
      {
        angle = M_PI * random();
        y = tan(angle);
        em = sq*y + mean;
      } while (em < 0.0 || em >= (en + 1.0)); // Reject.

      em = floor(em); // Trick for integer-valued distribution.
      t = 1.2 * sq * (1.0 + y*y)
              * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);

    } while (random() > t); // Reject. This happens about 1.5 times per deviate, on average.

    nb_success = (int32_t) rint(em);
  }


  // Undo the symmetry transformation.
  if (p != prob) nb_success = nb_drawings - nb_success;

  return nb_success;
}

double MersenneTwister::gaussian_random()
{
  double x1, x2;
  double r = 0;
  do
  {
    x1 = 2.0 * random() - 1.0;
    x2 = 2.0 * random() - 1.0;

    r = x1*x1 + x2*x2; // (x1,x2) must be in the unit circle
  }
  while ((r >= 1.0) || (r == 0));

  r = sqrt((-2.0 * log(r)) / r); // Box-muller transformation

  return x1 * r;
}

void MersenneTwister::multinomial_drawing(int32_t* destination, double* source, int32_t nb_drawings, int32_t nb_colors)
{
  //    This function generates a vector of random variates, each with the
  //    binomial distribution (source code from http://www.agner.org/random/).
  //    The multinomial distribution is the distribution you get when drawing
  //    balls from an urn with more than two colors, with replacement.

  //    Parameters:
  //    destination:    An output array to receive the number of balls of each
  //                    color. Must have space for at least 'nb_colors' elements.
  //    source:         An input array containing the probability or fraction
  //                    of each color in the urn. Must have 'nb_colors' elements.
  //                    All elements must be non-negative. The sum doesn't have
  //                    to be 1, but the sum must be positive.
  //    nb_drawings:    The number of balls drawn from the urn.
  //    nb_colors:      The number of possible colors.


  if (nb_drawings < 0 || nb_colors < 0)
  {
    printf("Negative parameter in multinomial function\n");
    exit(EXIT_FAILURE);
  }
  if (nb_colors == 0) return;

  // compute sum of probabilities
  double sum = 0.0;
  double p;
  for (int32_t i = 0 ; i < nb_colors ; i++)
  {
    p = source[i];
    if (p < 0)
    {
      printf("Negative parameter in multinomial function\n");
      exit(EXIT_FAILURE);
    }
    sum += p;
  }
  if (sum == 0 && nb_drawings > 0)
  {
    printf("Zero sum in multinomial function\n");
    exit(EXIT_FAILURE);
  }

  int32_t x;
  int32_t n = nb_drawings;
  for (int32_t i = 0 ; i < nb_colors - 1 ; i++)
  {
    // generate output by calling binomial (nb_colors-1) times
    p = source[i];
    if (sum <= p)
    {
      // this fixes two problems:
      // 1. prevent division by 0 when sum = 0
      // 2. prevent p/sum getting bigger than 1 in case of rounding errors
      x = n;
    }
    else
    {
      x = binomial_random(n, p/sum);
    }
    n = n - x;
    sum = sum - p;
    destination[i] = x;
  }

  // get the last one
  destination[nb_colors-1] = n;
}

void MersenneTwister::multinomial_roulette(int32_t* destination, double* source, int32_t nb_drawings, int32_t nb_colors)
{
  double rand_value;

  // Initialize output
  memset(destination, 0, nb_colors * sizeof(*destination));

  // Make a table of accumulative probabilities
  // Example : if source was { 0.25, 0.6, 0.15 }
  //           accumulative_probs will be { 0.25, 0.85, 1.0 }
  double* accumulative_probs = new double[nb_colors];
  accumulative_probs[0] = source[0];
  for (int32_t i = 1 ; i < nb_colors ; i++)
  {
    accumulative_probs[i] = accumulative_probs[i-1] + source[i];
  }

  if (accumulative_probs[nb_colors-1] - 1 > 0.0000000001)
  {
    accumulative_probs[nb_colors-1] = 1.0;
  }

  // Do the actual drawing
  for (int32_t i = 0 ; i < nb_drawings ; i++)
  {
    rand_value = random();

    int32_t j = 0;
    while (rand_value >= accumulative_probs[j]) j++;

    //~ printf("rand_value = %f   accumulative_probs[j] = %f\n", rand_value, accumulative_probs[j]);

    destination[j]++;
  }
}


// =================================================================
//                           Protected Methods
// =================================================================
double MersenneTwister::gammln(double X)
// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt).
// When the argument z is an integer, the gamma function is just the familiar factorial
// function, but offset by one, n! = gamma(n + 1).
{
  double x, y, tmp, ser;
  static double cof[6] = {  76.18009172947146,
                            -86.50532032941677,
                            24.01409824083091,
                            -1.231739572450155,
                            0.1208650973866179e-2,
                            -0.5395239384953e-5 };

  y = x = X;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;

  for (int8_t j = 0 ; j <= 5 ; j++)
  {
    ser += cof[j] / ++y;
  }

  return -tmp + log(2.5066282746310005 * ser / x);
}
} // namespace aevol
