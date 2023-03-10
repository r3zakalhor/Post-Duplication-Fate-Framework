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
//
// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MersenneTwister
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2003, Richard J. Wagner
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote
//      products derived from this software without specific prior written
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
//
//     When you use this, send an email to: matumoto@math.keio.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
// when you write.


#ifndef AEVOL_MERSENNE_TWISTER_H_
#define AEVOL_MERSENNE_TWISTER_H_


#include <inttypes.h>

#define MT_RAND_MAX         4294967295.0
#define MT_RAND_MAX_PLUS_1  4294967296.0

// That is  MT_RAND_MAX = OxFFFFFFFF
//          MT_RAND_MAX = OxFFFFFFFF + 1 = 2e32

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MersenneTwister object)

#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "AeTime.h"
#include <math.h>
#include <zlib.h>
#include "macros.h"

namespace aevol {

class MersenneTwister
{
  // Data
  public:
    enum { N = 624 };       // length of state vector
    enum { SAVE = N + 1 };  // length of array for save()

  protected:
    enum { M = 397 };  // period parameter

    uint32_t state[N];   // internal state
    uint32_t *pNext;     // next value to get from state
    int16_t  left;       // number of values left before reload needed


  //Methods
  public:
    // Constructors
    MersenneTwister(const uint32_t& oneSeed);  // initialize with a simple uint32_t
    MersenneTwister(gzFile backup_file);      // Load from a gz backup file
    MersenneTwister(const MersenneTwister & model);

    // Destructors
    virtual ~MersenneTwister();

    // Main generator
    inline uint32_t rand_next(); // Draw a 32-bit-long integer in [0, 2e32[

    // AEvol wrappers
    inline double   random();         // Double in [0, 1[ (uniform distribution)
    inline int8_t   random(int8_t max);   // ~
    inline int16_t  random(int16_t max);  // ~ > Integer in [0, max[ (uniform distribution)
    inline int32_t  random(int32_t max);  // ~
    int32_t         binomial_random(int32_t nb, double prob); // Binomial drawing of parameters (nb, prob)
    double          gaussian_random();
    void            multinomial_drawing(int32_t* destination, double* source, int32_t nb_drawings, int32_t colors);
    void            multinomial_roulette(int32_t* destination, double* source, int32_t nb_drawings, int32_t colors);
    // Multinomial drawing of parameters (nb, {source[0], source[1], ... source[colors-1]})
    // WARNING : The roulette implementation as it is, is 40 times slower that the other implementation

    //~ inline int32_t   norm_geometric_random(double p, int32_t max_nb_trials);


    // Re-seeding functions with same behavior as initializers
    inline void seed(const uint32_t oneSeed);

    // Saving and loading generator state
    inline void save(uint32_t* saveArray) const;  // to array of size SAVE
    inline void load(uint32_t *const loadArray);  // from such array
    inline void write_to_backup(gzFile backup_file) const;

  protected:
    MersenneTwister()
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };
    inline void initialize(const uint32_t oneSeed);
    inline void reload();
    inline uint32_t hiBit(const uint32_t& u) const { return u & 0x80000000UL; }
    inline uint32_t loBit(const uint32_t& u) const { return u & 0x00000001UL; }
    inline uint32_t loBits(const uint32_t& u) const { return u & 0x7fffffffUL; }
    inline uint32_t mixBits(const uint32_t& u, const uint32_t& v) const
      { return hiBit(u) | loBits(v); }
    inline uint32_t twist(const uint32_t& m, const uint32_t& s0, const uint32_t& s1) const
      { return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL); }
    static inline uint32_t hash(time_t t, clock_t c);
    double gammln(double X);
};



inline uint32_t MersenneTwister::rand_next()
{
  // Pull a 32-bit integer from the generator state
  // Every other access function simply transforms the numbers extracted here

  if(left == 0) reload();
  --left;

  register uint32_t s1;
  s1 = *pNext++;
  s1 ^= (s1 >> 11);
  s1 ^= (s1 <<  7) & 0x9d2c5680UL;
  s1 ^= (s1 << 15) & 0xefc60000UL;
  return (s1 ^ (s1 >> 18));
}

inline double MersenneTwister::random()  // Double in [0, 1) (uniform distribution)
{
  return ((double)rand_next()) / MT_RAND_MAX_PLUS_1;
}

inline int8_t MersenneTwister::random(int8_t max) // Integer in [0, max[ (uniform distribution)
{
  return (int8_t)(((double)max) * (((double)rand_next()) / MT_RAND_MAX_PLUS_1));
}

inline int16_t MersenneTwister::random(int16_t max) // Integer in [0, max[ (uniform distribution)
{
  return (int16_t)(((double)max) * (((double)rand_next()) / MT_RAND_MAX_PLUS_1));
}

inline int32_t MersenneTwister::random(int32_t max) // Integer in [0, max[ (uniform distribution)
{
  return (int32_t)(((double)max) * (((double)rand_next()) / MT_RAND_MAX_PLUS_1));
}

inline void MersenneTwister::seed(const uint32_t oneSeed)
{
  // Seed the generator with a simple uint32_t
  initialize(oneSeed);
  reload();
}

inline void MersenneTwister::initialize(const uint32_t seed)
{
  // Initialize generator state with seed
  // See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
  // In previous versions, most significant bits (MSBs) of the seed affect
  // only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
  register uint32_t *s = state;
  register uint32_t *r = state;
  register int16_t  i  = 1;
  *s++ = seed & 0xffffffffUL;
  for(; i < N; ++i)
  {
    *s++ = (1812433253UL * (*r ^ (*r >> 30)) + i) & 0xffffffffUL;
    r++;
  }
}


void MersenneTwister::reload()
{
  // Generate N new values in state
  // Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
  register uint32_t *p = state;
  register int16_t i;
  for(i = N - M; i--; ++p)
    *p = twist(p[M], p[0], p[1]);
  for(i = M; --i; ++p)
    *p = twist(p[M-N], p[0], p[1]);
  *p = twist(p[M-N], p[0], state[0]);

  left = N, pNext = state;
}


uint32_t MersenneTwister::hash(time_t t, clock_t c)
{
  // Get a uint32_t from t and c
  // Better than uint32_t(x) in case x is floating point in [0,1]
  // Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)

  static uint32_t differ = 0;  // guarantee time-based seeds will change

  uint32_t h1 = 0;
  unsigned char *p = (unsigned char *) &t;
  for(size_t i = 0; i < sizeof(t); ++i)
  {
    h1 *= UCHAR_MAX + 2U;
    h1 += p[i];
  }
  uint32_t h2 = 0;
  p = (unsigned char *) &c;
  for(size_t j = 0; j < sizeof(c); ++j)
  {
    h2 *= UCHAR_MAX + 2U;
    h2 += p[j];
  }
  return (h1 + differ++) ^ h2;
}


void MersenneTwister::save(uint32_t* saveArray) const
{
  register uint32_t *sa = saveArray;
  register const uint32_t *s = state;
  register int16_t i = N;
  for(; i--; *sa++ = *s++) {}
  *sa = left;
}


void MersenneTwister::load(uint32_t *const loadArray)
{
  register uint32_t *s = state;
  register uint32_t *la = loadArray;
  register int16_t i = N;
  for(; i--; *s++ = *la++) {}
  left = *la;
  pNext = &state[N-left];
}

void MersenneTwister::write_to_backup(gzFile backup_file) const
{
  uint32_t saveArray[SAVE];
  save(saveArray);
  gzwrite(backup_file, saveArray, SAVE * sizeof(saveArray[0]));
}
} // namespace aevol

#endif
