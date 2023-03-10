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


#ifndef AEVOL_UTILS_H_
#define AEVOL_UTILS_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cassert>
#include <cstdlib>

#include <string>
#include <iostream>


namespace aevol {
class JumpingMT;

class Utils {
 public :
  Utils() = delete;
  Utils(const Utils&) = delete;
  Utils(Utils&&) = delete;
  ~Utils() = delete;
  static inline int32_t mod(int32_t a, int32_t b);
  static inline int64_t mod(int64_t a, int64_t b);
  static inline int32_t min(int32_t a, int32_t b);
  static inline int32_t max(int32_t a, int32_t b);
  static inline void exchange(int32_t& a, int32_t& b);
  static void ApplyAutoregressiveStochasticProcess(double& value,
                                                   double sigma,
                                                   int16_t tau,
                                                   JumpingMT& prng);
  static int16_t hamming(const char* str1, const char* str2);

  static void ExitWithUsrMsg(const std::string& msg);
  static void ExitWithDevMsg(const std::string& msg,
                             const std::string& file, int line);
  static void PrintAevolVersion();
};

inline int32_t Utils::mod(int32_t a, int32_t b)
{

  assert(b > 0);

  while (a < 0)  a += b;
  while (a >= b) a -= b;

  return a;
  //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
}

inline int64_t Utils::mod(int64_t a, int64_t b)
{

  assert(b > 0);

  while (a < 0)  a += b;
  while (a >= b) a -= b;

  return a;
  //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
}

inline int32_t Utils::min(int32_t a, int32_t b)
{
  return ((a < b)? a : b);
}

inline int32_t Utils::max(int32_t a, int32_t b)
{
  return ((a > b)? a : b);
}

inline void Utils::exchange(int32_t& a, int32_t& b)
{
  int32_t tmp = a;
  a = b;
  b = tmp;
}
} // namespace aevol
#endif // AEVOL_UTILS_H_
