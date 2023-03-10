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




// ============================================================================
//                                   Includes
// ============================================================================
#include "Utils.h"

#include "JumpingMT.h"
#include "aevol_version.h"

using std::string;
using std::cout;
using std::endl;


namespace aevol {
//##############################################################################
//
//                           Class Utils
//
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================


/**
 * Compute value(t+1) with the provided parameters and prng
 *
 * value(t+1) = value(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
 */
void Utils::ApplyAutoregressiveStochasticProcess(double& value,
                                                 double sigma,
                                                 int16_t tau,
                                                 JumpingMT& prng) {
  value = value * (1.0 - 1.0/tau) +
          (sigma/tau) * sqrt(2*tau- 1.0) * prng.gaussian_random();
}

/**
 * Print an error message and exit (with error status)
 *
 * \param msg message to be printed
 * \param file should be __FILE__
 * \param line should be __LINE__
 */
void Utils::ExitWithDevMsg(const std::string& msg,
                           const std::string& file, int line) {
  cout << file << ":" << line << ": error: " << msg << endl;
  exit(EXIT_FAILURE);
}

/**
 * Print an error message and exit (with error status)
 *
 * \param msg message to be printed
 * \param bin the binary that caused the error
 */
void Utils::ExitWithUsrMsg(const string& msg) {
  cout << "error: " << msg << endl;
  exit(EXIT_FAILURE);
}

/**
 * Print aevol version number
 */
void Utils::PrintAevolVersion() {
  printf("aevol %s\n", aevol::version_string);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
