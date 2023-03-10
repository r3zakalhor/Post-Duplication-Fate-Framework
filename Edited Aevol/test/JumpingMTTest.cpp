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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================

#include <cstdio>
#include <cstdlib>

#include "gtest/gtest.h"

#include "JumpingMT.h"


// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================

using namespace aevol;


//############################################################################
//                                                                           #
//                         Class JumpingMTTest                          #
//                                                                           #
//############################################################################
class JumpingMTTest : public testing::Test {
 protected:
  virtual void SetUp(void);
  virtual void TearDown(void);

  JumpingMT* obj1;
  JumpingMT* obj2;
};
// ===========================================================================
//                               Static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================

// ===========================================================================
//                                  Destructors
// ===========================================================================

// ===========================================================================
//                                   Operators
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void JumpingMTTest::SetUp(void)
{
  obj1 = new JumpingMT(56873);
  obj2 = new JumpingMT(*obj1);
}

void JumpingMTTest::TearDown(void)
{
  delete obj1;
  delete obj2;
}

TEST_F(JumpingMTTest, test1)
{
  for (int i = 0 ; i < 10 ; i++)
  {
    obj1->random();
    obj2->random();
  }

  EXPECT_TRUE(obj1->random() == obj2->random());
  EXPECT_TRUE(obj1->random((int8_t) INT8_MAX) == obj2->random((int8_t) INT8_MAX));
  EXPECT_TRUE(obj1->random((int16_t) INT16_MAX) == obj2->random((int16_t) INT16_MAX));
  EXPECT_TRUE(obj1->random((int32_t) INT32_MAX) == obj2->random((int32_t) INT32_MAX));
  EXPECT_TRUE(obj1->random((int64_t) INT64_MAX) == obj2->random((int64_t) INT64_MAX));
  EXPECT_TRUE(obj1->binomial_random(435, 0.5) == obj2->binomial_random(435, 0.5));
  EXPECT_TRUE(obj1->gaussian_random() == obj2->gaussian_random());
}

TEST_F(JumpingMTTest, test2)
{
  EXPECT_TRUE(true);
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
