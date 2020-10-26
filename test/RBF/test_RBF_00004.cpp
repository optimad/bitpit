/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
/*! @file   test_RBF_00004.cpp
 *  @author Alessandro Alaia (alessandro.alaia@optimad.it)
 *  @brief  Unit test for class RadialFunction
*/

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <iostream>
# include <exception>

// Bitpit
# include "rbf.hpp"

namespace bitpit
{
namespace rbf
{
namespace testing
{
  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int test_rf_constructors( bitpit::rbf::eRBFType type )
  {
    // Scope variables
    int err = 0;

    // Test default constructor.
    std::cout << "    testing default constructor" << std::endl;
    try
    {
        auto rbf = RadialFunct<d, coeff_t>::New( type );
    }
    catch ( std::exception &e )
    {
      std::cout << "bitpit::rbf::testing::test_constructors: ** ERROR ** " << e.what() << std::endl;
      ++err;
    }

    return err;
  }

  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int test_rf_operators( bitpit::rbf::eRBFType type )
  {
    return 0;
  }

  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int run_unit_tests( bitpit::rbf::eRBFType type )
  {
    // Scope variables
    int err = 0;

    // Test constructor(s)
    std::cout << "   testing constructor(s)" << std::endl;
    err += test_rf_constructors<d, coeff_t>( type );

    // Test operator(s)
    std::cout << "   testing operator(s)" << std::endl;
    err += test_rf_operators<d, coeff_t>( type );

    // Test member(s)
    std::cout << "   testing method(s)" << std::endl;

    return err;

  }

  // ------------------------------------------------------------------------ //
  template< size_t d, class coeff_t >
  int rbf_test_suite()
  {
    // Scope variables
    int err = 0;

    // Output message
    std::cout << "- testing d = " << d << std::endl;

    // Run tests for various type of RBF functions
    std::cout << "  testing " + bitpit::rbf::getRBFTag( bitpit::rbf::eRBFType::kWendlandC2 ) + " rbf " << std::endl;
    err += run_unit_tests<d, coeff_t>( bitpit::rbf::eRBFType::kWendlandC2 );
    std::cout << "  testing " + bitpit::rbf::getRBFTag( bitpit::rbf::eRBFType::kMultiQuadric2 ) + " rbf " << std::endl;
    err += run_unit_tests<d, coeff_t>( bitpit::rbf::eRBFType::kMultiQuadric2 );
    std::cout << "  testing " + bitpit::rbf::getRBFTag( bitpit::rbf::eRBFType::kMultiQuadric3_2 ) + " rbf " << std::endl;
    err += run_unit_tests<d, coeff_t>( bitpit::rbf::eRBFType::kMultiQuadric3_2 );
    std::cout << "  testing " + bitpit::rbf::getRBFTag( bitpit::rbf::eRBFType::kMultiQuadric5_2 ) + " rbf " << std::endl;
    err += run_unit_tests<d, coeff_t>( bitpit::rbf::eRBFType::kMultiQuadric5_2 );


    // Output message
    std::cout << err << " unit test(s) failed" << std::endl;
  }
} //end namespace testings
} //end namespace rbf
} //end namespace bitpit

/*! @brief Main driver. */
int main()
{
  // Scope variables
  int err = 0;

  // Run test for various dimensions of the working space and data types

  // Float
  {
    using coeff_t = float;
    std::cout << "* Testing float" << std::endl;
    err += ( bitpit::rbf::testing::rbf_test_suite<1, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<2, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<3, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<5, coeff_t>() );

  }

  // Double
  {
    using coeff_t = double;
    std::cout << "* Testing double" << std::endl;
    err += ( bitpit::rbf::testing::rbf_test_suite<1, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<2, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<3, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<5, coeff_t>() );

  }

  // Long double
  {
    using coeff_t = long double;
    std::cout << "* Testing long double" << std::endl;
    err += ( bitpit::rbf::testing::rbf_test_suite<1, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<2, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<3, coeff_t>() );
    err += ( bitpit::rbf::testing::rbf_test_suite<5, coeff_t>() );
  }

  return err;
}
