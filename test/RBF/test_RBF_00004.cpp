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
# include <vector>
# include <random>

// Bitpit
# include "rbf.hpp"
# include "metaprogramming.hpp"

namespace bitpit
{
namespace rbf
{
namespace testing
{
  template< std::size_t Dim, class CoordT >
  bool compare( RF<Dim, CoordT> const &in1, RF<Dim, CoordT> const &in2 )
  {
    int err = 0;
    if ( in1.getType() != in2.getType() ) {
      std::cout << "**ERROR** Mismatching type: found "
                << bitpit::rbf::getRBFTag(in1.getType() )
                << " vs "
                << bitpit::rbf::getRBFTag(in2.getType() )
                << std::endl;
      ++err;
    }
    if ( in1.getRadius() != in2.getRadius() ) {
      std::cout << "**ERROR** Mimatching radius: found "
                << in1.getRadius()
                << " vs "
                << in2.getRadus()
                << std::endl;
      ++err;
    }
    if( in1.getCenter() != in2.getCenter() ) {
      std::cout << "**ERROR** Mimatching geometry kernel: found "
                << in1.getCenter()
                << " vs "
                << in2.getCenter()
                << std::endl;
      ++err;
    }
    if (in1.getNumberOfParameters() != in2.getNumberOfParameters() ) {
      std::cout << "**ERROR** Mimatching nr. of parameters: found "
                << in1.getNumberOfParameters()
                << " vs "
                << in2.getNumberOfParameters()
                << std::endl;
      ++err;
    }
    else {
      auto p1 = in1.getParameters(),
           p2 = in2.getParameters();
      for (int i = 0, n = in1.getNumberOfParameters(); i < n; ++i )
      {
        if ( p1[i] !=  p2[i] ) {
          std::cout << "**ERROR** Mimatching value for parameter #" << i << ": found "
                    << in1.getNumberOfParameters()
                    << " vs "
                    << in2.getNumberOfParameters()
                    << std::endl;
          ++err;
        }
      } //next i
    }

    return err;

  }

  template< std::size_t Dim, class CoordT >
  class EvaluationTester
  {
    private:
    using rf_t      = RF<Dim, CoordT>;
    using coord_t   = typename rf_t::coord_t;
    using point_t   = typename rf_t::point_t;

    public:
    EvaluationTester( bitpit::rbf::eRBFType type, std::size_t N = 100, coord_t e = std::numeric_limits<coord_t>::epsilon() ) :
      rf( RF<Dim, CoordT>::New( type ) )
      , test_values(N)
      , n_trials(N)
      , tol(e)
    {
      pars = new coord_t[ rf->getNumberOfParameters() ];
    }
    int operator()() const
    {
      int err = 0;
      init();
      test();
      err += check();
      init();
      test();
      err += check();
      cleanup();
      return err;
    }
    protected:
    typename rf_t::rf_funct_t initGenerator() const
    {
        auto p = rf->getParameters();
        switch( rf->getType() )
        {
          case( bitpit::rbf::eRBFType::kWendlandC2 ) :
            return &wendland_c2<coord_t>;
          case( bitpit::rbf::eRBFType::kHardy ) :
              return std::bind( &bitpit::rbf::generalized_multiquadrics<coord_t, 1, 2>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kMultiQuadric2 ) :
            return std::bind( &bitpit::rbf::generalized_multiquadrics<coord_t, 2, 1>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kMultiQuadric3_2 ) :
            return std::bind( &bitpit::rbf::generalized_multiquadrics<coord_t, 3, 2>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kMultiQuadric5_2 ) :
            return std::bind( &bitpit::rbf::generalized_multiquadrics<coord_t, 5, 2>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kGaussian ) :
            return std::bind( &bitpit::rbf::gaussian<coord_t>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kLinear ) :
            return std::bind( &bitpit::rbf::linear<coord_t>, std::placeholders::_1 );
          case( bitpit::rbf::eRBFType::kQuadratic ) :
            return std::bind( &bitpit::rbf::generalized_power<coord_t, 2>, std::placeholders::_1 );
          case( bitpit::rbf::eRBFType::kCubic ) :
            return std::bind( &bitpit::rbf::generalized_power<coord_t, 3>, std::placeholders::_1 );
          case( bitpit::rbf::eRBFType::kQuartic ) :
            return std::bind( &bitpit::rbf::generalized_power<coord_t, 4>, std::placeholders::_1 );
          case( bitpit::rbf::eRBFType::kThinPlateSpline ) :
            return std::bind( &bitpit::rbf::thin_plate_spline<coord_t>, std::placeholders::_1 );
          case( bitpit::rbf::eRBFType::kPolyharmonic2 ) :
            return std::bind( &bitpit::rbf::polyharmonic<coord_t, 2>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kPolyharmonic3 ) :
            return std::bind( &bitpit::rbf::polyharmonic<coord_t, 3>, std::placeholders::_1, p[0] );
          case( bitpit::rbf::eRBFType::kPolyharmonic4 ) :
            return std::bind( &bitpit::rbf::polyharmonic<coord_t, 4>, std::placeholders::_1, p[0] );
          default: throw std::runtime_error( "bitpit::rbf::testing: **ERROR** unsupported rbf type" );
        }
    }
    bool  init() const
    {
      // Scope variables
      std::random_device dev;
      std::default_random_engine eng(dev());
      std::uniform_real_distribution<coord_t> dist( (coord_t)0, (coord_t)1);
      coord_t r;
      point_t c;

      // Initialize the radial funct
      {
        for ( std::size_t j = 0; j < Dim; ++j ) {
          c[j] = dist(eng);
          rf->center[j] = c[j];
        }
        r = dist(eng) + (coord_t)1;
        rf->radius = r;
        if ( std::size_t n = rf->getNumberOfParameters() )
        {
          for ( std::size_t i = 0; i < n; ++i )
            pars[i] = (coord_t)1;//dist(eng);
          rf->setParameters( pars );
        }
      }

      // Generate a cloud of scattered points in a space of dimension Dim
      {
        test_values.resize(n_trials);
        for ( std::size_t i = 0; i < n_trials; ++i ) {
          auto &p = std::get<0>( test_values[i] );
          for ( std::size_t j = 0; j < Dim; ++j )
            p[j] = dist(eng);
        } //next i
      }

      // Generate expected values
      {
        auto f = initGenerator();
        for ( std::size_t i = 0; i < n_trials; ++i ) {
          auto &v = std::get<1>( test_values[i] );
          const auto &p = std::get<0>( test_values[i] );
          v = f( norm2(p - c)/r );
        }
      }

      return true;
    }
    bool  cleanup() const
    {
      delete rf;
      delete [] pars;
      std::vector<std::tuple<point_t, coord_t, coord_t>>(0).swap( test_values );
    }
    bool  test() const
    {
      for ( std::size_t i = 0; i < n_trials; ++i ) {
        auto &trial = test_values[i];
        auto &v = std::get<2>( trial );
        const auto &p = std::get<0>( trial );
        v = rf->operator()(p);
      } //next i
      return true;
    }
    int   check() const
    {
      int err = 0;
      for ( std::size_t i = 0; i < n_trials; ++i ) {
        const auto &trial = test_values[i];
        const auto &v = std::get<2>( trial );
        const auto &e = std::get<1>( trial );
        auto ee = std::abs(v-e);
        if ( ee > tol )
        {
          std::cout << "**ERROR** Wrong value from rf evaluation. Found "
                    << v << ", expecting " << e
                    << " (err: " << ee << ")"
                    << std::endl;
          ++err;
        }
      } //next i
      return err;
    }
    protected:
    mutable std::vector<std::tuple<point_t, coord_t, coord_t>>  test_values;
    mutable RF< Dim, CoordT >     *rf;
    mutable coord_t               *pars;
    coord_t                       tol;
    std::size_t                   n_trials;
  }; //end class EvaluationTester

  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int test_rf_constructors( bitpit::rbf::eRBFType type )
  {
    // Scope variables
    int err = 0;

    // Test default constructor.
    std::cout << "    testing constructor #1" << std::endl;
    try
    {
        auto rbf = RF<d, coeff_t>::New( type );
        delete rbf;
    }
    catch ( std::exception &e )
    {
      std::cout << "bitpit::rbf::testing::test_rf_constructors: ** ERROR ** " << e.what() << std::endl;
      ++err;
    }

    return err;
  }

  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int test_rf_operators( bitpit::rbf::eRBFType type )
  {
    // Scope variables
    int err = 0;

    // Run rf evaluations on random data
    std::cout << "    testing evaluation operator" << std::endl;
    try
    {
      EvaluationTester<d, coeff_t> tester(type, 100);
      err += tester();
    }
    catch ( std::exception &e )
    {
      std::cout << "**ERROR** " << e.what() << std::endl;
      ++err;
    }

    // Output message
    std::cout << "    test completed with " << err << " error(s)" << std::endl;
    return err;
  }

  // ------------------------------------------------------------------------ //
  template<size_t d, class coeff_t>
  int test_rf_setters_getters( bitpit::rbf::eRBFType type )
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
    err += ( test_rf_constructors<d, coeff_t>( type ) != 0 );

    // Test member(s)
    std::cout << "   testing setter(s)/getter(s)" << std::endl;
    err += ( test_rf_setters_getters<d, coeff_t>( type ) != 0 );

    // Test operator(s)
    std::cout << "   testing operator(s)" << std::endl;
    err += ( test_rf_operators<d, coeff_t>( type ) != 0 );

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
    for ( int i = bitpit::rbf::eRBFType::kUndefined+1, n = bitpit::rbf::eRBFType::kUserDefined; i < n; ++i ) {
      std::cout << "  testing " + bitpit::rbf::getRBFTag( static_cast<bitpit::rbf::eRBFType>(i) ) + " rbf " << std::endl;
      err += run_unit_tests<d, coeff_t>( static_cast<bitpit::rbf::eRBFType>(i) );
    } //next type



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
