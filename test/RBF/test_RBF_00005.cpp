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
/*! @file   test_RBF_00005.cpp
 *  @author Alessandro Alaia (alessandro.alaia@optimad.it)
 *  @brief  Unit test for class RBF
*/
# undef __FULL_OUT__

// ========================================================================= //
// INCLUDES                                                                  //
// ========================================================================= //

// Standard Template Library
# include <random>
# include <vector>
# include <iostream>
# include <string>
# include <exception>

// bitpit
# include "rbf.hpp"

namespace bitpit
{
namespace rbf
{
  namespace testing
  {
    template< class T >
    struct test_tol{ static const T value; };
    template<>
    const float test_tol<float>::value = 1.e-4F;
    template<>
    const double test_tol<double>::value = 1.e-6;
    template<>
    const long double test_tol<long double>::value = 1.e-8L;

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    class RBFTester
    {
      private:
      using coord_t = CoordT;
      using rbf_t   = bitpit::rbf::RFBasis< Dim, coord_t >;

      private:
      mutable rbf_t rbf;
      mutable std::vector<coord_t> data_v;
      mutable std::vector<typename rbf_t::point_t> data_p;
      mutable std::size_t nrbf, irbf;
      CoordT tol;

      public:
      RBFTester( std::size_t n, bitpit::rbf::eRBFType type, CoordT tolerance = test_tol<CoordT>::value ) :
        data_v()
        , data_p()
        , rbf( n, type )
        , nrbf(n)
        , irbf(-1)
        , tol( tolerance )
      {}
      private:
      bool init( std::size_t np ) const
      {

        // Scope variables
        std::random_device dev;
        std::default_random_engine eng(dev());
        std::uniform_real_distribution<coord_t> dist( (coord_t)0, (coord_t)1);
        std::vector<std::size_t> ids(np);
        std::size_t M = np;

        // Functor(s)
        auto sample_wo_replacement = [&ids, &M, &dist, &eng]( )->std::size_t
        {
          if ( M == 0 )
            throw std::runtime_error(
              "bitpit::rbf::testing::RBFTester: ** ERROR ** all ids in the range have been extracted!"
            );
          auto p = dist(eng);
          auto i = static_cast<std::size_t>( (coord_t)(M-1) * p );
          std::size_t extracted_id = ids.at(i);
          std::swap( ids.at(i), ids.at(M-1) );
          --M;
          return extracted_id;
        };//end sample_wo_replacement

        // General initializations
        std::iota(ids.begin(), ids.end(), 0 );

        // Generate scattered data points
        {
          data_p.resize(np);
          for ( std::size_t i = 0; i < np; ++i )
          {
            auto &p = data_p.at(i);
            for ( auto j = 0; j < Dim; ++j )
              p[j] = dist(eng);
          } //next i
        }

        // Assign the value of radius and center to each rf.
        {
          // Set the radius for all radial functions (random value)
          rbf.setRadius( (coord_t) dist(eng) + (coord_t)1 / (coord_t)5 );

          // Set the center for each radial basis function
          for ( std::size_t i = 0; i < nrbf; ++i )
          {
            auto &rf = *( rbf.at(i).second );
            auto id = sample_wo_replacement();
            rf.center = data_p.at(id);
          } //next i
        }

        // Choose a random rf within the basis to generate data
        {
          std::uniform_int_distribution<std::size_t> idist(0, nrbf-1);
          irbf = idist(eng);
        }

        // Generate data values with the chosen distribution
        {
          const auto &rf = *( rbf.at(irbf).second );
          data_v.resize(np);
          for ( std::size_t i = 0; i < np; ++i )
            data_v[i] = rf( data_p.at(i) );
          //next i
        }

        return true;
      }
      bool test( ) const
      {
        if ( !bitpit::rbf::computeRBFWeights( data_p, data_v, rbf ) )
          throw std::runtime_error(
            "bitpit::rbf::testing::RBFTester: ** ERROR ** Failed to compute RBF weights"
          );
        # ifdef __FULL_OUT__
        rbf.display( std::cout, 8 );
        # endif
        return true;
      }
      bool cleanup() const
      {
        /*nothing to be done*/
        return true;
      }
      int check() const
      {
        // Scope variables
        int err = 0;

        // Check outcome of interpolation
        for ( std::size_t i = 0, n = data_p.size(); i < n; ++i )
        {
          const auto &p = data_p.at(i);
          auto v = rbf(p);
          const auto &ve = data_v.at(i);
          auto e = std::abs( v - ve );
          if ( e > tol )
          {
            std::cout << "** ERROR ** Wrong value from interpolation. Found "
                      << v << ", expecting " << ve
                      << " (err: " << e << ")"
                      << std::endl;
            ++err;
          }
        } //next i

        return err;
      }
      public:
      int operator()( std::size_t np ) const
      {
        // Scope variables
        int err = 0;

        if ( !init(np) )  return err;
        if ( !test() )    return err;
        err = check();
        if ( !cleanup() ) return err;

        return err;
      }
    }; //end class RBFTester

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int test_rbf_constructors( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int test_rbf_operators( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int test_rbf_getterse_setters( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int test_rbf_interpolation( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      // Output message
      std::vector<std::size_t> tested_size = {{10, 20, 30}};
      for ( std::size_t nrbf : tested_size ) {
        std::cout << "      testing Nrbf = " << nrbf << " Ndata = " << nrbf << std::endl;
        try
        {
          RBFTester<Dim, CoordT> tester( nrbf, type );
          err += tester( nrbf );
        }
        catch ( std::exception &e )
        {
          std::cout << "** ERROR ** bitpit::rbf::test_rbf_interpolation: " << e.what() << std::endl;
          ++err;
        }
        std::cout << "      testin Nrbf = " << nrbf << " Ndata = " << 2*nrbf << std::endl;
        try
        {
          RBFTester<Dim, CoordT> tester( nrbf, type );
          err += tester( 2*nrbf );
        }
        catch ( std::exception &e )
        {
          std::cout << "** ERROR ** bitpit::rbf::test_rbf_interpolation: " << e.what() << std::endl;
          ++err;
        }
      } //next nrbf

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int rbf_basic_unit_tests( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      // Test constructor
      std::cout << "    - testing constructor(s)" << std::endl;
      err += bitpit::rbf::testing::test_rbf_constructors<Dim,CoordT>(type);

      // Test operators
      std::cout << "    - testing operator(s)" << std::endl;
      err += bitpit::rbf::testing::test_rbf_operators<Dim,CoordT>(type);

      // Test getters/setters
      std::cout << "    - testing getter/setter(s)" << std::endl;
      err += bitpit::rbf::testing::test_rbf_getterse_setters<Dim,CoordT>(type);

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int rbf_interp_unit_tests( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      // Test interpolation
      if ( std::is_same< CoordT, float >::value || std::is_same< CoordT, double >::value )
      {
        std::cout << "    - testing interpolation" << std::endl;
        err += bitpit::rbf::testing::test_rbf_interpolation<Dim,CoordT>(type);
      }

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int rbf_unit_tests( bitpit::rbf::eRBFType type )
    {
      // Scope variables
      int err = 0;

      // Output message
      std::cout << "    testing RBF type \"" << bitpit::rbf::getRBFTag( type ) << "\"" << std::endl;

      // Run tests
      err += rbf_basic_unit_tests<Dim,CoordT>(type);
      err += rbf_interp_unit_tests<Dim,CoordT>(type);

      // Output message
      std::cout << "    test completed with " << err << " error(s)" << std::endl;

      return err;
    }

    // --------------------------------------------------------------------- //
    template< std::size_t Dim, class CoordT >
    int rbf_test_suite( )
    {
      // Scope variables
      int err = 0;

      for ( int type = bitpit::rbf::eRBFType::kUndefined+1
          , n = bitpit::rbf::eRBFType::kUserDefined
          ; type < n
          ; ++type ) {
        err += ( rbf_unit_tests<Dim,CoordT>( static_cast< bitpit::rbf::eRBFType >(type) ) != 0 );
      } //next type

      return err;
    }

    // --------------------------------------------------------------------- //
    template< class CoordT >
    int rbf_test_dim( )
    {
      // Scope variables
      int err = 0;

      // Run test for various dimensions of the working space.
      std::cout << "  * Testing Dim = 1" << std::endl;
      err += rbf_test_suite<1, CoordT>();
      std::cout << "  * Testing Dim = 2" << std::endl;
      err += rbf_test_suite<2, CoordT>();
      std::cout << "  * Testing Dim = 3" << std::endl;
      err += rbf_test_suite<3, CoordT>();
      std::cout << "  * Testing Dim = 4" << std::endl;
      err += rbf_test_suite<4, CoordT>();
      std::cout << "  * Testing Dim = 5" << std::endl;
      err += rbf_test_suite<5, CoordT>();

      return err;
    }

  } //end namespace testing
} //end namespace rbf
} //end namespace bitpit

// Main driver
int main()
{
  // Scope variables
  int err = 0;

  // Test float
  std::cout << "Testing CoordT = float" << std::endl;
  err = bitpit::rbf::testing::rbf_test_dim<float>();
  std::cout << err << " unit test(s) failed"  << std::endl;

  // Test double
  std::cout << "Testing CoordT = double" << std::endl;
  err = bitpit::rbf::testing::rbf_test_dim<double>();
  std::cout << err << " unit test(s) failed"  << std::endl;

  // Test long double
  std::cout << "Testing CoordT = long double" << std::endl;
  err = bitpit::rbf::testing::rbf_test_dim<long double>();
  std::cout << err << " unit test(s) failed"  << std::endl;

  return err;


}
