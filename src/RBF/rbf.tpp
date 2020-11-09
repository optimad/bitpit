/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2019 OPTIMAD engineering Srl
*
*  --------------------------------------------------------------------------
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
/*! @file     rbf.tpp
 *  @author   Alessandro Alaia (mail to: alessandro.alaia@optimad.it)
 *  @brief    Template implementations for Radial Basis Functions.
*/
namespace bitpit
{
namespace rbf
{
  // ======================================================================== //
  // IMPLEMENTATION OF RADIAL FUNCTION SUPPORTED BY bitpit.                   //
  // ======================================================================== //

  // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT  wendland_c2( CoordT r )
  {
    return r > (CoordT)1 ?
      (CoordT)0 :
      std::pow( (CoordT)1 - r, 4 ) * ( (CoordT)4 * r + (CoordT)1 );
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, unsigned Alpha, unsigned Beta, typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) && (Beta > 0) >::type* >
  CoordT generalized_multiquadrics( CoordT r, CoordT c )
  {
    return CoordT(1)/std::pow( c*c + r*r, CoordT(Alpha)/CoordT(Beta) );
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT gaussian( CoordT r, CoordT c )
  {
    return c * std::exp( - (r*r) );
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT linear( CoordT r )
  {
    return r;
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, int Alpha, typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) >::type* >
  CoordT generalized_power( CoordT r )
  {
    return std::pow( r, Alpha );
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT    thin_plate_spline( CoordT r )
  {
    return r*r * std::log(r);
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, int Beta, typename std::enable_if< std::is_floating_point<CoordT>::value && (Beta >= 0) >::type* >
  CoordT polyharmonic( CoordT r, CoordT c )
  {
    return ( Beta % 2 == 0 ? (CoordT)-1 : (CoordT)1 ) * std::pow( c, 1+Beta) * std::pow(r, 2*Beta) * std::log(r);
  }

  // ======================================================================== //
  // IMPLEMENTATION OF CLASS RF                                               //
  // ======================================================================== //


   // Static member function(s) ============================================== //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RF<D,C>* RF<D,C>::New( eRBFType type )
  {
    switch( type )
    {
      default: {
        throw std::runtime_error(
          "rbf::RF::New: Undefined RF type."
        );
      } //end default
      case( bitpit::rbf::eRBFType::kWendlandC2 ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::wendland_c2<coord_t>
        );
        out->mHasCompactSupport = true;
        return out;
      } //end case kWendlandC2
      case( bitpit::rbf::eRBFType::kGaussian ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::gaussian<coord_t>
        );
        out->mHasCompactSupport = false;
        out->mParams[0] = 1;
        return out;
      }
      case( bitpit::rbf::eRBFType::kHardy ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::generalized_multiquadrics<coord_t, 1, 2>
        );
        out->mHasCompactSupport = false;
        out->mParams[0] = 1;
        return out;
      } //end case kHardy
      case( bitpit::rbf::eRBFType::kMultiQuadric2 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::generalized_multiquadrics<coord_t, 2, 1 >
        );
        out->mHasCompactSupport = false;
        out->mParams[0] = 1;
        return out;
      } //end case kMultiQuadrics2
      case( bitpit::rbf::eRBFType::kMultiQuadric3_2 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::generalized_multiquadrics<coord_t, 3, 2>
        );
        out->mHasCompactSupport = false;
        out->mParams[0] = 1;
        return out;
      } //end case kMultiQuadric3_2
      case( bitpit::rbf::eRBFType::kMultiQuadric5_2 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::generalized_multiquadrics<coord_t, 5, 2>
        );
        out->mHasCompactSupport = false;
        out->mParams[0] = 1;
        return out;
      } //end case kMultiQuadric3_2
      case( bitpit::rbf::eRBFType::kLinear ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::linear<coord_t>
        );
        out->mHasCompactSupport = false;
        return out;
      } //end case kLinear
      case( bitpit::rbf::eRBFType::kQuadratic ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::generalized_power<coord_t, 2>
        );
        out->mHasCompactSupport = false;
        return out;
      } //end case kQuadratic
      case( bitpit::rbf::eRBFType::kCubic ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::generalized_power<coord_t, 3>
        );
        out->mHasCompactSupport = false;
        return out;
      } //end case kCubic
      case( bitpit::rbf::eRBFType::kQuartic ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::generalized_power<coord_t, 4>
        );
        out->mHasCompactSupport = false;
        return out;
      } //end case kQuartic
      case( bitpit::rbf::eRBFType::kThinPlateSpline ): {
        auto out = new bitpit::rbf::RF<D, C>(
          type,
          &bitpit::rbf::thin_plate_spline<coord_t>
        );
        out->mHasCompactSupport = false;
        return out;
      } //end case kThinPlateSpline
      case( bitpit::rbf::eRBFType::kPolyharmonic2 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::polyharmonic<coord_t, 2>
        );
        out->mParams[0] = (coord_t) 1;
        out->mHasCompactSupport = false;
        return out;
      } //end case kPolyharmonic2
      case( bitpit::rbf::eRBFType::kPolyharmonic3 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::polyharmonic<coord_t, 3>
        );
        out->mParams[0] = (coord_t) 1;
        out->mHasCompactSupport = false;
        return out;
      } //end case kPolyharmonic3
      case( bitpit::rbf::eRBFType::kPolyharmonic4 ): {
        auto out = new bitpit::rbf::RFP<D, 1, C>(
          type,
          &bitpit::rbf::polyharmonic<coord_t, 4>
        );
        out->mParams[0] = (coord_t) 1;
        out->mHasCompactSupport = false;
        return out;
      } //end case kPolyharmonic3
    } //end switch

    return nullptr;
  }

  // Constructor(s) ========================================================== //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RF<D,C>::RF() :
    mFunct(nullptr)
    , radius()
    , center()
    , mType( bitpit::rbf::eRBFType::kUndefined )
  {
    setDefault();
  }

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RF<D,C>::RF( bitpit::rbf::eRBFType type, coord_t (*f)(coord_t) ) :
    mFunct(f)
    , mType(type)
  {
    setDefault();
  }

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RF<D,C>::RF( bitpit::rbf::eRBFType type, rf_funct_t const &f ) :
    mFunct(f)
    , mType(type)
  {
    setDefault();
  }

  // Operator(s) ============================================================ //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RF<D,C>::coord_t RF<D,C>::operator()( const point_t &coords ) const
  {
    return mFunct( norm2( coords - center )/radius );
  }

  // Getter(s)/Info ========================================================= //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  std::size_t	RF<D,C>::getNumberOfParameters() const { return 0; };

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  bool RF<D,C>::hasCompactSupport() const { return mHasCompactSupport; };

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  const typename RF<D,C>::coord_t* RF<D,C>::getParameters() const { return nullptr; };

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  eRBFType RF<D,C>::getType() const { return mType; }

  // ---------------------------------------------------------------------- //
  template< std::size_t D, class C >
  void RF<D,C>::display( std::ostream &out /*= std::cout*/, unsigned int indent /*= 0*/ ) const
  {
    std::string   s( indent, ' ' );
    out << s << "Type:        " << bitpit::rbf::getRBFTag( mType ) << '\n'
        << s << "radius:      " << radius << "\n"
        << s << "center:      [ ";
    for ( std::size_t i = 0; i < D; ++i )
      out << center[i] << ' ';
    out << "]\n";
  }

  // Setter(s) ============================================================== //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RF<D,C>::setDefault()
  {
    radius = (C) 1;
    center.fill( (C)0 );
  }

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RF<D,C>::setParameters( const coord_t * ) {};

  // ======================================================================== //
  // IMPLEMENTATION OF CLASS RFP                                              //
  // ======================================================================== //

  // Constructor(s) ========================================================= //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  RFP<D,N,C>::RFP() :
    base_t()
  {}

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  template<class ...Args>
  RFP<D,N,C>::RFP( bitpit::rbf::eRBFType type, coord_t(*f)( coord_t, Args ...args) ) :
    base_t( type, bindParameters(f, mParams) )
  { }

  // Getter(s)/Info ========================================================= //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  std::size_t RFP<D,N,C>::getNumberOfParameters() const { return N; };

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  const typename RFP<D,N,C>::coord_t* RFP<D,N,C>::getParameters() const { return mParams; }

  // Setter(s) ============================================================== //

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  template<class f_in_t, std::size_t... I>
  typename RFP<D,N,C>::rf_funct_t RFP<D,N,C>::doBind( f_in_t f, coord_t* data, bitpit::index_sequence<I...> )
  {
       return std::bind( f, std::placeholders::_1, std::cref(data[I]...) ); // A trick here
  }

  // ------------------------------------------------------------------------ //
  template< std::size_t D, std::size_t N, class C >
  template<class f_in_t>
  typename RFP<D,N,C>::rf_funct_t RFP<D,N,C>::bindParameters( f_in_t f, coord_t* data )
  {
       return doBind( f, data, bitpit::make_index_sequence<N>{} );
  }

  // ---------------------------------------------------------------------- //
  template< std::size_t D, std::size_t N, class C >
  void RFP<D,N,C>::setParameters( const coord_t *values )
  {
    std::copy( values, values + N, mParams );
  }

  // ---------------------------------------------------------------------- //
  template< std::size_t D, std::size_t N, class C >
  void RFP<D,N,C>::display( std::ostream &out /*= std::cout*/, unsigned int indent /*= 0*/ ) const
  {
    std::string   s( indent, ' ' );
    base_t::display(out, indent);
    out << s  << "# params:    " << N << '\n'
        << s  << "params:      [ ";
    for ( std::size_t i = 0; i < N; ++i )
      out << mParams[i] << ' ';
    out << "]\n";
  }

} //end namespace rbf
} //end namespace bitpit
