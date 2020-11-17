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
  CoordT wendland_c2( CoordT r )
  {
    return r > (CoordT)1 ?
      (CoordT)0 :
      std::pow( (CoordT)1 - r, 4 ) * ( (CoordT)4 * r + (CoordT)1 );
  }

  // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT wendland_c2_der1( CoordT r )
  {
    return r > (CoordT)1 ?
      (CoordT)0 :
      - (CoordT)4 * std::pow( (CoordT)1 - r, 3 ) * ( (CoordT)4 * r + (CoordT)1 )
      + (CoordT)4 * std::pow( (CoordT)1 - r, 4 );
  }
  
    // ------------------------------------------------------------------------ //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* >
  CoordT wendland_c2_der2( CoordT r )
  {
    return r > (CoordT)1 ?
      (CoordT)0 :
      + (CoordT)12 * std::pow( (CoordT)1 - r, 2 ) * ( (CoordT)4 * r + (CoordT)1 )
      - (CoordT)32 * std::pow( (CoordT)1 - r, 3 );
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

  // ----------------------------------------------------------------------- //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr >
  CoordT gaussian_der1( CoordT r, CoordT c )
  {
    return - (CoordT)2 * c * r * std::exp( - (r*r) );
  }
  
  // ----------------------------------------------------------------------- //
  template< class CoordT, typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr >
  CoordT gaussian_der2( CoordT r, CoordT c )
  {
    return - (CoordT)2 * c * std::exp( - (r*r) );
           + (CoordT)4 * c * r * r * std::exp( - (r*r) );
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
  // HELPER FUNCTIONS                                                         //
  // ======================================================================== //
  
  // ------------------------------------------------------------------------ //
  template<
    std::size_t Dim,
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* /*= nullptr*/
  >
  bool computeRBFWeights( const std::vector<typename RFBasis<Dim,CoordT>::point_t> &data_points, const std::vector<CoordT> data_values, RFBasis<Dim,CoordT> &rbf )
  {
    if ( data_points.size() != data_values.size() )
      throw std::runtime_error(
        "bitpit::rbf::computeRBFWeights: The nr. of data points and the nr. of data values mismatch!"
      );
    if ( data_points.size() < rbf.size() )
      throw std::runtime_error(
        "bitpit::rbf::computeRBFWeights: The nr. of data points must be at least the size of the radial function basis."
      );

    // Constant(s)
    const CoordT  rcond = (CoordT)-1;

    // Scope variables
    int n_rbf   = rbf.size();
    int n_data  = data_points.size();

    // Implementation with LAPACKE support -------------------------- //
    // Implementation node: dgesv suffer from numerical stability issues.
    // Better implementation in Eigen with QR factorization and full pivoting.
    # ifdef __RBF_USE_LAPACKE__
    // Scope variables
    std::vector<CoordT> mat( n_data * n_rbf );
    std::vector<CoordT> rhs( n_data );

    // Fill rhs
    for( std::size_t i_data = 0; i_data < n_data; ++i_data )
      rhs[i_data] = data_values[i_data];
    //next j

    // Fill coeff. matrix
    for( std::size_t i_data = 0, i_mat = 0; i_data < n_data; ++i_data ) {
        const auto &p = data_points.at(i_data);
        for( std::size_t i_rbf = 0; i_rbf < n_rbf; ++i_rbf ) {
            const auto &rf = *( rbf.at(i_rbf).second );
            mat[i_mat++] = rf( p );
        } //next i
    } //next j

    // Single precision
    lapack_int info;
    if ( LAPACKE_xgels_type_wrap<CoordT>::xgels_ptr )
    {
      info = LAPACKE_xgels_type_wrap<CoordT>::xgels_ptr( 
        LAPACK_ROW_MAJOR, //matrix layout
        'N',              // 'M' for col major, 'N' for row major 
        n_data,           //nr. of matrix row
        n_rbf,            //nr. of matrix cols
        1,                //nr. of rhs cols.
        mat.data(),       //ptr to matrix entries
        n_rbf,            //leading dim for matrix array
        rhs.data(),       //ptr. to rhs entries
        1                 //leading dim for rhs array
      );
    }
    else
    {
      throw std::runtime_error(
        "bitpit::rbf::computeRBFWeights: ** ERROR ** LAPACK does not support the required floating point precision."
      );
    }
    
    // Check for error(s)
    if( info > 0 )
      return false;
    
    // Set weights (if success)
    rbf.setWeights(rhs.cbegin(), rhs.cbegin() + n_rbf );
    # endif

    // Implementation with EIGEN support ---------------------------- //
    # ifdef __RBF_USE_EIGEN__
    
    // Typedef(s)
    using matrix_t  = Eigen::Matrix<CoordT, Eigen::Dynamic, Eigen::Dynamic>;
    using vec_t     = Eigen::Matrix<CoordT, Eigen::Dynamic, 1>;
    
    // Fill coeff. matrix
    matrix_t  mat( n_data, n_rbf );
    vec_t     rhs( n_data );
    
     // Fill rhs
    for( std::size_t i_data = 0; i_data < n_data; ++i_data )
      rhs(i_data) = data_values[i_data];
    //next j

    // Fill coeff. matrix
    for( std::size_t i_data = 0; i_data < n_data; ++i_data ) {
        const auto &p = data_points.at(i_data);
        for( std::size_t i_rbf = 0; i_rbf < n_rbf; ++i_rbf ) {
            const auto &rf = *( rbf.at(i_rbf).second );
            mat(i_data, i_rbf) = rf( p );
        } //next i
    } //next j
      
    // Solve with QR fact. full pivoting
    vec_t sol = mat.colPivHouseholderQr().solve(rhs);
    
    // Assign weights to RBF
    rbf.setWeights( sol.data(), sol.data() + sol.size() );
    
    # endif

    
    // Return success
    return true;
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

  // ------------------------------------------------------------------------ //
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

  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RF<D,C>::rf_funct_t RF<D,C>::getFirstDerivative() const
    {
      switch( mType )
      {
        default :
        {
          throw std::runtime_error( 
            "bitpit::rbf::RF::getFirstDerivative: "
            "** ERROR ** First derivative is not available for this radial function."
          );
        }
        case( bitpit::rbf::eRBFType::kWendlandC2 ):
        {
          // No additional parameters. Nothing to bind
          return &bitpit::rbf::wendland_c2_der1<coord_t>;
        }
        case( bitpit::rbf::eRBFType::kGaussian ):
        {
          auto this_as_rfp = dynamic_cast< RFP<D,1,C>* >( const_cast< RF<D,C>* >( this ) );
          return RFP<D,1,C>::bindParameters( &bitpit::rbf::gaussian_der1<coord_t>, this_as_rfp->mParams );
        }
       }
       return nullptr;
    }
    
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RF<D,C>::rf_funct_t RF<D,C>::getSecondDerivative() const
  {
    switch( mType )
    {
      default:
        throw std::runtime_error(
          "bitpit::rbf::RF::getSecondDerivative: "
          "** ERROR ** Second derivative is not available for this radial function."
        );
      case( bitpit::rbf::eRBFType::kWendlandC2 ):
      {
        // No additional parameters. Nothing to bind.
        return &bitpit::rbf::wendland_c2_der2<coord_t>;
      }
      case( bitpit::rbf::eRBFType::kGaussian ):
      {
        auto this_as_rfp = dynamic_cast< RFP<D,1,C>* >( const_cast< RF<D,C>* >( this ) );
        return RFP<D,1,C>::bindParameters( &bitpit::rbf::gaussian_der2<coord_t>, this_as_rfp->mParams );
      }
    }
    return nullptr;
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


  // ======================================================================== //
  // IMPLEMENTATION OF CLASS RFBasis                                          //
  // ======================================================================== //
      
  // Constructor(s) ========================================================= //
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RFBasis<D,C>::RFBasis() :
    base_t()
  {}
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  RFBasis<D,C>::RFBasis( std::size_t n, eRBFType type ) :
    base_t( n )
  {
    for ( auto &rf : *this )
      rf.second.reset( rf_t::New( type ) ); 
  }
  
  // Operator(s) ============================================================ //
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RFBasis<D,C>::coord_t RFBasis<D,C>::operator()( const point_t &coords ) const
  {
    coord_t out = (coord_t)0;
    for ( const auto &rf : *this )
      out += rf.first * rf.second->operator()(coords);
    return out;
  }
  
  // Getter(s)/Info ========================================================= //
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  std::vector<typename RFBasis<D,C>::coord_t> RFBasis<D,C>::collectWeights() const
  {
    std::vector<coord_t> out( base_t::size() );
    auto  i = out.begin(), 
          e = out.end();
    auto  j = base_t::cbegin();
    for ( ; i != e; ++i, ++j )
      (*i) = j->first;
    return out;
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  const typename RFBasis<D,C>::coord_t& RFBasis<D,C>::getWeight( std::size_t i ) const
  {
    return this->at(i).first;
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RFBasis<D,C>::coord_t& RFBasis<D,C>::getWeight( std::size_t i )
  {
    return const_cast< coord_t& >( const_cast<const self_t*>( this )->getWeight(i) );
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  const typename RFBasis<D,C>::rf_t& RFBasis<D,C>::getRadialFunction( std::size_t i ) const
  {
    return *( this->at(i).second );
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  typename RFBasis<D,C>::rf_t& RFBasis<D,C>::getRadialFunction( std::size_t i )
  {
    return const_cast< rf_t& >( const_cast< const self_t* >( this )->getRadialFunction(i) );
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RFBasis<D,C>::display( std::ostream &out /*= std::cout*/, unsigned int indent /*= 0*/ ) const
  {
    std::string s(indent, ' ');
    out << s << "# of RBF:   " << this->size() << '\n';
    std::size_t i = 0;
    for ( const auto &rf : *this )
    {
      out << s << "  #" << i++ << '\n'
          << s << "  weight: " << rf.first << '\n';
      rf.second->display( out, indent +2 );
    } //next rf
  }
  
  // Setter(s) ============================================================== //
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  std::size_t RFBasis<D,C>::add( std::unique_ptr<rf_t> rf, C weight /*= (C)1*/ )
  {
    base_t::push_back( std::make_pair( weight, std::move(rf) ) );
    return base_t::size()-1;
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RFBasis<D,C>::remove( std::size_t i )
  {
    base_t::erase( base_t::begin() + i );
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RFBasis<D,C>::setWeights( const std::vector<coord_t> &weights )
  {
    setWeights( weights.cbegin(), weights.cend() );
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  template< class IteratorType >
  void RFBasis<D,C>::setWeights( IteratorType first, IteratorType last )
  {
    if ( std::distance( first, last ) != base_t::size() )
      throw std::runtime_error(
        "bitpit::rbf::RFBasis::setWeights: ** ERROR** The size of the input range "
        "and the size of this basis mismatch!"          
      );
    auto j = base_t::begin();
    for ( ; first != last; ++first, ++j )
      j->first = (*first);
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RFBasis<D,C>::setRadius( coord_t radius )
  {
    for ( auto &rf : *this )
      rf.second->radius = radius;
  }
  
  // ------------------------------------------------------------------------ //
  template< std::size_t D, class C >
  void RFBasis<D,C>::reset( eRBFType type, std::size_t n /*= -1*/ )
  {
    base_t( n == -1 ? this->size() : n ).swap( *this );
    for ( auto &rf : *this )
      rf.second.reset( rf_t::New( type ) ); 
  }
} //end namespace rbf
} //end namespace bitpit
