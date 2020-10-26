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

#pragma once

// Standard Template Library
# include <vector>
# include <array>
# include <functional>
# include <exception>
# include <type_traits>

// Bitpit
# include <bitpit_operators.hpp>


namespace bitpit{

/*!
 * @enum RBFBasisFunction
 * @ingroup RBF
 * @brief Enum class defining types of RBF kernel functions that could be used in bitpit::RBF class
 */
enum class RBFBasisFunction {
    CUSTOM     = 0,  /**<Identify custom linked support function */
    WENDLANDC2 = 1,  /**< Compact support Wendland C2 function */
    LINEAR     = 2,  /**< Compact support linear function */
    GAUSS90    = 3,  /**< Non compact gaussian with 90% of reduction at unary radius */
    GAUSS95    = 4,  /**< Non compact gaussian with 95% of reduction at unary radius */
    GAUSS99    = 5,  /**< Non compact gaussian with 99% of reduction at unary radius */
    C1C0       = 6,  /**< Compact quadratic funct, C1 on r=0, C0 on r=1, 0 outside */
    C2C0       = 7,  /**< Compact cubic funct, C2 on r=0, C0 on r=1, 0 outside */
    C0C1       = 8,  /**< Compact quadratic funct, C0 on r=0, C1 on r=1, 0 outside */
    C1C1       = 9,  /**< Compact cubic funct, C1 on r=0, C1 on r=1, 0 outside */
    C2C1       = 10, /**< Compact biquadratic funct, C2 on r=0, C1 on r=1, 0 outside */
    C0C2       = 11, /**< Compact cubic funct, C0 on r=0, C2 on r=1, 0 outside */
    C1C2       = 12, /**< Compact biquadratic funct, C1 on r=0, C2 on r=1, 0 outside */
    C2C2       = 13, /**< Compact poly (degree 5) funct, C2 on r=0, C2 on r=1, 0 outside */
};

/*!
 * @enum RBFMode
 * @ingroup RBF
 * @brief Enum class defining behaviour of the bitpit::RBF class
 */
enum class RBFMode {
    INTERP = 1, /**< RBF class interpolate external field data */
    PARAM  = 2  /**< RBF class used as pure parameterizator*/
};

class RBFKernel{

private:
    int     m_fields;                               /**<Number of data fields defined on RBF nodes.*/
    RBFMode m_mode;                                 /**<Behaviour of RBF class (interpolation or parametrization).*/
    double  m_supportRadius;                        /**<Support radius of function used as Radiabl Basis Function.*/
    RBFBasisFunction m_typef;                       /**<Recognize type of RBF shape function actually in the class. */
    double  (*m_fPtr)(double);

    std::vector<double>                 m_error;    /**<Interpolation error of a field evaluated on each RBF node (auxiliary memeber used in Greedy algorithm).*/

    protected:
    std::vector<std::vector<double>>    m_value;    /**< displ value to be interpolated on RBF nodes */
    std::vector<std::vector<double>>    m_weight;   /**< weight of your RBF interpolation */
    std::vector<bool>                   m_activeNodes;   /**<Vector of active/inactive node (m_activeNodes[i] = true/false -> the i-th node is used/not used during RBF evaluation).*/
    int m_maxFields;                                /**< fix the maximum number of fields that can be added to your class*/
    int m_nodes;                                    /**<Number of RBF nodes.*/

public:
    ~RBFKernel();
    RBFKernel();
    RBFKernel(const RBFKernel & other);

    void                    setFunction(RBFBasisFunction);
    void                    setFunction(double (&funct)(double ));

    RBFBasisFunction        getFunctionType();
    int                     getDataCount();
    int                     getActiveCount();
    std::vector<int>        getActiveSet();

    bool                    isActive(int );

    bool                    activateNode(int );
    bool                    activateNode(const std::vector<int> &);
    void                    activateAllNodes();
    bool                    deactivateNode(int );
    bool                    deactivateNode(const std::vector<int> &);
    void                    deactivateAllNodes();

    void                    setSupportRadius(double);
    double                  getSupportRadius();

    void                    setMode(RBFMode mode);
    RBFMode                 getMode();

    void                    setDataToNode (int , const std::vector<double> &);
    void                    setDataToAllNodes(int , const std::vector<double> &);

    int                     addData();
    int                     addData(const std::vector<double> &);
    bool                    removeData(int);
    bool                    removeData(std::vector<int> &);
    void                    removeAllData();

    void                    fitDataToNodes();
    void                    fitDataToNodes(int);

    std::vector<double>     evalRBF(const std::array<double,3> &);
    std::vector<double>     evalRBF(int jnode);
    double                  evalBasis(double);

    int                     solve();
    int                     greedy(double);

protected:
    double                  evalError();
    int                     addGreedyPoint();
    int                     solveLSQ();
    void                    swap(RBFKernel & x) noexcept;

private:

    virtual double calcDist(int i, int j) = 0;
    virtual double calcDist(const std::array<double,3> & point, int j) = 0;

};

class RBF : public RBFKernel {

protected:
    std::vector<std::array<double,3>>   m_node;     /**< list of RBF nodes */

public:
    ~RBF();
    RBF(RBFBasisFunction = RBFBasisFunction::WENDLANDC2);
    RBF(const RBF & other);
    RBF & operator=(RBF other);

    int                     getTotalNodesCount();

    int                     addNode(const std::array<double,3> &);
    std::vector<int>        addNode(const std::vector<std::array<double,3>> &);
    bool                    removeNode(int);
    bool                    removeNode(std::vector<int> &);
    void                    removeAllNodes();

protected:
    void     swap(RBF & x) noexcept;

private:
    double calcDist(int i, int j);
    double calcDist(const std::array<double,3> & point, int j);
};

/*!
 * @ingroup  RBF
 * @brief Utility fuctions for RBF
 */
namespace rbf
{
    double                  wendlandc2(double);
    double                  linear(double);
    double                  gauss90(double);
    double                  gauss95(double);
    double                  gauss99(double);
    double                  c1c0(double);
    double                  c2c0(double);
    double                  c0c1(double);
    double                  c1c1(double);
    double                  c2c1(double);
    double                  c0c2(double);
    double                  c1c2(double);
    double                  c2c2(double);

  /*! @addtogroup RBF
   *  @{
  */

  /*! @brief C2-continuous Wendland'2 functions.
   *
   *  @tparam         CoordT      type of coeffs. (e.g. double, float, etc.)
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr >
  CoordT  wendland_c2( CoordT r )
  {
    return r > (CoordT)1 ?
      (CoordT)0 :
      std::pow( (CoordT)1 - r, 4 ) * ( (CoordT)4 * r + (CoordT)1 );
  }

  /*! @brief Generalized multiquadric functions.
   *  @tparam         CoordT      type of coeffs. (e.g. double, float, etc. )
   *  @tparam         Alha, Beta  expoenent values.
   *
   *  @param [in]     r           radial distance
   *  @param [in]     c           value for the bias coeff.
  */
  template<
    unsigned Alpha,
    unsigned Beta,
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr >
  CoordT generalized_multiquadrics( CoordT r, CoordT c )
  {
    return CoordT(1)/std::pow( c*c + r*r, CoordT(Alpha)/CoordT(Beta) );
  }

  // ================================================================ //
  // FORWARD DECLARATIONS                                             //
  // ================================================================ //
  template< size_t, class >
  class WendlandC2;
  template< size_t, unsigned, unsigned, class >
  class GeneralizedMultiquadrics;


  /*! @brief Enum for supported families of radial functions. */
  enum eRBFType
  {
    /*! @brief undefined type. */
    kUndefined                = 0,
    /*! @brief WendLand C2-continuous radial functions. */
    kWendlandC2               = 10,
    /*! @brief Generalized multiquadrics */
    kGeneralizedMultiQuadric  = 20,
    /*! @brief Generalized multiquadrics with \f$\alpha = 2, \beta = 1\f$) */
    kMultiQuadric2            = 21,
    /*! @brief Generalized multiquarics with \f$\alpha = 3, \beta = 2\f$ */
    kMultiQuadric3_2          = 22,
    /*! @brief Generalized multiquarics with \f$\alpha = 5, \beta = 2\f$ */
    kMultiQuadric5_2          = 23

  }; //end enum eRBFType

  // ================================================================ //
  // DEFINITION OF CLASS RadialFunct									                //
  // ================================================================ //
  /*! @brief Base class used to derive radial functions of different families.
   *
   *	@tparam 			Dim 			   nr. of dimension in the working space.
   *	@tparam 			CoordT 	     (default = double) type of coeffs., coordinates, e.g. coord_t = double, float.
   *									           (only scalar floating point type are allowed).
  */
  template<
    size_t 	Dim,
    class 	CoordT = double
  >
  class RadialFunct
  {
    // Static assertions ========================================== //
    static_assert(
      (Dim > 0),
      "**ERROR** RadialFunct<Dim,CoordT>: nr. of working dimensions, Dim, "
      "must be greater than 0"
    );
    static_assert(
      std::is_floating_point<CoordT>::value,
      "**ERROR** RadialFunct<Dim,CoordT>: CoordT must be a integral floating point type "
      "(i.e. float, double or long double)"
    );

    // Typedef(s) ================================================= //
    public:
    /*!	@brief Coeffs. type. */
    using coord_t 		= CoordT;
    /*!	@brief Point type in the working space. */
    using point_t 		= std::array<coord_t, Dim>;
    /*!	@brief Type of functor holding the actual implementation */
    using rf_funct_t  = std::function< coord_t( coord_t ) >;
    private:
    /*!	@brief Type of this object. */
    using self_t    	= RadialFunct<Dim, CoordT>;

    // Member variable(s) ========================================= //
    protected:
    /*!	@brief Function implementing the expression of the Radial Basis Function. */
    rf_funct_t        mFunct;
    /*! @brief Type of this radial function. */
    eRBFType          mType;
    public:
    /*!	@brief Radius of this Radial function. */
    coord_t           radius;
    /*!	@brief Center of the Radial Function (sometimes referred to as control point) */
    point_t           center;

    // Static member function(s) ================================== //
    template< class ... Args >
    static RadialFunct*		New( eRBFType type, Args &&...args )
    {
      switch( type )
      {
        default:
        {
          throw std::runtime_error(
            "rbf::RadialFunct::new undefined radial function type"
          );
        } //end default
        case( bitpit::rbf::eRBFType::kWendlandC2 ):
        {
          return new bitpit::rbf::WendlandC2<Dim, CoordT>( args... );
        } //end case kWendlandC2
        case( bitpit::rbf::eRBFType::kMultiQuadric2 ):
        {
          return new bitpit::rbf::GeneralizedMultiquadrics<Dim, 2, 1, CoordT>( args ... );
        } //end case kMultiQuadrics2
        case( bitpit::rbf::eRBFType::kMultiQuadric3_2 ):
        {
          return new bitpit::rbf::GeneralizedMultiquadrics<Dim, 3, 2, CoordT>( args ... );
        }
        case( bitpit::rbf::eRBFType::kMultiQuadric5_2 ):
        {
          return new bitpit::rbf::GeneralizedMultiquadrics<Dim, 5, 2, CoordT>( args ... );
        }
      } //end switch

      return nullptr;
    }

    // Member function(s) ========================================= //

    // Constructor(s) --------------------------------------------- //
    public:
    /*! @brief Default constructor.
     *
     *  Initialize a radial function of undefined type with default values
     *  for the parameters (see #setDefault )
    */
    RadialFunct() :
      mFunct(nullptr)
      , radius()
      , center()
      , mType( bitpit::rbf::eRBFType::kUndefined )
    { }
    /*! @brief Constructor #1.
     *
     *  Initialize a radial function with the specified radius and center.
     *
     *  @param [in]     r         radius of this radial function.
     *  @param [in]     c         center of this radial function.
    */
    RadialFunct( coord_t r, const point_t &c ) :
      mFunct(nullptr)
      , radius(r)
      , center(c)
      , mType( bitpit::rbf::eRBFType::kUndefined )
    { }

    /*! @brief Copy-constructor.
     *
     *  Initialize this by object by copy.
     *
     *  @param [in]     other     input object.
    */
    RadialFunct( const self_t &other ) :
      mFunct(nullptr)
      , radius( other.radius )
      , center( other.center )
      , mType( other.mType )
    { }
    /*! @brief Move-constructor.
     *
     *  Initialize this ojbect by moving the content of the input
     *  class.
     *
     *  @note Depending on the implementation of the derived class,
     *  the input object might be left in a undefined state after its
     *  content is moved into this.
    */
    RadialFunct( self_t &&other ) :
      mFunct(nullptr)
      , radius( std::move( other.radius ) )
      , center( std::move( other.center ) )
      , mType( std::move( other.mType ) )
    { }

    // Operator(s) =============================================== //
    public:
    /*!	@brief Copy assignement operator.
     *
     *	Copy the content of the input object into this.
     *
     *	@param [in]			rf			input object.
     *
     *	@result Returns reference to this object.
    */
    self_t&					operator=( const self_t &other )
    {
      radius = other.radius;
      center = other.center;
      mType  = other.mType;
      return *this;
    }
    /*!	@brief Move assignment operator.
     *
     *	Move the content of the input object into this.
     *
     *	@note Depending on the impementation of the derived classes,
     *	after moving its content the input object might be in a undefined state.
     *
     *	@param [in]			other 		input object.
     *
     *	@result Returns reference to this.
    */
    self_t&					        operator=( self_t &&other )
    {
      radius = std::move( other.radius );
      center = std::move( other.center );
      mType  = std::move( other.mType );
      return *this;
    }
    /*!	@brief Evaluate this RF function at the input point.
     *
     *	@param [in]			coords 		point coordinates.
    */
    coord_t 				        operator()( const point_t &coords ) const
    {
      //using bitpit::norm_2;
      //using bitpit::operator-;
      (*mFunct)( norm_2( coords - center )/radius );
    }

    // Getter(s)/Info --------------------------------------------- //
    public:
    /*!	@brief Returns the nr. of additional parameters of this RF.
     *
     *  The default behavior is to assume that the radial function
     *  does not depend on any additional parameter.
     *
     *  @note If a family of radial basis functions depends on additional
     *  parameters, the derived class implementing the family should
     *  override this method.
    */
    virtual size_t 			    getNumberOfParameters() const { return 0; };
    /*!	@brief Returns (true) if this RF has a compact support.
     *
     *  The default behavior is to assume that the radial function
     *  has compact support.
     *
     *  @note If a family of radial basis function does not have a compact support,
     *  the derived class implementing the radial function, should override this method.
    */
    virtual bool 			      hasCompactSupport() const { return true; };
    /*!	@brief Returns const pointer to the internal array of parameters
     *	for this Radial Function.
     *
     *  The default behavior is to assume that the radial function does not
     *  depend on any additional parameter.
     *
     *	@note If a family of radial basis function depends on additional parameters,
     *  the defrived class implementing the radial function, should overload this method.
    */
    virtual const coord_t*	getParameters() const { return nullptr; };
    /*! @brief Returns the type of this radial function. */
    eRBFType                getType() const { return mType; }

    // Setter(s) -------------------------------------------------- //
    public:
    /*! @brief Set default values for function parameters, ie.:
     *  * radius = 1          (assumes uniform behavior across the basis, and normalized space)
     *  * center = (0, 0, 0)  (radial function centered in the origin)
    */
    virtual void            setDefault()
    {
      radius = (CoordT)1;
      center.fill( (CoordT)0 );
    }
    /*!	@brief Set the value of the parameters of this Radial Function.
     *
     *  The default behavior is to assume that the radial function does not
     *  depend on any additional parameters.
     *
     *	@note If a family of radial basis function depends on additional parameters,
     *  the defrived class implementing the radial function, should overload this method.
    */
    virtual void            setParameters( const coord_t * ) {};

  }; //end class RadialFunct

  // ============================================================= //
  // DEFINITION OF CLASS WendlandC2                                //
  // ============================================================= //
  /*! @brief Definition of WendLand C2-continuous radial function.
   *
   *  C2-continuous Wendland's functions have the following expression:
   *
   *  \f$
   *  \begin{equation}
   *      \left\{
   *        \begin{aligned}
   *          &(1 - z )^4 ( 4 z + 1 ), \; \text{if} z < 0, \\
   *          &0, \; \text{otherwise}
   *        \end{aligned}
   *    \right.
   *  \end{equation}
   *  \f$
   *
   *  where \f$ z := \frac{\left| x - c \right|}{r}\f$ is the radial distance
   *  between a given point, \f$x\f$ and the center of the radial function, \f$c\f$
   *  rescaled by the radius of the function, \f$r\f$.
   *
   *  @tparam       Dim       nr . of dimensions in the working space.
   *  @tparam       CoordT    (default = double) type of coeffs. e.g. double, float, etc.
   *                          (only scalar floating point type are allowed).
  */
  template< size_t Dim, class CoordT = double >
  class WendlandC2 : public RadialFunct<Dim, CoordT>
  {
    // Typedef(s) ================================================ //
    private:
    /*! @brief Type of base class. */
    using base_t      = RadialFunct<Dim, CoordT>;
    /*! @Type of this class. */
    using self_t      = WendlandC2<Dim, CoordT>;
    public:
    /*! @brief Coeffs. type. */
    using coord_t     = typename base_t::coord_t;
    /*! @brief Point type. */
    using point_t     = typename base_t::point_t;

    // Member variable(s) ======================================== //
    protected:
    using base_t::mFunct;
    using base_t::mType;
    public:
    using base_t::radius;
    using base_t::center;

    // Constructor(s) ============================================ //
    public:
    /*! @brief Default constructor.
     *
     *  Initialize a  WendLand C2 function with default values for
     *  for function parameters.
    */
    WendlandC2() :
      base_t()
    {
      mFunct  = &bitpit::rbf::wendland_c2<coord_t>;
      mType   = bitpit::rbf::eRBFType::kWendlandC2;
    }
    /*! @brief Constructor #1.
     *
     *  Initialize a Wendland C2 function with the specified radius and center.
     *
     *  @param [in]     r       radius for this radial function.
     *  @param [in]     c       center for this radial functions.
    */
    WendlandC2( coord_t r, const point_t &c ) :
      base_t( r, c )
    {}

    /*! @brief Copy-constructor (defaulted). */
    WendlandC2( const self_t & ) = default;
    /*! @brief Move-constructor (defaulted). */
    WendlandC2( self_t && ) = default;

    // Operator(s) =============================================== //
    public:
    /*! @brief Copy-assignement operator (defaulted). */
    self_t&       operator=( const self_t & ) = default;
    /*! @brief Move-assignement operator (defaulted). */
    self_t&       operator=( self_t && ) = default;

    // Member function(s) ======================================== //

    // Getter(s)/Info ------------------------------------------- //
    public:
    // Wendland C2 radial functions do not have any additional parameters
    // and are compactly supported. The default behavior implemented in
    // the base class should suffice.

    // Setter(s) ------------------------------------------------- //
    public:
    // Wendland C2 radial functions do not have any additional parameters
    // and are compactly supported. The default behavior implemented in
    // the base class should suffice.

  }; //end class WendlandC2

  // ============================================================= //
  // DEFINITION OF CLASS WendlandC2                                //
  // ============================================================= //
  /*! @brief Definitino of generalized multiquadrics radial functinos.
   *
   *  Generalized multi-quadratics radial function have the following expression:
   *
   *  \f$
   *  \begin{equation}
   *    \frac{1}{ \left( c^2 + z^2 \right)^{\frac{\alpha}{\beta}} }
   *  \end{equation}
   *  \f$
   *
   *  for \f$\alpha, \beta > 0\f$.
   *  In the above equation, \f$z := \frac{ \left|x - c\right|}{r}\f$ is the radial
   *  distance of a given point, \f$x\f$ from the center of the radial function, \f$c\f$
   *  with \f$r\f$ being the radius of the function.
   *
   *  @tparam       Dim           nr. of dimensions in the working space (Dim>0)
   *  @tparam       Alpha, Beta   (>0) fractional expoents.
   *  @tparam       CoordT        (default = double) type of coefficients, e.g. double, float, etc.
   *                              (only scalar floating point types are allowed).
  */
  template< size_t Dim, unsigned Alpha, unsigned Beta, class CoordT = double >
  class GeneralizedMultiquadrics : public RadialFunct<Dim, CoordT>
  {
    // Static assertions ========================================= //
    static_assert(
      ( (Alpha > 0) && (Beta > 0) ),
      "**ERROR bitpit::rbf::GeneralizedMultiquadrics<Dim, Alpha, Beta, CoordT>: "
      "exponents Alpha and Beta must be positve"
    );

    // Typedef(s) ================================================ //
    private:
    /*! @brief Type of the base class. */
    using base_t    = RadialFunct<Dim, CoordT>;
    /*! @brief Type of this class. */
    using self_t    = GeneralizedMultiquadrics<Dim, Alpha, Beta, CoordT >;
    public:
    /*! @brief Coeff. type. */
    using coord_t   = typename base_t::coord_t;
    /*! @brief Point type. */
    using point_t   = typename base_t::point_t;

    // Member variable(s) ======================================== //
    protected:
    using base_t::mType;
    using base_t::mFunct;
    /*! @brief List of values for the additional parameters of this radial function. */
    coord_t     mParams[1];
    public:
    using base_t::center;
    using base_t::radius;

    // Typedef(s) ================================================ //
    public:
    /*! @brief Default constructor.
     *
     *  Initialize a generalized multiquadric radial function with default values
     *  for the parameters.
    */
    GeneralizedMultiquadrics() :
      base_t()
    {
      setDefault();
    }
    /*! @brief Contructor #1.
     *
     *  Initialize a generalized multiquadrics radial function with the specified
     *  values for the parameters.
     *
     *  @param [in]       r         radius of this function.
     *  @param [in]       c         center of this radial function.
     *  @param [in]       bias      bias coeff.
    */
    GeneralizedMultiquadrics( coord_t r, const point_t &c, coord_t bias = (coord_t) 0 ) :
      base_t( r, c )
    {
      mType = bitpit::rbf::eRBFType::kGeneralizedMultiQuadric;
      bindParameters();
      setDefault();
    }
    /*! @brief Copy constructor.
     *
     *  Initialize this radial function by copying the content of the input object.
     *
     *  @param [in]     other     object to be copied into this.
    */
    GeneralizedMultiquadrics( const self_t &other ) :
      base_t( other )
    {
      bindParameters();
      mParams[0] = other.mParams[0];
    }
    /*! @brief Move-constructor.
     *
     *  Initalize this radial function by moving the content of the input object.
     *
     *  @note After moving its content the input object might be in a undefined state,
     *
     *  @param [in]       other     object to be moved into this.
    */
    GeneralizedMultiquadrics( self_t &&other ) :
      base_t( std::move( other ) )
    {
      bindParameters();
      mParams[0] = other.mParams[0];
    }

    // Operator(s) =============================================== //
    public:
    /*! @brief Copy-assignment operator.
     *
     *  Copy the content of the input object into this.
     *
     *  @param [in]       other     object to be copied into this.
     *
     *  @result Returns reference to this object.
    */
    self_t&       operator=( const self_t &other )
    {
      base_t::operator=( other );
      mParams[0] = other.mParams[0];
      return *this;
    }
    /*! @brief Move-assignement operator.
     *
     *  Move the content of the input object into this.
     *
     *  @note After moving its content, the input object might be in a undefined
     *  state.
     *
     *  @param [in]       other       object to be moved into this.
     *
     *  @result Returns refernce to this object.
    */
    self_t&       operator=( self_t &&other )
    {
      base_t::operator=( std::move(other) );
      mParams[0] = other.mParams[0];
      return *this;
    }

    // Member function(s) ======================================== //

    // Getter(s)/Info -------------------------------------------- //
    public:

    /*! @brief Returns the nr. of additional parameters for this radial function. */
    virtual size_t      getNumberOfParameters() const override { return 1; };
    /*! @brief Returns (const) pointer to the list of additional parameters of this funciton. */
    virtual coord_t*    getParameters(){ return mParams; }
    /*! @brief Returns (false) since this function does not have a compact support. */
    virtual bool        hasCompactSupport() const override { return false; }

    // Setter(s) ------------------------------------------------- //
    private:
    /*! @brief Bind the expression of the radial function to the parameters of
     *  this instance.
    */
    void                bindParameters()
    {
      mFunct = std::bind(
        &bitpit::rbf::generalized_multiquadrics<Alpha,Beta,coord_t>,
        std::placeholders::_1,
        std::cref( mParams[0] )
      );
    }
    public:
    /*! @brief Set default values for the paraters, i.e.:
     *  * radius = 1
     *  * center = (0, 0, ..., 0)
     *  * c = 0
    */
    virtual void        setDefault() override
    {
      base_t::setDefault();
      mParams[0] = (coord_t) 0;
    }
    /*! @brief Set the value for the additional parameters of this function.*/
    virtual void        setParameters( const coord_t *values ) override
    {
      mParams[0] = values[0];
    }
  }; //end class GeneralizedMultiquadrics

  /*! @brief Typename for multiquadric RF with \f$\alpha = 2, \beta = 1\f$ */
  template< size_t Dim, class CoordT = double >
  using GeneralizedMultiquadrics2 = GeneralizedMultiquadrics<Dim, 2, 1, CoordT >;
  /*! @brief Typename for multiquadric RF with \f$\alpha = 3, \beta = 2\f$ */
  template< size_t Dim, class CoordT = double >
  using GeneralizedMultiquadrics3_2 = GeneralizedMultiquadrics<Dim, 3, 2, CoordT >;
  /*! @brief Typename for multiquadric RF with \f$\alpha = 5, \beta = 2\f$ */
  template< size_t Dim, class CoordT = double >
  using GeneralizedMultiquadrics5_2 = GeneralizedMultiquadrics<Dim, 5, 2, CoordT >;

  // ============================================================= //
  // EXPLICIT SPECIALIZATIONS                                      //
  // ============================================================= //
  extern template class WendlandC2<1, float>;
  extern template class WendlandC2<2, float>;
  extern template class WendlandC2<3, float>;
  extern template class WendlandC2<1, double>;
  extern template class WendlandC2<2, double>;
  extern template class WendlandC2<3, double>;
  extern template class WendlandC2<1, long double>;
  extern template class WendlandC2<2, long double>;
  extern template class WendlandC2<3, long double>;

  extern template class GeneralizedMultiquadrics<1, 2, 1, float>;
  extern template class GeneralizedMultiquadrics<2, 2, 1, float>;
  extern template class GeneralizedMultiquadrics<3, 2, 1, float>;
  extern template class GeneralizedMultiquadrics<1, 2, 1, double>;
  extern template class GeneralizedMultiquadrics<2, 2, 1, double>;
  extern template class GeneralizedMultiquadrics<3, 2, 1, double>;
  extern template class GeneralizedMultiquadrics<1, 2, 1, long double>;
  extern template class GeneralizedMultiquadrics<2, 2, 1, long double>;
  extern template class GeneralizedMultiquadrics<3, 2, 1, long double>;

  extern template class GeneralizedMultiquadrics<1, 3, 2, float>;
  extern template class GeneralizedMultiquadrics<2, 3, 2, float>;
  extern template class GeneralizedMultiquadrics<3, 3, 2, float>;
  extern template class GeneralizedMultiquadrics<1, 3, 2, double>;
  extern template class GeneralizedMultiquadrics<2, 3, 2, double>;
  extern template class GeneralizedMultiquadrics<3, 3, 2, double>;
  extern template class GeneralizedMultiquadrics<1, 3, 2, long double>;
  extern template class GeneralizedMultiquadrics<2, 3, 2, long double>;
  extern template class GeneralizedMultiquadrics<3, 3, 2, long double>;

  extern template class GeneralizedMultiquadrics<1, 5, 2, float>;
  extern template class GeneralizedMultiquadrics<2, 5, 2, float>;
  extern template class GeneralizedMultiquadrics<3, 5, 2, float>;
  extern template class GeneralizedMultiquadrics<1, 5, 2, double>;
  extern template class GeneralizedMultiquadrics<2, 5, 2, double>;
  extern template class GeneralizedMultiquadrics<3, 5, 2, double>;
  extern template class GeneralizedMultiquadrics<1, 5, 2, long double>;
  extern template class GeneralizedMultiquadrics<2, 5, 2, long double>;
  extern template class GeneralizedMultiquadrics<3, 5, 2, long double>;
  /*! @} */

} //end namespace rbf

} //end namespace bitpit
