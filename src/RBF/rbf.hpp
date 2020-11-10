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
# pragma once
# define __USE_DEPRECATED__
# define __RBF_USE_EIGEN__
# ifdef __RBF_USE_EIGEN__
# undef __RBF_USE_LAPACKE__
# else
# define __RBF_USE_LAPACKE__
# endif

// ========================================================================= //
// INCLUDES                                                                  //
// ========================================================================= //

// Standard Template Library
# include <cmath>
# include <vector>
# include <array>
# include <functional>
# include <exception>
# include <type_traits>
# include <memory>
# include <string>
# include <iostream>
# include <stdexcept>

// Bitpit
# include "bitpit_operators.hpp"
# include "metaprogramming.hpp"

// Eigen
# ifdef __RBF_USE_EIGEN__
# include <Eigen/Dense>
# endif

// Lapack
# ifdef __RBF_USE_LAPACKE__
# include "bitpit_private_lapacke.hpp"
# endif


namespace bitpit{

# ifdef __USE_DEPRECATED__
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

# endif
/*!
 * @ingroup  	RBF
 * @brief 		Classes and utility fuctions for R(adial)B(asis)F(unction)
 */
namespace rbf
{
# ifdef __USE_DEPRECATED__

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
# endif

  /*! @addtogroup RBF
   *  @{
  */

  // ======================================================================= //
  // LIBRARY OF SUPPORTED RADIAL FUNCTIONS                                   //
  // ======================================================================= //
  
  // ----------------------------------------------------------------------- //
  /*! @brief C2-continuous Wendland'2 functions.
   *
   *  C2-continuous Wendland's functions have the following expression up to dimension 3:
   *  \f$
   *  \begin{equation}
   *      f(r;c) := \left\{
   *        \begin{aligned}
   *          &(1 - z )^4 ( 4 z + 1 ), \; \text{if} z < 0, \\
   *          &0, \; \text{otherwise}
   *        \end{aligned}
   *    \right.
   *  \end{equation}
   *  \f$
   *  In the above expression, \f$r\f$ is the radial distance of a given point
   *  from the geometrical kernel (e.g. the center of the RF).
   *
   *  @tparam         CoordT      Type for coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT  wendland_c2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Generalized multiquadric functions.
   *
   *  The family of generalized multi-quadrics RF has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r;c) := \frac{1}{ \left( c^2 + z^2 \right)^{\frac{\alpha}{\beta}} }
   *  \end{equation}
   *  \f$
   *  for \f$\alpha, \beta > 0\f$.
   *  In the above expression, \f$r\f$ is the radial distance of a given point
   *  from the kernel of the RF (e.g. the center of the RF),
   *  and \f$c\f$ is a additional parameter (bias).
   *
   *  @tparam         CoordT      type of coeffs. Only scalar floating 
   *                              point types are supported (e.g. double, float, etc. )
   *  @tparam         Alha, Beta  values for the fractional exponent.
   *
   *  @param [in]     r           radial distance
   *  @param [in]     c           value of the bias coeff.
  */
  template<
    class CoordT,
    unsigned Alpha,
    unsigned Beta,
    typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) && (Beta > 0) >::type* = nullptr
  >
  CoordT generalized_multiquadrics( CoordT r, CoordT c );

  // ---------------------------------------------------------------- //
  /*! @brief Gaussian radial function.
   *
   *  A Guassian radial basis function has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r;c) := c exp( -r^2 );
   *  \end{equation}
   *  \f$
   *  where \f$c\f$ is a amplification/reduction factor, and \f$r\f$ is the
   *  radial distance from the geometry kernel (e.g. the center of the radial function).
   *
   *  @tparam         CoordT      type of coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
   *  @param [in]     c           amplification/reduction factor.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT gaussian( CoordT r, CoordT c );
  
  // ---------------------------------------------------------------- //
  /*! @brief Linear function.
   *
   *  Linear function with trivial expression (generally used to as a bias):
   *  \f$ f(r) := r\ f$
   *  where \f$r\f$ is the radial distance from the geometry kernel (e.g. the 
   *  center of the radial function).
   *
   *  @tparam         CoordT      type of coeffs.  Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT linear( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Generalized power function (also referred to as radial power)
   *
   *  Generalized power functions have the following expression:
   *  \f$ f(r) := r^\alpha \f$
   *  where \f$\alpha>0\f$ is the power exponent and \f$r\f$ is the radial distance
   *  from the geomtry kernel of the radial function (e.g. the center of the RF).
   *
   *  @tparam         CoordT      type of coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. )
   *  @tparam         Alpha       power expoenent.
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    int   Alpha,
    typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) >::type* = nullptr
  >
  CoordT generalized_power( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Thin plate splines
   *
   *  Thin plate spline has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r):= r^2 log(r)
   *  \end{equation}
   *  \f$
   *  where \f$r\f$ is the radial distance from the geometry kernel
   *  (e.g. the center of the RF).
   *
   *  @tparam         CoordT      type of coeffs. Only scalar floating point types are
   *                              supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT    thin_plate_spline( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Polyharmonic radial functions.
   *
   *  The family of polyharmonic functions has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r):= (-c)^{1+\beta) r^{2\beta} \, log(r), \; \beta \neq 0
   *  \end{equation}
   *  \f$
   *  where \f$r\f$ is the radial distance, \f$\beta \ge 0 \f$ is the power
   *  exponent, and \f$c\f$ is a amplification/reduction factor.
   *
   *  @note The special case \f$beta = 2\f$ correspnds to thin-plate-splines.
   *
   *  @tparam         CoordT      type of coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *  @tparam         Beta        power expoenent (Beta >= 0).
   *
   *  @param [in]     r           radial distance.
   *  @param [in]     c           amplification/reduction factor.
  */
  template<
    class CoordT,
    int   Beta,
    typename std::enable_if< std::is_floating_point<CoordT>::value && (Beta >= 0) >::type* = nullptr
  >
  CoordT polyharmonic( CoordT r, CoordT c );

  // ================================================================ //
  // FORWARD DECLARATIONS                                             //
  // ================================================================ //
  template< std::size_t, class >
  class RF;
  template< std::size_t, std::size_t, class >
  class RFP;
  template< std::size_t, class >
  class RFBasis;

  // ================================================================ //
  // DEFINITIONS                                                      //
  // ================================================================ //
  
  // ---------------------------------------------------------------- //
  /*! @brief Enum for supported families of radial functions. */
  enum eRBFType
  {
    /*! @brief undefined type. */
    // Leave it first to facilitate auto-loop through eRBFType
    kUndefined,
    /*! @brief WendLand C2-continuous radial functions. */
    kWendlandC2,
    /*! @brief Gaussian radial function. */
    kGaussian,
    /*! @brief Thin plate spline of order 1*/
    kThinPlateSpline,
    /*! @brief Polyharmonic of order 2*/
    kPolyharmonic2,
    /*! @brief Polyharmonic of order 3*/
    kPolyharmonic3,
    /*! @brief Polyharmonic of order 4*/
    kPolyharmonic4,
    /*! @brief Hardy's radial function (from the family of multiquadrics with \f$\alpha = 1, \beta = 2\f$)*/
    kHardy,
    /*! @brief Generalized multiquadrics with \f$\alpha = 2, \beta = 1\f$) */
    kMultiQuadric2,
    /*! @brief Generalized multiquarics with \f$\alpha = 3, \beta = 2\f$ */
    kMultiQuadric3_2,
    /*! @brief Generalized multiquarics with \f$\alpha = 5, \beta = 2\f$ */
    kMultiQuadric5_2,
    /*! @brief Linear */
    kLinear,
    /*! @brief Quadratic */
    kQuadratic,
    /*! @brief Cubic */
    kCubic,
    /*! @brief Quartic. */
    kQuartic,
    /*! @brief User-defined */
    // Leave it last to facilitate auto-loop through eRBFType.
    kUserDefined
  }; //end enum eRBFType

  // ================================================================ //
  // HELPER FUNCTIONS                                                 //
  // ================================================================ //
  
  // ---------------------------------------------------------------- //
  /*! @brief Helper function returning the tag associated to each type
   *  of RBF.
  */
  std::string   getRBFTag( eRBFType type );
  
  // ---------------------------------------------------------------- //
  /*! @brief Compute the weights of a basis of radial functions
   *  for a least square interpolation problem.
   *
   *  Given a set of scattered data in a Dim dimensional space, 
   *  this function computes the weights associated to each radial function
   *  which best fit the input data in a least square sense.
   *
   *  @tparam         Dim           nr. of woking dimensions.
   *  @tparam         CoordT        type of coordinates. Only scalar floating point types
   *                                are supported (e.g. float, double, etc. ).
   *
   *  @param [in]     data_points   scattered points in the Dim-dimensional space.
   *  @param [in]     data_values   value associated to each data point.
   *  @param [in,out] rbf           basis of radial functions.
   *
   *  @result Returns (true) on success.
  */
  # ifdef __RBF_USE_LAPACKE__
  template< class CoordT >
  struct LAPACKE_xgels_type_wrap
  {
    using xgels_signature_t = lapack_int (*)( lapack_int, char, lapack_int, lapack_int, lapack_int, CoordT*, lapack_int, CoordT*, lapack_int );
    static const xgels_signature_t xgels_ptr;
  }; //end struct LAPACKE_xgels_type_wrap
  template< class CoordT >
  const typename LAPACKE_xgels_type_wrap<CoordT>::xgels_signature_t LAPACKE_xgels_type_wrap<CoordT>::xgels_ptr = nullptr;
  template<>
  const LAPACKE_xgels_type_wrap<float>::xgels_signature_t LAPACKE_xgels_type_wrap<float>::xgels_ptr;
  template<>
  const LAPACKE_xgels_type_wrap<double>::xgels_signature_t LAPACKE_xgels_type_wrap<double>::xgels_ptr;
  # endif
  template<
    std::size_t Dim,
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  bool computeRBFWeights( const std::vector<typename RFBasis<Dim,CoordT>::point_t> &data_points, const std::vector<CoordT> data_values, RFBasis<Dim,CoordT> &rbf );
  
  // ================================================================ //
  // DEFINITION OF CLASS RadialFunct									                //
  // ================================================================ //
  /*! @brief Class used for radial function with no additional parameters.
   *
   *  This class holds a radial function (RF in short) which does not depend
   *  on any additional parameter beyond its radius and its center.
   *
   *	@tparam 			Dim 		    nr. of dimension in the working space.
   *	@tparam 			CoordT 	    (default = double) type for coefficients and coordinates.
   *                            Only scalar floating point type are supported (e.g. coord_t = double, float).
  */
  template<
    std::size_t Dim,
    class CoordT = double
  >
  class RF
  {
    // Static assertions ============================================ //
    static_assert(
      std::is_floating_point<CoordT>::value,
      "**ERROR** bitpit::rbf::RF<Dim,CoordT>: CoordT must be a scalar floating point type "
      ", e.g. float, double or long double"
    );

    // Typedef(s) =================================================== //
    public:
    /*!	@brief Coeffs. type. */
    using coord_t 		= CoordT;
    /*!	@brief Point type in the working space. */
    using point_t 		= std::array<coord_t, Dim>;
    /*!	@brief Type of functor holding the actual implementation */
    using rf_funct_t  	= std::function< coord_t( coord_t ) >;
    private:
    /*!	@brief Type of this object. */
    using self_t    	= RF<Dim, CoordT>;

    // Member variable(s) =========================================== //
    protected:
    /*!	@brief Functor implementing the expression of the Radial Function. */
    rf_funct_t mFunct;
    /*! @brief Type of this radial function. */
    eRBFType mType;
    /*! @brief Bool for compactly supported RFs.*/
    bool mHasCompactSupport;
    public:
    /*!	@brief Radius of this RF. */
    coord_t radius;
    /*!	@brief Center of this RF (sometimes referred to as control point, geometry kernel) */
    point_t center;

    // Static member function(s) ==================================== //
    public:
    /*! @brief Returns a pointer to a new instance of a radial function 
     *  of the specified type.
     *
     *  @param [in]     type        type of the RF.    */
    static self_t* New( eRBFType type );
    
    // Constructor(s) =============================================== //
    protected:
    /*! @brief Default constructor.
     *
     *  Initialize a radial function of undefined type.
    */
    RF();
    public:
    /*! @brief Constructor #1.
     *
     *  Constructs a radial function of the specified type using the
     *  expression provided as input argument.
     *
     *  @param [in]     type      type of this RF.
     *  @param [in]     f         pointer to the function implementing the expression
     *	 						              of this radial function.
    */
    RF( bitpit::rbf::eRBFType type, coord_t (*f)(coord_t) );
    /*! @brief Constructor #2.
     *
     *  Constructs a radial function of the specified type using the
     *  functor provided as input argument.
     *
     *  @param [in]     type      type of this RF.
     *  @param [in]     f         functor implementing the expression of this RF.
    */
    RF( bitpit::rbf::eRBFType type, rf_funct_t const &f );
    /*! @brief Copy-constructor (deleted) */
    RF( const self_t & ) = delete;
    /*! @brief Move-constructor (deleted) */
    RF( self_t && ) = delete;

    // Operator(s) ================================================= //
    public:
    /*!	@brief Copy assignement operator (deleted). */
    self_t& operator=( const self_t & ) = delete;
    /*!	@brief Move assignment operator (deleted). */
    self_t&	operator=( self_t &&other ) = delete;
    /*!	@brief Evaluation operator.
     *
     *  Evaluate this RF function at the input point.
     *
     *	@param [in]		coords 	   coordinates of the input point.
    */
    coord_t operator()( const point_t &coords ) const;
    
    // Getter(s)/Info ============================================== //
    public:
    /*!	@brief Returns the nr. of additional parameters for this RF.
     *
     *  By default, it is assumed that the radial function
     *  does not depend on any additional parameter.
     *
     *  If the RF depends on some other parameter, this method
     *  must be overridden.
    */
    virtual std::size_t	getNumberOfParameters() const;
    /*!	@brief Returns (true) if this RF is compactly supported. */
    bool hasCompactSupport() const;
    /*!	@brief Returns const pointer to the internal array of parameters
     *	of this RF.
     *
     *  The default behavior is to assume that the radial function does not
     *  depend on any additional parameter.
     *
     *  If the RF depends on additional parameters, this method must be overridden
     *  by the derived class.
    */
    virtual const coord_t* getParameters() const;
    /*! @brief Returns the type of this radial function. */
    eRBFType getType() const;
    /*! @brief Display info for this RBF (mostly meant for debugging purposes).
     *
     *  @param [in,out]   out       (default = std::cout) output stream.
     *  @param [in]       indent    (default = 0) indentation level.
    */
    virtual void display( std::ostream &out = std::cout, unsigned int indent = 0 ) const;
    
    // Setter(s) =================================================== //
    public:
    /*! @brief Set default values for function parameters, ie.:
     *  * radius = 1          (assumes uniform behavior across the basis, and normalized space)
     *  * center = (0, 0, 0)  (radial function centered in the origin)
    */
    virtual void setDefault();
    /*!	@brief Set the value of the parameters of this Radial Function.
     *
     *  The default behavior is to assume that the radial function does not
     *  depend on any additional parameters.
     *
     *  If the RF depends on additional parameters, this method must be overridden
     *  by the derived class.
    */
    virtual void setParameters( const coord_t * );

  }; //end class RF

  // =============================================================== //
  // DEFINITION OF CLASS RFP                                         //
  // =============================================================== //
  /*! @brief Definition of generalized RF.
   *
   *  @tparam       Dim           nr. of dimensions in the working space (Dim>0)
   *  @tparam       NParams       Nr of additional parameters for this RF (NParams>0)
   *  @tparam       CoordT        (default = double) type of coefficients. Only scalar floating
   *                              point types are supported (e.g. double, float, etc.)
   *
   *  This class stores a generic radial function which depends
   *  on a arbitrary nr. of parameters of type CoordT.
   *
   *  ** Usage **
   *  * Default usage *
   *  If you wish to use one of the radial function implemented in bitpit just create a new instance
   *  by calling:
   *  RF::New (specify the type of radial function you want to use).
   *
   *  * Passing the ownership of the additional parameters to RFP class. *
   *  If you wish to encapsulate a user-defined function (which depends on N parameters),
   *  you can invoke the RFP constructor passing the pointer to that function. The nr. of 
   *  additional parameters, which you function depends on, must be specified in the list of template arguments.
   *  For instance, assumung D=3, CoordT = double and a user-defined funtion which depends
   *  on 2 additional parameters:
   *
   *  CoordT = my_radial_funct( CoordT r, CoordT par1, CoordT par2 ) defined somewhere
   *  auto my_rf = new RFP<3, 2, double>(
   *    bitpit::rbf::eRBFType::kUserDefined,
   *    &my_radial_funct
   *  )
   *  will create a instance of the radial function (my_rf) which has exclusive ownership 
   *  over the additional parameters. The value of these parameters can be accesed/modified
   *  any time via #setParameters and #getParameters
   *
   *  In some cases user might want to keep ownership on some (all) of the additional parameters.
   *  You can keep the ownership of such parameters by binding yourself the desired parameters
   *  and leaving RFP the ownership of the remaining parameters.
   *  For instance, assuming D=3, CoordT = double, and a user-defined function
   *  whihc depends on 3 parameters:
   *
   *  CoordT = my_radial_funct_2( CoordT r, Coord_T par1, CoordT par2, CoordT par3 ) defined somewhere,
   *  double par1 = //somevalue,
   *         par3 = //somevalue;
   *  auto my_rf = new RFP<3, 1 //nr. of parameters that will be managed by RFP// , double>(
   *    bitpit::rbf::eRBFType::kUserDefined,
   *    std::bind( &my_radial_funct_2, std::placeholders::_1, std::cref(par1), std::placeholders::_2, std::cred(par3)
   *  );
   *
   *  In this case, my_rf will have ownership on 1 parameter (par2) while the ownership of the 2 
   *  other parameters is left to the user.
  */
  template<
    std::size_t Dim,
    std::size_t NParams,
    class CoordT = double
  >
  class RFP : public RF<Dim, CoordT>
  {
    // Static asserts ============================================== //
    static_assert(
      (NParams > 0),
      "bitpit::rbf::RFP<Dim, NParams, CoordT>: "
      "**ERROR** The nr. of additional parameters must be greater than 0. "
      "If the radial function does not depends on any parameter, use RF<Dim, CoordT>"
    );

    // Typedef(s) ================================================== //
    private:
    /*! @brief Type of the base class. */
    using base_t      = RF<Dim, CoordT>;
    /*! @brief Type of this class. */
    using self_t      = RFP<Dim, NParams, CoordT >;
    public:
    /*! @brief Type of functor holding the actual expression of the RF. */
    using rf_funct_t  = typename base_t::rf_funct_t;
    /*! @brief Coeff. type. */
    using coord_t     = typename base_t::coord_t;
    /*! @brief Point type. */
    using point_t     = typename base_t::point_t;

    // Friendships ================================================= //
    friend class RF<Dim, CoordT>;

    // Member variable(s) ========================================== //
    protected:
    using base_t::mType;
    using base_t::mFunct;
    /*! @brief Values of the additional parameters for this radial function. */
    coord_t mParams[NParams];
    public:
    using base_t::center;
    using base_t::radius;

    // Constructor(s) ============================================== //
    private:
    /*! @brief Default constructor.
     *
     *  Initialize a generalized radial function with default values for the parameters.
    */
    RFP();
    public:
    /*! @brief Contructor #1.
     *
     *  Initialize a generalized radial function of the specified type,
     *  from the input pointer to the function implementating the expression of the RF.
     *
     *  @param [in]       type      type of this RF.
     *  @param [in]       f         pointer to the function implementing the expression
     *                              of this RF.
    */
    template<class ...Args>
    RFP( bitpit::rbf::eRBFType type, coord_t(*f)( coord_t, Args ...args) );
    /*! @brief Copy constructor (deleted) */
    RFP( const self_t & ) = delete;
    /*! @brief Move-constructor (deleted). */
    RFP( self_t && ) = delete;

    // Operator(s) ================================================= //
    public:
    /*! @brief Copy-assignment operator (delete) */
    self_t& operator=( const self_t & ) = delete;
    /*! @brief Move-assignement operator (delete) */
    self_t& operator=( self_t &&other ) = delete;

    // Getter(s)/Info ============================================== //
    public:
    /*! @brief Returns the nr. of additional parameters for this radial function. */
    virtual std::size_t getNumberOfParameters() const override;
    /*! @brief Returns (const) pointer to the list of additional parameters of this funciton. */
    virtual const coord_t*  getParameters() const override;
    /*! @brief Display info for this RBF (mostly meant for debugging purposes).
     *
     *  @param [in,out]   out       (default = std::cout) output stream.
     *  @param [in]       indent    (default = 0) indentation level.
    */
    virtual void display( std::ostream &out = std::cout, unsigned int indent = 0 ) const override;

    // Setter(s) =================================================== //
    private:
    /*! @brief Bind the expression of the radial function to the parameters of
     *  this instance.
    */
    template<class f_in_t, std::size_t... I>
    static rf_funct_t       doBind( f_in_t f, coord_t* data, bitpit::index_sequence<I...> );
    template<class f_in_t>
    static rf_funct_t       bindParameters( f_in_t f, coord_t* data );public:
    /*! @brief Set the value for the additional parameters of this function.*/
    virtual void            setParameters( const coord_t *values ) override;
  }; //end class RFP

  // =============================================================== //
  // DEFINITION OF CLASS RFBasis                                     //
  // =============================================================== //
  /*! @brief Radial Basis Function.
   *
   *  Class storing a basis of (heterogeneous) Radial Functions.
   *
   *  @tparam     Dim     nr. of working dimensions.
   *  @tparam     CoordT  (default = double) type for coordinates/coefficients. Only scalar
   *                      floating point types are supported (e.g. float, double, etc.)
  */
  template<
    std::size_t Dim,
    class CoordT = double
  >
  class RFBasis : private std::vector< std::pair<CoordT, std::unique_ptr< RF<Dim,CoordT> > > >
  {
    // Static assertion(s) ========================================= //
    static_assert(
      std::is_floating_point<CoordT>::value,
      "**ERROR** bitpit::rbf::RFBasis<Dim,CoordT>: CoordT must be a scalar floating point type "
      ", e.g. float, double or long double"
    );
    
    // Typedef(s) ================================================== //
    private:
    /*! @brief Type of the base class. */
    using base_t    = std::vector< std::pair<CoordT, std::unique_ptr<RF<Dim,CoordT> > > >;
    /*! @brief Self type. */
    using self_t    = RFBasis<Dim, CoordT>;
    public:
    /*! @brief Type of Radial function. */
    using rf_t      = RF<Dim,CoordT>;
    /*! @brief Type of coordinate. */
    using coord_t   = typename rf_t::coord_t;
    /*! @brief Type of point in the Euclidean space. */
    using point_t   = typename rf_t::point_t;
    
    // Constructor(s) ============================================== //
    public:
    /*! @brief Default constructor.
     *
     *  Initialize a empty radial basis.
    */
    RFBasis();
    /*! @brief Constructor #1.
     *
     *  Initialize a radial basis function with N functions of the specified,
     *  type.
     *
     *  @param [in]     N       nr. of functions.
     *  @param [in]     type    type of radial basis.
    */
    RFBasis( std::size_t n, eRBFType type );
    /*! @brief Copy-constructor (deleted). */
    RFBasis( const self_t & ) = delete;
    /*! @brief Move-constructor (deleted). */
    RFBasis( self_t && ) = delete;
    
    // Operator(s) ================================================= //
    public:
    /*! @brief Copy-assignment operator (deleted). */
    self_t& operator=( const self_t & ) = delete;
    /*! @brief Move-assignment operator (deleted). */
    self_t& operator=( self_t && ) = delete;
    /*! @brief Evaluate the linear combination of the RFs at the input point. */
    coord_t operator()( const point_t &coords ) const;
    
    // Getter(s)/Info ============================================== //
    public:
    /*! @brief Returns the size of this basis */
    using base_t::size;
    /*! @brief Returns (const/non-const) reference to the i-th radial function in the basis. */
    using   base_t::operator[];
    using   base_t::at;
    /*! @brief Returns the collection of weights for this radial basis function. */
    std::vector<coord_t> collectWeights() const;
    /*! @brief Returns const/non-const reference to the weight of the i-th radial function. */
    const coord_t& getWeight( std::size_t i ) const;
    coord_t& getWeight( std::size_t i );
    /*! @brief Returns const/non-const reference to the i-th radial function. */
    const rf_t& getRadialFunction( std::size_t i ) const;
    rf_t& getRadialFunction( std::size_t i );
    /*! @brief Display info to the output stream provided as input.
     *  (Mostly meant for debugging purposes).
     *
     *  @param [in,out]   out     (default = std::cout) output stream.
     *  @param [in]       indent  (default = 0) indentation level.
    */
    void display( std::ostream &out = std::cout, unsigned int indent = 0 ) const;

    // Setter(s) =================================================== //
    public:
    /*! @brief Add a new radial function to the basis. 
     *
     *  @param [in]   rf      (unique) Pointer to the radial function.
     *  @param [in]   weight  (default = 1) weight to be assigned to the new function.
     *
     *  @result Returns the position in the basis where the new function
     *  has been added.
    */
    std::size_t add( std::unique_ptr<rf_t> rf, CoordT weight = (CoordT)1 );
    /*! @brief Remove the i-th function from the basis. */
    void remove( std::size_t i );
    /*! @brief Utility function to set the weights of each radial function. */
    void setWeights( const std::vector<coord_t> &weights );
    /*! @brief Overloading of #setWeights taking a range iterator as input. 
     *
     *  Set the weight for each radial function from the range [first, last).
     *  [first, last) must contain at least N values (N being the size of this basis).
     *
     *  @tparam       ItaratorType      iterator type. IteratorType must be at
     *                                  least a forward iterator.
     * 
     *  @param [in]   first, last       iterators pointing to the beginning/end of the 
     *                                  input range.
    */
    template< class IteratorType >
    void setWeights( IteratorType first, IteratorType last );
    /*! @brief Set the specified radius for all radial function in this basis. */
    void setRadius( coord_t radius );
  }; //end class RFBasis
  
  // =============================================================== //
  // EXPLICIT SPECIALIZATIONS                                        //
  // =============================================================== //
  extern template class RF<1, float>;
  extern template class RF<2, float>;
  extern template class RF<3, float>;
  extern template class RF<1, double>;
  extern template class RF<2, double>;
  extern template class RF<3, double>;
  extern template class RF<1, long double>;
  extern template class RF<2, long double>;
  extern template class RF<3, long double>;

  extern template class RFP<1, 1, float>;
  extern template class RFP<2, 1, float>;
  extern template class RFP<3, 1, float>;
  extern template class RFP<1, 1, double>;
  extern template class RFP<2, 1, double>;
  extern template class RFP<3, 1, double>;
  extern template class RFP<1, 1, long double>;
  extern template class RFP<2, 1, long double>;
  extern template class RFP<3, 1, long double>;

  extern template class RFP<1, 2, float>;
  extern template class RFP<2, 2, float>;
  extern template class RFP<3, 2, float>;
  extern template class RFP<1, 2, double>;
  extern template class RFP<2, 2, double>;
  extern template class RFP<3, 2, double>;
  extern template class RFP<1, 2, long double>;
  extern template class RFP<2, 2, long double>;
  extern template class RFP<3, 2, long double>;
  
  extern template class RFBasis<1, float>;
  extern template class RFBasis<2, float>;
  extern template class RFBasis<3, float>;
  extern template class RFBasis<1, double>;
  extern template class RFBasis<2, double>;
  extern template class RFBasis<3, double>;
  extern template class RFBasis<1, long double>;
  extern template class RFBasis<2, long double>;
  extern template class RFBasis<3, long double>;
  /*! @} */

} //end namespace rbf

} //end namespace bitpit


// =============================================================== //
// TEMPLATE IMPLEMENTATIONS                                        //
// =============================================================== //
# include "rbf.tpp"