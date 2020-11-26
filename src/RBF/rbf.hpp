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
   *          &(1 - r )^4 ( 4 r + 1 ), \; \text{if} r < 1, \\
   *          &0, \; \text{otherwise}
   *        \end{aligned}
   *    \right.
   *  \end{equation}
   *  \f$
   *  where \f$r\f$ is the radial distance of a given point
   *  from the geometrical kernel (i.e. the center of the RF).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT  wendland_c2( CoordT r );
  
  // ----------------------------------------------------------------------- //
  /*! @brief First derivative of C2-continuous Wendland's functions.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT wendland_c2_der1( CoordT r );
  
  // ----------------------------------------------------------------------- //
  /*! @brief Second derivative of C2-continuous Wendland's functions.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT wendland_c2_der2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Generalized multiquadric functions.
   *
   *  The family of generalized multi-quadrics RF has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r;c) := \left( r^2 + c^2 \right)^{\frac{\alpha}{\beta}}
   *  \end{equation}
   *  \f$
   *  for \f$\alpha, \beta > 0\f$.
   *  In the above expression, \f$r\f$ is the radial distance of a given point
   *  from the geometric kernel (i.e. the center of the RF), and \f$c\f$ is 
   *  a additional parameter (bias).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating 
   *                              point types are supported (e.g. double, float, etc. )
   *  @tparam         Alpha, Beta Values for the fractional exponent (Alpha, Beta > 0 ).
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
  /*! @brief First derivative of generalized multiquadric. 
   *
   *  @tparam         CoordT      Tyoe of RF coeffs. Only scalar floating 
   *                              point types are supported (e.g. double, float, etc. )
   *  @tparam         Alpha, Beta Values for the fractional exponent (Alpha, Beta > 0 ).
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
  CoordT generalized_multiquadrics_der1( CoordT r, CoordT c );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of generalized multiquadric. 
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating 
   *                              point types are supported (e.g. double, float, etc. )
   *  @tparam         Alpha, Beta Values for the fractional exponent (Alpha, Beta > 0 ).
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
  CoordT generalized_multiquadrics_der2( CoordT r, CoordT c );
  
  // ---------------------------------------------------------------- //
  /*! @brief Gaussian radial function.
   *
   *  A Guassian radial function has the following expression:
   *  A Guassian radial basis function has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r;c) := c exp( -r^2 );
   *  \end{equation}
   *  \f$
   *  where \f$c\f$ is a amplification/reduction factor, and \f$r\f$ is the
   *  radial distance from the geometry kernel (i.e. the center of the radial function).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
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
  
  // ----------------------------------------------------------------------- //
  /*! @brief First derivative of Gaussian function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
   *  @param [in]     c           amplification/reduction factor.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT gaussian_der1( CoordT r, CoordT c );
  
  // ----------------------------------------------------------------------- //
  /*! @brief Second derivative of Gaussian function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance.
   *  @param [in]     c           amplification/reduction factor.
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT gaussian_der2( CoordT r, CoordT c );
  
  // ---------------------------------------------------------------- //
  /*! @brief Linear function.
   *
   *  Linear function with trivial expression (generally used to as a bias):
   *  \f$ f(r) := r\ f$
   *  where \f$r\f$ is the radial distance from the geometry kernel (i.e. the 
   *  center of the radial function).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
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
  /*! @brief First derivative of linear function.
   *
   *  @tparam         CoordT      TYpe of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT linear_der1( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of linear function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT linear_der2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Constant function.
   *
   *  Constant function with trivial expression (generally used to as a bias):
   *  \f$ f(r) := 1\ f$
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT constant( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief First derivative of constant function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT constant_der1( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of constant function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT constant_der2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Generalized power function (radial power).
   *
   *  Generalized power functions have the following expression:
   *  \f$ f(r) := r^\alpha \f$
   *  where \f$\alpha>0\f$ is the power exponent and \f$r\f$ is the radial distance
   *  from the geomtry kernel (i.e. the center of the RF).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. )
   *  @tparam         Alpha       power exponent (>0).
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
  /*! @brief First derivative of generalized power function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. )
   *  @tparam         Alpha       power exponent (>0).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    int   Alpha,
    typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) >::type* = nullptr
  >
  CoordT generalized_power_der1( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of generalized power function.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
   *                              are supported (e.g. double, float, etc. )
   *  @tparam         Alpha       power exponent (>0).
   *
   *  @param [in]     r           radial distance.
  */
  template<
    class CoordT,
    int   Alpha,
    typename std::enable_if< std::is_floating_point<CoordT>::value && (Alpha > 0) >::type* = nullptr
  >
  CoordT generalized_power_der2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Thin plate splines.
   *
   *  Thin plate splines have the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r):= r^2 log(r)
   *  \end{equation}
   *  \f$
   *  where \f$r\f$ is the radial distance from the geometry kernel
   *  (i.e. the center of the RF).
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types are
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
  /*! @brief First derivative of thin plate splines.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types are
   *                              supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT    thin_plate_spline_der1( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of thin plate splines.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types are
   *                              supported (e.g. double, float, etc. ).
   *
   *  @param [in]     r           radial distance
  */
  template<
    class CoordT,
    typename std::enable_if< std::is_floating_point<CoordT>::value >::type* = nullptr
  >
  CoordT    thin_plate_spline_der2( CoordT r );
  
  // ---------------------------------------------------------------- //
  /*! @brief Polyharmonic radial functions.
   *
   *  The family of polyharmonic functions has the following expression:
   *  \f$
   *  \begin{equation}
   *    f(r):= (-c)^{1+\beta) r^{2\beta} \, log(r), \; \beta \neq 0
   *  \end{equation}
   *  \f$
   *  where \f$r\f$ is the radial distance, \f$\beta \ge 0 \f$is the power
   *  exponent, and \f$c\f$ is a amplification/reduction factor.
   *
   *  @note The special case \f$beta = 2\f$ correspnds to thin-plate-splines.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
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
  
  // ---------------------------------------------------------------- //
  /*! @brief First derivative of polyharmonic radial functions.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
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
  CoordT polyharmonic_der1( CoordT r, CoordT c );
  
  // ---------------------------------------------------------------- //
  /*! @brief Second derivative of polyharmonic radial functions.
   *
   *  @tparam         CoordT      Type of RF coeffs. Only scalar floating point types
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
  CoordT polyharmonic_der2( CoordT r, CoordT c );

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
    /*! @brief Generalized multiquadric with \f$\alpha = 1, \beta = 2\f$ */
    kMultiQuadric1_2,
    /*! @brief Generalized multiquadric with \f$\alpha = 2, \beta = 1\f$ */
    kMultiQuadric2,
    /*! @brief Generalized multiquadric with \f$\alpha = 3, \beta = 2\f$ */
    kMultiQuadric3_2,
    /*! @brief Generalized multiquadric with \f$\alpha = 5, \beta = 2\f$ */
    kMultiQuadric5_2,
    /*! @brief Linear */
    kLinear,
    /*! @brief Constant */
    kConstant,
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
  // DEFINITION OF CLASS RF         									                //
  // ================================================================ //
  /*! @brief Radial function with no additional parameters.
   *
   *  This class holds a radial function (in short, RF) which does not depend
   *  on any additional parameter beyond its radius and its center.
   *
   *	@tparam 			Dim 		    Nr. of dimension in the working space.
   *	@tparam 			CoordT 	    (default = double) Type of coefficients.
   *                            Only scalar floating point type are supported (e.g. double, float, etc.).
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
      "** ERROR ** bitpit::rbf::RF<Dim,CoordT>: CoordT must be a scalar floating point type "
      ", e.g. float, double or long double"
    );
    static_assert(
      (Dim > 0),
      "** ERROR ** bitpit::rbf::RF<Dim,CoordT>: Dim must be greater than 0."
    );

    // Typedef(s) =================================================== //
    public:
    /*!	@brief Type of coeffs. */
    using coord_t 		= CoordT;
    /*!	@brief Point type. */
    using point_t 		= std::array<coord_t, Dim>;
    /*!	@brief Type of functor holding the actual implementation of the RF.*/
    using rf_funct_t  	= std::function< coord_t( coord_t ) >;
    private:
    /*!	@brief Type of this object. */
    using self_t    	= RF<Dim, CoordT>;

    // Member variable(s) =========================================== //
    protected:
    /*!	@brief Functor implementing the actual expression of the Radial Function. */
    rf_funct_t mFunct;
    /*! @brief Type of this radial function. */
    eRBFType mType;
    /*! @brief Flag for compactly supported RFs.*/
    bool mHasCompactSupport;
    public:
    /*!	@brief Radius of this RF. */
    coord_t radius;
    /*!	@brief Center of this RF (also referred to as control point, geometry kernel) */
    point_t center;

    // Static member function(s) ==================================== //
    public:
    /*! @brief Returns a pointer to a new instance of a RF 
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
     *	@param [in]		  coords 	   coordinates of the input point.
    */
    virtual coord_t operator()( const point_t &coords ) const;
    
    // Getter(s)/Info ============================================== //
    public:
    /*!	@brief Returns the nr. of additional parameters for this RF. */
    virtual std::size_t	getNumberOfParameters() const;
    /*!	@brief Returns (true) if this RF is compactly supported. */
    bool hasCompactSupport() const;
    /*!	@brief Returns (const) pointer to the internal array storing
     *  the values of the additional parameters of this RF. */
    virtual const coord_t* getParameters() const;
    /*! @brief Returns the type of this radial function. */
    eRBFType getType() const;
    /*! @brief Display info to a output stream
     *  (mostly meant for debugging purposes).
     *
     *  @param [in,out]   out       (default = std::cout) output stream.
     *  @param [in]       indent    (default = 0) indentation level.
    */
    virtual void display( std::ostream &out = std::cout, unsigned int indent = 0 ) const;
    /*! @brief Returns reference to the actual functor implementing the expression of this
     *  radial function. */
    const rf_funct_t& getFunctor() const
    {
      return mFunct;
    }
    /*! @brief Returns a functor implementing the the first derivative 
     *  of this radial function w.r.t. the radial distance. */
    rf_funct_t getFirstDerivative() const;
    /*! @brief Returns a functor implementing the second derivative 
     *  of this radial function w.r.t. the radial distance. */
    rf_funct_t getSecondDerivative() const;
    
    // Setter(s) =================================================== //
    public:
    /*! @brief Set default values for function parameters, ie.:
     *  * radius = 1 
     *  * center = (0, 0, 0)
    */
    virtual void setDefault();
    /*!	@brief Set the value of the parameters of this Radial Function.
     *
     *  @param [in]       vals      const pointer to the array storing the values
     *                              of the parameter.
    */
    virtual void setParameters( const coord_t * );

  }; //end class RF

  // =============================================================== //
  // DEFINITION OF CLASS RFP                                         //
  // =============================================================== //
  /*! @brief Definition of generalized RF.
   *
   *  This class can be used to store a generic radial function
   *  which depends on a arbitrary nr. of parameters of type CoordT
   *  beyond its radius and center.
   *
   *  ** Usage **
   *  * Default usage *
   *  If you wish to use one of the radial function implemented in bitpit, just create a new instance
   *  by calling:
   *  RF::New (specify the type of radial function you want to use).
   *
   *  * Passing the ownership of the additional parameters to RFP class. *
   *  If you wish to encapsulate a given function (which depends on N additional parameters),
   *  you just need to invoke the RFP constructor passing the pointer to that function. 
   *  The nr. of additional parameters must be specified as a template argument of the RF class.
   *
   *  For instance, assumung Dim=3, CoordT = double and a user-defined funtion which depends
   *  on 2 additional parameters:
   *
   *  CoordT = my_radial_funct( CoordT r, CoordT par1, CoordT par2 ) { ... }
   *  auto my_rf = new RFP<3, 2, double>( bitpit::rbf::eRBFType::kUserDefined, &my_radial_funct )
   *  
   *  will create a radial function with the specified expression and automatically bind 
   *  the additional parameters to some inernal member variable which can be accesed/modified
   *  any time via #setParameters and #getParameters member functions.
   *
   *  * Retaining the ownership of the additional parameters *
   *  In some cases you might want to keep the ownership over some of the additional parameters.
   *  You can keep the ownership of such parameters by binding yourself the desired parameters
   *  to the function implementing the expression of the RF and leaving the ownership of the remaining parameters.
   *  to RFP.
   *
   *  For instance, assuming Dim=3, CoordT = double, and a user-defined function
   *  which depends on 3 parameters:
   *
   *  CoordT = my_radial_funct_2( CoordT r, Coord_T par1, CoordT par2, CoordT par3 ) {...},
   *  double par1 = //some value,
   *         par3 = //some value;
   *  auto my_rf = new RFP<3, 1 //nr. of parameters that will be managed by RFP// , double>(
   *    bitpit::rbf::eRBFType::kUserDefined,
   *    std::bind( &my_radial_funct_2, std::placeholders::_1, std::cref(par1), std::placeholders::_2, std::cred(par3)
   *  );
   *
   *  will create a instance of RFP with exclusive ownership on just 1 parameter (par2), while the ownership of the 2 
   *  other parameters is left to you.
   *
   *  @tparam     Dim       Nr. of working dimensions (Dim > 0)
   *  @tparam     NParams   Nr. of additional parameters which the RF function depends on.
   *  @tparam     CoordT    (default = double) type of coeff. Only scalar floating point
   *                        types are supported (e.g. double, float, etc.)
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
      "** ERROR ** The nr. of additional parameters must be greater than 0. "
      "If the radial function does not depends on any parameter, use RF<Dim, CoordT>"
    );
    static_assert(
      Dim > 0,
      "bitpit::rbf::RFP<Dim, NParams, CoordT>: "
      "** ERROR ** The nr. of dimensions, Dim, must be greater than 0"
    );
    static_assert(
      std::is_floating_point<CoordT>::value,
      "bitpit::rbf::RFP<Dim, NParams, CoordT>: "
      "** ERROR ** Only floating point types are supported for the coeff. type, CoordT"
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
    /*! @brief Type of this radial function. */
    using base_t::mType;
    /*! @brief Functor implementing the actual expression of the RF. */
    using base_t::mFunct;
    /*! @brief Values of the additional parameters for this radial function which are binded at construction time. */
    coord_t mParams[NParams];
    public:
    /*! @brief Center of this radial function. */
    using base_t::center;
    /*! @brief Radius of this radial function. */
    using base_t::radius;

    // Static method(s) ============================================ //
    private:
    /*! @brief Do the actual binding. */
    template<class f_in_t, std::size_t... I>
    static rf_funct_t       doBind( f_in_t f, coord_t* data, bitpit::index_sequence<I...> );
    public:
    /*! @brief Bind the expression of the radial function to the parameters of
     *  this instance. */
    template<class f_in_t>
    static rf_funct_t       bindParameters( f_in_t f, coord_t* data );
    
    // Constructor(s) ============================================== //
    private:
    /*! @brief Default constructor.
     *
     *  Initialize a generalized radial function with default values for the parameters
     *  (see #setDefault).
    */
    RFP();
    public:
    /*! @brief Contructor #1.
     *
     *  Initialize a generalized radial function of a specified type and given expression.
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
    /*! @brief Returns the nr. of additional parameters of this radial function. */
    virtual std::size_t getNumberOfParameters() const override;
    /*! @brief Returns (const) pointer to the list of additional parameters of this funciton. */
    virtual const coord_t* getParameters() const override;
    /*! @brief Display info to a output stream (mostly meant for debugging purposes).
     *
     *  @param [in,out]   out       (default = std::cout) output stream.
     *  @param [in]       indent    (default = 0) indentation level.
    */
    virtual void display( std::ostream &out = std::cout, unsigned int indent = 0 ) const override;

    // Setter(s) =================================================== //
    public:
    /*! @brief Set the value of the additional parameters for this function.*/
    virtual void  setParameters( const coord_t *values ) override;
  }; //end class RFP

  // =============================================================== //
  // DEFINITION OF CLASS RFBasis                                     //
  // =============================================================== //
  /*! @brief Radial Basis Function.
   *
   *  Class storing a basis of (heterogeneous) Radial Functions.
   *
   *  @tparam     Dim     nr. of working dimensions (Dim > 0).
   *  @tparam     CoordT  (default = double) type of coefficients. Only scalar
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
      Dim > 0,
      "bitpit::rbf::RFBasis<Dim,CoordT>: "
      "** ERROR ** Nr. of working dimensions, Dim, must be greater than 0"
    );
    static_assert(
      std::is_floating_point<CoordT>::value,
      "bitpit::rbf::RFBasis<Dim,CoordT>: "
      "** ERROR ** CoordT must be a scalar floating point type (e.g. float, double, etc. )"
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
     *  Initialize a radial basis function with N functions of the specified type.
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
    /*! @brief Reserve memory for the insertion of new radial functions.
     *
     *  If the specified reserve is smaller than the current capacity
     *  no action is taken.
     *
     *  @param [in]     n       memory to be reserved in terms of nr. of elements.
    */
    using base_t::reserve;
    /*! @brief Returns const/non-const iterator pointing to the first
     *  pair radial function-weight. */
    using base_t::cbegin;
    using base_t::begin;
    /*! @brief Returns const/non-const iterator pointing at the end
     *  of the list of radial functions (i.e. after the last radial function). */
    using base_t::cend;
    using base_t::end;
    /*! @brief Returns the size of this basis. */
    using base_t::size;
    /*! @brief Returns (const/non-const) reference to a pair consisting of the i-th radial function 
     *  and the associated weight. */
    using base_t::operator[];
    using base_t::at;
    /*! @brief Returns (true) if this radial function basis is empty (i.e. it does not contain any function). */
    using base_t::empty;
    /*! @brief Clear the content of this class and reset its state to default (empty basis). */
    using base_t::clear;
    /*! @brief Returns the collection of weights for this radial basis function. */
    std::vector<coord_t> collectWeights() const;
    /*! @brief Returns const/non-const reference to the weight of the i-th radial function. */
    const coord_t& getWeight( std::size_t i ) const;
    coord_t& getWeight( std::size_t i );
    /*! @brief Returns const/non-const reference to the i-th radial function. */
    const rf_t& getRadialFunction( std::size_t i ) const;
    rf_t& getRadialFunction( std::size_t i );
    /*! @brief Display info to a output stream
     *  (mostly meant for debugging purposes).
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
     *  @note All iterators, pointer, references to any radial function/weight
     *  might be invalidated after calling this method.
     *
     *  @result Returns the ID assigned to the new radial function. This ID can
     *  be used to access the radial function via #getRadialFunction and #getWeight
     *  methods.
    */
    std::size_t add( std::unique_ptr<rf_t> rf, coord_t weight = (coord_t)1 );
    /*! @brief Remove the i-th function from the basis.
     *  
     *  @note All iterators, pointers, references to any radial function/weight
     *  might be invalidated after calling this method.
     *
     *  @param [in]   i       ID of the radial function to be removed.
    */
    void remove( std::size_t i );
    /*! @brief Utility function which assign a constant weight to all radial function
     *  of this basis. */
    void setWeights( coord_t weight );
    /*! @brief Set the weight of each radial function from the input list of values. */
    void setWeights( std::vector<coord_t> const &weights );
    /*! @brief Set the weight of each radial function from the input range of values.
     *
     *  @tparam       ItaratorType      iterator type. IteratorType must be at
     *                                  least a forward iterator.
     * 
     *  @param [in]   first, last       iterators pointing to the beginning/end of the 
     *                                  input range.
    */
    template< class IteratorType >
    void setWeights( IteratorType first, IteratorType last );
    /*! @brief Set the radius of all radial functions in this basis to the specified value. */
    void setRadii( coord_t radius );
    /*! @brief Set the radius of each radial function from the input list of values. */
    void setRadii( std::vector<coord_t> const &radii );
    /*! @brief Set the radius of each radial function from the input range of values.
     *
     *  @tparam       ItaratorType      iterator type. IteratorType must be at
     *                                  least a forward iterator.
     * 
     *  @param [in]   first, last       iterators pointing to the beginning/end of the 
     *                                  input range.
    */
    template< class IteratorType >
    void setRadii( IteratorType first, IteratorType last );
    /*! @brief Set the center of all radial functions in this basis to the specified value. */
    void setCenters( point_t const &center );
    /*! @brief Set the center of all radial functions in this basis from the input list of points. */
    void setCenters( std::vector<point_t> const &centers );
    /*! @brief Set the center of all radial functions in this basis from the input range of points.
     *
     *  @tparam       ItaratorType      iterator type. IteratorType must be at
     *                                  least a forward iterator.
     * 
     *  @param [in]   first, last       iterators pointing to the beginning/end of the 
     *                                  input range.
    */
    template< class IteratorType >
    void setCenters( IteratorType first, IteratorType last );
    /*! @brief Reset the type of this basis of radial functions.
     *
     *  @param [in]   type    new type for this radial basis.
     *  @param [in]   n       (default = -1) new size of this radial basis.
     *                        If not specified the current size will be kept.
    */
    void reset( eRBFType type, std::size_t n = -1 );
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