/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#ifndef __BITPIT_RBF_HPP__
#define __BITPIT_RBF_HPP__

#include "bitpit_discretization.hpp"
#include <array>
#include <set>
#include <vector>

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
    COSINUS    = 14, /**< Compact cosinusoidal funct, value of 1 on r=0, 0 outside */
    THINPLATE  = 15, /**< Non compact thin plate funct */
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
    /*!
     * Linear Polynomial class. A linear multivariate polynomial in function of the 
     * variables defining the RBF nodes is used to regularize the interpolation 
     * problem through RBFs. It can be used only when RBFMode::INTERP is enabled. 
     * The LinearPolynomial class derives from ReconstructionPolynomial class
     * that is intrisecally defined to work in a three dimensional space.
     */
    class LinearPolynomial : public ReconstructionPolynomial
    {
private:
        int m_dim;                          /**< Polynomial variables space dimension */
        int m_fields;                       /**< Number of equations of the polynomial object */

    public:
        LinearPolynomial();
        void clear();
        void setDimension(int dim);
        void setDataCount(int fields);
        void evalBasis(const double *x, double *basis);
        void initialize();
    };

private:
    int     m_fields;                               /**<Number of data fields defined on RBF nodes.*/
    RBFMode m_mode;                                 /**<Behaviour of RBF class (interpolation or parametrization).*/
    std::vector<double>  m_supportRadii;            /**<Support radii of function used as Radial Basis Function.*/
    RBFBasisFunction m_typef;                       /**<Recognize type of RBF shape function actually in the class. */
    double  (*m_fPtr)(double);

    std::vector<double>                 m_error;    /**<Interpolation error of a field evaluated on each RBF node (auxiliary memeber used in Greedy algorithm).*/

protected:
    std::vector<std::vector<double>>    m_value;    /**< displ value to be interpolated on RBF nodes */
    std::vector<std::vector<double>>    m_weight;   /**< weight of your RBF interpolation */
    std::vector<bool>                   m_activeNodes;   /**<Vector of active/inactive node (m_activeNodes[i] = true/false -> the i-th node is used/not used during RBF evaluation).*/   
    int m_maxFields;                                /**< fix the maximum number of fields that can be added to your class*/
    int m_nodes;                                    /**<Number of RBF nodes.*/
    bool m_polyEnabled;                             /**< Enable/disable the use of the linear polynomial term in interpolation */
    LinearPolynomial m_polynomial;                  /**< Linear polynomial object */
    std::vector<int> m_polyActiveBasis;             /**< Active terms of linear polynomial, 0 is constant, i+1 the i-th system coordinate */

public:
    RBFKernel();
    RBFKernel(const RBFKernel & other);

    virtual ~RBFKernel() = default;

    void                    setFunction(RBFBasisFunction);
    void                    setFunction(double (&funct)(double ));

    RBFBasisFunction        getFunctionType();
    int                     getDataCount();
    int                     getActiveCount();
    std::vector<int>        getActiveSet();

    void                    enablePolynomial(bool enable = true);
    int                     getPolynomialDimension();
    int                     getPolynomialWeightsCount();

    bool                    isActive(int );

    bool                    activateNode(int );
    bool                    activateNode(const std::vector<int> &);
    void                    activateAllNodes();
    bool                    deactivateNode(int );
    bool                    deactivateNode(const std::vector<int> &);
    void                    deactivateAllNodes();

    void                    setSupportRadius(double);
    void                    setSupportRadius(const std::vector<double> &);
    double                  getSupportRadius();
    double                  getSupportRadius(int);

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
    double                  evalBasisPair(int i, int j);

    int                     solve();
    int                     greedy(double);

    const std::vector<std::vector<double>> & getWeights() const;

protected:
    double                  evalError();
    int                     addGreedyPoint();
    int                     solveLSQ();
    void                    swap(RBFKernel & x) noexcept;

private:

    virtual double calcDist(int i, int j) = 0;
    virtual double calcDist(const std::array<double,3> & point, int j)                  = 0;
    virtual std::vector<double> evalPolynomialBasis(int i)                              = 0;
    virtual std::vector<double> evalPolynomialBasis(const std::array<double,3> &point)  = 0;
    virtual void initializePolynomialActiveBasis()                                      = 0;
    virtual void initializePolynomial()                                                 = 0;
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
    double calcDist(int i, int j) override;
    double calcDist(const std::array<double,3> & point, int j) override;
    void initializePolynomialActiveBasis() override;
    void initializePolynomial() override;
    std::vector<double> evalPolynomialBasis(int i) override;
    std::vector<double> evalPolynomialBasis(const std::array<double, 3> &point) override;
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
    double                  cosinus(double);
    double                  thinplate(double);
} // namespace rbf
}

#endif
