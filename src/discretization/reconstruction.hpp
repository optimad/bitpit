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

#ifndef __BTPIT_RECONSTRUCTION_HPP__
#define __BTPIT_RECONSTRUCTION_HPP__

#include <array>
#include <iostream>
#include <memory>
#include <vector>

#include "bitpit_patchkernel.hpp"

namespace bitpit {

class ReconstructionPolynomial {

friend class ReconstructionKernel;
friend class ReconstructionAssembler;

public:
    ReconstructionPolynomial();
    ReconstructionPolynomial(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin, int nFields=1);

    ReconstructionPolynomial(const ReconstructionPolynomial &other);
    ReconstructionPolynomial(ReconstructionPolynomial &&other) = default;
    ReconstructionPolynomial & operator = (const ReconstructionPolynomial &other);
    ReconstructionPolynomial & operator=(ReconstructionPolynomial &&other) = default;

    void swap(ReconstructionPolynomial &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin, int nFields=1, bool release = true);
    void clear(bool release = true);

    uint8_t getDegree() const;
    uint8_t getDimensions() const;

    uint16_t getCoefficientCount() const;

    int getFieldCount() const;

    const double * getCoefficients() const;
    double * getCoefficients();
    const double * getCoefficients(int field) const;
    double * getCoefficients(int field);
    const double * getDegreeCoefficients(uint8_t degree, int field = 0) const;
    double * getDegreeCoefficients(uint8_t degree, int field = 0);

    void computeValue(const std::array<double, 3> &point, int field, double *values) const;
    void computeValues(const std::array<double, 3> &point, double *values) const;
    void computeValues(const std::array<double, 3> &point, int nFields, double *values) const;
    void computeValues(const std::array<double, 3> &point, int nFields, int offset, double *values) const;
    void computeValue(int degree, const std::array<double, 3> &point, int field, double *values) const;
    void computeValues(int degree, const std::array<double, 3> &point, double *values) const;
    void computeValues(int degree, const std::array<double, 3> &point, int nFields, double *values) const;
    void computeValues(int degree, const std::array<double, 3> &point, int nFields, int offset, double *values) const;

    void computeValueLimited(const std::array<double, 3> &point, const double *limiters, int field, double *value) const;
    void computeValuesLimited(const std::array<double, 3> &point, const double *limiters, double *values) const;
    void computeValuesLimited(const std::array<double, 3> &point, const double *limiters, int nFields, double *values) const;
    void computeValuesLimited(const std::array<double, 3> &point, const double *limiters, int nFields, int offset, double *values) const;
    void computeValueLimited(int degree, const std::array<double, 3> &point, const double *limiters, int field, double *values) const;
    void computeValuesLimited(int degree, const std::array<double, 3> &point, const double *limiters, double *values) const;
    void computeValuesLimited(int degree, const std::array<double, 3> &point, const double *limiters, int nFields, double *values) const;
    void computeValuesLimited(int degree, const std::array<double, 3> &point, const double *limiters, int nFields, int offset, double *values) const;

    void computeDerivative(const std::array<double, 3> &point, const std::array<double, 3> &direction, int field, double *derivative) const;
    void computeDerivatives(const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivatives) const;
    void computeDerivatives(const std::array<double, 3> &point, const std::array<double, 3> &direction, int nFields, double *derivatives) const;
    void computeDerivatives(const std::array<double, 3> &point, const std::array<double, 3> &direction, int nFields, int offset, double *derivatives) const;
    void computeDerivative(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, int field, double *derivative) const;
    void computeDerivatives(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivatives) const;
    void computeDerivatives(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, int nFields, double *derivatives) const;
    void computeDerivatives(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, int nFields, int offset, double *derivatives) const;

    void computeDerivativeLimited(const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int field, double *derivative) const;
    void computeDerivativesLimited(const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, double *derivatives) const;
    void computeDerivativesLimited(const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int nFields, double *derivatives) const;
    void computeDerivativesLimited(const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int nFields, int offset, double *derivatives) const;
    void computeDerivativeLimited(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int field, double *derivative) const;
    void computeDerivativesLimited(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, double *derivatives) const;
    void computeDerivativesLimited(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int nFields, double *derivatives) const;
    void computeDerivativesLimited(int degree, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, int nFields, int offset, double *derivatives) const;

    void computeGradient(const std::array<double, 3> &point, int field, std::array<double, 3> *gradient) const;
    void computeGradients(const std::array<double, 3> &point, std::array<double, 3> *gradients) const;
    void computeGradients(const std::array<double, 3> &point, int nFields, std::array<double, 3> *gradients) const;
    void computeGradients(const std::array<double, 3> &point, int nFields, int offset, std::array<double, 3> *gradients) const;
    void computeGradient(int degree, const std::array<double, 3> &point, int field, std::array<double, 3> *gradient) const;
    void computeGradients(int degree, const std::array<double, 3> &point, std::array<double, 3> *gradients) const;
    void computeGradients(int degree, const std::array<double, 3> &point, int nFields, std::array<double, 3> *gradients) const;
    void computeGradients(int degree, const std::array<double, 3> &point, int nFields, int offset, std::array<double, 3> *gradients) const;

    void computeGradientLimited(const std::array<double, 3> &point, const double *limiters, int field, std::array<double, 3> *gradient) const;
    void computeGradientsLimited(const std::array<double, 3> &point, const double *limiters, std::array<double, 3> *gradients) const;
    void computeGradientsLimited(const std::array<double, 3> &point, const double *limiters, int nFields, std::array<double, 3> *gradients) const;
    void computeGradientsLimited(const std::array<double, 3> &point, const double *limiters, int nFields, int offset, std::array<double, 3> *gradients) const;
    void computeGradientLimited(int degree, const std::array<double, 3> &point, const double *limiters, int field, std::array<double, 3> *gradient) const;
    void computeGradientsLimited(int degree, const std::array<double, 3> &point, const double *limiters, std::array<double, 3> *gradients) const;
    void computeGradientsLimited(int degree, const std::array<double, 3> &point, const double *limiters, int nFields, std::array<double, 3> *gradients) const;
    void computeGradientsLimited(int degree, const std::array<double, 3> &point, const double *limiters, int nFields, int offset, std::array<double, 3> *gradients) const;

    void display(std::ostream &out) const;

protected:
    static const uint8_t MAX_DEGREE;
    static const uint8_t MAX_DIMENSIONS;

    static uint16_t getCoefficientCount(uint8_t degree, uint8_t dimensions);
    static const std::vector<uint16_t> & getCoefficientsCount(uint8_t dimensions);
    static uint16_t countCoefficients(uint8_t degree, uint8_t dimensions);
    static const std::vector<uint16_t> & getDegreeCoefficientsCount(uint8_t dimensions);

    static uint16_t getDegreeCoefficientCount(uint8_t degree, uint8_t dimensions);
    static uint16_t countDegreeCoefficients(uint8_t degree, uint8_t dimensions);

    static void evalPointBasisValues(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin, const std::array<double, 3> &point, double *csi);
    static void evalPointBasisDerivatives(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *dcsi);

    static void evalCellBasisValues(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin, const Cell &cell, const std::array<double, 3> *vertexCoords, double *csi);

private:
    static const bool ENABLE_FAST_PATH_OPTIMIZATIONS;
    static const int MAX_STACK_WORKSPACE_SIZE;

    const static std::vector<std::vector<uint16_t>> m_countCoefficientCache;
    const static std::vector<std::vector<uint16_t>> m_countDegreeCoefficientCache;

    static std::vector<std::vector<uint16_t>> generateCountCoefficientCache();
    static std::vector<std::vector<uint16_t>> generateCountDegreeCoefficientCache();

    uint8_t m_degree;
    uint8_t m_dimensions;

    int m_nFields;

    uint16_t m_nCoeffs;
    std::unique_ptr<double[]> m_coeffs;

    std::array<double, 3> m_origin;

    std::size_t computeFieldCoefficientsOffset(uint8_t degree, int field) const;
    std::size_t getFieldCoefficientsStride() const;

};

class ReconstructionKernel {

public:
    ReconstructionKernel();
    ReconstructionKernel(uint8_t degree, uint8_t dimensions, int nEquations);

    ReconstructionKernel(const ReconstructionKernel &other);
    ReconstructionKernel(ReconstructionKernel &&other) = default;
    ReconstructionKernel & operator = (const ReconstructionKernel &other);
    ReconstructionKernel & operator=(ReconstructionKernel &&other) = default;

    void swap(ReconstructionKernel &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions, int nEquations, bool release = true);

    uint8_t getDegree() const;
    uint8_t getDimensions() const;

    uint16_t getCoefficientCount() const;

    int getEquationCount() const;

    const double * getPolynomialWeights() const;
    double * getPolynomialWeights();

    void assemblePolynomial(const std::array<double, 3> &origin, const double *values, ReconstructionPolynomial *polynomial) const;
    void assemblePolynomial(const std::array<double, 3> &origin, int nFields, const double **values, ReconstructionPolynomial *polynomial) const;
    void assemblePolynomial(uint8_t degree, const std::array<double, 3> &origin, const double *values, ReconstructionPolynomial *polynomial) const;
    void assemblePolynomial(uint8_t degree, const std::array<double, 3> &origin, int nFields, const double **values, ReconstructionPolynomial *polynomial) const;

    void updatePolynomial(const double *values, ReconstructionPolynomial *polynomial) const;
    void updatePolynomial(int nFields, const double **values, ReconstructionPolynomial *polynomial) const;
    void updatePolynomial(uint8_t degree, const double *values, ReconstructionPolynomial *polynomial) const;
    void updatePolynomial(uint8_t degree, int nFields, const double **values, ReconstructionPolynomial *polynomial) const;

    void computeValueWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, double *valueWeights) const;
    void computeValueWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, double *valueWeights) const;
    void computeValueLimitedWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiters, double *valueWeights) const;
    void computeValueLimitedWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiters, double *valueWeights) const;

    void computeDerivativeWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivativeWeights) const;
    void computeDerivativeWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivativeWeights) const;
    void computeDerivativeLimitedWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, double *derivativeWeights) const;
    void computeDerivativeLimitedWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, const double *limiters, double *derivativeWeights) const;

    void computeGradientWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, std::array<double, 3> *gradientWeights) const;
    void computeGradientWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, std::array<double, 3> *gradientWeights) const;
    void computeGradientLimitedWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiters, std::array<double, 3> *gradientWeights) const;
    void computeGradientLimitedWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiters, std::array<double, 3> *gradientWeights) const;

    void display(std::ostream &out, double tolerance = 1.e-10) const;

protected:
    void applyLimiter(uint8_t degree, const double *limiters, double *coeffs) const;

private:
    static const int MAX_STACK_WORKSPACE_SIZE;

    std::unique_ptr<double[]> m_weights;

    int m_nEquations;

    uint16_t m_nCoeffs;

    uint8_t m_degree;
    uint8_t m_dimensions;

};

class ReconstructionAssembler {

public:
    enum ReconstructionType {
        TYPE_CONSTRAINT,
        TYPE_LEAST_SQUARE
    };

    ReconstructionAssembler();
    ReconstructionAssembler(uint8_t degree, uint8_t dimensions);

    void swap(ReconstructionAssembler &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions, bool release = true);
    void clear(bool release = true);

    uint8_t getDegree() const;
    uint8_t getDimensions() const;

    uint16_t getCoefficientCount() const;

    int countConstraints() const;
    int countLeastSquares() const;
    int countEquations() const;

    void addPointValueEquation(ReconstructionType type, const std::array<double, 3> &origin, const std::array<double, 3> &point, double scaleFactor = 1.);
    void addPointDerivativeEquation(ReconstructionType type, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double scaleFactor = 1.);
    void addCellAverageEquation(ReconstructionType type, const Cell &cell, const std::array<double, 3> &origin, const std::array<double, 3> *vertexCoords, double scaleFactor = 1.);

    void assembleKernel(ReconstructionKernel *kernel) const;

    void updateKernel(ReconstructionKernel *kernel) const;

private:
    static double SVD_ZERO_THRESHOLD;

    uint8_t m_degree;
    uint8_t m_dimensions;

    uint16_t m_nCoeffs;

    std::vector<int> m_constraintsOrder;
    std::vector<int> m_leastSquaresOrder;

    std::vector<double> m_leastSquaresScaleFactors;

    std::vector<double> m_A;
    std::vector<double> m_C;

    mutable std::vector<double> m_sigma;
    mutable std::vector<double> m_U;
    mutable std::vector<double> m_S;
    mutable std::vector<double> m_Vt;
    mutable std::vector<double> m_SVDWorkspace;
    mutable std::vector<double> m_w;

    double * _addEquation(ReconstructionType type, double scaleFactor);

    void computePseudoInverse(int m, int n, double tolerance, double *A) const;

};

class Reconstruction : public ReconstructionKernel, public ReconstructionAssembler {

public:
    Reconstruction(uint8_t degree, uint8_t dimensions);

    void swap(Reconstruction &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions, bool release = true);
    void clear(bool release = true);

    using ReconstructionAssembler::getDegree;
    using ReconstructionAssembler::getDimensions;

    using ReconstructionAssembler::getCoefficientCount;

    void assemble();

};

}

#endif
