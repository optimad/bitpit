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

#include <cassert>
#include <limits>
#include <cblas.h>

#include "bitpit_IO.hpp"
#include "bitpit_LA.hpp"
#include "bitpit_operators.hpp"

#include "bitpit_private_lapacke.hpp"

#include "reconstruction.hpp"

namespace bitpit {

/*!
 * \class ReconstructionPolynomial
 * \ingroup discretization
 *
 * \brief The ReconstructionPolynomial class allows to apply a reconstruction
 * polinoymial previously assembled.
 */

/*!
 * Maximum allowed degree.
 */
const uint8_t ReconstructionPolynomial::MAX_DEGREE = 2;

/*!
 * Maximum allowed dimensions.
 */
const uint8_t ReconstructionPolynomial::MAX_DIMENSIONS = 3;

/*!
 * Enable fast calculation for some special cases.
 */
const bool ReconstructionPolynomial::ENABLE_FAST_PATH_OPTIMIZATIONS = true;

/*!
 * Maximum size of the workspace that will be allocated on the stack
 */
const int ReconstructionPolynomial::MAX_STACK_WORKSPACE_SIZE = 10;

/*!
 * Cache for the function that counts the coefficients.
 */
const std::vector<std::vector<uint16_t>> ReconstructionPolynomial::m_countCoefficientCache = generateCountCoefficientCache();

/*!
 * Cache for the function that counts the coefficients associated with a
 * specified degree.
 */
const std::vector<std::vector<uint16_t>> ReconstructionPolynomial::m_countDegreeCoefficientCache = generateCountDegreeCoefficientCache();

/*!
 * Generate ache for the function that counts the coefficients.
 */
std::vector<std::vector<uint16_t>> ReconstructionPolynomial::generateCountCoefficientCache()
{
    std::vector<std::vector<uint16_t>> cache(MAX_DIMENSIONS + 1, std::vector<uint16_t>(MAX_DEGREE + 1, 0));
    for (uint8_t dimensions = 0; dimensions <= MAX_DIMENSIONS; ++dimensions) {
        for (uint8_t degree = 0; degree <= MAX_DEGREE; ++degree) {
            cache[dimensions][degree] = countCoefficients(degree, dimensions);
        }
    }

    return cache;
}

/*!
 * Generate cache for the function that counts the coefficients associated with
 * a specified degree.
 */
std::vector<std::vector<uint16_t>> ReconstructionPolynomial::generateCountDegreeCoefficientCache()
{
    std::vector<std::vector<uint16_t>> cache(MAX_DIMENSIONS + 1, std::vector<uint16_t>(MAX_DEGREE + 1, 0));
    for (uint8_t dimensions = 0; dimensions <= MAX_DIMENSIONS; ++dimensions) {
        for (uint8_t degree = 0; degree <= MAX_DEGREE; ++degree) {
            cache[dimensions][degree] = countDegreeCoefficients(degree, dimensions);
        }
    }

    return cache;
}

/*!
 * Get the number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 *
 * \param degree is the degree
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 */
uint16_t ReconstructionPolynomial::getCoefficientCount(uint8_t degree, uint8_t dimensions)
{
    assert(degree <= ReconstructionPolynomial::MAX_DEGREE);
    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);

    return m_countCoefficientCache[dimensions][degree];
}

/*!
 * Get the number of coefficients of the reconstruction polynomial up to
 * the different degrees.
 *
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients of the reconstruction polynomial up to
 * the different degrees.
 */
const std::vector<uint16_t> & ReconstructionPolynomial::getCoefficientsCount(uint8_t dimensions)
{
    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);

    return m_countCoefficientCache[dimensions];
}

/*!
 * Count the number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 *
 * \param degree is the degree
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 */
uint16_t ReconstructionPolynomial::countCoefficients(uint8_t degree, uint8_t dimensions)
{
    uint16_t nCoeffs = 0;
    for (int i = 0; i <= degree; ++i) {
        nCoeffs += countDegreeCoefficients(i, dimensions);
    }

    return nCoeffs;
}

/*!
 * Get the number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 *
 * \param degree is the degree
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 */
uint16_t ReconstructionPolynomial::getDegreeCoefficientCount(uint8_t degree, uint8_t dimensions)
{
    assert(degree <= ReconstructionPolynomial::MAX_DEGREE);
    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);

    return m_countDegreeCoefficientCache[dimensions][degree];
}

/*!
 * Get the number of coefficients associates with each degree of the
 * reconstruction polynomial.
 *
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients associates with the each degree of the
 * reconstruction polynomial.
 */
const std::vector<uint16_t> & ReconstructionPolynomial::getDegreeCoefficientsCount(uint8_t dimensions)
{
    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);

    return m_countDegreeCoefficientCache[dimensions];
}

/*!
 * Count the number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 *
 * \param degree is the degree
 * \param dimensions is the number of space dimensions
 * \return The number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 */
uint16_t ReconstructionPolynomial::countDegreeCoefficients(uint8_t degree, uint8_t dimensions)
{
    if (dimensions == 0) {
        return 0;
    }

    return (utils::factorial(dimensions - 1 + degree) / utils::factorial(dimensions - 1) / utils::factorial(degree));
}

/*!
 * Evaluates the values of the basis for point reconstruction.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the basis will be evaluated
 * \param[out] csi on output will contain values of the basis for point
 * reconstruction.
 */
void ReconstructionPolynomial::evalPointBasisValues(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin,
                                                    const std::array<double, 3> &point, double *csi)
{
    // Set 0-th degree coefficients
    csi[0] = 1.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;
        const std::array<double, 3> distance = point - origin;

        // Set 1-st degree coefficients
        for (int i = 0; i < dimensions; ++i) {
            csi[offset++] = distance[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            for (int i = 0; i < dimensions; ++i) {
                csi[offset++] = 0.5 * distance[i] * distance[i];
            }

            if (dimensions >= 2) {
                csi[offset++] = distance[0] * distance[1];

                if (dimensions >= 3) {
                    csi[offset++] = distance[0] * distance[2];
                    csi[offset++] = distance[1] * distance[2];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == ReconstructionPolynomial::getCoefficientCount(degree, dimensions));
    }
}

/*!
 * Evaluates the derivatives of the basis for point reconstruction.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the basis will be evaluated
 * \param direction is the direction of the derivative
 * \param point is the point were the basis will be evaluated
 * \param[out] dcsi on output will contain values of the basis for point
 * derivative reconstruction.
 */
void ReconstructionPolynomial::evalPointBasisDerivatives(uint8_t degree, uint8_t dimensions,  const std::array<double, 3> &origin,
                                                         const std::array<double, 3> &point, const std::array<double, 3> &direction,
                                                         double *dcsi)
{
    // Set 0-th degree coefficients
    dcsi[0] = 0.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;

        // Set 1-st degree coefficients
        for (int i = 0; i < dimensions; ++i) {
            dcsi[offset++] = direction[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            const std::array<double, 3> distance = point - origin;

            for (int i = 0; i < dimensions; ++i) {
                dcsi[offset++] = distance[i] * direction[i];
            }

            if (dimensions >= 2) {
                dcsi[offset++] = distance[0] * direction[1] + distance[1] * direction[0];

                if (dimensions >= 3) {
                    dcsi[offset++] = distance[0] * direction[2] + distance[2] * direction[0];
                    dcsi[offset++] = distance[1] * direction[2] + distance[2] * direction[1];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == ReconstructionPolynomial::getCoefficientCount(degree, dimensions));
    }
}

/*!
 * Evaluates the values of the basis for cell reconstruction.
 *
 * The method works only for ElementType::Voxel and ElementType::Pixel.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param cell is the cell
 * \param origin is the point chosen as origin of the reconstruction
 * \param vertexCoords are the vertecx coordinates
 * \param[out] csi on output will contain the coefficients of the
 * cell average reconstruction
 */
void ReconstructionPolynomial::evalCellBasisValues(uint8_t degree, uint8_t dimensions, const std::array<double, 3> &origin,
                                                   const Cell &cell,  const std::array<double, 3> *vertexCoords,
                                                   double *csi)
{
    // Check if cell type is supported
    ElementType cellType = cell.getType();

    bool cellTypeSupported = false;
    if (cellType == ElementType::PIXEL) {
        cellTypeSupported = true;
    } else if (cellType == ElementType::VOXEL) {
        cellTypeSupported = true;
    }

    if (!cellTypeSupported) {
        throw std::runtime_error("Cell type not supported.");
    }

    // Set 0-th degree coefficients
    csi[0] = 1.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;
        const std::array<double, 3> distance = cell.evalCentroid(vertexCoords) - origin;

        // Set 1-st degree coefficients
        for (int i = 0; i < dimensions; ++i) {
            csi[offset++] = distance[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            double cellSize = cell.evalSize(vertexCoords);

            for (int i = 0; i < dimensions; ++i) {
                csi[offset++] = 0.5 * (distance[i] * distance[i] + cellSize * cellSize / 12.);
            }

            if (dimensions >= 2) {
                csi[offset++] = distance[0] * distance[1];

                if (dimensions >= 3) {
                    csi[offset++] = distance[0] * distance[2];
                    csi[offset++] = distance[1] * distance[2];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == ReconstructionPolynomial::getCoefficientCount(degree, dimensions));
    }
}

/*!
 * Constructor.
 */
ReconstructionPolynomial::ReconstructionPolynomial()
    : m_degree(0), m_dimensions(0), m_nFields(0), m_nCoeffs(0)
{
    initialize(0, 0, {{0., 0., 0.}}, 0);
}

/*!
 * Constructor.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param nFields is the number of fields associated with the polynomial
 */
ReconstructionPolynomial::ReconstructionPolynomial(uint8_t degree, uint8_t dimensions,
                                                   const std::array<double, 3> &origin, int nFields)
    : m_degree(0), m_dimensions(0), m_nFields(0), m_nCoeffs(0)
{
    initialize(degree, dimensions, origin, nFields);
}

/*!
    Copy constructor

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionPolynomial::ReconstructionPolynomial(const ReconstructionPolynomial &other)
    : ReconstructionPolynomial(other.getDegree(), other.getDimensions(), other.m_origin, other.m_nFields)
{
    if (m_nFields > 0) {
        int nWeights = m_nCoeffs * m_nFields;
        std::copy(other.m_coeffs.get(), other.m_coeffs.get() + nWeights, m_coeffs.get());
    }
}

/*!
    Copy-assignament operator.

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionPolynomial & ReconstructionPolynomial::operator=(const ReconstructionPolynomial &other)
{
    ReconstructionPolynomial tmp(other);
    swap(tmp);

    return *this;
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionPolynomial::swap(ReconstructionPolynomial &other) noexcept
{
    std::swap(other.m_degree, m_degree);
    std::swap(other.m_dimensions, m_dimensions);
    std::swap(other.m_nCoeffs, m_nCoeffs);
    std::swap(other.m_nFields, m_nFields);
    std::swap(other.m_coeffs, m_coeffs);
    std::swap(other.m_origin, m_origin);
}

/*!
 * Initialize the polynomial.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param origin is the point chosen as origin of the reconstruction
 * \param nFields is the number of fields associated with the polynomial
 * \param release if true, possible unneeded memory hold by the polynomial
 * will be released, otherwise the polynomial will be initialized but possible
 * unneeded memory will not be released
 */
void ReconstructionPolynomial::initialize(uint8_t degree, uint8_t dimensions,
                                          const std::array<double, 3> &origin,
                                          int nFields, bool release)
{
    assert(degree <= ReconstructionPolynomial::MAX_DEGREE);
    m_degree = degree;

    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);
    m_dimensions = dimensions;

    int currentStorageSize = m_nCoeffs;

    m_nCoeffs = ReconstructionPolynomial::getCoefficientCount(m_degree, m_dimensions);

    if (nFields > 0) {
        m_nFields = nFields;

        int storageSize = m_nCoeffs * m_nFields;

        bool reallocate;
        if (release) {
            reallocate = (currentStorageSize != storageSize);
        } else {
            reallocate = (currentStorageSize < storageSize);
        }

        if (reallocate) {
            m_coeffs = std::unique_ptr<double[]>(new double[storageSize]);
        }
    } else {
        clear(release);
    }

    m_origin = origin;
}

/*!
 * Clear the reconstruction.
 *
 * \param release if true, the memory hold by the polynomial will be released,
 * otherwise the polynomial will be cleared but its memory will not be released
 */
void ReconstructionPolynomial::clear(bool release)
{
    m_nFields = 0;

    if (release) {
        m_coeffs.reset();
    }
}

/*!
 * Get the degree of the polynomial.
 *
 * \return The degree of the polynomial.
 */
uint8_t ReconstructionPolynomial::getDegree() const
{
    return m_degree;
}

/*!
 * Get the number of space dimensions.
 *
 * \return The number of space dimensions.
 */
uint8_t ReconstructionPolynomial::getDimensions() const
{
    return m_dimensions;
}

/*!
 * Get the number of coefficients of the reconstruction polynomial.
 *
 * \return The number of coefficients of the reconstruction polynomial.
 */
uint16_t ReconstructionPolynomial::getCoefficientCount() const
{
    return m_nCoeffs;
}

/*!
 * Get the number of fields associated with the polynomial.
 *
 * \return The number of fields associated with the polynomial.
 */
int ReconstructionPolynomial::getFieldCount() const
{
    return m_nFields;
}

/*!
 * Returns a constant pointer to the coefficients of the polynomial.
 *
 * \result A constant pointer to the coefficients of the polynomial.
 */
const double * ReconstructionPolynomial::getCoefficients() const
{
    return m_coeffs.get();
}

/*!
 * Returns a pointer to the coefficients of the polynomial.
 *
 * \result A pointer to the coefficients of the polynomial.
 */
double * ReconstructionPolynomial::getCoefficients()
{
    return m_coeffs.get();
}

/*!
 * Returns a constant pointer to the coefficients of the polynomial.
 *
 * \param field is the field associated with the coefficients
 * \result A constant pointer to the coefficients of the polynomial.
 */
const double * ReconstructionPolynomial::getCoefficients(int field) const
{
    const double *coefficients = m_coeffs.get() + computeFieldCoefficientsOffset(0, field);

    return coefficients;
}

/*!
 * Returns a pointer to the coefficients of the polynomial.
 *
 * \param field is the field associated with the coefficients
 * \result A pointer to the coefficients of the polynomial.
 */
double * ReconstructionPolynomial::getCoefficients(int field)
{
    double *coefficients = m_coeffs.get() + computeFieldCoefficientsOffset(0, field);

    return coefficients;
}

/*!
 * Returns a constant pointer to the coefficients of the polynomial for the
 * specified field.
 *
 * \param degree is the requested degree
 * \param field is the field associated with the coefficients
 * \result A constant pointer to the coefficients of the polynomial for the
 * specified field.
 */
const double * ReconstructionPolynomial::getDegreeCoefficients(uint8_t degree, int field) const
{
    const double *coefficients = m_coeffs.get() + computeFieldCoefficientsOffset(degree, field);

    return coefficients;
}

/*!
 * Returns a pointer to the coefficients of the polynomial for the specified
 * field.
 *
 * \param degree is the requested degree
 * \param field is the field associated with the coefficients
 * \result A pointer to the coefficients of the polynomial for the specified
 * field.
 */
double * ReconstructionPolynomial::getDegreeCoefficients(uint8_t degree, int field)
{
    double *coefficients = m_coeffs.get() + computeFieldCoefficientsOffset(degree, field);

    return coefficients;
}

/*!
 * Gets the offset for accessing the coefficient of the given order for the
 * specified field.
 *
 * \param degree is the requested degree
 * \param field is the field associated with the coefficients
 * \result The offset for accessing the coefficient of the given order for the
 * specified field.
 */
std::size_t ReconstructionPolynomial::computeFieldCoefficientsOffset(uint8_t degree, int field) const
{
    std::size_t offset = field * getFieldCoefficientsStride() + degree;

    return offset;
}

/*!
 * Gets the offset for accessing the coefficient of the given order for the
 * specified field.
 *
 * \param degree is the requested degree
 * \param field is the field associated with the coefficients
 * \result The offset for accessing the coefficient of the given order for the
 * specified field.
 */
std::size_t ReconstructionPolynomial::getFieldCoefficientsStride() const
{
    return m_nCoeffs;
}

/*!
 * Computes the value of the polynomial at the specified point.
 *
 * \param point is the point were the polynomial will be evaluated
 * \param field is the field that will be evaluated
 * \param[out] value on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValue(const std::array<double, 3> &point,
                                            int field, double *value) const
{
    computeValues(m_degree, point, 1, field, value);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param point is the point were the polynomial will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(const std::array<double, 3> &point,
                                             double *values) const
{
    computeValues(m_degree, point, m_nFields, 0, values);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param point is the point were the polynomial will be evaluated
 * \param nFields is the number of fields that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(const std::array<double, 3> &point,
                                             int nFields, double *values) const
{
    computeValues(m_degree, point, nFields, 0, values);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param point is the point were the polynomial will be evaluated
 * \param nFields is the number of fields that will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(const std::array<double, 3> &point,
                                             int nFields, int offset, double *values) const
{
    computeValues(m_degree, point, nFields, offset, values);
}

/*!
 * Computes the value of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param field is the field that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValue(int degree, const std::array<double, 3> &point,
                                            int field, double *values) const
{
    computeValues(degree, point, 1, field, values);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(int degree, const std::array<double, 3> &point,
                                             double *values) const
{
    computeValues(degree, point, m_nFields, 0, values);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param nFields is the number of fields that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(int degree, const std::array<double, 3> &point,
                                             int nFields, double *values) const
{
    computeValues(degree, point, nFields, 0, values);
}

/*!
 * Computes the values of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param nFields is the number of fields that will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValues(int degree, const std::array<double, 3> &point,
                                             int nFields, int offset, double *values) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get values
    const std::size_t fieldValuesStride = 1;
    double *fieldValue = values;
    const double *fieldValueEnd = fieldValue + fieldValuesStride * nFields;

    // Get coefficients
    std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldValue = fieldCoeffs[0];

            fieldCoeffs += fieldCoeffsStride;
            fieldValue  += fieldValuesStride;
        } while (fieldValue != fieldValueEnd);

        return;
    }

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisValues(degree, m_dimensions, m_origin, point, csi);

    do {
        *fieldValue = fieldCoeffs[0] * csi[0];
        for (int i = 1; i < nCoeffs; ++i) {
            *fieldValue += fieldCoeffs[i] * csi[i];
        }

        fieldCoeffs += fieldCoeffsStride;
        fieldValue  += fieldValuesStride;
    } while (fieldValue != fieldValueEnd);
}

/*!
 * Computes the limited value of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param field is the field that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValueLimited(const std::array<double, 3> &point,
                                                   const double *limiters, int field,
                                                   double *values) const
{
    computeValuesLimited(m_degree, point, limiters, 1, field, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(const std::array<double, 3> &point,
                                                    const double *limiters,
                                                    double *values) const
{
    computeValuesLimited(m_degree, point, limiters, m_nFields, 0, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of fields that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(const std::array<double, 3> &point,
                                                    const double *limiters, int nFields,
                                                    double *values) const
{
    computeValuesLimited(m_degree, point, limiters, nFields, 0, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of fields that will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(const std::array<double, 3> &point,
                                                    const double *limiters, int nFields, int offset,
                                                    double *values) const
{
    computeValuesLimited(m_degree, point, limiters, nFields, offset, values);
}

/*!
 * Computes the limited value of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param field is the field that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValueLimited(int degree, const std::array<double, 3> &point,
                                                   const double *limiters, int field,
                                                   double *values) const
{
    computeValuesLimited(degree, point, limiters, 1, field, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(int degree, const std::array<double, 3> &point,
                                                    const double *limiters,
                                                    double *values) const
{
    computeValuesLimited(degree, point, limiters, m_nFields, 0, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of fields that will be evaluated
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(int degree, const std::array<double, 3> &point,
                                                    const double *limiters, int nFields,
                                                    double *values) const
{
    computeValuesLimited(degree, point, limiters, nFields, 0, values);
}

/*!
 * Computes the limited values of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the polynomial will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of fields that will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] values on output will contain the values of the polynomial
 */
void ReconstructionPolynomial::computeValuesLimited(int degree, const std::array<double, 3> &point,
                                                    const double *limiters, int nFields, int offset,
                                                    double *values) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get values
    std::size_t fieldValuesStride = 1;
    double *fieldValue = values;
    double *fieldValueEnd = fieldValue + fieldValuesStride * nFields;

    // Get coefficients
    std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldValue = fieldCoeffs[0];

            fieldCoeffs += fieldCoeffsStride;
            fieldValue  += fieldValuesStride;
        } while (fieldValue != fieldValueEnd);

        return;
    }

    // Get limiters
    std::size_t fieldLimitersStride = degree;
    const double *fieldLimiters = limiters;

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);
    const std::vector<uint16_t> &nDegreeCoeffs = ReconstructionPolynomial::getDegreeCoefficientsCount(m_dimensions);

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisValues(degree, m_dimensions, m_origin, point, csi);

    do {
        // Degree 0
        *fieldValue = fieldCoeffs[0] * csi[0];

        // Degrees greater than 0
        int coeffEnd = nDegreeCoeffs[0];
        for (int n = 1; n <= degree; ++n) {
            int coeffBegin = coeffEnd;
            coeffEnd = coeffBegin + nDegreeCoeffs[n];

            double fieldDegreeValue = 0;
            for (int i = coeffBegin; i < coeffEnd; ++i) {
                fieldDegreeValue += fieldCoeffs[i] * csi[i];
            }

            *fieldValue += fieldDegreeValue * fieldLimiters[n - 1];
        }

        fieldCoeffs   += fieldCoeffsStride;
        fieldValue    += fieldValuesStride;
        fieldLimiters += fieldLimitersStride;
    } while (fieldValue != fieldValueEnd);
}

/*!
 * Computes the derivative of the polynomial at the specified point.
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivative
 * \param field is the field for which derivative will be evaluated
 * \param[out] derivative on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivative(const std::array<double, 3> &point,
                                                 const std::array<double, 3> &direction,
                                                 int field, double *derivative) const
{
    computeDerivatives(m_degree, point, direction, 1, field, derivative);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  double *derivatives) const
{
    computeDerivatives(m_degree, point, direction, m_nFields, 0, derivatives);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  int nFields, double *derivatives) const
{
    computeDerivatives(m_degree, point, direction, nFields, 0, derivatives);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  int nFields, int offset, double *derivatives) const
{
    computeDerivatives(m_degree, point, direction, nFields, offset, derivatives);
}

/*!
 * Computes the derivative of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivative
 * \param field is the field for which derivative will be evaluated
 * \param[out] derivative on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivative(int degree, const std::array<double, 3> &point,
                                                 const std::array<double, 3> &direction,
                                                 int field, double *derivative) const
{
    computeDerivatives(degree, point, direction, 1, field, derivative);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(int degree, const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  double *derivatives) const
{
    computeDerivatives(degree, point, direction, m_nFields, 0, derivatives);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(int degree, const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  int nFields, double *derivatives) const
{
    computeDerivatives(degree, point, direction, nFields, 0, derivatives);
}

/*!
 * Computes the derivatives of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivatives(int degree, const std::array<double, 3> &point,
                                                  const std::array<double, 3> &direction,
                                                  int nFields, int offset, double *derivatives) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get derivatives
    std::size_t fieldDerivativeStride = 1;
    double *fieldDerivative = derivatives;
    double *fieldDerivativeEnd = fieldDerivative + fieldDerivativeStride * nFields;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldDerivative = 0.;

            fieldDerivative += fieldDerivativeStride;
        } while (fieldDerivative != fieldDerivativeEnd);

        return;
    }

    // Get coefficients
    std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Linear polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 1) {
        do {
            *fieldDerivative = fieldCoeffs[1] * direction[0];
            for (int d = 1; d < m_dimensions; ++d) {
                *fieldDerivative += fieldCoeffs[d + 1] * direction[d];
            }

            fieldCoeffs     += fieldCoeffsStride;
            fieldDerivative += fieldDerivativeStride;
        } while (fieldDerivative != fieldDerivativeEnd);

        return;
    }

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);
    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, direction, csi);

    do {
        *fieldDerivative = fieldCoeffs[0] * csi[0];
        for (int i = 1; i < nCoeffs; ++i) {
            *fieldDerivative += fieldCoeffs[i] * csi[i];
        }

        fieldCoeffs     += fieldCoeffsStride;
        fieldDerivative += fieldDerivativeStride;
    } while (fieldDerivative != fieldDerivativeEnd);
}

/*!
 * Computes the limited derivative of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivative
 * \param limiters are the scale factors of the limiters
 * \param field is the field for which derivative will be evaluated
 * \param[out] derivative on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativeLimited(const std::array<double, 3> &point,
                                                        const std::array<double, 3> &direction,
                                                        const double *limiters, int field,
                                                        double *derivative) const
{
    computeDerivativesLimited(m_degree, point, direction, limiters, 1, field, derivative);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters,
                                                         double *derivatives) const
{
    computeDerivativesLimited(m_degree, point, direction, limiters, m_nFields, 0, derivatives);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters, int nFields,
                                                         double *derivatives) const
{
    computeDerivativesLimited(m_degree, point, direction, limiters, nFields, 0, derivatives);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
  *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters, int nFields, int offset,
                                                         double *derivatives) const
{
    computeDerivativesLimited(m_degree, point, direction, limiters, nFields, offset, derivatives);
}

/*!
 * Computes the limited derivative of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivative
 * \param limiters are the scale factors of the limiters
 * \param field is the field for which derivative will be evaluated
 * \param[out] derivative on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativeLimited(int degree, const std::array<double, 3> &point,
                                                        const std::array<double, 3> &direction,
                                                        const double *limiters, int field,
                                                        double *derivative) const
{
    computeDerivativesLimited(degree, point, direction, limiters, 1, field, derivative);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(int degree, const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters,
                                                         double *derivatives) const
{
    computeDerivativesLimited(degree, point, direction, limiters, m_nFields, 0, derivatives);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(int degree, const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters, int nFields,
                                                         double *derivatives) const
{
    computeDerivativesLimited(degree, point, direction, limiters, nFields, 0, derivatives);
}

/*!
 * Computes the limited derivatives of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the derivatives will be evaluated
 * \param direction is the direction of the derivatives
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which derivatives will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] derivatives on output will contain the derivatives
 */
void ReconstructionPolynomial::computeDerivativesLimited(int degree, const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         const double *limiters, int nFields, int offset,
                                                         double *derivatives) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get derivatives
    std::size_t fieldDerivativeStride = 1;
    double *fieldDerivative = derivatives;
    double *fieldDerivativeEnd = fieldDerivative + fieldDerivativeStride * nFields;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldDerivative = 0.;

            fieldDerivative += fieldDerivativeStride;
        } while (fieldDerivative != fieldDerivativeEnd);

        return;
    }

    // Get coefficients
    std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Get limiters
    std::size_t fieldLimitersStride = degree;
    const double *fieldLimiters = limiters;

    // Linear polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 1) {
        do {
            *fieldDerivative = fieldCoeffs[1] * direction[0];
            for (int d = 1; d < m_dimensions; ++d) {
                *fieldDerivative += fieldCoeffs[d + 1] * direction[d];
            }

            *fieldDerivative *= fieldLimiters[0];

            fieldCoeffs     += fieldCoeffsStride;
            fieldDerivative += fieldDerivativeStride;
            fieldLimiters   += fieldLimitersStride;
        } while (fieldDerivative != fieldDerivativeEnd);

        return;
    }

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);
    const std::vector<uint16_t> &nDegreeCoeffs = ReconstructionPolynomial::getDegreeCoefficientsCount(m_dimensions);

    BITPIT_CREATE_WORKSPACE(dcsi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, direction, dcsi);

    do {
        *fieldDerivative = fieldCoeffs[0] * dcsi[0];

        int coeffEnd = nDegreeCoeffs[0];
        for (int n = 1; n <= degree; ++n) {
            int coeffBegin = coeffEnd;
            coeffEnd = coeffBegin + nDegreeCoeffs[n];

            double fieldDegreeDerivative = 0;
            for (int i = coeffBegin; i < coeffEnd; ++i) {
                fieldDegreeDerivative += fieldCoeffs[i] * dcsi[i];
            }

            *fieldDerivative += fieldLimiters[n - 1] * fieldDegreeDerivative;
        }

        fieldCoeffs     += fieldCoeffsStride;
        fieldDerivative += fieldDerivativeStride;
        fieldLimiters   += fieldLimitersStride;
    } while (fieldDerivative != fieldDerivativeEnd);
}

/*!
 * Computes the gradient of the polynomial at the specified point.
 *
 * \param point is the point were the gradient will be evaluated
 * \param field is the field for which gradient will be evaluated
 * \param[out] gradient on output will contain the gradient
 */
void ReconstructionPolynomial::computeGradient(const std::array<double, 3> &point,
                                               int field, std::array<double, 3> *gradient) const
{
    computeGradients(m_degree, point, 1, field, gradient);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param point is the point were the gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(const std::array<double, 3> &point,
                                                std::array<double, 3> *gradients) const
{
    computeGradients(m_degree, point, m_nFields, 0, gradients);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param point is the point were the gradients will be evaluated
 * \param nFields is the number of field for which gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(const std::array<double, 3> &point,
                                                int nFields, std::array<double, 3> *gradients) const
{
    computeGradients(m_degree, point, nFields, 0, gradients);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param point is the point were the gradients will be evaluated
 * \param nFields is the number of field for which gradients will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(const std::array<double, 3> &point,
                                                int nFields, int offset,
                                                std::array<double, 3> *gradients) const
{
    computeGradients(m_degree, point, nFields, offset, gradients);
}

/*!
 * Computes the gradient of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradient will be evaluated
 * \param field is the field for which gradient will be evaluated
 * \param[out] gradient on output will contain the gradient
 */
void ReconstructionPolynomial::computeGradient(int degree, const std::array<double, 3> &point,
                                               int field, std::array<double, 3> *gradient) const
{
    computeGradients(degree, point, 1, field, gradient);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(int degree, const std::array<double, 3> &point,
                                                std::array<double, 3> *gradients) const
{
    computeGradients(degree, point, m_nFields, 0, gradients);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param nFields is the number of field for which gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(int degree, const std::array<double, 3> &point,
                                                int nFields, std::array<double, 3> *gradients) const
{
    computeGradients(degree, point, nFields, 0, gradients);
}

/*!
 * Computes the gradients of the polynomial at the specified point.
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param nFields is the number of field for which gradients will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradients(int degree, const std::array<double, 3> &point,
                                                int nFields, int offset,
                                                std::array<double, 3> *gradients) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get gradients
    const std::size_t fieldGradientStride = 1;
    std::array<double, 3> *fieldGradient = gradients;
    const std::array<double, 3> *fieldGradientEnd = fieldGradient + fieldGradientStride * nFields;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldGradient = {{0., 0., 0.}};

            fieldGradient += fieldGradientStride;
        } while (fieldGradient != fieldGradientEnd);

        return;
    }

    // Get coefficients
    const std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Linear polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 1) {
        do {
            // Evaluate gradients
            for (int d = 0; d < m_dimensions; ++d) {
                (*fieldGradient)[d] = fieldCoeffs[1 + d];
            }

            // Explicitly zero unused components
            for (int d = m_dimensions; d < ReconstructionPolynomial::MAX_DIMENSIONS; ++d) {
                (*fieldGradient)[d] = 0.;
            }

            // Advance to the next field
            fieldCoeffs   += fieldCoeffsStride;
            fieldGradient += fieldGradientStride;
        } while (fieldGradient != fieldGradientEnd);

        return;
    }

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);

    BITPIT_CREATE_WORKSPACE(dcsi, double, static_cast<std::size_t>(m_dimensions * nCoeffs), 3 * MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{1., 0., 0.}}, dcsi);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{0., 1., 0.}}, dcsi + nCoeffs);
    if (m_dimensions == 3) {
        ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{0., 0., 1.}}, dcsi + 2 * nCoeffs);
    }

    do {
        // Evaluate gradient
        for (int d = 0; d < m_dimensions; ++d) {
            const double *dcsi_dimension = dcsi + d * nCoeffs;

            (*fieldGradient)[d] = fieldCoeffs[0] * dcsi_dimension[0];
            for (int i = 1; i < nCoeffs; ++i) {
                (*fieldGradient)[d] += fieldCoeffs[i] * dcsi_dimension[i];
            }
        }

        // Explicitly zero unused components
        for (int d = m_dimensions; d < ReconstructionPolynomial::MAX_DIMENSIONS; ++d) {
            (*fieldGradient)[d] = 0.;
        }

        // Advance to the next field
        fieldCoeffs   += fieldCoeffsStride;
        fieldGradient += fieldGradientStride;
    } while (fieldGradient != fieldGradientEnd);
}

/*!
 * Computes the limited gradient of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the gradient will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param field is the field for which gradient will be evaluated
 * \param[out] gradient on output will contain the gradient
 */
void ReconstructionPolynomial::computeGradientLimited(const std::array<double, 3> &point,
                                                      const double *limiters, int field,
                                                      std::array<double, 3> *gradient) const
{
    computeGradientsLimited(m_degree, point, limiters, 1, field, gradient);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(const std::array<double, 3> &point,
                                                       const double *limiters,
                                                       std::array<double, 3> *gradients) const
{
    computeGradientsLimited(m_degree, point, limiters, m_nFields, 0, gradients);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(const std::array<double, 3> &point,
                                                       const double *limiters, int nFields,
                                                       std::array<double, 3> *gradients) const
{
    computeGradientsLimited(m_degree, point, limiters, nFields, 0, gradients);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which gradients will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(const std::array<double, 3> &point,
                                                       const double *limiters, int nFields, int offset,
                                                       std::array<double, 3> *gradients) const
{
    computeGradientsLimited(m_degree, point, limiters, nFields, offset, gradients);
}

/*!
 * Computes the limited gradient of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param limiters are the scale factors of the limiters
 * \param point is the point were the gradient will be evaluated
 * \param field is the field for which gradient will be evaluated
 * \param[out] gradient on output will contain the gradient
 */
void ReconstructionPolynomial::computeGradientLimited(int degree, const std::array<double, 3> &point,
                                                      const double *limiters, int field,
                                                      std::array<double, 3> *gradient) const
{
    computeGradientsLimited(degree, point, limiters, 1, field, gradient);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(int degree, const std::array<double, 3> &point,
                                                       const double *limiters,
                                                       std::array<double, 3> *gradients) const
{
    computeGradientsLimited(degree, point, limiters, m_nFields, 0, gradients);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which gradients will be evaluated
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(int degree, const std::array<double, 3> &point,
                                                       const double *limiters, int nFields,
                                                       std::array<double, 3> *gradients) const
{
    computeGradientsLimited(degree, point, limiters, nFields, 0, gradients);
}

/*!
 * Computes the limited gradients of the polynomial at the specified point.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param point is the point were the gradients will be evaluated
 * \param limiters are the scale factors of the limiters
 * \param nFields is the number of field for which gradients will be evaluated
 * \param offset is the offset that will be used when reading field data
 * \param[out] gradients on output will contain the gradients
 */
void ReconstructionPolynomial::computeGradientsLimited(int degree, const std::array<double, 3> &point,
                                                       const double *limiters, int nFields, int offset,
                                                       std::array<double, 3> *gradients) const
{
    assert(degree <= getDegree());

    // Early return if there are no field to process
    if (nFields == 0) {
        return;
    }

    // Get gradients
    const std::size_t fieldGradientStride = 1;
    std::array<double, 3> *fieldGradient = gradients;
    const std::array<double, 3> *fieldGradientEnd = fieldGradient + fieldGradientStride * nFields;

    // Constant polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 0) {
        do {
            *fieldGradient = {{0., 0., 0.}};

            fieldGradient += fieldGradientStride;
        } while (fieldGradient != fieldGradientEnd);

        return;
    }

    // Get coefficients
    const std::size_t fieldCoeffsStride = getFieldCoefficientsStride();
    const double *fieldCoeffs = m_coeffs.get() + fieldCoeffsStride * offset;

    // Get limiters
    const std::size_t fieldLimitersStride = degree;
    const double *fieldLimiters = limiters;

    // Linear polynomial
    if (ENABLE_FAST_PATH_OPTIMIZATIONS && degree == 1) {
        do {
            // Evaluate gradients
            for (int d = 0; d < m_dimensions; ++d) {
                (*fieldGradient)[d] = fieldLimiters[0] * fieldCoeffs[1 + d];
            }

            // Explicitly zero unused components
            for (int d = m_dimensions; d < ReconstructionPolynomial::MAX_DIMENSIONS; ++d) {
                (*fieldGradient)[d] = 0.;
            }

            // Advance to the next field
            fieldCoeffs   += fieldCoeffsStride;
            fieldGradient += fieldGradientStride;
            fieldLimiters += fieldLimitersStride;
        } while (fieldGradient != fieldGradientEnd);

        return;
    }

    // Generic polynomial
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, m_dimensions);
    const std::vector<uint16_t> &nDegreeCoeffs = ReconstructionPolynomial::getDegreeCoefficientsCount(m_dimensions);

    BITPIT_CREATE_WORKSPACE(dcsi, double, static_cast<std::size_t>(m_dimensions * nCoeffs), 3 * MAX_STACK_WORKSPACE_SIZE);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{1., 0., 0.}}, dcsi);
    ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{0., 1., 0.}}, dcsi + nCoeffs);
    if (m_dimensions == 3) {
        ReconstructionPolynomial::evalPointBasisDerivatives(degree, m_dimensions, m_origin, point, {{0., 0., 1.}}, dcsi + 2 * nCoeffs);
    }

    do {
        // Evaluate gradients
        for (int d = 0; d < m_dimensions; ++d) {
            const double *dcsi_dimension = dcsi + d * nCoeffs;

            (*fieldGradient)[d] = fieldCoeffs[0] * dcsi_dimension[0];
        }

        int coeffEnd = nDegreeCoeffs[0];
        for (int n = 1; n <= degree; ++n) {
            double fieldsLimiter = fieldLimiters[n - 1];

            int coeffBegin = coeffEnd;
            coeffEnd = coeffBegin + nDegreeCoeffs[n];
            for (int i = coeffBegin; i < coeffEnd; ++i) {
                for (int d = 0; d < m_dimensions; ++d) {
                    const double *dcsi_dimension = dcsi + d * nCoeffs;

                    (*fieldGradient)[d] += fieldsLimiter * fieldCoeffs[i] * dcsi_dimension[i];
                }
            }
        }

        // Explicitly zero unused components
        for (int d = m_dimensions; d < ReconstructionPolynomial::MAX_DIMENSIONS; ++d) {
            (*fieldGradient)[d] = 0.;
        }

        // Advance to the next field
        fieldCoeffs   += fieldCoeffsStride;
        fieldGradient += fieldGradientStride;
        fieldLimiters += fieldLimitersStride;
    } while (fieldGradient != fieldGradientEnd);
}

/*!
 * Displays polynomial information to an output stream.
 *
 * \param out is the output stream
 */
void ReconstructionPolynomial::display(std::ostream &out) const
{
    uint8_t dimensions = getDimensions();

    for (int k = 0; k < m_nFields; ++k) {
        out << " field " << k << "\n";
        for (int i = 0; i <= m_degree; ++i) {
            out << "   degree = " << i << " : " ;

            uint16_t nDegreeCoeffs = ReconstructionPolynomial::getDegreeCoefficientCount(i, dimensions);
            const double *degreeCoeffs = getDegreeCoefficients(i, k);
            for (int n = 0; n < nDegreeCoeffs; ++n) {
                out << degreeCoeffs[n];
                if (n != nDegreeCoeffs - 1) {
                    out << " , ";
                }
            }

            out << "\n";
        }
    }
}

/*!
 * \class ReconstructionKernel
 * \ingroup discretization
 *
 * \brief The ReconstructionKernel class allows to evaluate the weight of
 * a reconstruction polinoymial previously assembled.
 */

/*!
 * Maximum size of the workspace that will be allocated on the stack
 */
const int ReconstructionKernel::MAX_STACK_WORKSPACE_SIZE = 10;

/*!
 * Constructor.
 */
ReconstructionKernel::ReconstructionKernel()
    : m_nEquations(0), m_nCoeffs(0), m_degree(0), m_dimensions(0)
{
    initialize(0, 0, 0);
}

/*!
 * Constructor.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param nEquations is the number of equations that defines the reconstruction
 */
ReconstructionKernel::ReconstructionKernel(uint8_t degree, uint8_t dimensions, int nEquations)
    : m_nEquations(0), m_nCoeffs(0), m_degree(0), m_dimensions(0)
{
    initialize(degree, dimensions, nEquations);
}

/*!
    Copy constructor

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionKernel::ReconstructionKernel(const ReconstructionKernel &other)
    : ReconstructionKernel(other.getDegree(), other.getDimensions(), other.m_nEquations)
{
    if (m_nEquations > 0) {
        int nWeights = m_nCoeffs * m_nEquations;
        std::copy(other.m_weights.get(), other.m_weights.get() + nWeights, m_weights.get());
    }
}

/*!
    Copy-assignament operator.

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionKernel & ReconstructionKernel::operator=(const ReconstructionKernel &other)
{
    ReconstructionKernel tmp(other);
    swap(tmp);

    return *this;
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionKernel::swap(ReconstructionKernel &other) noexcept
{
    std::swap(other.m_degree, m_degree);
    std::swap(other.m_dimensions, m_dimensions);
    std::swap(other.m_nCoeffs, m_nCoeffs);
    std::swap(other.m_nEquations, m_nEquations);
    std::swap(other.m_weights, m_weights);
}

/*!
 * Initialize the kernel.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param nEquations is the number of equations that defines the reconstruction
 * \param release if true, possible unneeded memory hold by the kernel
 * will be released, otherwise the kernel will be initialized but possible
 * unneeded memory will not be released
 */
void ReconstructionKernel::initialize(uint8_t degree, uint8_t dimensions, int nEquations, bool release)
{
    assert(degree <= ReconstructionPolynomial::MAX_DEGREE);
    m_degree = degree;

    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);
    m_dimensions = dimensions;

    int currentStorageSize = m_nCoeffs * m_nEquations;

    m_nCoeffs    = ReconstructionPolynomial::getCoefficientCount(m_degree, m_dimensions);
    m_nEquations = nEquations;

    int storageSize = m_nCoeffs * m_nEquations;

    bool reallocate;
    if (release) {
        reallocate = (currentStorageSize != storageSize);
    } else {
        reallocate = (currentStorageSize < storageSize);
    }

    if (reallocate) {
        m_weights = std::unique_ptr<double[]>(new double[storageSize]);
    }
}

/*!
 * Get the degree of the polynomial.
 *
 * \return The degree of the polynomial.
 */
uint8_t ReconstructionKernel::getDegree() const
{
    return m_degree;
}

/*!
 * Get the number of space dimensions.
 *
 * \return The number of space dimensions.
 */
uint8_t ReconstructionKernel::getDimensions() const
{
    return m_dimensions;
}

/*!
 * Get the number of coefficients of the reconstruction polynomial.
 *
 * \return The number of coefficients of the reconstruction polynomial.
 */
uint16_t ReconstructionKernel::getCoefficientCount() const
{
    return m_nCoeffs;
}

/*!
 * Get the number of equations associated with the reconstruction.
 *
 * \return The number of equations associated with the reconstruction.
 */
int ReconstructionKernel::getEquationCount() const
{
    return m_nEquations;
}

/*!
 * Returns a constant pointer to the internal weight that are used to evaluate
 * the polynomial coefficients.
 *
 * \result A constant pointer to the internal weight that are used to evaluate
 * the polynomial coefficients.
 */
const double * ReconstructionKernel::getPolynomialWeights() const
{
    return m_weights.get();
}

/*!
 * Returns a pointer to the internal weight that are used to evaluate the
 * polynomial coefficients.
 *
 * \result A pointer to the internal weight that are used to evaluate the
 * polynomial coefficients.
 */
double * ReconstructionKernel::getPolynomialWeights()
{
    return m_weights.get();
}

/*!
 * Assembles the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will be properly
 * initialized.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param values are the values associated with the equations
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::assemblePolynomial(const std::array<double, 3> &origin,
                                              const double *values,
                                              ReconstructionPolynomial *polynomial) const
{
    uint8_t degree = getDegree();
    uint8_t dimensions = getDimensions();

    // Initialize the polynomial
    polynomial->initialize(degree, dimensions, origin, 1, true);

    // Update the polynomial
    updatePolynomial(degree, values, polynomial);
}

/*!
 * Assembles the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will be properly
 * initialized and possible unneeded memory hold by the polynomial will be
 * released.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param nFields is the number of fields associated with the polynomial
 * \param values are the values associated with the equations, values are
 * grouped by equation, this means that the first element will contain
 * all the values associated with the first equation, the second element
 * will contain the values associated with the second equation, and so on
 * and so forth
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::assemblePolynomial(const std::array<double, 3> &origin,
                                              int nFields, const double **values,
                                              ReconstructionPolynomial *polynomial) const
{
    uint8_t degree = getDegree();
    uint8_t dimensions = getDimensions();

    // Initialize the polynomial
    polynomial->initialize(degree, dimensions, origin, nFields, true);

    // Update the polynomial
    updatePolynomial(degree, nFields, values, polynomial);
}

/*!
 * Assembles the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will be properly
 * initialized.
 *
 * The degree of the polynomial cannot be higher than the degree of the
 * kernel.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param values are the values associated with the equations
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::assemblePolynomial(uint8_t degree, const std::array<double, 3> &origin,
                                              const double *values,
                                              ReconstructionPolynomial *polynomial) const
{
    uint8_t dimensions = getDimensions();

    // Initialize the polynomial
    polynomial->initialize(degree, dimensions, origin, 1, true);

    // Update the polynomial
    updatePolynomial(degree, values, polynomial);
}

/*!
 * Assembles the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will be properly
 * initialized.
 *
 * The degree of the polynomial cannot be higher than the degree of the
 * kernel.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param values are the values associated with the equations, values are
 * grouped by equation, this means that the first element will contain
 * all the values associated with the first equation, the second element
 * will contain the values associated with the second equation, and so on
 * and so forth
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::assemblePolynomial(uint8_t degree, const std::array<double, 3> &origin,
                                              int nFields, const double **values,
                                              ReconstructionPolynomial *polynomial) const
{
    uint8_t dimensions = getDimensions();

    // Initialize the polynomial
    polynomial->initialize(degree, dimensions, origin, nFields, true);

    // Update the polynomial
    updatePolynomial(degree, nFields, values, polynomial);
}

/*!
 * Updates the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will not be
 * initialized.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param values are the values associated with the equations
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::updatePolynomial(const double *values,
                                            ReconstructionPolynomial *polynomial) const
{
    updatePolynomial(getDegree(), values, polynomial);
}

/*!
 * Updates the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will not be
 * initialized.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param values are the values associated with the equations, values are
 * grouped by equation, this means that the first element will contain
 * all the values associated with the first equation, the second element
 * will contain the values associated with the second equation, and so on
 * and so forth
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::updatePolynomial(int nFields, const double **values,
                                            ReconstructionPolynomial *polynomial) const
{
    updatePolynomial(getDegree(), nFields, values, polynomial);
}

/*!
 * Updates the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will not be
 * initialized.
 *
 * The degree of the polynomial cannot be higher than the degree of the
 * kernel.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param degree is the degree of the polynomial
 * \param values are the values associated with the equations
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::updatePolynomial(uint8_t degree, const double *values,
                                            ReconstructionPolynomial *polynomial) const
{
    assert(degree <= getDegree());

    uint8_t dimensions = getDimensions();

    int nEquations = getEquationCount();
    int nCoeffs    = ReconstructionPolynomial::getCoefficientCount(degree, dimensions);

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nEquations, nCoeffs, 1., getPolynomialWeights(), nEquations,
                values, 1, 0, polynomial->getCoefficients(), 1);
}

/*!
 * Updates the polynomial for the specified values.
 *
 * Before computing polynomial coefficients, the polynomial will not be
 * initialized.
 *
 * The degree of the polynomial cannot be higher than the degree of the
 * kernel.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param degree is the degree of the polynomial
 * \param nFields is the number of fields associated with the polynomial
 * \param values are the values associated with the equations, values are
 * grouped by equation, this means that the first element will contain
 * all the values associated with the first equation, the second element
 * will contain the values associated with the second equation, and so on
 * and so forth
 * \param[out] polynomial on output will contain the polynomial
 */
void ReconstructionKernel::updatePolynomial(uint8_t degree, int nFields, const double **values,
                                            ReconstructionPolynomial *polynomial) const
{
    assert(degree <= getDegree());

    uint8_t dimensions = getDimensions();

    int nEquations = getEquationCount();
    int nCoeffs    = ReconstructionPolynomial::getCoefficientCount(degree, dimensions);

    const double *weights = m_weights.get();
    double *coeffs = polynomial->getCoefficients();
    for (int k = 0; k < nFields; ++k) {
        double *fieldCoeffs = coeffs + polynomial->computeFieldCoefficientsOffset(0, k);

        const double *weight = weights;
        for (int i = 0; i < nCoeffs; ++i) {
            fieldCoeffs[i] = *weight * values[0][k];
            ++weight;

            for (int j = 1; j < nEquations; ++j) {
                fieldCoeffs[i] += *weight * values[j][k];
                ++weight;
            }
        }
    }
}

/*!
 * Computes the weights for evaluating the value at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] valueWeights on output will contain the weights of the
 * reconstruction
 */
void ReconstructionKernel::computeValueWeights(const std::array<double, 3> &origin,
                                               const std::array<double, 3> &point, double *valueWeights) const
{
    computeValueLimitedWeights(getDegree(), origin, point, nullptr, valueWeights);
}

/*!
 * Computes the weights for evaluating the value at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] valueWeights on output will contain the weights of the
 * reconstruction
 */
void ReconstructionKernel::computeValueWeights(uint8_t degree, const std::array<double, 3> &origin,
                                               const std::array<double, 3> &point, double *valueWeights) const
{
    computeValueLimitedWeights(degree, origin, point, nullptr, valueWeights);
}

/*!
 * Computes the weights for evaluating the limited value at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param limiters are the scale factors of the limiter
 * \param[out] valueWeights on output will contain the weights of the
 * reconstruction
 */
void ReconstructionKernel::computeValueLimitedWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point,
                                                      const double *limiters, double *valueWeights) const
{
    computeValueLimitedWeights(getDegree(), origin, point, limiters, valueWeights);
}

/*!
 * Computes the weights for evaluating the limited value at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param limiters are the scale factors of the limiter
 * \param[out] valueWeights on output will contain the weights of the
 * reconstruction
 */
void ReconstructionKernel::computeValueLimitedWeights(uint8_t degree, const std::array<double, 3> &origin, const std::array<double, 3> &point,
                                                      const double *limiters, double *valueWeights) const
{
    assert(degree <= getDegree());

    uint8_t dimensions = getDimensions();

    int nEquations = getEquationCount();
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, dimensions);

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);

    ReconstructionPolynomial::evalPointBasisValues(degree, dimensions, origin, point, csi);
    if (limiters) {
        applyLimiter(degree, limiters, csi);
    }

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                nEquations, nCoeffs, 1., getPolynomialWeights(), nEquations,
                csi, 1, 0, valueWeights, 1);
}

/*!
 * Computes the weights for evaluating the directional derivative at a given
 * point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param direction is the direction of the derivative
 * \param[out] derivativeWeights on output will contain the weights of the
 * derivative
 */
void ReconstructionKernel::computeDerivativeWeights(const std::array<double, 3> &origin,
                                                    const std::array<double, 3> &point, const std::array<double, 3> &direction,
                                                    double *derivativeWeights) const
{
    computeDerivativeLimitedWeights(getDegree(), origin, point, direction, nullptr, derivativeWeights);
}

/*!
 * Computes the weights for evaluating the directional derivative at a given
 * point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param direction is the direction of the derivative
 * \param[out] derivativeWeights on output will contain the weights of the
 * derivative
 */
void ReconstructionKernel::computeDerivativeWeights(uint8_t degree, const std::array<double, 3> &origin,
                                                    const std::array<double, 3> &point, const std::array<double, 3> &direction,
                                                    double *derivativeWeights) const
{
    computeDerivativeLimitedWeights(degree, origin, point, direction, nullptr, derivativeWeights);

}

/*!
 * Computes the weights for evaluating the limited directional derivative
 * at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param direction is the direction of the derivative
 * \param limiters are the scale factors of the limiters
 * \param[out] derivativeWeights on output will contain the weights of the
 * derivative
 */
void ReconstructionKernel::computeDerivativeLimitedWeights(const std::array<double, 3> &origin,
                                                           const std::array<double, 3> &point, const std::array<double, 3> &direction,
                                                           const double *limiters,
                                                           double *derivativeWeights) const
{
    computeDerivativeLimitedWeights(getDegree(), origin, point, direction, limiters, derivativeWeights);
}

/*!
 * Computes the weights for evaluating the limited directional derivative
 * at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param direction is the direction of the derivative
 * \param limiters are the scale factors of the limiters
 * \param[out] derivativeWeights on output will contain the weights of the
 * derivative
 */
void ReconstructionKernel::computeDerivativeLimitedWeights(uint8_t degree, const std::array<double, 3> &origin,
                                                           const std::array<double, 3> &point, const std::array<double, 3> &direction,
                                                           const double *limiters,
                                                           double *derivativeWeights) const
{
    assert(degree <= getDegree());

    uint8_t dimensions = getDimensions();

    int nEquations = getEquationCount();
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, dimensions);

    BITPIT_CREATE_WORKSPACE(dcsi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);

    ReconstructionPolynomial::evalPointBasisDerivatives(degree, dimensions, origin, point, direction, dcsi);
    if (limiters) {
        applyLimiter(degree, limiters, dcsi);
    }

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                nEquations, nCoeffs, 1., getPolynomialWeights(), nEquations,
                dcsi, 1, 0, derivativeWeights, 1);
}

/*!
 * Computes the weights for evaluating the gradient at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] gradientWeights on output will contain the weights of the
 * gradient
 */
void ReconstructionKernel::computeGradientWeights(const std::array<double, 3> &origin,
                                                  const std::array<double, 3> &point,
                                                  std::array<double, 3> *gradientWeights) const
{
    computeGradientLimitedWeights(getDegree(), origin, point, nullptr, gradientWeights);
}

/*!
 * Computes the weights for evaluating the gradient at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] gradientWeights on output will contain the weights of the
 * gradient
 */
void ReconstructionKernel::computeGradientWeights(uint8_t degree, const std::array<double, 3> &origin,
                                                  const std::array<double, 3> &point,
                                                  std::array<double, 3> *gradientWeights) const
{
    computeGradientLimitedWeights(degree, origin, point, nullptr, gradientWeights);
}

/*!
 * Computes the weights for evaluating the limited gradient at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] gradientWeights on output will contain the weights of the
 * gradient
 */
void ReconstructionKernel::computeGradientLimitedWeights(const std::array<double, 3> &origin,
                                                         const std::array<double, 3> &point, const double *limiters,
                                                         std::array<double, 3> *gradientWeights) const
{
    computeGradientLimitedWeights(getDegree(), origin, point, limiters, gradientWeights);
}

/*!
 * Computes the weights for evaluating the limited gradient at a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction will be evaluated
 * \param[out] gradientWeights on output will contain the weights of the
 * gradient
 */
void ReconstructionKernel::computeGradientLimitedWeights(uint8_t degree, const std::array<double, 3> &origin,
                                                         const std::array<double, 3> &point, const double *limiters,
                                                         std::array<double, 3> *gradientWeights) const
{
    assert(degree <= getDegree());

    uint8_t dimensions = getDimensions();

    int nEquations = getEquationCount();
    int nCoeffs = ReconstructionPolynomial::getCoefficientCount(degree, dimensions);

    BITPIT_CREATE_WORKSPACE(dcsi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);

    const double *polynomialWeights = getPolynomialWeights();

    std::array<double, 3> direction = {{0., 0., 0.}};
    for (int d = 0; d < dimensions; ++d) {
        direction[d] = 1.;

        ReconstructionPolynomial::evalPointBasisDerivatives(degree, dimensions, origin, point, direction, dcsi);
        if (limiters) {
            applyLimiter(degree, limiters, dcsi);
        }

        // Weights are stored in contiguous three-dimensional arrays, this
        // means we can access the weights of the current dimensions using
        // the dimension as offset and a stride of three elements
        cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                    nEquations, nCoeffs, 1., polynomialWeights, nEquations,
                    dcsi, 1, 0, gradientWeights->data() + d, 3);

        direction[d] = 0.;
    }

    // Explicitly zero unused components
    if (dimensions != ReconstructionPolynomial::MAX_DIMENSIONS) {
        for (int i = 0; i < nEquations; ++i) {
            for (int d = dimensions; d < ReconstructionPolynomial::MAX_DIMENSIONS; ++d) {
                gradientWeights[i][d] = 0.;
            }
        }
    }
}

/*!
 * Apply the limiter to the given coefficients.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param degree is the degree of the polynomial
 * \param limiters are the scale factors of the limiter
 * \param[in,out] coeffs are the coefficients to be limited
 */
void ReconstructionKernel::applyLimiter(uint8_t degree, const double *limiters, double *coeffs) const
{
    if (degree < 1) {
        return;
    }

    uint8_t dimensions = getDimensions();
    const std::vector<uint16_t> &nDegreeCoeffs = ReconstructionPolynomial::getDegreeCoefficientsCount(dimensions);

    int coeffEnd = nDegreeCoeffs[0];
    for (int n = 1; n <= degree; ++n) {
        double limiterScaleFactor = limiters[n - 1];

        int coeffBegin = coeffEnd;
        coeffEnd = coeffBegin + nDegreeCoeffs[n];
        for (int k = coeffBegin; k < coeffEnd; ++k) {
            coeffs[k] *= limiterScaleFactor;
        }
    }
}

/*!
 * Displays reconstruction polynomial information to an output stream.
 *
 * \param out is the output stream
 * \param tolerance is the tolerance below which the weight will not be
 * displayed
 */
void ReconstructionKernel::display(std::ostream &out, double tolerance) const
{
    int nCoeffs = getCoefficientCount();
    for (int i = 0; i < nCoeffs; ++i) {
        out << i << " ";
        for (int j = 0; j < m_nEquations; ++j) {
            double weigth = m_weights[linearalgebra::linearIndexColMajor(j, i, m_nEquations, m_nCoeffs)];
            if (std::abs(weigth) < tolerance) {
                continue;
            }

            out << "(" << j << "," << weigth << ") ";
        }

        out << std::endl;
    }
}

/*!
 * \class ReconstructionAssembler
 * \ingroup discretization
 *
 * \brief The ReconstructionAssembler class allows to define reconstruction
 * polinoymial.
 */

/*!
 * Is the threshold for which a singuler value is considered zero.
 */
double ReconstructionAssembler::SVD_ZERO_THRESHOLD = 1e-9;


/*!
 * Constructor.
 */
ReconstructionAssembler::ReconstructionAssembler()
{
    initialize(0, 0, false);
}

/*!
 * Constructor.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 */
ReconstructionAssembler::ReconstructionAssembler(uint8_t degree, uint8_t dimensions)
{
    initialize(degree, dimensions, false);
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionAssembler::swap(ReconstructionAssembler &other) noexcept
{
    std::swap(other.m_degree, m_degree);
    std::swap(other.m_dimensions, m_dimensions);
    std::swap(other.m_nCoeffs, m_nCoeffs);
    std::swap(other.m_constraintsOrder, m_constraintsOrder);
    std::swap(other.m_leastSquaresOrder, m_leastSquaresOrder);
    std::swap(other.m_leastSquaresScaleFactors, m_leastSquaresScaleFactors);
    std::swap(other.m_A, m_A);
    std::swap(other.m_C, m_C);
    std::swap(other.m_sigma, m_sigma);
    std::swap(other.m_U, m_U);
    std::swap(other.m_S, m_S);
    std::swap(other.m_Vt, m_Vt);
    std::swap(other.m_SVDWorkspace, m_SVDWorkspace);
    std::swap(other.m_w, m_w);
}

/*!
 * Initialize the assembler.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param release if true, possible unneeded memory hold by the assembler
 * will be released, otherwise the assembler will be initialized but possible
 * unneeded memory will not be released
 */
void ReconstructionAssembler::initialize(uint8_t degree, uint8_t dimensions, bool release)
{
    assert(degree <= ReconstructionPolynomial::MAX_DEGREE);
    m_degree = degree;

    assert(dimensions <= ReconstructionPolynomial::MAX_DIMENSIONS);
    m_dimensions = dimensions;

    m_nCoeffs = ReconstructionPolynomial::getCoefficientCount(m_degree, m_dimensions);

    clear(release);
}

/*!
 * Clear the reconstruction.
 *
 * \param release if true, the memory hold by the assembler will be released,
 * otherwise the assembler will be cleared but its memory will not be released
 */
void ReconstructionAssembler::clear(bool release)
{
    m_constraintsOrder.clear();
    m_leastSquaresOrder.clear();
    m_leastSquaresScaleFactors.clear();

    m_A.clear();
    m_C.clear();

    m_sigma.clear();
    m_U.clear();
    m_S.clear();
    m_Vt.clear();
    m_SVDWorkspace.resize(1);
    m_w.clear();

    if (release) {
        m_constraintsOrder.shrink_to_fit();
        m_leastSquaresOrder.shrink_to_fit();
        m_leastSquaresScaleFactors.shrink_to_fit();

        m_A.shrink_to_fit();
        m_C.shrink_to_fit();

        m_sigma.shrink_to_fit();
        m_U.shrink_to_fit();
        m_S.shrink_to_fit();
        m_Vt.shrink_to_fit();
        m_SVDWorkspace.shrink_to_fit();
        m_w.shrink_to_fit();
    }
}

/*!
 * Get the degree of the polynomial.
 *
 * \return The degree of the polynomial.
 */
uint8_t ReconstructionAssembler::getDegree() const
{
    return m_degree;
}

/*!
 * Get the number of space dimensions.
 *
 * \return The number of space dimensions.
 */
uint8_t ReconstructionAssembler::getDimensions() const
{
    return m_dimensions;
}

/*!
 * Get the number of coefficients of the reconstruction polynomial.
 *
 * \return The number of coefficients of the reconstruction polynomial.
 */
uint16_t ReconstructionAssembler::getCoefficientCount() const
{
    return m_nCoeffs;
}

/*!
 * Count the number of constrain-type equations added to the assembler.
 *
 * \return The number of constraint-type equations added to the assembler.
 */
int ReconstructionAssembler::countConstraints() const
{
    return m_constraintsOrder.size();
}

/*!
 * Count the number of least square-type equations added to the assembler.
 *
 * \return The number of least square-type equations added to the assembler.
 */
int ReconstructionAssembler::countLeastSquares() const
{
    return m_leastSquaresOrder.size();
}

/*!
 * Count the number of equations added to the assembler.
 *
 * \return The number of equations added to the assembler.
 */
int ReconstructionAssembler::countEquations() const
{
    int nConstraints  = countConstraints();
    int nLeastSquares = countLeastSquares();
    int nEquations    = nConstraints + nLeastSquares;

    return nEquations;
}

/*!
 * Add a point value equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the basis will be evaluated
 * \param scaleFactor the scale factor associated to the equation
 */
void ReconstructionAssembler::addPointValueEquation(ReconstructionType type,
                                                    const std::array<double, 3> &origin,
                                                    const std::array<double, 3> &point,
                                                    double scaleFactor)
{
    double *equationCoeffs = _addEquation(type, scaleFactor);
    ReconstructionPolynomial::evalPointBasisValues(getDegree(), getDimensions(), origin, point, equationCoeffs);
}

/*!
 * Add a point derivative equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the basis will be evaluated
 * \param direction is the direction of the derivative
 * \param scaleFactor the scale factor associated to the equation
 */
void ReconstructionAssembler::addPointDerivativeEquation(ReconstructionType type,
                                                         const std::array<double, 3> &origin,
                                                         const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction,
                                                         double scaleFactor)
{
    double *equationCoeffs = _addEquation(type, scaleFactor);
    ReconstructionPolynomial::evalPointBasisDerivatives(getDegree(), getDimensions(), origin, point, direction, equationCoeffs);
}

/*!
 * Add a cell average equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param cell is the cell
 * \param origin is the point chosen as origin of the reconstruction
 * \param vertexCoords are the vertecx coordinates
 * \param scaleFactor the scale factor associated to the equation
 */
void ReconstructionAssembler::addCellAverageEquation(ReconstructionType type,
                                                     const Cell &cell,
                                                     const std::array<double, 3> &origin,
                                                     const std::array<double, 3> *vertexCoords,
                                                     double scaleFactor)
{
    double *equationCoeffs = _addEquation(type, scaleFactor);
    ReconstructionPolynomial::evalCellBasisValues(getDegree(), getDimensions(), origin, cell, vertexCoords, equationCoeffs);
}

/*!
 * Internal function to add an equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param scaleFactor the scale factor associated to the equation
 */
double * ReconstructionAssembler::_addEquation(ReconstructionType type, double scaleFactor)
{
    // Update equation information
    int nEquations = countEquations();
    switch (type) {

    case TYPE_CONSTRAINT:
        m_constraintsOrder.emplace_back(nEquations);
        break;

    case TYPE_LEAST_SQUARE:
        m_leastSquaresOrder.emplace_back(nEquations);
        m_leastSquaresScaleFactors.emplace_back(scaleFactor);
        break;

    }

    // Prepare storage for equation coefficients
    int nCoeffs = getCoefficientCount();

    double *equationCoeffsStorage = nullptr;
    switch (type) {

    case TYPE_CONSTRAINT:
        m_C.resize(m_C.size() + nCoeffs);
        equationCoeffsStorage = m_C.data() + m_C.size() - nCoeffs;
        break;

    case TYPE_LEAST_SQUARE:
        m_A.resize(m_A.size() + nCoeffs);
        equationCoeffsStorage = m_A.data() + m_A.size() - nCoeffs;
        break;

    }

    return equationCoeffsStorage;
}

/*!
 * Assembles the reconstruction kernel.
 *
 * Before computing kernel weights, the kernel will be properly initialized
 * and possible unneeded memory hold by the kernel will be released.
 *
 * \param[out] kernel on output will contain the reconstrucion kernel
 */
void ReconstructionAssembler::assembleKernel(ReconstructionKernel *kernel) const
{
    // Initialize reconstruciton kernel
    uint8_t degree     = getDegree();
    uint8_t dimensions = getDimensions();

    int nEquations = countEquations();

    kernel->initialize(degree, dimensions, nEquations, true);

    // Update the kernel
    updateKernel(kernel);
}

/*!
 * Updates the reconstruction kernel.
 *
 * Before computing kernel weights, the kernel will not be initialized.
 *
 * \param[out] kernel on output will contain the reconstrucion kernel
 */
void ReconstructionAssembler::updateKernel(ReconstructionKernel *kernel) const
{
    // Get the number of equations
    int nConstraints  = countConstraints();
    int nLeastSquares = countLeastSquares();
    int nEquations    = countEquations();

    // Get the number of polynomial coefficients
    int nCoeffs = getCoefficientCount();

    // Evaluate normalized least square scale factors
    if (nLeastSquares > 0) {
        double maxLeastSquareScaleFactor = std::abs(m_leastSquaresScaleFactors[0]);
        for (int k = 1; k < nLeastSquares; ++k) {
            maxLeastSquareScaleFactor = std::max(std::abs(m_leastSquaresScaleFactors[k]), maxLeastSquareScaleFactor);
        }

        m_w.resize(nLeastSquares);
        for (int k = 0; k < nLeastSquares; ++k) {
            m_w[k] = m_leastSquaresScaleFactors[k] / maxLeastSquareScaleFactor;
        }
    }

    // The linear-constrained are introduced in the least-squares problem
    // through Lagrange multipliers. The resulting linear system is:
    //
    // | A^t A w  C^t | |x     | = |A^t w b|
    // | C        0   | |lambda|   |d      |
    //
    // with A and C the least-squares and costraints equations respectively,
    // and b and d their corresponding RHSs. x are the coefficients of
    // the polynomial and lambda the lagrange multipliers. w are the normalized
    // east square scale factors.
    //
    // This system is denoted by S:
    //
    //     |  x   |   |A^t w  0| |b|
    // |S| |      | = |        | | |
    //     |lambda|   |0      I| |d|
    //
    // The matrices S and S^-1 are symmetric and only the upper portions are
    // computed
    int nUnknowns = nCoeffs + nConstraints;

    m_S.resize(nUnknowns * nUnknowns);
    for (int i = 0; i < nCoeffs; ++i) {
        for (int j = i; j < nCoeffs; ++j) {
            // Compute A^t A on the fly
            double ATA_ij = 0.;
            for (int k = 0; k < nLeastSquares; ++k) {
                int A_ki_idx = linearalgebra::linearIndexRowMajor(k, i, nLeastSquares, nCoeffs);
                int A_kj_idx = linearalgebra::linearIndexRowMajor(k, j, nLeastSquares, nCoeffs);
                ATA_ij += m_A[A_ki_idx] * m_A[A_kj_idx] * m_w[k];
            }

            int l = linearalgebra::linearIndexColMajor(i, j, nUnknowns, nUnknowns);
            m_S[l] = ATA_ij;

            int m = linearalgebra::linearIndexColMajor(j, i, nUnknowns, nUnknowns);
            m_S[m] = ATA_ij;
        }

        for (int j = nCoeffs; j < nUnknowns; ++j) {
            int l     = linearalgebra::linearIndexColMajor(i, j, nUnknowns, nUnknowns);
            int C_idx = linearalgebra::linearIndexRowMajor(j - nCoeffs, i, nConstraints, nCoeffs);
            m_S[l] = m_C[C_idx];

            int m = linearalgebra::linearIndexColMajor(j, i, nUnknowns, nUnknowns);
            m_S[m] = m_S[l];
        }
    }

    // Compute inverse S matrix
    // Since S may me be rank-deficit (eg if not enough neighbours are available)
    // the pseudo-inverse is used. This corresponds of computing the least-norm
    // solution of the problem.
    computePseudoInverse(nUnknowns, nUnknowns, SVD_ZERO_THRESHOLD, m_S.data());

    // Weights needed to evaluate the polynomial coefficients come from the
    // following equation:
    //
    // |  x   |        |A^t w  0| |b|          |b|
    // |      | = S^-1 |        | | | = S^-1 Q | |
    // |lambda|        |0      I| |d|          |d|
    //
    // Since we are interested only in x (the polynomial coefficients) only
    // the first nCoeffs rows of the matrix S^-1 Q are computed. Those values
    // are the polynomial weights.
    //
    // Weigths are stored according the order in which the equations have been
    // added.
    double *weights = kernel->getPolynomialWeights();
    for (int i = 0; i < nCoeffs; ++i) {
        for (int j = 0; j < nEquations; ++j) {
            double value = 0;
            for (int k = 0; k < nUnknowns; ++k) {
                int l = linearalgebra::linearIndexColMajorSymmetric(i, k, nUnknowns, nUnknowns, 'U');
                if (k < nCoeffs && j < nLeastSquares) {
                    int A_jk_idx = linearalgebra::linearIndexRowMajor(j, k, nLeastSquares, nCoeffs);
                    value += m_S[l] * m_A[A_jk_idx] * m_w[j];
                } else if ((k - nCoeffs) == (j - nLeastSquares)) {
                    value += m_S[l];
                }
            }

            int equation;
            if (j < nLeastSquares) {
                equation = m_leastSquaresOrder[j];
            } else {
                equation = m_constraintsOrder[j - nLeastSquares];
            }

            int weightLineraIndex = linearalgebra::linearIndexColMajor(equation, i, nEquations, nCoeffs);
            weights[weightLineraIndex] = value;
        }
    }
}
/*!
 * Computes the pseudo inverse of a matrix using a singular value decomposition
 *
 * See "Solving Ill-Conditioned And Singular Linear Systems: A Tutorial On
 * Regularization", by Arnold Neumaier (see https://www.mat.univie.ac.at/~neum/ms/regtutorial.pdf).
 *
 * \param m number of columns
 * \param n number of rows
 * \param zeroThreshold is the threshold below which a singuler value is
 * considered zero
 * \param[in,out] A on input matrix in coumn-major ordering, on output its
 * pseudo-inverse
 */
void ReconstructionAssembler::computePseudoInverse(int m, int n, double zeroThreshold, double *A) const
{
    // Compute SVD
    //
    // A = U * Sigma * Vt (Equation 21)
    int k = std::min(m,n);

    m_sigma.resize(k);
    m_U.resize(m * k);
    m_Vt.resize(k * n);

    char jobU  = 'S';
    char jobVT = 'S';

    int workspaceSize = -1;

    int info;

    LAPACK_dgesvd(&jobU, &jobVT, &m, &n, A, &m, m_sigma.data(), m_U.data(), &m, m_Vt.data(), &k,
                  m_SVDWorkspace.data(), &workspaceSize, &info);

    workspaceSize = m_SVDWorkspace[0];
    if (workspaceSize > (int) m_SVDWorkspace.size()) {
        m_SVDWorkspace.resize(workspaceSize);
    }

    LAPACK_dgesvd(&jobU, &jobVT, &m, &n, A, &m, m_sigma.data(), m_U.data(), &m, m_Vt.data(), &k,
                  m_SVDWorkspace.data(), &workspaceSize, &info);

    if (info > 0) {
        log::cout() << "SVD failed in ReconstructionAssembler::computePseudoInverse()" <<std::endl;
        exit(1);
    }

    // Inv(A) = V * Sigma^ + *U^T (Equation 22)
    //
    // u = sigma^ + *U
    // and is stored in U
    for (int i = 0; i < k; ++i) {
       double sigma_plus = (m_sigma[i] > zeroThreshold) ? (1. / m_sigma[i]) : m_sigma[i];
       cblas_dscal(m, sigma_plus, &m_U[i*m], 1);
    }

    // Inv(A) = (Vt)^T * u^T
    cblas_dgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans, CBLAS_TRANSPOSE::CblasTrans,
                n, m, k, 1., m_Vt.data(), k, m_U.data(), m, 0., A, n);
}

/*!
 * \class Reconstruction
 * \ingroup discretization
 *
 * \brief The Reconstruction class allows to build and apply a polynomial
 * reconstructions.
 */

/*!
 * Default constructor.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 */
Reconstruction::Reconstruction(uint8_t degree, uint8_t dimensions)
    : ReconstructionAssembler(degree, dimensions)
{
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void Reconstruction::swap(Reconstruction &other) noexcept
{
    ReconstructionKernel::swap(other);
    ReconstructionAssembler::swap(other);
}

/*!
 * Initialize the reconstruction.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param release if true, possible unneeded memory hold by the reconstruction
 * will be released, otherwise the reconstruction will be initialized but
 * possible unneeded memory will not be released
 */
void Reconstruction::initialize(uint8_t degree, uint8_t dimensions, bool release)
{
    ReconstructionAssembler::initialize(degree, dimensions, release);
}

/*!
 * Clear the reconstruction.
 *
 * \param release if true, the memory hold by the reconstruction will be
 * released, otherwise the reconstruction will be cleared but its memory
 * will not be released
 */
void Reconstruction::clear(bool release)
{
    ReconstructionAssembler::clear(release);

    if (release) {
        ReconstructionKernel().swap(*this);
    }
}

/*!
 * Computes the weights to be used in order to calculate the coefficients
 * of the polynomial. The coefficients are such that the conditions decoded
 * in equations are enforced.
 */
void Reconstruction::assemble()
{
    if (ReconstructionKernel::getEquationCount() != ReconstructionAssembler::countEquations()) {
        ReconstructionAssembler::assembleKernel(this);
    } else {
        ReconstructionAssembler::updateKernel(this);
    }
}

}
