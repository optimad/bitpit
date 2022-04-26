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

#define __BITPIT_STENCIL_SRC__

#include <cassert>

#include "bitpit_operators.hpp"

#include "stencil.hpp"

// Explicit instantization
namespace bitpit {

template class DiscreteStencilWeightPool<double>;
template class DiscreteStencilWeightPool<std::array<double, 3>>;
template class DiscreteStencilWeightPool<std::vector<double>>;

template class DiscreteStencil<double>;
template class DiscreteStencil<std::array<double, 3>>;
template class DiscreteStencil<std::vector<double>>;

template class MPDiscreteStencil<double>;
template class MPDiscreteStencil<std::array<double, 3>>;
template class MPDiscreteStencil<std::vector<double>>;

/*!
 * Sum the specified value to the target.
 *
 * \param value is the value that will be summed
 * \param factor is the factor the value will be multiplied with
 * \param target on output will contain the original weight plus the value multiplied by the factor
 */
template<>
void DiscreteStencil<std::array<double, 3>>::rawSumValue(const std::array<double, 3> &value, double factor, std::array<double, 3> *target)
{
    for (int i = 0; i < 3; ++i) {
        (*target)[i] += factor * value[i];
    }
}

/*!
 * Sum the specified value to the target.
 *
 * The target will be resized to match the size of the value to be summed. If the value size is
 * greater that the target size, missing target elements will be initialized to zero before
 * summing the specified value.
 *
 * \param value is the value that will be summed
 * \param factor is the factor the value will be multiplied with
 * \param target on output will contain the original weight plus the value multiplied by the factor
 */
template<>
void DiscreteStencil<std::vector<double>>::rawSumValue(const std::vector<double> &value, double factor, std::vector<double> *target)
{
    std::size_t valueSize  = value.size();
    std::size_t targetSize = target->size();
    std::size_t commonSize = std::min(valueSize, targetSize);

    for (std::size_t i = 0; i < commonSize; ++i) {
        (*target)[i] += factor * value[i];
    }

    if (valueSize > targetSize) {
        target->insert(target->end(), value.cbegin() + commonSize, value.cend());

        if (factor != 1.) {
            auto targetBegin = target->begin();
            auto targetEnd   = target->end();
            for (auto itr = targetBegin + commonSize; itr != targetEnd; ++itr) {
                *itr *= factor;
            }
        }
    }
}

/**
 * Set the source value into the target.
 *
 * \param source is the value that will be set
 * \param[out] target on output will contain the source value
 */
template<>
void DiscreteStencil<std::array<double, 3>>::rawCopyValue(const std::array<double, 3> &source, std::array<double, 3> *target)
{
    std::copy_n(source.data(), source.size(), target->data());
}

/**
 * Copy the source value into the target.
 *
 * The target will be resized to match the size of the source.
 *
 * \param source is the value that will be set
 * \param[out] target on output will contain the source value
 */
template<>
void DiscreteStencil<std::vector<double>>::rawCopyValue(const std::vector<double> &source, std::vector<double> *target)
{
    std::size_t sourceSize = source.size();
    std::size_t targetSize = target->size();
    std::size_t commonSize = std::min(sourceSize, targetSize);

    std::copy_n(source.data(), commonSize, target->data());

    if (sourceSize < targetSize) {
        target->resize(sourceSize);
    } else if (sourceSize > targetSize) {
        target->insert(target->end(), source.begin() + commonSize, source.end());
    }
}

/*!
* Optimize the specified weight.
*
* \param pos is the position of the weight to check
* \param tolerance is the tolerance that will be used for the check
* \result Returns true if the whole weight is neglibile
*/
template<>
bool DiscreteStencil<std::array<double, 3>>::optimizeWeight(std::size_t pos, double tolerance)
{
    std::array<double, 3> &weight = m_weights[pos];
    int nItems = weight.size();

    int nNegligibleItems = 0;
    for (int k = 0; k < nItems; ++k) {
        if (std::abs(weight[k] - m_zero[k]) > tolerance) {
            continue;
        }

        weight[k] = m_zero[k];
        ++nNegligibleItems;
    }

    return (nNegligibleItems == nItems);
}

/*!
* Optimize the specified weight.
*
* \param pos is the position of the weight to check
* \param tolerance is the tolerance that will be used for the check
* \result Returns true if the whole weight is neglibile
*/
template<>
bool DiscreteStencil<std::vector<double>>::optimizeWeight(std::size_t pos, double tolerance)
{
    std::vector<double> &weight = m_weights[pos];
    int nItems = weight.size();

    int nNegligibleItems = 0;
    for (int k = 0; k < nItems; ++k) {
        if (std::abs(weight[k] - m_zero[k]) > tolerance) {
            continue;
        }

        weight[k] = m_zero[k];
        ++nNegligibleItems;
    }

    return (nNegligibleItems == nItems);
}

}

/*!
* The multiplication operator between a scalar stencil and a vector.
*
* \param stencil is the stencil
* \param vector is the vector
* \result The result fo the multiplication.
*/
bitpit::StencilVector operator*(const bitpit::StencilScalar &stencil, const std::array<double, 3> &vector)
{
    return (vector * stencil);
}

/*!
* The multiplication operator between a vector and a scalar stencil.
*
* \param vector is the vector
* \param stencil is the stencil
* \result The result fo the multiplication.
*/
bitpit::StencilVector operator*(const std::array<double, 3> &vector, const bitpit::StencilScalar &stencil)
{
    const std::size_t nItems = stencil.size();
    bitpit::StencilVector stencil_B;
    stencil_B.resize(nItems);
    for (std::size_t n = 0; n < nItems; ++n) {
        stencil_B.setPattern(n, stencil.getPattern(n));
        stencil_B.setWeight(n, ::operator*(stencil.getWeight(n), vector));
    }
    stencil_B.setConstant(::operator*(stencil.getConstant(), vector));

    return stencil_B;
}

/*!
* The dot procduct operator betwee a vector stencil and a vector.
*
* \param stencil is the vector stencil
* \param vector is the vector
* \result The result fo the dot product.
*/
bitpit::StencilScalar dotProduct(const bitpit::StencilVector &stencil, const bitpit::StencilVector::weight_type &vector)
{
    bitpit::StencilScalar stencil_dotProduct;
    dotProduct(stencil, vector, &stencil_dotProduct);

    return stencil_dotProduct;
}

/*!
* The dot procduct operator betwee a vector stencil and a vector.
*
* \param stencil is the vector stencil
* \param vector is the vector
* \param[out] stencil_dotProduct on output will contain the dot product
*/
void dotProduct(const bitpit::StencilVector &stencil, const bitpit::StencilVector::weight_type &vector, bitpit::StencilScalar *stencil_dotProduct)
{
    const std::size_t nItems = stencil.size();
    stencil_dotProduct->resize(nItems);

    const long *patternData = stencil.patternData();
    const bitpit::StencilVector::weight_type *weightData = stencil.weightData();
    long *patternData_dotProduct = stencil_dotProduct->patternData();
    bitpit::StencilScalar::weight_type *weightData_dotProduct = stencil_dotProduct->weightData();
    for (std::size_t n = 0; n < nItems; ++n) {
        patternData_dotProduct[n] = patternData[n];
        weightData_dotProduct[n]  = ::dotProduct(weightData[n], vector);
    }

    stencil_dotProduct->setConstant(::dotProduct(stencil.getConstant(), vector));
}

/*!
* Project the stencil along the specified direction.
*
* \param stencil is the vector stencil
* \param direction is the direction
* \result The projection of the stencil along the specified direction.
*/
bitpit::StencilVector project(const bitpit::StencilVector &stencil, const std::array<double, 3> &direction)
{
    bitpit::StencilVector stencil_projection;
    project(stencil, direction, &stencil_projection);

    return stencil_projection;
}

/*!
* Project the stencil along the specified direction.
*
* \param stencil is the vector stencil
* \param direction is the direction
* \param[out] stencil_projection on output will contain the projection
*/
void project(const bitpit::StencilVector &stencil, const std::array<double, 3> &direction, bitpit::StencilVector *stencil_projection)
{
    stencil_projection->initialize(stencil);
    const std::size_t nItems = stencil_projection->size();

    bitpit::StencilVector::weight_type *weightData_projection = stencil_projection->weightData();
    for (std::size_t n = 0; n < nItems; ++n) {
        bitpit::StencilVector::weight_type &weight = weightData_projection[n];
        weight = ::dotProduct(weight, direction) * direction;
    }

    stencil_projection->setConstant(::dotProduct(stencil_projection->getConstant(), direction) * direction);
}
