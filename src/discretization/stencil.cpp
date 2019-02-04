/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

template class DiscreteStencil<double>;
template class DiscreteStencil<std::array<double, 3>>;

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
    int nBuckets = stencil.getBucketCount();
    std::vector<int> sizes(nBuckets);
    for (int k = 0; k < nBuckets; ++k) {
        sizes[k] = stencil.size(k);
    }

    bitpit::StencilVector stencil_B(sizes);
    for (std::size_t k = 0; k < stencil_B.size(); ++k) {
        stencil_B.rawSetPattern(k, stencil.rawGetPattern(k));
        stencil_B.rawSetWeight(k, ::operator*(stencil.rawGetWeight(k), vector));
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
    int nBuckets = stencil.getBucketCount();
    std::vector<int> sizes(nBuckets);
    for (int k = 0; k < nBuckets; ++k) {
        sizes[k] = stencil.size(k);
    }

    stencil_dotProduct->initialize(sizes);

    for (std::size_t k = 0; k < stencil_dotProduct->size(); ++k) {
        stencil_dotProduct->rawSetPattern(k, stencil.rawGetPattern(k));
        stencil_dotProduct->rawSetWeight(k, ::dotProduct(stencil.rawGetWeight(k), vector));
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

    for (std::size_t k = 0; k < stencil_projection->size(); ++k) {
        bitpit::StencilVector::weight_type &weight = stencil_projection->rawGetWeight(k);
        weight = ::dotProduct(weight, direction) * direction;
    }

    bitpit::StencilVector::weight_type &constant = stencil_projection->getConstant();
    constant = ::dotProduct(constant, direction) * direction;
}
