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

#define __BITPIT_STENCIL_SOLVER_SRC__

#include "stencil_solver.hpp"

namespace bitpit {

// Explicit instantization
template class DiscretizationStencilSolverAssembler<StencilScalar>;
template class DiscretizationStencilSolverAssembler<StencilVector>;
template class DiscretizationStencilSolverAssembler<StencilBlock>;

template class DiscretizationStencilSolver<StencilScalar>;
template class DiscretizationStencilSolver<StencilVector>;
template class DiscretizationStencilSolver<StencilBlock>;

// Template specializations

/*!
 * Set block size.
 */
template<>
void DiscretizationStencilSolverAssembler<StencilScalar>::setBlockSize()
{
    setBlockSize(1);
}

/*!
 * Get the raw value of the specified element.
 *
 * \param element is the stencil element
 * \param item is the requested block item
 * \result The values of the specified weight.
 */
template<>
double DiscretizationStencilSolverAssembler<StencilScalar>::getRawValue(const StencilScalar::weight_type &element, int item) const
{
    BITPIT_UNUSED(item);

    return element;
}

/*!
 * Set block size.
 */
template<>
void DiscretizationStencilSolverAssembler<StencilVector>::setBlockSize()
{
    setBlockSize(sizeof(typename StencilVector::weight_type) / sizeof(typename StencilVector::weight_type::value_type));
}

/*!
 * Get the raw value of the specified element.
 *
 * \param element is the stencil element
 * \param item is the requested block item
 * \result The values of the specified weight.
 */
template<>
double DiscretizationStencilSolverAssembler<StencilVector>::getRawValue(const StencilVector::weight_type &element, int item) const
{
    return element[item];
}

/*!
 * Initialize block size.
 *
 * The block size is set equal to the square root of the weight size; if the
 * square root of the weight type is not an integer number, an exception is
 * throw.
 *
 * Block size is evaluated from the first weight. If the size of the other
 * weights don't match the size of the first one, an exception is thrown
 * only when bitpit is compiled in debug mode, otherwise the error is
 * silently ignored.
 */
template<>
void DiscretizationStencilSolverAssembler<StencilBlock>::setBlockSize()
{
    // Get block size
    //
    // Block size is evaluated from the first weight.
    //
    // The block size is set equal to the square root of the weight size;
    // if the square root of the weight type is not an integer number, an
    // exception is throw.
    m_blockSize = - 1;
    for (long i = 0; i < getRowCount(); ++i) {
        const StencilBlock &stencil = getRowStencil(i);
        std::size_t stencilSize = stencil.size();
        if (stencilSize == 0) {
            break;
        }

        const StencilBlock::weight_type *weightData = stencil.weightData();
        double blockSize = std::sqrt(weightData[0].size());
        if (blockSize != std::sqrt(weightData[0].size())) {
            throw std::runtime_error("Weights size should be a square.");
        }
        setBlockSize(static_cast<int>(blockSize));
        break;
    }

    if (m_blockSize == -1) {
        throw std::runtime_error("All weights should have a size greater than zero.");
    }

#ifdef DEBUG
    // Validate block size
    //
    // All weight sizes should match
    for (long i = 0; i < getRowCount(); ++i) {
        const StencilBlock &stencil = getRowStencil(i);
        const StencilBlock::weight_type *weightData = stencil.weightData();
        std::size_t stencilSize = stencil.size();

        for (std::size_t k = 0; k < stencilSize; ++k) {
            if (weightData[k].size() != m_blockSize)) {
                throw std::runtime_error("All stencils weights should have the same size.");
            }
        }

        if (stencil.getConstant().size() != m_blockSize) {
            throw std::runtime_error("The stencils constant should have the same size of the stencil weights.");
        }
    }
#endif
}

/*!
 * Get the raw value of the specified element.
 *
 * \param element is the stencil element
 * \param item is the requested block item
 * \result The values of the specified weight.
 */
template<>
double DiscretizationStencilSolverAssembler<StencilBlock>::getRawValue(const StencilBlock::weight_type &element, int item) const
{
    return element[item];
}

}
