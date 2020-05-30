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
 * Initialize block size.
 *
 * Block size can only be initialized once.
 */
template<>
void DiscretizationStencilSolverAssembler<StencilScalar>::initializeBlockSize()
{
    assert(m_blockSize == -1);

    m_blockSize = 1;
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
 * Initialize block size.
 *
 * Block size can only be initialized once.
 */
template<>
void DiscretizationStencilSolverAssembler<StencilVector>::initializeBlockSize()
{
    assert(m_blockSize == -1);

    m_blockSize = sizeof(typename StencilVector::weight_type) / sizeof(typename StencilVector::weight_type::value_type);
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
 * Vector can only be initialized once.
 *
 * \param stencils are the stencils
 */
template<>
void DiscretizationStencilSolverAssembler<StencilBlock>::initializeBlockSize()
{
    assert(m_blockSize == -1);

    std::size_t nStencils = m_stencils->size();

    // Set the block size equal to the size of the first weight
    m_blockSize = - 1;
    for (std::size_t i = 0; i < nStencils; ++i) {
        const StencilBlock &stencil = (*m_stencils)[i];
        const StencilBlock::weight_type *weightData = stencil.weightData();
        std::size_t stencilSize = stencil.size();

        for (std::size_t k = 0; k < stencilSize; ++k) {
            m_blockSize = weightData[k].size();
            break;
        }

        if (m_blockSize >= 0) {
            break;
        }
    }

    assert(m_blockSize == -1);

#ifdef DEBUG
    // All weight sizes should match
    for (std::size_t i = 0; i < nStencils; ++i) {
        const StencilBlock &stencil = (*m_stencils)[i];
        const StencilBlock::weight_type *weightData = stencil.weightData();
        std::size_t stencilSize = stencil.size();

        for (std::size_t k = 0; k < stencilSize; ++k) {
            assert((int) weightData[k].size() == m_blockSize);
        }

        assert((int) stencil.getConstant().size() == m_blockSize);
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
