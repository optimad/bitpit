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

template class DiscretizationStencilSolver<StencilScalar>;

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

}
