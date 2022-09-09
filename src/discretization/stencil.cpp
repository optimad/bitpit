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

}
