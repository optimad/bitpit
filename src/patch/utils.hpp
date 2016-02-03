/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#ifndef __BITP_MESH_UTILS_HPP__
#define __BITP_MESH_UTILS_HPP__

/*! \file */

#include <array>
#include <functional>
#include <vector>

namespace utils {

std::array<double, 3> cross_3D(std::array<double, 3> &x, std::array<double, 3> &y);
void normalize_3D(std::array<double, 3> &x);
void transpose_3D(std::array<std::array<double, 3>, 3> &A);

}

#endif
