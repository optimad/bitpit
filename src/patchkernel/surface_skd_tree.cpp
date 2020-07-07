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

#include "surface_skd_tree.hpp"

namespace bitpit {

/*!
* \class SurfaceSkdTree
*
* \brief The SurfaceSkdTree implements a Bounding Volume Hierarchy tree for
* surface patches.
*/

/*!
* Constructor.
*
* \param patch is the surface patch that will be use to build the tree
* \param includeGhosts if set to true (default value) the ghost cells are included in the tree
*/
SurfaceSkdTree::SurfaceSkdTree(const SurfaceKernel *patch, bool includeGhosts)
    : PatchSkdTree(patch,includeGhosts)
{
}

}
