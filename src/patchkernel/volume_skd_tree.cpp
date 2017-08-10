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

#include "volume_skd_tree.hpp"

namespace bitpit {

/*!
* \class VolumeSkdTree
*
* \brief The VolumeSkdTree implements a Bounding Volume Hierarchy tree for
* volume patches.
*/

/*!
* Constructor.
*
* \param patch is the volume patch that will be use to build the tree
*/
VolumeSkdTree::VolumeSkdTree(const VolumeKernel *patch)
    : PatchSkdTree(patch)
{
}

}
