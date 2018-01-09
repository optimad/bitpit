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

#include "pod_voloctree.hpp"

namespace bitpit {

/**
 * \class PODVolOctree
 * \ingroup POD
 *
 * \brief The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 *
 * The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 */

/**
 * Creates a new PODVolOctree object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the pod object. MPI_COMM_WORLD is the default value.
 */
PODVolOctree::PODVolOctree(MPI_Comm comm):PODKernel(comm)
# else
PODVolOctree::PODVolOctree():PODKernel()
# endif
{
}

/**
 * Destructor of PODVolOctree
 */
PODVolOctree::~PODVolOctree()
{
}

/**
 * Create the mesh for POD.
 *
 * return Pointer to instantiated mesh.
 */
VolumeKernel* PODVolOctree::createMesh()
{
    VolumeKernel *mesh = new VolOctree();

#if BITPIT_ENABLE_MPI
    mesh->setCommunicator(m_communicator);
#endif

    return mesh;
}

}
