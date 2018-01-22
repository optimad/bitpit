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

#include <cassert>

#include "volume_mapper.hpp"

namespace bitpit {

/**
 * \class VolumeMapper
 * \ingroup volumepatches
 *
 * \brief The VolumeMapper is the class to map two meshes.
 *
 * The VolumeMapper allows to map meshes of class VolumeKernel. The meshes are
 * defined as a reference mesh and a mapped mesh.
 *
 * The object can provide a direct mapper and a inverse mapper between only the
 * cells of the meshes.
 *
 * The information given by a mapper object is analogous to an adaptation info
 * to adapt the mapped mesh to the reference one (or vice versa in case of
 * inverse mapper).
 *
 * The two meshes have to be imperatively linked at declaration of the mapper
 * object.
 *
 * To compute the mapper the first time call initialize method. Then, if the
 * reference mesh OR the mapped mesh (one at a time) is adapted, the mapper
 * can be adapted together by passing the adaptation information to the prepare
 * and alter methods. A load-balancing procedure is not allowed, i.e. the mapper
 * has to be entirely recomputed after a load-balance of a mesh.
 *
 * Note that the VolumeMapper class is a pure virtual class, only objects of
 * its complete derived classes can be instantiated.
 *
 * The core method of the mapping procedure and the adapting mesh functions
 * have to be implemented in a derived class.
 */

/**
 * Constructor.
 *
 * \param[in] referencePatch is the reference mesh
 * \param[in] mappedPatch is the mapped mesh
 */
#if BITPIT_ENABLE_MPI
/**
 * \param[in] communicator is the MPI communicator
 */
VolumeMapper::VolumeMapper(VolumeKernel *referencePatch, VolumeKernel *mappedPatch, MPI_Comm communicator)
#else
VolumeMapper::VolumeMapper(VolumeKernel *referencePatch, VolumeKernel *mappedPatch)
#endif
    : m_referencePatch(referencePatch), m_mappedPatch(mappedPatch), m_mapping(1)
#if BITPIT_ENABLE_MPI
      , m_communicator(MPI_COMM_NULL)
#endif
{
#if BITPIT_ENABLE_MPI
    initializeCommunicator(communicator);
#endif
}

/**
 * Clear mapping members
 */
void VolumeMapper::clear()
{
    clearMapping();
    clearInverseMapping();
}

/**
 * Clear direct mapping
 */
void VolumeMapper::clearMapping()
{
    if (!m_mapping.getKernel()) {
        return;
    }

    m_mapping.unsetKernel(true);
}

/**
 * Clear inverse mapping
 */
void VolumeMapper::clearInverseMapping()
{
    if (!m_inverseMapping.getKernel()) {
        return;
    }

    m_inverseMapping.unsetKernel(false);
}

/**
 * Get direct mapping
 */
const bitpit::PiercedStorage<mapping::Info> & VolumeMapper::getMapping()
{
    return m_mapping;
}

/**
 * Get inverse mapping
 */
const bitpit::PiercedStorage<mapping::Info> & VolumeMapper::getInverseMapping()
{
    return m_inverseMapping;
}

/**
 * Initialize (compute for the first time) the mapper of the mapped mesh on the
 * reference mesh.
 *
 * The two meshes are set in the constructor of the derived object.
 *
 * \param[in] fillInverse if set to true the inverse mapped (reference mesh
 * to input mesh) will be filled
 */
void VolumeMapper::initialize(bool fillInverse)
{
    clear();

    _mapMeshes(fillInverse);
}

#if BITPIT_ENABLE_MPI
/**
 * Initialize the MPI communicator to be used for parallel communications.
 *
 * \param communicator is the communicator.
 */
void VolumeMapper::initializeCommunicator(MPI_Comm communicator)
{
    // Communication can be set just once
    if (isCommunicatorSet()) {
        throw std::runtime_error ("VolumeMapper communicator can be set just once");
    }

    // The communicator has to be valid
    if (communicator == MPI_COMM_NULL) {
        throw std::runtime_error ("VolumeMapper communicator is not valid");
    }

    // Create a copy of the user-specified communicator
    //
    // No library routine should use MPI_COMM_WORLD as the communicator;
    // instead, a duplicate of a user-specified communicator should always
    // be used.
    MPI_Comm_dup(communicator, &m_communicator);

    // Initialize parallel information
    MPI_Comm_size(m_communicator, &m_nProcs);
    MPI_Comm_rank(m_communicator, &m_rank);
}

/**
 * Returns the MPI communicator stored within LevelSetKernel.
 *
 * \result Return the MPI communicator.
 */
MPI_Comm VolumeMapper::getCommunicator() const
{
    return m_communicator;
}

/**
 * Checks if the communicator to be used for parallel communications has
 * already been set.
 *
 * \result Returns true if the communicator has been set, false otherwise.
 */
bool VolumeMapper::isCommunicatorSet() const
{
    return (getCommunicator() != MPI_COMM_NULL);
}

/**
 * Frees the MPI communicator associated to the mapper.
 */
void VolumeMapper::freeCommunicator()
{
    if (!isCommunicatorSet()) {
        return;
    }

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled) {
        return;
    }

    MPI_Comm_free(&m_communicator);
}
#endif

}
