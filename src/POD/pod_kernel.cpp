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

#include <cassert>
#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>

#if BITPIT_ENABLE_MPI
#	include <mpi.h>
#endif

#include <lapacke.h>
#include <bitpit_voloctree.hpp>

#include "pod_kernel.hpp"

namespace bitpit {

/**
 * \class PODKernel
 * \ingroup POD
 * 
 * \brief The PODKernel class provides an interface to manage the mesh dependent members and functions of a POD object.
 * 
 * PODKernel is the base class to manage the mesh dependent members and functions of a POD object.
 */

/**
 * Creates a new PODKernel object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the pod object. MPI_COMM_WORLD is the default value.
 */
PODKernel::PODKernel(MPI_Comm comm)
# else
PODKernel::PODKernel()
# endif
:m_cellsVolume(1)
{

#if BITPIT_ENABLE_MPI
    m_communicator = MPI_COMM_NULL;
#endif  

    m_meshPOD = nullptr;

# if BITPIT_ENABLE_MPI
    initializeCommunicator(comm);
    MPI_Comm_size(m_communicator, &m_nProcs);
    MPI_Comm_rank(m_communicator, &m_rank);
#else
    m_rank = 0;
    m_nProcs = 1;
#endif

    m_dirtymap = true;

}

/**
 * Destructor of PODKernel
 */
PODKernel::~PODKernel()
{
    clear();
}

/**
 * Cleaning PODKernel
 */
void PODKernel::clear()
{
    delete m_meshPOD;

# if BITPIT_ENABLE_MPI
    freeCommunicator();
# endif
}

/**
 * Set the pointer to the POD mesh.
 *
 * param[in] The pointer to the POD mesh.
 */
void PODKernel::setMesh(bitpit::VolumeKernel* mesh)
{
    m_meshPOD = mesh;
}

/**
 * Get a pointer to the POD mesh.
 * 
 * \return The const pointer to the POD mesh.
 */
bitpit::VolumeKernel* PODKernel::getMesh()
{
    return m_meshPOD;
}

/**
 * Read the mesh from file.
 *
 * \param[in] snap Snapshot filename.
 */
VolumeKernel* PODKernel::readMesh(const pod::SnapshotFile &snap)
{
    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;
    std::string filename = std::string(snap.directory) + "/" + std::string(snap.name) + ".mesh";
    IBinaryArchive binaryReader(filename, dumpBlock);

    VolumeKernel *mesh = createMesh();

    mesh->restore(binaryReader.getStream());

    binaryReader.close();

    return mesh;
}

/**
 * Restore the mesh from file.
 *
 * \param[in] snap Snapshot filename.
 */
void PODKernel::restoreMesh(const pod::SnapshotFile &snap)
{
    m_meshPOD = readMesh(snap);
    evalCellsVolume();
}

/**
 * Compute the volume of the cells of the mesh and store them in a member.
 */
void PODKernel::evalCellsVolume()
{
    m_cellsVolume.setStaticKernel(&m_meshPOD->getCells());
    for (Cell & cell : m_meshPOD->getCells()){
        long id = cell.getId();
        m_cellsVolume[id] = m_meshPOD->evalCellVolume(id);
    }
}

/**
 * Get the volume of a cell of the mesh.
 *
 * \param[in] id Cell id
 */
double PODKernel::getCellVolume(long id)
{
    return m_cellsVolume[id];
}

/**
 * Get the volume of a cell of the mesh.
 *
 * \param[in] rawIndex Cell raw index
 */
double PODKernel::getRawCellVolume(long rawIndex)
{
    return m_cellsVolume.rawAt(rawIndex);
}

/**
 * Compute the mapping of an input mesh on the pod mesh.
 * \param[in] mesh Pointer to input mesh.
 * \param[in] fillInv If true even the inverse mapping is computed.
 */
void PODKernel::computeMapper(VolumeKernel * mesh, bool fillInv)
{
    m_meshmap.mapMeshes(m_meshPOD, mesh, fillInv);
    setMapperDirty(false);
}

/**
 * Prepare the pre-computed mapping of an input mesh to the POD mesh by
 * the info given before an adaptation of the input mesh (internal method).
 * \param[in] info Info vector result of adaptation prepare of the input mesh
 */
void PODKernel::prepareMapper(const std::vector<adaption::Info> & info)
{
    m_meshmap.mappingAdaptionPreparare(info, false);
    setMapperDirty(true);
}

/**
 * Update the mapping pre-computed of an input mesh to the POD mesh by
 * the info given after an adaptation of the input mesh (internal method).
 * \param[in] info Info vector result of adaptation of the input mesh
 * \param[in] fillInv If true even the inverse mapping is computed.
 */
void PODKernel::updateMapper(const std::vector<adaption::Info> & info, bool fillInv)
{
    m_meshmap.mappingAdaptionUpdate(info, false, fillInv);
    setMapperDirty(false);
}

/**
 * Get the stored pre-computed mapping.
 * \return mesh mapper.
 */
MeshMapper & PODKernel::getMeshMapper()
{
    return m_meshmap;
}

/**
 * Clear the mapping info.
 */
void PODKernel::clearMapper()
{
    m_meshmap.clear();
    m_dirtymap = true;
}
/**
 * Set if the mapper has to be recomputed.
 * param[in] Dirty mapping flag
 */
void PODKernel::setMapperDirty(bool dirty)
{
    m_dirtymap = dirty;
}

/**
 * Get if the mapper has to be recomputed.
 * return Is the mapping dirty?
 */
bool PODKernel::isMapperDirty()
{
    return m_dirtymap;
}

#if BITPIT_ENABLE_MPI
/**
 * Initializes the MPI communicator to be used for parallel communications.
 *
 * \param communicator is the communicator.
 */
void PODKernel::initializeCommunicator(MPI_Comm communicator)
{
    // Communication can be set just once
    if (isCommunicatorSet())
        throw std::runtime_error ("PODKernel communicator can be set just once");

    // The communicator has to be valid
    if (communicator == MPI_COMM_NULL)
        throw std::runtime_error ("PODKernel communicator is not valid");

    // Create a copy of the user-specified communicator
    //
    // No library routine should use MPI_COMM_WORLD as the communicator;
    // instead, a duplicate of a user-specified communicator should always
    // be used.
    MPI_Comm_dup(communicator, &m_communicator);
}

/**
 * Returns the MPI communicator stored within LevelSetKernel.
 * @return MPI communicator.
 */
MPI_Comm PODKernel::getCommunicator() const 
{
    return m_communicator;
}

/**
 * Checks if the communicator to be used for parallel communications has
 * already been set.
 *
 * \result Returns true if the communicator has been set, false otherwise.
 */
bool PODKernel::isCommunicatorSet() const 
{

    return (getCommunicator() != MPI_COMM_NULL);
}

/**
 * Frees the MPI communicator associated to the patch
 */
void PODKernel::freeCommunicator()
{
    if (!isCommunicatorSet()) 
        return;

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled)
        return;

    MPI_Comm_free(&m_communicator);
}
#endif

}
