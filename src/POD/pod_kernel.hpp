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

#ifndef __BITPIT_POD_KERNEL_HPP__
#define __BITPIT_POD_KERNEL_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "pod_common.hpp"
#include "mesh_mapper.hpp"
#include "bitpit_patchkernel.hpp"
#include "bitpit_IO.hpp"

namespace bitpit {

class PODKernel : public VTKBaseStreamer {

    friend class POD;

    public:
# if BITPIT_ENABLE_MPI
    PODKernel(MPI_Comm comm = MPI_COMM_WORLD);
# else
    PODKernel();
# endif

    virtual ~PODKernel();

    PODKernel(PODKernel&& other) = default;

protected:
    VolumeKernel            *m_meshPOD;     /**< Pointer to POD mesh*/
    PiercedStorage<double>   m_cellsVolume; /**< Cells volume of POD mesh*/

#if BITPIT_ENABLE_MPI
    MPI_Comm                m_communicator; /**< MPI communicator */
#endif
    int                     m_rank;         /**< Local rank of process. */
    int                     m_nProcs;       /**< Number of processes. */

    MeshMapper              m_meshmap;      /**< Mapping object TO/FROM pod mesh.*/

    bool                    m_dirtymap;     /**< True if mapping has to be recomputed/updated [to be set by set method]. */

    void clear();

    void setMesh(VolumeKernel*);
    VolumeKernel* getMesh();

    VolumeKernel* readMesh(const pod::SnapshotFile &snap);
    void restoreMesh(const pod::SnapshotFile &snap);

    void evalCellsVolume();
    double getCellVolume(long id);
    double getRawCellVolume(long rawIndex);

    void    computeMapping(const VolumeKernel * mesh, bool fillInv = true);
    MeshMapper & getMeshMapper();

    void setMappingDirty(bool dirty = true);
    bool isMappingDirty();

    virtual VolumeKernel* createMesh() = 0;

    virtual pod::PODField mapPODFieldToPOD(const pod::PODField & field, const std::unordered_set<long> * targetCells) = 0;
    virtual void mapPODFieldFromPOD(pod::PODField & field, const std::unordered_set<long> * targetCells, const pod::PODField & mappedField) = 0;

    virtual PiercedStorage<double> mapFieldsToPOD(const PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds) = 0;
    virtual void mapFieldsFromPOD(PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const PiercedStorage<double> & mappedFields,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds) = 0;

    virtual PiercedStorage<bool> mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells) = 0;
    virtual void mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells, PiercedStorage<bool> & mappedField) = 0;

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

};

}

#endif
