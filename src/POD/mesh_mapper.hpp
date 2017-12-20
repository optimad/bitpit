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

#ifndef __BITPIT_MESH_MAPPER_HPP__
#define __BITPIT_MESH_MAPPER_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_voloctree.hpp"

namespace bitpit {

namespace mapping{

    typedef adaption::Info Info;

}

class MeshMapper {



public:
# if BITPIT_ENABLE_MPI
    MeshMapper(MPI_Comm comm = MPI_COMM_WORLD);
# else
    MeshMapper();
# endif

    ~MeshMapper();
    MeshMapper(MeshMapper&& other) = default;

    const bitpit::PiercedStorage<bitpit::adaption::Info> & getMapping();

    void mapMeshes(bitpit::VolumeKernel * meshReference, bitpit::VolumeKernel * meshMapped);

protected:

#if BITPIT_ENABLE_MPI
    MPI_Comm                m_communicator; /**< MPI communicator */
#endif
    int                     m_rank;         /**< Local rank of process. */
    int                     m_nProcs;       /**< Number of processes. */

    bitpit::PiercedStorage<bitpit::adaption::Info> m_mapper;  /**< Mapping info for each cell of reference mesh. */

    void _mapMeshes(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped);
    void _mapMeshesSamePartition(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped);

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

};

}
#endif
