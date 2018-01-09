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

    ~PODKernel();

    PODKernel(PODKernel&& other) = default;

protected:
    VolumeKernel    *m_meshPOD;     /**< Pointer to POD mesh*/

#if BITPIT_ENABLE_MPI
    MPI_Comm                m_communicator; /**< MPI communicator */
#endif
    int                     m_rank;         /**< Local rank of process. */
    int                     m_nProcs;       /**< Number of processes. */

    void clear();

    void setMesh(VolumeKernel*);
    VolumeKernel* getMesh();

    VolumeKernel* readMesh(const pod::SnapshotFile &snap);
    void restoreMesh(const pod::SnapshotFile &snap);

    virtual VolumeKernel* createMesh() = 0;

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

};

}

#endif
