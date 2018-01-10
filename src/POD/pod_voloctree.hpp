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

#ifndef __BITPIT_POD_VOLOCTREE_HPP__
#define __BITPIT_POD_VOLOCTREE_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "pod_kernel.hpp"
#include "bitpit_voloctree.hpp"

namespace bitpit {

class PODVolOctree: public PODKernel {

    friend class POD;

public:
# if BITPIT_ENABLE_MPI
    PODVolOctree(MPI_Comm comm = MPI_COMM_WORLD);
# else
    PODVolOctree();
# endif

    ~PODVolOctree() override;

    PODVolOctree(PODVolOctree&& other) = default;

private:

    VolumeKernel* createMesh() override;


    bitpit::PiercedStorage<bitpit::adaption::Info> mapMesh(bitpit::VolOctree * mesh);
    void mapMeshSamePartition(bitpit::VolOctree * mesh, bitpit::PiercedStorage<bitpit::adaption::Info> & mapper);

    pod::PODField mapPODFieldToPOD(const pod::PODField & field, const std::unordered_set<long> * targetCells);
    void mapPODFieldFromPOD(pod::PODField & field, const std::unordered_set<long> * targetCells, const pod::PODField & mappedField);

    PiercedStorage<double> mapFieldsToPOD(const PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds);
    void mapFieldsFromPOD(PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const PiercedStorage<double> & mappedFields,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds);

    PiercedStorage<bool> mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells);
    void mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells, PiercedStorage<bool> & mappedField);

    std::unordered_set<long> mapCellsToPOD(const std::unordered_set<long> * targetCells);

};

}
#endif
