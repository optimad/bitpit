/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

    /**
     * Default copy constructor.
     * \param[in] other Input PODVolOctree object
     */
    PODVolOctree(PODVolOctree&& other) = default;

private:

    std::unique_ptr<VolumeKernel> createMesh() override;

    std::unique_ptr<VolumeMapper> _computeMapper(const VolumeKernel * mesh,  bool fillInv) override;

    bitpit::PiercedStorage<bitpit::adaption::Info> mapMesh(bitpit::VolOctree * mesh);
    void mapMeshSamePartition(bitpit::VolOctree * mesh, bitpit::PiercedStorage<bitpit::adaption::Info> & mapper);

    pod::PODField mapPODFieldToPOD(const pod::PODField & field, const std::unordered_set<long> * targetCells) override;
    void mapPODFieldFromPOD(pod::PODField & field, const std::unordered_set<long> * targetCells, const pod::PODField & mappedField) override;

    PiercedStorage<double> mapFieldsToPOD(const PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds) override;
    void mapFieldsFromPOD(PiercedStorage<double> & fields, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells,
            const PiercedStorage<double> & mappedFields,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds) override;

    PiercedStorage<bool> mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells) override;
    void mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh, const std::unordered_set<long> * targetCells, PiercedStorage<bool> & mappedField) override;

    std::unordered_set<long> mapCellsToPOD(const std::unordered_set<long> * targetCells) override;

    void adaptMeshToMesh(const VolumeKernel* meshToAdapt, const VolumeKernel * meshReference) override;

# if BITPIT_ENABLE_MPI
    void communicatePODField(const pod::PODField & field, std::map<int, std::map<long, bool> > & dataBrec, std::map<int, std::map<long, std::vector<double> > > & dataSrec, std::map<int, std::map<long, std::vector<std::array<double,3> > > > & dataVrec, std::map<int, std::map<long, double> > & volrec);
    void communicatePODFieldFromPOD(const pod::PODField & field, std::map<int, std::map<long, bool> > & dataBrec, std::map<int, std::map<long, std::vector<double> > > & dataSrec, std::map<int, std::map<long, std::vector<std::array<double,3> > > > & dataVrec, std::map<int, std::map<long, double> > & volrec);
    void communicateBoolField(const PiercedStorage<bool> & field, std::map<int, std::map<long, bool> > & dataBrec);
    void communicateField(const PiercedStorage<double> & field, const VolumeKernel * mesh, std::map<int, std::map<long, std::vector<double> > > & datarec, std::map<int, std::map<long, double> > & volrec);
    void communicateFieldFromPOD(const PiercedStorage<double> & field, const VolumeKernel * mesh, std::map<int, std::map<long, std::vector<double> > > & datarec, std::map<int, std::map<long, double> > & volrec);
#endif

};

}
#endif
