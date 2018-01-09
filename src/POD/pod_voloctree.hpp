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

public:
# if BITPIT_ENABLE_MPI
    PODVolOctree(MPI_Comm comm = MPI_COMM_WORLD);
# else
    PODVolOctree();
# endif

    ~PODVolOctree();

    PODVolOctree(PODVolOctree&& other) = default;

protected:

    VolumeKernel* createMesh();

};

}
#endif
