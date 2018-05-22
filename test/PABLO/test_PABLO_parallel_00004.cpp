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

#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing generation of multi-layer halo.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    // Instantation of a 2D para_tree object
    ParaTree pablo(2);

    // Set 2:1 balance
    uint32_t idx=0;
    pablo.setBalance(idx,false);

    // Set the number of ghost layers
    pablo.setNofGhostLayers(3);

    // Compute the connectivity and write the octree
    pablo.computeConnectivity();
    pablo.write("PabloParallel001_iter0");

    // Refine globally two level and write the octree
    for (int iter=1; iter<4; iter++){
        pablo.adaptGlobalRefine();
        pablo.updateConnectivity();
        pablo.write("PabloParallel001_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    // Call loadBalance, the octree is now distributed over the processes
    pablo.loadBalance();

    log::cout() << "Nof ghost in " << pablo.getRank() << " = " << pablo.getNumGhosts() << std::endl;

    // Done
    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    // Initialize the logger
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing parallel octree dump and restore." << std::endl;

    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
