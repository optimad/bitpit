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

using namespace std;
using namespace bitpit;

/*!
* Subtest 001
*
* Testing parallel refinement of a 2D octree.
*/
int subtest_001()
{
    /**<Instantation of a 2D para_tree object.*/
    ParaTree pablo(2);

    /**<Set NO 2:1 balance for the octree.*/
    uint32_t idx=0;
    pablo.setBalance(idx,false);

    /**<Compute the connectivity and write the para_tree.*/
    pablo.updateConnectivity();
    pablo.write("PabloParallelEmptyPartitions_serial");

    /**<GLOBAL REFINE*/
    pablo.adaptGlobalRefine();
    pablo.setMarker(2,1);
    pablo.adapt();

    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo.loadBalance();

//    if(pablo.getRank()==1 || pablo.getRank()==2){
//        uint32_t nocts = pablo.getNumOctants();
//        for(uint32_t i = 0; i < nocts; ++i)
//            pablo.setMarker(i,-1);
//    }

    uint32_t nocts = pablo.getNumOctants();
    for(uint32_t i = 0; i < nocts; ++i){
        if(pablo.getLevel(i) > 1)
            pablo.setMarker(i,-1);
    }

    pablo.adapt();

    pablo.adaptGlobalCoarse();

    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo.loadBalance();

    /**<Compute the connectivity and write the para_tree.*/
    pablo.updateConnectivity();
    pablo.write("PabloParallelEmptyPartitions_parallel");

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
    int    rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing parallel octree refinement" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
