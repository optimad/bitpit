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

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"

using namespace std;
using namespace bitpit;

// =================================================================================== //
void testParallel01() {

    /**<Instantation and setup of a default (named bitpit) logfile.*/
    int nproc;
    int    rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm,&nproc);
    MPI_Comm_rank(comm,&rank);
#else
    nproc = 1;
    rank = 0;
#endif
    log::manager().initialize(log::SEPARATE, false, nproc, rank);
    log::cout() << fileVerbosity(log::NORMAL);
    log::cout() << consoleVerbosity(log::QUIET);

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

#if BITPIT_ENABLE_MPI==1
    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo.loadBalance();
#endif

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

#if BITPIT_ENABLE_MPI==1
    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo.loadBalance();
#endif

    /**<Compute the connectivity and write the para_tree.*/
    pablo.updateConnectivity();
    pablo.write("PabloParallelEmptyPartitions_parallel");

    return ;
}

// =================================================================================== //
int main( int argc, char *argv[] ) {

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc, &argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    /**<Calling Pablo Test routines*/
    testParallel01() ;

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
