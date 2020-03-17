/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
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

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_PABLO.hpp"

using namespace bitpit;

// =================================================================================== //
/*!
    \example PABLO_example_00011.cpp

    \brief 3D adaptive mesh refinement (AMR) using PABLO

    This example creates a 3D Octree mesh on the cube domain [0,1]x[0,1]x[0,1].

    The domain is refined globally one time and periodic conditions are imposed
    on the boundaries of the domain.
    Then, the face, vertex and edge neighbors of the cells are found and print.

    <b>To run</b>: ./PABLO_example_00011 \n

    Thanks to <a href="https://www.surrey.ac.uk/people/maxime-delorme">Maxime Delorme</a> for this example.
*/
// =================================================================================== //

/**
 * Run the example.
 */
void run()
{
    /**<Instantation of a 3D pablo uniform object.*/
    PabloUniform pablo11(3);

    /** Set Periodic boundary conditions */
    pablo11.setPeriodic(0);
    pablo11.setPeriodic(2);
    pablo11.setPeriodic(4);

    /**<Compute the connectivity and write the octree.*/
    pablo11.computeConnectivity();
    pablo11.write("pablo000011.0");

    /**<Refine globally one level and write the octree.*/
    pablo11.adaptGlobalRefine();
    pablo11.updateConnectivity();
    pablo11.write("pablo000011.1");

    /**<Extract and print face neighbors.*/
    uint32_t nOct = pablo11.getNumOctants();
    log::cout() << "Number of Octants : " << nOct << std::endl;
    log::cout() << "Extracting the four face neighbours of each Octant " << std::endl;
    for (uint32_t iOct = 0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << iOct << std::endl;
        for (uint8_t iFace=0; iFace<pablo11.getNfaces(); ++iFace) {
            pablo11.findNeighbours(iOct, iFace, 1, neigh, isGhost);
            log::cout() << " - For face " << (int)iFace << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh : neigh) {
                log::cout() << iNeigh << " ";
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print vertex neighbors .*/
    log::cout() << "Extracting the four vertex neighbours of each Octant " << std::endl;
    for (uint32_t iOct = 0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << iOct << std::endl;
        for (uint8_t iVertex=0; iVertex<pablo11.getNnodes(); ++iVertex) {
            pablo11.findNeighbours(iOct, iVertex, 3, neigh, isGhost);
            log::cout() << " - For vertex " << (int)iVertex << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh: neigh) {
                log::cout() << iNeigh << " ";
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print edge neighbors .*/
    log::cout() << "Extracting the four edge neighbours of each Octant " << std::endl;
    for (uint32_t iOct=0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << iOct << std::endl;
        for (uint8_t iEdge = 0; iEdge<pablo11.getNedges(); ++iEdge) {
            pablo11.findNeighbours(iOct, iEdge, 2, neigh, isGhost);
            log::cout() << " - For edge " << (int)iEdge << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh: neigh) {
                log::cout() << iNeigh << " ";
            }
            log::cout() << "]" << std::endl;
        }
    }
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;
#endif

    // Initialize the logger
    log::manager().initialize(log::SEPARATE, false, nProcs, rank);
    log::cout() << fileVerbosity(log::NORMAL);
    log::cout() << consoleVerbosity(log::QUIET);

    // Run the example
    try {
        run();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
