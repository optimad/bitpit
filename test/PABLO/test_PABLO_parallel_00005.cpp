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

#include <mpi.h>

#include "bitpit_PABLO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing adaptive mesh refinement (AMR) with full periodic boundary conditions.
*/
int subtest_001()
{
    /**<Instantation of a 3D pablo uniform object.*/
    PabloUniform tree(3);

    /** Set Periodic boundary conditions */
    tree.setPeriodic(0);
    tree.setPeriodic(2);
    tree.setPeriodic(4);

    /**<Refine globally one level and write the octree.*/
    tree.adaptGlobalRefine();

    /**<Refine globally one level and write the octree.*/
    tree.adaptGlobalRefine();

    /**<Partition the tree.*/
    tree.loadBalance();

    /**<Update tree connectivity.*/
    tree.updateConnectivity();

    /**<Tree infromation.*/
    uint32_t nOctants = tree.getNumOctants();

    /**<Write the tree.*/
    std::vector<double> octantIds(nOctants);
    for (uint32_t i = 0; i < nOctants; ++i) {
        octantIds[i] = i;
    }
    tree.writeTest("pablo_parallel_00005_01", octantIds);

    /**<Extract and print face neighsbors.*/
    log::cout() << "Number of Octants : " << nOctants << std::endl;
    log::cout() << "Extracting the four face neighsbours of each Octant " << std::endl;
    for (uint32_t i = 0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k=0; k<tree.getNfaces(); ++k) {
            tree.findNeighbours(i, k, 1, neighs, isGhost);
            log::cout() << " - For face " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print vertex neighsbors .*/
    log::cout() << "Extracting the four vertex neighsbours of each Octant " << std::endl;
    for (uint32_t i = 0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k=0; k<tree.getNnodes(); ++k) {
            tree.findNeighbours(i, k, 3, neighs, isGhost);
            log::cout() << " - For vertex " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print edge neighsbors .*/
    log::cout() << "Extracting the four edge neighsbours of each Octant " << std::endl;
    for (uint32_t i=0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k = 0; k<tree.getNedges(); ++k) {
            tree.findNeighbours(i, k, 2, neighs, isGhost);
            log::cout() << " - For edge " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    return 0;
}

/*!
* Subtest 002
*
* Testing adaptive mesh refinement (AMR) with partial periodic boundary
* conditions.
*/
int subtest_002()
{
    /**<Instantation of a 3D pablo uniform object.*/
    PabloUniform tree(3);

    /** Set Periodic boundary conditions */
    tree.setPeriodic(0);
    tree.setPeriodic(4);

    /**<Refine globally one level and write the octree.*/
    tree.adaptGlobalRefine();

    /**<Refine globally one level and write the octree.*/
    tree.adaptGlobalRefine();

    /**<Partition the tree.*/
    tree.loadBalance();

    /**<Update tree connectivity.*/
    tree.updateConnectivity();

    /**<Tree infromation.*/
    uint32_t nOctants = tree.getNumOctants();

    /**<Write the tree.*/
    std::vector<double> octantIds(nOctants);
    for (uint32_t i = 0; i < nOctants; ++i) {
        octantIds[i] = i;
    }
    tree.writeTest("pablo_parallel_00005_02", octantIds);

    /**<Extract and print face neighsbors.*/
    log::cout() << "Number of Octants : " << nOctants << std::endl;
    log::cout() << "Extracting the four face neighsbours of each Octant " << std::endl;
    for (uint32_t i = 0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k=0; k<tree.getNfaces(); ++k) {
            tree.findNeighbours(i, k, 1, neighs, isGhost);
            log::cout() << " - For face " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print vertex neighsbors .*/
    log::cout() << "Extracting the four vertex neighsbours of each Octant " << std::endl;
    for (uint32_t i = 0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k=0; k<tree.getNnodes(); ++k) {
            tree.findNeighbours(i, k, 3, neighs, isGhost);
            log::cout() << " - For vertex " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    /**<Extract and print edge neighsbors .*/
    log::cout() << "Extracting the four edge neighsbours of each Octant " << std::endl;
    for (uint32_t i=0; i < nOctants; ++i) {
        std::vector<uint32_t> neighs;
        std::vector<bool> isGhost;
        log::cout() << ". Octant index : " << i << std::endl;
        for (uint8_t k = 0; k<tree.getNedges(); ++k) {
            tree.findNeighbours(i, k, 2, neighs, isGhost);
            log::cout() << " - For edge " << (int)k << "; " << neighs.size() << " neighsbours: [ ";
            for (std::size_t n = 0; n < neighs.size(); ++n) {
                uint32_t neigh = neighs[n];
                if (!isGhost[n]) {
                    Octant *neighOctant = tree.getOctant(neigh);
                    log::cout() << neigh << " " << " : Center = " << tree.getCenter(neighOctant);
                } else {
                    Octant *neighOctant = tree.getGhostOctant(neigh);
                    log::cout() << neigh << " " << " (ghost) : Center = " << tree.getCenter(neighOctant);
                }
            }
            log::cout() << "]" << std::endl;
        }
    }

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initialize the logger
    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing parallel AMR with periodic boundary conditions." << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }

        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
