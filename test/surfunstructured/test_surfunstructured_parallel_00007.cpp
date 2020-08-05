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

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "bitpit_surfunstructured.hpp"

#include <mpi.h>

#include <exception>

using namespace bitpit;

// ========================================================================== //
// SUBTEST #001 Test surface partitioning                                     //
// ========================================================================== //
int subtest_001(int nProcs, int myRank)
{
    BITPIT_UNUSED(nProcs);

    //
    // Initialize patch
    //

    // Create the patch
    SurfUnstructured mesh(2, 3);
    mesh.setCommunicator(MPI_COMM_WORLD);

    // Read the file
    log::cout() << "Reading STL file..." << std::endl;

    if (myRank == 0) {
        const std::string &filename = "data/sphere.stl";
        mesh.importSTL(filename, true);
    }

    //  Build adjacencies
    mesh.deleteCoincidentVertices();
    mesh.buildAdjacencies();

    mesh.write("test00007_subtest_001_initial");

    //
    // First partitioning
    //
    log::cout() << "First partitioning..." << std::endl;

    std::unordered_map<long, int> firstPartitioningRanks;
    if (myRank == 0) {
        firstPartitioningRanks[ 2] = 1;
        firstPartitioningRanks[ 3] = 1;
        firstPartitioningRanks[ 4] = 1;
        firstPartitioningRanks[ 5] = 1;
        firstPartitioningRanks[10] = 1;
        firstPartitioningRanks[11] = 1;
        firstPartitioningRanks[12] = 1;
        firstPartitioningRanks[13] = 1;
        firstPartitioningRanks[25] = 1;
        firstPartitioningRanks[27] = 1;
        firstPartitioningRanks[37] = 1;
        firstPartitioningRanks[39] = 1;

        firstPartitioningRanks[16] = 1;
        firstPartitioningRanks[18] = 1;
        firstPartitioningRanks[21] = 1;
        firstPartitioningRanks[22] = 1;
        firstPartitioningRanks[24] = 1;
        firstPartitioningRanks[26] = 1;
        firstPartitioningRanks[36] = 1;
        firstPartitioningRanks[38] = 1;
        firstPartitioningRanks[40] = 1;
        firstPartitioningRanks[42] = 1;
        firstPartitioningRanks[44] = 1;
        firstPartitioningRanks[46] = 1;
    }

    mesh.partition(firstPartitioningRanks, false);

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00007_subtest_001_first_partitioning");

    mesh.exportSTL("test00007_subtest_001_first_partitioning."+std::to_string(myRank)+".stl", false, false);

    //
    // Delete and store ghosts
    //
    std::vector<long> ghostIds;
    std::vector<bitpit::Cell> ghostCells;
    ghostIds.reserve(mesh.getGhostCount());
    ghostCells.reserve(mesh.getGhostCount());
    const auto itbegin = mesh.ghostBegin();
    const auto itend = mesh.ghostEnd();
    for (auto it = itbegin; it != itend; it++){
        ghostIds.push_back(it->getId());
        ghostCells.push_back(*it);
    }
    mesh.deleteCells(ghostIds);

    //
    // Re-add ghost cells with different Ids
    //
    int owner = -1;
    if (myRank == 0){
        owner = 1;
    }
    else {
        owner = 0;
    }
    for (const Cell & ghost : ghostCells){
        mesh.addCell(ghost, owner, ghost.getId()+1000);
    }

    mesh.buildAdjacencies();

    //
    // Second false partitioning to reset partitioning information
    //
    log::cout() << "Second false partitioning..." << std::endl;

    std::unordered_map<long, int> secondPartitioningRanks;
    for (const Cell & cell : mesh.getCells()){
        if (cell.isInterior()){
            long cellId = cell.getId();
            secondPartitioningRanks[cellId] = myRank;
        }
    }

    mesh.partition(secondPartitioningRanks, false, false);

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00007_subtest_001_second_partitioning");

    return 0;
}

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    //======================================================================= //
    MPI_Init(&argc, &argv);

    // ====================================================================== //
    // Initialize the logger
    // ====================================================================== //
    int nProcs;
    int    myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    log::manager().initialize(log::COMBINED, true, nProcs, myRank);
    log::cout().setVisibility(log::GLOBAL);

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001(nProcs, myRank);
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
    MPI_Finalize();

    return status;
}
