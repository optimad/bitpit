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
    SurfUnstructured mesh(2, 3, MPI_COMM_WORLD);

    // Read the file
    log::cout() << "Reading STL file..." << std::endl;

    if (myRank == 0) {
        const std::string &filename = "data/sphere.stl";
        mesh.importSTL(filename, true);
    }

    //  Build adjacencies
    mesh.initializeAdjacencies();

    mesh.write("test00005_subtest_001_initial");

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

        firstPartitioningRanks[16] = 2;
        firstPartitioningRanks[18] = 2;
        firstPartitioningRanks[21] = 2;
        firstPartitioningRanks[22] = 2;
        firstPartitioningRanks[24] = 2;
        firstPartitioningRanks[26] = 2;
        firstPartitioningRanks[36] = 2;
        firstPartitioningRanks[38] = 2;
        firstPartitioningRanks[40] = 2;
        firstPartitioningRanks[42] = 2;
        firstPartitioningRanks[44] = 2;
        firstPartitioningRanks[46] = 2;
    }

    mesh.partition(firstPartitioningRanks, false);

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00005_subtest_001_first_partitioning");


    //
    // Second partitioning
    //
    log::cout() << "Second partitioning..." << std::endl;

    std::unordered_map<long, int> secondPartitioningRanks;
    if (myRank == 0) {
        secondPartitioningRanks[ 8] = 1;
        secondPartitioningRanks[ 9] = 1;
        secondPartitioningRanks[14] = 1;
        secondPartitioningRanks[15] = 1;
        secondPartitioningRanks[19] = 1;
        secondPartitioningRanks[43] = 1;
        secondPartitioningRanks[47] = 1;

        secondPartitioningRanks[17] = 2;
        secondPartitioningRanks[20] = 2;
        secondPartitioningRanks[23] = 2;
        secondPartitioningRanks[28] = 2;
        secondPartitioningRanks[30] = 2;
        secondPartitioningRanks[31] = 2;
        secondPartitioningRanks[34] = 2;
        secondPartitioningRanks[35] = 2;
    } else if (myRank == 1) {
        secondPartitioningRanks[ 2] = 0;
        secondPartitioningRanks[ 3] = 0;
        secondPartitioningRanks[ 4] = 0;
        secondPartitioningRanks[ 5] = 0;
        secondPartitioningRanks[37] = 0;

        secondPartitioningRanks[11] = 2;
        secondPartitioningRanks[12] = 2;
        secondPartitioningRanks[13] = 2;
        secondPartitioningRanks[25] = 2;
        secondPartitioningRanks[27] = 2;
    } else if (myRank == 2) {
        secondPartitioningRanks[36] = 0;

        secondPartitioningRanks[18] = 1;
        secondPartitioningRanks[22] = 1;
        secondPartitioningRanks[38] = 1;
        secondPartitioningRanks[40] = 1;
        secondPartitioningRanks[42] = 1;
        secondPartitioningRanks[46] = 1;
    }

    mesh.partition(secondPartitioningRanks, false);

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00005_subtest_001_second_partitioning");

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
