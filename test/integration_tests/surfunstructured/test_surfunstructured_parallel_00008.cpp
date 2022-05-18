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
#include "bitpit_CG.hpp"
#include "bitpit_surfunstructured.hpp"

#include <mpi.h>

#include <exception>

using namespace bitpit;

/*!
* Subtest 001
*
* Testing initialization of a partitioned patch.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    //
    // Initialize patch
    //

    // Create the patch
    SurfUnstructured mesh(2, MPI_COMM_WORLD);

    // Build the adjacencies
    mesh.initializeAdjacencies();

    // Build the interfaces
    mesh.initializeInterfaces();

    // Fill the patch
    if (rank == 0) {
        mesh.addVertex({{0., 0., 0.00000000}},  0);
        mesh.addVertex({{0., 1., 0.00000000}},  1);
        mesh.addVertex({{1., 0., 0.00000000}},  2);
        mesh.addVertex({{1., 1., 0.00000000}},  3);
        mesh.addVertex({{2., 0., 0.00000000}},  4);
        mesh.addVertex({{2., 1., 0.00000000}},  5);
        mesh.addVertex({{3., 0., 0.00000000}},  6);
        mesh.addVertex({{3., 1., 0.00000000}},  7);
        mesh.addVertex({{4., 0., 0.00000000}},  8);
        mesh.addVertex({{4., 1., 0.00000000}},  9);
        mesh.addVertex({{5., 0., 0.00000000}}, 10);
        mesh.addVertex({{5., 1., 0.00000000}}, 11);

        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 0,  3,  1}}), 0, -1, 0);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 0,  2,  3}}), 0, -1, 1);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 2,  5,  3}}), 0, -1, 2);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 2,  4,  5}}), 0, -1, 3);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 4,  7,  5}}), 0, -1, 4);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 4,  6,  7}}), 0, -1, 5);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 6,  9,  7}}), 0, -1, 6);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 6,  8,  9}}), 0, -1, 7);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 11,  9}}), 1,  0, 8);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 10, 11}}), 1,  0, 9);
    } else if (rank == 1) {
        mesh.addVertex({{3., 0., 0.00000000}},  6);
        mesh.addVertex({{3., 1., 0.00000000}},  7);
        mesh.addVertex({{4., 0., 0.00000000}},  8);
        mesh.addVertex({{4., 1., 0.00000000}},  9);
        mesh.addVertex({{5., 0., 0.00000000}}, 10);
        mesh.addVertex({{5., 1., 0.00000000}}, 11);
        mesh.addVertex({{6., 0., 0.00000000}}, 12);
        mesh.addVertex({{6., 1., 0.00000000}}, 13);
        mesh.addVertex({{7., 0., 0.00000000}}, 14);
        mesh.addVertex({{7., 1., 0.00000000}}, 15);
        mesh.addVertex({{8., 0., 0.00000000}}, 16);
        mesh.addVertex({{8., 1., 0.00000000}}, 17);

        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 6,  9,  7}}),     0,  0,  8);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 6,  8,  9}}),     0,  0,  9);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 11,  9}}),     1, -1, 10);
        mesh.addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 10, 11}}),     1, -1, 11);
        mesh.addCell(ElementType::PIXEL,    std::vector<long>({{10, 12, 11, 13}}), 1, -1, 12);
        mesh.addCell(ElementType::PIXEL,    std::vector<long>({{12, 14, 13, 15}}), 1, -1, 13);
        mesh.addCell(ElementType::PIXEL,    std::vector<long>({{14, 16, 15, 17}}), 2,  0, 14);
    } else if (rank == 2) {
        mesh.addVertex({{ 6., 0., 0.00000000}}, 12);
        mesh.addVertex({{ 6., 1., 0.00000000}}, 13);
        mesh.addVertex({{ 7., 0., 0.00000000}}, 14);
        mesh.addVertex({{ 7., 1., 0.00000000}}, 15);
        mesh.addVertex({{ 8., 0., 0.00000000}}, 16);
        mesh.addVertex({{ 8., 1., 0.00000000}}, 17);
        mesh.addVertex({{ 9., 0., 0.00000000}}, 18);
        mesh.addVertex({{ 9., 1., 0.00000000}}, 19);
        mesh.addVertex({{10., 0., 0.00000000}}, 20);
        mesh.addVertex({{10., 1., 0.00000000}}, 21);

        mesh.addCell(ElementType::PIXEL, std::vector<long>({{12, 14, 13, 15}}), 1,  0, 13);
        mesh.addCell(ElementType::PIXEL, std::vector<long>({{14, 16, 15, 17}}), 2, -1, 14);
        mesh.addCell(ElementType::PIXEL, std::vector<long>({{16, 18, 17, 19}}), 2, -1, 15);
        mesh.addCell(ElementType::PIXEL, std::vector<long>({{18, 20, 19, 21}}), 2, -1, 16);
    }

    mesh.update();

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00007_subtest_001_initial");


    //
    // First partitioning
    //
    log::cout() << "First partitioning..." << std::endl;

    std::unordered_map<long, int> firstPartitioningRanks;
    if (rank == 0) {
        firstPartitioningRanks[ 6] = 1;
        firstPartitioningRanks[ 7] = 1;
    } else if (rank == 1) {
        firstPartitioningRanks[13] = 2;
    }

    mesh.partition(firstPartitioningRanks, false);

    for (Cell & cell : mesh.getCells()){
        int nCellInterfaces = cell.getInterfaceCount();
        const long *cellInterfaces = cell.getInterfaces();

        std::cout << "#" << rank << " cell id " << cell.getId();
        for (int k = 0; k < nCellInterfaces; ++k) {
            std::cout << " , interface " << cellInterfaces[k];
        }
        std::cout << std::endl;
    }

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00007_subtest_001_first_partitioning");


    //
    // Second partitioning
    //
    log::cout() << "Second partitioning..." << std::endl;

    std::unordered_map<long, int> secondPartitioningRanks;
    if (rank == 0) {
        secondPartitioningRanks[ 2] = 1;
        secondPartitioningRanks[ 3] = 1;
        secondPartitioningRanks[ 4] = 2;
        secondPartitioningRanks[ 5] = 2;
    } else if (rank == 1) {
        secondPartitioningRanks[ 8] = 0;
        secondPartitioningRanks[ 9] = 0;
        secondPartitioningRanks[12] = 2;
    } else if (rank == 2) {
        secondPartitioningRanks[13] = 0;
        secondPartitioningRanks[14] = 1;
        secondPartitioningRanks[16] = 0;
    }

    mesh.partition(secondPartitioningRanks, false);

    for (Cell & cell : mesh.getCells()){
        int nCellInterfaces = cell.getInterfaceCount();
        const long *cellInterfaces = cell.getInterfaces();

        std::cout << "#" << rank << " cell id " << cell.getId();
        for (int k = 0; k < nCellInterfaces; ++k) {
            std::cout << " , interface " << cellInterfaces[k];
        }
        std::cout << std::endl;
    }

    // Write mesh
    log::cout() << "Writing mesh..." << std::endl;

    mesh.write("test00007_subtest_001_second_partitioning");

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    // Initialize the logger
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::MODE_COMBINE, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::VISIBILITY_GLOBAL);

    // Run the subtests
    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();

    return status;
}
