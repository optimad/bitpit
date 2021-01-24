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

#include <array>
#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_volunstructured.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing partitioning
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    // Create the patch
    std::unique_ptr<VolUnstructured> patch = std::unique_ptr<VolUnstructured>(new VolUnstructured(3, MPI_COMM_WORLD));
    patch->getVTK().setName("test_00003_partitioned_mesh");
    if (rank == 0) {
        patch->addVertex({{0.00000000, 0.00000000,  0.00000000}},  1);
        patch->addVertex({{1.00000000, 0.00000000,  0.00000000}},  2);
        patch->addVertex({{1.00000000, 1.00000000,  0.00000000}},  3);
        patch->addVertex({{0.00000000, 1.00000000, -0.75000000}},  4);
        patch->addVertex({{0.00000000, 0.00000000,  1.00000000}},  5);
        patch->addVertex({{1.00000000, 0.00000000,  1.00000000}},  6);
        patch->addVertex({{1.00000000, 1.00000000,  1.00000000}},  7);
        patch->addVertex({{0.00000000, 1.00000000,  1.00000000}},  8);
        patch->addVertex({{1.00000000, 0.00000000,  0.54678323}},  9);
        patch->addVertex({{0.50000000, 0.00000000,  1.00000000}}, 10);
        patch->addVertex({{0.00000000, 0.00000000,  0.54678323}}, 11);
        patch->addVertex({{1.00000000, 1.00000000,  0.54678323}}, 12);
        patch->addVertex({{0.50000000, 1.00000000,  1.00000000}}, 13);
        patch->addVertex({{0.00000000, 1.00000000,  0.54678323}}, 14);
        patch->addVertex({{1.00000000, 0.50000000,  1.00000000}}, 15);
        patch->addVertex({{0.00000000, 0.50000000,  1.00000000}}, 16);
        patch->addVertex({{0.51053620, 0.00000000,  0.34680184}}, 17);
        patch->addVertex({{0.36278402, 0.00000000,  0.68603230}}, 18);
        patch->addVertex({{0.69618860, 0.00000000,  0.73234294}}, 19);
        patch->addVertex({{0.51053620, 1.00000000,  0.34680184}}, 20);
        patch->addVertex({{0.36278402, 1.00000000,  0.68603230}}, 21);
        patch->addVertex({{0.69618860, 1.00000000,  0.73234294}}, 22);
        patch->addVertex({{1.00000000, 0.51053620,  0.34680184}}, 23);
        patch->addVertex({{1.00000000, 0.36278402,  0.68603230}}, 24);
        patch->addVertex({{1.00000000, 0.69618860,  0.73234294}}, 25);
        patch->addVertex({{0.00000000, 0.51053620,  0.34680184}}, 26);
        patch->addVertex({{0.00000000, 0.36278402,  0.68603230}}, 27);
        patch->addVertex({{0.00000000, 0.69618860,  0.73234294}}, 28);
        patch->addVertex({{0.50000000, 0.50000000,  1.00000000}}, 29);
        patch->addVertex({{0.75000000, 0.25000000,  1.00000000}}, 30);
        patch->addVertex({{0.25000000, 0.25000000,  1.00000000}}, 31);
        patch->addVertex({{0.00000000, 0.00000000, -1.50000000}}, 32);
        patch->addVertex({{1.00000000, 0.00000000, -1.00000000}}, 33);
        patch->addVertex({{1.00000000, 1.00000000, -0.50000000}}, 34);
        patch->addVertex({{0.00000000, 1.00000000, -1.00000000}}, 35);
        patch->addVertex({{0.00000000, 0.00000000, -2.25000000}}, 36);
        patch->addVertex({{1.00000000, 0.00000000, -1.66666666}}, 37);
        patch->addVertex({{1.00000000, 1.00000000, -2.12500000}}, 38);
        patch->addVertex({{0.00000000, 1.00000000, -2.00000000}}, 39);
        patch->addVertex({{0.00000000, 0.00000000, -3.00000000}}, 40);
        patch->addVertex({{1.00000000, 0.00000000, -3.00000000}}, 41);
        patch->addVertex({{1.00000000, 1.00000000, -3.00000000}}, 42);
        patch->addVertex({{0.00000000, 1.00000000, -3.00000000}}, 43);
        patch->addVertex({{0.00000000, 0.00000000, -4.00000000}}, 44);
        patch->addVertex({{1.00000000, 0.00000000, -4.00000000}}, 45);
        patch->addVertex({{1.00000000, 1.00000000, -4.00000000}}, 46);
        patch->addVertex({{0.00000000, 1.00000000, -4.00000000}}, 47);
        patch->addVertex({{0.50000000, 0.00000000, -4.00000000}}, 48);
        patch->addVertex({{1.00000000, 0.50000000, -4.00000000}}, 49);
        patch->addVertex({{0.50000000, 1.00000000, -4.00000000}}, 50);
        patch->addVertex({{0.33333333, 0.66666666, -4.00000000}}, 51);
        patch->addVertex({{0.00000000, 0.50000000, -4.00000000}}, 52);

        patch->addCell(ElementType::TETRA,      std::vector<long>({{29, 22, 25, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 27, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 28, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26, 23, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 25, 23}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{29, 21, 22, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 28, 27, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{23, 26, 29, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 23, 29, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 29, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 23, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 29, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{29, 25, 23, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 23,  3, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{31, 18, 27, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{22, 12, 25, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 11, 27}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 19, 24, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 14, 28}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{15, 22,  7, 25}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{15, 29, 22, 25}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{25, 12, 23, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{23,  4,  3,  2}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26, 18, 11, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 24, 23, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{21, 26, 14, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{28, 21, 13, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{30, 18, 19, 10}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 29, 18}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 18, 10}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{21,  8, 28, 13}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 29, 28}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 28,  8}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 22,  7}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 29, 22}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{15, 24, 29, 25}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{21, 13, 29, 22}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{28, 16, 27, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{31, 18,  5, 27}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{24, 30, 15, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{16, 31, 27, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 11,  5, 27}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{30, 19, 24,  6}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 1, 26, 11, 17}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{19,  9, 24,  6}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{14, 21,  8, 28}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 31,  5, 10}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{23, 12,  3, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{19, 30, 10,  6}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 7, 22, 12, 25}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{24, 30,  6, 15}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{16, 31,  5, 27}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{23,  9, 17,  2}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{26,  4, 14, 20}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 17, 24, 19}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{24, 17, 18, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{18, 24, 30, 19}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{30, 24, 18, 29}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 17, 26,  1}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 26, 17, 23}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 4,  2, 26,  1}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26,  2, 23}}));
        patch->addCell(ElementType::WEDGE,      std::vector<long>({{ 2, 4, 1, 33, 35, 32}}));
        patch->addCell(ElementType::WEDGE,      std::vector<long>({{ 4, 2, 3, 35, 33, 34}}));
        patch->addCell(ElementType::PYRAMID,    std::vector<long>({{36, 37, 38, 39, 33}}));
        patch->addCell(ElementType::PYRAMID,    std::vector<long>({{39, 38, 34, 35, 33}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{35, 36, 39, 33}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{32, 35, 33, 36}}));
        patch->addCell(ElementType::HEXAHEDRON, std::vector<long>({{42, 43, 40, 41, 38, 39, 36, 37}}));
        patch->addCell(ElementType::POLYHEDRON, std::vector<long>({{11,
                                                                     4, 42, 43, 40, 41,
                                                                     5, 52, 51, 50, 49, 48,
                                                                     3, 41, 40, 48,
                                                                     3, 42, 41, 49,
                                                                     3, 43, 42, 50,
                                                                     3, 40, 43, 52,
                                                                     3, 43, 51, 52,
                                                                     3, 43, 50, 51,
                                                                     3, 42, 49, 50,
                                                                     3, 41, 48, 49,
                                                                     3, 40, 52, 48
                                                                        }}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{52, 51, 47, 43}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{51, 50, 47, 43}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{44, 48, 52, 40}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{48, 45, 49, 41}}));
        patch->addCell(ElementType::TETRA,      std::vector<long>({{50, 49, 46, 42}}));
    }

    patch->initializeAdjacencies();
    patch->initializeInterfaces();

    //
    // Partition the patch
    //

    // Evaluate cell ranks
    log::cout() << "Evaluating cell ranks..." << std::endl;

    std::unordered_map<long, int> cellRanks;

    std::unordered_set<int> rrr;

    if (rank == 0) {
        int nProcs;
        MPI_Comm_size(patch->getCommunicator(), &nProcs);
        std::size_t nMaxCellsPerProc = std::ceil((double) patch->getInternalCellCount() / nProcs);

        std::size_t index = 0;
        for (auto itr = patch->internalCellBegin(); itr != patch->internalCellEnd(); ++itr) {
            int rank = std::floor((double) index / nMaxCellsPerProc);
            ++index;

            cellRanks[itr.getId()] = rank;
        }
    }

    // Partition
    log::cout() << "Partitioning the patch..." << std::endl;

    std::clock_t partitioningStartTime = clock();

    patch->partition(cellRanks, true, true);

    std::clock_t partitioningEndTime = clock();

    double partitioningElapsed = double(partitioningEndTime - partitioningStartTime) / CLOCKS_PER_SEC;

    log::cout() << "    Partition completed in " << partitioningElapsed << " seconds" << std::endl;

    // Show patch info
    log::cout() << "Cell count: " << patch->getCellCount() << std::endl;
    log::cout() << "Internal cell count: " << patch->getInternalCellCount() << std::endl;
    log::cout() << "Ghost cell count: " << patch->getGhostCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch->getVertexCount() << std::endl;

    patch->write();

    //
    // Test for cell data exchange
    //

    log::cout() << "Testing cell data exchange..." << std::endl;

    // Exchange data on cells
    DataCommunicator cellDataCommunicator(patch->getCommunicator());

    long cellData;
    size_t cellDataSize = sizeof(cellData);

    // Set and start the receives
    for (const auto &entry : patch->getGhostCellExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &cellList = entry.second;
        std::size_t listSize = cellList.size();

        cellDataCommunicator.setRecv(rank, listSize * cellDataSize);
        cellDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : patch->getGhostCellExchangeSources()) {
        const int rank = entry.first;
        const std::vector<long> &cellList = entry.second;
        std::size_t listSize = cellList.size();

        cellDataCommunicator.setSend(rank, listSize * cellDataSize);

        SendBuffer &buffer = cellDataCommunicator.getSendBuffer(rank);
        cellData = 0;
        for (std::size_t i = 0; i < listSize; ++i) {
            buffer << cellData;
            ++cellData;
        }

        cellDataCommunicator.startSend(rank);

    }

    // Receive data
    int nCompletedCellRecvs = 0;
    while (nCompletedCellRecvs < cellDataCommunicator.getRecvCount()) {
        int rank = cellDataCommunicator.waitAnyRecv();
        const std::vector<long> &cellList = patch->getGhostCellExchangeTargets(rank);
        std::size_t listSize = cellList.size();

        RecvBuffer &buffer = cellDataCommunicator.getRecvBuffer(rank);

        long expectedData = 0;
        for (std::size_t i = 0; i < listSize; ++i) {
            buffer >> cellData;


            if (cellData != expectedData) {
                log::cout() << "Error in cell data excahge. Expected data doen't match received data" << std::endl;
                return 1;
            }
            ++expectedData;
        }

        ++nCompletedCellRecvs;
    }

    cellDataCommunicator.waitAllSends();

    log::cout() << " Test completed successfully." << std::endl;

    //
    // Test for vertex data exchange
    //

    log::cout() << "Testing vertex data exchange..." << std::endl;

    // Exchange data on vertices
    DataCommunicator vertexDataCommunicator(patch->getCommunicator());

    long vertexData;
    size_t vertexDataSize = sizeof(vertexData);

    // Set and start the receives
    for (const auto &entry : patch->getGhostVertexExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &vertexList = entry.second;
        std::size_t listSize = vertexList.size();

        vertexDataCommunicator.setRecv(rank, listSize * vertexDataSize);

        vertexDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : patch->getGhostVertexExchangeSources()) {
        const int rank = entry.first;
        const std::vector<long> &vertexList = entry.second;
        std::size_t listSize = vertexList.size();

        vertexDataCommunicator.setSend(rank, listSize * vertexDataSize);

        SendBuffer &buffer = vertexDataCommunicator.getSendBuffer(rank);
        vertexData = 0;
        for (std::size_t i = 0; i < listSize; ++i) {
            buffer << vertexData;
            ++vertexData;
        }

        vertexDataCommunicator.startSend(rank);
    }

    // Receive data
    int nCompletedVertexRecvs = 0;
    while (nCompletedVertexRecvs < vertexDataCommunicator.getRecvCount()) {
        int rank = vertexDataCommunicator.waitAnyRecv();
        const std::vector<long> &vertexList = patch->getGhostVertexExchangeTargets(rank);
        std::size_t listSize = vertexList.size();

        RecvBuffer &buffer = vertexDataCommunicator.getRecvBuffer(rank);

        long expectedData = 0;
        for (std::size_t i = 0; i < listSize; ++i) {
            buffer >> vertexData;
            if (vertexData != expectedData) {
                log::cout() << "Error in vertex data excahge. Expected data doen't match received data" << std::endl;
                return 1;
            }
            ++expectedData;
        }

        ++nCompletedVertexRecvs;
    }

    vertexDataCommunicator.waitAllSends();

    log::cout() << " Test completed successfully." << std::endl;

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
    log::cout() << "Testing partitioning of unstructured patches" << std::endl;

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
