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

#include <ctime>
#include <chrono>

#include <bitpit_CG.hpp>
#include <bitpit_IO.hpp>
#include <bitpit_surfunstructured.hpp>

using namespace bitpit;

// Subtest 001
//
// Two-dimensional data exchange.
int subtest_001()
{
    int status = 0;
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    log::cout() << "** ================================================================= **" << std::endl;
    log::cout() << "** Subtest #001 - Data exchange                                      **" << std::endl;
    log::cout() << "** ================================================================= **" << std::endl;

    // Importing STL
    log::cout() << std::endl;
    log::cout() << "Importing STL..." << std::endl;

    std::unique_ptr<SurfUnstructured> surfaceMesh(new SurfUnstructured(2, 3, MPI_COMM_WORLD));

    int myRank = surfaceMesh->getRank();

    surfaceMesh->setExpert(true);
    if (myRank == 0) {
        surfaceMesh->importSTL("./data/buddha.stl");
        surfaceMesh->deleteCoincidentVertices();
    }
    surfaceMesh->initializeAdjacencies();
    surfaceMesh->getVTK().setName("skd_test_STL");
    surfaceMesh->write();

    {
        // Scope variables ------------------------------------------------------ //
        long nCells = surfaceMesh->getCellCount();

        // Evaluation of baricenter ----------------------------------------------//
        std::array<double, 3> baricenter = {{0, 0, 0}};
        for (const auto &cell : surfaceMesh->getCells()) {
            baricenter += surfaceMesh->evalCellCentroid(cell.getId());
        }
        baricenter = baricenter / ((double) nCells);

        // Partitioning ----------------------------------------------------------//
        log::cout() << "Mesh partitioning..." << std::endl;

        std::unordered_map<long, int> cellRanks;
        if (myRank == 0) {
            int nProcs;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

            for (const auto &cell : surfaceMesh->getCells()) {
                long cellId = cell.getId();

                int side_x = (surfaceMesh->evalCellCentroid(cellId)[0] > baricenter[0]) ? 0 : 1;
                int side_y = (surfaceMesh->evalCellCentroid(cellId)[1] > baricenter[1]) ? 0 : 1;
                int side_z = (surfaceMesh->evalCellCentroid(cellId)[2] > baricenter[2]) ? 0 : 1;

                int rank = -1;
                if (side_z == 0 && side_y == 0 && side_x == 0) {
                    rank = 0;
                } else if (side_z == 0 && side_y == 0 && side_x == 1) {
                    rank = 1;
                } else if (side_z == 0 && side_y == 1 && side_x == 0) {
                    rank = 2;
                } else if (side_z == 0 && side_y == 1 && side_x == 1) {
                    rank = 3;
                } else if (side_z == 1 && side_y == 0 && side_x == 0) {
                    rank = 4;
                } else if (side_z == 1 && side_y == 0 && side_x == 1) {
                    rank = 5;
                } else if (side_z == 1 && side_y == 1 && side_x == 0) {
                    rank = 6;
                } else if (side_z == 1 && side_y == 1 && side_x == 1) {
                    rank = 7;
                }
                rank = rank % nProcs;

                cellRanks.insert({cellId, rank});
            }
        }

        surfaceMesh->partition(cellRanks, false);

        // Write info ----------------------------------------------------------- //
        log::cout() << "    Number of vertices: " << surfaceMesh->getVertexCount() << std::endl;
        log::cout() << "    Number of elements : " << surfaceMesh->getCellCount() << std::endl;

        // Write mesh ----------------------------------------------------------- //
        log::cout() << "Writing mesh..." << std::endl;
        surfaceMesh->write();

    }

    //
    // Test for cell data exchange
    //

    log::cout() << "Testing cell data exchange..." << std::endl;

    // Exchange data on cells
    DataCommunicator cellDataCommunicator(surfaceMesh->getCommunicator());

    long cellData;
    size_t cellDataSize = sizeof(cellData);

    // Set and start the receives
    for (const auto &entry : surfaceMesh->getGhostCellExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &cellList = entry.second;
        std::size_t listSize = cellList.size();

        cellDataCommunicator.setRecv(rank, listSize * cellDataSize);
        cellDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : surfaceMesh->getGhostCellExchangeSources()) {
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
        const std::vector<long> &cellList = surfaceMesh->getGhostCellExchangeTargets(rank);
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
    DataCommunicator vertexDataCommunicator(surfaceMesh->getCommunicator());

    long vertexData;
    size_t vertexDataSize = sizeof(vertexData);

    // Set and start the receives
    for (const auto &entry : surfaceMesh->getGhostVertexExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &vertexList = entry.second;
        std::size_t listSize = vertexList.size();

        vertexDataCommunicator.setRecv(rank, listSize * vertexDataSize);

        vertexDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : surfaceMesh->getGhostVertexExchangeSources()) {
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
        const std::vector<long> &vertexList = surfaceMesh->getGhostVertexExchangeTargets(rank);
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

    return status;
}

// Subtest 002
//
// Three-dimensional data exchange.
int subtest_002()
{
    int status = 0;
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    log::cout() << "** ================================================================= **" << std::endl;
    log::cout() << "** Subtest #001 - Data exchange                                      **" << std::endl;
    log::cout() << "** ================================================================= **" << std::endl;

    // Importing STL
    log::cout() << std::endl;
    log::cout() << "Importing STL..." << std::endl;

    std::unique_ptr<SurfUnstructured> surfaceMesh(new SurfUnstructured(2, 3, MPI_COMM_WORLD));

    int myRank = surfaceMesh->getRank();

    surfaceMesh->setExpert(true);
    if (myRank == 0) {
        surfaceMesh->importSTL("./data/buddha.stl");
        surfaceMesh->deleteCoincidentVertices();
    }
    surfaceMesh->initializeAdjacencies();
    surfaceMesh->getVTK().setName("skd_test_STL");
    surfaceMesh->write();

    {
        // Scope variables ------------------------------------------------------ //
        long nCells = surfaceMesh->getCellCount();

        // Evaluation of baricenter ----------------------------------------------//
        std::array<double, 3> baricenter = {{0, 0, 0}};
        for (const auto &cell : surfaceMesh->getCells()) {
            baricenter += surfaceMesh->evalCellCentroid(cell.getId());
        }
        baricenter = baricenter / ((double) nCells);

        // Partitioning ----------------------------------------------------------//
        log::cout() << "Mesh partitioning..." << std::endl;

        std::unordered_map<long, int> cellRanks;
        if (myRank == 0) {
            int nProcs;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

            for (const auto &cell : surfaceMesh->getCells()) {
                long cellId = cell.getId();

                int side_x = (surfaceMesh->evalCellCentroid(cellId)[0] > baricenter[0]) ? 0 : 1;
                int side_y = (surfaceMesh->evalCellCentroid(cellId)[1] > baricenter[1]) ? 0 : 1;
                int side_z = (surfaceMesh->evalCellCentroid(cellId)[2] > baricenter[2]) ? 0 : 1;

                int rank = -1;
                if (side_z == 0 && side_y == 0 && side_x == 0) {
                    rank = 0;
                } else if (side_z == 0 && side_y == 0 && side_x == 1) {
                    rank = 1;
                } else if (side_z == 0 && side_y == 1 && side_x == 0) {
                    rank = 2;
                } else if (side_z == 0 && side_y == 1 && side_x == 1) {
                    rank = 3;
                } else if (side_z == 1 && side_y == 0 && side_x == 0) {
                    rank = 4;
                } else if (side_z == 1 && side_y == 0 && side_x == 1) {
                    rank = 5;
                } else if (side_z == 1 && side_y == 1 && side_x == 0) {
                    rank = 6;
                } else if (side_z == 1 && side_y == 1 && side_x == 1) {
                    rank = 7;
                }
                rank = rank % nProcs;

                cellRanks.insert({cellId, rank});
            }
        }

        surfaceMesh->partition(cellRanks, false);

        // Write info ----------------------------------------------------------- //
        log::cout() << "    Number of vertices: " << surfaceMesh->getVertexCount() << std::endl;
        log::cout() << "    Number of elements : " << surfaceMesh->getCellCount() << std::endl;

        // Write mesh ----------------------------------------------------------- //
        log::cout() << "Writing mesh..." << std::endl;
        surfaceMesh->write();

    }

    //
    // Test for cell data exchange
    //

    log::cout() << "Testing cell data exchange..." << std::endl;

    // Exchange data on cells
    DataCommunicator cellDataCommunicator(surfaceMesh->getCommunicator());

    long cellData;
    size_t cellDataSize = sizeof(cellData);

    // Set and start the receives
    for (const auto &entry : surfaceMesh->getGhostCellExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &cellList = entry.second;
        std::size_t listSize = cellList.size();

        cellDataCommunicator.setRecv(rank, listSize * cellDataSize);
        cellDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : surfaceMesh->getGhostCellExchangeSources()) {
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
        const std::vector<long> &cellList = surfaceMesh->getGhostCellExchangeTargets(rank);
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
    DataCommunicator vertexDataCommunicator(surfaceMesh->getCommunicator());

    long vertexData;
    size_t vertexDataSize = sizeof(vertexData);

    // Set and start the receives
    for (const auto &entry : surfaceMesh->getGhostVertexExchangeTargets()) {
        const int rank = entry.first;
        const std::vector<long> &vertexList = entry.second;
        std::size_t listSize = vertexList.size();

        vertexDataCommunicator.setRecv(rank, listSize * vertexDataSize);

        vertexDataCommunicator.startRecv(rank);
    }

    // Set the sends
    for (const auto &entry : surfaceMesh->getGhostVertexExchangeSources()) {
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
        const std::vector<long> &vertexList = surfaceMesh->getGhostVertexExchangeTargets(rank);
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

    return status;
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

    // ====================================================================== //
    // Initialize the logger
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);
#endif

    // Run the subtests
    int status = 0;
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }

        status = subtest_002();
        if (status != 0) {
            return (20 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
