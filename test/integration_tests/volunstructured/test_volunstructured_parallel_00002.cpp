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
* Fill the test mesh.
*/
void fillMesh(VolUnstructured *mesh)
{
    double radiusin  = 2.0;
    double radiusout = 5.0;

    double azimuthin  = 0.0;
    double azimuthout = 0.5 * BITPIT_PI;

    double heightbottom = -1.0;
    double heighttop    =  1.0;

    int nr = 4;
    int nt = 5;
    int nh = 7;

    double deltar = (radiusout - radiusin) / double(nr);
    double deltat = (azimuthout - azimuthin) / double(nt);
    double deltah = (heighttop - heightbottom) / double(nh);

    long nVertices = (nr + 1) * (nt + 1) * (nh + 1);
    long nCells    = nr * nt * nh;

    // Add vertices
    std::vector<std::array<double,3>> vertexCoords(nVertices);

    int counter = 0;
    for (int k = 0; k <= nh; ++k) {
        for (int j = 0; j <= nt; ++j) {
            for (int i = 0; i <= nr; ++i) {
                vertexCoords[counter][0] = (radiusin + i * deltar) * std::cos(azimuthin + j * deltat);
                vertexCoords[counter][1] = (radiusin + i * deltar) * std::sin(azimuthin + j * deltat);
                vertexCoords[counter][2] = (heightbottom + k * deltah);
                ++counter;
            }
        }
    }

    mesh->reserveVertices(nVertices);
    for (const std::array<double,3> &coords : vertexCoords) {
        mesh->addVertex(coords);
    }

    // Add cells
    mesh->reserveCells(nCells);

    ElementType cellType = ElementType::HEXAHEDRON;
    int nCellVertices = ReferenceElementInfo::getInfo(cellType).nVertices;

    for (int k = 0; k < nh; ++k) {
        for (int j = 0; j < nt; ++j) {
            for (int i = 0; i < nr; ++i) {
                std::unique_ptr<long[]> connectivity = std::unique_ptr<long[]>(new long[nCellVertices]);
                connectivity[0] = (nr + 1) * (nt + 1) * k + (nr + 1) * j + i;
                connectivity[1] = (nr + 1) * (nt + 1) * k + (nr + 1) * j + i + 1;
                connectivity[2] = (nr + 1) * (nt + 1) * k + (nr + 1) * (j + 1) + i + 1;
                connectivity[3] = (nr + 1) * (nt + 1) * k + (nr + 1) * (j + 1) + i;
                connectivity[4] = (nr + 1) * (nt + 1) * (k + 1) + (nr + 1) * j + i;
                connectivity[5] = (nr + 1) * (nt + 1) * (k + 1) + (nr + 1) * j + i + 1;
                connectivity[6] = (nr + 1) * (nt + 1) * (k + 1) + (nr + 1) * (j + 1) + i + 1;
                connectivity[7] = (nr + 1) * (nt + 1) * (k + 1) + (nr + 1) * (j + 1) + i;

                mesh->addCell(ElementType::HEXAHEDRON, std::move(connectivity));
            }
        }
    }
}

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
    log::cout() << "Creating patch..." << std::endl;

    std::unique_ptr<VolUnstructured> patch = std::unique_ptr<VolUnstructured>(new VolUnstructured(3, MPI_COMM_WORLD, 3));
    if (rank == 0) {
        fillMesh(patch.get());
    }
    patch->initializeAdjacencies();
    patch->initializeInterfaces();

    // Show patch info
    log::cout() << "Cell count: " << patch->getCellCount() << std::endl;
    log::cout() << "Internal cell count: " << patch->getInternalCellCount() << std::endl;
    log::cout() << "Ghost cell count: " << patch->getGhostCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch->getVertexCount() << std::endl;

    patch->write("serial_mesh");

    //
    // Partition the patch
    //

    // Evaluate cell ranks
    log::cout() << "Evaluating cell ranks..." << std::endl;

    std::unordered_map<long, int> cellRanks;
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

    patch->write("partitioned_mesh");

    //
    // Re-serialize the patch
    //

    // Evaluate cell ranks
    cellRanks.clear();
    for (auto itr = patch->internalCellBegin(); itr != patch->internalCellEnd(); ++itr) {
        cellRanks[itr.getId()] = 0;
    }

    // Re-serialize
    log::cout() << "Re-serializing the patch..." << std::endl;

    std::clock_t serializationStartTime = clock();

    patch->partition(cellRanks, true, true);

    std::clock_t serializationEndTime = clock();

    double serializationElapsed = double(serializationEndTime - serializationStartTime) / CLOCKS_PER_SEC;

    log::cout() << "    Serialization completed in " << serializationElapsed << " seconds" << std::endl;

    // Show patch info
    log::cout() << "Cell count: " << patch->getCellCount() << std::endl;
    log::cout() << "Internal cell count: " << patch->getInternalCellCount() << std::endl;
    log::cout() << "Ghost cell count: " << patch->getGhostCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch->getVertexCount() << std::endl;

    patch->write("reserialized_mesh");

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

    log::manager().initialize(log::MODE_COMBINE, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::VISIBILITY_GLOBAL);

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
