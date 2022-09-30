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
* Testing parallel creation of a volunstructured patch with 2 halo layers
* The test is meant to be run with 2 ranks
*
*/
int subtest_001()
{
    // Initialize the logger
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::MODE_COMBINE, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::VISIBILITY_GLOBAL);

    log::cout() << "Testing partitioning of unstructured patches" << std::endl;

    //
    // Create the patch
    //
    int haloLayers = 2;
    std::unique_ptr<VolUnstructured> patch = std::unique_ptr<VolUnstructured>(new VolUnstructured(2, MPI_COMM_WORLD, haloLayers));
    patch->getVTK().setName("test_00005_partitioned_mesh");
    patch->setVertexAutoIndexing(false);

    // Set mesh size
    double domainSize = 1.;
    double cellRankOffset = domainSize / nProcs;
    long nofCellPerDirection = 16;
    double cellSize = domainSize / nofCellPerDirection;

    // Interior vertices
    long vertexId = 0;
    for (int j = 0; j < nofCellPerDirection / nProcs + 1; ++j) {
        for (int i = 0; i < nofCellPerDirection + 1; ++i) {
            std::array<double, 3> vertex = {{ i * cellSize, rank * cellRankOffset + j * cellSize , 0.00000000}};
            patch->addVertex(vertex,  vertexId);
            ++vertexId;
        }
    }

    int nofInteriorVertices = vertexId;

    // Interior cells
    long baseVertexId = 0;
    long cellId = 0;
    for (int j = 0; j < nofCellPerDirection / nProcs; ++j) {
        for (int i = 0; i < nofCellPerDirection; ++i) {
            patch->addCell(ElementType::QUAD,
                    std::vector<long>({{baseVertexId, baseVertexId + 1, baseVertexId + nofCellPerDirection + 2, baseVertexId + nofCellPerDirection + 1}}));
            ++baseVertexId;
            ++cellId;
        }
        ++baseVertexId;
    }

    // Ghost vertices
    int rankDirection = (rank ? -1 : 1);
    for (int j = 0; j < haloLayers; ++j) {
        for (int i = 0; i < nofCellPerDirection + 1; ++i) {
            std::array<double, 3> vertex = {{ i * cellSize, cellRankOffset + rankDirection * (j + 1) * cellSize , 0.00000000}};
            patch->addVertex(vertex,  vertexId);
            ++vertexId;
        }
    }

    // Ghost cells
    int layer = 0;
    baseVertexId += rank * (nofCellPerDirection + 1);
    for (int j = 0; j < haloLayers; ++j) {
        for (int i = 0; i < nofCellPerDirection; ++i) {
            std::vector<long> connectivity({{baseVertexId, baseVertexId + 1,
                baseVertexId + rankDirection * nofCellPerDirection + (-2 * rank + 2), baseVertexId + rankDirection *  nofCellPerDirection + (-2 * rank + 1)}});
            if ( layer == 0 && rank == 1) {
                connectivity = std::vector<long>({{baseVertexId, baseVertexId + 1,
                    baseVertexId - nofInteriorVertices + 1, baseVertexId - nofInteriorVertices}});
            }
            patch->addCell(ElementType::QUAD, connectivity, - rank + 1, layer);
            ++baseVertexId;
            ++cellId;
        }
        ++baseVertexId;
        ++layer;
    }

    // Update the patch
    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    patch->update();

    //Write the patch
    patch->write();

    auto sourceLists = patch->getGhostExchangeSources();
    auto targetLists = patch->getGhostExchangeTargets();

    // Check results
    bool enableOutput = false;
    if (enableOutput) {
        for (auto rankSourceList : sourceLists) {
            log::cout() << "On rank " << rank << " interior cells for rank " << rankSourceList.first << " are " << rankSourceList.second << std::endl;
        }

        for (auto rankTargetList : targetLists) {
            log::cout() << "On rank " << rank << " ghost cells on rank " << rankTargetList.first << " are " << rankTargetList.second << std::endl;
        }
    }

    std::vector<long> rankSources, rankTargets;
    if (rank == 0) {
        rankSources = std::vector<long>({{96, 112, 97, 113, 98, 114, 99, 115, 100, 116, 101, 117, 102, 118, 103, 119, 104, 120, 105, 121, 106, 122, 107, 123, 108, 124, 109, 125, 110, 126, 111, 127}});
        rankTargets = std::vector<long>({{128, 144, 129, 145, 130, 146, 131, 147, 132, 148, 133, 149, 134, 150, 135, 151, 136, 152, 137, 153, 138, 154, 139, 155, 140, 156, 141, 157, 142, 158, 143,159}});
        for (long source : sourceLists.at(1)) {
            if (find(rankSources.begin(), rankSources.end(), source) == rankSources.end()) {
                return 1;
            }
        }
        for (long target : targetLists.at(1)) {
            if (find(rankTargets.begin(), rankTargets.end(), target) == rankTargets.end()) {
                return 1;
            }
        }
    } else {
        rankSources = std::vector<long>({{0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31}});
        rankTargets = std::vector<long>({{144, 128, 145, 129, 146, 130, 147, 131, 148, 132, 149, 133, 150, 134, 151, 135, 152, 136, 153, 137, 154, 138, 155, 139, 156, 140, 157, 141, 158, 142, 159, 143}});
        for (long source : sourceLists.at(0)) {
            if (find(rankSources.begin(), rankSources.end(), source) == rankSources.end()) {
                return 1;
            }
        }
        for (long target : targetLists.at(0)) {
            if (find(rankTargets.begin(), rankTargets.end(), target) == rankTargets.end()) {
                return 1;
            }
        }
    }

    log::cout() << "Test passed!" << std::endl;

    return 0;

}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    // Run the subtests
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

