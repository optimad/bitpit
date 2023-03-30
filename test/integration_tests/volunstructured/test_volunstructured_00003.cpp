/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2023 OPTIMAD engineering Srl
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
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_volunstructured.hpp"

#include <unordered_map>

using namespace bitpit;

/*!
* Subtest 001
*
* Testing evaluation of cell information of a 3D unstructured patch.
*/
int subtest_001()
{
    // Create the patch
    log::cout() << std::endl;
    log::cout() << "Creating 2D patch..." << std::endl;

#if BITPIT_ENABLE_MPI
    VolUnstructured patch_2D(2, MPI_COMM_NULL);
#else
    VolUnstructured patch_2D(2);
#endif
    patch_2D.getVTK().setName("test_volunstructured_00003_2D");

    patch_2D.setVertexAutoIndexing(false);

    patch_2D.addVertex({{0.00000000, 0.00000000, 0.00000000}},  1);
    patch_2D.addVertex({{0.00000000, 1.00000000, 0.00000000}},  2);
    patch_2D.addVertex({{1.00000000, 1.00000000, 0.00000000}},  3);
    patch_2D.addVertex({{1.00000000, 0.00000000, 0.00000000}},  4);
    patch_2D.addVertex({{0.00000000, 0.50000000, 0.00000000}},  5);
    patch_2D.addVertex({{1.00000000, 0.50000000, 0.00000000}},  6);

    patch_2D.addCell(ElementType::QUAD,     std::vector<long>({{ 1, 2, 3, 4}}));
    patch_2D.addCell(ElementType::PIXEL,    std::vector<long>({{ 1, 2, 4, 3}}));
    patch_2D.addCell(ElementType::TRIANGLE, std::vector<long>({{ 1, 2, 4}}));
    patch_2D.addCell(ElementType::POLYGON,  std::vector<long>({{ 6, 1, 5, 2, 3, 6, 4}}));

    // Write patch
    patch_2D.write();

    // Evaluate element information
    std::unordered_map<long, double> expectedAreas(patch_2D.getCellCount());
    expectedAreas[0] = 1.;
    expectedAreas[1] = 1.;
    expectedAreas[2] = 0.5;
    expectedAreas[3] = 1.;

    std::unordered_map<long, double> expectedSizes(patch_2D.getCellCount());
    expectedSizes[0] = 1.;
    expectedSizes[1] = 1.;
    expectedSizes[2] = 3. / (2. + std::sqrt(2.));
    expectedSizes[3] = 1.;

    std::vector<std::array<double, 3>> cellVertexCoordinates;

    log::cout() << std::endl;
    log::cout() << " Evaluate cell information..." << std::endl;
    log::cout() << std::endl;
    for (const Cell &cell : patch_2D.getCells()) {
        log::cout() << " Cell #" << cell.getId() << std::endl;

        long cellId = cell.getId();
        std::size_t nCellVertices = cell.getVertexCount();
        cellVertexCoordinates.resize(nCellVertices);
        patch_2D.getCellVertexCoordinates(cellId, cellVertexCoordinates.data());

        double area = cell.evalArea(cellVertexCoordinates.data());
        log::cout() << "   Area          = " << area << std::endl;
        log::cout() << "   Expected area = " << expectedAreas.at(cellId) << std::endl;
        if (!utils::DoubleFloatingEqual()(area, expectedAreas.at(cellId))) {
            log::cout() << "   Area doesn't match the expected value!" << std::endl;
            return 1;
        }

        double size = cell.evalSize(cellVertexCoordinates.data());
        log::cout() << "   Size          = " << size << std::endl;
        log::cout() << "   Expected size = " << expectedSizes.at(cellId) << std::endl;
        if (!utils::DoubleFloatingEqual()(size, expectedSizes.at(cellId))) {
            log::cout() << "   Size doesn't match the expected value!" << std::endl;
            return 1;
        }
    }

    return 0;
}

/*!
* Subtest 002
*
* Testing evaluation of cell information of a 3D unstructured patch.
*/
int subtest_002()
{
    // Create the patch
    log::cout() << std::endl;
    log::cout() << "Creating 3D patch..." << std::endl;

#if BITPIT_ENABLE_MPI
    VolUnstructured patch_3D(3, MPI_COMM_NULL);
#else
    VolUnstructured patch_3D(3);
#endif
    patch_3D.getVTK().setName("test_volunstructured_00003_3D");

    patch_3D.setVertexAutoIndexing(false);

    patch_3D.addVertex({{0.00000000, 0.00000000,  0.00000000}},  1);
    patch_3D.addVertex({{1.00000000, 0.00000000,  0.00000000}},  2);
    patch_3D.addVertex({{1.00000000, 1.00000000,  0.00000000}},  3);
    patch_3D.addVertex({{0.00000000, 1.00000000,  0.00000000}},  4);
    patch_3D.addVertex({{0.00000000, 0.00000000,  1.00000000}},  5);
    patch_3D.addVertex({{1.00000000, 0.00000000,  1.00000000}},  6);
    patch_3D.addVertex({{1.00000000, 1.00000000,  1.00000000}},  7);
    patch_3D.addVertex({{0.00000000, 1.00000000,  1.00000000}},  8);
    patch_3D.addVertex({{0.00000000, 0.00000000,  0.50000000}},  9);
    patch_3D.addVertex({{1.00000000, 0.00000000,  0.50000000}}, 10);
    patch_3D.addVertex({{1.00000000, 1.00000000,  0.50000000}}, 11);
    patch_3D.addVertex({{0.00000000, 1.00000000,  0.50000000}}, 12);

    patch_3D.addCell(ElementType::TETRA,      std::vector<long>({{1, 2, 4, 5}}));
    patch_3D.addCell(ElementType::WEDGE,      std::vector<long>({{5, 6, 8, 1, 2, 4}}));
    patch_3D.addCell(ElementType::PYRAMID,    std::vector<long>({{1, 2, 3, 4, 5}}));
    patch_3D.addCell(ElementType::VOXEL,      std::vector<long>({{1, 2, 4, 3, 5, 6, 8, 7}}));
    patch_3D.addCell(ElementType::HEXAHEDRON, std::vector<long>({{1, 2, 3, 4, 5, 6, 7, 8}}));
    patch_3D.addCell(ElementType::POLYHEDRON, std::vector<long>({{8,
                                                                   4, 1, 4, 3, 2,
                                                                   6, 2, 3, 11, 7, 6, 10,
                                                                   4, 5, 6, 7, 8,
                                                                   6, 1, 9, 5, 8, 12, 4,
                                                                   4, 1, 2, 10, 9,
                                                                   4, 5, 9, 10, 6,
                                                                   4, 3, 4, 12, 11,
                                                                   4, 11, 12, 8, 7
                                                                 }}));

    patch_3D.initializeAdjacencies();
    patch_3D.initializeInterfaces();

    // Write patch
    patch_3D.write();

    // Evaluate element information
    std::unordered_map<long, double> expectedVolumes(patch_3D.getCellCount());
    expectedVolumes[0] = 1. / 6.;
    expectedVolumes[1] = 1. / 2.;
    expectedVolumes[2] = 1. / 3.;
    expectedVolumes[3] = 1.;
    expectedVolumes[4] = 1.;
    expectedVolumes[5] = 1.;

    std::unordered_map<long, double> expectedSizes(patch_3D.getCellCount());
    expectedSizes[0] = 4. / (std::sqrt(2.) * (1. + std::sqrt(3.)));
    expectedSizes[1] = 1. / 2. / std::sqrt(2.);
    expectedSizes[2] = 1. / 3.;
    expectedSizes[3] = 1.;
    expectedSizes[4] = 1.;
    expectedSizes[5] = 1.;

    std::vector<std::array<double, 3>> cellVertexCoordinates;

    log::cout() << std::endl;
    log::cout() << " Evaluate cell information..." << std::endl;
    log::cout() << std::endl;
    for (const Cell &cell : patch_3D.getCells()) {
        log::cout() << " Cell #" << cell.getId() << std::endl;

        long cellId = cell.getId();
        std::size_t nCellVertices = cell.getVertexCount();
        cellVertexCoordinates.resize(nCellVertices);
        patch_3D.getCellVertexCoordinates(cellId, cellVertexCoordinates.data());

        double volume = cell.evalVolume(cellVertexCoordinates.data());
        log::cout() << "   Volume          = " << volume << std::endl;
        log::cout() << "   Expected volume = " << expectedVolumes.at(cellId) << std::endl;
        if (!utils::DoubleFloatingEqual()(volume, expectedVolumes.at(cellId))) {
            log::cout() << "   Volume doesn't match the expected value!" << std::endl;
            return 1;
        }

        double size = cell.evalSize(cellVertexCoordinates.data());
        log::cout() << "   Size            = " << size << std::endl;
        log::cout() << "   Expected size   = " << expectedSizes.at(cellId) << std::endl;
        if (!utils::DoubleFloatingEqual()(size, expectedSizes.at(cellId))) {
            log::cout() << "   Size doesn't match the expected value!" << std::endl;
            return 1;
        }
    }

    return 0;
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

    // Initialize the logger
    log::manager().initialize(log::MODE_COMBINE);

    // Run the subtests
    log::cout() << "Testing evaluation of cell information" << std::endl;

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

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
