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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_volcartesian.hpp"

using namespace bitpit;

struct ExpectedResults
{
    std::array<int, 3> nVertices;
    std::array<double, 3> vertexCoords_BL;
    std::array<double, 3> vertexCoords_TR;

    std::array<int, 3> nCells;
    std::array<double, 3> cellCoords_BL;
    std::array<double, 3> cellCoords_TR;
};

/*!
* Test a degenerate patch.
*
* \param dimension is the dimension of the patch
* \param expectedResults are the expected results
*/
int executeTest(int dimension, const ExpectedResults &expectedResults)
{
    log::cout() << std::endl;
    log::cout() << " Starting test..." << std::endl;
    log::cout() << std::endl;

    log::cout() << " Dimension: " << dimension << std::endl;

    std::array<double, 3> origin = {{10., 100., 1000.}};
    log::cout() << " Origin: " << origin << std::endl;

    std::array<double, 3> lengths = {{2., 2., 2.}};
    log::cout() << " Lengths: " << lengths << std::endl;

    log::cout() << std::endl;
    log::cout() << " Cell count: " << expectedResults.nCells << std::endl;

    log::cout() << std::endl;
    log::cout() << " Generating mesh... " << expectedResults.nCells << std::endl;
    VolCartesian patch(dimension, origin, lengths, expectedResults.nCells);

    log::cout() << std::endl;
    log::cout() << " Checking mesh... " << expectedResults.nCells << std::endl;

    for (int d = 0; d < 3; ++d) {
        int nCells = patch.getCellCount(d);
        log::cout() << "   - Cell count along direction " << d << ": " << nCells << std::endl;
        if (nCells != expectedResults.nCells[d]) {
            log::cout() << "       Expected cell count along direction " << d << ": " << expectedResults.nCells[d] << std::endl;
            throw std::runtime_error("The calculated number of cells differs from the expected one");
        }
    }

    for (int d = 0; d < 3; ++d) {
        int nVertices = patch.getVertexCount(d);
        log::cout() << "   - Vertex count along direction " << d << ": " << nVertices << std::endl;
        if (nVertices != expectedResults.nVertices[d]) {
            log::cout() << "       Expected vertex count along direction " << d << ": " << expectedResults.nVertices[d] << std::endl;
            throw std::runtime_error("The calculated number of vertices differs from the expected one");
        }
    }

# if defined(_MSC_VER) && BITPIT_ENABLE_DEBUG == 0
    std::array<double, 3> cellCoords_BL = patch.evalCellCentroid(patch.getCellCartesianId(patch.getCellLinearId({0,0,0})));
    //TODO MSVC puzzling error in optimized /O2 /O1 (Release versions): 
    // using VolCartesian::evalCellCentroid(long id) instead of VolCartesian::evalCellCentroid(array ijk)
    // get a crash with undefined error. 
#else 
    std::array<double, 3> cellCoords_BL = patch.evalCellCentroid(patch.getCellLinearId({0,0,0}));
#endif    
    log::cout() << "   - Bottom-left cell (BL cell) centroid: " << cellCoords_BL << std::endl;
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(cellCoords_BL[d], expectedResults.cellCoords_BL[d])) {
            log::cout() << "       Expected BL cell centroid " << d << ": " << expectedResults.cellCoords_BL << std::endl;
            throw std::runtime_error("The calculated centroid of the bottom-left cell differs from the expected one.");
        }
    }

    std::array<double, 3> cellCoords_TR = patch.evalCellCentroid(patch.getCellLinearId({expectedResults.nCells[0] - 1, expectedResults.nCells[1] - 1, expectedResults.nCells[2] - 1}));
    log::cout() << "   - Bottom-left cell (TR cell) centroid: " << cellCoords_TR << std::endl;
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(cellCoords_TR[d], expectedResults.cellCoords_TR[d])) {
            log::cout() << "       Expected TR cell centroid " << d << ": " << expectedResults.cellCoords_TR << std::endl;
            throw std::runtime_error("The calculated centroid of the bottom-left cell differs from the expected one.");
        }
    }

    std::array<double, 3> vertexCoords_BL = patch.evalVertexCoords(patch.getVertexLinearId({0,0,0}));
    log::cout() << "   - Bottom-left vertex (BL vertex) coordinates: " << vertexCoords_BL << std::endl;
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(vertexCoords_BL[d], expectedResults.vertexCoords_BL[d])) {
            log::cout() << "       Expected BL vertex coordinates " << d << ": " << expectedResults.vertexCoords_BL << std::endl;
            throw std::runtime_error("The calculated coordinates of the bottom-left vertex differs from the expected ones.");
        }
    }

    std::array<double, 3> vertexCoords_TR = patch.evalVertexCoords(patch.getVertexLinearId({expectedResults.nVertices[0] - 1, expectedResults.nVertices[1] - 1, expectedResults.nVertices[2] - 1}));
    log::cout() << "   - Bottom-left vertex (TR vertex) coordinates: " << vertexCoords_TR << std::endl;
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(vertexCoords_TR[d], expectedResults.vertexCoords_TR[d])) {
            log::cout() << "       Expected TR vertex coordinates " << d << ": " << expectedResults.vertexCoords_TR << std::endl;
            throw std::runtime_error("The calculated coordinates of the bottom-left vertex differs from the expected ones.");
        }
    }

    log::cout() << std::endl;
    log::cout() << " Test completed." << std::endl;

    return 0;
}

/*!
* Subtest 001
*
* Testing two-dimensional degenerate patches.
*/
int subtest_001()
{
    log::cout() << " Testing degenerate two-dimensional patches" << std::endl;
    log::cout() << std::endl;

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{0, 0, 0}};
        expectedResults.cellCoords_BL   = {{10, 100., 1000.}};
        expectedResults.cellCoords_TR   = {{10, 100., 1000.}};
        expectedResults.nVertices       = {{1, 1, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{10, 100., 1000.}};

        executeTest(2, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{10, 0, 0}};
        expectedResults.cellCoords_BL   = {{10.1, 100., 1000.}};
        expectedResults.cellCoords_TR   = {{11.9, 100., 1000.}};
        expectedResults.nVertices       = {{11, 1, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{12, 100., 1000.}};

        executeTest(2, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{0, 10, 0}};
        expectedResults.cellCoords_BL   = {{10., 100.1, 1000.}};
        expectedResults.cellCoords_TR   = {{10., 101.9, 1000.}};
        expectedResults.nVertices       = {{1, 11, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{10, 102., 1000.}};

        executeTest(2, expectedResults);
    }

    return 0;
}

/*!
* Subtest 002
*
* Testing three-dimensional degenerate patches.
*/
int subtest_002()
{
    log::cout() << " Testing degenerate three-dimensional patches" << std::endl;
    log::cout() << std::endl;

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{0, 0, 0}};
        expectedResults.cellCoords_BL   = {{10, 100., 1000.}};
        expectedResults.cellCoords_TR   = {{10, 100., 1000.}};
        expectedResults.nVertices       = {{1, 1, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{10, 100., 1000.}};

        executeTest(3, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{10, 0, 0}};
        expectedResults.cellCoords_BL   = {{10.1, 100., 1000.}};
        expectedResults.cellCoords_TR   = {{11.9, 100., 1000.}};
        expectedResults.nVertices       = {{11, 1, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{12, 100., 1000.}};

        executeTest(3, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{0, 10, 0}};
        expectedResults.cellCoords_BL   = {{10., 100.1, 1000.}};
        expectedResults.cellCoords_TR   = {{10., 101.9, 1000.}};
        expectedResults.nVertices       = {{1, 11, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{10, 102., 1000.}};

        executeTest(3, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{0, 0, 10}};
        expectedResults.cellCoords_BL   = {{10, 100., 1000.1}};
        expectedResults.cellCoords_TR   = {{10, 100., 1001.9}};
        expectedResults.nVertices       = {{1, 1, 11}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{10, 100., 1002.}};

        executeTest(3, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{10, 10, 0}};
        expectedResults.cellCoords_BL   = {{10.1, 100.1, 1000.}};
        expectedResults.cellCoords_TR   = {{11.9, 101.9, 1000.}};
        expectedResults.nVertices       = {{11, 11, 1}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{12, 102., 1000.}};

        executeTest(3, expectedResults);
    }

    {
        ExpectedResults expectedResults;
        expectedResults.nCells          = {{10, 0, 10}};
        expectedResults.cellCoords_BL   = {{10.1, 100., 1000.1}};
        expectedResults.cellCoords_TR   = {{11.9, 100., 1001.9}};
        expectedResults.nVertices       = {{11, 1, 11}};
        expectedResults.vertexCoords_BL = {{10, 100., 1000.}};
        expectedResults.vertexCoords_TR = {{12, 100., 1002.}};

        executeTest(3, expectedResults);
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
    log::cout() << "Testing basic features of Cartesian patches" << std::endl;

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
