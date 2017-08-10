/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
// Evaluation of the closest cell of a surface patch
int subtest_001()
{
    int status = 0;
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Subtest #001 - Evaluation of the closest cell of a surface patch  **" << endl;
    log::cout() << "** ================================================================= **" << endl;

    // Importing STL
    log::cout() << std::endl;
    log::cout() << "Importing STL..." << std::endl;

    std::unique_ptr<SurfUnstructured> surfaceMesh(new SurfUnstructured (0, 2, 3));
    surfaceMesh->setExpert(true);
    surfaceMesh->importSTL("./data/buddha.stl");
    surfaceMesh->deleteCoincidentVertices();
    surfaceMesh->buildAdjacencies();
    surfaceMesh->getVTK().setName("skd_test_STL");
    surfaceMesh->write();

    log::cout() << "    Number of vertices: " << surfaceMesh->getVertexCount() << std::endl;
    log::cout() << "    Number of elements : " << surfaceMesh->getCellCount() << std::endl;

    // Build skd-tree
    log::cout() << std::endl;
    log::cout() << "Building skd-tree..." << std::endl;

    int elapsed_initalization = 0;
    start = std::chrono::system_clock::now();

    SurfaceSkdTree searchTree(surfaceMesh.get());
    searchTree.build();

    end = std::chrono::system_clock::now();
    elapsed_initalization = chrono::duration_cast<chrono::milliseconds>(end-start).count();

    long treeDepth = searchTree.evalMaxDepth();
    log::cout() << "    Maximum tree depth................. " << treeDepth << std::endl;

    long expectedTreeDepth = searchTree.evalMaxDepth();
    log::cout() << "    Expected maximum tree depth........ " << expectedTreeDepth << std::endl;

    if (treeDepth != expectedTreeDepth) {
        log::cout() << std::endl;
        log::cout() << "    <<< Maximum tree depth doesn't match expected value >>>>" << std::endl;

        status = 1;
    }

    assert(status == 0);
    if (status != 0) {
        return status;
    }

    log::cout() << std::endl;
    log::cout() << "    Elapsed time for skd-tree initialization " << elapsed_initalization << " ms" << std::endl;

    // Evaluate the closest cells
    log::cout() << std::endl;
    log::cout() << "Evaluation of the distance..." << std::endl;

    std::vector<std::array<double, 3>> points;
    points.push_back({{0., 0., 0.}});
    points.push_back({{1., 0., 0.}});
    points.push_back({{0., 1., 0.}});
    points.push_back({{0., 0., 1.}});
    points.push_back({{0.11111, 0.33333, 0.66666}});

    std::vector<long> expectedIds;
    expectedIds.push_back(8264);
    expectedIds.push_back(171212);
    expectedIds.push_back(5934);
    expectedIds.push_back(8264);
    expectedIds.push_back(5934);

    int elapsed_distance_evaluation = 0.;
    start = std::chrono::system_clock::now();

    std::size_t nPoints = points.size();
    for (std::size_t k = 0; k < nPoints; ++k) {
        long cellId;
        double cellDistance;
        const std::array<double, 3> &point = points[k];
        searchTree.findPointClosestCell(point, &cellId, &cellDistance);

        log::cout() << std::endl;
        log::cout() << "    Point .................... (" << point[0] << ", " <<  point[1] << ", " << point[2] << ")" << std::endl;
        log::cout() << "    Closest cell ............. " << cellId << std::endl;
        log::cout() << "    Expected closest cell .... " << expectedIds[k] << std::endl;
        log::cout() << "    Distance ................. " << cellDistance << std::endl;

        if (cellId != expectedIds[k]) {
            log::cout() << std::endl;
            log::cout() << "    <<< Closest cell doesn't match expected value >>>>" << std::endl;

            status = 1;
        }

        assert(status == 0);
    }

    end = std::chrono::system_clock::now();
    elapsed_distance_evaluation = chrono::duration_cast<chrono::milliseconds>(end-start).count();
    log::cout() << std::endl;
    log::cout() << "    Elapsed time for distance evaluation " << elapsed_distance_evaluation << " ms" << std::endl;

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

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    int status = 0;
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
