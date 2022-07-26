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
#include <ctime>
#include <chrono>

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_levelset.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_volcartesian.hpp"

const int SPACE_DIMENSION = 3;

/*!
 * Generate segmentation.
 *
 * \result The generated segmentation.
 */
std::unique_ptr<bitpit::SurfUnstructured> generateSegmentation()
{
    // Initialize segmentation
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> segmentation(new bitpit::SurfUnstructured(0, SPACE_DIMENSION - 1, MPI_COMM_NULL));
#else
    std::unique_ptr<bitpit::SurfUnstructured> segmentation(new bitpit::SurfUnstructured(0, SPACE_DIMENSION - 1));
#endif

    // Create vertex list
    //
    // Use non-consecutive vertex ids to test if the levelset can handle them.
    segmentation->addVertex({{ -0.5, -0.5,  0.5}}, 0);
    segmentation->addVertex({{  0.5, -0.5,  0.5}}, 1);
    segmentation->addVertex({{  0.5, -0.5, -0.5}}, 2);
    segmentation->addVertex({{ -0.5, -0.5, -0.5}}, 3);
    segmentation->addVertex({{ -0.5,  0.5,  0.5}}, 4);
    segmentation->addVertex({{  0.5,  0.5,  0.5}}, 5);
    segmentation->addVertex({{  0.5,  0.5, -0.5}}, 6);
    segmentation->addVertex({{ -0.5,  0.5, -0.5}}, 7);
    segmentation->addVertex({{  0.0,  0.6, -0.5}}, 8);

    // Create simplex list
    //
    // Use non-consecutive cell ids to test if the levelset can handle them.
    std::vector<long> cellConnect(6);

    cellConnect[0] = 3;
    cellConnect[1] = 2;
    cellConnect[2] = 1;
    cellConnect[3] = 0;
    segmentation->addCell(bitpit::ElementType::QUAD, cellConnect);

    cellConnect[0] = 6;
    cellConnect[1] = 5;
    cellConnect[2] = 1;
    cellConnect[3] = 2;
    segmentation->addCell(bitpit::ElementType::QUAD, cellConnect);

    cellConnect[0] = 7;
    cellConnect[1] = 3;
    cellConnect[2] = 0;
    cellConnect[3] = 4;
    segmentation->addCell(bitpit::ElementType::QUAD, cellConnect);

    cellConnect[0] = 8;
    cellConnect[1] = 7;
    cellConnect[2] = 4;
    segmentation->addCell(bitpit::ElementType::TRIANGLE, cellConnect);
    cellConnect[0] = 8;
    cellConnect[1] = 5;
    cellConnect[2] = 6;
    segmentation->addCell(bitpit::ElementType::TRIANGLE, cellConnect);
    cellConnect[0] = 8;
    cellConnect[1] = 4;
    cellConnect[2] = 5;
    segmentation->addCell(bitpit::ElementType::TRIANGLE, cellConnect);

    cellConnect[0] = 5;
    cellConnect[1] = 4;
    cellConnect[2] = 1;
    segmentation->addCell(bitpit::ElementType::TRIANGLE, cellConnect);
    cellConnect[0] = 4;
    cellConnect[1] = 0;
    cellConnect[2] = 1;
    segmentation->addCell(bitpit::ElementType::TRIANGLE, cellConnect);

    cellConnect[ 0] = 5;
    cellConnect[ 1] = 8;
    cellConnect[ 2] = 6;
    cellConnect[ 3] = 2;
    cellConnect[ 4] = 3;
    cellConnect[ 5] = 7;
    segmentation->addCell(bitpit::ElementType::POLYGON, cellConnect);

    segmentation->initializeAdjacencies();

    return segmentation;
}

/*!
 * Generate the Cartesian mesh.
 *
 * \result The generated Cartesian mesh.
 */
std::unique_ptr<bitpit::VolCartesian> generateCartesianMesh(const bitpit::SurfUnstructured &segmentation)
{
    std::array<double, 3> segmentationMin;
    std::array<double, 3> segmentationMax;
    segmentation.getBoundingBox(segmentationMin, segmentationMax);

    std::array<double, 3> length = 3. * (segmentationMax - segmentationMin);
    std::array<double, 3> origin = -0.5 * length;

    std::array<int,3> nc = {{64, 64, 64}};

    std::unique_ptr<bitpit::VolCartesian> mesh(new bitpit::VolCartesian(SPACE_DIMENSION, origin, length, nc));

    return mesh;
}

/*!
 * Generate the Octree mesh.
 *
 * \result The generated Octree mesh.
 */
std::unique_ptr<bitpit::VolOctree> generateOctreeMesh(const bitpit::SurfUnstructured &segmentation)
{
    std::array<double, 3> segmentationMin;
    std::array<double, 3> segmentationMax;
    segmentation.getBoundingBox(segmentationMin, segmentationMax);

    double length = 0.;
    for (int i = 0; i < 3; ++i) {
        length = std::max(length, 3 * (segmentationMax[i] - segmentationMin[i]));
    };

    std::array<double, 3> origin = - 0.5 * std::array<double, 3>{{length, length, length}};

    double dh = length / 64;

#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::VolOctree> mesh(new bitpit::VolOctree(SPACE_DIMENSION, origin, length, dh, MPI_COMM_NULL));
#else
    std::unique_ptr<bitpit::VolOctree> mesh(new bitpit::VolOctree(SPACE_DIMENSION, origin, length, dh));
#endif

    return mesh;
}

/*!
* Subtest 001
*
* Testing dense and sparse storage on a Cartesian mesh in default memory mode.
*/
int subtest_001()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing dense and sparse storage on an Cartesian mesh in default memory mode" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_008");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolCartesian> mesh = generateCartesianMesh(*segmentation);
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->initializeAdjacencies();
    mesh->update();

    // Initialize test
    long testCellId0 = 137120;
    long testCellId1 = 189856;
    long testCellId2 = 233888;

    int objectId = 0;

    //
    // Levelset - Sparse storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetSparse(bitpit::LevelSetStorageType::SPARSE);
    levelsetSparse.setPropagateSign(true);
    levelsetSparse.setSizeNarrowBand(0.25);
    levelsetSparse.setMesh(mesh.get());
    levelsetSparse.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using sprase storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startSparse = std::chrono::system_clock::now();
    levelsetSparse.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endSparse = std::chrono::system_clock::now();
    int elapsedTimeSparse = std::chrono::duration_cast<std::chrono::milliseconds>(endSparse - startSparse).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeSparse << " ms" << std::endl;

    levelsetSparse.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_cartesian_default_sparse");
    mesh->write();

    double sparseValue0 = levelsetSparse.getObject(objectId).getValue(testCellId0);
    double sparseValue1 = levelsetSparse.getObject(objectId).getValue(testCellId1);
    double sparseValue2 = levelsetSparse.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId0 << " is equal to " << sparseValue0 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId1 << " is equal to " << sparseValue1 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId2 << " is equal to " << sparseValue2 << std::endl;

    //
    // Levelset - Dense storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetDense(bitpit::LevelSetStorageType::DENSE);
    levelsetDense.setPropagateSign(true);
    levelsetDense.setSizeNarrowBand(0.25);
    levelsetDense.setMesh(mesh.get());
    levelsetDense.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using dense storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startDense = std::chrono::system_clock::now();
    levelsetDense.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endDense = std::chrono::system_clock::now();
    int elapsedTimeDense = std::chrono::duration_cast<std::chrono::milliseconds>(endDense - startDense).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeDense << " ms" << std::endl;

    levelsetDense.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_cartesian_default_dense");
    mesh->write();

    double denseValue0 = levelsetDense.getObject(objectId).getValue(testCellId0);
    double denseValue1 = levelsetDense.getObject(objectId).getValue(testCellId1);
    double denseValue2 = levelsetDense.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId0 << " is equal to " << denseValue0 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId1 << " is equal to " << denseValue1 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId2 << " is equal to " << denseValue2 << std::endl;

    //
    // Comparison
    //

    bitpit::log::cout() << " Checking levelset values" << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue0, denseValue0)) {
        bitpit::log::cout() << "  - Value obtained for test cell #0 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #0 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue1, denseValue1)) {
        bitpit::log::cout() << "  - Value obtained for test cell #1 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #1 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue2, denseValue2)) {
        bitpit::log::cout() << "  - Value obtained for test cell #2 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #2 sparse storage matches the one obtained using dense storage." << std::endl;

    return 0;
}

/*!
* Subtest 003
*
* Testing dense and sparse storage on a Cartesian mesh in light memory mode.
*/
int subtest_002()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing dense and sparse storage on an Cartesian mesh in light memory mode" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_008");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolCartesian> mesh = generateCartesianMesh(*segmentation);
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    // Initialize test
    long testCellId0 = 137120;
    long testCellId1 = 189856;
    long testCellId2 = 233888;

    int objectId = 0;

    //
    // Levelset - Sparse storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetSparse(bitpit::LevelSetStorageType::SPARSE);
    levelsetSparse.setSizeNarrowBand(0.25);
    levelsetSparse.setMesh(mesh.get());
    levelsetSparse.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using sprase storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startSparse = std::chrono::system_clock::now();
    levelsetSparse.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endSparse = std::chrono::system_clock::now();
    int elapsedTimeSparse = std::chrono::duration_cast<std::chrono::milliseconds>(endSparse - startSparse).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeSparse << " ms" << std::endl;

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->initializeAdjacencies();
    mesh->update();

    levelsetSparse.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_cartesian_light_sparse");
    mesh->write();

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    double sparseValue0 = levelsetSparse.getObject(objectId).getValue(testCellId0);
    double sparseValue1 = levelsetSparse.getObject(objectId).getValue(testCellId1);
    double sparseValue2 = levelsetSparse.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId0 << " is equal to " << sparseValue0 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId1 << " is equal to " << sparseValue1 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId2 << " is equal to " << sparseValue2 << std::endl;

    //
    // Levelset - Dense storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetDense(bitpit::LevelSetStorageType::DENSE);
    levelsetDense.setSizeNarrowBand(0.25);
    levelsetDense.setMesh(mesh.get());
    levelsetDense.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using dense storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startDense = std::chrono::system_clock::now();
    levelsetDense.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endDense = std::chrono::system_clock::now();
    int elapsedTimeDense = std::chrono::duration_cast<std::chrono::milliseconds>(endDense - startDense).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeDense << " ms" << std::endl;

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->initializeAdjacencies();
    mesh->update();

    levelsetDense.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_cartesian_light_dense");
    mesh->write();

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    double denseValue0 = levelsetDense.getObject(objectId).getValue(testCellId0);
    double denseValue1 = levelsetDense.getObject(objectId).getValue(testCellId1);
    double denseValue2 = levelsetDense.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId0 << " is equal to " << denseValue0 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId1 << " is equal to " << denseValue1 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId2 << " is equal to " << denseValue2 << std::endl;

    //
    // Comparison
    //

    bitpit::log::cout() << " Checking levelset values" << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue0, denseValue0)) {
        bitpit::log::cout() << "  - Value obtained for test cell #0 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #0 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue1, denseValue1)) {
        bitpit::log::cout() << "  - Value obtained for test cell #1 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #1 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue2, denseValue2)) {
        bitpit::log::cout() << "  - Value obtained for test cell #2 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #2 sparse storage matches the one obtained using dense storage." << std::endl;

    return 0;
}

/*!
* Subtest 003
*
* Testing dense and sparse storage on an Octree mesh.
*/
int subtest_003()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing dense and sparse storage on an Octree mesh" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_008");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolOctree> mesh = generateOctreeMesh(*segmentation);
    mesh->initializeAdjacencies();
    mesh->update();

    // Initialize test
    long testCellId0 = 173202;
    long testCellId1 = 174512;
    long testCellId2 = 182672;

    int objectId = 0;

    //
    // Levelset - Sparse storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetSparse(bitpit::LevelSetStorageType::SPARSE);
    levelsetSparse.setPropagateSign(true);
    levelsetSparse.setSizeNarrowBand(0.25);
    levelsetSparse.setMesh(mesh.get());
    levelsetSparse.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using sprase storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startSparse = std::chrono::system_clock::now();
    levelsetSparse.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endSparse = std::chrono::system_clock::now();
    int elapsedTimeSparse = std::chrono::duration_cast<std::chrono::milliseconds>(endSparse - startSparse).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeSparse << " ms" << std::endl;

    levelsetSparse.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_octree_sparse");
    mesh->write();

    double sparseValue0 = levelsetSparse.getObject(objectId).getValue(testCellId0);
    double sparseValue1 = levelsetSparse.getObject(objectId).getValue(testCellId1);
    double sparseValue2 = levelsetSparse.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId0 << " is equal to " << sparseValue0 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId1 << " is equal to " << sparseValue1 << std::endl;
    bitpit::log::cout() << " Sparse storage mode: levelset on cell " << testCellId2 << " is equal to " << sparseValue2 << std::endl;

    //
    // Levelset - Dense storage
    //

    // Initialize levelset
    bitpit::LevelSet levelsetDense(bitpit::LevelSetStorageType::DENSE);
    levelsetDense.setPropagateSign(true);
    levelsetDense.setSizeNarrowBand(0.25);
    levelsetDense.setMesh(mesh.get());
    levelsetDense.addObject(segmentation.get(), BITPIT_PI, objectId);

    bitpit::log::cout() << "Computing levelset using dense storage... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> startDense = std::chrono::system_clock::now();
    levelsetDense.compute( ) ;
    std::chrono::time_point<std::chrono::system_clock> endDense = std::chrono::system_clock::now();
    int elapsedTimeDense = std::chrono::duration_cast<std::chrono::milliseconds>(endDense - startDense).count();
    bitpit::log::cout() << "Computation compreted in " << elapsedTimeDense << " ms" << std::endl;

    levelsetDense.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_008_octree_dense");
    mesh->write();

    double denseValue0 = levelsetDense.getObject(objectId).getValue(testCellId0);
    double denseValue1 = levelsetDense.getObject(objectId).getValue(testCellId1);
    double denseValue2 = levelsetDense.getObject(objectId).getValue(testCellId2);

    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId0 << " is equal to " << denseValue0 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId1 << " is equal to " << denseValue1 << std::endl;
    bitpit::log::cout() << " Dense storage mode: levelset on cell " << testCellId2 << " is equal to " << denseValue2 << std::endl;

    //
    // Comparison
    //

    bitpit::log::cout() << " Checking levelset values" << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue0, denseValue0)) {
        bitpit::log::cout() << "  - Value obtained for test cell #0 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #0 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue1, denseValue1)) {
        bitpit::log::cout() << "  - Value obtained for test cell #0 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #1 sparse storage matches the one obtained using dense storage." << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(sparseValue2, denseValue2)) {
        bitpit::log::cout() << "  - Value obtained for test cell #2 using sparse storage doesn't match the one obtained using dense storage." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained for test cell #2 sparse storage matches the one obtained using dense storage." << std::endl;

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
    bitpit::log::manager().initialize(bitpit::log::COMBINED);

    // Run the subtests
    bitpit::log::cout() << "Testing storage types" << std::endl;

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

        status = subtest_003();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        bitpit::log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
