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

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_levelset.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_volcartesian.hpp"

/*!
 * Generate segmentation.
 *
 * \result The generated segmentation.
 */
std::unique_ptr<bitpit::SurfUnstructured> generateSegmentation()
{
    const double R = 1.;
    const long N = 32;
    const double dtheta = 2. * BITPIT_PI / ((double) N);

    // Initialize segmentation
    int dimensions = 2;

#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> segmentation(new bitpit::SurfUnstructured(0, dimensions - 1, MPI_COMM_NULL));
#else
    std::unique_ptr<bitpit::SurfUnstructured> segmentation(new bitpit::SurfUnstructured(0, dimensions - 1));
#endif

    // Create vertex list
    //
    // Use non-consecutive vertex ids to test if the levelset can handle them.
    const long vertexIdOffset = 101;
    const long vertexIdStride = 2;

    std::array<double, 3> point;
    point[2] = 0.0;
    for (long i = 0; i < N; ++i) {
        double theta = i * dtheta;
        point[0] = R * cos(theta);
        point[1] = R * sin(theta);
        segmentation->addVertex(point, vertexIdOffset + vertexIdStride * i);
    }

    // Create simplex list
    //
    // Use non-consecutive cell ids to test if the levelset can handle them.
    const long cellIdOffset = 202;
    const long cellIdStride = 3;

    for (long i = 0; i < N; ++i) {
        std::vector<long> cellConnect(2, bitpit::Element::NULL_ID);
        cellConnect[0] = vertexIdOffset + vertexIdStride * i;
        cellConnect[1] = vertexIdOffset + vertexIdStride * ((i + 1) % N);
        segmentation->addCell(bitpit::ElementType::LINE, cellConnect, cellIdOffset + cellIdStride * i);
    }

    segmentation->initializeAdjacencies();

    return segmentation;
}

/*!
 * Generate mesh.
 *
 * \result The generated mesh.
 */
std::unique_ptr<bitpit::VolCartesian> generateMesh(const bitpit::SurfUnstructured &segmentation)
{
    int dimensions = 2;

    std::array<double, 3> meshMin, meshMax;
    segmentation.getBoundingBox(meshMin, meshMax);

    std::array<double, 3> delta = meshMax - meshMin;
    meshMin -= 0.1 * delta;
    meshMax += 0.1 * delta;
    delta = meshMax - meshMin;

    std::array<int,3> nc = {{64, 64, 0}};

    std::unique_ptr<bitpit::VolCartesian> mesh(new bitpit::VolCartesian(dimensions, meshMin, delta, nc));

    return mesh;
}

/*!
* Subtest 001
*
* Testing evaluation of levelset on Cartesian mesh in light mode.
*/
int subtest_001()
{
    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_006");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolCartesian> mesh = generateMesh(*segmentation);
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    // Get test cell
    long testCellId = 3671;

    // Initialize levelset
    bitpit::LevelSet levelset;
    int objectId = 0;

    // Compute levleset with the patch in normal memory mode
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->initializeAdjacencies();
    mesh->update();

    levelset.clear();
    levelset.setMesh(mesh.get());
    levelset.addObject(segmentation.get(), BITPIT_PI, objectId);

    levelset.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_006_normal");
    mesh->write();

    double normalValue = levelset.getObject(objectId).evalCellValue(testCellId, true);

    // Compute levleset with the patch in light memory mode
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    levelset.clear();
    levelset.setMesh(mesh.get());
    levelset.addObject(segmentation.get(), BITPIT_PI, objectId);

    double lightValue = levelset.getObject(objectId).evalCellValue(testCellId, true);

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    levelset.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_006_light");
    mesh->write();

    // Check if the light and normal values match
    bitpit::log::cout() << " Checking levelset values" << std::endl;

    if (!bitpit::utils::DoubleFloatingEqual()(normalValue, lightValue)) {
        bitpit::log::cout() << "  - Value obtained in light mode doesn't match the one obtained in normal mode." << std::endl;
        return 1;
    }
    bitpit::log::cout() << "  - Value obtained in light mode matches the one obtained in normal mode." << std::endl;

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
    bitpit::log::cout() << "Testing Cartesian mesh in light mode" << std::endl;

    int status;
    try {
        status = subtest_001();
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
