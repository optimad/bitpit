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

//Standard Template Library
# include <ctime>
# include <chrono>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"

/*!
 * Generate segmentation.
 *
 * \result The generated segmentation.
 */
std::unique_ptr<bitpit::SurfUnstructured> generateSegmentation()
{
    // Input geometry
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> segmentation( new bitpit::SurfUnstructured(2, MPI_COMM_NULL) );
#else
    std::unique_ptr<bitpit::SurfUnstructured> segmentation( new bitpit::SurfUnstructured(2) );
#endif

    segmentation->importSTL("./data/cube.stl", true);

    segmentation->deleteCoincidentVertices() ;
    segmentation->initializeAdjacencies() ;

    segmentation->getVTK().setName("geometry_002") ;
    segmentation->write() ;

    return segmentation;
}

/*!
 * Generate the Cartesian mesh.
 *
 * \result The generated Cartesian mesh.
 */
std::unique_ptr<bitpit::VolCartesian> generateCartesianMesh(const bitpit::SurfUnstructured &segmentation)
{
    int dimensions = 3;

    std::array<double, 3> meshMin, meshMax;
    segmentation.getBoundingBox(meshMin, meshMax);

    std::array<double, 3> delta = meshMax - meshMin;
    meshMin -= 0.1 * delta;
    meshMax += 0.1 * delta;
    delta = meshMax - meshMin;

    std::array<int, 3> nc = {{64, 64, 64}};

    std::unique_ptr<bitpit::VolCartesian> mesh(new bitpit::VolCartesian(dimensions, meshMin, delta, nc));

    return mesh;
}

/*!
 * Generate the Octree mesh.
 *
 * \result The generated Octree mesh.
 */
std::unique_ptr<bitpit::VolOctree> generateOctreeMesh(const bitpit::SurfUnstructured &segmentation)
{
    int dimensions = 3;

    std::array<double, 3> segmentationMin;
    std::array<double, 3> segmentationMax;
    segmentation.getBoundingBox(segmentationMin, segmentationMax);

    std::array<double, 3> delta = segmentationMax - segmentationMin;
    segmentationMin -= 0.1 * delta;
    segmentationMax += 0.1 * delta;
    delta = segmentationMax - segmentationMin;

    std::array<double, 3> origin = segmentationMin;

    double length = 0.;
    for (int i = 0; i < 3; ++i) {
        length = std::max(length, segmentationMax[i] - segmentationMin[i]);
    };

    double dh = length / 64;

#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::VolOctree> mesh(new bitpit::VolOctree(dimensions, origin, length, dh, MPI_COMM_NULL));
#else
    std::unique_ptr<bitpit::VolOctree> mesh(new bitpit::VolOctree(dimensions, origin, length, dh));
#endif

    return mesh;
}

/*!
* Subtest 001
*
* Testing basic features of a 3D levelset on a Cartesian mesh in default memory mode.
*/
int subtest_001()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing three-dimensional levelset on a Cartesian mesh in default memory mode" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_002");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolCartesian> mesh = generateCartesianMesh(*segmentation);
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->initializeAdjacencies();
    mesh->update();

    // Initialize levelset
    bitpit::log::cout() << " - Initializing levelset" << std::endl;

    int objectId = 0;

    bitpit::LevelSet levelset ;
    levelset.setPropagateSign(true);
    levelset.setMesh(mesh.get());
    levelset.addObject(segmentation.get(), BITPIT_PI, objectId);

    // Compute levelset
    bitpit::log::cout() << " - Evaluating levelset" << std::endl;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    levelset.compute( ) ;

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    // Write output
    bitpit::log::cout() << " - Writing output" << std::endl;

    levelset.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_002_cartesian_normal");
    mesh->write();

    return 0;
}

/*!
* Subtest 002
*
* Testing basic features of a 3D levelset on a Cartesian mesh in light memory mode.
*/
int subtest_002()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing three-dimensional levelset on a Cartesian mesh in light memory mode" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_002");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolCartesian> mesh = generateCartesianMesh(*segmentation);
    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_LIGHT);

    // Initialize levelset
    bitpit::log::cout() << " - Initializing levelset" << std::endl;

    int objectId = 0;

    bitpit::LevelSet levelset ;
    levelset.setPropagateSign(false);
    levelset.setMesh(mesh.get());
    levelset.addObject(segmentation.get(), BITPIT_PI, objectId);

    // Compute levelset
    bitpit::log::cout() << " - Evaluating levelset" << std::endl;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    levelset.compute( ) ;

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    // Write output
    bitpit::log::cout() << " - Writing output" << std::endl;

    mesh->switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh->getVTK().setName("levelset_002_cartesian_light");
    mesh->write();

    return 0;
}

/*!
* Subtest 003
*
* Testing basic features of a 3D levelset on an Octreee mesh.
*/
int subtest_003()
{
    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing three-dimensional levelset on an Octree mesh" << std::endl;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation();
    segmentation->getVTK().setName("geometry_002");
    segmentation->write();

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Create the mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;

    std::unique_ptr<bitpit::VolOctree> mesh = generateOctreeMesh(*segmentation);
    mesh->initializeAdjacencies();
    mesh->update();

    // Initialize levelset
    bitpit::log::cout() << " - Initializing levelset" << std::endl;

    int objectId = 0;

    bitpit::LevelSet levelset ;
    levelset.setPropagateSign(true);
    levelset.setMesh(mesh.get());
    levelset.addObject(segmentation.get(), BITPIT_PI, objectId);

    // Compute levelset
    bitpit::log::cout() << " - Evaluating levelset" << std::endl;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    levelset.compute( ) ;

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    // Write output
    bitpit::log::cout() << " - Writing output" << std::endl;

    levelset.getObject(objectId).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh->getVTK().setName("levelset_002_octree");
    mesh->write();

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
	bitpit::log::manager().initialize(bitpit::log::MODE_COMBINE);

	// Run the subtests
	bitpit::log::cout() << "Testing basic levelset features" << std::endl;

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
