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


# include <mpi.h>

//Standard Template Library
# include <ctime>
# include <chrono>

// bitpit
# include "bitpit_IO.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"

/*!
 * Generate segmentation.
 *
 * \param name is the name of the segmentation
 * \result The generated segmentation.
 */
std::unique_ptr<bitpit::SurfUnstructured> generateSegmentation(const std::string &name)
{
    // Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> segmentation( new bitpit::SurfUnstructured(2, MPI_COMM_NULL) );

    segmentation->importSTL("./data/cube.stl", true);

    segmentation->deleteCoincidentVertices();
    segmentation->initializeAdjacencies();

    segmentation->getVTK().setName(name);

    return segmentation;
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

    double dh = length / 16;

    std::unique_ptr<bitpit::VolOctree> mesh(new bitpit::VolOctree(dimensions, origin, length, dh, MPI_COMM_WORLD));

    return mesh;
}

/*!
* Subtest 001
*
* Testing basic features of a 3D levelset on an Octreee mesh.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    BITPIT_UNUSED(rank);

    bitpit::log::cout() << std::endl;
    bitpit::log::cout() << "Testing three-dimensional levelset on an Octree mesh" << std::endl;

    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    // Input geometry
    bitpit::log::cout() << " - Loading geometry" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation("geometry_002");
    if (rank == 0) {
        segmentation->write();
    }

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

    // Compute levelset in serial
    bitpit::log::cout() << " - Evaluating the levelset" << std::endl;

    start = std::chrono::system_clock::now();
    levelset.compute();
    end = std::chrono::system_clock::now();

    int elapsed_init = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    bitpit::log::cout() << " - Exporting serial levelset" << std::endl;
    mesh->getVTK().setName("levelset_parallel_001_octree_serial") ;
    mesh->write() ;

    // Partition the patch
    bitpit::log::cout() << " - Partitioning the patch" << std::endl;

    std::vector<bitpit::adaption::Info> partitioningData = mesh->partition(true) ;

    // Compute levelset in parallel
    start = std::chrono::system_clock::now();
    levelset.update(partitioningData) ;
    end = std::chrono::system_clock::now();

    int elapsed_part = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    bitpit::log::cout() << " - Exporting partitioned levelset" << std::endl;
    mesh->getVTK().setName("levelset_parallel_001_octree_partitioned") ;
    mesh->write() ;

    // Refine mesh and update levelset
    const bitpit::LevelSetObject &object0 = levelset.getObject(objectId);

    mesh->getVTK().setName("levelset_parallel_001_octree_refined") ;
    mesh->getVTK().setCounter() ;

    int elapsed_refi = 0;
    for (int i=0; i<3; ++i) {
        for (const bitpit::Cell &cell : mesh->getCells()) {
            long id = cell.getId() ;
            if (std::abs(object0.getValue(id)) < 100.) {
                mesh->markCellForRefinement(id) ;
            }
        }

        std::vector<bitpit::adaption::Info> adaptionData = mesh->update(true) ;

        start = std::chrono::system_clock::now();
        levelset.update(adaptionData) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

        bitpit::log::cout() << " - Exporting serial levelset" << std::endl;
        mesh->write();
    }

    // Write elapsed times
    bitpit::log::cout() << "elapsed time initialization " << elapsed_init << " ms" << std::endl;
    bitpit::log::cout() << "elapsed time partitioning   " << elapsed_part << " ms" << std::endl;
    bitpit::log::cout() << "elapsed time refinement     " << elapsed_refi << " ms" << std::endl;

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

    bitpit::log::manager().initialize(bitpit::log::MODE_COMBINE, true, nProcs, rank);
    bitpit::log::cout().setDefaultVisibility(bitpit::log::VISIBILITY_GLOBAL);

    // Run the subtests
    bitpit::log::cout() << "Testing basic levelset features" << std::endl;

    int status;
    try {
        status = subtest_001(rank);
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
