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
# define _USE_MATH_DEFINES

# include <cmath>
# include <ctime>
# include <chrono>

# include <array>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_levelset.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"

/*!
* Generate a 2D surface mesh.
*
* \param dimension is the dimension of the patch
*/
std::unique_ptr<bitpit::SurfUnstructured> generateSegmentation(int dimension)
{
    // Surface parameters
    const double R = 1.0;
    const long   N = 32;

    // Creat an empty patch
#if BITPIT_ENABLE_MPI
    auto segmentation = std::unique_ptr<bitpit::SurfUnstructured>(new bitpit::SurfUnstructured(dimension, MPI_COMM_NULL));
#else
    auto segmentation = std::unique_ptr<bitpit::SurfUnstructured>(new bitpit::SurfUnstructured(dimension));
#endif

    // Create vertex list
    // Use non-consecutive vertex ids to test if the levelset can handle them.
    const long vertexIdOffset = 101;
    const long vertexIdStride = 2;

    double dtheta = 2. * BITPIT_PI/((double) N);

    std::array<double,3> point;
    point[2] = 0.0;
    for (long i = 0; i < N; ++i) {
        double theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        segmentation->addVertex(point, vertexIdOffset + vertexIdStride * i);
    }

    // Create simplex list
    //
    // Use non-consecutive cell ids to test if the levelset can handle them.
    std::vector<long> connect(2, bitpit::Element::NULL_ID);

    const long cellIdOffset = 202;
    const long cellIdStride = 3;

    for (long i = 0; i < N; ++i) {
        connect[0] = vertexIdOffset + vertexIdStride * i;
        connect[1] = vertexIdOffset + vertexIdStride * ((i + 1) % N);
        segmentation->addCell(bitpit::ElementType::LINE, connect, cellIdOffset + cellIdStride * i);
    }

    return segmentation;
}

/*!
* Subtest 001
*
* Testing basic features of a 2D levelset on a Certesian mesh.
*/
int subtest_001()
{
    int dimensions(2) ;

    // Initialize segmentation
    bitpit::log::cout() << " - Generating segmentation" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation(dimensions - 1);
    segmentation->initializeAdjacencies();
    segmentation->getVTK().setName("geometry_cartesian_001") ;
    segmentation->write() ;

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Initialize mesh
    bitpit::log::cout() << " - Generating mesh" << std::endl;
    std::array<double,3> meshMin, meshMax, delta ;
    std::array<int,3> nc = {{64, 64, 0}} ;

    segmentation->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);
    mesh.update() ;
    mesh.initializeAdjacencies() ;
    mesh.initializeInterfaces() ;

    // mark cells within R=0.5
    std::unordered_set<long> mask;
    for( auto & cell : mesh.getCells() ){
        long id = cell.getId() ;
        std::array<double,3> center = mesh.evalCellCentroid(id);
        double r = norm2(center);
        if(r<=0.5){
            mask.insert(id);
        }
    }

    // Compute level set in narrow band
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int elapsed_seconds;

    bitpit::LevelSet levelset ;

    levelset.setMesh(&mesh) ;

    int id0 = levelset.addObject(std::move(segmentation),BITPIT_PI) ;
    int id1 = levelset.addObject(mask) ;
    std::vector<int> ids;
    ids.push_back(id0);
    ids.push_back(id1);

    start = std::chrono::system_clock::now();
    levelset.compute( ) ;
    end = std::chrono::system_clock::now();

    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    bitpit::log::cout() << " - Exporting data" << std::endl;

    levelset.getObject(id0).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id1).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh.getVTK().setName("levelset_cartesian_001") ;
    mesh.write() ;

    bitpit::log::cout() << " - Exported data" << std::endl;

    return 0;
}

/*!
* Subtest 002
*
* Testing basic features of a 2D levelset on an Octree mesh.
*/
int subtest_002()
{
    int dimensions(2) ;

    // Initialize segmentation
    bitpit::log::cout() << " - Generating segmentation" << std::endl;

    std::unique_ptr<bitpit::SurfUnstructured> segmentation = generateSegmentation(dimensions - 1);
    segmentation->initializeAdjacencies();
    segmentation->getVTK().setName("geometry_cartesian_001") ;
    segmentation->write() ;

    bitpit::log::cout() << "n. vertex: " << segmentation->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << segmentation->getCellCount() << std::endl;

    // Initialize mesh
    bitpit::log::cout() << " - Generating mesh" << std::endl;

    std::array<double, 3> segmentationMin;
    std::array<double, 3> segmentationMax;
    segmentation->getBoundingBox(segmentationMin, segmentationMax);

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
    bitpit::VolOctree mesh(dimensions, origin, length, dh, MPI_COMM_NULL);
#else
    bitpit::VolOctree mesh(dimensions, origin, length, dh);
#endif

    mesh.update() ;
    mesh.initializeAdjacencies() ;
    mesh.initializeInterfaces() ;

    // mark cells within R=0.5
    std::unordered_set<long> mask;
    for( auto & cell : mesh.getCells() ){
        long id = cell.getId() ;
        std::array<double,3> center = mesh.evalCellCentroid(id);
        double r = norm2(center);
        if(r<=0.5){
            mask.insert(id);
        }
    }

    // Compute level set in narrow band
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int elapsed_seconds;

    bitpit::LevelSet levelset ;

    levelset.setMesh(&mesh) ;

    int id0 = levelset.addObject(std::move(segmentation),BITPIT_PI) ;
    int id1 = levelset.addObject(mask) ;
    std::vector<int> ids;
    ids.push_back(id0);
    ids.push_back(id1);

    start = std::chrono::system_clock::now();
    levelset.compute( ) ;
    end = std::chrono::system_clock::now();

    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    bitpit::log::cout() << " - Exporting data" << std::endl;

    levelset.getObject(id0).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id1).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh.getVTK().setName("levelset_octree_001") ;
    mesh.write() ;

    bitpit::log::cout() << " - Exported data" << std::endl;

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
	} catch (const std::exception &exception) {
		bitpit::log::cout() << exception.what();
		exit(1);
	}

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
}
