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

/*!
* Generate a 2D surface mesh.
*
* \param[in,out] mesh on output it will contain the generated surface mesh
*/
void Generate2DSurfMesh( bitpit::SurfUnstructured &mesh )
{

    const double R = 1.0;
    const long N = 32;
    double theta;
    double dtheta = 2. * BITPIT_PI/((double) N);
    
    std::array<double,3> point;
    std::vector<long> connect(2, bitpit::Element::NULL_ID);

    // Create vertex list
    // Use non-consecutive vertex ids to test if the levelset can handle them.
    const long vertexIdOffset = 101;
    const long vertexIdStride = 2;

    point[2] = 0.0;
    for (long i = 0; i < N; ++i) {
        theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        mesh.addVertex(point, vertexIdOffset + vertexIdStride * i);
    }

    // Create simplex list
    // Use non-consecutive cell ids to test if the levelset can handle them.
    const long cellIdOffset = 202;
    const long cellIdStride = 3;

    for (long i = 0; i < N; ++i) {
        connect[0] = vertexIdOffset + vertexIdStride * i;
        connect[1] = vertexIdOffset + vertexIdStride * ((i + 1) % N);
        mesh.addCell(bitpit::ElementType::LINE, connect, cellIdOffset + cellIdStride * i);
    }

}

/*!
* Subtest 001
*
* Testing basic features of a 2D levelset.
*/
int subtest_001()
{
    int dimensions(2) ;

    // Input geometry
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(dimensions - 1,MPI_COMM_NULL) );
#else
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(dimensions - 1) );
#endif
    STL->initializeAdjacencies();

    bitpit::log::cout() << " - Loading dgf geometry" << std::endl;

    Generate2DSurfMesh( *(STL.get()) ) ;

    STL->getVTK().setName("geometry_001") ;
    STL->write() ;

    bitpit::log::cout() << "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL->getCellCount() << std::endl;

    // create cartesian mesh around geometry 
    bitpit::log::cout() << " - Setting mesh" << std::endl;
    std::array<double,3> meshMin, meshMax, delta ;
    std::array<int,3> nc = {{64, 64, 0}} ;

    STL->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);
    mesh.switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL) ;
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
    std::chrono::milliseconds elapsed_seconds;
    start = std::chrono::system_clock::now();

    bitpit::LevelSet levelset ;
    levelset.setNarrowBandSize(0) ;
    levelset.setMesh(&mesh) ;

    int id0 = levelset.addObject(std::move(STL),BITPIT_PI) ;
    int id1 = levelset.addObject(mask) ;

    bitpit::LevelSetObject *object0 = static_cast<bitpit::LevelSetObject *>(levelset.getObjectPtr(id0));
    bitpit::LevelSetObject *object1 = static_cast<bitpit::LevelSetObject *>(levelset.getObjectPtr(id1));

    object0->setCellBulkEvaluationMode(bitpit::LevelSetBulkEvaluationMode::NONE);
    object1->setCellBulkEvaluationMode(bitpit::LevelSetBulkEvaluationMode::SIGN_PROPAGATION);

    object0->enableFieldCellCache(bitpit::LevelSetField::VALUE, bitpit::LevelSetCacheMode::NARROW_BAND);
    object1->enableFieldCellCache(bitpit::LevelSetField::VALUE, bitpit::LevelSetCacheMode::FULL);

    end = std::chrono::system_clock::now();
    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds.count() << " ms" << std::endl;

    bitpit::log::cout() << " - Writing output" << std::endl;

    levelset.getObject(id0).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id1).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    mesh.getVTK().setName("levelset_001") ;
    mesh.write() ;

    return 0;
}


/*!
* Subtest 002
*
* Testing projection and interface intersection on surface feature of a 2D levelset on
* a Cartesian mesh in default memory mode.
*/
int subtest_002()
{
    int dimensions(2) ;

    // Input geometry
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(dimensions - 1,MPI_COMM_NULL) );
#else
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(dimensions - 1) );
#endif

    bitpit::log::cout() << " - Loading dgf geometry" << std::endl;

    STL->importDGF("./data/naca0012_coarse.dgf", true);

    STL->initializeAdjacencies();

    bitpit::log::cout() << "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL->getCellCount() << std::endl;

    // create cartesian mesh around geometry 
    bitpit::log::cout() << " - Setting mesh" << std::endl;
    std::array<double,3> meshMin, meshMax, delta ;
    std::array<int,3> nc = {{64, 64, 0}} ;

    STL->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);
    mesh.switchMemoryMode(bitpit::VolCartesian::MEMORY_NORMAL);
    mesh.initializeAdjacencies() ;
    mesh.initializeInterfaces() ;
    mesh.update() ;

    // Compute level set in narrow band
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int elapsed_seconds;
    start = std::chrono::system_clock::now();

    bitpit::LevelSet levelset ;
    levelset.setNarrowBandSize(0) ;
    levelset.setMesh(&mesh) ;

    int objectId = 0;
    double limitAngle = 85.0 * BITPIT_PI / 180.0;
    levelset.addObject(STL.get(), limitAngle, bitpit::LevelSetSurfaceSmoothing::HIGH_ORDER, objectId);

    bitpit::LevelSetObject *object = static_cast<bitpit::LevelSetObject *>(levelset.getObjectPtr(objectId));

    object->setCellBulkEvaluationMode(bitpit::LevelSetBulkEvaluationMode::SIGN_PROPAGATION);
    object->enableFieldCellCache(bitpit::LevelSetField::VALUE, bitpit::LevelSetCacheMode::FULL);

    // Compute projections on curved surface
    bitpit::LevelSetSegmentationBaseObject *levelSetSegmentation = dynamic_cast<bitpit::LevelSetSegmentationBaseObject *>(object);

    int cellId = 640;
    std::array<double, 3> point = mesh.evalCellCentroid(cellId);
    std::array<double, 3> projectionPoint;
    std::array<double, 3> projectionNormal;

    levelSetSegmentation->evalProjection(point, true, &projectionPoint, &projectionNormal);

    std::array<double, 3> projectionPointTarget  = {0.0021934049109738, -0.008062212041176, 0.0};
    std::array<double, 3> projectionNormalTarget = {-0.87574759289519, -0.48276925496379, 0.0};

    double pointDeviation  = norm2(projectionPoint - projectionPointTarget);
    double normalDeviation = norm2(projectionNormal - projectionNormalTarget);

    double distanceTolerance = mesh.getTol();
    bool projectionFound = (bitpit::utils::DoubleFloatingEqual()(pointDeviation, 0., distanceTolerance, distanceTolerance)
                         && bitpit::utils::DoubleFloatingEqual()(normalDeviation, 0., distanceTolerance, distanceTolerance));

    // Compute intersected interface
    int intrId = 744;

    std::array<std::array<double, 3>, 2> intersection;
    std::vector<std::array<double, 3>> polygon;
    bool intersectionFound = (levelSetSegmentation->isInterfaceIntersected(intrId, false, &intersection, &polygon) == bitpit::LevelSetIntersectionStatus::TRUE);

    std::array<double, 3> intersectionMean = 0.5 * (intersection[0] + intersection[1]);

    int size = polygon.size();
    std::unique_ptr<long[]> connectivity(new long[size + 1]);
    connectivity[0] = size;
    for (int i = 0; i < size; ++i) {
        connectivity[i + 1] = i;
    }

    bitpit::Element element(0, bitpit::ElementType::POLYGON, std::move(connectivity));
    std::array<double, 3> polygonMean = element.evalCentroid(polygon.data());

    std::array<double, 3> intersectionMeanTarget = {0.234654275, -0.059034230965178, 0.0};
    std::array<double, 3> polygonMeanTarget      = {0.234654275, -0.059897734876829, 0.0};

    double intersectionDeviation = norm2(intersectionMean - intersectionMeanTarget);
    double polygonDeviation      = norm2(polygonMean - polygonMeanTarget);

    intersectionFound = intersectionFound && (bitpit::utils::DoubleFloatingEqual()(intersectionDeviation, 0., distanceTolerance, distanceTolerance)
                       && bitpit::utils::DoubleFloatingEqual()(polygonDeviation, 0., distanceTolerance, distanceTolerance));

    end = std::chrono::system_clock::now();
    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    return !(projectionFound && intersectionFound);
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
