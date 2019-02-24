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

//Standard Template Library
# include <ctime>
# include <chrono>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_levelset.hpp"

/*!
* Subtest 001
*
* Testing basic features of a 3D levelset.
*/
int subtest_001()
{
    int dimensions(3) ;

    // Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(2, 3) );

    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;

    STL->importSTL("./data/cube.stl", true);

    STL->deleteCoincidentVertices() ;
    STL->buildAdjacencies() ;

    STL->getVTK().setName("geometry_002") ;
    STL->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL->getCellCount() << std::endl;


    // create cartesian mesh around geometry 
    bitpit::log::cout()<< " - Setting mesh" << std::endl;
    std::array<double,3> meshMin, meshMax, delta ;
    std::array<int,3> nc = {{64, 64, 64}} ;

    STL->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);
    mesh.update() ;
    mesh.buildAdjacencies() ;
    mesh.buildInterfaces() ;

    // Compute level set  in narrow band
    bitpit::LevelSet levelset ;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    levelset.setMesh(&mesh) ;
    int id0 = levelset.addObject( std::move(STL), BITPIT_PI/3. ) ;

    levelset.setPropagateSign(true) ;
    start = std::chrono::system_clock::now();
    levelset.compute( ) ;
    end = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout()<< "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    bitpit::log::cout()<< " - Exporting data" << std::endl;
    mesh.getVTK().setName("levelset_002") ;
    bitpit::LevelSetObject &object = levelset.getObject(id0);

    object.enableVTKOutput( bitpit::LevelSetWriteField::ALL);

    mesh.write() ;

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
	bitpit::log::cout() << "Testing basic levelset features" << std::endl;

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
