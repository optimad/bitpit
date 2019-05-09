
/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
# include "bitpit_CG.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"


/*!
* Subtest 001
*
* Testing levelset refinement.
*/
int subtest_001()
{
    uint8_t dimensions(2);

    // First Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (0,1,dimensions) );

    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;

    STL0->importDGF("./data/naca0012.dgf");

    STL0->deleteCoincidentVertices() ;
    STL0->buildAdjacencies() ;

    STL0->getVTK().setName("geometry_003_0") ;
    STL0->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL0->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL0->getCellCount() << std::endl;


    // Second Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL1( new bitpit::SurfUnstructured (1,dimensions) );

    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;

    STL1->importDGF("./data/square.dgf");

    STL1->deleteCoincidentVertices() ;
    STL1->buildAdjacencies() ;

    STL1->getVTK().setName("geometry_003_1") ;
    STL1->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL1->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL1->getCellCount() << std::endl;

    // Third Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL2( new bitpit::SurfUnstructured (1,dimensions) );

    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;

    STL2->importDGF("./data/rectangle.dgf");

    STL2->deleteCoincidentVertices() ;
    STL2->buildAdjacencies() ;

    STL2->getVTK().setName("geometry_003_2") ;
    STL2->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL2->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL2->getCellCount() << std::endl;

    // Create initial octree mesh
    bitpit::log::cout()<< " - Setting mesh" << std::endl;
    std::array<double,3>    mesh0, mesh1;
    std::array<double,3>    meshMin, meshMax, delta ;
    double h(0), dh ;

    STL0->getBoundingBox( meshMin, meshMax ) ;

    STL1->getBoundingBox( mesh0, mesh1 ) ;
    bitpit::CGElem::unionAABB( meshMin, meshMax, mesh0, mesh1, meshMin, meshMax ) ;

    STL2->getBoundingBox( mesh0, mesh1 ) ;
    bitpit::CGElem::unionAABB( meshMin, meshMax, mesh0, mesh1, meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    for( int i=0; i<3; ++i){
        h = std::max( h, meshMax[i]-meshMin[i] ) ;
    };

    dh = h / 16. ;
    bitpit::VolOctree mesh(dimensions, meshMin, h, dh );
    mesh.buildAdjacencies();
    mesh.buildInterfaces();
    mesh.update() ;

    // Set levelset configuration
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int elapsed_init, elapsed_refi(0);

    bitpit::LevelSet levelset;
    
    std::vector<bitpit::adaption::Info> mapper ;
    int id0, id1, id2, id3, id4, id5;

    levelset.setMesh(&mesh) ;
    id0 = levelset.addObject(std::move(STL0),BITPIT_PI) ;
    id1 = levelset.addObject(std::move(STL1),BITPIT_PI) ;
    id2 = levelset.addObject(std::move(STL2),BITPIT_PI/10.) ;

    id3 = levelset.addObject(bitpit::LevelSetBooleanOperation::UNION,id0,id1) ;
    id4 = levelset.addObject(bitpit::LevelSetBooleanOperation::SUBTRACTION,id3,id2) ;

    std::vector<int> ids;
    ids.push_back(id0);
    ids.push_back(id1);
    ids.push_back(id2);
    id5 = levelset.addObject(bitpit::LevelSetBooleanOperation::UNION,ids) ;

    ids.push_back(id3);
    ids.push_back(id4);
    ids.push_back(id5);

    levelset.getObject(id0).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id1).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id2).enableVTKOutput(bitpit::LevelSetWriteField::DEFAULT);
    levelset.getObject(id3).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    levelset.getObject(id4).enableVTKOutput(bitpit::LevelSetWriteField::DEFAULT);
    levelset.getObject(id5).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);

    mesh.getVTK().setName("levelset_003") ;
    mesh.getVTK().setCounter() ;

    levelset.setPropagateSign(true);

    // Compute and write level set on initial mesh
    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    bitpit::log::cout() << " - Exporting data" << std::endl;
    mesh.write() ;

    // Refine mesh, update levelset and write data
    const bitpit::LevelSetObject &object0 = levelset.getObject(id0);
    const bitpit::LevelSetObject &object1 = levelset.getObject(id1);
    const bitpit::LevelSetObject &object2 = levelset.getObject(id2);

    for( int i=0; i<10; ++i){

        for( auto & cell : mesh.getCells() ){
            const long &cellId = cell.getId() ;
            if( std::abs(object0.getLS(cellId)) < mesh.evalCellSize(cellId)  ){
                mesh.markCellForRefinement(cellId) ;
            }

            if( i<3) {
                if( std::abs(object1.getLS(cellId)) < mesh.evalCellSize(cellId)  ){
                    mesh.markCellForRefinement(cellId) ;
                }
            }

            if( i<6) {
                if( std::abs(object2.getLS(cellId)) < mesh.evalCellSize(cellId)  ){
                    mesh.markCellForRefinement(cellId) ;
                }
            }
        }

        mapper = mesh.update(true) ;
        start = std::chrono::system_clock::now();
        levelset.update(mapper) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

        mesh.write() ;
    }

    bitpit::log::cout() << "elapsed time initialization " << elapsed_init << " ms" << std::endl;
    bitpit::log::cout() << "elapsed time refinement     " << elapsed_refi << " ms" << std::endl;

    return 0;
};

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
	bitpit::log::cout() << "Testing levelset refinement" << std::endl;

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
