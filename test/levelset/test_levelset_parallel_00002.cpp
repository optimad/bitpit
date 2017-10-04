
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

/*!
 *	\date			10/jul/2014
 *	\authors		Alessandro Alaia
 *	\authors		Haysam Telib
 *	\authors		Edoardo Lombardi
 *	\version		0.1
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Level Set Class Demos
 */

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

//Standard Template Library
# include <ctime>
# include <chrono>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_CG.hpp"
# include "bitpit_IO.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

/*!
* Subtest 001
*
* Testing 2D levelset parallel refinement.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    uint8_t                 dimensions(2);

    // First Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (1,dimensions) );

    std::cout << " - Loading stl geometry" << std::endl;

    STL0->importDGF("./data/naca0012.dgf");

    STL0->deleteCoincidentVertices() ;
    STL0->buildAdjacencies() ;

    STL0->getVTK().setName("geometry_parallel_002_0") ;
    if(rank==0){
        STL0->write() ;
    }

    std::cout << "n. vertex: " << STL0->getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL0->getCellCount() << std::endl;


    // Second Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL1( new bitpit::SurfUnstructured (1,dimensions) );

    std::cout << " - Loading stl geometry" << std::endl;

    STL1->importDGF("./data/square.dgf");

    STL1->deleteCoincidentVertices() ;
    STL1->buildAdjacencies() ;

    STL1->getVTK().setName("geometry_parallel_002_1") ;
    if(rank==0){
        STL1->write() ;
    }

    std::cout << "n. vertex: " << STL1->getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL1->getCellCount() << std::endl;

    // ========================================================================== //
    // CREATE MESH                                                                //
    // ========================================================================== //
    std::cout << " - Setting mesh" << std::endl;
    std::array<double,3>    meshMin0, meshMax0;
    std::array<double,3>    meshMin1, meshMax1;
    std::array<double,3>    meshMin, meshMax, delta ;
    double                  h(0), dh ;

    STL0->getBoundingBox( meshMin0, meshMax0 ) ;
    STL1->getBoundingBox( meshMin1, meshMax1 ) ;
    bitpit::CGElem::unionAABB( meshMin0, meshMax0, meshMin1, meshMax1, meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    for( int i=0; i<3; ++i){
        h = max( h, meshMax[i]-meshMin[i] ) ;
    };

    dh = h / 16. ;
    bitpit::VolOctree    mesh(dimensions, meshMin, h, dh );
    mesh.update() ;


    // COMPUTE LEVEL SET in NARROW BAND
    std::chrono::time_point<std::chrono::system_clock>    start, end;
    int                                         elapsed_init, elapsed_refi(0);

    bitpit::LevelSet                levelset;

    std::vector<bitpit::adaption::Info> mapper ;
    int                             id0, id1, id2 ;
    std::vector<double>             LS0, LS1, LS2 ;
    std::vector<double>::iterator   it0, it1, it2 ;

    levelset.setMesh(&mesh) ;
    id0 = levelset.addObject(std::move(STL0),M_PI) ;
    id1 = levelset.addObject(std::move(STL1),M_PI) ;
    id2 = levelset.addObject(bitpit::LevelSetBooleanOperation::INTERSECTION,id0,id1) ;
    const bitpit::LevelSetObject &object0 = levelset.getObject(id0) ;
    const bitpit::LevelSetObject &object1 = levelset.getObject(id1) ;
    const bitpit::LevelSetObject &object2 = levelset.getObject(id2) ;

    mesh.getVTK().addData("ls0", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS0) ;
    mesh.getVTK().addData("ls1", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS1) ;
    mesh.getVTK().addData("ls2", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS2) ;
    mesh.getVTK().setName("levelset_parallel_002") ;

    levelset.setPropagateSign(true);

    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();


    // Export level set ------------------------------------------------------- //
    std::cout << " - Exporting data" << endl;
    LS0.resize(mesh.getCellCount() ) ;
    LS1.resize(mesh.getCellCount() ) ;
    LS2.resize(mesh.getCellCount() ) ;
    it0 = LS0.begin() ;
    it1 = LS1.begin() ;
    it2 = LS2.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &cellId = cell.getId();
        *it0 = object0.getLS(cellId);
        *it1 = object1.getLS(cellId);
        *it2 = object2.getLS(cellId);
        ++it0;
        ++it1;
        ++it2;
    };

    if(rank==0){
        mesh.write() ;
    }

    //Mesh Partitioning
    mapper = mesh.partition(MPI_COMM_WORLD, true) ;
    mesh.getVTK().setCounter() ;
    levelset.partition(mapper) ;

    //Refinement
    for( int i=0; i<10; ++i){

        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            if( std::abs(object0.getLS(id)) < mesh.evalCellSize(id) ){
                mesh.markCellForRefinement(id) ;
            }

            if( i<4){
                if( std::abs(object1.getLS(id)) < mesh.evalCellSize(id) ){
                    mesh.markCellForRefinement(id) ;
                }
            }
        }

        mapper = mesh.update(true) ;
        start = std::chrono::system_clock::now();
        levelset.update(mapper) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += chrono::duration_cast<chrono::milliseconds>(end-start).count();

        std::cout << "Refinement " << i << " Narrow Band id " << id0 << " " << object0.getSizeNarrowBand() << endl;
        std::cout << "Refinement " << i << " Narrow Band id " << id1 << " " << object1.getSizeNarrowBand() << endl;

        LS0.resize(mesh.getCellCount());
        LS1.resize(mesh.getCellCount());
        LS2.resize(mesh.getCellCount());
        it0 = LS0.begin() ;
        it1 = LS1.begin() ;
        it2 = LS2.begin() ;
        for( auto & cell : mesh.getCells() ){
            const long &cellId = cell.getId();
            *it0 = object0.getLS(cellId);
            *it1 = object1.getLS(cellId);
            *it2 = object2.getLS(cellId);
            ++it0 ;
            ++it1 ;
            ++it2 ;
        };
        mesh.write() ;
    }

    cout << "elapsed time initialization " << elapsed_init << " ms" << endl;
    cout << "elapsed time refinement     " << elapsed_refi << " ms" << endl;

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
    int    rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bitpit::log::manager().initialize(bitpit::log::COMBINED, true, nProcs, rank);
    bitpit::log::cout().setVisibility(bitpit::log::GLOBAL);

    // Run the subtests
    bitpit::log::cout() << "Testing levelset parallel refinement." << std::endl;

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

    MPI_Finalize();
}
