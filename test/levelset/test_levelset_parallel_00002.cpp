
/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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
# include "bitpit_levelset.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

/*!Demo for 2D level set of complex geometries on a Pablo octree mesh.
*/
int main( int argc, char *argv[]){

    uint8_t                 dimensions(2);

    // MPI information
    int nProcessors, rank ;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initialize logger
    bitpit::log::manager().initialize(bitpit::log::SEPARATE, false, nProcessors, rank) ;
    bitpit::log::cout().setVisibility( bitpit::log::GLOBAL ) ;

    // First Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (0,1,dimensions) );

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
    std::unique_ptr<bitpit::SurfUnstructured> STL1( new bitpit::SurfUnstructured (0,1,dimensions) );

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
    bitpit::VolOctree    mesh(1, dimensions, meshMin, h, dh );
    mesh.update() ;


    // COMPUTE LEVEL SET in NARROW BAND
    std::chrono::time_point<std::chrono::system_clock>    start, end;
    int                                         elapsed_init, elapsed_refi(0);

    bitpit::LevelSet                levelset;

    std::vector<bitpit::adaption::Info> mapper ;
    int                             id0, id1 ;
    std::vector<double>             LS, LS0, LS1 ;
    std::vector<double>::iterator   itLS, itLS0, itLS1 ;

    levelset.setMesh(&mesh) ;
    id0 = levelset.addObject(std::move(STL0),M_PI) ;
    id1 = levelset.addObject(std::move(STL1),M_PI) ;

    mesh.getVTK().addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS) ;
    mesh.getVTK().addData("ls0", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS0) ;
    mesh.getVTK().addData("ls1", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS1) ;
    mesh.getVTK().setName("levelset_parallel_002") ;

    levelset.setPropagateSign(true);

    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();


    // Export level set ------------------------------------------------------- //
    std::cout << " Narrow Band " << levelset.getSizeNarrowBand() << endl;
    std::cout << " - Exporting data" << endl;

    LS.resize(mesh.getCellCount() ) ;
    LS0.resize(mesh.getCellCount() ) ;
    LS1.resize(mesh.getCellCount() ) ;
    itLS = LS.begin() ;
    itLS0 = LS0.begin() ;
    itLS1 = LS1.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &cellId = cell.getId() ;
        *itLS = levelset.getLS(cellId) ;
        *itLS0 = levelset.getLS(id0,cellId) ;
        *itLS1 = levelset.getLS(id1,cellId) ;
        ++itLS ;
        ++itLS0 ;
        ++itLS1 ;
    };

    if(rank==0){
        mesh.write() ;
    }

    //Mesh Partitioning
    mapper = mesh.partition(MPI_COMM_WORLD, true) ;
    mesh.getVTK().setCounter() ;
    levelset.update(mapper) ;

    //Refinement
    for( int i=0; i<10; ++i){

        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            if( std::abs(levelset.getLS(id)) < mesh.evalCellSize(id) ){
                mesh.markCellForRefinement(id) ;
            }
        }

        mapper = mesh.update(true) ;
        start = std::chrono::system_clock::now();
        levelset.update(mapper) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += chrono::duration_cast<chrono::milliseconds>(end-start).count();

        std::cout << " Narrow Band " << i << " " << levelset.getSizeNarrowBand() << endl;

        LS.resize(mesh.getCellCount() ) ;
        LS0.resize(mesh.getCellCount() ) ;
        LS1.resize(mesh.getCellCount() ) ;
        itLS = LS.begin() ;
        itLS0 = LS0.begin() ;
        itLS1 = LS1.begin() ;
        for( auto & cell : mesh.getCells() ){
            const long &cellId = cell.getId() ;
            *itLS = levelset.getLS(cellId) ;
            *itLS0 = levelset.getLS(id0,cellId) ;
            *itLS1 = levelset.getLS(id1,cellId) ;
            ++itLS ;
            ++itLS0 ;
            ++itLS1 ;
        };
        mesh.write() ;
    }

    cout << "elapsed time initialization " << elapsed_init << " ms" << endl;
    cout << "elapsed time refinement     " << elapsed_refi << " ms" << endl;

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return 0;

};


