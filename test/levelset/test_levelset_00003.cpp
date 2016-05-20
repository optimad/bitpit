
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
# include "bitpit_levelset.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

/*!Demo for 2D level set of complex geometries on a Pablo octree mesh.
*/
int main( int argc, char *argv[]){

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc, &argv);
#endif
    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    uint8_t                 dimensions(2);


    // Input geometry
    bitpit::SurfUnstructured    STL(0,1,dimensions);

    std::cout << " - Loading stl geometry" << std::endl;

    STL.importDGF("./data/naca0012.dgf");

    STL.deleteCoincidentVertices() ;
    STL.buildAdjacencies() ;

    STL.getVTK().setName("geometry_003") ;
    STL.write() ;

    std::cout << "n. vertex: " << STL.getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL.getCellCount() << std::endl;



    // ========================================================================== //
    // CREATE MESH                                                                //
    // ========================================================================== //
    std::cout << " - Setting mesh" << std::endl;
    std::array<double,3>    meshMin, meshMax, delta ;
    double                  h(0), dh ;

    STL.getBoundingBox( meshMin, meshMax ) ;

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
    std::vector<double>             LS ;
    std::vector<double>::iterator   itLS ;

    levelset.setMesh(&mesh) ;
    levelset.addObject(&STL) ;


    mesh.getVTK().addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS) ;
    mesh.getVTK().setName("levelset_003") ;
    mesh.getVTK().setCounter() ;

    levelset.setPropagateSign(true);

    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();

    // Export level set ------------------------------------------------------- //
    std::cout << " - Exporting data" << endl;

    LS.resize(mesh.getCellCount() ) ;
    itLS = LS.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        *itLS = levelset.getLS(id) ;
        ++itLS ;
    };

    mesh.write() ;

    //Refinement
    for( int i=0; i<10; ++i){

        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            if( std::abs(levelset.getLS(id)) < 100. ){
                mesh.markCellForRefinement(id) ;
            }
        }

        mapper = mesh.update(true) ;
        start = std::chrono::system_clock::now();
        levelset.update(mapper) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += chrono::duration_cast<chrono::milliseconds>(end-start).count();

        LS.resize(mesh.getCellCount() ) ;
        itLS = LS.begin() ;
        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            *itLS = levelset.getLS(id) ;
            ++itLS ;
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


