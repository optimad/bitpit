
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

# include <mpi.h>

//Standard Template Library
# include <ctime>
# include <chrono>

// bitpit
# include "bitpit_IO.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

/*!
    Test for 3D levelset of complex geometries on a Pablo octree mesh.
*/
int main( int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    // MPI information
    int nProcessors, rank ;

    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initialize logger
    bitpit::log::manager().initialize(bitpit::log::SEPARATE, false, nProcessors, rank) ;
    bitpit::log::cout().setVisibility( bitpit::log::GLOBAL ) ;

    // Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(0) );

    std::cout << " - Loading stl geometry" << std::endl;

    STL->importSTL("./data/cube.stl", true);

    STL->deleteCoincidentVertices() ;
    STL->buildAdjacencies() ;

    STL->getVTK().setName("geometry_002") ;
    if (rank == 0) {
        STL->write() ;
    }

    std::cout << "n. vertex: " << STL->getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL->getCellCount() << std::endl;

    // Create mesh
    std::cout << " - Setting mesh" << std::endl;
    std::array<double,3>    meshMin, meshMax, delta ;
    double                  h(0), dh ;
    int                     dimensions(3);

    STL->getBoundingBox( meshMin, meshMax ) ;

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

    // Compute level set in narrow band
    std::chrono::time_point<std::chrono::system_clock>    start, end;
    int                                         elapsed_init, elapsed_refi(0);

    bitpit::LevelSet                levelset;

    std::vector<bitpit::adaption::Info> mapper ;
    std::vector<double>             LS ;
    std::vector<double>::iterator   itLS ;

    levelset.setMesh(&mesh) ;
    levelset.addObject(std::move(STL),M_PI) ;

    mesh.getVTK().addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS) ;
    mesh.getVTK().setName("levelset_parallel_001_initial") ;

    levelset.setPropagateSign(true);

    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();

    if (rank == 0) {
        cout << " - Exporting data" << endl;

        LS.resize(mesh.getCellCount() ) ;
        itLS = LS.begin() ;
        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            *itLS = levelset.getLS(id) ;
            ++itLS ;
        };

        mesh.write() ;
    }

    // Partition the mesh
    mesh.getVTK().setCounter() ;
    mesh.getVTK().setName("levelset_parallel_001") ;

    mapper = mesh.partition(MPI_COMM_WORLD, true) ;
    levelset.update(mapper) ;

    // Write mesh
    if (rank == 0) {
        cout << " - Exporting data" << endl;
    }

    LS.resize(mesh.getCellCount() ) ;
    itLS = LS.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        *itLS = levelset.getLS(id) ;
        ++itLS ;
    };

    mesh.write() ;

//    std::fstream    file ;
//    file.open("levelset_004.dump", std::ios::out | std::ios::binary );
//
//    levelset.dump(file);
//    geometry.dump(file) ;
//    file.close();

//    {// dump levelset to file and restor into new class
//
//        bitpit::LevelSetOctree          levelset2(mesh);
//        bitpit::LevelSetSegmentation    geometry2(0,&STL);
//
//        std::fstream    file ;
//        file.open("levelset_004.dump", std::ios::in | std::ios::binary );
//
//        levelset2.restore(file);
//        geometry2.restore(file);
//
//
        //Refinement
        for( int i=0; i<3; ++i){

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
//    }

    cout << "elapsed time initialization " << elapsed_init << " ms" << endl;
    cout << "elapsed time refinement     " << elapsed_refi << " ms" << endl;

    MPI_Finalize();
}
