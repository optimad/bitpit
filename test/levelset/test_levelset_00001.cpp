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
# define _USE_MATH_DEFINES

# include <cmath>
# include <ctime>
# include <chrono>

# include <array>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_IO.hpp"
# include "bitpit_levelset.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace bitpit;

void Generate2DSurfMesh(
    SurfUnstructured                   &mesh
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
const double                            R = 1.0;
const long                              N = 32;
double                                  theta;
double                                  dtheta = 2. * M_PI/((double) N);

// Counters
long                                    i;

// ========================================================================== //
// GENERATE MESH                                                              //
// ========================================================================== //
{
    // Local variables ------------------------------------------------------ //
    std::array<double, 3>                    point;
    std::vector<long>                        connect(2, Element::NULL_ID);

    // Create vertex list --------------------------------------------------- //
    point[2] = 0.0;
    for (i = 0; i < N; ++i) {
        theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        mesh.addVertex(point);
    } //next i

    // Create simplex list -------------------------------------------------- //
    for (i = 0; i < N; ++i) {
        connect[0] = i;
        connect[1] = (i+1) % N;
        mesh.addCell(ElementInfo::LINE, true, connect);
    } //next i
}

return; }

int main( int argc, char *argv[]){


#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc, &argv);
#endif
    int                    dimensions(2) ;

    // Input geometry
    std::unique_ptr<SurfUnstructured> STL( new SurfUnstructured(0,1,dimensions) );

    log::cout() << " - Loading dgf geometry" << std::endl;

    Generate2DSurfMesh( *(STL.get()) ) ;

    STL->getVTK().setName("geometry_001") ;
    STL->write() ;

    log::cout() << "n. vertex: " << STL->getVertexCount() << std::endl;
    log::cout() << "n. simplex: " << STL->getCellCount() << std::endl;



    // create cartesian mesh around geometry 
    log::cout() << " - Setting mesh" << std::endl;
    std::array<double,3>     meshMin, meshMax, delta ;
    std::array<int,3>        nc = {{64, 64, 0}} ;

    STL->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    VolCartesian mesh( 1, dimensions, meshMin, delta, nc);

    // Compute level set  in narrow band
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int                                      elapsed_seconds;

    LevelSet                levelset ;

    levelset.setMesh(&mesh) ;

    levelset.addObject(std::move(STL),M_PI) ;

    start = std::chrono::system_clock::now();
    levelset.compute( ) ;
    end = std::chrono::system_clock::now();

    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    log::cout() << " - Exporting data" << std::endl;
    mesh.update() ;
    std::vector<double> LS(mesh.getCellCount() ) ;
    std::vector<double>::iterator it = LS.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        *it = levelset.getLS(id) ;
        ++it ;
    };

    mesh.getVTK().addData("ls", VTKFieldType::SCALAR, VTKLocation::CELL, LS) ;
    mesh.getVTK().setName("levelset_001") ;
    mesh.write() ;

    log::cout() << " - Exported data" << std::endl;

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return 0;

};


