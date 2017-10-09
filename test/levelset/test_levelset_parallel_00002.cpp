
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
# include "bitpit_CG.hpp"
# include "bitpit_IO.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"

/*!
* Subtest 001
*
* Testing 2D levelset parallel refinement.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    uint8_t dimensions(2);

    // First Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (1,dimensions) );

    bitpit::log::cout() << " - Loading stl geometry" << std::endl;

    STL0->importDGF("./data/naca0012.dgf");

    STL0->deleteCoincidentVertices() ;
    STL0->buildAdjacencies() ;

    STL0->getVTK().setName("geometry_parallel_002_0") ;
    if(rank==0){
        STL0->write() ;
    }

    bitpit::log::cout() << "n. vertex: " << STL0->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL0->getCellCount() << std::endl;


    // Second Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL1( new bitpit::SurfUnstructured (1,dimensions) );

    bitpit::log::cout() << " - Loading stl geometry" << std::endl;

    STL1->importDGF("./data/square.dgf");

    STL1->deleteCoincidentVertices() ;
    STL1->buildAdjacencies() ;

    STL1->getVTK().setName("geometry_parallel_002_1") ;
    if(rank==0){
        STL1->write() ;
    }

    bitpit::log::cout() << "n. vertex: " << STL1->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL1->getCellCount() << std::endl;

    // Create mesh
    bitpit::log::cout() << " - Setting mesh" << std::endl;
    std::array<double,3> meshMin0, meshMax0;
    std::array<double,3> meshMin1, meshMax1;
    std::array<double,3> meshMin, meshMax, delta ;
    double h(0), dh ;

    STL0->getBoundingBox( meshMin0, meshMax0 ) ;
    STL1->getBoundingBox( meshMin1, meshMax1 ) ;
    bitpit::CGElem::unionAABB( meshMin0, meshMax0, meshMin1, meshMax1, meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    for( int i=0; i<3; ++i){
        h = std::max( h, meshMax[i]-meshMin[i] ) ;
    };

    dh = h / 16. ;
    bitpit::VolOctree    mesh(dimensions, meshMin, h, dh );
    mesh.update() ;


    // Configure levelset
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int elapsed_init, elapsed_part, elapsed_refi(0);

    bitpit::LevelSet levelset;

    std::vector<bitpit::adaption::Info> mapper ;
    int id0, id1, id2 ;
    std::vector<int> ids;

    levelset.setMesh(&mesh) ;
    id0 = levelset.addObject(std::move(STL0),M_PI) ;
    id1 = levelset.addObject(std::move(STL1),M_PI) ;
    id2 = levelset.addObject(bitpit::LevelSetBooleanOperation::INTERSECTION,id0,id1) ;
    ids.push_back(id0);
    ids.push_back(id1);
    ids.push_back(id2);

    bitpit::LevelSetObject &object0 = levelset.getObject(id0) ;
    bitpit::LevelSetObject &object1 = levelset.getObject(id1) ;
    bitpit::LevelSetObject &object2 = levelset.getObject(id2) ;

    object0.enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    object1.enableVTKOutput(bitpit::LevelSetWriteField::VALUE);
    object2.enableVTKOutput(bitpit::LevelSetWriteField::VALUE);

    levelset.setPropagateSign(true);

    // Compute levelset in narrow band
    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    if(rank==0){
        mesh.getVTK().setName("levelset_parallel_002_serial");
        mesh.write() ;
    }


    // Mesh Partitioning
    mapper = mesh.partition(MPI_COMM_WORLD, true) ;

    start = std::chrono::system_clock::now();
    levelset.partition(mapper) ;
    end = std::chrono::system_clock::now();

    elapsed_part = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    mesh.getVTK().setName("levelset_parallel_002_partitioned");
    mesh.write() ;

    // Refinement
    mesh.getVTK().setCounter() ;
    mesh.getVTK().setName("levelset_parallel_002_refinement");
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

        elapsed_refi += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

        mesh.write() ;
    }

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
