
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

#include <ctime>
#include <chrono>

#include <mpi.h>

#include "bitpit_IO.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_voloctree.hpp"
#include "bitpit_levelset.hpp"

using namespace std;

/*!
* Subtest 001
*
* Testing creation of a levelset from an existing tree.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    // Info
    bitpit::log::cout() << "Testing creating a levelset from an existing tree" << std::endl;

    // Input geometry
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(2, 3) );

    bitpit::log::cout() << " - Loading stl geometry" << std::endl;

    STL->importSTL("./data/cube.stl", true);

    STL->deleteCoincidentVertices();
    STL->buildAdjacencies();

    STL->getVTK().setName("geometry_003");
    if (rank == 0) {
        STL->write();
    }

    bitpit::log::cout() << "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL->getCellCount() << std::endl;

    // Create the octree
    bitpit::log::cout() << " - Creating the octree" << std::endl;

    std::array<double,3> meshMin, meshMax;
    STL->getBoundingBox(meshMin, meshMax);

    std::array<double,3> delta = meshMax -meshMin;
    meshMin -= 0.1 * delta;
    meshMax += 0.1 * delta;
    delta = meshMax -meshMin;

    double h = 0.;
    for (int i=0; i<3; ++i) {
        h = max(h, meshMax[i] - meshMin[i]);
    }

    int dimensions = 3;

    std::unique_ptr<bitpit::PabloUniform> octree = std::unique_ptr<bitpit::PabloUniform>(new bitpit::PabloUniform(meshMin[0], meshMin[1], meshMin[2], h, dimensions));
    octree->adaptGlobalRefine();
    octree->adaptGlobalRefine();
    octree->adaptGlobalRefine();
    octree->adaptGlobalRefine();

    octree->loadBalance();

    // Create the mesh
    bitpit::log::cout() << " - Creating the mesh" << std::endl;

    bitpit::VolOctree mesh(0, std::move(octree));
    mesh.update();

    mesh.getVTK().setCounter();
    mesh.getVTK().setName("levelset_parallel_003");

    // Compute level set in narrow band
    std::chrono::time_point<std::chrono::system_clock>    start, end;
    int                                         elapsed_init, elapsed_refi(0);

    bitpit::LevelSet                levelset;

    std::vector<bitpit::adaption::Info> mapper;
    std::vector<double>             LS;
    std::vector<double>::iterator   itLS;

    levelset.setMesh(&mesh);
    int id0 = levelset.addObject(std::move(STL),M_PI);
    const bitpit::LevelSetObject &object = levelset.getObject(id0);

    mesh.getVTK().addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS);

    levelset.setPropagateSign(true);

    start = std::chrono::system_clock::now();
    levelset.compute( );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();

    // Write mesh
    if (rank == 0) {
        bitpit::log::cout() << " - Exporting data" << endl;
    }

    LS.resize(mesh.getCellCount() );
    itLS = LS.begin();
    for (auto & cell : mesh.getCells() ){
        const long &id = cell.getId();
        *itLS = object.getLS(id);
        ++itLS;
    }

    mesh.write();

    // Refinement
    for (int i=0; i<3; ++i){

        for (auto & cell : mesh.getCells() ){
            const long &id = cell.getId();
            if( std::abs(object.getLS(id)) < 100. ){
                mesh.markCellForRefinement(id);
            }
        }

        mapper = mesh.update(true);
        start = std::chrono::system_clock::now();
        levelset.update(mapper);
        end = std::chrono::system_clock::now();

        elapsed_refi += chrono::duration_cast<chrono::milliseconds>(end-start).count();

        LS.resize(mesh.getCellCount() );
        itLS = LS.begin();
        for (auto & cell : mesh.getCells() ){
            const long &id = cell.getId();
            *itLS = object.getLS(id);
            ++itLS;
        }

        mesh.write();
    }

    bitpit::log::cout() << " Elapsed time initialization " << elapsed_init << " ms" << endl;
    bitpit::log::cout() << " Elapsed time refinement     " << elapsed_refi << " ms" << endl;

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
    bitpit::log::cout() << "Testing creation of a levelset from an existing tree" << std::endl;

    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        bitpit::log::cout() << exception.what();
    }

    MPI_Finalize();
}
