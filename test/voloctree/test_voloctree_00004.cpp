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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/*!
* Subtest 002
*
* Testing dump/restore of a 3D octree patch.
*
* \param patch_3D is the patch that will be created by the test
* \param patch_3D_restored is the patch that will be restored by the test
*/
int subtest_002()
{
    std::array<double, 3> origin = {{-8., -8., -8.}};
    double length = 16;
    double dh = 1;

    int archiveVersion = 1;

    log::cout() << "  >> 3D octree patch" << std::endl;

    // Create the patch
    log::cout() << "Creating 3D patch..." << std::endl;

    VolOctree *patch_3D = new VolOctree(3, origin, length, dh);
    patch_3D->getVTK().setName("octree_uniform_patch_3D");
	patch_3D->setCommunicator(MPI_COMM_WORLD);
    patch_3D->initializeAdjacencies();
    patch_3D->initializeInterfaces();
    patch_3D->update();


    // Refine the patch
    std::vector<uint32_t> refineList3D;
    refineList3D.push_back(2351);
    refineList3D.push_back(2365);
    refineList3D.push_back(2367);
    refineList3D.push_back(2422);
    refineList3D.push_back(2423);
    refineList3D.push_back(2431);
    refineList3D.push_back(2477);
    refineList3D.push_back(2479);
    refineList3D.push_back(2533);
    refineList3D.push_back(2534);
    refineList3D.push_back(2535);
    refineList3D.push_back(2541);
    refineList3D.push_back(2543);
    refineList3D.push_back(2559);
    refineList3D.push_back(2878);
    refineList3D.push_back(2988);
    refineList3D.push_back(2990);
    refineList3D.push_back(2997);
    refineList3D.push_back(3004);
    refineList3D.push_back(3005);
    refineList3D.push_back(3006);
    refineList3D.push_back(3375);
    refineList3D.push_back(3389);
    refineList3D.push_back(3391);
    refineList3D.push_back(3430);
    refineList3D.push_back(3431);
    refineList3D.push_back(3439);
    refineList3D.push_back(3445);
    refineList3D.push_back(3446);
    refineList3D.push_back(3447);
    refineList3D.push_back(3453);
    refineList3D.push_back(3455);
    refineList3D.push_back(3501);
    refineList3D.push_back(3503);
    refineList3D.push_back(3567);
    refineList3D.push_back(3886);
    refineList3D.push_back(3900);
    refineList3D.push_back(3902);
    refineList3D.push_back(4005);
    refineList3D.push_back(4012);
    refineList3D.push_back(4013);
    refineList3D.push_back(4014);

    for (long id : refineList3D) {
        patch_3D->markCellForRefinement(id);
    }

    for (int k = 0; k < 20; ++k) {
        long nCells = patch_3D->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 20; ++i) {
            long cellId = rand() % nCells;
            if (!patch_3D->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch_3D->findCellNeighs(cellId)) {
                patch_3D->markCellForRefinement(neighId);
            }
        }
        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch_3D->update();

        nCells = patch_3D->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }

    patch_3D->partition(MPI_COMM_WORLD);


    // Show patch info
    log::cout() << "Cell count:   " << patch_3D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_3D->getVertexCount() << std::endl;

    patch_3D->write();

    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = 0; i < nProcs; ++i) {
        if (i == rank) {
            for (PatchKernel::CellIterator itr = patch_3D->cellBegin(); itr != patch_3D->cellEnd(); ++itr) {
                for (int i = 0; i < itr->getVertexCount(); ++i) {
                    patch_3D->findCellVertexNeighs(itr.getId(), i);
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

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
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Testing dump/restore of octree patches" << std::endl;

    int status;
    try {
        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
