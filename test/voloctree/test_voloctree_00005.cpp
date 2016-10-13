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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    log::manager().initialize(log::COMBINED);
    log::cout() << "Testing creating a patch from an existing tree" << std::endl;

    //
    // Create the tree
    //
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    PabloUniform octree(x_0, y_0, z_0, l, 2);

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    // Refine and write the octree
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    std::vector<uint32_t> refineList;
    refineList.push_back(  7);
    refineList.push_back( 13);
    refineList.push_back( 15);
    refineList.push_back( 26);
    refineList.push_back( 27);
    refineList.push_back( 31);
    refineList.push_back( 37);
    refineList.push_back( 39);
    refineList.push_back( 49);
    refineList.push_back( 50);
    refineList.push_back( 51);
    refineList.push_back( 53);
    refineList.push_back( 55);
    refineList.push_back( 63);
    refineList.push_back( 78);
    refineList.push_back(100);
    refineList.push_back(102);
    refineList.push_back(105);
    refineList.push_back(108);
    refineList.push_back(109);
    refineList.push_back(110);
    refineList.push_back(135);
    refineList.push_back(141);
    refineList.push_back(143);
    refineList.push_back(146);
    refineList.push_back(147);
    refineList.push_back(151);
    refineList.push_back(153);
    refineList.push_back(154);
    refineList.push_back(155);
    refineList.push_back(157);
    refineList.push_back(159);
    refineList.push_back(165);
    refineList.push_back(167);
    refineList.push_back(183);
    refineList.push_back(198);
    refineList.push_back(204);
    refineList.push_back(206);
    refineList.push_back(225);
    refineList.push_back(228);
    refineList.push_back(229);
    refineList.push_back(230);

    for (uint32_t id : refineList) {
        octree.setMarker(id, 3);
    }
    octree.adapt(false);

#if BITPIT_ENABLE_MPI==1
    octree.loadBalance();
#endif

    //
    // Create the patch from the existing tree
    //

    VolOctree *patch_2D;
    std::unique_ptr<PabloUniform> treePointer = std::unique_ptr<PabloUniform>(&octree);

    // Create the patch
    patch_2D = new VolOctree(0, std::move(treePointer), &treePointer);
    patch_2D->getVTK().setName("octree_patch_from_tree_2D_initial");
    patch_2D->update();

    // Show patch info
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->write();

    // Refine the patch
    for (int k = 0; k < 5; ++k) {
        long nCells = patch_2D->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 200; ++i) {
            long cellId = rand() % nCells;
            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch_2D->findCellNeighs(cellId)) {
                patch_2D->markCellForRefinement(neighId);
            }
        }
        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch_2D->update();

        nCells = patch_2D->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }

    // Show patch info
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->getVTK().setName("octree_patch_from_tree_2D_refined");
    patch_2D->write();

    // Destroy the patch
    delete patch_2D;

    // The tree has now being adopted and can used again.
    patch_2D = new VolOctree(0, std::move(treePointer), &treePointer);
    patch_2D->getVTK().setName("octree_patch_from_adopted_tree_2D_initial");
    patch_2D->update();

    // Show patch info
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

}
