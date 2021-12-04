/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
* Subtest 001
*
* Testing refining a 2D patch up to the maximum level.
*/
int subtest_001()
{
    log::cout() << "  >> Checking maximum refinement for a 2D octree patch" << "\n";

    //
    // Create the patch
    //
    std::array<double, 3> origin = {{0., 0., 0.}};
    double length = 10;
    double dh = 10;

    log::cout() << "  >> 2D octree patch" << "\n";

#if BITPIT_ENABLE_MPI
    VolOctree patch_2D(2, origin, length, dh, MPI_COMM_NULL);
#else
    VolOctree patch_2D(2, origin, length, dh);
#endif
    patch_2D.getVTK().setName("octree_patch_at_maximum_level_2D");
    patch_2D.initializeAdjacencies();
    patch_2D.update();

    // Get the tree
    const PabloUniform &tree = patch_2D.getTree();

    // Refine the patch up to the maximum level
    int maximumLevel = patch_2D.getTree().getMaxLevel();
    for (int i = 0; i < maximumLevel; ++i) {
        patch_2D.markCellForRefinement(0);
        patch_2D.update();
    }

    const Octant *cellOctant_1 = tree.getOctant(patch_2D.getCellOctant(0).id);
    std::cout << " Level of cell 0 is " << static_cast<int>(cellOctant_1->getLevel()) << std::endl;
    if (cellOctant_1->getLevel() == maximumLevel) {
        std::cout << "   Cell has reached its maximum refinement level" << std::endl;
    } else {
        throw std::runtime_error("Cell has not reach its maximum refinement level.");
    }

    // Firther refinement should not alter the mesh
    patch_2D.markCellForRefinement(0);
    patch_2D.update();

    const Octant *cellOctant_2 = tree.getOctant(patch_2D.getCellOctant(0).id);
    std::cout << " Level of cell 0 is " << static_cast<int>(cellOctant_2->getLevel()) << std::endl;
    if (cellOctant_2->getLevel() == maximumLevel) {
        std::cout << "   Cell is still at its maximum refinement level" << std::endl;
    } else {
        throw std::runtime_error("Cell is not anymore at its maximum refinement level.");
    }

    // Write the patch
    patch_2D.write();

    return 0;
}

/*!
* Subtest 002
*
* Testing refining a 3D patch up to the maximum level.
*/
int subtest_002()
{
    log::cout() << "  >> Checking maximum refinement for a 3D octree patch" << "\n";

    //
    // Create the patch
    //
    std::array<double, 3> origin = {{0., 0., 0.}};
    double length = 10;
    double dh = 10;

    log::cout() << "  >> 3D octree patch" << "\n";

#if BITPIT_ENABLE_MPI
    VolOctree patch_3D(3, origin, length, dh, MPI_COMM_NULL);
#else
    VolOctree patch_3D(3, origin, length, dh);
#endif
    patch_3D.getVTK().setName("octree_patch_at_maximum_level_3D");
    patch_3D.initializeAdjacencies();
    patch_3D.update();

    // Get the tree
    const PabloUniform &tree = patch_3D.getTree();

    // Refine the patch up to the maximum level
    int maximumLevel = patch_3D.getTree().getMaxLevel();
    for (int i = 0; i < maximumLevel; ++i) {
        patch_3D.markCellForRefinement(0);
        patch_3D.update();
    }

    const Octant *cellOctant_1 = tree.getOctant(patch_3D.getCellOctant(0).id);
    std::cout << " Level of cell 0 is " << static_cast<int>(cellOctant_1->getLevel()) << std::endl;
    if (cellOctant_1->getLevel() == maximumLevel) {
        std::cout << "   Cell has reached its maximum refinement level" << std::endl;
    } else {
        throw std::runtime_error("Cell has not reach its maximum refinement level.");
    }

    // Firther refinement should not alter the mesh
    patch_3D.markCellForRefinement(0);
    patch_3D.update();

    const Octant *cellOctant_2 = tree.getOctant(patch_3D.getCellOctant(0).id);
    std::cout << " Level of cell 0 is " << static_cast<int>(cellOctant_2->getLevel()) << std::endl;
    if (cellOctant_2->getLevel() == maximumLevel) {
        std::cout << "   Cell is still at its maximum refinement level" << std::endl;
    } else {
        throw std::runtime_error("Cell is not anymore at its maximum refinement level.");
    }

    // Write the patch
    patch_3D.write();

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
    log::cout() << "Testing creating an octree patch from an existing tree" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }

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

