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

#include <array>

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/**
    Test for the 2D patch.

    \param origin is the origin of the domain
    \param length is the length of the domain
    \param dh is the maximum cell size
*/
void test_2D(const std::array<double, 3> &origin, double length, double dh)
{
    log::cout() << std::endl;
    log::cout() << "  :: 2D adaption test ::" << std::endl;

    log::cout() << std::endl;
    log::cout() << ">> Creating the patch" << std::endl;

    VolOctree *patch = new VolOctree(0, 2, origin, length, dh);
    patch->getVTK().setName("octree_adapted_patch_2D");
    patch->update();
    patch->write();

    log::cout() << std::endl;
    log::cout() << ">> Adapting patch" << std::endl;

    for (int k = 0; k < 5; ++k) {
        long nCells = patch->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 150; ++i) {
            long cellId = rand() % nCells;
            if (!patch->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch->findCellNeighs(cellId)) {
                for (auto coarseId : patch->findCellNeighs(neighId)) {
                    patch->markCellForCoarsening(coarseId);
                }
            }
        }

        for (int i = 0; i < 200; ++i) {
            long cellId = rand() % nCells;
            if (!patch->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch->findCellNeighs(cellId)) {
                patch->markCellForRefinement(neighId);
            }
        }

        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch->update();

        nCells = patch->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }
    patch->write();
}

/**
    Test for the 3D patch.

    \param origin is the origin of the domain
    \param length is the length of the domain
    \param dh is the maximum cell size
*/
void test_3D(const std::array<double, 3> &origin, double length, double dh)
{
    log::cout() << std::endl;
    log::cout() << "  :: 3D adaption test ::" << std::endl;

    log::cout() << std::endl;
    log::cout() << ">> Creating patch" << std::endl;

    VolOctree *patch = new VolOctree(0, 3, origin, length, dh);
    patch->getVTK().setName("octree_adapted_patch_3D");
    patch->update();
    patch->write();

    log::cout() << std::endl;
    log::cout() << ">> Adapting patch" << std::endl;

    for (int k = 0; k < 5; ++k) {
        long nCells = patch->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 150; ++i) {
            long cellId = rand() % nCells;
            if (!patch->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch->findCellNeighs(cellId)) {
                for (auto coarseId : patch->findCellNeighs(neighId)) {
                    patch->markCellForCoarsening(coarseId);
                }
            }
        }

        for (int i = 0; i < 200; ++i) {
            long cellId = rand() % nCells;
            if (!patch->getCells().exists(cellId)) {
                continue;
            }

            for (auto neighId : patch->findCellNeighs(cellId)) {
                patch->markCellForRefinement(neighId);
            }
        }

        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch->update();

        nCells = patch->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }
    patch->write();
}

/**
    Main proggram.
*/
int main(int argc, char *argv[])
{
    log::manager().initialize(log::COMBINED);
    log::cout() << "Testing adaption on octree patch" << std::endl;

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    std::array<double, 3> origin = {{0., 0., 0.}};
    double length = 20;
    double dh = 1.0;

    std::srand(1);

    test_2D(origin, length, dh);

    test_3D(origin, length, dh);

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

}
