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

#include "bitpit_octreepatch.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

	std::cout << "Testing adaption on octree patch" << std::endl;

#if ENABLE_MPI==1
	MPI::Init(argc,argv);
#endif

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 1.0;

	std::cout << std::endl;
	std::cout << "  :: 2D adaption test ::" << std::endl;

	std::cout << std::endl;
	std::cout << ">> Creating the patch" << std::endl;

	Patch *patch_2D = new OctreePatch(0, 2, origin, length, dh);
	patch_2D->setName("octree_adapted_patch_2D");
	patch_2D->update();
	patch_2D->writeMesh();

	std::cout << std::endl;
	std::cout << ">> Adapting patch" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch_2D->getCellCount();
		std::cout << std::endl;
		std::cout << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch_2D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_2D->extractCellNeighs(cellId)) {
				for (auto coarseId : patch_2D->extractCellNeighs(neighId)) {
					patch_2D->markCellForCoarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch_2D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_2D->extractCellNeighs(cellId)) {
				patch_2D->markCellForRefinement(neighId);
			}
		}

		std::cout << std::endl;
		std::cout << ">> Initial number of cells... " << nCells << std::endl;

		patch_2D->update();

		nCells = patch_2D->getCellCount();
		std::cout << ">> Final number of cells... " << nCells << std::endl;
	}
	patch_2D->writeMesh();

	std::cout << std::endl;
	std::cout << "  :: 3D adaption test ::" << std::endl;

	std::cout << std::endl;
	std::cout << ">> Creating patch" << std::endl;

	Patch *patch_3D = new OctreePatch(0, 3, origin, length, dh);
	patch_3D->setName("octree_adapted_patch_3D");
	patch_3D->update();
	patch_3D->writeMesh();

	std::cout << std::endl;
	std::cout << ">> Adapting patch" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch_3D->getCellCount();
		std::cout << std::endl;
		std::cout << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch_3D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_3D->extractCellNeighs(cellId)) {
				for (auto coarseId : patch_3D->extractCellNeighs(neighId)) {
					patch_3D->markCellForCoarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch_3D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_3D->extractCellNeighs(cellId)) {
				patch_3D->markCellForRefinement(neighId);
			}
		}

		std::cout << std::endl;
		std::cout << ">> Initial number of cells... " << nCells << std::endl;

		patch_3D->update();
		
		nCells = patch_3D->getCellCount();
		std::cout << ">> Final number of cells... " << nCells << std::endl;
	}
	patch_3D->writeMesh();

#if ENABLE_MPI==1
	MPI::Finalize();
#endif

}

