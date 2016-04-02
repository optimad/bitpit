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

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

	MPI_Init(&argc,&argv);

	int nProcs;
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setVisibility(log::GLOBAL);
	log::cout() << "Testing octree patch" << "\n";

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 1.0;

	// 2D Test
	log::cout() << "  >> 2D octree patch" << "\n";

	// Create the patch
	VolOctree *patch_2D = new VolOctree(0, 2, origin, length, dh);
	patch_2D->setCommunicator(MPI_COMM_WORLD);
	patch_2D->setName("octree_parallel_uniform_patch_2D");
	patch_2D->update();

	// Partition the patch
	patch_2D->partition(true);

	// Refine one processor
	if(rank==0){
		for(const auto& cell:patch_2D->getCells()){
			patch_2D->markCellForRefinement(cell.getId());
		}
	}
	patch_2D->update(true);

	// Repartiion the patch
	patch_2D->partition(true);

	// Write the patch
	patch_2D->write();

	delete patch_2D;

	// 3D Test
	log::cout() << "  >> 3D octree mesh" << "\n";

	// Create the patch
	VolOctree *patch_3D = new VolOctree(0, 3, origin, length, dh);
	patch_3D->setCommunicator(MPI_COMM_WORLD);
	patch_3D->setName("octree_parallel_uniform_patch_3D");
	patch_3D->update();

	// Partition the patch
	patch_3D->partition(true);

	// Refine one processor
	if(rank==0){
		for(const auto& cell:patch_3D->getCells()){
			patch_3D->markCellForRefinement(cell.getId());
		}
	}
	patch_3D->update(true);

	// Repartiion the patch
	patch_3D->partition(true);

	// Write the patch
	patch_3D->write();

	delete patch_3D;

	MPI_Finalize();

}
