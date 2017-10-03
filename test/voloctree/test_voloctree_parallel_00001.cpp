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
#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing basic features of a parallel 2D patch.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 1.0;

	log::cout() << "  >> 2D octree patch" << "\n";

	// Create the patch
	VolOctree *patch_2D = new VolOctree(2, origin, length, dh);
	patch_2D->setCommunicator(MPI_COMM_WORLD);
	patch_2D->getVTK().setName("octree_parallel_uniform_patch_2D");
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

	return 0;
}

/*!
* Subtest 002
*
* Testing basic features of a parallel 2D patch.
*
* \param rank is the rank of the process
*/
int subtest_002(int rank)
{
	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 1.0;

	log::cout() << "  >> 3D octree mesh" << "\n";

	// Create the patch
	VolOctree *patch_3D = new VolOctree(3, origin, length, dh);
	patch_3D->setCommunicator(MPI_COMM_WORLD);
	patch_3D->getVTK().setName("octree_parallel_uniform_patch_3D");
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
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setVisibility(log::GLOBAL);

	// Run the subtests
    log::cout() << "Testing basic features of parallel octree patches" << std::endl;

	int status;
	try {
		status = subtest_001(rank);
		if (status != 0) {
			return status;
		}

		status = subtest_002(rank);
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		log::cout() << exception.what();
	}

	MPI_Finalize();
}
