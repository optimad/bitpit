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
* Testing coarsening of a 2D patch that merges octants from different
* processors.
*/
int subtest_001()
{
	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 20;
	double dh = 2;

	log::cout() << "  >> 2D octree patch" << "\n";

	// Create the patch
	VolOctree *patch_2D_original = new VolOctree(2, origin, length, dh);
	patch_2D_original->setCommunicator(MPI_COMM_WORLD);
	patch_2D_original->getVTK().setName("octree_parallel_uniform_patch_2D");
	patch_2D_original->buildInterfaces();
	patch_2D_original->update();
	patch_2D_original->getVTK().setCounter(0);

	// Clone the patch
	std::unique_ptr<VolOctree> patch_2D = PatchKernel::clone(patch_2D_original);
	delete patch_2D_original;

	// Partition the patch
	patch_2D->partition(true);
	patch_2D->write();

	// Coarsening
	bool keepCoarsening_2D = true;
	while (keepCoarsening_2D) {
		for(const auto &cell : patch_2D->getCells()){
			if (!cell.isInterior()) {
				continue;
			}

			patch_2D->markCellForCoarsening(cell.getId());
		}

		// Update the patch
		const std::vector<adaption::Info> adaptionInfo = patch_2D->update(true);
		patch_2D->write();

		// Show adaption info
		log::cout() << " Current internal cells count: " << patch_2D->getInternalCount() << std::endl;
		for (const auto &adaptionData : adaptionInfo) {
			if (adaptionData.entity != adaption::ENTITY_CELL) {
				continue;
			}

			log::cout() << "   Adpation type:     " << adaptionData.type << std::endl;
			log::cout() << "   Adpation entity:   " << adaptionData.entity << std::endl;
			log::cout() << "   Adpation rank:     " << adaptionData.rank << std::endl;
			log::cout() << "   Adpation previous: " << adaptionData.previous << std::endl;
			log::cout() << "   Adpation current:  " << adaptionData.current << std::endl;
			log::cout() << std::endl;
		}

		// If one partition has only one cell stop
		bool local_keepCoargening = (adaptionInfo.size() > 0);
		MPI_Allreduce(&local_keepCoargening, &keepCoarsening_2D, 1,
					  MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
	}

	return 0;
}

/*!
* Subtest 002
*
* Testing coarsening of a 3D patch that merges octants from different
* processors.
*/
int subtest_002()
{
	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 20;
	double dh = 2;

	log::cout() << "  >> 3D octree mesh" << "\n";

	// Create the patch
	VolOctree *patch_3D_original = new VolOctree(3, origin, length, dh);
	patch_3D_original->setCommunicator(MPI_COMM_WORLD);
	patch_3D_original->getVTK().setName("octree_parallel_uniform_patch_3D");
	patch_3D_original->buildInterfaces();
	patch_3D_original->update();
	patch_3D_original->getVTK().setCounter(0);

	// Clone the patch
	std::unique_ptr<VolOctree> patch_3D = PatchKernel::clone(patch_3D_original);
	delete patch_3D_original;

	// Partition the patch
	patch_3D->partition(true);
		patch_3D->write();

	// Coarsening
	bool keepCoarsening_3D = true;
	while (keepCoarsening_3D) {
		for(const auto &cell : patch_3D->getCells()){
			if (!cell.isInterior()) {
				continue;
			}

			patch_3D->markCellForCoarsening(cell.getId());
		}

		// Update the patch
		const std::vector<adaption::Info> adaptionInfo = patch_3D->update(true);
		patch_3D->write();

		// Show adaption info
		log::cout() << " Current internal cells count: " << patch_3D->getInternalCount() << std::endl;
		for (const auto &adaptionData : adaptionInfo) {
			if (adaptionData.entity != adaption::ENTITY_CELL) {
				continue;
			}

			log::cout() << "   Adpation type:     " << adaptionData.type << std::endl;
			log::cout() << "   Adpation entity:   " << adaptionData.entity << std::endl;
			log::cout() << "   Adpation rank:     " << adaptionData.rank << std::endl;
			log::cout() << "   Adpation previous: " << adaptionData.previous << std::endl;
			log::cout() << "   Adpation current:  " << adaptionData.current << std::endl;
			log::cout() << std::endl;
		}

		// If one partition has only one cell stop
		bool local_keepCoargening = (adaptionInfo.size() > 0);
		MPI_Allreduce(&local_keepCoargening, &keepCoarsening_3D, 1,
					  MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
	}

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
    log::cout() << "Testing coarsening that merges octants from different processors" << std::endl;

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

	MPI_Finalize();
}
