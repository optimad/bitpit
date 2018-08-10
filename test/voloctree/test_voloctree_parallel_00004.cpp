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
* Testing multi-layer ghost cells halo for 2D patches.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
	BITPIT_UNUSED(rank);

	// Patch data
	//
	// The number of ghosts layers in the halo should high enough for every
	// partition to contain all the global octants.
	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 8;
	double dh = 1.0;
	int nGhostLayers = 16;

	log::cout() << "  >> 2D octree patch" << "\n";

	// Create the patch
	VolOctree *patch_2D = new VolOctree(2, origin, length, dh);
	patch_2D->getVTK().setName("octree_parallel_uniform_patch_2D");
	patch_2D->update();

	// Partition the patch
	patch_2D->partition(MPI_COMM_WORLD, true, true, nGhostLayers);

	// Get information on patch numbering
	PatchNumberingInfo patchNumberingInfo(patch_2D);

	// Refine the mesh to obtain a non-uniform mesh
	for(const Cell &cell : patch_2D->getCells()){
		long globalId = patchNumberingInfo.getCellGlobalId(cell.getId());
		if (((globalId + 2) % 4) == 0) {
			patch_2D->markCellForRefinement(cell.getId());
		}
	}
	patch_2D->update(true);

	// Repartiion the patch
	patch_2D->partition(true);

	// Evaluate the ghost layer the cells belongs to
	const PatchKernel::CellConstRange &cellWriteRange = patch_2D->getVTKCellWriteRange();
	std::vector<int> cellHaloLayer(cellWriteRange.evalSize());

	std::size_t flatIndex = 0;
	for (const Cell &cell : cellWriteRange) {
		cellHaloLayer[flatIndex++] = patch_2D->getCellHaloLayer(cell.getId());
	}

	patch_2D->getVTK().addData("haloLayer", VTKFieldType::SCALAR, VTKLocation::CELL, cellHaloLayer);

	// Write the patch
	patch_2D->write();

	// Check if the ghost layer has been performed correctly
	//
	// Every partition contains all the global octant, this means that the
	// volume of the partition is equal to the volume of the whole patch.
	double partitionVolume = 0;
	for(const Cell &cell : patch_2D->getCells()){
		partitionVolume += patch_2D->evalCellVolume(cell.getId());
	}

	double globalVolume = length * length;

	log::cout() << "  Checking halo layers" << "\n";
	log::cout() << "\n";
	log::cout() << "    >> Partition volume = " << partitionVolume << "\n";
	log::cout() << "    >> Global volume = " << globalVolume << "\n";

	int errorFlag = 0;
	if (!utils::DoubleFloatingEqual()(partitionVolume, globalVolume)) {
		errorFlag = 1;
	}

	delete patch_2D;

	return errorFlag;
}

/*!
* Subtest 002
*
* Testing multi-layer ghost cells halo for 3D patches.
*
* \param rank is the rank of the process
*/
int subtest_002(int rank)
{
	BITPIT_UNUSED(rank);

	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 8;
	double dh = 1.0;
	int nGhostLayers = 16;

	log::cout() << "  >> 3D octree mesh" << "\n";

	// Create the patch
	VolOctree *patch_3D = new VolOctree(3, origin, length, dh);
	patch_3D->getVTK().setName("octree_parallel_uniform_patch_3D");
	patch_3D->update();

	// Partition the patch
	patch_3D->partition(MPI_COMM_WORLD, true, true, nGhostLayers);

	// Get information on patch numbering
	PatchNumberingInfo patchNumberingInfo(patch_3D);

	// Refine the mesh to obtain a non-uniform mesh
	for(const Cell &cell : patch_3D->getCells()){
		long globalId = patchNumberingInfo.getCellGlobalId(cell.getId());
		if (((globalId + 2) % 8) == 0) {
			patch_3D->markCellForRefinement(cell.getId());
		}
	}
	patch_3D->update(true);

	// Repartiion the patch
	patch_3D->partition(true);

	// Evaluate the ghost layer the cells belongs to
	const PatchKernel::CellConstRange &cellWriteRange = patch_3D->getVTKCellWriteRange();
	std::vector<int> cellHaloLayer(cellWriteRange.evalSize());

	std::size_t flatIndex = 0;
	for (const Cell &cell : cellWriteRange) {
		cellHaloLayer[flatIndex++] = patch_3D->getCellHaloLayer(cell.getId());
	}

	patch_3D->getVTK().addData("haloLayer", VTKFieldType::SCALAR, VTKLocation::CELL, cellHaloLayer);

	// Write the patch
	patch_3D->write();

	// Check if the ghost layer has been performed correctly
	//
	// Every partition contains all the global octant, this means that the
	// volume of the partition is equal to the volume of the whole patch.
	double partitionVolume = 0;
	for(const Cell &cell : patch_3D->getCells()){
		partitionVolume += patch_3D->evalCellVolume(cell.getId());
	}

	double globalVolume = length * length * length;

	log::cout() << "  Checking halo layers" << "\n";
	log::cout() << "\n";
	log::cout() << "    >> Partition volume = " << partitionVolume << "\n";
	log::cout() << "    >> Global volume = " << globalVolume << "\n";

	int errorFlag = 0;
	if (!utils::DoubleFloatingEqual()(partitionVolume, globalVolume)) {
		errorFlag = 1;
	}

	delete patch_3D;

	return errorFlag;
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
		exit(1);
	}

	MPI_Finalize();
}
