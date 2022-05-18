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
#if BITPIT_ENABLE_MPI==1

#include <mpi.h>
#include <unordered_set>

#include "voloctree.hpp"

using namespace std;

namespace bitpit {

/*!
	Initialize tree partitioning.

	All octants are moved to the process identified by the rank zero in the
	communicator.
*/
void VolOctree::initializeTreePartitioning()
{
	// Move all the octants to the first processor
	//
	// Assigning a null weight to every octant but the last octant and then
	// doing a load balance will have the effect of moving all the octants
	// to the first process.
	std::size_t nOctants = m_tree->getNumOctants();
	std::vector<double> octantWeights(nOctants, 0.);
	if (nOctants > 0) {
		octantWeights[nOctants - 1] = 1.;
	}
	m_tree->loadBalance(m_partitioningOctantWeights.get());
}

/*!
	Initialize the size, expressed in number of layers, of the tree ghost
	cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void VolOctree::initializeTreeHaloSize(std::size_t haloSize)
{
	m_tree->setNofGhostLayers(haloSize);
}

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void VolOctree::setCommunicator(MPI_Comm communicator)
{
	PatchKernel::setCommunicator(communicator);

	m_tree->setComm(communicator);
}

/*!
	Gets the maximum allowed size, expressed in number of layers, of the ghost
	cells halo.

	\result The maximum allowed size, expressed in number of layers, of the
	ghost cells halo.
*/
std::size_t VolOctree::_getMaxHaloSize()
{
	return std::numeric_limits<uint32_t>::max();
}

/*!
	Internal function to set the size, expressed in number of layers, of the
	ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void VolOctree::_setHaloSize(std::size_t haloSize)
{
	initializeTreeHaloSize(haloSize);
}

/*!
	Prepares the patch for performing the partitioning.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated to a specified cell. If no weight
	is specified for a cell, the default weight will be used
	\param defaultWeight is the default weight that will assigned to the cells
	for which an explicit weight has not been defined
	\param trackPartitioning if set to true the function will return the
	changes that will be performed in the alter step
	\result If the partitioning is tracked, returns a vector of adaption::Info
	that can be used to discover what changes will be performed in the alter
	step, otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> VolOctree::_partitioningPrepare(const std::unordered_map<long, double> &cellWeights, double defaultWeight, bool trackPartitioning)
{
	// Eval partitioning weights
	computePartitioningOctantWeights(cellWeights, defaultWeight);

	// Generate partitioning information
	std::vector<adaption::Info> partitioningData;
	if (trackPartitioning) {
		int currentRank = getRank();
		PabloUniform::LoadBalanceRanges loadBalanceRanges = m_tree->evalLoadBalanceRanges(m_partitioningOctantWeights.get());
		for (const auto &entry : loadBalanceRanges.sendRanges) {
			int receiver = entry.first;
			if (receiver == currentRank) {
				continue;
			}

			adaption::Type adaptionType;
			if (loadBalanceRanges.sendAction == PabloUniform::LoadBalanceRanges::ACTION_DELETE) {
				adaptionType = adaption::TYPE_DELETION;
			} else {
				adaptionType = adaption::TYPE_PARTITION_SEND;
			}

			partitioningData.emplace_back();
			adaption::Info &partitioningInfo = partitioningData.back();
			partitioningInfo.entity = adaption::ENTITY_CELL;
			partitioningInfo.type   = adaptionType;
			partitioningInfo.rank   = receiver;

			const std::array<uint32_t, 2> &sendRange = entry.second;
			uint32_t beginTreeId = sendRange[0];
			uint32_t endTreeId   = sendRange[1];
			partitioningInfo.previous.reserve(endTreeId - beginTreeId);
			for (uint32_t treeId = beginTreeId; treeId < endTreeId; ++treeId) {
				OctantInfo octantInfo(treeId, true);
				partitioningInfo.previous.emplace_back();
				long &cellId = partitioningInfo.previous.back();
				cellId = getOctantId(octantInfo);
			}
		}
	}

	return partitioningData;
}

/*!
	Alter the patch performing the partitioning.

	Default implementation is a no-op function.

	\param trackPartitioning if set to true the function will return the changes
	done to the patch during the partitioning
	\result If the partitioning is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the adaption, otherwise an
	empty vector will be returned.
*/
std::vector<adaption::Info> VolOctree::_partitioningAlter(bool trackPartitioning)
{
	std::vector<adaption::Info> partitioningData;

	// Early return if the dimension of the tree is null
	if (m_tree->getDim() == 0) {
		return partitioningData;
	}

	// Updating the tree
	m_tree->loadBalance(m_partitioningOctantWeights.get());

	// Sync the patch
	partitioningData = sync(trackPartitioning);

	// The bounding box is frozen, it is not updated automatically
	setBoundingBox();

	return partitioningData;
}

/*!
	Cleanup patch data structured after the partitioning.
*/
void VolOctree::_partitioningCleanup()
{
	// Clear partitioning weights
	clearPartitioningOctantWeights();
}

/*!
	Finds the internal cells that will be ghost cells for the processes
	with the specified ranks. During data exchange, these cells will be
	the sources form which data will be read from.

	\param rank is the rank for which the information will be built
*/
std::vector<long> VolOctree::_findGhostCellExchangeSources(int rank)
{
	const std::vector<uint32_t> &ghostOctantSources = m_tree->getBordersPerProc().at(rank);
	std::size_t nGhostSources = ghostOctantSources.size();

	std::vector<long> ghostCellSources(nGhostSources);
	for (std::size_t k = 0; k < nGhostSources; ++k) {
		OctantInfo octantInfo(ghostOctantSources[k], true);
		ghostCellSources[k] = getOctantId(octantInfo);
	}

	return ghostCellSources;
}

/*!
	Compute partitioning weights for the octants.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated to a specified cell. If no weight
	is specified for a cell, a default weight will be used
	\param defaultWeight is the default weight that will assigned to the cells
	for which an explicit weight has not been defined
*/
void VolOctree::computePartitioningOctantWeights(const std::unordered_map<long, double> &cellWeights, double defaultWeight)
{
	bool useWeights = !cellWeights.empty();
	MPI_Allreduce(MPI_IN_PLACE, &useWeights, 1, MPI_C_BOOL, MPI_LOR, getCommunicator());
	if (!useWeights) {
		clearPartitioningOctantWeights();
		return;
	}

	CellConstIterator beginItr = internalCellConstBegin();
	CellConstIterator endItr   = internalCellConstEnd();

	std::size_t nOctants = m_tree->getNumOctants();
	m_partitioningOctantWeights = std::unique_ptr<std::vector<double>>(new std::vector<double>(nOctants, defaultWeight));
	for (CellConstIterator cellItr = beginItr; cellItr != endItr; ++cellItr) {
		long cellId = cellItr.getId();
		auto weightItr = cellWeights.find(cellId);
		if (weightItr != cellWeights.end()) {
			OctantInfo octantInfo = getCellOctant(cellId);
			uint32_t treeId = octantInfo.id;
			(*m_partitioningOctantWeights)[treeId] = weightItr->second;
		}
	}
}

/*!
	Compute partitioning weights for the octants.
*/
void VolOctree::clearPartitioningOctantWeights()
{
	m_partitioningOctantWeights.reset();
}

}

#endif
