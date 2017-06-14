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
#if BITPIT_ENABLE_MPI==1

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include <mpi.h>
#include <chrono>
#include <unordered_set>

#include "bitpit_SA.hpp"

#include "patch_kernel.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace chrono;

namespace bitpit {

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void PatchKernel::setCommunicator(MPI_Comm communicator)
{
	// Communication can be set just once
	if (isCommunicatorSet()) {
		throw std::runtime_error ("Patch communicator can be set just once");
	}

	// The communicator has to be valid
	if (communicator == MPI_COMM_NULL) {
		throw std::runtime_error ("Patch communicator is not valid");
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	MPI_Comm_dup(communicator, &m_communicator);

	// Get MPI information
	MPI_Comm_size(m_communicator, &m_nProcessors);
	MPI_Comm_rank(m_communicator, &m_rank);

	// Set parallel data for the VTK output
	if (m_nProcessors > 1) {
		m_vtk.setParallel(m_nProcessors, m_rank);
	}
}

/*!
	Checks if the communicator to be used for parallel communications has
	already been set.

	\result Returns true if the communicator has been set, false otherwise.
*/
bool PatchKernel::isCommunicatorSet() const
{
	return (getCommunicator() != MPI_COMM_NULL);
}

/*!
	Gets the MPI communicator associated to the patch

	\return The MPI communicator associated to the patch.
*/
const MPI_Comm & PatchKernel::getCommunicator() const
{
	return m_communicator;
}

/*!
	Frees the MPI communicator associated to the patch
*/
void PatchKernel::freeCommunicator()
{
	if (!isCommunicatorSet()) {
		return;
	}

	int finalizedCalled;
	MPI_Finalized(&finalizedCalled);
	if (finalizedCalled) {
		return;
	}

	MPI_Comm_free(&m_communicator);
}

/*!
	Gets the MPI rank associated to the patch

	\return The MPI rank associated to the patch.
*/
int PatchKernel::getRank() const
{
	return m_rank;
}

/*!
	Gets the MPI processors in the communicator associated to the patch

	\return The MPI processors in the communicator associated to the patch
*/
int PatchKernel::getProcessorCount() const
{
	return m_nProcessors;
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackChanges if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
const std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackChanges, bool squeezeStorage)
{
	setCommunicator(communicator);

	return partition(cellRanks, trackChanges, squeezeStorage);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning.
	\param trackChanges if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
const std::vector<adaption::Info> PatchKernel::partition(const std::vector<int> &cellRanks, bool trackChanges, bool squeezeStorage)
{
	std::vector<adaption::Info> adaptionData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check if the patch allow custom partition
	if (!isExpert()) {
		log::cout() << "The patch does not allow custom partition" << std::endl;

		return adaptionData;
	}

	// Build the send map
	std::unordered_map<int, std::vector<long>> sendMap;

	auto cellItr = cellBegin();
	for (int k = 0; k < getInternalCount(); ++k) {
		const int &rank = cellRanks[k];
		if (rank == getRank()) {
			cellItr++;
			continue;
		}

		sendMap[rank].push_back(cellItr->getId());

		cellItr++;
	}

	// Local sender-receiver pairs
	std::vector<int> localPairs;
	localPairs.reserve(2 * sendMap.size());
	for (const auto &entry : sendMap) {
		localPairs.push_back(getRank());
		localPairs.push_back(entry.first);
	}

	// Exchange the size of the communication
	int globalPairsSize = localPairs.size();
	std::vector<int> globalPairsSizes(getProcessorCount());
	MPI_Allgather(&globalPairsSize, 1, MPI_INT, globalPairsSizes.data(), 1, MPI_INT, getCommunicator());

	std::vector<int> globalPairsOffsets(getProcessorCount());
	globalPairsOffsets[0] = 0;
	for (int i = 1; i < getProcessorCount(); ++i) {
		globalPairsOffsets[i] = globalPairsOffsets[i-1] + globalPairsSizes[i-1];
	}

	// Global sender-receiver pairs
	std::vector<int> globalPairs;
	globalPairs.resize(globalPairsOffsets.back() + globalPairsSizes.back());

	MPI_Allgatherv(localPairs.data(), localPairs.size(), MPI_INT, globalPairs.data(),
				   globalPairsSizes.data(), globalPairsOffsets.data(), MPI_INT,
                   getCommunicator());

	// Exchange the cells
	std::vector<long> emptyCellList;
	for (size_t i = 0; i < globalPairs.size(); i += 2) {
		int srcRank = globalPairs[i];
		int dstRank = globalPairs[i+1];

		std::vector<long> *cellList;
		if (srcRank == getRank()) {
			cellList = &(sendMap[dstRank]);
		} else {
			cellList = &emptyCellList;
		}

		adaption::Info adaptionInfo = sendCells(srcRank, dstRank, *cellList);
		if (trackChanges && adaptionInfo.type != adaption::TYPE_NONE) {
			adaptionData.push_back(std::move(adaptionInfo));
		}
	}

	// Squeeze the storage
	if (squeezeStorage) {
		squeeze();
	}

	// Patch is now partitioned
	setPartitioned(true);

	return adaptionData;
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param communicator is the communicator that will be used
	\param trackChanges if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
const std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, bool trackChanges, bool squeezeStorage)
{
	setCommunicator(communicator);

	return partition(trackChanges, squeezeStorage);
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param trackChanges if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partition.
*/
const std::vector<adaption::Info> PatchKernel::partition(bool trackChanges, bool squeezeStorage)
{
	return balancePartition(trackChanges, squeezeStorage);
}

/*!
	Tries to balance the computational load among the processors redistributing
	the cells among the processors.

	\param trackChanges if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
const std::vector<adaption::Info> PatchKernel::balancePartition(bool trackChanges, bool squeezeStorage)
{
	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Balance patch
	const std::vector<adaption::Info> adaptionData = _balancePartition(trackChanges, squeezeStorage);

	// Update the bouding box
	updateBoundingBox();

	// Patch is now partitioned
	setPartitioned(true);

	// Done
	return adaptionData;
}

/*!
	Checks if the patch has been partitioned.

	\result Returns true if the patch has been partitioned, false otherwise.
*/
bool PatchKernel::isPartitioned() const
{
	return m_partitioned;
}

/*!
	Sets the partitioned flag.

	\param partitioned is the flag that will be set
*/
void PatchKernel::setPartitioned(bool partitioned)
{
	m_partitioned = partitioned;
}

/*!
	Internal function that tries to balance the computational load among the
	processors moving redistributing the cells among the processors.

	\param trackChanges if set to true the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<adaption::Info> PatchKernel::_balancePartition(bool trackChanges, bool squeezeStorage)
{
	BITPIT_UNUSED(trackChanges);
	BITPIT_UNUSED(squeezeStorage);

	log::cout() << "The patch does not implement a algortihm for balacing the partition" << std::endl;

	return std::vector<adaption::Info>();
}

/*!
	Gets the rank of the processor that owns the specified cell.

	\param id is the id of the requested cell
	\result The rank that owns the specified cell.
*/
int PatchKernel::getCellRank(const long &id) const
{
	const Cell &cell = getCell(id);
	if (cell.isInterior()) {
		return m_rank;
	} else {
		return m_ghostOwners.at(id);
	}
}

/*!
	Check if the processors associated to the specified rank is a neighbour.

	\param rank is the rank associated to the processor
	\result True is the processor is a neighbour, false otherwise.
*/
bool PatchKernel::isRankNeighbour(int rank)
{
	return (m_ghostExchangeTargets.count(rank) > 0);
}

/*!
	Get a list of neighbour ranks.

	\result A list of neighbour ranks.
*/
std::vector<int> PatchKernel::getNeighbourRanks()
{
	std::vector<int> neighRanks;
	neighRanks.reserve(m_ghostExchangeTargets.size());
	for (const auto &entry : m_ghostExchangeTargets) {
		neighRanks.push_back(entry.first);
	}

	return neighRanks;
}

/*!
	Gets a reference to the ghost targets needed for data exchange.

	\result A reference to the ghost targets needed for data exchange.
*/
std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeTargets()
{
	return m_ghostExchangeTargets;
}

/*!
	Gets a constant reference to the ghost targets needed for data exchange.

	\result A constant reference to the ghost targets needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeTargets() const
{
	return m_ghostExchangeTargets;
}

/*!
	Gets a reference to the ghost targets needed for data exchange for
	the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A reference to the ghost targets needed for data exchange for
	the specified rank.
*/
std::vector<long> & PatchKernel::getGhostExchangeTargets(int rank)
{
	return m_ghostExchangeTargets.at(rank);
}

/*!
	Gets a constant reference to the ghost targets needed for data
	exchange for the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A constant reference to the ghost targets needed for data
	exchange for the specified rank.
*/
const std::vector<long> & PatchKernel::getGhostExchangeTargets(int rank) const
{
	return m_ghostExchangeTargets.at(rank);
}

/*!
	Gets a reference to the ghost sources needed for data exchange.

	\result A reference to the ghost sources needed for data exchange.
*/
std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeSources()
{
	return m_ghostExchangeSources;
}

/*!
	Gets a constant reference to the ghost sources needed for data exchange.

	\result A constant reference to the ghost sources needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeSources() const
{
	return m_ghostExchangeSources;
}

/*!
	Gets a reference to the ghost sources needed for data exchange for
	the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A reference to the ghost sources needed for data exchange for
	the specified rank.
*/
std::vector<long> & PatchKernel::getGhostExchangeSources(int rank)
{
	return m_ghostExchangeSources.at(rank);
}

/*!
	Gets a constant reference to the ghost sources needed for data
	exchange for the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A constant reference to the ghost sources needed for data
	exchange for the specified rank.
*/
const std::vector<long> & PatchKernel::getGhostExchangeSources(int rank) const
{
	return m_ghostExchangeSources.at(rank);
}

/*!
	Sets the owner of the specified ghost.

	\param id is the id of the ghost cell
	\param rank is the rank of the processors that owns the ghost cell
	\param updateExchangeData if set to true exchange data will be updated
*/
void PatchKernel::setGhostOwner(int id, int rank, bool updateExchangeData)
{
	// Rebuild the exchange data information of the previous owner
	if (updateExchangeData) {
		if (m_ghostOwners.count(id) > 0) {
			removeGhostFromExchangeTargets(id);
		}
	}

	// Assign the owner to the cell
	m_ghostOwners[id] = rank;

	// Rebuild the exchange data information of the current owner
	if (updateExchangeData) {
		addGhostToExchangeTargets(id);
	}
}

/*!
	Unsets the owner of the specified ghost.

	\param id is the id of the ghost cell
	\param updateExchangeData if set to true exchange data will be updated
*/
void PatchKernel::unsetGhostOwner(int id, bool updateExchangeData)
{
	if (m_ghostOwners.count(id) <= 0) {
		return;
	}

	// Rebuild the exchange data information of the previous owner
	if (updateExchangeData) {
		removeGhostFromExchangeTargets(id);
	}

	// Remove the owner
	m_ghostOwners.erase(id);
}

/*!
	Clear the owners of all the ghosts.

	\param updateExchangeData if set to true exchange data will be updated
*/
void PatchKernel::clearGhostOwners(bool updateExchangeData)
{
	// Clear the owners
	m_ghostOwners.clear();

	// Clear exchange data
	if (updateExchangeData) {
		deleteGhostExchangeData();
	}
}

/*!
	Reset the ghost information needed for data exchange.
*/
void PatchKernel::deleteGhostExchangeData()
{
	m_ghostExchangeTargets.clear();
	m_ghostExchangeSources.clear();
}

/*!
	Reset the ghost information needed for data exchange for the specified rank.

	\param rank is the rank for which the information will be reset
*/
void PatchKernel::deleteGhostExchangeData(int rank)
{
	if (!isRankNeighbour(rank)) {
		return;
	}

	m_ghostExchangeTargets.erase(rank);
	m_ghostExchangeSources.erase(rank);
}

/*!
	Builds the ghost information needed for data exchange.
*/
void PatchKernel::buildGhostExchangeData()
{
	std::vector<long> ghosts;
	for (const auto &entry : m_ghostOwners) {
		long ghostId = entry.first;
		ghosts.push_back(ghostId);
	}

	deleteGhostExchangeData();
	addGhostsToExchangeTargets(ghosts);
}

/*!
	Builds the ghost information needed for data exchange for the specified
	rank.

	\param rank is the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeData(int rank)
{
	buildGhostExchangeData(std::vector<int>{rank});
}

/*!
	Builds the ghost information needed for data exchange for the specified
	list of ranks.

	\param ranks are the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeData(const std::vector<int> &ranks)
{
	// List of ghost to add
	std::unordered_set<int> buildRanks(ranks.begin(), ranks.end());

	std::vector<long> ghosts;
	for (const auto &entry : m_ghostOwners) {
		int ghostRank = entry.second;
		if (buildRanks.count(ghostRank) == 0) {
			continue;
		}

		long ghostId = entry.first;
		ghosts.push_back(ghostId);
	}

	// Build exchange data
	for (const int rank : ranks) {
		deleteGhostExchangeData(rank);
	}
	addGhostsToExchangeTargets(ghosts);
}

/*!
	Adds the specified ghosts to the exchange targets.

	No check will be perfomed to ensure that

	\param ghostIds are the ids of the ghosts that will be added
*/
void PatchKernel::addGhostsToExchangeTargets(const std::vector<long> &ghostIds)
{
	// Add the ghost to the targets
	std::unordered_set<int> ranks;
	for (const long ghostId : ghostIds) {
		// Rank of the ghost
		int rank = m_ghostOwners[ghostId];
		ranks.insert(rank);

		// Add the ghost to the targets
		m_ghostExchangeTargets[rank].push_back(ghostId);
	}

	// Sort the targets
	for (const int rank : ranks) {
		std::vector<long> &rankTargets = m_ghostExchangeTargets[rank];
		rankTargets.erase(std::unique(rankTargets.begin(), rankTargets.end()), rankTargets.end());
		std::sort(rankTargets.begin(), rankTargets.end(), CellPositionLess(*this));
	}

	// Add the sources
	addExchangeSources(ghostIds);
}

/*!
	Adds the specified ghost to the exchange list.

	\param ghostId is the id of the ghost that will be added
*/
void PatchKernel::addGhostToExchangeTargets(const long ghostId)
{
	addGhostsToExchangeTargets(std::vector<long>{ghostId});
}

/*!
	Removes the specified ghosts from the exchange list.

	\param ghostIds are the ids of the ghosts that will be removed
*/
void PatchKernel::removeGhostsFromExchangeTargets(const std::vector<long> &ghostIds)
{
	// Remove ghost from targets
	std::unordered_set<int> ranks;
	for (const long ghostId : ghostIds) {
		// Rank of the ghost
		int rank = m_ghostOwners[ghostId];
		ranks.insert(rank);

		// Remove targets
		std::vector<long> &ghostTargets = m_ghostExchangeTargets[rank];
		auto iterator = std::lower_bound(ghostTargets.begin(), ghostTargets.end(), rank, CellPositionLess(*this));
		ghostTargets.erase(iterator);
	}

	// Rebuild information of the sources
	for (const int rank : ranks) {
		m_ghostExchangeSources[rank].clear();
		addExchangeSources(m_ghostExchangeTargets[rank]);
	}
}

/*!
	Removes the specified ghost from the exchange list.

	\param ghostId id the id of the ghost that will be removed
*/
void PatchKernel::removeGhostFromExchangeTargets(const long ghostId)
{
	removeGhostsFromExchangeTargets(std::vector<long>{ghostId});
}

/*!
	Finds the internal cells that will be sources for the neighbour processors
	that owns the specified ghost cells and add those cells to the sources
	for that processor.

	\param ghostIds are the ids of the ghosts
*/
void PatchKernel::addExchangeSources(const std::vector<long> &ghostIds)
{
	// Get the sources
	std::unordered_map<int, std::unordered_set<long>> ghostSources;
	for (long ghostId : ghostIds) {
		// Owner of the ghost
		int rank = m_ghostOwners[ghostId];

		// The internal neighbourss will be sources for the rank
		for (long neighId : findCellNeighs(ghostId)) {
			if (m_ghostOwners.count(neighId) > 0) {
				continue;
			}

			ghostSources[rank].insert(neighId);
		}
	}

	// Add the sources
	for (auto entry : ghostSources) {
		int rank = entry.first;
		std::unordered_set<long> &updatedSources = entry.second;

		std::vector<long> &rankSources = m_ghostExchangeSources[rank];
		for (long rankSourceId : rankSources) {
			updatedSources.insert(rankSourceId);
		}
		rankSources = std::vector<long>(updatedSources.begin(), updatedSources.end());
		std::sort(rankSources.begin(), rankSources.end(), CellPositionLess(*this));
	}
}

/*!
    Sends the specified list of cells from process with rank sendRank (sender)
    to process with rank recvRank (receiver). If the rank the process currently
    hosting the mesh is neither the sender or the receiver, a notification is
    received in case ghost cells has changed owner.

    \param[in] sendRank sender rank
    \param[in] recvRank receiver rank
    \param[in] cellsToSend list of cells to be moved
 */
adaption::Info PatchKernel::sendCells(const int &sendRank, const int &recvRank, const std::vector<long> &cellsToSend)
{
	adaption::Info adaptionInfo;
	if (m_rank == sendRank) {
		adaptionInfo = sendCells_sender(recvRank, cellsToSend);
	} else if (m_rank == recvRank) {
		adaptionInfo = sendCells_receiver(sendRank);
	} else {
		adaptionInfo = sendCells_notified(sendRank, recvRank);
	}

	return adaptionInfo;
}

/*!
    Sends the given list of cells to the process with the specified rank.

    \param[in] recvRank is the receiver rank
    \param[in] cellsToSend is the list of cells to be sent
 */
adaption::Info PatchKernel::sendCells_sender(const int &recvRank, const std::vector<long> &cellsToSend)
{
    //
    // Initialize adaption info
    //
    adaption::Info adaptionInfo;
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_SEND;
    adaptionInfo.rank   = recvRank;

    //
    // Create a set with the cells to communicate
    //
    // For now the set will contain the cells explicitly marked for sending.
    // Later, a layer of surrounding cells will be added.
    //
    // Only internal cells can sent.
    std::vector<long> cellsToCommunicate;
    std::unordered_map<long, int> cellRankOnReceiver;

    cellsToCommunicate.reserve(cellsToSend.size());
    cellRankOnReceiver.reserve(cellsToSend.size());
    for (long cellId : cellsToSend) {
		const Cell &cell = m_cells[cellId];
		if (!cell.isInterior()) {
			throw std::runtime_error ("Only internal cells can sent.");
		}

		cellsToCommunicate.push_back(cellId);
        cellRankOnReceiver.insert({{cellId, recvRank}});
    }

    //
    // Create the notifications of ownership change
    //
    // Processors that have, among their ghost, cells that will be sent to
    // the receiver need to be notified about the ownership change. We will
    // communicate to the neighbours only the index in the exchange data
    // structure of the cells that have change ownership.
    //
    // We need to create the notification now that the set with the cells
    // to communicate contains only the cells explicitly marked for sending.
    std::unordered_map<int, std::vector<long>> ownershipNotifications;
    for (const auto &rankExchangeData : getGhostExchangeSources()) {
        int neighRank = rankExchangeData.first;
		if (neighRank == recvRank) {
			continue;
		}

        std::vector<long> &notificationList = ownershipNotifications[neighRank];

        auto &rankExchangeSources = rankExchangeData.second;
        int nRankExchangeSources = rankExchangeSources.size();
        for (long k = 0; k < nRankExchangeSources; ++k) {
            long cellId = rankExchangeSources[k];
			if (cellRankOnReceiver.count(cellId) == 0) {
                continue;
            }

            notificationList.push_back(k);
        }
	}

    // Find the frame of the cells explicitly marked for sending. These cells
    // are the cells marked for sending that have at least one neighbour that
    // will not be sent.
    std::unordered_set<long> cellsToSendFrame;
	for (long cellId : cellsToSend) {
		auto neighs = findCellNeighs(cellId);
		int nNeighs = neighs.size();
		for (int j = 0; j < nNeighs; ++j) {
			long neighId = neighs[j];
			if (cellRankOnReceiver.count(neighId) == 0) {
				cellsToSendFrame.insert(cellId);
				break;
			}
		}
	}

    // Along with the cells explicitly marked for sending, we need to send
    // also an halo of surrounfing cells. Those cells will be used by the
    // receiver to connect the cells it receives to the existing cells,
    // plus some of them will become ghost cells.
    //
    // Some cells on the halo may be already on the receiver (becuase they
    // are already ghosts owned by another processor). However we don't have
    // enough information to identify those duplicate cells. The receiver
    // needs to permorm a check to avoid inserting duplicate cells.
    //
    // Cells owned by receiver are already on the receiver, so there is no
    // need to send them.
    for (long cellId : cellsToSendFrame) {
        auto neighs = findCellNeighs(cellId);
        int nNeighs = neighs.size();
        for (int j = 0; j < nNeighs; ++j) {
            long neighId = neighs[j];
            if (cellRankOnReceiver.count(neighId) > 0) {
                continue;
            }

            int ownerRank;
            if (m_ghostOwners.count(neighId) == 0) {
                ownerRank = m_rank;
            } else {
                ownerRank = m_ghostOwners[neighId];
                if (ownerRank == recvRank) {
                    continue;
                }
            }

			cellsToCommunicate.push_back(neighId);
            cellRankOnReceiver.insert({{neighId, ownerRank}});
        }
    }

    // The cells that are exchange sources for the receiver are already on
    // the receiver (these cells are ghost owned by this processor). Mark
    // is as such to help the receiver to remove duplicate cells.
    std::unordered_map<long, long> senderGhostsToPromote;
    if (getGhostExchangeSources().count(recvRank) > 0) {
        const auto &exchangeSources = getGhostExchangeSources(recvRank);
        long nExchangeSources = exchangeSources.size();
        for (long k = 0; k < nExchangeSources; ++k) {
            long cellId = exchangeSources[k];
            if (cellRankOnReceiver.count(cellId) > 0) {
                senderGhostsToPromote.insert({{cellId, k}});
            }
        }
    }

    //
    // Create the list of vertices to send
    //
    std::unordered_set<long> vertexToCommunicate;
    for (const long &cellId : cellsToCommunicate) {
        const Cell &cell = m_cells[cellId];

        int nCellVertices = cell.getVertexCount();
        for (int j = 0; j < nCellVertices; ++j) {
            long vertexId = cell.getVertex(j);
            if (vertexToCommunicate.count(vertexId) > 0) {
                continue;
            }

            vertexToCommunicate.insert(vertexId);
        }
    }

    //
    // Send vertex data
    //
    OBinaryStream vertexBuffer;
    long vertexBufferCapacity = 0;

    // Fill buffer with vertex data
    vertexBufferCapacity += sizeof(long);
    for (long vertexId : vertexToCommunicate) {
        vertexBufferCapacity += m_vertices[vertexId].getBinarySize();
    }
    vertexBuffer.setCapacity(vertexBufferCapacity);

    vertexBuffer << (long) vertexToCommunicate.size();
    for (long vertexId : vertexToCommunicate) {
        vertexBuffer << m_vertices[vertexId];
    }

    if (vertexBufferCapacity != (long) vertexBuffer.capacity()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&vertexBufferCapacity, 1, MPI_LONG, recvRank, 10, m_communicator);
    MPI_Send(vertexBuffer.rawData(), vertexBuffer.capacity(), MPI_CHAR, recvRank, 11, m_communicator);

    //
    // Send cell data
    //
    OBinaryStream cellBuffer;
    long cellBufferCapacity = 0;

	// Fill the buffer with information on the cells that will be send to the
	// receiver, but are already there because they are ghosts owned by the
	// sender. These cells will be promoted by the receiver to internal cells.
    cellBufferCapacity += sizeof(long) + 2 * senderGhostsToPromote.size() * sizeof(long);
    cellBuffer.setCapacity(cellBufferCapacity);

    cellBuffer << (long) senderGhostsToPromote.size();
    for (const auto &entry : senderGhostsToPromote) {
        // Index of the cell in the exchange data structure
		long index = entry.second;
		cellBuffer << index;

		// Id of the cell on this processor
		long cellId = entry.first;
		cellBuffer << cellId;
	}

    // Fill the buffer with cell data
    cellBufferCapacity += sizeof(long);
    for (const long &cellId : cellsToCommunicate) {
        cellBufferCapacity += sizeof(int) + m_cells[cellId].getBinarySize();
    }
    cellBuffer.setCapacity(cellBufferCapacity);

    cellBuffer << (long) cellsToCommunicate.size();
    for (const long &cellId : cellsToCommunicate) {
        // Owner of the cell
        int cellOwner = cellRankOnReceiver[cellId];
        cellBuffer << cellOwner;

        // Cell data
        cellBuffer << m_cells[cellId];
    }

    if (cellBufferCapacity != (long) cellBuffer.capacity()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&cellBufferCapacity, 1, MPI_LONG, recvRank, 20, m_communicator);
    MPI_Send(cellBuffer.rawData(), cellBuffer.capacity(), MPI_CHAR, recvRank, 21, m_communicator);

    //
    // Send ownership notifications
    //
    for (const auto &ownershipNotification : ownershipNotifications) {
        int neighRank = ownershipNotification.first;
		const auto &notificationList = ownershipNotification.second;

        // Communicate the size of the notification
        long notificationBufferCapacity = 0;
        if (notificationList.size() > 0) {
            notificationBufferCapacity += sizeof(long) + notificationList.size() * sizeof(long);
        }
        MPI_Send(&notificationBufferCapacity, 1, MPI_LONG, neighRank, 30, m_communicator);

        // Fill the buffer and send the notification
        if (notificationBufferCapacity > 0) {
            OBinaryStream notificationBuffer(notificationBufferCapacity);

			notificationBuffer << (long) notificationList.size();
            for (long index : notificationList) {
                notificationBuffer << index;
            }

            MPI_Send(notificationBuffer.rawData(), notificationBuffer.capacity(), MPI_CHAR, neighRank, 31, m_communicator);
        }
    }

    //
    // Update cells that no longer belong to this processor
    //

    // Add the sned ids int the adaption info
    //
    // The ids will be sorted by the position of the cells, this is the
    // same order that will be used on the processor that has received the
    // octants. Since the order is the same, the two processors are able
    // to exchange cell data without any additional extra communication
    // (they already know the list of cells for which data is needed and
    // the order in which these data will be sent).
    for (long cellId : cellsToSend) {
        adaptionInfo.previous.push_back(cellId);
    }

    std::sort(adaptionInfo.previous.begin(), adaptionInfo.previous.end(), CellPositionLess(*this));

    // Delete sent cells or mark them as ghosts owned by the receiver.
    for (long cellId : cellsToSend) {
		// Check if a cell has to be delete or is a ghost owned by the receiver
		//
		// A cell will become a ghost if at least one of his neighbours is
		// an internal cell. If the processors is sending all its cells
		// there will be no ghosts cell.
		bool moveToGhosts = false;
		if (cellsToSendFrame.count(cellId) == 0) {
			moveToGhosts = false;
		} else {
			auto neighs = findCellNeighs(cellId);
			int nNeighs = neighs.size();
			for (int j = 0; j < nNeighs; ++j) {
				long neighId = neighs[j];
				if (m_ghostOwners.count(neighId) == 0) {
					moveToGhosts = true;
					break;
				}
			}
		}

        // Delete the cell or mark is as a ghost owned by the receiver.
        if (moveToGhosts) {
            moveInternal2Ghost(cellId);
            setGhostOwner(cellId, recvRank, false);
        } else {
			unsetGhostOwner(cellId, false);
            deleteCell(cellId, true, true);
        }
    }

    // Delete stale ghosts
    //
    // Loop over all the ghosts and keep only the cells that have at least
    // one internal neighbour.
    std::unordered_set<int> involvedRanks;
    involvedRanks.insert(recvRank);
    for (const auto &entry : ownershipNotifications) {
        involvedRanks.insert(entry.first);
    }

    auto itr = m_ghostOwners.cbegin();
    while (itr != m_ghostOwners.cend()) {
        long ghostId = itr->first;

        bool keep   = false;
        auto neighs = findCellNeighs(ghostId);
        int nNeighs = neighs.size();
        for (int j = 0; j < nNeighs; ++j) {
			long neighId = neighs[j];
			if (m_ghostOwners.count(neighId) == 0) {
                keep = true;
                break;
            }
        }

        auto nextItr = itr;
        nextItr++;
        if (!keep) {
            long ghostOwner = itr->second;
            involvedRanks.insert(ghostOwner);

            unsetGhostOwner(ghostId, false);
            deleteCell(ghostId, true, true);
        }
        itr = nextItr;
    }

    m_cells.flush();

	// Delete orphan vertices
	deleteOrphanVertices();

	// Rebuild ghost information
	std::vector<int> involvedRankList(involvedRanks.begin(), involvedRanks.end());
	buildGhostExchangeData(involvedRankList);

	// Return adaption info
    return adaptionInfo;
}

/*!
    Recevies a list of cells from the specified processor.

    \param[in] sendRank is the rank of the processors sending the cells
 */
adaption::Info PatchKernel::sendCells_receiver(const int &sendRank)
{
    //
    // Initialize adaption info
    //
    adaption::Info adaptionInfo;
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_RECV;
    adaptionInfo.rank   = sendRank;

    //
    // Add vertices
    //

    // Receive data
    long vertexBufferCapacity;
    MPI_Recv(&vertexBufferCapacity, 1, MPI_LONG, sendRank, 10, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream vertexBuffer(vertexBufferCapacity);
    MPI_Recv(vertexBuffer.rawData(), vertexBuffer.capacity(), MPI_CHAR, sendRank, 11, m_communicator, MPI_STATUS_IGNORE);

    // Build a kd-tree with the vertices on the ghosts cells
    //
    // These are the only vertices to check for duplicates when receiving the
    // list vertices from the sender.
    //
    // The kd-tree stores the pointer to the vertices. Ifwe try to store in the
    // kd-tree the pointers to the vertices of the patch, the first resize of
    // the vertex container would invalidate the pointer. Create a copy of the
    // vertices and store the pointer to that copy.
    long nGhostsMaxVertices = 0;
    for (const auto &entry : m_ghostOwners) {
        long ghostId = entry.first;
        const Cell &ghost = m_cells[ghostId];
        nGhostsMaxVertices += ghost.getVertexCount();
    }

    KdTree<3, Vertex, long> ghostVerticesTree(nGhostsMaxVertices);
    std::unordered_map<long, Vertex> ghostVertices(nGhostsMaxVertices);
    for (const auto &entry : m_ghostOwners) {
        long ghostId = entry.first;
        const Cell &ghost = m_cells[ghostId];
        int nGhostVertices = ghost.getVertexCount();
        for (int k = 0; k < nGhostVertices; ++k) {
            long vertexId = ghost.getVertex(k);
            if (ghostVertices.count(vertexId) == 0) {
                ghostVertices.insert({{vertexId, m_vertices[vertexId]}});
                ghostVerticesTree.insert(&ghostVertices.at(vertexId), vertexId);
            }
        }
    }

    // Receive vertices
    //
    // There are no duplicate in the received vertices, but some of them may
    // be already a local vertex of a ghost cell.
    long nRecvVertices;
    vertexBuffer >> nRecvVertices;

    std::unordered_map<long, long> recvVertexMap;
    recvVertexMap.reserve(nRecvVertices);
    for (long i = 0; i < nRecvVertices; ++i) {
        Vertex vertex;
        vertexBuffer >> vertex;
        long recvVertexId = vertex.getId();

        long localVertexId;
        if (ghostVerticesTree.exist(&vertex, localVertexId) < 0) {
            localVertexId = generateVertexId();
            addVertex(std::move(vertex), localVertexId);
        }

        recvVertexMap.insert({{recvVertexId, localVertexId}});
    }

    std::unordered_map<long, Vertex>().swap(ghostVertices);

    //
    // Add cells
    //

    // Receive data
    long cellBufferCapacity;
    MPI_Recv(&cellBufferCapacity, 1, MPI_LONG, sendRank, 20, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream cellBuffer(cellBufferCapacity);
    MPI_Recv(cellBuffer.rawData(), cellBuffer.capacity(), MPI_CHAR, sendRank, 21, m_communicator, MPI_STATUS_IGNORE);

    // Ghost owned by the receiver that will be promoted as internal cell
    //
    // Among all the received cells, some cells may be already here because
    // they are ghost cells owned by the sending processor. To identify those
    // cells, the sender is sending us the position in the ghost exchange data
    // structures of the duplicate cells owned by this processor. These cells
    // will be promoted to internal cells and the data received by the sender
    // will allow to properly set the adjacencies.
    //
    // Other cells may be already here, i.e. ghost cells owned by a processor
    // that is not the sender. However the sender doens't have the necessary
    // information to tell us which ghost cells are already here. The data we
    // receive now is only related to the ghost cells owned by the sender.
    // Cells owned by other processors need to be check when they are received.
    long nSenderGhostsToPromote;
    cellBuffer >> nSenderGhostsToPromote;

	unordered_map<long, long> senderGhostsToPromote;
	senderGhostsToPromote.reserve(nSenderGhostsToPromote);

	if (nSenderGhostsToPromote > 0) {
		const auto &sendExchangeTargets = getGhostExchangeTargets(sendRank);
		for (int n = 0; n < nSenderGhostsToPromote; ++n) {
			long index;
			cellBuffer >> index;
			long localCellId = sendExchangeTargets[index];

			long senderCellId;
			cellBuffer >> senderCellId;

			senderGhostsToPromote.insert({{senderCellId, localCellId}});
		}
	}

    // Receive cells
    //
    // The sender is sending also an halo surroudning the cells explicitly
    // marked for sending. This halo will be used for connecting the received
    // cells to the existing ones and to build the ghosts.
    //
    // Some cells on the halo can be duplicate. For some of these cells we
    // already know the apping between sender cell ids and local ids. For
    // all the other cells we need to check if they are duplicate. All halo
    // are inserted in the cell list, if they are duplicates they will be
    // deleted later.
    long nRecvCells;
    cellBuffer >> nRecvCells;

    std::unordered_set<int> involvedRanks;
	involvedRanks.insert(sendRank);

    std::unordered_map<long, long> recvCellMap;
    std::unordered_set<long> recvCells;
	recvCells.reserve(nRecvCells);

    std::unordered_map<long, FlatVector2D<long>> linkAdjacencies;

    m_cells.reserve(nRecvCells);
    for (long i = 0; i < nRecvCells; ++i) {
        // Cell data
        int recvCellOwner;
        cellBuffer >> recvCellOwner;

        Cell recvCell;
        cellBuffer >> recvCell;
        long senderCellId = recvCell.getId();

        // Set cell interior flag
        bool recvIsInterior = (recvCellOwner == m_rank);
        recvCell.setInterior(recvIsInterior);

        // Remap connectivity
        int nCellVertices = recvCell.getVertexCount();
        for (int j = 0; j < nCellVertices; ++j) {
            long senderVertexId = recvCell.getVertex(j);
            long localVertexId  = recvVertexMap.at(senderVertexId);

            recvCell.setVertex(j, localVertexId);
        }

        // Check if the cells is a duplicate
        long localCellId = Cell::NULL_ID;
        if (recvCellOwner != m_rank && getGhostExchangeTargets().count(recvCellOwner) > 0) {
            const auto &rankExchangeTargets = getGhostExchangeTargets(recvCellOwner);
            for (long ghostId : rankExchangeTargets) {
                const Cell &ghostCell = m_cells[ghostId];
                if (ghostCell.getType() != recvCell.getType()) {
                    continue;
                }

                bool cellsCoincide = true;
                int nGhostVertices = ghostCell.getVertexCount();
                for (int vertex = 0; vertex < nGhostVertices; ++vertex) {
                    long ghostVertexId = ghostCell.getVertex(vertex);
                    long recvVertexId  = recvCell.getVertex(vertex);
                    if (ghostVertexId != recvVertexId) {
                        cellsCoincide = false;
                        break;
                    }
                }

                if (cellsCoincide) {
                    localCellId = ghostId;
                    break;
                }
            }
        } else if (senderGhostsToPromote.count(senderCellId) > 0) {
            localCellId = senderGhostsToPromote[senderCellId];
        }

        // If the cell is not a duplicate add it in the cell data structure,
        // otherwise merge the connectivity of the duplicate cell to the
        // existing cell. This ensure that the received cell will be
        // properly connected to the received cells
        if (localCellId < 0) {
            // Add cell
            localCellId = generateCellId();
            addCell(std::move(recvCell), localCellId);
            if (!recvIsInterior) {
                setGhostOwner(localCellId, recvCellOwner, false);
            }

            recvCells.insert(localCellId);

            // Update adaption info
            adaptionInfo.current.push_back(localCellId);
        } else {
            // Check if the existing cells needs to become an internal cell
            Cell &localCell = m_cells[localCellId];
            if (recvIsInterior && !localCell.isInterior()) {
                unsetGhostOwner(localCellId, false);
                moveGhost2Internal(localCellId);
            }

            // Save the adjacencies of the received cell, this adjacencies
            // will link together the recevied cell to the existing ones.
            FlatVector2D<long> &recvAdjacencies = linkAdjacencies[localCellId];

            int nCellFaces = recvCell.getFaceCount();
            recvAdjacencies.reserve(nCellFaces);
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = recvCell.getAdjacencyCount(face);

                std::vector<long> faceAdjacencies;
                faceAdjacencies.reserve(nFaceAdjacencies);
                for (int k = 0; k < nFaceAdjacencies; ++k) {
                    faceAdjacencies.push_back(recvCell.getAdjacency(face, k));
                }

                recvAdjacencies.pushBack(faceAdjacencies);
            }
        }

        // Ranks involved in the communication
        if (recvCellOwner != m_rank) {
            involvedRanks.insert(recvCellOwner);
        }

        // Add the cell to the cell map
        recvCellMap.insert({{senderCellId, localCellId}});
    }

    // Remap adjacencies
    for (auto cellId : recvCells) {
        Cell &cell = m_cells[cellId];

        int nCellFaces = cell.getFaceCount();
        for (int face = 0; face < nCellFaces; ++face) {
            int nFaceAdjacencies = cell.getAdjacencyCount(face);
            for (int k = 0; k < nFaceAdjacencies; ++k) {
                long senderAdjacencyId = cell.getAdjacency(face, k);
                if (senderAdjacencyId < 0) {
                    continue;
                } else if (recvCellMap.count(senderAdjacencyId) == 0) {
					cell.deleteAdjacency(face, k);
                    continue;
                } else {
					long localAdjacencyId = recvCellMap.at(senderAdjacencyId);
					cell.setAdjacency(face, k, localAdjacencyId);
				}
            }
        }
    }

    // Link received cells with the current cells
    for (auto &entry : linkAdjacencies) {
        long cellId = entry.first;
        Cell &cell = m_cells[cellId];

        int nCellFaces = cell.getFaceCount();
        FlatVector2D<long> &cellLinkAdjacencies = entry.second;
        for (int face = 0; face < nCellFaces; ++face) {
            int nFaceLinkAdjacencies = cellLinkAdjacencies.getItemCount(face);
            for (int k = 0; k < nFaceLinkAdjacencies; ++k) {
                long senderAdjacencyId = cellLinkAdjacencies.getItem(face, k);
                long localAdjacencyId  = recvCellMap[senderAdjacencyId];
                if (cell.findAdjacency(face, localAdjacencyId) >= 0) {
                    continue;
                }

                cell.pushAdjacency(face, localAdjacencyId);
            }
        }
    }

    // Sort the ids in the adaption info
    //
    // The ids will be sorted by the position of the cells, this is the
    // same order that will be used on the processor that has received the
    // octants. Since the order is the same, the two processors are able
    // to exchange cell data without any additional extra communication
    // (they already know the list of cells for which data is needed and
    // the order in which these data will be sent).
    std::sort(adaptionInfo.current.begin(), adaptionInfo.current.end(), CellPositionLess(*this));

    // Rebuild ghost information
    std::vector<int> involvedRankList(involvedRanks.begin(), involvedRanks.end());
    buildGhostExchangeData(involvedRankList);

    // Return adaption info
    return adaptionInfo;
}

/*!
    Notifies the current processor of changes in ghost ownership after a
    cell send operation.

    \param[in] sendRank is the rank of the processor sending the cells
    \param[in] recvRank is the rank of the processor receiving the cells
 */
adaption::Info PatchKernel::sendCells_notified(const int &sendRank, const int &recvRank)
{
    adaption::Info adaptionInfo;

    // This processor will receive a notification only if it has ghosts owned
    // by the sender
    if (getGhostExchangeTargets().count(sendRank) == 0) {
        adaptionInfo.type = adaption::TYPE_NONE;
        return adaptionInfo;
    }

    // Receive ownership changes
    long bufferCapacity;
    MPI_Recv(&bufferCapacity, 1, MPI_LONG, sendRank, 30, m_communicator, MPI_STATUS_IGNORE);
    if (bufferCapacity == 0) {
        adaptionInfo.type = adaption::TYPE_NONE;
        return adaptionInfo;
    }

    IBinaryStream buffer(bufferCapacity);
    MPI_Recv(buffer.rawData(), buffer.capacity(), MPI_CHAR, sendRank, 31, m_communicator, MPI_STATUS_IGNORE);

    // Initialize adaption info
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_NOTICE;
    adaptionInfo.rank   = sendRank;

    // Apply ownership changes
    long nChanges;
    buffer >> nChanges;

    const auto exchangeTargets = getGhostExchangeTargets(sendRank);
    for (long k = 0; k < nChanges; ++k) {
        long index;
        buffer >> index;

        long cellId = exchangeTargets[index];
        setGhostOwner(cellId, recvRank, false);
    }

    // Rebuild ghost exchagne data
    buildGhostExchangeData(std::vector<int>{sendRank, recvRank});

	// Return adaption info
    return adaptionInfo;
}

}

#endif
