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
	Sets the size, expressed in number of layers, of the ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::setHaloSize(std::size_t haloSize)
{
	if (isPartitioned()) {
		throw std::runtime_error ("Halo size can only be set before partitionig the patch.");
	}

	std::size_t maxHaloSize = _getMaxHaloSize();
	if (haloSize > maxHaloSize) {
		throw std::runtime_error ("Halo size exceeds the maximum allowed value.");
	}

	m_haloSize = haloSize;

	_setHaloSize(haloSize);
}

/*!
	Gets the size, expressed in number of layers, of the ghost cells halo.

	\result The size, expressed in number of layers, of the ghost cells halo.
*/
std::size_t PatchKernel::getHaloSize() const
{
	return m_haloSize;
}

/*!
	Gets the maximum allowed size, expressed in number of layers, of the ghost
	cells halo.

	\result The maximum allowed size, expressed in number of layers, of the
	ghost cells halo.
*/
std::size_t PatchKernel::_getMaxHaloSize()
{
	return 1;
}

/*!
	Internal function to set the size, expressed in number of layers, of the
	ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::_setHaloSize(std::size_t haloSize)
{
	BITPIT_UNUSED(haloSize);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(cellRanks, trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning.
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	partitioningPrepare(cellRanks, false);

	partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	partitioningPrepare(false);

	partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(cellRanks, trackPartitioning);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning.
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(const std::vector<int> &cellRanks, bool trackPartitioning)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Build the list of ids to be sent
	auto cellItr = cellBegin();
	for (int k = 0; k < getInternalCount(); ++k) {
		const int &rank = cellRanks[k];
		if (rank == getRank()) {
			cellItr++;
			continue;
		}

		m_partitioningLocalSendList[rank].push_back(cellItr->getId());

		cellItr++;
	}

	// Local senders and receivers
	int nLocalExchanges = m_partitioningLocalSendList.size();

	std::vector<int> localSenders;
	std::vector<int> localReceivers;
	localSenders.reserve(nLocalExchanges);
	localReceivers.reserve(nLocalExchanges);
	for (const auto &entry : m_partitioningLocalSendList) {
		localSenders.push_back(getRank());
		localReceivers.push_back(entry.first);
	}

	// Prepare the communication for exchanging the sender/receiver pairs
	std::vector<int> globalExchangeSizes(getProcessorCount());
	MPI_Allgather(&nLocalExchanges, 1, MPI_INT, globalExchangeSizes.data(), 1, MPI_INT, getCommunicator());

	std::vector<int> globalExchangeOffsets(getProcessorCount());
	globalExchangeOffsets[0] = 0;
	for (int i = 1; i < getProcessorCount(); ++i) {
		globalExchangeOffsets[i] = globalExchangeOffsets[i-1] + globalExchangeSizes[i-1];
	}

	// Gather global information
	m_nPartitioningGlobalExchanges = globalExchangeOffsets.back() + globalExchangeSizes.back();

	m_partitioningGlobalSenders.resize(m_nPartitioningGlobalExchanges);
	m_partitioningGlobalReceivers.resize(m_nPartitioningGlobalExchanges);

	MPI_Allgatherv(localSenders.data(), localSenders.size(), MPI_INT, m_partitioningGlobalSenders.data(),
				   globalExchangeSizes.data(), globalExchangeOffsets.data(), MPI_INT, getCommunicator());

	MPI_Allgatherv(localReceivers.data(), localReceivers.size(), MPI_INT, m_partitioningGlobalReceivers.data(),
				   globalExchangeSizes.data(), globalExchangeOffsets.data(), MPI_INT, getCommunicator());

	// Build the information on the cells that will be sent
	if (trackPartitioning) {
		for (const auto &entry : m_partitioningLocalSendList) {
			int receiver = entry.first;
			const std::vector<long> &ids = entry.second;

			partitioningData.emplace_back();
			adaption::Info &partitioningInfo = partitioningData.back();
			partitioningInfo.entity   = adaption::ENTITY_CELL;
			partitioningInfo.type     = adaption::TYPE_PARTITION_SEND;
			partitioningInfo.rank     = receiver;
			partitioningInfo.previous = ids;
		}
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(trackPartitioning);
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(bool trackPartitioning)
{
	std::vector<adaption::Info> partitioningData;

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Reset partitioning information
	m_nPartitioningGlobalExchanges = 0;

	// Execute the partitioning preparation
	partitioningData = _partitioningPrepare(trackPartitioning);

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Alter the patch performing the partitioning.

	The actual modification of the patch takes place during this phase. After
	this phase the adapton is completed and the patch is in its final state.
	Optionally the patch can track the changes performed to the patch.

	\param trackPartitioning if set to true the function will return the changes
	done to the patch during the partitioning
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the partitioning
	\result If the partitioning is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the partitioning, otherwise
	an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_UNSUPPORTED || partitioningStatus == PARTITIONING_CLEAN) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_PREPARED) {
		throw std::runtime_error ("The prepare function has no been called.");
	}

	// Begin patch alteration
	beginAlteration();

	// Alter patch
	if (m_nPartitioningGlobalExchanges == 0) {
		partitioningData = _partitioningAlter(trackPartitioning);
	} else {
		std::vector<long> emptyCellList;
		for (int i = 0; i < m_nPartitioningGlobalExchanges; ++i) {
			int sender   = m_partitioningGlobalSenders[i];
			int receiver = m_partitioningGlobalReceivers[i];

			std::vector<long> *ids;
			if (sender == getRank()) {
				ids = &(m_partitioningLocalSendList[receiver]);
			} else {
				ids = &emptyCellList;
			}

			adaption::Info partitioningInfo = sendCells_any(sender, receiver, *ids);
			if (trackPartitioning && partitioningInfo.type != adaption::TYPE_NONE) {
				partitioningData.push_back(std::move(partitioningInfo));
			}
		}
	}

	// End patch alteration
	endAlteration(squeezeStorage);

	// The patch is now partitioned
	setPartitioned(true);

	// Update the status
	setPartitioningStatus(PARTITIONING_ALTERED);

	return partitioningData;
}

/*!
	Cleanup patch data structured after the partitioning.

	The patch will only clean-up the data structures needed during the
	partitioning.
*/
void PatchKernel::partitioningCleanup()
{
	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_UNSUPPORTED || partitioningStatus == PARTITIONING_CLEAN) {
		return;
	} else if (partitioningStatus == PARTITIONING_PREPARED) {
		throw std::runtime_error ("It is not yet possible to abort a partitioning.");
	} else if (partitioningStatus != PARTITIONING_ALTERED) {
		throw std::runtime_error ("The alter function has no been called.");
	}

	// Clean-up the partitioning
	_partitioningCleanup();

	if (m_nPartitioningGlobalExchanges != 0) {
		m_nPartitioningGlobalExchanges = 0;
		std::unordered_map<int, std::vector<long>>().swap(m_partitioningLocalSendList);
		std::vector<int>().swap(m_partitioningGlobalSenders);
		std::vector<int>().swap(m_partitioningGlobalReceivers);
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_CLEAN);
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
	Returns the current partitioning status.

	\param global if set to true the partitioning status will be
	\return The current partitioning status.
*/
PatchKernel::PartitioningStatus PatchKernel::getPartitioningStatus(bool global) const
{
	int partitioningStatus = static_cast<int>(m_partitioningStatus);

	if (global && isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &partitioningStatus, 1, MPI_INT, MPI_MAX, communicator);
	}

	return static_cast<PartitioningStatus>(partitioningStatus);
}

/*!
	Set the current partitioning status.

	\param status is the partitioning status that will be set
*/
void PatchKernel::setPartitioningStatus(PartitioningStatus status)
{
	m_partitioningStatus = status;
}


/*!
	Evaluate partitioning load unbalance index.

	\result Partitioning load unbalance index.
*/
double PatchKernel::evalPartitioningUnbalance()
{
	if (!isPartitioned()) {
		return 0.;
	}

	// Evaluate partition weight
	double localWeight = getInternalCount();

	// Evalaute global weights
	double totalWeight;
	MPI_Allreduce(&localWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, getCommunicator());

	double minimumWeight;
	MPI_Allreduce(&localWeight, &minimumWeight, 1, MPI_DOUBLE, MPI_MIN, getCommunicator());

	double maximumWeight;
	MPI_Allreduce(&localWeight, &maximumWeight, 1, MPI_DOUBLE, MPI_MAX, getCommunicator());

	// Evaluate the unbalance
	double unbalance = (maximumWeight - minimumWeight) / totalWeight;

	return unbalance;
}

/*!
	Prepares the patch for performing the partitioning.

	Default implementation is a no-op function.

	\param trackPartitioning if set to true the function will return the
	changes that will be performed in the alter step
	\result If the partitioning is tracked, returns a vector of adaption::Info
	that can be used to discover what changes will be performed in the alter
	step, otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningPrepare(bool trackPartitioning)
{
	BITPIT_UNUSED(trackPartitioning);

	return std::vector<adaption::Info>();
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
std::vector<adaption::Info> PatchKernel::_partitioningAlter(bool trackPartitioning)
{
	BITPIT_UNUSED(trackPartitioning);

	assert(false && "The patch needs to implement _partitioningAlter");

	return std::vector<adaption::Info>();
}

/*!
	Cleanup patch data structured after the partitioning.

	Default implementation is a no-op function.
*/
void PatchKernel::_partitioningCleanup()
{
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
	Gets the halo layer of the specified cell.

	\param id is the id of the requested cell
	\result The halo layer of the specified cell.
*/
int PatchKernel::getCellHaloLayer(const long &id) const
{
	const Cell &cell = getCell(id);
	if (cell.isInterior()) {
		return 0;
	} else {
		return -1;
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
	Gets a constant reference to the ghost targets needed for data exchange.

	\result A constant reference to the ghost targets needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeTargets() const
{
	return m_ghostExchangeTargets;
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
	Gets a constant reference to the ghost sources needed for data exchange.

	\result A constant reference to the ghost sources needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeSources() const
{
	return m_ghostExchangeSources;
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
	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::setGhostOwner(int id, int rank, bool updateExchangeInfo)
{
	// Rebuild the exchange info information of the previous owner
	if (updateExchangeInfo) {
		if (m_ghostOwners.count(id) > 0) {
			removeGhostFromExchangeInfo(id);
		}
	}

	// Assign the owner to the cell
	m_ghostOwners[id] = rank;

	// Rebuild the exchange info information of the current owner
	if (updateExchangeInfo) {
		addGhostToExchangeInfo(id);
	}
}

/*!
	Unsets the owner of the specified ghost.

	\param id is the id of the ghost cell
	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::unsetGhostOwner(int id, bool updateExchangeInfo)
{
	if (m_ghostOwners.count(id) <= 0) {
		return;
	}

	// Rebuild the exchange info information of the previous owner
	if (updateExchangeInfo) {
		removeGhostFromExchangeInfo(id);
	}

	// Remove the owner
	m_ghostOwners.erase(id);
}

/*!
	Clear the owners of all the ghosts.

	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::clearGhostOwners(bool updateExchangeInfo)
{
	// Clear the owners
	m_ghostOwners.clear();

	// Clear exchange info
	if (updateExchangeInfo) {
		deleteGhostExchangeInfo();
	}
}

/*!
	Reset the ghost information needed for data exchange.
*/
void PatchKernel::deleteGhostExchangeInfo()
{
	m_ghostExchangeTargets.clear();
	m_ghostExchangeSources.clear();
}

/*!
	Reset the ghost information needed for data exchange for the specified rank.

	\param rank is the rank for which the information will be reset
*/
void PatchKernel::deleteGhostExchangeInfo(int rank)
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
void PatchKernel::buildGhostExchangeInfo()
{
	std::vector<long> ghosts;
	for (const auto &entry : m_ghostOwners) {
		long ghostId = entry.first;
		ghosts.push_back(ghostId);
	}

	deleteGhostExchangeInfo();
	addGhostsToExchangeInfo(ghosts);
}

/*!
	Builds the ghost information needed for data exchange for the specified
	rank.

	\param rank is the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeInfo(int rank)
{
	buildGhostExchangeInfo(std::vector<int>{rank});
}

/*!
	Builds the ghost information needed for data exchange for the specified
	list of ranks.

	\param ranks are the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeInfo(const std::vector<int> &ranks)
{
	// Check if all structures needed are ready
	assert(getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

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

	// Build exchange info
	for (const int rank : ranks) {
		deleteGhostExchangeInfo(rank);
	}
	addGhostsToExchangeInfo(ghosts);
}

/*!
	Adds the specified ghosts to the exchange targets.

	No check will be perfomed to ensure that

	\param ghostIds are the ids of the ghosts that will be added
*/
void PatchKernel::addGhostsToExchangeInfo(const std::vector<long> &ghostIds)
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

	// Build the sources
	for (const int rank : ranks) {
		buildGhostExchangeSources(rank);
	}
}

/*!
	Adds the specified ghost to the exchange list.

	\param ghostId is the id of the ghost that will be added
*/
void PatchKernel::addGhostToExchangeInfo(const long ghostId)
{
	addGhostsToExchangeInfo(std::vector<long>{ghostId});
}

/*!
	Removes the specified ghosts from the exchange list.

	\param ghostIds are the ids of the ghosts that will be removed
*/
void PatchKernel::removeGhostsFromExchangeInfo(const std::vector<long> &ghostIds)
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
		assert(iterator != ghostTargets.end());
		ghostTargets.erase(iterator);
	}

	// Rebuild information of the sources
	for (const int rank : ranks) {
		buildGhostExchangeSources(rank);
	}
}

/*!
	Removes the specified ghost from the exchange list.

	\param ghostId id the id of the ghost that will be removed
*/
void PatchKernel::removeGhostFromExchangeInfo(const long ghostId)
{
	removeGhostsFromExchangeInfo(std::vector<long>{ghostId});
}

/*!
	Build the information about the internal cells that will be ghost cells
	for the processor with the specified rank.

	\param rank is the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeSources(int rank)
{
	buildGhostExchangeSources(std::vector<int>{rank});
}

/*!
	Build the information about the internal cells that will be ghost cells
	for the processors with the specified ranks.

	\param ranks are the rank for which the information will be built
*/
void PatchKernel::buildGhostExchangeSources(const std::vector<int> &ranks)
{
	for (int rank : ranks) {
		std::vector<long> &rankSources = m_ghostExchangeSources[rank];

		// Clear current sources
		rankSources.clear();

		// Update the source list
		rankSources = _findGhostExchangeSources(rank);

		// Sort the sources
		std::sort(rankSources.begin(), rankSources.end(), CellPositionLess(*this));
	}
}

/*!
	Finds the internal cells that will be ghost cells for the processors
	with the specified ranks. During data exchange, these cells will be
	the sources form which data will be read from.

	\param rank is the rank for which the information will be built
*/
std::vector<long> PatchKernel::_findGhostExchangeSources(int rank)
{
	// The internal neighbours of the ghosts will be sources for the rank
	std::vector<long> neighIds;
	std::unordered_set<long> exchangeSources;
	exchangeSources.reserve(m_ghostExchangeTargets[rank].size());
	for (long ghostId : m_ghostExchangeTargets[rank]) {
		neighIds.clear();
		findCellNeighs(ghostId, &neighIds);
		for (long neighId : neighIds) {
			if (m_ghostOwners.count(neighId) > 0) {
				continue;
			}

			exchangeSources.insert(neighId);
		}
	}

	return std::vector<long>(exchangeSources.begin(), exchangeSources.end());
}

/*!
    Sends the specified list of cells from process with rank sendRank (sender)
    to process with rank recvRank (receiver). If the rank the process currently
    hosting the mesh is neither the sender or the receiver, a notification is
    received in case ghost cells has changed owner.

    \param[in] sendRank sender rank
    \param[in] recvRank receiver rank
    \param[in] cellsToSend list of cells to be moved
    \param[in] squeezeStorage if set to true the vector that store patch information
    will be squeezed after the synchronization
 */
adaption::Info PatchKernel::sendCells(const int &sendRank, const int &recvRank,
                                      const std::vector<long> &cellsToSend,
                                      bool squeezeStorage)
{
	//
	// Pereare partitioning alteration
	//
	partitioningPrepare(false);

	//
	// Alter partitioning
	//

	// Begin patch alteration
	beginAlteration();

	// Send cells
	adaption::Info adaptionInfo = sendCells_any(sendRank, recvRank, cellsToSend);

	// End patch alteration
	endAlteration(squeezeStorage);

	// The patch is now partitioned
	setPartitioned(true);

	// Update the status
	setPartitioningStatus(PARTITIONING_ALTERED);

	//
	// Cleanup partitioning alteration
	//
	partitioningCleanup();

	return adaptionInfo;
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
adaption::Info PatchKernel::sendCells_any(const int &sendRank, const int &recvRank,
                                          const std::vector<long> &cellsToSend)
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
    std::vector<long> neighIds;

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
    // communicate to the neighbours only the index in the exchange info
    // structure of the cells that have change ownership.
    //
    // We need to create the notification now that the set with the cells
    // to communicate contains only the cells explicitly marked for sending.
    std::unordered_map<int, std::vector<long>> ownershipNotifications;
    for (const auto &rankExchangeInfo : getGhostExchangeSources()) {
        int neighRank = rankExchangeInfo.first;
		if (neighRank == recvRank) {
			continue;
		}

        std::vector<long> &notificationList = ownershipNotifications[neighRank];

        auto &rankExchangeSources = rankExchangeInfo.second;
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
		neighIds.clear();
		findCellNeighs(cellId, &neighIds);
		for (long neighId : neighIds) {
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
        neighIds.clear();
        findCellNeighs(cellId, &neighIds);
        for (long neighId : neighIds) {
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

        ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
        int nCellVertices = cellVertexIds.size();
        for (int j = 0; j < nCellVertices; ++j) {
            long vertexId = cellVertexIds[j];
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
    long vertexBufferSize = 0;

    // Fill buffer with vertex data
    vertexBufferSize += sizeof(long);
    for (long vertexId : vertexToCommunicate) {
        vertexBufferSize += m_vertices[vertexId].getBinarySize();
    }
    vertexBuffer.setSize(vertexBufferSize);

    vertexBuffer << (long) vertexToCommunicate.size();
    for (long vertexId : vertexToCommunicate) {
        vertexBuffer << m_vertices[vertexId];
    }

    if (vertexBufferSize != (long) vertexBuffer.getSize()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&vertexBufferSize, 1, MPI_LONG, recvRank, 10, m_communicator);
    MPI_Send(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, recvRank, 11, m_communicator);

    //
    // Send cell data
    //
    OBinaryStream cellBuffer;
    long cellBufferSize = 0;

	// Fill the buffer with information on the cells that will be send to the
	// receiver, but are already there because they are ghosts owned by the
	// sender. These cells will be promoted by the receiver to internal cells.
    cellBufferSize += sizeof(long) + 2 * senderGhostsToPromote.size() * sizeof(long);
    cellBuffer.setSize(cellBufferSize);

    cellBuffer << (long) senderGhostsToPromote.size();
    for (const auto &entry : senderGhostsToPromote) {
        // Index of the cell in the exchange info structure
		long index = entry.second;
		cellBuffer << index;

		// Id of the cell on this processor
		long cellId = entry.first;
		cellBuffer << cellId;
	}

    // Fill the buffer with cell data
    cellBufferSize += sizeof(long);
    for (const long &cellId : cellsToCommunicate) {
        cellBufferSize += sizeof(int) + m_cells[cellId].getBinarySize();
    }
    cellBuffer.setSize(cellBufferSize);

    cellBuffer << (long) cellsToCommunicate.size();
    for (const long &cellId : cellsToCommunicate) {
        // Owner of the cell
        int cellOwner = cellRankOnReceiver[cellId];
        cellBuffer << cellOwner;

        // Cell data
        cellBuffer << m_cells[cellId];
    }

    if (cellBufferSize != (long) cellBuffer.getSize()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&cellBufferSize, 1, MPI_LONG, recvRank, 20, m_communicator);
    MPI_Send(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, recvRank, 21, m_communicator);

    //
    // Send ownership notifications
    //
    for (const auto &ownershipNotification : ownershipNotifications) {
        int neighRank = ownershipNotification.first;
		const auto &notificationList = ownershipNotification.second;

        // Communicate the size of the notification
        long notificationBufferSize = 0;
        if (notificationList.size() > 0) {
            notificationBufferSize += sizeof(long) + notificationList.size() * sizeof(long);
        }
        MPI_Send(&notificationBufferSize, 1, MPI_LONG, neighRank, 30, m_communicator);

        // Fill the buffer and send the notification
        if (notificationBufferSize > 0) {
            OBinaryStream notificationBuffer(notificationBufferSize);

			notificationBuffer << (long) notificationList.size();
            for (long index : notificationList) {
                notificationBuffer << index;
            }

            MPI_Send(notificationBuffer.data(), notificationBuffer.getSize(), MPI_CHAR, neighRank, 31, m_communicator);
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
			neighIds.clear();
			findCellNeighs(cellId, &neighIds);
			for (long neighId : neighIds) {
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

        neighIds.clear();
        findCellNeighs(ghostId, &neighIds);
        bool keep = false;
        for (long neighId : neighIds) {
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
	buildGhostExchangeInfo(involvedRankList);

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
    long vertexBufferSize;
    MPI_Recv(&vertexBufferSize, 1, MPI_LONG, sendRank, 10, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream vertexBuffer(vertexBufferSize);
    MPI_Recv(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, sendRank, 11, m_communicator, MPI_STATUS_IGNORE);

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
        ConstProxyVector<long> ghostVertexIds = ghost.getVertexIds();
        int nGhostVertices = ghostVertexIds.size();
        for (int k = 0; k < nGhostVertices; ++k) {
            long vertexId = ghostVertexIds[k];
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

    std::unordered_map<long, long> vertexMap;
    vertexMap.reserve(nRecvVertices);
    for (long i = 0; i < nRecvVertices; ++i) {
        Vertex vertex;
        vertexBuffer >> vertex;
        long vertexId = vertex.getId();

        long localVertexId;
        if (ghostVerticesTree.exist(&vertex, localVertexId) < 0) {
            localVertexId = generateVertexId();
            addVertex(std::move(vertex), localVertexId);
        }

        vertexMap.insert({{vertexId, localVertexId}});
    }

    std::unordered_map<long, Vertex>().swap(ghostVertices);

    //
    // Add cells
    //

    // Receive data
    long cellBufferSize;
    MPI_Recv(&cellBufferSize, 1, MPI_LONG, sendRank, 20, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream cellBuffer(cellBufferSize);
    MPI_Recv(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, sendRank, 21, m_communicator, MPI_STATUS_IGNORE);

    // Ghost owned by the receiver that will be promoted as internal cell
    //
    // Among all the received cells, some cells may be already here because
    // they are ghost cells owned by the sending processor. To identify those
    // cells, the sender is sending us the position in the ghost exchange info
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
			long cellId = sendExchangeTargets[index];

			long senderCellId;
			cellBuffer >> senderCellId;

			senderGhostsToPromote.insert({{senderCellId, cellId}});
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
    long nReceivedCells;
    cellBuffer >> nReceivedCells;

    std::unordered_set<int> involvedRanks;
	involvedRanks.insert(sendRank);

    std::unordered_map<long, long> cellMap;
    std::unordered_set<long> receivedCells;
	receivedCells.reserve(nReceivedCells);

    std::unordered_map<long, FlatVector2D<long>> linkAdjacencies;

    m_cells.reserve(nReceivedCells);
    for (long i = 0; i < nReceivedCells; ++i) {
        // Cell data
        int cellOwner;
        cellBuffer >> cellOwner;

        Cell cell;
        cellBuffer >> cell;
        long cellOriginalId = cell.getId();
        ConstProxyVector<long> vertexIds = cell.getVertexIds();

        // Set cell interior flag
        bool isInterior = (cellOwner == m_rank);
        cell.setInterior(isInterior);

        // Remap connectivity
        cell.renumberVertices(vertexMap);

        // Check if the cells is a duplicate
        long cellId = Cell::NULL_ID;
        if (cellOwner != m_rank && getGhostExchangeTargets().count(cellOwner) > 0) {
            const auto &rankExchangeTargets = getGhostExchangeTargets(cellOwner);
            for (long ghostId : rankExchangeTargets) {
                const Cell &ghostCell = m_cells[ghostId];
                if (ghostCell.getType() != cell.getType()) {
                    continue;
                }

                bool cellsCoincide = true;
                ConstProxyVector<long> ghostVertexIds = ghostCell.getVertexIds();
                int nGhostVertices = ghostVertexIds.size();
                for (int vertex = 0; vertex < nGhostVertices; ++vertex) {
                    long ghostVertexId = ghostVertexIds[vertex];
                    long vertexId  = vertexIds[vertex];
                    if (ghostVertexId != vertexId) {
                        cellsCoincide = false;
                        break;
                    }
                }

                if (cellsCoincide) {
                    cellId = ghostId;
                    break;
                }
            }
        } else if (senderGhostsToPromote.count(cellOriginalId) > 0) {
            cellId = senderGhostsToPromote[cellOriginalId];
        }

        // If the cell is not a duplicate add it in the cell data structure,
        // otherwise merge the connectivity of the duplicate cell to the
        // existing cell. This ensure that the received cell will be
        // properly connected to the received cells
        if (cellId < 0) {
            // Add cell
            cellId = generateCellId();
            addCell(std::move(cell), cellId);
            if (!isInterior) {
                setGhostOwner(cellId, cellOwner, false);
            }

            receivedCells.insert(cellId);

            // Update adaption info
            adaptionInfo.current.push_back(cellId);
        } else {
            // Check if the existing cells needs to become an internal cell
            Cell &localCell = m_cells[cellId];
            if (isInterior && !localCell.isInterior()) {
                unsetGhostOwner(cellId, false);
                moveGhost2Internal(cellId);
            }

            // Save the adjacencies of the received cell, this adjacencies
            // will link together the recevied cell to the existing ones.
            FlatVector2D<long> &cellAdjacencies = linkAdjacencies[cellId];

            int nCellFaces = cell.getFaceCount();
            cellAdjacencies.reserve(nCellFaces);
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = cell.getAdjacencyCount(face);

                std::vector<long> faceAdjacencies;
                faceAdjacencies.reserve(nFaceAdjacencies);
                for (int k = 0; k < nFaceAdjacencies; ++k) {
                    faceAdjacencies.push_back(cell.getAdjacency(face, k));
                }

                cellAdjacencies.pushBack(faceAdjacencies);
            }
        }

        // Ranks involved in the communication
        if (cellOwner != m_rank) {
            involvedRanks.insert(cellOwner);
        }

        // Add the cell to the cell map
        cellMap.insert({{cellOriginalId, cellId}});
    }

    // Remap adjacencies
    for (auto cellId : receivedCells) {
        Cell &cell = m_cells[cellId];

        int nCellFaces = cell.getFaceCount();
        for (int face = 0; face < nCellFaces; ++face) {
            int nFaceAdjacencies = cell.getAdjacencyCount(face);
            for (int k = 0; k < nFaceAdjacencies; ++k) {
                long senderAdjacencyId = cell.getAdjacency(face, k);
                if (senderAdjacencyId < 0) {
                    continue;
                } else if (cellMap.count(senderAdjacencyId) == 0) {
					cell.deleteAdjacency(face, k);
                    continue;
                } else {
					long localAdjacencyId = cellMap.at(senderAdjacencyId);
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
                // We need to updated the adjacencies only if they are cells
                // that have been send.
                long senderAdjacencyId = cellLinkAdjacencies.getItem(face, k);
                if (cellMap.count(senderAdjacencyId) == 0) {
                    continue;
                }

                // If the send cell is already in the adjacency list there is
                // nothing to update.
                long localAdjacencyId  = cellMap.at(senderAdjacencyId);
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
    buildGhostExchangeInfo(involvedRankList);

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
    long bufferSize;
    MPI_Recv(&bufferSize, 1, MPI_LONG, sendRank, 30, m_communicator, MPI_STATUS_IGNORE);
    if (bufferSize == 0) {
        adaptionInfo.type = adaption::TYPE_NONE;
        return adaptionInfo;
    }

    IBinaryStream buffer(bufferSize);
    MPI_Recv(buffer.data(), buffer.getSize(), MPI_CHAR, sendRank, 31, m_communicator, MPI_STATUS_IGNORE);

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
    buildGhostExchangeInfo(std::vector<int>{sendRank, recvRank});

	// Return adaption info
    return adaptionInfo;
}

}

#endif
