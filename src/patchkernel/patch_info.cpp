/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#	include <mpi.h>
#endif

#include "patch_info.hpp"
#include "patch_kernel.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	\class PatchInfo

	\brief The PatchInfo class provides an interface for defining patch info.
*/

/*!
	Default constructor

	\param patch is patch from which the informations will be extracted
*/
PatchInfo::PatchInfo(PatchKernel const *patch)
	: m_patch(nullptr)
{
	setPatch(patch, false);
}

/*!
	Destructor
*/
PatchInfo::~PatchInfo()
{
}

/*!
	Sets the patch associated to the info.

	\param patch is patch from which the informations will be extracted
*/
void PatchInfo::setPatch(PatchKernel const *patch)
{
	setPatch(patch, true);
}

/*!
	Sets the patch associated to the info.

	\param patch is patch from which the informations will be extracted
	\param initialize if set to true the initializatino function will be called
*/
void PatchInfo::setPatch(PatchKernel const *patch, bool initialize)
{
	// Reset information
	if (m_patch != nullptr) {
		reset();
	}

	// Set the patch
	m_patch = patch;

	// Call initialization
	if (initialize) {
		_init();
	}
}

/*!
	Resets the information.
*/
void PatchInfo::reset()
{
	_reset();
}

/*!
	Extracts the information.
*/
void PatchInfo::extract()
{
	if (m_patch == nullptr) {
		return;
	}

	// Reset information
	reset();

	// Extract new information
	_extract();
}

/*!
	Updates the information.
*/
void PatchInfo::update()
{
	extract();
}

/*!
	\ingroup patchkernel
	\class PatchNumberingInfo

	\brief Numbering information about the patch.
*/

/*!
	Creates a new info.

	\param patch is patch from which the informations will be extracted
*/
PatchNumberingInfo::PatchNumberingInfo(PatchKernel const *patch)
	: PatchInfo(patch)
{
	PatchNumberingInfo::_init();

	extract();
}

/*!
	Internal function to initialize the information.
*/
void PatchNumberingInfo::_init()
{
}

/*!
	Internal function to reset the information.
*/
void PatchNumberingInfo::_reset()
{
	m_cellLocalToConsecutiveMap.clear();
#if BITPIT_ENABLE_MPI==1
	m_nGlobalInternals.clear();
#endif
}

/*!
	Internal function to extract information from the patch.
*/
void PatchNumberingInfo::_extract()
{
	long consecutiveId;

#if BITPIT_ENABLE_MPI==1
	// Initialize communications
	size_t exchangeDataSize = sizeof(consecutiveId);
	std::unique_ptr<DataCommunicator> dataCommunicator;

	if (m_patch->getProcessorCount() > 1) {
		// Create the data communicator
		dataCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(m_patch->getCommunicator()));

		// Set and start the receives
		for (const auto entry : m_patch->getGhostExchangeTargets()) {
			const int rank = entry.first;
			const auto &list = entry.second;

			dataCommunicator->setRecv(rank, list.size() * exchangeDataSize);
			dataCommunicator->startRecv(rank);
		}
	}
#endif

#if BITPIT_ENABLE_MPI==1
	// Get the internal count of all the partitions
	m_nGlobalInternals.resize(m_patch->getProcessorCount());
	if (m_patch->getProcessorCount() > 1) {
		long nLocalInternals = m_patch->getInternalCount();
		MPI_Allgather(&nLocalInternals, 1, MPI_LONG, m_nGlobalInternals.data(), 1, MPI_LONG, m_patch->getCommunicator());
	} else {
		m_nGlobalInternals[0] = m_patch->getInternalCount();
	}
#endif

	// Evalaute the consecutive id of the internal cells
	if (m_patch->getInternalCount() > 0) {
		std::map<long,long> nativeIds;
		for (auto itr = m_patch->internalConstBegin(); itr != m_patch->internalConstEnd(); ++itr) {
			long id = itr.getId();
			long nativeId = m_patch->_getCellNativeIndex(id);

			nativeIds.insert({nativeId, id});
		}

#if BITPIT_ENABLE_MPI==1
		consecutiveId = getCellGlobalCountOffset();
#else
		consecutiveId = 0;
#endif
		for (auto itr = nativeIds.begin(); itr != nativeIds.end(); ++itr) {
			m_cellLocalToConsecutiveMap.insert({itr->second, consecutiveId++});
		}

		std::map<long,long>().swap(nativeIds);
	}

#if BITPIT_ENABLE_MPI==1
	// Communicate the consecutive id of the ghost cells
	if (m_patch->getProcessorCount() > 1) {
		// Set and start the sends
		for (const auto entry : m_patch->getGhostExchangeSources()) {
			const int rank = entry.first;
			auto &list = entry.second;

			dataCommunicator->setSend(rank, list.size() * exchangeDataSize);
			SendBuffer &buffer = dataCommunicator->getSendBuffer(rank);
			for (long id : list) {
				buffer << m_cellLocalToConsecutiveMap.at(id);
			}
			dataCommunicator->startSend(rank);
		}

		// Receive the consecutive ids of the ghosts
		int nCompletedRecvs = 0;
		while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
			int rank = dataCommunicator->waitAnyRecv();
			const auto &list = m_patch->getGhostExchangeTargets(rank);

			RecvBuffer &buffer = dataCommunicator->getRecvBuffer(rank);
			for (long id : list) {
				buffer >> consecutiveId;
				m_cellLocalToConsecutiveMap.insert({id, consecutiveId});
			}

			++nCompletedRecvs;
		}

		// Wait for the sends to finish
		dataCommunicator->waitAllSends();
	}
#endif
}

/*!
	Return the consecutive id of the cell with the specified local id

	\param id is the local id of the cell
	\return The consecutive id of the specified cell.
*/
long PatchNumberingInfo::getCellConsecutiveId(long id) const
{
	return m_cellLocalToConsecutiveMap.at(id);
}

/*!
	Return the map between local indexes and consecutive indexes.

	\result The map between local indexes and consecutive indexes.
*/
const std::unordered_map<long, long> & PatchNumberingInfo::getCellConsecutiveMap() const
{
	return m_cellLocalToConsecutiveMap;
}


#if BITPIT_ENABLE_MPI==1
/*!
	Gets the global number of cells.

	\return The global number of cells.
*/
long PatchNumberingInfo::getCellGlobalCount() const
{
	long nGlobalCells = 0;
	for (long count : m_nGlobalInternals) {
		nGlobalCells += count;
	}

	return nGlobalCells;
}

/*!
	Return the offset of the global cell count for the current parition.

	\result The offset of the global cell count for the current parition.
*/
long PatchNumberingInfo::getCellGlobalCountOffset() const
{
	return getCellGlobalCountOffset(m_patch->getRank());
}

/*!
	Return the offset of the global cell count for the specified parition.

	\param rank is the rank
	\result The offset of the global cell count for the specified parition.
*/
long PatchNumberingInfo::getCellGlobalCountOffset(int rank) const
{
	long offset = 0;
	for (int i = 0; i < rank; ++i) {
		offset += m_nGlobalInternals[i];
	}

	return offset;
}

/*!
	Return the global id of the cell with the specified local id

	\param id is the local id of the cell
	\return The global id of the specified cell.
*/
long PatchNumberingInfo::getCellGlobalId(long id) const
{
    // Global ids and consecutive ids are the same
	return m_cellLocalToConsecutiveMap.at(id);
}

/*!
	Return the map between local indexes and global indexes.

	\result The map between local indexes and global indexes.
*/
const std::unordered_map<long, long> & PatchNumberingInfo::getCellGlobalMap() const
{
    // Global ids and consecutive ids are the same
	return m_cellLocalToConsecutiveMap;
}


/*!
	Gets the rank of the cell with the specified local id

	\param id is the local id of the cell
	\return The rank of the specified cell.
*/
int PatchNumberingInfo::getCellRankFromLocal(long id) const
{
	auto ghostOwnerItr = m_patch->m_ghostOwners.find(id);
	if (ghostOwnerItr != m_patch->m_ghostOwners.end()) {
		return ghostOwnerItr->second;
	} else {
		return m_patch->getRank();
	}
}

/*!
	Gets the rank of the cell with the specified consecutive id

	\param id is the consecutive id of the cell
	\return The rank of the specified cell.
*/
int PatchNumberingInfo::getCellRankFromConsecutive(long id) const
{
	long offset = 0;
	for (int k = 0; k < m_patch->getProcessorCount(); ++k) {
		offset += m_nGlobalInternals[k];
		if (id < offset) {
			return k;
		}
	}

	return -1;
}

/*!
	Gets the rank of the cell with the specified global id

	\param id is the global id of the cell
	\return The rank of the specified cell.
*/
int PatchNumberingInfo::getCellRankFromGlobal(long id) const
{
    // Global ids and consecutive ids are the same
    return getCellRankFromConsecutive(id);
}
#endif

}
