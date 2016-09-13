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

#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif

#include "patch_info.hpp"
#include "patch_kernel.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class PatchInfo

	\brief The PatchInfo class provides an interface for defining patch info.
*/

/*!
	@}
*/

#if BITPIT_ENABLE_MPI==1
/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class PatchGlobalInfo

	\brief Global information about the patch.
*/

/*!
	Creates a new info.
*/
PatchGlobalInfo::PatchGlobalInfo(PatchKernel const *patch)
{
	extract(patch);
}

/*!
	Extracts the global information from the patch.

	\param patch is patch from which the informations will be extracted
*/
void PatchGlobalInfo::extract(PatchKernel const *patch)
{
	m_patch = patch;

	long globalId;

	size_t exchangeDataSize = sizeof(globalId);
	std::unique_ptr<DataCommunicator> dataCommunicator;

	// Initialize the communication of the ghost global ids
	if (m_patch->getProcessorCount() > 1) {
		// Create the data communicator
		dataCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(m_patch->getCommunicator()));
		dataCommunicator->setTag(108);

		// Set and start the receives
		for (const auto entry : m_patch->getGhostExchangeTargets()) {
			const int rank = entry.first;
			const auto &list = entry.second;

			dataCommunicator->setRecv(rank, list.size() * exchangeDataSize);
			dataCommunicator->startRecv(rank);
		}
	}

	// Get the internal count of all the partitions
	m_nGlobalInternals.resize(m_patch->getProcessorCount());
	if (m_patch->getProcessorCount() > 1) {
		long nLocalInternals = m_patch->getInternalCount();
		MPI_Allgather(&nLocalInternals, 1, MPI_LONG, m_nGlobalInternals.data(), 1, MPI_LONG, m_patch->getCommunicator());
	} else {
		m_nGlobalInternals[0] = m_patch->getInternalCount();
	}

	// Evaluate the offset for the current partition
	long offset = 0;
	for (int i = 0; i < m_patch->getRank() - 1; ++i) {
		offset += m_nGlobalInternals[i];
	}

	// Evalaute the global id of the internal cells
	if (m_patch->getInternalCount() > 0) {
		auto cbeginInternals = m_patch->m_cells.cbegin();
		auto cendInternals   = ++m_patch->m_cells.getConstIterator(m_patch->m_lastInternalId);

		globalId = offset;
		for (auto cellItr = cbeginInternals; cellItr != cendInternals; ++cellItr) {
			m_cellLocalToGlobalMap.insert({cellItr.getId(), globalId++});
		}
	}

	// Communicate the global id of the ghost cells
	if (m_patch->getProcessorCount() > 1) {
		// Set and start the sends
		for (const auto entry : m_patch->getGhostExchangeSources()) {
			const int rank = entry.first;
			auto &list = entry.second;

			dataCommunicator->setSend(rank, list.size() * exchangeDataSize);
			SendBuffer &buffer = dataCommunicator->getSendBuffer(rank);
			for (long id : list) {
				buffer << m_cellLocalToGlobalMap.at(id);
			}
			dataCommunicator->startSend(rank);
		}

		// Receive the global ids of the ghosts
		int nCompletedRecvs = 0;
		while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
			int rank = dataCommunicator->waitAnyRecv();
			const auto &list = m_patch->getGhostExchangeTargets(rank);

			RecvBuffer &buffer = dataCommunicator->getRecvBuffer(rank);
			for (long id : list) {
				buffer >> globalId;
				m_cellLocalToGlobalMap.insert({id, globalId});
			}

			++nCompletedRecvs;
		}

		// Wait for the sends to finish
		dataCommunicator->waitAllSends();
	}
}

/*!
	Gets the rank of the cell with the specified local id

	\param id is the local id of the cell
	\return The rank of the specified cell.
*/
int PatchGlobalInfo::getCellRankFromLocal(long id)
{
	auto ghostOwnerItr = m_patch->m_ghostOwners.find(id);
	if (ghostOwnerItr != m_patch->m_ghostOwners.end()) {
		return ghostOwnerItr->second;
	} else {
		return m_patch->getRank();
	}
}

/*!
	Gets the rank of the cell with the specified global id

	\param id is the global id of the cell
	\return The rank of the specified cell.
*/
int PatchGlobalInfo::getCellRankFromGlobal(long id)
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
	Return the global id of the cell with the specified local id

	\param id is the local id of the cell
	\return The global id of the specified cell.
*/
long PatchGlobalInfo::getCellGlobalId(long id)
{
	return m_cellLocalToGlobalMap.at(id);
}

/*!
	Return the map between local indexes and global indexes.

	\result The map between local indexes and global indexes.
*/
const std::unordered_map<long, long> & PatchGlobalInfo::getCellGlobalMap()
{
	return m_cellLocalToGlobalMap;
}

/*!
	@}
*/
#endif

}
