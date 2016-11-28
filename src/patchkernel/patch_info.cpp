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
	\class PatchInfo

	\brief The PatchInfo class provides an interface for defining patch info.
*/

/*!
	Default constructor
*/
PatchInfo::PatchInfo(PatchKernel const *patch)
{
	setPatch(patch);
}

/*!
	Destructor
*/
PatchInfo::~PatchInfo()
{
}

/*!
	Sets the patch associated to the info.
*/
void PatchInfo::setPatch(PatchKernel const *patch)
{
	m_patch = patch;
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

	\param patch is patch from which the informations will be extracted
*/
void PatchInfo::update()
{
	extract();
}

#if BITPIT_ENABLE_MPI==1
/*!
	\ingroup patchkernel
	\class PatchGlobalInfo

	\brief Global information about the patch.
*/

/*!
	Creates a new info.

	\param patch is patch from which the informations will be extracted
*/
PatchGlobalInfo::PatchGlobalInfo(PatchKernel const *patch)
	: PatchInfo(patch)
{
	extract();
}

/*!
	Internal function to reset the information.
*/
void PatchGlobalInfo::_reset()
{
	m_cellLocalToGlobalMap.clear();
	m_nGlobalInternals.clear();
}

/*!
	Internal function to extract global information from the patch.
*/
void PatchGlobalInfo::_extract()
{
	long globalId;

	size_t exchangeDataSize = sizeof(globalId);
	std::unique_ptr<DataCommunicator> dataCommunicator;

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
	long offset = getCellGlobalOffset();

	// Evalaute the global id of the internal cells
	if (m_patch->getInternalCount() > 0) {
		std::map<long,long> nativeIds;
		for (auto itr = m_patch->internalConstBegin(); itr != m_patch->internalConstEnd(); ++itr) {
			long id = itr.getId();
			long nativeId = m_patch->_getCellNativeIndex(id);

			nativeIds.insert({nativeId, id});
		}

		globalId = offset;
		for (auto itr = nativeIds.begin(); itr != nativeIds.end(); ++itr) {
			m_cellLocalToGlobalMap.insert({itr->second, globalId++});
		}

		std::map<long,long>().swap(nativeIds);
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
int PatchGlobalInfo::getCellRankFromLocal(long id) const
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
int PatchGlobalInfo::getCellRankFromGlobal(long id) const
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
long PatchGlobalInfo::getCellGlobalId(long id) const
{
	return m_cellLocalToGlobalMap.at(id);
}

/*!
	Return the offset of the global numbering for the current parition.

	\result The offset of the global numbering for the current parition.
*/
long PatchGlobalInfo::getCellGlobalOffset() const
{
	return getCellGlobalOffset(m_patch->getRank());
}

/*!
	Return the offset of the global numbering for the specified parition.

	\param rank is the rank
	\result The offset of the global numbering for the current parition.
*/
long PatchGlobalInfo::getCellGlobalOffset(int rank) const
{
	long offset = 0;
	for (int i = 0; i < rank; ++i) {
		offset += m_nGlobalInternals[i];
	}

	return offset;
}

/*!
	Return the map between local indexes and global indexes.

	\result The map between local indexes and global indexes.
*/
const std::unordered_map<long, long> & PatchGlobalInfo::getCellGlobalMap() const
{
	return m_cellLocalToGlobalMap;
}
#endif

}
