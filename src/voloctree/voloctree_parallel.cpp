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

#include <mpi.h>
#include <unordered_set>

#include "voloctree.hpp"

using namespace std;

namespace bitpit {

/*!
	\ingroup voloctree
	@{
*/

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void VolOctree::setCommunicator(MPI_Comm communicator)
{
	PatchKernel::setCommunicator(communicator);

	m_tree.setComm(communicator);
}

/*!
	Updates the partition to optimize the load balance.

	\param trackChanges if set to true the changes to the patch will be
	tracked
	\result Returns all the changes applied to the patch.
*/
const std::vector<adaption::Info> VolOctree::_balancePartition(bool trackChanges)
{
	// Updating the tree
	log::cout() << ">> Load balancing...";

	m_tree.loadBalance();
	m_lastTreeOperation = OP_LOAD_BALANCE;

	log::cout() << " Done" << std::endl;

	// Sync the patch
	return sync(trackChanges);
}

/*!
	Rebuilds the ghost data needed for data exchange

	\param ghostTreeIds is the list of tree ids of the ghosts
*/
void VolOctree::rebuildGhostExchangeData(std::unordered_map<int, std::vector<uint32_t>> ghostTreeIds)
{
	const int MSG_TAG = 100;

	resetGhostExchangeData();
	std::unordered_map<short, std::unordered_map<long, long>> &ghostExchangeData = getGhostExchangeData();

	const auto &communicator = getCommunicator();
	for (int i = 0; i < getProcessorCount(); ++i) {
		long nRankGhosts;
		if (getRank() == i) {
			for (int j = 0; j < getProcessorCount(); ++j) {
				if (i == j) {
					continue;
				}

				// Send the number of ghosts octants for which we want to
				// retreive the id
				if (ghostTreeIds.count(j) != 0) {
					nRankGhosts = ghostTreeIds.at(j).size();
				} else {
					nRankGhosts = 0;
				}
				MPI_Send(&nRankGhosts, 1, MPI_LONG, j, MSG_TAG, communicator);

				// Retreive the id of the octants
				if (nRankGhosts > 0) {
					// List of tree ids for the considered rank
					const std::vector<uint32_t> &rankGhostTreeIds = ghostTreeIds.at(j);

					// List of global tree ids for the considered rank
					std::vector<uint64_t> rankGhostGlobalTreeIds;
					rankGhostGlobalTreeIds.reserve(nRankGhosts);
					for (const uint32_t &treeId : rankGhostTreeIds) {
						rankGhostGlobalTreeIds.push_back(m_tree.getGhostGlobalIdx(treeId));
					}

					// Send the list of global tree ids for which we want the
					// cells ids on the other processor
					MPI_Send(rankGhostGlobalTreeIds.data(), nRankGhosts, MPI_UINT64_T, j, MSG_TAG, communicator);

					// Receive the ids on the other processor
					std::vector<long> rankGhostRemoteIds(nRankGhosts);
					MPI_Recv(rankGhostRemoteIds.data(), nRankGhosts, MPI_LONG, j, MSG_TAG, communicator, MPI_STATUS_IGNORE);

					// Store the id of the ghosts
					for (int n = 0; n < nRankGhosts; ++n) {
						// Id of the ghost in the local processor
						OctantInfo ghostOctantInfo(rankGhostTreeIds[n], false);
						long ghostLocalId = getOctantId(ghostOctantInfo);

						// Store the id of the ghost in the cuttent patch and
						// the id of the same cell in the neighbour partition
						ghostExchangeData[j].insert({{rankGhostRemoteIds[n], ghostLocalId}});
					}
				}
			}
		} else {
			// Get the number of octants for which we want to retreive the id
			MPI_Recv(&nRankGhosts, 1, MPI_LONG, i, MSG_TAG, communicator, MPI_STATUS_IGNORE);

			// Send the id of the octants
			if (nRankGhosts > 0) {
				// Receive the requested octants
				std::vector<uint64_t> octantGlobalTreeIds(nRankGhosts);
				MPI_Recv(octantGlobalTreeIds.data(), nRankGhosts, MPI_UINT64_T, i, MSG_TAG, communicator, MPI_STATUS_IGNORE);

				// Retreive the ids
				std::vector<long> octantIds;
				octantIds.reserve(nRankGhosts);
				for (uint64_t globalTreeId : octantGlobalTreeIds) {
					// Local tree of the octant
					uint32_t treeId = m_tree.getLocalIdx(globalTreeId);

					// Local id of the cell associated to the octant
					OctantInfo octantInfo(treeId, true);
					long id = getOctantId(octantInfo);

					octantIds.push_back(id);
				}

				// Send the ids
				MPI_Send(octantIds.data(), nRankGhosts, MPI_LONG, i, MSG_TAG, communicator);
			}
		}
	}
}

/*!
	@}
*/

}

#endif
