/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2022 OPTIMAD engineering Srl
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

namespace bitpit {

/**
 * Confirm that the specified halo layer is the halo layer of the given cell.
 *
 * The check is performed looping the neighbours and searching a negihbour that belongs
 * to the previous layer. If that neguhbour is found, the halo layer is confirmed.
 *
 * \param cell is the cell that will be check
 * \param haloLayer is the halo layer the cell supposedly belongs to
 * \param excludeList is a list of cells that will not be considered when processing the
 * neighbours
 * \result Return true if the cell belongs to the specified halo layer, false otherwise.
 */
template<typename ExcludeList>
bool PatchKernel::confirmCellHaloLayer(const Cell &cell, int haloLayer, const ExcludeList &excludeList) const
{
	// Search if a face neighbour is in the previous halo layer
	//
	// We first search in the face neighbours only, because its faster than
	// searching among all the neighbours.
	//
	// Searching for neighbours in the previous layer works also for ghosts in
	// the first layer, because inner cells has a layer equal to minus one.
	const long *adjacencies = cell.getAdjacencies();
	int nCellAdjacencies = cell.getAdjacencyCount();
	for (int i = 0; i < nCellAdjacencies; ++i) {
		long neighId = adjacencies[i];
		if (excludeList.count(neighId) > 0) {
			continue;
		} else if (getCellHaloLayer(neighId) == (haloLayer - 1)) {
			return true;
		}
	}

	// Search if a generic neighbour is in the previous halo layer
	//
	// If we haven't found a cell in the previous layer among the face neighbours,
	// we perform a more expensive search among all the neighbours. a search among
	// all the neighbours.
	//
	// Searching for neighbours in the previous layer works also for ghosts in the
	// first layer, because inner cells has a layer equal to minus one.
	std::vector<long> neighIds;
	findCellNeighs(cell.getId(), &neighIds);
	for (long neighId : neighIds) {
		if (excludeList.count(neighId) > 0) {
			continue;
		} else if (getCellHaloLayer(neighId) == (haloLayer - 1)) {
			return true;
		}
	}

	// There are no neighbours in the previous halo layer
	return false;
}

}

#endif
