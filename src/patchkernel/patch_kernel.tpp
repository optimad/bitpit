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

#ifndef __BITPIT_PATCH_KERNEL_TPP__
#define __BITPIT_PATCH_KERNEL_TPP__

#include <stdexcept>

namespace bitpit {

/*!
	Creates a clone of the specified pach.

	\param original is the patch that will be cloned
	\result A clone of the pach.
*/
template<typename patch_t>
std::unique_ptr<patch_t> PatchKernel::clone(const patch_t *original)
{
	static_assert(std::is_base_of<PatchKernel, patch_t>::value, "Specified pointer is not derived from PatchKernel");

	patch_t *clone = static_cast<patch_t *>(original->clone().release());

	return std::unique_ptr<patch_t>(clone);
}

/*!
	Deletes a list of cells.

	\param ids are the ids of the cells to be deleted
 */
template<typename IdStorage>
bool PatchKernel::deleteCells(const IdStorage &ids)
{
	if (!isExpert()) {
		return false;
	}

	// Deleteing the last internal cell requires some additional work. If the
	// ids of cells to be deleted contains the last internal cell, we delete
	// that cell ater deleting all other cells. In this way we make sure to
	// deleting the last internal cells just once (after deleting the last
	// internal, another cells becomes the last one and  that cells may be
	// on the deletion list, and on and so forth). The same applies for the
	// first ghost.
	bool deleteLastInternalCell = false;
#if BITPIT_ENABLE_MPI==1
	bool deleteFirstGhostCell = false;
#endif
	for (long id : ids) {
		if (id == m_lastInternalCellId) {
			deleteLastInternalCell = true;
			continue;
		}
#if BITPIT_ENABLE_MPI==1
		else if (id == m_firstGhostCellId) {
			deleteFirstGhostCell = true;
			continue;
		}
#endif

		deleteCell(id);
	}

	if (deleteLastInternalCell) {
		deleteCell(m_lastInternalCellId);
	}

#if BITPIT_ENABLE_MPI==1
	if (deleteFirstGhostCell) {
		deleteCell(m_firstGhostCellId);
	}
#endif

	return true;
}

/*!
	Deletes a list of vertices.

	\param ids are the ids of the vertices to be deleted
*/
template<typename IdStorage>
bool PatchKernel::deleteVertices(const IdStorage &ids)
{
	if (!isExpert()) {
		return false;
	}

	// Deleting the last vertex requires some additional work. If the ids
	// of vertices to be deleted contain the last vertex, that vertex is
	// deleted after deleting all other vertices. In this way we make sure
	// to delete the last vertex just once (after deleting the last vertex,
	// another vertex becomes the last one and that vertex may be on the
	// deletion list as well, and on and so forth). The same applies for
	// the first ghost.
	bool deleteLastInternalVertex = false;
#if BITPIT_ENABLE_MPI==1
	bool deleteFirstGhostVertex = false;
#endif
	for (long id : ids) {
		if (id == m_lastInternalVertexId) {
			deleteLastInternalVertex = true;
			continue;
		}
#if BITPIT_ENABLE_MPI==1
		else if (id == m_firstGhostVertexId) {
			deleteFirstGhostVertex = true;
			continue;
		}
#endif

		deleteVertex(id);
	}

	if (deleteLastInternalVertex) {
		deleteVertex(m_lastInternalVertexId);
	}

#if BITPIT_ENABLE_MPI==1
	if (deleteFirstGhostVertex) {
		deleteVertex(m_firstGhostVertexId);
	}
#endif

	return true;
}

/*!
	Deletes a list of interfaces.

	\param ids are the ids of the interfaces to be deleted
*/
template<typename IdStorage>
bool PatchKernel::deleteInterfaces(const IdStorage &ids)
{
	if (!isExpert()) {
		return false;
	}

	// Deleting the last interface requires some additional work. If the ids
	// of interfaces to be deleted contain the last interface, that interface
	// is deleted after deleting all other interfaces. In this way we make
	// sure to delete the last interface just once (after deleting the last
	// interface, another interface becomes the last one and that interface
	// may be on the deletion list as well, and on and so forth).
	std::size_t lastId;
	if (!m_interfaces.empty()) {
		lastId = m_interfaces.back().getId();
	} else {
		lastId = Interface::NULL_ID;
	}

	bool deleteLast = false;
	for (long id : ids) {
		if (id == lastId) {
			deleteLast = true;
			continue;
		}

		deleteInterface(id);
	}

	if (deleteLast) {
		deleteInterface(lastId);
	}

	return true;
}

/*!
	Renumber the ids of the items in the specified container.

	\param container is the container
	\param offset is the offset that will be used
*/
template<typename item_t, typename id_t>
std::unordered_map<id_t, id_t> PatchKernel::consecutiveItemRenumbering(PiercedVector<item_t, id_t> &container, long offset)
{
	// Build renumber map
	std::unordered_map<id_t, id_t> renumberMap;

	id_t counter = offset;
	for(const item_t &item : container) {
		id_t originalId = item.getId();
		id_t finalId    = counter++;

		renumberMap.insert({originalId, finalId});
	}

	// Renumber items
	mappedItemRenumbering(container, renumberMap);

	return renumberMap;
}

/*!
	Renumber the ids of the items in the specified container.

	\param container is the container
	\param renumberMap is the map that will be used for the renumer
*/
template<typename item_t, typename id_t>
void PatchKernel::mappedItemRenumbering(PiercedVector<item_t, id_t> &container,
                                        const std::unordered_map<id_t, id_t> &renumberMap)
{
	// Find an unused id
	id_t unusedId = std::numeric_limits<id_t>::max();
	while (container.exists(unusedId)) {
		--unusedId;
		if (unusedId < 0) {
			break;
		}
	}

	// Renumber the items
	std::unordered_map<id_t, id_t> conflictMap;

	for(const item_t &item : container) {
		// Original id of the item
		//
		// To get the original id of the item we use the conflict map, once
		// we get the id of an item we can remove it from the conflict map.
		id_t currentId = item.getId();

		id_t originalId;
		auto conflictMapItr = conflictMap.find(currentId);
		if (conflictMapItr != conflictMap.end()) {
			originalId = conflictMapItr->second;
			conflictMap.erase(conflictMapItr);
		} else {
			originalId = currentId;
		}

		// Final id of the item
		id_t finalId = renumberMap.at(originalId);

		// Nothing else to do if the item has already the correct id, otherwise
		// we need to renumber the item.
		if (currentId == finalId) {
			continue;
		}

		// Check if the final id is already in the container.
		//
		// If there is a conflict, the item that holds the final id will
		// be temporarly renumbered with the unused id find before, the current
		// item will be renumberd with is final id and eventually the item with
		// the temporary id will be assoicated with the, now avilable, id of
		// the current item.
		id_t oldConflictId;
		id_t newConflicdId;
		id_t tmpConflictId;

		bool conflict = container.exists(finalId);
		if (conflict) {
			if (unusedId < 0) {
				throw std::runtime_error("Renumbering requires an unused id, but all ids are used.");
			}

			oldConflictId = finalId;
			newConflicdId = currentId;
			tmpConflictId = unusedId;
		}

		// Temporary renumber of the conflicting item
		//
		// There is no need to change the id of the element, it is sufficient
		// update the container. The id will be updated later.
		if (conflict) {
			container.updateId(oldConflictId, tmpConflictId);
		}

		// Renumber of the current element
		container[currentId].setId(finalId);
		container.updateId(currentId, finalId);

		// Renumber the conflicting element with the previous id of the item
		if (conflict) {
			container.updateId(tmpConflictId, newConflicdId);
			container[newConflicdId].setId(newConflicdId);

			// Update conflic map
			auto conflictMapItr = conflictMap.find(oldConflictId);
			if (conflictMapItr == conflictMap.end()) {
				conflictMap.insert({newConflicdId, oldConflictId});
			} else {
				conflictMap.insert({newConflicdId, conflictMapItr->second});
				conflictMap.erase(conflictMapItr);
			}
		}
	}
}

}

#endif
