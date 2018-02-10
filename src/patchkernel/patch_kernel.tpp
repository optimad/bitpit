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
