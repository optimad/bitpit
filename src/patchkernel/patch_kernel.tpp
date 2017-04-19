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
		id_t id_original = item.getId();
		id_t id_final    = counter++;

		renumberMap.insert({id_original, id_final});
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
		id_t id_current = item.getId();

		id_t id_original;
		auto conflictMapItr = conflictMap.find(id_current);
		if (conflictMapItr != conflictMap.end()) {
			id_original = conflictMapItr->second;
			conflictMap.erase(conflictMapItr);
		} else {
			id_original = id_current;
		}

		// Final id of the item
		id_t id_final = renumberMap.at(id_original);

		// Nothing else to do if the item has already the correct id, otherwise
		// we need to renumber the item.
		if (id_current == id_final) {
			continue;
		}

		// Check if the final id is already in the container.
		//
		// If there is a conflict, the item that holds the final id will
		// be temporarly renumbered with the unused id find before, the current
		// item will be renumberd with is final id and eventually the item with
		// the temporary id will be assoicated with the, now avilable, id of
		// the current item.
		id_t conflic_id_old;
		id_t conflic_id_new;
		id_t conflic_id_tmp;

		bool conflict = container.exists(id_final);
		if (conflict) {
			if (unusedId < 0) {
				throw std::runtime_error("Renumbering requires an unused id, but all ids are used.");
			}

			conflic_id_old = id_final;
			conflic_id_new = id_current;
			conflic_id_tmp = unusedId;
		}

		// Temporary renumber of the conflicting item
		//
		// There is no need to change the id of the element, it is sufficient
		// update the container. The id will be updated later.
		if (conflict) {
			container.updateId(conflic_id_old, conflic_id_tmp);
		}

		// Renumber of the current element
		container[id_current].setId(id_final);
		container.updateId(id_current, id_final);

		// Renumber the conflicting element with the previous id of the item
		if (conflict) {
			container.updateId(conflic_id_tmp, conflic_id_new);
			container[conflic_id_new].setId(conflic_id_new);

			// Update conflic map
			auto conflictMapItr = conflictMap.find(conflic_id_old);
			if (conflictMapItr == conflictMap.end()) {
				conflictMap.insert({conflic_id_new, conflic_id_old});
			} else {
				conflictMap.insert({conflic_id_new, conflictMapItr->second});
				conflictMap.erase(conflictMapItr);
			}
		}
	}
}

}

#endif
