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

#include <unordered_map>
#include <unordered_set>

#include "patch_kernel.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	\brief The namespace 'adaption' contains the routines and the data
	structures for handling patch adaption.
*/
namespace adaption
{

/*!
	\enum adaption::Marker

	\brief The Marker enum defines the type of adaption that has been
	requested.
*/

/*!
	\enum adaption::Type

	\brief The Type enum defines the type of adaption that has been
	performed.
*/

/*!
	\enum adaption::Entity

	\brief The Entity enum defines the type of entities on which the
	adaption has been performed.
*/

/*!
	\struct Info

	\brief The Info struct defines the infomation associated to an
	adaption.
*/

/*!
	\class InfoCollection

	\brief The InfoCollection class is a container that holds one or more
	adaption info items.
*/

/*!
	Default constructor.
*/
InfoCollection::InfoCollection()
{
	m_cachedTypes.insert(adaption::TYPE_RENUMBERING);
	m_cachedTypes.insert(adaption::TYPE_DELETION);
	m_cachedTypes.insert(adaption::TYPE_CREATION);
	m_cachedTypes.insert(adaption::TYPE_PARTITION_RECV);
	m_cachedTypes.insert(adaption::TYPE_PARTITION_SEND);
}


/*!
 * Get the size of the collection
 *
 * \return The number of entries in the collection.
 */
std::size_t InfoCollection::size() const
{
    return m_collection.size();
}

/*!
 * Get a read/write iterator to the first entry of the collection.
 *
 * \return The iterator pointing to the first entry of the collection.
 */
InfoCollection::iterator InfoCollection::begin() noexcept
{
    return iterator(m_collection.begin());
}

/*!
 * Get a read/write iterator that points one past the last enty of the collection.
 *
 * \return The iterator pointing to one past the last entry of the collection.
 */
InfoCollection::iterator InfoCollection::end() noexcept
{
    return iterator(m_collection.end());
}

/*!
 * Get a read-only iterator to the first enty of the collection.
 *
 * \return the const iterator pointing to the first entry of the collection.
 */
InfoCollection::const_iterator InfoCollection::begin() const noexcept
{
    return const_iterator(cbegin());
}

/*!
 * Get a read-only iterator that points one past the last enty of the collection.
 *
 * \return the const iterator pointing to one past the last entry of the collection.
 */
InfoCollection::const_iterator InfoCollection::end() const noexcept
{
    return const_iterator(cend());
}

/*!
 * Get a read-only iterator to the first enty of the collection.
 *
 * \return the const iterator pointing to the first entry of the collection.
 */
InfoCollection::const_iterator InfoCollection::cbegin() const noexcept
{
    return const_iterator(m_collection.cbegin());
}

/*!
 * Get a read-only iterator that points one past the last enty of the collection.
 *
 * \return the const iterator pointing to one past the last entry of the collection.
 */
InfoCollection::const_iterator InfoCollection::cend() const noexcept
{
    return const_iterator(m_collection.cend());
}

/*!
	Returns a constant reference to the vector used internally by the
	collection to store its owned adaption info.

	\result A constant reference to the vector used internally by the
	collection to store its owned adaption info.
*/
const std::vector<Info> & InfoCollection::data() const noexcept
{
	return m_collection;
}

/*!
	Returns a reference to the vector used internally by the collection to
	store its owned adaption info.

	\result A reference to the vector used internally by the collection to
	store its owned adaption info.
*/
std::vector<Info> & InfoCollection::data() noexcept
{
	return m_collection;
}

/*!
	Insert an empty adaption info.

	\result The index of the newly added adaption info.
*/
std::size_t InfoCollection::insert()
{
	m_collection.emplace_back();

	return m_collection.size() - 1;
}

/*!
	Insert an adaption info with the requested data.

	If the requested type is among the types that are cached and an adaption
	info with the requested data already exists, the existing value will be
	returned, otherwise a new adaption info will be created. If the requested
	type is not among the types that are cached a new adaption info will be
	created even if the collection already contains adaption info with the
	requested data.

	\param type is the type of adaption info
	\param entity is the entity associated to the adaption info
	\param rank is the rank associated to the adaption info
	\result The position inside the collection of the adaption info.
*/
std::size_t InfoCollection::insert(Type type, Entity entity, int rank)
{
	infoData_t infoData = infoData_t(type, entity, rank);
	bool useCache = (m_cachedTypes.count(type) > 0);

	if (useCache) {
		if (m_cache.count(infoData) > 0) {
			return m_cache.at(infoData);
		}
	}

	m_collection.emplace_back(type, entity, rank);

	std::size_t id = m_collection.size() - 1;
	if (useCache) {
		m_cache.insert({{infoData, id}});
	}
	return id;
}

/*!
	Copy the specified adaption in the collection.

	If the requested type is among the types that are cached and an adaption
	info with the requested data already exists, the existing value will be
	returned, otherwise a new adaption info will be created. If the requested
	type is not among the types that are cached a new adaption info will be
	created even if the collection already contains adaption info with the
	same data.

	\param info is the adaption info that will be moved in the collection
	\result The index of the newly added (or already existing) adaption
	info
*/
std::size_t InfoCollection::insert(const Info &other)
{
	return insert(Info(other));
}

/*!
	Move the specified adaption in the collection.

	If the requested type is among the types that are cached and an adaption
	info with the requested data already exists, the existing value will be
	returned, otherwise a new adaption info will be created. If the requested
	type is not among the types that are cached a new adaption info will be
	created even if the collection already contains adaption info with the
	same data.

	\param info is the adaption info that will be moved in the collection
	\result The index of the newly added (or already existing) adaption
	info
*/
std::size_t InfoCollection::insert(Info &&other)
{
	std::size_t id = insert(other.type, other.entity, other.rank);
	Info &info = m_collection[id];

	if (info.previous.empty()) {
		info.previous = std::move(other.previous);
	} else {
		appendIds(std::move(other.previous), true, &(info.previous));
	}

	if (info.current.empty()) {
		info.current = std::move(other.current);
	} else {
		appendIds(std::move(other.current), true, &(info.current));
	}

	return id;
}

/*!
	Erase the adaption info at the specified position in the collection.

	If the provided position is not valid, no actions will be performed.

	\param id is the index of the adaption info
	\result Returns true if the adaption info was deleted, otherwise it
	return false.
*/
bool InfoCollection::erase(std::size_t id)
{
    if (id >= size()) {
        return false;
    }

    const Info &info = m_collection[id];
    if (m_cachedTypes.count(info.type)) {
        infoData_t infoData = infoData_t(info.type, info.entity, info.rank);
        auto cacheItr = m_cache.find(infoData);
        if (cacheItr != m_cache.end()) {
            m_cache.erase(cacheItr);
        }
    }

    m_collection.erase(m_collection.begin() + id);

    return true;
}

/*!
	Erase all the adaption info of the specified type.

	\param type is the type of adaption info that wil be deleted
	\result Returns the number of adaption info that have been deleted
*/
std::size_t InfoCollection::erase(Type type)
{
    std::size_t nErasedInfo = 0;

    std::size_t id = 0;
    while (id < m_collection.size()) {
        const Info &info = at(id);
        if (info.type != type) {
            ++id;
            continue;
        }

        erase(id);
        ++nErasedInfo;
    }

    return nErasedInfo;
}

/*!
	Returns a reference to the adaption info at the specified position in
	the collection.

	\param id is the index of the adaption info
	\result Returns a reference to requested adaption info.
*/
Info & InfoCollection::at(std::size_t id)
{
	if (id >= m_collection.size()) {
		throw std::out_of_range("Requested adaption info is not in the collection");
	}

	return (*this)[id];
}

/*!
	Returns a constant reference to the adaption info at the specified position
	in the collection.

	\param id is the index of the adaption info
	\result Returns a constant reference to requested adaption info.
*/
const Info & InfoCollection::at(std::size_t id) const
{
	if (id >= m_collection.size()) {
		throw std::out_of_range("Requested adaption info is not in the collection");
	}

	return (*this)[id];
}

/*!
	Get a reference to the specified adaption info

	\param id is the index of the requested adaption info
	\result Returns a reference to requested adaption info.
*/
Info & InfoCollection::operator[](std::size_t id)
{
	return m_collection[id];
}

/*!
	Get a constant reference to the specified adaption info

	\param id is the index of the requested adaption info
	\result Returns a constant reference to requested adaption info
*/
const Info & InfoCollection::operator[](std::size_t id) const
{
	return m_collection[id];
}

/*!
	Dumps the collection of adaption info.

	Once the collection is dumped, the internal list of adaption info and
	the internal cache are deleted.

	\result The collection of adaption info.
*/
std::vector<Info> InfoCollection::dump()
{
	std::vector<Info> exportedCollection;
	exportedCollection.swap(m_collection);

	m_cache.clear();

	return exportedCollection;
}

/*!
	Append the ids contained in the source list into the destination list.

	Optionally, it is possible to append only the ids of the source list that
	are not already contained in the destination list.

	\param src is the list that contains the ids to be appended
	\param unique if set to true, only the ids of the source list that are not
	already contained in the destination list will be appended
	\param dst is the list where the ids will be copied to
*/
void InfoCollection::appendIds(std::vector<long> src, bool unique, std::vector<long> *dst)
{
	// Early return if the source is empty
	if (src.empty()) {
		return;
	}

	// Early return if the destination is empty
	if (dst->empty()) {
		*dst = std::move(src);
		return;
	}

	// Append the ids
	if (unique) {
		std::unordered_set<long> dstSet(dst->begin(), dst->end());
		for (long id : src) {
			if (dstSet.count(id) == 0) {
				dst->push_back(id);
			}
		}
	} else {
		dst->insert(dst->end(), src.begin(), src.end());
	}
}

/*!
	Remove the ids contained in the source list from the destination list.

	\param src is the list that contains the ids to be removed
	\param dst is the list where the ids will be copied to
	\result The number of ids that have been deleted.
*/
std::size_t InfoCollection::removeIds(std::unordered_set<long> src, std::vector<long> *dst)
{
	// Early return if the source is empty
	if (src.empty()) {
		return 0;
	}

	// Early return if the destination is empty
	if (dst->empty()) {
		return 0;
	}

	// Remove the ids
	std::size_t nDeletedIds = 0;
	for (auto itr = dst->begin(); itr != dst->end(); ++itr) {
		long vertexId = *itr;
		if (src.count(vertexId) > 0) {
			++nDeletedIds;
		} else if (nDeletedIds > 0) {
			*(itr - nDeletedIds) = *itr;
		}
	}

	if (nDeletedIds > 0) {
		dst->resize(dst->size() - nDeletedIds);
		dst->shrink_to_fit();
	}

	return nDeletedIds;
}

}

/*!
	\class FlatMapping
	\ingroup patchkernel

	\brief The FlatMapping class allows to generate a mapping between an
	id-base numeration to a continuous-index numeration.
*/

/*!
	Default constructor.
*/
FlatMapping::FlatMapping()
	: m_patch(nullptr)
{
}

/*!
	Creates a new flat mapping.

	\param patch is the patch from witch the flat numbering will be built
*/
FlatMapping::FlatMapping(PatchKernel *patch)
	: m_patch(patch)
{
}

/*!
	Gets the numbering associated to the flat mapping.

	\result The numbering associated to the flat mapping.
*/
const std::vector<long> & FlatMapping::getNumbering() const
{
	return m_numbering;
}

/*!
	Gets the mapping associated to the flat mapping.

	\result The mapping associated to the flat mapping.
*/
const std::vector<long> & FlatMapping::getMapping() const
{
	return m_mapping;
}

/*!
	\class CellFlatMapping
	\ingroup patchkernel

	\brief The CellFlatMapping class allows to generate a cell mapping
	between an id-base numeration to a continuous-index numeration.
*/

/*!
	Default constructor.
*/
CellFlatMapping::CellFlatMapping()
	: FlatMapping()
{
}

/*!
	Creates a new cell flat mapping.

	\param patch is the patch from witch the flat numbering will be built
*/
CellFlatMapping::CellFlatMapping(PatchKernel *patch)
	: FlatMapping(patch)
{
	m_numbering.reserve(m_patch->getCellCount());
	m_mapping.reserve(m_patch->getCellCount());

	long flatId = -1;
	for (const auto &cell : m_patch->getCells()) {
		flatId++;

		m_numbering.emplace_back();
		long &cellId = m_numbering.back();
		cellId = cell.getId();

		m_mapping.emplace_back();
		long &cellMapping = m_mapping.back();
		cellMapping = flatId;
	}
}

/*!
	Updates the cell flat mapping.

	\param adaptionData is adaption data that will be used to update
	the mapping
*/
void CellFlatMapping::update(const std::vector<adaption::Info> &adaptionData)
{
	// Previous number of cells
	long nPreviousCells = m_numbering.size();

	// Current number of cells
	long nCurrentCells = m_patch->getCellCount();

	// Find the first change
	long firstChangedFlatId = std::min(nCurrentCells, nPreviousCells);
	long firstChangedId     = m_patch->getCells().getSizeMarker(firstChangedFlatId, Element::NULL_ID);
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != adaption::ENTITY_CELL) {
			continue;
		}

		for (const auto &currentId : adaptionInfo.current) {
			long currentFlatId = m_patch->getCells().evalFlatIndex(currentId);
			if (currentFlatId < firstChangedFlatId) {
				firstChangedId     = currentId;
				firstChangedFlatId = currentFlatId;
			}
		}
	}

	// The mapping, until the flat id with the first change, is a no-op mapping
	m_mapping.resize(nCurrentCells);
	for (long flatId = 0; flatId < firstChangedFlatId; ++flatId) {
		m_mapping[flatId] = flatId;
	}

	// If there are no changes the mapping is already updated
	if (firstChangedId < 0) {
		return;
	}

	// List of cells that have been deleted
	std::unordered_set<long> removedIds;
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != adaption::ENTITY_CELL) {
			continue;
		}

		for (const auto &previousId : adaptionInfo.previous) {
			removedIds.insert(previousId);
		}
	}
	long nRemovedCells = removedIds.size();

	// Create a mapping between cell ids and previous flat ids
	std::unordered_map<long, long> previousFlatIds;
	std::unordered_map<long, long> removedFlatIds;
	previousFlatIds.reserve(nCurrentCells - firstChangedFlatId - nRemovedCells);
	removedFlatIds.reserve(nRemovedCells);
	for (long previousFlatId = firstChangedFlatId; previousFlatId < nPreviousCells; ++previousFlatId) {
		long previousId = m_numbering[previousFlatId];
		if (removedIds.count(previousId) > 0) {
			removedFlatIds.insert({{previousId, previousFlatId}});
		} else {
			previousFlatIds.insert({{previousId, previousFlatId}});
		}
	}

	// Add to the mapping the added cells
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != adaption::ENTITY_CELL) {
			continue;
		}

		// Ancestor flat index of the added cells
		long ancestorFlatId;
		if (adaptionInfo.previous.size() > 0) {
			ancestorFlatId = removedFlatIds.at(adaptionInfo.previous[0]);
		} else {
			ancestorFlatId = -1;
		}

		// Mapping of the added cells
		for (const auto &currentId : adaptionInfo.current) {
			previousFlatIds.insert({{currentId, ancestorFlatId}});
		}
	}

	std::unordered_map<long, long>().swap(removedFlatIds);

	// Update numbering and mapping past the flat id with the first change
	m_numbering.resize(nCurrentCells);
	auto cellIterator = m_patch->getCells().find(firstChangedId);
	for (long flatId = firstChangedFlatId; flatId < nCurrentCells; ++flatId) {
		long cellId = cellIterator->getId();

		m_numbering[flatId] = cellId;
		m_mapping[flatId]   = previousFlatIds[cellId];

		cellIterator++;
	}
}

}
