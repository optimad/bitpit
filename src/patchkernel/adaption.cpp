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

	\brief The Info struct defines the information associated to an
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
	m_cachedTypes.insert(adaption::TYPE_DELETION);
	m_cachedTypes.insert(adaption::TYPE_CREATION);
	m_cachedTypes.insert(adaption::TYPE_PARTITION_RECV);
	m_cachedTypes.insert(adaption::TYPE_PARTITION_SEND);
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
	Creates an empty adaption info.

	\result An empty adaption info.
*/
std::size_t InfoCollection::create()
{
	m_collection.emplace_back();

	return m_collection.size() - 1;
}

/*!
	Creates an adaption info with the requested data.

	If an adaption info with the requested data already exists, the existing
	value will be returned, otherwise a new	adaption info will be created.

	\param type is the type of adaption info
	\param entity is the entity associated to the adaption info
	\param rank is the rank associated to the adaption info
	\result The requested adaption info. If an adaption info with the
	requested data already exists, the existing value will be returned,
	otherwise a new	adaption info will be created.
*/
std::size_t InfoCollection::create(Type type, Entity entity, int rank)
{
	infoData_t infoData = infoData_t(type, entity, rank);
	bool useCache = (m_cachedTypes.count(type) > 0);

	if (useCache) {
		if (m_cache.count(infoData) > 0) {;
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
	Returns a reference to the adaption info at the specified position in
	the collection.

	\param n is the position of an adaption info in the collection.
	\result Returns a reference to requested adaption info.
*/
Info & InfoCollection::at(std::size_t n)
{
	if (n >= m_collection.size()) {
		throw std::out_of_range("Requested adaption info is not in the collection");
	}

	return (*this)[n];
}

/*!
	Returns a constant reference to the adaption info at the specified position
	in the collection.

	\param n is the position of an adaption info in the collection.
	\result Returns a constant reference to requested adaption info.
*/
const Info & InfoCollection::at(std::size_t n) const
{
	if (n >= m_collection.size()) {
		throw std::out_of_range("Requested adaption info is not in the collection");
	}

	return (*this)[n];
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
	Default destructor.
*/
FlatMapping::~FlatMapping()
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
	Default destructor.
*/
CellFlatMapping::~CellFlatMapping()
{
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
