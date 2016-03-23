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

#include <unordered_map>
#include <unordered_set>

#include "patch_kernel.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\struct Adaption

	\brief The Adaption struct defines the information associated to an
	adaption.
*/


/*!
	\enum Adaption::Type

	\brief The Type enum defines the type of adaption that has been
	performed.
*/

/*!
	\enum Adaption::Entity

	\brief The Entity enum defines the type of entities on which the
	adaption has been performed.
*/

/*!
	\struct Info

	\brief The Info struct defines the information associated to an
	adaption.
*/

/*!
	@}
*/

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class FlatMapping

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
	@}
*/

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class CellFlatMapping

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
void CellFlatMapping::update(const std::vector<Adaption::Info> adaptionData)
{
	// Previous number of cells
	long nPreviousCells = m_numbering.size();

	// Current number of cells
	long nCurrentCells = m_patch->getCellCount();

	// Resize the data structures
	m_mapping.resize(nCurrentCells);
	m_numbering.resize(nCurrentCells);

	// Initialize the mapping
	for (long flatId = 0; flatId < nCurrentCells; ++flatId) {
		m_mapping[flatId] = flatId;
	}

	// Find the first change
	long firstChangedFlatId = std::min(nCurrentCells, nPreviousCells);
	long firstChangedId     = m_patch->getCells().getSizeMarker(firstChangedFlatId, Element::NULL_ID);
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
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

	// If there are no changes the mapping is already updated
	if (firstChangedId < 0) {
		return;
	}

	// Build a map for the previous numbering
	std::unordered_set<long> previousIds;
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
			continue;
		}

		previousIds.insert(adaptionInfo.previous[0]);
	}

	std::unordered_map<long, long> previousNumberingMap;
	previousNumberingMap.reserve(nPreviousCells - firstChangedFlatId);
	for (long flatId = 0; flatId < nPreviousCells; ++flatId) {
		long previousId = m_numbering[flatId];
		if (previousIds.count(previousId) > 0) {
			previousNumberingMap.insert({{previousId, flatId}});
		}
	}

	// Update the current numbering
	//
	// We only need to update the numbering after the flat id with the first
	// change.
	auto cellIterator = m_patch->getCells().getConstIterator(firstChangedId);
	for (long flatId = firstChangedFlatId; flatId < nCurrentCells; ++flatId) {
		long cellId = cellIterator->getId();
		m_numbering[flatId] = cellId;
		cellIterator++;
	}

	// Update the mapping
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
			continue;
		}

		// Id of the ancestor
		long previousId = adaptionInfo.previous[0];

		// Flat id previously associated to the ancestor
		long previousFlatId;
		if (previousNumberingMap.count(previousId) != 0) {
			previousFlatId = previousNumberingMap.at(previousId);
		} else {
			previousFlatId = -1;
		}

		// Update the mapping for all the current cells
		for (const auto &currentId : adaptionInfo.current) {
			// Flat id associated to the current cell
			long currentFlatId = m_patch->getCells().evalFlatIndex(currentId);

			// Mapping between the two flat ids
			m_mapping[currentFlatId] = previousFlatId;
		}
	}
}

}
