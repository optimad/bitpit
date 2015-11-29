//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*! \file */

#include <unordered_map>

#include "patch.hpp"

namespace pman {

FlatMapping::FlatMapping()
	: m_patch(nullptr)
{
}

FlatMapping::FlatMapping(Patch *patch)
	: m_patch(patch)
{
}

FlatMapping::~FlatMapping()
{
}

const std::vector<long> & FlatMapping::get_numbering() const
{
	return m_numbering;
}

const std::vector<long> & FlatMapping::get_mapping() const
{
	return m_mapping;
}

CellFlatMapping::CellFlatMapping()
	: FlatMapping()
{
}

CellFlatMapping::CellFlatMapping(Patch *patch)
	: FlatMapping(patch)
{
	m_numbering.reserve(m_patch->get_cell_count());
	m_mapping.reserve(m_patch->get_cell_count());

	long flatId = -1;
	for (const auto &cell : m_patch->cells()) {
		flatId++;

		m_numbering.emplace_back();
		long &cellId = m_numbering.back();
		cellId = cell.get_id();

		m_mapping.emplace_back();
		long &cellMapping = m_mapping.back();
		cellMapping = flatId;
	}
}

CellFlatMapping::~CellFlatMapping()
{
}

void CellFlatMapping::update(const std::vector<Adaption::Info> adaptionData)
{
	// Previous number of cells
	long nPreviousCells = m_numbering.size();

	// Current number of cells
	long nCurrentCells = m_patch->get_cell_count();

	// Map for renumbering the elements
	std::unordered_map<long, long> backwardMap;
	for (auto &adaptionInfo : adaptionData) {
		if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
			continue;
		}

		if (adaptionInfo.type == Adaption::TYPE_REFINEMENT) {
			long previousId = adaptionInfo.previous[0];
			for (const auto &currentId : adaptionInfo.current) {
				backwardMap.insert({{currentId, previousId}});
			}
		} else if (adaptionInfo.type == Adaption::TYPE_COARSENING) {
			long previousId = adaptionInfo.previous[0];
			long currentId  = adaptionInfo.current[0];
			backwardMap.insert({{currentId, previousId}});
		} else if (adaptionInfo.type == Adaption::TYPE_CREATION) {
			long previousId = nPreviousCells;
			long currentId  = adaptionInfo.current[0];
			backwardMap.insert({{currentId, previousId}});
		}
	}

	// Update the mapping up to the first change
	m_mapping.resize(nCurrentCells);
	m_numbering.resize(nCurrentCells);

	long flatId = -1;
	long firstChangedFlatId = -1;
	auto cellIterator = m_patch->cells().cbegin();
	while (cellIterator != m_patch->cells().cend()) {
		flatId++;
		if (flatId == nCurrentCells) {
			break;
		}

		long cellId = cellIterator->get_id();
		if (cellId == m_numbering[flatId]) {
			m_mapping[flatId] = flatId;
			cellIterator++;
			continue;
		}

		firstChangedFlatId = flatId;
		break;
	}

	// If there are no changes the mapping is updated
	if (firstChangedFlatId < 0 && nCurrentCells <= nPreviousCells) {
		return;
	}

	// Build a map for the previous numbering
	//
	// We only need the flat ids after the flat id with the first change.
	std::unordered_map<long, long> previousNumberingMap;
	for (long n = firstChangedFlatId; n < nPreviousCells; ++n) {
		previousNumberingMap.insert({{m_numbering[n], n}});
	}

	// Continue the update of the mapping
	flatId = firstChangedFlatId - 1;
	while (cellIterator != m_patch->cells().cend()) {
		flatId++;
		long currentId = cellIterator->get_id();

		long previousId;
		if (backwardMap.count(currentId) != 0) {
			previousId = backwardMap.at(currentId);
		} else {
			previousId = currentId;
		}

		long previousFlatId;
		if (previousId == m_numbering[flatId]) {
			previousFlatId = flatId;
		} else {
			long previousId;
			if (backwardMap.count(currentId) != 0) {
				previousId = backwardMap.at(currentId);
			} else {
				previousId = currentId;
			}

			if (previousNumberingMap.count(previousId) != 0) {
				previousFlatId = previousNumberingMap.at(previousId);
			} else {
				previousFlatId = Element::NULL_ELEMENT_ID;
			}
		}

		m_numbering[flatId] = currentId;
		m_mapping[flatId]   = previousFlatId;

		cellIterator++;
	}
}

}
