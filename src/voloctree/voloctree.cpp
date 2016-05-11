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

#include "logger.hpp"

#include "voloctree.hpp"

namespace bitpit {

/*!
	\ingroup voloctree
	\class OctreeLevelInfo

	\brief The OctreeLevelInfo class defines the information associated to
	an octree level.
*/

/*!
	\ingroup voloctree
	@{
*/

/*!
	\class VolOctree

	\brief The VolOctree defines a Octree patch.

	VolOctree defines a Octree patch.
*/

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
*/
VolOctree::VolOctree(const int &id, const int &dimension,
				 std::array<double, 3> origin, double length, double dh )
	: VolumeKernel(id, dimension, false),
	  m_tree(origin[0], origin[1], origin[2], length, dimension,
	         PabloUniform::DEFAULT_MAX_LEVELS, PabloUniform::DEFAULT_LOG_FILE
#if BITPIT_ENABLE_MPI==1
	         , MPI_COMM_NULL
#endif
	        ),
	  m_lastTreeOperation(OP_INITIALIZATION)
{
	log::cout() << ">> Initializing Octree mesh\n";

	// Initialize local edges/vertex/faces association
	if (getDimension() == 3) {
		m_octantLocalFacesOnVertex.reserve(8);
		m_octantLocalFacesOnVertex.push_back({{0, 2, 4}});
		m_octantLocalFacesOnVertex.push_back({{1, 2, 4}});
		m_octantLocalFacesOnVertex.push_back({{0, 3, 4}});
		m_octantLocalFacesOnVertex.push_back({{1, 3, 4}});
		m_octantLocalFacesOnVertex.push_back({{0, 2, 5}});
		m_octantLocalFacesOnVertex.push_back({{1, 2, 5}});
		m_octantLocalFacesOnVertex.push_back({{0, 3, 5}});
		m_octantLocalFacesOnVertex.push_back({{1, 3, 5}});

		m_octantLocalEdgesOnVertex.reserve(8);
		m_octantLocalEdgesOnVertex.push_back({{0, 2,  4}});
		m_octantLocalEdgesOnVertex.push_back({{1, 2,  5}});
		m_octantLocalEdgesOnVertex.push_back({{0, 3,  6}});
		m_octantLocalEdgesOnVertex.push_back({{1, 3,  7}});
		m_octantLocalEdgesOnVertex.push_back({{4, 8, 10}});
		m_octantLocalEdgesOnVertex.push_back({{5, 9, 10}});
		m_octantLocalEdgesOnVertex.push_back({{6, 8, 11}});
		m_octantLocalEdgesOnVertex.push_back({{7, 9, 11}});

		m_octantLocalFacesOnEdge.reserve(12);
		m_octantLocalFacesOnEdge.push_back({{0, 4}});
		m_octantLocalFacesOnEdge.push_back({{1, 4}});
		m_octantLocalFacesOnEdge.push_back({{2, 4}});
		m_octantLocalFacesOnEdge.push_back({{3, 4}});
		m_octantLocalFacesOnEdge.push_back({{0, 2}});
		m_octantLocalFacesOnEdge.push_back({{1, 2}});
		m_octantLocalFacesOnEdge.push_back({{0, 3}});
		m_octantLocalFacesOnEdge.push_back({{1, 3}});
		m_octantLocalFacesOnEdge.push_back({{0, 5}});
		m_octantLocalFacesOnEdge.push_back({{1, 5}});
		m_octantLocalFacesOnEdge.push_back({{2, 5}});
		m_octantLocalFacesOnEdge.push_back({{3, 5}});
	} else {
		m_octantLocalFacesOnVertex.reserve(4);
		m_octantLocalFacesOnVertex.push_back({{0, 2}});
		m_octantLocalFacesOnVertex.push_back({{1, 2}});
		m_octantLocalFacesOnVertex.push_back({{0, 3}});
		m_octantLocalFacesOnVertex.push_back({{1, 3}});
	}

	// Inizializzazione dell'octree
	double initial_level = ceil(log2(std::max(1., length / dh)));

	m_tree.setMarker((uint32_t) 0, initial_level);

	// Info sull'octree
	initializeTreeGeometry();

	// Info of the cell type
	ElementInfo::Type cellType;
	if (isThreeDimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	m_cellTypeInfo = &ElementInfo::getElementInfo(cellType);

	// Info on the interface type
	ElementInfo::Type interfaceType;
	if (isThreeDimensional()) {
		interfaceType = ElementInfo::PIXEL;
	} else {
		interfaceType = ElementInfo::LINE;
	}

	m_interfaceTypeInfo = &ElementInfo::getElementInfo(interfaceType);

	// Info sulle interfacce
	for (int i = 0; i < dimension; i++) {
		for (int n = -1; n <= 1; n += 2) {
			std::array<double, 3> normal = {{0.0, 0.0, 0.0}};
			normal[i] = n;

			m_normals.push_back(normal);
		}
	}
}

/*!
	Destroys the patch.
*/
VolOctree::~VolOctree()
{

}

/*!
	Initializes octree geometry.
*/
void VolOctree::initializeTreeGeometry()
{
	int maxLevels = m_tree.getMaxLevel();
	double length = m_tree.getL();

	m_tree_dh.clear();
	m_tree_area.clear();
	m_tree_volume.clear();

	m_tree_dh.reserve(maxLevels);
	m_tree_area.reserve(maxLevels);
	m_tree_volume.reserve(maxLevels);
	for(int i = 0; i < maxLevels; i++) {
	    double levelLength = length / ((double) pow(2,i));

	    m_tree_dh.push_back(pow(levelLength, 1.));
	    m_tree_area.push_back(pow(levelLength, (double) (getDimension() - 1)));
	    m_tree_volume.push_back(pow(levelLength, (double) (getDimension())));
	}
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolOctree::evalCellVolume(const long &id)
{
	int level = getCellLevel(id);

	return m_tree_volume[level];
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> VolOctree::evalCellCentroid(const long &id)
{
	OctantInfo octantInfo = getCellOctant(id);

	Octant *octant;
	if (octantInfo.internal) {
		octant = m_tree.getOctant(octantInfo.id);
	} else {
		octant = m_tree.getGhostOctant(octantInfo.id);
	}

	return m_tree.getCenter(octant);
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolOctree::evalCellSize(const long &id)
{
	int level = getCellLevel(id);

	return m_tree_dh[level];
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double VolOctree::evalInterfaceArea(const long &id)
{
	const Interface &interface = getInterface(id);
	int owner = interface.getOwner();
	int level = getCellLevel(owner);

	return m_tree_area[level];
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> VolOctree::evalInterfaceNormal(const long &id)
{
	const Interface &interface = getInterface(id);
	int ownerFace = interface.getOwnerFace();

	return m_normals[ownerFace];
}

/*!
	Gets the octant of the cell with the specified id.

	\param id the id of the cell
	\result The octant info of the specified cell
*/
VolOctree::OctantInfo VolOctree::getCellOctant(const long &id) const
{
	OctantInfo octantInfo;
	octantInfo.internal = m_cells[id].isInterior();
	if (octantInfo.internal) {
		octantInfo.id = m_cellToOctant.at(id);
	} else {
		octantInfo.id = m_cellToGhost.at(id);
	}

	return octantInfo;
}

/*!
	\brief Gets a reference to the octree associated with the patch.

	\result A reference to the octree associated to the patch.
*/
PabloUniform & VolOctree::getTree()
{
	return m_tree;
}

/*!
	Gets the id of the specified octant.

	\param octantInfo the data of the octant
	\result The id of the specified octant
*/
long VolOctree::getOctantId(const OctantInfo &octantInfo) const
{
	std::unordered_map<uint32_t, long>::const_iterator octantItr;
	if (octantInfo.internal) {
		octantItr = m_octantToCell.find(octantInfo.id);
		if (octantItr == m_octantToCell.end()) {
			return Element::NULL_ID;
		}
	} else {
		octantItr = m_ghostToCell.find(octantInfo.id);
		if (octantItr == m_ghostToCell.end()) {
			return Element::NULL_ID;
		}
	}

	return octantItr->second;
}

/*!
	Gets the connectivity of the specified octant.

	\param octantInfo the data of the octant
	\result A reference to the octant's connectivity
*/
const std::vector<uint32_t> & VolOctree::getOctantConnect(const OctantInfo &octantInfo)
{
	if (octantInfo.internal) {
		return m_tree.getConnectivity()[octantInfo.id];
	} else {
		return m_tree.getGhostConnectivity()[octantInfo.id];
	}
}

/*!
	Evaluates a unique hash for the octant.

	\param octantInfo the data of the octant
	\result A unique hash for the octant.
*/
VolOctree::OctantHash VolOctree::evaluateOctantHash(const OctantInfo &octantInfo)
{
	uint8_t level   = m_tree.getLevel(octantInfo.id);
	uint64_t morton = m_tree.getMorton(octantInfo.id);

	OctantHash octantHash;
	octantHash |= morton;
	octantHash <<= 8;
	octantHash |= level;

	return octantHash;
}

/*!
	Gets the refinement level of the cell with the specified id.

	\param id is the id of the cell
	\result The refinement level of the specified cell.
*/
int VolOctree::getCellLevel(const long &id)
{
	OctantInfo octantInfo = getCellOctant(id);

	Octant* octant;
	if (octantInfo.internal) {
		octant = m_tree.getOctant(octantInfo.id);
	} else {
		octant = m_tree.getGhostOctant(octantInfo.id);
	}
	return m_tree.getLevel(octant);
}

/*!
	Updates the patch.

	\param trackAdaption if set to true the changes to the patch will be
	tracked
	\result Returns all the changes applied to the patch.
*/
const std::vector<adaption::Info> VolOctree::_updateAdaption(bool trackAdaption)
{

	// Updating the tree
	log::cout() << ">> Adapting tree...";

	bool buildMapping = (getCellCount() != 0);
	bool updated = m_tree.adapt(buildMapping);
	if (trackAdaption) {
		m_lastTreeOperation = OP_ADAPTION_MAPPED;
	} else {
		m_lastTreeOperation = OP_ADAPTION_UNMAPPED;
	}

	if (!updated) {
		log::cout() << " Already updated" << std::endl;
		return std::vector<adaption::Info>();
	}
	log::cout() << " Done" << std::endl;

	// Sync the patch
	return sync(trackAdaption);
}

/*!
	Syncronizes the patch with the underlying octree.

	\param trackChanges if set to true the changes to the patch will be
	tracked
	\result Returns all the changes applied to the patch.
*/
const std::vector<adaption::Info> VolOctree::sync(bool trackChanges)
{
	log::cout() << ">> Syncing patch..." << std::endl;

	// If the current mesh is empty we need to import all the octants
	bool importAll = (getCellCount() == 0);

	// Last operation on the tree
	TreeOperation lastTreeOperation = m_lastTreeOperation;
	if (lastTreeOperation == OP_ADAPTION_UNMAPPED && !importAll) {
		throw std::runtime_error ("Unable to sync the patch after an unmapped adaption");
	}

	// Info on the tree
	long nOctants = m_tree.getNumOctants();
	long nPreviousOctants = m_octantToCell.size();

	log::cout() << ">> Number of octants : " << nOctants << std::endl;

	// Info on the tree
	long nGhostsOctants = m_tree.getNumGhosts();
	long nPreviousGhosts = m_ghostToCell.size();

	// Evaluate tree conenctivity
	log::cout() << ">> Evaluating Octree connectivity...";

	m_tree.computeConnectivity();

	log::cout() << " Done" << std::endl;

	// Initialize tracking data
	adaption::InfoCollection adaptionData;

	// Extract information for transforming the patch
	//
	// If there are no cells in the mesh we need to import all
	// octants.
	log::cout() << ">> Extract information for transforming the patch...";

	std::vector<OctantInfo> addedOctants;
	std::vector<RenumberInfo> renumberedOctants;
	std::unordered_set<long> removedCells;

	addedOctants.reserve(nOctants + nGhostsOctants);
	renumberedOctants.reserve(nPreviousOctants + nPreviousGhosts);
	removedCells.reserve(nPreviousOctants + nPreviousGhosts);

	uint32_t treeId = 0;
	while (treeId < (uint32_t) nOctants) {
		// Octant mapping
		std::vector<uint32_t> mapper_octantMap;
		std::vector<bool> mapper_ghostFlag;
		std::vector<int> mapper_octantRank;
		if (!importAll) {
			m_tree.getMapping(treeId, mapper_octantMap, mapper_ghostFlag, mapper_octantRank);
		}

		// Adaption type
		adaption::Type adaptionType = adaption::TYPE_NONE;
		if (importAll) {
			adaptionType = adaption::TYPE_CREATION;
		} else if (lastTreeOperation == OP_ADAPTION_MAPPED) {
			bool isNewR = m_tree.getIsNewR(treeId);
			if (isNewR) {
				adaptionType = adaption::TYPE_REFINEMENT;
			} else {
				bool isNewC = m_tree.getIsNewC(treeId);
				if (isNewC) {
					adaptionType = adaption::TYPE_COARSENING;
				} else if (treeId != mapper_octantMap.front()) {
					adaptionType = adaption::TYPE_RENUMBERING;
				}
			}
#if BITPIT_ENABLE_MPI==1
		} else if (lastTreeOperation == OP_LOAD_BALANCE) {
			if (getRank() != mapper_octantRank.front()) {
				adaptionType = adaption::TYPE_PARTITION_RECV;
			} else if (treeId != mapper_octantMap.front()) {
				adaptionType = adaption::TYPE_RENUMBERING;
			}
#endif
		}

		// If the octant cell has not been modified we can skip to the next
		// octant.
		if (adaptionType == adaption::TYPE_NONE) {
			++treeId;
			continue;
		}

		// Re-numbered cells just need to be added to the proper list.
		//
		// Renumbered cells are not tracked, because the re-numbering
		// only happens inside VolOctree.
		if (adaptionType == adaption::TYPE_RENUMBERING) {
			uint32_t previousTreeId = mapper_octantMap.front();
			OctantInfo previousOctantInfo(previousTreeId, !mapper_ghostFlag.front());
			renumberedOctants.emplace_back(previousOctantInfo, treeId);

			// No more work needed, skip to the next octant
			++treeId;
			continue;
		}

		// Handle other kind of adaption
		//
		// New octants need to be imported into the patch,
		// whereas cells associated to previous octants
		// need to be removed.
		//
		// If the user want to track adaption, adaption
		// data needs to be filled.

		// Current tree ids that will be imported
		long nCurrentTreeIds;
		if (importAll) {
			nCurrentTreeIds = nOctants - treeId;
		} else if (adaptionType == adaption::TYPE_REFINEMENT) {
			nCurrentTreeIds = pow(2, getDimension());
		} else {
			nCurrentTreeIds = 1;
		}

		const long lastCurrentTreeId = treeId + nCurrentTreeIds;
		for (int currentTreeId = treeId; currentTreeId < lastCurrentTreeId; ++currentTreeId) {
			addedOctants.emplace_back(currentTreeId, true);
		}

		// Cells that will be removed
		//
		// Mark the cells associated to previous local octants for deletion.
		int nPreviousTreeIds = mapper_octantMap.size();
		for (int k = 0; k < nPreviousTreeIds; ++k) {
#if BITPIT_ENABLE_MPI==1
			// Only local cells can be deleted
			if (mapper_octantRank[k] != getRank()) {
				continue;
			}
#endif

			// Mark the cell for deletion
			OctantInfo previousOctantInfo(mapper_octantMap[k], !mapper_ghostFlag[k]);
			long previousCellId = getOctantId(previousOctantInfo);

			removedCells.insert(previousCellId);
		}

		// Adaption tracking
		//
		// The adaption info associated to the octants that has been received
		// from external partitions will contain the current octants sorted by
		// their tree id (we are looping over the octants in that order), this
		// is the same order that will be used on the processor that has sent
		// the octants. Since the order is the same, the two processors are able
		// to exchange cell data without any additional extra communication
		// (they already know the list of cells for which data is needed and
		// the order in which these data will be sent).
		if (trackChanges) {
			// Rank assocated to the adaption info
			int rank = -1;
#if BITPIT_ENABLE_MPI==1
			if (adaptionType == adaption::TYPE_PARTITION_RECV) {
				rank = mapper_octantRank[0];
			} else {
				rank = getRank();
			}
#endif

			// Get the adaption info
			std::size_t infoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, rank);
			adaption::Info &adaptionInfo = adaptionData[infoId];

			// Current status
			//
			// We don't know the id of the current status, because those
			// cells are not yet in the mesh. Store the trre id, and
			// make the translation later.
			//
			// WARNING: tree id are uint32_t wherase adaptionInfo stores
			//          id as unsigned long.
			adaptionInfo.current.reserve(nCurrentTreeIds);
			auto addedOctantsIter = addedOctants.cend() - nCurrentTreeIds;
			while (addedOctantsIter != addedOctants.cend()) {
				adaptionInfo.current.emplace_back();
				unsigned long &adaptionId = adaptionInfo.current.back();
				adaptionId = (*addedOctantsIter).id;

				addedOctantsIter++;
			}

			// Previous cells
			//
			// A coarsening can merge togheter cells of different processors.
			// However, since the coarsening is limited to one level, the
			// previous cells will always be internal or among the ghost of
			// the current processor.
			int nPreviousCellIds = mapper_octantMap.size();
			adaptionInfo.previous.reserve(nPreviousCellIds);
			for (int k = 0; k < nPreviousCellIds; ++k) {
				long previousCellId;
#if BITPIT_ENABLE_MPI==1
				if (mapper_octantRank[k] != getRank()) {
					previousCellId = Cell::NULL_ID;
				} else
#endif
				{
					OctantInfo previousOctantInfo(mapper_octantMap[k], !mapper_ghostFlag[k]);
					previousCellId = getOctantId(previousOctantInfo);
				}

				adaptionInfo.previous.emplace_back();
				unsigned long &adaptionId = adaptionInfo.previous.back();
				adaptionId = previousCellId;
			}
		}

		// Incremente tree id
		treeId += nCurrentTreeIds;
	}

	log::cout() << " Done" << std::endl;

#if BITPIT_ENABLE_MPI==1
	// Cells that have been send to other processors need to be removed
	std::unordered_map<int, std::array<uint32_t, 4>> sendOctants = m_tree.getSentIdx();
	for (const auto &rankEntry : sendOctants) {
		int rank = rankEntry.first;

		// If needed create the adaption info
		std::size_t adaptionInfoId = -1;
		if (trackChanges) {
			adaption::Type adaptionType;
			if (rank >= 0) {
				adaptionType = adaption::TYPE_PARTITION_SEND;
			} else {
				adaptionType = adaption::TYPE_DELETION;
			}

			adaptionInfoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, rank);
		}

		// Add the send cells to the list of cells to be removed
		//
		// The ids will be stored sorted by the id of the octants, this is the
		// same order that will be used on the processor that has received the
		// octants. Since the order is the same, the two processors are able
		// to exchange cell data without any additional extra communication
		// (they already know the list of cells for which data is needed and
		// the order in which these data will be sent).
		for (int k = 0; k < 2; ++k) {
			uint32_t beginTreeId = rankEntry.second[2 * k];
			uint32_t endTreeId   = rankEntry.second[2 * k + 1];
			for (uint32_t treeId = beginTreeId; treeId < endTreeId; ++treeId) {
				OctantInfo octantInfo(treeId, true);
				long cellId = getOctantId(octantInfo);
				removedCells.insert(cellId);

				if (trackChanges) {
					adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
					adaptionInfo.previous.emplace_back();
					unsigned long &sentId = adaptionInfo.previous.back();
					sentId = cellId;
				}
			}
		}
	}

	// Previous ghosts cells need to be removed
	if (nPreviousGhosts > 0) {
		clearGhostOwners(true);

		std::size_t adaptionInfoId = -1;
		if (trackChanges) {
			adaptionInfoId = adaptionData.create(adaption::TYPE_DELETION, adaption::ENTITY_CELL, getRank());
		}

		auto cellIterator = m_cellToGhost.cbegin();
		while (cellIterator != m_cellToGhost.cend()) {
			long id = cellIterator->first;
			removedCells.insert(id);

			// Adaption tracking
			if (trackChanges) {
				adaption::Info &deletedGhostsInfo = adaptionData[adaptionInfoId];
				deletedGhostsInfo.previous.emplace_back();
				unsigned long &adaptionId = deletedGhostsInfo.previous.back();
				adaptionId = id;
			}

			// Increment iterator
			cellIterator++;
		}
	}

	// New ghost octants need to be added
	for (uint32_t treeId = 0; treeId < (uint32_t) nGhostsOctants; ++treeId) {
		addedOctants.emplace_back(treeId, false);
	}
	addedOctants.shrink_to_fit();
#endif

	// Enable advanced editing
	setExpert(true);

	// Delete removed cells
	FaceInfoSet danglingFaces;
	if (removedCells.size() > 0) {
		log::cout() << ">> Removing non-existing cells...";

		// Track deleted interfaces
		if (trackChanges) {
			// List of unique interfaces that will be deleted
			std::unordered_set<long> removedInterfaces;
			for (const auto &cellId : removedCells) {
				const Cell &cell = m_cells.at(cellId);
				long nCellInterfaces = cell.getInterfaceCount();
				const long *interfaces = cell.getInterfaces();
				for (int k = 0; k < nCellInterfaces; ++k) {
					long interfaceId = interfaces[k];
					if (interfaceId >= 0) {
						removedInterfaces.insert(interfaceId);
					}
				}
			}

			// Rank assocated to the adaption info
			int rank = -1;
#if BITPIT_ENABLE_MPI==1
			rank = getRank();
#endif

			// Adaption info
			std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_DELETION, adaption::ENTITY_INTERFACE, rank);
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
			for (const long &interfaceId : removedInterfaces) {
				adaptionInfo.previous.emplace_back();
				unsigned long &deletedInterfaceId = adaptionInfo.previous.back();
				deletedInterfaceId = interfaceId;
			}
		}

		// Delete cells
		danglingFaces = removeCells(removedCells);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Cells removed: " <<  removedCells.size() << std::endl;
	}

	std::unordered_set<long>().swap(removedCells);

	// Reserve space in the octant-to-cell maps
	m_cellToOctant.reserve(nOctants);
	m_octantToCell.reserve(nOctants);

	// Remap renumbered cells
	if (renumberedOctants.size() > 0) {
		log::cout() << ">> Rebuilding octant-to-cell map for renumbered cells...";

		renumberOctants(renumberedOctants);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Cells renumbered: " <<  renumberedOctants.size() << std::endl;
	}

	std::vector<RenumberInfo>().swap(renumberedOctants);

	// Reset ghost maps
	m_cellToGhost.clear();
	m_cellToGhost.reserve(nGhostsOctants);

	m_ghostToCell.clear();
	m_ghostToCell.reserve(nGhostsOctants);

	// Import added octants
	std::vector<long> createdCells;
	if (addedOctants.size() > 0) {
		log::cout() << ">> Importing new octants...";

		createdCells = importOctants(addedOctants, danglingFaces);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Octants imported: " <<  addedOctants.size() << std::endl;
	}

	FaceInfoSet().swap(danglingFaces);

	// Rebuild the ghost information
#if BITPIT_ENABLE_MPI==1
	buildGhostExchangeData();
#endif

	// Disable advanced editing
	setExpert(false);

	// Track mesh adaption
	if (trackChanges) {
		// Map ids of the added cells
		for (auto &adaptionInfo : adaptionData.data()) {
			if (adaptionInfo.entity != adaption::ENTITY_CELL) {
				continue;
			}

			int nCurrentIds = adaptionInfo.current.size();
			for (int k = 0; k < nCurrentIds; ++k) {
				long cellId = m_octantToCell.at(adaptionInfo.current[k]);
				adaptionInfo.current[k] = cellId;
			}
		}

		// Track created ghosts cells
#if BITPIT_ENABLE_MPI==1
		if (nGhostsOctants > 0) {
			std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_CREATION, adaption::ENTITY_CELL, getRank());
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];

			adaptionInfo.current.reserve(nGhostsOctants);
			auto cellIterator = m_cellToGhost.cbegin();
			while (cellIterator != m_cellToGhost.cend()) {
				adaptionInfo.current.emplace_back();
				unsigned long &adaptionId = adaptionInfo.current.back();
				adaptionId = cellIterator->first;

				cellIterator++;
			}
		}
#endif

		// Track created interfaces
		if (createdCells.size() > 0) {
			// List of unique interfaces that have been created
			std::unordered_set<long> createdInterfaces;
			for (const auto &cellId : createdCells) {
				const Cell &cell = m_cells.at(cellId);
				long nCellInterfaces = cell.getInterfaceCount();
				const long *interfaces = cell.getInterfaces();
				for (int k = 0; k < nCellInterfaces; ++k) {
					long interfaceId = interfaces[k];
					if (interfaceId >= 0) {
						createdInterfaces.insert(interfaceId);
					}
				}
			}

			// Rank assocated to the adaption info
			int rank = -1;
#if BITPIT_ENABLE_MPI==1
			rank = getRank();
#endif

			// Adaption info
			std::size_t infoId = adaptionData.create(adaption::TYPE_CREATION, adaption::ENTITY_INTERFACE, rank);
			adaption::Info &adaptionInfo = adaptionData[infoId];
			for (const long &interfaceId : createdInterfaces) {
				adaptionInfo.current.emplace_back();
				unsigned long &createdInterfaceId = adaptionInfo.current.back();
				createdInterfaceId = interfaceId;
			}
		}
	}

	// Delete tree conenctivity
	m_tree.clearConnectivity();

	// Done
	return adaptionData.dump();
}

/*!
	Renumber a list of octants.

	\param renumberedOctants is the list of octant to renumber
*/
void VolOctree::renumberOctants(std::vector<RenumberInfo> &renumberedOctants)
{
	// Remove previous cell-to-tree associations
	std::vector<long> cellIds;
	cellIds.reserve(renumberedOctants.size());
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		const OctantInfo &octantInfo = renumberInfo.octantInfo;

		long cellId = getOctantId(octantInfo);
		cellIds.push_back(cellId);

		uint32_t previousTreeId = octantInfo.id;
		m_octantToCell.erase(previousTreeId);
	}

	// Creeate new cell-to-tree associations
	std::vector<long>::iterator cellIdsItr = cellIds.begin();
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		uint32_t treeId = renumberInfo.newTreeId;

		long cellId = (*cellIdsItr);
		cellIdsItr++;

		m_cellToOctant[cellId] = treeId;
		m_octantToCell[treeId] = cellId;
	}
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
std::vector<long> VolOctree::importOctants(std::vector<OctantInfo> &octantInfoList)
{
	FaceInfoSet danglingFaces;

	return importOctants(octantInfoList, danglingFaces);
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
	\param danglingVertices is the list of dangling vertices in the current mesh
*/
std::vector<long> VolOctree::importOctants(std::vector<OctantInfo> &octantInfoList,
                                           FaceInfoSet &danglingFaces)
{
	// Info of the cells
	const int &nCellVertices = m_cellTypeInfo->nVertices;
	const std::vector<std::vector<int>> &cellLocalFaceConnect = m_cellTypeInfo->faceConnect;

	// Info on the interfaces
	const int &nInterfaceVertices = m_interfaceTypeInfo->nVertices;

	// Add the vertex of the dangling faces to the vertex map
	std::unordered_map<uint32_t, long> vertexMap;
	for (auto &danglingFaceInfo : danglingFaces) {
		// List of faces with the vertx to be added
		//
		// We need to consider the dangling face itself and all interfaces
		// on that face.
		long danglingId = danglingFaceInfo.id;
		Cell &danglingCell = m_cells[danglingId];
		int danglingFace = danglingFaceInfo.face;

		int nInterfaces = danglingCell.getInterfaceCount(danglingFace);
		std::vector<FaceInfo> vertexSourceList;
		vertexSourceList.reserve(1 + nInterfaces);

		vertexSourceList.emplace_back(danglingFaceInfo.id, danglingFaceInfo.face);
		for (int k = 0; k < nInterfaces; ++k) {
			long interfaceId = danglingCell.getInterface(danglingFace, k);
			if (interfaceId < 0) {
				continue;
			}

			Interface &interface = m_interfaces[interfaceId];
			if (interface.getOwner() != danglingId) {
				vertexSourceList.emplace_back(interface.getOwner(), interface.getOwnerFace());
			} else {
				vertexSourceList.emplace_back(interface.getNeigh(), interface.getNeighFace());
			}
		}

		// Add the vertices to the map
		for (auto & vertexSource : vertexSourceList) {
			// Cell data
			Cell &cell = m_cells[vertexSource.id];
			const long *cellConnect = cell.getConnect();

			// Octant data
			OctantInfo octantInfo = getCellOctant(vertexSource.id);
			const std::vector<uint32_t> &octantTreeConnect = getOctantConnect(octantInfo);

			// List of vertices
			const std::vector<int> &localConnect = cellLocalFaceConnect[vertexSource.face];
			for (int k = 0; k < nInterfaceVertices; ++k) {
				uint32_t vertexTreeId = octantTreeConnect[localConnect[k]];
				if (vertexMap.count(vertexTreeId) == 0) {
					long vertexId = cellConnect[localConnect[k]];
					vertexMap.insert({{vertexTreeId, vertexId}});
				}
			}
		}
	}

	// Create the new vertices
	for (OctantInfo &octantInfo : octantInfoList) {
		const std::vector<uint32_t> &octantTreeConnect = getOctantConnect(octantInfo);
		for (int k = 0; k < nCellVertices; ++k) {
			uint32_t vertexTreeId = octantTreeConnect[k];
			if (vertexMap.count(vertexTreeId) == 0) {
				vertexMap[vertexTreeId] = addVertex(vertexTreeId);
			}
		}
	}

	// Add the cells
	std::vector<long> createdCells;
	createdCells.reserve(octantInfoList.size());

	for (OctantInfo &octantInfo : octantInfoList) {
		// Octant connectivity
		const std::vector<uint32_t> &octantTreeConnect = getOctantConnect(octantInfo);

		// Cell connectivity
		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int k = 0; k < nCellVertices; ++k) {
			uint32_t vertexTreeId = octantTreeConnect[k];
			cellConnect[k] = vertexMap.at(vertexTreeId);
		}

		// Add cell
		long cellId = addCell(octantInfo, cellConnect);
		createdCells.push_back(cellId);
	}

	// Build adjacencies
	updateAdjacencies(createdCells, false);
	updateInterfaces(createdCells, false);

	// Done
	return createdCells;
}

/*!
	Remove a list of octants from the patch.

	\param cellIds is the list of cells ids to remove
*/
VolOctree::FaceInfoSet VolOctree::removeCells(std::unordered_set<long> &cellIds)
{
	FaceInfoSet danglingFaces;

	// Info of the cells
	ElementInfo::Type cellType;
	if (isThreeDimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
	const std::vector<std::vector<int>> &cellLocalFaceConnect = cellTypeInfo.faceConnect;

	// Info on the faces
	ElementInfo::Type faceType;
	if (isThreeDimensional()) {
		faceType = ElementInfo::PIXEL;
	} else {
		faceType = ElementInfo::LINE;
	}

	const ElementInfo &faceTypeInfo = ElementInfo::getElementInfo(faceType);
	const int &nFaceVertices = faceTypeInfo.nVertices;

	// Delete the cells
	std::unordered_set<long> deadVertices;
	std::unordered_set<long> deadInterfaces;
	for (long cellId : cellIds) {
		Cell &cell = m_cells[cellId];

		// List vertices to remove
		//
		// For now, all cell vertices will be listed. Later, the vertex of
		// the dangling faces will be removed from the list.
		int nCellVertices = cell.getVertexCount();
		for (int k = 0; k < nCellVertices; ++k) {
			long vertexId = cell.getVertex(k);
			deadVertices.insert(vertexId);
		}

		// List of interfaces to delete
		//
		// All the interfaces will be deleted, this means that some of the
		// cells that are not deleted will have a face not donnected to
		// anything. This faces is called dangling face.
		//
		// The vertices of a dangling face has to be kept. Morover information
		// about adjacencies and interfaces of a dangling cell needs to be
		// updated (it it not necessary to update those information for the
		// cells that will be deleted, because they are going to be removed).
		int nCellInterfaces = cell.getInterfaceCount();
		const long *interfaces = cell.getInterfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];
			if (interfaceId < 0) {
				continue;
			}

			// Interfaces has to be considered just once
			if (deadInterfaces.count(interfaceId) > 0) {
				continue;
			}

			// Find if the owner or the neighbour of the interface will be
			// kept. An interface whose owner or neighbour will not be
			// deleted is called dangling.
			Interface &interface = m_interfaces[interfaceId];

			int danglingSide = -1;
			if (!interface.isBorder()) {
				if (cellIds.count(interface.getOwner()) == 0) {
					danglingSide = 0;
				} else if (cellIds.count(interface.getNeigh()) == 0) {
					danglingSide = 1;
				}
			}

			// Handle dangling faces
			if (danglingSide >= 0) {
				// Info on the dangling face
				long danglingCellId;
				long danglingNeighId;
				int danglingCellFace;
				if (danglingSide == 0) {
					danglingCellId   = interface.getOwner();
					danglingNeighId  = interface.getNeigh();
					danglingCellFace = interface.getOwnerFace();
				} else {
					danglingCellId   = interface.getNeigh();
					danglingNeighId  = interface.getOwner();
					danglingCellFace = interface.getNeighFace();
				}

				Cell &danglingCell = m_cells[danglingCellId];

				// Since the dangling cell will not be deleted, we have to
				// updated its interface and adjacency data structures.
				int cellInterfaceIndex = danglingCell.findInterface(danglingCellFace, interfaceId);
				danglingCell.deleteInterface(danglingCellFace, cellInterfaceIndex);

				int cellAdjacencyIndex = danglingCell.findAdjacency(danglingCellFace, danglingNeighId);
				danglingCell.deleteAdjacency(danglingCellFace, cellAdjacencyIndex);

				// Add the face to the dangling face list
				danglingFaces.insert({{danglingCellId, danglingCellFace}});
			}

			// Add the interface to the list of interfaces to delete
			deadInterfaces.insert(interfaceId);
		}

		// Delete cell
		deleteCell(cellId);
	}

	// Delete the interfaces
	for (long interfaceId : deadInterfaces) {
		VolumeKernel::deleteInterface(interfaceId, false, true);
	}

	// All the vertices belonging to the dangling cells has to be kept
	//
	// It's not enough to consider only the vertices on the dangling faces,
	// we have to consider the vertices of the whole cell. That's because
	// we may need to keep vertices on the edges of the cell, and those
	// vertices may not be on interfaces of the dangling face
	for (const FaceInfo &danglingFaceInfo : danglingFaces) {
		const Cell &cell = m_cells[danglingFaceInfo.id];

		// Vertices of the face
		const std::vector<int> &localConnect = cellLocalFaceConnect[danglingFaceInfo.face];
		for (int n = 0; n < nFaceVertices; ++n) {
			int localVertexId = localConnect[n];
			long vertexId = cell.getVertex(localVertexId);
			deadVertices.erase(vertexId);
		}

		// Vertices of all other interfaces of the cell
		int nCellInterfaces = cell.getInterfaceCount();
		const long *interfaces = cell.getInterfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];
			if (interfaceId < 0) {
				continue;
			}

			const Interface &interface = m_interfaces[interfaceId];
			for (int n = 0; n < nFaceVertices; ++n) {
				long vertexId = interface.getVertex(n);
				deadVertices.erase(vertexId);
			}
		}
	}

	// Delete the vertices
	for (long vertexId : deadVertices) {
		VolumeKernel::deleteVertex(vertexId, true);
	}

	// Done
	return danglingFaces;
}

/*!
	Creates a new patch vertex from the specified tree vertex.

	\param treeId is the id of the vertex in the tree
	\result The id of the newly created vertex.
*/
long VolOctree::addVertex(uint32_t treeId)
{
	// Vertex coordinates
	std::array<double, 3> nodeCoords = m_tree.getNodeCoordinates(treeId);

	// Create the vertex
	VertexIterator vertexIterator = VolumeKernel::addVertex(std::move(nodeCoords));

	// Done
	return vertexIterator->getId();
}

/*!
	Creates a new patch cell from the specified tree octant.

	\param octantInfo is the octant associated to the cell
	\param vertices are the vertices of the cell
	\param interfaces are the interfaces of the cell
	\param interfacesOwner is a flag that defines if an interface is owned
	by the cell
	\result The id of the newly created cell.
*/
long VolOctree::addCell(OctantInfo octantInfo, std::unique_ptr<long[]> &vertices)
{
	// Create the cell
	ElementInfo::Type cellType;
	if (isThreeDimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	CellIterator cellIterator = VolumeKernel::addCell(cellType, octantInfo.internal);
	Cell &cell = *cellIterator;
	long id = cell.getId();

	// Connectivity
	cell.setConnect(std::move(vertices));

	// If the cell is a ghost set its owner
#if BITPIT_ENABLE_MPI==1
	if (!octantInfo.internal) {
		uint64_t globalTreeId = m_tree.getGhostGlobalIdx(octantInfo.id);
		int rank = m_tree.getOwnerRank(globalTreeId);

		setGhostOwner(id, rank, false);
	}
#endif

	// Update cell to octant mapping
	if (octantInfo.internal) {
		m_cellToOctant.insert({{id, octantInfo.id}});
		m_octantToCell.insert({{octantInfo.id, id}});
	} else {
		m_cellToGhost.insert({{id, octantInfo.id}});
		m_ghostToCell.insert({{octantInfo.id, id}});
	}

	// Done
	return id;
}

/*!
	Deletes a cell from the patch.

	\param id is the id of the cell
*/
void VolOctree::deleteCell(long id)
{
	// Remove the information that link the cell to the octant
	bool interior = m_cells[id].isInterior();

	std::unordered_map<long, uint32_t, Element::IdHasher> *cellMap;
	if (interior) {
		cellMap = &m_cellToOctant;
	} else {
		cellMap = &m_cellToGhost;
	}

	std::unordered_map<long, uint32_t, Element::IdHasher>::const_iterator cellItr = cellMap->find(id);
	if (cellItr != cellMap->end()) {
		// Delete octant-to-cell entry
		std::unordered_map<uint32_t, long> *octantMap;
		if (interior) {
			octantMap = &m_octantToCell;
		} else {
			octantMap = &m_ghostToCell;
		}

		uint32_t treeId = cellItr->second;
		octantMap->erase(treeId);

		// Delete cell-to-octant entry
		cellMap->erase(cellItr);
	}

	// Delete the cell
	VolumeKernel::deleteCell(id, false, true);
}


/*!
	Build the adjacencies the cells.
*/
void VolOctree::updateAdjacencies(const std::vector<long> &cellIds, bool resetAdjacencies)
{
	// Face information
	int nCellFaces = 2 * getDimension();
	uint8_t oppositeFace[nCellFaces];
	m_tree.getOppface(oppositeFace);

	// Reset the adjacencies
	if (resetAdjacencies) {
		for (long cellId : cellIds) {
			m_cells[cellId].resetAdjacencies();
		}
	}

	// Sort the cells beased on their tree level
	int maxLevel = m_tree.getMaxDepth();
	size_t averageSize = cellIds.size() / maxLevel;
	std::vector<std::vector<long>> hierarchicalCellIds(maxLevel + 1);
	for (int level = 0; level <= maxLevel; ++level) {
		hierarchicalCellIds[level].reserve(averageSize);
	}

	for (long cellId : cellIds) {
		int cellLevel = getCellLevel(cellId);
		hierarchicalCellIds[cellLevel].push_back(cellId);
	}

	// Update the adjacencies
	FaceInfoSet processedFaces;
	processedFaces.reserve(cellIds.size() * getDimension());

	for (int level = 0; level <= maxLevel; ++level) {
		for (long cellId : hierarchicalCellIds[level]) {
			Cell &cell = m_cells[cellId];
			OctantInfo octantInfo = getCellOctant(cellId);
			for (int face = 0; face < nCellFaces; ++face) {
				FaceInfo currentFaceInfo(cellId, face);
				if (processedFaces.count(currentFaceInfo) > 0) {
					continue;
				}

				// Find cell neighbours
				std::vector<uint32_t> neighTreeIds;
				std::vector<bool> neighGhostFlags;
				if (octantInfo.internal) {
					m_tree.findNeighbours(octantInfo.id, face, 1, neighTreeIds, neighGhostFlags);
				} else {
					m_tree.findGhostNeighbours(octantInfo.id, face, 1, neighTreeIds, neighGhostFlags);
				}

				// Set the adjacencies
				//
				// Adjacencies will processed twice, once while processing the
				// current cell, and once while processing the neighbour cell.
				// However they will be set only once, because the function that
				// insert the adjacency in the cell will insert only unique
				// adjacencies.
				int nNeighs = neighTreeIds.size();
				for (int k = 0; k < nNeighs; ++k) {
					OctantInfo neighOctantInfo(neighTreeIds[k], !neighGhostFlags[k]);
					long neighId = getOctantId(neighOctantInfo);

					// Set cell data
					cell.pushAdjacency(face, neighId);

					// Set neighbour data
					int neighFace = oppositeFace[face];
					Cell &neigh = m_cells[neighId];
					neigh.pushAdjacency(neighFace, cellId);

					FaceInfo neighFaceInfo(neighId, neighFace);
					processedFaces.insert(neighFaceInfo);
				}
			}
		}
	}
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
bool VolOctree::_markCellForRefinement(const long &id)
{
	return set_marker(id, 1);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
bool VolOctree::_markCellForCoarsening(const long &id)
{
	return set_marker(id, -1);
}

/*!
	Set the marker on a cell.

	\param id is the id of the cell
	\param value is the value of the marker
*/
bool VolOctree::set_marker(const long &id, const int8_t &value)
{
	OctantInfo octantInfo = getCellOctant(id);
	if (!octantInfo.internal) {
		return false;
	}

	m_tree.setMarker(octantInfo.id, value);

	return true;
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
bool VolOctree::_enableCellBalancing(const long &id, bool enabled)
{
	OctantInfo octantInfo = getCellOctant(id);
	if (!octantInfo.internal) {
		return false;
	}

	m_tree.setBalance(octantInfo.id, enabled);

	return true;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] point is the point to be checked
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolOctree::isPointInside(const std::array<double, 3> &point)
{
	return (m_tree.getPointOwner(point) != nullptr);
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the idof the cell
	\param[in] point is the point to be checked
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolOctree::isPointInside(const long &id, const std::array<double, 3> &point)
{
	const Cell &cell = m_cells[id];

    int lowerLeftVertex  = 0;
	int upperRightVertex = pow(2, getDimension()) - 1;

	std::array<double, 3> lowerLeft  = getVertexCoords(cell.getVertex(lowerLeftVertex));
	std::array<double, 3> upperRight = getVertexCoords(cell.getVertex(upperRightVertex));

	const double EPS = getTol();
    for (int d = 0; d < 3; ++d){
        if (point[d] < lowerLeft[d] - EPS || point[d] > upperRight[d] + EPS) {
            return false;
        }
    }

    return true;
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] point is the point to be checked
	\result Returns the id of the cell the contains the point. If the point
	is not inside the patch, the function returns the id of the null element.
*/
long VolOctree::locatePoint(const std::array<double, 3> &point)
{
	Octant *octant = m_tree.getPointOwner(point);
	if (m_tree.getPointOwner(point) == nullptr) {
		return Element::NULL_ID;
	}

	OctantInfo octantInfo(m_tree.getIdx(octant), true);
	return getOctantId(octantInfo);
}

/*!
	Internal function to set the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void VolOctree::_setTol(double tolerance)
{
	m_tree.setTol(tolerance);

	VolumeKernel::_setTol(tolerance);
}

/*!
	Internal function to reset the tolerance for the geometrical checks.
*/
void VolOctree::_resetTol()
{
	m_tree.setTol();

	double tolerance = m_tree.getTol();
	VolumeKernel::_setTol(tolerance);
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
 */
void VolOctree::translate(std::array<double, 3> translation)
{
	m_tree.setOrigin(m_tree.getOrigin() + translation);

	VolumeKernel::translate(translation);
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor vector
 */
void VolOctree::scale(std::array<double, 3> scaling)
{
	bool uniformScaling = true;
	uniformScaling &= (fabs(scaling[0] - scaling[1]) > 1e-14);
	uniformScaling &= (fabs(scaling[0] - scaling[2]) > 1e-14);
	assert(uniformScaling);
	if (!uniformScaling) {
		log::cout() << "octree patch only allows uniform scaling)" << std::endl;
		return;
	}

	m_tree.setL(m_tree.getL() * scaling[0]);

	initializeTreeGeometry();

	VolumeKernel::scale(scaling);
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> VolOctree::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return neighs;
	}

	// Get edge neighbours
	int codimension = getDimension() - 1;
	neighs = findCellCodimensionNeighs(id, edge, codimension, blackList);

	// Add face neighbours
	for (int face : m_octantLocalFacesOnEdge[edge]) {
		for (auto &neigh : _findCellFaceNeighs(id, face, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a local vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> VolOctree::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;

	// Get vertex neighbours
	int codimension = getDimension();
	neighs = findCellCodimensionNeighs(id, vertex, codimension, blackList);

	// Add edge neighbours or face neighbours
	if (isThreeDimensional()) {
		for (int edge : m_octantLocalEdgesOnVertex[vertex]) {
			for (auto &neigh : _findCellEdgeNeighs(id, edge, blackList)) {
				utils::addToOrderedVector<long>(neigh, neighs);
			}
		}
	} else {
		for (int face : m_octantLocalFacesOnVertex[vertex]) {
			for (auto &neigh : _findCellFaceNeighs(id, face, blackList)) {
				utils::addToOrderedVector<long>(neigh, neighs);
			}
		}
	}

	return neighs;
}

/*!
	Finds the neighbours for the given co-dimension of the specified cell.

	Only the neighbours for the specified co-dimension are found, neighbours
	of higher co-dimensions are not inserted in the returned list.

	\param id is the id of the cell
	\param codimension is the co-dimension
	\param index is the local index of the entity (vertex, edge or face)
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours for the given codimension of the specified cell.
*/
std::vector<long> VolOctree::findCellCodimensionNeighs(const long &id, const int &index,
											           const int &codimension, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	int dimension = getDimension();
	if (codimension > dimension || codimension <= 0) {
		return neighs;
	}

	OctantInfo octantInfo = getCellOctant(id);

	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	if (octantInfo.internal) {
		m_tree.findNeighbours(octantInfo.id, index, codimension, neighTreeIds, neighGhostFlags);
	} else {
		m_tree.findGhostNeighbours(octantInfo.id, index, codimension, neighTreeIds);
		neighGhostFlags.resize(neighTreeIds.size(), false);
	}

	int nNeighs = neighTreeIds.size();
	for (int i = 0; i < nNeighs; ++i) {
		OctantInfo neighOctantInfo(neighTreeIds[i], !neighGhostFlags[i]);
		long neighId = getOctantId(neighOctantInfo);

		if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	@}
*/

}
