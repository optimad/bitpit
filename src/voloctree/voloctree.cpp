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
	  m_tree(origin[0], origin[1], origin[2], length, dimension),
	  m_lastTreeOperation(OP_INITIALIZATION)
{
	log::cout() << ">> Initializing Octree mesh\n";

	// Inizializzazione dell'octree
	double initial_level = ceil(log2(std::max(1., length / dh)));

	m_tree.setMarker((uint32_t) 0, initial_level);

	// Info sull'octree
	initializeTreeGeometry();

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
const std::vector<Adaption::Info> VolOctree::_updateAdaption(bool trackAdaption)
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
		return std::vector<Adaption::Info>();
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
const std::vector<Adaption::Info> VolOctree::sync(bool trackChanges)
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

	// Initialize intersections
	log::cout() << ">> Evaluating Octree intersections...";

	m_tree.computeIntersections();

	log::cout() << " Done" << std::endl;

	// Initialize tracking data
	std::vector<Adaption::Info> adaptionData;

	// Extract information for transforming the patch
	//
	// If there are no cells in the mesh we need to import all
	// octants.
	log::cout() << ">> Extract information for transforming the patch...";

	std::vector<OctantInfo> newOctants;
	std::unordered_map<uint32_t, long> renumberedOctants;
	std::vector<long> removedCells;

	newOctants.reserve(nOctants + nGhostsOctants);
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
		Adaption::Type adaptionType = Adaption::TYPE_NONE;
		if (importAll) {
			adaptionType = Adaption::TYPE_CREATION;
		} else if (lastTreeOperation == OP_ADAPTION_MAPPED) {
			bool isNewR = m_tree.getIsNewR(treeId);
			if (isNewR) {
				adaptionType = Adaption::TYPE_REFINEMENT;
			} else {
				bool isNewC = m_tree.getIsNewC(treeId);
				if (isNewC) {
					adaptionType = Adaption::TYPE_COARSENING;
				} else if (treeId != mapper_octantMap.front()) {
					adaptionType = Adaption::TYPE_RENUMBERING;
				}
			}
#if BITPIT_ENABLE_MPI==1
		} else if (lastTreeOperation == OP_LOAD_BALANCE) {
			if (getRank() != mapper_octantRank.front()) {
				adaptionType = Adaption::TYPE_PARTITION_RECV;
			} else if (treeId != mapper_octantMap.front()) {
				adaptionType = Adaption::TYPE_RENUMBERING;
			}
#endif
		}

		// If the octant cell has not been modified we can skip to the next
		// octant.
		if (adaptionType == Adaption::TYPE_NONE) {
			++treeId;
			continue;
		}

		// Re-numbered cells just need to be added to the proper list.
		//
		// Renumbered cells are not tracked, because the re-numbering
		// only happens inside VolOctree.
		if (adaptionType == Adaption::TYPE_RENUMBERING) {
			OctantInfo previousOctantInfo(mapper_octantMap.front(), !mapper_ghostFlag.front());
			long cellId = getOctantId(previousOctantInfo);

			renumberedOctants.insert({{treeId, cellId}});

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
		} else if (adaptionType == Adaption::TYPE_REFINEMENT) {
			nCurrentTreeIds = pow(2, getDimension());
		} else {
			nCurrentTreeIds = 1;
		}

		const long lastCurrentTreeId = treeId + nCurrentTreeIds;
		for (int currentTreeId = treeId; currentTreeId < lastCurrentTreeId; ++currentTreeId) {
			newOctants.emplace_back(currentTreeId, true);
		}

		// Cells that will be removed
		//
		// For refinement and coarsening the cells associated to the previous
		// octant ids needs to be removed.
		if (adaptionType == Adaption::TYPE_REFINEMENT || adaptionType == Adaption::TYPE_COARSENING) {
			int nPreviousTreeIds = mapper_octantMap.size();
			for (int k = 0; k < nPreviousTreeIds; ++k) {
				OctantInfo previousOctantInfo(mapper_octantMap[k], !mapper_ghostFlag[k]);

				removedCells.emplace_back();
				long &cellId = removedCells.back();
				cellId = getOctantId(previousOctantInfo);
			}
		}

		// Adaption tracking
		if (trackChanges) {
			adaptionData.emplace_back();
			Adaption::Info &adaptionInfo = adaptionData.back();

			// Type
			adaptionInfo.type = adaptionType;

			// Element
			adaptionInfo.entity = Adaption::ENTITY_CELL;

			// Current status
			//
			// We don't know the id of the current status, because those
			// cells are not yet in the mesh. Store the trre id, and
			// make the translation later.
			//
			// WARNING: tree id are uint32_t wherase adaptionInfo stores
			//          id as unsigned long.
			adaptionInfo.current.reserve(nCurrentTreeIds);
			auto newOctantsIter = newOctants.cend() - nCurrentTreeIds;
			while (newOctantsIter != newOctants.cend()) {
				adaptionInfo.current.emplace_back();
				unsigned long &adaptionId = adaptionInfo.current.back();
				adaptionId = (*newOctantsIter).id;

				newOctantsIter++;
			}

			// Previous cell and interface ids
			if (adaptionType != Adaption::TYPE_PARTITION_RECV) {
				int nPreviousCellIds = mapper_octantMap.size();
				adaptionInfo.previous.reserve(nPreviousCellIds);
				auto removedCellsIter = removedCells.cend() - nPreviousCellIds;
				while (removedCellsIter != removedCells.cend()) {
					const long &id = *removedCellsIter;

					adaptionInfo.previous.emplace_back();
					unsigned long &adaptionId = adaptionInfo.previous.back();
					adaptionId = id;

					removedCellsIter++;
				}
			}

			// Rank of received cells
			if (adaptionType == Adaption::TYPE_PARTITION_RECV) {
				adaptionInfo.rank = mapper_octantRank[0];
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
		if (trackChanges) {
			adaptionData.emplace_back();
			Adaption::Info &adaptionInfo = adaptionData.back();
			adaptionInfo.entity = Adaption::ENTITY_CELL;
			if (rank >= 0) {
				adaptionInfo.type = Adaption::TYPE_PARTITION_SEND;
				adaptionInfo.rank = rank;
			} else {
				adaptionInfo.type = Adaption::TYPE_DELETION;
			}
		}

		// Add the send cells to the list of cells to be removed
		for (int k = 0; k < 2; ++k) {
			uint32_t beginTreeId = rankEntry.second[2 * k];
			uint32_t endTreeId   = rankEntry.second[2 * k + 1];
			for (uint32_t treeId = beginTreeId; treeId < endTreeId; ++treeId) {
				OctantInfo octantInfo(treeId, true);

				removedCells.emplace_back();
				long &cellId = removedCells.back();
				cellId = getOctantId(octantInfo);

				if (trackChanges) {
					Adaption::Info &adaptionInfo = adaptionData.back();
					adaptionInfo.previous.emplace_back();
					unsigned long &sentId = adaptionInfo.previous.back();
					sentId = cellId;
				}
			}
		}
	}

	// Previous ghosts cells need to be removed
	if (nPreviousGhosts > 0) {
		long deletedGhostsInfoIdx = -1;
		if (trackChanges) {
			adaptionData.emplace_back();
			Adaption::Info &deletedGhostsInfo = adaptionData.back();
			deletedGhostsInfo.type   = Adaption::TYPE_DELETION;
			deletedGhostsInfo.entity = Adaption::ENTITY_CELL;
			deletedGhostsInfo.previous.reserve(nPreviousGhosts);

			deletedGhostsInfoIdx = adaptionData.size() - 1;
		}

		auto cellIterator = m_cellToGhost.cbegin();
		while (cellIterator != m_cellToGhost.cend()) {
			long id = cellIterator->first;

			removedCells.emplace_back();
			long &removedId = removedCells.back();
			removedId = id;

			// Adaption tracking
			if (trackChanges) {
				Adaption::Info &deletedGhostsInfo = adaptionData[deletedGhostsInfoIdx];
				deletedGhostsInfo.previous.emplace_back();
				unsigned long &adaptionId = deletedGhostsInfo.previous.back();
				adaptionId = id;
			}

			// Increment iterator
			cellIterator++;
		}
	}
	removedCells.shrink_to_fit();

	// New ghost octants need to be added
	std::unordered_map<int, std::vector<uint32_t>> ghostTreeIds;
	for (uint32_t treeId = 0; treeId < (uint32_t) nGhostsOctants; ++treeId) {
		newOctants.emplace_back(treeId, false);

		uint64_t globalTreeId = m_tree.getGhostGlobalIdx(treeId);
		int rank = m_tree.getOwnerRank(globalTreeId);
		ghostTreeIds[rank].push_back(treeId);
	}
	newOctants.shrink_to_fit();
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

			// Adaption info
			adaptionData.emplace_back();
			Adaption::Info &deletedInterfacesInfo = adaptionData.back();
			deletedInterfacesInfo.type   = Adaption::TYPE_DELETION;
			deletedInterfacesInfo.entity = Adaption::ENTITY_INTERFACE;
			for (const long &interfaceId : removedInterfaces) {
				deletedInterfacesInfo.previous.emplace_back();
				unsigned long &deletedInterfaceId = deletedInterfacesInfo.previous.back();
				deletedInterfaceId = interfaceId;
			}
		}

		// Delete cells
		danglingFaces = removeCells(removedCells);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Cells removed: " <<  removedCells.size() << std::endl;
	}

	std::vector<long>().swap(removedCells);

	// Reserve space in the octant-to-cell maps
	m_cellToOctant.reserve(nOctants);
	m_octantToCell.reserve(nOctants);

	// Remap renumbered cells
	if (renumberedOctants.size() > 0) {
		log::cout() << ">> Rebuilding octant-to-cell map for renumbered cells...";

		// Remap cells to the new tree ids
		auto cellIterator = renumberedOctants.begin();
		while (cellIterator != renumberedOctants.end()) {
			long cellId = cellIterator->second;

			uint32_t currentTreeId  = cellIterator->first;
			uint32_t previousTreeId = m_cellToOctant.at(cellId);

			m_cellToOctant[cellId] = currentTreeId;
			m_octantToCell[currentTreeId] = cellId;
			if (renumberedOctants.count(previousTreeId) == 0) {
				m_octantToCell.erase(previousTreeId);
			}

			cellIterator++;
		}

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Cells renumbered: " <<  renumberedOctants.size() << std::endl;
	}

	std::unordered_map<uint32_t, long>().swap(renumberedOctants);

	// Reset ghost maps
	m_cellToGhost.clear();
	m_cellToGhost.reserve(nGhostsOctants);

	m_ghostToCell.clear();
	m_ghostToCell.reserve(nGhostsOctants);

	// Import added octants
	std::vector<unsigned long> createdInterfaces;
	if (newOctants.size() > 0) {
		log::cout() << ">> Importing new octants...";

		createdInterfaces = importOctants(newOctants, danglingFaces);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Octants imported: " <<  newOctants.size() << std::endl;
	}

	FaceInfoSet().swap(danglingFaces);

	// Rebuild the ghost information
#if BITPIT_ENABLE_MPI==1
	rebuildGhostExchangeData(ghostTreeIds);
#endif

	// Disable advanced editing
	setExpert(false);

	// Track mesh adaption
	if (trackChanges) {
		// Map ids of the added cells
		for (auto &adaptionInfo : adaptionData) {
			if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
				continue;
			}

			int nCurrentIds = adaptionInfo.current.size();
			for (int k = 0; k < nCurrentIds; ++k) {
				long cellId = m_octantToCell.at(adaptionInfo.current[k]);
				adaptionInfo.current[k] = cellId;
			}
		}

		// Track created ghosts cells
		if (nGhostsOctants > 0) {
			adaptionData.emplace_back();
			Adaption::Info &createdGhostsInfo = adaptionData.back();
			createdGhostsInfo.type   = Adaption::TYPE_CREATION;
			createdGhostsInfo.entity = Adaption::ENTITY_CELL;

			createdGhostsInfo.previous.reserve(nGhostsOctants);
			auto cellIterator = m_cellToGhost.cbegin();
			while (cellIterator != m_cellToGhost.cend()) {
				createdGhostsInfo.previous.emplace_back();
				unsigned long &adaptionId = createdGhostsInfo.previous.back();
				adaptionId = cellIterator->first;

				cellIterator++;
			}
		}

		// Track created interfaces
		if (createdInterfaces.size() > 0) {
			adaptionData.emplace_back();
			Adaption::Info &createdInterfacesInfo = adaptionData.back();
			createdInterfacesInfo.type   = Adaption::TYPE_CREATION;
			createdInterfacesInfo.entity = Adaption::ENTITY_INTERFACE;

			createdInterfacesInfo.current.swap(createdInterfaces);
		}
	} else {
		adaptionData.emplace_back();
	}

	// Delete tree conenctivity
	m_tree.clearConnectivity();

	// Done
	return adaptionData;
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
std::vector<unsigned long> VolOctree::importOctants(std::vector<OctantInfo> &octantInfoList)
{
	FaceInfoSet danglingFaces;

	return importOctants(octantInfoList, danglingFaces);
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
	\param danglingFaces is the list of dangling faces in the current mesh
*/
std::vector<unsigned long> VolOctree::importOctants(std::vector<OctantInfo> &octantInfoList,
                                 FaceInfoSet &danglingFaces)
{
	// Info of the cells
	ElementInfo::Type cellType;
	if (isThreeDimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
	const int &nCellFaces = cellTypeInfo.nFaces;
	const int &nCellVertices = cellTypeInfo.nVertices;
	const std::vector<std::vector<int>> &cellLocalFaceConnect = cellTypeInfo.faceConnect;

	// Info on the interfaces
	ElementInfo::Type interfaceType;
	if (isThreeDimensional()) {
		interfaceType = ElementInfo::PIXEL;
	} else {
		interfaceType = ElementInfo::LINE;
	}

	const ElementInfo &interfaceTypeInfo = ElementInfo::getElementInfo(interfaceType);
	const int &nInterfaceVertices = interfaceTypeInfo.nVertices;

	uint32_t nIntersections = m_tree.getNumIntersections();

	std::vector<unsigned long> createdInterfaces;
	createdInterfaces.reserve(nCellFaces * octantInfoList.size());

	// Add the vertex of the dangling faces to the vertex map
	std::unordered_map<uint32_t, long> vertexMap;
	for (auto &danglingFaceInfo : danglingFaces) {
		// List of faces with the vertx to be added
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
				long vertexId = cellConnect[localConnect[k]];
				uint32_t vertexTreeId = octantTreeConnect[localConnect[k]];

				vertexMap.insert({{vertexTreeId, vertexId}});
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

	// Create the interfaces
	std::unordered_map<uint32_t, long> interfaceMap;
	std::unordered_map<uint32_t, std::vector<uint32_t>> octantInterfaces;
	std::unordered_map<uint32_t, std::vector<uint32_t>> octantTreeInterfaces;
	for (uint32_t interfaceTreeId = 0; interfaceTreeId < nIntersections; ++interfaceTreeId) {
		// Skip the interface is already inserted in the patch
		if (interfaceMap.count(interfaceTreeId) != 0) {
			continue;
		}

		// Info on the interface
		Intersection *treeInterface = m_tree.getIntersection(interfaceTreeId);

		int owner       = m_tree.getOut(treeInterface);
		int ownerFace   = m_tree.getFace(treeInterface);
		int neigh       = m_tree.getIn(treeInterface);
		bool isBoundary = m_tree.getBound(treeInterface);
		bool isGhost    = m_tree.getIsGhost(treeInterface);

		bool ownerIsGhost;
		bool neighIsGhost;
		if (isGhost) {
			ownerIsGhost = m_tree.getOutIsGhost(treeInterface);
			neighIsGhost = !ownerIsGhost;
		} else {
			ownerIsGhost = false;
			neighIsGhost = false;
		}

		// Decide if we need to build the interface
		bool buildInterface = false;

		OctantInfo ownerOctantInfo(owner, !ownerIsGhost);
		long ownerId = getOctantId(ownerOctantInfo);
		if (ownerId < 0) {
			octantTreeInterfaces[ownerOctantInfo.id].push_back(interfaceTreeId);
			buildInterface = true;
		}

		long neighId = Element::NULL_ID;
		if (!isBoundary) {
			OctantInfo neighOctantInfo(neigh, !neighIsGhost);
			neighId = getOctantId(neighOctantInfo);
			if (neighId < 0) {
				octantTreeInterfaces[neighOctantInfo.id].push_back(interfaceTreeId);
				buildInterface = true;
			}
		}

		if (!buildInterface) {
			continue;
		}

		// Interface owner and neighbour faces
		//
		// The owner or the neighbour may be not know yet, but the
		// corresponding faces are. So, set the faces also for
		// the unknown cells because they are needed when finding
		// the interfaces associated to a cell.
		std::array<FaceInfo, 2> interfaceFaces;
		interfaceFaces[0] = FaceInfo(ownerId, ownerFace);
		interfaceFaces[1] = FaceInfo(neighId, ownerFace + 1 - 2 * (ownerFace % 2));

		// Interface connectivity
		const std::vector<uint32_t> &octantTreeConnect = getOctantConnect(ownerOctantInfo);
		const std::vector<int> &localConnect = cellLocalFaceConnect[ownerFace];
		std::unique_ptr<long[]> interfaceConnect = std::unique_ptr<long[]>(new long[nInterfaceVertices]);
		for (int k = 0; k < nInterfaceVertices; ++k) {
			interfaceConnect[k] = vertexMap.at(octantTreeConnect[localConnect[k]]);
		}

		// Create the interface
		createdInterfaces.emplace_back();
		unsigned long &interfaceId = createdInterfaces.back();
		interfaceId = addInterface(interfaceTreeId, interfaceConnect, interfaceFaces);
		interfaceMap[interfaceTreeId] = interfaceId;

		// If the interface is on an dangling faces, the owner or
		// the neigbour of the interface already exists. We need
		// to identify that cell and set that information on the
		// interface.
		//
		// A boundary interface cannnot be on a dangling face.
		if (!isBoundary) {
			FaceInfo *guessDanglingInfo;
			if (ownerId < 0) {
				guessDanglingInfo = &interfaceFaces[1];
			} else {
				guessDanglingInfo = &interfaceFaces[0];
			}

			if (danglingFaces.count(*guessDanglingInfo) != 0) {
				Cell &danglingCell = m_cells[guessDanglingInfo->id];
				danglingCell.pushInterface(guessDanglingInfo->face, interfaceId);
			}
		}
	}

	// Add the cells
	std::vector<std::vector<long>> cellAdjacencies(nCellFaces, std::vector<long>());
	std::vector<std::vector<long>> cellInterfaces(nCellFaces, std::vector<long>());
	std::vector<std::vector<bool>> cellInterfacesOwner(nCellFaces, std::vector<bool>());
	for (OctantInfo &octantInfo : octantInfoList) {
		// Octant connectivity
		const std::vector<uint32_t> &octantTreeConnect = getOctantConnect(octantInfo);

		// Cell connectivity
		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int k = 0; k < nCellVertices; ++k) {
			uint32_t vertexTreeId = octantTreeConnect[k];
			cellConnect[k] = vertexMap.at(vertexTreeId);
		}

		// Cell interfaces and adjacencies
		for (int k = 0; k < nCellFaces; ++k) {
			cellAdjacencies[k].clear();
			cellInterfaces[k].clear();
			cellInterfacesOwner[k].clear();
		}

		for (uint32_t interfaceTreeId : octantTreeInterfaces[octantInfo.id]) {
			const long &interfaceId = interfaceMap.at(interfaceTreeId);
			const Interface &interface = m_interfaces[interfaceId];

			Intersection *treeInterface = m_tree.getIntersection(interfaceTreeId);
			uint32_t owner = m_tree.getOut(treeInterface);
			bool ownsInterface = (owner == octantInfo.id);

			int cellFace;
			if (ownsInterface) {
				cellFace = interface.getOwnerFace();
			} else {
				cellFace = interface.getNeighFace();
			}

			cellAdjacencies[cellFace].emplace_back();
			long &cellAdjacencyId = cellAdjacencies[cellFace].back();
			if (ownsInterface) {
				cellAdjacencyId = interface.getNeigh();
			} else {
				cellAdjacencyId = interface.getOwner();
			}

			cellInterfaces[cellFace].emplace_back();
			long &cellInterfaceId = cellInterfaces[cellFace].back();
			cellInterfaceId = interfaceId;

			cellInterfacesOwner[cellFace].push_back(ownsInterface);
		}

		// Add cell
		addCell(octantInfo, cellConnect, cellAdjacencies,
			   cellInterfaces, cellInterfacesOwner);
	}

	// Done
	createdInterfaces.shrink_to_fit();

	return createdInterfaces;
}

/*!
	Remove a list of octants from the patch.

	\param cellIds is the list of cells ids to remove
*/
VolOctree::FaceInfoSet VolOctree::removeCells(std::vector<long> &cellIds)
{
	// Delete cells
	//
	// To make deletion of interfaces more efficient, keep the list of
	// faces to remove in an ordered map. The order should be based upon
	// the position of the interfaces in the container, rather than upon
	// the id of the faces.
	typedef bitpit::PiercedVector<Interface>::positionLess InterfacePositionCompare;
	typedef std::map<long, int, InterfacePositionCompare> InterfaceOrderedMap;

	std::unordered_set<long> deadVertices;
	InterfaceOrderedMap deadInterfaces((InterfacePositionCompare(m_interfaces)));
	for (long cellId : cellIds) {
		Cell &cell = m_cells[cellId];

		// List dead vertices
		//
		// For now, all cell vertices will be listed. Later, the vertex of
		// the dangling faces will be removed from the list.
		int nCellVertices = cell.getVertexCount();
		const long *vertices = cell.getConnect();
		for (int k = 0; k < nCellVertices; ++k) {
			deadVertices.insert(vertices[k]);
		}

		// List dead interface and set the dangling status. An interface
		// is dangling if only the owner or the neighbour will be deleted.
		int nCellInterfaces = cell.getInterfaceCount();
		const long *interfaces = cell.getInterfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];
			if (interfaceId < 0) {
				continue;
			}

			int danglingSide = -1;
			if (deadInterfaces.count(interfaceId) == 0) {
				Interface &interface = m_interfaces[interfaceId];
				if (!interface.isBorder()) {
					if (interface.getOwner() == cellId) {
						danglingSide = 1;
					} else {
						danglingSide = 0;
					}
				}
			}
			deadInterfaces[interfaceId] = danglingSide;
		}

		// Delete cell
		deleteCell(cellId);
	}

	// Delete interfaces
	FaceInfoSet danglingFaces;

	for (auto it = deadInterfaces.begin(); it != deadInterfaces.end(); ++it) {
		long interfaceId  = it->first;
		long danglingSide = it->second;

		// Handle dangling interfaces
		if (danglingSide >= 0) {
			Interface &interface = m_interfaces[interfaceId];

			long danglingCellId;
			long danglingNeighId;
			long danglingCellFace;
			if (danglingSide == 0) {
				danglingCellId   = interface.getOwner();
				danglingNeighId  = interface.getNeigh();
				danglingCellFace = interface.getOwnerFace();
			} else {
				danglingCellId   = interface.getNeigh();
				danglingNeighId  = interface.getOwner();
				danglingCellFace = interface.getNeighFace();
			}

			// Remove interface vertices from dead vertices
			int nFaceVertices = interface.getVertexCount();
			const long *vertices = interface.getConnect();
			for (int k = 0; k < nFaceVertices; ++k) {
				deadVertices.erase(vertices[k]);
			}

			// Remove interface and adjacency from dangling cell
			Cell &danglingCell = m_cells[danglingCellId];

			int cellInterfaceIndex = danglingCell.findInterface(danglingCellFace, interfaceId);
			danglingCell.deleteInterface(danglingCellFace, cellInterfaceIndex);

			int cellAdjacencyIndex = danglingCell.findAdjacency(danglingCellFace, danglingNeighId);
			danglingCell.deleteAdjacency(danglingCellFace, cellAdjacencyIndex);

			// Add the associated cell face to the dangling faces list
			FaceInfo danglingFace(danglingCellId, danglingCellFace);
			danglingFaces.insert(danglingFace);
		}

		// Add the interface to the list of interfaces to delete
		VolumeKernel::deleteInterface(interfaceId, false, true);
	}

	// Delete vertices
	for (auto it = deadVertices.begin(); it != deadVertices.end(); ++it) {
		VolumeKernel::deleteVertex(*it, true);
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
	Creates a new patch interface from the specified tree intersection.

	\param treeId is the id of the intersection in the tree
	\param vertices are the vertices of the interface
	\param faces are the faces of the interface
	\result The id of the newly created interface.
*/
long VolOctree::addInterface(uint32_t treeId,
                                   std::unique_ptr<long[]> &vertices,
                                   std::array<FaceInfo, 2> &faces)
{
	BITPIT_UNUSED(treeId);

	// Info on the interfaces
	ElementInfo::Type interfaceType;
	if (isThreeDimensional()) {
		interfaceType = ElementInfo::PIXEL;
	} else {
		interfaceType = ElementInfo::LINE;
	}

	// Create the interface
	InterfaceIterator interfaceIterator = VolumeKernel::addInterface(interfaceType);
	Interface &interface = *interfaceIterator;

	// Connectivity
	interface.setConnect(std::move(vertices));

	// Owner and neighbour
	interface.setOwner(faces[0].id, faces[0].face);
	interface.setNeigh(faces[1].id, faces[1].face);

	// Done
	return interface.getId();
}

/*!
	Creates a new patch cell from the specified tree octant.

	\param octantInfo is the octant associated to the cell
	\param vertices are the vertices of the cell
	\param adjacencies are the adjacencies of the cell
	\param interfaces are the interfaces of the cell
	\param interfacesOwner is a flag that defines if an interface is owned
	by the cell
	\result The id of the newly created cell.
*/
long VolOctree::addCell(OctantInfo octantInfo,
                              std::unique_ptr<long[]> &vertices,
                              std::vector<std::vector<long>> &adjacencies,
                              std::vector<std::vector<long>> &interfaces,
                              std::vector<std::vector<bool>> &interfacesOwner)
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

	// Interfaces
	cell.setInterfaces(interfaces);

	// Adjacencies
	cell.setAdjacencies(adjacencies);

	// Update data of interfaces and neighbours
	int nCellFaces = cell.getFaceCount();
	for (int face = 0; face < nCellFaces; ++face) {
		int nNeighbours = adjacencies[face].size();
		for (int k = 0; k < nNeighbours; ++k) {
			long interfaceId = interfaces[face][k];
			Interface &interface = m_interfaces[interfaceId];
			bool ownsInterface = interfacesOwner[face][k];

			// Update data of interfaces
			if (ownsInterface) {
				interface.setOwner(id, face);
			} else {
				interface.setNeigh(id, face);
			}

			// Update data of neighbours
			long neighId = adjacencies[face][k];
			if (neighId > 0) {
				int neighFace;
				if (ownsInterface) {
					neighFace = interface.getNeighFace();
				} else {
					neighFace = interface.getOwnerFace();
				}

				Cell &neigh = m_cells.at(neighId);
				neigh.pushAdjacency(neighFace, id);
			}
		}
	}

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
	@}
*/

}
