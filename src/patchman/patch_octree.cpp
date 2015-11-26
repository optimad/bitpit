//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#ifndef DISABLE_OCTREE
#include "patch_octree.hpp"

#include <math.h>

#include "utils.hpp"

namespace pman {

/*!
	\class PatchOctree

	\brief The PatchOctree defines a Octree patch.

	PatchOctree defines a Octree patch.
*/

/*!
	Creates a new patch.
*/
PatchOctree::PatchOctree(const int &id, const int &dimension,
				 std::array<double, 3> origin, double length, double dh )
	: Patch(id, dimension)
{
	std::cout << ">> Initializing Octree mesh\n";

	// Info sull'octree
	int maxLevels;
	if (is_three_dimensional()) {
		maxLevels = 32;
	} else {
		maxLevels = 32;
	}

	m_tree_dh.reserve(maxLevels);
	m_tree_area.reserve(maxLevels);
	m_tree_volume.reserve(maxLevels);
	for(int i = 0; i < maxLevels; i++) {
	    double levelLength = length / ((double) pow(2,i));

	    m_tree_dh.push_back(pow(levelLength, 1.));
	    m_tree_area.push_back(pow(levelLength, (double) (get_dimension() - 1)));
	    m_tree_volume.push_back(pow(levelLength, (double) (get_dimension())));
	};

	// Inizializzazione dell'octree
	double initial_level = ceil(log2(max(1., length / dh)));
	if (is_three_dimensional()) {
		m_tree_3D = Class_Para_Tree<3>(origin[0], origin[1], origin[2], length);
		m_tree_3D.setMarker((uint32_t) 0, initial_level);
	} else {
		m_tree_2D = Class_Para_Tree<2>(origin[0], origin[1], origin[2], length);
		m_tree_2D.setMarker((uint32_t) 0, initial_level);
	}

	// Info sulle interfacce
	for (int i = 0; i < dimension; i++) {
		for (int n = -1; n <= 1; n += 2) {
			std::array<double, 3> normal = {0.0, 0.0, 0.0};
			normal[i] = n;

			m_normals.push_back(normal);
		}
	}
}

/*!
	Destroys the patch.
*/
PatchOctree::~PatchOctree()
{

}

/*!
	Gets a pointer to the the opposite normal.

	\param normal is a pointer to the normal
	\result A pointer to the opposite normal.
 */
std::array<double, 3> & PatchOctree::_get_opposite_normal(std::array<double, 3> &normal)
{
	std::vector<std::array<double, 3> >::iterator itr_current = std::find(m_normals.begin(), m_normals.end(), normal);
	if (itr_current == m_normals.end()) {
		 throw std::out_of_range ("Input normal is not among the stored normals");
	}

	int dimension = get_dimension();

	int id_current  = std::distance(m_normals.begin(), itr_current);
	int id_opposite = (id_current + dimension) % (2 * dimension);

	return m_normals[id_opposite];
}

/*!
	Gets the octant of the cell with the specified id.

	\param id the id of the cell
	\result The octant info of the specified cell
*/
PatchOctree::OctantInfo PatchOctree::get_cell_octant(const long &id) const
{
	bool internal = (m_cells[id].get_position_type() == Cell::INTERNAL);

	OctantInfo octantInfo;
	octantInfo.internal = internal;
	if (internal) {
		octantInfo.id = m_cell_to_octant.at(id);
	} else {
		octantInfo.id = m_cell_to_ghost.at(id);
	}

	return octantInfo;
}

/*!
	Gets the id of the specified octant.

	\param octantInfo the data of the octant
	\result The id of the specified octant
*/
long PatchOctree::get_octant_id(const OctantInfo &octantInfo) const
{
	std::unordered_map<uint32_t, long>::const_iterator octantItr;
	if (octantInfo.internal) {
		octantItr = m_octant_to_cell.find(octantInfo.id);
		if (octantItr == m_octant_to_cell.end()) {
			return Element::NULL_ELEMENT_ID;
		}
	} else {
		octantItr = m_ghost_to_cell.find(octantInfo.id);
		if (octantItr == m_ghost_to_cell.end()) {
			return Element::NULL_ELEMENT_ID;
		}
	}

	return octantItr->second;
}

/*!
	Gets the connectivity of the specified octant.

	\param octantInfo the data of the octant
	\result A reference to the octant's connectivity
*/
const std::vector<uint32_t> & PatchOctree::get_octant_connect(const OctantInfo &octantInfo)
{
	if (is_three_dimensional()) {
		bool isGhost = m_tree_3D.getIsGhost(octantInfo.id);
		if (!isGhost) {
			return m_tree_3D.getConnectivity()[octantInfo.id];
		} else {
			return m_tree_3D.getGhostConnectivity()[octantInfo.id];
		}
	} else {
		bool isGhost = m_tree_2D.getIsGhost(octantInfo.id);
		if (!isGhost) {
			return m_tree_2D.getConnectivity()[octantInfo.id];
		} else {
			return m_tree_2D.getGhostConnectivity()[octantInfo.id];
		}
	}
}

/*!
	Evaluates a unique hash for the octant.

	\param octantInfo the data of the octant
	\result A unique hash for the octant.
*/
PatchOctree::OctantHash PatchOctree::evaluate_octant_hash(const OctantInfo &octantInfo)
{
	uint8_t level;
	uint64_t morton;
	if (is_three_dimensional()) {
		Class_Octant<3> *octant = m_tree_3D.getOctant(octantInfo.id);
		level  = octant->getLevel();
		morton = octant->computeMorton();
	} else {
		Class_Octant<2> *octant = m_tree_2D.getOctant(octantInfo.id);
		level  = octant->getLevel();
		morton = octant->computeMorton();
	}

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
int PatchOctree::get_cell_level(const long &id)
{
	OctantInfo octantInfo = get_cell_octant(id);
	if (is_three_dimensional()) {
		Class_Octant<3> *octant;
		if (m_tree_3D.getIsGhost(octantInfo.id)) {
			octant = m_tree_3D.getGhostOctant(octantInfo.id);
		} else {
			octant = m_tree_3D.getOctant(octantInfo.id);
		}
		return m_tree_3D.getLevel(octant);
	} else {
		Class_Octant<2> *octant;
		if (m_tree_2D.getIsGhost(octantInfo.id)) {
			octant = m_tree_2D.getGhostOctant(octantInfo.id);
		} else {
			octant = m_tree_2D.getOctant(octantInfo.id);
		}
		return m_tree_2D.getLevel(octant);
	}
}

/*!
	Updates the patch.

	\result Returns true if the mesh was updated, false otherwise.
*/
const std::vector<Adaption::Info> PatchOctree::_update(bool trackAdaption)
{
	if (!is_dirty()) {
		return std::vector<Adaption::Info>();
	}

	// Updating the tree
	std::cout << ">> Adapting tree...";

	bool updated;
	if (is_three_dimensional()) {
		updated = m_tree_3D.adapt(true);
	} else {
		updated = m_tree_2D.adapt(true);
	}

	if (!updated) {
		std::cout << " Already updated" << std::endl;

		return std::vector<Adaption::Info>();
	}

	std::cout << " Done" << std::endl;

	// Info on the tree
	long nOctants;
	if (is_three_dimensional()) {
		nOctants = m_tree_3D.getNumOctants();
	} else {
		nOctants = m_tree_2D.getNumOctants();
	}
	long nPreviousOctants = m_octant_to_cell.size();

	std::cout << ">> Number of octants : " << nOctants << std::endl;

	// Info on the tree
	long nGhostsOctants;
	if (is_three_dimensional()) {
		nGhostsOctants = m_tree_3D.getNumGhosts();
	} else {
		nGhostsOctants = m_tree_2D.getNumGhosts();
	}
	long nPreviousGhosts = m_ghost_to_cell.size();

	// Evaluate tree conenctivity
	std::cout << ">> Evaluating Octree connectivity...";

	if (is_three_dimensional()) {
		m_tree_3D.computeConnectivity();
	} else {
		m_tree_2D.computeConnectivity();
	}

	std::cout << " Done" << std::endl;

	// Initialize intersections
	std::cout << ">> Evaluating Octree intersections...";

	if (is_three_dimensional()) {
		m_tree_3D.computeIntersections();
	} else {
		m_tree_2D.computeIntersections();
	}

	std::cout << " Done" << std::endl;

	// Initialize tracking data
	std::vector<Adaption::Info> adaptionData;

	// Extract information for transforming the patch
	//
	// If there are no cells in the mesh we need to import all
	// octants.
	std::cout << ">> Extract information for transforming the patch...";

	std::vector<OctantInfo> newOctants;
	std::unordered_map<uint32_t, long> renumberedOctants;
	std::vector<long> removedCells;
	std::unordered_set<long> removedInterfaces;

	newOctants.reserve(nOctants + nGhostsOctants);
	renumberedOctants.reserve(nPreviousOctants + nPreviousGhosts);
	removedCells.reserve(nPreviousOctants + nPreviousGhosts);

	uint32_t treeId = 0;
	bool allNew = (get_cell_count() == 0);
	while (treeId < nOctants) {
		// Octant mapping
		std::vector<uint32_t> mapper_octantMap;
		std::vector<bool> mapper_ghostFlag;
		if (is_three_dimensional()) {
			m_tree_3D.getMapping(treeId, mapper_octantMap, mapper_ghostFlag);
		} else {
			m_tree_2D.getMapping(treeId, mapper_octantMap, mapper_ghostFlag);
		}

		// Adaption type
		Adaption::Type adaptionType;
		if (allNew) {
			adaptionType = Adaption::TYPE_CREATION;
		} else {
			bool isNewR;
			if (is_three_dimensional()) {
				isNewR = m_tree_3D.getIsNewR(treeId);
			} else {
				isNewR = m_tree_2D.getIsNewR(treeId);
			}

			if (isNewR) {
				adaptionType = Adaption::TYPE_REFINEMENT;
			} else {
				bool isNewC;
				if (is_three_dimensional()) {
					isNewC = m_tree_3D.getIsNewC(treeId);
				} else {
					isNewC = m_tree_2D.getIsNewC(treeId);
				}

				if (isNewC) {
					adaptionType = Adaption::TYPE_COARSENING;
				} else if (treeId != mapper_octantMap.front()) {
					adaptionType = Adaption::TYPE_RENUMBERING;
				} else {
					++treeId;
					continue;
				}
			}
		}

		// Re-numbered cells just need to be added to the proper list.
		//
		// Renumbered cells are not tracked, because the re-numbering
		// only happens inside PatchMan.
		if (adaptionType == Adaption::TYPE_RENUMBERING) {
			uint32_t previousTreeId = mapper_octantMap.front();
			long cellId = m_octant_to_cell.at(previousTreeId);

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
		if (allNew) {
			nCurrentTreeIds = nOctants - treeId;
		} else if (adaptionType == Adaption::TYPE_REFINEMENT) {
			nCurrentTreeIds = pow(2, get_dimension());
		} else {
			nCurrentTreeIds = 1;
		}

		const long lastCurrentTreeId = treeId + nCurrentTreeIds;
		for (int currentTreeId = treeId; currentTreeId < lastCurrentTreeId; ++currentTreeId) {
			newOctants.emplace_back(currentTreeId, true);
		}

		// Previous cell and interfaces ids that will be removed
		if (adaptionType != Adaption::TYPE_CREATION) {
			auto mapperIter = mapper_octantMap.cbegin();
			while (mapperIter != mapper_octantMap.cend()) {
				const uint32_t &previousTreeId = *mapperIter;

				removedCells.emplace_back();
				long &cellId = removedCells.back();
				cellId = m_octant_to_cell.at(previousTreeId);

				mapperIter++;
			}
		}

		// Adaption tracking
		if (trackAdaption) {
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
			if (adaptionType != Adaption::TYPE_CREATION) {
				int nPreviousCellIds = mapper_octantMap.size();
				adaptionInfo.previous.reserve(nPreviousCellIds);
				auto removedCellsIter = removedCells.cend() - nPreviousCellIds;
				while (removedCellsIter != removedCells.cend()) {
					const long &id = *removedCellsIter;

					adaptionInfo.previous.emplace_back();
					unsigned long &adaptionId = adaptionInfo.previous.back();
					adaptionId = id;

					const Cell &cell = m_cells.at(id);
					long nCellInterfaces = cell.get_interface_count();
					const long *interfaces = cell.get_interfaces();
					for (int k = 0; k < nCellInterfaces; ++k) {
						removedInterfaces.insert(interfaces[k]);
					}

					removedCellsIter++;
				}
			}
		}

		// Incremente tree id
		treeId += nCurrentTreeIds;
	}

	std::cout << " Done" << std::endl;

	// Previous ghosts cells need to be removed
	if (nPreviousGhosts > 0) {
		long deletedGhostsInfoIdx;
		if (trackAdaption) {
			adaptionData.emplace_back();
			Adaption::Info &deletedGhostsInfo = adaptionData.back();
			deletedGhostsInfo.type   = Adaption::TYPE_DELETION;
			deletedGhostsInfo.entity = Adaption::ENTITY_CELL;
			deletedGhostsInfo.previous.reserve(nPreviousGhosts);

			deletedGhostsInfoIdx = adaptionData.size() - 1;
		}

		auto cellIterator = m_cell_to_ghost.cbegin();
		while (cellIterator != m_cell_to_ghost.cend()) {
			long id = cellIterator->first;

			removedCells.emplace_back();
			long &removedId = removedCells.back();
			removedId = id;

			// Adaption tracking
			if (trackAdaption) {
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
	for (uint32_t treeId = 0; treeId < nGhostsOctants; ++treeId) {
		newOctants.emplace_back(treeId, false);
	}
	newOctants.shrink_to_fit();

	// Delete removed cells
	FaceInfoSet danglingFaces;
	if (removedCells.size() > 0) {
		std::cout << ">> Removing non-existing cells...";

		danglingFaces = remove_cells(removedCells);

		std::cout << " Done" << std::endl;
		std::cout << ">> Cells removed: " <<  removedCells.size() << std::endl;
	}

	std::vector<long>().swap(removedCells);

	// Reserve space in the octant-to-cell maps
	m_cell_to_octant.reserve(nOctants);
	m_octant_to_cell.reserve(nOctants);

	// Remap renumbered cells
	if (renumberedOctants.size() > 0) {
		std::cout << ">> Rebuilding octant-to-cell map for renumbered cells...";

		// Remap cells to the new tree ids
		auto cellIterator = renumberedOctants.begin();
		while (cellIterator != renumberedOctants.end()) {
			long cellId = cellIterator->second;

			uint32_t currentTreeId  = cellIterator->first;
			uint32_t previousTreeId = m_cell_to_octant.at(cellId);

			m_cell_to_octant[cellId] = currentTreeId;
			m_octant_to_cell[currentTreeId] = cellId;
			if (renumberedOctants.count(previousTreeId) == 0) {
				m_octant_to_cell.erase(previousTreeId);
			}

			cellIterator++;
		}

		std::cout << " Done" << std::endl;
		std::cout << ">> Cells renumbered: " <<  renumberedOctants.size() << std::endl;
	}

	std::unordered_map<uint32_t, long>().swap(renumberedOctants);

	// Reset ghost maps
	m_cell_to_ghost.clear();
	m_cell_to_ghost.reserve(nGhostsOctants);

	m_ghost_to_cell.clear();
	m_ghost_to_cell.reserve(nGhostsOctants);

	// Import added octants
	std::vector<unsigned long> createdInterfaces;
	if (newOctants.size() > 0) {
		std::cout << ">> Importing new octants...";

		createdInterfaces = import_octants(newOctants, danglingFaces);

		std::cout << " Done" << endl;
		std::cout << ">> Octants imported: " <<  newOctants.size() << std::endl;
	}

	FaceInfoSet().swap(danglingFaces);

	// Track mesh adaption
	if (trackAdaption) {
		// Map ids of the added cells
		for (auto &adaptionInfo : adaptionData) {
			if (adaptionInfo.entity != Adaption::ENTITY_CELL) {
				continue;
			}

			int nCurrentIds = adaptionInfo.current.size();
			for (int k = 0; k < nCurrentIds; ++k) {
				long cellId = m_octant_to_cell.at(adaptionInfo.current[k]);
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
			auto cellIterator = m_cell_to_ghost.cbegin();
			while (cellIterator != m_cell_to_ghost.cend()) {
				createdGhostsInfo.previous.emplace_back();
				unsigned long &adaptionId = createdGhostsInfo.previous.back();
				adaptionId = cellIterator->first;

				cellIterator++;
			}
		}

		// Track deleted interfaces
		if (removedInterfaces.size() > 0) {
			adaptionData.emplace_back();
			Adaption::Info &deletedInterfacesInfo = adaptionData.back();
			deletedInterfacesInfo.type   = Adaption::TYPE_DELETION;
			deletedInterfacesInfo.entity = Adaption::ENTITY_INTERFACE;

			deletedInterfacesInfo.previous.reserve(removedInterfaces.size());
			auto interfaceIterator = removedInterfaces.begin();
			while (interfaceIterator != removedInterfaces.end()) {
				deletedInterfacesInfo.previous.emplace_back();
				unsigned long &adaptionId = deletedInterfacesInfo.previous.back();
				adaptionId = *interfaceIterator;

				interfaceIterator++;
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
	if (is_three_dimensional()) {
		m_tree_3D.clearConnectivity();
	} else {
		m_tree_2D.clearConnectivity();
	}

	// Done
	return adaptionData;
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
std::vector<unsigned long> PatchOctree::import_octants(std::vector<OctantInfo> &octantInfoList)
{
	FaceInfoSet danglingFaces;

	return import_octants(octantInfoList, danglingFaces);
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
std::vector<unsigned long> PatchOctree::import_octants(std::vector<OctantInfo> &octantInfoList,
                                 FaceInfoSet &danglingFaces)
{
	// Info of the cells
	Element::Type cellType;
	if (is_three_dimensional()) {
		cellType = Element::BRICK;
	} else {
		cellType = Element::RECTANGLE;
	}

	int nCellFaces;
	if (is_three_dimensional()) {
		nCellFaces = Element::get_face_count(Element::BRICK);
	} else {
		nCellFaces = Element::get_face_count(Element::RECTANGLE);
	}

	int nCellVertices;
	if (is_three_dimensional()) {
		nCellVertices = Element::get_vertex_count(Element::BRICK);
	} else {
		nCellVertices = Element::get_vertex_count(Element::RECTANGLE);
	}

	std::vector<std::vector<int>> cellLocalFaceConnect(nCellFaces);
	for (int k = 0; k < nCellFaces; ++k) {
		cellLocalFaceConnect[k] = Element::get_face_local_connect(cellType, k);
	}

	// Info on the interfaces
	int nInterfaceVertices;
	if (is_three_dimensional()) {
		nInterfaceVertices = Element::get_face_count(Element::RECTANGLE);
	} else {
		nInterfaceVertices = Element::get_face_count(Element::LINE);
	}

	uint32_t nIntersections;
	if (is_three_dimensional()) {
		nIntersections = m_tree_3D.getNumIntersections();
	} else {
		nIntersections = m_tree_2D.getNumIntersections();
	}

	std::vector<unsigned long> createdInterfaces;
	createdInterfaces.reserve(nCellFaces * octantInfoList.size());

	// Add the vertex of the dangling faces to the vertex map
	std::unordered_map<uint32_t, long> vertexMap;
	for (auto &danglingFaceInfo : danglingFaces) {
		// List of faces with the vertx to be added
		long danglingId = danglingFaceInfo.id;
		Cell &danglingCell = m_cells[danglingId];
		int danglingFace = danglingFaceInfo.face;

		int nInterfaces = danglingCell.get_interface_count(danglingFace);
		std::vector<FaceInfo> vertexSourceList(1 + nInterfaces);
		vertexSourceList[0] = danglingFaceInfo;
		for (int k = 0; k < nInterfaces; ++k) {
			FaceInfo vertexSource;
			long interfaceId = danglingCell.get_interface(danglingFace, k);

			Interface &interface = m_interfaces[interfaceId];
			if (interface.get_owner() != danglingId) {
				vertexSource.id   = interface.get_owner();
				vertexSource.face = interface.get_owner_face();
			} else {
				vertexSource.id   = interface.get_neigh();
				vertexSource.face = interface.get_neigh_face();
			}
			vertexSourceList[1 + k] = vertexSource;
		}

		// Add the vertices to the map
		for (auto & vertexSource : vertexSourceList) {
			// Cell data
			Cell &cell = m_cells[vertexSource.id];
			const long *cellConnect = cell.get_connect();

			// Octant data
			OctantInfo octantInfo = get_cell_octant(vertexSource.id);
			const std::vector<uint32_t> &octantTreeConnect = get_octant_connect(octantInfo);

			// List of vertices
			std::vector<int> &localConnect = cellLocalFaceConnect[vertexSource.face];
			for (int k = 0; k < nInterfaceVertices; ++k) {
				long vertexId = cellConnect[localConnect[k]];
				uint32_t vertexTreeId = octantTreeConnect[localConnect[k]];

				vertexMap.insert({{vertexTreeId, vertexId}});
			}
		}
	}

	// Create the new vertices
	for (OctantInfo &octantInfo : octantInfoList) {
		const std::vector<uint32_t> &octantTreeConnect = get_octant_connect(octantInfo);
		for (int k = 0; k < nCellVertices; ++k) {
			uint32_t vertexTreeId = octantTreeConnect[k];
			if (vertexMap.count(vertexTreeId) == 0) {
				vertexMap[vertexTreeId] = create_vertex(vertexTreeId);
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
		vector<uint32_t> cells;
		int owner;
		int ownerFace;
		bool isBoundary;
		bool isGhost;
		Class_Intersection<2> *treeInterface_2D;
		Class_Intersection<3> *treeInterface_3D;
		if (is_three_dimensional()) {
			treeInterface_3D = m_tree_3D.getIntersection(interfaceTreeId);

			cells      = m_tree_3D.getOwners(treeInterface_3D);
			owner      = m_tree_3D.getFiner(treeInterface_3D);
			ownerFace  = m_tree_3D.getFace(treeInterface_3D);
			isBoundary = m_tree_3D.getBound(treeInterface_3D);
			isGhost    = m_tree_3D.getIsGhost(treeInterface_3D);
		} else {
			treeInterface_2D = m_tree_2D.getIntersection(interfaceTreeId);

			cells      = m_tree_2D.getOwners(treeInterface_2D);
			owner      = m_tree_2D.getFiner(treeInterface_2D);
			ownerFace  = m_tree_2D.getFace(treeInterface_2D);
			isBoundary = m_tree_2D.getBound(treeInterface_2D);
			isGhost    = m_tree_2D.getIsGhost(treeInterface_2D);
		}

		// Decide if we need to build the interface
		bool createInterface = false;

		OctantInfo ownerOctantInfo(cells[owner], true);
		long ownerId = get_octant_id(ownerOctantInfo);
		if (ownerId == Element::NULL_ELEMENT_ID) {
			octantTreeInterfaces[ownerOctantInfo.id].push_back(interfaceTreeId);
			createInterface = true;
		}

		long neighId = Element::NULL_ELEMENT_ID;
		if (!isBoundary) {
			OctantInfo neighOctantInfo(cells[!owner], !isGhost);
			neighId = get_octant_id(neighOctantInfo);
			if (neighId == Element::NULL_ELEMENT_ID) {
				octantTreeInterfaces[neighOctantInfo.id].push_back(interfaceTreeId);
				createInterface = true;
			}
		}

		if (!createInterface) {
			continue;
		}

		// Interface owner and neighbour faces
		//
		// The owner or the neighbour may be not know yet, but the
		// corresponding faces are. So, set the faces also for
		// the unknown cells because they are needed when finding
		// the interfaces associated to a cell.
		std::array<FaceInfo, 2> interfaceFaces;

		interfaceFaces[0].id = ownerId;
		interfaceFaces[0].face = ownerFace;

		interfaceFaces[1].id = neighId;
		interfaceFaces[1].face = ownerFace + 1 - 2 * (ownerFace % 2);

		// Interface connectivity
		const std::vector<uint32_t> &octantTreeConnect = get_octant_connect(ownerOctantInfo);
		std::vector<int> &localConnect = cellLocalFaceConnect[ownerFace];
		std::unique_ptr<long[]> interfaceConnect = std::unique_ptr<long[]>(new long[nInterfaceVertices]);
		for (int k = 0; k < nInterfaceVertices; ++k) {
			interfaceConnect[k] = vertexMap.at(octantTreeConnect[localConnect[k]]);
		}

		// Create the interface
		createdInterfaces.emplace_back();
		unsigned long &interfaceId = createdInterfaces.back();
		interfaceId = create_interface(interfaceTreeId, interfaceConnect, interfaceFaces);
		interfaceMap[interfaceTreeId] = interfaceId;

		// If the interface is on an dangling faces, the owner or
		// the neigbour of the interface already exists. We need
		// to identify that cell and set that information on the
		// interface.
		//
		// A boundary interface cannnot be on a dangling face.
		if (!isBoundary) {
			FaceInfo guessDanglingInfo;
			if (ownerId == Element::NULL_ELEMENT_ID) {
				guessDanglingInfo = interfaceFaces[1];
			} else {
				guessDanglingInfo = interfaceFaces[0];
			}

			if (danglingFaces.count(guessDanglingInfo) != 0) {
				Cell &danglingCell = m_cells[guessDanglingInfo.id];
				danglingCell.push_interface(guessDanglingInfo.face, interfaceId);
			}
		}
	}

	// Add the cells
	std::vector<std::vector<long>> cellInterfaces(nCellFaces, std::vector<long>());
	std::vector<std::vector<bool>> interfaceOwnerFlags(nCellFaces, std::vector<bool>());
	for (OctantInfo &octantInfo : octantInfoList) {
		// Octant connectivity
		const std::vector<uint32_t> &octantTreeConnect = get_octant_connect(octantInfo);

		// Cell connectivity
		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int k = 0; k < nCellVertices; ++k) {
			uint32_t vertexTreeId = octantTreeConnect[k];
			cellConnect[k] = vertexMap.at(vertexTreeId);
		}

		// Add interfaces
		for (int k = 0; k < nCellFaces; ++k) {
			cellInterfaces[k].clear();
			interfaceOwnerFlags[k].clear();
		}

		for (uint32_t interfaceTreeId : octantTreeInterfaces[octantInfo.id]) {
			const long &interfaceId = interfaceMap.at(interfaceTreeId);
			const Interface &interface = m_interfaces[interfaceId];

			vector<uint32_t> cells;
			int owner;
			Class_Intersection<2> *treeInterface_2D;
			Class_Intersection<3> *treeInterface_3D;
			if (is_three_dimensional()) {
				treeInterface_3D = m_tree_3D.getIntersection(interfaceTreeId);

				cells      = m_tree_3D.getOwners(treeInterface_3D);
				owner      = m_tree_3D.getFiner(treeInterface_3D);
			} else {
				treeInterface_2D = m_tree_2D.getIntersection(interfaceTreeId);

				cells      = m_tree_2D.getOwners(treeInterface_2D);
				owner      = m_tree_2D.getFiner(treeInterface_2D);
			}
			bool ownerFlag = (cells[owner] == octantInfo.id);

			int cellFace;
			if (ownerFlag) {
				cellFace = interface.get_owner_face();
			} else {
				cellFace = interface.get_neigh_face();
			}

			cellInterfaces[cellFace].emplace_back();
			long &cellInterfaceId = cellInterfaces[cellFace].back();
			cellInterfaceId = interfaceId;

			interfaceOwnerFlags[cellFace].push_back(ownerFlag);
		}

		// Add cell
		create_cell(octantInfo.id, octantInfo.internal, cellConnect, cellInterfaces, interfaceOwnerFlags);
	}

	// Done
	createdInterfaces.shrink_to_fit();

	return createdInterfaces;
}

/*!
	Remove a list of octants from the patch.

	\param octantTreeIds is the list of octant ids to remove
*/
PatchOctree::FaceInfoSet PatchOctree::remove_cells(std::vector<long> &cellIds)
{
	// Delete cells
	//
	// To make deletion of interfaces more efficient, keep the list of
	// faces to remove in an ordered map. The order should be based upon
	// the position of the interfaces in the container, rather than upon
	// the id of the faces.
	typedef PiercedVector<Interface>::position_less InterfacePositionCompare;
	typedef std::map<long, int, InterfacePositionCompare> InterfaceOrderedMap;

	std::unordered_set<long> deadVertices;
	InterfaceOrderedMap deadInterfaces((InterfacePositionCompare(m_interfaces)));
	for (long cellId : cellIds) {
		Cell &cell = m_cells[cellId];

		// List dead vertices
		//
		// For now, all cell vertices will be listed. Later, the vertex of
		// the dangling faces will be removed from the list.
		int nCellVertices = cell.get_vertex_count();
		const long *vertices = cell.get_connect();
		for (int k = 0; k < nCellVertices; ++k) {
			deadVertices.insert(vertices[k]);
		}

		// List dead interface and set the dangling status. An interface
		// is dangling if only the owner or the neighbour will be deleted.
		int nCellInterfaces = cell.get_interface_count();
		const long *interfaces = cell.get_interfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];

			int danglingSide = -1;
			if (deadInterfaces.count(interfaceId) == 0) {
				Interface &interface = m_interfaces[interfaceId];
				if (interface.get_position_type() != Interface::BOUNDARY) {
					if (interface.get_owner() == cellId) {
						danglingSide = 1;
					} else {
						danglingSide = 0;
					}
				}
			}
			deadInterfaces[interfaceId] = danglingSide;
		}

		// Delete cell
		delete_cell(cellId);
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
			long danglingCellFace;
			if (danglingSide == 0) {
				danglingCellId   = interface.get_owner();
				danglingCellFace = interface.get_owner_face();
			} else {
				danglingCellId   = interface.get_neigh();
				danglingCellFace = interface.get_neigh_face();
			}

			// Remove interface vertices from dead vertices
			int nFaceVertices = interface.get_vertex_count();
			const long *vertices = interface.get_connect();
			for (int k = 0; k < nFaceVertices; ++k) {
				deadVertices.erase(vertices[k]);
			}

			// Remove interface from dangling cell
			Cell &danglingCell = m_cells[danglingCellId];

			int j = 0;
			while (danglingCell.get_interface(danglingCellFace, j) != interfaceId) {
				++j;
			}
			danglingCell.delete_interface(danglingCellFace, j);

			// Add the associated cell face to the dangling faces list
			FaceInfo danglingFace;
			danglingFace.id   = danglingCellId;
			danglingFace.face = danglingCellFace;

			danglingFaces.insert(danglingFace);
		}

		// Add the interface to the list of interfaces to delete
		Patch::delete_interface(interfaceId, true);
	}

	// Delete vertices
	for (auto it = deadVertices.begin(); it != deadVertices.end(); ++it) {
		Patch::delete_vertex(*it, true);
	}

	// Done
	return danglingFaces;
}

/*!
	Creates a new patch vertex from the specified tree vertex.

	\param treeId is the id of the vertex in the tree
*/
long PatchOctree::create_vertex(uint32_t treeId)
{
	// Create the vertex
	long id = Patch::create_vertex();
	Node &vertex = m_vertices[id];

	// Coordinate
	vector<double> nodeCoords;
	if (is_three_dimensional()) {
		nodeCoords = m_tree_3D.getNodeCoordinates(treeId);
	} else {
		nodeCoords = m_tree_2D.getNodeCoordinates(treeId);
	}

	std::array<double, 3> coords;
	std::copy_n(nodeCoords.data(), coords.size(), coords.data());
	vertex.set_coords(coords);

	// Done
	return id;
}

/*!
	Creates a new patch interface from the specified tree intersection.

	\param treeId is the id of the intersection in the tree
*/
long PatchOctree::create_interface(uint32_t treeId,
                                   std::unique_ptr<long[]> &vertices,
                                   std::array<FaceInfo, 2> &faces)
{
	// Create the interface
	long id = Patch::create_interface();
	Interface &interface = m_interfaces[id];

	// Info associate al tree
	int level;
	bool isGhost;
	bool isBoundary;
	int ownerFace;
	vector<double> faceCenter;
	if (is_three_dimensional()) {
		Class_Intersection<3> *treeInterface = m_tree_3D.getIntersection(treeId);

		level      = m_tree_3D.getLevel(treeInterface);
		isGhost    = m_tree_3D.getPbound(treeInterface);
		isBoundary = m_tree_3D.getBound(treeInterface);
		ownerFace  = m_tree_3D.getFace(treeInterface);
		faceCenter = m_tree_3D.getCenter(treeInterface);
	} else {
		Class_Intersection<2> *treeInterface = m_tree_2D.getIntersection(treeId);

		level      = m_tree_2D.getLevel(treeInterface);
		isGhost    = m_tree_2D.getPbound(treeInterface);
		isBoundary = m_tree_2D.getBound(treeInterface);
		ownerFace  = m_tree_2D.getFace(treeInterface);
		faceCenter = m_tree_2D.getCenter(treeInterface);
	}

	// Tipo
	if (is_three_dimensional()) {
		interface.set_type(Element::RECTANGLE);
	} else {
		interface.set_type(Element::LINE);
	}

	// Area
	interface.set_area(&m_tree_area[level]);

	// Normal
	interface.set_normal(&m_normals[ownerFace]);

	// Centroid
	std::array<double, 3> centroid;
	for (unsigned int k = 0; k < centroid.size(); k++) {
		centroid[k] = faceCenter[k];
	}

	interface.set_centroid(centroid);

	// Position
	if (isGhost) {
		interface.set_position_type(Interface::GHOST);
	} else if (isBoundary) {
		interface.set_position_type(Interface::BOUNDARY);
	} else {
		interface.set_position_type(Interface::INTERNAL);
	}

	// Connectivity
	interface.set_connect(std::move(vertices));

	// Owner and neighbour
	interface.set_owner(faces[0].id, faces[0].face);
	interface.set_neigh(faces[1].id, faces[1].face);

	// Done
	return id;
}

/*!
	Creates a new patch cell from the specified tree octant.

	\param treeId is the id of the octant in the tree
*/
long PatchOctree::create_cell(uint32_t treeId, bool internal,
                              std::unique_ptr<long[]> &vertices,
                              std::vector<std::vector<long>> &interfaces,
                              std::vector<std::vector<bool>> &ownerFlags)
{
	// Create the cell
	long id = Patch::create_cell(internal);
	Cell &cell = m_cells[id];

	// Octant info
	int octantLevel;
	vector<double> octantCentroid;
	if (is_three_dimensional()) {
		Class_Octant<3> *octant;
		if (internal) {
			octant = m_tree_3D.getOctant(treeId);
		} else {
			octant = m_tree_3D.getGhostOctant(treeId);
		}
		octantLevel    = m_tree_3D.getLevel(octant);
		octantCentroid = m_tree_3D.getCenter(octant);
	} else {
		Class_Octant<2> *octant;
		if (internal) {
			octant = m_tree_2D.getOctant(treeId);
		} else {
			octant = m_tree_2D.getGhostOctant(treeId);
		}
		octantLevel    = m_tree_2D.getLevel(octant);
		octantCentroid = m_tree_2D.getCenter(octant);
	}

	// Tipo
	if (is_three_dimensional()) {
		cell.set_type(Element::BRICK);
	} else {
		cell.set_type(Element::RECTANGLE);
	}

	// Position
	if (internal) {
		cell.set_position_type(Cell::INTERNAL);
	} else {
		cell.set_position_type(Cell::GHOST);
	}

	// Volume
	cell.set_volume(&m_tree_volume[octantLevel]);

	// Centroid
	std::array<double, 3> centroid;
	for (unsigned int k = 0; k < centroid.size(); k++) {
		centroid[k] = octantCentroid[k];
	}

	cell.set_centroid(centroid);

	// Connectivity
	cell.set_connect(std::move(vertices));

	// Interfaces
	int nCellFaces = interfaces.size();
	for (int face = 0; face < nCellFaces; ++face) {
		int nFaceInterfaces = interfaces[face].size();
		for (int k = 0; k < nFaceInterfaces; ++k) {
			long interfaceId   = interfaces[face][k];
			bool ownsInterface = ownerFlags[face][k];

			pman::Interface &interface = m_interfaces[interfaceId];
			if (ownsInterface) {
				interface.set_owner(id, face);
			} else {
				interface.set_neigh(id, face);
			}
		}
	}
	cell.initialize_interfaces(interfaces);

	// Update cell to octant mapping
	if (internal) {
		m_cell_to_octant.insert({{id,treeId}});
		m_octant_to_cell.insert({{treeId,id}});
	} else {
		m_cell_to_ghost.insert({{id,treeId}});
		m_ghost_to_cell.insert({{treeId,id}});
	}

	// Done
	return id;
}

/*!
	Deletes a cell from the patch.

	\param id is the id of the cell
*/
void PatchOctree::delete_cell(long id)
{
	// Remove the information that link the cell to the octant
	bool internal = (m_cells[id].get_position_type() == Cell::INTERNAL);

	std::unordered_map<long, uint32_t, Element::IdHasher> *cellMap;
	if (internal) {
		cellMap = &m_cell_to_octant;
	} else {
		cellMap = &m_cell_to_ghost;
	}

	std::unordered_map<long, uint32_t, Element::IdHasher>::const_iterator cellItr = cellMap->find(id);
	if (cellItr != cellMap->end()) {
		// Delete octant-to-cell entry
		std::unordered_map<uint32_t, long> *octantMap;
		if (internal) {
			octantMap = &m_octant_to_cell;
		} else {
			octantMap = &m_ghost_to_cell;
		}

		uint32_t treeId = cellItr->second;
		octantMap->erase(treeId);

		// Delete cell-to-octant entry
		cellMap->erase(cellItr);
	}

	// Delete the cell
	Patch::delete_cell(id, true);
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
bool PatchOctree::_mark_cell_for_refinement(const long &id)
{
	return set_marker(id, 1);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
bool PatchOctree::_mark_cell_for_coarsening(const long &id)
{
	return set_marker(id, -1);
}

/*!
	Set the marker on a cell.

	\param id is the id of the cell
	\param value is the value of the marker
*/
bool PatchOctree::set_marker(const long &id, const int8_t &value)
{
	OctantInfo octantInfo = get_cell_octant(id);
	if (!octantInfo.internal) {
		return false;
	}

	if (is_three_dimensional()) {
		m_tree_3D.setMarker(octantInfo.id, value);
	} else {
		m_tree_2D.setMarker(octantInfo.id, value);
	}

	return true;
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
bool PatchOctree::_enable_cell_balancing(const long &id, bool enabled)
{
	OctantInfo octantInfo = get_cell_octant(id);
	if (!octantInfo.internal) {
		return false;
	}

	if (is_three_dimensional()) {
		m_tree_3D.setBalance(octantInfo.id, enabled);
	} else {
		m_tree_2D.setBalance(octantInfo.id, enabled);
	}

	return true;
}

}

#endif
