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

	m_tree_dh.resize(maxLevels);
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
		const std::unordered_map<long, uint32_t, Element::IdHasher>::const_iterator cellItr = m_cell_to_octant.find(id);
		if (cellItr == m_cell_to_octant.end()) {
			octantInfo.exists = false;
			return octantInfo;
		}
		octantInfo.id = cellItr->second;
	} else {
		const std::unordered_map<long, uint32_t, Element::IdHasher>::const_iterator cellItr = m_cell_to_ghost.find(id);
		if (cellItr == m_cell_to_ghost.end()) {
			octantInfo.exists = false;
			return octantInfo;
		}
		octantInfo.id = cellItr->second;
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
bool PatchOctree::_update(vector<uint32_t> &cellMapping)
{
	if (!is_dirty()) {
		return false;
	}

	std::cout << ">> Updating the mesh\n";

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

		return updated;
	}

	std::cout << " Done" << std::endl;

	// Info on the tree
	long nOctants;
	if (is_three_dimensional()) {
		nOctants = m_tree_3D.getNumOctants();
	} else {
		nOctants = m_tree_2D.getNumOctants();
	}

	std::cout << ">> Number of octants : " << nOctants << std::endl;

	// Info on the tree
	long nGhostsOctants;
	if (is_three_dimensional()) {
		nGhostsOctants = m_tree_3D.getNumGhosts();
	} else {
		nGhostsOctants = m_tree_2D.getNumGhosts();
	}

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

	// Import octants
	std::cout << ">> Importing octants...";

	m_cells.clear();
	m_interfaces.clear();
	m_vertices.clear();

	m_cell_to_octant.clear();
	m_cell_to_octant.reserve(nOctants);

	m_cell_to_ghost.clear();
	m_cell_to_ghost.reserve(nGhostsOctants);

	m_octant_to_cell.clear();
	m_octant_to_cell.reserve(nOctants);

	m_ghost_to_cell.clear();
	m_ghost_to_cell.reserve(nGhostsOctants);

	std::vector<OctantInfo> addedOctants;

	for (uint32_t treeId = 0; treeId < nOctants; ++treeId) {
		OctantInfo octantInfo;
		octantInfo.id       = treeId;
		octantInfo.internal = true;

		addedOctants.push_back(octantInfo);
	}

	for (uint32_t treeId = 0; treeId < nGhostsOctants; ++treeId) {
		OctantInfo octantInfo;
		octantInfo.id       = treeId;
		octantInfo.internal = false;

		addedOctants.push_back(octantInfo);
	}

	import_octants(addedOctants);

	std::cout << " Done" << endl;
	std::cout << ">> Octants imported: " <<  addedOctants.size() << std::endl;

	// Delete tree conenctivity
	if (is_three_dimensional()) {
		m_tree_3D.clearConnectivity();
	} else {
		m_tree_2D.clearConnectivity();
	}

	return updated;
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
void PatchOctree::import_octants(std::vector<OctantInfo> &octantInfoList)
{
	FaceInfoSet danglingFaces;

	import_octants(octantInfoList, danglingFaces);
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
*/
void PatchOctree::import_octants(std::vector<OctantInfo> &octantInfoList,
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

	// Add the vertex of the dangling faces to the vertex map
	std::unordered_map<uint32_t, long> vertexMap;
	for (auto &danglingFaceInfo : danglingFaces) {
		// List of faces with the vertx to be added
		long danglingId = danglingFaceInfo.id;
		Cell &danglingCell = m_cells[danglingId];
		int danglingFace = danglingFaceInfo.face;

		std::vector<FaceInfo> vertexSourceList;
		vertexSourceList.push_back(danglingFaceInfo);
		for (int k = 0; k < danglingCell.get_interface_count(danglingFace); ++k) {
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
			vertexSourceList.push_back(vertexSource);
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
			std::unordered_map<uint32_t, long>::const_iterator vertexItr = vertexMap.find(vertexTreeId);
			if (vertexItr == vertexMap.end()) {
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
		std::unordered_map<uint32_t, long>::const_iterator interfaceItr = interfaceMap.find(interfaceTreeId);
		if (interfaceItr != interfaceMap.end()) {
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

		OctantInfo ownerOctantInfo;
		ownerOctantInfo.id       = cells[owner];
		ownerOctantInfo.internal = true;

		long ownerId = get_octant_id(ownerOctantInfo);
		if (ownerId == Element::NULL_ELEMENT_ID) {
			octantTreeInterfaces[ownerOctantInfo.id].push_back(interfaceTreeId);
			createInterface = true;
		}

		long neighId = Element::NULL_ELEMENT_ID;
		if (!isBoundary) {
			OctantInfo neighOctantInfo;
			neighOctantInfo.id       = cells[owner ? 0 : 1];
			neighOctantInfo.internal = !isGhost;

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
		long interfaceId = create_interface(interfaceTreeId, interfaceConnect, interfaceFaces);
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

			FaceInfoSet::const_iterator guessItr = danglingFaces.find(guessDanglingInfo);
			if (guessItr != danglingFaces.end()) {
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

			cellInterfaces[cellFace].push_back(interfaceId);
			interfaceOwnerFlags[cellFace].push_back(ownerFlag);
		}

		// Add cell
		create_cell(octantInfo.id, octantInfo.internal, cellConnect, cellInterfaces, interfaceOwnerFlags);
	}
}

/*!
	Creates a new patch vertex from the specified tree vertex.

	\param treeId is the id of the vertex in the tree
*/
long PatchOctree::create_vertex(uint32_t treeId)
{
	long id = m_vertices.size();

	// Create the vertex
	m_vertices.emplace(id);
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
	long id = m_interfaces.size();

	// Create the interface
	m_interfaces.emplace(id, this);
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
	long id = m_cells.size();

	// Create the cell
	if (internal) {
		m_cells.emplace(id, this);
	} else {
		m_cells.emplace_back(id, this);
	}
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
