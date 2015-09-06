//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\class PatchOctree

	\brief The PatchOctree defines a Octree patch.

	PatchOctree defines a Octree patch.
*/

#ifndef DISABLE_OCTREE
#include "patch_octree.hpp"

#include <math.h>

namespace pman {

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
	double initial_level = (int) log2(max(1., length / dh)) ;
	if (is_three_dimensional()) {
		m_tree_3D = Class_Para_Tree<3>(origin[0], origin[1], origin[2], length);
		m_tree_3D.setMarker((uint32_t) 0, initial_level);
	} else {
		m_tree_2D = Class_Para_Tree<2>(origin[0], origin[1], origin[2], length);
		m_tree_2D.setMarker((uint32_t) 0, initial_level);
	}

	// Info sulle interfacce
	m_normals = unique_ptr<CollapsedArray2D<double> >(new CollapsedArray2D<double>(2 * dimension, 2 * dimension * dimension));
	for (int n = -1; n <= 1; n += 2) {
		for (int i = 0; i < dimension; i++) {
			vector<double> normal(dimension, 0);
			normal[i] = 1.0;

			m_normals->push_back(dimension, normal.data());
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
double * PatchOctree::_get_opposite_normal(double *normal)
{
	int dimension = get_dimension();

	int id_current  = std::distance(m_normals->data(), normal) / dimension;
	int id_opposite = (id_current + dimension) % (2 * dimension);

	return m_normals->get(id_opposite);
}

/*!
	Gets the octant of the cell with the specified id.

	\param id the id of the cell
	\result The octant of the specified cell
*/
int PatchOctree::get_cell_octant(const int &id) const
{
	return m_cell_to_octant.at(id);
}

/*!
	Updates the patch.
*/
void PatchOctree::_update(const vector<uint32_t> &cellMapping)
{
	if (!is_dirty()) {
		return;
	}

	std::cout << ">> Updating the mesh\n";

	// Updating the tree
	if (is_three_dimensional()) {
		m_tree_3D.adapt();
	} else {
		m_tree_2D.adapt();
	}

	// Evaluate tree conenctivity
	if (is_three_dimensional()) {
		m_tree_3D.computeConnectivity();
	} else {
		m_tree_2D.computeConnectivity();
	}

	// Definition of mesh
	update_vertices();
	update_cells();
	update_interfaces();

	// Delete tree conenctivity
	if (is_three_dimensional()) {
		m_tree_3D.clearConnectivity();
	} else {
		m_tree_2D.clearConnectivity();
	}
}

/*!
	Updates the vertices.
*/
void PatchOctree::update_vertices()
{
	if (!is_dirty()) {
		return;
	}

	if (m_vertices.empty()) {
		import_vertices();
	} else {
		reload_vertices();
	}
}

/*!
	Imports the vertices from the tree.
*/
void PatchOctree::import_vertices()
{
	std::cout << ">> Importing the vertices\n";

	// Reset
	reset_vertices();

	// Defintion of the vertices
	int nVertices;
	if (is_three_dimensional()) {
		nVertices = m_tree_3D.getNumNodes();
	} else {
		nVertices = m_tree_2D.getNumNodes();
	}

	m_vertices.reserve(nVertices);
	for (int n = 0; n < nVertices; n++) {
		m_vertices.emplace_back(n);
		Node &vertex = m_vertices.back();

		// Coordinate
		vector<double> nodeCoords;
		if (is_three_dimensional()) {
			nodeCoords = m_tree_3D.getNodeCoordinates(n);
		} else {
			nodeCoords = m_tree_2D.getNodeCoordinates(n);
		}

		std::unique_ptr<double[]> coords = std::unique_ptr<double[]>(new double[get_dimension()]);
		coords[Node::COORD_X] = nodeCoords[Node::COORD_X];
		coords[Node::COORD_Y] = nodeCoords[Node::COORD_Y];
		if (is_three_dimensional()) {
			coords[Node::COORD_Z] = nodeCoords[Node::COORD_Z];
		}

		vertex.set_coords(std::move(coords));
	}
}

/*!
	Updates the vertices after a tree change.
*/
void PatchOctree::reload_vertices()
{
	std::cout << ">> Updating the vertices\n";

	reset_vertices();
	import_vertices();
}

/*!
	Updates the vertices.
*/
void PatchOctree::update_cells()
{
	if (!is_dirty()) {
		return;
	}

	if (m_cells.empty()) {
		import_cells();
	} else {
		reload_cells();
	}
}

/*!
	Imports the cells from the tree.
*/
void PatchOctree::import_cells()
{
	std::cout << ">> Importing the cells\n";

	// Reset
	reset_cells();

	// Info of the cells
	int nCellVertices;
	if (is_three_dimensional()) {
		nCellVertices = Element::get_vertex_count(Element::HEXAHEDRON);
	} else {
		nCellVertices = Element::get_vertex_count(Element::QUADRANGLE);
	}

	// Mapping della connettività
	vector<int> connectMap(nCellVertices);
	if (is_three_dimensional()) {
		connectMap[0] = 0;
		connectMap[1] = 1;
		connectMap[2] = 3;
		connectMap[3] = 2;
		connectMap[4] = 4;
		connectMap[5] = 5;
		connectMap[6] = 7;
		connectMap[7] = 6;
	} else {
		connectMap[0] = 0;
		connectMap[1] = 1;
		connectMap[2] = 3;
		connectMap[3] = 2;
	}

	// Importing the cells
	int m_nInternalCells;
	int m_nGhostsCells;
	if (is_three_dimensional()) {
		m_nInternalCells = m_tree_3D.getNumOctants();
		m_nGhostsCells   = m_tree_3D.getNumGhosts();
	} else {
		m_nInternalCells = m_tree_2D.getNumOctants();
		m_nGhostsCells   = m_tree_2D.getNumGhosts();
	}

	int nCells = m_nInternalCells + m_nGhostsCells;

	m_cell_to_octant.reserve(nCells);
	m_cells.reserve(nCells);
	for (int n = 0; n < nCells; n++) {
		m_cells.emplace_back(n, this);
		Cell &cell = m_cells.back();
		m_cell_to_octant[n] = n;

		// Distinguo tra celle interne e celle ghost
		bool isInternal = (n < m_nInternalCells);

		// Info associate al tree
		int octantId = n;
		if (!isInternal) {
			octantId -= m_nInternalCells;
		}

		int octantLevel;
		vector<double> octantCentroid;
		if (is_three_dimensional()) {
			Class_Octant<3> *octant;
			if (isInternal) {
				octant = m_tree_3D.getOctant(octantId) ;
			} else {
				octant = m_tree_3D.getGhostOctant(octantId) ;
			}

			octantLevel    = m_tree_3D.getLevel(octant);
			octantCentroid = m_tree_3D.getCenter(octant);
		} else {
			Class_Octant<2> *octant;
			if (isInternal) {
				octant = m_tree_2D.getOctant(octantId) ;
			} else {
				octant = m_tree_2D.getGhostOctant(octantId) ;
			}

			octantLevel    = m_tree_2D.getLevel(octant);
			octantCentroid = m_tree_2D.getCenter(octant);
		}

		// Tipo
		if (is_three_dimensional()) {
			cell.set_type(Element::HEXAHEDRON);
		} else {
			cell.set_type(Element::QUADRANGLE);
		}

		// Position
		if (isInternal) {
			cell.set_position_type(Cell::INTERNAL);
		} else {
			cell.set_position_type(Cell::GHOST);
		}

		// Volume
		cell.set_volume(&m_tree_volume[octantLevel]);

		// Centroide
		std::unique_ptr<double[]> centroid = std::unique_ptr<double[]>(new double[get_dimension()]);
		for (int k = 0; k < get_dimension(); k++) {
			centroid[k] = octantCentroid[k];
		}

		cell.set_centroid(std::move(centroid));

		// Connettività
		vector<uint32_t> octantConnect;
		if (is_three_dimensional()) {
			octantConnect = m_tree_3D.getOctantConnectivity(n) ;
		} else {
			octantConnect = m_tree_2D.getOctantConnectivity(n) ;
		}

		std::unique_ptr<int[]> connect = std::unique_ptr<int[]>(new int[nCellVertices]);
		for (int k = 0; k < nCellVertices; k++) {
			connect[k] = octantConnect[connectMap[k]];
		}

		cell.set_connect(std::move(connect));

		// Struttura dati per le interfacce
		// cell.initialize_interfaces(nInterfacesForSide);
	}
}

/*!
	Updates the cells after a tree change.
*/
void PatchOctree::reload_cells()
{
	std::cout << ">> Updating the cells\n";

	reset_cells();
	import_cells();
}


/*!
	Updates the interfaces.
*/
void PatchOctree::update_interfaces()
{
	if (!is_dirty()) {
		return;
	}

	if (m_interfaces.empty()) {
		import_interfaces();
	} else {
		update_interfaces();
	}
}

/*!
	Imports the interfaces from the tree.
*/
void PatchOctree::import_interfaces()
{
	std::cout << ">> Importing the interfaces\n";

	// Reset
	reset_interfaces();

	// Initialize intersections
	if (is_three_dimensional()) {
		m_tree_3D.computeIntersections();
	} else {
		m_tree_2D.computeIntersections();
	}

	// Importing the interfaces
	int nInterfaces;
	if (is_three_dimensional()) {
		nInterfaces = m_tree_3D.getNumIntersections();
	} else {
		nInterfaces = m_tree_2D.getNumIntersections();
	}

	vector<vector<int>> cellToInterfaceMap(m_cells.size(), vector<int>(0));

	m_interfaces.reserve(nInterfaces);
	for (int n = 0; n < nInterfaces; n++) {
		m_interfaces.emplace_back(n, this);
		Interface &interface = m_interfaces.back();

		// Info associate al tree
		int level;
		bool isGhost;
		bool isBoundary;
		vector<uint32_t> cells;
		int ownerFace;
		int ownerCell;
		if (is_three_dimensional()) {
			Class_Intersection<3> *treeInterface = m_tree_3D.getIntersection(n);

			level      = m_tree_3D.getLevel(treeInterface);
			isGhost    = m_tree_3D.getPbound(treeInterface);
			isBoundary = m_tree_3D.getBound(treeInterface);
			cells      = m_tree_3D.getOwners(treeInterface);
			ownerFace  = m_tree_3D.getFace(treeInterface);
			ownerCell  = m_tree_3D.getFiner(treeInterface) ? 0 : 1;
		} else {
			Class_Intersection<2> *treeInterface = m_tree_2D.getIntersection(n);

			level      = m_tree_2D.getLevel(treeInterface);
			isGhost    = m_tree_2D.getPbound(treeInterface);
			isBoundary = m_tree_2D.getBound(treeInterface);
			cells      = m_tree_2D.getOwners(treeInterface);
			ownerFace  = m_tree_2D.getFace(treeInterface);
			ownerCell  = m_tree_2D.getFiner(treeInterface) ? 0 : 1;
		}

		int ownerId = cells[ownerCell];
		int neighId = cells[(ownerCell + 1) % 2];
		if (isGhost) {
			neighId += m_nInternalCells;
		}

		// Area
		interface.set_area(&m_tree_area[level]);

		// Normal
		int normalIdx;
		switch (ownerFace) {

		case 0:
		case 1:
			normalIdx = Node::COORD_X;
			break;

		case 2:
		case 3:
			normalIdx = Node::COORD_Y;
			break;

		case 4:
		case 5:
			normalIdx = Node::COORD_Z;
			break;

		}

		if (ownerId == 1) {
			normalIdx += get_dimension();
		}

		interface.set_normal(m_normals->get(normalIdx));

		// Position
		if (isGhost) {
			interface.set_position_type(Interface::GHOST);
		} else if (isBoundary) {
			interface.set_position_type(Interface::BOUNDARY);
		} else {
			interface.set_position_type(Interface::INTERNAL);
		}

		// Owner
		Cell &owner = m_cells[ownerId];
		cellToInterfaceMap[ownerId].push_back(interface.get_id());

		interface.set_owner(owner.get_id(), ownerFace);

		// Neighbour
		if (isBoundary) {
			interface.unset_neigh();
		} else {
			Cell &neigh = m_cells[neighId];
			cellToInterfaceMap[neighId].push_back(interface.get_id());

			int neighFace = ownerFace + 1 - 2 * (ownerFace % 2);
			interface.set_neigh(neigh.get_id(), neighFace);
		}
	}

	// Add interface information to the cells
	int nCellFaces;
	if (is_three_dimensional()) {
		nCellFaces = Element::get_face_count(Element::HEXAHEDRON);
	} else {
		nCellFaces = Element::get_face_count(Element::QUADRANGLE);
	}

	vector<vector<int>> interfaceList(nCellFaces, vector<int>(0));
	for (unsigned int n = 0; n < m_cells.size(); n++) {
		for (int k = 0; k < nCellFaces; k++) {
			interfaceList[k].clear();
		}

		Cell &cell = m_cells[n];
		for (unsigned int k = 0; k < cellToInterfaceMap[n].size(); k++) {
			int interfaceId = cellToInterfaceMap[n][k];
			Interface &interface = m_interfaces[interfaceId];

			int face;
			int ownerId = interface.get_owner();
			Cell &owner = m_cells[ownerId];
			if (owner.get_id() == cell.get_id()) {
				face = interface.get_owner_face();
			} else {
				face = interface.get_neigh_face();
			}

			interfaceList[face].push_back(interface.get_id());
		}

		cell.initialize_interfaces(interfaceList);

		cellToInterfaceMap[n].clear();
		vector<int>().swap(cellToInterfaceMap[n]);
	}

	// Clean intersections
	// if (is_three_dimensional()) {
	//	m_tree_3D.clearIntersections();
	//} else {
	//	m_tree_2D.clearIntersections();
	//}
}

/*!
	Updates the interfaces after a tree change.
*/
void PatchOctree::reload_interfaces()
{
	std::cout << ">> Updating the interfaces\n";

	reset_interfaces();
	import_interfaces();

}

/*!
	Marks a cell for refinement.

	\param cell the cell to be refined
*/
void PatchOctree::_mark_for_refinement(Cell &cell)
{
	if (is_three_dimensional()) {
		m_tree_2D.setMarker(cell.get_id(), 1);
	} else {
		m_tree_3D.setMarker(cell.get_id(), 1);
	}
}

}

#endif
