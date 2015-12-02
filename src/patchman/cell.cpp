//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "cell.hpp"
#include "interface.hpp"
#include "patch.hpp"
#include "utils.hpp"

#include<iostream>

namespace pman {

/*!
	\class Cell

	\brief The Cell class defines the cells.

	Cell is class that defines the cells.
*/

/*!
	Default constructor.
*/
Cell::Cell()
	: Element()
{

}

/*!
	Creates a new cell.
*/
Cell::Cell(const long &id, ElementInfo::Type type, Patch *patch)
	: Element(id, type, patch)
{

}

/*!
	Sets if the cells belongs to the the interior domain.

	\param interior defines if the cells belongs to the the interior domain
*/
void Cell::set_interior(bool interior)
{
	m_interior = interior;
}

/*!
	Gets if the cells belongs to the the interior domain.

	\result Returns true if the cells belongs to the the interior domain,
	otherwise it returns false.
*/
bool Cell::is_interior() const
{
	return m_interior;
}

/*!
	Sets the volume of the cell.

	\param volume the volume of the cell
*/
void Cell::set_volume(double * const volume)
{
	m_volume = volume;
}

/*!
	Gets the volume of the cell.

	\return The volume of the cell
*/
const double & Cell::get_volume() const
{
    return *m_volume;
}  

/*!
	Initialize all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::initialize_interfaces(std::vector<std::vector<long>> &interfaces)
{
	m_interfaces = CollapsedVector2D<long>(interfaces);
}

/*!
	Initialize the data structure that holds the information about the
	interfaces.

	\param nInterfaces An array with the numbero of interfaces for each face
*/
void Cell::initialize_empty_interfaces(const int nInterfaces[])
{
	for (int k = 0; k < get_face_count(); ++k) {
		m_interfaces.push_back(nInterfaces[k], 0);
	}
}

/*!
	Sets the i-th interface associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the interface
	\param interface is the index of the interface 
*/
void Cell::set_interface(const int &face, const int &index, const long &interface)
{
	m_interfaces.set(face, index, interface);
}

/*!
	Add an interface to the given face of the cell.

	\param face is the face of the cell
	\param interface is the index of the interface that will be added
*/
void Cell::push_interface(const int &face, const long &interface)
{
	m_interfaces.push_back_in_sub_array(face, interface);
}

/*!
	Deletes the specified interface from the interfaces associate to the
	given face of the cell.

	\param face the face of the cell
	\param i is the index of the interface to delete
*/
void Cell::delete_interface(const int &face, const int &i)
{
	m_interfaces.erase(face, i);
}

/*!
	Unsets the interfaces associated to the cell.
*/
void Cell::unset_interfaces()
{
	m_interfaces.clear();
}

/*!
	Gets the total number of interfaces of the cell.

	\result The total number of interfaces of the cell.
*/
int Cell::get_interface_count() const
{
	return m_interfaces.sub_arrays_total_size();
}

/*!
	Gets the number of interfaces of the specified face of the cell.

	\param face the face of the cell
	\result The number of interfaces of the specified face of the cell.
*/
int Cell::get_interface_count(const int &face) const
{
	return m_interfaces.sub_array_size(face);
}

/*!
	Gets the i-th interface of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface.
*/
long Cell::get_interface(const int &face, const int &index) const
{
	return m_interfaces.get(face, index);
}

/*!
	Gets all the interfaces of the cell.

	\result The interfaces of the cell.
*/
const long * Cell::get_interfaces() const
{
	return m_interfaces.get(0);
}

/*!
	Gets the interfaces of the given face of the cell.

	\as get_interface(const int &face, const int &index) const

	\param face the face of the cell
	\result The requested interfaces
*/
const long * Cell::get_interfaces(const int &face) const
{
	return m_interfaces.get(face);
}

/*!
	Extracts the neighbours of all the faces of the cell.

	\result The neighbours of all the faces of the cell.
*/
std::vector<long> Cell::extract_face_neighs() const
{
	std::vector<long> neighs;
	for (int i = 0; i < get_face_count(); ++i) {
		std::vector<long> faceNeighs = extract_face_neighs(i);
		for (auto &faceNeigh : faceNeighs) {
			utils::add_to_ordered_vector<long>(faceNeigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts all the neighbours of the cell

	\result All the neighbours of the cell.
*/
std::vector<long> Cell::extract_neighs() const
{
	return extract_vertex_neighs();
}

/*!
	Extracts all the neighbours for the specified codimension.

	\param codimension the codimension for which the neighbours
	are requested. For a three-dimensional cell a codimension
	equal 1 will extract the face neighbours, a codimension equal
	2 will extract the edge negihbours and a codimension equal
	3 will extract the vertex neighbours. For a two-dimensional
	cell a codimension qual 1 will extract the face neighbours,
	and a codimension equal 2 will extract the vertex neighbours.
	\param complete controls if the list of neighbours should contain
	only the neighbours for the specified codimension, or should contain
	also the neighbours for lower codimensions.
	\result The neighbours for the specified codimension.
*/
std::vector<long> Cell::extract_neighs(int codimension, bool complete) const
{
	assert(codimension >= 1 && codimension <= get_dimension());

	if (codimension == 1) {
		return extract_face_neighs();
	} else if (codimension == get_dimension()) {
		return extract_vertex_neighs(complete);
	} else if (codimension == 2) {
		return extract_edge_neighs(complete);
	} else {
		return std::vector<long>();
	}
}

/*!
	Extracts the neighbours of the specified face.

	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the face.
*/
std::vector<long> Cell::extract_face_neighs(const int &face, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	for (int i = 0; i < get_interface_count(face); ++i) {
		long interfaceId = get_interface(face,i);
		Interface &interface = get_patch()->get_interface(interfaceId);
		if (interface.get_position_type() == Interface::BOUNDARY) {
			continue;
		}

		long neighId = interface.get_neigh();
		if (neighId == get_id()) {
			neighId = interface.get_owner();
		}

		if(std::find(blackList.begin(), blackList.end(), neighId) != blackList.end()) {
				continue;
		}

		// Add the cell to the negihbour list
		utils::add_to_ordered_vector<long>(neighId, neighs);
	}

	return neighs;
}

/*!
	Extracts the neighbours of all the edges of the cell.

	This function can be only used with three-dimensional cells.

	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified edge, or should
	contain also neighbours that share an entire face
	\result The neighbours of all the edges of the cell.
*/
std::vector<long> Cell::extract_edge_neighs(bool complete) const
{
	assert(is_three_dimensional());
	if (!is_three_dimensional()) {
		return std::vector<long>();
	}

	std::vector<long> blackList;
	if (!complete) {
		blackList = extract_face_neighs();
	}

	std::vector<long> neighs;
	for (int i = 0; i < get_edge_count(); ++i) {
		for (auto &neigh : extract_edge_neighs(i, blackList)) {
			utils::add_to_ordered_vector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified edge.

	This function can be only used with three-dimensional cells.

	\param vertex is an edge of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the edge.
*/
std::vector<long> Cell::extract_edge_neighs(const int &edge, const std::vector<long> &blackList) const
{
	assert(is_three_dimensional());
	if (!is_three_dimensional()) {
		return std::vector<long>();
	}

	std::vector<int> vertices;
	if (is_three_dimensional()) {
		vertices = get_edge_local_connect(edge);
	} else {
		vertices.push_back(edge);
	}

	return Cell::extract_vertex_neighs(vertices, blackList);
}

/*!
	Extracts the neighbours of all the vertices of the cell.

	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the cell.
*/
std::vector<long> Cell::extract_vertex_neighs(bool complete) const
{
	std::vector<long> blackList;
	if (!complete) {
		if (is_three_dimensional()) {
			blackList = extract_edge_neighs();
		} else {
			blackList = extract_face_neighs();
		}
	}

	std::vector<long> neighs;
	for (int i = 0; i < get_vertex_count(); ++i) {
		for (auto &neigh : extract_vertex_neighs(i, blackList)) {
			utils::add_to_ordered_vector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified vertex.

	Cells that has only a vertex in common are considered neighbours only
	if there are other cells "connecting" them.

	                  .-----.                   .-----.
	                  |     |                   |     |
	                V | A1  |                 V | A2  |
	            .-----+-----.             .-----+-----.
	            |     |                   |     |     |
	            | B1  |                   | B2  | C2  |
	            .-----.                   .-----.-----.

	For example, A1 and B1 are not neighbours (although they share the
	vertex V), whereas A2 and B2 are neighbours.

	\param vertex is a vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the vertex.
*/
std::vector<long> Cell::extract_vertex_neighs(const int &vertex, const std::vector<long> &blackList) const
{
	std::vector<int> vertices(1);
	vertices[0] = vertex;

	return Cell::extract_vertex_neighs(vertices, blackList);
}

/*!
	Extracts the neighbours that share the specified vertices.

	Cells that has only a vertex in common are considered neighbours only
	if there are other cells "connecting" them.

	                  .-----.                   .-----.
	                  |     |                   |     |
	                V | A1  |                 V | A2  |
	            .-----+-----.             .-----+-----.
	            |     |                   |     |     |
	            | B1  |                   | B2  | C2  |
	            .-----.                   .-----.-----.

	For example, A1 and B1 are not neighbours (although they share the
	vertex V), whereas A2 and B2 are neighbours.

	\param vertices is the list of vertices of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours that share the specified vertices.
*/
std::vector<long> Cell::extract_vertex_neighs(const std::vector<int> &vertices, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;

	int nVerticesToFound = vertices.size();

	const long *connectivity = get_connect();

	std::vector<long> alreadyScanned;
	std::vector<long> processingQueue;
	processingQueue.push_back(get_id());
	while (!processingQueue.empty()) {
		// Get a cell to scan and remove it form the list
		long cellId(processingQueue.back());
		processingQueue.pop_back();
		Cell &cell = get_patch()->get_cell(cellId);

		// Scan the interfaces of the cell
		const long *interfaces = cell.get_interfaces();
		for (int i = 0; i < cell.get_interface_count(); i++) {
			long interfaceId = interfaces[i];
			Interface &interface = get_patch()->get_interface(interfaceId);

			// Neighbour cell assocated to the interface
			//
			// Only consider the cells that are not
			long neighId = interface.get_neigh();
			if (neighId < 0 || neighId == cell.get_id()) {
				neighId = interface.get_owner();
			}

			if (neighId == get_id()) {
				continue;
			} else if(std::find(alreadyScanned.begin(), alreadyScanned.end(), neighId) != alreadyScanned.end()) {
				continue;
			}

			// Number of vertices owned by the interface
			int nCommonVertices = 0;
			const long *interfaceConnect = interface.get_connect();
			for (int k = 0; k < interface.get_vertex_count(); ++k) {
				for (int n = 0; n < nVerticesToFound; ++n) {
					if (interfaceConnect[k] == connectivity[vertices[n]]) {
						nCommonVertices++;
						break;
					}
				}

				if (nCommonVertices == nVerticesToFound) {
					break;
				}
			}

			// If the interface contains all the requested vertices,
			// add the neighbour cell of the interface to the list
			// of cells neighbours.
			if (nCommonVertices == nVerticesToFound) {
				if(std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
					utils::add_to_ordered_vector<long>(neighId, neighs);
				}
				processingQueue.push_back(neighId);
			}

			// The cell has been scanned
			alreadyScanned.push_back(neighId);
		}
	}

	return neighs;
}

}
