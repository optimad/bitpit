//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "cell.hpp"
#include "interface.hpp"
#include "patch.hpp"

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
Cell::Cell(const int &id)
	: Element(id)
{

}

/*!
	Creates a new cell.
*/
Cell::Cell(const int &id, Patch *patch)
	: Element(id, patch)
{

}

/*!
	Sets the position type of the interface.

	\param position the position type of the interface
*/
void Cell::set_position_type(PositionType positionType)
{
	m_positionType = positionType;
}

/*!
	Gets the position type of the interface.

	\result The position type of the nterface
*/
Cell::PositionType Cell::get_position_type() const
{
  return m_positionType;
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
double & Cell::get_volume() const
{
    return *m_volume;
}  

/*!
	Initialize all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::initialize_interfaces(std::vector<std::vector<int>> &interfaces)
{
	m_interfaces.reset();
	m_interfaces = std::unique_ptr<CollapsedArray2D<int> >(new CollapsedArray2D<int>(interfaces));
}

/*!
	Initialize the data structure that holds the information about the
	interfaces.

	\param nInterfaces An array with the numbero of interfaces for each face
*/
void Cell::initialize_empty_interfaces(const int nInterfaces[])
{
	m_interfaces = std::unique_ptr<CollapsedArray2D<int> >(new CollapsedArray2D<int>(get_face_count(), nInterfaces));
}

/*!
	Sets the i-th interface associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the interface
	\param interface A pointer to the interface
*/
void Cell::set_interface(const int &face, const int &index, const int &interface)
{
	m_interfaces->set(face, index, interface);
}

/*!
	Sets the interfaces associated to the specified face of the cell.

	\param face the face of the cell
	\param interface a pointer to the interfaces
*/
void Cell::set_interfaces(const int &face, int interfaces[])
{
	m_interfaces->set(face, interfaces);
}

/*!
	Unsets the interfaces associated to the cell.
*/
void Cell::unset_interfaces()
{
	m_interfaces.reset();
}

/*!
	Gets the total number of interfaces of the cell.

	\result The total number of interfaces of the cell.
*/
int Cell::get_interface_count() const
{
	return m_interfaces->data_size();
}

/*!
	Gets the number of interfaces of the specified face of the cell.

	\param face the face of the cell
	\result The number of interfaces of the specified face of the cell.
*/
int Cell::get_interface_count(const int &face) const
{
	return m_interfaces->sub_array_size(face);
}

/*!
	Gets the i-th interface of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface.
*/
int Cell::get_interface(const int &face, const int &index) const
{
	return m_interfaces->get(face, index);
}

/*!
	Gets all the interfaces of the cell.

	\result The interfaces of the cell.
*/
const int * Cell::get_interfaces() const
{
	return m_interfaces->get(0);
}

/*!
	Gets the interfaces of the given face of the cell.

	\as get_interface(const int &face, const int &index) const

	\param face the face of the cell
	\result The requested interfaces
*/
const int * Cell::get_interfaces(const int &face) const
{
	return m_interfaces->get(face);
}

/*!
	Extracts the neighbours of all the faces of the cell.

	\result The neighbours of all the faces of the cell.
*/
std::vector<int> Cell::extract_face_neighs() const
{
	std::vector<int> neighs;
	for (int i = 0; i < get_face_count(); ++i) {
		std::vector<int> faceNeighs = extract_face_neighs(i);
		for (auto &faceNeigh : faceNeighs) {
			add_id_to_ordered_list(faceNeigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified face.

	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the face.
*/
std::vector<int> Cell::extract_face_neighs(const int &face, const std::vector<int> &blackList) const
{
	std::vector<int> neighs;
	for (int i = 0; i < get_interface_count(face); ++i) {
		int interfaceId = get_interface(face,i);
		Interface &interface = get_patch()->get_interface(interfaceId);
		if (interface.get_position_type() == Interface::BOUNDARY) {
			continue;
		}

		int neighId = interface.get_neigh();
		if (neighId == get_id()) {
			neighId = interface.get_owner();
		}

		if(std::find(blackList.begin(), blackList.end(), neighId) != blackList.end()) {
				continue;
		}

		// Add the cell to the negihbour list
		add_id_to_ordered_list(neighId, neighs);
	}

	return neighs;
}

/*!
	Extracts the neighbours of all the vertices of the cell.

	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the cell.
*/
std::vector<int> Cell::extract_vertex_neighs(bool complete) const
{
	std::vector<int> blackList;
	if (!complete) {
		if (is_three_dimensional()) {
			blackList = extract_edge_neighs();
		} else {
			blackList = extract_face_neighs();
		}
	}

	std::vector<int> neighs;
	for (int i = 0; i < get_vertex_count(); ++i) {
		for (auto &neigh : extract_vertex_neighs(i, blackList)) {
			add_id_to_ordered_list(neigh, neighs);
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
std::vector<int> Cell::extract_vertex_neighs(const int &vertex, const std::vector<int> &blackList) const
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
std::vector<int> Cell::extract_vertex_neighs(const std::vector<int> &vertices, const std::vector<int> &blackList) const
{
	std::vector<int> neighs;

	int nVerticesToFound = vertices.size();

	std::vector<int> alreadyScanned;
	std::vector<int> processingQueue;
	processingQueue.push_back(get_id());
	while (!processingQueue.empty()) {
		// Get a cell to scan and remove it form the list
		int cellId(processingQueue.back());
		processingQueue.pop_back();
		Cell &cell = get_patch()->get_cell(cellId);

		// Scan the interfaces of the cell
		const int *interfaces = cell.get_interfaces();
		for (int i = 0; i < cell.get_interface_count(); i++) {
			int interfaceId = interfaces[i];
			Interface &interface = get_patch()->get_interface(interfaceId);

			// Neighbour cell assocated to the interface
			//
			// Only consider the cells that are not
			int neighId = interface.get_neigh();
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
			const int *interfaceConnect = interface.get_connect();
			for (int k = 0; k < interface.get_vertex_count(); ++k) {
				for (int n = 0; n < nVerticesToFound; ++n) {
					if (interfaceConnect[k] == get_connect()[vertices[n]]) {
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
					add_id_to_ordered_list(neighId, neighs);
				}
				processingQueue.push_back(neighId);
			}

			// The cell has been scanned
			alreadyScanned.push_back(neighId);
		}
	}

	return neighs;
}

/*!
	Sets the data of the cell.

	\param data a pointer to the data of the cell
*/
void Cell::set_data(std::unique_ptr<CellData> data)
{
	m_data = std::move(data);
}

/*!
	Gets the data of the cell.

	\return A pointer to the data of the cell
*/
CellData * Cell::get_data() const
{
    return m_data.get();
}

}
