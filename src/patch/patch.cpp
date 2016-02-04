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

#include <sstream>
#include <typeinfo>
#include <unordered_map>

#include "patch.hpp"
#include "utils.hpp"

namespace bitpit {

// Explicit instantiation of the bitpit::PiercedVectors
template class PiercedVector<Cell>;
template class PiercedVector<Interface>;
template class PiercedVector<Vertex>;

/*!
	\ingroup patch
	@{
*/

/*!
	\class Patch

	\brief The Patch class provides an interface for defining patches.

	Patch is the base class for defining patches like .
*/

/*!
	Creates a new patch.
*/
Patch::Patch(const int &id, const int &dimension)
	: m_dirty(true)
{
	set_id(id) ;
	set_dimension(dimension);

	std::ostringstream convert;
	convert << get_id();
	set_name(convert.str());
}

/*!
	Destroys the patch.
*/
Patch::~Patch()
{
	reset();
}

/*!
	Updates the mesh

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> Patch::update(bool trackAdaption)
{
	const std::vector<Adaption::Info> adaptionInfo = _update(trackAdaption);

	m_cells.flush();
	m_interfaces.flush();
	m_vertices.flush();

	set_dirty(false);

	return adaptionInfo;
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void Patch::mark_cell_for_refinement(const long &id)
{
	bool updated = _mark_cell_for_refinement(id);

	set_dirty(updated);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void Patch::mark_cell_for_coarsening(const long &id)
{
	bool updated = _mark_cell_for_coarsening(id);

	set_dirty(updated);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void Patch::enable_cell_balancing(const long &id, bool enabled)
{
	bool updated = _enable_cell_balancing(id, enabled);

	set_dirty(updated);
}

/*!
	Resest the patch.
*/
void Patch::reset()
{
	reset_vertices();
	reset_cells();
	reset_interfaces();
}

/*!
	Resest the vertices of the patch.
*/
void Patch::reset_vertices()
{
	m_vertices.clear();
	bitpit::PiercedVector<Vertex>().swap(m_vertices);

	for (auto &cell : m_cells) {
		cell.unset_connect();
	}
}

/*!
	Resest the cells of the patch.
*/
void Patch::reset_cells()
{
	m_cells.clear();
	bitpit::PiercedVector<Cell>().swap(m_cells);

	for (auto &interface : m_interfaces) {
		interface.unset_neigh();
		interface.unset_owner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void Patch::reset_interfaces()
{
	m_interfaces.clear();
	bitpit::PiercedVector<Interface>().swap(m_interfaces);

	for (auto &cell : m_cells) {
		cell.unset_interfaces();
	}
}

/*!
	Writes the mesh to filename specified in input.

	\param filename the filename where the mesh will be written to
*/
void Patch::write_mesh(std::string filename)
{
	bitpit::VTKUnstructuredGrid::setCodex(bitpit::VTKFormat::APPENDED);
	bitpit::VTKUnstructuredGrid::setNames(".", filename);
	bitpit::VTKUnstructuredGrid::write();
}

/*!
	Writes the mesh a filename with the same name of the mesh
*/
void Patch::write_mesh()
{
	write_mesh(get_name());
}

/*!
	Writes a field defined on the patch.

	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::write_field(std::string name, bitpit::VTKLocation location, const std::vector<double> &values)
{
	write_field(get_name(), name, location, values);
}

/*!
	Writes a field defined on the patch using the specified filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::write_field(std::string filename, std::string name, bitpit::VTKLocation location, const std::vector<double> &values)
{
	bitpit::VTKUnstructuredGrid::addData(name, bitpit::VTKFieldType::SCALAR, location);
	m_dataFields[name] = &values;
	m_dataLocations[name] = location;
	m_dataType[name] = bitpit::VTKFieldType::SCALAR;

	write_mesh(filename);

	bitpit::VTKUnstructuredGrid::removeData(name);
	m_dataFields.erase(name);
	m_dataLocations.erase(name);
	m_dataType.erase(name);
}

/*!
	Writes a field defined on the cells of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_cell_field(std::string name, const std::vector<double> &values)
{
	write_cell_field(get_name(), name, values);
}

/*!
	Writes a field defined on the cells of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_cell_field(std::string filename, std::string name, const std::vector<double> &values)
{
	write_field(filename, name, bitpit::VTKLocation::CELL, values);
}

/*!
	Writes a field defined on the vertices of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_vertex_field(std::string name, const std::vector<double> &values)
{
	write_vertex_field(get_name(), name, values);
}

/*!
	Writes a field defined on the vertices of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_vertex_field(std::string filename, std::string name, const std::vector<double> &values)
{
	write_field(filename, name, bitpit::VTKLocation::POINT, values);
}

/*!
	Flags the mesh for update.

	\param dirty if true, then mesh is informed that something in the mesh
	definition has changed and thus the current data structures are not
	valid anymore.
*/
void Patch::set_dirty(bool dirty)
{
	if (m_dirty == dirty) {
		return;
	}

	m_dirty = dirty;
}

/*!
	Returns true if the the mesh needs to update its data strucutres.

	\return This method returns true to indicate the mesh needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool Patch::is_dirty() const
{
	return m_dirty;
}

/*!
	Sets the ID of the patch.

	\param id the ID of the patch
*/
void Patch::set_id(int id)
{
	m_id = id;
}

/*!
	Gets the ID of the patch.

	\return The ID of the patch
*/
int Patch::get_id() const
{
	return m_id;
}

/*!
	Sets the dimension of the patch.

	\param dimension the dimension of the patch
*/
void Patch::set_dimension(int dimension)
{
	m_dimension = dimension;
}

/*!
	Gets the dimension of the patch.

	\return The dimension of the patch
*/
int Patch::get_dimension() const
{
	return m_dimension;
}

/*!
	Returns true if the mesh is a three-dimensional mesh.

	\return This method returns true to indicate the mesh is three-dimensional
*/
bool Patch::is_three_dimensional() const
{
	return (m_dimension == 3);
}

/*!
	Sets the name of the patch.

	\param id the name of the patch
*/
void Patch::set_name(std::string name)
{
	m_name = name;
}

/*!
	Gets the name of the patch.

	\return The name of the patch
*/
std::string Patch::get_name() const
{
	return m_name;
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long Patch::get_vertex_count() const
{
	return m_vertices.size();
}

/*!
	Gets the nodes owned by the patch.

	\return The nodes owned by the patch.
*/
bitpit::PiercedVector<Vertex> & Patch::vertices()
{
	return m_vertices;
}

/*!
	Gets a reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A reference to the vertex with the specified id.
*/
Vertex & Patch::get_vertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Gets a constant reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A constant reference to the vertex with the specified id.
*/
const Vertex & Patch::get_vertex(const long &id) const
{
	return m_vertices[id];
}

/*!
	Creates a new vertex with the specified id.

	\param id is the id of the new vertex
*/
long Patch::create_vertex(const long &id)
{
	m_vertices.reclaim(id);

	return id;
}

/*!
	Creates a new vertex.
*/
long Patch::create_vertex()
{
	long id;
	if (m_unusedVertexIds.empty()) {
		id = m_vertices.size();
	} else {
		id = m_unusedVertexIds.front();
		m_unusedVertexIds.pop_front();
	}

	return create_vertex(id);
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
*/
void Patch::delete_vertex(const long &id, bool delayed)
{
	m_vertices.erase(id, delayed);
	m_unusedVertexIds.push_back(id);
}

/*!
	Gets the coordinates of the specified vertex.

	\param is is the id of the vertex
	\result The coordinates of the specified vertex.
*/
const std::array<double, 3> & Patch::get_vertex_coords(const long &id) const
{
	return get_vertex(id).get_coords();
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long Patch::get_cell_count() const
{
	return m_cells.size();
}

/*!
	Gets the cells owned by the patch.

	\return The cells owned by the patch.
*/
bitpit::PiercedVector<Cell> & Patch::cells()
{
	return m_cells;
}

/*!
	Gets a reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A reference to the cell with the specified id.
*/
Cell & Patch::get_cell(const long &id)
{
	return m_cells[id];
}

/*!
	Gets a constant reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A constant reference to the cell with the specified id.
*/
const Cell & Patch::get_cell(const long &id) const
{
	return m_cells[id];
}

/*!
	Creates a new cell with the specified id.

	\param id is the id of the new cell
	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::create_cell(const long &id, bool internal, ElementInfo::Type type)
{
	bitpit::PiercedVector<Cell>::iterator iterator;
	if (internal) {
		iterator = m_cells.reclaim(id);
	} else {
		iterator = m_cells.reclaim_back(id);
	}

	Cell &cell = *iterator;
	cell.initialize(type);

	return id;
}

/*!
	Creates a new cell.

	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::create_cell(bool internal, ElementInfo::Type type)
{
	long id;
	if (m_unusedCellIds.empty()) {
		id = m_cells.size();
	} else {
		id = m_unusedCellIds.front();
		m_unusedCellIds.pop_front();
	}

	return create_cell(id, internal, type);
}

/*!
	Deletes a cell.

	\param id is the id of the cell
*/
void Patch::delete_cell(const long &id, bool delayed)
{
	m_cells.erase(id, delayed);
	m_unusedCellIds.push_back(id);
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\result The neighbours of all the faces of the specified cell.
*/
std::vector<long> Patch::extract_cell_face_neighs(const long &id) const
{
	std::vector<long> neighs;
	const Cell &cell = get_cell(id);
	for (int i = 0; i < cell.get_face_count(); ++i) {
		std::vector<long> faceNeighs = extract_cell_face_neighs(id, i);
		for (auto &neighId : faceNeighs) {
			bitpit::utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\result All the neighbours of the specified cell.
*/
std::vector<long> Patch::extract_cell_neighs(const long &id) const
{
	return extract_cell_vertex_neighs(id);
}

/*!
	Extracts all the neighbours of the specified cell for the given
	codimension.

	\param id is the id of the cell
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
std::vector<long> Patch::extract_cell_neighs(const long &id, int codimension, bool complete) const
{
	assert(codimension >= 1 && codimension <= get_dimension());

	if (codimension == 1) {
		return extract_cell_face_neighs(id);
	} else if (codimension == get_dimension()) {
		return extract_cell_vertex_neighs(id, complete);
	} else if (codimension == 2) {
		return extract_cell_edge_neighs(id, complete);
	} else {
		return std::vector<long>();
	}
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> Patch::extract_cell_face_neighs(const long &id, const int &face, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	const Cell &cell = get_cell(id);
	for (int i = 0; i < cell.get_interface_count(face); ++i) {
		long interfaceId = cell.get_interface(face, i);
		const Interface &interface = get_interface(interfaceId);
		if (interface.is_border()) {
			continue;
		}

		long neighId = interface.get_neigh();
		if (neighId == cell.get_id()) {
			neighId = interface.get_owner();
		}

		if (std::find(blackList.begin(), blackList.end(), neighId) != blackList.end()) {
			continue;
		}

		// Add the cell to the negihbour list
		bitpit::utils::addToOrderedVector<long>(neighId, neighs);
	}

	return neighs;
}

/*!
	Extracts the neighbours of all the edges of the specified cell.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified edge, or should
	contain also neighbours that share an entire face
	\result The neighbours of all the edges of the specified cell.
*/
std::vector<long> Patch::extract_cell_edge_neighs(const long &id, bool complete) const
{
	assert(is_three_dimensional());
	if (!is_three_dimensional()) {
		return std::vector<long>();
	}

	std::vector<long> blackList;
	if (!complete) {
		blackList = extract_cell_face_neighs(id);
	}

	std::vector<long> neighs;
	const Cell &cell = get_cell(id);
	for (int i = 0; i < cell.get_edge_count(); ++i) {
		for (auto &neigh : extract_cell_edge_neighs(id, i, blackList)) {
			bitpit::utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param vertex is an edge of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> Patch::extract_cell_edge_neighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	assert(is_three_dimensional());
	if (!is_three_dimensional()) {
		return std::vector<long>();
	}

	const Cell &cell = get_cell(id);
	std::vector<int> vertices = cell.get_edge_local_connect(edge);

	return extract_cell_vertex_neighs(id, vertices, blackList);
}

/*!
	Extracts the neighbours of all the vertices of the specified cell.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the specified cell.
*/
std::vector<long> Patch::extract_cell_vertex_neighs(const long &id, bool complete) const
{
	std::vector<long> blackList;
	if (!complete) {
		if (is_three_dimensional()) {
			blackList = extract_cell_edge_neighs(id);
		} else {
			blackList = extract_cell_face_neighs(id);
		}
	}

	std::vector<long> neighs;
	const Cell &cell = get_cell(id);
	for (int i = 0; i < cell.get_vertex_count(); ++i) {
		for (auto &neigh : extract_cell_vertex_neighs(id, i, blackList)) {
			bitpit::utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given vertex.

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

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> Patch::extract_cell_vertex_neighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	std::vector<int> vertexList(1);
	vertexList[0] = vertex;

	return extract_cell_vertex_neighs(id, vertexList, blackList);
}

/*!
	Extracts the neighbours of the specified cell for the given vertices.

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

	\param id is the id of the cell
	\param vertices is the list of vertices of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertices.
*/
std::vector<long> Patch::extract_cell_vertex_neighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;

	int nVerticesToFound = vertices.size();

	const Cell &cell = get_cell(id);
	const long *cellConnect = cell.get_connect();

	std::vector<long> alreadyScanned;
	std::vector<long> processingQueue;
	processingQueue.push_back(cell.get_id());
	while (!processingQueue.empty()) {
		// Get a cell to scan and remove it form the list
		long scanId(processingQueue.back());
		processingQueue.pop_back();
		const Cell &scanCell = get_cell(scanId);

		// Scan the interfaces of the cell
		const long *interfaces = scanCell.get_interfaces();
		for (int i = 0; i < scanCell.get_interface_count(); i++) {
			long interfaceId = interfaces[i];
			const Interface &interface = get_interface(interfaceId);

			// Neighbour cell assocated to the interface
			//
			// Only consider the cells that are not
			long neighId = interface.get_neigh();
			if (neighId < 0 || neighId == scanId) {
				neighId = interface.get_owner();
			}

			if (neighId == id) {
				continue;
			} else if (std::find(alreadyScanned.begin(), alreadyScanned.end(), neighId) != alreadyScanned.end()) {
				continue;
			}

			// Number of vertices owned by the interface
			int nCommonVertices = 0;
			const long *interfaceConnect = interface.get_connect();
			for (int k = 0; k < interface.get_vertex_count(); ++k) {
				for (int n = 0; n < nVerticesToFound; ++n) {
					if (interfaceConnect[k] == cellConnect[vertices[n]]) {
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
				if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
					bitpit::utils::addToOrderedVector<long>(neighId, neighs);
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
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
long Patch::get_interface_count() const
{
	return m_interfaces.size();
}

/*!
	Gets the interfaces owned by the patch.

	\return The interfaces owned by the patch.
*/
bitpit::PiercedVector<Interface> & Patch::interfaces()
{
	return m_interfaces;
}

/*!
	Gets a reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A reference to the interface with the specified id.
*/
Interface & Patch::get_interface(const long &id)
{
	return m_interfaces[id];
}

/*!
	Gets a constant reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A constant reference to the interface with the specified id.
*/
const Interface & Patch::get_interface(const long &id) const
{
	return m_interfaces[id];
}

/*!
	Creates a new interface with the specified id.

	\param id is the id of the new interface
*/
long Patch::create_interface(const long &id, ElementInfo::Type type)
{
	bitpit::PiercedVector<Interface>::iterator iterator = m_interfaces.reclaim(id);

	Interface &interface = *iterator;
	interface.initialize(type);

	return id;
}

/*!
	Creates a new interface.
*/
long Patch::create_interface(ElementInfo::Type type)
{
	long id;
	if (m_unusedInterfaceIds.empty()) {
		id = m_interfaces.size();
	} else {
		id = m_unusedInterfaceIds.front();
		m_unusedInterfaceIds.pop_front();
	}

	return create_interface(id, type);
}

/*!
	Deletes an interface.

	\param id is the id of the interface
*/
void Patch::delete_interface(const long &id, bool delayed)
{
	m_interfaces.erase(id, delayed);
	m_unusedInterfaceIds.push_back(id);
}

/*!
	Sorts the internal storage for cells, vertices and interfaces in
	ascending id order.
*/
void Patch::sort()
{
	m_vertices.sort();
	m_cells.sort();
	m_interfaces.sort();
}

/*!
	Requests the patch to compact the data structures and reduce its capacity
	to fit its size.

	The request is non-binding, and after the function call the patch can
	still occupy more memory than it actually needs.
*/
void Patch::squeeze()
{
	m_vertices.squeeze();
	m_cells.squeeze();
	m_interfaces.squeeze();
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> Patch::eval_cell_centroid(const long &id)
{
	Cell &cell = get_cell(id);

	return eval_element_centroid(cell);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> Patch::eval_interface_centroid(const long &id)
{
	Interface &interface = get_interface(id);

	return eval_element_centroid(interface);
}

/*!
	Evaluates the centroid of the specified element.

	\param element is the element
	\result The centroid of the specified element.
*/
std::array<double, 3> Patch::eval_element_centroid(const Element &element)
{
	const int nDimensions = 3;

	const long *elementConnect = element.get_connect();
	const ElementInfo &elementInfo = element.get_info();

	std::array<double, nDimensions> centroid = {{0., 0., 0.}};
	for (int i = 0; i < elementInfo.nVertices; ++i) {
		Vertex &vertex = get_vertex(elementConnect[i]);
		const std::array<double, nDimensions> &vertexCoords = vertex.get_coords();
		for (int k = 0; k < nDimensions; ++k) {
			centroid[k] += vertexCoords[k];
		}
	}

	for (int k = 0; k < nDimensions; ++k) {
		centroid[k] /= elementInfo.nVertices;
	}

	return centroid;
}

/*!
 *  Interface method for obtaining field meta Data
 *
 *  @param[in] name is the name of the field to be written
 *  @return Returns a bitpit::VTKFieldMetaData struct containing the metadata
 *  of the requested custom data.
 */
const bitpit::VTKFieldMetaData Patch::getMetaData(std::string name)
{
	if (name == "Points") {
		std::cout << "Numero di punti: " << 3 * m_vertices.size() << std::endl;
		return bitpit::VTKFieldMetaData(3 * m_vertices.size(), typeid(double));
	} else if (name == "offsets") {
		std::cout << "Offset size: " << m_cells.size() << std::endl;
		return bitpit::VTKFieldMetaData(m_cells.size(), typeid(int));
	} else if (name == "types") {
		std::cout << "Type size: " << m_cells.size() << std::endl;
		return bitpit::VTKFieldMetaData(m_cells.size(), typeid(bitpit::VTKElementType));
	} else if (name == "connectivity") {
		long connectSize = 0;
		for (Cell &cell : m_cells) {
			connectSize += cell.get_info().nVertices;
		}
		std::cout << "Connect size: " << connectSize << std::endl;
		return bitpit::VTKFieldMetaData(connectSize, typeid(long));
	} else if (m_dataFields.count(name) > 0) {
		long fieldSize = 0;

		if (m_dataLocations[name] == bitpit::VTKLocation::CELL) {
			fieldSize = m_cells.size();
		} else {
			fieldSize = m_vertices.size();
		}

		if (m_dataType[name] == bitpit::VTKFieldType::VECTOR) {
			fieldSize *= 3;
		}

		std::cout << "Field size: " << fieldSize << std::endl;

		return bitpit::VTKFieldMetaData(fieldSize, typeid(double));
	}
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] codex is the codex which must be used. Supported options
 *  are "ascii" or "appended". For "appended" type an unformatted binary
 *  stream must be used
 *  @param[in] name is the name of the data to be written. Either user
 *  data or grid data
 */
void Patch::flushData(std::fstream &stream, bitpit::VTKFormat format, std::string name)
{
	assert(format == bitpit::VTKFormat::APPENDED);

	static std::unordered_map<long, long> vertexMap;

	if (name == "Points") {
		long vertexId = 0;
		for (Vertex &vertex : m_vertices) {
			vertexMap[vertex.get_id()] = vertexId++;

			bitpit::genericIO::flushBINARY(stream, vertex.get_coords());
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (Cell &cell : m_cells) {
			offset += cell.get_info().nVertices;
			bitpit::genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "types") {
		for (Cell &cell : m_cells) {
			bitpit::VTKElementType VTKType;
			switch (cell.get_type())  {

			case ElementInfo::VERTEX:
				VTKType = bitpit::VTKElementType::VERTEX;
				break;

			case ElementInfo::LINE:
				VTKType = bitpit::VTKElementType::LINE;
				break;

			case ElementInfo::TRIANGLE:
				VTKType = bitpit::VTKElementType::TRIANGLE;
				break;

			case ElementInfo::PIXEL:
				VTKType = bitpit::VTKElementType::PIXEL;
				break;

			case ElementInfo::QUAD:
				VTKType = bitpit::VTKElementType::QUAD;
				break;

			case ElementInfo::TETRA:
				VTKType = bitpit::VTKElementType::TETRA;
				break;

			case ElementInfo::VOXEL:
				VTKType = bitpit::VTKElementType::VOXEL;
				break;

			case ElementInfo::HEXAHEDRON:
				VTKType = bitpit::VTKElementType::HEXAHEDRON;
				break;

			case ElementInfo::WEDGE:
				VTKType = bitpit::VTKElementType::WEDGE;
				break;

			case ElementInfo::PYRAMID:
				VTKType = bitpit::VTKElementType::PYRAMID;
				break;

			default:
				VTKType = bitpit::VTKElementType::UNDEFINED;
				break;

			}

			bitpit::genericIO::flushBINARY(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (Cell &cell : m_cells) {
			for (int i = 0; i < cell.get_info().nVertices; ++i) {
				bitpit::genericIO::flushBINARY(stream, vertexMap.at(cell.get_vertex(i)));
			}
		}

		vertexMap.clear();
		std::unordered_map<long, long>().swap(vertexMap);
	} else if (m_dataFields.count(name) > 0) {
		bitpit::genericIO::flushBINARY(stream, *(m_dataFields.at(name)));
	}
}

/*!
	@}
*/

}
