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
	: m_dirty(true), m_hasCustomTolerance(false)
{
	set_id(id) ;
	setDimension(dimension);

	std::ostringstream convert;
	convert << get_id();
	setName(convert.str());
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

	updateBoundingBox();

	setDirty(false);

	return adaptionInfo;
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void Patch::markCellForRefinement(const long &id)
{
	bool updated = _markCellForRefinement(id);

	setDirty(updated);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void Patch::markCellForCoarsening(const long &id)
{
	bool updated = _markCellForCoarsening(id);

	setDirty(updated);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void Patch::enableCellBalancing(const long &id, bool enabled)
{
	bool updated = _enableCellBalancing(id, enabled);

	setDirty(updated);
}

/*!
	Resest the patch.
*/
void Patch::reset()
{
	resetVertices();
	resetCells();
	resetInterfaces();
}

/*!
	Resest the vertices of the patch.
*/
void Patch::resetVertices()
{
	m_vertices.clear();
	PiercedVector<Vertex>().swap(m_vertices);

	for (auto &cell : m_cells) {
		cell.unsetConnect();
	}
}

/*!
	Resest the cells of the patch.
*/
void Patch::resetCells()
{
	m_cells.clear();
	PiercedVector<Cell>().swap(m_cells);

	for (auto &interface : m_interfaces) {
		interface.unsetNeigh();
		interface.unsetOwner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void Patch::resetInterfaces()
{
	m_interfaces.clear();
	PiercedVector<Interface>().swap(m_interfaces);

	for (auto &cell : m_cells) {
		cell.unsetInterfaces();
	}
}

/*!
	Writes the mesh to filename specified in input.

	\param filename the filename where the mesh will be written to
*/
void Patch::writeMesh(std::string filename)
{
	VTKUnstructuredGrid::setCodex(VTKFormat::APPENDED);
	VTKUnstructuredGrid::setNames(".", filename);
	VTKUnstructuredGrid::write();
}

/*!
	Writes the mesh a filename with the same name of the mesh
*/
void Patch::writeMesh()
{
	writeMesh(getName());
}

/*!
	Writes a field defined on the patch.

	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::writeField(std::string name, VTKLocation location, const std::vector<double> &values)
{
	writeField(getName(), name, location, values);
}

/*!
	Writes a field defined on the patch using the specified filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::writeField(std::string filename, std::string name, VTKLocation location, const std::vector<double> &values)
{
	VTKUnstructuredGrid::addData(name, VTKFieldType::SCALAR, location);
	m_dataFields[name] = &values;
	m_dataLocations[name] = location;
	m_dataType[name] = VTKFieldType::SCALAR;

	writeMesh(filename);

	VTKUnstructuredGrid::removeData(name);
	m_dataFields.erase(name);
	m_dataLocations.erase(name);
	m_dataType.erase(name);
}

/*!
	Writes a field defined on the cells of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeCellField(std::string name, const std::vector<double> &values)
{
	writeCellField(getName(), name, values);
}

/*!
	Writes a field defined on the cells of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeCellField(std::string filename, std::string name, const std::vector<double> &values)
{
	writeField(filename, name, VTKLocation::CELL, values);
}

/*!
	Writes a field defined on the vertices of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeVertexField(std::string name, const std::vector<double> &values)
{
	writeVertexField(getName(), name, values);
}

/*!
	Writes a field defined on the vertices of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeVertexField(std::string filename, std::string name, const std::vector<double> &values)
{
	writeField(filename, name, VTKLocation::POINT, values);
}

/*!
	Flags the mesh for update.

	\param dirty if true, then mesh is informed that something in the mesh
	definition has changed and thus the current data structures are not
	valid anymore.
*/
void Patch::setDirty(bool dirty)
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
bool Patch::isDirty() const
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
void Patch::setDimension(int dimension)
{
	m_dimension = dimension;
}

/*!
	Gets the dimension of the patch.

	\return The dimension of the patch
*/
int Patch::getDimension() const
{
	return m_dimension;
}

/*!
	Returns true if the mesh is a three-dimensional mesh.

	\return This method returns true to indicate the mesh is three-dimensional
*/
bool Patch::isThreeDimensional() const
{
	return (m_dimension == 3);
}

/*!
	Sets the name of the patch.

	\param id the name of the patch
*/
void Patch::setName(std::string name)
{
	m_name = name;
}

/*!
	Gets the name of the patch.

	\return The name of the patch
*/
std::string Patch::getName() const
{
	return m_name;
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long Patch::getVertexCount() const
{
	return m_vertices.size();
}

/*!
	Gets the nodes owned by the patch.

	\return The nodes owned by the patch.
*/
PiercedVector<Vertex> & Patch::vertices()
{
	return m_vertices;
}

/*!
	Gets a reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A reference to the vertex with the specified id.
*/
Vertex & Patch::getVertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Gets a constant reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A constant reference to the vertex with the specified id.
*/
const Vertex & Patch::getVertex(const long &id) const
{
	return m_vertices[id];
}

/*!
	Adds a new vertex with the specified id.

	\param id is the id of the new vertex
*/
long Patch::addVertex(const long &id)
{
	m_vertices.reclaim(id);

	return id;
}

/*!
	Adds a new vertex.
*/
long Patch::addVertex()
{
	long id;
	if (m_unusedVertexIds.empty()) {
		id = m_vertices.size();
	} else {
		id = m_unusedVertexIds.front();
		m_unusedVertexIds.pop_front();
	}

	return addVertex(id);
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
*/
void Patch::deleteVertex(const long &id, bool delayed)
{
	m_vertices.erase(id, delayed);
	m_unusedVertexIds.push_back(id);
}

/*!
	Gets the coordinates of the specified vertex.

	\param is is the id of the vertex
	\result The coordinates of the specified vertex.
*/
const std::array<double, 3> & Patch::getVertexCoords(const long &id) const
{
	return getVertex(id).getCoords();
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long Patch::getCellCount() const
{
	return m_cells.size();
}

/*!
	Gets the cells owned by the patch.

	\return The cells owned by the patch.
*/
PiercedVector<Cell> & Patch::cells()
{
	return m_cells;
}

/*!
	Gets a reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A reference to the cell with the specified id.
*/
Cell & Patch::getCell(const long &id)
{
	return m_cells[id];
}

/*!
	Gets a constant reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A constant reference to the cell with the specified id.
*/
const Cell & Patch::getCell(const long &id) const
{
	return m_cells[id];
}

/*!
	Adds a new cell with the specified id.

	\param id is the id of the new cell
	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::addCell(const long &id, bool internal, ElementInfo::Type type)
{
	PiercedVector<Cell>::iterator iterator;
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
	Adds a new cell.

	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::addCell(bool internal, ElementInfo::Type type)
{
	long id;
	if (m_unusedCellIds.empty()) {
		id = m_cells.size();
	} else {
		id = m_unusedCellIds.front();
		m_unusedCellIds.pop_front();
	}

	return addCell(id, internal, type);
}

/*!
	Deletes a cell.

	\param id is the id of the cell
*/
void Patch::deleteCell(const long &id, bool delayed)
{
	m_cells.erase(id, delayed);
	m_unusedCellIds.push_back(id);
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\result The neighbours of all the faces of the specified cell.
*/
std::vector<long> Patch::extractCellFaceNeighs(const long &id) const
{
	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getFaceCount(); ++i) {
		std::vector<long> faceNeighs = extractCellFaceNeighs(id, i);
		for (auto &neighId : faceNeighs) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\result All the neighbours of the specified cell.
*/
std::vector<long> Patch::extractCellNeighs(const long &id) const
{
	return extractCellVertexNeighs(id);
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
std::vector<long> Patch::extractCellNeighs(const long &id, int codimension, bool complete) const
{
	assert(codimension >= 1 && codimension <= getDimension());

	if (codimension == 1) {
		return extractCellFaceNeighs(id);
	} else if (codimension == getDimension()) {
		return extractCellVertexNeighs(id, complete);
	} else if (codimension == 2) {
		return extractCellEdgeNeighs(id, complete);
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
std::vector<long> Patch::extractCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getInterfaceCount(face); ++i) {
		long interfaceId = cell.getInterface(face, i);
		const Interface &interface = getInterface(interfaceId);
		if (interface.isBorder()) {
			continue;
		}

		long neighId = interface.getNeigh();
		if (neighId == cell.get_id()) {
			neighId = interface.getOwner();
		}

		if (std::find(blackList.begin(), blackList.end(), neighId) != blackList.end()) {
			continue;
		}

		// Add the cell to the negihbour list
		utils::addToOrderedVector<long>(neighId, neighs);
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
std::vector<long> Patch::extractCellEdgeNeighs(const long &id, bool complete) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return std::vector<long>();
	}

	std::vector<long> blackList;
	if (!complete) {
		blackList = extractCellFaceNeighs(id);
	}

	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getEdgeCount(); ++i) {
		for (auto &neigh : extractCellEdgeNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
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
std::vector<long> Patch::extractCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return std::vector<long>();
	}

	const Cell &cell = getCell(id);
	std::vector<int> vertices = cell.getEdgeLocalConnect(edge);

	return extractCellVertexNeighs(id, vertices, blackList);
}

/*!
	Extracts the neighbours of all the vertices of the specified cell.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the specified cell.
*/
std::vector<long> Patch::extractCellVertexNeighs(const long &id, bool complete) const
{
	std::vector<long> blackList;
	if (!complete) {
		if (isThreeDimensional()) {
			blackList = extractCellEdgeNeighs(id);
		} else {
			blackList = extractCellFaceNeighs(id);
		}
	}

	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getVertexCount(); ++i) {
		for (auto &neigh : extractCellVertexNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
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
std::vector<long> Patch::extractCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	std::vector<int> vertexList(1);
	vertexList[0] = vertex;

	return extractCellVertexNeighs(id, vertexList, blackList);
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
std::vector<long> Patch::extractCellVertexNeighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;

	int nVerticesToFound = vertices.size();

	const Cell &cell = getCell(id);
	const long *cellConnect = cell.getConnect();

	std::vector<long> alreadyScanned;
	std::vector<long> processingQueue;
	processingQueue.push_back(cell.get_id());
	while (!processingQueue.empty()) {
		// Get a cell to scan and remove it form the list
		long scanId(processingQueue.back());
		processingQueue.pop_back();
		const Cell &scanCell = getCell(scanId);

		// Scan the interfaces of the cell
		const long *interfaces = scanCell.getInterfaces();
		for (int i = 0; i < scanCell.getInterfaceCount(); i++) {
			long interfaceId = interfaces[i];
			const Interface &interface = getInterface(interfaceId);

			// Neighbour cell assocated to the interface
			//
			// Only consider the cells that are not
			long neighId = interface.getNeigh();
			if (neighId < 0 || neighId == scanId) {
				neighId = interface.getOwner();
			}

			if (neighId == id) {
				continue;
			} else if (std::find(alreadyScanned.begin(), alreadyScanned.end(), neighId) != alreadyScanned.end()) {
				continue;
			}

			// Number of vertices owned by the interface
			int nCommonVertices = 0;
			const long *interfaceConnect = interface.getConnect();
			for (int k = 0; k < interface.getVertexCount(); ++k) {
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
					utils::addToOrderedVector<long>(neighId, neighs);
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
long Patch::getInterfaceCount() const
{
	return m_interfaces.size();
}

/*!
	Gets the interfaces owned by the patch.

	\return The interfaces owned by the patch.
*/
PiercedVector<Interface> & Patch::interfaces()
{
	return m_interfaces;
}

/*!
	Gets a reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A reference to the interface with the specified id.
*/
Interface & Patch::getInterface(const long &id)
{
	return m_interfaces[id];
}

/*!
	Gets a constant reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A constant reference to the interface with the specified id.
*/
const Interface & Patch::getInterface(const long &id) const
{
	return m_interfaces[id];
}

/*!
	Adds a new interface with the specified id.

	\param id is the id of the new interface
*/
long Patch::addInterface(const long &id, ElementInfo::Type type)
{
	PiercedVector<Interface>::iterator iterator = m_interfaces.reclaim(id);

	Interface &interface = *iterator;
	interface.initialize(type);

	return id;
}

/*!
	Adds a new interface.
*/
long Patch::addInterface(ElementInfo::Type type)
{
	long id;
	if (m_unusedInterfaceIds.empty()) {
		id = m_interfaces.size();
	} else {
		id = m_unusedInterfaceIds.front();
		m_unusedInterfaceIds.pop_front();
	}

	return addInterface(id, type);
}

/*!
	Deletes an interface.

	\param id is the id of the interface
*/
void Patch::deleteInterface(const long &id, bool delayed)
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
std::array<double, 3> Patch::evalCellCentroid(const long &id)
{
	Cell &cell = getCell(id);

	return evalElementCentroid(cell);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> Patch::evalInterfaceCentroid(const long &id)
{
	Interface &interface = getInterface(id);

	return evalElementCentroid(interface);
}

/*!
	Evaluates the centroid of the specified element.

	\param element is the element
	\result The centroid of the specified element.
*/
std::array<double, 3> Patch::evalElementCentroid(const Element &element)
{
	const int nDimensions = 3;

	const long *elementConnect = element.getConnect();
	const ElementInfo &elementInfo = element.get_info();

	std::array<double, nDimensions> centroid = {{0., 0., 0.}};
	for (int i = 0; i < elementInfo.nVertices; ++i) {
		Vertex &vertex = getVertex(elementConnect[i]);
		const std::array<double, nDimensions> &vertexCoords = vertex.getCoords();
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
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool Patch::isPointInside(const double &x, const double &y, const double &z)
{
	return isPointInside({{x, y, z}});
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns the id of the cell the contains the point. If the point
	is not inside the patch, the function returns the id of the null element.
*/
long Patch::locatePoint(const double &x, const double &y, const double &z)
{
	return locatePoint({{x, y, z}});
}

/*!
	Updates the stored patch bounding box.
*/
void Patch::updateBoundingBox()
{
	evalBoundingBox(m_minPoint, m_maxPoint);
}

/*!
	Gets the previously stored patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void Patch::getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	minPoint = m_minPoint;
	maxPoint = m_maxPoint;
}

/*!
	Evalautes patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void Patch::evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	// Initialize bounding box
	for (int k = 0; k < 3; ++k) {
		minPoint[k] =   std::numeric_limits<double>::max();
		maxPoint[k] = - std::numeric_limits<double>::max();
	}

	// Compute bounding box limits
	PiercedVector<Vertex>::iterator e = m_vertices.end();
	for (PiercedVector<Vertex>::iterator i = m_vertices.begin(); i != e; ++i) {
		for (int k = 0; k < 3; ++k) {
			double value = i->getCoords()[k];

			minPoint[k] = std::min(minPoint[k], value);
			maxPoint[k] = std::max(maxPoint[k], value);
		}
	}
}

/*!
	Sort patch vertices on regular bins.

	\param[in] n_bins (default = 128) is the number of bins (on each space
	direction)
	\result Returns the bin index associated to each vertex.
*/
std::unordered_map<long, long> Patch::binSortVertex(int nBins)
{
	// ====================================================================== //
	// VARIABLES DECLARATION                                                  //
	// ====================================================================== //

	// Local variables
	double                              dx, dy, dz;

	// Counters
	long                                i, j, k;
	PiercedVector<Vertex>::iterator     V, E = m_vertices.end();

	// ====================================================================== //
	// ASSOCIATE EACH VERTEX WITH A BIN                                       //
	// ====================================================================== //

	// Bin's spacing
	dx = max(1.0e-12, m_maxPoint[0] - m_minPoint[0]) / ((double) nBins);
	dy = max(1.0e-12, m_maxPoint[1] - m_minPoint[1]) / ((double) nBins);
	dz = max(1.0e-12, m_maxPoint[2] - m_minPoint[2]) / ((double) nBins);

	// Loop over vertices
	std::unordered_map<long, long> bin_index;
	for (V = m_vertices.begin(); V != E; ++V) {
		i = std::min(nBins - 1L, long((V->getCoords()[0] - m_minPoint[0]) / dx));
		j = std::min(nBins - 1L, long((V->getCoords()[1] - m_minPoint[1]) / dy));
		k = std::min(nBins - 1L, long((V->getCoords()[2] - m_minPoint[2]) / dz));
		bin_index[V->get_id()] = nBins * nBins * k + nBins * j + i;
	}

	return bin_index;
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
*/
void Patch::translate(std::array<double, 3> translation)
{
	for (auto &vertex : m_vertices) {
		vertex.translate(translation[0], translation[1], translation[2]);
	}
}

/*!
	Translates the patch.

	\param[in] sx translation along x direction
	\param[in] sy translation along y direction
	\param[in] sy translation along z direction
*/
void Patch::translate(double sx, double sy, double sz)
{
	translate({{sx, sy, sz}});
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor vector
*/
void Patch::scale(std::array<double, 3> scaling)
{
	for (auto &vertex : m_vertices) {
		vertex.scale(scaling[0], scaling[1], scaling[2]);
	}
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor
*/
void Patch::scale(double scaling)
{
	scale({{scaling, scaling, scaling}});
}

/*!
	Scales the patch.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sy scaling factor along z direction
*/
void Patch::scale(double sx, double sy, double sz)
{
	scale({{sx, sy, sz}});
}

/*!
	Sets the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void Patch::setTol(double tolerance)
{
	_setTol(tolerance);

	m_hasCustomTolerance = true;
}

/*!
	Internal function to set the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void Patch::_setTol(double tolerance)
{
	m_tolerance = tolerance;
}

/*!
	Gets the tolerance for the geometrical checks.

	\result The tolerance fot the geometrical checks.
*/
double Patch::getTol() const
{
	return m_tolerance;
}

/*!
	Resets the tolerance for the geometrical checks.
*/
void Patch::resetTol()
{
	_resetTol();

	m_hasCustomTolerance = false;
}

/*!
	Internal function to reset the tolerance for the geometrical checks.
*/
void Patch::_resetTol()
{
	m_tolerance = 1;
	for (int k = 0; k < 3; ++k) {
		m_tolerance = std::max(m_maxPoint[0] - m_minPoint[0], m_tolerance);
	}
	m_tolerance *= 1e-14;
}

/*!
	Checks if the tolerance for the geometrical checks has been customized
	by the user.

	\result True if the tolerance was customized by the user, false otherwise.
*/
bool Patch::isTolCustomized() const
{
	return m_hasCustomTolerance;
}

/*!
 *  Interface method for obtaining field meta Data
 *
 *  @param[in] name is the name of the field to be written
 *  @return Returns a VTKFieldMetaData struct containing the metadata
 *  of the requested custom data.
 */
const VTKFieldMetaData Patch::getMetaData(std::string name)
{
	if (name == "Points") {
		return VTKFieldMetaData(3 * m_vertices.size(), typeid(double));
	} else if (name == "offsets") {
		return VTKFieldMetaData(m_cells.size(), typeid(int));
	} else if (name == "types") {
		return VTKFieldMetaData(m_cells.size(), typeid(VTKElementType));
	} else if (name == "connectivity") {
		long connectSize = 0;
		for (Cell &cell : m_cells) {
			connectSize += cell.get_info().nVertices;
		}

		return VTKFieldMetaData(connectSize, typeid(long));
	} else if (m_dataFields.count(name) > 0) {
		long fieldSize = 0;

		if (m_dataLocations[name] == VTKLocation::CELL) {
			fieldSize = m_cells.size();
		} else {
			fieldSize = m_vertices.size();
		}

		if (m_dataType[name] == VTKFieldType::VECTOR) {
			fieldSize *= 3;
		}

		return VTKFieldMetaData(fieldSize, typeid(double));
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
void Patch::flushData(std::fstream &stream, VTKFormat format, std::string name)
{
	assert(format == VTKFormat::APPENDED);

	static std::unordered_map<long, long> vertexMap;

	if (name == "Points") {
		long vertexId = 0;
		for (Vertex &vertex : m_vertices) {
			vertexMap[vertex.get_id()] = vertexId++;

			genericIO::flushBINARY(stream, vertex.getCoords());
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (Cell &cell : m_cells) {
			offset += cell.get_info().nVertices;
			genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "types") {
		for (Cell &cell : m_cells) {
			VTKElementType VTKType;
			switch (cell.getType())  {

			case ElementInfo::VERTEX:
				VTKType = VTKElementType::VERTEX;
				break;

			case ElementInfo::LINE:
				VTKType = VTKElementType::LINE;
				break;

			case ElementInfo::TRIANGLE:
				VTKType = VTKElementType::TRIANGLE;
				break;

			case ElementInfo::PIXEL:
				VTKType = VTKElementType::PIXEL;
				break;

			case ElementInfo::QUAD:
				VTKType = VTKElementType::QUAD;
				break;

			case ElementInfo::TETRA:
				VTKType = VTKElementType::TETRA;
				break;

			case ElementInfo::VOXEL:
				VTKType = VTKElementType::VOXEL;
				break;

			case ElementInfo::HEXAHEDRON:
				VTKType = VTKElementType::HEXAHEDRON;
				break;

			case ElementInfo::WEDGE:
				VTKType = VTKElementType::WEDGE;
				break;

			case ElementInfo::PYRAMID:
				VTKType = VTKElementType::PYRAMID;
				break;

			default:
				VTKType = VTKElementType::UNDEFINED;
				break;

			}

			genericIO::flushBINARY(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (Cell &cell : m_cells) {
			for (int i = 0; i < cell.get_info().nVertices; ++i) {
				genericIO::flushBINARY(stream, vertexMap.at(cell.getVertex(i)));
			}
		}

		vertexMap.clear();
		std::unordered_map<long, long>().swap(vertexMap);
	} else if (m_dataFields.count(name) > 0) {
		genericIO::flushBINARY(stream, *(m_dataFields.at(name)));
	}
}

/*!
	@}
*/

}
