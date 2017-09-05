/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#include <cassert>
#include <limits>
#include <set>

#include "bitpit_common.hpp"
#include "bitpit_operators.hpp"

#include "element.hpp"

/*!
	Input stream operator for class Element

	\param[in] buffer is the input stream
	\param[in] element is the element to be streamed
	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::Element &element)
{
	// Initialize the element
	bitpit::ElementType type;
	buffer >> type;

	long id;
	buffer >> id;

	element._initialize(id, type);

	// Set the connectivity
	int nVertices = element.getVertexCount();
	for (int i = 0; i < nVertices; ++i) {
	    buffer >> element.m_connect[i];
	}

	return buffer;
}

/*!
	Output stream operator for element

	\param[in] buffer is the output stream
	\param[in] element is the element to be streamed
	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const bitpit::Element &element)
{
	int nVertices = element.getVertexCount();
	buffer << element.getType();
	buffer << element.getId();
	for (int i = 0; i < nVertices; ++i) {
	    buffer << element.m_connect[i];
	}

	return buffer;
}

namespace bitpit {

/*!
	\class Tesselation
	\ingroup patchelements

	\brief The Tesselation class allows to tessalete polygons and polyhedrons.

	The Tesselation class allows to divde polygons and polyhedrons elements
	"regular" elements, i.e., elements that are associated to a reference
	element.
*/

/*!
	Constructor.
*/
Element::Tesselation::Tesselation()
	: m_nTiles(0)
{
}

/*!
	Import the specified vertex coordinates in the tesselation.

	\param coordinates are the coordinates of the vertex
	\result The id associated by the tesselation to the imported vertex
	coordinates.
*/
int Element::Tesselation::importVertexCoordinates(const std::array<double, 3> &coordinates)
{
    return importVertexCoordinates(std::array<double, 3>(coordinates));
}

/*!
	Import the specified vertex coordinates in the tesselation.

	\param coordinates are the coordinates of the vertex
	\result The id associated by the tesselation to the imported vertex
	coordinates.
*/
int Element::Tesselation::importVertexCoordinates(std::array<double, 3> &&coordinates)
{
    m_coordinates.push_back(std::move(coordinates));

    return (m_coordinates.size() - 1);
}

/*!
	Import the specified vertex coordinates in the tesselation.

	\param coordinates are the coordinates of the vertices
	\param nVertices is the number of vertices
	\result The ids associated by the tesselation to the imported vertex
	coordinates.
*/
std::vector<int> Element::Tesselation::importVertexCoordinates(const std::array<double, 3> *coordinates, int nVertices)
{
    int nStoredVertices = m_coordinates.size();
    m_coordinates.reserve(nStoredVertices + nVertices);

    std::vector<int> ids(nVertices);
    for (int k = 0; k < nVertices; ++k) {
        m_coordinates.push_back(coordinates[k]);
        ids[k] = nStoredVertices++;
    }

    return ids;
}

/*!
	Import the specified polygon in the tesselation.

	\param vertexIds are the ids of the polygon's vertices
*/
void Element::Tesselation::importPolygon(const std::vector<int> &vertexIds)
{
	int nVertices = vertexIds.size();
	if (nVertices == 3 || nVertices == 4) {
		m_nTiles++;
		m_types.push_back(nVertices == 3 ? ElementType::TRIANGLE : ElementType::QUAD);
		m_connects.push_back(vertexIds);

		return;
	}

	// Add the centroid
	std::array<double, 3> centroid = {{0., 0., 0.}};
	for (int k = 0; k < nVertices; ++k) {
		centroid += m_coordinates[vertexIds[k]];
	}
	centroid = centroid / double(nVertices);

	int centroidId = importVertexCoordinates(std::move(centroid));

	// Decompose the polygon in triangles
	//
	// Each triangle is composed by the two vertices of a side and by the
	// centroid.
	ElementType tileType = ElementType::TRIANGLE;
	int nTileVertices = ReferenceElementInfo::getInfo(tileType).nVertices;

	int nSides = nVertices;
	m_types.resize(m_nTiles + nSides, tileType);
	m_connects.resize(m_nTiles + nSides, std::vector<int>(nTileVertices));
	for (int i = 0; i < nSides; ++i) {
		m_nTiles++;
		for (int k = 0; k < nTileVertices; ++k) {
		m_connects[i][k] = (i + k) % nVertices;
		}
		m_connects[i][nTileVertices] = centroidId;
	}
}

/*!
	Import the specified polygon in the tesselation.

	\param vertexIds are the ids of the polygon's vertices
	\param faceVertexIds are the ids of the polygon's face vertices
*/
void Element::Tesselation::importPolyhedron(const std::vector<int> &vertexIds, const std::vector<std::vector<int>> &faceVertexIds)
{
	int nFaces    = faceVertexIds.size();
	int nVertices = vertexIds.size();

	// Generate the tesselation of the surface
	for (int i = 0; i < nFaces; ++i) {
		importPolygon(faceVertexIds[i]);
	}

	// Add the centroid of the element to the tesselation
	std::array<double, 3> centroid = {{0., 0., 0.}};
	for (int k = 0; k < nVertices; ++k) {
		centroid += m_coordinates[vertexIds[k]];
	}
	centroid = centroid / double(nVertices);

	int centroidTesselationId = importVertexCoordinates(std::move(centroid));

	// Decompose the polyhedron in prisms and pyramids.
	//
	// Each prism/pyramid has a tile of the surface tesselation as the
	// base and the centroid of the element as the apex. Since we have
	// already generated the surface tesselation, we can "convert" the
	// two-dimensional tiles of that tesselation in three-dimensional
	// tiles to obtain the volume tesselation.
	for (int i = 0; i < m_nTiles; ++i) {
		// Change the tile type
		ElementType &tileType = m_types[i];
		if (tileType == ElementType::TRIANGLE) {
			tileType = ElementType::TETRA;
		} else if (tileType == ElementType::QUAD) {
			tileType = ElementType::PYRAMID;
		} else {
			BITPIT_UNREACHABLE("Unsupported tile");
			throw std::runtime_error ("Unsupported tile");
		}

		// Fix the order of the connectivity the match the order of the
		// three-dimensional element
		if (tileType == ElementType::TETRA) {
			std::swap(m_connects[i][0], m_connects[i][2]);
		} else if (tileType == ElementType::PYRAMID) {
			std::swap(m_connects[i][1], m_connects[i][3]);
		}

		// Add the centroid to the tile connectivity
		m_connects[i].push_back(centroidTesselationId);
	}
}

/*!
	Get the tiles contained in the tesselation.

	\result The tiles contained in the tesselation.
*/
int Element::Tesselation::getTileCount() const
{
    return m_nTiles;
}

/*!
	Get the type of the specified tile.

	\param tile is the tile
	\result The type of the specified tile.
*/
ElementType Element::Tesselation::getTileType(int tile) const
{
    return m_types[tile];
}

/*!
	Get the coordinates of the vertices of the specified tile.

	\param tile is the tile
	\result The coordinates of the vertices of the specified tile..
*/
std::vector<std::array<double, 3>> Element::Tesselation::getTileVertexCoordinates(int tile) const
{
    const ElementType tileType = getTileType(tile);
    const int nTileVertices = ReferenceElementInfo::getInfo(tileType).nVertices;
    const std::vector<int> &tileConnect = m_connects[tile];

    std::vector<std::array<double, 3>> coordinates(nTileVertices);
    for (int i = 0; i < nTileVertices; ++i) {
        coordinates[i] = m_coordinates[tileConnect[i]];
    }

    return coordinates;
}

/*!
	\class Element
	\ingroup patchelements

	\brief The Element class provides an interface for defining elements.

	Element is the base class for defining elements like cells and
	intefaces.
*/

const long Element::NULL_ID = std::numeric_limits<long>::min();

/*!
	Default constructor.
*/
Element::Element()
{
	_initialize(NULL_ID, ElementType::UNDEFINED);
}

/*!
	Creates a new element.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
*/
Element::Element(long id, ElementType type, int connectSize)
{
	_initialize(id, type, connectSize);
}

/*!
	Creates a new element.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
Element::Element(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage)
{
	_initialize(id, type, std::move(connectStorage));
}

/*!
	Copy constructor
*/
Element::Element(const Element &other)
{
	_initialize(other.m_id, other.m_type);

	if (other.m_connect) {
		int connectSize = getInfo().nVertices;
		std::copy(other.m_connect.get(), other.m_connect.get() + connectSize, m_connect.get());
	}
}

/*!
	Copy-assignament operator.
*/
Element & Element::operator=(const Element &other)
{
	Element tmp(other);
	swap(tmp);

	return *this;
}

/**
* Exchanges the content of the element by the content the specified other
* element.
*
* \param other is another element whose content is swapped with that of this
* element
*/
void Element::swap(Element &other) noexcept
{
	std::swap(other.m_id, m_id);
	std::swap(other.m_type, m_type);
	std::swap(other.m_connect, m_connect);
}

/*!
	Initializes the data structures of the element.

	\param id the id of the element
	\param type the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
*/
void Element::initialize(long id, ElementType type, int connectSize)
{
	_initialize(id, type, connectSize);
}

/*!
	Initializes the data structures of the element.

	\param type the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
void Element::initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage)
{
	_initialize(id, type, std::move(connectStorage));
}

/*!
	Internal function to initialize the data structures of the element.

	\param id the id of the element
	\param type the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
*/
void Element::_initialize(long id, ElementType type, int connectSize)
{
	// Get previous connect size
	int previousConnectSize = 0;
	if (m_connect) {
		if (hasInfo()) {
			previousConnectSize = getInfo().nVertices;
		}
	}

	// Initialize connectivity storage
	if (ReferenceElementInfo::hasInfo(type)) {
		connectSize = ReferenceElementInfo::getInfo(type).nVertices;
	}

	std::unique_ptr<long[]> connectStorage;
	if (connectSize != previousConnectSize) {
		connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	} else {
		connectStorage = std::move(m_connect);
	}

	// Initialize element
	_initialize(id, type, std::move(connectStorage));
}

/*!
	Internal function to initialize the data structures of the element.

	\param id is the ID of the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
void Element::_initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage)
{
	// Set the id
	setId(id);

	// Set type
	setType(type);

	// Initialize connectivity
	setConnect(std::move(connectStorage));
}

/*!
	Sets the ID of the element.

	\param id the ID of the element
*/
void Element::setId(const long &id)
{
	m_id = id;
}

/*!
	Gets the ID of the element.

	\return The ID of the element
*/
long Element::getId() const
{
	return m_id;
}

/*!
	Check if the element is associated to a reference element.

	\result Returns true if the element is associated to a reference element,
	false otherwise.
*/
bool Element::hasInfo() const
{
	return ReferenceElementInfo::hasInfo(m_type);
}

/*!
	Gets the basic information of the element.

	\result A constant reference to the basic information of the element.
*/
const ReferenceElementInfo & Element::getInfo() const
{
	return ReferenceElementInfo::getInfo(m_type);
}

/*!
	Sets the element type.

	\param type the element type
*/
void Element::setType(ElementType type)
{
	m_type = type;
}

/*!
	Gets the element type.

	\result The element type
*/
ElementType Element::getType() const
{
	return m_type;
}

/*!
	Sets the vertex connectivity of the element.

	\param connect a pointer to the connectivity of the element
*/
void Element::setConnect(std::unique_ptr<long[]> &&connect)
{
	m_connect = std::move(connect);
}

/*!
	Unsets the vertex connectivity of the element.
*/
void Element::unsetConnect()
{
	m_connect.reset(nullptr);
}

/*!
	Gets the vertex connectivity of the element.

	\result A constant pointer to the connectivity of the element
*/
const long * Element::getConnect() const
{
	return m_connect.get();
}

/*!
	Gets the vertex connectivity of the element.

	\result A pointer to the connectivity of the element
*/
long * Element::getConnect()
{
	return m_connect.get();
}

/*!
	Given the id of a vertex, evaluates thhe local index of that vertex within
	the element. If the specified vertex does not exist in the element
	connectivity list, a negative number is returned.

	\param vertexId is the vertex id
	\result The local index of the vertex if the element contains the vertex,
	a negative number otherwise.
*/
int Element::findVertex(long vertexId) const
{
	ConstProxyVector<long> cellVertexIds = getVertexIds();

	int localVertexId = std::distance(cellVertexIds.begin(), std::find(cellVertexIds.begin(), cellVertexIds.end(), vertexId));
	if (localVertexId >= (int) cellVertexIds.size()) {
		return -1;
	}

	return localVertexId;
}

/*!
	Gets the size of the connectivity of the element.

	\result The size of the connectivity of the element.
*/
int Element::getConnectSize() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");

	default:
		return getVertexCount();

	}
}

/*!
	Gets the number of faces of the element.

	\result The number of vertices of the element
*/
int Element::getFaceCount() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return getInfo().nFaces;

	}
}

/*!
	Gets the face type of the specified face of the element.

	\result The face type of specified face of the element
*/
ElementType Element::getFaceType(const int &face) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return ElementType::UNDEFINED;

	default:
		return getInfo().face_type[face];

	}
}

/*!
	Gets the number of vertices of the specified face.

	\param face is the face for which the number of vertices is requested
	\result The number of vertices of the specified face.
*/
int Element::getFaceVertexCount(const int &face) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
    {
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
    }

	default:
    {
		ElementType faceType = getFaceType(face);

		return ReferenceElementInfo::getInfo(faceType).nVertices;
    }

	}
}

/*!
	Gets the local connectivity of the specified face of the element.

	\param face is the face for which the connectivity is reqested
	\result The local connectivity of the specified face of the element.
*/
ConstProxyVector<int> Element::getFaceLocalConnect(const int &face) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		const std::vector<int> &localFaceConnect = getInfo().faceConnect[face];

		return ConstProxyVector<int>(localFaceConnect.data(), localFaceConnect.size());
	}

	}
}

/*!
	Gets the connectivity of the specified face of the element.

	\param face is the face for which the connectivity is reqested
	\result The connectivity of the specified face of the element.
*/
ConstProxyVector<long> Element::getFaceConnect(int face) const
{
	const long *connectivity = getConnect();

	ConstProxyVector<int> localFaceConnect = getFaceLocalConnect(face);
	int faceConnectSize = localFaceConnect.size();

	std::vector<long> faceConnect(faceConnectSize);
	for (int k = 0; k < faceConnectSize; ++k) {
		int localVertexId = localFaceConnect[k];
		long vertexId = connectivity[localVertexId];
		faceConnect[k] = vertexId;
	}

	return ConstProxyVector<long>(std::move(faceConnect));
}

/*!
	Gets the number of edges of the element.

	\result The number of edges of the element
*/
int Element::getEdgeCount() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return getInfo().nEdges;

	}
}

/*!
	Gets the local connectivity of the specified edge of the element.

	\param edge is the edge for which the connectivity is reqested
	\result The local connectivity of the specified edge of the element.
*/
ConstProxyVector<int> Element::getEdgeLocalConnect(const int &edge) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		const std::vector<int> &localEdgeConnect = getInfo().edgeConnect[edge];

		return ConstProxyVector<int>(localEdgeConnect.data(), localEdgeConnect.size());
	}

	}
}

/*!
	Gets the connectivity of the specified edge of the element.

	\param edge is the edge for which the connectivity is reqested
	\result The connectivity of the specified edge of the element.
*/
ConstProxyVector<long> Element::getEdgeConnect(int edge) const
{
	const long *connectivity = getConnect();

	ConstProxyVector<int> localEdgeConnect = getEdgeLocalConnect(edge);
	int nEdgeVertices = localEdgeConnect.size();

	std::vector<long> edgeConnect(nEdgeVertices);
	for (int k = 0; k < nEdgeVertices; ++k) {
		int localVertexId = localEdgeConnect[k];
		edgeConnect[k] = connectivity[localVertexId];
	}

	return ConstProxyVector<long>(std::move(edgeConnect));
}

/*!
	Gets the dimension of the element.

	\param type the type of the element
	\return The dimension of the element
*/
int Element::getDimension(ElementType type)
{
	switch (type) {

	case ElementType::POLYGON:
		return 2;

	case ElementType::POLYHEDRON:
		return 3;

	case ElementType::UNDEFINED:
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return ReferenceElementInfo::getInfo(type).dimension;

	}
}

/*!
	Gets the dimension of the element.

	\return The dimension of the element
*/
int Element::getDimension() const
{
	return getDimension(m_type);
}

/*!
	Returns true if the element is a three-dimensional element.

	\param type the type of the element
	\return Returns true if the element is a three-dimensional element,
	false otherwise.
*/
bool Element::isThreeDimensional(ElementType type)
{
	return (getDimension(type) == 3);
}

/*!
	Returns true if the element is a three-dimensional element.

	\return Returns true if the element is a three-dimensional element,
	false otherwise.
*/
bool Element::isThreeDimensional() const
{
	return (getDimension() == 3);
}

/*!
	Gets the number of vertices of the element.

	\result The number of vertices of the element
*/
int Element::getVertexCount() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return getInfo().nVertices;

	}
}

/*!
	Gets the list of the vertex ids.

	\result The list of the vertex ids.
*/
ConstProxyVector<long> Element::getVertexIds() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		return ConstProxyVector<long>(getConnect(), getVertexCount());
	}

	}
}

/*!
	Gets the list of vertex ids for the specified face of the element.

	\param face is the face for which the vertex ids is reqested
	\result The list of vertex ids for the specified face of the element.
*/
ConstProxyVector<long> Element::getFaceVertexIds(int face) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		ConstProxyVector<int> localFaceVertexIds = getFaceLocalConnect(face);
		int nFaceVertices = localFaceVertexIds.size();

		ConstProxyVector<long> cellVertexIds = getVertexIds();

		std::vector<long> faceVertexIds(nFaceVertices);
		for (int k = 0; k < nFaceVertices; ++k) {
			int localVertexId = localFaceVertexIds[k];
			long vertexId = cellVertexIds[localVertexId];
			faceVertexIds[k] = vertexId;
		}

		return ConstProxyVector<long>(std::move(faceVertexIds));
	}

	}
}

/*!
	Renumber the vertices of a cell.

	\param map is the map that will be used for the renumbering
*/
void Element::renumberVertices(const std::unordered_map<long, long> &map)
{
	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
    {
		int nVertices = getVertexCount();
		long *connectivity = getConnect();
		for (int k = 0; k < nVertices; ++k) {
			connectivity[k] = map.at(connectivity[k]);
		}

		break;
	}

	}
}

/*!
	Evaluates the characteristics size of the element.

	\param coordinates are the coordinate of the vertices
	\result The characteristics size of the element.
*/
double Element::evalSize(const std::array<double, 3> *coordinates) const
{
	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		return 0.;
	}

	default:
    {
		const ReferenceElementInfo &referenceInfo = static_cast<const ReferenceElementInfo &>(getInfo());

		return referenceInfo.evalSize(coordinates);
	}

	}
}


/*!
	Evaluates the volume of the element.

	\param coordinates are the coordinate of the vertices
	\result The volume of the element.
*/
double Element::evalVolume(const std::array<double, 3> *coordinates) const
{
	assert(getDimension() == 3);

	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		const Reference3DElementInfo &referenceInfo = static_cast<const Reference3DElementInfo &>(getInfo());

		return referenceInfo.evalVolume(coordinates);

		break;
	}

	}
}

/*!
	Evaluates the area of the element.

	\param coordinates are the coordinate of the vertices
	\result The area of the specified element.
*/
double Element::evalArea(const std::array<double, 3> *coordinates) const
{
	assert(getDimension() == 2);

	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		const Reference2DElementInfo &referenceInfo = static_cast<const Reference2DElementInfo &>(getInfo());

		return referenceInfo.evalArea(coordinates);
	}

	}
}

/*!
	Evaluates the length of the element.

	\param coordinates are the coordinate of the vertices
	\result The length of the element.
*/
double Element::evalLength(const std::array<double, 3> *coordinates) const
{
	assert(getDimension() == 1);

	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		return 0.;
	}

	default:
	{
		const Reference1DElementInfo &referenceInfo = static_cast<const Reference1DElementInfo &>(getInfo());

		return referenceInfo.evalLength(coordinates);
	}

	}
}

/*!
	Evaluates the normal of an element.

	\param coordinates are the coordinate of the vertices
	\param orientation is a vector carring the additional information needed
	to un-ambigously define a normal to the element (e.g., when evaluating
	the normal of a one-dimensional element, this versor is perpendicular to
	the plane where the normal should lie)
	\param point are the element reference coordinates of the point where the
	normal should be evaluated
	\result The normal of the element.
*/
std::array<double, 3> Element::evalNormal(const std::array<double, 3> *coordinates,
										  const std::array<double, 3> &orientation,
										  const std::array<double, 3> &point) const
{
	assert(getDimension() != 3);

	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");
	}

	default:
	{
		int dimension = getDimension();
		if (dimension == 2) {
			const Reference2DElementInfo &referenceInfo = static_cast<const Reference2DElementInfo &>(getInfo());

			return referenceInfo.evalNormal(coordinates, point);
		} else if (dimension == 1) {
			const Reference1DElementInfo &referenceInfo = static_cast<const Reference1DElementInfo &>(getInfo());

			return referenceInfo.evalNormal(coordinates, orientation, point);
		} else {
			return orientation;
		}
	}

	}
}

/*!
	Generate a tesselation for the element.

	\param coordinates are the coordinate of the vertices
	\result A tesselation for the element.
*/
Element::Tesselation Element::generateTesselation(const std::array<double, 3> *coordinates) const
{
	Tesselation tesselation;

	// Add the coordinates of the vertices to the tesselation
	int nVertices = getVertexCount();
	std::vector<int> vertexTesselationIds = tesselation.importVertexCoordinates(coordinates, nVertices);

	// Generate the tesselation
	ElementType type = getType();
	switch(type) {

	case ElementType::POLYGON:
	{
		tesselation.importPolygon(vertexTesselationIds);

		break;
	}

	case ElementType::POLYHEDRON:
	{
		int nFaces = getFaceCount();
		std::vector<std::vector<int>> faceTesselationIds(nFaces);
		for (int i = 0; i < nFaces; ++i) {
			ConstProxyVector<int> localFaceConnect = getFaceLocalConnect(i);
			int nFaceVertices = getFaceVertexCount(i);

			faceTesselationIds[i].resize(nFaceVertices);
			for (int k = 0; k < nFaceVertices; ++k) {
				faceTesselationIds[i][k] = vertexTesselationIds[localFaceConnect[k]];
			}
		}

		tesselation.importPolyhedron(vertexTesselationIds, faceTesselationIds);

		break;
	}

	default:
	{
		assert(ReferenceElementInfo::hasInfo(type));

		tesselation.m_nTiles = 1;
		tesselation.m_types.push_back(type);
		tesselation.m_connects.push_back(vertexTesselationIds);

		break;
	}

	}

	return tesselation;
}

/*!
	Evaluates the connectivity of all edges.

	This function does not use the information of the reference element, so it
	is slow and should only be used for polyhedral elements.

	\param nRequestedEdges is the number of edges to extract
	\result The connectivity of all edges.
*/
std::vector<ConstProxyVector<long>> Element::evalEdgeConnects(int nRequestedEdges) const
{
	if (nRequestedEdges == -1) {
		nRequestedEdges = getEdgeCount();
	}
	assert(nRequestedEdges <= getEdgeCount());

	std::set<std::pair<long, long>> edgeSet;
	std::vector<ConstProxyVector<long>> edgeConnects(nRequestedEdges);

	int nFaces = getFaceCount();
	for (int i = 0; i < nFaces; ++i) {
		ConstProxyVector<long> faceVertexIds = getFaceVertexIds(i);
		int nFaceVertices = faceVertexIds.size();
		for (int k = 0; k < nFaceVertices; ++k) {
			long vertex_A = faceVertexIds[k];
			long vertex_B = faceVertexIds[(k + 1) % nFaceVertices];
			if (vertex_A > vertex_B) {
				std::swap(vertex_A, vertex_B);
			}

			std::pair<long, long> edgePair = std::pair<long, long>(vertex_A, vertex_B);
			std::pair<std::set<std::pair<long, long>>::iterator, bool> insertResult = edgeSet.insert(edgePair);
			if (insertResult.second) {
				std::vector<long> edgeConnect(2, edgePair.first);
				edgeConnect[1] = edgePair.second;
				edgeConnects[edgeSet.size() - 1] = ConstProxyVector<long>(std::move(edgeConnect));

				if (edgeSet.size() == (std::size_t) nRequestedEdges) {
					return edgeConnects;
				}
			}
		}
	}

	assert(edgeSet.size() == nRequestedEdges);

	return edgeConnects;
}

/*!
        Returns the buffer size required to communicate cell data

        \result buffer size (in bytes)
*/
unsigned int Element::getBinarySize()
{
	return (sizeof(ElementType) + (getVertexCount() + 1) * sizeof(long));
}

// Explicit instantiation of the Element containers
template class PiercedVector<Element>;

}
