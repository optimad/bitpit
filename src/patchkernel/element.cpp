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
#include "bitpit_CG.hpp"
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

	int connectSize;
	if (bitpit::ReferenceElementInfo::hasInfo(type)) {
		element._initialize(id, type);
		connectSize = element.getConnectSize();
	} else {
		buffer >> connectSize;
		element._initialize(id, type, connectSize);
	}

	// Set connectivity
	for (int i = 0; i < connectSize; ++i) {
	    buffer >> element.m_connect[i];
	}

	// Set PID
	int pid;
	buffer >> pid;

	element.setPID(pid);

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
	buffer << element.getType();
	buffer << element.getId();

	int connectSize = element.getConnectSize();
	if (!bitpit::ReferenceElementInfo::hasInfo(element.m_type)) {
		buffer << connectSize;
	}

	for (int i = 0; i < connectSize; ++i) {
	    buffer << element.m_connect[i];
	}

	buffer << element.getPID();

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
	int nSideVertices = ReferenceElementInfo::getInfo(ElementType::LINE).nVertices;

	int nSides = nVertices;
	m_types.resize(m_nTiles + nSides, tileType);
	m_connects.resize(m_nTiles + nSides, std::vector<int>(nTileVertices));
	for (int i = 0; i < nSides; ++i) {
		m_nTiles++;
		for (int k = 0; k < nSideVertices; ++k) {
			m_connects[m_nTiles - 1][k] = vertexIds[(i + k) % nVertices];
		}
		m_connects[m_nTiles - 1][nSideVertices] = centroidId;
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
	int nInitialTiles = m_nTiles;
	for (int i = 0; i < nFaces; ++i) {
		importPolygon(faceVertexIds[i]);
	}
	int nFinalTiles = m_nTiles;

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
	for (int i = nInitialTiles; i < nFinalTiles; ++i) {
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

	\param other is another element whose content is copied in this element
*/
Element::Element(const Element &other)
{
	int connectSize;
	if (other.m_connect) {
		connectSize = other.getConnectSize();
	} else {
		connectSize = 0;
	}

	_initialize(other.m_id, other.m_type, connectSize);

	m_pid = other.m_pid;

	if (other.m_connect) {
		std::copy(other.m_connect.get(), other.m_connect.get() + connectSize, m_connect.get());
	}
}

/*!
	Copy-assignament operator.

	\param other is another element whose content is copied in this element
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
	std::swap(other.m_pid, m_pid);
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

	\param id the id of the element
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

	// Initialize PID
	setPID(0);

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
	Sets the PID associated to the element.

	\param pid is the PID associated to the element.
*/
void Element::setPID(int pid)
{
	m_pid = pid;
}

/*!
	Gets the PID associated to the element.

	\result The PID associated to the element.
*/
int Element::getPID() const
{
	return m_pid;
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
	Gets the size of the connectivity of the element.

	\result The size of the connectivity of the element.
*/
int Element::getConnectSize() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
		return 1 + getVertexCount();

	case (ElementType::POLYHEDRON):
		return getFaceStreamSize();

	default:
		assert(m_type != ElementType::UNDEFINED);

		return getVertexCount();

	}
}

/*!
	Checks if the connectivity of this element and the connectivity of the
	other element are the same.

	\param other is the other element
	\result True if the connectivity of this element and the connectivity of
	the other element are the same, false otherwise.
*/
bool Element::hasSameConnect(const Element &other) const
{
    int cellConnectSize = getConnectSize();
    if (other.getConnectSize() != cellConnectSize) {
        return false;
    }

    const long *cellConnect  = getConnect();
    const long *otherConnect = other.getConnect();
    for (int k = 0; k < cellConnectSize; ++k) {
        if (cellConnect[k] != otherConnect[k]) {
            return false;
        }
    }

    return true;
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
		return getConnect()[0];

	default:
		assert(m_type != ElementType::UNDEFINED);

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
	{
		return ElementType::LINE;
	}

	case (ElementType::POLYHEDRON):
	{
		int nFaceVertices = getFaceVertexCount(face);
		switch (nFaceVertices) {

		case 3:
			return ElementType::TRIANGLE;

		case 4:
			return ElementType::QUAD;

		default:
			return ElementType::POLYGON;

		}
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		return getInfo().face_type[face];
	}

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

	case (ElementType::POLYHEDRON):
	{
		int facePos = getFaceStreamPosition(face);
		const long *connectivity = getConnect();

		return connectivity[facePos];
	}

	default:
    {
		assert(m_type != ElementType::UNDEFINED);

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
	{
		int connectSize = getConnectSize();

		int faceConnectSize = getFaceVertexCount(face);
		std::vector<int> localFaceConnect(faceConnectSize);
		for (int i = 0; i < faceConnectSize; ++i) {
			int localVertexId = (face + i) % connectSize;

			localFaceConnect[i] = localVertexId;
		}

		return ConstProxyVector<int>(std::move(localFaceConnect));
	}

	case (ElementType::POLYHEDRON):
	{
		ElementType faceType = getFaceType(face);
		ConstProxyVector<long> faceVertexIds = getFaceVertexIds(face);
		int nFaceVertices = faceVertexIds.size();
		bool faceHasReferenceInfo = ReferenceElementInfo::hasInfo(faceType);

		int faceConnectSize  = nFaceVertices;
		int localVertexOffset = 0;
		if (!faceHasReferenceInfo) {
			++faceConnectSize;
			++localVertexOffset;
		}

		std::vector<int> localFaceConnect(faceConnectSize);
		if (!faceHasReferenceInfo) {
			localFaceConnect[0] = nFaceVertices;
		}

		// The list of the vertices in sorted, therfore we can search the
		// vertices using a binary search.
		ConstProxyVector<long> vertexIds = getVertexIds();
		for (int k = 0; k < nFaceVertices; ++k) {
			int vertexId = faceVertexIds[k];
			auto localVertexIdItr = std::lower_bound(vertexIds.begin(), vertexIds.end(), vertexId);
			assert(localVertexIdItr != vertexIds.end() && *localVertexIdItr == vertexId);
			localFaceConnect[localVertexOffset + k] = std::distance(vertexIds.begin(), localVertexIdItr);
		}

		return ConstProxyVector<int>(std::move(localFaceConnect));
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

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

	switch (m_type) {

	case (ElementType::POLYGON):
	{
		int connectSize = getConnectSize();

		int facePos         = 1 + face;
		int faceConnectSize = getFaceVertexCount(face);
		if (facePos + faceConnectSize <= connectSize) {
			return ConstProxyVector<long>(connectivity + facePos, faceConnectSize);
		} else {
			std::vector<long> faceConnect(faceConnectSize);
			for (int i = 0; i < faceConnectSize; ++i) {
				int position = facePos + i;
				if (position >= connectSize) {
					position = position % connectSize + 1;
				}

				faceConnect[i] = connectivity[position];
			}

			return ConstProxyVector<long>(std::move(faceConnect));
		}
	}

	case (ElementType::POLYHEDRON):
	{
		ElementType faceType = getFaceType(face);

		int facePos          = getFaceStreamPosition(face);
		int faceConnectSize  = connectivity[facePos];
		int faceConnectBegin = facePos + 1;
		if (!ReferenceElementInfo::hasInfo(faceType)) {
			faceConnectSize++;
			faceConnectBegin--;
		}

		return ConstProxyVector<long>(connectivity + faceConnectBegin, faceConnectSize);
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		// If we are here, the element has a reference element, therefore we
		// can retrieve the local face connectivity directly form the info
		// associated the to referenc element.
		const std::vector<int> &localFaceConnect = getInfo().faceConnect[face];
		int faceConnectSize = localFaceConnect.size();

		std::vector<long> faceConnect(faceConnectSize);
		for (int k = 0; k < faceConnectSize; ++k) {
			int localVertexId = localFaceConnect[k];
			long vertexId = connectivity[localVertexId];
			faceConnect[k] = vertexId;
		}

		return ConstProxyVector<long>(std::move(faceConnect));
	}

	}
}

/*!
	Gets the number of edges of the element.

	\result The number of edges of the element
*/
int Element::getEdgeCount() const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	{
		return getVertexCount();
	}

	case (ElementType::POLYHEDRON):
	{
		int nVertices = getVertexCount();
		int nFaces    = getFaceCount();
		int nEdges    = nVertices + nFaces - 2;

		return nEdges;
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		return getInfo().nEdges;
	}

	}
}

/*!
	Gets the type of the specified edge of the element.

	\result The type of specified edge of the element
*/
ElementType Element::getEdgeType(const int &edge) const
{
	BITPIT_UNUSED(edge);

	int dimension = getDimension();
	switch (dimension) {

	case 0:
		return ElementType::UNDEFINED;

	case 1:
	case 2:
		return ElementType::VERTEX;

	default:
		return ElementType::LINE;

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
	{
		int nEdgeVertices = ReferenceElementInfo::getInfo(ElementType::VERTEX).nVertices;

		std::vector<int> localEdgeConnect(nEdgeVertices);
		for (int k = 0; k < nEdgeVertices; ++k) {
			localEdgeConnect[k] = k;
		}

		return ConstProxyVector<int>(std::move(localEdgeConnect));
	}

	case (ElementType::POLYHEDRON):
	{
		ConstProxyVector<long> edgeVertexIds = getEdgeVertexIds(edge);
		int nEdgeVertices = edgeVertexIds.size();

		// The list of the vertices in sorted, therfore we can search the
		// vertices using a binary search.
		ConstProxyVector<long> vertexIds = getVertexIds();

		std::vector<int> localEdgeConnect(nEdgeVertices);
		for (int k = 0; k < nEdgeVertices; ++k) {
			int vertexId = edgeVertexIds[k];
			auto localVertexIdItr = std::lower_bound(vertexIds.begin(), vertexIds.end(), vertexId);
			assert(localVertexIdItr != vertexIds.end() && *localVertexIdItr == vertexId);
			localEdgeConnect[k] = std::distance(vertexIds.begin(), localVertexIdItr);
		}

		return ConstProxyVector<int>(std::move(localEdgeConnect));
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

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

	switch (m_type) {

	case (ElementType::POLYGON):
	{
		return ConstProxyVector<long>(connectivity + 1 + edge, 1);
	}

	case (ElementType::POLYHEDRON):
	{
		std::vector<ConstProxyVector<long>> edgeConnects = evalEdgeConnects(edge + 1);

		return edgeConnects[edge];
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		ConstProxyVector<int> localEdgeConnect = getEdgeLocalConnect(edge);
		int nEdgeVertices = localEdgeConnect.size();

		std::vector<long> edgeConnect(nEdgeVertices);
		for (int k = 0; k < nEdgeVertices; ++k) {
			int localVertexId = localEdgeConnect[k];
			edgeConnect[k] = connectivity[localVertexId];
		}

		return ConstProxyVector<long>(std::move(edgeConnect));
	}

	}
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

	default:
		assert(type != ElementType::UNDEFINED);

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
	{
		return getConnect()[0];
	}

	case (ElementType::POLYHEDRON):
	{
		return getVertexIds().size();
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		return getInfo().nVertices;
	}

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
	{
		return ConstProxyVector<long>(getConnect() + 1, getVertexCount());
	}

	case (ElementType::POLYHEDRON):
	{
		// The list of the vertices has to be sorted, this will make it easy
		// to extract the connectivity of the face.
		int nFaces = getFaceCount();
		const long *connectivity = getConnect();

		std::set<long> vertexIds;
		for (int i = 0; i < nFaces - 1; ++i) {
			int facePos = getFaceStreamPosition(i);

			int beginVertexPos = facePos + 1;
			int endVertexPos   = facePos + 1 + connectivity[facePos];
			vertexIds.insert(connectivity + beginVertexPos, connectivity + endVertexPos);
		}

		return ConstProxyVector<long>(std::vector<long>(vertexIds.begin(), vertexIds.end()));
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		return ConstProxyVector<long>(getConnect(), getVertexCount());
	}

	}
}

/*!
	Gets the vertex id of the specified local vertex.

	If more than one vertex is needed, the function getVertexIds may be a
	better choice. This is specially true for polygons and polyhedra, where
	the whole list of vertex ids has to be evaluated at each function call.

	\param vertex is the local index of the vertex
	\result The id of the specified vertex.
*/
long Element::getVertexId(int vertex) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	{
		ConstProxyVector<long> vertexIds = getVertexIds();

		return vertexIds[vertex];
	}

	default:
	{
		assert(m_type != ElementType::UNDEFINED);

		const long *connectivity = getConnect();

		return connectivity[vertex];
	}

	}
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
	Gets the list of vertex ids for the specified face of the element.

	\param face is the face for which the vertex ids is reqested
	\result The list of vertex ids for the specified face of the element.
*/
ConstProxyVector<long> Element::getFaceVertexIds(int face) const
{
	ConstProxyVector<long> vertexIds = getFaceConnect(face);
	if  (m_type == ElementType::POLYHEDRON) {
		ElementType faceType = getFaceType(face);
		if (faceType == ElementType::POLYGON) {
			vertexIds.set(vertexIds.data() + 1, vertexIds.size() - 1);
		}
	}

	return vertexIds;
}

/*!
	Gets the vertex id of the specified local vertex in the given face of
	the element.

	\param face is the face for which the vertex id is reqested
	\param vertex is the local index of the vertex
	\result The vertex id of the specified local vertex in the given face of
	the element.
*/
long Element::getFaceVertexId(int face, int vertex) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	{
		ConstProxyVector<long> faceVertexIds = getFaceVertexIds(face);

		return faceVertexIds[vertex];
	}

	default:
	{
		ConstProxyVector<long> cellVertexIds = getVertexIds();
		ConstProxyVector<int> faceLocalVertexIds = getFaceLocalVertexIds(face);

		return cellVertexIds[faceLocalVertexIds[vertex]];
	}

	}
}

/*!
	Gets the list of local vertex ids for the specified face of the element.

	\param face is the face for which the vertex ids is reqested
	\result The list of local vertex ids for the specified face of the element.
*/
ConstProxyVector<int> Element::getFaceLocalVertexIds(int face) const
{
	ConstProxyVector<int> localVertexIds = getFaceLocalConnect(face);
	if  (m_type == ElementType::POLYHEDRON) {
		ElementType faceType = getFaceType(face);
		if (faceType == ElementType::POLYGON) {
			localVertexIds.set(localVertexIds.data() + 1, localVertexIds.size() - 1);
		}
	}

	return localVertexIds;
}

/*!
	Gets the list of vertex ids for the specified edge of the element.

	\param edge is the edge for which the vertex ids is reqested
	\result The list of vertex ids for the specified edge of the element.
*/
ConstProxyVector<long> Element::getEdgeVertexIds(int edge) const
{
	return getEdgeConnect(edge);
}

/*!
	Gets the vertex id of the specified local vertex in the given edge of
	the element.

	\param edge is the edge for which the vertex id is reqested
	\param vertex is the local index of the vertex
	\result The vertex id of the specified local vertex in the given edge of
	the element.
*/
long Element::getEdgeVertexId(int edge, int vertex) const
{
	ConstProxyVector<long> edgeVertexIds = getEdgeVertexIds(edge);

	return edgeVertexIds[vertex];
}

/*!
	Gets the list of local vertex ids for the specified edge of the element.

	\param edge is the edge for which the vertex ids is reqested
	\result The list of local vertex ids for the specified edge of the element.
*/
ConstProxyVector<int> Element::getEdgeLocalVertexIds(int edge) const
{
	return getEdgeLocalConnect(edge);
}

/*!
	Renumber the vertices of a cell.

	\param map is the map that will be used for the renumbering
*/
void Element::renumberVertices(const std::unordered_map<long, long> &map)
{
	switch (m_type) {

	case ElementType::POLYGON:
    {
		int nVertices = getVertexCount();
		long *connectivity = getConnect();
		for (int k = 1; k < nVertices + 1; ++k) {
			connectivity[k] = map.at(connectivity[k]);
		}

		break;
	}

	case ElementType::POLYHEDRON:
    {
		int nFaces = getFaceCount();
		long *connectivity = getConnect();

		for (int i = 0; i < nFaces; ++i) {
			int facePos = getFaceStreamPosition(i);

			int beginVertexPos = facePos + 1;
			int endVertexPos   = facePos + 1 + connectivity[facePos];
			for (int k = beginVertexPos; k < endVertexPos; ++k) {
				connectivity[k] = map.at(connectivity[k]);
			}
		}

		break;
	}

	default:
    {
		assert(m_type != ElementType::UNDEFINED);

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
	switch (m_type) {

	case ElementType::POLYHEDRON:
	{
		Tesselation tesselation = generateTesselation(coordinates);
		int nTiles = tesselation.getTileCount();

		double volume = 0.;
		for (int i = 0; i < nTiles; ++i) {
			const ElementType &tileType = tesselation.getTileType(i);
			const std::vector<std::array<double, 3>> tileCoordinates = tesselation.getTileVertexCoordinates(i);
			const Reference3DElementInfo &referenceInfo = static_cast<const Reference3DElementInfo &>(ReferenceElementInfo::getInfo(tileType));

			volume += referenceInfo.evalVolume(tileCoordinates.data());
		}

		return volume;
	}

	default:
	{
		assert(getDimension() == 3);

		const Reference3DElementInfo &referenceInfo = static_cast<const Reference3DElementInfo &>(getInfo());

		return referenceInfo.evalVolume(coordinates);
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
	switch (m_type) {

	case ElementType::POLYGON:
	{
		Tesselation tesselation = generateTesselation(coordinates);
		int nTiles = tesselation.getTileCount();

		double area = 0.;
		for (int i = 0; i < nTiles; ++i) {
			const ElementType &tileType = tesselation.getTileType(i);
			const std::vector<std::array<double, 3>> tileCoordinates = tesselation.getTileVertexCoordinates(i);
			const Reference2DElementInfo &referenceInfo = static_cast<const Reference2DElementInfo &>(ReferenceElementInfo::getInfo(tileType));

			area += referenceInfo.evalArea(tileCoordinates.data());
		}

		return area;
	}

	default:
	{
		assert(getDimension() == 2);

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
	switch (m_type) {

	case ElementType::POLYGON:
	case ElementType::POLYHEDRON:
	case ElementType::UNDEFINED:
	{
		return 0.;
	}

	default:
	{
		assert(getDimension() == 1);

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
	switch (m_type) {

	case ElementType::POLYGON:
	{
		int dimension = getDimension();

		Tesselation tesselation = generateTesselation(coordinates);
		int nTiles = tesselation.getTileCount();

		double surfaceArea = 0.;
		std::array<double, 3> normal = {{0., 0., 0.}};
		for (int i = 0; i < nTiles; ++i) {
			const ElementType &tileType = tesselation.getTileType(i);
			const std::vector<std::array<double, 3>> tileCoordinates = tesselation.getTileVertexCoordinates(i);
			const Reference2DElementInfo &referenceInfo = static_cast<const Reference2DElementInfo &>(ReferenceElementInfo::getInfo(tileType));

			double tileArea = referenceInfo.evalArea(tileCoordinates.data());
			std::array<double, 3> tileNormal = {{0., 0., 0.,}};
			if (dimension == 2) {
				const Reference2DElementInfo &referenceInfo = static_cast<const Reference2DElementInfo &>(ReferenceElementInfo::getInfo(tileType));
				tileNormal = referenceInfo.evalNormal(tileCoordinates.data(), point);
			} else if (dimension == 1) {
				const Reference1DElementInfo &referenceInfo = static_cast<const Reference1DElementInfo &>(ReferenceElementInfo::getInfo(tileType));
				tileNormal = referenceInfo.evalNormal(tileCoordinates.data(), orientation, point);
			} else if (dimension == 0) {
				tileNormal = orientation;
			}

			normal      += tileArea * tileNormal;
			surfaceArea += tileArea;
		}
		normal = (1. / surfaceArea) * normal;

		return normal;
	}

	default:
	{
		assert(getDimension() != 3);

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
	Evaluates the distance between the element and the specified point.

	\param[in] point is the point
	\param coordinates are the coordinate of the vertices
	\result The distance between the element and the specified point.
*/
double Element::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *coordinates) const
{
	switch (m_type) {

	case ElementType::POLYGON:
	{
		return CGElem::distancePointPolygon(point, getVertexCount(), coordinates);
	}

	case ElementType::POLYHEDRON:
	{
		double distance = std::numeric_limits<double>::max();

		int nFaces = getFaceCount();
		std::vector<std::array<double, 3>> faceCoordinates;
		for (int i = 0; i < nFaces; ++i) {
			ElementType faceType = getFaceType(i);
			ConstProxyVector<int> faceVertexIds = getFaceLocalVertexIds(i);
			int nFaceVertices = faceVertexIds.size();
			faceCoordinates.resize(nFaceVertices);
			for (int k = 0; k < nFaceVertices; ++k) {
				faceCoordinates[k] = coordinates[faceVertexIds[k]];
			}

			double faceDistance;
			bool faceHasReferenceInfo = ReferenceElementInfo::hasInfo(faceType);
			if (faceHasReferenceInfo) {
				faceDistance = ReferenceElementInfo::getInfo(faceType).evalPointDistance(point, faceCoordinates.data());
			} else {
				faceDistance = CGElem::distancePointPolygon(point, nFaceVertices, faceCoordinates.data());
			}

			distance = std::min(faceDistance, distance);
		}

		return distance;
	}

	default:
	{
		assert(ReferenceElementInfo::hasInfo(m_type));

		return getInfo().evalPointDistance(point, coordinates);
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
			ConstProxyVector<int> localVertexIds = getFaceLocalVertexIds(i);
			int nFaceVertices = localVertexIds.size();
			faceTesselationIds[i].resize(nFaceVertices);
			for (int k = 0; k < nFaceVertices; ++k) {
				faceTesselationIds[i][k] = vertexTesselationIds[localVertexIds[k]];
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
	Gets the size of the face stream that describes the element.

	\result The size of the face stream that describes the element.
*/
int Element::getFaceStreamSize() const
{
	int nFaces = getFaceCount();
	int size   = 1 + (getFaceStreamPosition(nFaces) - 1);

	return size;
}

/*!
	Gets the face stream that describes the element.

	\result The face stream that describes the element.
*/
std::vector<long> Element::getFaceStream() const
{
	int nFaces = getFaceCount();
	int faceStreamSize = getFaceStreamSize();
	std::vector<long> faceStream(faceStreamSize);

	int pos = 0;
	faceStream[pos] = nFaces;
	for (int i = 0; i < nFaces; ++i) {
		ConstProxyVector<long> faceVertexIds = getFaceVertexIds(i);
		int nFaceVertices = faceVertexIds.size();

		++pos;
		faceStream[pos] = getFaceVertexCount(i);
		for (int k = 0; k < nFaceVertices; ++k) {
			++pos;
			faceStream[pos] = faceVertexIds[k];
		}
	}

	return faceStream;
}

/*!
	Renumber the vertices of the specified face stream according to the
	given map.

	\param map is the vertex map
	\param faceStream is the face stream to be renumbered
*/
void Element::renumberFaceStream(const PiercedStorage<long, long> &map, std::vector<long> *faceStream)
{
	int pos = 0;
	int nFaces = (*faceStream)[pos];
	for (int i = 0; i < nFaces; ++i) {
		++pos;
		int nFaceVertices = (*faceStream)[pos];
		for (int k = 0; k < nFaceVertices; ++k) {
			++pos;
			(*faceStream)[pos] = map.at((*faceStream)[pos]);
		}
	}
}

/*!
	Gets the position of the specified face in the face stream.

	\param face is the face
	\result The position of the specified face in the face stream.
*/
int Element::getFaceStreamPosition(int face) const
{
	const long *connectivity = getConnect();

	int position = 1;
	for (int i = 0; i < face; ++i) {
		position += 1 + connectivity[position];
	}

	return position;
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

	assert((int) edgeSet.size() == nRequestedEdges);

	return edgeConnects;
}

/*!
        Returns the buffer size required to communicate cell data

        \result buffer size (in bytes)
*/
unsigned int Element::getBinarySize() const
{
	unsigned int binarySize = sizeof(m_type) + sizeof(m_id) + getConnectSize() * sizeof(long) + sizeof(m_pid);
	if (bitpit::ReferenceElementInfo::hasInfo(m_type)) {
		binarySize += sizeof(int);
	}

	return binarySize;
}

/*!
	\class ElementHalfFace
	\ingroup patchelements

	\brief The ElementHalfFace class defines element half-faces.

	ElementHalfFace is the class that defines element half-faces. Each face
	can be seen as two half-faces: one belonging to an element and the other
	belonging to the neighbouring element. A half-face is identify by its
	vertices and by the winding order of the vertices.
*/

/*!
	Constructor.

	\param element is a reference to the element the owns the face
	\param face if the local face of the element
	\param winding is the winding order of the vertices
*/
ElementHalfFace::ElementHalfFace(Element &element, int face, Winding winding)
    : m_element(element), m_face(face),
      m_connect(m_element.getFaceConnect(m_face)),
      m_winding(winding)
{
	// Find the vertex with the lowest id, this will be used as first vertex
	// when iterating over the connectivity.
	std::size_t connectSize = m_connect.size();

	m_connectBegin = 0;
	if (connectSize > 0) {
		long smallestVertexId = m_connect[m_connectBegin];
		for (std::size_t i = 1; i < connectSize; ++i) {
			long vertexId = m_connect[i];
			if (vertexId < smallestVertexId) {
				m_connectBegin = i;
				smallestVertexId = vertexId;
			}
		}
	}
}

/*!
	Get the element the face belongs to.

	\result Returns the element the face belongs to.
*/
Element & ElementHalfFace::getElement() const
{
    return m_element;
}

/*!
	Get the local face index.

	\result Returns the local face index.
*/
int ElementHalfFace::getFace() const
{
    return m_face;
}

/*!
	Get the connectivity of the face

	\result Returns the connectivity of the face.
*/
const ConstProxyVector<long> & ElementHalfFace::getConnect() const
{
    return m_connect;
}

/*!
	Get vertex winding order.

	\result Vertex winding order.
*/
ElementHalfFace::Winding ElementHalfFace::getWinding() const
{
	return m_winding;
}

/*!
	Set vertex winding order.

	\param winding is the vertex winding order
*/
void ElementHalfFace::setWinding(Winding winding)
{
	m_winding = winding;
}

/*!
	Comparison operator for the ElementHalfFace class.

	\param other is the object to be compared with
	\result Returns true if the two half-faces define the same half-face,
	false otherwise.
*/
bool ElementHalfFace::operator==(const bitpit::ElementHalfFace &other) const
{
	const bitpit::ConstProxyVector<long> &connect1 = m_connect;
	const bitpit::ConstProxyVector<long> &connect2 = other.m_connect;

	std::size_t connectSize = connect1.size();
	if (connectSize != connect2.size()) {
		return false;
	}

	int winding = m_winding * other.m_winding;
	std::size_t offset = other.m_connectBegin - winding * m_connectBegin + 2 * connectSize;
	for (std::size_t i = 0; i < connectSize; ++i) {
		std::size_t k1 = i;
		std::size_t k2 = (offset + winding * i) % connectSize;
		if (connect1[k1] != connect2[k2]) {
			return false;
		}
	}

	return true;
}

/*!
	Comparison operator for the ElementHalfFace class.

	\param other is the object to be compared with
	\result Returns true if the two half-faces define different half-faces,
	false otherwise.
*/
bool ElementHalfFace::operator!=(const bitpit::ElementHalfFace &other) const
{
	return !((*this) == other);
}

/*!
	\class ElementHalfFace::Hasher
	\ingroup patchelements

	\brief The ElementHalfFace::Hasher class allows to create hashes for the
	face items.
*/

/*!
	Generete the hash for the specified face item.

	\param item is the face item
	\result The hash of the specified face item.
*/
std::size_t ElementHalfFace::Hasher::operator()(const ElementHalfFace &item) const
{
	const ConstProxyVector<long> &connectivity = item.m_connect;
	std::size_t connectivitySize = connectivity.size();

	std::size_t hash = connectivitySize;
	if (item.m_winding == ElementHalfFace::WINDING_NATURAL) {
		for (std::size_t i = 0; i < connectivitySize; ++i) {
			std::size_t k = (item.m_connectBegin + i) % connectivitySize;
			utils::hashing::hash_combine(hash, connectivity[k]);
		}
	} else {
		for (std::size_t i = connectivitySize; i > 0; --i) {
			std::size_t k = (item.m_connectBegin + i) % connectivitySize;
			utils::hashing::hash_combine(hash, connectivity[k]);
		}
	}

	return hash;
}

// Explicit instantiation of the Element containers
template class PiercedVector<Element>;

}
