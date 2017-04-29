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

#include "bitpit_common.hpp"

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

	\param id the id of the element
	\param type is the type of the element
*/
Element::Element(long id, ElementType type)
{
	_initialize(id, type);
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
*/
void Element::initialize(long id, ElementType type)
{
	_initialize(id, type);
}

/*!
	Internal function to initialize the data structures of the element.

	\param id the id of the element
	\param type the type of the element
*/
void Element::_initialize(long id, ElementType type)
{
	// Set the id
	setId(id);

	// Get previous connect size
	int previousConnectSize;
	if (type != ElementType::UNDEFINED) {
		if (m_connect) {
			assert(getType() != ElementType::UNDEFINED);
			previousConnectSize = getInfo().nVertices;
		} else {
			previousConnectSize = 0;
		}
	}

	// Set element type
	setType(type);

	// Set connectivity
	if (type != ElementType::UNDEFINED) {
		int connectSize = getInfo().nVertices;
		if (connectSize != previousConnectSize) {
			setConnect(std::unique_ptr<long[]>(new long[connectSize]));
		}
	} else {
		unsetConnect();
	}
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
	int nVertices = getVertexCount();
	const long *connectivity = getConnect();

	int localVertexId = std::find(connectivity, connectivity + nVertices, vertexId) - connectivity;
	if (localVertexId >= nVertices) {
		return -1;
	}

	return localVertexId;
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
const std::vector<int> & Element::getFaceLocalConnect(const int &face) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");

	default:
		return getInfo().faceConnect[face];

	}
}

/*!
	Gets the connectivity of the specified face of the element.

	\param face is the face for which the connectivity is reqested
	\result The connectivity of the specified face of the element.
*/
std::vector<long> Element::getFaceConnect(int face) const
{
	const std::vector<int> &localFaceConnect = getFaceLocalConnect(face);
	int nFaceVertices = localFaceConnect.size();

	std::vector<long> faceConnect(nFaceVertices);
	for (int k = 0; k < nFaceVertices; ++k) {
		int localVertexId = localFaceConnect[k];
		long vertexId = getVertex(localVertexId);
		faceConnect[k] = vertexId;
	}

	return faceConnect;
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
const std::vector<int> & Element::getEdgeLocalConnect(const int &edge) const
{
	switch (m_type) {

	case (ElementType::POLYGON):
	case (ElementType::POLYHEDRON):
	case (ElementType::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");

	default:
		return getInfo().edgeConnect[edge];

	}
}

/*!
	Gets the connectivity of the specified edge of the element.

	\param edge is the edge for which the connectivity is reqested
	\result The connectivity of the specified edge of the element.
*/
std::vector<long> Element::getEdgeConnect(int edge) const
{
	const std::vector<int> &localEdgeConnect = getEdgeLocalConnect(edge);
	int nEdgeVertices = localEdgeConnect.size();

	std::vector<long> edgeConnect(nEdgeVertices);
	for (int k = 0; k < nEdgeVertices; ++k) {
		int localVertexId = localEdgeConnect[k];
		long vertexId = getVertex(localVertexId);
		edgeConnect[k] = vertexId;
	}

	return edgeConnect;
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
	Sets the vertex with the specified local index.

	\param index is the local index of the vertex
	\param vertex is the id of the vertex.
*/
void Element::setVertex(const int &index, const long &vertex)
{
	m_connect[index] = vertex;
}

/*!
	Gets the vertex with the specified local index.

	\param vertex is the local index of the vertex
	\result The id of the specified vertex.
*/
long Element::getVertex(const int &vertex) const
{
	return m_connect[vertex];
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
