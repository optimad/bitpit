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

#include "element.hpp"

#include <assert.h>
#include <limits>

/*!
	Input stream operator for class Element

	\param[in] buffer is the input stream
	\param[in] element is the element to be streamed
	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::Element &element)
{
	buffer >> element.m_type;
	element.initialize(element.m_type);
	int nVertices = element.get_vertex_count();
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
	int nVertices = element.get_vertex_count();
	buffer << element.get_type();
	for (int i = 0; i < nVertices; ++i) {
	    buffer << element.m_connect[i];
	}

	return buffer;
}

namespace bitpit {

/*!
	\ingroup patch
	@{
*/

/*!
	\class ElementInfo

	\brief The ElementInfo class allow to define elements information.

	Element is a struct that hold the basic geometrical information of
	an element.

	The local numbering scheme of element vertices is shown below.

	\image html common_elements.png
*/

/*!
	\enum ElementInfo::Type

	This enum defines the element types that can be used.

	\var ElementInfo::Type ElementInfo::UNDEFINED
	An undefined element.

	\var ElementInfo::Type ElementInfo::VERTEX
	A point.

	\var ElementInfo::Type ElementInfo::LINE
	A line.

	\var ElementInfo::Type ElementInfo::TRIANGLE
	A triangle.

	\var ElementInfo::Type ElementInfo::PIXEL
	A rectangle.

	\var ElementInfo::Type ElementInfo::QUAD
	A quadrangle.

	\var ElementInfo::Type ElementInfo::POLYGON
	A polygon.

	\var ElementInfo::Type ElementInfo::TETRA
	A tetrahedron.

	\var ElementInfo::Type ElementInfo::VOXEL
	A brick.

	\var ElementInfo::Type ElementInfo::HEXAHEDRON
	A hexahedron.

	\var ElementInfo::Type ElementInfo::PYRAMID
	A pyramid.

	\var ElementInfo::Type ElementInfo::WEDGE
	A wedge (triangular prism).

	\var ElementInfo::Type ElementInfo::POLYHEDRON
	A polyhedron.
*/

const ElementInfo ElementInfo::undefinedInfo  = ElementInfo(ElementInfo::UNDEFINED);
const ElementInfo ElementInfo::vertexInfo     = ElementInfo(ElementInfo::VERTEX);
const ElementInfo ElementInfo::lineInfo       = ElementInfo(ElementInfo::LINE);
const ElementInfo ElementInfo::triangleInfo   = ElementInfo(ElementInfo::TRIANGLE);
const ElementInfo ElementInfo::pixelInfo      = ElementInfo(ElementInfo::PIXEL);
const ElementInfo ElementInfo::quadInfo       = ElementInfo(ElementInfo::QUAD);
const ElementInfo ElementInfo::tetraInfo      = ElementInfo(ElementInfo::TETRA);
const ElementInfo ElementInfo::voxelInfo      = ElementInfo(ElementInfo::VOXEL);
const ElementInfo ElementInfo::hexahedronInfo = ElementInfo(ElementInfo::HEXAHEDRON);
const ElementInfo ElementInfo::pyramidInfo    = ElementInfo(ElementInfo::PYRAMID);
const ElementInfo ElementInfo::wedgeInfo      = ElementInfo(ElementInfo::WEDGE);

/*!
	Default constructor
*/
ElementInfo::ElementInfo()
	: type(UNDEFINED)
{
}

/*!
	Creates a new set of element information.

	\param type is the type of element
*/
ElementInfo::ElementInfo(ElementInfo::Type type)
{
	switch (type) {

	case (VERTEX):
		initializeVertexInfo();
		break;

	case (LINE):
		initializeLineInfo();
		break;

	case (TRIANGLE):
		initializeTriangleInfo();
		break;

	case (PIXEL):
		initializePixelInfo();
		break;

	case (QUAD):
		initializeQuadInfo();
		break;

	case (TETRA):
		initializeTetraInfo();
		break;

	case (VOXEL):
		initializeVoxelInfo();
		break;

	case (HEXAHEDRON):
		initializeHexahedronInfo();
		break;

	case (PYRAMID):
		initializePyramidInfo();
		break;

	case (WEDGE):
		initializeWedgeInfo();
		break;

	default:
		initializeUndefinedInfo();
		break;

	}
}

/*!
	Gets the information for the specified element type.

	\param type is the type of element
	\result The information for the specified element type.
*/
const ElementInfo & ElementInfo::get_element_info(ElementInfo::Type type)
{
	switch (type) {

	case (VERTEX):
		return vertexInfo;

	case (LINE):
		return lineInfo;

	case (TRIANGLE):
		return triangleInfo;

	case (PIXEL):
		return pixelInfo;

	case (QUAD):
		return quadInfo;

	case (TETRA):
		return tetraInfo;

	case (VOXEL):
		return voxelInfo;

	case (HEXAHEDRON):
		return hexahedronInfo;

	case (PYRAMID):
		return pyramidInfo;

	case (WEDGE):
		return wedgeInfo;

	default:
		assert(false);
		return undefinedInfo;

	}
}

/*!
	Initializes the information for the undefined element.
*/
void ElementInfo::initializeUndefinedInfo()
{
	type      = UNDEFINED;
	dimension = -1;

	nVertices = -1;
	nEdges    = -1;
	nFaces    = -1;
}

/*!
	Initializes the information for the point element.
*/
void ElementInfo::initializeVertexInfo()
{
	type      = VERTEX;
	dimension = 0;

	// Vertices data
	nVertices = 1;

	// Edge data
	nEdges = 1;

	edge_type = std::vector<Type>(nEdges);
	edge_type[0] = VERTEX;

	edge_connect = std::vector<std::vector<int>>(nEdges);
	edge_connect[0] = std::vector<int>(nVertices);
	edge_connect[0][0] = 0;

	// Face data
	nFaces = 1;

	face_type = std::vector<Type>(nFaces);
	face_type[0] = VERTEX;

	face_connect = std::vector<std::vector<int>>(nFaces);
	face_connect[0] = std::vector<int>(nVertices);
	face_connect[0][0] = 0;
}

/*!
	Initializes the information for the line element.
*/
void ElementInfo::initializeLineInfo()
{
	ElementInfo vertexInfo(VERTEX);

	type      = LINE;
	dimension = 1;

	// Vertices data
	nVertices = 2;

	// Edge data
	nEdges = 2;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edge_connect[k]    = std::vector<int>(vertexInfo.nVertices);
		edge_connect[k][0] = k;
	}

	// Face data
	nFaces = 2;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = VERTEX;
		face_connect[k]    = std::vector<int>(vertexInfo.nVertices);
		face_connect[k][0] = k;
	}
}

/*!
	Initializes the information for the triangle element.
*/
void ElementInfo::initializeTriangleInfo()
{
	ElementInfo vertexInfo(VERTEX);
	ElementInfo lineInfo(LINE);

	type      = TRIANGLE;
	dimension = 2;

	// Vertices data
	nVertices = 3;

	// Edge data
	nEdges = 3;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edge_connect[k]    = std::vector<int>(vertexInfo.nVertices);
		edge_connect[k][0] = k;
	}

	// Face data
	nFaces = 3;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = LINE;
		face_connect[k]    = std::vector<int>(lineInfo.nVertices);
		face_connect[k][0] = k;
		face_connect[k][1] = (k + 1) % nVertices;
	}
}

/*!
	Initializes the information for the quadrangle element.
*/
void ElementInfo::initializePixelInfo()
{
	ElementInfo vertexInfo(VERTEX);
	ElementInfo lineInfo(LINE);

	type      = QUAD;
	dimension = 2;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 4;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edge_connect[k]    = std::vector<int>(vertexInfo.nVertices);
		edge_connect[k][0] = k;
	}

	// Face data
	nFaces = 4;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = LINE;
		face_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	face_connect[0][0] = 2;
	face_connect[0][1] = 0;

	face_connect[1][0] = 1;
	face_connect[1][1] = 3;

	face_connect[2][0] = 0;
	face_connect[2][1] = 1;

	face_connect[3][0] = 3;
	face_connect[3][1] = 2;
}

/*!
	Initializes the information for the quadrangle element.
*/
void ElementInfo::initializeQuadInfo()
{
	ElementInfo vertexInfo(VERTEX);
	ElementInfo lineInfo(LINE);

	type      = QUAD;
	dimension = 2;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 4;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edge_connect[k]    = std::vector<int>(vertexInfo.nVertices);
		edge_connect[k][0] = k;
	}

	// Face data
	nFaces = 4;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = LINE;
		face_connect[k]    = std::vector<int>(lineInfo.nVertices);
		face_connect[k][0] = k;
		face_connect[k][1] = (k + 1) % nVertices;
	}
}

/*!
	Initializes the information for the tetrahedron element.
*/
void ElementInfo::initializeTetraInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo triangleInfo(TRIANGLE);

	type      = TETRA;
	dimension = 3;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 6;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edge_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edge_connect[0][0] = 0;
	edge_connect[0][1] = 1;

	edge_connect[1][0] = 1;
	edge_connect[1][1] = 2;

	edge_connect[2][0] = 2;
	edge_connect[2][1] = 0;

	edge_connect[3][0] = 3;
	edge_connect[3][1] = 0;

	edge_connect[4][0] = 3;
	edge_connect[4][1] = 1;

	edge_connect[5][0] = 3;
	edge_connect[5][1] = 2;

	// Face data
	nFaces = 4;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = TRIANGLE;
		face_connect[k] = std::vector<int>(triangleInfo.nVertices);
	}

	face_connect[0][0] = 1;
	face_connect[0][1] = 0;
	face_connect[0][2] = 2;

	face_connect[1][0] = 0;
	face_connect[1][1] = 3;
	face_connect[1][2] = 2;

	face_connect[2][0] = 3;
	face_connect[2][1] = 1;
	face_connect[2][2] = 2;

	face_connect[3][0] = 0;
	face_connect[3][1] = 1;
	face_connect[3][2] = 3;
}

/*!
	Initializes the information for the brick element.
*/
void ElementInfo::initializeVoxelInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo pixelInfo(PIXEL);

	type      = VOXEL;
	dimension = 3;

	// Vertices data
	nVertices = 8;

	// Edge data
	nEdges = 12;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edge_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edge_connect[0][0] = 1;
	edge_connect[0][1] = 0;

	edge_connect[1][0] = 1;
	edge_connect[1][1] = 2;

	edge_connect[2][0] = 2;
	edge_connect[2][1] = 3;

	edge_connect[3][0] = 3;
	edge_connect[3][1] = 0;

	edge_connect[4][0] = 4;
	edge_connect[4][1] = 5;

	edge_connect[5][0] = 5;
	edge_connect[5][1] = 6;

	edge_connect[6][0] = 6;
	edge_connect[6][1] = 7;

	edge_connect[7][0] = 7;
	edge_connect[7][1] = 4;

	edge_connect[8][0] = 0;
	edge_connect[8][1] = 4;

	edge_connect[9][0] = 1;
	edge_connect[9][1] = 5;

	edge_connect[10][0] = 2;
	edge_connect[10][1] = 6;

	edge_connect[11][0] = 3;
	edge_connect[11][1] = 7;

	// Face data
	nFaces = 6;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = PIXEL;
		face_connect[k] = std::vector<int>(pixelInfo.nVertices);
	}

	face_connect[0][0] = 2;
	face_connect[0][1] = 0;
	face_connect[0][2] = 4;
	face_connect[0][3] = 6;

	face_connect[1][0] = 1;
	face_connect[1][1] = 3;
	face_connect[1][2] = 7;
	face_connect[1][3] = 5;

	face_connect[2][0] = 0;
	face_connect[2][1] = 1;
	face_connect[2][2] = 5;
	face_connect[2][3] = 4;

	face_connect[3][0] = 3;
	face_connect[3][1] = 2;
	face_connect[3][2] = 6;
	face_connect[3][3] = 7;

	face_connect[4][0] = 2;
	face_connect[4][1] = 3;
	face_connect[4][2] = 1;
	face_connect[4][3] = 0;

	face_connect[5][0] = 4;
	face_connect[5][1] = 5;
	face_connect[5][2] = 7;
	face_connect[5][3] = 6;
}

/*!
	Initializes the information for the hexahedron element.
*/
void ElementInfo::initializeHexahedronInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo quadInfo(QUAD);

	type      = HEXAHEDRON;
	dimension = 3;

	// Vertices data
	nVertices = 8;

	// Edge data
	nEdges = 12;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edge_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edge_connect[0][0] = 1;
	edge_connect[0][1] = 0;

	edge_connect[1][0] = 1;
	edge_connect[1][1] = 2;

	edge_connect[2][0] = 2;
	edge_connect[2][1] = 3;

	edge_connect[3][0] = 3;
	edge_connect[3][1] = 0;

	edge_connect[4][0] = 4;
	edge_connect[4][1] = 5;

	edge_connect[5][0] = 5;
	edge_connect[5][1] = 6;

	edge_connect[6][0] = 6;
	edge_connect[6][1] = 7;

	edge_connect[7][0] = 7;
	edge_connect[7][1] = 4;

	edge_connect[8][0] = 0;
	edge_connect[8][1] = 4;

	edge_connect[9][0] = 1;
	edge_connect[9][1] = 5;

	edge_connect[10][0] = 2;
	edge_connect[10][1] = 6;

	edge_connect[11][0] = 3;
	edge_connect[11][1] = 7;

	// Face data
	nFaces = 6;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = QUAD;
		face_connect[k] = std::vector<int>(quadInfo.nVertices);
	}

	face_connect[0][0] = 1;
	face_connect[0][1] = 0;
	face_connect[0][2] = 3;
	face_connect[0][3] = 2;

	face_connect[1][0] = 4;
	face_connect[1][1] = 5;
	face_connect[1][2] = 6;
	face_connect[1][3] = 7;

	face_connect[2][0] = 7;
	face_connect[2][1] = 3;
	face_connect[2][2] = 0;
	face_connect[2][3] = 4;

	face_connect[3][0] = 5;
	face_connect[3][1] = 1;
	face_connect[3][2] = 2;
	face_connect[3][3] = 6;

	face_connect[4][0] = 4;
	face_connect[4][1] = 0;
	face_connect[4][2] = 1;
	face_connect[4][3] = 5;

	face_connect[5][0] = 6;
	face_connect[5][1] = 2;
	face_connect[5][2] = 3;
	face_connect[5][3] = 7;
}

/*!
	Initializes the information for the pyramid element.
*/
void ElementInfo::initializePyramidInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo triangleInfo(TRIANGLE);
	ElementInfo quadInfo(QUAD);

	type      = PYRAMID;
	dimension = 3;

	// Vertices data
	nVertices = 5;

	// Edge data
	nEdges = 8;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edge_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edge_connect[0][0] = 0;
	edge_connect[0][1] = 1;

	edge_connect[1][0] = 1;
	edge_connect[1][1] = 2;

	edge_connect[2][0] = 2;
	edge_connect[2][1] = 3;

	edge_connect[3][0] = 3;
	edge_connect[3][1] = 0;

	edge_connect[4][0] = 4;
	edge_connect[4][1] = 0;

	edge_connect[5][0] = 4;
	edge_connect[5][1] = 1;

	edge_connect[6][0] = 4;
	edge_connect[6][1] = 2;

	edge_connect[7][0] = 4;
	edge_connect[7][1] = 3;

	// Face data
	nFaces = 5;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		if (k == 0) {
			face_type[k]    = QUAD;
			face_connect[k] = std::vector<int>(quadInfo.nVertices);
		} else {
			face_type[k]    = TRIANGLE;
			face_connect[k] = std::vector<int>(triangleInfo.nVertices);
		}
	}

	face_connect[0][0] = 0;
	face_connect[0][1] = 3;
	face_connect[0][2] = 2;
	face_connect[0][3] = 1;

	face_connect[1][0] = 3;
	face_connect[1][1] = 0;
	face_connect[1][2] = 4;

	face_connect[2][0] = 0;
	face_connect[2][1] = 1;
	face_connect[2][2] = 4;

	face_connect[3][0] = 1;
	face_connect[3][1] = 2;
	face_connect[3][2] = 4;

	face_connect[4][0] = 2;
	face_connect[4][1] = 3;
	face_connect[4][2] = 4;
}

/*!
	Initializes the information for the wedge (triangular prism) element.
*/
void ElementInfo::initializeWedgeInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo triangleInfo(TRIANGLE);
	ElementInfo quadInfo(QUAD);

	type      = WEDGE;
	dimension = 3;

	// Vertices data
	nVertices = 6;

	// Edge data
	nEdges = 9;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edge_connect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edge_connect[0][0] = 1;
	edge_connect[0][1] = 0;

	edge_connect[1][0] = 1;
	edge_connect[1][1] = 2;

	edge_connect[2][0] = 2;
	edge_connect[2][1] = 0;

	edge_connect[3][0] = 3;
	edge_connect[3][1] = 4;

	edge_connect[4][0] = 4;
	edge_connect[4][1] = 5;

	edge_connect[5][0] = 5;
	edge_connect[5][1] = 3;

	edge_connect[6][0] = 3;
	edge_connect[6][1] = 0;

	edge_connect[7][0] = 4;
	edge_connect[7][1] = 1;

	edge_connect[8][0] = 5;
	edge_connect[8][1] = 2;

	// Face data
	nFaces = 5;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		if (k == 0 || k == 1) {
			face_type[k]    = TRIANGLE;
			face_connect[k] = std::vector<int>(triangleInfo.nVertices);
		} else {
			face_type[k]    = QUAD;
			face_connect[k] = std::vector<int>(quadInfo.nVertices);
		}
	}

	face_connect[0][0] = 1;
	face_connect[0][1] = 0;
	face_connect[0][2] = 2;

	face_connect[1][0] = 3;
	face_connect[1][1] = 4;
	face_connect[1][2] = 5;

	face_connect[2][0] = 3;
	face_connect[2][1] = 0;
	face_connect[2][2] = 1;
	face_connect[2][3] = 4;

	face_connect[3][0] = 4;
	face_connect[3][1] = 1;
	face_connect[3][2] = 2;
	face_connect[3][3] = 5;

	face_connect[4][0] = 5;
	face_connect[4][1] = 2;
	face_connect[4][2] = 0;
	face_connect[4][3] = 3;
}


/*!
	\class Element

	\brief The Element class provides an interface for defining elements.

	Element is the base class for defining elements like cells and
	intefaces.
*/

const long Element::NULL_ELEMENT_ID = std::numeric_limits<long>::min();

/*!
	Default constructor.
*/
Element::Element()
{
	initialize(ElementInfo::UNDEFINED);

	set_id(NULL_ELEMENT_ID);
}

/*!
	Creates a new element.
*/
Element::Element(const long &id, ElementInfo::Type type)
{
	initialize(type);

	set_id(id);
}

/*!
	Initializes the data structures of the element.

	\param type the type of the element
*/
void Element::initialize(ElementInfo::Type type)
{
	set_type(type);

	if (get_type() != ElementInfo::UNDEFINED) {
		const int &nVertices = get_info().nVertices;
		set_connect(std::unique_ptr<long[]>(new long[nVertices]));
	} else {
		unset_connect();
	}
}

/*!
	Sets the ID of the element.

	\param id the ID of the element
*/
void Element::set_id(const long &id)
{
	m_id = id;
}

/*!
	Gets the ID of the element.

	\return The ID of the element
*/
long Element::get_id() const
{
	return m_id;
}

/*!
	Gets the basic information of the element.

	\result A constant reference to the basic information of the element.
*/
const ElementInfo & Element::get_info() const
{
	return ElementInfo::get_element_info(m_type);
}

/*!
	Sets the element type.

	\param type the element type
*/
void Element::set_type(ElementInfo::Type type)
{
	m_type = type;
}

/*!
	Gets the element type.

	\result The element type
*/
ElementInfo::Type Element::get_type() const
{
	return m_type;
}

/*!
	Sets the vertex connectivity of the element.

	\param connect a pointer to the connectivity of the element
*/
void Element::set_connect(std::unique_ptr<long[]> connect)
{
	m_connect = std::move(connect);
}

/*!
	Unsets the vertex connectivity of the element.
*/
void Element::unset_connect()
{
	m_connect.reset();
}

/*!
	Gets the vertex connectivity of the element.

	\result A constant pointer to the connectivity of the element
*/
const long * Element::get_connect() const
{
	return m_connect.get();
}

/*!
	Gets the vertex connectivity of the element.

	\result A pointer to the connectivity of the element
*/
long * Element::get_connect()
{
	return m_connect.get();
}

/*!
	Gets the number of faces of the element.

	\result The number of vertices of the element
*/
int Element::get_face_count() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return -1;

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.nFaces;

	}
}

/*!
	Gets the face type of the specified face of the element.

	\result The face type of specified face of the element
*/
ElementInfo::Type Element::get_face_type(const int &face) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return ElementInfo::UNDEFINED;

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.face_type[face];

	}
}

/*!
	Gets the local connectivity of the specified face of the element.

	\param face is the face for which the connectiviy is reqested
	\result The local connectivity of the specified face of the element.
*/
std::vector<int> Element::get_face_local_connect(const int &face) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return std::vector<int>();

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.face_connect[face];

	}
}

/*!
	Gets the number of edges of the element.

	\result The number of edges of the element
*/
int Element::get_edge_count() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return -1;

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.nEdges;

	}
}

/*!
	Gets the local connectivity of the specified edge of the element.

	\param edge is the edge for which the connectiviy is reqested
	\result The local connectivity of the specified edge of the element.
*/
std::vector<int> Element::get_edge_local_connect(const int &edge) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return std::vector<int>();

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.edge_connect[edge];

	}
}

/*!
	Gets the dimension of the element.

	\return The dimension of the element
*/
int Element::get_dimension() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
		return 2;

	case (ElementInfo::POLYHEDRON):
		return 3;

	case (ElementInfo::UNDEFINED):
		assert(false);
		return -1;

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.dimension;

	}
}

/*!
	Returns true if the element is a three-dimensional element.

	\return Returns true if the element is a three-dimensional element,
	false otherwise.
*/
bool Element::is_three_dimensional() const
{
	return (get_dimension() == 3);
}

/*!
	Gets the number of vertices of the element.

	\result The number of vertices of the element
*/
int Element::get_vertex_count() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		assert(false);
		return -1;

	default:
		const ElementInfo &elementInfo = ElementInfo::get_element_info(m_type);
		return elementInfo.nVertices;

	}
}

/*!
	Sets the vertex with the specified local index.

	\param index is the local index of the vertex
	\param vertex is the id of the vertex.
*/
void Element::set_vertex(const int &index, const long &vertex)
{
	m_connect[index] = vertex;
}

/*!
	Gets the vertex with the specified local index.

	\param vertex is the local index of the vertex
	\result The id of the specified vertex.
*/
long Element::get_vertex(const int &vertex) const
{
	return m_connect[vertex];
}

/*!
        Returns the buffer size required to communicate cell data

        \result buffer size (in bytes)
*/
unsigned int Element::get_binary_size()
{
    return (sizeof(ElementInfo::Type) + get_vertex_count() * sizeof(long));
}

/*!
	@}
*/

}
