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
	buffer >> element.m_type;
	buffer >> element.m_id;
	element._initialize(element.m_type);
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
	\class ElementInfo
	\ingroup patchelements

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
const ElementInfo & ElementInfo::getElementInfo(ElementInfo::Type type)
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
		BITPIT_UNREACHABLE("Unsupported element");
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

	edgeConnect = std::vector<std::vector<int>>(nEdges);
	edgeConnect[0] = std::vector<int>(nVertices);
	edgeConnect[0][0] = 0;

	// Face data
	nFaces = 1;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	face_type[0] = VERTEX;

	facesInfo[0] = this;

	faceConnect = std::vector<std::vector<int>>(nFaces);
	faceConnect[0] = std::vector<int>(nVertices);
	faceConnect[0][0] = 0;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
		edgeConnect[k][0] = k;
	}

	// Face data
	nFaces = 2;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = VERTEX;
		facesInfo[k]      = &vertexInfo;
		faceConnect[k]    = std::vector<int>(vertexInfo.nVertices);
		faceConnect[k][0] = k;
	}

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
		edgeConnect[k][0] = k;
	}

	// Face data
	nFaces = 3;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = LINE;
		facesInfo[k]      = &lineInfo;
		faceConnect[k]    = std::vector<int>(lineInfo.nVertices);
		faceConnect[k][0] = k;
		faceConnect[k][1] = (k + 1) % nVertices;
	}

	initializeFaceEdges(facesInfo);
}

/*!
	Initializes the information for the quadrangle element.
*/
void ElementInfo::initializePixelInfo()
{
	ElementInfo vertexInfo(VERTEX);
	ElementInfo lineInfo(LINE);

	type      = PIXEL;
	dimension = 2;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 4;

	edge_type = std::vector<Type>(nEdges);
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
		edgeConnect[k][0] = k;
	}

	// Face data
	nFaces = 4;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = LINE;
		facesInfo[k]   = &lineInfo;
		faceConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	faceConnect[0][0] = 2;
	faceConnect[0][1] = 0;

	faceConnect[1][0] = 1;
	faceConnect[1][1] = 3;

	faceConnect[2][0] = 0;
	faceConnect[2][1] = 1;

	faceConnect[3][0] = 3;
	faceConnect[3][1] = 2;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = VERTEX;
		edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
		edgeConnect[k][0] = k;
	}

	// Face data
	nFaces = 4;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = LINE;
		facesInfo[k]      = &lineInfo;
		faceConnect[k]    = std::vector<int>(lineInfo.nVertices);
		faceConnect[k][0] = k;
		faceConnect[k][1] = (k + 1) % nVertices;
	}

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edgeConnect[0][0] = 0;
	edgeConnect[0][1] = 1;

	edgeConnect[1][0] = 1;
	edgeConnect[1][1] = 2;

	edgeConnect[2][0] = 2;
	edgeConnect[2][1] = 0;

	edgeConnect[3][0] = 3;
	edgeConnect[3][1] = 0;

	edgeConnect[4][0] = 3;
	edgeConnect[4][1] = 1;

	edgeConnect[5][0] = 3;
	edgeConnect[5][1] = 2;

	// Face data
	nFaces = 4;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = TRIANGLE;
		facesInfo[k]   = &triangleInfo;
		faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
	}

	faceConnect[0][0] = 1;
	faceConnect[0][1] = 0;
	faceConnect[0][2] = 2;

	faceConnect[1][0] = 0;
	faceConnect[1][1] = 3;
	faceConnect[1][2] = 2;

	faceConnect[2][0] = 3;
	faceConnect[2][1] = 1;
	faceConnect[2][2] = 2;

	faceConnect[3][0] = 0;
	faceConnect[3][1] = 1;
	faceConnect[3][2] = 3;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edgeConnect[0][0] = 0;
	edgeConnect[0][1] = 2;

	edgeConnect[1][0] = 1;
	edgeConnect[1][1] = 3;

	edgeConnect[2][0] = 0;
	edgeConnect[2][1] = 1;

	edgeConnect[3][0] = 2;
	edgeConnect[3][1] = 3;

	edgeConnect[4][0] = 0;
	edgeConnect[4][1] = 4;

	edgeConnect[5][0] = 1;
	edgeConnect[5][1] = 5;

	edgeConnect[6][0] = 2;
	edgeConnect[6][1] = 6;

	edgeConnect[7][0] = 3;
	edgeConnect[7][1] = 7;

	edgeConnect[8][0] = 4;
	edgeConnect[8][1] = 6;

	edgeConnect[9][0] = 5;
	edgeConnect[9][1] = 7;

	edgeConnect[10][0] = 4;
	edgeConnect[10][1] = 5;

	edgeConnect[11][0] = 6;
	edgeConnect[11][1] = 7;

	// Face data
	nFaces = 6;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = PIXEL;
		facesInfo[k]   = &pixelInfo;
		faceConnect[k] = std::vector<int>(pixelInfo.nVertices);
	}

	faceConnect[0][0] = 2;
	faceConnect[0][1] = 0;
	faceConnect[0][2] = 6;
	faceConnect[0][3] = 4;

	faceConnect[1][0] = 1;
	faceConnect[1][1] = 3;
	faceConnect[1][2] = 5;
	faceConnect[1][3] = 7;

	faceConnect[2][0] = 0;
	faceConnect[2][1] = 1;
	faceConnect[2][2] = 4;
	faceConnect[2][3] = 5;

	faceConnect[3][0] = 3;
	faceConnect[3][1] = 2;
	faceConnect[3][2] = 7;
	faceConnect[3][3] = 6;

	faceConnect[4][0] = 2;
	faceConnect[4][1] = 3;
	faceConnect[4][2] = 0;
	faceConnect[4][3] = 1;

	faceConnect[5][0] = 4;
	faceConnect[5][1] = 5;
	faceConnect[5][2] = 6;
	faceConnect[5][3] = 7;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edgeConnect[0][0] = 1;
	edgeConnect[0][1] = 0;

	edgeConnect[1][0] = 1;
	edgeConnect[1][1] = 2;

	edgeConnect[2][0] = 2;
	edgeConnect[2][1] = 3;

	edgeConnect[3][0] = 3;
	edgeConnect[3][1] = 0;

	edgeConnect[4][0] = 4;
	edgeConnect[4][1] = 5;

	edgeConnect[5][0] = 5;
	edgeConnect[5][1] = 6;

	edgeConnect[6][0] = 6;
	edgeConnect[6][1] = 7;

	edgeConnect[7][0] = 7;
	edgeConnect[7][1] = 4;

	edgeConnect[8][0] = 0;
	edgeConnect[8][1] = 4;

	edgeConnect[9][0] = 1;
	edgeConnect[9][1] = 5;

	edgeConnect[10][0] = 2;
	edgeConnect[10][1] = 6;

	edgeConnect[11][0] = 3;
	edgeConnect[11][1] = 7;

	// Face data
	nFaces = 6;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]    = QUAD;
		facesInfo[k]   = &quadInfo;
		faceConnect[k] = std::vector<int>(quadInfo.nVertices);
	}

	faceConnect[0][0] = 1;
	faceConnect[0][1] = 0;
	faceConnect[0][2] = 3;
	faceConnect[0][3] = 2;

	faceConnect[1][0] = 4;
	faceConnect[1][1] = 5;
	faceConnect[1][2] = 6;
	faceConnect[1][3] = 7;

	faceConnect[2][0] = 7;
	faceConnect[2][1] = 3;
	faceConnect[2][2] = 0;
	faceConnect[2][3] = 4;

	faceConnect[3][0] = 5;
	faceConnect[3][1] = 1;
	faceConnect[3][2] = 2;
	faceConnect[3][3] = 6;

	faceConnect[4][0] = 4;
	faceConnect[4][1] = 0;
	faceConnect[4][2] = 1;
	faceConnect[4][3] = 5;

	faceConnect[5][0] = 6;
	faceConnect[5][1] = 2;
	faceConnect[5][2] = 3;
	faceConnect[5][3] = 7;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edgeConnect[0][0] = 0;
	edgeConnect[0][1] = 1;

	edgeConnect[1][0] = 1;
	edgeConnect[1][1] = 2;

	edgeConnect[2][0] = 2;
	edgeConnect[2][1] = 3;

	edgeConnect[3][0] = 3;
	edgeConnect[3][1] = 0;

	edgeConnect[4][0] = 4;
	edgeConnect[4][1] = 0;

	edgeConnect[5][0] = 4;
	edgeConnect[5][1] = 1;

	edgeConnect[6][0] = 4;
	edgeConnect[6][1] = 2;

	edgeConnect[7][0] = 4;
	edgeConnect[7][1] = 3;

	// Face data
	nFaces = 5;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		if (k == 0) {
			face_type[k]    = QUAD;
			facesInfo[k]   = &quadInfo;
			faceConnect[k] = std::vector<int>(quadInfo.nVertices);
		} else {
			face_type[k]    = TRIANGLE;
			facesInfo[k]   = &triangleInfo;
			faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
		}
	}

	faceConnect[0][0] = 0;
	faceConnect[0][1] = 3;
	faceConnect[0][2] = 2;
	faceConnect[0][3] = 1;

	faceConnect[1][0] = 3;
	faceConnect[1][1] = 0;
	faceConnect[1][2] = 4;

	faceConnect[2][0] = 0;
	faceConnect[2][1] = 1;
	faceConnect[2][2] = 4;

	faceConnect[3][0] = 1;
	faceConnect[3][1] = 2;
	faceConnect[3][2] = 4;

	faceConnect[4][0] = 2;
	faceConnect[4][1] = 3;
	faceConnect[4][2] = 4;

	initializeFaceEdges(facesInfo);
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
	edgeConnect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]    = LINE;
		edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
	}

	edgeConnect[0][0] = 1;
	edgeConnect[0][1] = 0;

	edgeConnect[1][0] = 1;
	edgeConnect[1][1] = 2;

	edgeConnect[2][0] = 2;
	edgeConnect[2][1] = 0;

	edgeConnect[3][0] = 3;
	edgeConnect[3][1] = 4;

	edgeConnect[4][0] = 4;
	edgeConnect[4][1] = 5;

	edgeConnect[5][0] = 5;
	edgeConnect[5][1] = 3;

	edgeConnect[6][0] = 3;
	edgeConnect[6][1] = 0;

	edgeConnect[7][0] = 4;
	edgeConnect[7][1] = 1;

	edgeConnect[8][0] = 5;
	edgeConnect[8][1] = 2;

	// Face data
	nFaces = 5;

	std::vector<ElementInfo *> facesInfo(nFaces);

	face_type = std::vector<Type>(nFaces);
	faceConnect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		if (k == 0 || k == 1) {
			face_type[k]    = TRIANGLE;
			facesInfo[k]   = &triangleInfo;
			faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
		} else {
			face_type[k]    = QUAD;
			facesInfo[k]   = &quadInfo;
			faceConnect[k] = std::vector<int>(quadInfo.nVertices);
		}
	}

	faceConnect[0][0] = 1;
	faceConnect[0][1] = 0;
	faceConnect[0][2] = 2;

	faceConnect[1][0] = 3;
	faceConnect[1][1] = 4;
	faceConnect[1][2] = 5;

	faceConnect[2][0] = 3;
	faceConnect[2][1] = 0;
	faceConnect[2][2] = 1;
	faceConnect[2][3] = 4;

	faceConnect[3][0] = 4;
	faceConnect[3][1] = 1;
	faceConnect[3][2] = 2;
	faceConnect[3][3] = 5;

	faceConnect[4][0] = 5;
	faceConnect[4][1] = 2;
	faceConnect[4][2] = 0;
	faceConnect[4][3] = 3;

	initializeFaceEdges(facesInfo);
}

/*!
	Initializes the list of edges associated to the faces
*/
void ElementInfo::initializeFaceEdges(const std::vector<ElementInfo *> &facesInfo)
{
	faceEdges = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		const ElementInfo &faceInfo = *(facesInfo[k]);

		int nFaceEdges = faceInfo.nFaces;
		for (int i = 0; i < nFaceEdges; ++i) {
			// Number of vertices of the edge associated to the face
			std::size_t nFaceEdgeVertices = faceInfo.faceConnect[i].size();

			// Connectivity of the edge associated to the face
			const std::vector<int> &localFaceEdgeConnect = faceInfo.faceConnect[i];

			std::vector<int> faceEdgeConnect(nFaceEdgeVertices);
			for (std::size_t n = 0; n < nFaceEdgeVertices; ++n) {
				int localVertexId = localFaceEdgeConnect[n];
				int vertexId      = faceConnect[k][localVertexId];

				faceEdgeConnect[n] = vertexId;
			}

			// Search the edge that has the same connectivity of the face edge
			for (int j = 0; j < nEdges; ++j) {
				// If face edge and the guess edge have a different number of
				// vertices, the two edge cannot be the same.
				std::size_t nGuessEdgeVertices = edgeConnect[j].size();
				if (nGuessEdgeVertices != nFaceEdgeVertices) {
					continue;
				}

				// If the connecitivity of the face edge and the one of the
				// guess edge are the same, the two edges coincides.
				const std::vector<int> commonVertices = utils::intersectionVector(faceEdgeConnect, edgeConnect[j]);
				if (commonVertices.size() == nFaceEdgeVertices) {
					faceEdges[k].push_back(j);
				}
			}
		}

		assert(faceEdges[k].size() == nFaceEdges);
	}
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
	_initialize(ElementInfo::UNDEFINED);

	setId(NULL_ID);
}

/*!
	Creates a new element.
*/
Element::Element(const long &id, ElementInfo::Type type)
{
	_initialize(type);

	setId(id);
}

/*!
	Copy constructor
*/
Element::Element(const Element &other)
{
	*this = other;
}

/*!
	Copy-assignament operator.
*/
Element & Element::operator=(const Element& other)
{
	m_id   = other.m_id;
	m_type = other.m_type;

	if (other.m_connect) {
		int nVertices = other.getVertexCount();
		m_connect = std::unique_ptr<long[]>(new long[nVertices]);
		std::copy(other.m_connect.get(), other.m_connect.get() + nVertices, m_connect.get());
	}

	return (*this);
}

/*!
	Initializes the data structures of the element.

	\param type the type of the element
*/
void Element::initialize(ElementInfo::Type type)
{
	_initialize(type);
}

/*!
	Internal function to initialize the data structures of the element.

	\param type the type of the element
*/
void Element::_initialize(ElementInfo::Type type)
{
	setType(type);

	if (getType() != ElementInfo::UNDEFINED) {
		const int &nVertices = getInfo().nVertices;
		setConnect(std::unique_ptr<long[]>(new long[nVertices]));
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
	Gets the basic information of the element.

	\result A constant reference to the basic information of the element.
*/
const ElementInfo & Element::getInfo() const
{
	return ElementInfo::getElementInfo(m_type);
}

/*!
	Sets the element type.

	\param type the element type
*/
void Element::setType(ElementInfo::Type type)
{
	m_type = type;
}

/*!
	Gets the element type.

	\result The element type
*/
ElementInfo::Type Element::getType() const
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
	Gets the number of faces of the element.

	\result The number of vertices of the element
*/
int Element::getFaceCount() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
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
ElementInfo::Type Element::getFaceType(const int &face) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return ElementInfo::UNDEFINED;

	default:
		return getInfo().face_type[face];

	}
}

/*!
	Gets the local connectivity of the specified face of the element.

	\param face is the face for which the connectiviy is reqested
	\result The local connectivity of the specified face of the element.
*/
const std::vector<int> & Element::getFaceLocalConnect(const int &face) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");

	default:
		return getInfo().faceConnect[face];

	}
}

/*!
	Gets the connectivity of the specified face of the element.

	\param face is the face for which the connectiviy is reqested
	\result The connectivity of the specified face of the element.
*/
std::vector<long> Element::getFaceConnect(int face) const
{
	const std::vector<int> &localFaceConnect = getFaceLocalConnect(face);
	int nFaceVertices = localFaceConnect.size();

	std::vector<long> faceConnect;
	faceConnect.resize(nFaceVertices);
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

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return getInfo().nEdges;

	}
}

/*!
	Gets the local connectivity of the specified edge of the element.

	\param edge is the edge for which the connectiviy is reqested
	\result The local connectivity of the specified edge of the element.
*/
const std::vector<int> & Element::getEdgeLocalConnect(const int &edge) const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		throw std::runtime_error ("Unsupported element");

	default:
		return getInfo().edgeConnect[edge];

	}
}

/*!
	Gets the connectivity of the specified edge of the element.

	\param edge is the edge for which the connectiviy is reqested
	\result The connectivity of the specified edge of the element.
*/
std::vector<long> Element::getEdgeConnect(int edge) const
{
	const std::vector<int> &localEdgeConnect = getEdgeLocalConnect(edge);
	int nEdgeVertices = localEdgeConnect.size();

	std::vector<long> edgeConnect;
	edgeConnect.resize(nEdgeVertices);
	for (int k = 0; k < nEdgeVertices; ++k) {
		int localVertexId = localEdgeConnect[k];
		long vertexId = getVertex(localVertexId);
		edgeConnect[k] = vertexId;
	}

	return edgeConnect;
}

/*!
	Gets the dimension of the element.

	\return The dimension of the element
*/
int Element::getDimension() const
{
	switch (m_type) {

	case (ElementInfo::POLYGON):
		return 2;

	case (ElementInfo::POLYHEDRON):
		return 3;

	case (ElementInfo::UNDEFINED):
		BITPIT_UNREACHABLE("Unsupported element");
		return -1;

	default:
		return getInfo().dimension;

	}
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

	case (ElementInfo::POLYGON):
	case (ElementInfo::POLYHEDRON):
	case (ElementInfo::UNDEFINED):
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
        Returns the buffer size required to communicate cell data

        \result buffer size (in bytes)
*/
unsigned int Element::getBinarySize()
{
	return (sizeof(ElementInfo::Type) + (getVertexCount() + 1) * sizeof(long));
}

// Explicit instantiation of the Element containers
template class PiercedVector<Element>;

}
