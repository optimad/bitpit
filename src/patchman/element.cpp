//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "element.hpp"
#include "patch.hpp"

#include <limits>

namespace pman {

/*!
	\class ElementInfo

	\brief The ElementInfo class allow to define elements information.

	Element is a struct that hold the basic geometrical information of
	an element.
*/

/*!
	\enum ElementInfo::Type

	This enum defines the element types that can be used.

	\var ElementInfo::Type ElementInfo::UNDEFINED
	An undefined element.

	\var ElementInfo::Type ElementInfo::POINT
	A point.

	\var ElementInfo::Type ElementInfo::LINE
	A line.

	\var ElementInfo::Type ElementInfo::TRIANGLE
	A triangle.

	\var ElementInfo::Type ElementInfo::RECTANGLE
	A rectangle.

	\var ElementInfo::Type ElementInfo::QUADRANGLE
	A quadrangle.

	\var ElementInfo::Type ElementInfo::POLYGON
	A polygon.

	\var ElementInfo::Type ElementInfo::TETRAHEDRON
	A tetrahedron.

	\var ElementInfo::Type ElementInfo::BRICK
	A brick.

	\var ElementInfo::Type ElementInfo::HEXAHEDRON
	A hexahedron.

	\var ElementInfo::Type ElementInfo::PYRAMID
	A pyramid.

	\var ElementInfo::Type ElementInfo::PRISM
	A prism.

	\var ElementInfo::Type ElementInfo::POLYHEDRON
	A polyhedron.
*/

const ElementInfo ElementInfo::undefinedInfo   = ElementInfo(ElementInfo::UNDEFINED);
const ElementInfo ElementInfo::pointInfo       = ElementInfo(ElementInfo::POINT);
const ElementInfo ElementInfo::lineInfo        = ElementInfo(ElementInfo::LINE);
const ElementInfo ElementInfo::triangleInfo    = ElementInfo(ElementInfo::TRIANGLE);
const ElementInfo ElementInfo::rectangleInfo   = ElementInfo(ElementInfo::RECTANGLE);
const ElementInfo ElementInfo::quadrangleInfo  = ElementInfo(ElementInfo::QUADRANGLE);
const ElementInfo ElementInfo::tetrahedronInfo = ElementInfo(ElementInfo::TETRAHEDRON);
const ElementInfo ElementInfo::brickInfo       = ElementInfo(ElementInfo::BRICK);
const ElementInfo ElementInfo::hexahedronInfo  = ElementInfo(ElementInfo::HEXAHEDRON);
const ElementInfo ElementInfo::pyramidInfo     = ElementInfo(ElementInfo::PYRAMID);
const ElementInfo ElementInfo::prismInfo       = ElementInfo(ElementInfo::PRISM);

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

	case (POINT):
		initializePointInfo();
		break;

	case (LINE):
		initializeLineInfo();
		break;

	case (TRIANGLE):
		initializeTriangleInfo();
		break;

	case (RECTANGLE):
		initializeRectangleInfo();
		break;

	case (QUADRANGLE):
		initializeQuadrangleInfo();
		break;

	case (TETRAHEDRON):
		initializeTetrahedronInfo();
		break;

	case (BRICK):
		initializeBrickInfo();
		break;

	case (HEXAHEDRON):
		initializeHexahedronInfo();
		break;

	case (PYRAMID):
		initializePyramidInfo();
		break;

	case (PRISM):
		initializePrismInfo();
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

	case (POINT):
		return pointInfo;

	case (LINE):
		return lineInfo;

	case (TRIANGLE):
		return triangleInfo;

	case (RECTANGLE):
		return rectangleInfo;

	case (QUADRANGLE):
		return quadrangleInfo;

	case (TETRAHEDRON):
		return tetrahedronInfo;

	case (BRICK):
		return brickInfo;

	case (HEXAHEDRON):
		return hexahedronInfo;

	case (PYRAMID):
		return pyramidInfo;

	case (PRISM):
		return prismInfo;

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
void ElementInfo::initializePointInfo()
{
	type      = POINT;
	dimension = 0;

	// Vertices data
	nVertices = 1;

	// Edge data
	nEdges = 1;

	edge_type = std::vector<Type>(nEdges);
	edge_type[0] = POINT;

	edge_connect = std::vector<std::vector<int>>(nEdges);
	edge_connect[0] = std::vector<int>(nVertices);
	edge_connect[0][0] = 0;

	// Face data
	nFaces = 1;

	face_type = std::vector<Type>(nFaces);
	face_type[0] = POINT;

	face_connect = std::vector<std::vector<int>>(nFaces);
	face_connect[0] = std::vector<int>(nVertices);
	face_connect[0][0] = 0;
}

/*!
	Initializes the information for the line element.
*/
void ElementInfo::initializeLineInfo()
{
	ElementInfo pointInfo(POINT);

	type      = LINE;
	dimension = 1;

	// Vertices data
	nVertices = 2;

	// Edge data
	nEdges = 2;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = POINT;
		edge_connect[k]    = std::vector<int>(pointInfo.nVertices);
		edge_connect[k][0] = k;
	}

	// Face data
	nFaces = 2;

	face_type = std::vector<Type>(nFaces);
	face_connect = std::vector<std::vector<int>>(nFaces);
	for (int k = 0; k < nFaces; ++k) {
		face_type[k]       = POINT;
		face_connect[k]    = std::vector<int>(pointInfo.nVertices);
		face_connect[k][0] = k;
	}
}

/*!
	Initializes the information for the triangle element.
*/
void ElementInfo::initializeTriangleInfo()
{
	ElementInfo pointInfo(POINT);
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
		edge_type[k]       = POINT;
		edge_connect[k]    = std::vector<int>(pointInfo.nVertices);
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
		face_connect[k][1] = (k + 1) % lineInfo.nVertices;
	}
}

/*!
	Initializes the information for the quadrangle element.
*/
void ElementInfo::initializeRectangleInfo()
{
	ElementInfo pointInfo(POINT);
	ElementInfo lineInfo(LINE);

	type      = QUADRANGLE;
	dimension = 2;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 4;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = POINT;
		edge_connect[k]    = std::vector<int>(pointInfo.nVertices);
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
void ElementInfo::initializeQuadrangleInfo()
{
	ElementInfo pointInfo(POINT);
	ElementInfo lineInfo(LINE);

	type      = QUADRANGLE;
	dimension = 2;

	// Vertices data
	nVertices = 4;

	// Edge data
	nEdges = 4;

	edge_type = std::vector<Type>(nEdges);
	edge_connect = std::vector<std::vector<int>>(nEdges);
	for (int k = 0; k < nEdges; ++k) {
		edge_type[k]       = POINT;
		edge_connect[k]    = std::vector<int>(pointInfo.nVertices);
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
		face_connect[k][1] = (k + 1) % lineInfo.nVertices;
	}
}

/*!
	Initializes the information for the tetrahedron element.
*/
void ElementInfo::initializeTetrahedronInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo triangleInfo(TRIANGLE);

	type      = TETRAHEDRON;
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
void ElementInfo::initializeBrickInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo rectangleInfo(RECTANGLE);

	type      = BRICK;
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
		face_type[k]    = RECTANGLE;
		face_connect[k] = std::vector<int>(rectangleInfo.nVertices);
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
	ElementInfo quadrangleInfo(QUADRANGLE);

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
		face_type[k]    = QUADRANGLE;
		face_connect[k] = std::vector<int>(quadrangleInfo.nVertices);
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
	ElementInfo quadrangleInfo(QUADRANGLE);

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
			face_type[k]    = QUADRANGLE;
			face_connect[k] = std::vector<int>(quadrangleInfo.nVertices);
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
	Initializes the information for the prism element.
*/
void ElementInfo::initializePrismInfo()
{
	ElementInfo lineInfo(LINE);
	ElementInfo triangleInfo(TRIANGLE);
	ElementInfo quadrangleInfo(QUADRANGLE);

	type      = PRISM;
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
			face_type[k]    = QUADRANGLE;
			face_connect[k] = std::vector<int>(quadrangleInfo.nVertices);
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

	set_patch(NULL);
	set_id(NULL_ELEMENT_ID);
}

/*!
	Creates a new element.
*/
Element::Element(const long &id, Patch *patch)
{
	initialize(ElementInfo::UNDEFINED);

	set_patch(patch);
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
	Sets the patch that owns the element.

	\param patch the patch that owns the element
*/
void Element::set_patch(Patch *patch)
{
	m_patch = patch;
}

/*!
	Gets the patch that owns the element.

	\return The patch that owns the element.
*/
Patch * Element::get_patch() const
{
	return m_patch;
}

/*!
	Gets the dimension of the owner patch.

	\return The dimension of the owner patch.
*/
int Element::get_patch_dimension() const
{
	return m_patch->get_dimension();
}

/*!
	Returns true if the owner patch is a three-dimensional patch.

	\return This method returns true if the owner patch is
	        three-dimensional, false otherwise.
*/
bool Element::is_patch_three_dimensional() const
{
	return m_patch->is_three_dimensional();
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

	\result A pointer to the connectivity of the element
*/
const long * Element::get_connect() const
{
	return m_connect.get();
}

/*!
	Sets the centroid of the element.

	\param centroid the centroid of the element.
*/
void Element::set_centroid(std::array<double, 3> &centroid)
{
	m_centroid = centroid;
}

/*!
	Gets the centroid of the element.

	\return The centroid of the element.
*/
const std::array<double, 3> & Element::get_centroid() const
{
	return m_centroid;
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
int Element::get_vertex(const int &vertex) const
{
	return m_connect[vertex];
}

/*!
	Evaluates the characteristic length of the element.

	\result The characteristic of the element.
*/
double Element::eval_length() const
{
	double length;
	if (m_type == ElementInfo::POINT) {
		length = 0.;
	} else if (m_type == ElementInfo::LINE) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);

		const std::array<double, 3> &coords_A = node_A.get_coords();
		const std::array<double, 3> &coords_B = node_B.get_coords();

		length = 0.;
		for (int k = 0; k < get_patch_dimension(); ++k) {
			length += pow(coords_B[k] - coords_A[k], 2);
		}
		length = pow(length, 0.5);
	} else if (m_type == ElementInfo::RECTANGLE) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);
		Node &node_C = m_patch->get_vertex(m_connect[3]);

		const std::array<double, 3> &coords_A = node_A.get_coords();
		const std::array<double, 3> &coords_B = node_B.get_coords();
		const std::array<double, 3> &coords_C = node_C.get_coords();

		double length_x = 0.0;
		double length_y = 0.0;
		for (unsigned int k = 0; k < coords_A.size(); ++k) {
			length_x += pow(coords_B[k] - coords_A[k], 2);
			length_y += pow(coords_C[k] - coords_A[k], 2);
		}
		length = pow(std::min({length_x, length_y}), 0.5);
	} else if (m_type == ElementInfo::BRICK) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);
		Node &node_C = m_patch->get_vertex(m_connect[3]);
		Node &node_D = m_patch->get_vertex(m_connect[4]);

		const std::array<double, 3> &coords_A = node_A.get_coords();
		const std::array<double, 3> &coords_B = node_B.get_coords();
		const std::array<double, 3> &coords_C = node_C.get_coords();
		const std::array<double, 3> &coords_D = node_D.get_coords();

		double length_x = 0.0;
		double length_y = 0.0;
		double length_z = 0.0;
		for (unsigned int k = 0; k < coords_A.size(); ++k) {
			length_x += pow(coords_B[k] - coords_A[k], 2);
			length_y += pow(coords_C[k] - coords_A[k], 2);
			length_z += pow(coords_D[k] - coords_A[k], 2);
		}
		length = pow(std::min({length_x, length_y, length_z}), 0.5);
	} else {
		length = 0.0;
	}

	return length;
}

}
