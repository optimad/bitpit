//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\class Element

	\brief The Element class provides an interface for defining elements.

	Element is the base class for defining elements like cells and
	intefaces.
*/

/*!
	\enum Element::Type

	This enum defines the element types that can be used.

	\var Element::Type Element::POINT
	A point.

	\var Element::Type Element::LINE
	A line.

	\var Element::Type Element::TRIANGLE
	A triangle.

	\var Element::Type Element::QUADRANGLE
	A quadrangle.

	\var Element::Type Element::POLYGON
	A polygon.

	\var Element::Type Element::TETRAHEDRON
	A tetrahedron.

	\var Element::Type Element::HEXAHEDRON
	A hexahedron.

	\var Element::Type Element::PYRAMID
	A pyramid.

	\var Element::Type Element::PRISM
	A prism.

	\var Element::Type Element::POLYHEDRON
	A polyhedron.
*/

#include "element.hpp"

#include <limits>

namespace pman {

const int Element::NULL_ELEMENT_ID = std::numeric_limits<int>::min();

/*!
	Default constructor.
*/
Element::Element()
{
	set_patch(NULL);
	set_id(NULL_ELEMENT_ID);
}


/*!
	Creates a new element.
*/
Element::Element(const int &id)
{
	set_patch(NULL);
	set_id(id);
}

/*!
	Creates a new element.
*/
Element::Element(const int &id, Patch *patch)
{
	set_patch(patch);
	set_id(id);
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
	Sets the ID of the element.

	\param id the ID of the element
*/
void Element::set_id(const int &id)
{
	m_id = id;
}

/*!
	Gets the ID of the element.

	\return The ID of the element
*/
int Element::get_id() const
{
	return m_id;
}

/*!
	Sets the local ID of the element.

	\param id the local ID of the element

*/
void Element::set_local_id(int id)
{
	m_local_id = id;
}

/*!
	Gets the local ID of the element.

	\return The local ID of the element
*/
int Element::get_local_id() const
{
	return m_local_id;
}

/*!
	Sets the element type.

	\param type the element type
*/
void Element::set_type(Element::Type type)
{
	m_type = type;
}

/*!
	Gets the element type.

	\result The element type
*/
Element::Type Element::get_type() const
{
	return m_type;
}

/*!
	Sets the vertex connectivity of the element.

	\param connect a pointer to the connectivity of the element
*/
void Element::set_connect(std::unique_ptr<int[]> connect)
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
int* Element::get_connect() const
{
	return m_connect.get();
}

/*!
	Gets the number of faces of the element.

	\result The number of vertices of the element
*/
int Element::get_face_count() const
{
	return get_face_count(m_type);
}

/*!
	Gets the number of faces of the given type of element.

	\param type the type of the element
	\result The number of faces of the element
*/
int Element::get_face_count(Element::Type type)
{
	switch (type) {

	case (Type::POINT):
	    return 1;

	case (Type::LINE):
	    return 2;

	case (Type::TRIANGLE):
	    return 3;

	case (Type::QUADRANGLE):
	    return 4;

	case (Type::POLYGON):
	    return -1;

	case (Type::TETRAHEDRON):
	    return 4;

	case (Type::HEXAHEDRON):
	    return 6;

	case (Type::PYRAMID):
	    return 5;

	case (Type::PRISM):
	    return 5;

	case (Type::POLYHEDRON):
	    return -1;

	default:
	    return -1;

	}

	return -1;
}

/*!
	Gets the number of vertices of the element.

	\result The number of vertices of the element
*/
int Element::get_vertex_count() const
{
	return get_vertex_count(m_type);
}

/*!
	Gets the number of vertices of the given type of element.

	\param type the type of the element
	\result The number of vertices of the element
*/
int Element::get_vertex_count(Element::Type type)
{
	switch (type) {

	case (Type::POINT):
	    return 1;

	case (Type::LINE):
	    return 2;

	case (Type::TRIANGLE):
	    return 3;

	case (Type::QUADRANGLE):
	    return 4;

	case (Type::POLYGON):
	    return -1;

	case (Type::TETRAHEDRON):
	    return 4;

	case (Type::HEXAHEDRON):
	    return 8;

	case (Type::PYRAMID):
	    return 5;

	case (Type::PRISM):
	    return 6;

	case (Type::POLYHEDRON):
	    return -1;

	default:
	    return -1;

	}

	return -1;
}

}
