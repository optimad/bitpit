//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "element.hpp"
#include "patch.hpp"

#include <limits>

namespace pman {

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

	\var Element::Type Element::RECTANGLE
	A rectangle.

	\var Element::Type Element::QUADRANGLE
	A quadrangle.

	\var Element::Type Element::POLYGON
	A polygon.

	\var Element::Type Element::TETRAHEDRON
	A tetrahedron.

	\var Element::Type Element::BRICK
	A brick.

	\var Element::Type Element::HEXAHEDRON
	A hexahedron.

	\var Element::Type Element::PYRAMID
	A pyramid.

	\var Element::Type Element::PRISM
	A prism.

	\var Element::Type Element::POLYHEDRON
	A polyhedron.
*/

const long Element::NULL_ELEMENT_ID = std::numeric_limits<long>::min();

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
Element::Element(const long &id)
{
	set_patch(NULL);
	set_id(id);
}

/*!
	Creates a new element.
*/
Element::Element(const long &id, Patch *patch)
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
	Sets the local ID of the element.

	\param id the local ID of the element

*/
void Element::set_local_id(long id)
{
	m_local_id = id;
}

/*!
	Gets the local ID of the element.

	\return The local ID of the element
*/
long Element::get_local_id() const
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

	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
	    return 4;

	case (Type::POLYGON):
	    return -1;

	case (Type::TETRAHEDRON):
	    return 4;

	case (Type::BRICK):
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
	Gets the face type of the specified face of the element.

	\result The face type of specified face of the element
*/
Element::Type Element::get_face_type(const int &face) const
{
	return get_face_type(m_type, face);
}

/*!
	Gets the face type of the specified face of the given type of element.

	\param type the type of the element
	\result The face type of specified face of the element
*/
Element::Type Element::get_face_type(Element::Type type, const int &face)
{
	switch (type) {

	case (Type::POINT):
		return Type::POINT;

	case (Type::LINE):
		return Type::POINT;

	case (Type::TRIANGLE):
	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
		return Type::LINE;

	case (Type::POLYGON):
		return Type::UNDEFINED;

	case (Type::TETRAHEDRON):
		return Type::TRIANGLE;

	case (Type::BRICK):
		return Type::RECTANGLE;

	case (Type::HEXAHEDRON):
		return Type::QUADRANGLE;

	case (Type::PYRAMID):
		if (face == 0) {
			return Type::QUADRANGLE;
		} else {
			return Type::TRIANGLE;
		}

	case (Type::PRISM):
		if (face == 0 || face == 1) {
			return Type::TRIANGLE;
		} else {
			return Type::QUADRANGLE;
		}

	case (Type::POLYHEDRON):
		return Type::UNDEFINED;

	default:
		return Type::UNDEFINED;

	}

	return Type::UNDEFINED;
}

/*!
	Gets the local connectivity of the specified face of the element.

	\param face is the face for which the connectiviy is reqested
	\result The local connectivity of the specified face of the element.
*/
std::vector<int> Element::get_face_local_connect(const int &face) const
{
	return get_face_local_connect(m_type, face);
}

/*!
	Gets the local connectivity of the specified face of the given type
	of element.

	\param type the type of the element
	\param face is the face for which the connectiviy is reqested
	\result The local connectivity of the specified face of the element.
*/
std::vector<int> Element::get_face_local_connect(Element::Type type, const int &face)
{
	Element::Type faceType = get_face_type(type, face);
	std::vector<int> faceConnect(get_vertex_count(faceType));

	switch (type) {

	case (Type::POINT):
	    faceConnect[0] = 0;

	    break;

	case (Type::LINE):
		if (face == 0) {
			faceConnect[0] = 0;
		} else {
			faceConnect[1] = 1;
		}

	    break;

	case (Type::TRIANGLE):
		if (face == 0) {
			faceConnect[0] = 0;
			faceConnect[1] = 1;
		} else if (face == 1) {
			faceConnect[0] = 1;
			faceConnect[1] = 2;
		} else if (face == 2) {
			faceConnect[0] = 2;
			faceConnect[1] = 0;
		}

	    break;

	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
		if (face == 0) {
			faceConnect[0] = 2;
			faceConnect[1] = 0;
		} else if (face == 1) {
			faceConnect[0] = 1;
			faceConnect[1] = 3;
		} else if (face == 2) {
			faceConnect[0] = 0;
			faceConnect[1] = 1;
		} else if (face == 3) {
			faceConnect[0] = 3;
			faceConnect[1] = 2;
		}

	    break;

	case (Type::POLYGON):
	    break;

	case (Type::TETRAHEDRON):
		if (face == 0) {
			faceConnect[0] = 1;
			faceConnect[1] = 0;
			faceConnect[2] = 2;
		} else if (face == 1) {
			faceConnect[0] = 0;
			faceConnect[1] = 3;
			faceConnect[2] = 2;
		} else if (face == 2) {
			faceConnect[0] = 3;
			faceConnect[1] = 1;
			faceConnect[2] = 2;
		} else if (face == 3) {
			faceConnect[0] = 0;
			faceConnect[1] = 1;
			faceConnect[2] = 3;
		}

	    break;

	case (Type::BRICK):
	case (Type::HEXAHEDRON):
		if (face == 0) {
			faceConnect[0] = 2;
			faceConnect[1] = 0;
			faceConnect[2] = 4;
			faceConnect[3] = 6;
		} else if (face == 1) {
			faceConnect[0] = 1;
			faceConnect[1] = 3;
			faceConnect[2] = 7;
			faceConnect[3] = 5;
		} else if (face == 2) {
			faceConnect[0] = 0;
			faceConnect[1] = 1;
			faceConnect[2] = 5;
			faceConnect[3] = 4;
		} else if (face == 3) {
			faceConnect[0] = 3;
			faceConnect[1] = 2;
			faceConnect[2] = 6;
			faceConnect[3] = 7;
		} else if (face == 4) {
			faceConnect[0] = 2;
			faceConnect[1] = 3;
			faceConnect[2] = 1;
			faceConnect[3] = 0;
		} else if (face == 5) {
			faceConnect[0] = 4;
			faceConnect[1] = 5;
			faceConnect[2] = 7;
			faceConnect[3] = 6;
		}

	    break;

	case (Type::PYRAMID):
		if (face == 0) {
			faceConnect[0] = 0;
			faceConnect[1] = 3;
			faceConnect[2] = 2;
			faceConnect[3] = 1;
		} else if (face == 1) {
			faceConnect[0] = 3;
			faceConnect[1] = 0;
			faceConnect[2] = 4;
		} else if (face == 2) {
			faceConnect[0] = 0;
			faceConnect[1] = 1;
			faceConnect[2] = 4;
		} else if (face == 3) {
			faceConnect[0] = 1;
			faceConnect[1] = 2;
			faceConnect[2] = 4;
		} else if (face == 4) {
			faceConnect[0] = 2;
			faceConnect[1] = 3;
			faceConnect[2] = 4;
		}

	    break;

	case (Type::PRISM):
		if (face == 0) {
			faceConnect[0] = 1;
			faceConnect[1] = 0;
			faceConnect[2] = 2;
		} else if (face == 1) {
			faceConnect[0] = 3;
			faceConnect[1] = 4;
			faceConnect[2] = 5;
		} else if (face == 2) {
			faceConnect[0] = 3;
			faceConnect[1] = 0;
			faceConnect[2] = 1;
			faceConnect[3] = 4;
		} else if (face == 3) {
			faceConnect[0] = 4;
			faceConnect[1] = 1;
			faceConnect[2] = 2;
			faceConnect[3] = 5;
		} else if (face == 4) {
			faceConnect[0] = 5;
			faceConnect[1] = 2;
			faceConnect[2] = 0;
			faceConnect[3] = 3;
		}

	    break;

	case (Type::POLYHEDRON):

	    break;

	default:

	    break;

	}

	return faceConnect;
}

/*!
	Gets the number of edges of the element.

	\result The number of edges of the element
*/
int Element::get_edge_count() const
{
	return get_edge_count(m_type);
}

/*!
	Gets the number of edges of the given type of element.

	\param type the type of the element
	\result The number of edges of the element
*/
int Element::get_edge_count(Element::Type type)
{
	switch (type) {

	case (Type::POINT):
		return 1;

	case (Type::LINE):
		return 2;

	case (Type::TRIANGLE):
		return 3;

	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
		return 4;

	case (Type::POLYGON):
		return -1;

	case (Type::TETRAHEDRON):
		return 6;

	case (Type::BRICK):
	case (Type::HEXAHEDRON):
		return 12;

	case (Type::PYRAMID):
		return 8;

	case (Type::PRISM):
		return 9;

	case (Type::POLYHEDRON):
		return -1;

	default:
		return -1;

	}

	return -1;
}

/*!
	Gets the local connectivity of the specified edge of the element.

	\param edge is the edge for which the connectiviy is reqested
	\result The local connectivity of the specified edge of the element.
*/
std::vector<int> Element::get_edge_local_connect(const int &edge) const
{
	return get_edge_local_connect(m_type, edge);
}

/*!
	Gets the local connectivity of the specified edge of the given type
	of element.

	\param type the type of the element
	\param edge is the edge for which the connectiviy is reqested
	\result The local connectivity of the specified edge of the element.
*/
std::vector<int> Element::get_edge_local_connect(Element::Type type, const int &edge)
{
	int nEdgeVertices = std::max(get_dimension(type) - 1, 1);
	std::vector<int> edgeConnect(nEdgeVertices);

	switch (type) {

	case (Type::POINT):
	    edgeConnect[0] = 0;

	    break;

	case (Type::LINE):
	case (Type::TRIANGLE):
	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
	case (Type::POLYGON):
		edgeConnect[0] = edge;

		break;

	case (Type::TETRAHEDRON):
		if (edge == 0) {
			edgeConnect[0] = 0;
			edgeConnect[1] = 1;
		} else if (edge == 1) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 2;
		} else if (edge == 2) {
			edgeConnect[0] = 2;
			edgeConnect[1] = 0;
		} else if (edge == 3) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 0;
		} else if (edge == 2) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 1;
		} else if (edge == 3) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 2;
		}

	    break;

	case (Type::BRICK):
	case (Type::HEXAHEDRON):
		if (edge == 0) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 0;
		} else if (edge == 1) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 2;
		} else if (edge == 2) {
			edgeConnect[0] = 2;
			edgeConnect[1] = 3;
		} else if (edge == 3) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 0;
		} else if (edge == 4) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 5;
		} else if (edge == 5) {
			edgeConnect[0] = 5;
			edgeConnect[1] = 6;
		} else if (edge == 6) {
			edgeConnect[0] = 6;
			edgeConnect[1] = 7;
		} else if (edge == 7) {
			edgeConnect[0] = 7;
			edgeConnect[1] = 4;
		} else if (edge == 8) {
			edgeConnect[0] = 0;
			edgeConnect[1] = 4;
		} else if (edge == 9) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 5;
		} else if (edge == 10) {
			edgeConnect[0] = 2;
			edgeConnect[1] = 6;
		} else if (edge == 11) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 7;
		}

	    break;

	case (Type::PYRAMID):
		if (edge == 0) {
			edgeConnect[0] = 0;
			edgeConnect[1] = 1;
		} else if (edge == 1) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 2;
		} else if (edge == 2) {
			edgeConnect[0] = 2;
			edgeConnect[1] = 3;
		} else if (edge == 3) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 0;
		} else if (edge == 4) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 0;
		} else if (edge == 5) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 1;
		} else if (edge == 6) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 2;
		} else if (edge == 7) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 3;
		}

	    break;

	case (Type::PRISM):
		if (edge == 0) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 0;
		} else if (edge == 1) {
			edgeConnect[0] = 1;
			edgeConnect[1] = 2;
		} else if (edge == 2) {
			edgeConnect[0] = 2;
			edgeConnect[1] = 0;
		} else if (edge == 3) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 4;
		} else if (edge == 4) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 5;
		} else if (edge == 5) {
			edgeConnect[0] = 5;
			edgeConnect[1] = 3;
		} else if (edge == 6) {
			edgeConnect[0] = 3;
			edgeConnect[1] = 0;
		} else if (edge == 7) {
			edgeConnect[0] = 4;
			edgeConnect[1] = 1;
		} else if (edge == 8) {
			edgeConnect[0] = 5;
			edgeConnect[1] = 2;
		}

	    break;

	case (Type::POLYHEDRON):

	    break;

	default:

	    break;

	}

	return edgeConnect;
}

/*!
	Gets the dimension of the element.

	\return The dimension of the element
*/
int Element::get_dimension() const
{
	return get_dimension(m_type);
}

/*!
	Gets the dimension of the given type of element.

	\param type the type of the element
	\return The dimension of the element
*/
int Element::get_dimension(Element::Type type)
{
	switch (type) {

	case (Type::POINT):
	    return 0;

	case (Type::LINE):
	    return 1;

	case (Type::TRIANGLE):
	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
	case (Type::POLYGON):
	    return 2;

	case (Type::TETRAHEDRON):
	case (Type::BRICK):
	case (Type::HEXAHEDRON):
	case (Type::PYRAMID):
	case (Type::PRISM):
	case (Type::POLYHEDRON):
	    return 3;

	default:
	    return -1;

	}

	return -1;
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
	Returns true if the given type of element is a three-dimensional element.

	\param type the type of the element
	\return Returns true if the element is a three-dimensional element,
	false otherwise.
*/
bool Element::is_three_dimensional(Element::Type type)
{
	return (get_dimension(type) == 3);
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

	case (Type::RECTANGLE):
	case (Type::QUADRANGLE):
	    return 4;

	case (Type::POLYGON):
	    return -1;

	case (Type::TETRAHEDRON):
	    return 4;

	case (Type::BRICK):
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
	Evaluates the minimum length of the element.

	\result The minimum length of the element.
*/
double Element::eval_min_length() const
{
	double length;
	if (m_type == Type::POINT) {
		length = 0.;
	} else if (m_type == Type::LINE) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);

		const std::array<double, 3> &coords_A = node_A.get_coords();
		const std::array<double, 3> &coords_B = node_B.get_coords();

		length = 0.;
		for (int k = 0; k < get_patch_dimension(); ++k) {
			length += pow(coords_B[k] - coords_A[k], 2);
		}
		length = pow(length, 0.5);
	} else if (m_type == Type::RECTANGLE) {
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
	} else if (m_type == Type::BRICK) {
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

/*!
	Evaluates the cross product of two arrays.

	\param x the first a three-dimensional array
	\param y the second a three-dimensional array
	\param cross the array used to store the cross product.
*/
void Element::cross_3D(std::array<double, 3> x, std::array<double, 3> y, std::array<double, 3> cross)
{
	cross[0] = x[1] * y[2] - x[2] * y[1];
	cross[1] = x[2] * y[0] - x[0] * y[2];
	cross[2] = x[0] * y[1] - x[1] * y[0];
}

/*!
	Normalizes an array.

	\param x the three-dimensional array to be normalized
*/
void Element::normalize_3D(std::array<double, 3> &x)
{
	double module = 0.;
	for (int k = 0; k < 3; ++k) {
		module += x[k] * x[k];
	}
	module = pow(module, 0.5);

	for (int k = 0; k < 3; ++k) {
		x[k] /= module;
	}
}

/*!
	Evaluates the transpose of a matrix.

	\param A the 3x3 matrix
*/
void Element::transpose_3D(std::array<std::array<double, 3>, 3> &A)
{
	for (int i = 1; i < 3; i++) {
		for (int j = i; j < 3; j++) {
			double tmp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = tmp;
		}
	}
}


/*!
	Adds an id to an ordered list of unique ids.

	\param id is the id to be added
	\param list is the ordered list of uniqe ids
	\result Returns true is the id was added to the list, false otherwise.
*/
bool Element::add_id_to_ordered_list(const long &id, std::vector<long> &list)
{
	if (list.empty()) {
		list.push_back(id);
		return true;
	}

	std::vector<long>::iterator itr = lower_bound(list.begin(), list.end(), id);
	if (itr == list.end() || *itr != id) {
		list.insert(itr, id);
		return true;
	} else {
		return false;
	}
}

}
