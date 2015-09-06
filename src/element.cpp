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

#include "element.hpp"
#include "patch.hpp"

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
const int * Element::get_connect() const
{
	return m_connect.get();
}

/*!
	Sets the centroid of the element.

	\param centroid the centroid of the element.
*/
void Element::set_centroid(std::unique_ptr<double[]> centroid)
{
	m_centroid = std::move(centroid);
}

/*!
	Gets the centroid of the element.

	\return The centroid of the element.
*/
const double * Element::get_centroid() const
{
	return m_centroid.get();
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

		const double *coords_A = node_A.get_coords();
		const double *coords_B = node_B.get_coords();

		length = 0.;
		for (int k = 0; k < get_patch_dimension(); ++k) {
			length += pow(coords_B[k] - coords_A[k], 2);
		}
		length = pow(length, 0.5);
	} else if (m_type == Type::RECTANGLE) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);
		Node &node_C = m_patch->get_vertex(m_connect[3]);

		const double *coords_A = node_A.get_coords();
		const double *coords_B = node_B.get_coords();
		const double *coords_C = node_C.get_coords();

		double length_x = 0.0;
		double length_y = 0.0;
		for (int k = 0; k < get_patch_dimension(); ++k) {
			length_x += pow(coords_B[k] - coords_A[k], 2);
			length_y += pow(coords_C[k] - coords_A[k], 2);
		}
		length = pow(std::min({length_x, length_y}), 0.5);
	} else if (m_type == Type::BRICK) {
		Node &node_A = m_patch->get_vertex(m_connect[0]);
		Node &node_B = m_patch->get_vertex(m_connect[1]);
		Node &node_C = m_patch->get_vertex(m_connect[3]);
		Node &node_D = m_patch->get_vertex(m_connect[4]);

		const double *coords_A = node_A.get_coords();
		const double *coords_B = node_B.get_coords();
		const double *coords_C = node_C.get_coords();
		const double *coords_D = node_C.get_coords();

		double length_x = 0.0;
		double length_y = 0.0;
		double length_z = 0.0;
		for (int k = 0; k < get_patch_dimension(); ++k) {
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
void Element::cross(double x[], double y[], double cross[])
{
	cross[0] = x[1] * y[2] - x[2] * y[1];
	cross[1] = x[2] * y[0] - x[0] * y[2];
	cross[2] = x[0] * y[1] - x[1] * y[0];
}

/*!
	Normalizes an array.

	\param x the array to be normalized
	\param size the size of the array, if no size is specified the default
	            size of 3 will be used
*/
void Element::normalize(double x[], int size)
{
	double module = 0.;
	for (int k = 0; k < size; ++k) {
		module += x[k] * x[k];
	}
	module = pow(module, 0.5);

	for (int k = 0; k < size; ++k) {
		x[k] /= module;
	}
}

/*!
	Evaluates the transpose of a matrix.

	\param A the matrix
	\param nRows the number of rows of the matrix
	\param nCols the number of columns of the matrix
*/
void Element::transpose(double **A, const int &nRows, const int &nCols)
{
	for (int i = 1; i < nRows; i++) {
		for (int j = i; j < nCols; j++) {
			double tmp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = tmp;
		}
	}
}

}
