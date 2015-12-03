//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "vertex.hpp"

#include <limits>

namespace pman {

/*!
	\class Vertex

	\brief The Vertex class defines the vertexs.

	Vertex is class that defines the vertexs.
*/

const long Vertex::NULL_VERTEX_ID = std::numeric_limits<long>::min();

/*!
	Default constructor.
*/
Vertex::Vertex()
{
	set_id(NULL_VERTEX_ID);
}

/*!
	Creates a new element.
*/
Vertex::Vertex(const long &id)
{
	set_id(id);
}

/*!
	Sets the ID of the vertex.

	\param id the ID of the vertex
*/
void Vertex::set_id(const long &id)
{
	m_id = id;
}

/*!
	Gets the ID of the vertex.

	\return The ID of the vertex
*/
long Vertex::get_id() const
{
	return m_id;
}

/*!
	Sets the coordinates of the vertex.

	\param coords are the coordinates of the vertex
*/
void Vertex::set_coords(std::array<double, 3> &coords)
{
	m_coords = coords;
}

/*!
	Gets the coordinates of the vertex.

	\return A pointer to the coordinates of the vertex
*/
const std::array<double, 3> & Vertex::get_coords() const
{
	return m_coords;
}

}
