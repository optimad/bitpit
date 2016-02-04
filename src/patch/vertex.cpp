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

#include <limits>

#include "vertex.hpp"

/*!
	Output stream operator for class Vertex. Stream vertex coordinates from
	object to communication buffer.

	\param[in] out_stream is the output memory stream
	\param[in] vertex is the vertex to be streamed
	\result updated output stream
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &out_stream, const bitpit::Vertex &vertex)
{
	out_stream << vertex[0]
		     << vertex[1]
		     << vertex[2];

	return out_stream;
}

/*!
	Input stream operator fro class Vertex. Stream vertex coordinates from
	communication buffer to object

	\param[in] in_stream is the input memory stream
	\param[in] vertex is the vertex to be streamed
	\result updated output stream
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &in_stream, bitpit::Vertex &vertex)
{
	in_stream >> vertex[0];
	in_stream >> vertex[1];
	in_stream >> vertex[2];

	return in_stream;
}


namespace bitpit {

/*!
	\ingroup patch
	@{
*/

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
	Create vertex from input coordinates

	\param[in] id is the id of the vertex
	\param[in] coords are the vertex coordinates
*/
Vertex::Vertex(const long &id, std::array<double, 3> &coords)
{
	set_id(id);
	setCoords(coords);
}

/*!
	Comparison operator

	\param[in] other is the object to be compared with
	\result Returns the boolean result of comparison of the values of the
	arguments, which are not modified.
*/
bool Vertex::operator==(const Vertex &other)
{
	return (m_coords == other.m_coords);
}

/*!
	Get a reference to the specified coordinate

	\param coord_id is the index of the requested coordinate
	\result Returns a reference to requested coordinate
*/
double & Vertex::operator[](int coord_id)
{
	return m_coords[coord_id];
}

/*!
	Get a constants reference to the specified coordinate

	\param coord_id is the index of the requested coordinate
	\result Returns a constant reference to requested coordinate
*/
const double & Vertex::operator[](int coord_id) const
{
	return m_coords[coord_id];
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
void Vertex::setCoords(std::array<double, 3> &coords)
{
	m_coords = coords;
}

/*!
	Gets the coordinates of the vertex.

	\return A pointer to the coordinates of the vertex
*/
const std::array<double, 3> & Vertex::getCoords() const
{
	return m_coords;
}

// Explicit instantiation of the Vertex containers
template class PiercedVector<Vertex>;
template class PositionalPiercedVector<Vertex>;

/*!
	@}
*/

}
