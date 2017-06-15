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

#include <limits>

#include "bitpit_common.hpp"

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
    out_stream << vertex.m_id;
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
    in_stream >> vertex.m_id;
	in_stream >> vertex[0];
	in_stream >> vertex[1];
	in_stream >> vertex[2];

	return in_stream;
}


namespace bitpit {

/*!
	\class Vertex
	\ingroup patchelements

	\brief The Vertex class defines the vertexs.

	Vertex is class that defines the vertexs.
*/

const long Vertex::NULL_ID = std::numeric_limits<long>::min();

/*!
	Default constructor.
*/
Vertex::Vertex()
{
	setId(NULL_ID);
}

/*!
	Creates a new element.
*/
Vertex::Vertex(const long &id)
{
	setId(id);
}

/*!
	Create vertex from input coordinates

	\param[in] id is the id of the vertex
	\param[in] coords are the vertex coordinates
*/
Vertex::Vertex(const long &id, const std::array<double, 3> &coords)
{
	setId(id);
	setCoords(coords);
}

/**
* Exchanges the content of the vertex by the content the specified other
* vertex.
*
* \param other is another vertex whose content is swapped with that of this
* vertex.
*/
void Vertex::swap(Vertex &other) noexcept
{
	std::swap(other.m_id, m_id);
	std::swap(other.m_coords, m_coords);
}

/*!
	Comparison operator

	\param[in] other is the object to be compared with
	\result Returns the boolean result of comparison of the values of the
	arguments, which are not modified.
*/
bool Vertex::operator==(const Vertex &other)
{
	for (int i = 0; i < 3; ++i) {
		if (!utils::DoubleFloatingEqual()(m_coords[i], other.m_coords[i])) {
			return false;
		}
	}

	return true;
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
void Vertex::setId(const long &id)
{
	m_id = id;
}

/*!
	Gets the ID of the vertex.

	\return The ID of the vertex
*/
long Vertex::getId() const
{
	return m_id;
}

/*!
	Sets the coordinates of the vertex.

	\param coords are the coordinates of the vertex
*/
void Vertex::setCoords(const std::array<double, 3> &coords)
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

/*!
	Translates the vertex.

	\param[in] translation is the translation vector
 */
void Vertex::translate(const std::array<double, 3> &translation)
{
	m_coords[0] += translation[0];
	m_coords[1] += translation[1];
	m_coords[2] += translation[2];
}

/*!
	Translates the vertex.

	\param[in] sx translation along x direction
	\param[in] sy translation along y direction
	\param[in] sz translation along z direction
 */
void Vertex::translate(const double &sx, const double &sy, const double &sz)
{
	translate({{sx, sy, sz}});
}

/*!
	Scales the vertex.

	\param[in] scaling is the scaling vector
	\param[in] center is the center of the scaling
 */
void Vertex::scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center)
{
	m_coords[0] = center[0] + scaling[0] * (m_coords[0] - center[0]);
	m_coords[1] = center[1] + scaling[1] * (m_coords[1] - center[1]);
	m_coords[2] = center[2] + scaling[2] * (m_coords[2] - center[2]);
}

/*!
	Scales the vertex.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sz scaling factor along z direction
	\param[in] cx is the x coordinate scaling center
	\param[in] cy is the y coordinate scaling center
	\param[in] cz is the z coordinate scaling center
 */
void Vertex::scale(const double &sx, const double &sy, const double &sz,
                   const double &cx, const double &cy, const double &cz)
{
	scale({{sx, sy, sz}}, {{cx, cy, cz}});
}

/*!
        Display vertex content to output stream in a nicely formatted form

        \param[in,out] out      output stream
        \param[in]     padding  (default = 0) number of trailing spaces
*/
void Vertex::display(std::ostream &out, unsigned int padding)
{
    // Local variables
    std::string         indent(padding, ' ');

    // Display vertex content
    out << indent << "id:     " << m_id << std::endl;
    out << indent << "coords: ";
    for (int i = 0; i < 3; ++i) {
        out << m_coords[i] << " ";
    } //next i
    out << std::endl;

    return;
}

/*!
	Displays vertex information to an output stream

	\param[in] out is the output stream
	\param[in] indent is the number of trailing spaces to prepend when
	writing the information
*/
void Vertex::display(std::ostream &out, unsigned short int indent) const
{
	std::string t_s = std::string(indent, ' ');

	// General info ----------------------------------------------------- //
	out << t_s << "ID:             " << getId() << std::endl;

	// Coordinates ------------------------------------------------------ //
	out << t_s << "coordinates :   (";
	out << m_coords[0] << ", ";
	out << m_coords[1] << ", ";
	out << m_coords[2] << ")" << std::endl;
}

/*!
        Returns the buffer size required to communicate vertex data

        \result buffer size (in bytes)
*/
unsigned int Vertex::getBinarySize()
{
    return (sizeof(m_id) + m_coords.size() * sizeof(double));
}

// Explicit instantiation of the Vertex containers
template class PiercedVector<Vertex>;

}
