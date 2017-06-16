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

#ifndef __BITPIT_VERTEX_HPP__
#define __BITPIT_VERTEX_HPP__

#include <array>
#include <memory>

#include "bitpit_containers.hpp"

namespace bitpit {
	class Vertex;
}

bitpit::OBinaryStream& operator<< (bitpit::OBinaryStream &out, const bitpit::Vertex &vertex);
bitpit::IBinaryStream& operator>> (bitpit::IBinaryStream &in, bitpit::Vertex &vertex);

namespace bitpit {

class Vertex {

friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream&, const Vertex &);
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream&, Vertex &);

public:
	enum Coordinate {
		COORD_X =0,
		COORD_Y,
		COORD_Z,
	};

	Vertex();
	Vertex(const long &id);
	Vertex(const long &id, const std::array<double, 3> &coords);

	Vertex(const Vertex &other) = default;
	Vertex(Vertex &&other) = default;
	Vertex& operator = (const Vertex &other) = default;
	Vertex& operator=(Vertex &&other) = default;

	void swap(Vertex &other) noexcept;

	void initialize(long id, const std::array<double, 3> &coords);

	bool operator==(const Vertex &other);

	double & operator[](int coord_id);
	const double & operator[](int coord_id) const;

	void setId(const long &id);
	long getId() const;

	void setCoords(const std::array<double, 3> &coords);
	const std::array<double, 3> & getCoords() const;

	void translate(const std::array<double, 3> &translation);
	void translate(const double &sx, const double &sy, const double &sz);
	void scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center);
	void scale(const double &sx, const double &sy, const double &sz,
	           const double &cx, const double &cy, const double &cz);
        void display(std::ostream &, unsigned int padding = 0);

	static const long NULL_ID;

    unsigned int getBinarySize();

	void display(std::ostream &out, unsigned short int indent) const;

private:
	long m_id;
	std::array<double, 3> m_coords;

	void _initialize(long id, const std::array<double, 3> &coords);

};

extern template class PiercedVector<Vertex>;

}

#endif
