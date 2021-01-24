/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
	class PatchKernel;
}

bitpit::OBinaryStream& operator<< (bitpit::OBinaryStream &out, const bitpit::Vertex &vertex);
bitpit::IBinaryStream& operator>> (bitpit::IBinaryStream &in, bitpit::Vertex &vertex);

namespace bitpit {

class Vertex {

friend class PatchKernel;

friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream&, const Vertex &);
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream&, Vertex &);

public:
	enum Coordinate {
		COORD_X =0,
		COORD_Y,
		COORD_Z,
	};

	Vertex();
	Vertex(long id, bool interior = true);
	Vertex(long id, const std::array<double, 3> &coords, bool interior = true);

	Vertex(const Vertex &other) = default;
	Vertex(Vertex &&other) = default;
	Vertex& operator = (const Vertex &other) = default;
	Vertex& operator=(Vertex &&other) = default;

	void swap(Vertex &other) noexcept;

	void initialize(long id, const std::array<double, 3> &coords, bool interior);

	bool isInterior() const;

	bool operator==(const Vertex &other) const;

	double & operator[](int coord_id);
	double  operator[](int coord_id) const;

	void setId(long id);
	long getId() const;

	void setCoords(const std::array<double, 3> &coords);
	const std::array<double, 3> & getCoords() const;

	void translate(const std::array<double, 3> &translation);
	void translate(double sx, double sy, double sz);
	void scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center);
	void scale(double sx, double sy, double sz, double cx, double cy, double cz);

	static const long NULL_ID;

	unsigned int getBinarySize() const;

	void display(std::ostream &out, unsigned short int indent) const;

protected:
	void setInterior(bool interior);

private:
	std::array<double, 3> m_coords;
	long m_id;
	bool m_interior;

	void _initialize(long id, const std::array<double, 3> &coords, bool interior);

};

extern template class PiercedVector<Vertex>;

}

#endif
