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

#ifndef __BITPIT_CELL_HPP__
#define __BITPIT_CELL_HPP__

/*! \file */

#include "element.hpp"

#include <memory>

#include <bitpit_containers.hpp>

namespace bitpit {
	class Cell;
}

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, bitpit::Cell& cell);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const bitpit::Cell& cell);

namespace bitpit {

/*!
	\ingroup patch
	@{
*/

class Cell : public Element {

friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream& buf, const Cell& cell);
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream& buf, Cell& cell);

public:
	Cell();
	Cell(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Cell(Cell&& other) = default;
	Cell& operator=(Cell&& other) = default;

	void initialize(ElementInfo::Type type, int nInterfacesPerFace = 0);

	void set_interior(bool interior);
	bool is_interior() const;
	
	void initialize_interfaces(std::vector<std::vector<long>> &interfaces);
	void initialize_empty_interfaces(const std::vector<int> interfaceCount);
	void set_interface(const int &face, const int &index, const long &interface);
	void push_interface(const int &face, const long &interface);
	void delete_interface(const int &face, const int &i);
	void unset_interfaces();
	int get_interface_count() const;
	int get_interface_count(const int &face) const;
	long get_interface(const int &face, const int &index = 0) const;
	const long * get_interfaces() const;
	const long * get_interfaces(const int &face) const;

	void display(std::ostream &out, unsigned short int indent);

	unsigned int get_binary_size( );

protected:

private:
	bool m_interior;

	bitpit::CollapsedVector2D<long> m_interfaces;

	Cell(const Cell &other) = delete;
	Cell& operator = (const Cell &other) = delete;

};

extern template class bitpit::PiercedVector<Cell>;
extern template class bitpit::PositionalPiercedVector<Cell>;

/*!
	@}
*/

}

#endif
