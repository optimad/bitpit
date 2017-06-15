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

#ifndef __BITPIT_CELL_HPP__
#define __BITPIT_CELL_HPP__

#include <memory>

#include "bitpit_containers.hpp"

#include "element.hpp"

namespace bitpit {
	class Cell;
	class PatchKernel;
}

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, bitpit::Cell& cell);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const bitpit::Cell& cell);

namespace bitpit {

class Cell : public Element {

friend class PatchKernel;

friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream& buf, const Cell& cell);
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream& buf, Cell& cell);

public:
	Cell();
	Cell(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED,
	     bool interior = true, bool storeNeighbourhood = true);

	Cell(const Cell &other);
	Cell(Cell&& other) = default;
	Cell& operator = (const Cell &other);
	Cell& operator=(Cell&& other) = default;

	void swap(Cell &other) noexcept;

	void initialize(ElementInfo::Type type, bool interior, bool storeNeighbourhood = true);

	bool isInterior() const;

	void setPID(int pid);
	int getPID() const;
	
	void deleteInterfaces();
	void resetInterfaces(bool storeInterfaces = true);
	void setInterfaces(std::vector<std::vector<long>> &interfaces);
	void setInterface(const int &face, const int &index, const long &interface);
	void pushInterface(const int &face, const long &interface);
	void deleteInterface(const int &face, const int &i);
	int getInterfaceCount() const;
	int getInterfaceCount(const int &face) const;
	long getInterface(const int &face, const int &index = 0) const;
	const long * getInterfaces() const;
	long * getInterfaces();
	const long * getInterfaces(const int &face) const;
	long * getInterfaces(const int &face);
	int findInterface(const int &face, const int &interface);
	int findInterface(const int &interface);

	void deleteAdjacencies();
	void resetAdjacencies(bool storeAdjacencies = true);
	void setAdjacencies(std::vector<std::vector<long>> &adjacencies);
	void setAdjacency(const int &face, const int &index, const long &adjacencies);
	void pushAdjacency(const int &face, const long &adjacency);
	void deleteAdjacency(const int &face, const int &i);
	int getAdjacencyCount() const;
	int getAdjacencyCount(const int &face) const;
	long getAdjacency(const int &face, const int &index = 0) const;
	const long * getAdjacencies() const;
	long * getAdjacencies();
	const long * getAdjacencies(const int &face) const;
	long * getAdjacencies(const int &face);
	int findAdjacency(const int &face, const int &adjacency);
	int findAdjacency(const int &adjacency);
        int findVertex(const long &vertex) const;

	bool isFaceBorder(int face) const;

	void display(std::ostream &out, unsigned short int indent) const;

	unsigned int getBinarySize( );

protected:
	void setInterior(bool interior);

private:
	bool m_interior;
	int m_pid;

	bitpit::FlatVector2D<long> m_interfaces;
	bitpit::FlatVector2D<long> m_adjacencies;

	void _initialize(bool interior, bool initializeNeighbourhood, bool storeNeighbourhood);

};

extern template class PiercedVector<Cell>;

}

#endif
