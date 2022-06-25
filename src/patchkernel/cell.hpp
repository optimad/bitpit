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

#ifndef __BITPIT_CELL_HPP__
#define __BITPIT_CELL_HPP__

#include <memory>

#include "bitpit_containers.hpp"

#include "element.hpp"

namespace bitpit {

class Cell;
class PatchKernel;

IBinaryStream & operator>>(IBinaryStream &buf, Cell& cell);
OBinaryStream & operator<<(OBinaryStream &buf, const Cell& cell);

class Cell : public Element {

friend class PatchKernel;

friend OBinaryStream& (operator<<) (OBinaryStream& buf, const Cell& cell);
friend IBinaryStream& (operator>>) (IBinaryStream& buf, Cell& cell);

public:
	Cell();
	BITPIT_DEPRECATED(Cell(long id, ElementType type, bool interior, bool storeNeighbourhood));
	BITPIT_DEPRECATED(Cell(long id, ElementType type, int connectSize, bool interior, bool storeNeighbourhood));
	BITPIT_DEPRECATED(Cell(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeNeighbourhood));
	Cell(long id, ElementType type, bool interior = true, bool storeInterfaces = true, bool storeAjacencies = true);
	Cell(long id, ElementType type, int connectSize, bool interior = true, bool storeInterfaces = true, bool storeAjacencies = true);
	Cell(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior = true, bool storeInterfaces = true, bool storeAjacencies = true);

	void swap(Cell &other) noexcept;

	BITPIT_DEPRECATED(void initialize(long id, ElementType type, bool interior, bool storeNeighbourhood ));
	BITPIT_DEPRECATED(void initialize(long id, ElementType type, int connectSize, bool interior, bool storeNeighbourhood));
	BITPIT_DEPRECATED(void initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeNeighbourhood));
	void initialize(long id, ElementType type, bool interior, bool storeInterfaces = true, bool storeAjacencies = true);
	void initialize(long id, ElementType type, int connectSize, bool interior, bool storeInterfaces = true, bool storeAjacencies = true);
	void initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeInterfaces = true, bool storeAjacencies = true);

	bool isInterior() const;
	
	void deleteInterfaces();
	void resetInterfaces(bool storeInterfaces = true);
	void setInterfaces(const std::vector<std::vector<long>> &interfaces);
	void setInterfaces(FlatVector2D<long> &&interfaces);
	void setInterface(int face, int index, long interface);
	bool pushInterface(int face, long interface);
	void deleteInterface(int face, int i);
	int getInterfaceCount() const;
	int getInterfaceCount(int face) const;
	long getInterface(int face, int index = 0) const;
	const long * getInterfaces() const;
	long * getInterfaces();
	const long * getInterfaces(int face) const;
	long * getInterfaces(int face);
	int findInterface(int face, int interface);
	int findInterface(int interface);

	void deleteAdjacencies();
	void resetAdjacencies(bool storeAdjacencies = true);
	void setAdjacencies(const std::vector<std::vector<long>> &adjacencies);
	void setAdjacencies(FlatVector2D<long> &&adjacencies);
	void setAdjacency(int face, int index, long adjacencies);
	bool pushAdjacency(int face, long adjacency);
	void deleteAdjacency(int face, int i);
	int getAdjacencyCount() const;
	int getAdjacencyCount(int face) const;
	long getAdjacency(int face, int index = 0) const;
	const long * getAdjacencies() const;
	long * getAdjacencies();
	const long * getAdjacencies(int face) const;
	long * getAdjacencies(int face);
	int findAdjacency(int face, int adjacency);
	int findAdjacency(int adjacency);

	bool isFaceBorder(int face) const;

	void display(std::ostream &out, unsigned short int indent) const;

	unsigned int getBinarySize() const;

protected:
	void setInterior(bool interior);

private:
	bitpit::FlatVector2D<long> m_interfaces;
	bitpit::FlatVector2D<long> m_adjacencies;

	bool m_interior;

	bitpit::FlatVector2D<long> createNeighbourhoodStorage(bool storeNeighbourhood);

	void _initialize(bool interior, bool initializeInterfaces, bool storeInterfaces, bool initializeAdjacency, bool storeAdjacencies);

};

template<typename QualifiedCell>
class QualifiedCellHalfEdge : public ElementHalfEdge<QualifiedCell> {

public:
	typedef typename ElementHalfEdge<QualifiedCell>::Winding Winding;

	QualifiedCellHalfEdge(QualifiedCell &cell, int edge, Winding winding = Winding::WINDING_NATURAL);

	QualifiedCell & getCell() const;

};

template<typename QualifiedCell>
class QualifiedCellHalfFace : public ElementHalfFace<QualifiedCell> {

public:
	typedef typename ElementHalfFace<QualifiedCell>::Winding Winding;

	QualifiedCellHalfFace(QualifiedCell &cell, int face, Winding winding = Winding::WINDING_NATURAL);

	QualifiedCell & getCell() const;

};

extern template class PiercedVector<Cell>;

extern template class QualifiedCellHalfEdge<Cell>;
extern template class QualifiedCellHalfEdge<const Cell>;

extern template class QualifiedCellHalfFace<Cell>;
extern template class QualifiedCellHalfFace<const Cell>;

typedef QualifiedCellHalfEdge<Cell> CellHalfEdge;
typedef QualifiedCellHalfEdge<const Cell> ConstCellHalfEdge;

typedef QualifiedCellHalfFace<Cell> CellHalfFace;
typedef QualifiedCellHalfFace<const Cell> ConstCellHalfFace;

}

// Include template implementations
#include "cell.tpp"

#endif
