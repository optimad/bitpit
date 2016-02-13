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

#include<iostream>

#include "bitpit_common.hpp"

#include "cell.hpp"
#include "interface.hpp"


/*!
	Input stream operator for class Cell. Stream cell data from memory
	input stream to container.

	\param[in] buffer is the input stream from memory
	\param[in] cell is the cell object
	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::Cell &cell)
{
	// Write connectivity data ---------------------------------------------- //
	bitpit::Element &element(cell);
	buffer >> element;

	// Write interface data ------------------------------------------------- //
	buffer >> cell.m_interfaces;

	// Write adjacencies data ----------------------------------------------- //
	buffer >> cell.m_adjacencies;

	return buffer;
}

/*!
	Output stream operator for class Cell. Stream cell data from container
	to output stream.

	\param[in] buffer is the output stream from memory
	\param[in] cell is the cell object
	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const bitpit::Cell &cell)
{
	// Write connectivity data ---------------------------------------------- //
	const bitpit::Element &element(cell);
	buffer << element;

	// Write interface data ------------------------------------------------- //
	buffer << cell.m_interfaces;

	// Write adjacencies data ----------------------------------------------- //
	buffer << cell.m_adjacencies;

	return buffer;
}

namespace bitpit {

/*!
	\ingroup patch
	@{
*/

/*!
	\class Cell

	\brief The Cell class defines the cells.

	Cell is class that defines the cells.
*/

/*!
	Default constructor.
*/
Cell::Cell()
	: Element()
{

}

/*!
	Creates a new cell.
*/
Cell::Cell(const long &id, ElementInfo::Type type, bool interior, bool storeNeighbourhood)
	: Element(id, type)
{
	_initialize(interior, storeNeighbourhood);
}

/*!
        Copy-constructor
*/
Cell::Cell(const Cell &other)
{
	*this = other;
}

/*!
	Copy assignament operator
*/
Cell & Cell::operator=(const Cell& other)
{
	Element::operator=(other);
	m_interior    = other.m_interior;
	m_interfaces  = other.m_interfaces;
	m_adjacencies = other.m_adjacencies;

	return (*this);
}

/*!
	Initializes the data structures of the cell.

	\param type is the type of the element
	\param interior if true the cell is flagged as interior
*/
void Cell::initialize(ElementInfo::Type type, bool interior, bool storeNeighbourhood)
{
	Element::initialize(type);

	_initialize(interior, storeNeighbourhood);
}

/*!
	Internal function to initialize the data structures of the cell.
*/
void Cell::_initialize(bool interior, bool storeNeighbourhood)
{
	setInterior(interior);

	resetInterfaces(storeNeighbourhood);
	resetAdjacencies(storeNeighbourhood);
}

/*!
	Sets if the cells belongs to the the interior domain.

	\param interior defines if the cells belongs to the the interior domain
*/
void Cell::setInterior(bool interior)
{
	m_interior = interior;
}

/*!
	Gets if the cells belongs to the the interior domain.

	\result Returns true if the cells belongs to the the interior domain,
	otherwise it returns false.
*/
bool Cell::isInterior() const
{
	return m_interior;
}

/*!
	Resets the interfaces of the cell.

	If the interfaces are stored, there will always be at least one
	interface entry for each face. If a face is not linked with an
	interface, its entry needs to be set to the placeholder value
	NULL_ID. When multiple interfaces are linked to a face, all
	entries must point to valid interfaces.

	The interface data structure can be prepared to store the interfaces
	only if the cell type is known.

	\param storeInterfaces if true the interface data structure will
	prepared to store the interfaces, otherwise the data structure will
	be cleared and it will not be possible to store interfaces.
*/
void Cell::resetInterfaces(bool storeInterfaces)
{
	if (!storeInterfaces || getType() == ElementInfo::UNDEFINED) {
		m_interfaces.clear();
	} else {
		std::vector<int> interfaceCount(getFaceCount(), 1);
		m_interfaces = CollapsedVector2D<long>(interfaceCount, NULL_ID);
	}
}

/*!
	Sets all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::setInterfaces(std::vector<std::vector<long>> &interfaces)
{
	m_interfaces = CollapsedVector2D<long>(interfaces);

	// The interface vector must have as many elements as faces
	if (getType() != ElementInfo::UNDEFINED) {
		int delta = m_interfaces.size() - getFaceCount();
		if (delta > 0) {
			for (int i = 0; i < delta; ++i) {
				m_interfaces.pop_back();
			}
		} else if (delta < 0) {
			for (int i = 0; i < delta; ++i) {
				m_interfaces.push_back(1, NULL_ID);
			}
		}
	}

	// Check that there is at least one interfaces for each face
	for (int i = 0; i < m_interfaces.size(); ++i) {
		if (m_interfaces.sub_array_size(i) == 0) {
			m_interfaces.push_back_in_sub_array(i, NULL_ID);
		}
	}
}

/*!
	Sets the i-th interface associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the interface
	\param interface is the index of the interface 
*/
void Cell::setInterface(const int &face, const int &index, const long &interface)
{
	m_interfaces.set(face, index, interface);
}

/*!
	Add an interface to the given face of the cell.

	\param face is the face of the cell
	\param interface is the index of the interface that will be added
*/
void Cell::pushInterface(const int &face, const long &interface)
{
	// If there is only one interface stored for the specified face and
	// that interface is negative, we need to overwrite the value and
	// not add another interface.
	if (m_interfaces.sub_array_size(face) == 1) {
		if (m_interfaces.get(face, 0) < 0) {
			m_interfaces.set(face, 0, interface);
			return;
		}
	}

	// There are multiple adjacency, we need to add the interface.
	m_interfaces.push_back_in_sub_array(face, interface);
}

/*!
	Deletes the specified interface from the interfaces associate to the
	given face of the cell.

	\param face the face of the cell
	\param i is the index of the interface to delete
*/
void Cell::deleteInterface(const int &face, const int &i)
{
	// If there is only one interface stored for the specified face, we
	// need to overwrite the value and not delete the interface.
	if (m_interfaces.sub_array_size(face) == 1) {
		m_interfaces.set(face, 0, NULL_ID);
		return;
	}

	// There are multiple interfaces, we need to delete the interface.
	m_interfaces.erase(face, i);
}

/*!
	Gets the total number of interfaces of the cell.

	The placeholder interface ids on faces not acutally linked to a
	real interfaces will be counted as well.

	\result The total number of interfaces of the cell.
*/
int Cell::getInterfaceCount() const
{
	return m_interfaces.sub_arrays_total_size();
}

/*!
	Gets the number of interfaces of the specified face of the cell.

	The placeholder interface id on a face not acutally linked to a
	real interfaces will be counted as well.

	\param face the face of the cell
	\result The number of interfaces of the specified face of the cell.
*/
int Cell::getInterfaceCount(const int &face) const
{
	return m_interfaces.sub_array_size(face);
}

/*!
	Gets the i-th interface of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface.
*/
long Cell::getInterface(const int &face, const int &index) const
{
	return m_interfaces.get(face, index);
}

/*!
	Gets all the interfaces of the cell.

	\result The interfaces of the cell.
*/
const long * Cell::getInterfaces() const
{
	return m_interfaces.get(0);
}

/*!
	Gets the interfaces of the given face of the cell.

	\param face the face of the cell
	\result The requested interfaces
*/
const long * Cell::getInterfaces(const int &face) const
{
	return m_interfaces.get(face);
}

/*!
	Resets the adjacencies of the cell.

	If the adjacencies are stored, there will always be at least one
	adjacency entry for each face. If a face is not linked with a
	neighbour, its entry needs to be set to the placeholder value
	NULL_ID. When multiple neighbours are linked to a face, all
	entries must point to valid cells.

	The adjacency data structure can be prepared to store the neighbours
	only if the cell type is known.

	\param storeAdjacencies if true the adjacency data structure will
	prepared to store the neighbours, otherwise the data structure will
	be cleared and it will not be possible to store neighbours.
*/
void Cell::resetAdjacencies(bool storeAdjacencies)
{
	if (!storeAdjacencies || getType() == ElementInfo::UNDEFINED) {
		m_adjacencies.clear();
	} else {
		std::vector<int> adjacencyCount(getFaceCount(), 1);
		m_adjacencies = CollapsedVector2D<long>(adjacencyCount, NULL_ID);
	}
}

/*!
	Sets all the adjacencies of the cell.

	\param adjacencies the list of all adjacencies associated to the cell
*/
void Cell::setAdjacencies(std::vector<std::vector<long>> &adjacencies)
{
	m_adjacencies = CollapsedVector2D<long>(adjacencies);

	// The adjacency vector must have as many elements as faces
	if (getType() != ElementInfo::UNDEFINED) {
		int delta = m_adjacencies.size() - getFaceCount();
		if (delta > 0) {
			for (int i = 0; i < delta; ++i) {
				m_adjacencies.pop_back();
			}
		} else if (delta < 0) {
			for (int i = 0; i < delta; ++i) {
				m_adjacencies.push_back(1, NULL_ID);
			}
		}
	}

	// Check that there is at least one adjacency for each face
	for (int i = 0; i < m_adjacencies.size(); ++i) {
		if (m_adjacencies.sub_array_size(i) == 0) {
			m_adjacencies.push_back_in_sub_array(i, NULL_ID);
		}
	}
}

/*!
	Sets the i-th adjacency associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the adjacency
	\param adjacency is the index of the adjacency
*/
void Cell::setAdjacency(const int &face, const int &index, const long &adjacency)
{
	m_adjacencies.set(face, index, adjacency);
}

/*!
	Add an adjacency to the given face of the cell.

	\param face is the face of the cell
	\param adjacency is the index of the adjacency that will be added
*/
void Cell::pushAdjacency(const int &face, const long &adjacency)
{
	// If there is only one adjacency stored for the specified face and
	// that adjacency is negative, we need to overwrite the value and
	// not add another adjacency.
	if (m_adjacencies.sub_array_size(face) == 1) {
		if (m_adjacencies.get(face, 0) < 0) {
			m_adjacencies.set(face, 0, adjacency);
			return;
		}
	}

	// There are multiple adjacency, we need to add the adjacency.
	m_adjacencies.push_back_in_sub_array(face, adjacency);
}

/*!
	Deletes the specified adjacency from the adjacencies associate to the
	given face of the cell.

	\param face the face of the cell
	\param i is the index of the adjacency to delete
*/
void Cell::deleteAdjacency(const int &face, const int &i)
{
	// If there is only one adjacency stored for the specified face, we
	// need to overwrite that value and not delete the adjacency.
	if (m_adjacencies.sub_array_size(face) == 1) {
		m_adjacencies.set(face, 0, NULL_ID);
		return;
	}

	// There are multiple adjacencies, we need to delete the adjacency.
	m_adjacencies.erase(face, i);
}

/*!
	Gets the total number of adjacencies of the cell.

	The placeholder neighbour ids on faces not acutally linked to a
	real neighbour will be counted as well.

	\result The total number of adjacencies of the cell.
*/
int Cell::getAdjacencyCount() const
{
	return m_adjacencies.sub_arrays_total_size();
}

/*!
	Gets the number of adjacencies of the specified face of the cell.

	The placeholder neighbour id on a face not acutally linked to a
	real neighbours will be counted as well.

	\param face the face of the cell
	\result The number of adjacencies of the specified face of the cell.
*/
int Cell::getAdjacencyCount(const int &face) const
{
	return m_adjacencies.sub_array_size(face);
}

/*!
	Gets the i-th adjacency of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the adjacency to retreive
	\result The requested adjacency.
*/
long Cell::getAdjacency(const int &face, const int &index) const
{
	return m_adjacencies.get(face, index);
}

/*!
	Gets all the adjacencies of the cell.

	\result The adjacencies of the cell.
*/
const long * Cell::getAdjacencies() const
{
	return m_adjacencies.get(0);
}

/*!
	Gets the adjacencies of the given face of the cell.

	\as getAdjacency(const int &face, const int &index) const

	\param face the face of the cell
	\result The requested adjacencies
*/
const long * Cell::getAdjacencies(const int &face) const
{
	return m_adjacencies.get(face);
}

/*!
	Displays the cell information to an output stream

	\param[in] out is the output stream
	\param[in] indent is the number of trailing spaces to prepend when
	writing the information
*/
void Cell::display(std::ostream &out, unsigned short int indent)
{
	// ====================================================================== //
	// VARIABLES DECLARATION                                                  //
	// ====================================================================== //

	// Local variables
	std::string                 t_s(indent, ' ');

	// Counters
	int                         i, j;
	int                         nn;

	// ====================================================================== //
	// DISPLAY INFOS                                                          //
	// ====================================================================== //
	if (getType() == ElementInfo::UNDEFINED) {
	    out << t_s << "cell type:    (unknown)" << std::endl;
	    return;
	}

	// Scope variables -------------------------------------------------- //
	int                         nv = getVertexCount();
	int                         nf = getFaceCount();

	// General info ----------------------------------------------------- //
	out << t_s << "cell type:    " << getType() << std::endl;
	out << t_s << "ID:           " << get_id() << std::endl;
	out << t_s << "is ghost:     ";
	if (m_interior)     { out << "(false)"; }
	else                { out << "(true)"; }
	out << std::endl;

	// Connectivity infos --------------------------------------------------- //
	out << t_s << "connectivity: [ ";
	for (i = 0; i < nv-1; ++i) {
		out << getVertex(i) << ", ";
	} //next i
	out << getVertex(nv-1) << " ]" << std::endl;

	// neighbors infos ------------------------------------------------------ //
        if (m_interfaces.size() > 0) {
            out << t_s << "neighbors:    [ ";
            for (i = 0; i < nf-1; ++i) {
                nn = getInterfaceCount(i);
                out << "[ ";
                for (j = 0; j < nn-1; ++j) {
                    out << getInterface(i, j) << ", ";
                } //next j
                out << getInterface(i, nn-1) << " ], ";
            } //next i
            nn = getInterfaceCount(nf-1);
            out << "[ ";
            for (j = 0; j < nn-1; ++j) {
                out << getInterface(nf-1, j) << ", ";
            } //next j
            out << getInterface(nf-1, nn-1) << " ]";
            out << " ]" << std::endl;
        }
}

/*!
	Get the size of the buffer required to communicate cell.

	\result Returns the buffer size (in bytes).
*/
unsigned int Cell::getBinarySize()
{
    return (Element::getBinarySize() + m_interfaces.get_binary_size() + m_adjacencies.get_binary_size());
}

// Explicit instantiation of the Cell containers
template class PiercedVector<Cell>;
template class PositionalPiercedVector<Cell>;

/*!
	@}
*/

}
