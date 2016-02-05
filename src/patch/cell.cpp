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
Cell::Cell(const long &id, ElementInfo::Type type)
	: Element(id, type)
{

}

/*!
	Initializes the data structures of the cell.

	\param type is the type of the element
	\param nInterfacesPerFace is the number of interfaces per face that
	will be used to initialize the interface's data structure. If this
	parameter is less or equal to zero, the interface's data structure
	will not be initialized.
*/
void Cell::initialize(ElementInfo::Type type, int nInterfacesPerFace)
{
	Element::initialize(type);

	if (type != ElementInfo::UNDEFINED && nInterfacesPerFace >= 1) {
		const ElementInfo &elementInfo = get_info();
		std::vector<int> interfaceCount(elementInfo.nFaces, nInterfacesPerFace);
		initializeEmptyInterfaces(interfaceCount);
	}
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
	Initialize all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::initializeInterfaces(std::vector<std::vector<long>> &interfaces)
{
	m_interfaces = CollapsedVector2D<long>(interfaces);
}

/*!
	Initialize the data structure that holds the information about the
	interfaces.

	\param nInterfaces An array with the number of interfaces for each face
*/
void Cell::initializeEmptyInterfaces(const std::vector<int> interfaceCount)
{
	m_interfaces = CollapsedVector2D<long>(interfaceCount, NULL_ELEMENT_ID);
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
	m_interfaces.erase(face, i);
}

/*!
	Unsets the interfaces associated to the cell.
*/
void Cell::unsetInterfaces()
{
	m_interfaces.clear();
}

/*!
	Gets the total number of interfaces of the cell.

	\result The total number of interfaces of the cell.
*/
int Cell::getInterfaceCount() const
{
	return m_interfaces.sub_arrays_total_size();
}

/*!
	Gets the number of interfaces of the specified face of the cell.

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

	\as getInterface(const int &face, const int &index) const

	\param face the face of the cell
	\result The requested interfaces
*/
const long * Cell::getInterfaces(const int &face) const
{
	return m_interfaces.get(face);
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
	if (getType() != ElementInfo::UNDEFINED) {
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
	out << t_s << "neighbors:   [ ";
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

/*!
	Get the size of the buffer required to communicate cell.

	\result Returns the buffer size (in bytes).
*/
unsigned int Cell::getBinarySize()
{
    return (Element::getBinarySize() + m_interfaces.get_binary_size());
}

// Explicit instantiation of the Cell containers
template class PiercedVector<Cell>;
template class PositionalPiercedVector<Cell>;

/*!
	@}
*/

}
