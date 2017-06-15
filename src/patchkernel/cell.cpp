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

#include<iostream>

#include "bitpit_common.hpp"

#include "vertex.hpp"
#include "cell.hpp"

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
	\class Cell
	\ingroup patchelements

	\brief The Cell class defines the cells.

	Cell is class that defines the cells.
*/

/*!
	Default constructor.
*/
Cell::Cell()
	: Element(), m_interior(true), m_pid(0),
      m_interfaces(false), m_adjacencies(false)
{

}

/*!
	Creates a new cell.
*/
Cell::Cell(const long &id, ElementInfo::Type type, bool interior, bool storeNeighbourhood)
	: Element(id, type), m_pid(0),
      m_interfaces(storeNeighbourhood ? ElementInfo(type).nFaces : 0, 1, NULL_ID),
      m_adjacencies(storeNeighbourhood ? ElementInfo(type).nFaces : 0, 1, NULL_ID)
{
	_initialize(interior, false, false);
}

/*!
        Copy-constructor
*/
Cell::Cell(const Cell &other)
	: Element(other)
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
	m_pid         = other.m_pid;
	m_interfaces  = other.m_interfaces;
	m_adjacencies = other.m_adjacencies;

	return (*this);
}

/**
* Exchanges the content of the cell by the content the specified other cell.
*
* \param other is another cell whose content is swapped with that of this cell.
*/
void Cell::swap(Cell &other) noexcept
{
	Element::swap(other);

	std::swap(other.m_interior, m_interior);
	std::swap(other.m_pid, m_pid);

	other.m_interfaces.swap(m_interfaces);
	other.m_adjacencies.swap(m_adjacencies);
}

/*!
	Initializes the data structures of the cell.

	\param type is the type of the element
	\param interior if true the cell is flagged as interior
	\param storeNeighbourhood if true the structures to store adjacencies
	and interfaces will be initialized
*/
void Cell::initialize(ElementInfo::Type type, bool interior, bool storeNeighbourhood)
{
	Element::initialize(type);

	_initialize(interior, true, storeNeighbourhood);
}

/*!
	Internal function to initialize the data structures of the cell.
*/
void Cell::_initialize(bool interior, bool initializeNeighbourhood, bool storeNeighbourhood)
{
	setInterior(interior);

	// Neighbourhood
	if (initializeNeighbourhood) {
		// To reduce memory fragmentation, destroy the interfaces/adjacencies
		// before resetting the adjacencies
		m_interfaces.destroy();
		m_adjacencies.destroy();

		resetInterfaces(storeNeighbourhood);
		resetAdjacencies(storeNeighbourhood);
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
	Sets the PID associated to the cell.

	\param pid is the PID associated to the cell.
*/
void Cell::setPID(int pid)
{
	m_pid = pid;
}

/*!
	Gets the PID associated to the cell.

	\result The PID associated to the cell.
*/
int Cell::getPID() const
{
	return m_pid;
}

/*!
	Deletes the interfaces of the cell.
*/
void Cell::deleteInterfaces()
{
	m_interfaces.clear();
	m_interfaces.shrinkToFit();
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
		deleteInterfaces();
	} else {
		m_interfaces.initialize(getFaceCount(), 1, NULL_ID);
	}
}

/*!
	Sets all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::setInterfaces(std::vector<std::vector<long>> &interfaces)
{
	m_interfaces.initialize(interfaces);

	// The interface vector must have as many elements as faces
	if (getType() != ElementInfo::UNDEFINED) {
		int delta = m_interfaces.size() - getFaceCount();
		if (delta > 0) {
			for (int i = 0; i < delta; ++i) {
				m_interfaces.popBack();
			}
		} else if (delta < 0) {
			for (int i = 0; i < delta; ++i) {
				m_interfaces.pushBack(1, NULL_ID);
			}
		}
	}

	// Check that there is at least one interfaces for each face
	for (int i = 0; i < m_interfaces.size(); ++i) {
		if (m_interfaces.getItemCount(i) == 0) {
			m_interfaces.pushBackItem(i, NULL_ID);
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
	m_interfaces.setItem(face, index, interface);
}

/*!
	Add an interface to the given face of the cell.

	\param face is the face of the cell
	\param interface is the index of the interface that will be added
*/
void Cell::pushInterface(const int &face, const long &interface)
{
	// Do not push an existing interface
	if (findInterface(face, interface) >= 0) {
		return;
	}

	// If there is only one interface stored for the specified face and
	// that interface is negative, we need to overwrite the value and
	// not add another interface.
	if (m_interfaces.getItemCount(face) == 1) {
		if (m_interfaces.getItem(face, 0) < 0) {
			m_interfaces.setItem(face, 0, interface);
			return;
		}
	}

	// There are multiple adjacency, we need to add the interface.
	m_interfaces.pushBackItem(face, interface);
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
	if (m_interfaces.getItemCount(face) == 1) {
		m_interfaces.setItem(face, 0, NULL_ID);
		return;
	}

	// There are multiple interfaces, we need to delete the interface.
	m_interfaces.eraseItem(face, i);
}

/*!
	Gets the total number of interfaces of the cell.

	The placeholder interface ids of the faces not acutally linked to a
	real interfaces will be counted as well.

	\result The total number of interfaces of the cell.
*/
int Cell::getInterfaceCount() const
{
	return m_interfaces.getItemCount();
}

/*!
	Gets the number of interfaces of the specified face of the cell.

	The placeholder interface ids of the faces not acutally linked to a
	real interfaces will be counted as well.

	\param face the face of the cell
	\result The number of interfaces of the specified face of the cell.
*/
int Cell::getInterfaceCount(const int &face) const
{
	return m_interfaces.getItemCount(face);
}

/*!
	Gets the i-th interface of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface.
*/
long Cell::getInterface(const int &face, const int &index) const
{
	return m_interfaces.getItem(face, index);
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
	Gets a pointer to the interfaces of the cell.

	\result A pointer to the interfaces of the cell.
*/
long * Cell::getInterfaces()
{
	return m_interfaces.get(0);
}

/*!
	Gets a pointer to the interfaces of the given face of the cell.

	\param face the face of the cell
	\result A pointer to the interfaces of the given face of the cell.
*/
long * Cell::getInterfaces(const int &face)
{
	return m_interfaces.get(face);
}

/*!
	Find the specified interface among the interfaces of the given face.

	The function returns the position in the interface face list of the
	specified interface. If the face doesn't contain the specified interface
	the function returns a dummy position equal to -1.

	\param face is the face of the cell
	\param interface is the interface to look for
	\result The position in the interface face list of the specfied interface.
*/
int Cell::findInterface(const int &face, const int &interface)
{
	int nFaceInterfaces = getInterfaceCount(face);
	for (int i = 0; i < nFaceInterfaces; i++) {
		if (getInterface(face, i) == interface) {
			return i;
		}
	}

	return -1;
}

/*!
	Find the specified interface among all the interfaces of the cell.

	The function returns the first position in the interface cell list of the
	specified interface. If the cell doesn't contain the specified interface
	the function returns a dummy position equal to -1.

	\param interface is the interface to look for
	\result The position in the interface cell list of the specfied interface.
*/
int Cell::findInterface(const int &interface)
{
	int nCellInterfaces = getInterfaceCount();
	const long *interfaces = getInterfaces();
	for (int i = 0; i < nCellInterfaces; i++) {
		if (interfaces[i] == interface) {
			return i;
		}
	}

	return -1;
}

/*!
	Deletes the adjacencies of the cell.
*/
void Cell::deleteAdjacencies()
{
	m_adjacencies.clear();
	m_adjacencies.shrinkToFit();
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
		deleteAdjacencies();
	} else {
		m_adjacencies.initialize(getFaceCount(), 1, NULL_ID);
	}
}

/*!
	Sets all the adjacencies of the cell.

	\param adjacencies the list of all adjacencies associated to the cell
*/
void Cell::setAdjacencies(std::vector<std::vector<long>> &adjacencies)
{
	m_adjacencies.initialize(adjacencies);

	// The adjacency vector must have as many elements as faces
	if (getType() != ElementInfo::UNDEFINED) {
		int delta = m_adjacencies.size() - getFaceCount();
		if (delta > 0) {
			for (int i = 0; i < delta; ++i) {
				m_adjacencies.popBack();
			}
		} else if (delta < 0) {
			for (int i = 0; i < delta; ++i) {
				m_adjacencies.pushBack(1, NULL_ID);
			}
		}
	}

	// Check that there is at least one adjacency for each face
	for (int i = 0; i < m_adjacencies.size(); ++i) {
		if (m_adjacencies.getItemCount(i) == 0) {
			m_adjacencies.pushBackItem(i, NULL_ID);
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
	m_adjacencies.setItem(face, index, adjacency);
}

/*!
	Add an adjacency to the given face of the cell.

	\param face is the face of the cell
	\param adjacency is the index of the adjacency that will be added
*/
void Cell::pushAdjacency(const int &face, const long &adjacency)
{
	// Do not push an existing adjacency
	if (findAdjacency(face, adjacency) >= 0) {
		return;
	}

	// If there is only one adjacency stored for the specified face and
	// that adjacency is negative, we need to overwrite the value and
	// not add another adjacency.
	if (m_adjacencies.getItemCount(face) == 1) {
		if (m_adjacencies.getItem(face, 0) < 0) {
			m_adjacencies.setItem(face, 0, adjacency);
			return;
		}
	}

	// There are multiple adjacency, we need to add the adjacency.
	m_adjacencies.pushBackItem(face, adjacency);
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
	if (m_adjacencies.getItemCount(face) == 1) {
		m_adjacencies.setItem(face, 0, NULL_ID);
		return;
	}

	// There are multiple adjacencies, we need to delete the adjacency.
	m_adjacencies.eraseItem(face, i);
}

/*!
	Gets the total number of adjacencies of the cell.

	The placeholder neighbour ids of the faces not acutally linked to a
	real neighbour will be counted as well.

	\result The total number of adjacencies of the cell.
*/
int Cell::getAdjacencyCount() const
{
	return m_adjacencies.getItemCount();
}

/*!
	Gets the number of adjacencies of the specified face of the cell.

	The placeholder neighbour ids of the faces not acutally linked to a
	real neighbours will be counted as well.

	\param face the face of the cell
	\result The number of adjacencies of the specified face of the cell.
*/
int Cell::getAdjacencyCount(const int &face) const
{
	return m_adjacencies.getItemCount(face);
}

/*!
	Gets the i-th adjacency of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the adjacency to retreive
	\result The requested adjacency.
*/
long Cell::getAdjacency(const int &face, const int &index) const
{
	return m_adjacencies.getItem(face, index);
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

	\param face the face of the cell
	\result The requested adjacencies
*/
const long * Cell::getAdjacencies(const int &face) const
{
	return m_adjacencies.get(face);
}

/*!
	Gets a pointer to the adjacencies of the cell.

	\result A pointer to the adjacencies of the cell.
*/
long * Cell::getAdjacencies()
{
	return m_adjacencies.get(0);
}

/*!
	Gets a pointer to the adjacencies of the given face of the cell.

	\param face the face of the cell
	\result A pointer to the adjacencies of the given face of the cell.
*/
long * Cell::getAdjacencies(const int &face)
{
	return m_adjacencies.get(face);
}

/*!
	Find the specified adjacency among the adjacencies of the given face.

	The function returns the position in the adjacency face list of the
	specified adjacency. If the face doesn't contain the specified adjacency
	the function returns a dummy position equal to -1.

	\param face is the face of the cell
	\param adjacency is the adjacency to look for
	\result The position in the adjacency face list of the specfied adjacency.
*/
int Cell::findAdjacency(const int &face, const int &adjacency)
{
	int nFaceAdjacencies = getAdjacencyCount(face);
	for (int i = 0; i < nFaceAdjacencies; i++) {
		if (getAdjacency(face, i) == adjacency) {
			return i;
		}
	}

	return -1;
}

/*!
        Given the global index of a vertex, returns the local index of that
        vertex within the cell (if exists). If the specified vertex does not
        exist in the cell connectivity list, returns Vertex::NULL_ID

        \param[in] vertex vertex global ID
        \result local index of the vertex with specified ID within the cell.
*/
int Cell::findVertex(const long &vertex) const
{
    int         loc_id;
    int         n_vert = getVertexCount();
    const long* c_ = getConnect();

    loc_id = std::find(c_, c_ + n_vert, vertex) - c_;
    if (loc_id >= n_vert) return(Vertex::NULL_ID);
    return (loc_id);
}

/*!
	Find the specified adjacency among all the adjacencies of the cell.

	The function returns the first position in the adjacency cell list of the
	specified adjacency. If the cell doesn't contain the specified adjacency
	the function returns a dummy position equal to -1.

	\param adjacency is the adjacency to look for
	\result The position in the adjacency cell list of the specfied adjacency.
*/
int Cell::findAdjacency(const int &adjacency)
{
	int nCellAdjacencies = getAdjacencyCount();
	const long *adjacencies = getAdjacencies();
	for (int i = 0; i < nCellAdjacencies; i++) {
		if (adjacencies[i] == adjacency) {
			return i;
		}
	}

	return -1;
}

/*!
	Checks if the specified face is a border.

	\return True if the face is a border, false otherwise.
*/
bool Cell::isFaceBorder(int face) const
{
	assert(m_adjacencies.size() > 0);

	if (m_adjacencies.size() == 0) {
		return true;
	} else {
		return (getAdjacency(face, 0) < 0);
	}
}

/*!
	Displays the cell information to an output stream

	\param[in] out is the output stream
	\param[in] indent is the number of trailing spaces to prepend when
	writing the information
*/
void Cell::display(std::ostream &out, unsigned short int indent) const
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
	out << t_s << "ID:           " << getId() << std::endl;
	out << t_s << "is ghost:     ";
	if (m_interior)     { out << "(false)"; }
	else                { out << "(true)"; }
	out << std::endl;

	// Connectivity infos --------------------------------------------------- //
	out << t_s << "connectivity: [ ";
	for (i = 0; i < nv-1; ++i) {
            if (getVertex(i) == Vertex::NULL_ID)   out << "n.a. ";
            else                                   out << getVertex(i) << ", ";
	} //next i
	if (getVertex(nv-1) == Vertex::NULL_ID)    out << "n.a. ";
	else                                       out << getVertex(nv-1) << " ]" << std::endl;

	// neighbors infos ------------------------------------------------------ //
        if (m_adjacencies.size() > 0) {
            out << t_s << "neighbors:    [ ";
            for (i = 0; i < nf-1; ++i) {
                nn = getAdjacencyCount(i);
                out << "[ ";
                for (j = 0; j < nn-1; ++j) {
                    if (getAdjacency(i,j) == Element::NULL_ID)  out << "n.a. ";
                    else                                        out << getAdjacency(i, j) << ", ";
                } //next j
                if (getAdjacency(i, nn-1) == Element::NULL_ID)  out << "n.a. ], ";
                else                                            out << getAdjacency(i, nn-1) << " ], ";
            } //next i
            nn = getAdjacencyCount(nf-1);
            out << "[ ";
            for (j = 0; j < nn-1; ++j) {
                if (getAdjacency(nf-1, j) == Element::NULL_ID)  out << "n.a. ";
                else                                            out << getAdjacency(nf-1, j) << ", ";
            } //next j
            if (getAdjacency(nf-1, nn-1) == Element::NULL_ID)   out << "n.a. ]";
            else                                                out << getAdjacency(nf-1, nn-1) << " ]";
            out << " ]" << std::endl;
        }

        // interface infos ------------------------------------------------------ //
        if (m_interfaces.size() > 0) {
            out << t_s << "interfaces:   [ ";
            for (i = 0; i < nf-1; ++i) {
                nn = getInterfaceCount(i);
                out << "[ ";
                for (j = 0; j < nn-1; ++j) {
                    if (getInterface(i,j) == Element::NULL_ID)  out << "n.a. ";
                    else                                        out << getInterface(i, j) << ", ";
                } //next j
                if (getInterface(i, nn-1) == Element::NULL_ID)  out << "n.a. ], ";
                else                                            out << getInterface(i, nn-1) << " ], ";
            } //next i
            nn = getInterfaceCount(nf-1);
            out << "[ ";
            for (j = 0; j < nn-1; ++j) {
                if (getInterface(nf-1, j) == Element::NULL_ID)  out << "n.a. ";
                else                                            out << getInterface(nf-1, j) << ", ";
            } //next j
            if (getInterface(nf-1, nn-1) == Element::NULL_ID)   out << "n.a. ]";
            else                                                out << getInterface(nf-1, nn-1) << " ]";
            out << " ]" << std::endl;
        }

}

/*!
	Get the size of the buffer required to communicate cell.

	\result Returns the buffer size (in bytes).
*/
unsigned int Cell::getBinarySize()
{
    return (Element::getBinarySize() + m_interfaces.getBinarySize() + m_adjacencies.getBinarySize());
}

// Explicit instantiation of the Cell containers
template class PiercedVector<Cell>;

}
