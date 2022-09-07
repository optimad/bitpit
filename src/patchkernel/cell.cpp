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

#include<iostream>

#include "bitpit_common.hpp"

#include "vertex.hpp"
#include "cell.hpp"

namespace bitpit {

/*!
	Input stream operator for class Cell. Stream cell data from memory
	input stream to container.

	\param[in] buffer is the input stream from memory
	\param[in] cell is the cell object
	\result Returns the same input stream received in input.
*/
IBinaryStream& operator>>(IBinaryStream &buffer, Cell &cell)
{
	// Write connectivity data ---------------------------------------------- //
	Element &element(cell);
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
OBinaryStream& operator<<(OBinaryStream  &buffer, const Cell &cell)
{
	// Write connectivity data ---------------------------------------------- //
	const Element &element(cell);
	buffer << element;

	// Write interface data ------------------------------------------------- //
	buffer << cell.m_interfaces;

	// Write adjacencies data ----------------------------------------------- //
	buffer << cell.m_adjacencies;

	return buffer;
}

/*!
	\class Cell
	\ingroup patchelements

	\brief The Cell class defines the cells.

	Cell is class that defines the cells.
*/

/*!
	Create an empty storage for interfaces or adjacencies.

	This function is called in the initializer list of the constructor, care
	must be taken to ensure that all members needed by the function are
	already initialized at the time the function is called.

	\param storeNeighbourhood defines if the cell should store neighbourhood
	information
*/
bitpit::FlatVector2D<long> Cell::createNeighbourhoodStorage(bool storeNeighbourhood)
{
	ElementType type = getType();
	if (!storeNeighbourhood || type == ElementType::UNDEFINED) {
		return bitpit::FlatVector2D<long>(false);
	}

	int nFaces = getFaceCount();
	if (nFaces <= 0) {
		return bitpit::FlatVector2D<long>(false);
	}

	return bitpit::FlatVector2D<long>(nFaces, 0);
}

/*!
	Default constructor.
*/
Cell::Cell()
	: Element(),
      m_interfaces(createNeighbourhoodStorage(false)),
      m_adjacencies(createNeighbourhoodStorage(false)),
      m_interior(true)
{

}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param interior defines if the cell is interior or ghost
	\param storeNeighbourhood defines if the cell should store neighbourhood
	information
*/
Cell::Cell(long id, ElementType type, bool interior, bool storeNeighbourhood)
	: Cell(id, type, interior, storeNeighbourhood, storeNeighbourhood)
{
}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
	\param interior defines if the cell is interior or ghost
	\param storeNeighbourhood defines if the cell should store neighbourhood
	information
*/
Cell::Cell(long id, ElementType type, int connectSize, bool interior, bool storeNeighbourhood)
	: Cell(id, type, connectSize, interior, storeNeighbourhood, storeNeighbourhood)
{
}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param interior defines if the cell is interior or ghost
	\param storeNeighbourhood defines if the cell should store neighbourhood
	information
*/
Cell::Cell(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeNeighbourhood)
	: Cell(id, type, std::move(connectStorage), interior, storeNeighbourhood, storeNeighbourhood)
{
}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param interior defines if the cell is interior or ghost
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
Cell::Cell(long id, ElementType type, bool interior, bool storeInterfaces, bool storeAdjacencies)
	: Element(id, type),
      m_interfaces(createNeighbourhoodStorage(storeInterfaces)),
      m_adjacencies(createNeighbourhoodStorage(storeAdjacencies))
{
	_initialize(interior, false, false, false, false);
}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
	\param interior defines if the cell is interior or ghost
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
Cell::Cell(long id, ElementType type, int connectSize, bool interior, bool storeInterfaces, bool storeAdjacencies)
	: Element(id, type, connectSize),
      m_interfaces(createNeighbourhoodStorage(storeInterfaces)),
      m_adjacencies(createNeighbourhoodStorage(storeAdjacencies))
{
	_initialize(interior, false, false, false, false);
}

/*!
	Creates a new cell.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param interior defines if the cell is interior or ghost
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
Cell::Cell(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeInterfaces, bool storeAdjacencies)
	: Element(id, type, std::move(connectStorage)),
      m_interfaces(createNeighbourhoodStorage(storeInterfaces)),
      m_adjacencies(createNeighbourhoodStorage(storeAdjacencies))
{
	_initialize(interior, false, false, false, false);
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

	other.m_interfaces.swap(m_interfaces);
	other.m_adjacencies.swap(m_adjacencies);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param interior if true the cell is flagged as interior
	\param storeNeighbourhood if true the structures to store adjacencies
	and interfaces will be initialized
*/
void Cell::initialize(long id, ElementType type, bool interior, bool storeNeighbourhood)
{
	initialize(id, type, interior, storeNeighbourhood, storeNeighbourhood);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
	\param interior if true the cell is flagged as interior
	\param storeNeighbourhood defines if the structures to store adjacencies
	and interfaces will be initialized
*/
void Cell::initialize(long id, ElementType type, int connectSize, bool interior, bool storeNeighbourhood)
{
	initialize(id, type, connectSize, interior, storeNeighbourhood, storeNeighbourhood);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param interior if true the cell is flagged as interior
	\param storeNeighbourhood defines if the cell should store neighbourhood
	information
*/
void Cell::initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeNeighbourhood)
{
	initialize(id, type, std::move(connectStorage), interior, storeNeighbourhood, storeNeighbourhood);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param interior if true the cell is flagged as interior
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
void Cell::initialize(long id, ElementType type, bool interior, bool storeInterfaces, bool storeAdjacencies)
{
	Element::initialize(id, type);

	_initialize(interior, true, storeInterfaces, true, storeAdjacencies);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
	\param interior if true the cell is flagged as interior
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
void Cell::initialize(long id, ElementType type, int connectSize, bool interior, bool storeInterfaces, bool storeAdjacencies)
{
	Element::initialize(id, type, connectSize);

	_initialize(interior, true, storeInterfaces, true, storeAdjacencies);
}

/*!
	Initializes the data structures of the cell.

	\param id is the id of the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param interior if true the cell is flagged as interior
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
void Cell::initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage, bool interior, bool storeInterfaces, bool storeAdjacencies)
{
	Element::initialize(id, type, std::move(connectStorage));

	_initialize(interior, true, storeInterfaces, true, storeAdjacencies);
}

/*!
	Internal function to initialize the data structures of the cell.

	\param interior if true the cell is flagged as interior
	\param initializeInterfaces defines if interface information will
	be initialized
	\param storeInterfaces defines if the cell should initialize the storage
	for storing interface information
	\param initializeAdjacency defines if adjacency information will
	be initialized
	\param storeAdjacencies defines if the cell should initialize the storage
	for storing adjacency information
*/
void Cell::_initialize(bool interior, bool initializeInterfaces, bool storeInterfaces, bool initializeAdjacency, bool storeAdjacencies)
{
	setInterior(interior);

	// Interface information
	if (initializeInterfaces) {
		resetInterfaces(storeInterfaces);
	}

	// Interface information
	if (initializeAdjacency) {
		resetAdjacencies(storeAdjacencies);
	}
}

/*!
	Sets if the cell belongs to the the interior domain.

	\param interior defines if the cell belongs to the the interior domain
*/
void Cell::setInterior(bool interior)
{
	m_interior = interior;
}

/*!
	Gets if the cell belongs to the the interior domain.

	\result Returns true if the cell belongs to the the interior domain,
	otherwise it returns false.
*/
bool Cell::isInterior() const
{
	return m_interior;
}

/*!
	Deletes the interfaces of the cell.
*/
void Cell::deleteInterfaces()
{
	resetInterfaces(false);
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
	// Early return if no interfaces should be stored
	ElementType type = getType();
	if (!storeInterfaces || type == ElementType::UNDEFINED) {
		m_interfaces.destroy();
		return;
	}

	// Early return if not element defines no faces
	int nFaces = getFaceCount();
	if (nFaces <= 0) {
		m_interfaces.destroy();
		return;
	}

	// Initialize the interface storage
	m_interfaces.initialize(nFaces, 0);
}

/*!
	Sets all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::setInterfaces(const std::vector<std::vector<long>> &interfaces)
{
	if (getType() == ElementType::UNDEFINED) {
	    return;
	}

	assert((int) interfaces.size() == getFaceCount());
	assert((int) m_interfaces.size() == getFaceCount());
	m_interfaces.initialize(interfaces);
}

/*!
	Sets all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::setInterfaces(FlatVector2D<long> &&interfaces)
{
	if (getType() == ElementType::UNDEFINED) {
	    return;
	}

	assert((int) m_interfaces.size() == getFaceCount());
	assert((int) interfaces.size() == getFaceCount());
	m_interfaces.swap(interfaces);
}

/*!
	Sets the i-th interface associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the interface
	\param interface is the index of the interface 
*/
void Cell::setInterface(int face, int index, long interface)
{
	m_interfaces.setItem(face, index, interface);
}

/*!
	Add an interface to the given face of the cell.

	Before adding the interface a check is performed to avoid adding duplicate
	adjacencies. The interface will not be added if the face already contains
	it.

	\param face is the face of the cell
	\param interface is the index of the interface that will be added
	\result Returns true if the interface has been added, false if the
	interface has not been added because the face already contained it.
*/
bool Cell::pushInterface(int face, long interface)
{
	// Do not push an existing interface
	if (findInterface(face, interface) >= 0) {
		return false;
	}

	// Add the interface
	m_interfaces.pushBackItem(face, interface);

	// Done
	return true;
}

/*!
	Deletes the specified interface from the interfaces associate to the
	given face of the cell.

	\param face the face of the cell
	\param i is the index of the interface to delete
*/
void Cell::deleteInterface(int face, int i)
{
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
int Cell::getInterfaceCount(int face) const
{
	return m_interfaces.getItemCount(face);
}

/*!
	Gets the i-th interface of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface.
*/
long Cell::getInterface(int face, int index) const
{
	return m_interfaces.getItem(face, index);
}

/*!
	Gets all the interfaces of the cell.

	\result The interfaces of the cell.
*/
const long * Cell::getInterfaces() const
{
	if (m_interfaces.empty()) {
		return nullptr;
	}

	return m_interfaces.get(0);
}

/*!
	Gets the interfaces of the given face of the cell.

	\param face the face of the cell
	\result The requested interfaces
*/
const long * Cell::getInterfaces(int face) const
{
	if (m_interfaces.empty()) {
		return nullptr;
	}

	return m_interfaces.get(face);
}

/*!
	Gets a pointer to the interfaces of the cell.

	\result A pointer to the interfaces of the cell.
*/
long * Cell::getInterfaces()
{
	if (m_interfaces.empty()) {
		return nullptr;
	}

	return m_interfaces.get(0);
}

/*!
	Gets a pointer to the interfaces of the given face of the cell.

	\param face the face of the cell
	\result A pointer to the interfaces of the given face of the cell.
*/
long * Cell::getInterfaces(int face)
{
	if (m_interfaces.empty()) {
		return nullptr;
	}

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
int Cell::findInterface(int face, int interface)
{
	long *faceInterfaces = getInterfaces(face);
	int nFaceInterfaces = getInterfaceCount(face);
	for (int i = 0; i < nFaceInterfaces; i++) {
		if (faceInterfaces[i] == interface) {
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
int Cell::findInterface(int interface)
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
	resetAdjacencies(false);
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
	// Early return if no adjacencies should be stored
	ElementType type = getType();
	if (!storeAdjacencies || type == ElementType::UNDEFINED) {
		m_adjacencies.destroy();
		return;
	}

	// Early return if not element defines no faces
	int nFaces = getFaceCount();
	if (nFaces <= 0) {
		m_adjacencies.destroy();
		return;
	}

	// Initialize the interface storage
	m_adjacencies.initialize(nFaces, 0);
}

/*!
	Sets all the adjacencies of the cell.

	\param adjacencies the list of all adjacencies associated to the cell
*/
void Cell::setAdjacencies(const std::vector<std::vector<long>> &adjacencies)
{
	if (getType() == ElementType::UNDEFINED) {
	    return;
	}

	assert((int) m_adjacencies.size() == getFaceCount());
	assert((int) adjacencies.size() == getFaceCount());
	m_adjacencies.initialize(adjacencies);
}

/*!
	Sets all the adjacencies of the cell.

	\param adjacencies the list of all adjacencies associated to the cell
*/
void Cell::setAdjacencies(FlatVector2D<long> &&adjacencies)
{
	if (getType() == ElementType::UNDEFINED) {
	    return;
	}

	assert((int) m_adjacencies.size() == getFaceCount());
	assert((int) adjacencies.size() == getFaceCount());
	m_adjacencies.swap(adjacencies);
}

/*!
	Sets the i-th adjacency associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the adjacency
	\param adjacency is the index of the adjacency
*/
void Cell::setAdjacency(int face, int index, long adjacency)
{
	m_adjacencies.setItem(face, index, adjacency);
}

/*!
	Add an adjacency to the given face of the cell.

	Before adding the adjacency a check is performed to avoid adding duplicate
	adjacencies. The adjacency will not be added if the face already contains
	it.

	\param face is the face of the cell
	\param adjacency is the index of the adjacency that will be added
	\result Returns true if the adjacency has been added, false if the
	adjacency has not been added because the face already contained it.
*/
bool Cell::pushAdjacency(int face, long adjacency)
{
	// Do not push an existing adjacency
	if (findAdjacency(face, adjacency) >= 0) {
		return false;
	}

	// Add the adjacency
	m_adjacencies.pushBackItem(face, adjacency);

	// Done
	return true;
}

/*!
	Deletes the specified adjacency from the adjacencies associate to the
	given face of the cell.

	\param face the face of the cell
	\param i is the index of the adjacency to delete
*/
void Cell::deleteAdjacency(int face, int i)
{
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
int Cell::getAdjacencyCount(int face) const
{
	return m_adjacencies.getItemCount(face);
}

/*!
	Gets the i-th adjacency of the specified face of the cell.

	\param face the face of the cell
	\param index the index of the adjacency to retreive
	\result The requested adjacency.
*/
long Cell::getAdjacency(int face, int index) const
{
	return m_adjacencies.getItem(face, index);
}

/*!
	Gets all the adjacencies of the cell.

	\result The adjacencies of the cell.
*/
const long * Cell::getAdjacencies() const
{
	if (m_adjacencies.empty()) {
		return nullptr;
	}

	return m_adjacencies.get(0);
}

/*!
	Gets the adjacencies of the given face of the cell.

	\param face the face of the cell
	\result The requested adjacencies
*/
const long * Cell::getAdjacencies(int face) const
{
	if (m_adjacencies.empty()) {
		return nullptr;
	}

	return m_adjacencies.get(face);
}

/*!
	Gets a pointer to the adjacencies of the cell.

	\result A pointer to the adjacencies of the cell.
*/
long * Cell::getAdjacencies()
{
	if (m_adjacencies.empty()) {
		return nullptr;
	}

	return m_adjacencies.get(0);
}

/*!
	Gets a pointer to the adjacencies of the given face of the cell.

	\param face the face of the cell
	\result A pointer to the adjacencies of the given face of the cell.
*/
long * Cell::getAdjacencies(int face)
{
	if (m_adjacencies.empty()) {
		return nullptr;
	}

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
int Cell::findAdjacency(int face, int adjacency)
{
	long *faceAdjacencies = getAdjacencies(face);
	int nFaceAdjacencies = getAdjacencyCount(face);
	for (int i = 0; i < nFaceAdjacencies; i++) {
		if (faceAdjacencies[i] == adjacency) {
			return i;
		}
	}

	return -1;
}

/*!
	Find the specified adjacency among all the adjacencies of the cell.

	The function returns the first position in the adjacency cell list of the
	specified adjacency. If the cell doesn't contain the specified adjacency
	the function returns a dummy position equal to -1.

	\param adjacency is the adjacency to look for
	\result The position in the adjacency cell list of the specfied adjacency.
*/
int Cell::findAdjacency(int adjacency)
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
	return (m_adjacencies.getItemCount(face) == 0);
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
	if (getType() == ElementType::UNDEFINED) {
	    out << t_s << "cell type:    (unknown)" << std::endl;
	    return;
	}

	// Scope variables -------------------------------------------------- //
	int                         nv = getVertexCount();
	int                         nf = getFaceCount();

	ConstProxyVector<long>      cellVertexIds = getVertexIds();

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
            if (cellVertexIds[i] == Vertex::NULL_ID)   out << "n.a. ";
            else                                   out << cellVertexIds[i] << ", ";
	} //next i
	if (cellVertexIds[nv-1] == Vertex::NULL_ID)    out << "n.a. ";
	else                                       out << cellVertexIds[nv-1] << " ]" << std::endl;

	// neighbors infos ------------------------------------------------------ //
        if (m_adjacencies.size() > 0) {
            out << t_s << "neighbors:    [ ";
            for (i = 0; i < nf; ++i) {
                nn = getAdjacencyCount(i);
                out << "[ ";
                if (nn == 0) {
                    out << "n.a. ";
                } else {
                    for (j = 0; j < nn; ++j) {
                        out << getAdjacency(i, j);
                        if (j != (nn - 1)) {
                            out << ",";
                        }
                        out << " ";
                    }
                }
                out << "]";
                if (i != (nf - 1)) {
                    out << ",";
                }
                out << " ";
            }
            out << "]" << std::endl;
        }

        // interface infos ------------------------------------------------------ //
        if (m_interfaces.size() > 0) {
            out << t_s << "interfaces:   [ ";
            for (i = 0; i < nf; ++i) {
                nn = getInterfaceCount(i);
                out << "[ ";
                if (nn == 0) {
                    out << "n.a. ";
                } else {
                    for (j = 0; j < nn; ++j) {
                        out << getInterface(i, j);
                        if (j != (nn - 1)) {
                            out << ",";
                        }
                        out << " ";
                    }
                }
                out << "]";
                if (i != (nf - 1)) {
                    out << ",";
                }
                out << " ";
            }
            out << "]" << std::endl;
        }

}

/*!
	Get the size of the buffer required to communicate cell.

	\result Returns the buffer size (in bytes).
*/
unsigned int Cell::getBinarySize() const
{
    return (Element::getBinarySize() + m_interfaces.getBinarySize() + m_adjacencies.getBinarySize());
}

// Explicit instantiation of the Cell containers
template class PiercedVector<Cell>;

template class QualifiedCellHalfEdge<Cell>;
template class QualifiedCellHalfEdge<const Cell>;

template class QualifiedCellHalfFace<Cell>;
template class QualifiedCellHalfFace<const Cell>;

}
