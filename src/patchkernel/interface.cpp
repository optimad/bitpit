/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include <cmath>

#include "bitpit_common.hpp"
#include "bitpit_operators.hpp"

#include "cell.hpp"
#include "interface.hpp"

namespace bitpit {

/*!
	\class Interface
	\ingroup patchelements

	\brief The Interface class defines the interfaces among cells.

	Interface is class that defines the interfaces among cells.
*/

/*!
	Default constructor.
*/
Interface::Interface()
	: Element()
{
	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Creates a new interface.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
*/
Interface::Interface(const long &id, ElementType type)
	: Element(id, type)
{
	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Creates a new interface.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
*/
Interface::Interface(const long &id, ElementType type, int connectSize)
	: Element(id, type, connectSize)
{
	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Creates a new interface.

	\param id is the id that will be assigned to the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
Interface::Interface(const long &id, ElementType type, std::unique_ptr<long[]> &&connectStorage)
	: Element(id, type, std::move(connectStorage))
{
	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/**
* Exchanges the content of the interface by the content the specified other
* interface.
*
* \param other is another interface whose content is swapped with that of this
* interface
*/
void Interface::swap(Interface &other) noexcept
{
	Element::swap(other);

	std::swap(other.m_owner, m_owner);
	std::swap(other.m_ownerFace, m_ownerFace);
	std::swap(other.m_neigh, m_neigh);
	std::swap(other.m_neighFace, m_neighFace);
}

/*!
	Initializes the data structures of the interface.

	\param id is the id of the element
	\param type is the type of the element
*/
void Interface::initialize(long id, ElementType type)
{
	Element::initialize(id, type);

	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Initializes the data structures of the interface.

	\param id is the id of the element
	\param type is the type of the element
	\param connectSize is the size of the connectivity, this is only used
	if the element is not associated to a reference element
*/
void Interface::initialize(long id, ElementType type, int connectSize)
{
	Element::initialize(id, type, connectSize);

	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Initializes the data structures of the interface.

	\param id is the id of the element
	\param type is the type of the element
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
void Interface::initialize(long id, ElementType type, std::unique_ptr<long[]> &&connectStorage)
{
	Element::initialize(id, type, std::move(connectStorage));

	_initialize(NULL_ID, -1, NULL_ID, -1);
}

/*!
	Internal function to initialize the data structures of the interface.

	\param owner is the cell on the "owner" side of the interface
	\param ownerFace is the owner face on the interface
	\param neigh is the cell on the "neighbour" side of the interface
	\param neighFace is the neighbour face on the interface
*/
void Interface::_initialize(long owner, long ownerFace, long neigh, long neighFace)
{
	m_owner     = owner;
	m_ownerFace = ownerFace;

	m_neigh     = neigh;
	m_neighFace = neighFace;
}

/*!
	Evaluates the rotation matrix from the Cartesian coordinate system
	to a coordinate system build starting from the specified versor.

	Evaluates the rotation matrix that needs to be applied to the
	Cartesian coordinate system to make it coincide with the
	coordinates system defined starting from the specified versor.
	The axes of the coordinate system are defined as follows:

	  - the x axis is aligned with the versor;

	  - the y axis is normal to the plane where the axis x-versor
	    and z-Cartesian lay, or, if this two vectors are aligned, to the
	    plane where the axis x-versor and x-Cartesian lay;

	  - the z axis is obtained evaluating the cross product of the axis
	    x-versor and y-versor.

	\param versor is the versor that defines the coordinate system
	\result The rotation matrix from the Cartesian coordinate system
	to a coordinate system build starting from the specified versor.
*/
std::array<std::array<double, 3>, 3> Interface::evalRotationFromCartesian(std::array<double, 3> &versor)
{
	// The rotation matrix has in its rows the versors that define
	// the interface coordinate system.
	//
	//                | [x_int] |
	//          R =   | [y_int] |
	//                | [z_int] |
	//
	std::array<std::array<double, 3>, 3> R;

	// x-interface axis
	for (int k = 0; k < 3; ++k) {
		R[0][k] = versor[k];
	}

	// y-interface axis
	if (std::abs(std::abs(versor[2]) - 1.) > 1e-8) {
		std::array<double, 3> z = {{0.0, 0.0, 1.0}};
		R[1] = crossProduct(z, R[0]);
	} else {
		std::array<double, 3> x = {{1.0, 0.0, 0.0}};
		R[1] = crossProduct(x, R[0]);
	}
	R[1] = R[1] / norm2(R[1]);

	// z-interface axis
	R[2] = crossProduct(R[0], R[1]);
	R[2] = R[2] / norm2(R[2]);

	return R;
}

/*!
	Evaluates the rotation matrix from the coordinate system build
	starting from the specified versor to the Cartesian coordinate
	system.

	Evaluates the rotation matrix that needs to be applied to the
	coordinates system defined starting from the specified versor
	to make it coincide with the Cartesian coordinates system.
	This matrix can be evaluated as the inverse of the rotation
	matrix from the Cartesian coordinate system to the versor
	coordinate system.

	\param versor is the three-dimensional versor that defines the
	coordinate system
	\result The rotation matrix from the coordinate system build
	starting from the specified versor to the Cartesian coordinate
	system.
*/
std::array<std::array<double, 3>, 3> Interface::evalRotationToCartesian(std::array<double, 3> &versor)
{
	return evalRotationInverse(evalRotationFromCartesian(versor));
}

/*!
	Evaluates the inverse of the specified rotation matrix.

	Rotation matrices are orthogonal matrices, hence the inverse of a rotation
	matrix is equal to its transpose.

	\param R the rotation matrix to transpose
	\result The inverse of the rotation matrix.
*/
std::array<std::array<double, 3>, 3> Interface::evalRotationInverse(const std::array<std::array<double, 3>, 3> &R)
{
	return evalRotationTranspose(R);
}

/*!
	Evaluates the transpose of the specified rotation matrix.

	\param R the rotation matrix to transpose
	\result The transposed rotation matrix.
*/
std::array<std::array<double, 3>, 3> Interface::evalRotationTranspose(const std::array<std::array<double, 3>, 3> &R)
{
	std::array<std::array<double, 3>, 3> R_transposed;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			R_transposed[j][i] = R[i][j];
		}
	}

	return R_transposed;
}

/*!
	Checks whether the interface is a border.

	\result Returns true if the interface is a border, false otherwise.
*/
bool Interface::isBorder() const
{
	return (m_neigh < 0);
}

/*!
	Sets the owner of the interface.

	\param owner the owner of the interface
	\param onwerFace the owner's face adjacent to the interface
*/
void Interface::setOwner(const long &owner, const int &onwerFace)
{
	m_owner     = owner;
	m_ownerFace = onwerFace;
}

/*!
	Deletes the owner of the interface.
*/
void Interface::unsetOwner()
{
	m_owner     = Element::NULL_ID;
	m_ownerFace = -1;
}

/*!
	Gets the owner of the interface.

	\result The owner of the nterface
*/
long Interface::getOwner() const
{
  return m_owner;
}

/*!
	Gets the face of the owner associated with the interface.

	\result The face of the owner associated with the interface
*/
int Interface::getOwnerFace() const
{
  return m_ownerFace;
}

/*!
	Sets the neighbour of the interface.

	\param neigh the neighbour of the interface
	\param neighFace the neighbour's face adjacent to the interface
*/
void Interface::setNeigh(const long &neigh, const int &neighFace)
{
	m_neigh     = neigh;
	m_neighFace = neighFace;
}

/*!
	Deletes the neighbour of the interface.
*/
void Interface::unsetNeigh()
{
	m_neigh     = Element::NULL_ID;
	m_neighFace = -1;
}

/*!
	Gets the neighbour of the interface.

	\result The neighbour of the nterface
*/
long Interface::getNeigh() const
{
  return m_neigh;
}

/*!
	Gets the face of the neighbour associated with the interface.

	\result The face of the neighbour associated with the interface
*/
int Interface::getNeighFace() const
{
  return m_neighFace;
}

/*!
	Gets both the owner and the neighbour of the interface.

	\result An array containing the owner and the neighbour of the
	        interface.
*/
std::array<long, 2> Interface::getOwnerNeigh() const
{
	std::array<long, 2> cells;
	cells[0] = m_owner;
	cells[1] = m_neigh;

	return cells;
}

/*!
	Displays interface information to an output stream

	\param[in] out is the output stream
	\param[in] indent is the number of trailing spaces to prepend when
	writing the information
*/
void Interface::display(std::ostream &out, unsigned short int indent) const
{
	std::string t_s = std::string(indent, ' ');

	// If the type is unknown there are no information to display
	if (getType() == ElementType::UNDEFINED) {
	    out << t_s << "interface type:    (unknown)" << std::endl;
	    return;
	}

	// General info ----------------------------------------------------- //
	out << t_s << "interface type: " << getType() << std::endl;
	out << t_s << "ID:             " << getId() << std::endl;
	out << t_s << "is border:      ";
	if (getNeigh() >= 0)  { out << "(false)"; }
	else                  { out << "(true)"; }
	out << std::endl;

	// Connectivity infos --------------------------------------------------- //
	int nVertices = getVertexCount();
	ConstProxyVector<long> cellVertexIds = getVertexIds();
	out << t_s << "connectivity: [ ";
	for (int i = 0; i < nVertices - 1; ++i) {
		out << cellVertexIds[i] << ", ";
	} //next i
	out << cellVertexIds[nVertices - 1] << " ]" << std::endl;

	// Onwer infos ---------------------------------------------------------- //
	out << t_s << "onwer ID:    " << getOwner() << std::endl;
	out << t_s << "owner face:  " << getOwnerFace() << std::endl;

	// Onwer infos ---------------------------------------------------------- //
	if (getNeigh() >= 0) {
		out << t_s << "neighbour ID:    " << getNeigh() << std::endl;
		out << t_s << "neighbour face:  " << getNeighFace() << std::endl;
	}
}

/*!
	\class InterfaceHalfEdge
	\ingroup patchelements

	\brief The InterfaceHalfEdge class defines interface half-edges.

	InterfaceHalfEdge is the class that defines interface half-edges. Each
	edge can be seen as two half-edges: one belonging to an interface and
	the other belonging to the neighbouring interface. A half-edge is
	identify by its vertices and by the winding order of the vertices.
*/

/*!
	Constructor.

	\param interface is a reference to the interface the owns the edge
	\param edge if the local edge of the interface
	\param winding is the winding order of the vertices
*/
InterfaceHalfEdge::InterfaceHalfEdge(Interface &interface, int edge, Winding winding)
    : ElementHalfFace(interface, edge, winding)
{
}

/*!
	Get the interface the edge belongs to.

	\result Returns the interface the edge belongs to.
*/
Interface & InterfaceHalfEdge::getInterface() const
{
	return static_cast<Interface &>(getElement());
}

// Explicit instantiation of the Interface containers
template class PiercedVector<Interface>;

}
