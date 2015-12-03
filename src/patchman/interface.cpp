//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <math.h>

#include "cell.hpp"
#include "interface.hpp"
#include "utils.hpp"

namespace pman {

/*!
	\class Interface

	\brief The Interface class defines the interfaces among cells.

	Interface is class that defines the interfaces among cells.
*/



/*!
	\enum Interface::Position

	This enum defines the element types that can be used.

	\var Interface::Type Interface::INTERNAL
	The interface is internal to the domain.

	\var Interface::Type Interface::BOUNDARY
	The interface is on a boundary.

	\var Interface::Type Interface::GHOST
	The interface is a ghost.
*/

/*!
	\enum Interface::Position

	This enum defines the element types that can be used.

	\var Interface::Type Interface::LEFT
	Left side of the interface.

	\var Interface::Type Interface::RIGHT
	Right side of the interface.
*/

/*!
	Default constructor.
*/
Interface::Interface()
	: Element(), m_normal(NULL)
{

}

/*!
	Creates a new interface.
*/
Interface::Interface(const long &id, ElementInfo::Type type)
	: Element(id, type), m_normal(NULL)
{

}


/*!
	Sets the area of the interface.

	\param area the area of the interface
*/
void Interface::set_area(double *area)
{
    m_area = area;
}

/*!
	Gets the area of the interface.

	\result The area of the interface
*/
const double & Interface::get_area() const
{
    return *m_area;
}

/*!
	Sets the normal of the interface.

	\param normal a pointer to the normal of the interface
*/
void Interface::set_normal(std::array<double, 3> *normal)
{
	m_normal = normal;
}

/*!
	Gets the normal of the interface.

	\result A pointer to the normal of the interface
*/
const std::array<double, 3> & Interface::get_normal() const
{
	return *m_normal;
}

/*!
	Evaluates the rotation matrix from the Cartesian coordinate system
	to the interface coordinate system.

	\result The rotation matrix from the Cartesian coordinate system
	to the interface coordinate system.
*/
std::array<std::array<double, 3>, 3> Interface::eval_rotation_from_cartesian()
{
	return eval_rotation_from_cartesian(*m_normal);
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
std::array<std::array<double, 3>, 3> Interface::eval_rotation_from_cartesian(std::array<double, 3> &versor)
{
	// The rotation matrix has in its rows the versors that define
	// the interface coordinate system.
	//
	//                | [x_int] |
	//          R =   | [y_int] |
	//                | [z_int] |
	//
	std::array<std::array<double, 3>, 3> R;
	for (int i = 0; i < 3; ++i) {
		R[i].fill(0.0);
	}

	// x-interface axis
	for (int k = 0; k < 3; ++k) {
		R[0][k] = versor[k];
	}

	// y-interface axis
	if (fabs(versor[2] - 1.) > 1e-8) {
		std::array<double, 3> z = {0.0, 0.0, 1.0};
		R[1] = utils::cross_3D(z, R[0]);
	} else {
		std::array<double, 3> x = {1.0, 0.0, 0.0};
		R[1] = utils::cross_3D(x, R[0]);
	}
	utils::normalize_3D(R[1]);

	// z-interface axis
	R[2] = utils::cross_3D(R[0], R[1]);
	utils::normalize_3D(R[2]);

	return R;
}

/*!
	Evaluates the rotation matrix from the interface coordinate system
	to the Cartesian coordinate system.

	Evaluates the rotation matrix that needs to be applied to the
	coordinates system defined on the interface to make it coincide
	with the Cartesian coordinates system.

	\result The rotation matrix from the interface coordinate system
	to the Cartesian coordinate system.
*/
std::array<std::array<double, 3>, 3> Interface::eval_rotation_to_cartesian()
{
	return eval_rotation_to_cartesian(*m_normal);
}

/*!
	Evaluates the rotation matrix from the coordinate system build
	starting from the specified versor to the Cartesian coordinate
	system.

	Evaluates the rotation matrix that needs to be applied to the
	coordinates system defined starting from the specified versor
	to make it coincide with the Cartesian coordinates system.
	This matrix can be evaluated as the transpose of the rotation
	matrix from the Cartesian coordinate system to the versor
	coordinate system.

	\param versor is the three-dimensional versor that defines the
	coordinate system
	\result The rotation matrix from the coordinate system build
	starting from the specified versor to the Cartesian coordinate
	system.
*/
std::array<std::array<double, 3>, 3> Interface::eval_rotation_to_cartesian(std::array<double, 3> &versor)
{
	std::array<std::array<double, 3>, 3> R = eval_rotation_from_cartesian(versor);
	utils::transpose_3D(R);

	return R;
}

/*!
	Evaluates the transpose of the specified rotation matrix.

	\param R the rotation matrix to transpose
	\result The transposed rotation matrix.
*/
std::array<std::array<double, 3>, 3> Interface::eval_rotation_transpose(const std::array<std::array<double, 3>, 3> &R)
{
	std::array<std::array<double, 3>, 3> R_prime = R;
	utils::transpose_3D(R_prime);

	return R_prime;
}

/*!
	Sets the position type of the interface.

	\param positionType the position type of the interface
*/
void Interface::set_position_type(PositionType positionType)
{
	m_positionType = positionType;
}

/*!
	Gets the position type of the interface.

	\result The position type of the nterface
*/
Interface::PositionType Interface::get_position_type() const
{
	return m_positionType;
}

/*!
	Sets the owner of the interface.

	\param owner the owner of the interface
*/
void Interface::set_owner(const long &owner, const int &onwerFace)
{
	m_owner     = owner;
	m_ownerFace = onwerFace;
}

/*!
	Deletes the owner of the interface.
*/
void Interface::unset_owner()
{
	m_owner     = Element::NULL_ELEMENT_ID;
	m_ownerFace = -1;
}

/*!
	Gets the owner of the interface.

	\result The owner of the nterface
*/
long Interface::get_owner() const
{
  return m_owner;
}

/*!
	Gets the face of the owner associated with the interface.

	\result The face of the owner associated with the interface
*/
int Interface::get_owner_face() const
{
  return m_ownerFace;
}

/*!
	Sets the neighbour of the interface.

	\param neigh the neighbour of the interface
*/
void Interface::set_neigh(const long &neigh, const int &onwerFace)
{
	m_neigh     = neigh;
	m_neighFace = onwerFace;
}

/*!
	Deletes the neighbour of the interface.
*/
void Interface::unset_neigh()
{
	m_neigh     = Element::NULL_ELEMENT_ID;
	m_neighFace = -1;
}

/*!
	Gets the neighbour of the interface.

	\result The neighbour of the nterface
*/
long Interface::get_neigh() const
{
  return m_neigh;
}

/*!
	Gets the face of the neighbour associated with the interface.

	\result The face of the neighbour associated with the interface
*/
int Interface::get_neigh_face() const
{
  return m_neighFace;
}

/*!
	Gets both the owner and the neighbour of the interface.

	\result An array containing the owner and the neighbour of the
	        interface.
*/
std::array<long, 2> Interface::get_owner_neigh() const
{
	std::array<long, 2> cells;
	cells[0] = m_owner;
	cells[1] = m_neigh;

	return cells;
}

}
