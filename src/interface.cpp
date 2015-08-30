//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

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

#include "cell.hpp"
#include "interface.hpp"

namespace pman {

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
Interface::Interface(const int &id)
	: Element(id), m_normal(NULL)
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
double Interface::get_area() const
{
    return *m_area;
}

/*!
	Sets the normal of the interface.

	\param normal a pointer to the normal of the interface
*/
void Interface::set_normal(double normal[])
{
	m_normal = normal;
}

/*!
	Gets the normal of the interface.

	\result A pointer to the normal of the interface
*/
double * Interface::get_normal() const
{
	return m_normal;
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
void Interface::set_owner(Cell *owner, const int &onwerFace)
{
	m_owner = owner;
	m_ownerFace = onwerFace;
}

/*!
	Deletes the owner of the interface.
*/
void Interface::unset_owner()
{
	m_owner = NULL;
	m_ownerFace = -1;
}

/*!
	Gets the owner of the interface.

	\result The owner of the nterface
*/
Cell * Interface::get_owner() const
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
void Interface::set_neigh(Cell * neigh, const int &onwerFace)
{
	m_neigh = neigh;
	m_neighFace = onwerFace;
}

/*!
	Deletes the neighbour of the interface.
*/
void Interface::unset_neigh()
{
	m_neigh = NULL;
	m_neighFace = -1;
}

/*!
	Gets the neighbour of the interface.

	\result The neighbour of the nterface
*/
Cell * Interface::get_neigh() const
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
	Sets the data of the interface.

	\param data a pointer to the data of the interface
*/
void Interface::set_data(std::unique_ptr<InterfaceData> data)
{
	m_data = std::move(data);
}

/*!
	Gets the data of the interface.

	\return A pointer to the data of the interface
*/
InterfaceData * Interface::get_data() const
{
    return m_data.get();
}

}
