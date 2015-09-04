//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\class Cell

	\brief The Cell class defines the cells.

	Cell is class that defines the cells.
*/

#include "cell.hpp"
#include "interface.hpp"

#include<iostream>

namespace pman {

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
Cell::Cell(const int &id)
	: Element(id)
{

}

/*!
	Creates a new cell.
*/
Cell::Cell(const int &id, Patch *patch)
	: Element(id, patch)
{

}

/*!
	Sets the position type of the interface.

	\param position the position type of the interface
*/
void Cell::set_position_type(PositionType positionType)
{
	m_positionType = positionType;
}

/*!
	Gets the position type of the interface.

	\result The position type of the nterface
*/
Cell::PositionType Cell::get_position_type() const
{
  return m_positionType;
}

/*!
	Sets the volume of the cell.

	\param volume the volume of the cell
*/
void Cell::set_volume(double * const volume)
{
	m_volume = volume;
}

/*!
	Gets the volume of the cell.

	\return The volume of the cell
*/
double & Cell::get_volume() const
{
    return *m_volume;
}  

/*!
	Sets the centroid of the cell.

	\param centroid the centroid of the cell
*/
void Cell::set_centroid(double * const centroid)
{
	m_centroid = centroid;
}

/*!
	Gets the centroid of the cell.

	\return The centroid of the cell
*/
double * Cell::get_centroid() const
{
	return m_centroid;
}

/*!
	Initialize all the interfaces of the cell.

	\param interfaces the list of all interfaces associated to the cell
*/
void Cell::initialize_interfaces(std::vector<std::vector<int>> &interfaces)
{
	m_interfaces.reset();
	m_interfaces = std::unique_ptr<CollapsedArrayArray<int> >(new CollapsedArrayArray<int>(interfaces));
}

/*!
	Initialize the data structure that holds the information about the
	interfaces.

	\param nInterfaces An array with the numbero of interfaces for each face
*/
void Cell::initialize_empty_interfaces(const int nInterfaces[])
{
	m_interfaces = std::unique_ptr<CollapsedArrayArray<int> >(new CollapsedArrayArray<int>(get_face_count(), nInterfaces));
}

/*!
	Sets the i-th interface associated the the given face of the cell.

	\param face the face of the cell
	\param index the index of the interface
	\param interface A pointer to the interface
*/
void Cell::set_interface(const int &face, const int &index, const int &interface)
{
	m_interfaces->set(face, index, interface);
}

/*!
	Sets the interfaces associated to the given face of the cell.

	\param face the face of the cell
	\param interface A pointer to the interface
*/
void Cell::set_interfaces(const int &face, int interfaces[])
{
	m_interfaces->set(face, interfaces);
}

/*!
	Unsets the interfaces associated to the cell.

	\param face the face of the cell
	\param interface A pointer to the interface
*/
void Cell::unset_interfaces()
{
	m_interfaces.reset();
}

/*!
	Gets the number of interfaces of the element for the given face.

	\result The number of interfaces on the face  elements
*/
int Cell::get_interface_count(const int &face) const
{
	return m_interfaces->sub_array_size(face);
}

/*!
	Gets the i-th interface of the given face of the cell.

	\as get_interfaces(const int &face)

	\param face the face of the cell
	\param index the index of the interface to retreive
	\result The requested interface
*/
int Cell::get_interface(const int &face, const int &index) const
{
	return m_interfaces->get(face, index);
}

/*!
	Gets the interfaces of the given face of the cell.

	\as get_interface(const int &face, const int &index) const

	\param face the face of the cell
	\result The requested interfaces
*/
int * Cell::get_interfaces(const int &face) const
{
	return m_interfaces->get(face);
}

/*!
	Sets the data of the cell.

	\param data a pointer to the data of the cell
*/
void Cell::set_data(std::unique_ptr<CellData> data)
{
	m_data = std::move(data);
}

/*!
	Gets the data of the cell.

	\return A pointer to the data of the cell
*/
CellData * Cell::get_data() const
{
    return m_data.get();
}

}
