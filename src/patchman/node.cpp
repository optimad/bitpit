//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "node.hpp"

#include <limits>

namespace pman {

/*!
    \ingroup    PatchMan
    @{
	\class Node

	\brief The Node class defines the nodes.

	Node is class that defines the nodes.
*/

const long Node::NULL_NODE_ID = std::numeric_limits<long>::min();

/*!
	Default constructor.
*/
Node::Node()
{
	set_id(NULL_NODE_ID);
}

/*!
	Creates a new element.
*/
Node::Node(const long &id)
{
	set_id(id);
}

/*!
	Sets the ID of the node.

	\param id the ID of the node
*/
void Node::set_id(const long &id)
{
	m_id = id;
}

/*!
	Gets the ID of the node.

	\return The ID of the node
*/
long Node::get_id() const
{
	return m_id;
}

/*!
	Sets the coordinates of the node.

	\param coords are the coordinates of the node
*/
void Node::set_coords(std::array<double, 3> &coords)
{
	m_coords = coords;
}

/*!
	Gets the coordinates of the node.

	\return A pointer to the coordinates of the node
*/
const std::array<double, 3> & Node::get_coords() const
{
	return m_coords;
}

}

/*! @} */
