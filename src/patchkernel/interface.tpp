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

#ifndef __BITPIT_INTERFACE_TPP__
#define __BITPIT_INTERFACE_TPP__

namespace bitpit {

/*!
	\class QualifiedInterfaceHalfEdge
	\ingroup patchelements

	\brief The QualifiedInterfaceHalfEdge class defines interface half-edges.

	QualifiedInterfaceHalfEdge is the class that defines interface half-edges.
	Each edge can be seen as two half-edges: one belonging to an interface
	and the other belonging to the neighbouring interface. A half-edge is
	identify by its vertices and by the winding order of the vertices.
*/

/*!
	Constructor.

	\param interface is a reference to the interface the owns the edge
	\param edge if the local edge of the interface
	\param winding is the winding order of the vertices
*/
template<class QualifiedInterface>
QualifiedInterfaceHalfEdge<QualifiedInterface>::QualifiedInterfaceHalfEdge(QualifiedInterface &interface, int edge, Winding winding)
    : ElementHalfFace<QualifiedInterface>(interface, edge, winding)
{
}

/*!
	Get the element the half-item belongs to.

	\result Returns the element the half-item belongs to.
*/
template<class QualifiedInterface>
QualifiedInterface & QualifiedInterfaceHalfEdge<QualifiedInterface>::getInterface() const
{
    return static_cast<QualifiedInterface &>(this->getElement());
}

}

#endif
