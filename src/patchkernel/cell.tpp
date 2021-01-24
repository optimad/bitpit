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

#ifndef __BITPIT_CELL_TPP__
#define __BITPIT_CELL_TPP__

namespace bitpit {

/*!
	\class CellHalfEdge
	\ingroup patchelements

	\brief The CellHalfEdge class defines cell half-edges.

	CellHalfEdge is the class that defines cell half-edges.
*/

/*!
	Constructor.

	\param cell is a reference to the cell the owns the edge
	\param edge if the local edge of the cell
	\param winding is the winding order of the vertices
*/
template<class QualifiedCell>
QualifiedCellHalfEdge<QualifiedCell>::QualifiedCellHalfEdge(QualifiedCell &cell, int edge, Winding winding)
    : ElementHalfEdge<QualifiedCell>(cell, edge, winding)
{
}

/*!
	Get the element the half-item belongs to.

	\result Returns the element the half-item belongs to.
*/
template<class QualifiedCell>
QualifiedCell & QualifiedCellHalfEdge<QualifiedCell>::getCell() const
{
    return static_cast<QualifiedCell &>(this->getElement());
}

/*!
	\class QualifiedCellHalfFace
	\ingroup patchelements

	\brief The QualifiedCellHalfFace class defines cell half-faces.

	QualifiedCellHalfFace is the class that defines cell half-faces. Each
	face can be seen as two half-faces: one belonging to a cell and the other
	belonging to the neighbouring cell. A half-face is identify by its vertices
	and by the winding order of the vertices.
*/

/*!
	Constructor.

	\param cell is a reference to the cell the owns the face
	\param face if the local face of the cell
	\param winding is the winding order of the vertices
*/
template<class QualifiedCell>
QualifiedCellHalfFace<QualifiedCell>::QualifiedCellHalfFace(QualifiedCell &cell, int face, Winding winding)
    : ElementHalfFace<QualifiedCell>(cell, face, winding)
{
}

/*!
	Get the element the half-item belongs to.

	\result Returns the element the half-item belongs to.
*/
template<class QualifiedCell>
QualifiedCell & QualifiedCellHalfFace<QualifiedCell>::getCell() const
{
    return static_cast<QualifiedCell &>(this->getElement());
}

}

#endif
