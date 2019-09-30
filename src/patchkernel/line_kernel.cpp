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

#include "line_kernel.hpp"

namespace bitpit {


#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(MPI_Comm communicator, std::size_t haloSize, bool expert)
	: PatchKernel(communicator, haloSize, expert)
#else
/*!
	Creates a patch.

	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(bool expert)
	: PatchKernel(expert)
#endif
{
	initialize();
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert)
	: PatchKernel(dimension, communicator, haloSize, expert)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(int dimension, bool expert)
	: PatchKernel(dimension, expert)
#endif
{
    initialize();
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert)
	: PatchKernel(id, dimension, communicator, haloSize, expert)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
LineKernel::LineKernel(int id, int dimension, bool expert)
	: PatchKernel(id, dimension, expert)
#endif
{
}

/*!
    Initialize the patch
*/
void LineKernel::initialize()
{
    // Nothing to do
}

/*!
 * Evaluates the characteristic size of the specified cell.
 *
 * \param id is the id of the cell
 * \result The characteristic size of the specified cell.
*/
double LineKernel::evalCellSize(long id) const
{
    return evalCellLength(id);
}

/*!
 * Evaluate the length of the specified cell.
 *
 * If cell is of type ElementType::VERTEX, the function returns 0.0
 *
 * \param id is the id of the cell
 * \result The length of the specified cell.
*/
double LineKernel::evalCellLength(long id) const
{
    const Cell &cell = m_cells[id];
    switch (cell.getType()) {

    case ElementType::LINE:
    {
        ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
        const Vertex &vertex_0 = getVertex(cellVertexIds[0]);
        const Vertex &vertex_1 = getVertex(cellVertexIds[1]);
        double length = norm2(vertex_1.getCoords() - vertex_0.getCoords());

        return length;
    }

    default:
    {
        return 0.;
    }

    }
}

/*!
 * Evaluate the normal of the specified cell.
 *
 * If cell is of type ElementType::VERTEX or ElementType::LINE, returns 0.0
 *
 * \param id is the id of the cell
 * \param orientation is a vector carring the additional information needed
 * to un-ambigously define a normal to the element (e.g., when evaluating
 * the normal of a one-dimensional element, this versor is perpendicular to
 * the plane where the normal should lie)
 * \result The normal of the specified cell.
*/
std::array<double, 3> LineKernel::evalCellNormal(long id, const std::array<double, 3> &orientation) const
{
    const Cell &cell = m_cells[id];
    switch (cell.getType()) {

    case ElementType::LINE:
    {
        ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
        const Vertex &vertex_0 = getVertex(cellVertexIds[0]);
        const Vertex &vertex_1 = getVertex(cellVertexIds[1]);

        std::array<double, 3> normal = vertex_1.getCoords() - vertex_0.getCoords();
        normal = crossProduct(normal, orientation);
        normal = normal / norm2(normal);

        return normal;
    }

    default:
    {
        return {{0., 0., 0.}};
    }

    }
}

}
