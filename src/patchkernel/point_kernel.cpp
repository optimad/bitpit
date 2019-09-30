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

#include "point_kernel.hpp"

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
	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(MPI_Comm communicator, bool expert)
	: PatchKernel(communicator, 0, expert)
#else
/*!
	Creates a patch.

	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(bool expert)
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
	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(int dimension, MPI_Comm communicator, bool expert)
	: PatchKernel(dimension, communicator, 0, expert)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(int dimension, bool expert)
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
	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(int id, int dimension, MPI_Comm communicator, bool expert)
	: PatchKernel(id, dimension, communicator, 0, expert)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PointKernel::PointKernel(int id, int dimension, bool expert)
	: PatchKernel(id, dimension, expert)
#endif
{
}

/*!
    Initialize the patch
*/
void PointKernel::initialize()
{
    // Nothing to do
}

/*!
 * Evaluates the characteristic size of the specified cell.
 *
 * \param id is the id of the cell
 * \result The characteristic size of the specified cell.
*/
double PointKernel::evalCellSize(long id) const
{
    BITPIT_UNUSED(id);

    return 0.;
}

/*!
 * Evaluates the distance between the specified points.
 *
 * \param id1 is the id of the first point
 * \param id2 is the id of the second point
 * \result The distance between the specified points.
*/
double PointKernel::evalPointsDistance(long id1, long id2) const
{
    const Cell &cell1 = m_cells[id1];
    const Cell &cell2 = m_cells[id2];

    ConstProxyVector<long> cellVertexIds1 = cell1.getVertexIds();
    const Vertex &vertex1 = getVertex(cellVertexIds1[0]);

    ConstProxyVector<long> cellVertexIds2 = cell2.getVertexIds();
    const Vertex &vertex2 = getVertex(cellVertexIds2[0]);

    double distance = norm2(vertex2.getCoords() - vertex1.getCoords());

    return distance;
}

/*!
 * Evaluates the direction of the line that connects the specified points.
 *
 * \param id1 is the id of the first point
 * \param id2 is the id of the second point
 * \result The direction of the line that connects the specified points.
*/
std::array<double, 3> PointKernel::evalPointsDirection(long id1, long id2) const
{
    const Cell &cell1 = m_cells[id1];
    const Cell &cell2 = m_cells[id2];

    ConstProxyVector<long> cellVertexIds1 = cell1.getVertexIds();
    const Vertex &vertex1 = getVertex(cellVertexIds1[0]);

    ConstProxyVector<long> cellVertexIds2 = cell2.getVertexIds();
    const Vertex &vertex2 = getVertex(cellVertexIds2[0]);

    std::array<double, 3> direction = vertex2.getCoords() - vertex1.getCoords();

    double distance = norm2(direction);
    if (!utils::DoubleFloatingEqual()(distance, 0., getTol())) {
        direction /= distance;
    }

    return direction;
}

}
