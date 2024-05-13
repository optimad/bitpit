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

#include "bitpit_CG.hpp"
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
	\param adaptionMode is the adaption mode that will be used for the patch
	\param partitioningMode is the partitioning mode that will be used for the
	patch
*/
LineKernel::LineKernel(MPI_Comm communicator, std::size_t haloSize,
                       AdaptionMode adaptionMode, PartitioningMode partitioningMode)
	: PatchKernel(communicator, haloSize, adaptionMode, partitioningMode)
#else
/*!
	Creates a patch.

	\param adaptionMode is the adaption mode that will be used for the patch
*/
LineKernel::LineKernel(AdaptionMode adaptionMode)
	: PatchKernel(adaptionMode)
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
	\param adaptionMode is the adaption mode that will be used for the patch
	\param partitioningMode is the partitioning mode that will be used for the
	patch
*/
LineKernel::LineKernel(int dimension, MPI_Comm communicator, std::size_t haloSize,
                       AdaptionMode adaptionMode, PartitioningMode partitioningMode)
	: PatchKernel(dimension, communicator, haloSize, adaptionMode, partitioningMode)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
LineKernel::LineKernel(int dimension, AdaptionMode adaptionMode)
	: PatchKernel(dimension, adaptionMode)
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
	\param adaptionMode is the adaption mode that will be used for the patch
	\param partitioningMode is the partitioning mode that will be used for the
	patch
*/
LineKernel::LineKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize,
                       AdaptionMode adaptionMode, PartitioningMode partitioningMode)
	: PatchKernel(id, dimension, communicator, haloSize, adaptionMode, partitioningMode)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
LineKernel::LineKernel(int id, int dimension, AdaptionMode adaptionMode)
	: PatchKernel(id, dimension, adaptionMode)
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
    Get the codimension of the patch in the volume space.

    \result The codimension of the patch in the volume space.
*/
int LineKernel::getVolumeCodimension() const
{
    return 2;
}

/*!
    Get the codimension of the patch in the surface space.

    \result The codimension of the patch in the surface space.
*/
int LineKernel::getSurfaceCodimension() const
{
    return 1;
}

/*!
    Get the codimension of the patch in the line space.

    \result The codimension of the patch in the line space.
*/
int LineKernel::getLineCodimension() const
{
    return 0;
}

/*!
    Get the codimension of the patch in the point space.

    \result The codimension of the patch in the point space.
*/
int LineKernel::getPointCodimension() const
{
    return -1;
}

/*!
    Extracts the external envelope and appends it to the given patch.

    The external envelope is composed by all the free faces of the patch.

    \param[in,out] envelope is the patch to which the external envelope
    will be appended
*/
void LineKernel::extractEnvelope(PointKernel &envelope) const
{
    PatchKernel::extractEnvelope(envelope);
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

/*!
 *
 * Evaluates the baricentric coordinates of the specified cell.
 *
 * If cell is not of type ElementType::LINE, the function returns 0.0
 *
 * \param[in] id is the id of the cell
 * \param[in] point are the coordinates of point
 * \param[out] lambda on output will contain the barycentric coordinates of the projection point
 */
void LineKernel::evalBarycentricCoordinates(long id, const std::array<double, 3> &point, double *lambda) const
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                   *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE BARYCENTRIC COORDINATES
    // ====================================================================== //
    switch (cell_->getType()) {

    case ElementType::LINE:
    {
        ConstProxyVector<long> vertexIds = cell_->getVertexIds();

        std::array<double,3> point0 = getVertexCoords(vertexIds[0]);
        std::array<double,3> point1 = getVertexCoords(vertexIds[1]);

        std::array<double, 3> projectionPoint = CGElem::projectPointSegment(point, point0, point1, lambda);
        BITPIT_UNUSED(projectionPoint);
        return;
    }

    default:
    {
        int nVertices = cell_->getVertexCount();
        for (int i = 0; i < nVertices; ++i) {
            lambda[i] = 0.0;
        }
        return;
    }

    }
}

}
