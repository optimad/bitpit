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

#define _USE_MATH_DEFINES

#include <cmath>
#include <limits>
#include <set>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#if BITPIT_ENABLE_MPI==1
#include "bitpit_communications.hpp"
#endif

#include "surface_kernel.hpp"

namespace bitpit {

/*!
	\class SurfaceKernel
	\ingroup surfacepatches

	\brief The SurfaceKernel class provides an interface for defining
	surface patches.

	SurfaceKernel is the base class for defining surface patches.
*/

const std::map<ElementType, unsigned short>  SurfaceKernel::m_selectionTypes({
                                {ElementType::QUAD,     SurfaceKernel::SELECT_QUAD},
                                {ElementType::TRIANGLE, SurfaceKernel::SELECT_TRIANGLE}
                                });
const unsigned short SurfaceKernel::SELECT_TRIANGLE = 1;
const unsigned short SurfaceKernel::SELECT_QUAD     = 2;
const unsigned short SurfaceKernel::SELECT_ALL      = 3;

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
*/
SurfaceKernel::SurfaceKernel(MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param adaptionMode is the adaption mode that will be used for the patch
*/
SurfaceKernel::SurfaceKernel(AdaptionMode adaptionMode)
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
*/
SurfaceKernel::SurfaceKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(dimension, communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
SurfaceKernel::SurfaceKernel(int dimension, AdaptionMode adaptionMode)
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
*/
SurfaceKernel::SurfaceKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(id, dimension, communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
SurfaceKernel::SurfaceKernel(int id, int dimension, AdaptionMode adaptionMode)
	: PatchKernel(id, dimension, adaptionMode)
#endif
{
}

/*!
	Initialize the patch
*/
void SurfaceKernel::initialize()
{
    // Nothing to do
}

/*!
	Get the codimension of the patch in the volume space.

	\result The codimension of the patch in the volume space.
*/
int SurfaceKernel::getVolumeCodimension() const
{
	return 1;
}

/*!
	Get the codimension of the patch in the surface space.

	\result The codimension of the patch in the surface space.
*/
int SurfaceKernel::getSurfaceCodimension() const
{
	return 0;
}

/*!
	Get the codimension of the patch in the line space.

	\result The codimension of the patch in the line space.
*/
int SurfaceKernel::getLineCodimension() const
{
	return -1;
}

/*!
	Get the codimension of the patch in the point space.

	\result The codimension of the patch in the point space.
*/
int SurfaceKernel::getPointCodimension() const
{
	return -2;
}

/*!
	Extracts the external envelope and appends it to the given patch.

	The external envelope is composed by all the free faces of the patch.

	\param[in,out] envelope is the patch to which the external envelope
	will be appended
*/
void SurfaceKernel::extractEnvelope(LineKernel &envelope) const
{
	PatchKernel::extractEnvelope(envelope);
}

/*!
        Evaluates the characteristic size of the specified cell.

        \param id is the id of the cell
        \result The characteristic size of the specified cell.
*/
double SurfaceKernel::evalCellSize(long id) const
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                   *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE CELL SIZE                                                      //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)) return 0.0;
    if (cell_->getType() == ElementType::LINE) {
        return evalCellArea(id);
    }
    return(sqrt(evalCellArea(id))); 
}

/*!
 * Evaluate facet area for a cell with specified ID. If cell is of type
 * ElementType::VERTEX or ElementType::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet area
*/
double SurfaceKernel::evalCellArea(long id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                   *cell_ = &m_cells[id];
    ConstProxyVector<long>       cellVertexIds = cell_->getVertexIds();

    // Counters
    // none

    // ====================================================================== //
    // EVALUATE FACET AREA                                                    //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)) return 0.0;
    if (cell_->getType() == ElementType::LINE) {
        return ( norm2(m_vertices[cellVertexIds[0]].getCoords()
                     - m_vertices[cellVertexIds[1]].getCoords()) );
    }
    std::array<double, 3>       d1, d2;
    if (cell_->getType() == ElementType::TRIANGLE) {
        d1 = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        d2 = m_vertices[cellVertexIds[2]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        return (0.5*norm2(crossProduct(d1, d2)));
    }
    else {
        int                     nvert = cellVertexIds.size();
        double                  coeff = 0.25;
        double                  area = 0.0;
        for (int i = 0; i < nvert; ++i) {
            int prevVertex = getFacetOrderedLocalVertex(*cell_, (nvert + i - 1) % nvert);
            int vertex     = getFacetOrderedLocalVertex(*cell_, i);
            int nextVertex = getFacetOrderedLocalVertex(*cell_, (i + 1) % nvert);

            const std::array<double, 3> &prevVertexCoords = m_vertices[cellVertexIds[prevVertex]].getCoords();
            const std::array<double, 3> &vertexCoords     = m_vertices[cellVertexIds[vertex]].getCoords();
            const std::array<double, 3> &nextVertexCoords = m_vertices[cellVertexIds[nextVertex]].getCoords();

            d1 = nextVertexCoords - vertexCoords;
            d2 = prevVertexCoords - vertexCoords;
            area += coeff*norm2(crossProduct(d1, d2));
        } //next i
        return(area);
    }
}

/*!
 *  Evaluate the length of the edge with specified local index for the cell
 *  with specified Iid.
 *
 *  If the cell is of type ElementType::VERTEX or ElementType::LINE returns 0.
 * 
 *  \param[in] cellId is the cell id
 *  \param[in] edgeId is the edge local index
 *  \result The edge length.
*/
double SurfaceKernel::evalEdgeLength(long cellId, int edgeId) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                           &cell = m_cells[cellId];
    double                               edge_length;

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    switch (cell.getType()) {

    case ElementType::LINE:
    case ElementType::VERTEX:
    case ElementType::UNDEFINED:
    {
        edge_length = 0.;
        break;
    }

    default:
    {
        ConstProxyVector<long> faceVertexIds = cell.getFaceVertexIds(edgeId);
        const Vertex &vertex_0 = m_vertices[faceVertexIds[0]];
        const Vertex &vertex_1 = m_vertices[faceVertexIds[1]];
        edge_length = norm2(vertex_0.getCoords() - vertex_1.getCoords());
        break;
    }

    }

    return edge_length;
}

/*!
 *  Evaluate the minimal edge length for e cell with specified ID.
 *  If the cell is of type ElementType::VERTEX or ElementType::LINE
 *  returns 0.0.
 * 
 *  \param[in] id cell id
 *  \param[in,out] edge_id on output stores the local index of edge with minimal edge length
 * 
 *  \result minimal edge length
*/
double SurfaceKernel::evalMinEdgeLength(long id, int &edge_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                          *cell_ = &m_cells[id];

    // Counters
    int                                 i;
    int                                 n_faces;

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::LINE)
     || (cell_->getType() == ElementType::VERTEX)
     || (cell_->getType() == ElementType::UNDEFINED)) return 0.0;

    double edge_length = std::numeric_limits<double>::max(), tmp;
    n_faces = cell_->getFaceCount();
    for (i = 0; i < n_faces; ++i) {
        tmp = evalEdgeLength(id, i);
        if (tmp < edge_length) {
            edge_length = tmp;
            edge_id = i;
        }
    } //next i

    return(edge_length);
}

/*!
 *  Evaluate the maximal edge length for e cell with specified ID.
 *  If the cell is of type ElementType::VERTEX or ElementType::LINE
 *  returns 0.0.
 * 
 *  \param[in] id cell id
 *  \param[in,out] edge_id on output stores the local inde xof edge with maximal edge length
 * 
 *  \result maximal edge length
*/
double SurfaceKernel::evalMaxEdgeLength(long id, int &edge_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                          *cell_ = &m_cells[id];

    // Counters
    int                                 i;
    int                                 n_faces;

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::LINE)
     || (cell_->getType() == ElementType::VERTEX)
     || (cell_->getType() == ElementType::UNDEFINED)) return 0.0;

    double edge_length = std::numeric_limits<double>::min(), tmp;
    n_faces = cell_->getFaceCount();
    for (i = 0; i < n_faces; ++i) {
        tmp = evalEdgeLength(id, i);
        if (tmp > edge_length) {
            edge_length = tmp;
            edge_id = i;
        }
    } //next i

    return(edge_length);
}

/*!
 * Evaluate the angle at specified vertex for a cell with specified id.
 *
 * The function will return zero if the cell is of type ElementType::UNDEFINED
 * or ElementType::VERTEX. It will also return zero if the cell is of type
 * ElementType::LINE and the codimension of the surface with respect to the
 * space in which it is embedded is not one.
 * 
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
*/
double SurfaceKernel::evalAngleAtVertex(long id, int vertex) const
{
    const Cell &cell = m_cells.at(id);
    ElementType cellType = cell.getType();
    if (cellType == ElementType::UNDEFINED) {
        return 0.;
    } else if (cellType == ElementType::VERTEX) {
        return 0.;
    } else if (cellType == ElementType::LINE) {
        return BITPIT_PI;
    }

    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();

    int prevVertex = getFacetOrderedLocalVertex(cell, (vertex - 1 + nCellVertices) % nCellVertices);
    int nextVertex = getFacetOrderedLocalVertex(cell, (vertex + 1) % nCellVertices);

    const std::array<double, 3> &prevVertexCoords = m_vertices[cellVertexIds[prevVertex]].getCoords();
    const std::array<double, 3> &vertexCoords     = m_vertices[cellVertexIds[vertex]].getCoords();
    const std::array<double, 3> &nextVertexCoords = m_vertices[cellVertexIds[nextVertex]].getCoords();

    std::array<double, 3> d1 = prevVertexCoords - vertexCoords;
    d1 = d1 / norm2(d1);

    std::array<double, 3> d2 = nextVertexCoords - vertexCoords;
    d2 = d2 / norm2(d2);

    double angle = std::acos(std::min(1., std::max(-1., dotProduct(d1, d2))));

    return angle;
}

/*!
 * Evaluate the minimal angle at vertex for a cell with specified ID.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, a returns zero.
 *
 * \param[in] id cell id
 * \param[in,out] vertex_id on output stores the local index of vertex with minimal angle
 * 
 * \result minimal angle at vertex
*/
double SurfaceKernel::evalMinAngleAtVertex(long id, int &vertex_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell          *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE MINIMAL ANGLE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)
     || (cell_->getType() == ElementType::LINE)) return 0.0;

    double angle = std::numeric_limits<double>::max(), tmp;
    int n_vert = cell_->getVertexCount();
    for (int i = 0; i < n_vert; ++i) {
        tmp = evalAngleAtVertex(id, i);
        if (tmp < angle) {
            angle = tmp;
            vertex_id = i;
        }
    } //next i

    return(angle);
}

/*!
 * Evaluate the maximal angle at vertex for a cell with specified ID.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, a returns zero.
 *
 * \param[in] id cell id
 * \param[in,out] vertex_id on output stores the local index of vertex with minimal angle
 * 
 * \result maximal angle at vertex
*/
double SurfaceKernel::evalMaxAngleAtVertex(long id, int &vertex_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell          *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE MINIMAL ANGLE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)
     || (cell_->getType() == ElementType::LINE)) return 0.0;

    double angle = std::numeric_limits<double>::min(), tmp;
    int n_vert = cell_->getVertexCount();
    for (int i = 0; i < n_vert; ++i) {
        tmp = evalAngleAtVertex(id, i);
        if (tmp > angle) {
            angle = tmp;
            vertex_id = i;
        }
    } //next i

    return(angle);
}

/*!
 * Evalute the aspect ratio for a cell with specified ID. The aspect ratio
 * is defined as the ratio between the longest and the shortest edge.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * \param[in,out] edge_id on output stores the index of the shortest edge
 * 
 * \result cell aspect ratio
*/
double SurfaceKernel::evalAspectRatio(long id, int &edge_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                   *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // EVALUATE ASPECT RATIO                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)
     || (cell_->getType() == ElementType::LINE)) return 0.0;

    double                      l_edge;
    double                      m_edge = std::numeric_limits<double>::max();
    double                      M_edge = std::numeric_limits<double>::min();
    int                         nfaces = cell_->getFaceCount();
    for (int i = 0; i < nfaces; ++i) {
        l_edge = evalEdgeLength(id, i);
        if (l_edge < m_edge) {
            m_edge = l_edge;
            edge_id = i;
        }
        M_edge = std::max(M_edge, l_edge);
    } //next i

    return (M_edge/m_edge);

}

/*!
 * Evaluate facet normal for a cell with specified ID.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * \param orientation is a vector carring the additional information needed
 * to un-ambigously define a normal to the element (e.g., when evaluating
 * the normal of a one-dimensional element, this versor is perpendicular to
 * the plane where the normal should lie)
 *
 * \result facet normal
*/
std::array<double, 3> SurfaceKernel::evalFacetNormal(long id, const std::array<double, 3> &orientation) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    std::array<double, 3>        normal = {{0.0, 0.0, 0.0}};
    const Cell                   *cell_ = &m_cells[id];
    ConstProxyVector<long>       cellVertexIds = cell_->getVertexIds();

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE NORMAL                                                         //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)) return normal;
    
    if (cell_->getType() == ElementType::LINE) {
        normal = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        normal = crossProduct(normal, orientation);
    }

    if (cell_->getType() == ElementType::TRIANGLE) {
        std::array<double, 3>           d1, d2;
        d1 = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        d2 = m_vertices[cellVertexIds[2]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        normal = crossProduct(d1, d2);
    }
    else {
        std::array<double, 3>           d1, d2;
        int                             i, nvert = cellVertexIds.size();
        double                          coeff = 1.0/double(nvert);
        for (i = 0; i < nvert; ++i) {
            int prevVertex = getFacetOrderedLocalVertex(*cell_, (nvert + i - 1) % nvert);
            int vertex     = getFacetOrderedLocalVertex(*cell_, i);
            int nextVertex = getFacetOrderedLocalVertex(*cell_, (i + 1) % nvert);

            const std::array<double, 3> &prevVertexCoords = m_vertices[cellVertexIds[prevVertex]].getCoords();
            const std::array<double, 3> &vertexCoords     = m_vertices[cellVertexIds[vertex]].getCoords();
            const std::array<double, 3> &nextVertexCoords = m_vertices[cellVertexIds[nextVertex]].getCoords();

            d1 = nextVertexCoords - vertexCoords;
            d2 = prevVertexCoords - vertexCoords;
            normal += coeff*crossProduct(d1, d2);
        } //next i
    }
    normal = normal/norm2(normal);

    return(normal);

}

/*!
 * Evaluate the normal unit vector to the edge with specified local id belonging to
 * the surface facet with specified global id. Edge normal is computed as the arithmetic
 * average of the normals to each facet incident to the edge. In case adjacencies
 * are not built, the edge normal will be the same as the facet normal.
 * 
 * \param[in] id cell global ID
 * \param[in] edge_id edge local ID on the specified cell
*/
std::array<double, 3> SurfaceKernel::evalEdgeNormal(long id, int edge_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    std::array<double, 3>               normal = evalFacetNormal(id);
    const Cell                          *cell_ = &m_cells[id];
    int                                 n_adj = cell_->getAdjacencyCount(edge_id);
    const long                          *adjacencies = cell_->getAdjacencies(edge_id);

    // ====================================================================== //
    // COMPUTE EDGE NORMAL                                                    //
    // ====================================================================== //
    if (n_adj > 0) {
        for (int i = 0; i < n_adj; ++i) {
            normal += evalFacetNormal(adjacencies[i]);
        } //next i
        normal = normal/norm2(normal);
    }
    return(normal);
}

/*!
 * Evaluate the normal unit vector of the specified local vertex.
 *
 * Vertex normals are evaluated as a weighted average of the normals of the
 * cells that define the one-ring of the vertex. For each cell, the weight
 * is the angle beteen the two side that share the specfied vertex.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \result The normal unit vector of the specified local vertex.
*/
std::array<double, 3> SurfaceKernel::evalVertexNormal(long id, int vertex) const
{
    std::array<double, 3> normal;
    std::vector<long> vertexNeighs = findCellVertexNeighs(id, vertex);
    evalVertexNormals(id, vertex, vertexNeighs.size(), vertexNeighs.data(), std::numeric_limits<double>::max(), &normal, nullptr);

    return normal;
}

/*!
 * Evaluate the normal unit vector of the specified local vertex.
 *
 * Vertex normals are evaluated as a weighted average of the normals of the
 * cells that define the one-ring of the vertex. For each cell, the weight
 * is the angle beteen the two side that share the specfied vertex.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \param[in] nVertexNeighs is the number of vertex neighbours
 * \param[in] vertexNeighs are the neighbours of the vertex
 * \result The normal unit vector of the specified local vertex.
*/
std::array<double, 3> SurfaceKernel::evalVertexNormal(long id, int vertex, std::size_t nVertexNeighs, const long *vertexNeighs) const
{
    std::array<double, 3> normal;
    evalVertexNormals(id, vertex, nVertexNeighs, vertexNeighs, std::numeric_limits<double>::max(), &normal, nullptr);

    return normal;
}

/*!
 * Evaluate the limited normal unit vector of the specified local vertex.
 *
 * Vertex normals are evaluated as a weighted average of the normals of the
 * cells that define the one-ring of the vertex. For each cell, the weight
 * is the angle beteen the two side that share the specfied vertex. Only the
 * normals whose angle with respect to the considered cell is less that the
 * specified limit are considered.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \param[in] limit is the maximum allowed misalignment between the normal
 * of the reference cell and the normal of facets used for evaualting the
 * vertex normal
 * \result The normal unit vector of the specified local vertex.
 */
std::array<double, 3> SurfaceKernel::evalLimitedVertexNormal(long id, int vertex, double limit) const
{
    std::array<double, 3> normal;
    std::vector<long> vertexNeighs = findCellVertexNeighs(id, vertex);
    evalVertexNormals(id, vertex, vertexNeighs.size(), vertexNeighs.data(), limit, nullptr, &normal);

    return normal;
}

/*!
 * Evaluate the limited normal unit vector of the specified local vertex.
 *
 * Vertex normals are evaluated as a weighted average of the normals of the
 * cells that define the one-ring of the vertex. For each cell, the weight
 * is the angle beteen the two side that share the specfied vertex. Only the
 * normals whose angle with respect to the considered cell is less that the
 * specified limit are considered.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \param[in] nVertexNeighs is the number of vertex neighbours
 * \param[in] vertexNeighs are the neighbours of the vertex
 * \param[in] limit is the maximum allowed misalignment between the normal
 * of the reference cell and the normal of facets used for evaualting the
 * vertex normal
 * \result The normal unit vector of the specified local vertex.
 */
std::array<double, 3> SurfaceKernel::evalLimitedVertexNormal(long id, int vertex, std::size_t nVertexNeighs, const long *vertexNeighs, double limit) const
{
    std::array<double, 3> normal;
    evalVertexNormals(id, vertex, nVertexNeighs, vertexNeighs, limit, nullptr, &normal);

    return normal;
}

/*!
 * Evaluate unlimited and limited normal unit vectors at the specified local
 * vertex.
 *
 * Vertex normals are evaluated as a weighted average of the normals of the
 * cells that define the one-ring of the vertex. For each cell, the weight
 * is the angle beteen the two side that share the specfied vertex. When
 * evaluating the limited normal, only the normals whose angle with respect
 * to the considered cell is less that the specified limit are considered.
 * Whereas, when evaluating the unlimited normal, all cells in the one-ring
 * are considered.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \param[in] nVertexNeighs is the number of vertex neighbours
 * \param[in] vertexNeighs are the neighbours of the vertex
 * \param[in] limit is the maximum allowed misalignment between the normal
 * of the reference cell and the normal of facets used for evaualting the
 * vertex normal
 * \param[out] unlimitedNormal if a valid pointer is provided, on output will
 * contain the unlimited normal
 * \param[out] limitedNormal if a valid pointer is provided, on output will
 * contain the limited normal
*/
void SurfaceKernel::evalVertexNormals(long id, int vertex, std::size_t nVertexNeighs, const long *vertexNeighs, double limit,
                                      std::array<double, 3> *unlimitedNormal, std::array<double, 3> *limitedNormal) const
{
    // Early return if no calculation is needed
    if (!unlimitedNormal && !limitedNormal) {
        return;
    }

    // Get vertex id
    const Cell &cell = getCell(id);
    long vertexId = cell.getVertexId(vertex);

    // Get cell information
    double cellVertexAngle = evalAngleAtVertex(id, vertex);
    std::array<double, 3> cellNormal = evalFacetNormal(id);

    // Initialize unlimited normal
    if (unlimitedNormal) {
        *unlimitedNormal = cellVertexAngle * cellNormal;
    }

    // Initialize limited normal
    if (limitedNormal) {
        *limitedNormal = cellVertexAngle * cellNormal;
    }

    // Add contribution of neighbouring cells
    for (std::size_t i = 0; i < nVertexNeighs; ++i) {
        // Get neighbour information
        long neighId = vertexNeighs[i];
        const Cell &neigh = getCell(neighId);
        std::array<double, 3> neighNormal = evalFacetNormal(neighId);
        double neighVertexAngle = evalAngleAtVertex(neighId, neigh.findVertex(vertexId));

        // Add contribution to unlimited normal
        if (unlimitedNormal) {
            *unlimitedNormal += neighVertexAngle * neighNormal;
        }

        // Add contribution to limited normal
        //
        // Only the negihbours whose normal has a misalignment less then the
        // specified limit are considered.
        if (limitedNormal) {
            // The argument of the acos function has to be in the range [-1, 1].
            // Rounding errors may lead to a dot product slightly outside this
            // range. Since the arguments of the dot product are unit vectors,
            // we can safetly clamp the dot product result to be between -1 and
            // 1.
            double misalignment = std::acos(std::min(1.0, std::max(-1.0, dotProduct(neighNormal, cellNormal))));
            if (misalignment <= std::abs(limit)) {
                *limitedNormal += neighVertexAngle * neighNormal;
            }
        }
    }

    // Normalize the unlimited normal
    if (unlimitedNormal) {
        *unlimitedNormal /= norm2(*unlimitedNormal);
    }

    // Normalize the unlimited normal
    if (limitedNormal) {
        *limitedNormal /= norm2(*limitedNormal);
    }
}

/*!
 * Check if the factes have a consistent orientation.
 *
 * \return Returns true if the factes have a consistent orientation, false otherwise.
 */
bool SurfaceKernel::isCellOrientationConsistent() const
{
    CellConstIterator internalBegin = internalCellConstBegin();
    CellConstIterator internalEnd   = internalCellConstEnd();

    std::size_t nUnprocessed = getInternalCellCount();
    PiercedStorage<bool> alreadyProcessed(1, &(getCells()));
    alreadyProcessed.fill(false);

    bool consistent = true;
    std::vector<std::size_t> processRawList;
    for (CellConstIterator seedItr = internalBegin; seedItr != internalEnd; ++seedItr) {
        // Get seed information
        long seedRawId = seedItr.getRawIndex();
        if (alreadyProcessed.rawAt(seedRawId)) {
            continue;
        }

        // Compare seed orientation with the orientation of its neighbours
        processRawList.assign(1, seedRawId);
        while (!processRawList.empty()) {
            // Get cell to process
            std::size_t cellRawId = processRawList.back();
            processRawList.pop_back();

            // Skip cells already processed
            if (alreadyProcessed.rawAt(cellRawId)) {
                continue;
            }

            // Process cell
            CellConstIterator cellItr = getCells().rawFind(cellRawId);
            const Cell &cell = *cellItr;
            int nCellFaces = cell.getFaceCount();
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = cell.getAdjacencyCount(face);
                const long *faceAdjacencies = cell.getAdjacencies(face);
                for (int i = 0; i < nFaceAdjacencies; ++i) {
                    // Neighbout information
                    long neighId = faceAdjacencies[i];
                    CellConstIterator neighItr = getCells().find(neighId);
                    const Cell &neigh = getCell(neighId);

                    // Check neighbour consistency
                    int neighFace = findAdjoinNeighFace(cell, face, neigh);
                    bool isNeighConsistent = haveSameOrientation(cell, face, neigh, neighFace);
                    if (!isNeighConsistent) {
                        consistent = false;
                        break;
                    }

                    // Add neighbour to the process list
                    //
                    // Only internal cells needs to be processed.
                    if (neigh.isInterior()) {
                        std::size_t neighRawId = neighItr.getRawIndex();
                        if (!alreadyProcessed.rawAt(neighRawId)) {
                            processRawList.push_back(neighRawId);
                        }
                    }
                }

                // Stop processing the cell is the orientation is not consistent
                if (!consistent) {
                    break;
                }
            }

            // Cell processing has been completed
            alreadyProcessed.rawSet(cellRawId, true);
            --nUnprocessed;
            if (nUnprocessed == 0) {
                break;
            }

            // Stop processing the patch is the orientation is not consistent
            if (!consistent) {
                break;
            }
        }

        // Stop processing the patch is the orientation is not consistent
        if (!consistent) {
            break;
        }

        // Check if all cells have been processed
        if (nUnprocessed == 0) {
            break;
        }
    }

#if BITPIT_ENABLE_MPI==1
    // Comunicate consistency flag among the partitions
    if (isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &consistent, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
    }
#endif

    return consistent;
}

/*!
 * Adjust the orientation of all facets of the local partition according to the
 * orientation of the first facet stored.
 *
 * The routine checks whether the surface is orientable; if the surface is not
 * orientable false is returned (with some normals flipped until the
 * orientablity condition has been violated), otherwise true is returned (with
 * all facet orientations coherent with the reference facet).
 *
 * \return Returns true if the surface is orientable, false otherwise.
 */
bool SurfaceKernel::adjustCellOrientation()
{
    // Early return if the patch is empty.
    if (empty()) {
        return false;
    }

    // Get the first facet
    long firstFacetId = Cell::NULL_ID;
    if (getCellCount() > 0) {
        firstFacetId = getCells().begin()->getId();
    }

#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        int patchRank = getRank();

        int firstFacetRank = std::numeric_limits<int>::max();
        if (firstFacetId >= 0) {
            firstFacetRank = patchRank;
        }

        MPI_Allreduce(MPI_IN_PLACE, &firstFacetRank, 1, MPI_INT, MPI_MIN, getCommunicator());
        if (patchRank != firstFacetRank) {
            firstFacetId = Cell::NULL_ID;
        }
    }
#endif

    // Adjust cell orientation according to the orientation of the first facet
    return adjustCellOrientation(firstFacetId);
}


/*!
 * Adjust the orientation of all facets of the local partition according to the
 * orientation of a given facet.
 *
 * The routine checks whether the surface is orientable; if the surface is not
 * orientable false is returned (with some normals flipped until the
 * orientablity condition has been violated), otherwise true is returned (with
 * all facet orientations coherent with the reference facet).
 *
 * \param[in] seed is the id of reference facet. For parallel computaions, if
 * a partition does not know how to orient its facets, Element::NULL_ID should
 * be passed for these partitions
 * \param[in] invert controls if the orientation should be the same of the
 * seed of should be inverted with respect to the seed
 * \return Returns true if the surface is orientable, false otherwise.
 */
bool SurfaceKernel::adjustCellOrientation(long seed, bool invert)
{
#if BITPIT_ENABLE_MPI==1
    // Check if the communications are needed
    bool ghostExchangeNeeded = isPartitioned();

    // Initialize ghost communications
    std::unordered_set<long> flippedCells;
    std::unique_ptr<DataCommunicator> ghostCommunicator;
    if (ghostExchangeNeeded) {
        ghostCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(getCommunicator()));

        for (auto &entry : getGhostCellExchangeSources()) {
            const int rank = entry.first;
            const std::vector<long> &sources = entry.second;
            ghostCommunicator->setSend(rank, 2 * sources.size() * sizeof(unsigned char));
        }

        for (auto &entry : getGhostCellExchangeTargets()) {
            const int rank = entry.first;
            const std::vector<long> &targets = entry.second;
            ghostCommunicator->setRecv(rank, 2 * targets.size() * sizeof(unsigned char));
        }
    }
#endif

    // Initialize list of cells to process
    std::vector<std::size_t> processRawList;
    if (seed != Element::NULL_ID) {
        assert(getCells().exists(seed));
        processRawList.push_back(getCells().find(seed).getRawIndex());
        if (invert) {
            flipCellOrientation(seed);
#if BITPIT_ENABLE_MPI==1
            if (ghostExchangeNeeded) {
                flippedCells.insert(seed);
            }
#endif
        }
    }

    // Orient the surface
    PiercedStorage<bool> alreadyProcessed(1, &(getCells()));
    alreadyProcessed.fill(false);

    bool orientable = true;
    while (true) {
        while (!processRawList.empty()) {
            // Get cell to process
            std::size_t cellRawId = processRawList.back();
            processRawList.pop_back();

            // Skip cells already processed
            if (alreadyProcessed.rawAt(cellRawId)) {
                continue;
            }

            // Process cell
            CellConstIterator cellItr = getCells().rawFind(cellRawId);
            const Cell &cell = *cellItr;
            int nCellFaces = cell.getFaceCount();
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = cell.getAdjacencyCount(face);
                const long *faceAdjacencies = cell.getAdjacencies(face);
                for (int i = 0; i < nFaceAdjacencies; ++i) {
                    // Neighbout information
                    long neighId = faceAdjacencies[i];
                    CellConstIterator neighItr = getCells().find(neighId);
                    std::size_t neighRawId = neighItr.getRawIndex();
                    const Cell &neigh = *neighItr;

                    // Skip ghost neighbours
                    if (!neigh.isInterior()) {
                        continue;
                    }

                    // Check neighbour orientation
                    //
                    // If the neighour has not been processed yet and its orientation
                    // is wrong we flip the orientation of the facet.If the neighour
                    // has already been processed and its orientation is wrong the
                    // surface is not orientable.
                    int neighFace = findAdjoinNeighFace(cell, face, neigh);
                    bool isNeighOriented = haveSameOrientation(cell, face, neigh, neighFace);
                    if (!isNeighOriented) {
                        if (!alreadyProcessed.rawAt(neighRawId)) {
                            flipCellOrientation(neighId);
#if BITPIT_ENABLE_MPI==1
                            if (ghostExchangeNeeded) {
                                flippedCells.insert(neighId);
                            }
#endif
                        } else {
                            orientable = false;
                            break;
                        }
                    }

                    // Add neighbour to the process list
                    if (!alreadyProcessed.rawAt(neighRawId)) {
                        processRawList.push_back(neighRawId);
                    }
                }

                // Stop processing the cell if the patch is not orientable
                if (!orientable) {
                    break;
                }
            }

            // Stop processing if the patch is not orientable
            if (!orientable) {
                break;
            }

            // Cell processing has been completed
            alreadyProcessed.rawSet(cellRawId, true);
        }

        // If the surface is globally not orientable we can exit
#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &orientable, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
        }
#endif

        if (!orientable) {
            break;
        }

#if BITPIT_ENABLE_MPI==1
        // Add ghost to the process list
        //
        // A ghost should be added to the process list if:
        //  - its orientation was flipped by the owner;
        //  - it has been processed by the owner but is hass't yet processed by
        //    the current process.
        if (ghostExchangeNeeded) {
            // Exchange ghost data
            ghostCommunicator->startAllRecvs();

            for(auto &entry : getGhostCellExchangeSources()) {
                int rank = entry.first;
                SendBuffer &buffer = ghostCommunicator->getSendBuffer(rank);
                for (long cellId : entry.second) {
                    std::size_t cellRawId = getCells().find(cellId).getRawIndex();

                    buffer << static_cast<unsigned char>(alreadyProcessed.rawAt(cellRawId));
                    if (flippedCells.count(cellId) > 0) {
                        buffer << static_cast<unsigned char>(1);
                    } else {
                        buffer << static_cast<unsigned char>(0);
                    }
                }
                ghostCommunicator->startSend(rank);
            }

            int nPendingRecvs = ghostCommunicator->getRecvCount();
            while (nPendingRecvs != 0) {
                int rank = ghostCommunicator->waitAnyRecv();
                RecvBuffer buffer = ghostCommunicator->getRecvBuffer(rank);

                for (long cellId : getGhostCellExchangeTargets(rank)) {
                    unsigned char processed;
                    buffer >> processed ;

                    unsigned char ghostFlipped;
                    buffer >> ghostFlipped;

                    if (ghostFlipped == 1 || processed == 1) {
                        std::size_t cellRawId = getCells().find(cellId).getRawIndex();

                        if (ghostFlipped == 1) {
                            flipCellOrientation(cellId);
                            processRawList.push_back(cellRawId);
                            alreadyProcessed.rawSet(cellRawId, false);
                        } else if (processed == 1) {
                            if (!alreadyProcessed.rawAt(cellRawId)) {
                                processRawList.push_back(cellRawId);
                            }
                        }
                    }
                }

                --nPendingRecvs;
            }

            ghostCommunicator->waitAllSends();

            // Clear list of flippedCells cells
            flippedCells.clear();
        }
#endif

        // Detect if the orientation is completed
        bool completed = processRawList.empty();
#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &completed, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
        }
#endif

        if (completed) {
            break;
        }
    }

    return orientable;
}

/*!
 * Check if the specified facets have the same orientation.
 *
 * The functions assumes that the specified facets share (at least partially)
 * a face, however it's up to the caller make sure this assumption is valid.
 *
 *  This guarantees  that is possible  are either coincident
 * (conformal faces), or lay on the same plane and one face contains the
 * other (non-conformal faces).
 *
 * If the faces share at least two vertices, it is possible to detect their
 * relative orientation only checking the order of the vertices, otherwise
 * is is necessary to evaluate the normals.
 *
 * \param[in] cell_A is the first cell
 * \param[in] face_A is the face on the first cell
 * \param[in] cell_B is the second cell
 * \param[in] face_B is the face on the second cell
 * \result Returns true if the two facets have the same orientation, false
 * otherwise.
*/
bool SurfaceKernel::haveSameOrientation(const Cell &cell_A, int face_A, const Cell &cell_B, int face_B) const
{
    //
    // Cell information
    //
    long cellId_A = cell_A.getId();
    ConstProxyVector<long> cellVertexIds_A = cell_A.getVertexIds();
    std::size_t nCellVertices_A = cellVertexIds_A.size();

    long cellId_B = cell_B.getId();
    ConstProxyVector<long> cellVertexIds_B = cell_B.getVertexIds();
    std::size_t nCellVertices_B = cellVertexIds_B.size();

    assert(cell_A.getDimension() == cell_B.getDimension());
    int cellDimension = cell_A.getDimension();
    assert(cellDimension < 3);

    //
    // Check relative orientation using vertices
    //

    // If the facets share at least two consecutive vertices, it is possible
    // to detect relative orientation checking the order in which these two
    // vertices appear while traversing the list of vertices of the facets.
    // If the vertices appear in the same order, the facets have opposite
    // orientation, if the vertices appear in reversed order, the facets have
    // the same orientation.
    for (std::size_t i = 0; i < nCellVertices_A; ++i) {
        long vertexId_A = cellVertexIds_A[getFacetOrderedLocalVertex(cell_A, i)];
        for (std::size_t j = 0; j < nCellVertices_B; ++j) {
            long vertexId_B = cellVertexIds_B[getFacetOrderedLocalVertex(cell_B, j)];
            if (vertexId_A == vertexId_B) {
                if (cellDimension == 2) {
                    long previousVertexId_A = cellVertexIds_A[getFacetOrderedLocalVertex(cell_A, (i - 1 + nCellVertices_A) % nCellVertices_A)];
                    long previousVertexId_B = cellVertexIds_B[getFacetOrderedLocalVertex(cell_B, (j - 1 + nCellVertices_B) % nCellVertices_B)];

                    long nextVertexId_A = cellVertexIds_A[getFacetOrderedLocalVertex(cell_A, (i + 1) % nCellVertices_A)];
                    long nextVertexId_B = cellVertexIds_B[getFacetOrderedLocalVertex(cell_B, (j + 1) % nCellVertices_B)];

                    if (nextVertexId_A == nextVertexId_B || previousVertexId_A == previousVertexId_B) {
                        return false;
                    } else if (nextVertexId_A == previousVertexId_B || previousVertexId_A == nextVertexId_B) {
                        return true;
                    }
                } else if (cellDimension == 1) {
                    return (i != j);
                } else {
                    throw std::runtime_error("Unable to identify if the factes has the same orientation.");
                }
            }
        }
    }

    //
    // Check relative orientation using facet normals
    //

    // The faces doesn't have enough common vertices to check the orientation
    // using topological checks, we need a geometric check. Starting from the
    // first vertex of the face, we identify three non-collinear vertices and
    // we detect their winding. The two facets have the same orientation if
    // they have the same winding.
    //
    // To detect the winding of the vertices, we evaluate the normal of the
    // plane defined by the vertices and then we calculate the projection of
    // that vector along a direction (called winding direction) which is
    // chosen equal for both the cells. The sign of the projection defines
    // the winding of the vertices. We don't care if the winding of a single
    // cell is clockwise or counter-clockwise, all we care about is if the
    // two cells have the same winding.

    // Pseudo normal directions
    //
    // For each face, we identify three non-collinear vertices and we evaluate
    // the normal direction of the plane defined by those three vertices.
    std::size_t nMaxCellVertices = std::max(nCellVertices_A, nCellVertices_B);
    BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nMaxCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);

    ConstProxyVector<int> faceLocalVertexIds_A = cell_A.getFaceLocalVertexIds(face_A);
    ConstProxyVector<int> faceLocalVertexIds_B = cell_B.getFaceLocalVertexIds(face_B);

    std::array<double, 3> pseudoNormal_A;
    std::array<double, 3> pseudoNormal_B;
    for (int i = 0; i < 2; ++i) {
        std::array<double, 3> *pseudoNormal;
        long firstFaceVertexLocalId;
        std::size_t nCellVertices;
        if (i == 0) {
            pseudoNormal = &pseudoNormal_A;
            firstFaceVertexLocalId = faceLocalVertexIds_A[0];
            nCellVertices = nCellVertices_A;
            getCellVertexCoordinates(cellId_A, vertexCoordinates);
        } else {
            pseudoNormal = &pseudoNormal_B;
            firstFaceVertexLocalId = faceLocalVertexIds_B[0];
            nCellVertices = nCellVertices_B;
            getCellVertexCoordinates(cellId_B, vertexCoordinates);
        }
        assert(nCellVertices >= 3);

        const std::array<double, 3> &V_0 = vertexCoordinates[(firstFaceVertexLocalId + 0) % nCellVertices];
        const std::array<double, 3> &V_1 = vertexCoordinates[(firstFaceVertexLocalId + 1) % nCellVertices];
        for (std::size_t k = 2; k < nCellVertices; ++k) {
            const std::array<double, 3> &V_2 = vertexCoordinates[(firstFaceVertexLocalId + k) % nCellVertices];

            *pseudoNormal = crossProduct(V_0 - V_1, V_2 - V_1);
            if (!utils::DoubleFloatingEqual()(norm2(*pseudoNormal), 0.)) {
                break;
            } else if (k == (nCellVertices - 1)) {
                throw std::runtime_error("Unable to identify the non-collinear vertices.");
            }
        }
    }

    // Direction for evaluating the winding
    //
    // This direction has to be the same for both the cells and should not be
    // perpendicular to the pseudo normals.
    std::array<double, 3> windingDirection = {{0., 0., 1.}};
    for (int d = 2; d <= 0; --d) {
        windingDirection = {{0., 0., 0.}};
        windingDirection[d] = 1.;

        double pseudNormalProjection_A = dotProduct(pseudoNormal_A, windingDirection);
        double pseudNormalProjection_B = dotProduct(pseudoNormal_B, windingDirection);
        if (!utils::DoubleFloatingEqual()(pseudNormalProjection_A, 0.) && !utils::DoubleFloatingEqual()(pseudNormalProjection_B, 0.)) {
            break;
        } else if (d == 0) {
            throw std::runtime_error("Unable to identify the relative orientation of the cells.");
        }
    }

    // Identify the winding of the cells
    int cellWinding_A = (dotProduct(pseudoNormal_A, windingDirection) > 0.);
    int cellWinding_B = (dotProduct(pseudoNormal_B, windingDirection) > 0.);

    return (cellWinding_A == cellWinding_B);
}

/*!
 * Flips the orientation of a cell.
 *
 * \param[in] id is the id of cell that will be flipped
 */
void SurfaceKernel::flipCellOrientation(long id)
{
    // Cell information
    Cell &cell = getCell(id);

    //
    // Flip connectivity
    //
    // In order to flip cell orientation, we need to reverse the order of
    // the vertices:
    //
    //  Vertices             0     1     2     3     4
    //
    //  Flipped vertices     4     3     2     1     0
    //
    //  This means we have to invert the connectivity of the cell.
    long *connectivity = cell.getConnect();

    if (cell.getType() != ElementType::PIXEL) {
        int connectBegin = 0;
        if (cell.getType() == ElementType::POLYGON) {
            ++connectBegin;
        }

        int nCellVertices = cell.getVertexCount();
        for (int i = 0; i < (int) std::floor(nCellVertices / 2); ++i) {
            std::swap(connectivity[connectBegin + i], connectivity[connectBegin + nCellVertices - i - 1]);
        }
    } else {
        std::swap(connectivity[0], connectivity[1]);
        std::swap(connectivity[2], connectivity[3]);
    }

    //
    // Flip adjacencies and interfaces
    //
    // Adjacencies and interfaces are associated with the faces of the cells,
    // since we have flipped the connectivity we need to flip also adjacencies
    // and interfaces.
    //
    // For two-dimensional elements faces are defined between two vertices:
    //
    //  Vertices             0     1     2     3     4
    //  Faces                   0     1     2     3     4
    //
    //  Flipped vertices     4     3     2     1     0
    //  Flipped faces           3     2     1     0     4
    //
    // We have reversed the order of the vertices, this means that we have
    // to invert the order of all the faces but the lust one. Pixel elements
    // should be handled separately.
    //
    // For one-dimensional elements faces are defined on the vertices itself:
    //
    //  Vertices             0     1
    //  Faces                0     1
    //
    //  Flipped vertices     1     0
    //  Flipped faces        1     0
    //
    // We have reversed the order of the vertices, this means that we have
    // to invert the order of all the faces.
    int nCellAdjacencies = cell.getAdjacencyCount();
    int nCellInterfaces = cell.getInterfaceCount();
    int nCellFaces = cell.getFaceCount();

    FlatVector2D<long> flippedCellAdjacencies;
    FlatVector2D<long> flippedCellInterfaces;
    flippedCellAdjacencies.reserve(nCellFaces, nCellAdjacencies);
    flippedCellInterfaces.reserve(nCellFaces, nCellInterfaces);

    if (cell.getType() != ElementType::PIXEL) {
        int flipOffset;
        if (cell.getDimension() == 1) {
            flipOffset = 0;
        } else {
            flipOffset = 1;
        }


        for (int i = 0; i < nCellFaces - flipOffset; ++i) {
            int nFaceAdjacencies = cell.getAdjacencyCount(nCellFaces - i - 1 - flipOffset);
            const long *faceAdjacencies = cell.getAdjacencies(nCellFaces - i - 1 - flipOffset);
            flippedCellAdjacencies.pushBack(nFaceAdjacencies, faceAdjacencies);

            int nFaceInterfaces = cell.getInterfaceCount(nCellFaces - i - 1 - flipOffset);
            const long *faceInterfaces = cell.getInterfaces(nCellFaces - i - 1 - flipOffset);
            flippedCellInterfaces.pushBack(nFaceInterfaces, faceInterfaces);
        }

        if (flipOffset == 1) {
            int nLastFaceAdjacencies = cell.getAdjacencyCount(nCellFaces - 1);
            const long *lastFaceAdjacencies = cell.getAdjacencies(nCellFaces - 1);
            flippedCellAdjacencies.pushBack(nLastFaceAdjacencies, lastFaceAdjacencies);

            int nLastFaceInterfaces = cell.getInterfaceCount(nCellFaces - 1);
            const long *lastFaceInterfaces = cell.getInterfaces(nCellFaces - 1);
            flippedCellInterfaces.pushBack(nLastFaceInterfaces, lastFaceInterfaces);
        }
    } else {
        std::array<int, 4> orderedFaces = {{1, 0, 2, 3}};
        for (int face : orderedFaces) {
            int nFaceAdjacencies = cell.getAdjacencyCount(face);
            const long *faceAdjacencies = cell.getAdjacencies(face);
            flippedCellAdjacencies.pushBack(nFaceAdjacencies, faceAdjacencies);

            int nFaceInterfaces = cell.getInterfaceCount(face);
            const long *faceInterfaces = cell.getInterfaces(face);
            flippedCellInterfaces.pushBack(nFaceInterfaces, faceInterfaces);
        }
    }

    // Set flipped adjacencies
    if (nCellAdjacencies > 0) {
        cell.setAdjacencies(std::move(flippedCellAdjacencies));
    }

    // Set flipped interfaces
    //
    // After setting the flipped interfaces, we also need to update the face
    // associated with the interfaces.
    if (nCellInterfaces > 0) {
        cell.setInterfaces(std::move(flippedCellInterfaces));

        const long *cellInterfaces = cell.getInterfaces();
        int nCellInterfaces = cell.getInterfaceCount();
        for (int i = 0; i < nCellInterfaces; ++i) {
            long interfaceId = cellInterfaces[i];
            Interface &interface = getInterface(interfaceId);

            long owner = interface.getOwner();
            if (owner == id) {
                int ownerFace = interface.getOwnerFace();
                if (ownerFace != nCellFaces - 1) {
                    ownerFace = nCellFaces - 2 - ownerFace;
                    interface.setOwner(id, ownerFace);
                }
            } else {
                int neighFace = interface.getNeighFace();
                if (neighFace != nCellFaces - 1) {
                    neighFace = nCellFaces - 2 - neighFace;
                    interface.setNeigh(id,neighFace);
                }
            }
        }
    }
}

/*!
 * Display quality stats for the mesh currently stored.
 * 
 * \param[in,out] out output stream where stats will be printed out
 * \param[in] padding (default = 0) number of trailing spaces
*/
void SurfaceKernel::displayQualityStats(std::ostream& out, unsigned int padding) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    std::string         indent(padding, ' ');

    // Counters
    // none

    // ====================================================================== //
    // DISPLAY STATS                                                          //
    // ====================================================================== //

    // Distribution for aspect ratio ---------------------------------------- //
    out << indent << "Aspect ratio distribution ---------------" << std::endl;
    {
        // Scope variables
        long                            count;
        std::vector<double>             bins, hist;

        // Compute histogram
        hist = computeHistogram(&SurfaceKernel::evalAspectRatio, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "AR", out, padding + 2);
    }

    // Distribution of min. angle ------------------------------------------- //
    out << indent << "Min. angle distribution -----------------" << std::endl;
    {
        // Scope variables
        long                            count;
        std::vector<double>             bins, hist;

        // Compute histogram
        hist = computeHistogram(&SurfaceKernel::evalMinAngleAtVertex, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "min. angle", out, padding + 2);
    }

    // Distribution of min. angle ------------------------------------------- //
    out << indent << "Max. angle distribution -----------------" << std::endl;
    {
        // Scope variables
        long                            count;
        std::vector<double>             bins, hist;

        // Compute histogram
        hist = computeHistogram(&SurfaceKernel::evalMaxAngleAtVertex, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "max. angle", out, padding + 2);
    }
    
    return;
}

/*!
 * Compute the histogram showing the distribution of aspect ratio among elements
 *
 * \param[in] funct_ pointer to member function to be used to compute stats
 * \param[in,out] bins bins used to construct histogram. If an empty vector is passed
 * as input, uniform bins will be generated in the interval [1.0, 5.0)
 * \param[in] count population size used to build histogram (depends on the selection mask
 * used). For instant, if the selection mask is set to the value SelectionMask::SELECT_QUAD |
 * SelectionMask::SELECT_TRIANGLE, only quad and tria element will be considered and count
 * will returns the number of such elements.
 * \param[in] n_intervals (default = 8), number of intervals to be used for histogram construction.
 * If bins is a non-empty vector, the number of intervals will be set equal to
 * bins.size()-1
 * \param[in] mask (default = SelectionMask::ALL) selection mask for element type
 *
 * \result on output returns a vector of dimensions n_internvals+2, storing the
 * histogram for the distribution of aspect ratio among cells.
 * hist[0] stores the % of element having aspect ratio below bins[0]
 * ...
 * hist[i] stores the % of element having aspect ratio between [bins[i], bins[i+1])
 * ...
 * hist[n_internvals+1] stores the % of element having aspect ratio above
 * bins[n_intervals].
*/
std::vector<double> SurfaceKernel::computeHistogram(eval_f_ funct_, std::vector<double> &bins, long &count, int n_intervals, unsigned short mask) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    double                      m_, M_;
    std::vector<double>         hist;

    // ====================================================================== //
    // BUILD BINS FOR HISTOGRAM                                               //
    // ====================================================================== //

    // Resize bin data structure
    if (bins.size() < 2) {
        double          db;
        bins.resize(n_intervals+1);
        m_ = 1.0;
        M_ = 3.0;
        db  = (M_ - m_)/double(n_intervals);
        for (int i = 0; i < n_intervals+1; ++i) {
            bins[i] = m_ + db * double(i);
        } //next i
    }
    else {
        n_intervals = bins.size()-1;
    }

    // Resize histogram data structure
    hist.resize(n_intervals+2, 0.0);

    // ====================================================================== //
    // COMPUTE HISTOGRAM                                                      //
    // ====================================================================== //

    // Scope variables
    long                id;
    int                 i;
    int                 dummy;
    double              ar;

    // Compute histogram
    count = 0;
    for (auto &cell_ : m_cells) {
        id = cell_.getId();
        if (compareSelectedTypes(mask, cell_.getType())) {
            ar = (this->*funct_)(id, dummy);
            i = 0;
            while ((bins[i] - ar < 1.0e-5) && (i < n_intervals+1)) ++i;
            ++hist[i];
            ++count;
        }
    } //next cell_

    // Normalize histogram
    hist = 100.0 * hist/double(count);

    return(hist);
}

/*!
 * Display histogram in a nicely formatted form.
 * 
 * \param[in] count population size used for histogram computation
 * \param[in] bins bins used for histogram computation
 * \param[in] hist histogram value
 * \param[in] stats_name statistcs name
 * \param[in,out] out output stream where stats will be printed out
 * \param[in] padding (default = 0) number of trailing spaces
*/
void SurfaceKernel::displayHistogram(
    long                         count,
    const std::vector<double>   &bins,
    const std::vector<double>   &hist,
    const std::string           &stats_name,
    std::ostream                &out,
    unsigned int                 padding
) const {
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                 nbins = bins.size();
    size_t              n_fill;
    std::string         indent(padding, ' ');
    std::stringstream   ss;

    // Counters
    int                 i;

    // ====================================================================== //
    // DISPLAY HISTOGRAM                                                      //
    // ====================================================================== //

    // Set stream properties
    ss << std::setprecision(3);

    // Display histograms
    out << indent << "poll size: " << count << std::endl;
    out << indent;
    ss << hist[0];
    n_fill = ss.str().length();
    ss << std::string(std::max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << "%, with          " << stats_name << " < ";
    ss.str("");
    ss << bins[0];
    out << ss.str() << std::endl;
    ss.str("");
    for (i = 0; i < nbins-1; ++i) {
        out << indent;
        ss << hist[i+1];
        n_fill = ss.str().length();
        ss << std::string(std::max(size_t(0),  6 - n_fill), ' ');
        out << ss.str() << "%, with ";
        ss.str("");
        ss << bins[i];
        n_fill = ss.str().length();
        ss << std::string(std::max(size_t(0), 6 - n_fill), ' ');
        out << ss.str() << " < " << stats_name << " < ";
        ss.str("");
        ss << bins[i+1];
        out << ss.str() << std::endl;
        ss.str("");
    } //next i
    out << indent;
    ss << hist[i+1];
    n_fill = ss.str().length();
    ss << std::string(std::max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << "%, with ";
    ss.str("");
    ss << bins[i];
    n_fill = ss.str().length();
    ss << std::string(std::max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << " < " << stats_name << std::endl;
    ss.str("");

    return;
}

/*!
 * Given and element type and the selection mask, returns true if the element type
 * is accounted for in the selection mask.
 * 
 * \param[in] mask_ selection mask_
 * \param[in] type_ element type
*/
bool SurfaceKernel::compareSelectedTypes(unsigned short mask_, ElementType type_) const
{
    unsigned short       masked = m_selectionTypes.at(type_);
    return ( (mask_ & masked) == masked );
}

/*!
	Get the counter-clockwise ordered list of vertex for the specified facet.

	\param facet is the facet
	\result The counter-clockwise ordered list of vertex for the specified facet.
*/
ConstProxyVector<long> SurfaceKernel::getFacetOrderedVertexIds(const Cell &facet) const
{
    ConstProxyVector<long> vertexIds = facet.getVertexIds();
    if (areFacetVerticesOrdered(facet)) {
        return vertexIds;
    }

    std::size_t nVertices = vertexIds.size();
    ConstProxyVector<long> orderedVertexIds(ConstProxyVector<long>::INTERNAL_STORAGE, nVertices);
    ConstProxyVector<long>::storage_pointer orderedVertexIdsStorage = orderedVertexIds.storedData();
    for (std::size_t k = 0; k < nVertices; ++k) {
        int vertex = getFacetOrderedLocalVertex(facet, k);
        orderedVertexIdsStorage[k] = vertexIds[vertex];
    }

    return orderedVertexIds;
}

/*!
 * Check if the vertices of the specified facet are counter-clockwise ordered.
 *
 * \param[in] facet is the facet
 * \result Return true if the vertices of the specified facet are counter-clockwise
 * ordered, false otherwise.
*/
bool SurfaceKernel::areFacetVerticesOrdered(const Cell &facet) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 1) {
        return true;
    }

    // Check if vertices are ordered
    ElementType facetType = facet.getType();
    switch (facetType) {

    case (ElementType::POLYGON):
    {
        return true;
    }

    default:
    {
        assert(facetType != ElementType::UNDEFINED);
        const Reference2DElementInfo &facetInfo = static_cast<const Reference2DElementInfo &>(facet.getInfo());
        return facetInfo.areVerticesCCWOrdered();
    }

    }
}

/*!
 * Get the local index of the vertex occupying the n-th position in the
 * counter-clockwise ordered list of vertex ids.
 *
 * \param[in] facet is the facet
 * \param[in] n is the requested position
 * \result The local index of the vertex occupying the n-th position in the
 * counter-clockwise ordered list of vertex ids.
*/
int SurfaceKernel::getFacetOrderedLocalVertex(const Cell &facet, std::size_t n) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 1) {
        return n;
    }

    // Get the request index
    ElementType facetType = facet.getType();
    switch (facetType) {

    case (ElementType::POLYGON):
    {
        return n;
    }

    default:
    {
        assert(facetType != ElementType::UNDEFINED);
        const Reference2DElementInfo &facetInfo = static_cast<const Reference2DElementInfo &>(facet.getInfo());
        return facetInfo.getCCWOrderedVertex(n);
    }

    }
}

/*!
	Get the counter-clockwise ordered list of edge for the specified facet.

	\param facet is the facet
	\result The counter-clockwise ordered list of edge for the specified facet.
*/
ConstProxyVector<long> SurfaceKernel::getFacetOrderedEdgeIds(const Cell &facet) const
{
    std::size_t nEdges = facet.getFaceCount();
    ConstProxyVector<long> orderedEdgeIds(ConstProxyVector<long>::INTERNAL_STORAGE, nEdges);
    ConstProxyVector<long>::storage_pointer orderedEdgeIdsStorage = orderedEdgeIds.storedData();
    for (std::size_t k = 0; k < nEdges; ++k) {
        orderedEdgeIdsStorage[k] = getFacetOrderedLocalEdge(facet, k);
    }

    return orderedEdgeIds;
}

/*!
 * Check if the edges of the specified facet are counter-clockwise ordered.
 *
 * \param[in] facet is the facet
 * \result Return true if the edges of the specified facet are counter-clockwise
 * ordered, false otherwise.
*/
bool SurfaceKernel::areFacetEdgesOrdered(const Cell &facet) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 1) {
        return true;
    }

    // Check if edges are ordered
    ElementType facetType = facet.getType();
    switch (facetType) {

    case (ElementType::POLYGON):
    {
        return true;
    }

    default:
    {
        assert(facetType != ElementType::UNDEFINED);
        const Reference2DElementInfo &facetInfo = static_cast<const Reference2DElementInfo &>(facet.getInfo());
        return facetInfo.areFacesCCWOrdered();
    }

    }
}

/*!
 * Get the local index of the edge occupying the n-th position in the
 * counter-clockwise ordered list of edge ids.
 *
 * \param[in] facet is the facet
 * \param[in] n is the requested position
 * \result The local index of the edge occupying the n-th position in the
 * counter-clockwise ordered list of edge ids.
*/
int SurfaceKernel::getFacetOrderedLocalEdge(const Cell &facet, std::size_t n) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 1) {
        return n;
    }

    // Get the request index
    ElementType facetType = facet.getType();
    switch (facetType) {

    case (ElementType::POLYGON):
    {
        return n;
    }

    default:
    {
        assert(facetType != ElementType::UNDEFINED);
        const Reference2DElementInfo &facetInfo = static_cast<const Reference2DElementInfo &>(facet.getInfo());
        return facetInfo.getCCWOrderedFace(n);
    }

    }
}

}
