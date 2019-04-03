/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#include <bitpit_communications.hpp>
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

/*!
	Creates a new patch.

	\param expert if true, the expert mode will be enabled
*/
SurfaceKernel::SurfaceKernel(bool expert)
	: PatchKernel(expert)
{
	initialize();
}

/*!
	Creates a new patch.

	\param patch_dim is the dimension of the patch
	\param space_dim is the dimension of the space
	\param expert if true, the expert mode will be enabled
*/
SurfaceKernel::SurfaceKernel(const int &patch_dim, const int& space_dim, bool expert)
	: PatchKernel(patch_dim, expert)
{
    initialize();

    // Set the sapce dimension
    setSpaceDimension(space_dim);
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param patch_dim is the dimension of the patch
	\param space_dim is the dimension of the space
	\param expert if true, the expert mode will be enabled
*/
SurfaceKernel::SurfaceKernel(const int &id, const int &patch_dim, const int& space_dim, bool expert)
	: PatchKernel(id, patch_dim, expert)
{
    m_spaceDim = space_dim;
}

/*!
	Destroys the patch.
*/
SurfaceKernel::~SurfaceKernel()
{

}

/*!
	Initialize the patch
*/
void SurfaceKernel::initialize()
{
    // Space dimension
    m_spaceDim = -1;
}

/*!
	Sets the dimension of the working space.

	\param dimension the dimension of the working patch
*/
void SurfaceKernel::setSpaceDimension(int dimension)
{
    // If the dimension was already assigned, reset the patch
    if (m_spaceDim > 0 && m_spaceDim != dimension) {
        reset();
    }

    // Set the dimension
    m_spaceDim = dimension;
}

/*!
 * Returns the number of dimensions of the working space (set at patch construction)
 *
 * \result The number of dimensions.
 */
int SurfaceKernel::getSpaceDimension(void) const
{
    return(m_spaceDim);
}

/*!
        Evaluates the characteristic size of the specified cell.

        \param id is the id of the cell
        \result The characteristic size of the specified cell.
*/
double SurfaceKernel::evalCellSize(const long &id) const
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE CELL SIZE                                                      //
    // ====================================================================== //
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
double SurfaceKernel::evalCellArea(const long &id) const
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
    array<double, 3>            d1, d2;
    if (cell_->getType() == ElementType::TRIANGLE) {
        d1 = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        d2 = m_vertices[cellVertexIds[2]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        return (0.5*norm2(crossProduct(d1, d2)));
    }
    else {
        int                     nvert = cell_->getVertexCount();
        int                     next, prev;
        double                  coeff = 0.25;
        double                  area = 0.0;
        for (int i = 0; i < nvert; ++i) {
            next = (i + 1) % nvert;
            prev = (nvert + i - 1) % nvert;
            d1 = m_vertices[cellVertexIds[next]].getCoords() - m_vertices[cellVertexIds[i]].getCoords();
            d2 = m_vertices[cellVertexIds[prev]].getCoords() - m_vertices[cellVertexIds[i]].getCoords();
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
double SurfaceKernel::evalEdgeLength(const long &cellId, const int &edgeId) const
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
double SurfaceKernel::evalMinEdgeLength(const long &id, int &edge_id) const
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
double SurfaceKernel::evalMaxEdgeLength(const long &id, int &edge_id) const
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
 * Evaluate the angle at specified vertex for a cell with specified ID.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, a returns zero.
 * 
 * \param[in] id cell global ID
 * \param[in] vertex_id vertex local ID
*/
double SurfaceKernel::evalAngleAtVertex(const long &id, const int &vertex_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    const Cell                   *cell_ = &m_cells[id];
    ConstProxyVector<long>       cellVertexIds = cell_->getVertexIds();

    // Counters

    // ====================================================================== //
    // EVALUATE ANGLE AT SPECIFIED VERTEX                                     //
    // ====================================================================== //
    if ((cell_->getType() == ElementType::UNDEFINED)
     || (cell_->getType() == ElementType::VERTEX)) return 0.0;
    if (cell_->getType() == ElementType::LINE) {
        if (m_spaceDim - getDimension() == 1)   return BITPIT_PI;
        else                                    return 0.0;
    }

    int                          n_vert = cell_->getVertexCount();
    int                          prev = (n_vert + vertex_id - 1) % n_vert;
    int                          next = (vertex_id + 1) % n_vert;
    double                       angle;
    array<double, 3>             d1, d2;

    d1 = m_vertices[cellVertexIds[next]].getCoords() - m_vertices[cellVertexIds[vertex_id]].getCoords();
    d2 = m_vertices[cellVertexIds[prev]].getCoords() - m_vertices[cellVertexIds[vertex_id]].getCoords();
    d1 = d1/norm2(d1);
    d2 = d2/norm2(d2);
    angle = acos( min(1.0, max(-1.0, dotProduct(d1, d2) ) ) );

    return(angle);
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
double SurfaceKernel::evalMinAngleAtVertex(const long&id, int &vertex_id) const
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
double SurfaceKernel::evalMaxAngleAtVertex(const long&id, int &vertex_id) const
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
double SurfaceKernel::evalAspectRatio(const long &id, int &edge_id) const
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
        M_edge = max(M_edge, l_edge);
    } //next i

    return (M_edge/m_edge);

}

/*!
 * Evaluate facet normal for a cell with specified ID.
 * If cell is of type ElementType::VERTEX or ElementType::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet normal
*/
array<double, 3> SurfaceKernel::evalFacetNormal(const long &id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    array<double, 3>             normal = {{0.0, 0.0, 0.0}};
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
        if (m_spaceDim - getDimension() == 1) {
            std::array<double, 3>       z = {{0.0, 0.0, 1.0}};
            normal = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
            normal = crossProduct(normal, z);
        }
        else return normal;
    }

    if (cell_->getType() == ElementType::TRIANGLE) {
        array<double, 3>                d1, d2;
        d1 = m_vertices[cellVertexIds[1]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        d2 = m_vertices[cellVertexIds[2]].getCoords() - m_vertices[cellVertexIds[0]].getCoords();
        normal = crossProduct(d1, d2);
    }
    else {
        array<double, 3>                d1, d2;
        int                             next, prev, i, nvert = cell_->getVertexCount();
        double                          coeff = 1.0/double(nvert);
        for (i = 0; i < nvert; ++i) {
            next = (i+1) % nvert;
            prev = (nvert + i - 1) % nvert;
            d1 = m_vertices[cellVertexIds[next]].getCoords() - m_vertices[cellVertexIds[i]].getCoords();
            d2 = m_vertices[cellVertexIds[prev]].getCoords() - m_vertices[cellVertexIds[i]].getCoords();
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
array<double, 3> SurfaceKernel::evalEdgeNormal(const long &id, const int &edge_id) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    array<double, 3>                    normal = evalFacetNormal(id);
    const Cell                          *cell_ = &m_cells[id];

    // Counters
    int                                 n_adj = cell_->getAdjacencyCount(edge_id);
    
    // ====================================================================== //
    // COMPUTE EDGE NORMAL                                                    //
    // ====================================================================== //
    if (cell_->getAdjacency(edge_id, 0) != Element::NULL_ID) {
        for (int i = 0; i < n_adj; ++i) {
            normal += evalFacetNormal(cell_->getAdjacency(edge_id, i));
        } //next i
        normal = normal/norm2(normal);
    }
    return(normal);
}

/*!
 * Evaluate the normal unit vector of the specified local vertex.
 *
 * Vertex normal is evaluated as a weighted average of the normals of the
 * 1 ring of the vertex. The weights used in the average are the normalized
 * angles of incident facets at vertex.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \result The normal unit vector of the specified local vertex.
*/
std::array<double, 3> SurfaceKernel::evalVertexNormal(const long &id, const int &vertex) const
{
    return evalLimitedVertexNormal(id, vertex, std::numeric_limits<double>::max());
}

/*!
 * Evaluate the limited normal unit vector of the specified local vertex.
 *
 * Vertex normal is evaluated as a weighted average of the normals of the
 * 1 ring of the vertex. Only the normal whose angle with respect to the
 * considered cell is less that the specified limit are considered. The
 * weights used in the average are the angles of incident facets at vertex.
 *
 * \param[in] id is the cell id
 * \param[in] vertex is the local vertex id
 * \param[in] limit is the maximum allowed misalignment between the normal
 * of the reference cell and the normal of facets used for evaualting the
 * vertex normal
 * \result The normal unit vector of the specified local vertex.
*/
std::array<double, 3> SurfaceKernel::evalLimitedVertexNormal(const long &id, const int &vertex, const double &limit) const
{
    // Get vertex id
    const Cell &cell = m_cells[id];
    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    long vertexId = cellVertexIds[vertex];

    // Cell normal
    std::array<double, 3> cellNormal = evalFacetNormal(id);

    // Cell angle at vertex
    double cellVertexAngle = evalAngleAtVertex(id, vertex);

    // Compute cell vertex neighs
    std::vector<long> vertexNeighs = findCellVertexNeighs(id, vertex);

    // Compute non-normalized normal
    std::array<double, 3> normal = cellVertexAngle * cellNormal;
    for (long facetId : vertexNeighs) {
        // Eval facet normal
        std::array<double, 3> facetNormal = evalFacetNormal(facetId);

        // Discard facets with a misalignment greater than the specified limit
        //
        // The argument of the acos function has to be in the range [-1, 1].
        // Rounding errors may lead to a dot product slightly outside this
        // range. Since the arguments of the dot product are unit vectors,
        // we can safetly clamp the dot product result to be between -1 and
        // 1.
        double misalignment = std::acos(std::min(1.0, std::max(-1.0, dotProduct(facetNormal, cellNormal))));
        if (misalignment > std::abs(limit)) {
            continue;
        }

        // Eval facet angle used as normalization weight
        const Cell &facet = m_cells[facetId];
        double facetVertexAngle = evalAngleAtVertex(facetId, facet.findVertex(vertexId));

        // Add facet contribution to normal
        normal += facetVertexAngle * facetNormal;
    }

    // Normalization
    normal = normal / norm2(normal);

    return normal;
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
    long id = Element::NULL_ID;

#if BITPIT_ENABLE_MPI==1
    if (getRank() == 0) {
        id = getCells().begin()->getId();
    }
#else
    id = getCells().begin()->getId();
#endif

    return adjustCellOrientation(id);
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
bool SurfaceKernel::adjustCellOrientation(const long &seed, const bool &invert)
{
#if BITPIT_ENABLE_MPI==1
    // Check if the communications are needed
    bool ghostExchangeNeeded = isPartitioned();

    // Initialize ghost communications
    std::unordered_set<long> flipped;
    std::unique_ptr<DataCommunicator> ghostComm;
    if (ghostExchangeNeeded) {
        ghostComm = std::unique_ptr<DataCommunicator>(new DataCommunicator(getCommunicator()));

        for (auto &entry : getGhostExchangeSources()) {
            const int rank = entry.first;
            const std::vector<long> sources = entry.second;
            ghostComm->setSend(rank, sources.size() * sizeof(bool));
        }

        for (auto &entry : getGhostExchangeTargets()) {
            const int rank = entry.first;
            const std::vector<long> targets = entry.second;
            ghostComm->setRecv(rank, targets.size() * sizeof(bool));
        }
    }
#endif

    // Initialize the seed
    std::set<long> toVisit;
    if (seed != Element::NULL_ID) {
        assert(getCells().exists(seed));
        toVisit.insert(seed);
        if (invert) {
            flipCellOrientation(seed);
#if BITPIT_ENABLE_MPI==1
            if (ghostExchangeNeeded) {
                flipped.insert(seed);
            }
#endif
        }
    }

    // Orient the surface
    std::unordered_set<long> visited;

    bool orientable = true;
    bool completed = false;
    while (!completed) {
        while (orientable && !toVisit.empty()) {
            auto toVisitBegin = toVisit.begin();
            long cellId = *toVisitBegin;
            toVisit.erase(toVisitBegin);

            visited.insert(cellId);

            const auto &cell = getCell(cellId);
            const long *adjacencyIds = cell.getAdjacencies();
            const long *interfaceIds = cell.getInterfaces();
            const int nInterfaces = cell.getInterfaceCount();

            for (int i = 0; i < nInterfaces; ++i) {
                const long neighId = adjacencyIds[i];

#if BITPIT_ENABLE_MPI==1
                if (ghostExchangeNeeded) {
                    const auto &neigh = getCell(neighId);
                    if (!neigh.isInterior()) {
                        continue;
                    }
                }
#endif

                bool isNeighOriented = sameOrientationAtInterface(interfaceIds[i]);
                if (visited.count(neighId) == 0) {
                    if (!isNeighOriented) {
                        flipCellOrientation(neighId);
#if BITPIT_ENABLE_MPI==1
                        if (ghostExchangeNeeded) {
                            flipped.insert(neighId);
                        }
#endif
                    }

                    toVisit.insert(neighId);
                } else if (!isNeighOriented) {
                    orientable = false;
                    break;
                }
            }
        }

#if BITPIT_ENABLE_MPI==1
        // Comunicate the orientable flag among the partitions
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &orientable, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
        }
#endif

        // If the surface is not orientable we can exit
        if (!orientable) {
            break;
        }

#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            // Add seeds from other partitions
            ghostComm->startAllRecvs();

            for(auto &entry : getGhostExchangeSources()) {
                int rank = entry.first;
                SendBuffer &buffer = ghostComm->getSendBuffer(rank);
                for (long cellId : entry.second) {
                    if (flipped.count(cellId) > 0) {
                        buffer << true;
                    } else {
                        buffer << false;
                    }
                }
                ghostComm->startSend(rank);
            }

            int nPendingRecvs = ghostComm->getRecvCount();
            while (nPendingRecvs != 0) {
                int rank = ghostComm->waitAnyRecv();
                RecvBuffer buffer = ghostComm->getRecvBuffer(rank);

                for (long cellId : getGhostExchangeTargets(rank)) {
                    bool ghostFlipped;
                    buffer >> ghostFlipped;
                    if (ghostFlipped) {
                        toVisit.insert(cellId);
                        flipCellOrientation(cellId);
                    }
                }

                --nPendingRecvs;
            }

            ghostComm->waitAllSends();

            // Clear list of flipped cells
            flipped.clear();
        }
#endif

        // Detect if the orientation is completed
        completed = (toVisit.size() == 0);
#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &completed, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
        }
#endif
    }

    return orientable;
}

/*!
 * Check if the cells sharing an interface have the same orientation.
 *
 * Two cells have the same orientation, if the two cell connectivity vectors
 * cycle through the common interface with opposite order.
 *
 * \param[in] id is the id of the interface
 * \return Returns true if the cells sharing an interface have the same
 * orientation.
 */
bool SurfaceKernel::sameOrientationAtInterface(const long &id)
{
    const auto &face = getInterface(id);
    if (face.isBorder()) {
        return true;
    }

    const long ownerId = face.getOwner();
    const auto &owner = getCell(ownerId);

    const long neighId = face.getNeigh();
    const auto &neigh = getCell(neighId);

    if (face.getType() == ElementType::VERTEX) {
        ConstProxyVector<long> ownerVertexIds = owner.getVertexIds();
        ConstProxyVector<long> neighVertexIds = neigh.getVertexIds();
        if (ownerVertexIds[0] == neighVertexIds[0] || ownerVertexIds[1] == neighVertexIds[1]) {
            return false;
        }
    } else if (face.getType() == ElementType::LINE) {
        const int ownerFace = face.getOwnerFace();
        ConstProxyVector<long> ownerVertexIds = owner.getFaceVertexIds(ownerFace);

        const int neighFace = face.getNeighFace();
        ConstProxyVector<long> neighVertexIds = neigh.getFaceVertexIds(neighFace);

        if (ownerVertexIds[0] != neighVertexIds[1] || ownerVertexIds[1] != neighVertexIds[0]) {
            return false;
        }
    } else {
        throw std::runtime_error ("Type of element not supported in SurfaceKernel");
    }

    return true;
}

/*!
 * Flips the orientation of a cell.
 *
 * \param[in] id is the id of cell that will be flipped
 */
void SurfaceKernel::flipCellOrientation(const long &id)
{
    int top, end;
    auto &cell = getCell(id);

    //
    // In order to flip the orientation of cell,
    // the connectivity will be inverted
    //
    int nCellVertices = cell.getVertexCount();
    const long *cellConnect = cell.getConnect();
    std::unique_ptr<long[]> newCellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
    for (int j = 0; j < nCellVertices; ++j) {
        newCellConnect[j] = cellConnect[j];
    }

    top = 0;
    end = nCellVertices - 1;
    for (int j = 0; j < (int) std::floor(nCellVertices/2); ++j) {
        std::swap(newCellConnect[top], newCellConnect[end]);
        ++top;
        --end;
    }

    //
    // The numbering of the adjacencies and interfaces must be 
    // changed according the new connectivity. This implies that
    // all numbering must be inverted, but the last entry of the
    // adjacencies and interfaces. Since both adjacencies and
    // interfaces are grouped by faces, in a second step each
    // ordering within a face must be inverted
    //
    int faceCount = cell.getFaceCount();
    std::vector<std::vector<long>> newAdjacency, newInterface;
    newAdjacency.resize(faceCount);
    newInterface.resize(faceCount);

    // Copy original ordering in newAdjacency and newInterface
    //
    for (int f = 0; f < faceCount; ++f) {
        // Adjacencies
        std::vector<long> &faceAdjacencies = newAdjacency[f];

        int nCellAdjacencies = cell.getAdjacencyCount(f);
        const long *cellAdjacencies = cell.getAdjacencies(f);

        faceAdjacencies.resize(nCellAdjacencies);
        for (int j = 0; j < nCellAdjacencies; ++j) {
            faceAdjacencies[j] = cellAdjacencies[j];
        }

        // Interfaces
        std::vector<long> &faceInterfaces = newInterface[f];

        int nCellInterfaces = cell.getInterfaceCount(f);
        const long *cellInterfaces = cell.getInterfaces(f);

        faceInterfaces.resize(nCellInterfaces);
        for (int j = 0; j < nCellInterfaces; ++j) {
            faceInterfaces[j] = cellInterfaces[j];
        }
    }

    //
    // Invert entire face based ordering but tha last entry
    //
    top = 0;
    end = faceCount - 2;
    for (int f = 0; f < std::max(1, (int) std::floor((faceCount - 1) / 2)); ++f) {
        std::swap(newAdjacency[top], newAdjacency[end]);
        std::swap(newInterface[top], newInterface[end]);
        ++top;
        --end;
    }

    //
    // Invert all enties within one face
    //
    for (int f = 0; f < faceCount; ++f) {
        // Adjacencies
        std::vector<long> &faceAdjacencies =newAdjacency[f];
        int nFaceAdjacencies = faceAdjacencies.size();

        top = 0;
        end = nFaceAdjacencies - 1;
        for (int j = 0; j<(int) std::floor(nFaceAdjacencies / 2); ++j) {
            std::swap(faceAdjacencies[top], faceAdjacencies[end]);
            ++top;
            --end;
        }

        // Interfaces
        std::vector<long> &faceInterfaces =newInterface[f];
        int nFaceInterfaces = faceInterfaces.size();

        top = 0;
        end = nFaceInterfaces - 1;
        for (int j = 0; j<(int) std::floor(nFaceInterfaces / 2); ++j) {
            std::swap(faceInterfaces[top], faceInterfaces[end]);
            ++top;
            --end;
        }
    }

    cell.setConnect(std::move(newCellConnect));

    if (newAdjacency[0].size() != 0) {
        cell.setAdjacencies(newAdjacency);
    }

    if (newInterface[0].size() != 0) {
        cell.setInterfaces(newInterface);

        //
        // For the interfaces, we need to correct either
        // ownerFace or neighFace, respectivly.
        //
        const long* intPtr = cell.getInterfaces(); 
        int interfaceCount = cell.getInterfaceCount();
        for (int i=0; i<interfaceCount; ++i) {
            long interfaceId = intPtr[i];
            if (interfaceId < 0) {
                continue;
            }

            Interface &interface=getInterface(interfaceId);
            long owner = interface.getOwner();
            if (owner == id) {
                int ownerFace = interface.getOwnerFace();
                if (ownerFace != faceCount - 1) {
                    ownerFace = faceCount - 2 - ownerFace;
                    interface.setOwner(id,ownerFace);
                }
            } else {
                int neighFace = interface.getNeighFace();
                if (neighFace != faceCount - 1) {
                    neighFace = faceCount - 2 - neighFace;
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
void SurfaceKernel::displayQualityStats(ostream& out, unsigned int padding) const
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    string              indent(padding, ' ');

    // Counters
    // none

    // ====================================================================== //
    // DISPLAY STATS                                                          //
    // ====================================================================== //

    // Distribution for aspect ratio ---------------------------------------- //
    out << indent << "Aspect ratio distribution ---------------" << endl;
    {
        // Scope variables
        long                            count;
        vector<double>                  bins, hist;

        // Compute histogram
        hist = computeHistogram(&SurfaceKernel::evalAspectRatio, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "AR", out, padding + 2);
    }

    // Distribution of min. angle ------------------------------------------- //
    out << indent << "Min. angle distribution -----------------" << endl;
    {
        // Scope variables
        long                            count;
        vector<double>                  bins, hist;

        // Compute histogram
        hist = computeHistogram(&SurfaceKernel::evalMinAngleAtVertex, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "min. angle", out, padding + 2);
    }

    // Distribution of min. angle ------------------------------------------- //
    out << indent << "Max. angle distribution -----------------" << endl;
    {
        // Scope variables
        long                            count;
        vector<double>                  bins, hist;

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
    vector<double>              hist;

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
    const long                  &count,
    const vector<double>        &bins,
    const vector<double>        &hist,
    const string                &stats_name,
    ostream                     &out,
    unsigned int                 padding
) const {
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                 nbins = bins.size();
    size_t              n_fill;
    string              indent(padding, ' ');
    stringstream        ss;

    // Counters
    int                 i;

    // ====================================================================== //
    // DISPLAY HISTOGRAM                                                      //
    // ====================================================================== //

    // Set stream properties
    ss << setprecision(3);

    // Display histograms
    out << indent << "poll size: " << count << endl;
    out << indent;
    ss << hist[0];
    n_fill = ss.str().length();
    ss << string(max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << "%, with          " << stats_name << " < ";
    ss.str("");
    ss << bins[0];
    out << ss.str() << endl;
    ss.str("");
    for (i = 0; i < nbins-1; ++i) {
        out << indent;
        ss << hist[i+1];
        n_fill = ss.str().length();
        ss << string(max(size_t(0),  6 - n_fill), ' ');
        out << ss.str() << "%, with ";
        ss.str("");
        ss << bins[i];
        n_fill = ss.str().length();
        ss << string(max(size_t(0), 6 - n_fill), ' ');
        out << ss.str() << " < " << stats_name << " < ";
        ss.str("");
        ss << bins[i+1];
        out << ss.str() << endl;
        ss.str("");
    } //next i
    out << indent;
    ss << hist[i+1];
    n_fill = ss.str().length();
    ss << string(max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << "%, with ";
    ss.str("");
    ss << bins[i];
    n_fill = ss.str().length();
    ss << string(max(size_t(0), 6 - n_fill), ' ');
    out << ss.str() << " < " << stats_name << endl;
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
bool SurfaceKernel::compareSelectedTypes(const unsigned short &mask_, const ElementType &type_) const
{
    unsigned short       masked = m_selectionTypes.at(type_);
    return ( (mask_ & masked) == masked );
}

}
