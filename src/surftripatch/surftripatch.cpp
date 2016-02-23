/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include "bitpit_common.hpp"

#include "surftripatch.hpp"

namespace bitpit {

const std::map<ElementInfo::Type, unsigned short>  SurfTriPatch::m_selectionTypes({{ElementInfo::QUAD, SurfTriPatch::SELECT_QUAD},
                                {ElementInfo::TRIANGLE, SurfTriPatch::SELECT_TRIANGLE}});
const unsigned short SurfTriPatch::SELECT_TRIANGLE = 1;
const unsigned short SurfTriPatch::SELECT_QUAD     = 2;
const unsigned short SurfTriPatch::SELECT_ALL      = 3;

/*!
	\ingroup surftripatch
	@{
*/

/*!
	\class SurfTriPatch

	\brief The SurfTriPatch class defines an unstructured surface
	triangulation.

	SurfTriPatch defines an unstructured surface triangulation.
*/

/*!
	Creates a new patch.

	\param id is the id of the patch
*/
SurfTriPatch::SurfTriPatch(const int &id)
	: SurfaceKernel(id, 2, true)
{

}

/*!
 * Enables or disables expert mode.
 *
 * When expert mode is enabled, it will be possible to change the
 * patch using low level functions (e.g., it will be possible to
 * add individual cells, add vertices, delete cells, ...).
 *
 * \param expert if true, the expert mode will be enabled
 */
void SurfTriPatch::setExpert(bool expert)
{
	SurfaceKernel::setExpert(expert);
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double SurfTriPatch::evalCellVolume(const long &id)
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double SurfTriPatch::evalCellSize(const long &id)
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
    return(sqrt(evalFacetArea(id))); 
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double SurfTriPatch::evalInterfaceArea(const long &id)
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> SurfTriPatch::evalInterfaceNormal(const long &id)
{
	BITPIT_UNUSED(id);

	return {{0., 0., 0.}};
}

/*!
	Updates the patch.

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> SurfTriPatch::_update(bool trackAdaption)
{
	if (!isDirty()) {
		return std::vector<Adaption::Info>();
	}

	std::cout << ">> Updating surface triangulation mesh\n";

	// Adaption info
	std::vector<Adaption::Info> adaptionData;
	if (trackAdaption) {

	}

	// Done
	return adaptionData;
}

/*!
	Marks a cell for refinement.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id is the id of the cell that needs to be refined
*/
bool SurfTriPatch::_markCellForRefinement(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Marks a cell for coarsening.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id the cell to be refined
*/
bool SurfTriPatch::_markCellForCoarsening(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Enables cell balancing.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
bool SurfTriPatch::_enableCellBalancing(const long &id, bool enabled)
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(enabled);

	return false;
}

/*!
 * Checks if the specified point is inside the patch.
 *
 * \param[in] point is the point to be checked
 * \result Returns true if the point is inside the patch, false otherwise.
 */
bool SurfTriPatch::isPointInside(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

/*!
 * Locates the cell the contains the point.
 *
 * If the point is not inside the patch, the function returns the id of the
 * null element.
 *
 * \param[in] point is the point to be checked
 * \result Returns the linear id of the cell the contains the point. If the
 * point is not inside the patch, the function returns the id of the null
 * element.
 */
long SurfTriPatch::locatePoint(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

/*!
 * Fill adjacencies info for each cell.
*/
void SurfTriPatch::buildAdjacencies(void)
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                                         n_faces;
    unordered_map<long, vector<long> >          V2S;

    // Counters
    int                                         i, j, k, m;
    long int                                    vertex_idx;
    long int                                    simplex_idx;
    PiercedVector< Cell >::iterator             c_, e_ = cellEnd();
    Cell                                        *n_;

    // ====================================================================== //
    // RESET ADJACENCY INFO                                                   //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        vector<int>                                     interfs;

        // Allocate memory for adjacencies ---------------------------------- //
        for (c_ = cellBegin(); c_ != e_; ++c_) {
            c_->resetAdjacencies();
        } //next c_
    }

    // ====================================================================== //
    // BUILD VERTEX->SIMPLEX CONNECTIVITY                                     //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        unordered_map< long, vector< long > >::iterator it_;

        // Build vertex->simplex connectivity ------------------------------- //
        for (c_ = cellBegin(); c_ != e_; ++c_) {
            n_faces = c_->getFaceCount();
            for (i = 0; i < n_faces; ++i) {
                vertex_idx = c_->getVertex(i);
                simplex_idx = c_->get_id();
                V2S[vertex_idx].push_back(simplex_idx);
            } //next i
        } //next c_

    }

    // ====================================================================== //
    // BUILD ADJACENCY                                                        //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        bool                    check;
        int                     n_candidates;
        int                     m_faces;
        int                     n_vertex_on_face;
        long                    candidate_idx;
        vector< int >           face_loc_connectivity;
        vector< long >          face_connectivity;
        vector< long >          candidate_list;

        // Build adjacencies ------------------------------------------------ //
        for (c_ = cellBegin(); c_ != e_; ++c_) {

            simplex_idx = c_->get_id();
            n_faces = c_->getFaceCount();

            for (m = 0; m < n_faces; m++) {

                // Build face connectivity
                face_loc_connectivity = c_->getFaceLocalConnect(m);
                n_vertex_on_face = face_loc_connectivity.size();
                face_connectivity.resize(n_vertex_on_face, -1);
                for (i = 0; i < n_vertex_on_face; ++i) {
                    face_connectivity[i] = c_->getVertex(face_loc_connectivity[i]);
                } //next i

                // Build list of candidates for adjacency test
                vertex_idx = face_connectivity[0];
                candidate_list = V2S[vertex_idx];
                j = 1;
                n_candidates = candidate_list.size();
                while ( ( n_candidates > 0 ) && ( j < n_vertex_on_face ) ) {
                    vertex_idx = face_connectivity[j];
                    candidate_list = utils::intersectionVector(candidate_list, V2S[vertex_idx]);
                    j++;
                    n_candidates = candidate_list.size();
                } //next j
                utils::eraseValue(candidate_list, simplex_idx);
                --n_candidates;

                // Update adjacencies
                for (j = 0; j < n_candidates; ++j) {
                    candidate_idx = candidate_list[j];
                    if ( (candidate_idx > simplex_idx) ) {
                        n_ = &(m_cells[candidate_idx]);
                        m_faces = n_->getFaceCount();
                        k = 0;
                        check = false;
                        while (!check && (k < m_faces)) {
                            check = isSameFace(simplex_idx, m, candidate_idx, k);
                            ++k;
                        } //next
                        --k;
                        if (check) {
                            c_->pushAdjacency(m, candidate_idx);
                            n_->pushAdjacency(k, simplex_idx);
                        }
                    }
                } //next j

            } //next m
        } //next c_

    }

    return;
}

/*!
 * Update adjacencies info for cells with specified ID.
 * 
 * \param[in] cell_ids list of cell ids
*/
void SurfTriPatch::updateAdjacencies(const std::vector<long> &cell_ids)
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                                         n_faces;
    unordered_map<long, vector<long> >          V2S;

    // Counters
    int                                         i, j, k, m;
    long int                                    vertex_idx;
    long int                                    simplex_idx;
    std::vector<long>::const_iterator           i_, e_ = cell_ids.cend();
    Cell                                        *c_, *n_;

    // ====================================================================== //
    // RESET ADJACENCY INFO                                                   //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        vector<int>                                     interfs;

        // Allocate memory for adjacencies ---------------------------------- //
        for (i_ = cell_ids.cbegin(); i_ != e_; ++i_) {
            m_cells[*i_].resetAdjacencies();
        } //next i_
    }

    // ====================================================================== //
    // BUILD VERTEX->SIMPLEX CONNECTIVITY                                     //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        CellIterator                    j_, k_ = cellEnd();

        // Build vertex->simplex connectivity ------------------------------- //
        for (j_ = cellBegin(); j_ != k_; ++j_) {
            n_faces = j_->getFaceCount();
            for (i = 0; i < n_faces; ++i) {
                vertex_idx = j_->getVertex(i);
                simplex_idx = j_->get_id();
                V2S[vertex_idx].push_back(simplex_idx);
            } //next i
        } //next c_

    }

    // ====================================================================== //
    // UPDATE ADJACENCY                                                       //
    // ====================================================================== //
    {
        // Scope variables -------------------------------------------------- //
        bool                    check;
        int                     n_candidates;
        int                     m_faces;
        int                     n_vertex_on_face;
        long                    candidate_idx;
        vector< int >           face_loc_connectivity;
        vector< long >          face_connectivity;
        vector< long >          candidate_list;

        // Build adjacencies ------------------------------------------------ //
        for (i_ = cell_ids.cbegin(); i_ != e_; ++i_) {

            c_ = &m_cells[*i_];
            simplex_idx = c_->get_id();
            n_faces = c_->getFaceCount();

            for (m = 0; m < n_faces; m++) {

                // Build face connectivity
                face_loc_connectivity = c_->getFaceLocalConnect(m);
                n_vertex_on_face = face_loc_connectivity.size();
                face_connectivity.resize(n_vertex_on_face, -1);
                for (i = 0; i < n_vertex_on_face; ++i) {
                    face_connectivity[i] = c_->getVertex(face_loc_connectivity[i]);
                } //next i

                // Build list of candidates for adjacency test
                vertex_idx = face_connectivity[0];
                candidate_list = V2S[vertex_idx];
                j = 1;
                n_candidates = candidate_list.size();
                while ( ( n_candidates > 0 ) && ( j < n_vertex_on_face ) ) {
                    vertex_idx = face_connectivity[j];
                    candidate_list = utils::intersectionVector(candidate_list, V2S[vertex_idx]);
                    j++;
                    n_candidates = candidate_list.size();
                } //next j
                utils::eraseValue(candidate_list, simplex_idx);
                --n_candidates;

                // Update adjacencies
                for (j = 0; j < n_candidates; ++j) {
                    candidate_idx = candidate_list[j];
//                     if ( (candidate_idx > simplex_idx) ) {
                        n_ = &(m_cells[candidate_idx]);
                        m_faces = n_->getFaceCount();
                        k = 0;
                        check = false;
                        while (!check && (k < m_faces)) {
                            check = isSameFace(simplex_idx, m, candidate_idx, k);
                            ++k;
                        } //next
                        --k;
                        if (check) {
                            c_->pushAdjacency(m, candidate_idx);
                            n_->pushAdjacency(k, simplex_idx);
                        }
//                     }
                } //next j

            } //next m
        } //next c_

    }

    return;
}

/*!
 *  Evaluate the length of the edge with specified local index
 *  for e cell with specified ID.
 *  If the cell is of type ElementInfo::VERTEX or ElementInfo::LINE
 *  returns 0.0.
 * 
 *  \param[in] id cell id
 *  \param[in] edge_id edge local index
 * 
 *  \result minimal edge length
*/
double SurfTriPatch::evalEdgeLength(const long &id, const int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 nv_on_edge = 2;

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::LINE)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::UNDEFINED)) return 0.0;

    double edge_length = 0.0;
    vector<int> face_loc_connect(2, Vertex::NULL_ID);
    face_loc_connect = cell_->getFaceLocalConnect(edge_id);
    face_loc_connect[0] = cell_->getVertex(face_loc_connect[0]);
    face_loc_connect[1] = cell_->getVertex(face_loc_connect[1]);
    edge_length = norm2(m_vertices[face_loc_connect[0]].getCoords() - m_vertices[face_loc_connect[1]].getCoords());

    return(edge_length);
}

/*!
 *  Evaluate the minimal edge length for e cell with specified ID.
 *  If the cell is of type ElementInfo::VERTEX or ElementInfo::LINE
 *  returns 0.0.
 * 
 *  \param[in] id cell id
 *  \param[in,out] edge_id on output stores the local index of edge with minimal edge length
 * 
 *  \result minimal edge length
*/
double SurfTriPatch::evalMinEdgeLength(const long &id, int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 i, j;
    int                                 n_faces;

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::LINE)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::UNDEFINED)) return 0.0;

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
 *  If the cell is of type ElementInfo::VERTEX or ElementInfo::LINE
 *  returns 0.0.
 * 
 *  \param[in] id cell id
 *  \param[in,out] edge_id on output stores the local inde xof edge with maximal edge length
 * 
 *  \result maximal edge length
*/
double SurfTriPatch::evalMaxEdgeLength(const long &id, int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 i, j;
    int                                 n_faces;

    // ====================================================================== //
    // COMPUTE MIN EDGE SIZE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::LINE)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::UNDEFINED)) return 0.0;

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
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, a returns zero.
 * 
 * \param[in] id cell global ID
 * \param[in] vertex_id vertex local ID
*/
double SurfTriPatch::evalAngleAtVertex(const long &id, const int &vertex_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                        *cell_ = &m_cells[id];

    // Counters

    // ====================================================================== //
    // EVALUATE ANGLE AT SPECIFIED VERTEX                                     //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return 0.0;

    int                          n_vert = cell_->getVertexCount();
    int                          prev = (n_vert + vertex_id - 1) % n_vert;
    int                          next = (vertex_id + 1) % n_vert;
    double                       angle;
    array<double, 3>             d1, d2;

    d1 = m_vertices[cell_->getVertex(next)].getCoords() - m_vertices[cell_->getVertex(vertex_id)].getCoords();
    d2 = m_vertices[cell_->getVertex(prev)].getCoords() - m_vertices[cell_->getVertex(vertex_id)].getCoords();
    d1 = d1/norm2(d1);
    d2 = d2/norm2(d2);
    angle = acos( min(1.0, max(-1.0, dotProduct(d1, d2) ) ) );

    return(angle);
}

/*!
 * Evaluate the minimal angle at vertex for a cell with specified ID.
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, a returns zero.
 *
 * \param[in] id cell id
 * \param[in,out] vertex_id on output stores the local index of vertex with minimal angle
 * 
 * \result minimal angle at vertex
*/
double SurfTriPatch::evalMinAngleAtVertex(const long&id, int &vertex_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell               *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE MINIMAL ANGLE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return 0.0;

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
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, a returns zero.
 *
 * \param[in] id cell id
 * \param[in,out] vertex_id on output stores the local index of vertex with minimal angle
 * 
 * \result maximal angle at vertex
*/
double SurfTriPatch::evalMaxAngleAtVertex(const long&id, int &vertex_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell               *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE MINIMAL ANGLE                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return 0.0;

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
 * Evaluate facet normal for a cell with specified ID.
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet normal
*/
array<double, 3> SurfTriPatch::evalFacetNormal(const long &id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    array<double, 3>             normal{0.0, 0.0, 0.0};
    Cell                        *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE NORMAL                                                         //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return normal;

    if (cell_->getType() == ElementInfo::TRIANGLE) {
        array<double, 3>                d1, d2;
        d1 = m_vertices[cell_->getVertex(1)].getCoords() - m_vertices[cell_->getVertex(0)].getCoords();
        d2 = m_vertices[cell_->getVertex(2)].getCoords() - m_vertices[cell_->getVertex(0)].getCoords();
        normal = crossProduct(d1, d2);
    }
    else {
        array<double, 3>                d1, d2;
        int                             next, prev, i, nvert = cell_->getVertexCount();
        double                          coeff = 1.0/double(nvert);
        for (i = 0; i < nvert; ++i) {
            next = (i+1) % nvert;
            prev = (nvert + i - 1) % nvert;
            d1 = m_vertices[cell_->getVertex(next)].getCoords() - m_vertices[cell_->getVertex(i)].getCoords();
            d2 = m_vertices[cell_->getVertex(prev)].getCoords() - m_vertices[cell_->getVertex(i)].getCoords();
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
array<double, 3> SurfTriPatch::evalEdgeNormal(const long &id, const int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    array<double, 3>                    normal = evalFacetNormal(id);
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 n_adj = cell_->getAdjacencyCount(edge_id);
    
    // ====================================================================== //
    // COMPUTE EDGE NORMAL                                                    //
    // ====================================================================== //
    if (cell_->getAdjacency(edge_id) != Element::NULL_ID) {
        for (int i = 0; i < n_adj; ++i) {
            normal += evalFacetNormal(cell_->getAdjacency(edge_id, i));
        } //next i
        normal = normal/norm2(normal);
    }
    return(normal);
}

/*!
 * Evaluate the normal unit vector at the vertex with specified local ID belonging to
 * the surface facet with specified global id. Vertex normal are computed as the arithmeic
 * mean of the normals of the edges incident to the vertex. In case adjacencies are
 * not build the vertex normal will be the same as the facet normal.
 * 
 * \param[in] id cell global ID
 * \param[in] vert_id vertex local ID on the specified cell
*/
array<double, 3> SurfTriPatch::evalVertexNormal(const long &id, const int &vert_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];
    long                                vert_ID = cell_->getVertex(vert_id);
    array<double, 3>                    normal;
    vector<long>                        ring1;

    // Counters
    vector<long>::const_iterator        i_, e_;
    
    // ====================================================================== //
    // COMPUTE VERTEX NORMAL                                                  //
    // ====================================================================== //

    // Compute 1-ring of vertex
    ring1 = findCellVertexOneRing(id, vert_id);

    // Compute angles of incident facet at vertex
    int i, n_ring1 = ring1.size();
    double disk_angle = 0.0;
    vector<double> angles(ring1.size(), 0.);
    for (i = 0; i < n_ring1; ++i) {
        cell_ = &m_cells[ring1[i]];
        angles[i] = evalAngleAtVertex(cell_->get_id(), cell_->findVertex(vert_ID));
    } //next i_
    sum(angles, disk_angle);
    angles = angles/disk_angle;
    e_ = ring1.cend();
    normal.fill(0.0);
    i = 0;
    for (i_ = ring1.cbegin(); i_ != e_; ++i_) {
        normal += angles[i]*evalFacetNormal(*i_);
        ++i;
    } //next i_

    // Average of incident element
//     e_ = ring1.cend();
//     normal.fill(0.0);
//     for (i_ = ring1.cbegin(); i_ != e_; ++i_) {
//         normal += evalFacetNormal(*i_);
//     } //next i_

    // Average of edge normal
//     int                                 loc_id, n_vert;
//     e_ = ring1.cend();
//     normal.fill(0.0);
//     for (i_ = ring1.cbegin(); i_ != e_; ++i_) {
//         cell_ = &m_cells[*i_];
//         n_vert = cell_->getVertexCount();
//         loc_id = cell_->findVertex(vert_ID);
//         normal += evalEdgeNormal(*i_, loc_id)/double(cell_->getAdjacencyCount(loc_id)
//                                                   + cell_->getAdjacency(loc_id) != Element::NULL_ID);
//         loc_id = (n_vert + loc_id - 1) % n_vert;
//         normal += evalEdgeNormal(*i_, loc_id)/double(cell_->getAdjacencyCount(loc_id)
//                                                   + cell_->getAdjacency(loc_id) != Element::NULL_ID);
//     } //next i_

    // Normalization
    normal = normal/norm2(normal);

    return(normal);
}

/*!
 * Evalute the aspect ratio for a cell with specified ID. The aspect ratio
 * is defined as the ratio between the longest and the shortest edge.
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * \param[in,out] edge_id on output stores the index of the shortest edge
 * 
 * \result cell aspect ratio
*/
double SurfTriPatch::evalAspectRatio(const long &id, int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                        *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // EVALUATE ASPECT RATIO                                                  //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return 0.0;

    double                      l_edge;
    double                      m_edge = std::numeric_limits<double>::max();
    double                      M_edge = std::numeric_limits<double>::min();
    double                      ar = 1.0;
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
 * Compute cell center coordinates. Cell center coordinates are computed
 * as the arithmetic average of cell's vertex coordinates.
 * 
 * \param[in] id cell id
 * 
 * \result returns cell's center coordinates
*/
std::array<double, 3> SurfTriPatch::evalCellCentroid(const long &id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];
    array<double, 3>                    C;

    // Counters
    int                                 i, n = cell_->getVertexCount();

    // ====================================================================== //
    // COMPUTE CELL CENTER COORDINATES                                        //
    // ====================================================================== //
    C.fill(.0);
    for (i = 0; i < n; ++i) {
        C += m_vertices[cell_->getVertex(i)].getCoords();
    } //next i
    C = C/double(n);

    return(C);
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
void SurfTriPatch::displayHistogram(
    const long                  &count,
    const vector<double>        &bins,
    const vector<double>        &hist,
    const string                &stats_name,
    ostream                     &out,
    unsigned int                 padding
) {
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
 * Display quality stats for the mesh currently stored.
 * 
 * \param[in,out] out output stream where stats will be printed out
 * \param[in] padding (default = 0) number of trailing spaces
*/
void SurfTriPatch::displayQualityStats(
    ostream                     &out,
    unsigned int                 padding
) {
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
        hist = computeHistogram(&SurfTriPatch::evalAspectRatio, bins, count, 8,
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
        hist = computeHistogram(&SurfTriPatch::evalMinAngleAtVertex, bins, count, 8,
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
        hist = computeHistogram(&SurfTriPatch::evalMaxAngleAtVertex, bins, count, 8,
                                 SELECT_TRIANGLE | SELECT_QUAD);

        // Display histogram
        displayHistogram(count, bins, hist, "max. angle", out, padding + 2);
    }
    
    return;
}

/*!
 * Evaluate facet area for a cell with specified ID. If cell is of type
 * ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet area
*/
double SurfTriPatch::evalFacetArea(const long &id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                        *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // EVALUATE FACET AREA                                                    //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)
     || (cell_->getType() == ElementInfo::LINE)) return 0.0;

    array<double, 3>            d1, d2;
    if (cell_->getType() == ElementInfo::TRIANGLE) {
        d1 = m_vertices[cell_->getVertex(1)].getCoords() - m_vertices[cell_->getVertex(0)].getCoords();
        d2 = m_vertices[cell_->getVertex(2)].getCoords() - m_vertices[cell_->getVertex(0)].getCoords();
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
            d1 = m_vertices[cell_->getVertex(next)].getCoords() - m_vertices[cell_->getVertex(i)].getCoords();
            d2 = m_vertices[cell_->getVertex(prev)].getCoords() - m_vertices[cell_->getVertex(i)].getCoords();
            area += coeff*norm2(crossProduct(d1, d2));
        } //next i
        return(area);
    }
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
 * \param[in] n_internvals (default = 8), number of intervals to be used for histogram construction.
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
vector<double> SurfTriPatch::computeHistogram(
    eval_f_                      funct_,
    vector<double>              &bins,
    long                        &count,
    int                          n_intervals,
    unsigned short               mask
) {
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
        id = cell_.get_id();
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

//TODO: Aggiungere un metodo in SurfTriPatch per aggiungere più vertici.
/*!
 * Extract the edge network from surface mesh. If adjacencies are not built
 * edges shared by more than 1 element are counted twice. Edges are appended
 * to the content of the input SurfTriPatch
 * 
 * \param[in,out] net on output stores the network of edges
*/
void SurfTriPatch::extractEdgeNetwork(SurfTriPatch &net)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    bool                                        check;
    int                                         n_faces, n_adj, n_vert;
    long                                        id;
    vector<int>                                 face_loc_connect;
    vector<long>                                face_connect;

    // Counters
    int                                         i, j;
    vector<int>::const_iterator                 i_;
    vector<long>::iterator                      j_;
    VertexIterator                              v_, ve_ = vertexEnd();
    CellIterator                                c_, ce_ = cellEnd();

    // ====================================================================== //
    // INITIALIZE DATA STRUCTURE                                              //
    // ====================================================================== //
    net.reserveCells(net.getCellCount() + countFaces());
    net.reserveVertices(net.getVertexCount() + m_nVertices);

    // ====================================================================== //
    // ADD VERTEX TO net                                                      //
    // ====================================================================== //
    for (v_ = vertexBegin(); v_ != ve_; ++v_) {
        net.addVertex(v_->getCoords(), v_->get_id());
    } //next v_

    // ====================================================================== //
    // ADD EDGES                                                              //
    // ====================================================================== //
    for (c_ = cellBegin(); c_ != ce_; ++c_) {
        id = c_->get_id();
        n_faces = c_->getFaceCount();
        for (i = 0; i < n_faces; ++i) {
            check = true;
            n_adj = c_->getAdjacencyCount(i);
            for (j = 0; j < n_adj; ++j) {
                check = check && (id > c_->getAdjacency(i, j));
            } //next j
            if (check) {
                face_loc_connect = c_->getFaceLocalConnect(i);
                n_vert = face_loc_connect.size();
                face_connect.resize(n_vert);
                j_ = face_connect.begin();
                for (i_ = face_loc_connect.cbegin(); i_ != face_loc_connect.cend(); ++i_) {
                    *j_ = c_->getVertex(*i_);
                    ++j_;
                } //next i_
                net.addCell(c_->getFaceType(i), true, face_connect);
            }
        } //next i
    } //next c_

    return;
}

/*!
 * Given and element type and the selection mask, returns true if the element type
 * is accounted for in the selection mask.
 * 
 * \param[in] mask_ selection mask_
 * \param[in] type_ element type
*/
bool SurfTriPatch::compareSelectedTypes(const unsigned short &mask_, const ElementInfo::Type &type_)
{
    unsigned short       masked = m_selectionTypes.at(type_);
    return ( (mask_ & masked) == masked );
}

//TODO: normals??
//TODO: error flag on output
//TODO: import a specified solid (ascii format only)
/*!
 * Import surface tasselation from S.T.L. file. STL facet are added at to the
 * present mesh, i.e. current mesh content is not discarded. Howver no checks
 * are performed to ensure that no duplicated vertices or cells are created.
 * 
 * \param[in] stl_name name of stl file
 * \param[in] isBinary flag for binary (true), of ASCII (false) stl file
 * 
 * \result on output returns an error flag for I/O error
*/
unsigned short SurfTriPatch::importSTL(const string &stl_name, const bool &isBinary)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Parameters
    const unordered_map<size_t, ElementInfo::Type> ele_type{
        {0, ElementInfo::UNDEFINED},
        {1, ElementInfo::VERTEX},
        {2, ElementInfo::LINE},
        {3, ElementInfo::TRIANGLE},
        {4, ElementInfo::QUAD}
    };

    // Local variables
    int                                         n_v;
    int                                         nVertex = 0;
    int                                         nSimplex = 0;
    vector<array<double, 3>>                    vertexList;
    vector<array<double, 3>>                    normalList;
    vector<vector<int>>                         connectivityList;
    vector<long>                                connect;
    unordered_map<long, long>                   vertexMap;
    STLObj                                      STL(stl_name, isBinary);
    

    // Counters
    int                                         i;
    long                                        v_counter = 0, c_counter = 0;
    VertexIterator                              i_;
    CellIterator                                j_;
    vector<array<double, 3>>::const_iterator    v_, ve_;
    vector<vector<int>>::const_iterator         c_, ce_;
    vector<int>::const_iterator                 w_, we_;

    // ====================================================================== //
    // READ MESH DATA FROM THE STL FILE                                       //
    // ====================================================================== //
    STL.load(nVertex, nSimplex, vertexList, normalList, connectivityList);

    // ====================================================================== //
    // PREPARE MESH FOR DATA IMPORT                                           //
    // ====================================================================== //
    reserveVertices(m_nVertices + nVertex);
    reserveCells(m_nInternals + m_nGhosts + nSimplex);

    // ====================================================================== //
    // ADD VERTICES TO MESH                                                   //
    // ====================================================================== //
    ve_ = vertexList.cend();
    for (v_ = vertexList.cbegin(); v_ != ve_; ++v_) {
        i_ = addVertex(*v_);
        vertexMap[v_counter] = i_->get_id();
        ++v_counter;
    } //next v_

    // ====================================================================== //
    // ADD CELLS TO MESH                                                      //
    // ====================================================================== //
    ce_ = connectivityList.cend();
    for (c_ = connectivityList.cbegin(); c_ != ce_; ++c_) {

        // Remap STL connectivity
        n_v = c_->size();
        connect.resize(n_v, Vertex::NULL_ID);
        we_ = c_->cend();
        i = 0;
        for (w_ = c_->cbegin(); w_ < we_; ++w_) {
            connect[i] = vertexMap[*w_];
            ++i;
        } //next w_

        // Add cell
        addCell(ele_type.at(n_v), true, connect);
    } //next c_

    return 0;
}

//TODO: normals??
//TODO: error flag on output
//TODO: conversion of quad into tria
/*!
 * Export surface tasselation in a STL format. No check is perfomed on element type
 * therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation.
 * 
 * \param[in] stl_name name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 * \param[in] exportInternalsOnly flag for exporting only internal cells (true),
 * or internal+ghost cells (false).
 * 
 * \result on output returns an error flag for I/O error.
*/
unsigned short SurfTriPatch::exportSTL(const string &stl_name, const bool &isBinary, bool exportInternalsOnly)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                                         n_v;
    int                                         nVertex;
    int                                         nSimplex;
    vector<array<double, 3>>                    vertexList;
    vector<array<double, 3>>                    normalList;
    vector<vector<int>>                         connectivityList;
    unordered_map<long, long>                   vertexMap;

    // Counters
    int                                         v_count ,j;
    vector<array<double, 3>>::iterator          i_;
    vector<vector<int>>::iterator               j_;
    vector<int>::iterator                       k_, ke_;
    VertexIterator                              v_, ve_;
    CellIterator                                c_, cb_, ce_;

    // ====================================================================== //
    // INITIALIZE DATA STRUCTURE                                              //
    // ====================================================================== //
    nSimplex = m_nInternals;
    if (!exportInternalsOnly) nSimplex += m_nGhosts;
    vertexList.resize(m_nVertices);
    normalList.resize(nSimplex);
    connectivityList.resize(nSimplex, vector<int>(3, 0));

    // ====================================================================== //
    // CREATE VERTEX LIST                                                     //
    // ====================================================================== //
    i_ = vertexList.begin();
    ve_ = vertexEnd();
    v_count = 0;
    for (v_ = vertexBegin(); v_ != ve_; ++v_) {

        // Store vertex coordinates
        *i_ = v_->getCoords();
        vertexMap[v_->get_id()] = v_count;

        // Update counters
        ++v_count;
        ++i_;

    } //next v_
    nVertex = m_nVertices;

    // ====================================================================== //
    // CREATE CONNECTIVITY                                                    //
    // ====================================================================== //
    if (exportInternalsOnly) {
        cb_ = internalBegin();
        ce_ = internalEnd();
    }
    else {
        cb_ = cellBegin();
        ce_ = cellEnd();
    }
    i_ = normalList.begin();
    j_ = connectivityList.begin();
    for (c_ = cb_; c_ != ce_; ++c_) {

        // Build normals
        *i_ = std::move(evalFacetNormal(c_->get_id()));
        
        // Build connectivity
        n_v = min(3, c_->getVertexCount());
        ke_ = j_->end();
        j = 0;
        for (k_ = j_->begin(); k_ != ke_; ++k_) {
            *k_ = vertexMap[c_->getVertex(j)];
            ++j;
        } //next k_

        // Update counters
        ++j_;
        ++i_;
    } //next c_

    // ====================================================================== //
    // EXPORT STL DATA                                                        //
    // ====================================================================== //
    STLObj                                      STL(stl_name, isBinary);
    STL.save("", nVertex, nSimplex, vertexList, normalList, connectivityList);

    return 0;
}

/*!
	@}
*/

}
