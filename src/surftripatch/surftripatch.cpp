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
	: Patch(id, 2, true)
{

}

/*!
	Destroys the patch.
*/
SurfTriPatch::~SurfTriPatch()
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
	Patch::setExpert(expert);
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
	BITPIT_UNUSED(id);

	return 0;
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
    Evaluate the length of the edge with specified local index
    for e cell with specified ID.
    If the cell is of type ElementInfo::SEGMENT or ElementInfo::POINT
    returns 0.0.

    \param[in] id cell id
    \param[in] edge_id edge local index

    \result minimal edge length
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
    Evaluate the minimal edge length for e cell with specified ID.
    If the cell is of type ElementInfo::SEGMENT or ElementInfo::POINT
    returns 0.0.

    \param[in] id cell id

    \result minimal edge length
*/
double SurfTriPatch::evalMinEdgeLength(const long &id)
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

    double edge_length = std::numeric_limits<double>::max();
    vector<int> face_loc_connect(2, Vertex::NULL_ID);
    n_faces = cell_->getFaceCount();
    for (i = 0; i < n_faces; ++i) {
        face_loc_connect = cell_->getFaceLocalConnect(i);
        face_loc_connect[0] = cell_->getVertex(face_loc_connect[0]);
        face_loc_connect[1] = cell_->getVertex(face_loc_connect[1]);
        edge_length = min( edge_length, norm2(m_vertices[face_loc_connect[0]].getCoords() - m_vertices[face_loc_connect[1]].getCoords()) );
    } //next i

    return(edge_length);
}

/*!
    Evaluate the maximal edge length for e cell with specified ID.
    If the cell is of type ElementInfo::SEGMENT or ElementInfo::POINT
    returns 0.0.

    \param[in] id cell id

    \result maximal edge length
*/
double SurfTriPatch::evalMaxEdgeLength(const long &id)
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

    double edge_length = std::numeric_limits<double>::min();
    vector<int> face_loc_connect(2, Vertex::NULL_ID);
    n_faces = cell_->getFaceCount();
    for (i = 0; i < n_faces; ++i) {
        face_loc_connect = cell_->getFaceLocalConnect(i);
        face_loc_connect[0] = cell_->getVertex(face_loc_connect[0]);
        face_loc_connect[1] = cell_->getVertex(face_loc_connect[1]);
        edge_length = max( edge_length, norm2(m_vertices[face_loc_connect[0]].getCoords() - m_vertices[face_loc_connect[1]].getCoords()) );
    } //next i

    return(edge_length);
}

/*!
	@}
*/

}
