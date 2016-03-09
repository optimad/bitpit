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
	\class SurfUnstructured

	\brief The SurfUnstructured class defines an unstructured surface
	triangulation.

	SurfUnstructured defines an unstructured surface triangulation.
*/

/*!
	Creates a new patch.

	\param id is the id of the patch
*/
SurfUnstructured::SurfUnstructured(const int &id, int patch_dim, int space_dim)
	: SurfaceKernel(id, patch_dim, space_dim, true)
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
void SurfUnstructured::setExpert(bool expert)
{
	SurfaceKernel::setExpert(expert);
}

/*!
	Updates the patch.

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> SurfUnstructured::_update(bool trackAdaption)
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
bool SurfUnstructured::_markCellForRefinement(const long &id)
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
bool SurfUnstructured::_markCellForCoarsening(const long &id)
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
bool SurfUnstructured::_enableCellBalancing(const long &id, bool enabled)
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
bool SurfUnstructured::isPointInside(const std::array<double, 3> &point)
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
long SurfUnstructured::locatePoint(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

/*!
 * Fill adjacencies info for each cell.
*/
void SurfUnstructured::buildAdjacencies(void)
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
                simplex_idx = c_->getId();
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

            simplex_idx = c_->getId();
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
void SurfUnstructured::updateAdjacencies(const std::vector<long> &cell_ids)
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
                simplex_idx = j_->getId();
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
            simplex_idx = c_->getId();
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

//TODO: Aggiungere un metodo in SurfUnstructured per aggiungere più vertici.
/*!
 * Extract the edge network from surface mesh. If adjacencies are not built
 * edges shared by more than 1 element are counted twice. Edges are appended
 * to the content of the input SurfUnstructured
 * 
 * \param[in,out] net on output stores the network of edges
*/
void SurfUnstructured::extractEdgeNetwork(SurfUnstructured &net)
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
    net.reserveVertices(net.getVertexCount() + getVertexCount());

    // ====================================================================== //
    // ADD VERTEX TO net                                                      //
    // ====================================================================== //
    for (v_ = vertexBegin(); v_ != ve_; ++v_) {
        net.addVertex(v_->getCoords(), v_->getId());
    } //next v_

    // ====================================================================== //
    // ADD EDGES                                                              //
    // ====================================================================== //
    for (c_ = cellBegin(); c_ != ce_; ++c_) {
        id = c_->getId();
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
unsigned short SurfUnstructured::importSTL(const string &stl_name, const bool &isBinary)
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
    reserveVertices(getVertexCount() + nVertex);
    reserveCells(m_nInternals + m_nGhosts + nSimplex);

    // ====================================================================== //
    // ADD VERTICES TO MESH                                                   //
    // ====================================================================== //
    ve_ = vertexList.cend();
    for (v_ = vertexList.cbegin(); v_ != ve_; ++v_) {
        i_ = addVertex(*v_);
        vertexMap[v_counter] = i_->getId();
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
unsigned short SurfUnstructured::exportSTL(const string &stl_name, const bool &isBinary, bool exportInternalsOnly)
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
    vertexList.resize(getVertexCount());
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
        vertexMap[v_->getId()] = v_count;

        // Update counters
        ++v_count;
        ++i_;

    } //next v_
    nVertex = getVertexCount();

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
        *i_ = std::move(evalFacetNormal(c_->getId()));
        
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
