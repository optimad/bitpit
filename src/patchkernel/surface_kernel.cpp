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

#include "surface_kernel.hpp"

namespace bitpit {


/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class SurfaceKernel

	\brief The SurfaceKernel class provides an interface for defining
	surface patches.

	SurfaceKernel is the base class for defining surface patches.
*/

const std::map<ElementInfo::Type, unsigned short>  SurfaceKernel::m_selectionTypes({
                                {ElementInfo::QUAD,     SurfaceKernel::SELECT_QUAD},
                                {ElementInfo::TRIANGLE, SurfaceKernel::SELECT_TRIANGLE}
                                });
const unsigned short SurfaceKernel::SELECT_TRIANGLE = 1;
const unsigned short SurfaceKernel::SELECT_QUAD     = 2;
const unsigned short SurfaceKernel::SELECT_ALL      = 3;

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
 * Returns the number of dimensions of the working space (set at patch construction)
 */
int SurfaceKernel::getSpaceDimensions(void) const
{
    return(m_spaceDim);
}

/*!
        Evaluates the characteristic size of the specified cell.

        \param id is the id of the cell
        \result The characteristic size of the specified cell.
*/
double SurfaceKernel::evalCellSize(const long &id)
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
 * ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet area
*/
double SurfaceKernel::evalCellArea(const long &id)
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
     || (cell_->getType() == ElementInfo::VERTEX)) return 0.0;
    if (cell_->getType() == ElementInfo::LINE) {
        return ( norm2(m_vertices[cell_->getVertex(0)].getCoords()
                     - m_vertices[cell_->getVertex(1)].getCoords()) );
    }
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
 * Compute cell center coordinates. Cell center coordinates are computed
 * as the arithmetic average of cell's vertex coordinates.
 * 
 * \param[in] id cell id
 * 
 * \result returns cell's center coordinates
*/
std::array<double, 3> SurfaceKernel::evalCellCentroid(const long &id)
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
double SurfaceKernel::evalEdgeLength(const long &id, const int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    // none

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
double SurfaceKernel::evalMinEdgeLength(const long &id, int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 i;
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
double SurfaceKernel::evalMaxEdgeLength(const long &id, int &edge_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    Cell                                *cell_ = &m_cells[id];

    // Counters
    int                                 i;
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
double SurfaceKernel::evalAngleAtVertex(const long &id, const int &vertex_id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    double const                _pi_ = 3.14159265358979;
    Cell                        *cell_ = &m_cells[id];

    // Counters

    // ====================================================================== //
    // EVALUATE ANGLE AT SPECIFIED VERTEX                                     //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)) return 0.0;
    if (cell_->getType() == ElementInfo::LINE) {
        if (m_spaceDim - getDimension() == 1)   return _pi_;
        else                                    return 0.0;
    }

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
double SurfaceKernel::evalMinAngleAtVertex(const long&id, int &vertex_id)
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
double SurfaceKernel::evalMaxAngleAtVertex(const long&id, int &vertex_id)
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
 * Evalute the aspect ratio for a cell with specified ID. The aspect ratio
 * is defined as the ratio between the longest and the shortest edge.
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * \param[in,out] edge_id on output stores the index of the shortest edge
 * 
 * \result cell aspect ratio
*/
double SurfaceKernel::evalAspectRatio(const long &id, int &edge_id)
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
 * If cell is of type ElementInfo::VERTEX or ElementInfo::LINE, returns 0.0
 * 
 * \param[in] id cell ID
 * 
 * \result facet normal
*/
array<double, 3> SurfaceKernel::evalFacetNormal(const long &id)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    array<double, 3>             normal = {{0.0, 0.0, 0.0}};
    Cell                        *cell_ = &m_cells[id];

    // Counters
    // none

    // ====================================================================== //
    // COMPUTE NORMAL                                                         //
    // ====================================================================== //
    if ((cell_->getType() == ElementInfo::UNDEFINED)
     || (cell_->getType() == ElementInfo::VERTEX)) return normal;
    
    if (cell_->getType() == ElementInfo::LINE) {
        if (m_spaceDim - getDimension() == 1) {
            std::array<double, 3>       z = {{0.0, 0.0, 1.0}};
            normal = m_vertices[cell_->getVertex(1)].getCoords() - m_vertices[cell_->getVertex(0)].getCoords();
            normal = crossProduct(normal, z);
        }
        else return normal;
    }

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
array<double, 3> SurfaceKernel::evalEdgeNormal(const long &id, const int &edge_id)
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
    if (cell_->getAdjacency(edge_id, 0) != Element::NULL_ID) {
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
array<double, 3> SurfaceKernel::evalVertexNormal(const long &id, const int &vert_id)
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
        angles[i] = evalAngleAtVertex(cell_->getId(), cell_->findVertex(vert_ID));
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
 * Display quality stats for the mesh currently stored.
 * 
 * \param[in,out] out output stream where stats will be printed out
 * \param[in] padding (default = 0) number of trailing spaces
*/
void SurfaceKernel::displayQualityStats(ostream& out, unsigned int padding)
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
vector<double> SurfaceKernel::computeHistogram(eval_f_ funct_, vector<double> &bins, long &count, int n_intervals, unsigned short mask)
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
 * Given and element type and the selection mask, returns true if the element type
 * is accounted for in the selection mask.
 * 
 * \param[in] mask_ selection mask_
 * \param[in] type_ element type
*/
bool SurfaceKernel::compareSelectedTypes(const unsigned short &mask_, const ElementInfo::Type &type_)
{
    unsigned short       masked = m_selectionTypes.at(type_);
    return ( (mask_ & masked) == masked );
}

/*!
	@}
*/

}
