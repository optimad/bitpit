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

#include "bitpit_common.hpp"

#include "lineunstructured.hpp"

namespace bitpit {

/*!
    \class LineUnstructured
    \ingroup linepatches

    \brief The LineUnstructured class defines an unstructured line
    tasselation.

    LineUnstructured defines an unstructured line tasselation.
*/

#if BITPIT_ENABLE_MPI==1
/*!
    Creates an uninitialized partitioned patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
LineUnstructured::LineUnstructured(MPI_Comm communicator)
    : LineKernel(communicator, 1, true)
#else
/*!
    Creates an uninitialized serial patch.
*/
LineUnstructured::LineUnstructured()
    : LineKernel(true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param dimension is the dimension of the patch
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
LineUnstructured::LineUnstructured(int dimension, MPI_Comm communicator)
    : LineKernel(PatchManager::AUTOMATIC_ID, dimension, communicator, 1, true)
#else
/*!
    Creates a patch.

    \param dimension is the dimension of the patch
*/
LineUnstructured::LineUnstructured(int dimension)
    : LineKernel(PatchManager::AUTOMATIC_ID, dimension, true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param id is the id of the patch
    \param dimension is the dimension of the patch
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
LineUnstructured::LineUnstructured(int id, int dimension, MPI_Comm communicator)
    : LineKernel(id, dimension, communicator, 1, true)
#else
/*!
    Creates a patch.

    \param id is the id of the patch
    \param dimension is the dimension of the patch
*/
LineUnstructured::LineUnstructured(int id, int dimension)
    : LineKernel(id, dimension, true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch restoring the patch saved in the specified stream.

    The number of processes in the communicator should be equal to the number
    of processes of the communicator used when dumping the patch.

    \param stream is the stream to read from
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
LineUnstructured::LineUnstructured(std::istream &stream, MPI_Comm communicator)
    : LineKernel(communicator, 1, false)
#else
/*!
    Creates a patch restoring the patch saved in the specified stream.

    \param stream is the stream to read from
*/
LineUnstructured::LineUnstructured(std::istream &stream)
    : LineKernel(false)
#endif
{
    // Restore the patch
    restore(stream);
}

/*!
    Creates a clone of the pach.

    \result A clone of the pach.
*/
std::unique_ptr<PatchKernel> LineUnstructured::clone() const
{
    return std::unique_ptr<LineUnstructured>(new LineUnstructured(*this));
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
void LineUnstructured::setExpert(bool expert)
{
    LineKernel::setExpert(expert);
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int LineUnstructured::_getDumpVersion() const
{
    const int DUMP_VERSION = 1;

    return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void LineUnstructured::_dump(std::ostream &stream) const
{
#if BITPIT_ENABLE_MPI==1
    // Dump works only for serial calculations
    if (getProcessorCount() != 1) {
        throw std::runtime_error ("Dump of lineunstructured is implemented only for serial calculations.");
    }
#endif

    // Save the vertices
    dumpVertices(stream);

    // Save the cells
    dumpCells(stream);

    // Save the interfaces
    dumpInterfaces(stream);
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void LineUnstructured::_restore(std::istream &stream)
{
#if BITPIT_ENABLE_MPI==1
    // Restore works only for serial calculations
    if (getProcessorCount() != 1) {
        throw std::runtime_error ("Restore of lineunstructured is implemented only for serial calculations.");
    }
#endif

    // Restore the vertices
    restoreVertices(stream);

    // Restore the cells
    restoreCells(stream);

    // Restore the interfaces
    restoreInterfaces(stream);
}

/*!
 * Locates the cell the contains the point.
 *
 * If the point is not inside the patch, the function returns the id of the
 * null element.
 *
 * NOTE: this function is not implemented yet.
 *
 * \param[in] point is the point to be checked
 * \result Returns the linear id of the cell the contains the point. If the
 * point is not inside the patch, the function returns the id of the null
 * element.
 */
long LineUnstructured::locatePoint(const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

    throw std::runtime_error ("The function 'locatePoint' is not implemented yet");

    return false;
}

/*!
 * Import line tasselation from DGF file.
 * 
 * \param[in] dgf_name name of dgf file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 * 
 * \result on output returns an error flag for I/O error.
*/
unsigned short LineUnstructured::importDGF(const std::string &dgf_name, int PIDOffset, bool PIDSquash)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    DGFObj                                                      dgf_in(dgf_name);
    int                                                         nV = 0, nS = 0;
    long                                                        vcount, idx;
    std::vector<std::array<double, 3>>                          vertex_list;
    std::vector<std::vector<int>>                               simplex_list;
    std::vector<int>                                            simplex_PID;
    std::vector<long>                                           vertex_map;
    std::vector<long>                                           connect;

    // Counters
    std::vector<std::array<double, 3>>::const_iterator          v_, ve_;
    std::vector<std::vector<int>>::iterator                     c_, ce_;
    std::vector<int>::iterator                                  i_, ie_;
    std::vector<long>::iterator                                 j_, je_;

    // ====================================================================== //
    // IMPORT DATA                                                            //
    // ====================================================================== //

    // Read vertices and cells from DGF file
    dgf_in.load(nV, nS, vertex_list, simplex_list, simplex_PID);

    // Add vertices
    ve_ = vertex_list.cend();
    vcount = 0;
    vertex_map.resize(nV);
    for (v_ = vertex_list.cbegin(); v_ != ve_; ++v_) {
        idx = addVertex(*v_)->getId();
        vertex_map[vcount] = idx;
        ++vcount;
    } //next v_

    // Update connectivity infos
    ce_ = simplex_list.end();
    for (c_ = simplex_list.begin(); c_ != ce_; ++c_) {
        ie_ = c_->end();
        for (i_ = c_->begin(); i_ != ie_; ++i_) {
            *i_ = vertex_map[*i_];
        } //next i_
    } //next c_

    // Add cells
    int k;
    for (c_ = simplex_list.begin(), k = 0; c_ != ce_; ++c_, ++k) {
        // Create cell
        i_ = c_->begin();
        connect.resize(c_->size(), Vertex::NULL_ID);
        je_ = connect.end();
        for (j_ = connect.begin(); j_ != je_; ++j_) {
            *j_ = *i_;
            ++i_;
        } //next j_
        CellIterator cellIterator = addCell(getDGFFacetType(c_->size()), connect);

        // Set cell PID
        int cellPID = PIDOffset;
        if (!PIDSquash) {
            cellPID += simplex_PID[k];
        }

        cellIterator->setPID(cellPID);
    } //next c_

    return 0;
}

/*!
 * Export line tasselation to DGF file
 * 
 * \param[in] dgf_name name of dgf file
 * 
 * \result on output returns an error flag for I/O error
*/
unsigned short LineUnstructured::exportDGF(const std::string &dgf_name)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    DGFObj                                                      dgf_in(dgf_name);
    int                                                         nV = getVertexCount(), nS = getCellCount();
    int                                                         v, nv;
    long                                                        vcount, ccount, idx;
    std::vector<std::array<double, 3>>                          vertex_list(nV);
    std::vector<std::vector<int>>                               simplex_list(nS);
    std::unordered_map<long, long>                              vertex_map;

    // Counters
    VertexIterator                                              v_, ve_;
    CellIterator                                                c_, ce_;

    // ====================================================================== //
    // EXPORT DATA                                                            //
    // ====================================================================== //

    // Create vertex list
    ve_ = vertexEnd();
    vcount = 0;
    for (v_ = vertexBegin(); v_ != ve_; ++v_) {
        idx = v_->getId();
        vertex_list[vcount] = v_->getCoords();
        vertex_map[idx] = vcount;
        ++vcount;
    } //next v_

    // Add cells
    ce_ = cellEnd();
    ccount = 0;
    for (c_ = cellBegin(); c_ != ce_; ++c_) {
        ConstProxyVector<long> cellVertexIds = c_->getVertexIds();
        nv = cellVertexIds.size();
        simplex_list[ccount].resize(nv);
        for (v = 0; v < nv; ++v) {
            simplex_list[ccount][v] = vertex_map[cellVertexIds[v]];
        } //next v
        ++ccount;
    } //next c_

    // Read vertices and cells from DGF file
    dgf_in.save(nV, nS, vertex_list, simplex_list);

    return 0;
}

/*!
* Get the element type of a facet with the specified number of vertices.
*
* \param[in] nFacetVertices is the number of the vertices of the facet
*
* \result The element type of a facet with the specified number of vertices.
*/
ElementType LineUnstructured::getDGFFacetType(int nFacetVertices)
{
    switch(nFacetVertices) {

    case 1:
        return ElementType::VERTEX;

    case 2:
        return ElementType::LINE;

    default:
        return ElementType::UNDEFINED;

    }
}

}
