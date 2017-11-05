/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include "surfunstructured.hpp"

namespace bitpit {

/*!
	\class SurfUnstructured
	\ingroup surfacepatches

	\brief The SurfUnstructured class defines an unstructured surface
	triangulation.

	SurfUnstructured defines an unstructured surface triangulation.
*/

/*!
	Creates an uninitialized patch.
*/
SurfUnstructured::SurfUnstructured()
	: SurfaceKernel(true)
{
#if BITPIT_ENABLE_MPI==1
	// This patch supports partitioning
	setPartitioningStatus(PARTITIONING_CLEAN);
#endif
}

/*!
	Creates a new patch.

	\param patch_dim is the dimension of the patch
	\param space_dim is the dimension of the space
*/
SurfUnstructured::SurfUnstructured(int patch_dim, int space_dim)
	: SurfaceKernel(PatchManager::AUTOMATIC_ID, patch_dim, space_dim, true)
{
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param patch_dim is the dimension of the patch
	\param space_dim is the dimension of the space
*/
SurfUnstructured::SurfUnstructured(const int &id, int patch_dim, int space_dim)
	: SurfaceKernel(id, patch_dim, space_dim, true)
{
#if BITPIT_ENABLE_MPI==1
	// This patch supports partitioning
	setPartitioningStatus(PARTITIONING_CLEAN);
#endif
}

/*!
	Creates a new patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
SurfUnstructured::SurfUnstructured(std::istream &stream)
	: SurfaceKernel(false)
{
#if BITPIT_ENABLE_MPI==1
	// This patch supports partitioning
	setPartitioningStatus(PARTITIONING_CLEAN);
#endif

	// Restore the patch
	restore(stream);
}

/*!
	Creates a clone of the pach.

	\result A clone of the pach.
*/
std::unique_ptr<PatchKernel> SurfUnstructured::clone() const
{
	return std::unique_ptr<SurfUnstructured>(new SurfUnstructured(*this));
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
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int SurfUnstructured::_getDumpVersion() const
{
	const int DUMP_VERSION = 1;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void SurfUnstructured::_dump(std::ostream &stream) const
{
#if BITPIT_ENABLE_MPI==1
	// Dump works only for serial calculations
	if (getProcessorCount() != 1) {
		throw std::runtime_error ("Dump of surfunstructured is implemented only for serial calculations.");
	}
#endif

	// Space dimension
	utils::binary::write(stream, getSpaceDimension());

	// Save the vertices
	utils::binary::write(stream, getVertexCount());

	for (const Vertex &vertex : m_vertices) {
		utils::binary::write(stream, vertex.getId());

		std::array<double, 3> coords = vertex.getCoords();
		utils::binary::write(stream, coords[0]);
		utils::binary::write(stream, coords[1]);
		utils::binary::write(stream, coords[2]);
	}

	// Save the cells
	utils::binary::write(stream, getInternalCount());
	utils::binary::write(stream, getGhostCount());

	for (const Cell &cell: m_cells) {
		const ReferenceElementInfo &cellInfo = cell.getInfo();

		utils::binary::write(stream, cell.getId());
		utils::binary::write(stream, cell.getPID());
		utils::binary::write(stream, cellInfo.type);

		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int i = 0; i < nCellVertices; ++i) {
			utils::binary::write(stream, cellVertexIds[i]);
		}
	}
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void SurfUnstructured::_restore(std::istream &stream)
{
#if BITPIT_ENABLE_MPI==1
	// Restore works only for serial calculations
	if (getProcessorCount() != 1) {
		throw std::runtime_error ("Restore of surfunstructured is implemented only for serial calculations.");
	}
#endif

	// Space dimension
	int spaceDimension;
	utils::binary::read(stream, spaceDimension);
	setSpaceDimension(spaceDimension);

	// Restore the vertices
	long nVertices;
	utils::binary::read(stream, nVertices);

	reserveVertices(nVertices);
	for (long i = 0; i < nVertices; ++i) {
		long id;
		utils::binary::read(stream, id);

		std::array<double, 3> coords;
		utils::binary::read(stream, coords[0]);
		utils::binary::read(stream, coords[1]);
		utils::binary::read(stream, coords[2]);

		addVertex(coords, id);
	}

	// Restore the cells
	long nInternals;
	utils::binary::read(stream, nInternals);

	long nGhosts;
	utils::binary::read(stream, nGhosts);

	long nCells = nInternals + nGhosts;

	reserveCells(nCells);
	for (long i = 0; i < nCells; ++i) {
		long id;
		utils::binary::read(stream, id);

		int PID;
		utils::binary::read(stream, PID);

		ElementType type;
		utils::binary::read(stream, type);
		const ReferenceElementInfo &cellInfo = ReferenceElementInfo::getInfo(type);

		int nCellVertices = cellInfo.nVertices;
		std::vector<long> connect(nCellVertices, Vertex::NULL_ID);
		for (int k = 0; k < nCellVertices; ++k) {
			utils::binary::read(stream, connect[k]);
		}

		CellIterator cellIterator = addCell(type, true, connect, id);
		cellIterator->setPID(PID);
	}

	// Build ghost exchange data
#if BITPIT_ENABLE_MPI==1
	if (getProcessorCount()) {
		buildGhostExchangeData();
	}
#endif
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
long SurfUnstructured::locatePoint(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	throw std::runtime_error ("The function 'locatePoint' is not implemented yet");

	return false;
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
    int                                         n_faces, n_adj;
    long                                        id;

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
        ConstProxyVector<long> cellVertexIds = c_->getVertexIds();
        for (i = 0; i < n_faces; ++i) {
            check = true;
            n_adj = c_->getAdjacencyCount(i);
            for (j = 0; j < n_adj; ++j) {
                check = check && (id > c_->getAdjacency(i, j));
            } //next j
            if (check) {
                // Get edge type
                ElementType edgeType = c_->getFaceType(i);

                // Get edge connect
                ConstProxyVector<long> faceConnect = c_->getFaceConnect(i);
                int faceConnectSize = faceConnect.size();

                std::unique_ptr<long[]> edgeConnect = std::unique_ptr<long[]>(new long[faceConnectSize]);
                for (int k = 0; k < faceConnectSize; ++k) {
                    edgeConnect[k] = faceConnect[k];
                }

                // Add edge
                net.addCell(edgeType, true, std::move(edgeConnect));
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
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 * 
 * \param[in] stl_name name of stl file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 * 
 * \result on output returns an error flag for I/O error
*/
unsigned short SurfUnstructured::importSTL(const string &stl_name,
                                           int PIDOffset, bool PIDSquash)
{
    // Create STL object
    STLObj STL(stl_name);

    // Import stl object
    importSTL(STL, PIDOffset, PIDSquash);

    return 0;
}

/*!
 * Import surface tasselation from S.T.L. file. STL facet are added at to the
 * present mesh, i.e. current mesh content is not discarded. Howver no checks
 * are performed to ensure that no duplicated vertices or cells are created.
 *
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 *
 * \param[in] stl_name name of stl file
 * \param[in] isBinary flag for binary (true), of ASCII (false) stl file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 *
 * \result on output returns an error flag for I/O error
*/
unsigned short SurfUnstructured::importSTL(const string &stl_name, const bool &isBinary,
                                           int PIDOffset, bool PIDSquash)
{
    // Create STL object
    STLObj STL(stl_name, isBinary);

    // Import stl object
    importSTL(STL, PIDOffset, PIDSquash);

    return 0;
}

/*!
 * Import surface tasselation from S.T.L. file. STL facet are added at to the
 * present mesh, i.e. current mesh content is not discarded. Howver no checks
 * are performed to ensure that no duplicated vertices or cells are created.
 *
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 *
 * \param[in] STL is the STL object to import
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 *
 * \result on output returns an error flag for I/O error
*/
unsigned short SurfUnstructured::importSTL(STLObj &STL, int PIDOffset, bool PIDSquash)
{
    // ====================================================================== //
    // OPEN STL FILE                                                          //
    // ====================================================================== //
    STL.open("in");
    if (STL.err != 0) {
        return STL.err;
    }

    // ====================================================================== //
    // LOAD ALL SOLID FROM THE STL FILE                                       //
    // ====================================================================== //
    int pid = PIDOffset;
    if (!PIDSquash) {
        --pid;
    }

    while (true) {
        // ====================================================================== //
        // LOAD SOLID FROM THE STL FILE                                           //
        // ====================================================================== //
        int nVertex = 0;
        int nSimplex = 0;
        std::vector<std::array<double, 3>> vertexList;
        std::vector<std::array<double, 3>> normalList;
        std::vector<std::array<int, 3>> connectivityList;

        STL.loadSolid(nVertex, nSimplex, vertexList, normalList, connectivityList);
        if (nVertex == 0) {
            break;
        }

        // ====================================================================== //
        // PID OF THE SOLID                                                       //
        // ====================================================================== //
        if (!PIDSquash) {
            ++pid;
        }

        // ====================================================================== //
        // PREPARE MESH FOR DATA IMPORT                                           //
        // ====================================================================== //
        reserveVertices(getVertexCount() + nVertex);
        reserveCells(m_nInternals + m_nGhosts + nSimplex);

        // ====================================================================== //
        // ADD VERTICES TO MESH                                                   //
        // ====================================================================== //
        vector<array<double, 3>>::const_iterator v_, ve_;

        std::vector<long> vertexMap;
        vertexMap.reserve(nVertex);

        long v_counter = 0;
        ve_ = vertexList.cend();
        for (v_ = vertexList.cbegin(); v_ != ve_; ++v_) {
            VertexIterator i_ = addVertex(*v_);
            vertexMap[v_counter] = i_->getId();
            ++v_counter;
        } //next v_

        // ====================================================================== //
        // ADD CELLS TO MESH                                                      //
        // ====================================================================== //
        vector<array<int,3>>::const_iterator c_, ce_;
        array<int,3>::const_iterator w_, we_;

        ce_ = connectivityList.cend();
        for (c_ = connectivityList.cbegin(); c_ != ce_; ++c_) {
            // Remap STL connectivity
            int n_v = c_->size();
            std::vector<long> connect(n_v, Vertex::NULL_ID);
            we_ = c_->cend();
            int i = 0;
            for (w_ = c_->cbegin(); w_ < we_; ++w_) {
                connect[i] = vertexMap[*w_];
                ++i;
            } //next w_

            // Add cell
            CellIterator cellIterator = addCell(getSTLFacetType(n_v), true, connect);
            cellIterator->setPID(pid);
        } //next c_

        // ====================================================================== //
        // Multi-body STL files are supported only in ASCII mode                        //
        // ====================================================================== //
        if (STL.stl_type) {
            break;
        }
    }

    // ====================================================================== //
    // CLOSE STL FILE                                                         //
    // ====================================================================== //
    STL.close("in");

    return 0;
}

/*!
 * Export surface tasselation in a STL format. No check is perfomed on element type
 * therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation.
 *
 * \param[in] stl_name name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 * \param[in] exportInternalsOnly flag for exporting only internal cells (true), or
 * internal and ghost cells (false).
 *
 * \result on output returns an error flag for I/O error.
 */
unsigned short SurfUnstructured::exportSTL(const string &stl_name, const bool &isBinary, bool exportInternalsOnly)
{
    return exportSTLSingle(stl_name, isBinary, exportInternalsOnly);
}

/*!
 * Export surface tasselation in a STL format. No check is perfomed on element type
 * therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation. Overloading supporting the the ascii multi-solid mode export.
 *
 * \param[in] stl_name name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 * \param[in] isMulti flag to write in ASCII multi-solid mode (true) or not (false).
 * If true, isBinary flag will be ignored.
 * \param[in] exportInternalsOnly flag for exporting only internal cells (true), or
 * internal+ghost cells (false).
 *
 * \result on output returns an error flag for I/O error.
 */
unsigned short SurfUnstructured::exportSTL(const string &stl_name, const bool &isBinary, const bool &isMulti, bool exportInternalsOnly)
{
    unsigned short flag = 0;
    if (isMulti) {
        flag = exportSTLMulti(stl_name, exportInternalsOnly);
    } else {
        flag = exportSTLSingle(stl_name, isBinary, exportInternalsOnly);
    }

    return flag;
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
unsigned short SurfUnstructured::exportSTLSingle(const string &stl_name, const bool &isBinary, bool exportInternalsOnly)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    int                                         nVertex;
    int                                         nSimplex;
    vector<array<double, 3>>                    vertexList;
    vector<array<double, 3>>                    normalList;
    vector<array<int,3>>                        connectivityList;
    unordered_map<long, long>                   vertexMap;
    array<int,3>                                dummyIntArray;

    // Counters
    int                                         v_count ,j;
    vector<array<double, 3>>::iterator          i_;
    vector<array<int,3>>::iterator              j_;
    array<int,3>::iterator                      k_, ke_;
    VertexIterator                              v_, ve_;
    CellIterator                                c_, cb_, ce_;

    // ====================================================================== //
    // INITIALIZE DATA STRUCTURE                                              //
    // ====================================================================== //
    dummyIntArray.fill(0) ;

    nSimplex = m_nInternals;
    if (!exportInternalsOnly) nSimplex += m_nGhosts;
    vertexList.resize(getVertexCount());
    normalList.resize(nSimplex);
    connectivityList.resize(nSimplex, dummyIntArray);

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
        ConstProxyVector<long> cellVertexIds = c_->getVertexIds();

        ke_ = j_->end();
        j = 0;
        for (k_ = j_->begin(); k_ != ke_; ++k_) {
            *k_ = vertexMap[cellVertexIds[j]];
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
 * Export surface tasselation in a STL Multi Solid format, in ascii mode only. Binary is not supported for STL Multisolid.
 * No check is perfomed on element type therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation. If available, ghost cells will be written in a stand-alone solid.
 *
 * \param[in] stl_name name of the stl file
 * \param[in] exportInternalsOnly flag for exporting only internal cells (true),
 * or internal+ghost cells (false).
 *
 * \result on output returns an error flag for I/O error 0-done, >0 errors.
 */
unsigned short SurfUnstructured::exportSTLMulti(const string &stl_name, bool exportInternalsOnly)
{
    int                                         nTotVertex;
    vector<array<double, 3>>                    totVertexList;
    unordered_map<long, long>                   vertexMap;

    int                                         nLocSimplex;
    vector<array<double, 3>>                    normalLocList;
    vector<array<int,3>>                        connectivityLocList;

    STLObj STL(stl_name, false);

    STL.open("out");
    if (STL.err != 0) {
        return STL.err;
    }

    // Create the vertex map
    nTotVertex = getVertexCount();
    totVertexList.resize(nTotVertex);
    unsigned int count_v = 0;
    for(const Vertex &v: getVertices()){
        totVertexList[count_v] = v.getCoords();
        vertexMap[v.getId()] = count_v;
        ++count_v;
    }

    // Export the internal cells
    for(int pid : getInternalPIDs()){
        std::vector<long> cells = getInternalsByPID(pid);

        nLocSimplex = (int) cells.size();
        connectivityLocList.resize(nLocSimplex);
        normalLocList.resize(nLocSimplex);

        vector<array<int,3>>::iterator itC = connectivityLocList.begin();
        vector<array<double, 3>>::iterator itN = normalLocList.begin();

        // Fill local connectivity and normals structures
        for (auto id : cells) {
            // Fill connectivity
            const Cell &cell = getCell(id);
            for (int iloc = 0; iloc<3; ++iloc) {
                (*itC)[iloc] = vertexMap[cell.getConnect()[iloc]];
            }

            // Fill normal
            *itN = std::move(evalFacetNormal(cell.getId()));

            // Increment  iterators
            ++itC;
            ++itN;
        }

        // Write the solid associated to the current PID
        STL.saveSolid(std::to_string(pid), nTotVertex, nLocSimplex, totVertexList, normalLocList, connectivityLocList);
    }

    // Export ghost cells
    long nGhosts = getGhostCount();
    if (!exportInternalsOnly && nGhosts > 0) {
        nLocSimplex = nGhosts;
        connectivityLocList.resize(nLocSimplex);
        normalLocList.resize(nLocSimplex);

        vector<array<int,3>>::iterator itC = connectivityLocList.begin();
        vector<array<double, 3>>::iterator itN = normalLocList.begin();

        // Fill local connectivity and normals structures
		CellConstIterator endItr = ghostConstEnd();
        for (CellConstIterator itr = ghostConstBegin(); itr != endItr; ++itr) {
            // Fill connectivity
            for (int iloc = 0; iloc<3; ++iloc) {
                (*itC)[iloc] = vertexMap[itr->getConnect()[iloc]];
            }

            // Fill normal
            *itN = std::move(evalFacetNormal(itr.getId()));

            // Increment iterators
            ++itC;
            ++itN;
        }

        // Write the solid associated to the ghosts
        STL.saveSolid("ghosts", nTotVertex, nLocSimplex, totVertexList, normalLocList, connectivityLocList);
    }

    STL.close("out");

    return 0;
}

/*!
* Get the element type of a facet with the specified number of vertices.
*
* \param[in] nFacetVertices is the number of the vertices of the facet
*
* \result The element type of a facet with the specified number of vertices.
*/
ElementType SurfUnstructured::getSTLFacetType(int nFacetVertices)
{
    switch(nFacetVertices) {

    case 1:
        return ElementType::VERTEX;

    case 2:
        return ElementType::LINE;

    case 3:
        return ElementType::TRIANGLE;

    case 4:
        return ElementType::QUAD;

    default:
        return ElementType::UNDEFINED;

    }
}

/*!
 * Import surface tasselation from DGF file.
 * 
 * \param[in] dgf_name name of dgf file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 * 
 * \result on output returns an error flag for I/O error.
*/
unsigned short SurfUnstructured::importDGF(const string &dgf_name, int PIDOffset, bool PIDSquash)
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
        CellIterator cellIterator = addCell(getDGFFacetType(c_->size()), true, connect);

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
 * Export surface tasselation to DGF file
 * 
 * \param[in] dgf_name name of dgf file
 * 
 * \result on output returns an error flag for I/O error
*/
unsigned short SurfUnstructured::exportDGF(const string &dgf_name)
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
ElementType SurfUnstructured::getDGFFacetType(int nFacetVertices)
{
    switch(nFacetVertices) {

    case 1:
        return ElementType::VERTEX;

    case 2:
        return ElementType::LINE;

    case 3:
        return ElementType::TRIANGLE;

    case 4:
        return ElementType::QUAD;

    default:
        return ElementType::UNDEFINED;

    }
}

}
