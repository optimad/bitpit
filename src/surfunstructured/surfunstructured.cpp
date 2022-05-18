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
#include "bitpit_IO.hpp"

#include "surfunstructured.hpp"

namespace bitpit {

/*!
	\class SurfUnstructured
	\ingroup surfacepatches

	\brief The SurfUnstructured class defines an unstructured surface
	triangulation.

	SurfUnstructured defines an unstructured surface triangulation.
*/

#if BITPIT_ENABLE_MPI==1
/*!
	Creates an uninitialized partitioned patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
SurfUnstructured::SurfUnstructured(MPI_Comm communicator, std::size_t haloSize)
	: SurfaceKernel(communicator, haloSize, true)
#else
/*!
	Creates an uninitialized serial patch.
*/
SurfUnstructured::SurfUnstructured()
	: SurfaceKernel(true)
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
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
SurfUnstructured::SurfUnstructured(int dimension, MPI_Comm communicator, std::size_t haloSize)
	: SurfaceKernel(PatchManager::AUTOMATIC_ID, dimension, communicator, haloSize, true)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
*/
SurfUnstructured::SurfUnstructured(int dimension)
	: SurfaceKernel(PatchManager::AUTOMATIC_ID, dimension, true)
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
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
SurfUnstructured::SurfUnstructured(int id, int dimension, MPI_Comm communicator, std::size_t haloSize)
	: SurfaceKernel(id, dimension, communicator, haloSize, true)
#else
/*!
	Creates a patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
*/
SurfUnstructured::SurfUnstructured(int id, int dimension)
	: SurfaceKernel(id, dimension, true)
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
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
SurfUnstructured::SurfUnstructured(std::istream &stream, MPI_Comm communicator, std::size_t haloSize)
	: SurfaceKernel(communicator, haloSize, false)
#else
/*!
	Creates a patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
SurfUnstructured::SurfUnstructured(std::istream &stream)
	: SurfaceKernel(false)
#endif
{
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
	const int DUMP_VERSION = 5;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void SurfUnstructured::_dump(std::ostream &stream) const
{
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
void SurfUnstructured::_restore(std::istream &stream)
{
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
long SurfUnstructured::locatePoint(const std::array<double, 3> &point) const
{
	BITPIT_UNUSED(point);

	throw std::runtime_error ("The function 'locatePoint' is not implemented yet");

	return Cell::NULL_ID;
}

//TODO: Aggiungere un metodo in SurfUnstructured per aggiungere più vertici.
/*!
 * Extract the edge network from surface mesh. If adjacencies are not built
 * edges shared by more than 1 element are counted twice. Edges are appended
 * to the content of the input SurfUnstructured
 * 
 * \param[in,out] net on output stores the network of edges
*/
void SurfUnstructured::extractEdgeNetwork(LineUnstructured &net)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    bool                                        check;
    int                                         n_faces, n_adj;
    long                                        id;
    const long                                  *adjacencies;

    // Counters
    int                                         i, j;
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
            adjacencies = c_->getAdjacencies(i);
            for (j = 0; j < n_adj; ++j) {
                check = check && (id > adjacencies[j]);
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
                net.addCell(edgeType, std::move(edgeConnect));
            }
        } //next i
    } //next c_

    return;
}

//TODO: normals??
//TODO: error flag on output
//TODO: import a specified solid (ascii format only)
/*!
 * Import surface tasselation from STL file.
 *
 * A separate set of vertices will be created for with each facet.
 *
 * STL facets are added to the present mesh, i.e. current mesh content is not
 * discarded. However, factes are not joined to the cells of the current mesh.
 * After importing the STL, the resulting mesh can contain duplicate vertices
 * and cells.
 *
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 *
 * \param[in] filename name of stl file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 *
 * \result on output returns an error flag for I/O error
*/
int SurfUnstructured::importSTL(const std::string &filename,
                                int PIDOffset, bool PIDSquash)
{
    return importSTL(filename, STLReader::FormatUnknown, false, PIDOffset, PIDSquash);
}

/*!
 * Import surface tasselation from STL file.
 *
 * A separate set of vertices will be created for with each facet.
 *
 * STL facets are added to the present mesh, i.e. current mesh content is not
 * discarded. However, factes are not joined to the cells of the current mesh.
 * After importing the STL, the resulting mesh can contain duplicate vertices
 * and cells.
 *
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 *
 * \param[in] filename name of stl file
 * \param[in] isBinary flag for binary (true), of ASCII (false) stl file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 * \param[in,out] PIDNames are the names of the PIDs, on output the names
 * of the newly imported PIDs will be added to the list
 *
 * \result on output returns an error flag for I/O error
*/
int SurfUnstructured::importSTL(const std::string &filename, bool isBinary,
                                int PIDOffset, bool PIDSquash,
                                std::unordered_map<int, std::string> *PIDNames)
{
    STLReader::Format format;
    if (isBinary) {
        format = STLReader::FormatBinary;
    } else {
        format = STLReader::FormatASCII;
    }

    return importSTL(filename, format, false, PIDOffset, PIDSquash, PIDNames);
}

/*!
 * Import surface tasselation from STL file.
 *
 * It is possible to control if STL factes that share the same vertices will
 * be joined together or if a separate set of vertices will be created for
 * with each facet. When faces are joined together, vertices with a distance
 * less than (10 * machine epsilon) are considered conincident and will be
 * merged together.
 *
 * STL facets are added to the present mesh, i.e. current mesh content is not
 * discarded. However, factes are not joined to the cells of the current mesh.
 * After importing the STL, the resulting mesh can contain duplicate vertices
 * and cells.
 *
 * If the input file is a multi-solid ASCII file, all solids will be loaded
 * and a different PID will be assigned to the PID of the different solids.
 * Solids will not be joined together, only facets whithin a solid can be
 * joined.
 *
 * \param[in] filename name of stl file
 * \param[in] format is the format of stl file
 * \param[in] joinFacets if set to true, facets sharing the same vertices will
 * be joined together, otherwise a separate set of vertices will be created for
 * each facet. In any case, factes will not be joined to the cells of the
 * current mesh
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file (false) or if the same PID will be assigned to all cells (true).
 * \param[in,out] PIDNames are the names of the PIDs, on output the names
 * of the newly imported PIDs will be added to the list
 *
 * \result on output returns an error flag for I/O error
*/
int SurfUnstructured::importSTL(const std::string &filename, STLReader::Format format,
                                bool joinFacets, int PIDOffset, bool PIDSquash,
                                std::unordered_map<int, std::string> *PIDNames)
{
    int readerError;

    // Initialize reader
    STLReader reader(filename, format);
    if (format == STLReader::FormatUnknown) {
        format = reader.getFormat();
    }

    // Begin reding the STL file
    readerError = reader.readBegin();
    if (readerError != 0) {
        return readerError;
    }

    // Initialize cache needed for joinin the facets
    //
    // The cache will contain the ids of the vertices that have been added
    // to the patch. It uses a special comparator that allows to compare
    // the coordinates of a candidate vertex with the coordinates of the
    // vertices in the cache.
    std::array<double, 3> *candidateVertexCoords;
    Vertex::Less vertexLess(10 * std::numeric_limits<double>::epsilon());

    auto vertexCoordsLess = [this, &vertexLess, &candidateVertexCoords](const long &id_1, const long &id_2)
    {
        const std::array<double, 3> &coords_1 = (id_1 >= 0) ? this->getVertex(id_1).getCoords() : *candidateVertexCoords;
        const std::array<double, 3> &coords_2 = (id_2 >= 0) ? this->getVertex(id_2).getCoords() : *candidateVertexCoords;

        return vertexLess(coords_1, coords_2);
    };

    std::set<long, decltype(vertexCoordsLess)> vertexCache(vertexCoordsLess);

    // Read all the solids in the STL file
    int pid = PIDOffset;

    ElementType facetType = ElementType::TRIANGLE;
    const int nFacetVertices = ReferenceElementInfo::getInfo(ElementType::TRIANGLE).nVertices;

    while (true) {
        // Read header
        std::size_t nFacets;
        std::string name;
        readerError = reader.readHeader(&name, &nFacets);
        if (readerError != 0) {
            if (readerError == -2) {
                break;
            } else {
                return readerError;
            }
        }

        // Generate patch cells from STL facets
        reserveCells(getCellCount() + nFacets);

        std::size_t nEstimatedVertices;
        if (joinFacets) {
            // The number of facets and the number of vertices of a triangulation
            // are related by the inequality nFacets <= 2 * nVertices - 4, where
            // the equality holds for a closed triangulation. To limit memory usage
            // we estimate the number of vertices assuming a closed triangulation
            // (this gives us the minimum number of nodes the triangulation could
            // possibly have).
            nEstimatedVertices = 0.5 * nFacets + 2;
        } else {
            nEstimatedVertices = nFacetVertices * nFacets;
        }
        reserveVertices(getVertexCount() + nEstimatedVertices);

        vertexCache.clear();

        for (std::size_t n = 0; n < nFacets; ++n) {
            // Read facet data
            std::array<Vertex, 3> faceVertices;
            std::array<double, 3> facetNormal;
            readerError = reader.readFacet(&(faceVertices[0].getCoords()), &(faceVertices[1].getCoords()), &(faceVertices[2].getCoords()), &facetNormal);
            if (readerError != 0) {
                return readerError;
            }

            // Add vertices
            std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[nFacetVertices]);
            if (joinFacets) {
                for (int i = 0; i < nFacetVertices; ++i) {
                    // The candidate vertex is set equal to the face vertex
                    // that needs to be added, then the cache is checked: if
                    // a vertex with the same coordinates is already in the
                    // cache, the vertex will not be added and the existing
                    // one will be used.
                    //
                    // The cache contains the ids of the vertices that have
                    // been added. It uses a special comparator that allows
                    // to compare the candidate vertex with the ones in the
                    // cache, in order to do that the null vertex should be
                    // passed to the find function.
                    candidateVertexCoords = &(faceVertices[i].getCoords());

                    long vertexId;
                    auto vertexCacheItr = vertexCache.find(Vertex::NULL_ID);
                    if (vertexCacheItr == vertexCache.end()) {
                        VertexIterator vertexItr = addVertex(*candidateVertexCoords);
                        vertexId = vertexItr.getId();
                        vertexCache.insert(vertexId);
                    } else {
                        vertexId = *vertexCacheItr;
                    }
                    connectStorage[i] = vertexId;
                }
            } else {
                for (int i = 0; i < nFacetVertices; ++i) {
                    VertexIterator vertexItr = addVertex(faceVertices[i].getCoords());
                    long vertexId = vertexItr.getId();
                    connectStorage[i] = vertexId;
                }
            }

            // Add cell
            CellIterator cellIterator = addCell(facetType, std::move(connectStorage));
            cellIterator->setPID(pid);
        }

        // Read footer
        readerError = reader.readFooter(name);
        if (readerError != 0) {
            return readerError;
        }

        // Assign PID name
        if (!PIDSquash) {
            if (PIDNames && !name.empty()) {
                PIDNames->insert({{pid, name}});
            }
            ++pid;
        }

        // Multi-body STL files are supported only in ASCII mode
        if (format == STLReader::FormatBinary) {
            break;
        }
    }

    // End reading
    readerError = reader.readEnd();
    if (readerError != 0) {
        return readerError;
    }

    return 0;
}

/*!
 * Export surface tasselation in a STL format. No check is perfomed on element type
 * therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation.
 *
 * \param[in] filename name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 *
 * \result on output returns an error flag for I/O error.
 */
int SurfUnstructured::exportSTL(const std::string &filename, bool isBinary)
{
    return exportSTLSingle(filename, isBinary);
}

/*!
 * Export surface tasselation in a STL format. No check is perfomed on element type
 * therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation. Overloading supporting the the ascii multi-solid mode export.
 *
 * \param[in] filename name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 * \param[in] isMulti flag to write in ASCII multi-solid mode (true) or not (false).
 * If true, isBinary flag will be ignored.
 * \param[in,out] PIDNames are the names of the PIDs, if a PIDs has no name
 * its number will be used
 * \result on output returns an error flag for I/O error.
 */
int SurfUnstructured::exportSTL(const std::string &filename, bool isBinary, bool isMulti,
                                std::unordered_map<int, std::string> *PIDNames)
{
    int flag = 0;
    if (isMulti) {
        flag = exportSTLMulti(filename, PIDNames);
    } else {
        flag = exportSTLSingle(filename, isBinary);
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
 * \param[in] filename name of the stl file
 * \param[in] isBinary flag for binary (true) or ASCII (false) file
 * \result on output returns an error flag for I/O error.
*/
int SurfUnstructured::exportSTLSingle(const std::string &filename, bool isBinary)
{
    int writerError;

    // Initialize writer
    STLReader::Format format;
    if (isBinary) {
        format = STLReader::FormatBinary;
    } else {
        format = STLReader::FormatASCII;
    }

    STLWriter writer(filename, format);

    // Solid name
    //
    // An empty solid name will be used.
    const std::string name = "";

    // Detect if this is the master writer
    bool isMasterWriter;
#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        isMasterWriter = (getRank() == 0);
    } else {
#else
    {
#endif
        isMasterWriter = true;
    }

    // Detect if the file will be written incrementally
    bool incrementalWrite;
#if BITPIT_ENABLE_MPI==1
        incrementalWrite = true;
#else
        incrementalWrite = false;
#endif

    // Count the number of facets
    long nFacets = getInternalCellCount();
#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &nFacets, 1, MPI_LONG, MPI_SUM, getCommunicator());
    }
#endif

    // Write header
    if (isMasterWriter) {
        // Begin writing
        writerError = writer.writeBegin(STLWriter::WriteOverwrite, incrementalWrite);
        if (writerError != 0) {
            writer.writeEnd();

            return writerError;
        }

        // Write header section
        writerError = writer.writeHeader(name, nFacets);
        if (writerError != 0) {
            writer.writeEnd();

            return writerError;
        }

#if BITPIT_ENABLE_MPI==1
        // Close the file to let other process write into it
        writerError = writer.writeEnd();
        if (writerError != 0) {
            return writerError;
        }
#endif
    }

    // Write facet data
#if BITPIT_ENABLE_MPI==1
    int nRanks = getProcessorCount();
    for (int i = 0; i < nRanks; ++i) {
        if (i == getRank()) {
            // Begin writing
            writerError = writer.writeBegin(STLWriter::WriteAppend, incrementalWrite);
            if (writerError != 0) {
                writer.writeEnd();
                return writerError;
            }

#else
    {
#endif
            // Write facet section
            CellConstIterator cellBegin = internalCellConstBegin();
            CellConstIterator cellEnd   = internalCellConstEnd();
            for (CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
                // Get cell
                const Cell &cell = *cellItr;

                // Get vertex coordinates
                ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
                assert(cellVertexIds.size() == 3);

                const std::array<double, 3> &coords_0 = getVertex(cellVertexIds[0]).getCoords();
                const std::array<double, 3> &coords_1 = getVertex(cellVertexIds[1]).getCoords();
                const std::array<double, 3> &coords_2 = getVertex(cellVertexIds[2]).getCoords();

                // Evaluate normal
                const std::array<double, 3> normal = evalFacetNormal(cell.getId());

                // Write data
                writerError = writer.writeFacet(coords_0, coords_1, coords_2, normal);
                if (writerError != 0) {
                    writer.writeEnd();
                    return writerError;
                }
            }

#if BITPIT_ENABLE_MPI==1
            // Close the file to let other process write into it
            writerError = writer.writeEnd();
            if (writerError != 0) {
                return writerError;
            }
        }

        if (isPartitioned()) {
            MPI_Barrier(getCommunicator());
        }
#endif
    }

    // Write footer
    if (isMasterWriter) {
#if BITPIT_ENABLE_MPI==1
        // Begin writning
        //
        // The file was previously closed to let other process write into it.
        writerError = writer.writeBegin(STLWriter::WriteAppend, incrementalWrite);
        if (writerError != 0) {
            writer.writeEnd();
            return writerError;
        }
#endif

        // Write footer section
        writerError = writer.writeFooter(name);
        if (writerError != 0) {
            writer.writeEnd();

            return writerError;
        }

        // End writing
        writerError = writer.writeEnd();
        if (writerError != 0) {
            return writerError;
        }
    }

    return 0;
}

/*!
 * Export surface tasselation in a STL Multi Solid format, in ascii mode only. Binary is not supported for STL Multisolid.
 * No check is perfomed on element type therefore tasselation containing vertex, line or quad elements will produce
 * ill-formed stl triangulation. If available, ghost cells will be written in a stand-alone solid.
 *
 * \param[in] filename name of the stl file
 * \param[in,out] PIDNames are the names of the PIDs, if a PIDs has no name
 * its number will be used
 * \result on output returns an error flag for I/O error 0-done, >0 errors.
 */
int SurfUnstructured::exportSTLMulti(const std::string &filename, std::unordered_map<int, std::string> *PIDNames)
{
    int writerError;

    // Initialize writer
    STLWriter writer(filename, STLReader::FormatASCII);

    // Detect if this is the master writer
    bool isMasterWriter;
#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        isMasterWriter = (getRank() == 0);
    } else {
#else
    {
#endif
        isMasterWriter = true;
    }

    // Detect if the file will be written incrementally
    bool incrementalWrite = true;

#if BITPIT_ENABLE_MPI==1
    // Count number of ranks
    int nRanks = getProcessorCount();
#endif

    // Get the global list of PID
    std::set<int> cellPIDs;
#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        std::set<int> internalPIDs = getInternalCellPIDs();
        std::vector<int> localPIDs(internalPIDs.begin(), internalPIDs.end());
        int nLocalPids = localPIDs.size();

        std::vector<int> gatherPIDCount(nRanks);
        MPI_Allgather(&nLocalPids, 1, MPI_INT, gatherPIDCount.data(), 1, MPI_INT, getCommunicator());

        std::vector<int> gatherPIDDispls(nRanks, 0);
        for (int i = 1; i < nRanks; ++i) {
            gatherPIDDispls[i] = gatherPIDDispls[i - 1] + gatherPIDCount[i - 1];
        }

        int gatherPIDsSize = gatherPIDDispls.back() + gatherPIDCount.back();

        std::vector<int> globalPIDs(gatherPIDsSize);
        MPI_Allgatherv(localPIDs.data(), nLocalPids, MPI_INT, globalPIDs.data(),
                       gatherPIDCount.data(), gatherPIDDispls.data(), MPI_INT, getCommunicator());

        cellPIDs = std::set<int>(globalPIDs.begin(), globalPIDs.end());
    } else {
#else
    {
#endif
        cellPIDs = getInternalCellPIDs();
    }

    // Export cells
    bool firstPID = true;
    for (int pid : cellPIDs) {
        // Cells associated with the PID
        std::vector<long> cells = getInternalCellsByPID(pid);

        // Count the number of facets
        long nFacets = cells.size();
#if BITPIT_ENABLE_MPI==1
        if (isPartitioned()) {
            MPI_Allreduce(MPI_IN_PLACE, &nFacets, 1, MPI_LONG, MPI_SUM, getCommunicator());
        }
#endif

        // Solid name
        std::string name;
        if (PIDNames && PIDNames->count(pid) > 0) {
            name = PIDNames->at(pid);
        } else {
            name = std::to_string(pid);
        }

        // Write header
        if (isMasterWriter) {
            // Begin writing
            STLWriter::WriteMode writeMode;
            if (firstPID) {
                writeMode = STLWriter::WriteOverwrite;
            } else {
                writeMode = STLWriter::WriteAppend;
            }

            writerError = writer.writeBegin(writeMode, incrementalWrite);
            if (writerError != 0) {
                writer.writeEnd();

                return writerError;
            }

            // Write header section
            writerError = writer.writeHeader(name, nFacets);
            if (writerError != 0) {
                writer.writeEnd();

                return writerError;
            }

#if BITPIT_ENABLE_MPI==1
            // Close the file to let other process write into it
            writerError = writer.writeEnd();
            if (writerError != 0) {
                return writerError;
            }
#endif
        }

        // Write facet data
#if BITPIT_ENABLE_MPI==1
        for (int i = 0; i < nRanks; ++i) {
            if (i == getRank()) {
                // Begin writing
                writerError = writer.writeBegin(STLWriter::WriteAppend, incrementalWrite);
                if (writerError != 0) {
                    writer.writeEnd();
                    return writerError;
                }

#else
        {
#endif
                // Write facet section
                for (long cellId : cells) {
                    // Get cell
                    const Cell &cell = getCell(cellId);

                    // Get vertex coordinates
                    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
                    assert(cellVertexIds.size() == 3);

                    const std::array<double, 3> &coords_0 = getVertex(cellVertexIds[0]).getCoords();
                    const std::array<double, 3> &coords_1 = getVertex(cellVertexIds[1]).getCoords();
                    const std::array<double, 3> &coords_2 = getVertex(cellVertexIds[2]).getCoords();

                    // Evaluate normal
                    const std::array<double, 3> normal = evalFacetNormal(cell.getId());

                    // Write data
                    writerError = writer.writeFacet(coords_0, coords_1, coords_2, normal);
                    if (writerError != 0) {
                        writer.writeEnd();
                        return writerError;
                    }
                }

    #if BITPIT_ENABLE_MPI==1
                // Close the file to let other process write into it
                writerError = writer.writeEnd();
                if (writerError != 0) {
                    return writerError;
                }
            }

            if (isPartitioned()) {
                MPI_Barrier(getCommunicator());
            }
#endif
        }

        // Write footer
        if (isMasterWriter) {
#if BITPIT_ENABLE_MPI==1
            // Begin writing
            //
            // The file was previously closed to let other process write into it.
            writerError = writer.writeBegin(STLWriter::WriteAppend, incrementalWrite);
            if (writerError != 0) {
                writer.writeEnd();
                return writerError;
            }
#endif

            // Write footer section
            writerError = writer.writeFooter(name);
            if (writerError != 0) {
                writer.writeEnd();

                return writerError;
            }

            // End writing
            writerError = writer.writeEnd();
            if (writerError != 0) {
                return writerError;
            }
        }

        // The first PID has been written
        firstPID = false;
    }

    return 0;
}

/*!
 * Import surface tasselation from DGF file.
 *
 * A separate set of vertices will be created for with each facet.
 *
 * DGF facets are added to the present mesh, i.e. current mesh content is not
 * discarded. However, factes are not joined to the cells of the current mesh.
 * After importing the DGF, the resulting mesh can contain duplicate vertices
 * and cells.
 *
 * \param[in] filename name of dgf file
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 *
 * \result on output returns an error flag for I/O error.
*/
int SurfUnstructured::importDGF(const std::string &filename, int PIDOffset, bool PIDSquash)
{
    return importDGF(filename, false, PIDOffset, PIDSquash);
}

/*!
 * Import surface tasselation from DGF file.
 *
 * It is possible to control if DGF factes that share the same vertices will
 * be joined together or if a separate set of vertices will be created for
 * with each facet. When faces are joined together, vertices with a distance
 * less than (10 * machine epsilon) are considered conincident and will be
 * merged together.
 *
 * DGF facets are added to the present mesh, i.e. current mesh content is not
 * discarded. However, factes are not joined to the cells of the current mesh.
 * After importing the DGF, the resulting mesh can contain duplicate vertices
 * and cells.
 * 
 * \param[in] filename name of dgf file
 * \param[in] joinFacets if set to true, facets sharing the same vertices will
 * be joined together, otherwise a separate set of vertices will be created for
 * each facet. In any case, factes will not be joined to the cells of the
 * current mesh
 * \param[in] PIDOffset is the offset for the PID numbering
 * \param[in] PIDSquash controls if the PID of the cells will be read from
 * the file or if the same PID will be assigned to all cells
 * 
 * \result on output returns an error flag for I/O error.
*/
int SurfUnstructured::importDGF(const std::string &filename, bool joinFacets, int PIDOffset, bool PIDSquash)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    DGFObj                                                      dgf_in(filename);
    int                                                         nV = 0, nS = 0;
    std::vector<std::array<double, 3>>                          vertex_list;
    std::vector<std::vector<int>>                               simplex_list;
    std::vector<int>                                            simplex_PID;
    std::vector<long>                                           vertex_map;
    std::vector<long>                                           connect;

    // Counters
    std::vector<std::vector<int>>::iterator                     c_, ce_;
    std::vector<int>::iterator                                  i_, ie_;
    std::vector<long>::iterator                                 j_, je_;

    // ====================================================================== //
    // IMPORT DATA                                                            //
    // ====================================================================== //

    // Read vertices and cells from DGF file
    dgf_in.load(nV, nS, vertex_list, simplex_list, simplex_PID);

    // Add vertices
    vertex_map.resize(nV);
    if (joinFacets) {
        // Vertices added to the patch are also added into a cache. Before
        // inserting a new vertex, the cache is checked: if a vertex with the
        // same coordinates is already in the cache, the new vertex will not
        // be added and the existing one will be used.
        //
        // The cache contains the position of the vertices in the list of
        // vertices read from the DGF file.
        Vertex::Less vertexLess(10 * std::numeric_limits<double>::epsilon());

        auto vertexCoordsLess = [&vertexLess, &vertex_list](const int &i, const int &j)
        {
            return vertexLess(vertex_list[i], vertex_list[j]);
        };

        std::set<int, decltype(vertexCoordsLess)> vertexCache(vertexCoordsLess);

        for (int i = 0; i < nV; ++i) {
            long vertexId;
            auto vertexCacheItr = vertexCache.find(i);
            if (vertexCacheItr == vertexCache.end()) {
                VertexIterator vertexItr = addVertex(vertex_list[i]);
                vertexId = vertexItr.getId();
                vertexCache.insert(i);
            } else {
                vertexId = vertex_map[*vertexCacheItr];
            }
            vertex_map[i] = vertexId;
        }
    } else {
        for (int i = 0; i < nV; ++i) {
            VertexIterator vertexItr = addVertex(vertex_list[i]);
            long vertexId = vertexItr.getId();
            vertex_map[i] = vertexId;
        }
    }

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
 * Export surface tasselation to DGF file
 * 
 * \param[in] filename name of dgf file
 * 
 * \result on output returns an error flag for I/O error
*/
int SurfUnstructured::exportDGF(const std::string &filename)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    DGFObj                                                      dgf_in(filename);
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
