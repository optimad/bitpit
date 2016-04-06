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

#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif

#include "bitpit_SA.hpp"

#include "patch_kernel.hpp"
#include "utils.hpp"

namespace bitpit {


/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class IndexGenerator

	\brief The IndexGenerator class allows to generate unique ids.
*/

/*!
	Creates a new generator.
*/
IndexGenerator::IndexGenerator()
	: m_id(-1)
{

}

/*!
	Generates a unique index.

	If the trash is empty a new index is generated, otherwise an index taken
	from the trash is recycled.

	\return A new unique index.
*/
long IndexGenerator::generateId()
{
	// If the trash is empty generate a new id
	if (m_trash.empty()) {
		assert(m_id < std::numeric_limits<long>::max());
		return ++m_id;
	}

	// If there are ids in the trash recycle te first id in the list
	long id = m_trash.front();
	m_trash.pop_front();

	return id;
}

/*!
	Gets the last assigned id.

	\return The last assigned index.
*/
long IndexGenerator::getLastId()
{
	return m_id;
}

/*!
	Trashes an index.

	A trashed index is an index no more used that can be recycled.

	\param id is the index that will be trashed
*/
void IndexGenerator::trashId(const long &id)
{
	m_trash.push_back(id);
}

/*!
	Reset the generator.
*/
void IndexGenerator::reset()
{
	m_id = -1;
	m_trash.clear();
}

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class PatchKernel

	\brief The PatchKernel class provides an interface for defining patches.

	PatchKernel is the base class for defining patches.
*/

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(const int &id, const int &dimension, bool expert)
	: m_nInternals(0), m_nGhosts(0),
	  m_lastInternalId(Element::NULL_ID),
	  m_firstGhostId(Element::NULL_ID),
	  m_boxFrozen(false), m_boxDirty(true),
	  m_adaptionDirty(true), m_expert(expert), m_hasCustomTolerance(false),
	  m_rank(0), m_nProcessors(1)
#if BITPIT_ENABLE_MPI==1
	  , m_communicator(MPI_COMM_NULL)
#endif
{
	setId(id) ;
	setDimension(dimension);

	std::ostringstream convert;
	convert << getId();
	setName(convert.str());

	// Initialize the geometrical tolerance to a default value
	_setTol(DEFAULT_TOLERANCE);

	// Initializes the bounding box
	setBoundingBoxFrozen(false);
	clearBoundingBox();

	// Set VTK codex
	VTKUnstructuredGrid::setCodex(VTKFormat::APPENDED);

	// Add VTK basic patch data
	VTKUnstructuredGrid::addData("cellIndex", VTKFieldType::SCALAR, VTKLocation::CELL);
	VTKUnstructuredGrid::addData("vertexIndex", VTKFieldType::SCALAR, VTKLocation::POINT);
#if BITPIT_ENABLE_MPI==1
	VTKUnstructuredGrid::addData("rank", VTKFieldType::SCALAR, VTKLocation::CELL);
#endif
}

/*!
	Destroys the patch.
*/
PatchKernel::~PatchKernel()
{
	reset();
}

/*!
	Updates the patch

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> PatchKernel::update(bool trackAdaption)
{
	const std::vector<Adaption::Info> adaptionInfo = updateAdaption(trackAdaption);

	updateBoundingBox();

	return adaptionInfo;
}

/*!
	Updates the adaption

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> PatchKernel::updateAdaption(bool trackAdaption)
{
	std::vector<Adaption::Info> adaptionInfo;
	if (!isAdaptionDirty(true)) {
		return adaptionInfo;
	}

	adaptionInfo = _updateAdaption(trackAdaption);

	m_cells.flush();
	m_interfaces.flush();
	m_vertices.flush();

	setAdaptionDirty(false);

	return adaptionInfo;
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void PatchKernel::markCellForRefinement(const long &id)
{
	bool updated = _markCellForRefinement(id);

	setAdaptionDirty(updated);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void PatchKernel::markCellForCoarsening(const long &id)
{
	bool updated = _markCellForCoarsening(id);

	setAdaptionDirty(updated);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void PatchKernel::enableCellBalancing(const long &id, bool enabled)
{
	bool updated = _enableCellBalancing(id, enabled);

	setAdaptionDirty(updated);
}

/*!
	Resest the patch.
*/
void PatchKernel::reset()
{
	resetVertices();
	resetCells();
	resetInterfaces();
}

/*!
	Resest the vertices of the patch.
*/
void PatchKernel::resetVertices()
{
	m_vertices.clear();
	PiercedVector<Vertex>().swap(m_vertices);
	m_vertexIdGenerator.reset();

	for (auto &cell : m_cells) {
		cell.unsetConnect();
	}
}

/*!
	Resest the cells of the patch.
*/
void PatchKernel::resetCells()
{
	m_cells.clear();
	PiercedVector<Cell>().swap(m_cells);
	m_cellIdGenerator.reset();
	m_nInternals = 0;
	m_nGhosts = 0;

	for (auto &interface : m_interfaces) {
		interface.unsetNeigh();
		interface.unsetOwner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void PatchKernel::resetInterfaces()
{
	m_interfaces.clear();
	PiercedVector<Interface>().swap(m_interfaces);
	m_interfaceIdGenerator.reset();

	for (auto &cell : m_cells) {
		cell.resetInterfaces();
	}
}

/*!
    Reserve memory for vertex storage.

    If the reserve size is smaller than the number of vertices currently stored
    within the patch no action will be taken.

    If instead, the reserve size is greater than the current number of vertices,
    reserve might cause re-location of the internal container into memory,
    potentially invalidating pointers and iterators to vertex entities.

    \param[in] nVertices size of memory reserve (in terms of number of
    vertices).
*/
bool PatchKernel::reserveVertices(size_t nVertices)
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.reserve(nVertices);

	return true;
}

/*!
	Reserve memory for cell storage.

	If the reserve size is smaller than the number of cells currently stored
	within the patch no action will be taken.

	If instead, the reserve size is greater than the current number of cells,
	reserve might cause re-location of the internal container into memory,
	potentially invalidating pointers and iterators to cell entities.

	\param[in] nCells is size of memory reserve (in terms of number of cells).
*/
bool PatchKernel::reserveCells(size_t nCells)
{
	if (!isExpert()) {
		return false;
	}

	m_cells.reserve(nCells);

	return true;
}

/*!
	Reserve memory for interface storage.

	If the reserve size is smaller than the number of interfaces currently
	stored within the patch no action will be taken.

	If instead, the reserve size is greater than the current number of
	interfaces, reserve might cause re-location of the internal container
	into memory, potentially invalidating pointers and iterators to cell
	entities.

	\param[in] nInterfaces is size of memory reserve (in terms of number of
	interfaces).
*/
bool PatchKernel::reserveInterfaces(size_t nInterfaces)
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.reserve(nInterfaces);

	return true;
}

/*!
	Writes the patch to filename specified in input.

	\param filename the filename where the patch will be written to
*/
void PatchKernel::write(std::string filename)
{
	std::string oldFilename = VTKUnstructuredGrid::getName();

	VTKUnstructuredGrid::setName(filename);
	write();
	VTKUnstructuredGrid::setName(oldFilename);
}

/*!
	Writes the patch a filename with the same name of the patch
*/
void PatchKernel::write()
{
	VTKUnstructuredGrid::write();
}

/*!
	Flags the patch for adaption update.

	\param dirty if true, then patch is informed that the patch needs to
	adapt after a refinement, coarsening, ... and thus the current data
	structures are not valid anymore.
*/
void PatchKernel::setAdaptionDirty(bool dirty)
{
	if (m_adaptionDirty == dirty) {
		return;
	}

	m_adaptionDirty = dirty;
}

/*!
	Returns true if the the patch needs to update after an adaption.

	\return This method returns true to indicate the patch needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool PatchKernel::isAdaptionDirty(bool global) const
{
	bool isDirty = m_adaptionDirty;
#if BITPIT_ENABLE_MPI==1
	if (global && isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(&m_adaptionDirty, &isDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif

	return isDirty;
}


/*!
	Returns true if the the patch needs to update its data strucutres.

	\return This method returns true to indicate the patch needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool PatchKernel::isDirty(bool global) const
{
	return (isAdaptionDirty(global) || isBoundingBoxDirty(global));
}

/*!
	Enables or disables expert mode.

	When expert mode is enabled, it will be possible to change the
	patch using low level functions (e.g., it will be possible to
	add individual cells, add vertices, delete cells, ...).

	\param expert if true, the expert mode will be enabled
*/
void PatchKernel::setExpert(bool expert)
{
	if (isExpert() == expert) {
		return;
	}

	m_expert = expert;
}

/*!
	Checks if the expert mode is enabled.

	When expert mode is enabled, it will be possible to change the
	patch using low level functions (e.g., it will be possible to
	add individual cells, add vertices, delete cells, ...).

	\return This method returns true when the expert is enabled,
	otherwise it returns false.
*/
bool PatchKernel::isExpert() const
{
	return m_expert;
}

/*!
	Sets the ID of the patch.

	\param id the ID of the patch
*/
void PatchKernel::setId(int id)
{
	m_id = id;
}

/*!
	Gets the ID of the patch.

	\return The ID of the patch
*/
int PatchKernel::getId() const
{
	return m_id;
}

/*!
	Sets the dimension of the patch.

	\param dimension the dimension of the patch
*/
void PatchKernel::setDimension(int dimension)
{
	m_dimension = dimension;
}

/*!
	Gets the dimension of the patch.

	\return The dimension of the patch
*/
int PatchKernel::getDimension() const
{
	return m_dimension;
}

/*!
	Returns true if the patch is a three-dimensional patch.

	\return This method returns true to indicate the patch is
	three-dimensional
*/
bool PatchKernel::isThreeDimensional() const
{
	return (m_dimension == 3);
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long PatchKernel::getVertexCount() const
{
	return m_vertices.size();
}

/*!
	Gets the nodes owned by the patch.

	\return The nodes owned by the patch.
*/
PiercedVector<Vertex> & PatchKernel::getVertices()
{
	return m_vertices;
}

/*!
	Gets a reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A reference to the vertex with the specified id.
*/
Vertex & PatchKernel::getVertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Gets a constant reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A constant reference to the vertex with the specified id.
*/
const Vertex & PatchKernel::getVertex(const long &id) const
{
	return m_vertices[id];
}

/*!
	Returns an iterator pointing to the specified vertex.

	\result An iterator to the specified vertex.
*/
PatchKernel::VertexIterator PatchKernel::getVertexIterator(const long &id)
{
	return m_vertices.getIterator(id);
}

/*!
	Returns iterator pointing to the first vertex.

	\result An iterator to the first vertex.
*/
PatchKernel::VertexIterator PatchKernel::vertexBegin()
{
	return m_vertices.begin();
}

/*!
	Returns iterator pointing to last vertex.

	\result An iterator to the last vertex.
*/
PatchKernel::VertexIterator PatchKernel::vertexEnd()
{
	return m_vertices.end();
}

/*!
	Generates a new unique id for the vertices.

	\result A new unique id for the vertices.
*/
long PatchKernel::generateVertexId()
{
	if (!isExpert()) {
		return Vertex::NULL_ID;
	}

	return m_vertexIdGenerator.generateId();
}

/*!
	Creates a new vertex with the specified id.

	\param coords are the coordinates of the vertex
	\param id is the id of the new vertex
	\return An iterator pointing to the newly created vertex.
*/
PatchKernel::VertexIterator PatchKernel::createVertex(const std::array<double, 3> &coords, long id)
{
	if (id < 0) {
		id = generateVertexId();
	}

	// Add the vertex
	PiercedVector<Vertex>::iterator iterator = m_vertices.reclaim(id);
    iterator->setId(id);
	iterator->setCoords(coords);

	// Update the bounding box
	addPointToBoundingBox(iterator->getCoords());

	return iterator;
}

/*!
	Adds a new vertex with the specified coordinates.

	\param coords are the coordinates of the vertex
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(const std::array<double, 3> &coords, const long &id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	return createVertex(coords, id);
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(const Vertex &source, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = createVertex(source.getCoords(), id);
	Vertex &vertex = (*iterator);
	id = vertex.getId();
	vertex = source;
	vertex.setId(id);

	return iterator;
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(Vertex &&source, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	if (id < 0) {
		id = source.getId();
	}

	VertexIterator iterator = createVertex(source.getCoords(), std::max(source.getId(), id));
	Vertex &vertex = (*iterator);
	id = vertex.getId();
	vertex = std::move(source);
	vertex.setId(id);

	return iterator;
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
	\param delayed is true a delayed delete will be performed
*/
bool PatchKernel::deleteVertex(const long &id, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

    // Update the bounding box
	removePointFromBoundingBox(m_vertices[id].getCoords(), delayed);

	// Delete the vertex
	m_vertices.erase(id, delayed);
	m_vertexIdGenerator.trashId(id);

    // If there are no more vertices reset them
    if (m_vertices.size() == 0) {
        resetVertices();
    }

	return true;
}

/*!
	Deletes a list of vertices.

	\param ids are the ids of the vertices to be deleted
	\param delayed is true a delayed delete will be performed
*/
bool PatchKernel::deleteVertices(const std::vector<long> &ids, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteVertex(*i, true);
	}

	if (!delayed) {
		m_vertices.flush();
		updateBoundingBox();
	}

	return true;
}

/*!
	Counts free vertices within the patch.

	A free vertex is a vertex on a free face.

	\return The number of free vertices.
*/
long PatchKernel::countFreeVertices() const
{
	std::unordered_set<long> freeVertices;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				continue;
			}

			std::vector<int> faceLocalConnect = cell.getFaceLocalConnect(i);
			for (unsigned int j = 0; j < faceLocalConnect.size(); ++j) {
				freeVertices.insert(cell.getVertex(faceLocalConnect[j]));
			}
		}
	}

        return freeVertices.size();
}

/*!
	Count orphan vertices in the patch.

	An orphan vertex is a vertex not linked by any cells.

	\result The number of orphan vertices.
*/
long PatchKernel::countOrphanVertices() const
{
	std::unordered_set<long> usedVertices;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int i = 0; i < nCellVertices; ++i) {
			usedVertices.insert(cell.getVertex(i));
		}
	}

	return (getVertexCount() - usedVertices.size());
}

/*!
	Find orphan vertices in the patch.

	An orphan vertex is a vertex not linked by any cells.

	\result The list of orphan vertice.
*/
std::vector<long> PatchKernel::findOrphanVertices()
{
	// Add all the vertices to the list
	std::unordered_set<long> vertexSet;
	for (const Vertex &vertex : m_vertices) {
		vertexSet.insert(vertex.getId());
	}

	// Remove used vertices
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int i = 0; i < nCellVertices; ++i) {
			vertexSet.erase(cell.getVertex(i));
		}
	}

	// Build a list
	std::vector<long> vertexList;
	vertexList.reserve(vertexSet.size());
	for (const long &id : vertexSet) {
		vertexList.emplace_back();
		long &lastId = vertexList.back();
		lastId = id;
	}

	return vertexList;
}

/*!
	Remove orphan vertices
*/
bool PatchKernel::deleteOrphanVertices()
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> list = findOrphanVertices();
	deleteVertices(list);

	return true;
}

/*!
	Find and collapse coincident vertices. Cell connectivity is
	automatically updated.

	\param[in] nBins (default = 128) is the number of bins used by
	bin-sorting algorithm to sort tasselation vertices
	\result The list of the of the collapsed vertices.
*/
std::vector<long> PatchKernel::collapseCoincidentVertices(int nBins)
{
	std::vector<long> collapsedVertices;
	if (!isExpert()) {
		return collapsedVertices;
	}

	std::vector<std::vector<std::array<long, 2> > > bins;

	// ====================================================================== //
	// INITIALIZE LOCAL VARIABLES                                             //
	// ====================================================================== //

	// Random number generator
	srand(1223145611);

	// Resize variables
	bins.resize(nBins * nBins * nBins);

	// ====================================================================== //
	// SORT VERTICES ON BINS                                                  //
	// ====================================================================== //

	// Sort vertices
	std::unordered_map<long, long> bin_index = binSortVertex(nBins);

	// Sort cells
	std::array<long, 2> binEntry;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; ++j) {
			binEntry[0] = cell.getId();
			binEntry[1] = j;

			long vertexId = cell.getVertex(j);
			long binId = bin_index[vertexId];
			bins[binId].push_back(binEntry);
		}
	}

	// Free memory
	bin_index.clear();

	// ====================================================================== //
	// COLLAPSE DOUBLE VERTICES                                               //
	// ====================================================================== //
	long collapsedVertexId;
	std::vector<bool> flag(getVertexCount(), false);
	for (auto &bin : bins) {
		int nBinCells = bin.size();
		if (nBinCells > 0) {
			// Randomize vertex insertion
			std::vector<int> list;
			utils::extractWithoutReplacement(nBinCells, nBinCells - 1, list);

			// Vertex insertion
			KdTree<3, Vertex, long> kd(nBinCells);
			for (int j = 0; j < nBinCells; ++j) {
				long cellId = bin[list[j]][0];
				Cell &cell  = m_cells[cellId];

				long k         = bin[list[j]][1];
				long vertexId  = cell.getVertex(k);
				Vertex &vertex = m_vertices[vertexId];

				if (kd.exist(&vertex, collapsedVertexId) >= 0) {
					cell.setVertex(k, collapsedVertexId);
					if (!flag[vertexId]) {
						flag[vertexId] = true;
						collapsedVertices.push_back(vertexId);
					}
				} else {
					flag[vertexId] = true;
					kd.insert(&vertex, vertexId);
				}
			}
		}
	}

	return collapsedVertices;
}

/*!
	Remove coincident vertices from the patch.

	\param[in] nBins (default = 128) is the number of bins used by bin
	sorting algotrithm to sort patch vertices.
*/
bool PatchKernel::deleteCoincidentVertices(int nBins)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> verticesToDelete = collapseCoincidentVertices(nBins);
	deleteVertices(verticesToDelete);

	return true;
}

/*!
	Gets the coordinates of the specified vertex.

	\param id is the id of the vertex
	\result The coordinates of the specified vertex.
*/
const std::array<double, 3> & PatchKernel::getVertexCoords(const long &id) const
{
	return getVertex(id).getCoords();
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long PatchKernel::getCellCount() const
{
	return m_cells.size();
}

/*!
	Gets the number of internal cells in the patch.

	\return The number of internal cells in the patch
*/
long PatchKernel::getInternalCount() const
{
	return m_nInternals;
}

/*!
	Gets the number of ghost cells in the patch.

	\return The number of ghost cells in the patch
*/
long PatchKernel::getGhostCount() const
{
	return m_nGhosts;
}

/*!
	Gets the cells owned by the patch.

	\return The cells owned by the patch.
*/
PiercedVector<Cell> & PatchKernel::getCells()
{
	return m_cells;
}

/*!
	Gets a reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A reference to the cell with the specified id.
*/
Cell & PatchKernel::getCell(const long &id)
{
	return m_cells[id];
}

/*!
	Gets a constant reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A constant reference to the cell with the specified id.
*/
const Cell & PatchKernel::getCell(const long &id) const
{
	return m_cells[id];
}

/*!
	Gets the element type for the cell with the specified id.

	\param id is the id of the requested cell
	\return The element type for the cell with the specified id.
*/
ElementInfo::Type PatchKernel::getCellType(const long &id) const
{
	return m_cells[id].getType();
}

/*!
	Gets a reference to the last internal cell.

	\return A reference to the last internal cell.
*/
Cell & PatchKernel::getLastInternal()
{
	return m_cells[m_lastInternalId];
}

/*!
	Gets a constant reference to the last internal cell.

	\return A constant reference to the last internal cell.
*/
const Cell & PatchKernel::getLastInternal() const
{
	return m_cells[m_lastInternalId];
}

/*!
	Gets a reference to the first ghost cell.

	\return A reference to the first ghost cell.
*/
Cell & PatchKernel::getFirstGhost()
{
	return m_cells[m_firstGhostId];
}

/*!
	Gets a constant reference to the first ghost cell.

	\return A constant reference to the first ghost cell.
*/
const Cell & PatchKernel::getFirstGhost() const
{
	return m_cells[m_firstGhostId];
}

/*!
	Returns an iterator pointing to the specified cell.

	\result An iterator to the specified cell.
*/
PatchKernel::CellIterator PatchKernel::getCellIterator(const long &id)
{
	return m_cells.getIterator(id);
}

/*!
	Returns iterator pointing to the first cell.

	\result An iterator to the first cell.
*/
PatchKernel::CellIterator PatchKernel::cellBegin()
{
	return m_cells.begin();
}

/*!
	Returns iterator pointing to last cell.

	\result An iterator to the last cell.
*/
PatchKernel::CellIterator PatchKernel::cellEnd()
{
	return m_cells.end();
}

/*!
	Returns iterator pointing to the first internal cell.

	\result An iterator to the first internal cell.
*/
PatchKernel::CellIterator PatchKernel::internalBegin()
{
	return m_cells.begin();
}

/*!
	Returns iterator pointing to the end of the list of internal cells.

	\result An iterator to the end of the list of internal cells.
*/
PatchKernel::CellIterator PatchKernel::internalEnd()
{
	return ++m_cells.getIterator(m_lastInternalId);
}

/*!
    Returns iterator to the first ghost cells within the cell list.

    \result An iterator to the first ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostBegin()
{
    return m_cells.getIterator(m_firstGhostId);
}

/*!
	Returns iterator to the end of the list of ghost cells.

	\result An iterator to the end of the list of ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostEnd()
{
	return m_cells.end();
}

/*!
	Generates a new unique id for the cells.

	\result A new unique id for the cells.
*/
long PatchKernel::generateCellId()
{
	if (!isExpert()) {
		return Element::NULL_ID;
	}

	return m_cellIdGenerator.generateId();
}

/*!
	Creates a new cell with the specified id.

	\param type is the type of the cell
	\param id is the id of the new cell
	\param interior is true if the cell is an interior cell, false otherwise
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::createCell(ElementInfo::Type type, bool interior, long id)
{
	if (id < 0) {
		id = generateCellId();
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(type);
	if (cellTypeInfo.dimension > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator;
	if (interior) {
		// Create an internal cell
		//
		// If there are ghosts cells, the internal cell should be inserted
		// before the first ghost cell.
		if (m_firstGhostId < 0) {
			iterator = m_cells.reclaim(id);
		} else {
			iterator = m_cells.reclaimBefore(m_firstGhostId, id);
		}
		m_nInternals++;

		// Update the id of the last internal cell
		if (m_lastInternalId < 0) {
			m_lastInternalId = id;
		} else if (m_cells.rawIndex(m_lastInternalId) < m_cells.rawIndex(id)) {
			m_lastInternalId = id;
		}
	} else {
		// Create a ghost cell
		//
		// If there are internal cells, the ghost cell should be inserted
		// after the last internal cell.
		if (m_lastInternalId < 0) {
			iterator = m_cells.reclaim(id);
		} else {
			iterator = m_cells.reclaimAfter(m_lastInternalId, id);
		}
		m_nGhosts++;

		// Update the id of the first ghost cell
		if (m_firstGhostId < 0) {
			m_firstGhostId = id;
		} else if (m_cells.rawIndex(m_firstGhostId) > m_cells.rawIndex(id)) {
			m_firstGhostId = id;
		}
	}
	iterator->setId(id);

	return iterator;
}

/*!
	Adds a new cell with the specified id.

	\param type is the type of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementInfo::Type type, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	return createCell(type, true, id);
}

/*!
	Adds a new cell with the specified id and type.

	\param type is the type of the cell
	\param interior is true if the cell is the interior of the patch,
	false otherwise
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementInfo::Type type, bool interior, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = createCell(type, interior, id);
	Cell &cell = (*iterator);
	cell.initialize(type, interior);

	return iterator;
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connect is the connectivity of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementInfo::Type type, bool interior,
                                   std::unique_ptr<long[]> &&connect, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = addCell(type, interior, id);
	Cell &cell = (*iterator);
	cell.setConnect(std::move(connect));

	return iterator;
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connect is the connectivity of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementInfo::Type type, bool interior,
								   const std::vector<long> &connect, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}


	// Add the cell
	CellIterator iterator = addCell(type, interior, id);

	// Set the connectivity
	Cell &cell = (*iterator);
	int nCellVertices = cell.getVertexCount();
	std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
	std::copy(connect.data(), connect.data() + nCellVertices, cellConnect.get());
	cell.setConnect(std::move(cellConnect));

	return iterator;
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = createCell(source.getType(), source.isInterior(), id);
	Cell &cell = (*iterator);
	id = cell.getId();
	cell = source;
	cell.setId(id);

	return iterator;
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (id < 0) {
		id = source.getId();
	}

	CellIterator iterator = createCell(source.getType(), source.isInterior(), id);
	Cell &cell = (*iterator);
	id = cell.getId();
	cell = std::move(source);
	cell.setId(id);

	return iterator;
}

/*!
	Deletes a cell.

	\param id is the id of the cell
	\param updateNeighs if true the neighbour data will be updated after
	removing the cell
	\param delayed is true a delayed delete will be performed
*/
bool PatchKernel::deleteCell(const long &id, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	// Update neighbours
	if (updateNeighs) {
		const Cell &cell = m_cells[id];
		int nCellFaces = m_cells[id].getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			// Update adjacency of the neighbours
			int nFaceAdjacencies = cell.getAdjacencyCount(i);
			for (int k = 0; k < nFaceAdjacencies; ++k) {
				long neighId = cell.getAdjacency(i,k);
				if (neighId >= 0) {
					Cell &neigh = m_cells[neighId];

					int neighFace, adjacencyId;
					findFaceNeighCell(neighId, id, neighFace, adjacencyId);
					if (neighFace >= 0) {
						neigh.deleteAdjacency(neighFace, adjacencyId);
					}
				}
			} //next k

			// Update interface
			int nFaceInterfaces = cell.getInterfaceCount(i);
			for (int k = 0; k < nFaceInterfaces; ++k) {
				long interfaceId = cell.getInterface(i,k);
				if (interfaceId >= 0) {
					Interface &interface = m_interfaces[interfaceId];
					if (interface.getOwner() == id) {
							interface.unsetOwner();
					} else {
							interface.unsetNeigh();
					}
				}
			} //next k
		}
	}

	// Delete cell
	bool isInternal = m_cells.at(id).isInterior();
	m_cells.erase(id, delayed);
	m_cellIdGenerator.trashId(id);
	if (isInternal) {
		m_nInternals--;
		if (m_nInternals == 0) {
			m_lastInternalId = Element::NULL_ID;
		} else if (id == m_lastInternalId) {
			m_lastInternalId = m_cells.getSizeMarker(m_nInternals - 1, Element::NULL_ID);
		}
	} else {
		m_nGhosts--;
		if (id == m_firstGhostId) {
			if (m_nGhosts == 0) {
				m_firstGhostId = Element::NULL_ID;
			} else if (m_nInternals == 0) {
				m_firstGhostId = m_cells.getSizeMarker(m_nInternals, Element::NULL_ID);
			} else {
				CellIterator first_ghost_iterator = ++m_cells.getIterator(m_lastInternalId);
				m_firstGhostId = first_ghost_iterator->getId();
			}
		}
	}

	// If there are no more cells reset them
	if (m_cells.size() == 0) {
        resetCells();
    }

	return true;
}

/*!
	Deletes a list of cells.

	\param ids are the ids of the cells to be deleted
	\param updateNeighs if true the neighbour data will be updated after
	removing the cell
	\param delayed is true a delayed delete will be performed
 */
bool PatchKernel::deleteCells(const std::vector<long> &ids, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteCell(*i, updateNeighs, true);
	}

	if (!delayed) {
		m_cells.flush();
	}

	return true;
}

/*!
	Sets the internal flag of a cell.

	\param[in] id is the index of the cell
	\param[in] isInternal is the internal flag that will be set
*/
bool PatchKernel::setCellInternal(const long &id, bool isInternal)
{
	if (!isExpert()) {
		return false;
	}

	if (m_cells[id].isInterior() == isInternal) {
		return true;
	} else if (isInternal) {
		moveGhost2Internal(id);
	} else {
		moveInternal2Ghost(id);
	}

	return true;
}

/*!
	Converts an internal cell to a ghost cell.

	\param[in] id is the index of the cell
*/
PatchKernel::CellIterator PatchKernel::moveInternal2Ghost(const long &id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the element with the last internal cell
	if (id != m_lastInternalId) {
		m_cells.swap(id, m_lastInternalId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.getIterator(id);

	// Update the interior flag
	iterator->setInterior(false);

	// Update cell counters
	--m_nInternals;
	++m_nGhosts;

	// Update the last internal and first ghost markers
	m_firstGhostId = id;
	if (m_nInternals == 0) {
		m_lastInternalId = Element::NULL_ID;
	} else {
		m_lastInternalId = m_cells.getSizeMarker(m_nInternals - 1, Element::NULL_ID);
	}

	// Return the iterator to the new position
	return iterator;
}

/*!
	Converts a ghost cell to an internal cell.

	\param[in] id is the index of the cell
*/
PatchKernel::CellIterator PatchKernel::moveGhost2Internal(const long &id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the cell with the first ghost
	if (id != m_firstGhostId) {
		m_cells.swap(id, m_firstGhostId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.getIterator(id);

	// Update the interior flag
	iterator->setInterior(true);

	// Update cell counters
	++m_nInternals;
	--m_nGhosts;

	// Update the last internal and first ghost markers
	m_lastInternalId = id;
	if (m_nGhosts == 0) {
		m_firstGhostId = Element::NULL_ID;
	} else {
		CellIterator firstGhostIterator = iterator;
		++firstGhostIterator;
		m_firstGhostId = firstGhostIterator->getId();
	}

	// Return the iterator to the new position
	return iterator;
}

/*!
	Counts free cells within the patch.

	A cell is free if contains at least one free face.

	\return The number of free cells.
*/
long PatchKernel::countFreeCells() const
{
	double nFreeCells = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nFreeCells;
				break;
			}
		}
	}

	return nFreeCells;
}

/*!
	Counts orphan cells within the patch.

	A cell is orphan if not adjacent to any cell in the patch (neither
	along an edge, nor at vertex)

	\return The number of orphan cells.
*/
long PatchKernel::countOrphanCells() const
{
	// Compute vertex valence
	std::unordered_map<long, short> vertexValence;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; j++) {
			vertexValence[cell.getVertex(j)] += 1;
		}
	}

	// Loop over cells
	long nOrphanCells = 0;
	for (const Cell &cell : m_cells) {
		long isIsolated = true;
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; j++) {
			long vertexId = cell.getVertex(j);
			if (vertexValence[vertexId] > 1) {
				isIsolated = false;
				break;
			}
		}

		if (isIsolated) {
			++nOrphanCells;
		}
        }

	return nOrphanCells;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\result All the neighbours of the specified cell.
*/
std::vector<long> PatchKernel::findCellNeighs(const long &id) const
{
	return findCellVertexNeighs(id);
}

/*!
	Extracts all the neighbours of the specified cell for the given
	codimension.

	\param id is the id of the cell
	\param codimension the codimension for which the neighbours
	are requested. For a three-dimensional cell a codimension
	equal 1 will extract the face neighbours, a codimension equal
	2 will extract the edge negihbours and a codimension equal
	3 will extract the vertex neighbours. For a two-dimensional
	cell a codimension qual 1 will extract the face neighbours,
	and a codimension equal 2 will extract the vertex neighbours.
	\param complete controls if the list of neighbours should contain
	only the neighbours for the specified codimension, or should contain
	also the neighbours for lower codimensions.
	\result The neighbours for the specified codimension.
*/
std::vector<long> PatchKernel::findCellNeighs(const long &id, int codimension, bool complete) const
{
	assert(codimension >= 1 && codimension <= getDimension());

	if (codimension == 1) {
		return findCellFaceNeighs(id);
	} else if (codimension == getDimension()) {
		return findCellVertexNeighs(id, complete);
	} else if (codimension == 2) {
		return findCellEdgeNeighs(id, complete);
	} else {
		return std::vector<long>();
	}
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\result The neighbours of all the faces of the specified cell.
*/
std::vector<long> PatchKernel::findCellFaceNeighs(const long &id) const
{
	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the face
	// count is evaluated using the ElementInfo associated to the cell.
	int nCellFaces;
	if (m_cells.size() == 0) {
		ElementInfo::Type cellType = getCellType(id);
		const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
		nCellFaces = cellTypeInfo.nFaces;
	} else {
		const Cell &cell = getCell(id);
		nCellFaces = cell.getFaceCount();
	}

	// Get the neighbours
	std::vector<long> neighs;
	for (int i = 0; i < nCellFaces; ++i) {
		std::vector<long> faceNeighs = _findCellFaceNeighs(id, i, std::vector<long>());
		for (auto &neighId : faceNeighs) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> PatchKernel::findCellFaceNeighs(const long &id, const int &face) const
{
	return _findCellFaceNeighs(id, face, std::vector<long>());
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> PatchKernel::_findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getAdjacencyCount(face); ++i) {
		long neighId = cell.getAdjacency(face, i);
		if (neighId < 0) {
			continue;
		}

		if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of all the edges of the specified cell.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified edge, or should
	contain also neighbours that share an entire face
	\result The neighbours of all the edges of the specified cell.
*/
std::vector<long> PatchKernel::findCellEdgeNeighs(const long &id, bool complete) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return std::vector<long>();
	}

	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the edge
	// count is evaluated using the ElementInfo associated to the cell.
	int nCellEdges;
	if (m_cells.size() == 0) {
		ElementInfo::Type cellType = getCellType(id);
		const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
		nCellEdges = cellTypeInfo.nEdges;
	} else {
		const Cell &cell = getCell(id);
		nCellEdges = cell.getEdgeCount();
	}

	// Get the neighbours
	std::vector<long> blackList;
	if (!complete) {
		blackList = findCellFaceNeighs(id);
	}

	std::vector<long> neighs;
	for (int i = 0; i < nCellEdges; ++i) {
		for (auto &neigh : _findCellEdgeNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> PatchKernel::findCellEdgeNeighs(const long &id, const int &edge) const
{
	return _findCellEdgeNeighs(id, edge, std::vector<long>());
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> PatchKernel::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return neighs;
	}

	const Cell &cell = getCell(id);
	std::vector<int> edgeVertices = cell.getEdgeLocalConnect(edge);
	std::size_t nEdgeVertices = edgeVertices.size();
	if (nEdgeVertices < 2) {
		return neighs;
	}

	// The neighbours of the edge are the cells that share all the edge vertices
	std::vector<long> firstVertexNeighs = _findCellVertexNeighs(id, edgeVertices[0], blackList);
	for (std::size_t k = 1; k < nEdgeVertices; ++k) {
		std::vector<long> vertexNeighs = _findCellVertexNeighs(id, edgeVertices[k], blackList);
		for (const long &neighId : vertexNeighs) {
			if (utils::findInOrderedVector<long>(neighId, firstVertexNeighs) != firstVertexNeighs.end()) {
				if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
					utils::addToOrderedVector<long>(neighId, neighs);
				}
			}
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of all the vertices of the specified cell.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the specified cell.
*/
std::vector<long> PatchKernel::findCellVertexNeighs(const long &id, bool complete) const
{
	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the vertex
	// count is evaluated using the ElementInfo associated to the cell.
	int nCellVertices;
	if (m_cells.size() == 0) {
		ElementInfo::Type cellType = getCellType(id);
		const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
		nCellVertices = cellTypeInfo.nVertices;
	} else {
		const Cell &cell = getCell(id);
		nCellVertices = cell.getVertexCount();
	}

	// Get the neighbours
	std::vector<long> blackList;
	if (!complete) {
		if (isThreeDimensional()) {
			blackList = findCellEdgeNeighs(id);
		} else {
			blackList = findCellFaceNeighs(id);
		}
	}

	std::vector<long> neighs;
	for (int i = 0; i < nCellVertices; ++i) {
		for (auto &neigh : _findCellVertexNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> PatchKernel::findCellVertexNeighs(const long &id, const int &vertex) const
{
	return _findCellVertexNeighs(id, vertex, std::vector<long>());
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	Cells that has only a vertex in common are considered neighbours only
	if there are other cells "connecting" them.

	                  .-----.                   .-----.
	                  |     |                   |     |
	                V | A1  |                 V | A2  |
	            .-----+-----.             .-----+-----.
	            |     |                   |     |     |
	            | B1  |                   | B2  | C2  |
	            .-----.                   .-----.-----.

	For example, A1 and B1 are not neighbours (although they share the
	vertex V), whereas A2 and B2 are neighbours.

	\param id is the id of the cell
	\param vertex is a local vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> PatchKernel::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	const Cell &cell = getCell(id);
	long vertexId = cell.getVertex(vertex);

	std::vector<long> neighs;
	std::unordered_set<long> scanQueue;
	std::unordered_set<long> alreadyScan;
	scanQueue.insert(cell.getId());
	while (!scanQueue.empty()) {
		// Pop a cell to process
		long scanId = *(scanQueue.begin());
		const Cell &scanCell = getCell(scanId);

		scanQueue.erase(scanId);
		alreadyScan.insert(scanId);

		// Info on the cell
		const ElementInfo &cellTypeInfo = scanCell.getInfo();
		const std::vector<std::vector<int>> &cellLocalFaceConnect = cellTypeInfo.faceConnect;
		const long *scanCellConnect = scanCell.getConnect();

		// Find the faces that share the vertex
		std::vector<long> faceList;
		for (int i = 0; i < scanCell.getFaceCount(); ++i) {
			// Info on the face
			ElementInfo::Type faceType = scanCell.getFaceType(i);
			const ElementInfo &faceTypeInfo = ElementInfo::getElementInfo(faceType);
			const std::vector<int> &faceLocalConnect = cellLocalFaceConnect[i];

			// Check if the face shares the vertex
			for (int k = 0; k < faceTypeInfo.nVertices; ++k) {
				long faceVertexId = scanCellConnect[faceLocalConnect[k]];
				if (faceVertexId == vertexId) {
					faceList.push_back(i);
					break;
				}
			}
		}

		// If there are no faces that share the vertices go to the next cell
		if (faceList.empty()) {
			continue;
		}

		// Add the current cell to the neighoburs
		if (scanId != id && std::find(blackList.begin(), blackList.end(), scanId) == blackList.end()) {
			utils::addToOrderedVector<long>(scanId, neighs);
		}

		// Add the neighbours of the faces to the scan list
		for (const long &face : faceList) {
			int nFaceNeighs = scanCell.getAdjacencyCount(face);
			for (int k = 0; k < nFaceNeighs; ++k) {
				long neighId = scanCell.getAdjacency(face, k);
				if (neighId >= 0 && alreadyScan.count(neighId) == 0) {
					scanQueue.insert(neighId);
				}
			}
		}
	}

	return neighs;
}

/*!
	Finds the one-ring of the specified vertex of the cell.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The one-ring of the specified vertex of the cell.
*/
std::vector<long> PatchKernel::findCellVertexOneRing(const long &id, const int &vertex) const
{
	std::vector<long> oneRing = findCellVertexNeighs(id, vertex);
	utils::addToOrderedVector<long>(id, oneRing);

	return oneRing;
}

/*!
        Stores the local index of the face shared by cell_idx and neigh_idx
        into face_loc_idx.
        If cell cell_idx and neigh_idx do not share any face, -1 is stored in
        face_loc_idx.

        \param[in] cell_idx cell index
        \param[in] neigh_idx neighbour index
        \param[in,out] face_loc_idx on output stores the local index (on cell
        cell_idx) shared by cell_idx and neigh_idx. If cells cell_idx and neigh_idx
        do not share any face, -1 is stored in face_loc_idx.
        \param[in,out] intf_loc_idx on output stores the index of the adjacency (on face
        face_loc_idx of cell cell_idx). If cells cell_idx and neigh_idx do not share
        any face, -1 is stored into intf_loc_idx.
*/
void PatchKernel::findFaceNeighCell(const long &cell_idx, const long &neigh_idx, int &face_loc_idx, int &intf_loc_idx)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //
    bool                        loop_continue = true;
    int                         n_faces, n_adj;
    int                         j, k;
    Cell                       &cell_ = m_cells[cell_idx];

    // ====================================================================== //
    // LOOP OVER ADJACENCIES                                                  //
    // ====================================================================== //
    n_faces = cell_.getFaceCount();
    j = 0;
    while ( loop_continue && (j < n_faces) ) {
        n_adj = cell_.getAdjacencyCount(j);
        k = 0;
        while ( loop_continue && (k < n_adj) ) {
            loop_continue = ( cell_.getAdjacency( j, k ) != neigh_idx );
            ++k;
        } //next k
        ++j;
    } //next j

    if ( loop_continue) { face_loc_idx = intf_loc_idx = -1; }
    else                { face_loc_idx = --j; intf_loc_idx = --k; }

    return;
}

/*!
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
long PatchKernel::getInterfaceCount() const
{
	return m_interfaces.size();
}

/*!
	Gets the interfaces owned by the patch.

	\return The interfaces owned by the patch.
*/
PiercedVector<Interface> & PatchKernel::getInterfaces()
{
	return m_interfaces;
}

/*!
	Gets a reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A reference to the interface with the specified id.
*/
Interface & PatchKernel::getInterface(const long &id)
{
	return m_interfaces[id];
}

/*!
	Gets a constant reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A constant reference to the interface with the specified id.
*/
const Interface & PatchKernel::getInterface(const long &id) const
{
	return m_interfaces[id];
}

/*!
	Gets the element type for the interface with the specified id.

	\param id is the id of the requested interface
	\return The element type for the interface with the specified id.
*/
ElementInfo::Type PatchKernel::getInterfaceType(const long &id) const
{
	return m_interfaces[id].getType();
}

/*!
	Returns an iterator pointing to the specified interface.

	\result An iterator to the specified interface.
*/
PatchKernel::InterfaceIterator PatchKernel::getInterfaceIterator(const long &id)
{
	return m_interfaces.getIterator(id);
}

/*!
	Returns iterator pointing to the first interface.

	\result An iterator to the first interface.
*/
PatchKernel::InterfaceIterator PatchKernel::interfaceBegin()
{
	return m_interfaces.begin();
}

/*!
	Returns iterator pointing to last interface.

	\result An iterator to the last interface.
*/
PatchKernel::InterfaceIterator PatchKernel::interfaceEnd()
{
	return m_interfaces.end();
}

/*!
 * Generates a new unique id for the interfaces.
 *
 * \result A new unique id for the interfaces.
 */
long PatchKernel::generateInterfaceId()
{
	if (!isExpert()) {
		return Element::NULL_ID;
	}

	return m_interfaceIdGenerator.generateId();
}

/*!
	Creates a new interface with the specified id.

	\param type is the type of the interface
	\param id is the id of the new interface
	\return An iterator pointing to the newly created interface.
*/
PatchKernel::InterfaceIterator PatchKernel::createInterface(ElementInfo::Type type, long id)
{
	if (id < 0) {
		id = generateInterfaceId();
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(type);
	if (cellTypeInfo.dimension > (getDimension() - 1)) {
		return interfaceEnd();
	}

	PiercedVector<Interface>::iterator iterator = m_interfaces.reclaim(id);
    iterator->setId(id);

	return iterator;
}

/*!
	Adds a new interface with the specified id.

	\param type is the type of the interface
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(ElementInfo::Type type, const long &id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = createInterface(type, id);
	Interface &interface = (*iterator);
	interface.initialize(type);

	return iterator;
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(const Interface &source, long id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = createInterface(source.getType(), id);
	Interface &interface = (*iterator);
	id = interface.getId();
	interface = source;
	interface.setId(id);

	return iterator;
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the interface
	\return An iterator pointing to the added interface.

*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(Interface &&source, long id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	if (id < 0) {
		id = source.getId();
	}

	InterfaceIterator iterator = createInterface(source.getType(), id);
	Interface &interface = (*iterator);
	id = interface.getId();
	interface = std::move(source);
	interface.setId(id);

	return iterator;
}

/*!
	Deletes an interface.

	\param id is the id of the interface
	\param updateNeighs if true the neighbour data will be updated after
	removing the interface
	\param delayed is true a delayed delete will be performed
*/
bool PatchKernel::deleteInterface(const long &id, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	// Update neighbours
	if (updateNeighs) {
		Interface &interface = m_interfaces[id];

		// Update owner
		long ownerId = interface.getOwner();
		Cell &owner = m_cells[ownerId];
		int ownerFace = interface.getOwnerFace();

		int ownerInterfaceId = 0;
		while (owner.getInterface(ownerFace, ownerInterfaceId) != id) {
			++ownerInterfaceId;
		}
		owner.deleteInterface(ownerFace, ownerInterfaceId);

		// Update neighbour
		long neighId = interface.getNeigh();
		if (neighId >= 0) {
			Cell &neigh = m_cells[neighId];
			int neighFace = interface.getNeighFace();

			int neighInterfaceId = 0;
			while (neigh.getInterface(neighFace, neighInterfaceId) != id) {
				++neighInterfaceId;
			}
			neigh.deleteInterface(neighFace, neighInterfaceId);
		}
	}

	// Delete interface
	m_interfaces.erase(id, delayed);
	m_interfaceIdGenerator.trashId(id);

    // If there are no more interfaces reset them
    if (m_interfaces.size() == 0) {
        resetInterfaces();
    }

	return true;
}

/*!
	Deletes a list of interfaces.

	\param ids are the ids of the interfaces to be deleted
	\param updateNeighs if true the neighbour data will be updated after
	removing the interface
	\param delayed is true a delayed delete will be performed
*/
bool PatchKernel::deleteInterfaces(const std::vector<long> &ids, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteInterface(*i, updateNeighs, true);
	}

	if (!delayed) {
		m_cells.flush();
	}

	return true;
}

/*!
	Counts free interfaces within the patch.

	An interface is free if belongs to just one cell.

	\result The number of free interfaces.
*/
long PatchKernel::countFreeInterfaces() const
{
	long nFreeInterfaces = 0;
	for (const Interface &interface : m_interfaces) {
		if (interface.getNeigh() < 0) {
			++nFreeInterfaces;
		}
        }

	return nFreeInterfaces;
}

/*!
	Counts orphan interfaces within the patch.

	An interface is orphan if not linked to any cell in the patch.

	\return The number of orphan interfaces.
*/
long PatchKernel::countOrphanInterfaces() const
{
	long nOrphanInterfaces = 0;
	for (const Interface &interface : m_interfaces) {
		if (interface.getOwner() < 0 && interface.getNeigh() < 0) {
			++nOrphanInterfaces;
		}
        }

	return nOrphanInterfaces;
}

/*!
	Count faces within the patch.

	\result The total number of faces in the patch.
*/
long PatchKernel::countFaces() const
{
	double nFaces = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				nFaces += 1. / (cell.getAdjacencyCount(i) + 1);
			} else {
				nFaces += 1.;
			}
		}
	}

	return ((long) round(nFaces));
}

/*!
	Counts free faces within the patch.

	A face is free if a cell has no adjacent along that faces.

	\result The number of free faces.
*/
long PatchKernel::countFreeFaces() const
{
	double nFreeFaces = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nFreeFaces;
			}
		}
	}

	return nFreeFaces;
}

/*!
	Sorts internal vertex storage in ascending id order.
*/
bool PatchKernel::sortVertices()
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.sort();

	return true;
}

/*!
	Sorts internal cell storage in ascending id order.
*/
bool PatchKernel::sortCells()
{
	if (!isExpert()) {
		return false;
	}

	m_cells.sort();

	return true;
}

/*!
	Sorts internal interface storage in ascending id order.
*/
bool PatchKernel::sortInterfaces()
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.sort();

	return true;
}

/*!
	Sorts internal storage for cells, vertices and interfaces in
	ascending id order.
*/
bool PatchKernel::sort()
{
	bool status = sortVertices();
	status |= sortCells();
	status |= sortInterfaces();

	return status;
}

/*!
	Requests the patch to compact the vertex data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the vertex
	data structure can still occupy more memory than it actually needs.
*/
bool PatchKernel::squeezeVertices()
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.squeeze();

	return true;
}

/*!
	Requests the patch to compact the cell data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the cell
	data structure can still occupy more memory than it actually needs.
*/
bool PatchKernel::squeezeCells()
{
	if (!isExpert()) {
		return false;
	}

	m_cells.squeeze();

	return true;
}

/*!
	Requests the patch to compact the interface data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the interface
	data structure can still occupy more memory than it actually needs.
*/
bool PatchKernel::squeezeInterfaces()
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.squeeze();

	return true;
}

/*!
	Requests the patch to compact the data structures and reduce its
	capacity to fit its size.

	The request is non-binding, and after the function call the patch
	can still occupy more memory than it actually needs.
*/
bool PatchKernel::squeeze()
{
	bool status = squeezeVertices();
	status |= squeezeCells();
	status |= squeezeInterfaces();

	return status;
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> PatchKernel::evalCellCentroid(const long &id)
{
	Cell &cell = getCell(id);

	return evalElementCentroid(cell);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> PatchKernel::evalInterfaceCentroid(const long &id)
{
	Interface &interface = getInterface(id);

	return evalElementCentroid(interface);
}

/*!
	Evaluates the centroid of the specified element.

	\param element is the element
	\result The centroid of the specified element.
*/
std::array<double, 3> PatchKernel::evalElementCentroid(const Element &element)
{
	const int nDimensions = 3;

	const long *elementConnect = element.getConnect();
	const ElementInfo &elementInfo = element.getInfo();

	std::array<double, nDimensions> centroid = {{0., 0., 0.}};
	for (int i = 0; i < elementInfo.nVertices; ++i) {
		Vertex &vertex = getVertex(elementConnect[i]);
		const std::array<double, nDimensions> &vertexCoords = vertex.getCoords();
		for (int k = 0; k < nDimensions; ++k) {
			centroid[k] += vertexCoords[k];
		}
	}

	for (int k = 0; k < nDimensions; ++k) {
		centroid[k] /= elementInfo.nVertices;
	}

	return centroid;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool PatchKernel::isPointInside(const double &x, const double &y, const double &z)
{
	return isPointInside({{x, y, z}});
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns the id of the cell the contains the point. If the point
	is not inside the patch, the function returns the id of the null element.
*/
long PatchKernel::locatePoint(const double &x, const double &y, const double &z)
{
	return locatePoint({{x, y, z}});
}

/*!
 * Check whether the i-th face on cell "cell_1" is the same as the j-th face
 * on cell "cell_2".
 * 
 * \param[in] cell_1 global ID of the 1st cell
 * \param[in] i      local index of face to be checked on cell_1
 * \param[in] cell_2 global ID of the 2nd cell
 * \param[in] j      local index of face to be checked on cell_2
 * 
 * \result returns true if face (cell_1, i) and face (cell_2, j) are the same.
*/
bool PatchKernel::isSameFace(
    const long                  &cell_1,
    const int                   &i,
    const long                  &cell_2,
    const int                   &j
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                            check = false;
std::vector<int>                face_loc_connect_A, face_loc_connect_B;
Cell                            *cell_1_ = &m_cells[cell_1], *cell_2_ = &m_cells[cell_2];

// Counters
size_t                          k;

// ========================================================================== //
// CHECK FOR COINCIDENT FACES                                                 //
// ========================================================================== //
face_loc_connect_A = cell_1_->getFaceLocalConnect(i);
face_loc_connect_B = cell_2_->getFaceLocalConnect(j);
if (face_loc_connect_A.size() == face_loc_connect_B.size()) {
    for (k = 0; k < face_loc_connect_A.size(); ++k) {
        face_loc_connect_A[k] = cell_1_->getVertex(face_loc_connect_A[k]);
    } //next i
    for (k = 0; k < face_loc_connect_B.size(); ++k) {
        face_loc_connect_B[k] = cell_2_->getVertex(face_loc_connect_B[k]);
    } //next i
    std::sort(face_loc_connect_A.begin(), face_loc_connect_A.end());
    std::sort(face_loc_connect_B.begin(), face_loc_connect_B.end());
    check = (face_loc_connect_A == face_loc_connect_B);
}

return(check);
    
};

/*!
	Clears the bounding box.

	The box will be cleared also if it declared frozen.
*/
void PatchKernel::clearBoundingBox()
{
	for (int k = 0; k < 3; ++k) {
		m_boxMinPoint[k]   =   std::numeric_limits<double>::max();
		m_boxMinCounter[k] = 0;

		m_boxMaxPoint[k]   = - std::numeric_limits<double>::max();
		m_boxMaxCounter[k] = 0;
	}

	setBoundingBoxDirty(getCellCount() > 0);
}

/*!
	Sets the bounding box.

	The box will be set also if it declared frozen.

	\param minPoint the minimum point of the patch
	\param maxPoint the maximum point of the patch
*/
void PatchKernel::setBoundingBox(const std::array<double, 3> &minPoint, const std::array<double, 3> &maxPoint)
{
	m_boxMinPoint = minPoint;
	m_boxMaxPoint = maxPoint;

	setBoundingBoxDirty(false);

	// Update geometrical tolerance
	if (!isTolCustomized()) {
		resetTol();
	}
}

/*!
	Gets the previously stored patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void PatchKernel::getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	minPoint = m_boxMinPoint;
	maxPoint = m_boxMaxPoint;
}

/*!
	Checks if the bounding box is frozen.

	\result Returns true if the bounding box is frozen, false otherwise.
*/
bool PatchKernel::isBoundingBoxFrozen() const
{
	return m_boxFrozen;
}

/*!
	Sets the bounding box as frozen.

	When the bounding box is frozen it won't be updated on insertion/deletion
	of vertices, neither when the function to update the bounding box is
	called. The only way to change a frozen bounding box is the usage of the
	functions that explicitly sets the bounding box.

	\param frozen controls if the bounding box will be set as frozen
*/
void PatchKernel::setBoundingBoxFrozen(bool frozen)
{
	m_boxFrozen = frozen;
}

/*!
	Checks if the bounding box is dirty.

	\result Returns true if the bounding box is dirty, false otherwise.
*/
bool PatchKernel::isBoundingBoxDirty(bool global) const
{
	bool isDirty = m_boxDirty;
#if BITPIT_ENABLE_MPI==1
	if (global && isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(&m_boxDirty, &isDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif

	return isDirty;
}

/*!
	Sets if the bounding box is dirty.

	\param dirty controls if the bounding box will be set as dirty
*/
void PatchKernel::setBoundingBoxDirty(bool dirty)
{
	m_boxDirty = dirty;
}

/*!
	Updates the stored patch bounding box.
*/
void PatchKernel::updateBoundingBox(bool forcedUpdated)
{
	if (isBoundingBoxFrozen()) {
		return;
	}

	// Check if the bounding box is dirty
	if (!isBoundingBoxDirty() && !forcedUpdated) {
		return;
	}

	// Initialize bounding box
	clearBoundingBox();

	// Unset the dirty flag in order to be able to update the bounding box
	setBoundingBoxDirty(false);

	// Compute bounding box
	for (const auto &vertex : m_vertices) {
		addPointToBoundingBox(vertex.getCoords());
	}
}

/*!
	Update the bounding adding the specified point.

	The bounding box is not updated if it's set as frozen, or if it's in a
	dirty state.

	\param point is the a new point that will be added to the bounding box
*/
void PatchKernel::addPointToBoundingBox(const std::array<double, 3> &point)
{
	if (isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		return;
	}

	bool boxUpdated = false;
	for (size_t k = 0; k < point.size(); ++k) {
		double value = point[k];

		// Update maximum value
		if (value > (m_boxMaxPoint[k] + getTol())) {
			m_boxMaxPoint[k]   = value;
			m_boxMaxCounter[k] = 1;

			boxUpdated = true;
		} else if (std::abs(value - m_boxMaxPoint[k]) <= getTol()) {
			++m_boxMaxCounter[k];
		}

		// Update minimum value
		if (value < (m_boxMinPoint[k] - getTol())) {
			m_boxMinPoint[k]   = value;
			m_boxMinCounter[k] = 1;

			boxUpdated = true;
		} else if (std::abs(value - m_boxMinPoint[k]) <= getTol()) {
			++m_boxMinCounter[k];
		}
	}

	// Update geometrical tolerance
	if (boxUpdated && !isTolCustomized()) {
		resetTol();
	}
}

/*!
	Update the bounding removing the specified point.

	The bounding box is not updated if it's set as frozen, or if it's in a
	dirty state.

	\param point is the point that will be removed from to the bounding box
	\param delayed if true a delayed update ofthe bounding box will
	be performed. This means that, if the bounding box requires an update,
	the update will not be performed and only a flag telling that the
	bounding box needs an update will set.
*/
void PatchKernel::removePointFromBoundingBox(const std::array<double, 3> &point, bool delayed)
{
	if (isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		return;
	}

	for (size_t k = 0; k < point.size(); ++k) {
		double value = point[k];

		// Check if maximum value is still valid
		assert(value <= (m_boxMaxPoint[k] + getTol()));
		if (value > (m_boxMaxPoint[k] - getTol())) {
			setBoundingBoxDirty(true);
		} else if (std::abs(value - m_boxMaxPoint[k]) <= getTol()) {
			--m_boxMaxCounter[k];
			if (m_boxMaxCounter[k] == 0) {
				setBoundingBoxDirty(true);
			}
		}

		// Update minimum value
		assert(value >= (m_boxMinPoint[k] - getTol()));
		if (value < (m_boxMinPoint[k] + getTol())) {
			setBoundingBoxDirty(true);
		} else if (std::abs(value - m_boxMinPoint[k]) <= getTol()) {
			--m_boxMinCounter[k];
			if (m_boxMinCounter[k] == 0) {
				setBoundingBoxDirty(true);
			}
		}
	}

	// Bounding box update
	if (!delayed) {
		updateBoundingBox();
	}
}

/*!
	Sort patch vertices on regular bins.

	\param[in] nBins (default = 128) is the number of bins (on each space
	direction)
	\result Returns the bin index associated to each vertex.
*/
std::unordered_map<long, long> PatchKernel::binSortVertex(int nBins)
{
	// ====================================================================== //
	// VARIABLES DECLARATION                                                  //
	// ====================================================================== //

	// Local variables
	double                              dx, dy, dz;

	// Counters
	long                                i, j, k;
	PiercedVector<Vertex>::iterator     V, E = m_vertices.end();

	// ====================================================================== //
	// ASSOCIATE EACH VERTEX WITH A BIN                                       //
	// ====================================================================== //

	// Update bounding box
	updateBoundingBox();

	// Bin's spacing
	dx = max(1.0e-12, m_boxMaxPoint[0] - m_boxMinPoint[0]) / ((double) nBins);
	dy = max(1.0e-12, m_boxMaxPoint[1] - m_boxMinPoint[1]) / ((double) nBins);
	dz = max(1.0e-12, m_boxMaxPoint[2] - m_boxMinPoint[2]) / ((double) nBins);

	// Loop over vertices
	std::unordered_map<long, long> bin_index;
	for (V = m_vertices.begin(); V != E; ++V) {
		i = std::min(nBins - 1L, long((V->getCoords()[0] - m_boxMinPoint[0]) / dx));
		j = std::min(nBins - 1L, long((V->getCoords()[1] - m_boxMinPoint[1]) / dy));
		k = std::min(nBins - 1L, long((V->getCoords()[2] - m_boxMinPoint[2]) / dz));
		bin_index[V->getId()] = nBins * nBins * k + nBins * j + i;
	}

	return bin_index;
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
*/
void PatchKernel::translate(std::array<double, 3> translation)
{
	// Translate the patch
	for (auto &vertex : m_vertices) {
		vertex.translate(translation);
	}

	// Update the bounding box
	if (!isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		m_boxMinPoint += translation;
		m_boxMaxPoint += translation;
	}
}

/*!
	Translates the patch.

	\param[in] sx translation along x direction
	\param[in] sy translation along y direction
	\param[in] sz translation along z direction
*/
void PatchKernel::translate(double sx, double sy, double sz)
{
	translate({{sx, sy, sz}});
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor vector
*/
void PatchKernel::scale(std::array<double, 3> scaling)
{
	// Scale the patch
	for (auto &vertex : m_vertices) {
		vertex.scale(scaling, m_boxMinPoint);
	}

	// Update the bounding box
	if (!isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		for (int k = 0; k < 3; ++k) {
			m_boxMaxPoint[k] = m_boxMinPoint[k] + scaling[k] * (m_boxMaxPoint[k] - m_boxMinPoint[k]);
		}
	}
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor
*/
void PatchKernel::scale(double scaling)
{
	scale({{scaling, scaling, scaling}});
}

/*!
	Scales the patch.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sz scaling factor along z direction
*/
void PatchKernel::scale(double sx, double sy, double sz)
{
	scale({{sx, sy, sz}});
}

/*!
	Sets the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void PatchKernel::setTol(double tolerance)
{
	_setTol(tolerance);

	m_hasCustomTolerance = true;
}

/*!
	Internal function to set the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void PatchKernel::_setTol(double tolerance)
{
	m_tolerance = tolerance;
}

/*!
	Gets the tolerance for the geometrical checks.

	\result The tolerance fot the geometrical checks.
*/
double PatchKernel::getTol() const
{
	return m_tolerance;
}

/*!
	Resets the tolerance for the geometrical checks.
*/
void PatchKernel::resetTol()
{
	_resetTol();

	m_hasCustomTolerance = false;
}

/*!
	Internal function to reset the tolerance for the geometrical checks.
*/
void PatchKernel::_resetTol()
{
	m_tolerance = 1;
	for (int k = 0; k < 3; ++k) {
		m_tolerance = std::max(m_boxMaxPoint[k] - m_boxMinPoint[k], m_tolerance);
	}
	m_tolerance *= DEFAULT_TOLERANCE;
}

/*!
	Checks if the tolerance for the geometrical checks has been customized
	by the user.

	\result True if the tolerance was customized by the user, false otherwise.
*/
bool PatchKernel::isTolCustomized() const
{
	return m_hasCustomTolerance;
}

/*!
	Extracts the external envelope and appends it to the given patch.

	The external envelope is composed by all the free faces of the patch.

	\param[in,out] envelope is the patch to which the external envelope
	will be appended
*/
void PatchKernel::extractEnvelope(PatchKernel &envelope) const
{

	// ====================================================================== //
	// RESIZE DATA STRUCTURES                                                 //
	// ====================================================================== //
	envelope.reserveVertices(envelope.getVertexCount() + countFreeVertices());
	envelope.reserveCells(envelope.getCellCount() + countFreeFaces());

	// ====================================================================== //
	// LOOP OVER CELLS                                                        //
	// ====================================================================== //
	std::unordered_map<long, long> vertexMap;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				continue;
			}

			// Add face vertices to the envelope and get face
			// connectivity in the envelope
			std::vector<int> faceLocalConnect = cell.getFaceLocalConnect(i);
			int nFaceVertices = faceLocalConnect.size();

			std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);
			for (int j = 0; j < nFaceVertices; ++j) {
				long vertexId = cell.getVertex(faceLocalConnect[j]);

				// If the vertex is not yet in the envelope
				// add it.
				if (vertexMap.count(vertexId) == 0) {
					const Vertex &vertex = getVertex(vertexId);
					VertexIterator envelopeVertex = envelope.addVertex(vertex);
					vertexMap[vertexId] = envelopeVertex->getId();
				}

				// Update face ace connectivity in the envelope
				faceEnvelopeConnect[j] = vertexMap.at(vertexId);
			}

			// Add face to envelope
			ElementInfo::Type faceType = cell.getFaceType(i);
			envelope.addCell(faceType, true, std::move(faceEnvelopeConnect));
		}
	}
}

/*!
	Display patch statistics.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void PatchKernel::displayTopologyStats(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');

	// ====================================================================== //
	// VERTEX STATS                                                           //
	// ====================================================================== //
	out << indent<< "Vertices --------------------------------"     << endl;
	out << indent<< "  # vertices        " << getVertexCount()      << endl;
	out << indent<< "  # orphan vertices " << countOrphanVertices() << endl;
	out << indent<< "  # free vertices   " << countFreeVertices()   << endl;
        //out << indent<< "  # free vertices   " << countDoubleVertices()   << endl;

	// ====================================================================== //
	// FACE STATS                                                             //
	// ====================================================================== //
	out << indent<< "Faces -----------------------------------"     << endl;
	out << indent<< "  # faces           " << countFaces()          << endl;
	out << indent<< "  # free faces      " << countFreeFaces()      << endl;

	// ====================================================================== //
	// CELLS STATS                                                            //
	// ====================================================================== //
	out << indent<< "Cells -----------------------------------"     << endl;
	out << indent<< "  # cells           " << getCellCount()        << endl;
	out << indent<< "  # orphan cells    " << countOrphanCells()    << endl;
	out << indent<< "  # free cells      " << countFreeCells()      << endl;
        //out << indent<< "  # free vertices   " << countDoubleCells()   << endl;
}

/*!
	Display all the vertices currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void PatchKernel::displayVertices(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Vertex &vertex : m_vertices) {
		out << indent << "vertex: " << std::endl;
		vertex.display(out, padding + 2);
	}
}

/*!
	Display all the cells currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void PatchKernel::displayCells(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Cell &cell : m_cells) {
		out << indent << "cell: " << std::endl;
		cell.display(out, padding + 2);
	}
}

/*!
	Display all the interfaces currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void PatchKernel::displayInterfaces(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Interface &interface : m_interfaces) {
		out << indent << "interface: " << std::endl;
		interface.display(out, padding + 2);
	}
}

/*!
 *  Interface method for obtaining field meta Data
 *
 *  @param[in] name is the name of the field to be written
 *  @return Returns a VTKFieldMetaData struct containing the metadata
 *  of the requested custom data.
 */
const VTKFieldMetaData PatchKernel::getMetaData(std::string name)
{
	if (name == "Points") {
		return VTKFieldMetaData(3 * m_vertices.size(), typeid(double));
	} else if (name == "offsets") {
		return VTKFieldMetaData(m_cells.size(), typeid(int));
	} else if (name == "types") {
		return VTKFieldMetaData(m_cells.size(), typeid(VTKElementType));
	} else if (name == "connectivity") {
		long connectSize = 0;
		for (Cell &cell : m_cells) {
			connectSize += cell.getInfo().nVertices;
		}

		return VTKFieldMetaData(connectSize, typeid(long));

	} else if (name == "cellIndex") {
		return VTKFieldMetaData(m_cells.size(), typeid(long));
	} else if (name == "vertexIndex") {
		return VTKFieldMetaData(m_vertices.size(), typeid(long));
#if BITPIT_ENABLE_MPI==1
	} else if (name == "rank") {
		return VTKFieldMetaData(m_cells.size(), typeid(int));
#endif
	}

	// This code should never be reached
	assert(false);
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] format is the format which must be used. Supported options
 *  are "ascii" or "appended". For "appended" type an unformatted binary
 *  stream must be used
 *  @param[in] name is the name of the data to be written. Either user
 *  data or patch data
 */
void PatchKernel::flushData(std::fstream &stream, VTKFormat format, std::string name)
{
	assert(format == VTKFormat::APPENDED);

	static std::unordered_map<long, long> vertexMap;

	if (name == "Points") {
		long vertexId = 0;
		for (Vertex &vertex : m_vertices) {
			vertexMap[vertex.getId()] = vertexId++;

			genericIO::flushBINARY(stream, vertex.getCoords());
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (Cell &cell : m_cells) {
			offset += cell.getInfo().nVertices;
			genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "types") {
		for (Cell &cell : m_cells) {
			VTKElementType VTKType;
			switch (cell.getType())  {

			case ElementInfo::VERTEX:
				VTKType = VTKElementType::VERTEX;
				break;

			case ElementInfo::LINE:
				VTKType = VTKElementType::LINE;
				break;

			case ElementInfo::TRIANGLE:
				VTKType = VTKElementType::TRIANGLE;
				break;

			case ElementInfo::PIXEL:
				VTKType = VTKElementType::PIXEL;
				break;

			case ElementInfo::QUAD:
				VTKType = VTKElementType::QUAD;
				break;

			case ElementInfo::TETRA:
				VTKType = VTKElementType::TETRA;
				break;

			case ElementInfo::VOXEL:
				VTKType = VTKElementType::VOXEL;
				break;

			case ElementInfo::HEXAHEDRON:
				VTKType = VTKElementType::HEXAHEDRON;
				break;

			case ElementInfo::WEDGE:
				VTKType = VTKElementType::WEDGE;
				break;

			case ElementInfo::PYRAMID:
				VTKType = VTKElementType::PYRAMID;
				break;

			default:
				VTKType = VTKElementType::UNDEFINED;
				break;

			}

			genericIO::flushBINARY(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (Cell &cell : m_cells) {
			for (int i = 0; i < cell.getInfo().nVertices; ++i) {
				genericIO::flushBINARY(stream, vertexMap.at(cell.getVertex(i)));
			}
		}

		vertexMap.clear();
		std::unordered_map<long, long>().swap(vertexMap);
	} else if (name == "cellIndex") {
		for (const Cell &cell : m_cells) {
			genericIO::flushBINARY(stream, cell.getId());
		}
	} else if (name == "vertexIndex") {
		for (const Vertex &vertex : m_vertices) {
			genericIO::flushBINARY(stream, vertex.getId());
		}
#if BITPIT_ENABLE_MPI==1
	} else if (name == "rank") {
		std::unordered_map<long, int> ghostRank;
		ghostRank.reserve(getGhostCount());

		for (const auto &neighghostInfo : m_ghost2id) {
			int rank = neighghostInfo.first;
			for (const auto &neighGhost : neighghostInfo.second) {
				ghostRank[neighGhost.second] = rank;
			}
		}

		for (Cell &cell : m_cells) {
			if (cell.isInterior()) {
				genericIO::flushBINARY(stream, m_rank);
			} else {
				genericIO::flushBINARY(stream, ghostRank.at(cell.getId()));
			}
		}
#endif
	}
}

/*!
	@}
*/

}
