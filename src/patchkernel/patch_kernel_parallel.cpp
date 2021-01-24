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
#if BITPIT_ENABLE_MPI==1

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include <mpi.h>
#include <chrono>
#include <unordered_set>

#include "bitpit_communications.hpp"
#include "bitpit_SA.hpp"

#include "patch_kernel.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace chrono;

namespace bitpit {

/*!
	Default cell weight used for patch partitioning
*/
const int PatchKernel::DEFAULT_PARTITIONING_WEIGTH = 1.;

/*!
	Initialize the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void PatchKernel::initializeCommunicator(MPI_Comm communicator)
{
	// Early return if setting a null communicator
	if (communicator == MPI_COMM_NULL) {
		m_communicator = communicator;

		m_nProcessors = 1;
		m_rank        = 0;

		return;
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	MPI_Comm_dup(communicator, &m_communicator);

	// Get MPI information
	MPI_Comm_size(m_communicator, &m_nProcessors);
	MPI_Comm_rank(m_communicator, &m_rank);

	// Set parallel data for the VTK output
	if (m_nProcessors > 1) {
		m_vtk.setParallel(m_nProcessors, m_rank);
	}
}

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void PatchKernel::setCommunicator(MPI_Comm communicator)
{
	// Communication can be set just once
	if (isCommunicatorSet()) {
		throw std::runtime_error ("Patch communicator can be set just once");
	}

	// The communicator has to be valid
	if (communicator == MPI_COMM_NULL) {
		throw std::runtime_error ("Patch communicator is not valid");
	}

	// Set the communicator
	initializeCommunicator(communicator);
}

/*!
	Checks if the communicator to be used for parallel communications has
	already been set.

	\result Returns true if the communicator has been set, false otherwise.
*/
bool PatchKernel::isCommunicatorSet() const
{
	return (getCommunicator() != MPI_COMM_NULL);
}

/*!
	Gets the MPI communicator associated with the patch

	\return The MPI communicator associated with the patch.
*/
const MPI_Comm & PatchKernel::getCommunicator() const
{
	return m_communicator;
}

/*!
	Frees the MPI communicator associated with the patch.
*/
void PatchKernel::freeCommunicator()
{
	if (!isCommunicatorSet()) {
		return;
	}

	int finalizedCalled;
	MPI_Finalized(&finalizedCalled);
	if (finalizedCalled) {
		return;
	}

	MPI_Comm_free(&m_communicator);
}

/*!
	Gets the MPI rank associated with the patch.

	\return The MPI rank associated with the patch.
*/
int PatchKernel::getRank() const
{
	return m_rank;
}

/*!
	Count the MPI processes in the communicator associated with the patch.

	\return The number of MPI processes in the communicator associated with
	the patch.
*/
int PatchKernel::getProcessorCount() const
{
	return m_nProcessors;
}

/*!
	Check if the patch is distributed among different processes.

	\return Return true if the patch is distributed among different processes,
	false otherwise.
*/
bool PatchKernel::isDistributed() const
{
	return (getOwner() < 0);
}

/*!
	If the path is NOT distributed among different processes, returns the
	process that owns the patch, otherwise returns a negative number.

	\return If the path is NOT distributed among different processes, returns
	the process that owns the patch, otherwise returns a negative number.
*/
int PatchKernel::getOwner() const
{
	return m_owner;
}

/*!
	Updates the owner of the patch.

	If the path is NOT distributed among different processes, the owner is set
	to the process that owns the cells, otherwise the owner is set to a negative
	number.
*/
void PatchKernel::updateOwner()
{
	long nInternalCells = getInternalCellCount();
	long nGlobalInternalCells = nInternalCells;
	if (isPartitioned()) {
		MPI_Allreduce(MPI_IN_PLACE, &nGlobalInternalCells, 1, MPI_LONG, MPI_SUM, getCommunicator());
	}

	if (nInternalCells == nGlobalInternalCells) {
		m_owner = getRank();
	} else {
		m_owner = -1;
	}

	if (isPartitioned()) {
		MPI_Allreduce(MPI_IN_PLACE, &m_owner, 1, MPI_INT, MPI_MAX, getCommunicator());
	}
}

/*!
	Initialize the size, expressed in number of layers, of the ghost cells halo.

	No checks will be perfomered on the validity of the halo size.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::initializeHaloSize(std::size_t haloSize)
{
	m_haloSize = haloSize;
}

/*!
	Sets the size, expressed in number of layers, of the ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::setHaloSize(std::size_t haloSize)
{
	if (isPartitioned()) {
		throw std::runtime_error ("Halo size can only be set before partitionig the patch.");
	}

	if (m_haloSize == haloSize) {
		return;
	}

	std::size_t maxHaloSize = _getMaxHaloSize();
	if (haloSize > maxHaloSize) {
		throw std::runtime_error ("Halo size exceeds the maximum allowed value.");
	}

	initializeHaloSize(haloSize);

	_setHaloSize(haloSize);

	if (isPartitioned()) {
		updatePartitioningInfo(true);
	}
}

/*!
	Gets the size, expressed in number of layers, of the ghost cells halo.

	\result The size, expressed in number of layers, of the ghost cells halo.
*/
std::size_t PatchKernel::getHaloSize() const
{
	return m_haloSize;
}

/*!
	Gets the maximum allowed size, expressed in number of layers, of the ghost
	cells halo.

	\result The maximum allowed size, expressed in number of layers, of the
	ghost cells halo.
*/
std::size_t PatchKernel::_getMaxHaloSize()
{
	return 1;
}

/*!
	Internal function to set the size, expressed in number of layers, of the
	ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::_setHaloSize(std::size_t haloSize)
{
	BITPIT_UNUSED(haloSize);
}

/*!
	Converts an internal vertex to a ghost vertex.

	\param[in] id is the index of the vertex
	\param[in] ownerRank is the owner of the vertex
*/
PatchKernel::VertexIterator PatchKernel::internalVertex2GhostVertex(long id, int ownerRank)
{
	if (!isExpert()) {
		return m_vertices.end();
	}

	// Swap the vertex with the last internal vertex
	if (id != m_lastInternalVertexId) {
		m_vertices.swap(id, m_lastInternalVertexId);
	}

	// Get the iterator pointing to the updated position of the vertex
	VertexIterator iterator = m_vertices.find(id);

	// Update the interior flag
	iterator->setInterior(false);

	// Update vertex counters
	--m_nInternalVertices;
	++m_nGhostVertices;

	// Update the last internal and first ghost markers
	m_firstGhostVertexId = id;
	if (m_nInternalVertices == 0) {
		m_lastInternalVertexId = Vertex::NULL_ID;
	} else {
		m_lastInternalVertexId = m_vertices.getSizeMarker(m_nInternalVertices - 1, Vertex::NULL_ID);
	}

	// Set ghost owner
	setGhostVertexOwner(id, ownerRank);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Converts a ghost vertex to an internal vertex.

	\param[in] id is the index of the vertex
*/
PatchKernel::VertexIterator PatchKernel::ghostVertex2InternalVertex(long id)
{
	if (!isExpert()) {
		return m_vertices.end();
	}

	// Swap the vertex with the first ghost
	if (id != m_firstGhostVertexId) {
		m_vertices.swap(id, m_firstGhostVertexId);
	}

	// Get the iterator pointing to the updated position of the vertex
	VertexIterator iterator = m_vertices.find(id);

	// Update the interior flag
	iterator->setInterior(true);

	// Update vertex counters
	++m_nInternalVertices;
	--m_nGhostVertices;

	// Update the last internal and first ghost markers
	m_lastInternalVertexId = id;
	if (m_nGhostVertices == 0) {
		m_firstGhostVertexId = Vertex::NULL_ID;
	} else {
		VertexIterator firstGhostVertexIterator = iterator;
		++firstGhostVertexIterator;
		m_firstGhostVertexId = firstGhostVertexIterator->getId();
	}

	// Unset ghost owner
	unsetGhostVertexOwner(id);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Gets the number of ghost vertices in the patch.

	\return The number of ghost vertices in the patch
*/
long PatchKernel::getGhostVertexCount() const
{
	return m_nGhostVertices;
}

/*!
	Gets a reference to the first ghost vertex.

	\return A reference to the first ghost vertex.
*/
Vertex & PatchKernel::getFirstGhostVertex()
{
	return m_vertices[m_firstGhostVertexId];
}

/*!
	Gets a constant reference to the first ghost vertex.

	\return A constant reference to the first ghost vertex.
*/
const Vertex & PatchKernel::getFirstGhostVertex() const
{
	return m_vertices[m_firstGhostVertexId];
}

/*!
	Restore the vertex with the specified id.

	The kernel should already contain the vertex, only the contents of the
	vertex will be updated.

	\param coords are the coordinates of the vertex
	\param rank is the rank that owns the vertex that will be restored
	\param id is the id of the vertex that will be restored
	\return An iterator pointing to the restored vertex.
*/
PatchKernel::VertexIterator PatchKernel::restoreVertex(const std::array<double, 3> &coords, int rank, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = m_vertices.find(id);
	if (iterator == m_vertices.end()) {
		throw std::runtime_error("Unable to restore the specified vertex: the kernel doesn't contain an entry for that vertex.");
	}

	// There is not need to set the id of the vertex as assigned, because
	// also the index generator will be restored.
	if (rank == getRank()) {
		_restoreInternalVertex(iterator, coords);
	} else {
		_restoreGhostVertex(iterator, coords, rank);
	}

	return iterator;
}

/*!
	Internal function to restore a ghost vertex.

	The kernel should already contain the vertex, only the contents of the
	vertex will be updated.

	\param iterator is an iterator pointing to the vertex to restore
	\param coords are the coordinates of the vertex
	\param rank is the rank that owns the vertex that will be restored
*/
void PatchKernel::_restoreGhostVertex(const VertexIterator &iterator, const std::array<double, 3> &coords, int rank)
{
	// Restore vertex
	Vertex &vertex = *iterator;
	vertex.initialize(iterator.getId(), coords, false);
	m_nGhostVertices++;

	// Update the bounding box
	addPointToBoundingBox(vertex.getCoords());

	// Set owner
	setGhostVertexOwner(vertex.getId(), rank);
}

/*!
	Internal function to delete a ghost vertex.

	\param id is the id of the vertex
*/
void PatchKernel::_deleteGhostVertex(long id)
{
	// Unset ghost owner
	unsetGhostVertexOwner(id);

	// Update the bounding box
	const Vertex &vertex = m_vertices[id];
	removePointFromBoundingBox(vertex.getCoords());

	// Delete vertex
	m_vertices.erase(id, true);
	m_nGhostVertices--;
	if (id == m_firstGhostVertexId) {
		updateFirstGhostVertexId();
	}
}

/*!
    Returns iterator to the first ghost vertex within the vertex list.

    \result An iterator to the first ghost vertex.
*/
PatchKernel::VertexIterator PatchKernel::ghostVertexBegin()
{
	if (m_nGhostVertices > 0) {
		return m_vertices.find(m_firstGhostVertexId);
	} else {
		return m_vertices.end();
	}
}

/*!
	Returns iterator to the end of the list of ghost vertices.

	\result An iterator to the end of the list of ghost vertex.
*/
PatchKernel::VertexIterator PatchKernel::ghostVertexEnd()
{
	return m_vertices.end();
}

/*!
    Returns a constant iterator to the first ghost vertices within the vertex
    list.

    \result A constant iterator to the first ghost vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::ghostVertexConstBegin() const
{
	if (m_nGhostVertices > 0) {
		return m_vertices.find(m_firstGhostVertexId);
	} else {
		return m_vertices.cend();
	}
}

/*!
	Returns a constant iterator to the end of the list of ghost vertices.

	\result A constant iterator to the end of the list of ghost vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::ghostVertexConstEnd() const
{
	return m_vertices.cend();
}

/*!
	Updates the id of the first ghost vertex.
*/
void PatchKernel::updateFirstGhostVertexId()
{
	if (m_nGhostVertices == 0) {
		m_firstGhostVertexId = Vertex::NULL_ID;
	} else if (m_nInternalVertices == 0) {
		VertexIterator firstGhostVertexItr = vertexBegin();
		m_firstGhostVertexId = firstGhostVertexItr->getId();
	} else {
		m_firstGhostVertexId = m_vertices.getSizeMarker(m_nInternalVertices, Vertex::NULL_ID);
	}
}

/*!
	Converts an internal cell to a ghost cell.

	\param[in] id is the index of the cell
	\param[in] ownerRank is the owner of the cell
*/
PatchKernel::CellIterator PatchKernel::internalCell2GhostCell(long id, int ownerRank)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the element with the last internal cell
	if (id != m_lastInternalCellId) {
		m_cells.swap(id, m_lastInternalCellId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.find(id);

	// Update the interior flag
	iterator->setInterior(false);

	// Update cell counters
	--m_nInternalCells;
	++m_nGhostCells;

	// Update the last internal and first ghost markers
	m_firstGhostCellId = id;
	if (m_nInternalCells == 0) {
		m_lastInternalCellId = Cell::NULL_ID;
	} else {
		m_lastInternalCellId = m_cells.getSizeMarker(m_nInternalCells - 1, Cell::NULL_ID);
	}

	// Set ghost owner
	setGhostCellOwner(id, ownerRank);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Converts a ghost cell to an internal cell.

	\param[in] id is the index of the cell
*/
PatchKernel::CellIterator PatchKernel::ghostCell2InternalCell(long id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the cell with the first ghost
	if (id != m_firstGhostCellId) {
		m_cells.swap(id, m_firstGhostCellId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.find(id);

	// Update the interior flag
	iterator->setInterior(true);

	// Update cell counters
	++m_nInternalCells;
	--m_nGhostCells;

	// Update the last internal and first ghost markers
	m_lastInternalCellId = id;
	if (m_nGhostCells == 0) {
		m_firstGhostCellId = Cell::NULL_ID;
	} else {
		CellIterator firstGhostCellIterator = iterator;
		++firstGhostCellIterator;
		m_firstGhostCellId = firstGhostCellIterator->getId();
	}

	// Unset ghost owner
	unsetGhostCellOwner(id);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Gets the number of ghost cells in the patch.

	\return The number of ghost cells in the patch.
*/
long PatchKernel::getGhostCellCount() const
{
	return m_nGhostCells;
}

/*!
	Gets the number of ghost cells in the patch.

	\return The number of ghost cells in the patch.
*/
long PatchKernel::getGhostCount() const
{
	return getGhostCellCount();
}

/*!
	Gets a reference to the first ghost cell.

	\return A reference to the first ghost cell.
*/
Cell & PatchKernel::getFirstGhostCell()
{
	return m_cells[m_firstGhostCellId];
}

/*!
	Gets a reference to the first ghost cell.

	\return A reference to the first ghost cell.
*/
Cell & PatchKernel::getFirstGhost()
{
	return getFirstGhostCell();
}

/*!
	Gets a constant reference to the first ghost cell.

	\return A constant reference to the first ghost cell.
*/
const Cell & PatchKernel::getFirstGhostCell() const
{
	return m_cells[m_firstGhostCellId];
}

/*!
	Gets a constant reference to the first ghost cell.

	\return A constant reference to the first ghost cell.
*/
const Cell & PatchKernel::getFirstGhost() const
{
	return getFirstGhostCell();
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, int rank, long id)
{
	Cell cell = source;
	cell.setId(id);

	return addCell(std::move(cell), rank, id);
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, int rank, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	int connectSize = source.getConnectSize();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	if (!source.hasInfo()){
		std::copy(source.getConnect(), source.getConnect() + connectSize, connectStorage.get());
	}

	CellIterator iterator = addCell(source.getType(), std::move(connectStorage), rank, id);

	Cell &cell = (*iterator);
	id = cell.getId();
	cell = std::move(source);
	cell.setId(id);

	return iterator;
}

/*!
	Adds a new cell with the specified id and type.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param type is the type of the cell
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, int rank, long id)
{
	std::unique_ptr<long[]> connectStorage;
	if (ReferenceElementInfo::hasInfo(type)) {
		int connectSize = ReferenceElementInfo::getInfo(type).nVertices;
		connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	} else {
		connectStorage = std::unique_ptr<long[]>(nullptr);
	}

	return addCell(type, std::move(connectStorage), rank, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param type is the type of the cell
	\param connectivity is the connectivity of the cell
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, const std::vector<long> &connectivity,
											   int rank, long id)
{
	int connectSize = connectivity.size();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	std::copy(connectivity.data(), connectivity.data() + connectSize, connectStorage.get());

	return addCell(type, std::move(connectStorage), rank, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
											   int rank, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (id < 0) {
		id = m_cellIdGenerator.generate();
	} else {
		m_cellIdGenerator.setAssigned(id);
	}

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	CellIterator iterator;
	if (rank == getRank()) {
		iterator = _addInternalCell(type, std::move(connectStorage), id);
	} else {
		iterator = _addGhostCell(type, std::move(connectStorage), rank, id);
	}

	return iterator;
}

/*!
	Internal function to add a ghost cell.

	It is not possible to create a new cell with an id already assigned to an
	existing cell of the patch or with an invalid id. If this happens, an
	exception is thrown. Ids are considered valid if they are greater or equal
	than zero.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::_addGhostCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												 int rank, long id)
{
	// Create the cell
	//
	// If there are internal cells, the ghost cell should be inserted
	// after the last internal cell.
	CellIterator iterator;
	if (m_lastInternalCellId < 0) {
		iterator = m_cells.emreclaim(id, id, type, std::move(connectStorage), false, true);
	} else {
		iterator = m_cells.emreclaimAfter(m_lastInternalCellId, id, id, type, std::move(connectStorage), false, true);
	}
	m_nGhostCells++;

	// Update the id of the first ghost cell
	if (m_firstGhostCellId < 0) {
		m_firstGhostCellId = id;
	} else if (m_cells.rawIndex(m_firstGhostCellId) > m_cells.rawIndex(id)) {
		m_firstGhostCellId = id;
	}

	// Set owner
	setGhostCellOwner(id, rank);

	// Set the alteration flags of the cell
	setAddedCellAlterationFlags(id);

	return iterator;
}

/*!
	Restore the cell with the specified id.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be restored
	\param id is the id of the cell that will be restored
	\return An iterator pointing to the restored cell.
*/
PatchKernel::CellIterator PatchKernel::restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												   int rank, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	CellIterator iterator = m_cells.find(id);
	if (iterator == m_cells.end()) {
		throw std::runtime_error("Unable to restore the specified cell: the kernel doesn't contain an entry for that cell.");
	}

	// There is not need to set the id of the cell as assigned, because
	// also the index generator will be restored.
	if (rank == getRank()) {
		_restoreInternalCell(iterator, type, std::move(connectStorage));
	} else {
		_restoreGhostCell(iterator, type, std::move(connectStorage), rank);
	}

	return iterator;
}

/*!
	Internal function to restore a ghost cell.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param iterator is an iterator pointing to the cell to restore
	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be restored
*/
void PatchKernel::_restoreGhostCell(const CellIterator &iterator, ElementType type,
								std::unique_ptr<long[]> &&connectStorage, int rank)
{
	// Restore cell
	long cellId = iterator.getId();
	Cell &cell = *iterator;
	cell.initialize(iterator.getId(), type, std::move(connectStorage), false, true);
	m_nGhostCells++;

	// Set owner
	setGhostCellOwner(cellId, rank);

	// Set the alteration flags of the cell
	setRestoredCellAlterationFlags(cellId);
}

/*!
	Internal function to delete a ghost cell.

	\param id is the id of the cell
*/
void PatchKernel::_deleteGhostCell(long id)
{
	// Unset ghost owner
	unsetGhostCellOwner(id);

	// Set the alteration flags of the cell
	setDeletedCellAlterationFlags(id);

	// Delete cell
	m_cells.erase(id, true);
	m_nGhostCells--;
	if (id == m_firstGhostCellId) {
		updateFirstGhostCellId();
	}
}

/*!
    Returns iterator to the first ghost cells within the cell list.

    \result An iterator to the first ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostCellBegin()
{
	if (m_nGhostCells > 0) {
		return m_cells.find(m_firstGhostCellId);
	} else {
		return m_cells.end();
	}
}

/*!
    Returns iterator to the first ghost cells within the cell list.

    \result An iterator to the first ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostBegin()
{
	return ghostCellBegin();
}

/*!
	Returns iterator to the end of the list of ghost cells.

	\result An iterator to the end of the list of ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostCellEnd()
{
	return m_cells.end();
}

/*!
	Returns iterator to the end of the list of ghost cells.

	\result An iterator to the end of the list of ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostEnd()
{
	return ghostCellEnd();
}

/*!
    Returns a constant iterator to the first ghost cells within the cell list.

    \result A constant iterator to the first ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostCellConstBegin() const
{
	if (m_nGhostCells > 0) {
		return m_cells.find(m_firstGhostCellId);
	} else {
		return m_cells.cend();
	}
}

/*!
    Returns a constant iterator to the first ghost cells within the cell list.

    \result A constant iterator to the first ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostConstBegin() const
{
	return ghostCellConstBegin();
}

/*!
	Returns a constant iterator to the end of the list of ghost cells.

	\result A constant iterator to the end of the list of ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostCellConstEnd() const
{
	return m_cells.cend();
}

/*!
	Returns a constant iterator to the end of the list of ghost cells.

	\result A constant iterator to the end of the list of ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostConstEnd() const
{
	return ghostCellConstEnd();
}

/*!
	Updates the id of the first ghost cell.
*/
void PatchKernel::updateFirstGhostCellId()
{
	if (m_nGhostCells == 0) {
		m_firstGhostCellId = Cell::NULL_ID;
	} else if (m_nInternalCells == 0) {
		CellIterator firstghostCellItr = cellBegin();
		m_firstGhostCellId = firstghostCellItr->getId();
	} else {
		m_firstGhostCellId = m_cells.getSizeMarker(m_nInternalCells, Cell::NULL_ID);
	}
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(cellRanks, trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, bool squeezeStorage)
{
	partitioningPrepare(cellRanks, false);

	std::vector<adaption::Info> partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	std::unordered_map<long, double> dummyCellWeights;

	return partition(dummyCellWeights, trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(bool trackPartitioning, bool squeezeStorage)
{
	std::unordered_map<long, double> dummyCellWeights;

	partitioningPrepare(dummyCellWeights, false);

	std::vector<adaption::Info> partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	\param communicator is the communicator that will be used
	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, a weight equal to one is used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(cellWeights, trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, a weight equal to one is used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, bool squeezeStorage)
{
	partitioningPrepare(cellWeights, false);

	std::vector<adaption::Info> partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(cellRanks, trackPartitioning);
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning)
{
	// Patch needs to support partitioning
	if (!isPartitioningSupported()) {
		throw std::runtime_error ("The patch does not support partitioning!");
	}

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Fill partitioning ranks
	int patchRank = getRank();

	std::set<int> recvRanks;
	m_partitioningOutgoings.clear();
	for (auto &entry : cellRanks) {
		int recvRank = entry.second;
		if (recvRank == patchRank) {
			continue;
		}

		long cellId = entry.first;
		if (m_ghostCellOwners.count(cellId) > 0) {
			continue;
		}

		m_partitioningOutgoings.insert(entry);
		recvRanks.insert(recvRank);
	}

	// Identify exchange entries
	int nRanks = getProcessorCount();

	int nExchanges = recvRanks.size();
	std::vector<std::pair<int, int>> exchanges;
	exchanges.reserve(nExchanges);
	for (int recvRank : recvRanks) {
		exchanges.emplace_back(patchRank, recvRank);
	}

	int exchangesGatherCount = 2 * nExchanges;
	std::vector<int> exchangeGatherCount(nRanks);
	MPI_Allgather(&exchangesGatherCount, 1, MPI_INT, exchangeGatherCount.data(), 1, MPI_INT, m_communicator);

	std::vector<int> exchangesGatherDispls(nRanks, 0);
	for (int i = 1; i < nRanks; ++i) {
		exchangesGatherDispls[i] = exchangesGatherDispls[i - 1] + exchangeGatherCount[i - 1];
	}

	int nGlobalExchanges = nExchanges;
	MPI_Allreduce(MPI_IN_PLACE, &nGlobalExchanges, 1, MPI_INT, MPI_SUM, m_communicator);

	m_partitioningGlobalExchanges.resize(nGlobalExchanges);
	MPI_Allgatherv(exchanges.data(), exchangesGatherCount, MPI_INT, m_partitioningGlobalExchanges.data(),
	               exchangeGatherCount.data(), exchangesGatherDispls.data(), MPI_INT, m_communicator);

	std::sort(m_partitioningGlobalExchanges.begin(), m_partitioningGlobalExchanges.end(), greater<std::pair<int,int>>());

	// Get global list of ranks that will send data
	std::unordered_set<int> globalSendRanks;
	for (const std::pair<int, int> &entry : m_partitioningGlobalExchanges) {
		globalSendRanks.insert(entry.first);
	}

	// Get global list of ranks that will receive data
	std::unordered_set<int> globalRecvRanks;
	for (const std::pair<int, int> &entry : m_partitioningGlobalExchanges) {
		globalRecvRanks.insert(entry.second);
	}

	// Identify if this is a serialization or a normal partitioning
	//
	// We are serializing the patch if all the processes are sending all
	// their cells to the same rank.
	m_partitioningSerialization = (globalRecvRanks.size() == 1);
	if (m_partitioningSerialization) {
		int receiverRank = *(globalRecvRanks.begin());
		if (patchRank != receiverRank) {
			if (m_partitioningOutgoings.size() != (std::size_t) getInternalCellCount()) {
				m_partitioningSerialization = false;
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, &m_partitioningSerialization, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
	}

	// Build the information on the cells that will be sent
	//
	// Only internal cells are tracked.
	std::vector<adaption::Info> partitioningData;
	if (trackPartitioning) {
		for (const std::pair<int, int> &entry : m_partitioningGlobalExchanges) {
			int sendRank = entry.first;
			if (sendRank != patchRank) {
				continue;
			}

			int recvRank = entry.second;

			std::vector<long> previous;
			for (const auto &entry : m_partitioningOutgoings) {
				int cellRank = entry.second;
				if (cellRank != recvRank) {
					continue;
				}

				long cellId = entry.first;
				previous.push_back(cellId);
			}

			if (!previous.empty()) {
				partitioningData.emplace_back();
				adaption::Info &partitioningInfo = partitioningData.back();
				partitioningInfo.entity   = adaption::ENTITY_CELL;
				partitioningInfo.type     = adaption::TYPE_PARTITION_SEND;
				partitioningInfo.rank     = recvRank;
				partitioningInfo.previous = std::move(previous);
			}
		}
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	std::unordered_map<long, double> dummyCellWeights;

	return partitioningPrepare(dummyCellWeights, trackPartitioning);
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(bool trackPartitioning)
{
	std::unordered_map<long, double> dummyCellWeights;

	return partitioningPrepare(dummyCellWeights, trackPartitioning);
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param communicator is the communicator that will be used
	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, a weight equal to one is used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(cellWeights, trackPartitioning);
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, a weight equal to one is used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(const std::unordered_map<long, double> &cellWeights, bool trackPartitioning)
{
	// Patch needs to support partitioning
	if (!isPartitioningSupported()) {
		throw std::runtime_error ("The patch does not support partitioning!");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Execute the partitioning preparation
	std::vector<adaption::Info> partitioningData = _partitioningPrepare(cellWeights, DEFAULT_PARTITIONING_WEIGTH, trackPartitioning);

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Alter the patch performing the partitioning.

	The actual modification of the patch takes place during this phase. After
	this phase the adaption is completed and the patch is in its final state.
	Optionally the patch can track the changes performed to the patch.

	\param trackPartitioning if set to true the function will return the changes
	done to the patch during the partitioning
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the partitioning
	\result If the partitioning is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the partitioning, otherwise
	an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;
	if (!isPartitioningSupported()) {
		return partitioningData;
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_CLEAN) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_PREPARED) {
		throw std::runtime_error ("The prepare function has no been called.");
	}

	// Partition the patch
	partitioningData = _partitioningAlter(trackPartitioning);

	// Finalize patch alterations
	finalizeAlterations(squeezeStorage);

	// Update the status
	setPartitioningStatus(PARTITIONING_ALTERED);

	return partitioningData;
}

/*!
	Cleanup patch data structured after the partitioning.

	The patch will only clean-up the data structures needed during the
	partitioning.
*/
void PatchKernel::partitioningCleanup()
{
	if (!isPartitioningSupported()) {
		return;
	}

	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_CLEAN) {
		return;
	} else if (partitioningStatus == PARTITIONING_PREPARED) {
		throw std::runtime_error ("It is not yet possible to abort a partitioning.");
	} else if (partitioningStatus != PARTITIONING_ALTERED) {
		throw std::runtime_error ("The alter function has no been called.");
	}

	// Clean-up the partitioning
	_partitioningCleanup();

	if (m_partitioningOutgoings.size() > 0) {
		std::unordered_map<long, int>().swap(m_partitioningOutgoings);
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_CLEAN);
}

/*!
	Checks if the patch has been partitioned.

	\result Returns true if the patch has been partitioned, false otherwise.
*/
bool PatchKernel::isPartitioned() const
{
	return isCommunicatorSet();
}

/*!
	Checks if the patch supports partitioning.

	\return Returns true if the patch supports partitioning, false otherwise.
*/
bool PatchKernel::isPartitioningSupported() const
{
    return (getPartitioningStatus() != PARTITIONING_UNSUPPORTED);
}

/*!
	Returns the current partitioning status.

	\param global if set to true the partitioning status will be
	\return The current partitioning status.
*/
PatchKernel::PartitioningStatus PatchKernel::getPartitioningStatus(bool global) const
{
	int partitioningStatus = static_cast<int>(m_partitioningStatus);

	if (global && isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &partitioningStatus, 1, MPI_INT, MPI_MAX, communicator);
	}

	return static_cast<PartitioningStatus>(partitioningStatus);
}

/*!
	Set the current partitioning status.

	\param status is the partitioning status that will be set
*/
void PatchKernel::setPartitioningStatus(PartitioningStatus status)
{
	m_partitioningStatus = status;
}


/*!
	Evaluate partitioning load unbalance index.

	This index measures the performance lost to imbalanced load or, conversely,
	the performance that could be reclaimed by balancing the load.

	Unbalance index is evaluate almost as in Equation 1 of paper "Quantifying
	the Effectiveness of Load Balance Algorithms", Olga Pearce, Todd Gamblin†,
	Bronis R. de Supinski†, Martin Schulz†, Nancy M. Amato, Department of
	Computer Science and Engineering, Texas A&M University, College Station,
	TX, USA. The difference is that the imbalance factor evaluate by this
	function is not a percentage, i.e., it is not multiplied by 100.

	\result Partitioning load unbalance index.
*/
double PatchKernel::evalPartitioningUnbalance() const
{
	std::unordered_map<long, double> dummyCellWeights;

	return evalPartitioningUnbalance(dummyCellWeights);
}

/*!
	Evaluate partitioning load unbalance index.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, a weight equal to one is used.
	\result Partitioning load unbalance index.
*/
double PatchKernel::evalPartitioningUnbalance(const std::unordered_map<long, double> &cellWeights) const
{
	if (!isPartitioned()) {
		return 0.;
	}

	// Evaluate partition weight
	double partitionWeight;
	if (!cellWeights.empty()) {
		partitionWeight = 0.;
		for (auto cellItr = internalCellConstBegin(); cellItr != internalCellConstEnd(); ++cellItr) {
			long cellId = cellItr.getId();

			double cellWeight;
			auto weightItr = cellWeights.find(cellId);
			if (weightItr != cellWeights.end()) {
				cellWeight = weightItr->second;
			} else {
				cellWeight = DEFAULT_PARTITIONING_WEIGTH;
			}

			partitionWeight += cellWeight;
		}
	} else {
		partitionWeight = DEFAULT_PARTITIONING_WEIGTH * getInternalCellCount();
	}

	// Evalaute global weights
	double totalPartitionWeight;
	MPI_Allreduce(&partitionWeight, &totalPartitionWeight, 1, MPI_DOUBLE, MPI_SUM, getCommunicator());

	double maximumPartitionWeight;
	MPI_Allreduce(&partitionWeight, &maximumPartitionWeight, 1, MPI_DOUBLE, MPI_MAX, getCommunicator());

	double meanPartitionWeight = totalPartitionWeight / getProcessorCount();

	// Evaluate the unbalance
	double unbalance = (maximumPartitionWeight / meanPartitionWeight - 1.);

	return unbalance;
}

/*!
	Prepares the patch for performing the partitioning.

	Default implementation is a no-op function and can only be used to
	partition an empty patch.

	\param cellWeights are the weights of the cells, the weight represents the
	relative computational cost associated with a specified cell. If no weight
	is specified for a cell, the default weight will be used
	\param defaultWeight is the default weight that will assigned to the cells
	for which an explicit weight has not been defined
	\param trackPartitioning if set to true the function will return the
	changes that will be performed in the alter step
	\result If the partitioning is tracked, returns a vector of adaption::Info
	that can be used to discover what changes will be performed in the alter
	step, otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningPrepare(const std::unordered_map<long, double> &cellWeights, double defaultWeight, bool trackPartitioning)
{
	BITPIT_UNUSED(cellWeights);
	BITPIT_UNUSED(defaultWeight);
	BITPIT_UNUSED(trackPartitioning);

	if (m_nInternalCells > 0) {
		throw std::runtime_error ("For this patch automatic partitioning can be used only to initialize partitioning on an empty patch.");
	}

	return std::vector<adaption::Info>();
}

/*!
    Alter the patch performing the partitioning.

    \param trackPartitioning if set to true the function will return the changes
    done to the patch during the partitioning
    \result If the partitioning is tracked, returns a vector of adaption::Info
    with all the changes done to the patch during the partitioning, otherwise
    an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter(bool trackPartitioning)
{
    // Adjacencies need to be up-to-date
    updateAdjacencies();

    // Ghost final ghost owners
    //
    // If we are serializing the patch, we need to delete all the ghosts cells.
    //
    // If we are not serializing the patch, some of the cells that will be send
    // during partitioning may be ghosts on other processes. Therefore, we
    // need to identify the owners of the ghosts after the partitioning.
    std::unordered_map<long, int> ghostCellOwnershipChanges;
    if (!m_partitioningSerialization) {
        if (isPartitioned()) {
            ghostCellOwnershipChanges = _partitioningAlter_evalGhostCellOwnershipChanges();
        }
    } else {
        _partitioningAlter_deleteGhosts();
    }

    // Communicate patch data structures
    int patchRank = getRank();

    m_partitioningCellsTag = communications::tags().generate(m_communicator);
    m_partitioningVerticesTag = communications::tags().generate(m_communicator);

    std::unordered_set<int> batchSendRanks;
    std::unordered_set<int> batchRecvRanks;

    std::vector<adaption::Info> partitioningData;
    std::vector<adaption::Info> rankPartitioningData;
    while (!m_partitioningGlobalExchanges.empty()) {
        // Add exchanges to the batch
        //
        // Inside a batch we can mix send and receives, but a process can be
        // either a sender or a receiver, it cannot be both.
        batchSendRanks.clear();
        batchRecvRanks.clear();
        while (!m_partitioningGlobalExchanges.empty()) {
            const std::pair<int, int> &entry = m_partitioningGlobalExchanges.back();

            int entrySendRank = entry.first;
            if (batchRecvRanks.count(entrySendRank)) {
                break;
            }

            int entryRecvRank = entry.second;
            if (batchSendRanks.count(entryRecvRank)) {
                break;
            }

            batchSendRanks.insert(entrySendRank);
            batchRecvRanks.insert(entryRecvRank);
            m_partitioningGlobalExchanges.pop_back();
        }

        // Detect the role of the current process
        bool isSender   = (batchSendRanks.count(patchRank) > 0);
        bool isReceiver = (batchRecvRanks.count(patchRank) > 0);

        // Perform sends/receives
        if (isSender) {
            rankPartitioningData = _partitioningAlter_sendCells(batchRecvRanks, trackPartitioning, &ghostCellOwnershipChanges);
        } else if (isReceiver) {
            rankPartitioningData = _partitioningAlter_receiveCells(batchSendRanks, trackPartitioning, &ghostCellOwnershipChanges);
        }

        if (trackPartitioning) {
            for (adaption::Info &rankPartitioningInfo : rankPartitioningData) {
                partitioningData.emplace_back(std::move(rankPartitioningInfo));
            }
        }

        // Update ghost ownership
        for (int sendRank : batchSendRanks) {
            _partitioningAlter_applyGhostCellOwnershipChanges(sendRank, &ghostCellOwnershipChanges);
        }
    }

    assert(ghostCellOwnershipChanges.size() == 0);

    communications::tags().trash(m_partitioningCellsTag, m_communicator);
    communications::tags().trash(m_partitioningVerticesTag, m_communicator);

    return partitioningData;
}

/*!
    Get the ghost that will change ownership after partitioning.

    Some of the cells that are send during a partitioning may be ghosts on
    other processes. We need to find out the final ghost owners after the
    partitioning.

    This function uses the ghost exchange information.

    \result The ghosts that will change ownership after the partitioning.
*/
std::unordered_map<long, int> PatchKernel::_partitioningAlter_evalGhostCellOwnershipChanges()
{
    int patchRank = getRank();

    // Initialize communications
    //
    // The communications will exchange the ranks on ghosts cells.
    DataCommunicator notificationCommunicator(getCommunicator());
    size_t notificationDataSize = sizeof(patchRank);

    // Set and start the receives
    for (const auto &entry : getGhostCellExchangeTargets()) {
        const int rank = entry.first;
        const auto &targetCells = entry.second;

        notificationCommunicator.setRecv(rank, targetCells.size() * notificationDataSize);
        notificationCommunicator.startRecv(rank);
    }

    // Set and start the sends
    //
    // Data buffer is filled with the ranks that will own source cells after
    // the partitioning.
    for (const auto &entry : getGhostCellExchangeSources()) {
        const int rank = entry.first;
        auto &sourceCells = entry.second;

        notificationCommunicator.setSend(rank, sourceCells.size() * notificationDataSize);
        SendBuffer &buffer = notificationCommunicator.getSendBuffer(rank);
        for (long id : sourceCells) {
            int finalOwner;
            if (m_partitioningOutgoings.count(id) > 0) {
                finalOwner = m_partitioningOutgoings.at(id);
            } else {
                finalOwner = patchRank;
            }

            buffer << finalOwner;
        }
        notificationCommunicator.startSend(rank);
    }

    // Receive the final owners of of the ghosts
    //
    // Data buffer contains the ranks that will own ghost cells after the
    // partitioning.
    std::unordered_map<long, int> ghostCellOwnershipChanges;

    int nCompletedRecvs = 0;
    while (nCompletedRecvs < notificationCommunicator.getRecvCount()) {
        int rank = notificationCommunicator.waitAnyRecv();
        const auto &list = getGhostCellExchangeTargets(rank);

        RecvBuffer &buffer = notificationCommunicator.getRecvBuffer(rank);
        for (long id : list) {
            int finalOwner;
            buffer >> finalOwner;

            if (finalOwner != m_ghostCellOwners.at(id)) {
                ghostCellOwnershipChanges[id] = finalOwner;
            }
        }

        ++nCompletedRecvs;
    }

    // Wait for the sends to finish
    notificationCommunicator.waitAllSends();

    return ghostCellOwnershipChanges;
}

/*!
    Apply changes of ghost ownership for ghosts previously owned by the
    specified sender process.

    \param sendRank is the rank of the process sending the cells
    \param[in,out] ghostCellOwnershipChanges are the ghosts that will change
    ownership after the partitioning, on output the list will be updated
    removing the ghosts that are no longer on the partition (i.e., ghosts
    deleted or promoted to internal cell) and the ghosts whose ownership
    have been updated
*/
void PatchKernel::_partitioningAlter_applyGhostCellOwnershipChanges(int sendRank,
                                                                std::unordered_map<long, int> *ghostCellOwnershipChanges)
{
    for (auto itr = ghostCellOwnershipChanges->begin(); itr != ghostCellOwnershipChanges->end();) {
        // Consider only ghosts previously owned by the sender
        long ghostCellId = itr->first;
        int previousGhostCellOwner = m_ghostCellOwners.at(ghostCellId);
        if (previousGhostCellOwner != sendRank) {
            ++itr;
            continue;
        }

        // Update ghost owner
        int finalGhostCellOwner = itr->second;
        setGhostCellOwner(ghostCellId, finalGhostCellOwner);
        itr = ghostCellOwnershipChanges->erase(itr);
    }
}

/*!
    Delete ghosts.
*/
void PatchKernel::_partitioningAlter_deleteGhosts()
{
    // Delete ghost cells
    std::vector<long> cellsDeleteList;
    cellsDeleteList.reserve(m_ghostCellOwners.size());
    for (const auto &entry : m_ghostCellOwners) {
        long cellId = entry.first;
        cellsDeleteList.emplace_back(cellId);
    }

    deleteCells(cellsDeleteList);

    // Prune stale adjacencies
    pruneStaleAdjacencies();

    // Prune stale interfaces
    //
    // Interfaces will be re-created after the partitioning.
    pruneStaleInterfaces();

    // Delete vertices no longer used
    deleteOrphanVertices();

    // Convert all ghost vertices to internal vertices
    while (getGhostVertexCount() > 0) {
        ghostVertex2InternalVertex(m_firstGhostVertexId);
    }
}

/*!
    Sends the given list of cells to the process with the specified rank.

    \param recvRanks are the receiver ranks
    \param trackPartitioning if set to true the function will return the changes
    done to the patch during the partitioning
    \param[in,out] ghostCellOwnershipChanges are the ghosts that will change
    ownership after the partitioning, on output the list will be updated
    removing the ghosts that are no longer on the partition (i.e., ghosts
    deleted or promoted to internal cell) and the ghosts whose ownership
    have been updated
    \result If the partitioning is tracked, returns a vector of adaption::Info
    with all the changes done to the patch during the partitioning, otherwise
    an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter_sendCells(const std::unordered_set<int> &recvRanks, bool trackPartitioning,
                                                                      std::unordered_map<long, int> *ghostCellOwnershipChanges)
{
    // Initialize adaption data
    std::vector<adaption::Info> partitioningData;
    if (trackPartitioning) {
        partitioningData.reserve(recvRanks.size());
    }

    // Detect if the rank is sending all its cells
    std::size_t nOutgoingsOverall = m_partitioningOutgoings.size();
    bool sendingAllCells = (nOutgoingsOverall == (std::size_t) getInternalCellCount());

    //
    // Send data to the receivers
    //
    const double ACTIVATE_MEMORY_LIMIT_THRESHOLD = 0.15;

    int patchDimension = getDimension();

    MPI_Request verticesSizeRequest = MPI_REQUEST_NULL;
    MPI_Request cellsSizeRequest    = MPI_REQUEST_NULL;

    MPI_Request verticesRequest = MPI_REQUEST_NULL;
    MPI_Request cellsRequest    = MPI_REQUEST_NULL;

    OBinaryStream verticesBuffer;
    OBinaryStream cellsBuffer;

    std::vector<long> neighsList;
    std::unordered_set<long> frontierNeighs;
    std::unordered_set<long> frontierVertices;
    std::unordered_set<long> frontierFaceVertices;
    std::unordered_set<CellHalfEdge, CellHalfEdge::Hasher> frontierEdges;

    std::unordered_set<long> outgoingCells;
    std::unordered_set<long> frameCells;
    std::unordered_set<long> haloCells;

    std::unordered_set<long> frameVertices;
    std::unordered_set<long> haloVertices;

    std::vector<long> cellSendList;
    std::unordered_set<long> vertexSendList;

    std::unordered_set<long> ghostCellsOverall;
    std::vector<long> frameCellsOverall;
    std::size_t frameCellIndexOverall = 0;

    for (int recvRank : recvRanks) {
        //
        // Fill the list of cells explicitly marked for sending
        //
        outgoingCells.clear();
        for (const auto &rankEntry : m_partitioningOutgoings) {
            int rank = rankEntry.second;
            if (rank != recvRank) {
                continue;
            }

            long cellId = rankEntry.first;
            outgoingCells.insert(cellId);
        }

        //
        // Identify frame and halo
        //
        // Along with the cells explicitly marked for sending, we need to send
        // also a halo of surrounding cells. Those cells will be used by the
        // receiver to connect the cells it receives to the existing cells and
        // to generate the ghost cells.
        //
        // The faces that tell apart cells explicitly marked for sending from
        // halo cells are called "frontier faces". We can identify frontier
        // faces looping through the outgoing cells an looking for faces with
        // a neighbour not explicitly marked for being sent to the receiver
        // rank.
        //
        // Frame is made up by cells explicitly marked for being sent to the
        // receiver rank that have a face, an edge or a vertex on a frontier
        // face.
        //
        // Halo is made up by cells not explicitly marked for being sent to
        // the receiver that have a face, an edge or a vertex on a frontier
        // face.
        //
        // We also identify cells that, at the end of the partitionig, will
        // be ghosts for the current rank. These are all the outgonig cells
        // that have a face/edge/vertex on the inner frontier. A face/edge/
        // vertex is on the inner frontier if its on the frontier and one of
        // the cells that contain that item is not an outgoing cell (i.e.,
        // is not explicitly marked for being sent to any rank).
        //
        // There are no halo nor frame cells if we are serializing a patch.
        //
        // NOTE: for three-dimensional unstructured non-conformal meshes we
        // need to explicitly search cells that have an edge on the frontier,
        // searching for cells that have a vertex on the frontier is not
        // enough (for unstructured meshes we may have a cell that touches the
        // frontier only through an edge).
        if (!m_partitioningSerialization) {
            haloCells.clear();
            frameCells.clear();

            frontierVertices.clear();
            frontierEdges.clear();
            for (long cellId : outgoingCells) {
                Cell &cell = m_cells.at(cellId);

                const int nFaces = cell.getFaceCount();
                for (int face = 0; face < nFaces; ++face) {
                    int nFaceAdjacencies = cell.getAdjacencyCount(face);
                    if (nFaceAdjacencies == 0) {
                        continue;
                    }

                    const long *faceAdjacencies = cell.getAdjacencies(face);
                    for (int i = 0; i < nFaceAdjacencies; ++i) {
                        // Check if the we are on the frontier
                        long neighId = faceAdjacencies[i];
                        if (outgoingCells.count(neighId) > 0) {
                            continue;
                        }

                        // Select the cell with only one adjacency
                        long frontierCellId;
                        long frontierFace;
                        Cell *frontierCell;
                        if (nFaceAdjacencies == 1) {
                            frontierCellId = cellId;
                            frontierFace = face;
                            frontierCell = &cell;
                        } else {
                            frontierCellId = neighId;
                            frontierFace = findAdjoinNeighFace(cellId, face, neighId);
                            frontierCell = &(m_cells.at(neighId));

                            assert(frontierFace >= 0);
                        }

                        // Clear neighbour list
                        frontierNeighs.clear();

                        //
                        // Add frontier face neighbours
                        //

                        // Check if the face is on the inner frontier
                        //
                        // To be on the inner frontier, the face needs to be
                        // on a cell not explicitly marked for being sent to
                        // any rank. We are iterating on the cells explicitly
                        // marked for being sent to the receiver rank, hence
                        // the current cell is definitely an outgoing cell.
                        // Therefore, the face will be on the inner frontier
                        // if the neighbour is not an outgoing cell.
                        bool innerFrontierFace = (m_partitioningOutgoings.count(neighId) == 0);

                        // Add the neighbours to the list
                        //
                        // If the face is on the inner frontier, the outgoing
                        // cell is the current cell (neighbour cannot be an
                        // outgoing cell, otherwise the face will not be on
                        // the inner frontier).
                        frontierNeighs.insert(cellId);
                        frontierNeighs.insert(neighId);

                        if (innerFrontierFace) {
                            ghostCellsOverall.insert(cellId);
                        }

                        //
                        // Add frontier vertex neighbours
                        //
                        ConstProxyVector<int> frontierFaceLocalVerticesIds = frontierCell->getFaceLocalVertexIds(frontierFace);
                        int nFrontierFaceVertices = frontierFaceLocalVerticesIds.size();

                        frontierFaceVertices.clear();
                        for (int k = 0; k < nFrontierFaceVertices; ++k) {
                            long vertexId = frontierCell->getFaceVertexId(frontierFace, k);
                            frontierFaceVertices.insert(vertexId);

                            // Avoid adding duplicate vertices
                            auto frontierVertexItr = frontierVertices.find(vertexId);
                            if (frontierVertexItr != frontierVertices.end()) {
                                continue;
                            }

                            // Add vertex to the list of frontier vertices
                            frontierVertices.insert(vertexId);

                            // Find vertex neighbours
                            neighsList.clear();
                            int vertex = frontierFaceLocalVerticesIds[k];
                            findCellVertexNeighs(frontierCellId, vertex, &neighsList);

                            // Check if the vertex is on the inner frontier
                            bool innerFrontierVertex = (m_partitioningOutgoings.count(frontierCellId) == 0);
                            if (!innerFrontierVertex) {
                                for (long frontierNeighId : neighsList) {
                                    if (m_partitioningOutgoings.count(frontierNeighId) == 0) {
                                        innerFrontierVertex = true;
                                    }
                                }
                            }

                            // Add neighbours to the list
                            for (long frontierNeighId : neighsList) {
                                frontierNeighs.insert(frontierNeighId);
                                if (innerFrontierVertex) {
                                    if (m_partitioningOutgoings.count(frontierNeighId) > 0) {
                                        ghostCellsOverall.insert(frontierNeighId);
                                    }
                                }
                            }
                        }

                        //
                        // Add frontier edge neighbours
                        //
                        if (patchDimension == 3) {
                            int nFrontierCellEdges = frontierCell->getEdgeCount();
                            for (int edge = 0; edge < nFrontierCellEdges; ++edge) {
                                // Edge information
                                ElementType edgeType = frontierCell->getEdgeType(edge);
                                int nEdgeVertices = frontierCell->getEdgeVertexCount(edge);

                                CellHalfEdge frontierEdge = CellHalfEdge(*frontierCell, edge, CellHalfEdge::WINDING_NATURAL);

                                // Discard edges that do not belong to the
                                // current frontier face
                                ConstProxyVector<long> edgeVertexIds = Cell::getVertexIds(edgeType, frontierEdge.getConnect().data());

                                bool isFrontierFaceEdge = true;
                                for (int k = 0; k < nEdgeVertices; ++k) {
                                    long edgeVertexId = edgeVertexIds[k];
                                    if (frontierFaceVertices.count(edgeVertexId) == 0) {
                                        isFrontierFaceEdge = false;
                                        break;
                                    }
                                }

                                if (!isFrontierFaceEdge) {
                                    continue;
                                }

                                // Avoid adding duplicate edges
                                frontierEdge.setWinding(CellHalfEdge::WINDING_REVERSE);
                                auto frontierEdgeItr = frontierEdges.find(frontierEdge);
                                if (frontierEdgeItr != frontierEdges.end()) {
                                    continue;
                                } else {
                                    frontierEdge.setWinding(CellHalfEdge::WINDING_NATURAL);
                                    frontierEdgeItr = frontierEdges.find(frontierEdge);
                                    if (frontierEdgeItr != frontierEdges.end()) {
                                        continue;
                                    }
                                }

                                // Add edge to the list of frontier edges
                                frontierEdges.insert(std::move(frontierEdge));

                                // Find edge neighbours
                                neighsList.clear();
                                findCellEdgeNeighs(frontierCellId, edge, &neighsList);

                                // Check if the vertex is on the inner frontier
                                bool innerFrontierEdge = (m_partitioningOutgoings.count(frontierCellId) == 0);
                                if (!innerFrontierEdge) {
                                    for (long frontierNeighId : neighsList) {
                                        if (m_partitioningOutgoings.count(frontierNeighId) == 0) {
                                            innerFrontierEdge = true;
                                        }
                                    }
                                }

                                // Add neighbours to the list
                                for (long frontierNeighId : neighsList) {
                                    frontierNeighs.insert(frontierNeighId);
                                    if (innerFrontierEdge) {
                                        if (m_partitioningOutgoings.count(frontierNeighId) > 0) {
                                            ghostCellsOverall.insert(frontierNeighId);
                                        }
                                    }
                                }
                            }
                        }

                        // Tell apart frame and halo cells
                        //
                        // Frame cells on inner frontier will be ghost cells
                        // for the current rank.
                        for (long frontierNeighId : frontierNeighs) {
                            if (outgoingCells.count(frontierNeighId) > 0) {
                                frameCells.insert(frontierNeighId);
                            } else {
                                haloCells.insert(frontierNeighId);
                            }
                        }
                    }
                }
            }
        }

        // Order cells to send
        //
        // To allow the receiver to efficently store the cells, first we
        // send the cells that are internal for the receiver, then the
        // halo cells.
        std::size_t nOutgoingCells = outgoingCells.size();
        std::size_t nHaloCells     = haloCells.size();

        cellSendList.resize(nOutgoingCells + nHaloCells);
        if (!sendingAllCells) {
            frameCellsOverall.resize(frameCellsOverall.size() + frameCells.size());
        }

        std::size_t outgoingIndex = 0;
        for (long cellId : outgoingCells) {
            cellSendList[outgoingIndex] = cellId;
            ++outgoingIndex;

            if (!sendingAllCells && frameCells.count(cellId) > 0) {
                frameCellsOverall[frameCellIndexOverall] = cellId;
                ++frameCellIndexOverall;
            }
        }

        std::size_t haloIndex = outgoingIndex;
        for (long cellId : haloCells) {
            cellSendList[haloIndex] = cellId;
            ++haloIndex;
        }

        //
        // Communicate cell buffer size
        //

        // Wait for previous comunications to finish
        if (cellsSizeRequest != MPI_REQUEST_NULL) {
            MPI_Wait(&cellsSizeRequest, MPI_STATUS_IGNORE);
        }

        // Start the communication
        long cellsBufferSize = 2 * sizeof(long) + 2 * nHaloCells * sizeof(int) + 2 * (nOutgoingCells + nHaloCells) * sizeof(bool);
        for (const long cellId : cellSendList) {
            cellsBufferSize += m_cells.at(cellId).getBinarySize();
        }

        MPI_Isend(&cellsBufferSize, 1, MPI_LONG, recvRank, m_partitioningCellsTag, m_communicator, &cellsSizeRequest);

        //
        // Create the list of vertices to send
        //
        // On serialization there are no frame nor halo cells, however the
        // vertices on border faces are considered frame vertices, becuase
        // they will be used to link the cells on the recevier (they will
        // be duplicate vertices on the receiver).
        //
        frameVertices.clear();
        haloVertices.clear();
        vertexSendList.clear();
        for (const long cellId : cellSendList) {
            const Cell &cell = m_cells.at(cellId);

            ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
            int nCellVertices = cellVertexIds.size();
            for (int j = 0; j < nCellVertices; ++j) {
                long vertexId = cellVertexIds[j];

                if (frameCells.count(cellId) > 0) {
                    frameVertices.insert(vertexId);
                } else if (haloCells.count(cellId) > 0) {
                    haloVertices.insert(vertexId);
                }

                vertexSendList.insert(vertexId);
            }
        }

        if (m_partitioningSerialization) {
            for (const long cellId : cellSendList) {
                const Cell &cell = m_cells.at(cellId);

                int nCellFaces = cell.getFaceCount();
                for (int face = 0; face < nCellFaces; ++face) {
                    if (!cell.isFaceBorder(face)) {
                        continue;
                    }

                    int nFaceVertices = cell.getFaceVertexCount(face);
                    for (int k = 0; k < nFaceVertices; ++k) {
                        long vertexId = cell.getFaceVertexId(face, k);
                        frameVertices.insert(vertexId);
                    }
                }
            }
        }

        //
        // Communicate vertex buffer size
        //

        // Wait for previous comunications to finish
        if (verticesSizeRequest != MPI_REQUEST_NULL) {
            MPI_Wait(&verticesSizeRequest, MPI_STATUS_IGNORE);
        }

        // Start the communication
        long verticesBufferSize = sizeof(long) + 2 * vertexSendList.size() * sizeof(bool);
        for (long vertexId : vertexSendList) {
            verticesBufferSize += m_vertices[vertexId].getBinarySize();
        }

        MPI_Isend(&verticesBufferSize, 1, MPI_LONG, recvRank, m_partitioningVerticesTag, m_communicator, &verticesSizeRequest);

        //
        // Send vertex data
        //

        // Wait for previous comunications to finish
        if (verticesRequest != MPI_REQUEST_NULL) {
            MPI_Wait(&verticesRequest, MPI_STATUS_IGNORE);
        }

        // Fill buffer with vertex data
        verticesBuffer.setSize(verticesBufferSize);
        verticesBuffer.seekg(0);

        verticesBuffer << (long) vertexSendList.size();
        for (long vertexId : vertexSendList) {
            // Cell information
            bool isFrame = (frameVertices.count(vertexId) > 0);
            bool isHalo  = (haloVertices.count(vertexId) > 0);

            verticesBuffer << isFrame;
            verticesBuffer << isHalo;

            // Certex data
            verticesBuffer << m_vertices[vertexId];
        }

        if (verticesBufferSize != (long) verticesBuffer.getSize()) {
            throw std::runtime_error ("Cell buffer size does not match calculated size");
        }

        MPI_Isend(verticesBuffer.data(), verticesBuffer.getSize(), MPI_CHAR, recvRank, m_partitioningVerticesTag, m_communicator, &verticesRequest);

        //
        // Send cell data
        //

        // Wait for previous comunications to finish
        if (cellsRequest != MPI_REQUEST_NULL) {
            MPI_Wait(&cellsRequest, MPI_STATUS_IGNORE);
        }

        // Fill the buffer with cell data
        cellsBuffer.setSize(cellsBufferSize);
        cellsBuffer.seekg(0);

        cellsBuffer << (long) nOutgoingCells;
        cellsBuffer << (long) nHaloCells;
        for (long cellId : cellSendList) {
            // Cell information
            bool isFrame = (frameCells.count(cellId) > 0);
            bool isHalo  = (haloCells.count(cellId) > 0);

            cellsBuffer << isFrame;
            cellsBuffer << isHalo;

            // Cell owner on receiver
            //
            // This is only needed if the cell is on the halo, in the other
            // case the owner is always the receiver itself.
            if (isHalo) {
                int cellOwnerOnReceiver;
                if (m_partitioningOutgoings.count(cellId) > 0) {
                    cellOwnerOnReceiver = m_partitioningOutgoings.at(cellId);
                } else if (m_ghostCellOwners.count(cellId) > 0) {
                    cellOwnerOnReceiver = m_ghostCellOwners.at(cellId);
                } else {
                    cellOwnerOnReceiver = m_rank;
                }

                cellsBuffer << cellOwnerOnReceiver;

                int ghostCellOwnershipChange;
                if (ghostCellOwnershipChanges->count(cellId) > 0) {
                    ghostCellOwnershipChange = ghostCellOwnershipChanges->at(cellId);
                } else {
                    ghostCellOwnershipChange = -1;
                }

                cellsBuffer << ghostCellOwnershipChange;
            }

            // Cell data
            const Cell &cell = m_cells.at(cellId);

            cellsBuffer << cell;
        }

        if (cellsBufferSize != (long) cellsBuffer.getSize()) {
            throw std::runtime_error ("Cell buffer size does not match calculated size");
        }

        MPI_Isend(cellsBuffer.data(), cellsBuffer.getSize(), MPI_CHAR, recvRank, m_partitioningCellsTag, m_communicator, &cellsRequest);

        //
        // Update adaption info
        //
        if (trackPartitioning) {
            // Update partition
            //
            // The ids of the cells send will be stored accordingly to the send
            // order, this is the same order that will be used on the process
            // that has received the cell. Since the order is the same, the
            // two processes are able to exchange cell data without additional
            // communications (they already know the list of cells for which
            // data is needed and the order in which these data will be sent).
            partitioningData.emplace_back();
            adaption::Info &partitioningInfo = partitioningData.back();
            partitioningInfo.entity = adaption::ENTITY_CELL;
            partitioningInfo.type   = adaption::TYPE_PARTITION_SEND;
            partitioningInfo.rank   = recvRank;

            partitioningInfo.previous.resize(nOutgoingCells);
            for (std::size_t i = 0; i < nOutgoingCells; ++i) {
                long cellId = cellSendList[i];
                partitioningInfo.previous[i] = cellId;
            }
        }

        // Delete outgoing cells not in the frame
        std::vector<long> deleteList;
        deleteList.clear();
        for (std::size_t i = 0; i < nOutgoingCells; ++i) {
            long cellId = cellSendList[i];
            if (frameCells.count(cellId) == 0) {
                deleteList.push_back(cellId);
            }
        }

        deleteCells(deleteList);

        // Prune cell adjacencies and interfaces
        //
        // At this stage we cannot fully update adjacencies and interfaces, but
        // we need to remove stale adjacencies and interfaces.
        pruneStaleAdjacencies();

        pruneStaleInterfaces();

        // If we are sending many cells try to reduced the used memory
        bool keepMemoryLimited = (nOutgoingCells > ACTIVATE_MEMORY_LIMIT_THRESHOLD * getInternalCellCount());
        if (keepMemoryLimited) {
            // Squeeze cells
            squeezeCells();

            // Squeeze interfaces
            squeezeInterfaces();

            // Delete orphan vertices
            deleteOrphanVertices();
            squeezeVertices();
        }
    }

    //
    // Update ghost and frame cells
    //

    // If the process is sending all its cells we can just clear the patch.
    if (!sendingAllCells) {
        std::vector<long> deleteList;
        std::vector<long> neighIds;

        // Delete frame cells or move them into the ghosts
        deleteList.clear();
        for (long cellId : frameCellsOverall) {
            bool moveToGhostCells = (ghostCellsOverall.count(cellId) > 0);
            if (moveToGhostCells) {
                int cellOwner = m_partitioningOutgoings.at(cellId);
                internalCell2GhostCell(cellId, cellOwner);
            } else {
                deleteList.push_back(cellId);
                ghostCellOwnershipChanges->erase(cellId);
            }
        }

        deleteCells(deleteList);

        // Prune cell adjacencies
        //
        // At this stage we cannot fully update the adjacencies, but we need
        // to remove stale adjacencies.
        pruneStaleAdjacencies();

        // Update cell adjacencies
        updateAdjacencies();

        // Delete stale ghosts
        //
        // Loop over all the ghosts and keep only the cells that have at least
        // one internal neighbour that is still on this process.
        //
        // Stale ghosts have to be deleted after processing frame cells,
        // because some of the frame cells moved into ghosts may be stale
        // ghosts.
        deleteList.clear();
        for (const auto &entry : m_ghostCellOwners) {
            long cellId = entry.first;
            bool keep = false;

            // First do a cheap search among the face neighbours.
            const Cell &cell = m_cells.at(cellId);
            const long *adjacencies = cell.getAdjacencies();
            int nCellAdjacencies = cell.getAdjacencyCount();
            for (int i = 0; i < nCellAdjacencies; ++i) {
                long neighId = adjacencies[i];
                if (m_ghostCellOwners.count(neighId) > 0) {
                    continue;
                } else if (m_partitioningOutgoings.count(neighId) > 0) {
                    continue;
                }

                keep = true;
                break;
            }

            // If we haven't find a neighbour that is still on this
            // process, do a search among all the neighbours.
            if (!keep) {
                neighIds.clear();
                findCellNeighs(cellId, &neighIds);
                for (long neighId : neighIds) {
                    if (m_ghostCellOwners.count(neighId) > 0) {
                        continue;
                    } else if (m_partitioningOutgoings.count(neighId) > 0) {
                        continue;
                    }

                    keep = true;
                    break;
                }
            }

            // Add the cell to the delete list
            if (!keep) {
                deleteList.push_back(cellId);
                ghostCellOwnershipChanges->erase(cellId);
            }
        }

        deleteCells(deleteList);

        // Prune cell adjacencies and interfaces
        //
        // At this stage we cannot fully update adjacencies and interfaces, but
        // we need to remove stale adjacencies and interfaces.
        pruneStaleAdjacencies();

        pruneStaleInterfaces();

        // Delete orphan vertices
        deleteOrphanVertices();
    } else {
        // The process has sent all its cells, the patch is now empty
        reset();
        ghostCellOwnershipChanges->clear();
    }

    // Wait for previous communications to finish
    if (cellsSizeRequest != MPI_REQUEST_NULL) {
        MPI_Wait(&cellsSizeRequest, MPI_STATUS_IGNORE);
    }

    if (verticesSizeRequest != MPI_REQUEST_NULL) {
        MPI_Wait(&verticesSizeRequest, MPI_STATUS_IGNORE);
    }

    if (verticesRequest != MPI_REQUEST_NULL) {
        MPI_Wait(&verticesRequest, MPI_STATUS_IGNORE);
    }

    if (cellsRequest != MPI_REQUEST_NULL) {
        MPI_Wait(&cellsRequest, MPI_STATUS_IGNORE);
    }

    // Return adaption info
    return partitioningData;
}

/*!
    Receives a list of cells from the specified process.

    \param sendRanks are the rank of the processes sending the cells
    \param trackPartitioning if set to true the function will return the changes
    done to the patch during the partitioning
    \param[in,out] ghostCellOwnershipChanges are the ghosts that will change
    ownership after the partitioning, on output the list will be updated
    removing the ghosts that are no longer on the partition (i.e., ghosts
    deleted or promoted to internal cell) and the ghosts whose ownership
    have been updated
    \result If the partitioning is tracked, returns a vector of adaption::Info
    with all the changes done to the patch during the partitioning, otherwise
    an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter_receiveCells(const std::unordered_set<int> &sendRanks, bool trackPartitioning,
                                                                         std::unordered_map<long, int> *ghostCellOwnershipChanges)
{
    //
    // Start receiving buffer sizes
    //
    int nSendRanks = sendRanks.size();

    std::vector<MPI_Request> cellsSizeRequests(nSendRanks, MPI_REQUEST_NULL);
    std::vector<MPI_Request> verticesSizeRequests(nSendRanks, MPI_REQUEST_NULL);

    std::vector<long> cellsBufferSizes(nSendRanks);
    std::vector<long> verticesBufferSizes(nSendRanks);

    std::unordered_map<int, int> sendRankIndexes;
    sendRankIndexes.reserve(nSendRanks);

    for (int sendRank : sendRanks) {
        int sendRankIndex = sendRankIndexes.size();
        sendRankIndexes[sendRank] = sendRankIndex;

        // Cell buffer size
        long &cellsBufferSize = cellsBufferSizes[sendRankIndex];
        MPI_Request &cellsSizeRequest = cellsSizeRequests[sendRankIndex];
        MPI_Irecv(&cellsBufferSize, 1, MPI_LONG, sendRank, m_partitioningCellsTag, m_communicator, &cellsSizeRequest);

        // Vertex buffer size
        long &verticesBufferSize = verticesBufferSizes[sendRankIndex];
        MPI_Request &verticesSizeRequest = verticesSizeRequests[sendRankIndex];
        MPI_Irecv(&verticesBufferSize, 1, MPI_LONG, sendRank, m_partitioningVerticesTag, m_communicator, &verticesSizeRequest);
    }

    // Initialize data structured for data exchange
    std::vector<MPI_Request> cellsRequests(nSendRanks, MPI_REQUEST_NULL);
    std::vector<MPI_Request> verticesRequests(nSendRanks, MPI_REQUEST_NULL);

    std::vector<std::unique_ptr<IBinaryStream>> cellsBuffers(nSendRanks);
    std::vector<std::unique_ptr<IBinaryStream>> verticesBuffers(nSendRanks);

    // Fill partitioning info
    std::vector<adaption::Info> partitioningData;
    if (trackPartitioning) {
        partitioningData.resize(nSendRanks);
        for (int sendRank : sendRanks) {
            int sendRankIndex = sendRankIndexes.at(sendRank);

            adaption::Info &partitioningInfo = partitioningData[sendRankIndex];
            partitioningInfo.entity = adaption::ENTITY_CELL;
            partitioningInfo.type   = adaption::TYPE_PARTITION_RECV;
            partitioningInfo.rank   = sendRank;
        }
    }

    // Mark border interfaces of ghost cells as dangling
    //
    // We may received cells that connect to the existing mesh through one
    // of the faces that are now borders. Marking those border interfaces as
    // dangling allows to delete them and create new internal interfaces.
    for (auto itr = ghostCellConstBegin(); itr != ghostCellConstEnd(); ++itr) {
        const Cell &cell = *itr;
        const long *interfaces = cell.getInterfaces();
        const int nCellInterfaces = cell.getInterfaceCount();

        setCellAlterationFlags(cell.getId(), FLAG_INTERFACES_DIRTY);

        for (int k = 0; k < nCellInterfaces; ++k) {
            long interfaceId = interfaces[k];
            const Interface &interface = getInterface(interfaceId);
            if (interface.isBorder()) {
                setInterfaceAlterationFlags(interfaceId, FLAG_DANGLING);
            }
        }
    }

    // Receive data
    int patchRank = getRank();

    int nSendCompleted = 0;
    std::vector<int> sendCompletedIndexes(nSendRanks);
    std::vector<MPI_Status> sendCompletedStatuses(nSendRanks);

    std::vector<long> neighIds;

    std::unordered_set<long> duplicateCellsCandidates;
    KdTree<3, Vertex, long> duplicateVerticesCandidatesTree;

    std::unordered_map<long, long> verticesMap;

    std::unordered_set<long> validReceivedAdjacencies;
    std::vector<long> addedCells;
    std::unordered_map<long, FlatVector2D<long>> duplicateCellsReceivedAdjacencies;
    std::unordered_map<long, long> cellsMap;

    std::unordered_set<int> awaitingSendRanks = sendRanks;
    while (!awaitingSendRanks.empty()) {
        //
        // Duplicate cell and vertex candidates
        //

        // Duplicate cell candidates
        //
        // The sender is sending a halo of cells surroudning the cells explicitly
        // marked for sending. This halo will be used for connecting the received
        // cells to the existing ones and to build the ghosts.
        //
        // We may receive duplicate cells for ghost targets and ghost sources. The
        // data structures that contains the list of these cells may not be updated
        // so we need to build the list on the fly. THe list will contain ghost
        // cells (target) and their neighbours (sources).
        duplicateCellsCandidates.clear();
        for (auto &entry : m_ghostCellOwners) {
            long ghostCellId = entry.first;

            duplicateCellsCandidates.insert(ghostCellId);

            findCellNeighs(ghostCellId, &neighIds);
            for (long neighId : neighIds) {
                duplicateCellsCandidates.insert(neighId);
            }
        }

        // Duplicate vertex candidates
        //
        // During normal partitioning, vertices of the duplicate cells candidates
        // are candidates for duplicate vertices.
        //
        // During mesh serialization, we will not receive duplicate cells, but in
        // order to link the recevied cells with the current one we still need to
        // properly identify duplicate vertices. Since we have delete all ghost
        // information, all the vertices of border faces need to be added to the
        // duplicate vertex candidates.
        //
        // Duplicate vertex candidates are stored in a kd-tree. The kd-tree stores
        // a pointer to the vertices. If we try to store in the kd-tree a pointers
        // to the vertices of the patch, the first resize of the vertex container
        // would invalidate the pointer. Create a copy of the vertices and store
        // the pointer to that copy.
        std::unordered_map<long, Vertex> duplicateVerticesCandidates;
        if (!m_partitioningSerialization) {
            for (long cellId : duplicateCellsCandidates) {
                const Cell &cell = m_cells.at(cellId);
                ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
                int nCellVertices = cellVertexIds.size();
                for (int k = 0; k < nCellVertices; ++k) {
                    long vertexId = cellVertexIds[k];
                    duplicateVerticesCandidates.insert({vertexId, m_vertices.at(vertexId)});
                }
            }
        } else {
            for (const Cell &cell : m_cells) {
                int nCellFaces = cell.getFaceCount();
                for (int face = 0; face < nCellFaces; ++face) {
                    if (!cell.isFaceBorder(face)) {
                        continue;
                    }

                    int nFaceVertices = cell.getFaceVertexCount(face);
                    for (int k = 0; k < nFaceVertices; ++k) {
                        long vertexId = cell.getFaceVertexId(face, k);
                        duplicateVerticesCandidates.insert({vertexId, m_vertices.at(vertexId)});
                    }
                }
            }
        }

        duplicateVerticesCandidatesTree.clear();
        for (auto &entry : duplicateVerticesCandidates) {
            long vertexId = entry.first;
            Vertex &vertex = entry.second;

            duplicateVerticesCandidatesTree.insert(&vertex, vertexId);
        }

        //
        // Identify rank that is sending data
        //
        // We need to start the receives and identify a rank for which both
        // cell and vertex recevies are active.
        int sendRank = -1;
        while (sendRank < 0) {
            // Start receiving vertices
            MPI_Waitsome(nSendRanks, verticesSizeRequests.data(), &nSendCompleted, sendCompletedIndexes.data(), sendCompletedStatuses.data());
            for (int i = 0; i < nSendCompleted; ++i) {
                int sourceIndex = sendCompletedIndexes[i];
                int sourceRank = sendCompletedStatuses[i].MPI_SOURCE;

                std::size_t bufferSize = verticesBufferSizes[sourceIndex];
                std::unique_ptr<IBinaryStream> buffer = std::unique_ptr<IBinaryStream>(new IBinaryStream(bufferSize));
                MPI_Request *request = verticesRequests.data() + sourceIndex;
                MPI_Irecv(buffer->data(), buffer->getSize(), MPI_CHAR, sourceRank, m_partitioningVerticesTag, m_communicator, request);
                verticesBuffers[sourceIndex] = std::move(buffer);
            }

            // Start receiving cells
            MPI_Waitsome(nSendRanks, cellsSizeRequests.data(), &nSendCompleted, sendCompletedIndexes.data(), sendCompletedStatuses.data());
            for (int i = 0; i < nSendCompleted; ++i) {
                int sourceIndex = sendCompletedIndexes[i];
                int sourceRank  = sendCompletedStatuses[i].MPI_SOURCE;

                std::size_t bufferSize = cellsBufferSizes[sourceIndex];
                std::unique_ptr<IBinaryStream> buffer = std::unique_ptr<IBinaryStream>(new IBinaryStream(bufferSize));
                MPI_Request *request = cellsRequests.data() + sourceIndex;
                MPI_Irecv(buffer->data(), buffer->getSize(), MPI_CHAR, sourceRank, m_partitioningCellsTag, m_communicator, request);
                cellsBuffers[sourceIndex] = std::move(buffer);
            }

            // Search for a rank that started sending both vertices and cells.
            // If there aren't any, wait for other request to complete.
            for (int awaitingSendRank : awaitingSendRanks) {
                int sourceIndex = sendRankIndexes.at(awaitingSendRank);
                if (verticesRequests[sourceIndex] == MPI_REQUEST_NULL) {
                    continue;
                } else if (cellsRequests[sourceIndex] == MPI_REQUEST_NULL) {
                    continue;
                }

                sendRank = awaitingSendRank;
                break;
            }
        }

        int sendRankIndex = sendRankIndexes.at(sendRank);

        //
        // Process vertices
        //

        // Wait until vertex communication is completed
        MPI_Request *verticesRequest = verticesRequests.data() + sendRankIndex;
        MPI_Wait(verticesRequest, MPI_STATUS_IGNORE);

        // Add vertices
        IBinaryStream &verticesBuffer = *(verticesBuffers[sendRankIndex]);

        long nRecvVertices;
        verticesBuffer >> nRecvVertices;

        reserveVertices(getVertexCount() + nRecvVertices);

        verticesMap.clear();
        for (long i = 0; i < nRecvVertices; ++i) {
            // Cell data
            bool isFrame;
            verticesBuffer >> isFrame;

            bool isHalo;
            verticesBuffer >> isHalo;

            Vertex vertex;
            verticesBuffer >> vertex;
            long originalVertexId = vertex.getId();

            // Set vertex interior flag
            //
            // All new vertices will be temporarly added as internal vertices,
            // is needed they will be converted to ghost vertices when updating
            // ghost information.
            vertex.setInterior(true);

            // Check if the vertex is a duplicate
            //
            // Only frame and halo vertices may be a duplicate.
            //
            // If the vertex is a duplicate the function that check its existance
            // will return the current id of the vertex.
            long vertexId;

            bool isDuplicate = (isFrame || isHalo);
            if (isDuplicate) {
                isDuplicate = (duplicateVerticesCandidatesTree.exist(&vertex, vertexId) >= 0);
            }

            // Add the vertex
            //
            // If the id of the received vertex is already assigned, let the
            // patch generate a new id. Otherwise, keep the id of the received
            // vertex.
            if (!isDuplicate) {
                if (m_vertexIdGenerator.isAssigned(originalVertexId)) {
                    vertex.setId(Vertex::NULL_ID);
                }

                VertexIterator vertexIterator = addVertex(std::move(vertex));
                vertexId = vertexIterator.getId();
            }

            if (originalVertexId != vertexId) {
                verticesMap.insert({{originalVertexId, vertexId}});
            }
        }

        // Cleanup
        verticesBuffers[sendRankIndex].reset();
        std::unordered_map<long, Vertex>().swap(duplicateVerticesCandidates);

        //
        // Process cells
        //

        // Wait until cell communication is completed
        MPI_Request *cellsRequest = cellsRequests.data() + sendRankIndex;
        MPI_Wait(cellsRequest, MPI_STATUS_IGNORE);

        // Receive cell data
        //
        // Internal cells will be sent first.
        IBinaryStream &cellsBuffer = *(cellsBuffers[sendRankIndex]);

        long nReceivedInternalCells;
        cellsBuffer >> nReceivedInternalCells;

        long nReceivedHalo;
        cellsBuffer >> nReceivedHalo;

        long nReceivedCells = nReceivedInternalCells + nReceivedHalo;

        reserveCells(getCellCount() + nReceivedCells);

        if (trackPartitioning) {
            // Only internal cells are tracked.
            partitioningData[sendRankIndex].current.reserve(nReceivedInternalCells);
        }

        validReceivedAdjacencies.clear();
        validReceivedAdjacencies.reserve(nReceivedCells);

        addedCells.clear();
        addedCells.reserve(nReceivedCells);

        duplicateCellsReceivedAdjacencies.clear();

        cellsMap.clear();

        for (long i = 0; i < nReceivedCells; ++i) {
            // Cell data
            bool isFrame;
            cellsBuffer >> isFrame;

            bool isHalo;
            cellsBuffer >> isHalo;

            int cellOwner;
            int ghostCellOwnershipChange;
            if (isHalo) {
                cellsBuffer >> cellOwner;
                cellsBuffer >> ghostCellOwnershipChange;
            } else {
                cellOwner = patchRank;
                ghostCellOwnershipChange = -1;
            }

            Cell cell;
            cellsBuffer >> cell;

            long cellOriginalId = cell.getId();

            // Set cell interior flag
            bool isInterior = (cellOwner == patchRank);
            cell.setInterior(isInterior);

            // Remap connectivity
            if (!verticesMap.empty()) {
                cell.renumberVertices(verticesMap);
            }

            // Check if the cells is a duplicate
            //
            // The received cell may be one of the current ghosts.
            long cellId = Cell::NULL_ID;
            if (isHalo || isFrame) {
                for (long candidateId : duplicateCellsCandidates) {
                    const Cell &candidateCell = m_cells.at(candidateId);
                    if (cell.hasSameConnect(candidateCell)) {
                        cellId = candidateId;
                        break;
                    }
                }
            }

            // If the cell is not a duplicate add it in the cell data structure,
            // otherwise merge the connectivity of the duplicate cell to the
            // existing cell. This ensure that the received cell will be
            // properly connected to the received cells
            bool isTracked = false;
            if (cellId < 0) {
                // Reset the interfaces of the cell
                if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
                    cell.resetInterfaces();
                }

                // Add cell
                //
                // If the id of the received cell is already assigned, let the
                // patch generate a new id. Otherwise, keep the id of the received
                // cell.
                if (m_cellIdGenerator.isAssigned(cellOriginalId)){
                    cell.setId(Cell::NULL_ID);
                }

                CellIterator cellIterator = addCell(std::move(cell), cellOwner);
                cellId = cellIterator.getId();
                addedCells.push_back(cellId);

                // Setup ghost ownership changes
                if (!cellIterator->isInterior()) {
                    if (ghostCellOwnershipChange >= 0) {
                        ghostCellOwnershipChanges->insert({cellId, ghostCellOwnershipChange});
                    }
                }

                // Interior cells are tacked
                isTracked = (trackPartitioning && isInterior);
            } else {
                // Check if the existing cells needs to become an internal cell
                const Cell &localCell = m_cells[cellId];
                if (isInterior && !localCell.isInterior()) {
                    ghostCell2InternalCell(cellId);
                    ghostCellOwnershipChanges->erase(cellId);
                    isTracked = trackPartitioning;
                }

                // Save the adjacencies of the received cell, this adjacencies
                // will link together the recevied cell to the existing ones.
                FlatVector2D<long> &cellAdjacencies = duplicateCellsReceivedAdjacencies[cellId];

                int nCellFaces = cell.getFaceCount();
                int nCellAdjacencies = cell.getAdjacencyCount();
                cellAdjacencies.reserve(nCellFaces, nCellAdjacencies);
                for (int face = 0; face < nCellFaces; ++face) {
                    int nFaceAdjacencies = cell.getAdjacencyCount(face);
                    const long *faceAdjacencies = cell.getAdjacencies(face);
                    cellAdjacencies.pushBack(nFaceAdjacencies, faceAdjacencies);
                }
            }

            // The interfaces of the cell need to be updated
            setCellAlterationFlags(cellId, FLAG_INTERFACES_DIRTY);

            // Add the cell to the cell map
            if (cellOriginalId != cellId) {
                cellsMap.insert({{cellOriginalId, cellId}});
            }

            // Add original cell id to the list of valid received adjacencies
            validReceivedAdjacencies.insert(cellOriginalId);

            // Update tracking information
            //
            // Only new internal cells and internal cells that were
            // previously ghosts are tracked.
            //
            // The ids of the cells send will be stored accordingly to the
            // receive order, this is the same order that will be used on the
            // process that has sent the cell. Since the order is the same,
            // the two processes are able to exchange cell data without
            // additional communications (they already know the list of cells
            // for which data is needed and the order in which these data will
            // be sent).
            if (isTracked) {
                partitioningData[sendRankIndex].current.emplace_back(cellId);
            }
        }

        // Cleanup
        cellsBuffers[sendRankIndex].reset();

        // Update adjacencies amonge added cells
        for (long cellId : addedCells) {
            Cell &cell = m_cells.at(cellId);

            // Remove stale adjacencies
            int nCellFaces = cell.getFaceCount();
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = cell.getAdjacencyCount(face);
                const long *faceAdjacencies = cell.getAdjacencies(face);

                int k = 0;
                while (k < nFaceAdjacencies) {
                    long senderAdjacencyId = faceAdjacencies[k];
                    if (validReceivedAdjacencies.count(senderAdjacencyId) == 0) {
                        cell.deleteAdjacency(face, k);
                        --nFaceAdjacencies;
                    } else {
                        ++k;
                    }
                }
            }

            // Remap adjacencies
            if (!cellsMap.empty()) {
                int nCellAdjacencies = cell.getAdjacencyCount();
                long *cellAdjacencies = cell.getAdjacencies();
                for (int k = 0; k < nCellAdjacencies; ++k) {
                    long &cellAdjacencyId = cellAdjacencies[k];
                    auto cellsMapItr = cellsMap.find(cellAdjacencyId);
                    if (cellsMapItr != cellsMap.end()) {
                        cellAdjacencyId = cellsMapItr->second;
                    }
                }
            }

            // Cell adjacencies are up-to-date
            if (duplicateCellsReceivedAdjacencies.count(cellId)) {
                unsetCellAlterationFlags(cellId, FLAG_ADJACENCIES_DIRTY);
            }
        }

        // Link received cells with the initial cells
        //
        // If we are serializing the patch, we don't have enough information to
        // link the recevied cells with the initial cells (ghost cells have been
        // deleted before receiving the cells). Cells that need to be linked will
        // have some faces without any adjacency, therefore to link those cells
        // we rebuild the adjacencies of all the cells that havefaces with no
        // adjacencis. Authentic border cells will have their adjacencies rebuilt
        // also if this may not be needed, but this is still the faster way to
        // link the received cells.
        if (!m_partitioningSerialization) {
            for (auto &entry : duplicateCellsReceivedAdjacencies) {
                long cellId = entry.first;
                Cell &cell = m_cells.at(cellId);

                int nCellFaces = cell.getFaceCount();
                FlatVector2D<long> &cellReceivedAdjacencies = entry.second;
                for (int face = 0; face < nCellFaces; ++face) {
                    int nFaceLinkAdjacencies = cellReceivedAdjacencies.getItemCount(face);
                    for (int k = 0; k < nFaceLinkAdjacencies; ++k) {
                        // We need to updated the adjacencies only if they are cells
                        // that have been send.
                        long receivedAdjacencyId = cellReceivedAdjacencies.getItem(face, k);
                        if (validReceivedAdjacencies.count(receivedAdjacencyId) == 0) {
                            continue;
                        }

                        // If the send cell is already in the adjacency list there is
                        // nothing to update.
                        long localAdjacencyId;
                        auto ajacenciyCellMapItr = cellsMap.find(receivedAdjacencyId);
                        if (ajacenciyCellMapItr != cellsMap.end()) {
                            localAdjacencyId = ajacenciyCellMapItr->second;
                        } else {
                            localAdjacencyId = receivedAdjacencyId;
                        }

                        if (cell.findAdjacency(face, localAdjacencyId) >= 0) {
                            continue;
                        }

                        cell.pushAdjacency(face, localAdjacencyId);
                    }
                }

                unsetCellAlterationFlags(cellId, FLAG_ADJACENCIES_DIRTY);
            }
        } else {
            for (const Cell &cell : m_cells) {
                int nCellFaces = cell.getFaceCount();
                for (int face = 0; face < nCellFaces; ++face) {
                    if (cell.isFaceBorder(face)) {
                        long cellId = cell.getId();
                        setCellAlterationFlags(cellId, FLAG_ADJACENCIES_DIRTY);
                        break;
                    }
                }
            }
        }

        // Receive is now complete
        awaitingSendRanks.erase(sendRank);
    }

    // Update adjacencies
    if (getAdjacenciesBuildStrategy() == ADJACENCIES_AUTOMATIC) {
        updateAdjacencies();
    }

    // Update interfaces
    if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
        updateInterfaces();
    }

    // Return adaption data
    return partitioningData;
}

/*!
	Cleanup patch data structured after the partitioning.

	Default implementation is a no-op function.
*/
void PatchKernel::_partitioningCleanup()
{
}

/*!
	Gets the rank of the process that owns the specified cell.

	\param id is the id of the requested cell
	\result The rank that owns the specified cell.
*/
int PatchKernel::getCellRank(long id) const
{
	auto cellOwner = m_ghostCellOwners.find(id);
	if (cellOwner == m_ghostCellOwners.end()) {
		return m_rank;
	} else {
		return cellOwner->second;
	}
}

/*!
	Gets the halo layer of the specified cell.

	\param id is the id of the requested cell
	\result The halo layer of the specified cell.
*/
int PatchKernel::getCellHaloLayer(long id) const
{
	const Cell &cell = getCell(id);
	if (cell.isInterior()) {
		return 0;
	} else {
		return -1;
	}
}

/*!
	Gets the rank of the process that owns the specified vertex.

	\param id is the id of the requested vertex
	\result The rank that owns the specified vertex.
*/
int PatchKernel::getVertexRank(long id) const
{
	auto vertexOwner = m_ghostVertexOwners.find(id);
	if (vertexOwner == m_ghostVertexOwners.end()) {
		return m_rank;
	} else {
		return vertexOwner->second;
	}
}

/*!
	Check if the processes associated with the specified rank is a neighbour.

	\param rank is the rank associated with the process
	\result True is the process is a neighbour, false otherwise.
*/
bool PatchKernel::isRankNeighbour(int rank)
{
	return (m_ghostCellExchangeTargets.count(rank) > 0);
}

/*!
	Get a list of neighbour ranks.

	\result A list of neighbour ranks.
*/
std::vector<int> PatchKernel::getNeighbourRanks()
{
	std::vector<int> neighRanks;
	neighRanks.reserve(m_ghostCellExchangeTargets.size());
	for (const auto &entry : m_ghostCellExchangeTargets) {
		neighRanks.push_back(entry.first);
	}

	return neighRanks;
}

/*!
	Gets a constant reference to the vertices that define the "targets" for the
	exchange of data on ghost vertices. For each process, the corresponding list
	of "targets" is returned.

	During data exchange, each partition will send the data of its "source"
	vertices to the corresponding "target" vertices (i.e., ghost vertices) on
	other partitions.

	Source and target vertices are ordered using a geometrical criterion (the
	same criterion is used on all the partitions), in this way the n-th source
	on a partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the vertices that define the "targets" for
	the exchange of data on ghost vertices.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostVertexExchangeTargets() const
{
	return m_ghostVertexExchangeTargets;
}

/*!
	Gets a constant reference to the vertices that define the "targets" for the
	exchange of data on ghost vertices with the specified process.

	During data exchange, each partition will send the data of its "source"
	vertices to the corresponding "target" vertices (i.e., ghost vertices) on
	other partitions.

	Source and target vertices are ordered using a geometrical criterion (the
	same criterion is used on all the partitions), in this way the n-th source
	on a partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the vertices that define the "targets" for
	the exchange of data on ghost vertices with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostVertexExchangeTargets(int rank) const
{
	return m_ghostVertexExchangeTargets.at(rank);
}

/*!
	Gets a constant reference to the vertices that define the "sources" for the
	exchange of data on ghost vertices. For each process, the corresponding list
	of "sources" is returned.

	Sources are internal vertices (i.e. vertices owned by the current partition)
	that are ghost vertices on other partitions. When exchanging data on ghost
	vertices, these vertices will be the sources form which data will be read
	from.

	Source and target vertices are ordered using a geometrical criterion (the
	same criterion is used on all the partitions), in this way the n-th source
	on a partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be send to
	\result A constant reference to the vertices that define the "sources" for
	the exchange of data on ghost vertices with the specified process.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostVertexExchangeSources() const
{
	return m_ghostVertexExchangeSources;
}

/*!
	Gets a constant reference to the vertices that define the "sources" for the
	exchange of data on ghost vertices with the specified process.

	Sources are internal vertices (i.e. vertices owned by the current partition)
	that are ghost vertices on other partitions. When exchanging data on ghost
	vertices, these vertices will be the sources form which data will be read
	from.

	Source and target vertices are ordered using a geometrical criterion (the
	same criterion is used on all the partitions), in this way the n-th source
	on a partition will correspond to the n-th target on the other partition.

	\param rank is the rank to which data will be send
	\result A constant reference to the vertices that define the "sources" for
	the exchange of data on ghost vertices with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostVertexExchangeSources(int rank) const
{
	return m_ghostVertexExchangeSources.at(rank);
}

/*!
	Gets a constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells. For each process, the corresponding list
	of "targets" is returned.

	During data exchange, each partition will send the data of its "source"
	cells to the corresponding "target" cells (i.e., ghost cells) on other
	partitions.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostCellExchangeTargets() const
{
	return m_ghostCellExchangeTargets;
}

/*!
	Gets a constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells. For each process, the corresponding list
	of "targets" is returned.

	During data exchange, each partition will send the data of its "source"
	cells to the corresponding "target" cells (i.e., ghost cells) on other
	partitions.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeTargets() const
{
	return getGhostCellExchangeTargets();
}

/*!
	Gets a constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells with the specified process.

	During data exchange, each partition will send the data of its "source"
	cells to the corresponding "target" cells (i.e., ghost cells) on other
	partitions.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostCellExchangeTargets(int rank) const
{
	return m_ghostCellExchangeTargets.at(rank);
}

/*!
	Gets a constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells with the specified process.

	During data exchange, each partition will send the data of its "source"
	cells to the corresponding "target" cells (i.e., ghost cells) on other
	partitions.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank data will be received from
	\result A constant reference to the cells that define the "targets" for the
	exchange of data on ghost cells with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostExchangeTargets(int rank) const
{
	return getGhostCellExchangeTargets(rank);
}

/*!
	Gets a constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells. For each process, the corresponding list
	of "sources" is returned.

	Sources are internal cells (i.e. cells owned by the current partition) that
	are ghost cells on other partitions. When exchanging data on ghost cells,
	these cells will be the sources form which data will be read from.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\result A constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostCellExchangeSources() const
{
	return m_ghostCellExchangeSources;
}

/*!
	Gets a constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells. For each process, the corresponding list
	of "sources" is returned.

	Sources are internal cells (i.e. cells owned by the current partition) that
	are ghost cells on other partitions. When exchanging data on ghost cells,
	these cells will be the sources form which data will be read from.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\result A constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeSources() const
{
	return getGhostCellExchangeSources();
}

/*!
	Gets a constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells with the specified process.

	Sources are internal cells (i.e. cells owned by the current partition) that
	are ghost cells on other partitions. When exchanging data on ghost cells,
	these cells will be the sources form which data will be read from.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank to which data will be send
	\result A constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostCellExchangeSources(int rank) const
{
	return m_ghostCellExchangeSources.at(rank);
}

/*!
	Gets a constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells with the specified process.

	Sources are internal cells (i.e. cells owned by the current partition) that
	are ghost cells on other partitions. When exchanging data on ghost cells,
	these cells will be the sources form which data will be read from.

	Source and target cells are ordered using a geometrical criterion (the same
	criterion is used on all the partitions), in this way the n-th source on a
	partition will correspond to the n-th target on the other partition.

	\param rank is the rank to which data will be send
	\result A constant reference to the cells that define the "sources" for the
	exchange of data on ghost cells with the specified process.
*/
const std::vector<long> & PatchKernel::getGhostExchangeSources(int rank) const
{
	return getGhostCellExchangeSources(rank);
}

/*!
	Sets the owner of the specified ghost vertex.

	\param id is the id of the ghost vertex
	\param rank is the rank of the process that owns the ghost vertex
*/
void PatchKernel::setGhostVertexOwner(int id, int rank)
{
	auto ghostVertexOwnerItr = m_ghostVertexOwners.find(id);
	if (ghostVertexOwnerItr != m_ghostVertexOwners.end()) {
		ghostVertexOwnerItr->second = rank;
	} else {
		m_ghostVertexOwners.insert({id, rank});
	}

	setPartitioningInfoDirty(true);
}

/*!
	Unsets the owner of the specified ghost vertex.

	\param id is the id of the ghost vertex
*/
void PatchKernel::unsetGhostVertexOwner(int id)
{
	auto ghostVertexOwnerItr = m_ghostVertexOwners.find(id);
	if (ghostVertexOwnerItr == m_ghostVertexOwners.end()) {
		return;
	}

	m_ghostVertexOwners.erase(ghostVertexOwnerItr);

	setPartitioningInfoDirty(true);
}

/*!
	Clear the owners of all the ghost vertices.

	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::clearGhostVertexOwners()
{
	m_ghostVertexOwners.clear();

	setPartitioningInfoDirty(true);
}

/*!
	Sets the owner of the specified ghost.

	\param id is the id of the ghost cell
	\param rank is the rank of the process that owns the ghost cell
*/
void PatchKernel::setGhostCellOwner(int id, int rank)
{
	auto ghostCellOwnerItr = m_ghostCellOwners.find(id);
	if (ghostCellOwnerItr != m_ghostCellOwners.end()) {
		ghostCellOwnerItr->second = rank;
	} else {
		m_ghostCellOwners.insert({id, rank});
	}

	setPartitioningInfoDirty(true);
}

/*!
	Unsets the owner of the specified ghost.

	\param id is the id of the ghost cell
*/
void PatchKernel::unsetGhostCellOwner(int id)
{
	auto ghostCellOwnerItr = m_ghostCellOwners.find(id);
	if (ghostCellOwnerItr == m_ghostCellOwners.end()) {
		return;
	}

	m_ghostCellOwners.erase(ghostCellOwnerItr);

	setPartitioningInfoDirty(true);
}

/*!
	Clear the owners of all the ghosts.

	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::clearGhostCellOwners()
{
	m_ghostCellOwners.clear();

	setPartitioningInfoDirty(true);
}

/*!
	Checks if the partitioning information are dirty.

	\param global if set to true, the dirty status will be evaluated globally
	across all the partitions
	\result Returns true if the partitioning information are dirty, false
	otherwise.
*/
bool PatchKernel::arePartitioningInfoDirty(bool global) const
{
	if (!isPartitioned()) {
		return false;
	}

	bool partitioningInfoDirty = m_partitioningInfoDirty;
	if (global) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &partitioningInfoDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}

	return partitioningInfoDirty;
}

/*!
	Sets if the partitioning information are dirty.

	\param dirty controls if the partitioning information will be set as dirty
*/
void PatchKernel::setPartitioningInfoDirty(bool dirty)
{
	if (dirty && !isPartitioned()) {
		return;
	}

	m_partitioningInfoDirty = dirty;
}

/*!
	Update the information needed for ghost data exchange.

	\param forcedUpdated if set to true, ghost exchange information will be
	updated also if they are not marked as dirty
*/
void PatchKernel::updatePartitioningInfo(bool forcedUpdated)
{
	// Check if ghost information are dirty
	if (!forcedUpdated && !arePartitioningInfoDirty()) {
		return;
	}

	// Cell exchange data
	updateGhostCellExchangeInfo();

	// Vertex exchange data
	updateGhostVertexExchangeInfo();

	// Update patch owner
	updateOwner();

	// Partitioning information are now up-to-date
	setPartitioningInfoDirty(false);
}

/*!
	Update the information needed for exchanging data on ghost vertices.

	Since vertices are shared among elements, data exchange does not involve
	only the ranks identified for the target cells, but also all other ranks
	that have the source cells of this process among their ghosts.

	For example, in a patch like this:

	       V1    V2    V3    V4
	       +-----+-----+-----+            Vx identifies the vertex x
	       |  0  |  1  |  2  |
	       +-----+-----+-----+             x identifies the cell x
	       V5    V6    V7    V8

	partitioned as the following:

	       +-----+-----+
	       | I0  | G1  |        Rank 1
	       +-----+-----+
	                                      Gx identifies the internal cell x
	       +-----+-----+-----+
	       | G0  | I1  | G2  |  Rank 2
	       +-----+-----+-----+
	                                      Ix identifies the ghost cell x
	             +-----+-----+
	             | I1  | G2  |  Rank 3
	             +-----+-----+

	the rank 0 has to send data for vertices V2 and V6 both to rank 2 and
	rank 3.
*/
void PatchKernel::updateGhostVertexExchangeInfo()
{
	// Patch information
	int patchRank = getRank();

	// Idenfity owners of target and source vertices
	std::unordered_map<long, int> exchangeVertexOwners = evaluateExchangeVertexOwners();

	//
	// Clear ghosts
	//

	// List of vertices that are no more ghosts.
	//
	// Previous ghost vertices will be converted to internal vertices.
	std::vector<long> previousGhosts;
	previousGhosts.reserve(getGhostVertexCount());
	for (VertexConstIterator vertexItr = ghostVertexBegin(); vertexItr != ghostVertexEnd(); ++vertexItr) {
		long vertexId = vertexItr.getId();
		if (exchangeVertexOwners.count(vertexId) != 0) {
			continue;
		}

		previousGhosts.push_back(vertexId);
	}

	// Convert previous ghosts to internal vertices
	for (long vertexId : previousGhosts) {
		ghostVertex2InternalVertex(vertexId);
	}

	//
	// Initialize ghost
	//

	// List of new ghost vertices
	//
	// Vertices need to be sorted with their order in the storage.
	std::map<std::size_t, long> newGhosts;
	for (const auto &entry : exchangeVertexOwners) {
		int vertexOwner = entry.second;
		if (vertexOwner == patchRank) {
			continue;
		}

		long vertexId = entry.first;
		VertexIterator vertexItr = m_vertices.find(vertexId);
		if (!vertexItr->isInterior()) {
			continue;
		}

		newGhosts.insert({vertexItr.getRawIndex(), vertexId});
	}

	// Set the owner of the existing ghosts
	for (VertexConstIterator vertexItr = ghostVertexBegin(); vertexItr != ghostVertexEnd(); ++vertexItr) {
		long vertexId = vertexItr->getId();
		int vertexOwner = exchangeVertexOwners.at(vertexId);
		setGhostVertexOwner(vertexId, vertexOwner);
	}

	// Create new ghosts
	//
	// The list of ghosts is processed backwards to reduce the number of
	// vertices that need to be swapped (there is a high probability that
	// most of the ghost vertices are already at the end of the storage).
	//
	// If a vertex is already a ghost, we still need to update its owner.
	for (auto iter = newGhosts.rbegin(); iter != newGhosts.rend(); ++iter) {
		long vertexId = iter->second;
		int vertexOwner = exchangeVertexOwners.at(vertexId);
		internalVertex2GhostVertex(vertexId, vertexOwner);
	}

	//
	// Identify exchange targets
	//
	// See function Doxygen documentation for an explanation of which ranks are
	// involved in vertex data exchange.

	// Clear targets
	m_ghostVertexExchangeTargets.clear();

	// Update targets
	for (const auto &entry : m_ghostVertexOwners) {
		int ghostVertexRank = entry.second;
		long ghostVertexId = entry.first;
		m_ghostVertexExchangeTargets[ghostVertexRank].push_back(ghostVertexId);
	}

	// Sort the targets
	for (auto &entry : m_ghostVertexExchangeTargets) {
		std::vector<long> &rankTargets = entry.second;
		std::sort(rankTargets.begin(), rankTargets.end(), VertexPositionLess(*this));
	}

	//
	// Identify exchange sources
	//
	// See function Doxygen documentation for an explanation of which ranks are
	// involved in vertex data exchange.

	// Compute source distribution
	//
	// For each source cell, we need to know the ranks that have that cell
	// among their ghosts.
	std::unordered_map<long, std::unordered_set<int>> sourcesDistribution;
	for (const auto &entry : getGhostCellExchangeSources()) {
		const int rank = entry.first;
		const std::vector<long> &cellList = entry.second;
		for (long cellId : cellList) {
			sourcesDistribution[cellId].insert(rank);
		}
	}

	// Initialize data communicator
	DataCommunicator sourceDataCommunicator(getCommunicator());

	// Prepare the sends
	for (const auto &entry : getGhostCellExchangeSources()) {
		const int rank = entry.first;
		const std::vector<long> &cellList = entry.second;

		// Set the sends
		std::size_t bufferSize = 0;
		for (long cellId : cellList) {
			bufferSize += sizeof(int) + sizeof(int) * sourcesDistribution[cellId].size();
		}
		sourceDataCommunicator.setSend(rank, bufferSize);

		// Fill send buffer
		SendBuffer &buffer = sourceDataCommunicator.getSendBuffer(rank);
		for (long cellId : cellList) {
			std::unordered_set<int> &cellRankDistribution = sourcesDistribution.at(cellId);
			int cellRankDistributionSize = cellRankDistribution.size();

			buffer << static_cast<int>(cellRankDistributionSize);
			for (int sourceRank : cellRankDistribution) {
				buffer << sourceRank;
			}
		}
	}

	// Discover the receives
	sourceDataCommunicator.discoverRecvs();

	// Start the communications
	sourceDataCommunicator.startAllRecvs();
	sourceDataCommunicator.startAllSends();

	// Clear the sources
	m_ghostVertexExchangeSources.clear();

	// Initialize list of unique sources
	std::unordered_map<long, std::unordered_set<int>> uniqueSources;

	// Identifiy sources for ranks defined by source cells
	for (auto &entry : m_ghostCellExchangeSources) {
		int rank = entry.first;

		// Identify unique sources
		for (long cellId : entry.second) {
			const Cell &cell = getCell(cellId);
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			for (long vertexId : cellVertexIds) {
				// Ghost vertices are not sources.
				if (m_ghostVertexOwners.count(vertexId) > 0) {
					continue;
				}

				// Add the source to the list
				uniqueSources[rank].insert(vertexId);
			}
		}
	}

	// Identifiy sources for ranks defined by target cells
	//
	// For each ghost we receive the list of ranks that have that cell among
	// their ghosts, then we loop thorugh the vertices of the ghost cell and
	// we add the vertices owned by this process to the sources of the all
	// the ranks that contain the ghost cell.
	std::vector<int> cellRankDistribution;

	int nCompletedCellRecvs = 0;
	while (nCompletedCellRecvs < sourceDataCommunicator.getRecvCount()) {
		int rank = sourceDataCommunicator.waitAnyRecv();
		const std::vector<long> &cellList = getGhostCellExchangeTargets(rank);

		RecvBuffer &buffer = sourceDataCommunicator.getRecvBuffer(rank);
		for (long cellId : cellList) {
			int cellRankDistributionSize;
			buffer >> cellRankDistributionSize;
			cellRankDistribution.resize(cellRankDistributionSize);
			for (int i = 0; i < cellRankDistributionSize; ++i) {
				buffer >> cellRankDistribution[i];
			}

			const Cell &cell = getCell(cellId);
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			for (long vertexId : cellVertexIds) {
				// Ghost vertices are not sources.
				if (m_ghostVertexOwners.count(vertexId) > 0) {
					continue;
				}

				// Add the source to the list
				for (int cellRank : cellRankDistribution) {
					if (cellRank == patchRank) {
						continue;
					}

					uniqueSources[cellRank].insert(vertexId);
				}
			}
		}

		++nCompletedCellRecvs;
	}

	// Store and sort sources
	for (const auto &entry : uniqueSources) {
		int rank = entry.first;
		const std::unordered_set<int> &rankUniqueSources = entry.second;

		// Store sources
		std::vector<long> &rankSources = m_ghostVertexExchangeSources[rank];
		rankSources.assign(rankUniqueSources.begin(), rankUniqueSources.end());

		// Sort sources
		std::sort(rankSources.begin(), rankSources.end(), VertexPositionLess(*this));
	}
}

/*!
	Evaluate owners of exchange (target and source) vertices.

	Vertices that are on a partition border, belong to the partition with the
	lowest rank.

	We don't have information on the neighbours of the last layer of ghost
	cells. To identify the owner of the vertices on the border faces of the
	last layer of ghosts cells we need to exchange some information among
	the partitions. Each partition will use the local information to define
	the owner of the vertices of the source and target cells. This tentative
	values is then communicated among the partitions (the communications is
	done using the cell exchange information that have already been built).
	At this point each partition knows which are all the tentative owners
	of the vertices, the correct owner is the one with the lowest rank.
*/
std::unordered_map<long, int> PatchKernel::evaluateExchangeVertexOwners() const
{
	// Patch information
	int patchRank = getRank();

	// Initialize communicator
	DataCommunicator vertexOwnerCommunicator(getCommunicator());

	// Start receives
	for (const auto &entry : m_ghostCellExchangeTargets) {
		int rank = entry.first;
		const std::vector<long> &cellIds = entry.second;

		std::size_t bufferSize = 0;
		for (long cellId : cellIds) {
			const Cell &cell = getCell(cellId);
			bufferSize += sizeof(int) * cell.getVertexCount();
		}
		vertexOwnerCommunicator.setRecv(rank, bufferSize);
		vertexOwnerCommunicator.startRecv(rank);
	}

	// Initialize owners using local information
	std::unordered_map<long, int> exchangeVertexOwners;

	for (const auto &entry : m_ghostCellExchangeSources) {
		for (long cellId : entry.second) {
			const Cell &cell = getCell(cellId);
			for (long vertexId : cell.getVertexIds()) {
				auto exchangeVertexOwnerItr = exchangeVertexOwners.find(vertexId);
				if (exchangeVertexOwnerItr == exchangeVertexOwners.end()) {
					exchangeVertexOwners.insert({vertexId, patchRank});
				}
			}
		}
	}

	for (const auto &entry : m_ghostCellExchangeTargets) {
		int targetCellOwner = entry.first;
		for (long cellId : entry.second) {
			const Cell &cell = getCell(cellId);
			for (long vertexId : cell.getVertexIds()) {
				auto exchangeVertexOwnerItr = exchangeVertexOwners.find(vertexId);
				if (exchangeVertexOwnerItr == exchangeVertexOwners.end()) {
					exchangeVertexOwners.insert({vertexId, targetCellOwner});
				} else {
					int &currentVertexOwner = exchangeVertexOwnerItr->second;
					if (targetCellOwner < currentVertexOwner) {
						currentVertexOwner = targetCellOwner;
					}
				}
			}
		}
	}

	// Send local vertex owners to the neighbouring partitions
	for (const auto &entry : m_ghostCellExchangeSources) {
		int rank = entry.first;
		const std::vector<long> &cellIds = entry.second;

		std::size_t bufferSize = 0;
		for (long cellId : cellIds) {
			const Cell &cell = getCell(cellId);
			bufferSize += sizeof(int) * cell.getVertexCount();
		}
		vertexOwnerCommunicator.setSend(rank, bufferSize);

		SendBuffer &buffer = vertexOwnerCommunicator.getSendBuffer(rank);
		for (long cellId : cellIds) {
			const Cell &cell = getCell(cellId);
			for (long vertexId : cell.getVertexIds()) {
				buffer << exchangeVertexOwners.at(vertexId);
			}
		}
		vertexOwnerCommunicator.startSend(rank);
	}

	// Receive vertex owners from the neighbouring partitions
	int nCompletedRecvs = 0;
	while (nCompletedRecvs < vertexOwnerCommunicator.getRecvCount()) {
		int rank = vertexOwnerCommunicator.waitAnyRecv();
		const std::vector<long> &cellIds = m_ghostCellExchangeTargets.at(rank);
		RecvBuffer &buffer = vertexOwnerCommunicator.getRecvBuffer(rank);

		for (long cellId : cellIds) {
			const Cell &cell = getCell(cellId);
			for (long vertexId : cell.getVertexIds()) {
				int remoteVertexOwner;
				buffer >> remoteVertexOwner;

				int &currentVertexOwner = exchangeVertexOwners.at(vertexId);
				if (remoteVertexOwner < currentVertexOwner) {
					currentVertexOwner = remoteVertexOwner;
				}
			}
		}

		++nCompletedRecvs;
	}

	// Wait until all sends are complete
	vertexOwnerCommunicator.waitAllSends();

	return exchangeVertexOwners;
}

/*!
	Update the information needed for exchanging data on ghost cells.
*/
void PatchKernel::updateGhostCellExchangeInfo()
{
	// Clear targets
	m_ghostCellExchangeTargets.clear();

	// Update targets
	for (const auto &entry : m_ghostCellOwners) {
		int ghostRank = entry.second;
		long ghostCellId = entry.first;
		m_ghostCellExchangeTargets[ghostRank].push_back(ghostCellId);
	}

	// Sort the targets
	for (auto &entry : m_ghostCellExchangeTargets) {
		std::vector<long> &rankTargets = entry.second;
		std::sort(rankTargets.begin(), rankTargets.end(), CellPositionLess(*this));
	}

	// Clear the sources
	m_ghostCellExchangeSources.clear();

	// Build the sources
	for (auto &entry : m_ghostCellExchangeTargets) {
		int rank = entry.first;

		// Generate the source list
		std::vector<long> rankSources = _findGhostCellExchangeSources(rank);
		if (rankSources.empty()) {
			m_ghostCellExchangeSources.erase(rank);
			continue;
		}

		// Sort the sources
		std::sort(rankSources.begin(), rankSources.end(), CellPositionLess(*this));

		// Store list
		m_ghostCellExchangeSources[rank] = std::move(rankSources);
	}
}

/*!
	Finds the internal cells that will be ghost cells for the process
	with the specified ranks. During data exchange, these cells will be
	the sources form which data will be read from.

	\param rank is the rank for which the information will be built
*/
std::vector<long> PatchKernel::_findGhostCellExchangeSources(int rank)
{
	// Get targets for the specified rank
	//
	// If there are no targets, there will be no soruces either.
	auto ghostExchangeTargetsItr = m_ghostCellExchangeTargets.find(rank);
	if (ghostExchangeTargetsItr == m_ghostCellExchangeTargets.end()) {
		return std::vector<long>();
	}

	std::vector<long> &rankTargets = ghostExchangeTargetsItr->second;

	// The internal neighbours of the ghosts will be sources for the rank
	assert(getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	std::vector<long> neighIds;
	std::unordered_set<long> exchangeSources;
	exchangeSources.reserve(rankTargets.size());
	for (long ghostCellId : rankTargets) {
		neighIds.clear();
		findCellNeighs(ghostCellId, &neighIds);
		for (long neighId : neighIds) {
			if (m_ghostCellOwners.count(neighId) > 0) {
				continue;
			}

			exchangeSources.insert(neighId);
		}
	}

	return std::vector<long>(exchangeSources.begin(), exchangeSources.end());
}

}

#endif
