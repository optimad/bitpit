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
#if BITPIT_ENABLE_METIS==1
#include <metis.h>
#endif
#include <mpi.h>

#include <chrono>
#include <unordered_set>
#include <unordered_map>

#include "bitpit_communications.hpp"

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

	A patch is distributed among different processes if it doesn't have an owner.

	Setting the appropriate function argument, this function can be called also
	when the patch is not up-to-date. If dirty patches are allowed and the patch
	is actually dirty, the function will evaluate the owner on-the-fly. Otherwise
	the function will return the owner evaluated during the last update. In any
	case, if dirt patches are allowed, the function is a collective function and
	needs to be called by all processes (otherwise a deadlock will occur).

	\return Return true if the patch is distributed among different processes,
	false otherwise.
*/
bool PatchKernel::isDistributed(bool allowDirty) const
{
	return (getOwner(allowDirty) < 0);
}

/*!
	If the path is NOT distributed among different processes, returns the
	process that owns the patch, otherwise returns a negative number.

	Setting the appropriate function argument, this function can be called also
	when the patch is not up-to-date. If dirty patches are allowed and the patch
	is actually dirty, the function will evaluate the owner on-the-fly. Otherwise
	the function will return the owner evaluated during the last update. In any
	case, if dirt patches are allowed, the function is a collective function and
	needs to be called by all processes (otherwise a deadlock will occur).

	\param allowDirty if set to true, the function will evaluate the owner of a dirty
	patch on on-the-fly; otherwise the function will return the owner evaluated during
	the last updated, even if the patch is currently dirty. If dirty patch are allowed,
	the function is a collective function and needs to be called by all processes
	(otherwise a deadlock will occur)
	\return If the path is NOT distributed among different processes, returns
	the process that owns the patch, otherwise returns a negative number.
*/
int PatchKernel::getOwner(bool allowDirty) const
{
	if (allowDirty) {
		return evalOwner();
	}

	assert(!arePartitioningInfoDirty(false));

	return m_owner;
}

/*!
	Evaluate the owner of the patch.

	This function can be called also when the patch is not up-to-date. If the patch
	is up-to-date, the function will return the same result of PatchKernel::getOwner().

	If the path is NOT distributed among different processes, the owner is set
	to the process that owns the cells, otherwise the owner is set to a negative
	number.
*/
int PatchKernel::evalOwner() const
{
	long nInternalCells = getInternalCellCount();
	long nGlobalInternalCells = nInternalCells;
	if (isPartitioned()) {
		MPI_Allreduce(MPI_IN_PLACE, &nGlobalInternalCells, 1, MPI_LONG, MPI_SUM, getCommunicator());
	}

	int owner = -1;
	if (nGlobalInternalCells > 0) {
		if (nInternalCells == nGlobalInternalCells) {
			owner = getRank();
		}

		if (isPartitioned()) {
			MPI_Allreduce(MPI_IN_PLACE, &owner, 1, MPI_INT, MPI_MAX, getCommunicator());
		}
	}

	return owner;
}

/*!
	Updates the owner of the patch.

	If the path is NOT distributed among different processes, the owner is set
	to the process that owns the cells, otherwise the owner is set to a negative
	number.
*/
void PatchKernel::updateOwner()
{
	m_owner = evalOwner();
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

	if (isDistributed()) {
		throw std::runtime_error ("Halo size can only be set when the patch is all on a single process");
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
	\param[in] owner is the rank of the process that owns the ghost cell
*/
PatchKernel::VertexIterator PatchKernel::internalVertex2GhostVertex(long id, int owner)
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
	setGhostVertexInfo(id, owner);

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
	unsetGhostVertexInfo(id);

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
	\param owner is the rank that owns the vertex that will be restored
	\param id is the id of the vertex that will be restored
	\return An iterator pointing to the restored vertex.
*/
PatchKernel::VertexIterator PatchKernel::restoreVertex(const std::array<double, 3> &coords, int owner, long id)
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
	if (owner == getRank()) {
		_restoreInternalVertex(iterator, coords);
	} else {
		_restoreGhostVertex(iterator, coords, owner);
	}

	return iterator;
}

/*!
	Internal function to restore a ghost vertex.

	The kernel should already contain the vertex, only the contents of the
	vertex will be updated.

	\param iterator is an iterator pointing to the vertex to restore
	\param coords are the coordinates of the vertex
	\param owner is the rank that owns the vertex that will be restored
*/
void PatchKernel::_restoreGhostVertex(const VertexIterator &iterator, const std::array<double, 3> &coords, int owner)
{
	// Restore the vertex
	//
	// There is no need to set the id of the vertex as assigned, because
	// also the index generator will be restored.
	Vertex &vertex = *iterator;
	vertex.initialize(iterator.getId(), coords, false);
	m_nGhostVertices++;

	// Update the bounding box
	addPointToBoundingBox(vertex.getCoords());

	// Set owner
	setGhostVertexInfo(vertex.getId(), owner);
}

/*!
	Internal function to delete a ghost vertex.

	\param id is the id of the vertex
*/
void PatchKernel::_deleteGhostVertex(long id)
{
	// Unset ghost owner
	unsetGhostVertexInfo(id);

	// Update the bounding box
	const Vertex &vertex = m_vertices[id];
	removePointFromBoundingBox(vertex.getCoords());

	// Delete vertex
	m_vertices.erase(id, true);
	m_nGhostVertices--;
	if (id == m_firstGhostVertexId) {
		updateFirstGhostVertexId();
	}

    // Vertex id is no longer used
    if (m_vertexIdGenerator) {
        m_vertexIdGenerator->trash(id);
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
	\param[in] owner is the rank of the process that owns the ghost cell
	\param[in] haloLayer is the halo layer the ghost cell belongs to
*/
PatchKernel::CellIterator PatchKernel::internalCell2GhostCell(long id, int owner, int haloLayer)
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

	// Set ghost information
	setGhostCellInfo(id, owner, haloLayer);

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

	// Unset ghost information
	unsetGhostCellInfo(id);

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
	\param owner is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, int owner, long id)
{
	return addCell(source, owner, 0, id);
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param owner is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, int owner, long id)
{
	return addCell(source, owner, 0, id);
}

/*!
	Adds a new cell with the specified id and type.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param type is the type of the cell
	\param owner is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, int owner, long id)
{
	return addCell(type, owner, 0, id);
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
	\param owner is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, const std::vector<long> &connectivity,
											   int owner, long id)
{
	return addCell(type, connectivity, owner, 0, id);
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
	\param owner is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
											   int owner, long id)
{
	return addCell(type, std::move(connectStorage), owner, 0, id);
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, int owner, int haloLayer, long id)
{
	Cell cell = source;
	cell.setId(id);

	return addCell(std::move(cell), owner, haloLayer, id);
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, int owner, int haloLayer, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	int connectSize = source.getConnectSize();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	if (!source.hasInfo()){
		std::copy(source.getConnect(), source.getConnect() + connectSize, connectStorage.get());
	}

	CellIterator iterator = addCell(source.getType(), std::move(connectStorage), owner, haloLayer, id);

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
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, int owner, int haloLayer, long id)
{
	std::unique_ptr<long[]> connectStorage;
	if (ReferenceElementInfo::hasInfo(type)) {
		int connectSize = ReferenceElementInfo::getInfo(type).nVertices;
		connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	} else {
		connectStorage = std::unique_ptr<long[]>(nullptr);
	}

	return addCell(type, std::move(connectStorage), owner, haloLayer, id);
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
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, const std::vector<long> &connectivity,
											   int owner, int haloLayer, long id)
{
	int connectSize = connectivity.size();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	std::copy(connectivity.data(), connectivity.data() + connectSize, connectStorage.get());

	return addCell(type, std::move(connectStorage), owner, haloLayer, id);
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
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
											   int owner, int haloLayer, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	CellIterator iterator;
	if (owner == getRank()) {
		iterator = _addInternalCell(type, std::move(connectStorage), id);
	} else {
		iterator = _addGhostCell(type, std::move(connectStorage), owner, haloLayer, id);
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
	\param owner is the rank that owns the cell that will be added
	\param haloLayer is the halo layer the ghost cell belongs to
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::_addGhostCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
                                                     int owner, int haloLayer, long id)
{
	// Get the id of the cell
	if (m_cellIdGenerator) {
		if (id < 0) {
			id = m_cellIdGenerator->generate();
		} else {
			m_cellIdGenerator->setAssigned(id);
		}
	} else if (id < 0) {
		throw std::runtime_error("No valid id has been provided for the cell.");
	}

	// Create the cell
	//
	// If there are internal cells, the ghost cell should be inserted
	// after the last internal cell.
	bool storeInterfaces  = (getInterfacesBuildStrategy() != INTERFACES_NONE);
	bool storeAdjacencies = storeInterfaces || (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	CellIterator iterator;
	if (m_lastInternalCellId < 0) {
		iterator = m_cells.emreclaim(id, id, type, std::move(connectStorage), false, storeInterfaces, storeAdjacencies);
	} else {
		iterator = m_cells.emreclaimAfter(m_lastInternalCellId, id, id, type, std::move(connectStorage), false, storeInterfaces, storeAdjacencies);
	}
	m_nGhostCells++;

	// Update the id of the first ghost cell
	if (m_firstGhostCellId < 0) {
		m_firstGhostCellId = id;
	} else if (m_cells.rawIndex(m_firstGhostCellId) > m_cells.rawIndex(id)) {
		m_firstGhostCellId = id;
	}

	// Set ghost information
	setGhostCellInfo(id, owner, haloLayer);

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
	\param owner is the rank that owns the cell that will be restored
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
	\param id is the id of the cell that will be restored
	\return An iterator pointing to the restored cell.
*/
PatchKernel::CellIterator PatchKernel::restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												   int owner, int haloLayer, long id)
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
	if (owner == getRank()) {
		_restoreInternalCell(iterator, type, std::move(connectStorage));
	} else {
		_restoreGhostCell(iterator, type, std::move(connectStorage), owner, haloLayer);
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
	\param owner is the rank that owns the cell that will be restored
	\param haloLayer is the halo layer the cell belongs to, this argument is only relevant
	if the cell is a ghost
*/
void PatchKernel::_restoreGhostCell(const CellIterator &iterator, ElementType type,
                                    std::unique_ptr<long[]> &&connectStorage,
                                    int owner, int haloLayer)
{
	// Restore the cell
	//
	// There is no need to set the id of the cell as assigned, because
	// also the index generator will be restored.
	bool storeInterfaces  = (getInterfacesBuildStrategy() != INTERFACES_NONE);
	bool storeAdjacencies = storeInterfaces || (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	long cellId = iterator.getId();
	Cell &cell = *iterator;
	cell.initialize(iterator.getId(), type, std::move(connectStorage), false, storeInterfaces, storeAdjacencies);
	m_nGhostCells++;

	// Set ghost information
	setGhostCellInfo(cellId, owner, haloLayer);

	// Set the alteration flags of the cell
	setRestoredCellAlterationFlags(cellId);
}

/*!
	Internal function to delete a ghost cell.

	\param id is the id of the cell
*/
void PatchKernel::_deleteGhostCell(long id)
{
	// Unset ghost information
	unsetGhostCellInfo(id);

	// Set the alteration flags of the cell
	setDeletedCellAlterationFlags(id);

	// Delete cell
	m_cells.erase(id, true);
	m_nGhostCells--;
	if (id == m_firstGhostCellId) {
		updateFirstGhostCellId();
	}

	// Cell id is no longer used
	if (m_cellIdGenerator) {
		m_cellIdGenerator->trash(id);
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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	Optionally, this funciton can track the changes performed to the patch. See
	PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
	for the documentation about the tracking information returned by this
	function.

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

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information returned by this function.

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

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information returned by this function.

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

	std::vector<adaption::Info> partitioningData = _partitioningPrepare(cellRanks, trackPartitioning);

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Partitions the patch among the processes. The partitioning is done using
	a criteria that tries to balance the load among the processes.

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information returned by this function.

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

	Tracking will only list changes that involve data exchange among internal
	data structured (e.g., partitioning preparation will not track ghost cells
	nor cell deletion).

	Information available on the sender side for tracking purposes are the
	following:
	 - internal cells that will be send;
	 - internal vertices that will be send.

	No information about tracking are provided on the receiver side.

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

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information returned by this function.

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

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information returned by this function.

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

	Optionally, this funciton can track the changes performed to the patch.

	Tracking information will not contain changes that involve data exchange
	among ghost data structured (e.g., ghost cells exchanges are not tracked,
	only deletion and creation of ghost cells will be tracked).

	Multiple adaption information can be associated to a single element, for
	example a cell can be sent to another processor and then "recreated" to
	become a ghost cell.

	Information available on the sender side for tracking purposes are the
	following:
	 - internal cells that have been sent;
	 - internal vertices that have been sent;
	 - new ghost cells that have been created (some of the internal cells that
	   have been sent may have become ghosts cells);
	 - new ghost vertices that have been created (some of the internal vertices
	   that have been sent may have become ghosts vertices);
	 - ghost cells that have been deleted;
	 - ghost vertices that have been deleted;
	 - new interfaces that have been created;
	 - interfaces that have been deleted.

	Information available on the receiver side for tracking purposes are the
	following:
	 - internal cells that have been received;
	 - internal vertices that have been received;
	 - new ghost cells that have been created;
	 - new ghost vertices that have been created;
	 - ghost cells that have been deleted (some ghost cells may have been
	   replaced by internal cells that have just been received);
	 - ghost vertices that have been deleted (some ghost vertices may have been
	   replaced by internal vertices that have just been received);
	 - new interfaces that have been created;
	 - interfaces that have been deleted.

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
	mergeAdaptionInfo(_partitioningAlter(trackPartitioning), partitioningData);

	// Finalize patch alterations
	mergeAdaptionInfo(finalizeAlterations(trackPartitioning, squeezeStorage), partitioningData);

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
	Checks if the patch is partitioned.

	A patch is considered partitioned if its MPI communicator spans multiple
	processes.

	\result Returns true if the patch is partitioned, false otherwise.
*/
bool PatchKernel::isPartitioned() const
{
	return (getProcessorCount() > 1);
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
		CellConstIterator beginItr = internalCellConstBegin();
		CellConstIterator endItr   = internalCellConstEnd();

		partitionWeight = 0.;
		for (CellConstIterator cellItr = beginItr; cellItr != endItr; ++cellItr) {
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

	See PatchKernel::partitioningPrepare(bool trackPartitioning) for the
	documentation about the tracking information that should be returned
	when re-implmeneting by this function.
*/
#if BITPIT_ENABLE_METIS==1
/*!
	This function uses METIS to evaluate patch partitioning, hence it can only be
	used when the patch is all on a single process.

	METIS expects weights to be integer numbers, non-integer weights are rounded to
	the nearest integer and weights less than 1 are set to 1.
*/
#else
/*!
	Default implementation is a no-op function and can only be used if the patch is
	empty or if the patch is not distributed and a single partition is requested.
*/
#endif
/*!
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
	// Default partitioning can only be used when the patch is all on a single process
	if (isDistributed()) {
		throw std::runtime_error("Default partitioning can only be used when the patch is all on a single process");
	}

	// Early return if the mesh is empty
	if (empty(true)) {
		return std::vector<adaption::Info>();
	}

	// Early return if there is a single partition
	if (getProcessorCount() == 1) {
		return std::vector<adaption::Info>();
	}

	// Evaluate partitioning
#if BITPIT_ENABLE_METIS==1
	std::unordered_map<long, int> cellRanks;
	if (getRank() == getOwner()) {
		// Initialize map between patch vertices and METIS nodes
		//
		// METIS node numbering needs to be contiguous and needs to start from zero.
		std::unordered_map<long, idx_t> vertexToNodeMap;

		idx_t nodeId = 0;
		for (const Vertex &vertex : m_vertices) {
			vertexToNodeMap[vertex.getId()] = nodeId;
			++nodeId;
		}

		// Create METIS mesh
		idx_t nElements = getInternalCellCount();
		idx_t nNodes = getInternalVertexCount();

		idx_t connectSize = 0;
		for (const Cell &cell : m_cells){
			connectSize += static_cast<idx_t>(cell.getVertexCount());
		}

		std::vector<idx_t> connectStorage(connectSize);
		std::vector<idx_t> connectRanges(nElements + 1);
		std::vector<idx_t> elementWeights(nElements, static_cast<idx_t>(std::max(defaultWeight, 1.)));

		idx_t i = 0;
		idx_t j = 0;
		connectRanges[0] = 0;
		for (const Cell &cell : m_cells) {
			// Element connectivity
			for (long vertex : cell.getVertexIds()) {
				connectStorage[j] = vertexToNodeMap[vertex];
				++j;
			}
			connectRanges[i + 1] = j;

			// Element weight
			auto weightItr = cellWeights.find(cell.getId());
			if (weightItr != cellWeights.end()) {
				elementWeights[i] = static_cast<idx_t>(std::max(std::round(weightItr->second), 1.));
			}

			// Increase counters
			i++;
		}

		// Get the number of common nodes two elements should share to be considered neighbours
		idx_t nCommonNodes = getDimension();

		// Get the number of partitions the mesh should be split into
		idx_t nPartitions = getProcessorCount();

		// Evaluate partitioning
		idx_t objectiveValue;

		std::vector<idx_t> elementRanks(nElements);
		std::vector<idx_t> nodeRanks(nNodes);

		int status = METIS_PartMeshDual(&nElements, &nNodes, connectRanges.data(), connectStorage.data(),
										elementWeights.data(), nullptr, &nCommonNodes, &nPartitions, nullptr,
										nullptr, &objectiveValue, elementRanks.data(), nodeRanks.data());

		if (status != METIS_OK) {
			throw std::runtime_error("Error during partitioning (METIS error " + std::to_string(status) + ")");
		}

		// Fill the cell ranks
		int metisId = 0;
		for (const Cell &cell : m_cells) {
			int metisRank = elementRanks[metisId];
			if (metisRank != getRank()) {
				cellRanks[cell.getId()] = elementRanks[metisId];
			}
			++metisId;
		}
	}

	return _partitioningPrepare(cellRanks, trackPartitioning);
#else
	BITPIT_UNUSED(cellWeights);
	BITPIT_UNUSED(defaultWeight);
	BITPIT_UNUSED(trackPartitioning);

	throw std::runtime_error("METIS library is required for automatic patch partitioning.");
#endif
}

/*!
	Partitions the patch among the processes. Each cell will be assigned
	to a specific process according to the specified input.

	See PatchKernel::partitioningPrepare(bool trackPartitioningge) for the
	documentation about the tracking information returned by this function.

	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::_partitioningPrepare(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning)
{
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
		if (m_ghostCellInfo.count(cellId) > 0) {
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

	// Track changes that involve data exchange
	std::vector<adaption::Info> partitioningData;
	if (trackPartitioning) {
		std::vector<long> cellsToSend;
		for (const std::pair<int, int> &entry : m_partitioningGlobalExchanges) {
			// Get send/receive ranks
			int sendRank = entry.first;
			if (sendRank != patchRank) {
				continue;
			}

			// Get list of cells to be sent
			int recvRank = entry.second;

			cellsToSend.clear();
			for (const auto &entry : m_partitioningOutgoings) {
				int cellRank = entry.second;
				if (cellRank != recvRank) {
					continue;
				}

				long cellId = entry.first;
				cellsToSend.push_back(cellId);
			}

			if (cellsToSend.empty()) {
				continue;
			}

			// Fill tracking data structures
			partitioningData.emplace_back();
			adaption::Info &partitioningVertexInfo = partitioningData.back();
			partitioningVertexInfo.entity   = adaption::ENTITY_VERTEX;
			partitioningVertexInfo.type     = adaption::TYPE_PARTITION_SEND;
			partitioningVertexInfo.rank     = recvRank;
			partitioningVertexInfo.previous = getOrderedCellsVertices(cellsToSend, true, false);

			partitioningData.emplace_back();
			adaption::Info &partitioningCellInfo = partitioningData.back();
			partitioningCellInfo.entity   = adaption::ENTITY_CELL;
			partitioningCellInfo.type     = adaption::TYPE_PARTITION_SEND;
			partitioningCellInfo.rank     = recvRank;
			partitioningCellInfo.previous = std::move(cellsToSend);
		}
	}

	return partitioningData;
}

/*!
    Alter the patch performing the partitioning.

    See PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
    for the documentation about the tracking information returned by this function.

    \param trackPartitioning if set to true the function will return the changes
    done to the patch during the partitioning
    \result If the partitioning is tracked, returns a vector of adaption::Info
    with all the changes done to the patch during the partitioning, otherwise
    an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter(bool trackPartitioning)
{
    std::vector<adaption::Info> partitioningData;

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
        mergeAdaptionInfo(_partitioningAlter_deleteGhosts(trackPartitioning), partitioningData);
    }

    // Communicate patch data structures
    int patchRank = getRank();

    m_partitioningCellsTag = communications::tags().generate(m_communicator);
    m_partitioningVerticesTag = communications::tags().generate(m_communicator);

    std::unordered_set<int> batchSendRanks;
    std::unordered_set<int> batchRecvRanks;

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

            if (finalOwner != m_ghostCellInfo.at(id).owner) {
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
        const GhostCellInfo &ghostCellInfo = m_ghostCellInfo.at(ghostCellId);
        int previousGhostCellOwner = ghostCellInfo.owner;
        if (previousGhostCellOwner != sendRank) {
            ++itr;
            continue;
        }

        // Update ghost owner
        int finalGhostCellOwner = itr->second;
        setGhostCellInfo(ghostCellId, finalGhostCellOwner, ghostCellInfo.haloLayer);
        itr = ghostCellOwnershipChanges->erase(itr);
    }
}

/*!
    Delete ghosts.

    See PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
    for the documentation about the tracking information returned by this function.

    \param trackPartitioning if set to true the function will return the changes
    done to the patch during the partitioning
    \result If the partitioning is tracked, returns a vector of adaption::Info
    with all the changes done to the patch during the partitioning, otherwise
    an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter_deleteGhosts(bool trackPartitioning)
{
    std::vector<adaption::Info> partitioningData;
    if (getGhostCellCount() == 0) {
        return partitioningData;
    }

    // Track changes
    if (trackPartitioning) {
        partitioningData.emplace_back();
        adaption::Info &partitioningVertexInfo = partitioningData.back();
        partitioningVertexInfo.entity= adaption::ENTITY_VERTEX;
        partitioningVertexInfo.type  = adaption::TYPE_DELETION;
        partitioningVertexInfo.current.reserve(getGhostVertexCount());
        for (VertexConstIterator itr = ghostVertexConstBegin(); itr != ghostVertexConstEnd(); ++itr) {
            partitioningVertexInfo.current.push_back(itr.getId());
        }

        partitioningData.emplace_back();
        adaption::Info &partitioningCellInfo = partitioningData.back();
        partitioningCellInfo.entity = adaption::ENTITY_CELL;
        partitioningCellInfo.type   = adaption::TYPE_DELETION;
        partitioningCellInfo.current.reserve(getGhostCellCount());
        for (CellConstIterator itr = ghostCellConstBegin(); itr != ghostCellConstEnd(); ++itr) {
            partitioningCellInfo.current.push_back(itr.getId());
        }
    }

    // Delete ghost cells
    std::vector<long> cellsDeleteList;
    cellsDeleteList.reserve(m_ghostCellInfo.size());
    for (const auto &entry : m_ghostCellInfo) {
        long cellId = entry.first;
        cellsDeleteList.emplace_back(cellId);
    }

    deleteCells(cellsDeleteList);

    // Prune stale adjacencies
    pruneStaleAdjacencies();

    // Prune stale interfaces
    //
    // Interfaces will be re-created after the partitioning.
    mergeAdaptionInfo(pruneStaleInterfaces(trackPartitioning), partitioningData);

    // Delete vertices no longer used
    deleteOrphanVertices();

    // Convert all ghost vertices to internal vertices
    while (getGhostVertexCount() > 0) {
        ghostVertex2InternalVertex(m_firstGhostVertexId);
    }

    return partitioningData;
}

/*!
    Sends the given list of cells to the process with the specified rank.

    See PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
    for the documentation about the tracking information returned by this function.

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
    std::unordered_set<ConstCellHalfEdge, ConstCellHalfEdge::Hasher> frontierEdges;

    std::unordered_set<long> outgoingCells;
    std::unordered_set<long> frameCells;
    std::unordered_map<long, std::size_t> haloCells;
    std::vector<long> firstHaloLayerCells;

    std::unordered_set<long> frameVertices;
    std::unordered_set<long> haloVertices;

    std::vector<long> cellSendList;
    std::unordered_set<long> vertexSendList;

    std::vector<long> frameCellsOverall;
    std::unordered_set<long> ghostHaloCellsOverall;
    std::unordered_set<long> innerFrontierCellsOverall;

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
        // The first layer of halo cells is made up by cells not explicitly
        // marked for being sent to the receiver that have a face, an edge
        // or a vertex on a frontier face. Cells on subsequent halo layers
        // are identifies as cells not explicitly marked for being sent
        // that are neighbours of the the previous halo layer. At this
        // stage we are only building the first layer of halo cells.
        //
        // We also identify outgoing cells that have a face/edge/vertex on
        // the inner frontier. A face/edge/vertex is on the inner frontier
        // if its on the frontier and one of the cells that contain that
        // item is not an outgoing cell (i.e., is not explicitly marked for
        // being sent to any rank). Cells on inner frontier will belong to
        // the first layer of ghost cells for the current rank.
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
                        const Cell *frontierCell;
                        if (nFaceAdjacencies == 1) {
                            frontierCellId = cellId;
                            frontierFace = face;
                            frontierCell = &cell;
                        } else {
                            const Cell &neigh = getCell(neighId);

                            frontierCellId = neighId;
                            frontierFace = findAdjoinNeighFace(cell, face, neigh);
                            frontierCell = &neigh;

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
                            innerFrontierCellsOverall.insert(cellId);
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
                                        break;
                                    }
                                }
                            }

                            // Add neighbours to the list
                            for (long frontierNeighId : neighsList) {
                                frontierNeighs.insert(frontierNeighId);
                                if (innerFrontierVertex) {
                                    if (m_partitioningOutgoings.count(frontierNeighId) > 0) {
                                        innerFrontierCellsOverall.insert(frontierNeighId);
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
                                int nEdgeVertices = frontierCell->getEdgeVertexCount(edge);

                                ConstCellHalfEdge frontierEdge = ConstCellHalfEdge(*frontierCell, edge, ConstCellHalfEdge::WINDING_NATURAL);

                                // Discard edges that do not belong to the
                                // current frontier face
                                ConstProxyVector<long> edgeVertexIds = frontierEdge.getVertexIds();

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
                                frontierEdge.setWinding(ConstCellHalfEdge::WINDING_REVERSE);
                                auto frontierEdgeItr = frontierEdges.find(frontierEdge);
                                if (frontierEdgeItr != frontierEdges.end()) {
                                    continue;
                                } else {
                                    frontierEdge.setWinding(ConstCellHalfEdge::WINDING_NATURAL);
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
                                            break;
                                        }
                                    }
                                }

                                // Add neighbours to the list
                                for (long frontierNeighId : neighsList) {
                                    frontierNeighs.insert(frontierNeighId);
                                    if (innerFrontierEdge) {
                                        if (m_partitioningOutgoings.count(frontierNeighId) > 0) {
                                            innerFrontierCellsOverall.insert(frontierNeighId);
                                        }
                                    }
                                }
                            }
                        }

                        // Tell apart frame and halo cells
                        for (long frontierNeighId : frontierNeighs) {
                            if (outgoingCells.count(frontierNeighId) > 0) {
                                frameCells.insert(frontierNeighId);
                            } else {
                                haloCells.insert({frontierNeighId, 0});
                            }
                        }
                    }
                }
            }

            // Identify cells on subsequent frame layers
            //
            // Only cells sent to the receiver rank can be frame cells.
            if (m_haloSize > 1) {
                auto frameSelector = [&outgoingCells](long cellId) {
                    return (outgoingCells.count(cellId) != 0);
                };

                auto frameBuilder = [&frameCells](long cellId, int layer) {
                    BITPIT_UNUSED(layer);

                    frameCells.insert(cellId);

                    return false;
                };

                processCellsNeighbours(frameCells, m_haloSize - 1, frameSelector, frameBuilder);
            }

            // Identify cells on subsequent halo layers
            //
            // The selector should discard cells in the frame, this ensures that
            // the halo will not extend to the cells being sent to the receiver
            // rank.
            if (m_haloSize > 1) {
                firstHaloLayerCells.clear();
                firstHaloLayerCells.reserve(haloCells.size());
                for (const auto &haloEntry : haloCells) {
                    firstHaloLayerCells.push_back(haloEntry.first);
                }

                auto haloSelector = [&frameCells](long cellId) {
                    return (frameCells.count(cellId) == 0);
                };

                auto haloBuilder = [&haloCells](long cellId, int layer) {
                    haloCells.insert({cellId, layer + 1});

                    return false;
                };

                processCellsNeighbours(firstHaloLayerCells, m_haloSize - 1, haloSelector, haloBuilder);
            }
        }

        // Order cells to send
        //
        // To allow the receiver to efficently store the cells, first we
        // send the cells that are internal for the receiver, then the
        // halo cells.
        std::size_t nOutgoingCells = outgoingCells.size();
        std::size_t nHaloCells     = haloCells.size();

        cellSendList.assign(outgoingCells.begin(), outgoingCells.end());

        std::size_t haloIndex = cellSendList.size();
        cellSendList.resize(cellSendList.size() + nHaloCells);
        for (const auto &haloEntry : haloCells) {
            cellSendList[haloIndex] = haloEntry.first;
            ++haloIndex;
        }

        // Update list of overall frame cells
        frameCellsOverall.insert(frameCellsOverall.end(), frameCells.begin(), frameCells.end());

        // Update list of halo cells
        for (const auto &haloEntry : haloCells) {
            long cellId = haloEntry.first;
            const Cell &cell = getCell(cellId);
            if (cell.isInterior()) {
                continue;
            }

            ghostHaloCellsOverall.insert(cellId);
        }

        // Track changes that involve data exchange
        if (trackPartitioning && !cellSendList.empty()) {
            // Track cells that have been sent to other processes
            //
            // The ids of the cells send will be stored accordingly to the send
            // order, this is the same order that will be used on the process
            // that has received the cell. Since the order is the same, the
            // two processes are able to exchange cell data without additional
            // communications (they already know the list of cells for which
            // data is needed and the order in which these data will be sent).
            partitioningData.emplace_back();
            adaption::Info &partitioningVertexInfo = partitioningData.back();
            partitioningVertexInfo.entity  = adaption::ENTITY_VERTEX;
            partitioningVertexInfo.type    = adaption::TYPE_PARTITION_SEND;
            partitioningVertexInfo.rank    = recvRank;
            partitioningVertexInfo.current = getOrderedCellsVertices(cellSendList, true, false);

            partitioningData.emplace_back();
            adaption::Info &partitioningCellInfo = partitioningData.back();
            partitioningCellInfo.entity   = adaption::ENTITY_CELL;
            partitioningCellInfo.type     = adaption::TYPE_PARTITION_SEND;
            partitioningCellInfo.rank     = recvRank;
            partitioningCellInfo.previous = std::vector<long>(cellSendList.begin(), cellSendList.begin() + nOutgoingCells);
        }

        //
        // Communicate cell buffer size
        //

        // Wait for previous comunications to finish
        if (cellsSizeRequest != MPI_REQUEST_NULL) {
            MPI_Wait(&cellsSizeRequest, MPI_STATUS_IGNORE);
        }

        // Start the communication
        long cellsBufferSize = 2 * sizeof(long) + 3 * nHaloCells * sizeof(int) + 2 * (nOutgoingCells + nHaloCells) * sizeof(bool);
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
            // Vertex information
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
                } else if (m_ghostCellInfo.count(cellId) > 0) {
                    cellOwnerOnReceiver = m_ghostCellInfo.at(cellId).owner;
                } else {
                    cellOwnerOnReceiver = m_rank;
                }

                cellsBuffer << cellOwnerOnReceiver;

                int cellHaloLayerOnReceiver;
                if (cellOwnerOnReceiver != recvRank) {
                    cellHaloLayerOnReceiver = haloCells.at(cellId);
                } else {
                    cellHaloLayerOnReceiver = -1;
                }

                cellsBuffer << cellHaloLayerOnReceiver;

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

        // Delete outgoing cells not in the frame
        //
        // Frame cells cannot be deleted just now, because they may be ghost cells
        // for other processes.
        std::vector<long> deleteList;

        deleteList.reserve(nOutgoingCells - frameCells.size());
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

        mergeAdaptionInfo(pruneStaleInterfaces(trackPartitioning), partitioningData);

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
        // Convert inner cells into ghosts
        //
        // All cells on the inner frontier should become ghost cells. These
        // cells belongs to the first halo layer and defines the seeds for
        // growing the subsequent halo layers.
        std::vector<long> trackedCreatedGhostCells;
        for (long cellId : innerFrontierCellsOverall) {
            const Cell &cell = getCell(cellId);
            std::size_t cellHaloLayer = 0;
            if (!confirmCellHaloLayer(cell, cellHaloLayer, m_partitioningOutgoings)) {
                continue;
            }

            int cellOwner = m_partitioningOutgoings.at(cellId);
            internalCell2GhostCell(cellId, cellOwner, cellHaloLayer);
            if (trackPartitioning) {
                trackedCreatedGhostCells.push_back(cellId);
            }
        }

        auto ghostFrameSelector = [this](long cellId) {
            return (m_partitioningOutgoings.count(cellId) != 0);
        };

        auto ghostFrameBuilder = [this, trackPartitioning, &trackedCreatedGhostCells](long cellId, int layer) {
            int cellOwner = m_partitioningOutgoings.at(cellId);
            internalCell2GhostCell(cellId, cellOwner, layer + 1);
            if (trackPartitioning) {
                trackedCreatedGhostCells.push_back(cellId);
            }

            return false;
        };

        processCellsNeighbours(innerFrontierCellsOverall, m_haloSize - 1, ghostFrameSelector, ghostFrameBuilder);
        if (!trackedCreatedGhostCells.empty()) {
            partitioningData.emplace_back();
            adaption::Info &cellCreationInfo = partitioningData.back();
            cellCreationInfo.entity  = adaption::ENTITY_CELL;
            cellCreationInfo.type    = adaption::TYPE_CREATION;
            cellCreationInfo.rank    = getRank();
            cellCreationInfo.current = getOrderedCellsVertices(trackedCreatedGhostCells, false, true);

            std::vector<long> deletedGhostVertices = getOrderedCellsVertices(trackedCreatedGhostCells, false, true);
            if (!deletedGhostVertices.empty()) {
                partitioningData.emplace_back();
                adaption::Info &vertexCreationInfo = partitioningData.back();
                vertexCreationInfo.entity  = adaption::ENTITY_VERTEX;
                vertexCreationInfo.type    = adaption::TYPE_CREATION;
                vertexCreationInfo.rank    = getRank();
                vertexCreationInfo.current = std::move(deletedGhostVertices);
            }
        }

        // Delete frame cells that are not ghosts
        //
        // Now that the new ghosts have been created, we can delete all the
        // frontier cells that are not ghosts.
        std::vector<long> deleteList;
        for (long cellId : frameCellsOverall) {
            const Cell &cell = getCell(cellId);
            if (!cell.isInterior()) {
                continue;
            }

            deleteList.push_back(cellId);
            ghostCellOwnershipChanges->erase(cellId);
        }

        deleteCells(deleteList);

        // Prune cell adjacencies
        //
        // At this stage we cannot fully update the adjacencies, but we need
        // to remove stale adjacencies.
        pruneStaleAdjacencies();

        // Delete or update ghosts
        //
        // Some of the halo cells may be stale ghosts. Ghosts in the first halo
        // layer are considered stale if they have has no inner neighbours. Ghost
        // on subsequent layers are considered stale if they have no neighbours
        // in the previous halo layer.
        //
        // A layer at a time, stale ghost cells are found and deleted, valid ghosts
        // are identified because we need to re-compute the layer layer associated
        // with those ghosts.
        std::vector<long> neighIds;
        std::vector<long> updateList;

        std::vector<long> haloDeleteList;
        std::vector<long> trackedDeletedGhostCells;
        for (int haloLayer = 0; haloLayer < static_cast<int>(m_haloSize); ++haloLayer) {
            haloDeleteList.clear();
            for (long cellId : ghostHaloCellsOverall) {
                // Consider only cells in the current halo layer
                const GhostCellInfo &ghostInfo = m_ghostCellInfo.at(cellId);
                if (ghostInfo.haloLayer != haloLayer) {
                    continue;
                }

                // Cells for which the halo layer is not confirmed will be delete, cells
                // for which the halo layer is confirmed will be updated.
                const Cell &cell = getCell(cellId);
                if (confirmCellHaloLayer(cell, haloLayer, m_partitioningOutgoings)) {
                    updateList.push_back(cellId);
                } else {
                    haloDeleteList.push_back(cellId);
                }
            }

            // Delete stale ghost in this layer
            for (long cellId : haloDeleteList) {
                ghostCellOwnershipChanges->erase(cellId);
                ghostHaloCellsOverall.erase(cellId);
            }
            deleteCells(haloDeleteList);

            // Track deleted cells
            if (trackPartitioning) {
                trackedDeletedGhostCells.insert(trackedDeletedGhostCells.end(), haloDeleteList.begin(), haloDeleteList.end());
            }

            // Prune stale adjacencies
            pruneStaleAdjacencies();
        }

        // Track ghost cell deletion
        if (trackPartitioning && !trackedDeletedGhostCells.empty()) {
            partitioningData.emplace_back();
            adaption::Info &cellDeletionInfo = partitioningData.back();
            cellDeletionInfo.entity  = adaption::ENTITY_CELL;
            cellDeletionInfo.type    = adaption::TYPE_DELETION;
            cellDeletionInfo.rank    = getRank();
            cellDeletionInfo.current = std::move(trackedDeletedGhostCells);
            trackedDeletedGhostCells.clear();
        }

        // Compute layer associated with ghosts in the halo of the outgoing cells
        //
        // Some cells have been deleted because they are now on a different processor,
        // we need to update the halo layer associated with the ghost cells left in
        // the halo of the outgoing cells.
        for (long ghostId : updateList) {
            computeCellHaloLayer(ghostId);
        }

        // Prune stale interfaces
        mergeAdaptionInfo(pruneStaleInterfaces(trackPartitioning), partitioningData);

        // Identify orphan vertices
        std::vector<long> orphanVertices = findOrphanVertices();

        // Track ghost vertices deletion
        //
        // Only ghost vertices need to be tracked, all orphan internal vertex
        // have already been tracked among the vertices that have been send.
        if (trackPartitioning && !orphanVertices.empty()) {
            partitioningData.emplace_back();
            adaption::Info &vertexDeletionInfo = partitioningData.back();
            vertexDeletionInfo.entity = adaption::ENTITY_VERTEX;
            vertexDeletionInfo.type   = adaption::TYPE_DELETION;
            vertexDeletionInfo.rank   = getRank();
            for (long vertexId : orphanVertices) {
                const Vertex &vertex = getVertex(vertexId);
                if (vertex.isInterior()) {
                    continue;
                }

                vertexDeletionInfo.current.push_back(vertexId);
            }
        }

        // Delete orphan vertices
        deleteVertices(orphanVertices);
    } else {
        // All ghost cells will be deleted
        if (trackPartitioning) {
            partitioningData.emplace_back();
            adaption::Info &partitioningVertexInfo = partitioningData.back();
            partitioningVertexInfo.entity= adaption::ENTITY_VERTEX;
            partitioningVertexInfo.type  = adaption::TYPE_DELETION;
            partitioningVertexInfo.current.reserve(getGhostVertexCount());
            for (VertexConstIterator itr = ghostVertexConstBegin(); itr != ghostVertexConstEnd(); ++itr) {
                partitioningVertexInfo.current.push_back(itr.getId());
            }

            partitioningData.emplace_back();
            adaption::Info &partitioningCellInfo = partitioningData.back();
            partitioningCellInfo.entity = adaption::ENTITY_CELL;
            partitioningCellInfo.type   = adaption::TYPE_DELETION;
            partitioningCellInfo.current.reserve(getGhostCellCount());
            for (CellConstIterator itr = ghostCellConstBegin(); itr != ghostCellConstEnd(); ++itr) {
                partitioningCellInfo.current.push_back(itr.getId());
            }
        }

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

    See PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
    for the documentation about the tracking information returned by this function.

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

    // Initialize partitioning info
    std::vector<adaption::Info> partitioningData;

    // Mark border interfaces of ghost cells as dangling
    //
    // We may received cells that connect to the existing mesh through one
    // of the faces that are now borders. Marking those border interfaces as
    // dangling allows to delete them and create new internal interfaces.
    if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
        CellConstIterator endItr = ghostCellConstEnd();
        for (CellConstIterator itr = ghostCellConstBegin(); itr != endItr; ++itr) {
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
    }

    // Receive data
    int patchRank = getRank();

    int nSendCompleted = 0;
    std::vector<int> sendCompletedIndexes(nSendRanks);
    std::vector<MPI_Status> sendCompletedStatuses(nSendRanks);

    std::vector<long> neighIds;

    std::unordered_set<long> duplicateCellsCandidates;

    std::array<double, 3> *candidateVertexCoords;
    Vertex::Less vertexLess(10 * std::numeric_limits<double>::epsilon());
    auto vertexCoordsLess = [this, &vertexLess, &candidateVertexCoords](const long &id_1, const long &id_2)
    {
        const std::array<double, 3> &coords_1 = (id_1 >= 0) ? this->getVertex(id_1).getCoords() : *candidateVertexCoords;
        const std::array<double, 3> &coords_2 = (id_2 >= 0) ? this->getVertex(id_2).getCoords() : *candidateVertexCoords;

        return vertexLess(coords_1, coords_2);
    };
    std::set<long, decltype(vertexCoordsLess)> duplicateVerticesCandidates(vertexCoordsLess);

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
        std::vector<long> firstHaloLayerCells;
        for (const auto &ghostOwnerEntry : m_ghostCellInfo) {
            long ghostId = ghostOwnerEntry.first;

            duplicateCellsCandidates.insert(ghostId);
            if (getCellHaloLayer(ghostId) == 0) {
                firstHaloLayerCells.push_back(ghostId);
            }
        }

        auto duplicateCandidatesSelector = [this](long cellId) {
            BITPIT_UNUSED(cellId);

            return (m_ghostCellInfo.count(cellId) == 0);
        };

        auto duplicateCandidatesBuilder = [&duplicateCellsCandidates](long cellId, int layer) {
            BITPIT_UNUSED(layer);

            duplicateCellsCandidates.insert({cellId});

            return false;
        };

        processCellsNeighbours(firstHaloLayerCells, m_haloSize, duplicateCandidatesSelector, duplicateCandidatesBuilder);

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
        duplicateVerticesCandidates.clear();
        if (!m_partitioningSerialization) {
            for (long cellId : duplicateCellsCandidates) {
                const Cell &cell = m_cells.at(cellId);
                ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
                for (long vertexId : cellVertexIds) {
                    duplicateVerticesCandidates.insert(vertexId);
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
                        duplicateVerticesCandidates.insert(vertexId);
                    }
                }
            }
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
            // Vertex data
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
            //
            // If the vertex is a duplicate, it will be possible to obtain the
            // id of the id of the coincident vertex already in the patch.
            long vertexId;

            bool isDuplicate = (isFrame || isHalo);
            if (isDuplicate) {
                // The container that stores the list of possible duplicate
                // vertices uses a customized comparator that allows to
                // compare the coordinates of the stored vertices with the
                // coordinates of an external vertex. To check if the external
                // vertex is among the vertex in the container, the search
                // should be performed using the NULL_ID as the key.
                candidateVertexCoords = &(vertex.getCoords());
                auto candidateVertexItr = duplicateVerticesCandidates.find(Cell::NULL_ID);
                isDuplicate = (candidateVertexItr != duplicateVerticesCandidates.end());
                if (isDuplicate) {
                    vertexId = *candidateVertexItr;
                }
            }

            // Add the vertex
            //
            // If the id of the received vertex is already assigned, let the
            // patch generate a new id. Otherwise, keep the id of the received
            // vertex.
            if (!isDuplicate) {
                if (m_vertexIdGenerator) {
                    if (m_vertexIdGenerator->isAssigned(originalVertexId)) {
                        vertex.setId(Vertex::NULL_ID);
                    }
                } else if (m_vertices.exists(originalVertexId)) {
                    throw std::runtime_error("A vertex with the same id of the received vertex already exists.");
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

        std::vector<long> trackedReceivedInteriorCells;
        std::vector<long> trackedCreatedGhostCells;
        if (trackPartitioning) {
            trackedReceivedInteriorCells.reserve(nReceivedInternalCells);
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
            int cellHaloLayer;
            int ghostCellOwnershipChange;
            if (isHalo) {
                cellsBuffer >> cellOwner;
                cellsBuffer >> cellHaloLayer;
                cellsBuffer >> ghostCellOwnershipChange;
            } else {
                cellOwner = patchRank;
                cellHaloLayer = -1;
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

            bool cellAlreadyExists = (cellId >= 0);

            // If the cell is not a duplicate add it in the cell data structure,
            // otherwise merge the connectivity of the duplicate cell to the
            // existing cell. This ensure that the received cell will be
            // properly connected to the received cells
            CellIterator cellIterator;
            if (!cellAlreadyExists) {
                // Reset the interfaces of the cell
                if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
                    cell.resetInterfaces();
                }

                // Add cell
                //
                // If the id of the received cell is already assigned, let the
                // patch generate a new id. Otherwise, keep the id of the received
                // cell.

                if (m_cellIdGenerator) {
                    if (m_cellIdGenerator->isAssigned(cellOriginalId)) {
                        cell.setId(Cell::NULL_ID);
                    }
                } else if (m_cells.exists(cellOriginalId)) {
                    throw std::runtime_error("A cell with the same id of the received cell already exists.");
                }

                cellIterator = addCell(std::move(cell), cellOwner, cellHaloLayer);
                cellId = cellIterator.getId();
                addedCells.push_back(cellId);

                // Setup ghost ownership changes
                if (!cellIterator->isInterior()) {
                    if (ghostCellOwnershipChange >= 0) {
                        ghostCellOwnershipChanges->insert({cellId, ghostCellOwnershipChange});
                    }
                }
            } else {
                // Check if the existing cells needs to become an internal cell
                cellIterator = getCellIterator(cellId);
                if (isInterior && !cellIterator->isInterior()) {
                    ghostCell2InternalCell(cellId);
                    ghostCellOwnershipChanges->erase(cellId);
                }

                // Update the halo layer associated with ghost cells
                if (!isInterior) {
                    const GhostCellInfo &ghostInfo = m_ghostCellInfo.at(cellId);
                    if (ghostInfo.haloLayer > cellHaloLayer) {
                        setGhostCellInfo(cellId, ghostInfo.owner, cellHaloLayer);
                    }
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
            if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
                setCellAlterationFlags(cellId, FLAG_INTERFACES_DIRTY);
            }

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
            if (trackPartitioning) {
                if (cellIterator->isInterior()) {
                    trackedReceivedInteriorCells.emplace_back(cellId);
                } else if (!cellAlreadyExists) {
                    trackedCreatedGhostCells.emplace_back(cellId);
                }
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

        // Track changes
        if (trackPartitioning) {
            if (!trackedReceivedInteriorCells.empty()) {
                partitioningData.emplace_back();
                adaption::Info &vertexRecvInfo = partitioningData.back();
                vertexRecvInfo.entity  = adaption::ENTITY_VERTEX;
                vertexRecvInfo.type    = adaption::TYPE_PARTITION_RECV;
                vertexRecvInfo.rank    = sendRank;
                vertexRecvInfo.current = getOrderedCellsVertices(trackedReceivedInteriorCells, true, false);

                partitioningData.emplace_back();
                adaption::Info &cellRecvInfo = partitioningData.back();
                cellRecvInfo.entity = adaption::ENTITY_CELL;
                cellRecvInfo.type   = adaption::TYPE_PARTITION_RECV;
                cellRecvInfo.rank   = sendRank;
                cellRecvInfo.current = std::move(trackedReceivedInteriorCells);
                trackedReceivedInteriorCells.clear();
            }

            if (!trackedCreatedGhostCells.empty()) {
                partitioningData.emplace_back();
                adaption::Info &vertexCreationInfo = partitioningData.back();
                vertexCreationInfo.entity  = adaption::ENTITY_VERTEX;
                vertexCreationInfo.type    = adaption::TYPE_CREATION;
                vertexCreationInfo.rank    = patchRank;
                vertexCreationInfo.current = getOrderedCellsVertices(trackedCreatedGhostCells, false, true);


                partitioningData.emplace_back();
                adaption::Info &cellCreationInfo = partitioningData.back();
                cellCreationInfo.entity  = adaption::ENTITY_CELL;
                cellCreationInfo.type    = adaption::TYPE_CREATION;
                cellCreationInfo.rank    = patchRank;
                cellCreationInfo.current = std::move(trackedCreatedGhostCells);
                trackedCreatedGhostCells.clear();
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
        mergeAdaptionInfo(updateInterfaces(false, trackPartitioning), partitioningData);
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
	auto cellInfoItr = m_ghostCellInfo.find(id);
	if (cellInfoItr == m_ghostCellInfo.end()) {
		return m_rank;
	} else {
		return cellInfoItr->second.owner;
	}
}

/*!
	Gets the rank of the process that owns the specified cell.

	\param id is the id of the requested cell
	\result The rank that owns the specified cell.
*/
int PatchKernel::getCellOwner(long id) const
{
	auto cellInfoItr = m_ghostCellInfo.find(id);
	if (cellInfoItr == m_ghostCellInfo.end()) {
		return m_rank;
	} else {
		return cellInfoItr->second.owner;
	}
}

/*!
	Gets the halo layer of the specified cell.

	\param id is the id of the requested cell
	\result The halo layer of the specified cell.
*/
int PatchKernel::getCellHaloLayer(long id) const
{
	auto cellInfoItr = m_ghostCellInfo.find(id);
	if (cellInfoItr == m_ghostCellInfo.end()) {
		return -1;
	} else {
		return cellInfoItr->second.haloLayer;
	}
}

/*!
	Gets the rank of the process that owns the specified vertex.

	\param id is the id of the requested vertex
	\result The rank that owns the specified vertex.
*/
int PatchKernel::getVertexRank(long id) const
{
	auto vertexInfoItr = m_ghostVertexInfo.find(id);
	if (vertexInfoItr == m_ghostVertexInfo.end()) {
		return m_rank;
	} else {
		return vertexInfoItr->second.owner;
	}
}

/*!
	Gets the rank of the process that owns the specified vertex.

	\param id is the id of the requested vertex
	\result The rank that owns the specified vertex.
*/
int PatchKernel::getVertexOwner(long id) const
{
	auto vertexInfoItr = m_ghostVertexInfo.find(id);
	if (vertexInfoItr == m_ghostVertexInfo.end()) {
		return m_rank;
	} else {
		return vertexInfoItr->second.owner;
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
	\param owner is the rank of the process that owns the ghost vertex
*/
void PatchKernel::setGhostVertexInfo(long id, int owner)
{
	auto ghostVertexInfoItr = m_ghostVertexInfo.find(id);
	if (ghostVertexInfoItr != m_ghostVertexInfo.end()) {
		ghostVertexInfoItr->second.owner = owner;
	} else {
		GhostVertexInfo ghostVertexInfo;
		ghostVertexInfo.owner = owner;
		m_ghostVertexInfo.insert({id, std::move(ghostVertexInfo)});
	}

	setPartitioningInfoDirty(true);
}

/*!
	Unsets the owner of the specified ghost vertex.

	\param id is the id of the ghost vertex
*/
void PatchKernel::unsetGhostVertexInfo(long id)
{
	auto ghostVertexInfoItr = m_ghostVertexInfo.find(id);
	if (ghostVertexInfoItr == m_ghostVertexInfo.end()) {
		return;
	}

	m_ghostVertexInfo.erase(ghostVertexInfoItr);

	setPartitioningInfoDirty(true);
}

/*!
	Clear the owners of all the ghost vertices.

	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::clearGhostVerticesInfo()
{
	m_ghostVertexInfo.clear();

	setPartitioningInfoDirty(true);
}

/*!
	Sets the information of the specified ghost.

	\param id is the id of the ghost cell
	\param owner is the rank of the process that owns the ghost cell
	\param haloLayer is the halo layer the ghost cell belongs to
*/
void PatchKernel::setGhostCellInfo(long id, int owner, int haloLayer)
{
	auto ghostCellInfoItr = m_ghostCellInfo.find(id);
	if (ghostCellInfoItr != m_ghostCellInfo.end()) {
		ghostCellInfoItr->second.owner     = owner;
		ghostCellInfoItr->second.haloLayer = haloLayer;
	} else {
		GhostCellInfo ghostCellInfo;
		ghostCellInfo.owner     = owner;
		ghostCellInfo.haloLayer = haloLayer;
		m_ghostCellInfo.insert({id, std::move(ghostCellInfo)});
	}

	setPartitioningInfoDirty(true);
}

/*!
	Unsets the information of the specified ghost.

	\param id is the id of the ghost cell
*/
void PatchKernel::unsetGhostCellInfo(long id)
{
	auto ghostCellInfoItr = m_ghostCellInfo.find(id);
	if (ghostCellInfoItr == m_ghostCellInfo.end()) {
		return;
	}

	m_ghostCellInfo.erase(ghostCellInfoItr);

	setPartitioningInfoDirty(true);
}

/*!
	Compute the halo layer associated with the specified ghost cell.

	The neighbours are processed one layer at the time until an internal cell
	is found. If an internal cell cannot be found in the specified halo size,
	an error is thrown.

	\param id is the id of the cell
*/
void PatchKernel::computeCellHaloLayer(int id)
{
	// Early return if the cell is interior
	if (getCell(id).isInterior()) {
		unsetGhostCellInfo(id);
		return;
	}

	// Process the neighbours until we find an internal cell
	bool haloLayerIdentified = false;

	std::array<long, 1> layerUpdateSeed = {{ id }};

	auto layerUpdateSelector = [](long cellId) {
		BITPIT_UNUSED(cellId);

		return true;
	};

	auto layerUpdateProcessor = [this, id, &haloLayerIdentified](long cellId, int layer) {
		const Cell &cell = getCell(cellId);
		if (cell.isInterior()) {
			const GhostCellInfo &cellGhostInfo = m_ghostCellInfo.at(id);
			if (cellGhostInfo.haloLayer != layer) {
				setGhostCellInfo(id, cellGhostInfo.owner, layer);
			}

			haloLayerIdentified = true;

			return true;
		}

		return false;
	};

	processCellsNeighbours(layerUpdateSeed, m_haloSize, layerUpdateSelector, layerUpdateProcessor);

	if (!haloLayerIdentified) {
		throw std::runtime_error ("Unable to identify the halo layer of the cell.");
	}
}

/*!
	Clear the information of all the ghosts.
*/
void PatchKernel::clearGhostCellsInfo()
{
	m_ghostCellInfo.clear();

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
	VertexConstIterator previousBeginItr = ghostVertexConstBegin();
	VertexConstIterator previousEndItr   = ghostVertexConstEnd();

	std::vector<long> previousGhosts;
	previousGhosts.reserve(getGhostVertexCount());
	for (VertexConstIterator vertexItr = previousBeginItr; vertexItr != previousEndItr; ++vertexItr) {
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
	VertexConstIterator beginItr = ghostVertexConstEnd();
	VertexConstIterator endItr   = ghostVertexConstEnd();

	for (VertexConstIterator vertexItr = beginItr; vertexItr != endItr; ++vertexItr) {
		long vertexId = vertexItr->getId();
		int vertexOwner = exchangeVertexOwners.at(vertexId);
		setGhostVertexInfo(vertexId, vertexOwner);
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
	for (const auto &entry : m_ghostVertexInfo) {
		int ghostVertexOwner = entry.second.owner;
		long ghostVertexId = entry.first;
		m_ghostVertexExchangeTargets[ghostVertexOwner].push_back(ghostVertexId);
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
				if (m_ghostVertexInfo.count(vertexId) > 0) {
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
				if (m_ghostVertexInfo.count(vertexId) > 0) {
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
	for (const auto &entry : m_ghostCellInfo) {
		int ghostOwner = entry.second.owner;
		long ghostCellId = entry.first;
		m_ghostCellExchangeTargets[ghostOwner].push_back(ghostCellId);
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
	// If there are no targets, there will be no sources either.
	auto ghostExchangeTargetsItr = m_ghostCellExchangeTargets.find(rank);
	if (ghostExchangeTargetsItr == m_ghostCellExchangeTargets.end()) {
		return std::vector<long>();
	}

	const std::vector<long> &rankTargets = ghostExchangeTargetsItr->second;

	// Generate sources from the targets
	std::vector<long> exchangeSources;

	auto sourceSelector = [](long cellId) {
		BITPIT_UNUSED(cellId);

		return true;
	};

	auto sourceBuilder = [this, &exchangeSources](long cellId, int layer) {
		BITPIT_UNUSED(layer);

		if (m_ghostCellInfo.count(cellId) == 0) {
			exchangeSources.push_back(cellId);
		}

		return false;
	};

	processCellsNeighbours(rankTargets, m_haloSize, sourceSelector, sourceBuilder);

	return exchangeSources;
}

}

#endif
