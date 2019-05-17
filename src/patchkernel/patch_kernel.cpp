/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_SA.hpp"

#include "patch_info.hpp"
#include "patch_kernel.hpp"
#include "patch_manager.hpp"

namespace bitpit {

/*!
	\class PatchKernel
	\ingroup patchkernel

	\brief The PatchKernel class provides an interface for defining patches.

	PatchKernel is the base class for defining patches.
*/

/*!
	Creates a new patch.

	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(bool expert)
	: m_expert(expert)
{
	// Initialize the patch
	initialize();

	// Register the patch
	patch::manager().registerPatch(this);
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(const int &dimension, bool expert)
	: m_expert(expert)
{
	// Initialize the patch
	initialize();

	// Register the patch
	patch::manager().registerPatch(this);

	// Set the dimension
	setDimension(dimension);
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(const int &id, const int &dimension, bool expert)
	: m_expert(expert)
{
	// Initialize the patch
	initialize();

	// Register the patch
	patch::manager().registerPatch(this, id);

	// Set the dimension
	//
	// Here we can only call the base function, if needed, every derived patch
	// has to call its derived function.
	PatchKernel::setDimension(dimension);
}

/*!
	Copy constructor

	\param other is another patch whose content is copied into this
*/
PatchKernel::PatchKernel(const PatchKernel &other)
    : VTKBaseStreamer(other),
      m_vertices(other.m_vertices),
      m_cells(other.m_cells),
      m_interfaces(other.m_interfaces),
      m_vertexIdGenerator(other.m_vertexIdGenerator),
      m_interfaceIdGenerator(other.m_interfaceIdGenerator),
      m_cellIdGenerator(other.m_cellIdGenerator),
      m_nInternals(other.m_nInternals),
#if BITPIT_ENABLE_MPI==1
      m_nGhosts(other.m_nGhosts),
#endif
      m_lastInternalId(other.m_lastInternalId),
#if BITPIT_ENABLE_MPI==1
      m_firstGhostId(other.m_firstGhostId),
#endif
      m_vtk(other.m_vtk),
      m_vtkWriteTarget(other.m_vtkWriteTarget),
      m_vtkVertexMap(other.m_vtkVertexMap),
      m_boxFrozen(other.m_boxFrozen),
      m_boxDirty(other.m_boxDirty),
      m_boxMinPoint(other.m_boxMinPoint),
      m_boxMaxPoint(other.m_boxMaxPoint),
      m_boxMinCounter(other.m_boxMinCounter),
      m_boxMaxCounter(other.m_boxMaxCounter),
      m_adjacenciesBuildStrategy(other.m_adjacenciesBuildStrategy),
      m_interfacesBuildStrategy(other.m_interfacesBuildStrategy),
      m_spawnStatus(other.m_spawnStatus),
      m_adaptionStatus(other.m_adaptionStatus),
      m_expert(other.m_expert),
      m_dimension(other.m_dimension),
      m_hasCustomTolerance(other.m_hasCustomTolerance),
      m_tolerance(other.m_tolerance),
      m_rank(other.m_rank),
      m_nProcessors(other.m_nProcessors)
#if BITPIT_ENABLE_MPI==1
      , m_communicator(MPI_COMM_NULL),
      m_partitioned(other.m_partitioned),
      m_partitioningStatus(other.m_partitioningStatus),
      m_haloSize(other.m_haloSize),
      m_ghostOwners(other.m_ghostOwners),
      m_ghostExchangeTargets(other.m_ghostExchangeTargets),
      m_ghostExchangeSources(other.m_ghostExchangeSources)
#endif
{
	// Register the patch
	patch::manager().registerPatch(this);

	// Update the VTK streamer
	//
	// The pointer to VTK streamers are copied, if there are pointer to the
	// original object they have to be replace with a pointer to this object.
	std::vector<std::string> streamedGeomFields;
	streamedGeomFields.reserve(m_vtk.getGeomDataCount());
	for (auto itr = m_vtk.getGeomDataBegin(); itr != m_vtk.getGeomDataEnd(); ++itr) {
		const VTKField &field = *itr;
		if (&field.getStreamer() != &other) {
			continue;
		}

		streamedGeomFields.push_back(field.getName());
	}

	for (const std::string &name : streamedGeomFields) {
		const VTKField &field = *(m_vtk.findGeomData(name));
		VTKField updatedField(field);
		updatedField.setStreamer(*this);

		m_vtk.setGeomData(std::move(updatedField));
	}

	std::vector<std::string> streamedDataFields;
	streamedDataFields.reserve(m_vtk.getDataCount());
	for (auto itr = m_vtk.getDataBegin(); itr != m_vtk.getDataEnd(); ++itr) {
		const VTKField &field = *itr;
		if (&field.getStreamer() != &other) {
			continue;
		}

		streamedDataFields.push_back(field.getName());
	}

	for (const std::string &name : streamedDataFields) {
		const VTKField &field = *(m_vtk.findData(name));
		VTKField updatedField(field);
		updatedField.setStreamer(*this);

		m_vtk.removeData(field.getName());
		m_vtk.addData(std::move(updatedField));
	}

#if BITPIT_ENABLE_MPI==1
	// Set the communicator
	MPI_Comm communicator = other.getCommunicator();
	if (communicator != MPI_COMM_NULL) {
		setCommunicator(communicator);
	}
#endif
}

/*!
	Initialize the patch
*/
void PatchKernel::initialize()
{
	// Id
	m_id = PatchManager::AUTOMATIC_ID;

	// Cell count
	m_nInternals = 0;
#if BITPIT_ENABLE_MPI==1
	m_nGhosts    = 0;
#endif

	m_lastInternalId = Cell::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostId   = Cell::NULL_ID;
#endif

	// Dimension
	m_dimension = -1;

	// Set adjacencies build strategy
	setAdjacenciesBuildStrategy(ADJACENCIES_NONE);

	// Set interfaces build strategy
	setInterfacesBuildStrategy(INTERFACES_NONE);

	// Set the spawn as unneeded
	//
	// Specific implementation will set the appropriate status during their
	// initialization.
	setSpawnStatus(SPAWN_UNNEEDED);

	// Set the adaption as unsupported
	//
	// Specific implementation will set the appropriate status during their
	// initialization.
	setAdaptionStatus(ADAPTION_UNSUPPORTED);

	// Parallel information
	m_rank        = 0;
	m_nProcessors = 1;
#if BITPIT_ENABLE_MPI==1

	// Patch is not partitioned
	m_communicator = MPI_COMM_NULL;
	m_haloSize = 0;
	setPartitioned(false);

	// Set the partitioning as unsupported
	//
	// Specific implementation will set the appropriate status during their
	// initialization.
	setPartitioningStatus(PARTITIONING_UNSUPPORTED);
#endif

	// Initialize the geometrical tolerance to a default value
	m_hasCustomTolerance = false;
	_setTol(DEFAULT_TOLERANCE);

	// Initializes the bounding box
	setBoundingBoxFrozen(false);
	clearBoundingBox();

	// Set VTK write target
	m_vtkWriteTarget = WRITE_TARGET_CELLS_ALL;

	// Set VTK information
	std::ostringstream convert;
	convert << getId();

	m_vtk.setName(convert.str());
	m_vtk.setCodex(VTKFormat::APPENDED);

	// Set VTK Geom Data
	m_vtk.setGeomData<double>(VTKUnstructuredField::POINTS, this);
	m_vtk.setGeomData<int>(VTKUnstructuredField::OFFSETS, this);
	m_vtk.setGeomData<int>(VTKUnstructuredField::TYPES, this);
	m_vtk.setGeomData<long>(VTKUnstructuredField::CONNECTIVITY, this);
	m_vtk.setGeomData<long>(VTKUnstructuredField::FACE_STREAMS, this);
	m_vtk.setGeomData<int>(VTKUnstructuredField::FACE_OFFSETS, this);

	// Add VTK basic patch data
	m_vtk.addData<long>("cellIndex", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("PID", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<long>("vertexIndex", VTKFieldType::SCALAR, VTKLocation::POINT, this);
#if BITPIT_ENABLE_MPI==1
	m_vtk.addData<long>("cellGlobalIndex", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("rank", VTKFieldType::SCALAR, VTKLocation::CELL, this);
#endif
}

/*!
	Destroys the patch.
*/
PatchKernel::~PatchKernel()
{
#if BITPIT_ENABLE_MPI==1
	freeCommunicator();
#endif

	try {
		patch::manager().unregisterPatch(this);
	} catch (const std::runtime_error &e) {
		log::cout() << "Unable to unregister the patch" << std::endl;
		log::cout() << " Error message: " << e.what() << endl;
	}
}

/*!
	Commit all pending changes.

	\param trackAdaption if set to true the changes to the patch will be
	tracked
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the update
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
std::vector<adaption::Info> PatchKernel::update(bool trackAdaption, bool squeezeStorage)
{
	std::vector<adaption::Info> updateInfo;

	// Check if there are pending changes
	bool spawnNeeed       = (getSpawnStatus() == SPAWN_NEEDED);
	bool adaptionDirty    = (getAdaptionStatus(true) == ADAPTION_DIRTY);
	bool boundingBoxDirty = isBoundingBoxDirty();

	bool pendingChanges = (spawnNeeed || adaptionDirty || boundingBoxDirty);
	if (!pendingChanges) {
		return updateInfo;
	}

	// Spawn
	if (spawnNeeed) {
		mergeAdaptionInfo(spawn(trackAdaption), updateInfo);
	}

	// Adaption
	if (adaptionDirty) {
		mergeAdaptionInfo(adaption(trackAdaption, squeezeStorage), updateInfo);
	}

	// Update bounding box
	//
	// Previous updates may already have updated the bounding box, before
	// doing the update re-check if it is still dirty.
	boundingBoxDirty = isBoundingBoxDirty();
	if (boundingBoxDirty) {
		updateBoundingBox();
	}

	return updateInfo;
}

/*!
	Generates the patch.

	\param trackSpawn if set to true the changes to the patch will be tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the spawn.
*/
std::vector<adaption::Info> PatchKernel::spawn(bool trackSpawn)
{
	std::vector<adaption::Info> spawnInfo;

#if BITPIT_ENABLE_MPI==1
	// This is a collevtive operation and should be called by all processes
	if (isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Barrier(communicator);
	}
#endif

	// Check spawn status
	SpawnStatus spawnStatus = getSpawnStatus();
	if (spawnStatus == SPAWN_UNNEEDED || spawnStatus == SPAWN_DONE) {
		return spawnInfo;
	}

	// Begin patch alteration
	beginAlteration();

	// Alter patch
	spawnInfo = _spawn(trackSpawn);

	// End patch alteration
	endAlteration(true);

	// Spwan is done
	setSpawnStatus(SPAWN_DONE);

	// Done
	return spawnInfo;
}

/*!
	Execute patch adaption.

	\param trackAdaption if set to true the changes to the patch will be
	tracked
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the update
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
std::vector<adaption::Info> PatchKernel::adaption(bool trackAdaption, bool squeezeStorage)
{
	std::vector<adaption::Info> adaptionInfo;

	// Check adaption status
	AdaptionStatus adaptionStatus = getAdaptionStatus(true);
	if (adaptionStatus == ADAPTION_UNSUPPORTED || adaptionStatus == ADAPTION_CLEAN) {
		return adaptionInfo;
	} else if (adaptionStatus != ADAPTION_DIRTY) {
		throw std::runtime_error ("An adaption is already in progress.");
	}

	adaptionPrepare(false);

	adaptionInfo = adaptionAlter(trackAdaption, squeezeStorage);

	adaptionCleanup();

	return adaptionInfo;
}

/*!
	Prepares the patch for performing the adaption.

	The patch will prepare its data structured for the adaption and optionally
	track the changes that will be performed to the patch. During this phase
	no changes will be performed to the patch.

	\param trackAdaption if set to true the function will return the changes
	that will be performed in the alter step
	\result If the adaption is tracked, returns a vector of adaption::Info that
	can be used to discover what changes will be performed in the alter step,
	otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::adaptionPrepare(bool trackAdaption)
{
	std::vector<adaption::Info> adaptionInfo;

	// Check adaption status
	AdaptionStatus adaptionStatus = getAdaptionStatus(true);
	if (adaptionStatus == ADAPTION_UNSUPPORTED || adaptionStatus == ADAPTION_CLEAN) {
		return adaptionInfo;
	} else if (adaptionStatus != ADAPTION_DIRTY) {
		throw std::runtime_error ("An adaption is already in progress.");
	}

	// Execute the adaption preparation
	adaptionInfo = _adaptionPrepare(trackAdaption);

	// Update the status
	setAdaptionStatus(ADAPTION_PREPARED);

	return adaptionInfo;
}

/*!
	Alter the patch performing the adpation.

	The actual modification of the patch takes place during this phase. After
	this phase the adapton is completed and the patch is in its final state.
	Optionally the patch can track the changes performed to the patch.

	\param trackAdaption if set to true the function will return the changes
	done to the patch during the adaption
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the adaption
	\result If the adaption is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the adaption, otherwise an
	empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::adaptionAlter(bool trackAdaption, bool squeezeStorage)
{
	std::vector<adaption::Info> adaptionInfo;

	// Check adaption status
	AdaptionStatus adaptionStatus = getAdaptionStatus();
	if (adaptionStatus == ADAPTION_UNSUPPORTED || adaptionStatus == ADAPTION_CLEAN) {
		return adaptionInfo;
	} else if (adaptionStatus != ADAPTION_PREPARED) {
		throw std::runtime_error ("The prepare function has not been called.");
	}

	// Begin patch alteration
	beginAlteration();

	// Alter patch
	adaptionInfo = _adaptionAlter(trackAdaption);

	// End patch alteration
	endAlteration(squeezeStorage);

	// Update the status
	setAdaptionStatus(ADAPTION_ALTERED);

	return adaptionInfo;
}

/*!
	Cleanup patch data structured after the adaption.

	The patch will only clean-up the data structures needed for adaption.
*/
void PatchKernel::adaptionCleanup()
{
	AdaptionStatus adaptionStatus = getAdaptionStatus();
	if (adaptionStatus == ADAPTION_UNSUPPORTED || adaptionStatus == ADAPTION_CLEAN) {
		return;
	} else if (adaptionStatus == ADAPTION_PREPARED) {
		throw std::runtime_error ("It is not yet possible to abort an adaption.");
	} else if (adaptionStatus != ADAPTION_ALTERED) {
		throw std::runtime_error ("The alter function has not been called.");
	}

	// Complete the adaption
	_adaptionCleanup();

	// Update the status
	setAdaptionStatus(ADAPTION_CLEAN);
}

/*!
	Make the adaption markers set by the user consistent with the internal
	criteria defined by the patch.
*/
void PatchKernel::settleAdaptionMarkers()
{
	// Nothing to do
}

/*!
	Begin patch alteration.
*/
void PatchKernel::beginAlteration()
{
}

/*!
	End patch alteration.

	\param squeezeStorage if set to true patch data structures will be
	squeezed after the adaption
*/
void PatchKernel::endAlteration(bool squeezeStorage)
{
	// Flush data structures
	m_cells.flush();
	m_interfaces.flush();
	m_vertices.flush();

	// Squeeze data structures
	if (squeezeStorage) {
		squeeze();
	}

	// Update geometric information
	updateBoundingBox();

#if BITPIT_ENABLE_MPI==1
	// Update information for ghost data exchange
	//
	// If we are partitioning the patch, the partition flag is not set yet (it
	// will be set after the call to this function). We call the update of
	// information for ghost data exchange unconditionally, if there are no
	// ghosts the function will exit without doing anything.
	if (getProcessorCount() > 1 && getAdjacenciesBuildStrategy() != ADJACENCIES_NONE) {
		updateGhostExchangeInfo();
	}
#endif

	// Synchronize storage
	m_cells.sync();
	m_interfaces.sync();
	m_vertices.sync();
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void PatchKernel::markCellForRefinement(const long &id)
{
	bool updated = _markCellForRefinement(id);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void PatchKernel::markCellForCoarsening(const long &id)
{
	bool updated = _markCellForCoarsening(id);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
}

/*!
	Resets the adaption marker of the specified cell.

	\param id is the id of the cell
*/
void PatchKernel::resetCellAdaptionMarker(const long &id)
{
	bool updated = _resetCellAdaptionMarker(id);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
}

/*!
	Returns the adaption marker of the specified cell.

	The marker only defines the type of adaption requested for the cell, it
	is not guaranteed that the adaption will effectively perfrom the requestd
	action (i.e., the requested marker may not be consistent with the internal
	criteria defined by the patch).

	\param id is the id of the cell
	\return The adaption marker of the cell.
*/
adaption::Marker PatchKernel::getCellAdaptionMarker(const long &id)
{
	return _getCellAdaptionMarker(id);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void PatchKernel::enableCellBalancing(const long &id, bool enabled)
{
	bool updated = _enableCellBalancing(id, enabled);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
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
#if BITPIT_ENABLE_MPI==1
	m_nGhosts = 0;
#endif
	m_lastInternalId = Cell::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostId = Cell::NULL_ID;
#endif

#if BITPIT_ENABLE_MPI==1
	clearGhostOwners();
#endif

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
	// Early return if no adjacencies have been built
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Clear interfaces
	m_interfaces.clear();
	PiercedVector<Interface>().swap(m_interfaces);
	m_interfaceIdGenerator.reset();

	for (auto &cell : m_cells) {
		cell.resetInterfaces();
	}

	// Set interface build strategy
	setInterfacesBuildStrategy(INTERFACES_NONE);
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
    \param mode is the VTK file mode that will be used for writing the patch
*/
void PatchKernel::write(const std::string &filename, VTKWriteMode mode)
{
	std::string oldFilename = m_vtk.getName();

	m_vtk.setName(filename);
	write(mode);
	m_vtk.setName(oldFilename);
}

/*!
	Writes the patch a filename with the same name of the patch

    \param mode is the VTK file mode that will be used for writing the patch
*/
void PatchKernel::write(VTKWriteMode mode)
{
	// Get VTK cell count
	long vtkCellCount = 0;
	if (m_vtkWriteTarget == WRITE_TARGET_CELLS_ALL) {
		vtkCellCount = getCellCount();
#if BITPIT_ENABLE_MPI==1
	} else if (m_vtkWriteTarget == WRITE_TARGET_CELLS_INTERNAL) {
		vtkCellCount = getInternalCount();
#endif
	}

	// Set the dimensinos of the mesh
	PiercedStorage<long, long> vertexWriteFlag(1, &m_vertices);
	vertexWriteFlag.fill(false);

	bool vtkFaceStreamNeeded = false;
	for (const Cell &cell : m_cells) {
		if (cell.getDimension() > 2 && !cell.hasInfo()) {
			vtkFaceStreamNeeded = true;
			break;
		}
	}

	long vtkConnectSize    = 0;
	long vtkFaceStreamSize = 0;
	for (const Cell &cell : getVTKCellWriteRange()) {
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		const int nCellVertices = cellVertexIds.size();
		for (int k = 0; k < nCellVertices; ++k) {
			long vertexId = cellVertexIds[k];
			vertexWriteFlag.at(vertexId) = true;
		}

		vtkConnectSize += nCellVertices;
		if (vtkFaceStreamNeeded) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				vtkFaceStreamSize += 1;
			} else {
				vtkFaceStreamSize += cell.getFaceStreamSize();
			}
		}
	}

	int vtkVertexCount = 0;
	m_vtkVertexMap.unsetKernel(true);
	m_vtkVertexMap.setStaticKernel(&m_vertices);
	for (VertexConstIterator itr = vertexConstBegin(); itr != vertexConstEnd(); ++itr) {
		std::size_t vertexRawId = itr.getRawIndex();
		if (vertexWriteFlag.rawAt(vertexRawId)) {
			m_vtkVertexMap.rawAt(vertexRawId) = vtkVertexCount++;
		} else {
			m_vtkVertexMap.rawAt(vertexRawId) = Vertex::NULL_ID;
		}
	}

	m_vtk.setDimensions(vtkCellCount, vtkVertexCount, vtkConnectSize, vtkFaceStreamSize);

	// Write the mesh
	m_vtk.write(mode);
}

/*!
	Returns the current spawn status.

	\return The current spawn status.
*/
PatchKernel::SpawnStatus PatchKernel::getSpawnStatus() const
{
	// There is no need to check the spawn status globally because the spawn
	// status will always be the same on all the processors

	return m_spawnStatus;
}

/*!
	Set the current spawn status.

	\param status is the spawn status that will be set
*/
void PatchKernel::setSpawnStatus(SpawnStatus status)
{
	m_spawnStatus = status;
}

/*!
	Returns the current adaption status.

	\param global if set to true the adaption status will be
	\return The current adaption status.
*/
PatchKernel::AdaptionStatus PatchKernel::getAdaptionStatus(bool global) const
{
	int adaptionStatus = static_cast<int>(m_adaptionStatus);

#if BITPIT_ENABLE_MPI==1
	if (global && isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &adaptionStatus, 1, MPI_INT, MPI_MAX, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif

	return static_cast<AdaptionStatus>(adaptionStatus);
}

/*!
	Set the current adaption status.

	\param status is the adaption status that will be set
*/
void PatchKernel::setAdaptionStatus(AdaptionStatus status)
{
	m_adaptionStatus = status;
}

/*!
	Returns true if the the patch needs to update its data strucutres.

	\return This method returns true to indicate the patch needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool PatchKernel::isDirty(bool global) const
{
	bool spawnNeeed       = (getSpawnStatus() == SPAWN_NEEDED);
	bool adaptionDirty    = (getAdaptionStatus(global) == ADAPTION_DIRTY);
	bool boundingBoxDirty = isBoundingBoxDirty(global);

	return (spawnNeeed || adaptionDirty || boundingBoxDirty);
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
	if (dimension == m_dimension) {
		return;
	}

	// If the dimension was already assigned, reset the patch
	if (m_dimension > 0) {
		reset();
	}

	// Set the dimension
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
	Gets a constant reference to the vertices owned by the patch.

	\return A constant reference to the vertices owned by the patch.
*/
const PiercedVector<Vertex> & PatchKernel::getVertices() const
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
	return m_vertices.find(id);
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
	Returns a constant iterator pointing to the specified vertex.

	\result A constant iterator to the specified vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::getVertexConstIterator(const long &id) const
{
	return m_vertices.find(id);
}

/*!
	Returns a constant iterator pointing to the first vertex.

	\result A constant iterator to the first vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::vertexConstBegin() const
{
	return m_vertices.cbegin();
}

/*!
	Returns a constant iterator pointing to last vertex.

	\result A constant iterator to the last vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::vertexConstEnd() const
{
	return m_vertices.cend();
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

	return m_vertexIdGenerator.generate();
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, a new unique id will be generated
	for the vertex
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(const Vertex &source, long id)
{
	if (id < 0) {
		id = generateVertexId();
	}

	Vertex vertex = source;

	return addVertex(std::move(vertex), id);
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(Vertex &&source, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	VertexIterator iterator = addVertex(source.getCoords(), id);

	Vertex &vertex = (*iterator);
	id = vertex.getId();
	vertex = std::move(source);
	vertex.setId(id);

	return iterator;
}

/*!
	Adds a new vertex with the specified coordinates.

	\param coords are the coordinates of the vertex
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, a new unique id will be generated
	for the vertex
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(const std::array<double, 3> &coords, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	if (id < 0) {
		id = generateVertexId();
	}

	// Add the vertex
	PiercedVector<Vertex>::iterator iterator = m_vertices.emreclaim(id, id, coords);

	// Update the bounding box
	addPointToBoundingBox(iterator->getCoords());

	return iterator;
}

/*!
	Resore the vertex with the specified id.

	The kernel should already contain the vertex, only the contents of the
	vertex will be updated.

	\param coords are the coordinates of the vertex
	\param id is the id of the vertex to restore
	\return An iterator pointing to the restored vertex.
*/
PatchKernel::VertexIterator PatchKernel::restoreVertex(const std::array<double, 3> &&coords, const long &id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = m_vertices.find(id);
	if (iterator == m_vertices.end()) {
		throw std::runtime_error("Unable to restore the specified vertex: the kernel doesn't contain an entry for that vertex.");
	}

	Vertex &vertex = *iterator;
	vertex.initialize(id, std::move(coords));

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
	m_vertexIdGenerator.trash(id);

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

			int nFaceVertices = cell.getFaceVertexCount(i);
			for (int k = 0; k < nFaceVertices; ++k) {
				long faceVertexId = cell.getFaceVertexId(i, k);
				freeVertices.insert(faceVertexId);
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
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int i = 0; i < nCellVertices; ++i) {
			usedVertices.insert(cellVertexIds[i]);
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
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int i = 0; i < nCellVertices; ++i) {
			vertexSet.erase(cellVertexIds[i]);
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
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int j = 0; j < nCellVertices; ++j) {
			binEntry[0] = cell.getId();
			binEntry[1] = j;

			long vertexId = cellVertexIds[j];
			long binId = bin_index[vertexId];
			bins[binId].push_back(binEntry);
		}
	}

	// Free memory
	std::unordered_map<long, long>().swap(bin_index);

	// ====================================================================== //
	// COLLAPSE DOUBLE VERTICES                                               //
	// ====================================================================== //
	std::size_t nCollapsedVertices = 0;
	std::unordered_map<long, long> vertexMap;
	vertexMap.reserve(getVertexCount());

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
				ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

				long k        = bin[list[j]][1];
				long vertexId = cellVertexIds[k];
				if (vertexMap.count(vertexId) > 0) {
					continue;
				}

				Vertex &vertex = m_vertices[vertexId];

				long collapsedVertexId;
				if (kd.exist(&vertex, collapsedVertexId) < 0) {
					collapsedVertexId = vertexId;
					kd.insert(&vertex, vertexId);
				} else {
					++nCollapsedVertices;
				}

				vertexMap.insert({vertexId, collapsedVertexId});
			}
		}
	}

	for (Cell &cell : m_cells) {
		cell.renumberVertices(vertexMap);
	}

	// Create the list of collapsed vertices
	collapsedVertices.resize(nCollapsedVertices);

	size_t k = 0;
	for (auto vertexItr = m_vertices.begin(); vertexItr != m_vertices.end(); ++vertexItr) {
		// Skip orphan vertices or unique vertices
		long vertexId = vertexItr.getId();
		auto vertexMapItr = vertexMap.find(vertexId);
		if (vertexMapItr == vertexMap.end() || vertexMap.at(vertexId) == vertexId) {
			continue;
		}

		// Add vertex id to the list of collpased vertices
		collapsedVertices[k++] = vertexId;
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
	Return true if the patch is emtpy.

	\return Return true if the patch is emtpy.
*/
bool PatchKernel::empty() const
{
	bool isEmpty = (getCellCount() == 0);
#if BITPIT_ENABLE_MPI==1
	if (isPartitioned()) {
		MPI_Allreduce(MPI_IN_PLACE, &isEmpty, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
	}
#endif

	return isEmpty;
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
	Gets the cells owned by the patch.

	\return The cells owned by the patch.
*/
PiercedVector<Cell> & PatchKernel::getCells()
{
	return m_cells;
}

/*!
	Gets a constant reference of the cells owned by the patch.

	\return A constant reference of the cells owned by the patch.
*/
const PiercedVector<Cell> & PatchKernel::getCells() const
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
ElementType PatchKernel::getCellType(const long &id) const
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
	Returns an iterator pointing to the specified cell.

	\result An iterator to the specified cell.
*/
PatchKernel::CellIterator PatchKernel::getCellIterator(const long &id)
{
	return m_cells.find(id);
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
	if (m_nInternals > 0) {
		return ++m_cells.find(m_lastInternalId);
	} else {
		return m_cells.end();
	}
}

/*!
	Returns a constant iterator pointing to the specified cell.

	\result A constant iterator to the specified cell.
*/
PatchKernel::CellConstIterator PatchKernel::getCellConstIterator(const long &id) const
{
	return m_cells.find(id);
}

/*!
	Returns a constant iterator pointing to the first cell.

	\result A constant iterator to the first cell.
*/
PatchKernel::CellConstIterator PatchKernel::cellConstBegin() const
{
	return m_cells.cbegin();
}

/*!
	Returns a constant iterator pointing to last cell.

	\result A constant iterator to the last cell.
*/
PatchKernel::CellConstIterator PatchKernel::cellConstEnd() const
{
	return m_cells.cend();
}

/*!
	Returns a constant iterator pointing to the first internal cell.

	\result A constant iterator to the first internal cell.
*/
PatchKernel::CellConstIterator PatchKernel::internalConstBegin() const
{
	return m_cells.cbegin();
}

/*!
	Returns a constant iterator pointing to the end of the list of internal cells.

	\result A constant iterator to the end of the list of internal cells.
*/
PatchKernel::CellConstIterator PatchKernel::internalConstEnd() const
{
	if (m_nInternals > 0) {
		return ++m_cells.find(m_lastInternalId);
	} else {
		return m_cells.cend();
	}
}

/*!
	Generates a new unique id for the cells.

	\result A new unique id for the cells.
*/
long PatchKernel::generateCellId()
{
	if (!isExpert()) {
		return Cell::NULL_ID;
	}

	return m_cellIdGenerator.generate();
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, long id)
{
	if (id < 0) {
		id = generateCellId();
	}

	Cell cell = source;

	return addCell(std::move(cell), id);
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	int connectSize = source.getConnectSize();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	if (!source.hasInfo()){
		std::copy(source.getConnect(), source.getConnect() + connectSize, connectStorage.get());
	}

	CellIterator iterator = addCell(source.getType(), std::move(connectStorage), id);

	Cell &cell = (*iterator);
	id = cell.getId();
	cell = std::move(source);
	cell.setId(id);

	return iterator;
}

/*!
	Adds a new cell with the specified id and type.

	\param type is the type of the cell
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, long id)
{
	std::unique_ptr<long[]> connectStorage;
	if (ReferenceElementInfo::hasInfo(type)) {
		int connectSize = ReferenceElementInfo::getInfo(type).nVertices;
		connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	} else {
		connectStorage = std::unique_ptr<long[]>(nullptr);
	}

	return addCell(type, std::move(connectStorage), id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param connectivity is the connectivity of the cell
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, const std::vector<long> &connectivity,
											   long id)
{
	int connectSize = connectivity.size();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	std::copy(connectivity.data(), connectivity.data() + connectSize, connectStorage.get());

	return addCell(type, std::move(connectStorage), id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
											   long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (id < 0) {
		id = generateCellId();
	}

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator = _addInternal(type, std::move(connectStorage), id);

	return iterator;
}

/*!
	Adds a new cell with the specified id and type.

	\param type is the type of the cell
	\param interior is true if the cell is the interior of the patch,
	false otherwise
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, bool interior, long id)
{
	if (!interior) {
		throw std::runtime_error("This legacy function can only be used to add internal cells.");
	}

	return addCell(type, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connectivity is the connectivity of the cell
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, bool interior,
								               const std::vector<long> &connectivity, long id)
{
	if (!interior) {
		throw std::runtime_error("This legacy function can only be used to add internal cells.");
	}

	return addCell(type, connectivity, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, bool interior,
                                               std::unique_ptr<long[]> &&connectStorage, long id)
{
	if (!interior) {
		throw std::runtime_error("This legacy function can only be used to add internal cells.");
	}

	return addCell(type, std::move(connectStorage), id);
}

/*!
	Internal function to add an internal cell.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::_addInternal(ElementType type, std::unique_ptr<long[]> &&connectStorage,
													long id)
{
	// Get the id of the cell before which the new cell should be inserted
#if BITPIT_ENABLE_MPI==1
	//
	// If there are ghosts cells, the internal cell should be inserted
	// before the first ghost cell.
#endif
	long referenceId;
#if BITPIT_ENABLE_MPI==1
	referenceId = m_firstGhostId;
#else
	referenceId = Cell::NULL_ID;
#endif

	// Create the cell
	PiercedVector<Cell>::iterator iterator;
	if (referenceId == Cell::NULL_ID) {
		iterator = m_cells.emreclaim(id, id, type, std::move(connectStorage), true, true);
	} else {
		iterator = m_cells.emreclaimBefore(referenceId, id, id, type, std::move(connectStorage), true, true);
	}
	m_nInternals++;

	// Update the id of the last internal cell
	if (m_lastInternalId < 0) {
		m_lastInternalId = id;
	} else if (m_cells.rawIndex(m_lastInternalId) < m_cells.rawIndex(id)) {
		m_lastInternalId = id;
	}

	return iterator;
}

/*!
	Resore the cell with the specified id.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param type is the type of the cell
	\param connectivity is the connectivity of the cell
	\param id is the id of the cell that will be restored
	\return An iterator pointing to the restored cell.
*/
PatchKernel::CellIterator PatchKernel::restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												   const long &id)
{
	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator = m_cells.find(id);
	if (iterator == m_cells.end()) {
		throw std::runtime_error("Unable to restore the specified cell: the kernel doesn't contain an entry for that cell.");
	}

	_restoreInternal(iterator, type, std::move(connectStorage));

	return iterator;
}

/*!
	Internal function to restore an internal cell.

	\param iterator is an iterator pointing to the cell to restore
	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
void PatchKernel::_restoreInternal(CellIterator iterator, ElementType type,
								   std::unique_ptr<long[]> &&connectStorage)
{
	Cell &cell = *iterator;
	cell.initialize(iterator.getId(), type, std::move(connectStorage), true, true);
	m_nInternals++;
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

	// Get cell information
	const Cell &cell = m_cells[id];

	// Update neighbours
	if (updateNeighs) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			// Update adjacency of the neighbours
			int nFaceAdjacencies = cell.getAdjacencyCount(i);
			const long *faceAdjacencies = cell.getAdjacencies(i);
			for (int k = 0; k < nFaceAdjacencies; ++k) {
				long neighId = faceAdjacencies[k];

				int neighFace, adjacencyId;
				findFaceNeighCell(neighId, id, &neighFace, &adjacencyId);
				if (neighFace >= 0) {
					Cell &neigh = m_cells[neighId];
					neigh.deleteAdjacency(neighFace, adjacencyId);
				}
			} //next k

			// Update interface
			int nFaceInterfaces = cell.getInterfaceCount(i);
			const long *faceInterfaces = cell.getInterfaces(i);
			for (int k = 0; k < nFaceInterfaces; ++k) {
				long interfaceId = faceInterfaces[k];
				Interface &interface = m_interfaces[interfaceId];
				if (interface.getOwner() == id) {
						interface.unsetOwner();
				} else {
						interface.unsetNeigh();
				}
			} //next k
		}
	}

	// Delete cell
#if BITPIT_ENABLE_MPI==1
	bool isInternal = cell.isInterior();
	if (isInternal) {
		_deleteInternal(id, delayed);
	} else {
		_deleteGhost(id, delayed);
	}
#else
	_deleteInternal(id, delayed);
#endif

	// Cell id is no longer used
	m_cellIdGenerator.trash(id);

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
	Internal function to delete an internal cell.

	\param id is the id of the cell
	\param delayed is true a delayed delete will be performed
*/
void PatchKernel::_deleteInternal(long id, bool delayed)
{
	m_cells.erase(id, delayed);
	m_nInternals--;
	if (id == m_lastInternalId) {
		updateLastInternalId();
	}
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
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int j = 0; j < nCellVertices; j++) {
			long vertexId = cellVertexIds[j];
			vertexValence[vertexId] += 1;
		}
	}

	// Loop over cells
	long nOrphanCells = 0;
	for (const Cell &cell : m_cells) {
		long isIsolated = true;
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int j = 0; j < nCellVertices; j++) {
			long vertexId = cellVertexIds[j];
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
	Get the native index of a cell.

	\param id is the id of the cell
	\result The native index of a cell.
*/
long PatchKernel::_getCellNativeIndex(long id) const
{
	return id;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\result All the neighbours of the specified cell.
*/
std::vector<long> PatchKernel::findCellNeighs(const long &id) const
{
	std::vector<long> neighs;
	findCellNeighs(id, &neighs);

	return neighs;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\param[in,out] neighs is the vector were the neighbours will be stored.
	The vector is not cleared before adding the neighbours, it is extended
	by appending all the neighbours found by this function
*/
void PatchKernel::findCellNeighs(const long &id, std::vector<long> *neighs) const
{
	std::vector<long> blackList;
	_findCellNeighs(id, blackList, neighs);
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
	std::vector<long> neighs;
	findCellNeighs(id, codimension, complete, &neighs);

	return neighs;
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
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given codimension. The vector is not cleared before dding
	the neighbours, it is extended by appending all the neighbours found by
	this function
*/
void PatchKernel::findCellNeighs(const long &id, int codimension, bool complete, std::vector<long> *neighs) const
{
	assert(codimension >= 1 && codimension <= getDimension());

	if (codimension == 1) {
		findCellFaceNeighs(id, neighs);
	} else if (codimension == getDimension()) {
		findCellVertexNeighs(id, complete, neighs);
	} else if (codimension == 2) {
		findCellEdgeNeighs(id, complete, neighs);
	}
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\result The neighbours of all the faces of the specified cell.
*/
std::vector<long> PatchKernel::findCellFaceNeighs(const long &id) const
{
	std::vector<long> neighs;
	findCellFaceNeighs(id, &neighs);

	return neighs;
}

/*!
	Extracts all the neighbours of the specified cell

	This implementation can NOT handle hanging nodes.

	\param id is the id of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellNeighs(const long &id, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the vertex
	// count is evaluated using the ReferenceElementInfo associated to the cell.
	int nCellVertices;
	if (m_cells.size() == 0) {
		ElementType cellType = getCellType(id);
		const ReferenceElementInfo &cellTypeInfo = ReferenceElementInfo::getInfo(cellType);
		nCellVertices = cellTypeInfo.nVertices;
	} else {
		const Cell &cell = getCell(id);
		nCellVertices = cell.getVertexCount();
	}

	for (int i = 0; i < nCellVertices; ++i) {
		_findCellVertexNeighs(id, i, blackList, neighs);
	}
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\param[in,out] neighs is the vector were the neighbours of all the faces
	of the specified cell will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::findCellFaceNeighs(const long &id, std::vector<long> *neighs) const
{
	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the face
	// count is evaluated using the ReferenceElementInfo associated to the cell.
	int nCellFaces;
	if (m_cells.size() == 0) {
		ElementType cellType = getCellType(id);
		const ReferenceElementInfo &cellTypeInfo = ReferenceElementInfo::getInfo(cellType);
		nCellFaces = cellTypeInfo.nFaces;
	} else {
		const Cell &cell = getCell(id);
		nCellFaces = cell.getFaceCount();
	}

	// Get the neighbours
	for (int i = 0; i < nCellFaces; ++i) {
		_findCellFaceNeighs(id, i, std::vector<long>(), neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> PatchKernel::findCellFaceNeighs(const long &id, const int &face) const
{
	std::vector<long> neighs;
	_findCellFaceNeighs(id, face, std::vector<long>(), &neighs);

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given face will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::findCellFaceNeighs(const long &id, const int &face, std::vector<long> *neighs) const
{
	_findCellFaceNeighs(id, face, std::vector<long>(), neighs);
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given face will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	const Cell &cell = getCell(id);

	int nFaceAdjacencies = cell.getAdjacencyCount(face);
	const long *faceAdjacencies = cell.getAdjacencies(face);
	for (int k = 0; k < nFaceAdjacencies; ++k) {
		long neighId = faceAdjacencies[k];
		if (utils::findInOrderedVector<long>(neighId, blackList) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
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
	std::vector<long> neighs;
	findCellEdgeNeighs(id, complete, &neighs);

	return neighs;
}

/*!
	Extracts the neighbours of all the edges of the specified cell.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified edge, or should
	contain also neighbours that share an entire face
	\param[in,out] neighs is the vector were the neighbours of all the edges
	of the specified cell will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::findCellEdgeNeighs(const long &id, bool complete, std::vector<long> *neighs) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return;
	}

	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the edge
	// count is evaluated using the ReferenceElementInfo associated to the cell.
	int nCellEdges;
	if (m_cells.size() == 0) {
		ElementType cellType = getCellType(id);
		const ReferenceElementInfo &cellTypeInfo = ReferenceElementInfo::getInfo(cellType);
		nCellEdges = cellTypeInfo.nEdges;
	} else {
		const Cell &cell = getCell(id);
		nCellEdges = cell.getEdgeCount();
	}

	// Get the neighbours
	std::vector<long> blackList;
	if (!complete) {
		findCellFaceNeighs(id, &blackList);
	}

	for (int i = 0; i < nCellEdges; ++i) {
		_findCellEdgeNeighs(id, i, blackList, neighs);
	}
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
	std::vector<long> neighs;
	findCellEdgeNeighs(id, edge, &neighs);

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given edge will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::findCellEdgeNeighs(const long &id, const int &edge, std::vector<long> *neighs) const
{
	_findCellEdgeNeighs(id, edge, std::vector<long>(), neighs);
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given edge will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return;
	}

	const Cell &cell = getCell(id);
	ConstProxyVector<int> edgeVertices = cell.getEdgeLocalConnect(edge);
	std::size_t nEdgeVertices = edgeVertices.size();
	if (nEdgeVertices < 2) {
		return;
	}

	// The neighbours of the edge are the cells that share all the edge vertices
	const int GUESS_NEIGHS_COUNT = 3;

	std::vector<long> firstVertexNeighs;
	firstVertexNeighs.reserve(GUESS_NEIGHS_COUNT);
	_findCellVertexNeighs(id, edgeVertices[0], blackList, &firstVertexNeighs);

	std::vector<long> secondVertexNeighs;
	secondVertexNeighs.reserve(GUESS_NEIGHS_COUNT);
	for (std::size_t k = 1; k < nEdgeVertices; ++k) {
		secondVertexNeighs.clear();
		_findCellVertexNeighs(id, edgeVertices[k], blackList, &secondVertexNeighs);

		for (long neighId : secondVertexNeighs) {
			if (utils::findInOrderedVector<long>(neighId, firstVertexNeighs) != firstVertexNeighs.end()) {
				utils::addToOrderedVector<long>(neighId, *neighs);
			}
		}
	}
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
	std::vector<long> neighs;
	findCellVertexNeighs(id, complete, &neighs);

	return neighs;
}

/*!
	Extracts the neighbours of all the vertices of the specified cell.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\param[in,out] neighs is the vector were the neighbours of all the
	vertices of the specified cell will be stored. The vector is not
	cleared before adding the neighbours, it is extended by appending all
	the neighbours found by this function
*/
void PatchKernel::findCellVertexNeighs(const long &id, bool complete, std::vector<long> *neighs) const
{
	// Some patches can work (at least partially) without initializing the
	// cell list. To handle those patches, if there are no cells the vertex
	// count is evaluated using the ReferenceElementInfo associated to the cell.
	int nCellVertices;
	if (m_cells.size() == 0) {
		ElementType cellType = getCellType(id);
		const ReferenceElementInfo &cellTypeInfo = ReferenceElementInfo::getInfo(cellType);
		nCellVertices = cellTypeInfo.nVertices;
	} else {
		const Cell &cell = getCell(id);
		nCellVertices = cell.getVertexCount();
	}

	// Get the neighbours
	std::vector<long> blackList;
	if (!complete) {
		if (isThreeDimensional()) {
			findCellEdgeNeighs(id, true, &blackList);
		} else {
			findCellFaceNeighs(id, true, &blackList);
		}
	}

	for (int i = 0; i < nCellVertices; ++i) {
		_findCellVertexNeighs(id, i, blackList, neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> PatchKernel::findCellVertexNeighs(const long &id, const int &vertex) const
{
	std::vector<long> neighs;
	findCellVertexNeighs(id, vertex, &neighs);

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::findCellVertexNeighs(const long &id, const int &vertex, std::vector<long> *neighs) const
{
	_findCellVertexNeighs(id, vertex, std::vector<long>(), neighs);
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

	This implementation can NOT handle hanging nodes.

	\param id is the id of the cell
	\param vertex is a local vertex of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	const Cell &cell = getCell(id);
	long vertexId = cell.getVertexId(vertex);

	// Since we are processing a small number of cells it's more efficient to
	// store the list of already processed cells in a vector instead of using
	// a set or an unordered_set. To speed-up the lookup, the vector is kept
	// sorted.
	const int GUESS_NEIGHS_COUNT = 8;

	std::vector<long> alreadyProcessed;
	alreadyProcessed.reserve(GUESS_NEIGHS_COUNT);

	std::vector<long> scanQueue;
	scanQueue.reserve(GUESS_NEIGHS_COUNT);
	scanQueue.push_back(id);

	ConstProxyVector<long> scanCellVertexIds;
	while (!scanQueue.empty()) {
		// Pop a cell to process
		long scanCellId = scanQueue.back();
		const Cell &scanCell = getCell(scanCellId);
		scanQueue.pop_back();
		utils::addToOrderedVector<long>(scanCell.getId(), alreadyProcessed);

		// Get vertex list
		if (scanCell.hasInfo()) {
			scanCellVertexIds = scanCell.getVertexIds();
		}

		// Use face adjacencies to find vertex negihbours
		int nFaces = scanCell.getFaceCount();
		for (int i = 0; i < nFaces; ++i) {
			// Discard faces with no neighbours
            if (scanCell.isFaceBorder(i)) {
				continue;
			}

			// Discard faces that don't own the vertex
			//
			// We use differents algorithm depending on whether the cell has a
			// reference element or not. If the cell has a reference element,
			// accessing the local connectivity of the face is cheap. If the
			// cell has no reference element, it is better to avoid using the
			// local connectivity of the face.
			bool faceOwnsVertex = false;
			if (scanCell.hasInfo()) {
				ConstProxyVector<int> faceLocalVertexIds = scanCell.getFaceLocalVertexIds(i);
				int nFaceVertices = faceLocalVertexIds.size();
				for (int k = 0; k < nFaceVertices; ++k) {
					long faceVertexId = scanCellVertexIds[faceLocalVertexIds[k]];
					if (faceVertexId == vertexId) {
						faceOwnsVertex = true;
						break;
					}
				}
			} else {
				ConstProxyVector<long> faceVertexIds = scanCell.getFaceVertexIds(i);
				int nFaceVertices = faceVertexIds.size();
				for (int k = 0; k < nFaceVertices; ++k) {
					long faceVertexId = faceVertexIds[k];
					if (faceVertexId == vertexId) {
						faceOwnsVertex = true;
						break;
					}
				}
			}

			if (!faceOwnsVertex) {
				continue;
			}

			// Loop through the adjacencies
			//
			// Non-manifold patches may have faces with multiple adjacencies.
			int nFaceAdjacencies = scanCell.getAdjacencyCount(i);
			const long *faceAdjacencies = scanCell.getAdjacencies(i);
			for (int k = 0; k < nFaceAdjacencies; ++k) {
				long faceNeighId = faceAdjacencies[k];

				// Discard neighbours that have already been processed
				if (utils::findInOrderedVector<long>(faceNeighId, alreadyProcessed) != alreadyProcessed.end()) {
					continue;
				}

				// Update list of vertex neighbours
				if (utils::findInOrderedVector<long>(faceNeighId, blackList) == blackList.end()) {
					utils::addToOrderedVector<long>(faceNeighId, *neighs);
				}

				// Update scan list
				scanQueue.push_back(faceNeighId);
			}
		}
	}
}

/*!
	Finds the one-ring of the specified vertex of the cell.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The one-ring of the specified vertex of the cell.
*/
std::vector<long> PatchKernel::findCellVertexOneRing(const long &id, const int &vertex) const
{
	std::vector<long> ring;
	findCellVertexOneRing(id, vertex, &ring);

	return ring;
}

/*!
	Finds the one-ring of the specified vertex of the cell.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param[in,out] ring is the vector were the one-ring of the specified vertex
	of the cell will be stored. The vector is not cleared before adding the
	neighbours, it is extended by appending all the neighbours found by this
	function
*/
void PatchKernel::findCellVertexOneRing(const long &id, const int &vertex, std::vector<long> *ring) const
{
	findCellVertexNeighs(id, vertex, ring);
	utils::addToOrderedVector<long>(id, *ring);
}

/*!
        Find the face and the adjacency that are shared by the specified cells.

        If the two cells don't share a face, the output arguments are set to
        a negative value.

        \param[in] cellId is the id of the cell
        \param[in] neighId is the id of the neighbour
        \param[out] cellFace is the face of the cell that is shared with the
        neighbour. If the cell doesn't share any face with the neighbour, the
        argument is set to a negative number
        \param[out] cellAdjacencyId is the adjacency of the cell that is shared
        with the neighbour. If the cell doesn't share any adjacency with the
        neighbour, the argument is set to a negative number
        \return Return true if the two cells share a face, false otherwise.
*/
bool PatchKernel::findFaceNeighCell(long cellId, long neighId, int *cellFace, int *cellAdjacencyId) const
{
    const Cell &cell = getCell(cellId);

    int nFaces = cell.getFaceCount();
    for (int i = 0; i < nFaces; ++i) {
        int nFaceAdjacencies = cell.getAdjacencyCount(i);
        const long *faceAdjacencies = cell.getAdjacencies(i);
        for (int k = 0; k < nFaceAdjacencies; ++k) {
            long faceNeighId = faceAdjacencies[k];
            if (faceNeighId < 0) {
                continue;
            } else if (faceNeighId == neighId) {
                *cellFace        = i;
                *cellAdjacencyId = k;
                return true;
            }
        }
    }

    *cellFace        = -1;
    *cellAdjacencyId = -1;

    return false;
}

/*!
 * Get the PIDs of the internal cells.
 *
 * \return The PIDs of the internal cells.
 */
std::set<int> PatchKernel::getInternalPIDs()
{
	std::set<int> list;
	CellConstIterator endItr = internalConstEnd();
	for (CellConstIterator itr = internalConstBegin(); itr != endItr; ++itr) {
		list.insert(itr->getPID());
	}

	return list;
}

/*!
 * Get all the internal cells which belongs to the specified PID.
 *
 * \param pid is the PID
 * \return All the internal cells which belongs to the specified PID.
 */
std::vector<long> PatchKernel::getInternalsByPID(int pid)
{
	std::vector<long> cells;
	CellConstIterator endItr = internalConstEnd();
	for (CellConstIterator itr = internalConstBegin(); itr != endItr; ++itr) {
		if (itr->getPID() == pid){
			cells.push_back(itr.getId());
		}
	}

	return cells;
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
	Gets a constant reference to the interfaces owned by the patch.

	\return A constant reference to the interfaces owned by the patch.
*/
const PiercedVector<Interface> & PatchKernel::getInterfaces() const
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
ElementType PatchKernel::getInterfaceType(const long &id) const
{
	return m_interfaces[id].getType();
}

/*!
	Returns an iterator pointing to the specified interface.

	\result An iterator to the specified interface.
*/
PatchKernel::InterfaceIterator PatchKernel::getInterfaceIterator(const long &id)
{
	return m_interfaces.find(id);
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
	Returns a constant iterator pointing to the specified interface.

	\result A constant iterator to the specified interface.
*/
PatchKernel::InterfaceConstIterator PatchKernel::getInterfaceConstIterator(const long &id) const
{
	return m_interfaces.find(id);
}

/*!
	Returnsa a constant iterator pointing to the first interface.

	\result A constant iterator to the first interface.
*/
PatchKernel::InterfaceConstIterator PatchKernel::interfaceConstBegin() const
{
	return m_interfaces.cbegin();
}

/*!
	Returns a constant iterator pointing to last interface.

	\result A constant iterator to the last interface.
*/
PatchKernel::InterfaceConstIterator PatchKernel::interfaceConstEnd() const
{
	return m_interfaces.cend();
}

/*!
 * Generates a new unique id for the interfaces.
 *
 * \result A new unique id for the interfaces.
 */
long PatchKernel::generateInterfaceId()
{
	if (!isExpert()) {
		return Interface::NULL_ID;
	}

	return m_interfaceIdGenerator.generate();
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the newly created interface.
	If a negative id value is specified, a new unique id will be generated
	for the interface
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(const Interface &source, long id)
{
	if (id < 0) {
		id = generateInterfaceId();
	}

	Interface interface = source;

	return addInterface(std::move(interface), id);
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the newly created interface.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added interface.

*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(Interface &&source, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	int connectSize = source.getConnectSize();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);

	InterfaceIterator iterator = addInterface(source.getType(), std::move(connectStorage), id);

	Interface &interface = (*iterator);
	id = interface.getId();
	interface = std::move(source);
	interface.setId(id);

	return iterator;
}

/*!
	Adds a new interface with the specified id.

	\param type is the type of the interface
	\param id is the id that will be assigned to the newly created interface.
	If a negative id value is specified, a new unique id will be generated
	for the interface
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(ElementType type, long id)
{
	int connectSize = ReferenceElementInfo::getInfo(type).nVertices;
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);

	return addInterface(type, std::move(connectStorage), id);
}

/*!
	Adds a new interface with the specified id.

	\param type is the type of the interface
	\param connectivity is the connectivity of the interface
	\param id is the id of the new interface. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(ElementType type,
														 const std::vector<long> &connectivity,
														 long id)
{
	int connectSize = connectivity.size();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	std::copy(connectivity.data(), connectivity.data() + connectSize, connectStorage.get());

	return addInterface(type, std::move(connectStorage), id);
}

/*!
	Adds a new interface with the specified id.

	\param type is the type of the interface
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(ElementType type,
														 std::unique_ptr<long[]> &&connectStorage,
														 long id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	if (id < 0) {
		id = generateInterfaceId();
	}

	if (Interface::getDimension(type) > (getDimension() - 1)) {
		return interfaceEnd();
	}

	PiercedVector<Interface>::iterator iterator = m_interfaces.emreclaim(id, id, type, std::move(connectStorage));

	return iterator;
}

/*!
	Resore the interface with the specified id.

	The kernel should already contain the interface, only the contents of the
	interface will be updated.

	\param type is the type of the interface
	\param connectivity is the connectivity of the interface
	\param id is the id of the interface to restore
	\return An iterator pointing to the restored interface.
*/
PatchKernel::InterfaceIterator PatchKernel::restoreInterface(ElementType type,
															 std::unique_ptr<long[]> &&connectStorage,
															 const long &id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	if (Interface::getDimension(type) > (getDimension() - 1)) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = m_interfaces.find(id);
	if (iterator == m_interfaces.end()) {
		throw std::runtime_error("Unable to restore the specified interface: the kernel doesn't contain an entry for that interface.");
	}

	Interface &interface = *iterator;
	interface.initialize(id, type, std::move(connectStorage));

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
		//
		// If the owner has been deleted before the interface, it may be null.
		long ownerId = interface.getOwner();
		if (ownerId >= 0) {
			Cell &owner = m_cells[ownerId];
			int ownerFace = interface.getOwnerFace();
			const long *ownerFaceInterfaces = owner.getInterfaces(ownerFace);

			int ownerInterfaceId = 0;
			while (ownerFaceInterfaces[ownerInterfaceId] != id) {
				++ownerInterfaceId;
			}
			owner.deleteInterface(ownerFace, ownerInterfaceId);
		}

		// Update neighbour
		long neighId = interface.getNeigh();
		if (neighId >= 0) {
			Cell &neigh = m_cells[neighId];
			int neighFace = interface.getNeighFace();
			const long *neighFaceInterfaces = neigh.getInterfaces(neighFace);

			int neighInterfaceId = 0;
			while (neighFaceInterfaces[neighInterfaceId] != id) {
				++neighInterfaceId;
			}
			neigh.deleteInterface(neighFace, neighInterfaceId);
		}
	}

	// Delete interface
	m_interfaces.erase(id, delayed);
	m_interfaceIdGenerator.trash(id);

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
		m_interfaces.flush();
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

	An interface is considered orphan if it has no owner nor neighbour or
	if it's on a border face of a ghost cell.

	\return The number of orphan interfaces.
*/
long PatchKernel::countOrphanInterfaces() const
{
	long nOrphanInterfaces = 0;
	for (InterfaceConstIterator itr = interfaceConstBegin(); itr != interfaceConstEnd(); ++itr) {
		const long interfaceId = itr.getId();
		if (isInterfaceOrphan(interfaceId)) {
			++nOrphanInterfaces;
		}
	}

	return nOrphanInterfaces;
}

/*!
	Find orphan interfaces in the patch.

	An interface is considered orphan if it has no owner nor neighbour or
	if it's on a border face of a ghost cell.

	\result The list of orphan interfaces.
*/
std::vector<long> PatchKernel::findOrphanInterfaces() const
{
	std::vector<long> orphanInterfaces;
	for (InterfaceConstIterator itr = interfaceConstBegin(); itr != interfaceConstEnd(); ++itr) {
		const long interfaceId = itr.getId();
		if (isInterfaceOrphan(interfaceId)) {
			orphanInterfaces.push_back(interfaceId);
		}
	}

	return orphanInterfaces;
}

/*!
	Remove orphan interfaces
*/
bool PatchKernel::deleteOrphanInterfaces()
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> list = findOrphanInterfaces();
	deleteInterfaces(list);

	return true;
}

/*!
	Check if the interfaces is orphan.

	An interface is considered orphan if it has no owner nor neighbour or
	if it's on a border face of a ghost cell.

	\result Returns true if the interfaces is orphan, false otherwise.
*/
bool PatchKernel::isInterfaceOrphan(long id) const
{
	const Interface &interface = getInterface(id);

	// Interface is not orphan if it has an owner or a neighbour
	long ownerId = interface.getOwner();
	long neighId = interface.getNeigh();
	if (ownerId >= 0 && neighId >= 0) {
		return false;
	}

	// Interface is not orphan if it's on a border face of an internal cell
	if (ownerId >= 0 && neighId < 0) {
		const Cell &owner = getCell(ownerId);
		if (owner.isInterior()) {
			return false;
		}
	}

	return true;
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
 *  Write the vertices to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PatchKernel::dumpVertices(std::ostream &stream) const
{
	// Dump kernel
	m_vertices.dumpKernel(stream);

	// Dump vertices
	for (const Vertex &vertex : m_vertices) {
		utils::binary::write(stream, vertex.getId());

		std::array<double, 3> coords = vertex.getCoords();
		utils::binary::write(stream, coords[0]);
		utils::binary::write(stream, coords[1]);
		utils::binary::write(stream, coords[2]);
	}
}

/*!
 *  Restore the vertices from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void PatchKernel::restoreVertices(std::istream &stream)
{
	// Restore kernel
	m_vertices.restoreKernel(stream);

	// Enable advanced editing
	bool originalExpertStatus = isExpert();
	setExpert(true);

	// Restore vertices
	long nVertices = m_vertices.size();
	for (long i = 0; i < nVertices; ++i) {
		long id;
		utils::binary::read(stream, id);

		std::array<double, 3> coords;
		utils::binary::read(stream, coords[0]);
		utils::binary::read(stream, coords[1]);
		utils::binary::read(stream, coords[2]);

		restoreVertex(std::move(coords), id);
	}

	// Set original advanced editing status
	setExpert(originalExpertStatus);
}

/*!
 *  Write the cells to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PatchKernel::dumpCells(std::ostream &stream) const
{
	// Dump kernel
	m_cells.dumpKernel(stream);

	// Dump cells
	for (const Cell &cell: m_cells) {
		utils::binary::write(stream, cell.getId());
		utils::binary::write(stream, cell.getPID());
		utils::binary::write(stream, cell.getType());
#if BITPIT_ENABLE_MPI==1
		utils::binary::write(stream, getCellRank(cell.getId()));
#else
		int dummyRank = 0;
		utils::binary::write(stream, dummyRank);
#endif

		int cellConnectSize = cell.getConnectSize();
		utils::binary::write(stream, cellConnectSize);

		const long *cellConnect = cell.getConnect();
		for (int i = 0; i < cellConnectSize; ++i) {
			utils::binary::write(stream, cellConnect[i]);
		}
	}

	// Dump ghost/internal subdivision
#if BITPIT_ENABLE_MPI==1
	utils::binary::write(stream, m_firstGhostId);
	utils::binary::write(stream, m_lastInternalId);
#else
	utils::binary::write(stream, Cell::NULL_ID);
	utils::binary::write(stream, Cell::NULL_ID);
#endif
}

/*!
 *  Restore the cells from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void PatchKernel::restoreCells(std::istream &stream)
{
	// Restore kernel
	m_cells.restoreKernel(stream);

	// Enable advanced editing
	bool originalExpertStatus = isExpert();
	setExpert(true);

	// Restore cells
	long nCells = m_cells.size();
	for (long i = 0; i < nCells; ++i) {
		long id;
		utils::binary::read(stream, id);

		int PID;
		utils::binary::read(stream, PID);

		ElementType type;
		utils::binary::read(stream, type);

#if BITPIT_ENABLE_MPI==1
		int rank;
		utils::binary::read(stream, rank);
#else
		int dummyRank;
		utils::binary::read(stream, dummyRank);
#endif

		int cellConnectSize;
		utils::binary::read(stream, cellConnectSize);

		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[cellConnectSize]);
		for (int k = 0; k < cellConnectSize; ++k) {
			utils::binary::read(stream, cellConnect[k]);
		}

		CellIterator iterator;
#if BITPIT_ENABLE_MPI==1
		iterator = restoreCell(type, std::move(cellConnect), rank, id);
#else
		iterator = restoreCell(type, std::move(cellConnect), id);
#endif
		iterator->setPID(PID);
	}

	// Restore ghost/internal subdivision
#if BITPIT_ENABLE_MPI==1
	utils::binary::read(stream, m_firstGhostId);
	utils::binary::read(stream, m_lastInternalId);
#else
	long dummyFirstGhostId;
	long dummyLastInternalId;
	utils::binary::read(stream, dummyFirstGhostId);
	utils::binary::read(stream, dummyLastInternalId);
#endif

	// Restore adjacencies
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_AUTOMATIC) {
		buildAdjacencies();
	}

#if BITPIT_ENABLE_MPI==1
	// Update information for ghost data exchange
	if (isPartitioned()) {
		updateGhostExchangeInfo();
	}
#endif

	// Set original advanced editing status
	setExpert(originalExpertStatus);
}

/*!
 *  Write the interfaces to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PatchKernel::dumpInterfaces(std::ostream &stream) const
{
	// Early return if the interfaces are not build
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Dump kernel
	m_interfaces.dumpKernel(stream);

	// Dump interfaces
	for (const Interface &interface : getInterfaces()) {
		utils::binary::write(stream, interface.getId());

		utils::binary::write(stream, interface.getOwner());
		utils::binary::write(stream, interface.getOwnerFace());

		long neighId = interface.getNeigh();
		utils::binary::write(stream, interface.getNeigh());
		if (neighId >= 0) {
			utils::binary::write(stream, interface.getNeighFace());
		}

		utils::binary::write(stream, interface.getPID());
	}
}

/*!
 *  Restore the interfaces from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void PatchKernel::restoreInterfaces(std::istream &stream)
{
	// Early return if the interfaces are not build
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Restore kernel
	m_interfaces.restoreKernel(stream);

	// Enable advanced editing
	bool originalExpertStatus = isExpert();
	setExpert(true);

	// Restore interfaces
	long nInterfaces = m_interfaces.size();
	for (long n = 0; n < nInterfaces; ++n) {
		long interfaceId;
		utils::binary::read(stream, interfaceId);

		long ownerId;
		int ownerFace;
		utils::binary::read(stream, ownerId);
		utils::binary::read(stream, ownerFace);
		Cell *owner = &(m_cells.at(ownerId));

		long neighId;
		int neighFace;
		utils::binary::read(stream, neighId);
		Cell *neigh;
		if (neighId >= 0) {
			utils::binary::read(stream, neighFace);
			neigh = &(m_cells.at(neighId));
		} else {
			neighFace = -1;
			neigh = nullptr;
		}

		int pid;
		utils::binary::read(stream, pid);

		InterfaceIterator interfaceIterator = buildCellInterface(owner, ownerFace, neigh, neighFace, interfaceId);
		interfaceIterator->setPID(pid);
	}

	// Set original advanced editing status
	setExpert(originalExpertStatus);
}

/*!
	Updates the id of the last internal cell.
*/
void PatchKernel::updateLastInternalId()
{
	if (m_nInternals == 0) {
		m_lastInternalId = Cell::NULL_ID;
	} else {
		m_lastInternalId = m_cells.getSizeMarker(m_nInternals - 1, Cell::NULL_ID);
	}
}

/*!
	Generates the patch.

	Default implementation is a no-op function.

	\param trackSpawn if set to true the changes to the patch will be tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the spawn.
*/
std::vector<adaption::Info> PatchKernel::_spawn(bool trackSpawn)
{
	BITPIT_UNUSED(trackSpawn);

	assert(false && "The patch needs to implement _spawn");

	return std::vector<adaption::Info>();
}

/*!
	Prepares the patch for performing the adaption.

	Default implementation is a no-op function.

	\param trackAdaption if set to true the function will return the changes
	that will be performed in the alter step
	\result If the adaption is tracked, returns a vector of adaption::Info that
	can be used to discover what changes will be performed in the alter step,
	otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_adaptionPrepare(bool trackAdaption)
{
	BITPIT_UNUSED(trackAdaption);

	return std::vector<adaption::Info>();
}

/*!
	Alter the patch performing the adpation.

	Default implementation is a no-op function.

	\param trackAdaption if set to true the function will return the changes
	done to the patch during the adaption
	\result If the adaption is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the adaption, otherwise an
	empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_adaptionAlter(bool trackAdaption)
{
	BITPIT_UNUSED(trackAdaption);

	assert(false && "The patch needs to implement _adaptionAlter");

	return std::vector<adaption::Info>();
}

/*!
	Cleanup patch data structured after the adaption.

	Default implementation is a no-op function.
*/
void PatchKernel::_adaptionCleanup()
{
}

/*!
	Marks a cell for refinement.

	Default implementation is a no-op function.

	\param id the cell to be refined
	\result Returns true if the marker was properly set, false otherwise.
*/
bool PatchKernel::_markCellForRefinement(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Marks a cell for coarsening.

	Default implementation is a no-op function.

	\param id the cell to be refined
	\result Returns true if the marker was properly set, false otherwise.
*/
bool PatchKernel::_markCellForCoarsening(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Resets the adaption marker of the specified cell.

	Default implementation is a no-op function.

	\param id the cell to be refined
	\result Returns true if the marker was properly reset, false otherwise.
*/
bool PatchKernel::_resetCellAdaptionMarker(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Returns the adaption marker of the specified cell.

	Default implementation always return an undefined marker.

	\param id is the id of the cell
	\return The adaption marker of the cell.
*/
adaption::Marker PatchKernel::_getCellAdaptionMarker(const long &id)
{
	BITPIT_UNUSED(id);

	return adaption::MARKER_UNDEFINED;
}

/*!
	Enables cell balancing.

	Default implementation is a no-op function.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
	\result Returns true if the falg was properly set, false otherwise.
*/
bool PatchKernel::_enableCellBalancing(const long &id, bool enabled)
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(enabled);

	return false;
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
	Sorts cell storage in ascending id order.

	Internal cells and ghost cells are sordered separately.
*/
bool PatchKernel::sortCells()
{
	if (!isExpert()) {
		return false;
	}

	// Sort internal cells
	if (m_nInternals > 0) {
		m_cells.sortBefore(m_lastInternalId, true);
		updateLastInternalId();
	}

#if BITPIT_ENABLE_MPI==1
	// Sort ghost cells
	if (m_nGhosts > 0) {
		m_cells.sortAfter(m_firstGhostId, true);
		updateFirstGhostId();
	}
#endif

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
std::array<double, 3> PatchKernel::evalCellCentroid(const long &id) const
{
	const Cell &cell = getCell(id);

	return evalElementCentroid(cell);
}

/*!
	Evaluates the bounding box of the specified cell.

	\param id is the id of the cell
	\param[out] minPoint is the minimum point of the bounding box
	\param[out] maxPoint is the maximum point of the bounding box
*/
void PatchKernel::evalCellBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const
{
	const Cell &cell = getCell(id);

	return evalElementBoundingBox(cell, minPoint, maxPoint);
}

/*!
	Get vertex coordinates of the specified cell.

	\param id is the id of the cell
	\param[in,out] externalStorage is an optional external storage that can
	be used to store vertex coordinates of standard elements (i.e., elements
	associated to a reference element)
	\result Vertex coordinates of the cell.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getCellVertexCoordinates(long id, std::array<double, 3> *externalStorage) const
{
	const Cell &cell = getCell(id);

	return getElementVertexCoordinates(cell, externalStorage);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> PatchKernel::evalInterfaceCentroid(const long &id) const
{
	const Interface &interface = getInterface(id);

	return evalElementCentroid(interface);
}

/*!
	Evaluates the bounding box of the specified interface.

	\param id is the id of the interface
	\param[out] minPoint is the minimum point of the bounding box
	\param[out] maxPoint is the maximum point of the bounding box
*/
void PatchKernel::evalInterfaceBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const
{
	const Interface &interface = getInterface(id);

	return evalElementBoundingBox(interface, minPoint, maxPoint);
}


/*!
	Get vertex coordinates of the specified interface.

	\param id is the id of the interface
	\param[in,out] externalStorage is an optional external storage that can
	be used to store vertex coordinates of standard elements (i.e., elements
	associated to a reference element)
	\result Vertex coordinates of the interface.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getInterfaceVertexCoordinates(long id, std::array<double, 3> *externalStorage) const
{
	const Interface &interface = getInterface(id);

	return getElementVertexCoordinates(interface, externalStorage);
}

/*!
	Evaluates the centroid of the specified element.

	Element centroid is computed as the arithmetic average of element
	vertex coordinates.

	\param element is the element
	\result The centroid of the specified element.
*/
std::array<double, 3> PatchKernel::evalElementCentroid(const Element &element) const
{
	ConstProxyVector<long> cellVertexIds = element.getVertexIds();
	int nElementVertices = cellVertexIds.size();

	std::array<double, 3> centroid = {{0., 0., 0.}};
	for (int i = 0; i < nElementVertices; ++i) {
		const Vertex &vertex = getVertex(cellVertexIds[i]);
		centroid += vertex.getCoords();
	}
	centroid /= (double) nElementVertices;

	return centroid;
}

/*!
	Evaluates the bounding box of the specified element.

	\param element is the element
	\param[out] minPoint is the minimum point of the bounding box
	\param[out] maxPoint is the maximum point of the bounding box
*/
void PatchKernel::evalElementBoundingBox(const Element &element, std::array<double,3> *minCoord, std::array<double,3> *maxCoord) const
{
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	const int nElementVertices = elementVertexIds.size();

	*minCoord = getVertexCoords(elementVertexIds[0]);
	*maxCoord = *minCoord;
	for (int i = 1; i < nElementVertices; ++i) {
		const std::array<double, 3> &vertexCoord = getVertexCoords(elementVertexIds[i]);
		for (int d = 0; d < 3; ++d) {
			(*minCoord)[d] = std::min(vertexCoord[d], (*minCoord)[d]);
			(*maxCoord)[d] = std::max(vertexCoord[d], (*maxCoord)[d]);
		}
	}
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
 * Check whether the face "face_A" on cell "cell_A" is the same as the face
 * "face_B" on cell "cell_B".
 * 
 * \param[in] cellId_A is the index of the first cell
 * \param[in] face_A is the face on the first cell
 * \param[in] cellId_B is the index of the second cell
 * \param[in] face_B is the face on the second cell
 * \result Returns true if the two faces are the same.
*/
bool PatchKernel::isSameFace(long cellId_A, int face_A, long cellId_B, int face_B)
{
	const Cell &cell_A = m_cells[cellId_A];
	const Cell &cell_B = m_cells[cellId_B];
	if (cell_A.getFaceType(face_A) != cell_B.getFaceType(face_B)) {
		return false;
	}

	int nFaceVertices = cell_A.getFaceVertexCount(face_A);

	std::vector<long> faceVertexIds_A(nFaceVertices);
	std::vector<long> faceVertexIds_B(nFaceVertices);
	for (int k = 0; k < nFaceVertices; ++k) {
		faceVertexIds_A[k] = cell_A.getFaceVertexId(face_A, k);
		faceVertexIds_B[k] = cell_B.getFaceVertexId(face_B, k);
	}

	std::sort(faceVertexIds_A.begin(), faceVertexIds_A.end());
	std::sort(faceVertexIds_B.begin(), faceVertexIds_B.end());

	return (faceVertexIds_A == faceVertexIds_B);
}

/*!
	Returns the current adjacencies build strategy.

	\return The current adjacencies build strategy.
*/
PatchKernel::AdjacenciesBuildStrategy PatchKernel::getAdjacenciesBuildStrategy() const
{
	return m_adjacenciesBuildStrategy;
}

/*!
	Set the current adjacencies build strategy.

	\param status is the adjacencies build strategy that will be set
*/
void PatchKernel::setAdjacenciesBuildStrategy(AdjacenciesBuildStrategy status)
{
	m_adjacenciesBuildStrategy = status;
}

/*!
	Fill adjacencies info for each cell.
*/
void PatchKernel::buildAdjacencies()
{
	// Reset adjacencies
	if (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE) {
		clearAdjacencies();
	}

	// Update the adjacencies
	updateAdjacencies(m_cells.getIds(false));
}

/*!
	Update the adjacencies of the specified list of cells and of their
	neighbours.

	This implementation can NOT handle hanging nodes.

	\param[in] cellIds is the list of cell ids
*/
void PatchKernel::updateAdjacencies(const std::vector<long> &cellIds)
{
	// The adjacencies are found looking for matching half-faces.
	//
	// On a three-dimensional patch each internal face is shared between two,
	// and only two, half-faces. On lower dimension patches one internal face
	// may be shared among multiple half-faces (this happend with non-manifold
	// patches).
	//
	// If faces can be shared only between two half-faces, there is a one-to-one
	// correspondence between the half-faces. Give one half-faces, its matching
	// can be found looking for an half-faces with the same vertices but with
	// reverse vertex winding order.
	//
	// If faces can be shared among more than two half-faces, it's not anymore
	// possible to define a unique winding order for the vertices of the faces.
	// There is no more a one-to-one correspondence between pairs of half-face.
	// In this case, when looking for half-face correspondence, it is necessary
	// to look for half-faces with reverse vertex winding order and also for
	// half-faces with the same vertex winding order.
	//
	// If we are updating the adjacencies of only part of the cells, we need
	// to populate the half-face list with the faces of the cells that will
	// not be updated (among those faces may find matches for the faces of
	// the updated cells). If faces can be shared only between two half-faces
	// we need to insert only the border faces of the non-updated cells,
	// otherwise we need to add all the faces of the non-updated cells.

	// Define matching windings
	bool multipleMatchesAllowed = (getDimension() < 3);

	std::vector<CellHalfFace::Winding> matchingWindings;
	matchingWindings.push_back(CellHalfFace::WINDING_REVERSE);
	if (multipleMatchesAllowed) {
		matchingWindings.push_back(CellHalfFace::WINDING_NATURAL);
	}

	// Initialize half-faces list
	std::unordered_multiset<CellHalfFace, CellHalfFace::Hasher> halfFaces;
	if (multipleMatchesAllowed) {
		halfFaces.reserve(4 * getCellCount());
	} else {
		halfFaces.reserve(2 * cellIds.size());
	}

	// Populate list with faces of non-updated cells
	if ((long) cellIds.size() != getCellCount()) {
		std::unordered_set<long> updateCellSet(cellIds.begin(), cellIds.end());

		for (Cell &cell : m_cells) {
			if (updateCellSet.count(cell.getId()) > 0) {
				continue;
			}

			int nCellFaces = cell.getFaceCount();
			for (int face = 0; face < nCellFaces; face++) {
				if (multipleMatchesAllowed || cell.isFaceBorder(face)) {
					halfFaces.emplace(cell, face, CellHalfFace::WINDING_NATURAL);
				}
			}
		}
	}

	// Update the adjacencies
	for (long cellId : cellIds) {
		Cell &cell = m_cells[cellId];

		const int nCellFaces = cell.getFaceCount();
		for (int face = 0; face < nCellFaces; face++) {
			// Generate the half-face
			CellHalfFace halfFace(cell, face);

			// Find matching half-faces
			int nAdjacencies = 0;
			for (CellHalfFace::Winding winding : matchingWindings) {
				// Set winding order
				halfFace.setWinding(winding);

				// Search for matching half-faces
				//
				// Each match is adjacent to the current half-face.
				auto matchRange = halfFaces.equal_range(halfFace);
				for (auto matchItr = matchRange.first; matchItr != matchRange.second; ++matchItr) {
					const CellHalfFace &neighHalfFace = *matchItr;

					// Generate the adjacency
					Cell &neigh    = neighHalfFace.getCell();
					long neighId   = neigh.getId();
					long neighFace = neighHalfFace.getFace();

					cell.pushAdjacency(face, neighId);
					neigh.pushAdjacency(neighFace, cellId);

					++nAdjacencies;
				}

				// Remove the matching half-faces from the list
				//
				// It is not possible to remove the matching half-faces
				// if multiple matches are allowed.
				if (!multipleMatchesAllowed) {
					halfFaces.erase(matchRange.first, matchRange.second);
				}
			}

			// Add the current half-face to the list.
			//
			// If no multiple matches are allowed, we need to add the current
			// half-face to the list only if no matchings were found.
			if (multipleMatchesAllowed || nAdjacencies == 0) {
				halfFace.setWinding(CellHalfFace::WINDING_NATURAL);
				halfFaces.insert(std::move(halfFace));
			}
		}
	}

	// Set adjacencies build strategy
	setAdjacenciesBuildStrategy(ADJACENCIES_AUTOMATIC);
}

/*!
	Clear the adjacencies.
*/
void PatchKernel::clearAdjacencies()
{
	// Early return if no adjacencies have been built
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_NONE) {
		return;
	}

	// Reset adjacencies
	for (Cell &cell : m_cells) {
		cell.resetAdjacencies();
	}

	// Set adjacencies status
	setAdjacenciesBuildStrategy(ADJACENCIES_NONE);
}

/*!
	Returns the current interfaces build strategy.

	\return The current interfaces build strategy.
*/
PatchKernel::InterfacesBuildStrategy PatchKernel::getInterfacesBuildStrategy() const
{
	return m_interfacesBuildStrategy;
}

/*!
	Set the current interfaces build strategy.

	\param status is the interfaces build strategy that will be set
*/
void PatchKernel::setInterfacesBuildStrategy(InterfacesBuildStrategy status)
{
	m_interfacesBuildStrategy = status;
}

/*!
	Build interfaces among the cells.
*/
void PatchKernel::buildInterfaces()
{
	// Reset interfaces
	if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
		clearInterfaces();
	}

	// Update interfaces
	updateInterfaces(m_cells.getIds(false));
}

/*!
	Update the interfaces of the specified list of cells and of their
	neighbours.

	\param[in] cellIds is the list of cell ids
*/
void PatchKernel::updateInterfaces(const std::vector<long> &cellIds)
{
	// Enable advanced editing
	bool originalExpertStatus = isExpert();
	setExpert(true);

	//
	// Update interfaces
	//
	// Adjacencies and interfaces are paired: the i-th adjacency correspondes
	// to the i-th interface. Moreover if we loop through the adjacencies of
	// a face, the adjacencies that have an interface are always listed first.
	// This meas that, to update the interfaces, we can count the interfaces
	// already associated to a face and loop only on the adjacencies which
	// have an index past the one of the last interface.
	//
	// On border faces of internal cells we need to build an interface, also
	// if there are no adjacencies.
	for (long cellId : cellIds) {
		Cell &cell = m_cells[cellId];
		const int nCellFaces = cell.getFaceCount();
		for (int face = 0; face < nCellFaces; face++) {
			bool isFaceBorder = cell.isFaceBorder(face);
			if (!isFaceBorder) {
				// Find the range of adjacencies that need an interface
				int updateEnd   = cell.getAdjacencyCount(face);
				int updateBegin = cell.getInterfaceCount(face);
				if (updateBegin == updateEnd) {
					continue;
				}

				// Build an interface for every adjacency
				const long *faceAdjacencies = cell.getAdjacencies(face);
				for (int k = updateBegin; k < updateEnd; ++k) {
					long neighId = faceAdjacencies[k];
					Cell *neigh  = &m_cells[neighId];

					int neighFace = findAdjoinNeighFace(cellId, neighId);

					buildCellInterface(&cell, face, neigh, neighFace);
				}
			} else if (cell.isInterior()) {
				// Internal borderes need an interface
				buildCellInterface(&cell, face, nullptr, -1);
			}
		}
	}

	// Set interfaces build strategy
	setInterfacesBuildStrategy(INTERFACES_AUTOMATIC);

	// Set original advanced editing status
	setExpert(originalExpertStatus);
}

/*!
	Clear the interfaces of the patch.

	This function just calls 'resetInterfaces', its only purpose is to make
	the API for handling the interfaces as much as possible similar to the
	API for handling the adjacencied.
*/
void PatchKernel::clearInterfaces()
{
	resetInterfaces();
}

/*!
	Given two cells, build the interface between them.

	After creating the interface, the data structures of the cells are update
	to take into account the newly created interface.

	\param cell_1 is the first cell
	\param face_1 is the face of the first cell
	\param cell_2 is the second cell
	\param face_2 is the face of the second cell
	\param interfaceId is the id of the interface that will be built. If a
	negative id value is specified, a new interface will be created, otherwise
	the existing interface will be overwritten with the new data.
 */
PatchKernel::InterfaceIterator PatchKernel::buildCellInterface(Cell *cell_1, int face_1, Cell *cell_2, int face_2, long interfaceId)
{
	// Validate the input
	assert(cell_1);
	assert(face_1 >= 0);

	if (cell_2) {
		assert(face_2 >= 0);
	} else {
        assert(face_2 < 0);
	}

	// Get the id of the first cell
	long id_1 = cell_1->getId();

	// Get the id of the second cell
	long id_2 = Cell::NULL_ID;
	if (cell_2) {
		id_2 = cell_2->getId();
	}

	// Owner and neighbour of the interface
	//
	// The interface is owned by the cell that has only one adjacency, i.e.,
	// by the cell that owns the smallest of the two faces. If the faces
	// of both cells have the same size, the interface is owned by the cell
	// with the "lower fuzzy positioning". It is not necessary to have a
	// precise comparison, it's only necassary to define a repetible order
	// between the two cells. It is therefore possible to use the "fuzzy"
	// cell comparison.
	bool cellOwnsInterface = true;
	if (cell_2) {
		if (cell_1->getAdjacencyCount(face_1) > 1) {
			cellOwnsInterface = false;
		} else if (cell_2->getAdjacencyCount(face_2) == 1) {
			assert(cell_1->getAdjacencyCount(face_1) == 1);
			cellOwnsInterface = CellFuzzyPositionLess(*this)(id_1, id_2);
		}
	}

	Cell *intrOwner;
	Cell *intrNeigh;
	long intrOwnerId;
	long intrNeighId;
	int intrOwnerFace;
	int intrNeighFace;
	if (cellOwnsInterface) {
		intrOwnerId   = id_1;
		intrOwner     = cell_1;
		intrOwnerFace = face_1;

		if (cell_2) {
			intrNeighId   = id_2;
			intrNeigh     = cell_2;
			intrNeighFace = face_2;
		} else {
			intrNeighId   = Cell::NULL_ID;
			intrNeigh     = nullptr;
			intrNeighFace = -1;
		}
	} else {
		intrOwnerId   = id_2;
		intrOwner     = cell_2;
		intrOwnerFace = face_2;

		intrNeighId   = id_1;
		intrNeigh     = cell_1;
		intrNeighFace = face_1;
	}

	// Connectivity of the interface
	ConstProxyVector<long> faceConnect = intrOwner->getFaceConnect(intrOwnerFace);

	int nInterfaceVertices = faceConnect.size();
	std::unique_ptr<long[]> interfaceConnect = std::unique_ptr<long[]>(new long[nInterfaceVertices]);
	for (int k = 0; k < nInterfaceVertices; ++k) {
		interfaceConnect[k] = faceConnect[k];
	}

	// Create the interface
	Interface *interface;
	InterfaceIterator interfaceIterator;

	ElementType interfaceType = intrOwner->getFaceType(intrOwnerFace);
	if (interfaceId < 0) {
		interfaceIterator = addInterface(interfaceType, std::move(interfaceConnect), interfaceId);
		interface = &(*interfaceIterator);
		interfaceId = interface->getId();
	} else {
		interfaceIterator = restoreInterface(interfaceType, std::move(interfaceConnect), interfaceId);
		interface = &(*interfaceIterator);
	}

	// Set owner and neighbour
	interface->setOwner(intrOwnerId, intrOwnerFace);
	if (intrNeighId >= 0) {
		interface->setNeigh(intrNeighId, intrNeighFace);
	}

	// Update owner and neighbour cell data
	//
	// Adjacencies and interfaces are paired: the i-th adjacency correspondes
	// to the i-th interface. Moreover if we loop through the adjacencies of
	// a face, the adjacencies that have an interface are always listed first.
	// When the interfaces are fully built, each adjacency is associated to an
	// interface, however this will not be true during the generation of the
	// interfaces. If the adjacencies already associated to an interface are
	// listed first, it will be easy to detect which interfaces are still to
	// be built: the missing interfaces are associated to the adjacencies which
	// have an index past the last interface.
	//
	// The above only matters if the neighbour cell exists, if there is no
	// neighbour there will be only one interface on that face, so there are
	// no ordering issues.
	intrOwner->pushInterface(intrOwnerFace, interfaceId);
	if (intrNeigh) {
		intrNeigh->pushInterface(intrNeighFace, interfaceId);

		// Fix adjacency order on the owner cell
		int ownerInterfaceIndex   = intrOwner->getInterfaceCount(intrOwnerFace) - 1;
		long ownerPairedAdjacency = intrOwner->getAdjacency(intrOwnerFace, ownerInterfaceIndex);
		if (ownerPairedAdjacency != intrNeighId) {
			int ownerPairedAdjacencyIndex = intrOwner->findAdjacency(intrOwnerFace, intrNeighId);
			intrOwner->setAdjacency(intrOwnerFace, ownerInterfaceIndex, intrNeighId);
			intrOwner->setAdjacency(intrOwnerFace, ownerPairedAdjacencyIndex, ownerPairedAdjacency);
		}

		// Fix adjacency order on the neighbour cell
		int neighInterfaceIndex   = intrNeigh->getInterfaceCount(intrNeighFace) - 1;
		long neighPairedAdjacency = intrNeigh->getAdjacency(intrNeighFace, neighInterfaceIndex);
		if (neighPairedAdjacency != intrOwnerId) {
			int neighPairedAdjacencyIndex = intrNeigh->findAdjacency(intrNeighFace, intrOwnerId);
			intrNeigh->setAdjacency(intrNeighFace, neighInterfaceIndex, intrOwnerId);
			intrNeigh->setAdjacency(intrNeighFace, neighPairedAdjacencyIndex, neighPairedAdjacency);
		}
	}

	return interfaceIterator;
}

/*!
	Given a cell and one of it's neighbours, finds the faces of the neighbour
	that adjoins the specified cell.

	The function doesn't check if the two cells are really neighbours.

	\param cellId is the id of the cell
	\param neighId is the id of a neighbour of the cell
	\result The face of the neighbour which adjoin the specified cell
 */
int PatchKernel::findAdjoinNeighFace(const long &cellId, const long &neighId) const
{
	const Cell &neigh = m_cells[neighId];
	const int nNeighFaces = neigh.getFaceCount();
	for (int face = 0; face < nNeighFaces; face++) {
		int nFaceAdjacencies = neigh.getAdjacencyCount(face);
		const long *faceAdjacencies = neigh.getAdjacencies(face);
		for (int k = 0; k < nFaceAdjacencies; ++k) {
			long geussId = faceAdjacencies[k];
			if (geussId == cellId) {
				return face;
			}
		}
	}

	return -1;
}

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
void PatchKernel::getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const
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
		MPI_Allreduce(const_cast<bool *>(&m_boxDirty), &isDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
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
	Get the coordinates of the specified element

	\param element is the element
	\param[in,out] externalStorage is an optional external storage that can
	be used to store vertex coordinates of standard elements (i.e., elements
	associated to a reference element)
	\result The coordinates of the element.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getElementVertexCoordinates(const Element &element, std::array<double, 3> *externalStorage) const
{
	std::unique_ptr<std::vector<std::array<double, 3>>> localStorage;
	std::array<double, 3> *storage;

	// Get the vertices ids
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	const int nElementVertices = elementVertexIds.size();

	// Store coordinates in storage
	bool useExternalStorage = (externalStorage != nullptr) && element.hasInfo();
	if (useExternalStorage) {
		storage = externalStorage;
	} else {
		localStorage = std::unique_ptr<std::vector<std::array<double, 3>>>(new std::vector<std::array<double, 3>>(nElementVertices));
		storage = localStorage->data();
	}

	for (int i = 0; i < nElementVertices; ++i) {
		storage[i] = getVertex(elementVertexIds[i]).getCoords();
	}

	// Build the proxy vector with the coordinates
	ConstProxyVector<std::array<double, 3>> vertexCoordinates;
	if (useExternalStorage) {
		vertexCoordinates.set(externalStorage, nElementVertices);
	} else {
		vertexCoordinates.set(std::move(*localStorage));
		localStorage.release();
	}

	return vertexCoordinates;
}

/*!
	Sort patch vertices on regular bins.

	\param[in] nBins (default = 128) is the number of bins (on each space
	direction)
	\result Returns the bin index associated to each vertex.
*/
std::unordered_map<long, long> PatchKernel::binSortVertex(int nBins)
{
	return PatchKernel::binSortVertex(m_vertices, nBins);
}

/*!
    Sort specified vertices on regular bins.

    \param[in] vertices are the vertices to be sorted
    \param[in] nBins (default = 128) is the number of bins (on each space
    direction)
    \result Returns the bin index associated to each vertex.
*/
std::unordered_map<long, long> PatchKernel::binSortVertex(const PiercedVector<Vertex> &vertices, int nBins)
{
    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variables
    double                              dx, dy, dz;

    // Counters
    long                                i, j, k;

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
    PiercedVector<Vertex>::const_iterator E = vertices.cend();
    for (PiercedVector<Vertex>::const_iterator V = vertices.cbegin(); V != E; ++V) {
        i = std::min(nBins - 1L, long((V->getCoords()[0] - m_boxMinPoint[0]) / dx));
        j = std::min(nBins - 1L, long((V->getCoords()[1] - m_boxMinPoint[1]) / dy));
        k = std::min(nBins - 1L, long((V->getCoords()[2] - m_boxMinPoint[2]) / dz));
        bin_index[V->getId()] = nBins * nBins * k + nBins * j + i;
    }

    return bin_index;
}

/*!
    Group vertices on regular bins.

    \param[in] vertices are the vertices to be sorted
    \param[in] nBins is the number of bins (on each space direction)
    \result Returns the vertices grouped into bins.
*/
std::unordered_map<long, std::vector<long>> PatchKernel::binGroupVertices(int nBins)
{
	return PatchKernel::binGroupVertices(m_vertices, nBins);
}

/*!
    Group specified vertices on regular bins.

    \param[in] vertices are the vertices to be sorted
    \param[in] nBins is the number of bins (on each space direction)
    \result Returns the vertices grouped into bins.
*/
std::unordered_map<long, std::vector<long>> PatchKernel::binGroupVertices(const PiercedVector<Vertex> &vertices, int nBins)
{
    // Update bounding box
    updateBoundingBox();

    // Bin's spacing
    double dx = max(1.0e-12, m_boxMaxPoint[0] - m_boxMinPoint[0]) / ((double) nBins);
    double dy = max(1.0e-12, m_boxMaxPoint[1] - m_boxMinPoint[1]) / ((double) nBins);
    double dz = max(1.0e-12, m_boxMaxPoint[2] - m_boxMinPoint[2]) / ((double) nBins);

    // Identify bins of vertices
    std::unordered_map<long, std::vector<long>> bins;
    PiercedVector<Vertex>::const_iterator E = vertices.cend();
    for (PiercedVector<Vertex>::const_iterator V = vertices.cbegin(); V != E; ++V) {
        int i = std::min(nBins - 1L, long((V->getCoords()[0] - m_boxMinPoint[0]) / dx));
        int j = std::min(nBins - 1L, long((V->getCoords()[1] - m_boxMinPoint[1]) / dy));
        int k = std::min(nBins - 1L, long((V->getCoords()[2] - m_boxMinPoint[2]) / dz));

        long binId = nBins * nBins * k + nBins * j + i;
        bins[binId].emplace_back(V->getId());
    }

    return bins;
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
			ConstProxyVector<long> faceConnect = cell.getFaceConnect(i);
			int nFaceVertices = faceConnect.size();

			std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);
			for (int j = 0; j < nFaceVertices; ++j) {
				long vertexId = faceConnect[j];

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
			ElementType faceType = cell.getFaceType(i);
			envelope.addCell(faceType, std::move(faceEnvelopeConnect));
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
	Get the VTK object.

	\result The VTK object.
*/
VTKUnstructuredGrid & PatchKernel::getVTK()
{
	return m_vtk;
}

/*!
	Get the VTK write target.

	\result The VTK write target.
*/
PatchKernel::WriteTarget PatchKernel::getVTKWriteTarget() const
{
	return m_vtkWriteTarget;
}

/*!
	Set the VTK write target.

	\param writeTarget is the VTK write target.
*/
void PatchKernel::setVTKWriteTarget(WriteTarget writeTarget)
{
	m_vtkWriteTarget = writeTarget;
}

/*!
	Get the VTK cell write range.

	\result The VTK cell write range.
*/
const PatchKernel::CellConstRange PatchKernel::getVTKCellWriteRange() const
{
	if (m_vtkWriteTarget == WRITE_TARGET_CELLS_ALL) {
		return CellConstRange(cellConstBegin(), cellConstEnd());
#if BITPIT_ENABLE_MPI==1
	} else if (m_vtkWriteTarget == WRITE_TARGET_CELLS_INTERNAL) {
		return CellConstRange(internalConstBegin(), internalConstEnd());
#endif
	} else {
		return CellConstRange(cellConstEnd(), cellConstEnd());
	}
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] name is the name of the data to be written. Either user
 *  data or patch data
 *  @param[in] format is the format which must be used. Supported options
 *  are "ascii" or "appended". For "appended" type an unformatted binary
 *  stream must be used
 */
void PatchKernel::flushData(std::fstream &stream, const std::string &name, VTKFormat format)
{
	assert(format == VTKFormat::APPENDED);
	BITPIT_UNUSED(format);

	if (name == "Points") {
		for (VertexConstIterator itr = vertexConstBegin(); itr != vertexConstEnd(); ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap.rawAt(vertexRawId);
			if (vertexVTKId != Vertex::NULL_ID) {
				const Vertex &vertex = m_vertices.rawAt(vertexRawId);
				genericIO::flushBINARY(stream, vertex.getCoords());
			}
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (const Cell &cell : getVTKCellWriteRange()) {
			offset += cell.getVertexCount();
			genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "types") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			VTKElementType VTKType;
			switch (cell.getType())  {

			case ElementType::VERTEX:
				VTKType = VTKElementType::VERTEX;
				break;

			case ElementType::LINE:
				VTKType = VTKElementType::LINE;
				break;

			case ElementType::TRIANGLE:
				VTKType = VTKElementType::TRIANGLE;
				break;

			case ElementType::PIXEL:
				VTKType = VTKElementType::PIXEL;
				break;

			case ElementType::QUAD:
				VTKType = VTKElementType::QUAD;
				break;

			case ElementType::POLYGON:
				VTKType = VTKElementType::POLYGON;
				break;

			case ElementType::TETRA:
				VTKType = VTKElementType::TETRA;
				break;

			case ElementType::VOXEL:
				VTKType = VTKElementType::VOXEL;
				break;

			case ElementType::HEXAHEDRON:
				VTKType = VTKElementType::HEXAHEDRON;
				break;

			case ElementType::WEDGE:
				VTKType = VTKElementType::WEDGE;
				break;

			case ElementType::PYRAMID:
				VTKType = VTKElementType::PYRAMID;
				break;

			case ElementType::POLYHEDRON:
				VTKType = VTKElementType::POLYHEDRON;
				break;

			default:
				VTKType = VTKElementType::UNDEFINED;
				break;

			}

			genericIO::flushBINARY(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			const int nCellVertices = cell.getVertexCount();
			for (int k = 0; k < nCellVertices; ++k) {
				long vertexId = cellVertexIds[k];
				long vtkVertexId = m_vtkVertexMap.at(vertexId);
				genericIO::flushBINARY(stream, vtkVertexId);
			}
		}
	} else if (name == "faces") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				genericIO::flushBINARY(stream, (long) 0);
			} else {
				std::vector<long> faceStream = cell.getFaceStream();
				Cell::renumberFaceStream(m_vtkVertexMap, &faceStream);
				int faceStreamSize = faceStream.size();
				for (int k = 0; k < faceStreamSize; ++k) {
					genericIO::flushBINARY(stream, faceStream[k]);
				}
			}
		}
	} else if (name == "faceoffsets") {
		int offset = 0;
		for (const Cell &cell : getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				offset += 1;
			} else {
				offset += cell.getFaceStreamSize();
			}

			genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "cellIndex") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			genericIO::flushBINARY(stream, cell.getId());
		}
	} else if (name == "PID") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			genericIO::flushBINARY(stream, cell.getPID());
		}
	} else if (name == "vertexIndex") {
		for (VertexConstIterator itr = vertexConstBegin(); itr != vertexConstEnd(); ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap.rawAt(vertexRawId);
			if (vertexVTKId != Vertex::NULL_ID) {
				std::size_t vertexId = itr.getId();
				genericIO::flushBINARY(stream, vertexId);
			}
		}
#if BITPIT_ENABLE_MPI==1
	} else if (name == "cellGlobalIndex") {
		PatchNumberingInfo numberingInfo(this);
		for (const Cell &cell : getVTKCellWriteRange()) {
			genericIO::flushBINARY(stream, numberingInfo.getCellGlobalId(cell.getId()));
		}
	} else if (name == "rank") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			genericIO::flushBINARY(stream, getCellRank(cell.getId()));
		}
#endif
	}
}

/*!
 *  Renumbers vertices ID consecutively, starting from a given offset.
 *
 *  \param[in] offset is the starting id
 */
void PatchKernel::consecutiveRenumberVertices(long offset)
{
	// Renumber vertices
	std::unordered_map<long, long > map = consecutiveItemRenumbering(m_vertices, offset);
	
	// Renumber cell connectivity
	for(Cell &cell : m_cells) {
		cell.renumberVertices(map);
	}

	// Renumber interface connectivity
	for(Interface &interface : getInterfaces()) {
		interface.renumberVertices(map);
	}
}	

/*!
 *  Renumbers cells consecutively, starting from a given offset.
 *
 *  \param[in] offset is the starting id
 */
void PatchKernel::consecutiveRenumberCells(long offset)
{
	// Renumber cells
	std::unordered_map<long, long > map = consecutiveItemRenumbering(m_cells, offset);
	
	// Renumber cell adjacencies
	for (auto &cell: m_cells) {
		long *adjacencies = cell.getAdjacencies();
		int nCellAdjacencies = cell.getAdjacencyCount();
		for (int i = 0; i < nCellAdjacencies; ++i) {
			long &neighId = adjacencies[i];
			neighId = map[neighId];
		}
	}
	
	// Renumber interface owner and neighbour
	for (Interface &interface: m_interfaces) {
		long ownerId = interface.getOwner();
		int ownerFace = interface.getOwnerFace();
		interface.setOwner(map.at(ownerId), ownerFace);

		long neighId = interface.getNeigh();
		if (neighId >= 0) {
			int neighFace = interface.getNeighFace();
			interface.setNeigh(map.at(neighId), neighFace);
		}
	}

	// Renumber last internal and first ghost markers
	if (m_lastInternalId >= 0) {
		m_lastInternalId = map.at(m_lastInternalId);
	}

#if BITPIT_ENABLE_MPI==1
	if (m_firstGhostId >= 0) {
		m_firstGhostId = map.at(m_firstGhostId);
	}
#endif

#if BITPIT_ENABLE_MPI==1
	// Update information for ghost data exchange
	if (isPartitioned()) {
		updateGhostExchangeInfo();
	}
#endif
}	

/*!
 *  Renumbers interfaces consecutively, starting from a given offset.
 *
 *  \param[in] offset is the starting id
 */
void PatchKernel::consecutiveRenumberInterfaces(long offset)
{
	// Renumber interfaces
	std::unordered_map<long, long > map = consecutiveItemRenumbering(m_interfaces, offset);
	
	// Renumber cell interfaces
	for (Cell &cell: m_cells) {
		long *interfaces = cell.getInterfaces();
		int nCellInterfaces = cell.getInterfaceCount();
		for (int i = 0; i < nCellInterfaces; ++i) {
			long &interfaceId = interfaces[i];
			interfaceId = map.at(interfaceId);
		}
	}
}

/*!
 *  Renumbering vertices, cells and interfaces consecutively, starting from
 *  given offsets.
 *
 *  \param[in] vertexOffset is the starting id of the vertices
 *  \param[in] cellOffset is the starting id of the cells
 *  \param[in] interfaceOffset is the starting id of the interfaces
 */
void PatchKernel::consecutiveRenumber(long vertexOffset, long cellOffset, long interfaceOffset)
{
	consecutiveRenumberVertices(vertexOffset);
	consecutiveRenumberCells(cellOffset);
	consecutiveRenumberInterfaces(interfaceOffset);
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int PatchKernel::getDumpVersion() const
{
	const int KERNEL_DUMP_VERSION = 9;

	return (KERNEL_DUMP_VERSION + _getDumpVersion());
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PatchKernel::dump(std::ostream &stream) const
{
	// Version
	utils::binary::write(stream, getDumpVersion());

	// Generic information
	utils::binary::write(stream, m_id);
	utils::binary::write(stream, m_dimension);
	utils::binary::write(stream, m_vtk.getName());
#if BITPIT_ENABLE_MPI==1
	utils::binary::write(stream, isPartitioned());
	utils::binary::write(stream, m_haloSize);
#else
	utils::binary::write(stream, false);
	utils::binary::write(stream, 0);
#endif

	// Spawn status
	utils::binary::write(stream, m_spawnStatus);

	// Adaption status
	utils::binary::write(stream, m_adaptionStatus);

	// Partition status
#if BITPIT_ENABLE_MPI==1
	utils::binary::write(stream, m_partitioningStatus);
#else
	utils::binary::write(stream, PARTITIONING_UNSUPPORTED);
#endif

	// Adjacencies build strategy
	utils::binary::write(stream, m_adjacenciesBuildStrategy);

	// Dump interfaces build strategy
	utils::binary::write(stream, m_interfacesBuildStrategy);

	// VTK data
	utils::binary::write(stream, m_vtkWriteTarget);

	// Specific dump
	_dump(stream);

	// Geometric tolerance
	utils::binary::write(stream, (int) m_hasCustomTolerance);
	if (m_hasCustomTolerance) {
		utils::binary::write(stream, m_tolerance);
	}

	// Index generators
	m_vertexIdGenerator.dump(stream);
	m_cellIdGenerator.dump(stream);
	m_interfaceIdGenerator.dump(stream);
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 *  \param reregister is true the patch will be unregistered and then
 *  registered again using the id found in the binary archive
 */
void PatchKernel::restore(std::istream &stream, bool reregister)
{
	// Reset the patch
	reset();

	// Version
	int version;
	utils::binary::read(stream, version);
	if (version != getDumpVersion()) {
		throw std::runtime_error ("The version of the file does not match the required version");
	}

	// Id
	int id;
	utils::binary::read(stream, id);
	if (reregister) {
		patch::manager().unregisterPatch(this);
		patch::manager().registerPatch(this, id);
	}

	// Dimension
	int dimension;
	utils::binary::read(stream, dimension);
	setDimension(dimension);

	// Name
	std::string name;
	utils::binary::read(stream, name);
	m_vtk.setName(name);

	// Partioned flag
	bool partitioned;
	utils::binary::read(stream, partitioned);
#if BITPIT_ENABLE_MPI==1
	setPartitioned(partitioned);
#endif

	// Halo size
#if BITPIT_ENABLE_MPI==1
	utils::binary::read(stream, m_haloSize);
#else
	int dummyHaloSize;
	utils::binary::read(stream, dummyHaloSize);
#endif

	// Spawn status
	utils::binary::read(stream, m_spawnStatus);

	// Adaption status
	utils::binary::read(stream, m_adaptionStatus);

	// Partition status
#if BITPIT_ENABLE_MPI==1
	utils::binary::read(stream, m_partitioningStatus);
#else
	PartitioningStatus dummyPartitioningStatus;
	utils::binary::read(stream, dummyPartitioningStatus);
#endif

	// Adjacencies build strategy
	utils::binary::read(stream, m_adjacenciesBuildStrategy);

	// Interfaces status
	utils::binary::read(stream, m_interfacesBuildStrategy);

	// VTK data
	utils::binary::read(stream, m_vtkWriteTarget);

	// Specific restore
	_restore(stream);

#if BITPIT_ENABLE_MPI==1
	// Update information for ghost data exchange
	if (isPartitioned()) {
		updateGhostExchangeInfo();
	}
#endif

	// Geometric tolerance
	int hasCustomTolerance;
	utils::binary::read(stream, hasCustomTolerance);
	if (hasCustomTolerance) {
		double tolerance;
		utils::binary::read(stream, tolerance);
		setTol(tolerance);
	} else {
		resetTol();
	}

	// Index generators
	m_vertexIdGenerator.restore(stream);
	m_cellIdGenerator.restore(stream);
	m_interfaceIdGenerator.restore(stream);
}

/*!
 *  Merge the specified adaption info.
 *
 *  \param[in] source is the source adaption info
 *  \param[in,out] destination is the destination adaption info
 */
void PatchKernel::mergeAdaptionInfo(std::vector<adaption::Info> &&source, std::vector<adaption::Info> &destination)
{
	if (source.empty()) {
		return;
	} else if (destination.empty()) {
		destination.swap(source);
		return;
	}

	throw std::runtime_error ("Unable to merge the adaption info.");
}

}
