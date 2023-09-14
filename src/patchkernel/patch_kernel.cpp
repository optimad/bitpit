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

#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif

#include "bitpit_CG.hpp"
#include "bitpit_common.hpp"

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

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	Patches that are filled automatically (e.g. VolOctree) will initialize
	the cells only on the process identified by the rank zero in the
	communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(MPI_Comm communicator, std::size_t haloSize, bool expert)
#else
/*!
	Creates a patch.

	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(bool expert)
#endif
	: m_expert(expert)
{
	// Initialize the patch
#if BITPIT_ENABLE_MPI==1
	initialize(communicator, haloSize);
#else
	initialize();
#endif

	// Register the patch
	patch::manager().registerPatch(this);
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	Patches that are filled automatically (e.g. VolOctree) will initialize
	the cells only on the process identified by the rank zero in the
	communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(int dimension, bool expert)
#endif
	: m_expert(expert)
{
	// Initialize the patch
#if BITPIT_ENABLE_MPI==1
	initialize(communicator, haloSize);
#else
	initialize();
#endif

	// Register the patch
	patch::manager().registerPatch(this);

	// Set the dimension
	setDimension(dimension);
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	Patches that are filled automatically (e.g. VolOctree) will initialize
	the cells only on the process identified by the rank zero in the
	communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
PatchKernel::PatchKernel(int id, int dimension, bool expert)
#endif
	: m_expert(expert)
{
	// Initialize the patch
#if BITPIT_ENABLE_MPI==1
	initialize(communicator, haloSize);
#else
	initialize();
#endif

	// Register the patch
	patch::manager().registerPatch(this, id);

	// Set the dimension
	//
	// Here we can only call the base function, if needed, every derived patch
	// has to call its derived function.
	PatchKernel::setDimension(dimension);
}

/*!
	Copy constructor.

	\param other is another patch whose content is copied into this
*/
PatchKernel::PatchKernel(const PatchKernel &other)
    : VTKBaseStreamer(other),
      m_vertices(other.m_vertices),
      m_cells(other.m_cells),
      m_interfaces(other.m_interfaces),
      m_alteredCells(other.m_alteredCells),
      m_alteredInterfaces(other.m_alteredInterfaces),
      m_nInternalVertices(other.m_nInternalVertices),
#if BITPIT_ENABLE_MPI==1
      m_nGhostVertices(other.m_nGhostVertices),
#endif
      m_lastInternalVertexId(other.m_lastInternalVertexId),
#if BITPIT_ENABLE_MPI==1
      m_firstGhostVertexId(other.m_firstGhostVertexId),
#endif
      m_nInternalCells(other.m_nInternalCells),
#if BITPIT_ENABLE_MPI==1
      m_nGhostCells(other.m_nGhostCells),
#endif
      m_lastInternalCellId(other.m_lastInternalCellId),
#if BITPIT_ENABLE_MPI==1
      m_firstGhostCellId(other.m_firstGhostCellId),
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
      m_toleranceCustom(other.m_toleranceCustom),
      m_tolerance(other.m_tolerance)
#if BITPIT_ENABLE_MPI==1
      , m_partitioningStatus(other.m_partitioningStatus),
      m_owner(other.m_owner),
      m_haloSize(other.m_haloSize),
      m_partitioningCellsTag(other.m_partitioningCellsTag),
      m_partitioningVerticesTag(other.m_partitioningVerticesTag),
      m_partitioningSerialization(other.m_partitioningSerialization),
      m_partitioningOutgoings(other.m_partitioningOutgoings),
      m_partitioningGlobalExchanges(other.m_partitioningGlobalExchanges),
      m_partitioningInfoDirty(other.m_partitioningInfoDirty),
      m_ghostVertexInfo(other.m_ghostVertexInfo),
      m_ghostVertexExchangeTargets(other.m_ghostVertexExchangeTargets),
      m_ghostVertexExchangeSources(other.m_ghostVertexExchangeSources),
      m_ghostCellInfo(other.m_ghostCellInfo),
      m_ghostCellExchangeTargets(other.m_ghostCellExchangeTargets),
      m_ghostCellExchangeSources(other.m_ghostCellExchangeSources)
#endif
{
#if BITPIT_ENABLE_MPI==1
	// Initialize the communicator
	initializeCommunicator(other.getCommunicator());
#else
	// Initialize serial communicator
	initializeSerialCommunicator();
#endif

	// Create index generators
	importVertexIndexGenerator(other);
	importInterfaceIndexGenerator(other);
	importCellIndexGenerator(other);

	// Register the patch
	patch::manager().registerPatch(this);

	// Update the VTK streamer
	replaceVTKStreamer(&other, this);
}

/*!
	Move constructor.

	\param other is another patch whose content is moved into this
*/
PatchKernel::PatchKernel(PatchKernel &&other)
    : VTKBaseStreamer(std::move(other)),
      m_vertices(std::move(other.m_vertices)),
      m_cells(std::move(other.m_cells)),
      m_interfaces(std::move(other.m_interfaces)),
      m_alteredCells(std::move(other.m_alteredCells)),
      m_alteredInterfaces(std::move(other.m_alteredInterfaces)),
      m_vertexIdGenerator(std::move(other.m_vertexIdGenerator)),
      m_interfaceIdGenerator(std::move(other.m_interfaceIdGenerator)),
      m_cellIdGenerator(std::move(other.m_cellIdGenerator)),
      m_nInternalVertices(std::move(other.m_nInternalVertices)),
#if BITPIT_ENABLE_MPI==1
      m_nGhostVertices(std::move(other.m_nGhostVertices)),
#endif
      m_lastInternalVertexId(std::move(other.m_lastInternalVertexId)),
#if BITPIT_ENABLE_MPI==1
      m_firstGhostVertexId(std::move(other.m_firstGhostVertexId)),
#endif
      m_nInternalCells(std::move(other.m_nInternalCells)),
#if BITPIT_ENABLE_MPI==1
      m_nGhostCells(std::move(other.m_nGhostCells)),
#endif
      m_lastInternalCellId(std::move(other.m_lastInternalCellId)),
#if BITPIT_ENABLE_MPI==1
      m_firstGhostCellId(std::move(other.m_firstGhostCellId)),
#endif
      m_vtk(std::move(other.m_vtk)),
      m_vtkWriteTarget(std::move(other.m_vtkWriteTarget)),
      m_vtkVertexMap(std::move(other.m_vtkVertexMap)),
      m_boxFrozen(std::move(other.m_boxFrozen)),
      m_boxDirty(std::move(other.m_boxDirty)),
      m_boxMinPoint(std::move(other.m_boxMinPoint)),
      m_boxMaxPoint(std::move(other.m_boxMaxPoint)),
      m_boxMinCounter(std::move(other.m_boxMinCounter)),
      m_boxMaxCounter(std::move(other.m_boxMaxCounter)),
      m_adjacenciesBuildStrategy(std::move(other.m_adjacenciesBuildStrategy)),
      m_interfacesBuildStrategy(std::move(other.m_interfacesBuildStrategy)),
      m_spawnStatus(std::move(other.m_spawnStatus)),
      m_adaptionStatus(std::move(other.m_adaptionStatus)),
      m_expert(std::move(other.m_expert)),
      m_id(std::move(other.m_id)),
      m_dimension(std::move(other.m_dimension)),
      m_toleranceCustom(std::move(other.m_toleranceCustom)),
      m_tolerance(std::move(other.m_tolerance)),
      m_rank(std::move(other.m_rank)),
      m_nProcessors(std::move(other.m_nProcessors))
#if BITPIT_ENABLE_MPI==1
      , m_communicator(std::move(MPI_COMM_NULL)),
      m_partitioningStatus(std::move(other.m_partitioningStatus)),
      m_owner(std::move(other.m_owner)),
      m_haloSize(std::move(other.m_haloSize)),
      m_partitioningCellsTag(std::move(other.m_partitioningCellsTag)),
      m_partitioningVerticesTag(std::move(other.m_partitioningVerticesTag)),
      m_partitioningSerialization(std::move(other.m_partitioningSerialization)),
      m_partitioningOutgoings(std::move(other.m_partitioningOutgoings)),
      m_partitioningGlobalExchanges(std::move(other.m_partitioningGlobalExchanges)),
      m_partitioningInfoDirty(std::move(other.m_partitioningInfoDirty)),
      m_ghostVertexInfo(std::move(other.m_ghostVertexInfo)),
      m_ghostVertexExchangeTargets(std::move(other.m_ghostVertexExchangeTargets)),
      m_ghostVertexExchangeSources(std::move(other.m_ghostVertexExchangeSources)),
      m_ghostCellInfo(std::move(other.m_ghostCellInfo)),
      m_ghostCellExchangeTargets(std::move(other.m_ghostCellExchangeTargets)),
      m_ghostCellExchangeSources(std::move(other.m_ghostCellExchangeSources))
#endif
{
	// Handle patch regstration
	patch::manager().unregisterPatch(&other);

	patch::manager().registerPatch(this, m_id);
	patch::manager().registerPatch(&other);

	// Update the VTK streamer
	//
	// Pointers to VTK streamers has been copied, we need to replace all the
	// pointer to the other object with pointer to this object.
	replaceVTKStreamer(&other, this);

#if BITPIT_ENABLE_MPI==1
	// Handle the communication
	std::swap(m_communicator, other.m_communicator);
#endif
}

/**
	Move assignment operator.

	\param other is another patch whose content is copied into this
*/
PatchKernel & PatchKernel::operator=(PatchKernel &&other)
{
	VTKBaseStreamer::operator=(std::move(other));
	m_vertices = std::move(other.m_vertices);
	m_cells = std::move(other.m_cells);
	m_interfaces = std::move(other.m_interfaces);
	m_alteredCells = std::move(other.m_alteredCells);
	m_alteredInterfaces = std::move(other.m_alteredInterfaces);
	m_vertexIdGenerator = std::move(other.m_vertexIdGenerator);
	m_interfaceIdGenerator = std::move(other.m_interfaceIdGenerator);
	m_cellIdGenerator = std::move(other.m_cellIdGenerator);
	m_nInternalVertices = std::move(other.m_nInternalVertices);
#if BITPIT_ENABLE_MPI==1
	m_nGhostVertices = std::move(other.m_nGhostVertices);
#endif
	m_lastInternalVertexId = std::move(other.m_lastInternalVertexId);
#if BITPIT_ENABLE_MPI==1
	m_firstGhostVertexId = std::move(other.m_firstGhostVertexId);
#endif
	m_nInternalCells = std::move(other.m_nInternalCells);
#if BITPIT_ENABLE_MPI==1
	m_nGhostCells = std::move(other.m_nGhostCells);
#endif
	m_lastInternalCellId = std::move(other.m_lastInternalCellId);
#if BITPIT_ENABLE_MPI==1
	m_firstGhostCellId = std::move(other.m_firstGhostCellId);
#endif
	m_vtk = std::move(other.m_vtk);
	m_vtkWriteTarget = std::move(other.m_vtkWriteTarget);
	m_vtkVertexMap = std::move(other.m_vtkVertexMap);
	m_boxFrozen = std::move(other.m_boxFrozen);
	m_boxDirty = std::move(other.m_boxDirty);
	m_boxMinPoint = std::move(other.m_boxMinPoint);
	m_boxMaxPoint = std::move(other.m_boxMaxPoint);
	m_boxMinCounter = std::move(other.m_boxMinCounter);
	m_boxMaxCounter = std::move(other.m_boxMaxCounter);
	m_adjacenciesBuildStrategy = std::move(other.m_adjacenciesBuildStrategy);
	m_interfacesBuildStrategy = std::move(other.m_interfacesBuildStrategy);
	m_spawnStatus = std::move(other.m_spawnStatus);
	m_adaptionStatus = std::move(other.m_adaptionStatus);
	m_expert = std::move(other.m_expert);
	m_id = std::move(other.m_id);
	m_dimension = std::move(other.m_dimension);
	m_toleranceCustom = std::move(other.m_toleranceCustom);
	m_tolerance = std::move(other.m_tolerance);
	m_rank = std::move(other.m_rank);
	m_nProcessors = std::move(other.m_nProcessors);
#if BITPIT_ENABLE_MPI==1
	m_communicator = std::move(MPI_COMM_NULL);
	m_partitioningStatus = std::move(other.m_partitioningStatus);
	m_owner = std::move(other.m_owner);
	m_haloSize = std::move(other.m_haloSize);
	m_partitioningCellsTag = std::move(other.m_partitioningCellsTag);
	m_partitioningVerticesTag = std::move(other.m_partitioningVerticesTag);
	m_partitioningSerialization = std::move(other.m_partitioningSerialization);
	m_partitioningOutgoings = std::move(other.m_partitioningOutgoings);
	m_partitioningGlobalExchanges = std::move(other.m_partitioningGlobalExchanges);
	m_partitioningInfoDirty = std::move(other.m_partitioningInfoDirty);
	m_ghostVertexInfo = std::move(other.m_ghostVertexInfo);
	m_ghostVertexExchangeTargets = std::move(other.m_ghostVertexExchangeTargets);
	m_ghostVertexExchangeSources = std::move(other.m_ghostVertexExchangeSources);
	m_ghostCellInfo = std::move(other.m_ghostCellInfo);
	m_ghostCellExchangeTargets = std::move(other.m_ghostCellExchangeTargets);
	m_ghostCellExchangeSources = std::move(other.m_ghostCellExchangeSources);
#endif

	// Handle patch regstration
	patch::manager().unregisterPatch(this);
	patch::manager().unregisterPatch(&other);

	patch::manager().registerPatch(this, m_id);
	patch::manager().registerPatch(&other);

	// Update the VTK streamer
	//
	// Pointers to VTK streamers has been moved, we need to replace all the
	// pointer to the other object with pointer to this object.
	replaceVTKStreamer(&other, this);

#if BITPIT_ENABLE_MPI==1
	// Handle the communication
	std::swap(m_communicator, other.m_communicator);
#endif

	return *this;
}

/*!
	Initialize the patch
*/
#if BITPIT_ENABLE_MPI==1
/*!
	\param communicator is the communicator to be used for exchanging data
	among the processes
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::initialize(MPI_Comm communicator, std::size_t haloSize)
#else
void PatchKernel::initialize()
#endif
{
	// Id
	m_id = PatchManager::AUTOMATIC_ID;

	// Vertex count
	m_nInternalVertices = 0;
#if BITPIT_ENABLE_MPI==1
	m_nGhostVertices = 0;
#endif

	m_lastInternalVertexId = Vertex::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostVertexId = Vertex::NULL_ID;
#endif

	// Cell count
	m_nInternalCells = 0;
#if BITPIT_ENABLE_MPI==1
	m_nGhostCells    = 0;
#endif

	m_lastInternalCellId = Cell::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostCellId   = Cell::NULL_ID;
#endif

	// Dimension
	m_dimension = -1;

	// Index generators
	setVertexAutoIndexing(true);
	setInterfaceAutoIndexing(true);
	setCellAutoIndexing(true);

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

#if BITPIT_ENABLE_MPI==1
	// Initialize communicator
	initializeCommunicator(communicator);

	// Set halo size
	initializeHaloSize(haloSize);

	// Mark patch as partioned
	if (isPartitioned()) {
		setPartitioningStatus(PARTITIONING_CLEAN);
	} else {
		setPartitioningStatus(PARTITIONING_UNSUPPORTED);
	}

	// Initialize partitioning tags
	m_partitioningCellsTag    = -1;
	m_partitioningVerticesTag = -1;

	// Update partitioning information
	if (isPartitioned()) {
		updatePartitioningInfo(true);
	}
#else
	// Initialize serial communicator
	initializeSerialCommunicator();
#endif

	// Initialize the geometrical tolerance
	resetTol();

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
	m_vtk.setGeomData<long>(VTKUnstructuredField::OFFSETS, this);
	m_vtk.setGeomData<int>(VTKUnstructuredField::TYPES, this);
	m_vtk.setGeomData<long>(VTKUnstructuredField::CONNECTIVITY, this);
	m_vtk.setGeomData<long>(VTKUnstructuredField::FACE_STREAMS, this);
	m_vtk.setGeomData<long>(VTKUnstructuredField::FACE_OFFSETS, this);

	// Add VTK basic patch data
	m_vtk.addData<long>("cellIndex", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("PID", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<long>("vertexIndex", VTKFieldType::SCALAR, VTKLocation::POINT, this);
#if BITPIT_ENABLE_MPI==1
	m_vtk.addData<long>("cellGlobalIndex", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("cellRank", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("cellHaloLayer", VTKFieldType::SCALAR, VTKLocation::CELL, this);
	m_vtk.addData<int>("vertexRank", VTKFieldType::SCALAR, VTKLocation::POINT, this);
#endif
}

#if BITPIT_ENABLE_MPI!=1
/*!
	Initialize a dummy communicator to be used when MPI support is disabled.
*/
void PatchKernel::initializeSerialCommunicator()
{
	m_rank        = 0;
	m_nProcessors = 1;
}
#endif

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
		log::cout() << " Error message: " << e.what() << std::endl;
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

	// Early return if the patch is not dirty
	//
	// If we need to squeeze the storage we need to perform the update also
	// if the patch is not dirty.
	if (!squeezeStorage && !isDirty(true)) {
		return updateInfo;
	}

	// Finalize alterations
	finalizeAlterations(squeezeStorage);

	// Spawn
	bool spawnNeeed = (getSpawnStatus() == SPAWN_NEEDED);
	if (spawnNeeed) {
		mergeAdaptionInfo(spawn(trackAdaption), updateInfo);
	}

	// Adaption
	bool adaptionDirty = (getAdaptionStatus(true) == ADAPTION_DIRTY);
	if (adaptionDirty) {
		mergeAdaptionInfo(adaption(trackAdaption, squeezeStorage), updateInfo);
	}

	return updateInfo;
}

/*!
	Simulate the adaption of the specified cell.

	\param id is the id of the cell
	\param marker is the adaption marker of the simulated update
	\param[out] virtualCells are the virtual cells that would be the outcome
	of the update
	\param[out] virtualVertices are the virtual vertices that would be the
	outcome of the update
*/
void PatchKernel::simulateCellUpdate(const long id, adaption::Marker marker, std::vector<Cell> *virtualCells, PiercedVector<Vertex, long> *virtualVertices) const
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(marker);
	BITPIT_UNUSED(virtualCells);
	BITPIT_UNUSED(virtualVertices);

	throw std::runtime_error ("This function has not been implemented for the specified patch.");
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
	if (isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Barrier(communicator);
	}
#endif

	// Check spawn status
	SpawnStatus spawnStatus = getSpawnStatus();
	if (spawnStatus == SPAWN_UNNEEDED || spawnStatus == SPAWN_DONE) {
		return spawnInfo;
	}

	// Spawn the patch
	spawnInfo = _spawn(trackSpawn);

	// Finalize patch alterations
	finalizeAlterations(true);

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
	Alter the patch performing the adaption.

	The actual modification of the patch takes place during this phase. After
	this phase the adaption is completed and the patch is in its final state.
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

	// Adapt the patch
	adaptionInfo = _adaptionAlter(trackAdaption);

	// Finalize patch alterations
	finalizeAlterations(squeezeStorage);

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
	Finalize patch alterations.

	\param squeezeStorage if set to true patch data structures will be
	squeezed
*/
void PatchKernel::finalizeAlterations(bool squeezeStorage)
{
	// Flush vertex data structures
	m_vertices.flush();

	// Update bounding box
	bool boundingBoxDirty = isBoundingBoxDirty();
	if (boundingBoxDirty) {
		updateBoundingBox();
	}

	// Flush cell data structures
	m_cells.flush();

	// Update adjacencies
	bool adjacenciesDirty = areAdjacenciesDirty();
	if (adjacenciesDirty) {
		updateAdjacencies();
	}

	// Update interfaces
	bool interfacesDirty = areInterfacesDirty();
	if (interfacesDirty) {
		updateInterfaces();
	}

	// Flush interfaces data structures
	m_interfaces.flush();

#if BITPIT_ENABLE_MPI==1
	// Update partitioning information
	bool partitioningInfoDirty = arePartitioningInfoDirty();
	if (partitioningInfoDirty) {
		updatePartitioningInfo(true);
	}
#endif

	// Clear alteration flags
	m_alteredCells.clear();
	m_alteredInterfaces.clear();

	// Squeeze the patch
	if (squeezeStorage) {
		squeeze();
	}

	// Synchronize storage
	m_cells.sync();
	m_interfaces.sync();
	m_vertices.sync();
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void PatchKernel::markCellForRefinement(long id)
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
void PatchKernel::markCellForCoarsening(long id)
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
void PatchKernel::resetCellAdaptionMarker(long id)
{
	bool updated = _resetCellAdaptionMarker(id);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
}

/*!
	Returns the adaption marker of the specified cell.

	The marker only defines the type of adaption requested for the cell, it
	is not guaranteed that the adaption will effectively perform the requested
	action (i.e., the requested marker may not be consistent with the internal
	criteria defined by the patch).

	\param id is the id of the cell
	\return The adaption marker of the cell.
*/
adaption::Marker PatchKernel::getCellAdaptionMarker(long id)
{
	return _getCellAdaptionMarker(id);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void PatchKernel::enableCellBalancing(long id, bool enabled)
{
	bool updated = _enableCellBalancing(id, enabled);

	if (updated) {
		setAdaptionStatus(ADAPTION_DIRTY);
	}
}

/*!
	Reset the patch.
*/
void PatchKernel::reset()
{
	resetVertices();
	resetCells();
	resetInterfaces();
}

/*!
	Reset the vertices of the patch.
*/
void PatchKernel::resetVertices()
{
	m_vertices.clear();
	if (m_vertexIdGenerator) {
		m_vertexIdGenerator->reset();
	}
	m_nInternalVertices = 0;
#if BITPIT_ENABLE_MPI==1
	m_nGhostVertices = 0;
#endif
	m_lastInternalVertexId = Vertex::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostVertexId = Vertex::NULL_ID;
#endif

	for (auto &cell : m_cells) {
		cell.unsetConnect();
	}
}

/*!
	Reset the cells of the patch.
*/
void PatchKernel::resetCells()
{
	m_cells.clear();
	if (m_cellIdGenerator) {
		m_cellIdGenerator->reset();
	}
	m_nInternalCells = 0;
#if BITPIT_ENABLE_MPI==1
	m_nGhostCells = 0;
#endif
	m_lastInternalCellId = Cell::NULL_ID;
#if BITPIT_ENABLE_MPI==1
	m_firstGhostCellId = Cell::NULL_ID;
#endif

#if BITPIT_ENABLE_MPI==1
	clearGhostCellsInfo();
	clearGhostVerticesInfo();
#endif

	resetAdjacencies();

	for (auto &interface : m_interfaces) {
		interface.unsetNeigh();
		interface.unsetOwner();
	}

	m_alteredCells.clear();
}

/*!
	Reset the interfaces of the patch.

	This function doesn't change the build strategy, it only resets the
	existing interface.
*/
void PatchKernel::resetInterfaces()
{
	// Early return if no interfaces have been built
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Prune stale interfaces
	pruneStaleInterfaces();

	// Reset the interfaces
	_resetInterfaces(false);

	// All remaining interfaces will be deleted
	setInterfaceAlterationFlags(FLAG_DELETED);

	// Mark cell interfaces as dirty
	setCellAlterationFlags(FLAG_INTERFACES_DIRTY);
}

/*!
	Internal function to reset the interfaces of the patch.

	This function doesn't change the alteration flags.

	\param release if it's true the memory hold by the interfaces will be
	released, otherwise the interfaces will be reset but their memory will
	not be relased
*/
void PatchKernel::_resetInterfaces(bool release)
{
	for (auto &cell : m_cells) {
		cell.resetInterfaces(!release);
	}

	m_interfaces.clear(release);
	if (m_interfaceIdGenerator) {
		m_interfaceIdGenerator->reset();
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
	Writes the patch a filename with the same name of the patch.

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
		vtkCellCount = getInternalCellCount();
#endif
	}

	// Set the dimensions of the mesh
	//
	// Even if there is just a single process that needs to write the VTK face stream, than
	// all the processes need to write the face stream as well (even the processes whose
	// local cells don't require the face steam).
	PiercedStorage<bool, long> vertexWriteFlag(1, &m_vertices);
	vertexWriteFlag.fill(false);

	bool vtkFaceStreamNeeded = false;
	for (const Cell &cell : m_cells) {
		if (cell.getDimension() > 2 && !cell.hasInfo()) {
			vtkFaceStreamNeeded = true;
			break;
		}
	}

#if BITPIT_ENABLE_MPI==1
	if (isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &vtkFaceStreamNeeded, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#endif

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
	VertexConstIterator endItr = vertexConstEnd();
	for (VertexConstIterator itr = vertexConstBegin(); itr != endItr; ++itr) {
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
	// status will always be the same on all the processes.

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
	Checks if the patch supports adaption.

	\return Returns true if the patch supports adaption, false otherwise.
*/
bool PatchKernel::isAdaptionSupported() const
{
    return (getAdaptionStatus() != ADAPTION_UNSUPPORTED);
}

/*!
	Returns the current adaption status.

	\param global if set to true, the adaption status will be evaluated
	globally across all the partitions
	\return The current adaption status.
*/
PatchKernel::AdaptionStatus PatchKernel::getAdaptionStatus(bool global) const
{
	int adaptionStatus = static_cast<int>(m_adaptionStatus);

#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
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
	Returns true if the the patch needs to update its data structures.

	\param global if set to true, the dirty status will be evaluated globally
	across all the partitions
	\return This method returns true to indicate the patch needs to update
	its data structures. Otherwise, it returns false.
*/
bool PatchKernel::isDirty(bool global) const
{
#if BITPIT_ENABLE_MPI==0
	BITPIT_UNUSED(global);
#endif

	bool isDirty = false;

	if (!isDirty) {
		isDirty |= areAdjacenciesDirty(false);
	}

	if (!isDirty) {
		isDirty |= !m_alteredCells.empty();
	}

	if (!isDirty) {
		isDirty |= areInterfacesDirty(false);
		assert(isDirty || m_alteredInterfaces.empty());
	}

	if (!isDirty) {
		isDirty |= (getSpawnStatus() == SPAWN_NEEDED);
	}

	if (!isDirty) {
		isDirty |= (getAdaptionStatus(false) == ADAPTION_DIRTY);
	}

	if (!isDirty) {
		isDirty |= isBoundingBoxDirty(false);
	}

#if BITPIT_ENABLE_MPI==1
	if (!isDirty) {
		isDirty |= arePartitioningInfoDirty(false);
	}
#endif

	if (!isDirty) {
		isDirty |= !m_vertices.isSynced();
	}

	if (!isDirty) {
		isDirty |= !m_cells.isSynced();
	}

	if (!isDirty) {
		isDirty |= !m_interfaces.isSynced();
	}

#if BITPIT_ENABLE_MPI==1
	// Get the global flag
	if (global && isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &isDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#endif

	return isDirty;
}

/*!
	Enables or disables expert mode.

	When expert mode is enabled, it will be possible to change the patch using
	low level functions (e.g., it will be possible to add individual cells,
	add vertices, delete cells, ...).

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

	When expert mode is enabled, it will be possible to change the patch using
	low level functions (e.g., it will be possible to add individual cells,
	add vertices, delete cells, ...).

	\return This method returns true when the expert is enabled,
	otherwise it returns false.
*/
bool PatchKernel::isExpert() const
{
	return m_expert;
}

/*!
	Sets the id of the patch.

	\param id the id of the patch
*/
void PatchKernel::setId(int id)
{
    if (id == m_id) {
        return;
    }

    patch::manager().unregisterPatch(this);
    patch::manager().registerPatch(this, id);
}

/*!
	Internal function to set the id of the patch.

	\param id the id of the patch
*/
void PatchKernel::_setId(int id)
{
	m_id = id;
}

/*!
	Gets the id associated with the patch.

	\return The id associated with the patch.
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

	\return The dimension of the patch.
*/
int PatchKernel::getDimension() const
{
	return m_dimension;
}

/*!
	Returns true if the patch is a three-dimensional patch.

	\return This method returns true to indicate the patch is three-dimensional.
*/
bool PatchKernel::isThreeDimensional() const
{
	return (m_dimension == 3);
}

/*!
	Returns true if auto-indexing is enabled for vertices.

	When auto-indexing is disabled, ids for newly added vertices has to be
	manually specified.

	\return Returns true if auto-indexing is enabled for vertices.
*/
bool PatchKernel::isVertexAutoIndexingEnabled() const
{
	return static_cast<bool>(m_vertexIdGenerator);
}

/*!
	Enables or disables auto-indexing for vertices.

	When auto-indexing is disabled, ids for newly added vertices has to be
	manually specified.

	\param enabled if set to true the auto-indexing will be enabled
*/
void PatchKernel::setVertexAutoIndexing(bool enabled)
{
	if (isVertexAutoIndexingEnabled() == enabled) {
		return;
	}

	if (enabled) {
		createVertexIndexGenerator(true);
	} else {
		m_vertexIdGenerator.reset();
	}
}

/*!
	Dump vertex auto indexing.

 *  \param stream is the stream to write to
*/
void PatchKernel::dumpVertexAutoIndexing(std::ostream &stream) const
{
	m_vertexIdGenerator->dump(stream);
}

/*!
	Restores vertex auto indexing.

	\param stream is the stream to read from
*/
void PatchKernel::restoreVertexAutoIndexing(std::istream &stream)
{
	createVertexIndexGenerator(false);
	m_vertexIdGenerator->restore(stream);
}

/*!
	Create the vertex index generator.

	If the index generator is already created, the existing index generator will be reset.

	\param populate if set to true, the index generator will be populated with current
	vertex ids
*/
void PatchKernel::createVertexIndexGenerator(bool populate)
{
	// Create or reset index generator
	if (!m_vertexIdGenerator) {
		m_vertexIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>());
	} else {
		m_vertexIdGenerator->reset();
	}

	// Populate index generator with current ids
	if (populate) {
		VertexConstIterator beginItr = m_vertices.cbegin();
		VertexConstIterator endItr   = m_vertices.cend();
		for (VertexConstIterator itr = beginItr; itr != endItr; ++itr) {
			m_vertexIdGenerator->setAssigned(itr.getId());
		}
	}
}

/*!
	Import the vertex index generator form the specified source patch.

	If the source patch doesn't define an index generator, the intexe generator of the current
	patch will be deleted.

	\param source is the patch from with the index generator will be imported from
*/
void PatchKernel::importVertexIndexGenerator(const PatchKernel &source)
{
	if (source.m_vertexIdGenerator) {
		m_vertexIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>(*(source.m_vertexIdGenerator)));
	} else {
		m_vertexIdGenerator.reset();
	}
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
	Gets the number of internal vertices in the patch.

	\return The number of internal vertices in the patch.
*/
long PatchKernel::getInternalVertexCount() const
{
	return m_nInternalVertices;
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
Vertex & PatchKernel::getVertex(long id)
{
	return m_vertices[id];
}

/*!
	Gets a constant reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A constant reference to the vertex with the specified id.
*/
const Vertex & PatchKernel::getVertex(long id) const
{
	return m_vertices[id];
}

/*!
	Gets a reference to the last internal vertex.

	\return A reference to the last internal vertex.
*/
Vertex & PatchKernel::getLastInternalVertex()
{
	return m_vertices[m_lastInternalVertexId];
}

/*!
	Gets a constant reference to the last internal vertex.

	\return A constant reference to the last internal vertex.
*/
const Vertex & PatchKernel::getLastInternalVertex() const
{
	return m_vertices[m_lastInternalVertexId];
}

/*!
	Returns an iterator pointing to the specified vertex.

	\result An iterator to the specified vertex.
*/
PatchKernel::VertexIterator PatchKernel::getVertexIterator(long id)
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
	Returns iterator pointing to the first internal vertex.

	\result An iterator to the first internal vertex.
*/
PatchKernel::VertexIterator PatchKernel::internalVertexBegin()
{
	return m_vertices.begin();
}

/*!
	Returns iterator pointing to the end of the list of internal vertices.

	\result An iterator to the end of the list of internal vertices.
*/
PatchKernel::VertexIterator PatchKernel::internalVertexEnd()
{
	if (m_nInternalVertices > 0) {
		return ++m_vertices.find(m_lastInternalVertexId);
	} else {
		return m_vertices.end();
	}
}

/*!
	Returns a constant iterator pointing to the specified vertex.

	\result A constant iterator to the specified vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::getVertexConstIterator(long id) const
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
	Returns a constant iterator pointing to the first internal vertex.

	\result A constant iterator to the first internal  vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::internalVertexConstBegin() const
{
	return m_vertices.cbegin();
}

/*!
	Returns a constant iterator pointing to last internal vertex.

	\result A constant iterator to the last internal vertex.
*/
PatchKernel::VertexConstIterator PatchKernel::internalVertexConstEnd() const
{
	if (m_nInternalVertices > 0) {
		return ++m_vertices.find(m_lastInternalVertexId);
	} else {
		return m_vertices.end();
	}
}

/*!
	Adds the specified vertex to the patch.

	All new vertices will be temporarily added as internal vertices, is needed
	they will be converted to ghost vertices when updating ghost information.

	If valid, the specified id will we assigned to the newly created vertex,
	otherwise a new unique id will be generated for the vertex. However, it
	is not possible to create a new vertex with an id already assigned to an
	existing vertex of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the vertex that will be added
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, a new unique id will be generated
	for the vertex
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(const Vertex &source, long id)
{
	Vertex vertex = source;
	vertex.setId(id);

	return addVertex(std::move(vertex), id);
}

/*!
	Adds the specified vertex to the patch.

	All new vertices will be temporarily added as internal vertices, is needed
	they will be converted to ghost vertices when updating ghost information.

	If valid, the specified id will we assigned to the newly created vertex,
	otherwise a new unique id will be generated for the vertex. However, it
	is not possible to create a new vertex with an id already assigned to an
	existing vertex of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the vertex that will be added
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added vertex.
*/
PatchKernel::VertexIterator PatchKernel::addVertex(Vertex &&source, long id)
{
	// Get the id
	if (id < 0) {
		id = source.getId();
	}

	// Add a dummy vertex
	//
	// It is not possible to directly add the source into the storage. First a
	// dummy vertex is created and then that vertex is replaced with the source.
	//
	// The corrdinates of the dummy vertex shuold match the coordinates of the
	// source, because the bounding box of the patch is updated every time a
	// new vertex is added.
	VertexIterator iterator = _addInternalVertex(source.getCoords(), id);

	// Replace the newly created vertex with the source
	//
	// Before replacing the newly created vertex we need to set the id
	// of the source to the id that has been assigned to the newly created
	// vertex, we also need to mark the vertex as internal.
	source.setId(iterator->getId());
	source.setInterior(true);

	Vertex &vertex = (*iterator);
	vertex = std::move(source);

	return iterator;
}

/*!
	Adds a new vertex with the specified coordinates.

	All new vertices will be temporarily added as internal vertices, is needed
	they will be converted to ghost vertices when updating ghost information.

	If valid, the specified id will we assigned to the newly created vertex,
	otherwise a new unique id will be generated for the vertex. However, it
	is not possible to create a new vertex with an id already assigned to an
	existing vertex of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	// Add the vertex
	VertexIterator iterator = _addInternalVertex(coords, id);

	return iterator;
}

/*!
	Internal function to add an internal vertex.

	It is not possible to create a new vertex with an id already assigned to an
	existing vertex of the patch or with an invalid id. If this happens, an
	exception is thrown. Ids are considered valid if they are greater or equal
	than zero.

	\param coords are the coordinates of the vertex
	\param id is the id that will be assigned to the newly created vertex.
	If a negative id value is specified, a new unique id will be generated
	for the vertex
	\return An iterator pointing to the newly created vertex.
*/
PatchKernel::VertexIterator PatchKernel::_addInternalVertex(const std::array<double, 3> &coords, long id)
{
	// Get the id
	if (m_vertexIdGenerator) {
		if (id < 0) {
			id = m_vertexIdGenerator->generate();
		} else {
			m_vertexIdGenerator->setAssigned(id);
		}
	} else if (id < 0) {
		throw std::runtime_error("No valid id has been provided for the vertex.");
	}

	// Get the id of the vertex before which the new vertex should be inserted
#if BITPIT_ENABLE_MPI==1
	//
	// If there are ghosts vertices, the internal vertex should be inserted
	// before the first ghost vertex.
#endif
	long referenceId;
#if BITPIT_ENABLE_MPI==1
	referenceId = m_firstGhostVertexId;
#else
	referenceId = Vertex::NULL_ID;
#endif

	// Create the vertex
	VertexIterator iterator;
	if (referenceId == Vertex::NULL_ID) {
		iterator = m_vertices.emreclaim(id, id, coords, true);
	} else {
		iterator = m_vertices.emreclaimBefore(referenceId, id, id, coords, true);
	}
	m_nInternalVertices++;

	// Update the id of the last internal vertex
	if (m_lastInternalVertexId < 0) {
		m_lastInternalVertexId = id;
	} else if (m_vertices.rawIndex(m_lastInternalVertexId) < m_vertices.rawIndex(id)) {
		m_lastInternalVertexId = id;
	}

	// Update the bounding box
	addPointToBoundingBox(iterator->getCoords());

#if BITPIT_ENABLE_MPI==1
	// Set partitioning information as dirty
	setPartitioningInfoDirty(true);
#endif

	return iterator;
}

#if BITPIT_ENABLE_MPI==0
/*!
	Restore the vertex with the specified id.

	The kernel should already contain the vertex, only the contents of the
	vertex will be updated.

	\param coords are the coordinates of the vertex
	\param id is the id of the vertex to restore
	\return An iterator pointing to the restored vertex.
*/
PatchKernel::VertexIterator PatchKernel::restoreVertex(const std::array<double, 3> &coords, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = m_vertices.find(id);
	if (iterator == m_vertices.end()) {
		throw std::runtime_error("Unable to restore the specified vertex: the kernel doesn't contain an entry for that vertex.");
	}

	_restoreInternalVertex(iterator, coords);

	return iterator;
}
#endif

/*!
	Internal function to restore an internal vertex.

	\param iterator is an iterator pointing to the vertex to restore
	\param coords are the coordinates of the vertex
*/
void PatchKernel::_restoreInternalVertex(const VertexIterator &iterator, const std::array<double, 3> &coords)
{
	// Restore the vertex
	//
	// There is no need to set the id of the vertex as assigned, because
	// also the index generator will be restored.
	Vertex &vertex = *iterator;
	vertex.initialize(iterator.getId(), coords, true);
	m_nInternalVertices++;

	// Update the bounding box
	addPointToBoundingBox(vertex.getCoords());
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
*/
bool PatchKernel::deleteVertex(long id)
{
	if (!isExpert()) {
		return false;
	}

	// Delete the vertex
#if BITPIT_ENABLE_MPI==1
	const Vertex &vertex = m_vertices[id];
	bool isInternal = vertex.isInterior();
	if (isInternal) {
		_deleteInternalVertex(id);
	} else {
		_deleteGhostVertex(id);
	}
#else
	_deleteInternalVertex(id);
#endif

	return true;
}

/*!
	Internal function to delete an internal vertex.

	\param id is the id of the vertex
*/
void PatchKernel::_deleteInternalVertex(long id)
{
	// Update the bounding box
	const Vertex &vertex = m_vertices[id];
	removePointFromBoundingBox(vertex.getCoords());

	// Delete vertex
	m_vertices.erase(id, true);
	m_nInternalVertices--;
	if (id == m_lastInternalVertexId) {
		updateLastInternalVertexId();
	}

	// Vertex id is no longer used
	if (m_vertexIdGenerator) {
		m_vertexIdGenerator->trash(id);
	}
}

/*!
	Counts free vertices within the patch.

	A free vertex is a vertex on a free face.

	\return The number of free vertices.
*/
long PatchKernel::countFreeVertices() const
{
	return countBorderVertices();
}

/*!
	Counts border vertices within the patch.

	A border vertex is a vertex on a border face.

	\return The number of border vertices.
*/
long PatchKernel::countBorderVertices() const
{
	std::unordered_set<long> borderVertices;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				continue;
			}

			int nFaceVertices = cell.getFaceVertexCount(i);
			for (int k = 0; k < nFaceVertices; ++k) {
				long faceVertexId = cell.getFaceVertexId(i, k);
				borderVertices.insert(faceVertexId);
			}
		}
	}

	return borderVertices.size();
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
	VertexConstIterator beginItr = m_vertices.cbegin();
	VertexConstIterator endItr   = m_vertices.cend();

	// Detect used vertices
	PiercedStorage<bool, long> vertexUsedFlag(1, &m_vertices);
	vertexUsedFlag.fill(false);
	for (const Cell &cell : m_cells) {
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		int nCellVertices = cellVertexIds.size();
		for (int j = 0; j < nCellVertices; ++j) {
			long vertexId = cellVertexIds[j];
			vertexUsedFlag[vertexId] = true;
		}
	}

	// Count the orphan vertices
	std::size_t nOrhpanVertices = 0;
	for (VertexConstIterator itr = beginItr; itr != endItr; ++itr) {
		std::size_t vertexRawIndex = itr.getRawIndex();
		if (!vertexUsedFlag.rawAt(vertexRawIndex)) {
			++nOrhpanVertices;
		}
	}

	// Build a list of unused vertices
	std::vector<long> orhpanVertices(nOrhpanVertices);

	std::size_t orphanVertexIndex = 0;
	for (VertexConstIterator itr = beginItr; itr != endItr; ++itr) {
		std::size_t vertexRawIndex = itr.getRawIndex();
		if (!vertexUsedFlag.rawAt(vertexRawIndex)) {
			orhpanVertices[orphanVertexIndex] = itr.getId();
			++orphanVertexIndex;
		}
	}

	return orhpanVertices;
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
	updateBoundingBox();

	return true;
}

/*!
	Find and collapse coincident vertices. Cell connectivity is
	automatically updated.

	\result The list of the of the collapsed vertices.
*/
std::vector<long> PatchKernel::collapseCoincidentVertices()
{
	std::vector<long> collapsedVertices;
	if (!isExpert()) {
		return collapsedVertices;
	}

	long nVertices = getVertexCount();
	if (nVertices == 0) {
		return collapsedVertices;
	}

	// Update bounding box
	updateBoundingBox();

	//
	// Collapse double vertices
	//
	Vertex::Less vertexLess(10 * std::numeric_limits<double>::epsilon());
	auto rawVertexLess = [this, &vertexLess](const std::size_t &i, const std::size_t &j)
	{
		return vertexLess(this->m_vertices.rawAt(i), this->m_vertices.rawAt(j));
	};
	std::set<std::size_t, decltype(rawVertexLess)> vertexTree(rawVertexLess);

	std::unordered_map<long, long> vertexMap;
	for (VertexConstIterator vertexItr = m_vertices.cbegin(); vertexItr != m_vertices.cend(); ++vertexItr) {
		std::size_t vertexRawId = vertexItr.getRawIndex();
		auto vertexTreeItr = vertexTree.find(vertexRawId);
		if (vertexTreeItr == vertexTree.end()) {
			vertexTree.insert(vertexRawId);
		} else {
			long vertexId = vertexItr.getId();
			long vertexCoincidentId = m_vertices.rawFind(*vertexTreeItr).getId();
			vertexMap.insert({vertexId, vertexCoincidentId});
		}
	}

	// Collapse vertices
	if (!vertexMap.empty()) {
		// Update cells
		//
		// If the adjacencies are currently up-to-date, they should kept up-to-date.
		AdjacenciesBuildStrategy adjacenciesBuildStrategy = getAdjacenciesBuildStrategy();

		bool keepAdjacenciesUpToDate;
		if (adjacenciesBuildStrategy != ADJACENCIES_NONE) {
			keepAdjacenciesUpToDate = (!areAdjacenciesDirty());
		} else {
			keepAdjacenciesUpToDate = false;
		}

		for (Cell &cell : m_cells) {
			// Renumber cell vertices
			int nRenumberedVertices = cell.renumberVertices(vertexMap);

			// Mark adjacencies are dirty
			//
			// If some vertices have been renumbered, the adjacencies of the cells
			// are now dirty.
			if ((adjacenciesBuildStrategy != ADJACENCIES_NONE) && (nRenumberedVertices > 0)) {
				setCellAlterationFlags(cell.getId(), FLAG_ADJACENCIES_DIRTY);
			}
		}

		if (keepAdjacenciesUpToDate) {
			updateAdjacencies();
		}

		// Update interface
		for (Interface &interface : m_interfaces) {
			// Renumber interface vertices
			interface.renumberVertices(vertexMap);
		}

		// Create the list of collapsed vertices
		collapsedVertices.resize(vertexMap.size());

		std::size_t k = 0;
		for (const auto &entry : vertexMap) {
			collapsedVertices[k++] = entry.first;
		}
	}

	return collapsedVertices;
}

/*!
	Remove coincident vertices from the patch.
*/
bool PatchKernel::deleteCoincidentVertices()
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> verticesToDelete = collapseCoincidentVertices();
	deleteVertices(verticesToDelete);

	return true;
}

/*!
	Gets the coordinates of the specified vertex.

	\param id is the id of the vertex
	\result The coordinates of the specified vertex.
*/
const std::array<double, 3> & PatchKernel::getVertexCoords(long id) const
{
	return getVertex(id).getCoords();
}

/*!
	Gets the coordinates of the specified vertices.

	\param nVertices is the number of vertices
	\param ids are the ids of the requested vertices
	\param[out] coordinates on output will contain the vertex coordinates
*/
void PatchKernel::getVertexCoords(std::size_t nVertices, const long *ids, std::unique_ptr<std::array<double, 3>[]> *coordinates) const
{
	*coordinates = std::unique_ptr<std::array<double, 3>[]>(new std::array<double, 3>[nVertices]);
	getVertexCoords(nVertices, ids, coordinates->get());
}

/*!
	Gets the coordinates of the specified vertices.

	\param nVertices is the number of vertices
	\param ids are the ids of the requested vertices
	\param[out] coordinates on output will contain the vertex coordinates, it
	is up to the caller to ensure that the storage has enough space for all
	the vertex coordinates
*/
void PatchKernel::getVertexCoords(std::size_t nVertices, const long *ids, std::array<double, 3> *coordinates) const
{
	for (std::size_t i = 0; i < nVertices; ++i) {
		coordinates[i] = getVertex(ids[i]).getCoords();
	}
}

/*!
	Return true if the patch is empty.

	\param global if set to true, the empty status will be evaluated globally
	across all the partitions
	\return Return true if the patch is empty.
*/
bool PatchKernel::empty(bool global) const
{
	bool isEmpty = (getCellCount() == 0);
#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
		MPI_Allreduce(MPI_IN_PLACE, &isEmpty, 1, MPI_C_BOOL, MPI_LAND, getCommunicator());
	}
#else
	BITPIT_UNUSED(global);
#endif

	return isEmpty;
}

/*!
	Updates the id of the last internal vertex.
*/
void PatchKernel::updateLastInternalVertexId()
{
	if (m_nInternalVertices == 0) {
		m_lastInternalVertexId = Vertex::NULL_ID;
		return;
	}

	VertexIterator lastInternalVertexItr;
#if BITPIT_ENABLE_MPI==1
	if (m_nGhostVertices == 0) {
		lastInternalVertexItr = --m_vertices.end();
		m_lastInternalVertexId = lastInternalVertexItr->getId();
	} else {
		m_lastInternalVertexId = m_vertices.getSizeMarker(m_nInternalVertices - 1, Vertex::NULL_ID);
	}
#else
	lastInternalVertexItr = --m_vertices.end();
	m_lastInternalVertexId = lastInternalVertexItr->getId();
#endif
}

/*!
	Returns true if auto-indexing is enabled for cells.

	When auto-indexing is disabled, cell ids for newly added cells should
	be provided by the user.

	\return Returns true if auto-indexing is enabled for cells.
*/
bool PatchKernel::isCellAutoIndexingEnabled() const
{
	return static_cast<bool>(m_cellIdGenerator);
}

/*!
	Enables or disables auto-indexing for cells.

	When auto-indexing is disabled, ids for newly added cells has to be
	manually specified.

	\param enabled if set to true the auto-indexing will be enabled
*/
void PatchKernel::setCellAutoIndexing(bool enabled)
{
	if (isCellAutoIndexingEnabled() == enabled) {
		return;
	}

	if (enabled) {
		createCellIndexGenerator(true);
	} else {
		m_cellIdGenerator.reset();
	}
}

/*!
	Dump cell auto indexing.

 *  \param stream is the stream to write to
*/
void PatchKernel::dumpCellAutoIndexing(std::ostream &stream) const
{
	m_cellIdGenerator->dump(stream);
}

/*!
	Restores cell auto indexing.

	\param stream is the stream to read from
*/
void PatchKernel::restoreCellAutoIndexing(std::istream &stream)
{
	createCellIndexGenerator(false);
	m_cellIdGenerator->restore(stream);
}

/*!
	Create the cell index generator.

	If the index generator is already created, the existing index generator will be reset.

	\param populate if set to true, the index generator will be populated with current
	cell ids
*/
void PatchKernel::createCellIndexGenerator(bool populate)
{
	// Create or reset index generator
	if (!m_cellIdGenerator) {
		m_cellIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>());
	} else {
		m_cellIdGenerator->reset();
	}

	// Populate index generator with current ids
	if (populate) {
		CellConstIterator beginItr = m_cells.cbegin();
		CellConstIterator endItr   = m_cells.cend();
		for (CellConstIterator itr = beginItr; itr != endItr; ++itr) {
			m_cellIdGenerator->setAssigned(itr.getId());
		}
	}
}

/*!
	Import the cell index generator form the specified source patch.

	If the source patch doesn't define an index generator, the intexe generator of the current
	patch will be deleted.

	\param source is the patch from with the index generator will be imported from
*/
void PatchKernel::importCellIndexGenerator(const PatchKernel &source)
{
	if (source.m_cellIdGenerator) {
		m_cellIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>(*(source.m_cellIdGenerator)));
	} else {
		m_cellIdGenerator.reset();
	}
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch.
*/
long PatchKernel::getCellCount() const
{
	return m_cells.size();
}

/*!
	Gets the number of internal cells in the patch.

	\return The number of internal cells in the patch.
*/
long PatchKernel::getInternalCellCount() const
{
	return m_nInternalCells;
}

/*!
	Gets the number of internal cells in the patch.

	\return The number of internal cells in the patch.
*/
long PatchKernel::getInternalCount() const
{
	return getInternalCellCount();
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
Cell & PatchKernel::getCell(long id)
{
	return m_cells[id];
}

/*!
	Gets a constant reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A constant reference to the cell with the specified id.
*/
const Cell & PatchKernel::getCell(long id) const
{
	return m_cells[id];
}

/*!
	Gets the element type for the cell with the specified id.

	\param id is the id of the requested cell
	\return The element type for the cell with the specified id.
*/
ElementType PatchKernel::getCellType(long id) const
{
	return m_cells[id].getType();
}

/*!
	Gets a reference to the last internal cell.

	\return A reference to the last internal cell.
*/
Cell & PatchKernel::getLastInternalCell()
{
	return m_cells[m_lastInternalCellId];
}

/*!
	Gets a constant reference to the last internal cell.

	\return A constant reference to the last internal cell.
*/
const Cell & PatchKernel::getLastInternalCell() const
{
	return m_cells[m_lastInternalCellId];
}

/*!
	Returns an iterator pointing to the specified cell.

	\result An iterator to the specified cell.
*/
PatchKernel::CellIterator PatchKernel::getCellIterator(long id)
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
PatchKernel::CellIterator PatchKernel::internalCellBegin()
{
	return m_cells.begin();
}

/*!
	Returns iterator pointing to the first internal cell.

	\result An iterator to the first internal cell.
*/
PatchKernel::CellIterator PatchKernel::internalBegin()
{
	return internalCellBegin();
}

/*!
	Returns iterator pointing to the end of the list of internal cells.

	\result An iterator to the end of the list of internal cells.
*/
PatchKernel::CellIterator PatchKernel::internalCellEnd()
{
	if (m_nInternalCells > 0) {
		return ++m_cells.find(m_lastInternalCellId);
	} else {
		return m_cells.end();
	}
}

/*!
	Returns iterator pointing to the end of the list of internal cells.

	\result An iterator to the end of the list of internal cells.
*/
PatchKernel::CellIterator PatchKernel::internalEnd()
{
	return internalCellEnd();
}

/*!
	Returns a constant iterator pointing to the specified cell.

	\result A constant iterator to the specified cell.
*/
PatchKernel::CellConstIterator PatchKernel::getCellConstIterator(long id) const
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
PatchKernel::CellConstIterator PatchKernel::internalCellConstBegin() const
{
	return m_cells.cbegin();
}

/*!
	Returns a constant iterator pointing to the first internal cell.

	\result A constant iterator to the first internal cell.
*/
PatchKernel::CellConstIterator PatchKernel::internalConstBegin() const
{
	return internalCellConstBegin();
}

/*!
	Returns a constant iterator pointing to the end of the list of internal
	cells.

	\result A constant iterator to the end of the list of internal cells.
*/
PatchKernel::CellConstIterator PatchKernel::internalCellConstEnd() const
{
	if (m_nInternalCells > 0) {
		return ++m_cells.find(m_lastInternalCellId);
	} else {
		return m_cells.cend();
	}
}

/*!
	Returns a constant iterator pointing to the end of the list of internal
	cells.

	\result A constant iterator to the end of the list of internal cells.
*/
PatchKernel::CellConstIterator PatchKernel::internalConstEnd() const
{
	return internalCellConstEnd();
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, long id)
{
	Cell cell = source;
	cell.setId(id);

	return addCell(std::move(cell), id);
}

/*!
	Adds the specified cell to the patch.

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, long id)
{
	// Get the id
	if (id < 0) {
		id = source.getId();
	}

	// Add a dummy cell
	//
	// It is not possible to directly add the source into the storage. First a
	// dummy cell is created and then that cell is replaced with the source.
	std::unique_ptr<long[]> dummyConnectStorage = std::unique_ptr<long[]>(nullptr);

	CellIterator iterator = _addInternalCell(ElementType::UNDEFINED, std::move(dummyConnectStorage), id);

	// Replace the newly created cell with the source
	//
	// Before replacing the newly created cell we need to set the id of the
	// source to the id that has been assigned to the newly created cell.
	source.setId(iterator->getId());

	Cell &cell = (*iterator);
	cell = std::move(source);

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

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	If valid, the specified id will we assigned to the newly created cell,
	otherwise a new unique id will be generated for the cell. However, it
	is not possible to create a new cell with an id already assigned to an
	existing cell of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	CellIterator iterator = _addInternalCell(type, std::move(connectStorage), id);

	return iterator;
}

/*!
	Internal function to add an internal cell.

	It is not possible to create a new cell with an id already assigned to an
	existing cell of the patch or with an invalid id. If this happens, an
	exception is thrown. Ids are considered valid if they are greater or equal
	than zero.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::_addInternalCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
													long id)
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

	// Get the id of the cell before which the new cell should be inserted
#if BITPIT_ENABLE_MPI==1
	//
	// If there are ghosts cells, the internal cell should be inserted
	// before the first ghost cell.
#endif
	long referenceId;
#if BITPIT_ENABLE_MPI==1
	referenceId = m_firstGhostCellId;
#else
	referenceId = Cell::NULL_ID;
#endif

	// Create the cell
	bool storeInterfaces  = (getInterfacesBuildStrategy() != INTERFACES_NONE);
	bool storeAdjacencies = storeInterfaces || (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	CellIterator iterator;
	if (referenceId == Cell::NULL_ID) {
		iterator = m_cells.emreclaim(id, id, type, std::move(connectStorage), true, storeInterfaces, storeAdjacencies);
	} else {
		iterator = m_cells.emreclaimBefore(referenceId, id, id, type, std::move(connectStorage), true, storeInterfaces, storeAdjacencies);
	}
	m_nInternalCells++;

	// Update the id of the last internal cell
	if (m_lastInternalCellId < 0) {
		m_lastInternalCellId = id;
	} else if (m_cells.rawIndex(m_lastInternalCellId) < m_cells.rawIndex(id)) {
		m_lastInternalCellId = id;
	}

	// Set the alteration flags of the cell
	setAddedCellAlterationFlags(id);

#if BITPIT_ENABLE_MPI==1
	// Set partitioning information as dirty
	setPartitioningInfoDirty(true);
#endif

	return iterator;
}

/*!
	Set the alteration flags for an added cell.

	Only alteration flags needed for tracking the status of the patch are
	added (for example, there is no explicit flag to tag newly added cells).

	\param id is the id of the cell
*/
void PatchKernel::setAddedCellAlterationFlags(long id)
{
	AlterationFlags flags = FLAG_NONE;
	if (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE) {
		flags |= FLAG_ADJACENCIES_DIRTY;
	}
	if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
		flags |= FLAG_INTERFACES_DIRTY;
	}

	setCellAlterationFlags(id, flags);
}

#if BITPIT_ENABLE_MPI==0
/*!
	Restore the cell with the specified id.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id of the cell that will be restored
	\return An iterator pointing to the restored cell.
*/
PatchKernel::CellIterator PatchKernel::restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												   long id)
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

	_restoreInternalCell(iterator, type, std::move(connectStorage));

	return iterator;
}
#endif

/*!
	Internal function to restore an internal cell.

	\param iterator is an iterator pointing to the cell to restore
	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
*/
void PatchKernel::_restoreInternalCell(const CellIterator &iterator, ElementType type,
								   std::unique_ptr<long[]> &&connectStorage)
{
	// Restore the cell
	//
	// There is no need to set the id of the cell as assigned, because
	// also the index generator will be restored.
	bool storeInterfaces  = (getInterfacesBuildStrategy() != INTERFACES_NONE);
	bool storeAdjacencies = storeInterfaces || (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	long cellId = iterator.getId();
	Cell &cell = *iterator;
	cell.initialize(cellId, type, std::move(connectStorage), true, storeInterfaces, storeAdjacencies);
	m_nInternalCells++;

	// Set the alteration flags of the cell
	setRestoredCellAlterationFlags(cellId);
}

/*!
	Set the alteration flags for a restored cell.

	\param id is the id of the cell
*/
void PatchKernel::setRestoredCellAlterationFlags(long id)
{
	AlterationFlags flags = FLAG_NONE;
	if (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE) {
		flags |= FLAG_ADJACENCIES_DIRTY;
	}
	if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
		flags |= FLAG_INTERFACES_DIRTY;
	}

	setCellAlterationFlags(id, flags);
}

/*!
	Deletes a cell.

	\param id is the id of the cell
*/
bool PatchKernel::deleteCell(long id)
{
	if (!isExpert()) {
		return false;
	}

	// Delete cell
#if BITPIT_ENABLE_MPI==1
	const Cell &cell = m_cells[id];
	bool isInternalCell = cell.isInterior();
	if (isInternalCell) {
		_deleteInternalCell(id);
	} else {
		_deleteGhostCell(id);
	}
#else
	_deleteInternalCell(id);
#endif

	return true;
}

/*!
	Internal function to delete an internal cell.

	\param id is the id of the cell
*/
void PatchKernel::_deleteInternalCell(long id)
{
	// Set the alteration flags of the cell
	setDeletedCellAlterationFlags(id);

	// Delete cell
	m_cells.erase(id, true);
	m_nInternalCells--;
	if (id == m_lastInternalCellId) {
		updateLastInternalCellId();
	}

	// Cell id is no longer used
	if (m_cellIdGenerator) {
		m_cellIdGenerator->trash(id);
	}
}

/*!
	Set the alteration flags for a deleted cell.

	Only alteration flags needed for tracking the status of the patch are set.

	\param id is the id of the cell
*/
void PatchKernel::setDeletedCellAlterationFlags(long id)
{
	const Cell &cell = getCell(id);

	// Set the alteration flags of the cell
	resetCellAlterationFlags(id, FLAG_DELETED);

	// Set the alteration flags of the adjacent cells
	const int nCellAdjacencies = cell.getAdjacencyCount();
	const long *cellAdjacencies = cell.getAdjacencies();
	for (int k = 0; k < nCellAdjacencies; ++k) {
		long adjacencyId = cellAdjacencies[k];
		if (!testCellAlterationFlags(adjacencyId, FLAG_DELETED)) {
			AlterationFlags flags = FLAG_DANGLING;
			if (getAdjacenciesBuildStrategy() != ADJACENCIES_NONE) {
				flags |= FLAG_ADJACENCIES_DIRTY;
			}
			if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
				flags |= FLAG_INTERFACES_DIRTY;
			}

			setCellAlterationFlags(adjacencyId, flags);
		}
	}

	// Set the alteration flags of the interfaces
	const int nCellInterfaces = cell.getInterfaceCount();
	const long *cellInterfaces = cell.getInterfaces();
	for (int k = 0; k < nCellInterfaces; ++k) {
		long interfaceId = cellInterfaces[k];
		if (!testInterfaceAlterationFlags(interfaceId, FLAG_DELETED)) {
			setInterfaceAlterationFlags(interfaceId, FLAG_DANGLING);
		}
	}
}

/*!
	Counts free cells within the patch.

	A cell is free if contains at least one free face.

	\return The number of free cells.
*/
long PatchKernel::countFreeCells() const
{
	return countBorderCells();
}

/*!
	Counts border cells within the patch.

	A cell is border if contains at least one border face.

	\return The number of border cells.
*/
long PatchKernel::countBorderCells() const
{
	double nBorderCells = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nBorderCells;
				break;
			}
		}
	}

	return nBorderCells;
}

/*!
	Counts orphan cells within the patch.

	A cell is orphan if not adjacent to any cell in the patch (neither
	along an edge, nor at vertex)

	\return The number of orphan cells.
*/
long PatchKernel::countOrphanCells() const
{
	return findOrphanCells().size();
}

/*!
	Finds orphan cells within the patch.

	A cell is orphan if not adjacent to any cell in the patch (neither
	along an edge, nor at vertex)

	\return The number of orphan cells.
*/
std::vector<long> PatchKernel::findOrphanCells() const
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
	std::vector<long> orphanCells;
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
			orphanCells.push_back(cell.getId());
		}
	}

	return orphanCells;
}

/*!
	Finds duplicate cells within the patch.

	A cell is a duplicate if there is at least one other cell with exactly
	the same vertices.

	\return The number of duplicate cells.
*/
long PatchKernel::countDuplicateCells() const
{
	return findDuplicateCells().size();
}

/*!
	Finds duplicate cells within the patch.

	A cell is a duplicate if there is at least one other cell with exactly
	the same vertices.

	\return The list of duplicate cells.
*/
std::vector<long> PatchKernel::findDuplicateCells() const
{
	// Define the hasher and the predicate to be used for lists of vertex ids
	//
	// These operators should be able to identify that two lists of vertex ids
	// are the same regardless of the order in which the ids are listed.
	auto hasher = [](const ConstProxyVector<long> &ids)
	{
		std::size_t hash = std::hash<long>{}(0);
		for (long id : ids) {
			hash = hash + std::hash<long>{}(id);
		}

		return hash;
	};

	std::unordered_map<long, std::size_t > counters;
	auto predicate = [&counters](const ConstProxyVector<long> &ids_A, const ConstProxyVector<long> &ids_B)
	{
		counters.clear();
		for (long id : ids_A) {
			counters[id]++;
		}
		for (long id : ids_B) {
			counters[id]--;
		}

		for (const auto &entry : counters) {
			if (entry.second != 0) {
				return false;
			}
		}

		return true;
	};

	// Detect if there are cells that share the same list of vertices
	//
	// For each cell, the list of vertex ids is added into a set. If a collision
	// is detected, we have found a duplicate cell.
	std::vector<long> duplicateCells;
	std::unordered_set<ConstProxyVector<long>, decltype(hasher), decltype(predicate)> vertexStash(getCellCount(), hasher, predicate);
	for (const Cell &cell : m_cells) {
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		auto vertexStashItr = vertexStash.find(cellVertexIds);
		if (vertexStashItr == vertexStash.end()) {
			vertexStash.insert(std::move(cellVertexIds));
		} else {
			duplicateCells.push_back(cell.getId());
		}
	}

	return duplicateCells;
}

/*!
	Updates the id of the last internal cell.
*/
void PatchKernel::updateLastInternalCellId()
{
	if (m_nInternalCells == 0) {
		m_lastInternalCellId = Cell::NULL_ID;
		return;
	}

	CellIterator lastInternalCellItr;
#if BITPIT_ENABLE_MPI==1
	if (m_nGhostCells == 0) {
		lastInternalCellItr = --m_cells.end();
		m_lastInternalCellId = lastInternalCellItr->getId();
	} else {
		m_lastInternalCellId = m_cells.getSizeMarker(m_nInternalCells - 1, Cell::NULL_ID);
	}
#else
	lastInternalCellItr = --m_cells.end();
	m_lastInternalCellId = lastInternalCellItr->getId();
#endif
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
std::vector<long> PatchKernel::findCellNeighs(long id) const
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
void PatchKernel::findCellNeighs(long id, std::vector<long> *neighs) const
{
	_findCellNeighs(id, nullptr, neighs);
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
std::vector<long> PatchKernel::findCellNeighs(long id, int codimension, bool complete) const
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
void PatchKernel::findCellNeighs(long id, int codimension, bool complete, std::vector<long> *neighs) const
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
std::vector<long> PatchKernel::findCellFaceNeighs(long id) const
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
	The blacklist has to be a pointer to a unique list of ordered cell ids
	or a null pointer if no cells should be excluded from the search
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellNeighs(long id, const std::vector<long> *blackList, std::vector<long> *neighs) const
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
void PatchKernel::findCellFaceNeighs(long id, std::vector<long> *neighs) const
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
		_findCellFaceNeighs(id, i, nullptr, neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> PatchKernel::findCellFaceNeighs(long id, int face) const
{
	std::vector<long> neighs;
	_findCellFaceNeighs(id, face, nullptr, &neighs);

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
void PatchKernel::findCellFaceNeighs(long id, int face, std::vector<long> *neighs) const
{
	_findCellFaceNeighs(id, face, nullptr, neighs);
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a pointer to a unique list of ordered cell ids
	or a null pointer if no cells should be excluded from the search
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given face will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellFaceNeighs(long id, int face, const std::vector<long> *blackList, std::vector<long> *neighs) const
{
	const Cell &cell = getCell(id);

	int nFaceAdjacencies = cell.getAdjacencyCount(face);
	const long *faceAdjacencies = cell.getAdjacencies(face);
	for (int k = 0; k < nFaceAdjacencies; ++k) {
		long neighId = faceAdjacencies[k];
		if (!blackList || utils::findInOrderedVector<long>(neighId, *blackList) == blackList->end()) {
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
std::vector<long> PatchKernel::findCellEdgeNeighs(long id, bool complete) const
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
void PatchKernel::findCellEdgeNeighs(long id, bool complete, std::vector<long> *neighs) const
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
	std::unique_ptr<std::vector<long>> blackList;
	if (!complete) {
		blackList = std::unique_ptr<std::vector<long>>(new std::vector<long>());
		findCellFaceNeighs(id, blackList.get());
	}

	for (int i = 0; i < nCellEdges; ++i) {
		_findCellEdgeNeighs(id, i, blackList.get(), neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> PatchKernel::findCellEdgeNeighs(long id, int edge) const
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
void PatchKernel::findCellEdgeNeighs(long id, int edge, std::vector<long> *neighs) const
{
	_findCellEdgeNeighs(id, edge, nullptr, neighs);
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a pointer to a unique list of ordered cell ids
	or a null pointer if no cells should be excluded from the search
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given edge will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellEdgeNeighs(long id, int edge, const std::vector<long> *blackList, std::vector<long> *neighs) const
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

	// Find candidates neighbours
	//
	// These are the negihbours of the first vertex of the edge.
	const int GUESS_NEIGHS_COUNT = 3;

	std::vector<long> candidateIds;
	candidateIds.reserve(GUESS_NEIGHS_COUNT);
	_findCellVertexNeighs(id, edgeVertices[0], blackList, &candidateIds);

	// Discard candidates that doesn't contain the last vertex of the edge
	long lastEdgeVertexId = cell.getEdgeVertexId(edge, nEdgeVertices - 1);

	ConstProxyVector<long> candidateVertexIds;
	for (long candidateId : candidateIds) {
		const Cell &candidateNeigh = m_cells.at(candidateId);
		candidateVertexIds = candidateNeigh.getVertexIds();
		std::size_t nCandidateVertices = candidateVertexIds.size();

		bool isEdgeNeighbour = false;
		for (std::size_t k = 0; k < nCandidateVertices; ++k) {
			if (candidateVertexIds[k] == lastEdgeVertexId) {
				isEdgeNeighbour = true;
				break;
			}
		}

		if (isEdgeNeighbour) {
			utils::addToOrderedVector<long>(candidateId, *neighs);
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
std::vector<long> PatchKernel::findCellVertexNeighs(long id, bool complete) const
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
void PatchKernel::findCellVertexNeighs(long id, bool complete, std::vector<long> *neighs) const
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
	std::unique_ptr<std::vector<long>> blackList;
	if (!complete) {
		blackList = std::unique_ptr<std::vector<long>>(new std::vector<long>());
		if (isThreeDimensional()) {
			findCellEdgeNeighs(id, true, blackList.get());
		} else {
			findCellFaceNeighs(id, true, blackList.get());
		}
	}

	for (int i = 0; i < nCellVertices; ++i) {
		_findCellVertexNeighs(id, i, blackList.get(), neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> PatchKernel::findCellVertexNeighs(long id, int vertex) const
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
void PatchKernel::findCellVertexNeighs(long id, int vertex, std::vector<long> *neighs) const
{
	_findCellVertexNeighs(id, vertex, nullptr, neighs);
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
	The blacklist has to be a pointer to a unique list of ordered cell ids
	or a null pointer if no cells should be excluded from the search
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void PatchKernel::_findCellVertexNeighs(long id, int vertex, const std::vector<long> *blackList, std::vector<long> *neighs) const
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

	while (!scanQueue.empty()) {
		// Pop a cell to process
		long scanCellId = scanQueue.back();
		const Cell &scanCell = getCell(scanCellId);
		scanQueue.pop_back();
		utils::addToOrderedVector<long>(scanCell.getId(), alreadyProcessed);

		// Get cell information
		const ReferenceElementInfo *scanCellInfo = nullptr;
		const long *scanCellConnectivity = nullptr;
		if (scanCell.hasInfo()) {
			scanCellInfo = &(scanCell.getInfo());
			scanCellConnectivity = scanCell.getConnect();
		}

		// Use face adjacencies to find vertex negihbours
		int nFaces;
		if (scanCellInfo) {
			nFaces = scanCellInfo->nFaces;
		} else {
			nFaces = scanCell.getFaceCount();
		}

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
			//
			// Moreover, to optimize the search for cells that has a reference
			// element, we work directly on the connectivity inseatd of on the
			// list of vertex ids. That's because for reference elements the
			// connectivity is just a list of vertices.
			bool faceOwnsVertex = false;
			if (scanCellInfo) {
				const ReferenceElementInfo &faceInfo = ReferenceElementInfo::getInfo(scanCellInfo->faceTypeStorage[i]);
				const int *faceLocalConnectivity = scanCellInfo->faceConnectStorage[i].data();
				const int nFaceVertices = faceInfo.nVertices;
				for (int k = 0; k < nFaceVertices; ++k) {
					long faceVertexId = scanCellConnectivity[faceLocalConnectivity[k]];
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
				if (!blackList || utils::findInOrderedVector<long>(faceNeighId, *blackList) == blackList->end()) {
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
std::vector<long> PatchKernel::findCellVertexOneRing(long id, int vertex) const
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
void PatchKernel::findCellVertexOneRing(long id, int vertex, std::vector<long> *ring) const
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
std::set<int> PatchKernel::getInternalCellPIDs()
{
	std::set<int> list;
	CellConstIterator endItr = internalCellConstEnd();
	for (CellConstIterator itr = internalCellConstBegin(); itr != endItr; ++itr) {
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
std::vector<long> PatchKernel::getInternalCellsByPID(int pid)
{
	std::vector<long> cells;
	CellConstIterator endItr = internalCellConstEnd();
	for (CellConstIterator itr = internalCellConstBegin(); itr != endItr; ++itr) {
		if (itr->getPID() == pid){
			cells.push_back(itr.getId());
		}
	}

	return cells;
}

/*!
 * Find the cells that share the specified vertex.
 *
 * This function doesn't support coincident vertices with different ids.
 *
 * \param vertexId is the index of the vertex
 * \result The cells that share the vertex.
 */
std::vector<long> PatchKernel::findVertexOneRing(long vertexId) const
{
    std::vector<long> ring;
    findVertexOneRing(vertexId, &ring);

    return ring;
}

/*!
 * Find the cells that share the specified vertex.
 *
 * This function doesn't support coincident vertices with different ids.
 *
 * \param vertexId is the index of the vertex
 * \param[in,out] ring is the vector were the one-ring of the specified vertex
 * will be stored. The vector is not cleared before adding the neighbours, it
 * is extended by appending all the neighbours found by this function.
 */
void PatchKernel::findVertexOneRing(long vertexId, std::vector<long> *ring) const
{
    // Find local id of the vertex
    //
    // Coincident vertices with different ids are not supported, this case is
    // special because a cell may contain a point with the vertex coordinates,
    // but may not contain the requested vertex (because it contains one of
    // the other coincident vertices).
    const std::array<double, 3> &coords = getVertexCoords(vertexId);
    long cellId = locatePoint(coords);
    if (cellId == bitpit::Element::NULL_ID) {
        return;
    }

    int vertexLocalId = getCell(cellId).findVertex(vertexId);
    assert(vertexLocalId >= 0);

    // Find vertex one-ring
    findCellVertexOneRing(cellId, vertexLocalId, ring);
}

/*!
	Returns true if auto-indexing is enabled for interfaces.

	When auto-indexing is disabled, ids for newly added interfaces has to be
	manually specified.

	\return Returns true if auto-indexing is enabled for interfaces.
*/
bool PatchKernel::isInterfaceAutoIndexingEnabled() const
{
	return static_cast<bool>(m_interfaceIdGenerator);
}

/*!
	Enables or disables auto-indexing for interfaces.

	When auto-indexing is disabled, ids for newly added interfaces has to be
	manually specified.

	Auto-indexing cannot be disabled if interfaces build strategy is set
	to "automatic".

	\param enabled if set to true the auto-indexing will be enabled
*/
void PatchKernel::setInterfaceAutoIndexing(bool enabled)
{
	if (isInterfaceAutoIndexingEnabled() == enabled) {
		return;
	}

	if (enabled) {
		createInterfaceIndexGenerator(true);
	} else {
		if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
			throw std::runtime_error("Auto-indexing cannot be disabled if interfaces build strategy is set to automatic.");
		}

		m_interfaceIdGenerator.reset();
	}
}

/*!
	Dump interface auto indexing.

 *  \param stream is the stream to write to
*/
void PatchKernel::dumpInterfaceAutoIndexing(std::ostream &stream) const
{
	m_interfaceIdGenerator->dump(stream);
}

/*!
	Restores interface auto indexing.

	\param stream is the stream to read from
*/
void PatchKernel::restoreInterfaceAutoIndexing(std::istream &stream)
{
	createInterfaceIndexGenerator(false);
	m_interfaceIdGenerator->restore(stream);
}

/*!
	Create the interface index generator.

	If the index generator is already created, the existing index generator will be reset.

	\param populate if set to true, the index generator will be populated with current
	interface ids
*/
void PatchKernel::createInterfaceIndexGenerator(bool populate)
{
	// Create or reset index generator
	if (!m_interfaceIdGenerator) {
		m_interfaceIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>());
	} else {
		m_interfaceIdGenerator->reset();
	}

	// Populate index generator with current ids
	if (populate) {
		InterfaceConstIterator beginItr = m_interfaces.cbegin();
		InterfaceConstIterator endItr   = m_interfaces.cend();
		for (InterfaceConstIterator itr = beginItr; itr != endItr; ++itr) {
			m_interfaceIdGenerator->setAssigned(itr.getId());
		}
	}
}

/*!
	Import the interface index generator form the specified source patch.

	If the source patch doesn't define an index generator, the intexe generator of the current
	patch will be deleted.

	\param source is the patch from with the index generator will be imported from
*/
void PatchKernel::importInterfaceIndexGenerator(const PatchKernel &source)
{
	if (source.m_interfaceIdGenerator) {
		m_interfaceIdGenerator = std::unique_ptr<IndexGenerator<long>>(new IndexGenerator<long>(*(source.m_interfaceIdGenerator)));
	} else {
		m_interfaceIdGenerator.reset();
	}
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
Interface & PatchKernel::getInterface(long id)
{
	return m_interfaces[id];
}

/*!
	Gets a constant reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A constant reference to the interface with the specified id.
*/
const Interface & PatchKernel::getInterface(long id) const
{
	return m_interfaces[id];
}

/*!
	Gets the element type for the interface with the specified id.

	\param id is the id of the requested interface
	\return The element type for the interface with the specified id.
*/
ElementType PatchKernel::getInterfaceType(long id) const
{
	return m_interfaces[id].getType();
}

/*!
	Returns an iterator pointing to the specified interface.

	\result An iterator to the specified interface.
*/
PatchKernel::InterfaceIterator PatchKernel::getInterfaceIterator(long id)
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
PatchKernel::InterfaceConstIterator PatchKernel::getInterfaceConstIterator(long id) const
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
	Adds the specified interface to the patch.

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the newly created interface.
	If a negative id value is specified, a new unique id will be generated
	for the interface
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(const Interface &source, long id)
{
	Interface interface = source;
	interface.setId(id);

	return addInterface(std::move(interface), id);
}

/*!
	Adds the specified interface to the patch.

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the newly created interface.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added interface.

*/
PatchKernel::InterfaceIterator PatchKernel::addInterface(Interface &&source, long id)
{
	// Get the id
	if (id < 0) {
		id = source.getId();
	}

	// Add a dummy interface
	//
	// It is not possible to directly add the source into the storage. First
	// a dummy interface is created and then that interface is replaced with
	// the source.
	std::unique_ptr<long[]> dummyConnectStorage = std::unique_ptr<long[]>(nullptr);

	InterfaceIterator iterator = _addInterface(ElementType::UNDEFINED, std::move(dummyConnectStorage), id);

	// Replace the newly created interface with the source
	//
	// Before replacing the newly created interface we need to set the id
	// of the source to the id that has been assigned to the newly created
	// interface.
	source.setId(iterator->getId());

	Interface &interface = (*iterator);
	interface = std::move(source);

	return iterator;
}

/*!
	Adds a new interface with the specified id.

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

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

	if (Interface::getDimension(type) > (getDimension() - 1)) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = _addInterface(type, std::move(connectStorage), id);

	return iterator;
}

/*!
	Internal function to add a new interface with the specified id.

	If valid, the specified id will we assigned to the newly created interface,
	otherwise a new unique id will be generated for the interface. However, it
	is not possible to create a new interface with an id already assigned to an
	existing interface of the patch. If this happens, an exception is thrown.
	Ids are considered valid if they are greater or equal than zero.

	\param type is the type of the interface
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
PatchKernel::InterfaceIterator PatchKernel::_addInterface(ElementType type,
														  std::unique_ptr<long[]> &&connectStorage,
														  long id)
{
	// Get the id
	if (m_interfaceIdGenerator) {
		if (id < 0) {
			id = m_interfaceIdGenerator->generate();
		} else {
			m_interfaceIdGenerator->setAssigned(id);
		}
	} else if (id < 0) {
		throw std::runtime_error("No valid id has been provided for the interface.");
	}

	// Create the interface
	PiercedVector<Interface>::iterator iterator = m_interfaces.emreclaim(id, id, type, std::move(connectStorage));

	// Set the alteration flags
	setAddedInterfaceAlterationFlags(id);

	return iterator;
}

/*!
	Set the alteration flags for an added interface.

	Only alteration flags needed for tracking the status of the patch are set
	(for example, there is no explicit flag to tag newly added interfaces).

	\param id is the id of the interface
*/
void PatchKernel::setAddedInterfaceAlterationFlags(long id)
{
	BITPIT_UNUSED(id);
}

/*!
	Resore the interface with the specified id.

	The kernel should already contain the interface, only the contents of the
	interface will be updated.

	\param type is the type of the interface
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param id is the id of the interface to restore
	\return An iterator pointing to the restored interface.
*/
PatchKernel::InterfaceIterator PatchKernel::restoreInterface(ElementType type,
															 std::unique_ptr<long[]> &&connectStorage,
															 long id)
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

	_restoreInterface(iterator, type, std::move(connectStorage));

	return iterator;
}

/*!
	Internal function to restore the interface with the specified id.

	The kernel should already contain the interface, only the contents of the
	interface will be updated.

	\param type is the type of the interface
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\return An iterator pointing to the restored interface.
*/
void PatchKernel::_restoreInterface(const InterfaceIterator &iterator, ElementType type,
									std::unique_ptr<long[]> &&connectStorage)
{
	// Restore the interface
	//
	// There is no need to set the id of the interfaces as assigned, because
	// also the index generator will be restored.
	long interfaceId = iterator.getId();
	Interface &interface = *iterator;
	interface.initialize(interfaceId, type, std::move(connectStorage));

	// Set the alteration flags
	setRestoredInterfaceAlterationFlags(interfaceId);
}

/*!
	Set the alteration flags for a restored interface.

	\param id is the id of the interface
*/
void PatchKernel::setRestoredInterfaceAlterationFlags(long id)
{
	BITPIT_UNUSED(id);
}

/*!
	Deletes an interface.

	\param id is the id of the interface
*/
bool PatchKernel::deleteInterface(long id)
{
	if (!isExpert()) {
		return false;
	}

	_deleteInterface(id);

	return true;
}

/*!
	Internal function to delete an interface.

	\param id is the id of the interface
*/
void PatchKernel::_deleteInterface(long id)
{
	// Set the alteration flags
	setDeletedInterfaceAlterationFlags(id);

	// Delete interface
	m_interfaces.erase(id, true);

	// Interface id is no longer used
	if (m_interfaceIdGenerator) {
		m_interfaceIdGenerator->trash(id);
	}
}

/*!
	Set the alteration flags for a deleted interface.

	Only alteration flags needed for tracking the status of the patch are set.

	\param id is the id of the interface
*/
void PatchKernel::setDeletedInterfaceAlterationFlags(long id)
{
	resetInterfaceAlterationFlags(id, FLAG_DELETED);
}

/*!
	Counts free interfaces within the patch.

	An interface is free if belongs to just one cell.

	\result The number of free interfaces.
*/
long PatchKernel::countFreeInterfaces() const
{
	return countBorderInterfaces();
}

/*!
	Counts border interfaces within the patch.

	An interface is border if belongs to just one cell.

	\result The number of border interfaces.
*/
long PatchKernel::countBorderInterfaces() const
{
	long nBorderInterfaces = 0;
	for (const Interface &interface : m_interfaces) {
		if (interface.getNeigh() < 0) {
			++nBorderInterfaces;
		}
        }

	return nBorderInterfaces;
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
	InterfaceConstIterator endItr = interfaceConstEnd();
	for (InterfaceConstIterator itr = interfaceConstBegin(); itr != endItr; ++itr) {
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
	InterfaceConstIterator endItr = interfaceConstEnd();
	for (InterfaceConstIterator itr = interfaceConstBegin(); itr != endItr; ++itr) {
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
	return countBorderFaces();
}

/*!
	Counts border faces within the patch.

	A face is border if a cell has no adjacent along that faces.

	\result The number of border faces.
*/
long PatchKernel::countBorderFaces() const
{
	double nBorderFaces = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nBorderFaces;
			}
		}
	}

	return nBorderFaces;
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
#if BITPIT_ENABLE_MPI==1
		utils::binary::write(stream, getVertexOwner(vertex.getId()));
#else
		int dummyOwner = 0;
		utils::binary::write(stream, dummyOwner);
#endif

		const std::array<double, 3> &coords = vertex.getCoords();
		utils::binary::write(stream, coords[0]);
		utils::binary::write(stream, coords[1]);
		utils::binary::write(stream, coords[2]);
	}

	// Dump ghost/internal subdivision
#if BITPIT_ENABLE_MPI==1
	utils::binary::write(stream, m_firstGhostVertexId);
	utils::binary::write(stream, m_lastInternalVertexId);
#else
	utils::binary::write(stream, Vertex::NULL_ID);
	utils::binary::write(stream, Vertex::NULL_ID);
#endif
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

#if BITPIT_ENABLE_MPI==1
		int owner;
		utils::binary::read(stream, owner);
#else
		int dummyOwner;
		utils::binary::read(stream, dummyOwner);
#endif

		std::array<double, 3> coords;
		utils::binary::read(stream, coords[0]);
		utils::binary::read(stream, coords[1]);
		utils::binary::read(stream, coords[2]);

#if BITPIT_ENABLE_MPI==1
		restoreVertex(std::move(coords), owner, id);
#else
		restoreVertex(std::move(coords), id);
#endif
	}

	// Restore ghost/internal subdivision
#if BITPIT_ENABLE_MPI==1
	utils::binary::read(stream, m_firstGhostVertexId);
	utils::binary::read(stream, m_lastInternalVertexId);
#else
	long dummyFirstGhostVertexId;
	long dummyLastInternalVertexId;
	utils::binary::read(stream, dummyFirstGhostVertexId);
	utils::binary::read(stream, dummyLastInternalVertexId);
#endif

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
		utils::binary::write(stream, getCellOwner(cell.getId()));
		utils::binary::write(stream, getCellHaloLayer(cell.getId()));
#else
		int dummyOwner = 0;
		utils::binary::write(stream, dummyOwner);

		int dummyHaloLayer = 0;
		utils::binary::write(stream, dummyHaloLayer);
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
	utils::binary::write(stream, m_firstGhostCellId);
	utils::binary::write(stream, m_lastInternalCellId);
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
		int owner;
		utils::binary::read(stream, owner);

		int haloLayer;
		utils::binary::read(stream, haloLayer);
#else
		int dummyOwner;
		utils::binary::read(stream, dummyOwner);

		int dummyHaloLayer;
		utils::binary::read(stream, dummyHaloLayer);
#endif

		int cellConnectSize;
		utils::binary::read(stream, cellConnectSize);

		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[cellConnectSize]);
		for (int k = 0; k < cellConnectSize; ++k) {
			utils::binary::read(stream, cellConnect[k]);
		}

		CellIterator iterator;
#if BITPIT_ENABLE_MPI==1
		iterator = restoreCell(type, std::move(cellConnect), owner, haloLayer, id);
#else
		iterator = restoreCell(type, std::move(cellConnect), id);
#endif
		iterator->setPID(PID);
	}

	// Restore ghost/internal subdivision
#if BITPIT_ENABLE_MPI==1
	utils::binary::read(stream, m_firstGhostCellId);
	utils::binary::read(stream, m_lastInternalCellId);
#else
	long dummyFirstGhostCellId;
	long dummyLastInternalCellId;
	utils::binary::read(stream, dummyFirstGhostCellId);
	utils::binary::read(stream, dummyLastInternalCellId);
#endif

	// Update adjacencies
	updateAdjacencies();

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

	// Interfaces need up-to-date adjacencies
	updateAdjacencies();

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

	// Interfaces are now updated
	unsetCellAlterationFlags(FLAG_INTERFACES_DIRTY);
	m_alteredInterfaces.clear();

	// Set original advanced editing status
	setExpert(originalExpertStatus);
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
	Alter the patch performing the adaption.

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
bool PatchKernel::_markCellForRefinement(long id)
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
bool PatchKernel::_markCellForCoarsening(long id)
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
bool PatchKernel::_resetCellAdaptionMarker(long id)
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
adaption::Marker PatchKernel::_getCellAdaptionMarker(long id)
{
	BITPIT_UNUSED(id);

	return adaption::MARKER_UNDEFINED;
}

/*!
	Enables cell balancing.

	Default implementation is a no-op function.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
	\result Returns true if the flag was properly set, false otherwise.
*/
bool PatchKernel::_enableCellBalancing(long id, bool enabled)
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

	// Sort internal vertices
	if (m_nInternalVertices > 0) {
		m_vertices.sortBefore(m_lastInternalVertexId, true);
		updateLastInternalVertexId();
	}

#if BITPIT_ENABLE_MPI==1
	// Sort ghost cells
	if (m_nGhostVertices > 0) {
		m_vertices.sortAfter(m_firstGhostVertexId, true);
		updateFirstGhostVertexId();
	}
#endif

	// Synchronize storage
	m_vertices.sync();

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
	if (m_nInternalCells > 0) {
		m_cells.sortBefore(m_lastInternalCellId, true);
		updateLastInternalCellId();
	}

#if BITPIT_ENABLE_MPI==1
	// Sort ghost cells
	if (m_nGhostCells > 0) {
		m_cells.sortAfter(m_firstGhostCellId, true);
		updateFirstGhostCellId();
	}
#endif

	// Synchronize storage
	m_cells.sync();

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

	m_interfaces.sync();

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

	m_vertices.sync();

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

	m_cells.sync();

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

	m_interfaces.sync();

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
	bool status = true;

	// Alteration flags
	if (m_alteredCells.empty()) {
		AlterationFlagsStorage().swap(m_alteredCells);
	}

	if (m_alteredInterfaces.empty()) {
		AlterationFlagsStorage().swap(m_alteredInterfaces);
	}

	// Vertices
	status |= squeezeVertices();

	// Cells
	status |= squeezeCells();

	// Interfaces
	status |= squeezeInterfaces();

	return status;
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> PatchKernel::evalCellCentroid(long id) const
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

	evalElementBoundingBox(cell, minPoint, maxPoint);
}

/*!
	Get vertex coordinates of the specified cell.

	\param id is the id of the cell
	\result Vertex coordinates of the cell.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getCellVertexCoordinates(long id) const
{
	// Get the vertices ids
	const Cell &cell = getCell(id);
	ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
	const int nCellVertices = cellVertexIds.size();

	// Build the proxy vector with the coordinates
	ConstProxyVector<std::array<double, 3>> vertexCoordinates(ConstProxyVector<std::array<double, 3>>::INTERNAL_STORAGE, nCellVertices);
	ConstProxyVector<std::array<double, 3>>::storage_pointer storage = vertexCoordinates.storedData();
	getVertexCoords(cellVertexIds.size(), cellVertexIds.data(), storage);

	return vertexCoordinates;
}

/*!
	Get vertex coordinates of the specified cell.

	\param id is the id of the cell
	\param[out] coordinates on output will contain the vertex coordinates
*/
void PatchKernel::getCellVertexCoordinates(long id, std::unique_ptr<std::array<double, 3>[]> *coordinates) const
{
	const Cell &cell = getCell(id);

	getElementVertexCoordinates(cell, coordinates);
}

/*!
	Get vertex coordinates of the specified cell.

	\param id is the id of the cell
	\param[out] coordinates on output will contain the vertex coordinates, it
	is up to the caller to ensure that the storage has enough space for all
	the vertex coordinates
*/
void PatchKernel::getCellVertexCoordinates(long id, std::array<double, 3> *coordinates) const
{
	const Cell &cell = getCell(id);

	getElementVertexCoordinates(cell, coordinates);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> PatchKernel::evalInterfaceCentroid(long id) const
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

	evalElementBoundingBox(interface, minPoint, maxPoint);
}

/*!
	Get vertex coordinates of the specified interface.

	\param id is the id of the interface
	\result Vertex coordinates of the interface.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getInterfaceVertexCoordinates(long id) const
{
	// Get the vertices ids
	const Interface &interface = getInterface(id);
	ConstProxyVector<long> interfaceVertexIds = interface.getVertexIds();
	const int nInterfaceVertices = interfaceVertexIds.size();

	// Build the proxy vector with the coordinates
	ConstProxyVector<std::array<double, 3>> vertexCoordinates(ConstProxyVector<std::array<double, 3>>::INTERNAL_STORAGE, nInterfaceVertices);
	ConstProxyVector<std::array<double, 3>>::storage_pointer storage = vertexCoordinates.storedData();
	getVertexCoords(interfaceVertexIds.size(), interfaceVertexIds.data(), storage);

	return vertexCoordinates;
}

/*!
	Get vertex coordinates of the specified interface.

	\param id is the id of the interface
	\param[out] coordinates on output will contain the vertex coordinates
*/
void PatchKernel::getInterfaceVertexCoordinates(long id, std::unique_ptr<std::array<double, 3>[]> *coordinates) const
{
	const Interface &interface = getInterface(id);

	getElementVertexCoordinates(interface, coordinates);
}

/*!
	Get vertex coordinates of the specified interface.

	\param id is the id of the interface
	\param[out] coordinates on output will contain the vertex coordinates, it
	is up to the caller to ensure that the storage has enough space for all
	the vertex coordinates
*/
void PatchKernel::getInterfaceVertexCoordinates(long id, std::array<double, 3> *coordinates) const
{
	const Interface &interface = getInterface(id);

	getElementVertexCoordinates(interface, coordinates);
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
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	std::size_t nCellVertices = elementVertexIds.size();
	BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
	getVertexCoords(nCellVertices, elementVertexIds.data(), vertexCoordinates);

	return element.evalCentroid(vertexCoordinates);
}

/*!
	Evaluates the bounding box of the specified element.

	\param element is the element
	\param[out] minPoint is the minimum point of the bounding box
	\param[out] maxPoint is the maximum point of the bounding box
*/
void PatchKernel::evalElementBoundingBox(const Element &element, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const
{
	ElementType elementType = element.getType();
	switch (elementType)
	{

	case ElementType::VERTEX:
	{
		const long *elementConnect = element.getConnect();
		*minPoint = getVertexCoords(elementConnect[0]);
		*maxPoint = *minPoint;
		break;
	}

	case ElementType::PIXEL:
	{
		const long *elementConnect = element.getConnect();

		const std::array<double, 3> &vertexCoord_0 = getVertexCoords(elementConnect[0]);
		const std::array<double, 3> &vertexCoord_3 = getVertexCoords(elementConnect[3]);

		for (int d = 0; d < 3; ++d) {
			(*minPoint)[d] = std::min(vertexCoord_0[d], vertexCoord_3[d]);
			(*maxPoint)[d] = std::max(vertexCoord_0[d], vertexCoord_3[d]);
		}

		break;
	}

	case ElementType::VOXEL:
	{
		const long *elementConnect = element.getConnect();

		const std::array<double, 3> &vertexCoord_0 = getVertexCoords(elementConnect[0]);
		const std::array<double, 3> &vertexCoord_7 = getVertexCoords(elementConnect[7]);

		for (int d = 0; d < 3; ++d) {
			(*minPoint)[d] = std::min(vertexCoord_0[d], vertexCoord_7[d]);
			(*maxPoint)[d] = std::max(vertexCoord_0[d], vertexCoord_7[d]);
		}

		break;
	}

	default:
	{
		ConstProxyVector<long> elementVertexIds = element.getVertexIds();
		const int nElementVertices = elementVertexIds.size();

		*minPoint = getVertexCoords(elementVertexIds[0]);
		*maxPoint = *minPoint;
		for (int i = 1; i < nElementVertices; ++i) {
			const std::array<double, 3> &vertexCoord = getVertexCoords(elementVertexIds[i]);
			for (int d = 0; d < 3; ++d) {
				(*minPoint)[d] = std::min(vertexCoord[d], (*minPoint)[d]);
				(*maxPoint)[d] = std::max(vertexCoord[d], (*maxPoint)[d]);
			}
		}

		break;
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
long PatchKernel::locatePoint(double x, double y, double z) const
{
	return locatePoint({{x, y, z}});
}

/*!
 * Check whether the face "face_A" on cell "cell_A" is the same as the face
 * "face_B" on cell "cell_B".
 * 
 * \param[in] cell_A is the first cell
 * \param[in] face_A is the face on the first cell
 * \param[in] cell_B is the the second cell
 * \param[in] face_B is the face on the second cell
 * \result Returns true if the two faces are the same.
*/
bool PatchKernel::isSameFace(const Cell &cell_A, int face_A, const Cell &cell_B, int face_B) const
{
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
	Checks if the adjacencies are dirty.

	\param global if set to true, the dirty status will be evaluated globally
	across all the partitions
	\result Returns true if the adjacencies are dirty, false otherwise.
*/
bool PatchKernel::areAdjacenciesDirty(bool global) const
{
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_NONE) {
		return false;
	}

	bool areDirty = false;
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
			areDirty = true;
			break;
		}
	}

#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &areDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif

	return areDirty;
}

/*!
	Fill adjacencies info for each cell.

	If adjacencies are already built, all adjacencies will be deleted and
	they will be re-generated from scratch.
*/
void PatchKernel::buildAdjacencies()
{
    initializeAdjacencies(ADJACENCIES_AUTOMATIC);
}

/*!
	Initialize the adjacencies using the specified build strategy.

	If the current strategy doesn't match the requested strategy, all
	adjacencies will be deleted and they will be re-generated from scratch.

	\param strategy is the build strategy that will be used
*/
void PatchKernel::initializeAdjacencies(AdjacenciesBuildStrategy strategy)
{
	// Get current build stategy
	AdjacenciesBuildStrategy currentStrategy = getAdjacenciesBuildStrategy();

	// Early return if we don't need adjacencies
	if (strategy == ADJACENCIES_NONE) {
		if (currentStrategy != ADJACENCIES_NONE) {
			destroyAdjacencies();
		}

		return;
	}

	// Update the build strategy
	if (currentStrategy != strategy) {
		setAdjacenciesBuildStrategy(strategy);
	}

	// Reset adjacencies
	resetAdjacencies();

	// Update the adjacencies
	updateAdjacencies();
}

/*!
	Update the adjacencies of the patch.

	\param forcedUpdated if set to true, bounding box information will be
	updated also if they are not marked as dirty
*/
void PatchKernel::updateAdjacencies(bool forcedUpdated)
{
	// Early return if adjacencies are not built
	AdjacenciesBuildStrategy currentStrategy = getAdjacenciesBuildStrategy();
	if (currentStrategy == ADJACENCIES_NONE) {
		return;
	}

	// Check if the adjacencies are dirty
	bool adjacenciesDirty = areAdjacenciesDirty();
	if (!adjacenciesDirty && !forcedUpdated) {
		return;
	}

	// Update adjacencies
	if (adjacenciesDirty) {
		// Prune stale adjacencies
		pruneStaleAdjacencies();

		// Update adjacencies
		_updateAdjacencies();

		// Adjacencies are now updated
		unsetCellAlterationFlags(FLAG_ADJACENCIES_DIRTY);
	} else {
		initializeAdjacencies(currentStrategy);
	}
}

/*!
	Destroy the adjacencies.

	After deleting the adjacencies, this function changes the build strategy
	to "None".
*/
void PatchKernel::destroyAdjacencies()
{
	// Early return if no adjacencies have been built
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_NONE) {
		return;
	}

	// Destroy the adjacencies
	_resetAdjacencies(true);

	// Clear list of cells with dirty adjacencies
	unsetCellAlterationFlags(FLAG_ADJACENCIES_DIRTY);

	// Set adjacencies status
	setAdjacenciesBuildStrategy(ADJACENCIES_NONE);
}

/*!
	Reset the adjacencies.

	This function doesn't change the build strategy, it only resets the
	existing adjacencies.
*/
void PatchKernel::resetAdjacencies()
{
	// Early return if no adjacencies have been built
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_NONE) {
		return;
	}

	// Reset adjacencies
	_resetAdjacencies(false);

	// All cells have now dirty adjacencies
	setCellAlterationFlags(FLAG_ADJACENCIES_DIRTY);
}

/*!
	Internal function to reset the adjacencies.

	This function doesn't change the alteration flags.

	\param release if it's true the memory hold by the adjacencies will be
	released, otherwise the adjacencies will be reset but their memory will
	not be relased
*/
void PatchKernel::_resetAdjacencies(bool release)
{
	for (Cell &cell : m_cells) {
		cell.resetAdjacencies(!release);
	}
}

/*!
	Prune stale adjacencies.

	The list of cells to process and the list of stale adjacencies are filled
	during cell deletion.
*/
void PatchKernel::pruneStaleAdjacencies()
{
	// Early return if adjacencies are not built
	AdjacenciesBuildStrategy currentStrategy = getAdjacenciesBuildStrategy();
	if (currentStrategy == ADJACENCIES_NONE) {
		return;
	}

	// Remove stale adjacencies
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
			continue;
		}

		long cellId = entry.first;
		Cell &cell = m_cells.at(cellId);
		if (cell.getAdjacencyCount() == 0) {
			continue;
		}

		int nCellFaces = cell.getFaceCount();
		for (int face = nCellFaces - 1; face >= 0; --face) {
			long *faceAdjacencies = cell.getAdjacencies(face);
			int nFaceAdjacencies = cell.getAdjacencyCount(face);
			for (int i = nFaceAdjacencies - 1; i >= 0; --i) {
				long adjacency = faceAdjacencies[i];
				if (!testCellAlterationFlags(adjacency, FLAG_DELETED)) {
					continue;
				}

				cell.deleteAdjacency(face, i);
			}
		}
	}
}

/*!
	Internal function to update the adjacencies of the patch.

	In addition to the cells whose adjacencies are marked as dirty, also the
	adjacencies of their neighbours will be updated.

	This implementation can NOT handle hanging nodes.
*/
void PatchKernel::_updateAdjacencies()
{
	int dimension = getDimension();

	// Check if multiple adjacencies are allowed
	//
	// On a three-dimensional patch each internal face is shared between two,
	// and only two, half-faces. On lower dimension patches one internal face
	// may be shared among multiple half-faces (this happens with non-manifold
	// patches).
	bool multipleMatchesAllowed = (dimension < 3);

	// Define matching windings
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
	std::vector<CellHalfFace::Winding> matchingWindings;
	matchingWindings.push_back(CellHalfFace::WINDING_REVERSE);
	if (multipleMatchesAllowed) {
		matchingWindings.push_back(CellHalfFace::WINDING_NATURAL);
	}

	// Count cells with dirty adjacencies
	long nDirtyAdjacenciesCells = 0;
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
			continue;
		}

		++nDirtyAdjacenciesCells;
	}

	// Get the vertices of cells with dirty adjacencies
	//
	// The list is only needed if multiple matches are allowed and there are
	// cells whose adjacencies don't need to be updated.
	std::unique_ptr<PiercedStorage<bool, long>> dirtyAdjacenciesVertices;
	if (multipleMatchesAllowed && (nDirtyAdjacenciesCells != getCellCount())) {
		dirtyAdjacenciesVertices = std::unique_ptr<PiercedStorage<bool, long>>(new PiercedStorage<bool, long>(1, &m_vertices));
		dirtyAdjacenciesVertices->fill(false);
		for (const auto &entry : m_alteredCells) {
			AlterationFlags cellAlterationFlags = entry.second;
			if (!testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
				continue;
			}

			long cellId = entry.first;
			const Cell &cell = m_cells.at(cellId);
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			for (long vertexId : cellVertexIds) {
				dirtyAdjacenciesVertices->at(vertexId) = true;
			}
		}
	}

	// List the cells that needs to be processed
	//
	// The list should contain the cells whose adjacencies need to be updated
	// and their neighbours candidates.
	//
	// Cells whose adjacencies needs to be updated are cells having the dirty
	// adjaceny flag set.
	//
	// Neighbours candidates are all border cells and, if multiple matches are
	// allowed, cells that share vertices with a cell with dirty adjacencies.
	std::vector<Cell *> processList;
	processList.reserve(nDirtyAdjacenciesCells);

	std::size_t nMaxHalfFaces = 0;
	for (Cell &cell : m_cells) {
		// Get cell information
		long cellId = cell.getId();
		int nCellFaces = cell.getFaceCount();

		// Add cells with dirty adjacencies
		if (testCellAlterationFlags(cellId, FLAG_ADJACENCIES_DIRTY)) {
			nMaxHalfFaces += nCellFaces;
			processList.push_back(&cell);
			continue;
		}

		// Add border cells
		bool isBorderCell = false;
		for (int face = 0; face < nCellFaces; face++) {
			if (cell.isFaceBorder(face)) {
				isBorderCell = true;
				continue;
			}
		}

		if (isBorderCell) {
			nMaxHalfFaces += nCellFaces;
			processList.push_back(&cell);
			continue;
		}

		// Add cell that share vertices with a cell with dirty adjacencies
		//
		// For performace reasons, the check is a bit rough. All cells that
		// share a number of vertices at least equal to the dimension af the
		// patch will be added to the process list. This will add also cells
		// that are not neighbour candidates, but it's faster to add some false
		// positives, rather trying to filter them out.
		if (multipleMatchesAllowed) {
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

			int nCelldirtyAdjacenciesVertices = 0;
			bool futureNeighbourCandidate = false;
			for (long vertexId : cellVertexIds) {
				if (dirtyAdjacenciesVertices->at(vertexId)) {
					++nCelldirtyAdjacenciesVertices;
					if (nCelldirtyAdjacenciesVertices >= dimension) {
						futureNeighbourCandidate = true;
						break;
					}
				}
			}

			if (futureNeighbourCandidate) {
				nMaxHalfFaces += nCellFaces;
				processList.push_back(&cell);
			}
		}
	}

	// Create the adjacencies
	std::unordered_set<CellHalfFace, CellHalfFace::Hasher> halfFaces;
	halfFaces.reserve(0.5 * nMaxHalfFaces);

	std::vector<std::pair<Cell *, int>> matchingAdjacencies;
	for (Cell *cell : processList) {
		long cellId = cell->getId();
		const int nCellFaces = cell->getFaceCount();
		bool areCellAdjacenciesDirty = testCellAlterationFlags(cellId, FLAG_ADJACENCIES_DIRTY);
		for (int face = 0; face < nCellFaces; face++) {
			// Generate the half-face
			CellHalfFace halfFace(*cell, face);

			// Find matching half-face
			auto matchingHalfFaceItr = halfFaces.end();
			for (CellHalfFace::Winding winding : matchingWindings) {
				// Set winding order
				halfFace.setWinding(winding);

				// Search for matching half-face
				matchingHalfFaceItr = halfFaces.find(halfFace);
				if (matchingHalfFaceItr != halfFaces.end()) {
					break;
				}
			}

			// Early return if no matching half-face has been found
			//
			// If no matching half-face has been found, we need to add the
			// current half-face to the list because it may match the face
			// of a cell that will be processed later.
			if (matchingHalfFaceItr == halfFaces.end()) {
				halfFace.setWinding(CellHalfFace::WINDING_NATURAL);
				halfFaces.insert(std::move(halfFace));
				continue;
			}

			// Get matching half-face information
			const CellHalfFace &matchingHalfFace = *matchingHalfFaceItr;
			Cell &matchingCell  = matchingHalfFace.getCell();
			long matchingCellId = matchingCell.getId();
			long matchingFace   = matchingHalfFace.getFace();

			// Only process cells with dirty adjacencies
			//
			// In order to limit the half faces that need to be stored, we are
			// processing at the same time cells with dirty adjacencies and
			// their neighbour candidates. This means we will also find matches
			// between faces of cell whose adjacensies are already up-to-date.
			bool areMatchingCellAdjacenciesDirty = testCellAlterationFlags(matchingCellId, FLAG_ADJACENCIES_DIRTY);
			if (!areCellAdjacenciesDirty && !areMatchingCellAdjacenciesDirty) {
				continue;
			}

			// Identifty matching adjacencies
			//
			// If multiple matches are allowed, in addition to the entry
			// related to the matching half face, we also need to add the
			// entries associated to the adjacencies already identified
			// for the half face.
			matchingAdjacencies.clear();
			matchingAdjacencies.emplace_back(std::make_pair<Cell *, int>(&matchingCell, matchingFace));
			if (multipleMatchesAllowed) {
				const int nMachingFaceNeighs = matchingCell.getAdjacencyCount(matchingFace);
				const long *machingFaceNeighs = matchingCell.getAdjacencies(matchingFace);
				for (int k = 0; k < nMachingFaceNeighs; ++k) {
					long neighId = machingFaceNeighs[k];
					if (neighId == cellId) {
						continue;
					}

					Cell &neigh = m_cells.at(neighId);
					int neighFace = findAdjoinNeighFace(matchingCell, matchingFace, neigh);
					matchingAdjacencies.emplace_back(std::make_pair<Cell *, int>(&neigh, std::move(neighFace)));
				}
			} else {
				// When cell adjacencies are makred dirty this means that some
				// adjacencies are dirty not that all adjacencies are dirty.
				// The matchin cell may have no adjaceny associated to the
				// matching face or a signle adjacency that should match
				// the current cell.
				assert(matchingCell.getAdjacencyCount(matchingFace) == 0 || (matchingCell.getAdjacencyCount(matchingFace) == 1 && (*(matchingCell.getAdjacencies(matchingFace)) == cellId)));
			}

			// Create adjacencies
			for (const std::pair<Cell *, int> &matchingAdjacency : matchingAdjacencies) {
				Cell *adjacentCell = matchingAdjacency.first;
				long adjacentCellId = adjacentCell->getId();
				int adjacentFace = matchingAdjacency.second;

				cell->pushAdjacency(face, adjacentCellId);
				adjacentCell->pushAdjacency(adjacentFace, cellId);
			}

			// Remove the matching half-face from the list
			//
			// If multiple matches are allowed, we need to keep the half
			// face, because that entry may match the face of a cell that
			// will be processed later.
			if (!multipleMatchesAllowed) {
				halfFaces.erase(matchingHalfFaceItr);
			}
		}
	}
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
	if (status == INTERFACES_AUTOMATIC) {
		if (!isInterfaceAutoIndexingEnabled()) {
			throw std::runtime_error("Automatic build strategy requires auto-indexing.");
		}
	}

	m_interfacesBuildStrategy = status;
}

/*!
	Checks if the interfaces are dirty.

	\param global if set to true, the dirty status will be evaluated globally
	across all the partitions
	\result Returns true if the interfaces are dirty, false otherwise.
*/
bool PatchKernel::areInterfacesDirty(bool global) const
{
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return false;
	}

	bool areDirty = !m_alteredInterfaces.empty();
	if (!areDirty) {
		for (const auto &entry : m_alteredCells) {
			AlterationFlags cellAlterationFlags = entry.second;
			if (testAlterationFlags(cellAlterationFlags, FLAG_INTERFACES_DIRTY)) {
				areDirty = true;
				break;
			}
		}
	}

#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &areDirty, 1, MPI_C_BOOL, MPI_LOR, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif

	return areDirty;
}

/*!
	Build interfaces among the cells.

	If interfaces are already built, all interfaces will be deleted and
	they will be re-generated from scratch.
*/
void PatchKernel::buildInterfaces()
{
    initializeInterfaces(INTERFACES_AUTOMATIC);
}

/*!
	Initialize the interfaces using the specified build strategy.

	If the current strategy doesn't match the requested strategy, all
	interfaces will be deleted and they will be re-generated from scratch.

	Adjacencies are a mandatory requirement for building interfaces, if the
	adjacencies are not yet initialized an exception is thrown.

	\param strategy is the build strategy that will be used
*/
void PatchKernel::initializeInterfaces(InterfacesBuildStrategy strategy)
{
	// Interfaces need adjacencies
	if (getAdjacenciesBuildStrategy() == ADJACENCIES_NONE) {
		throw std::runtime_error ("Adjacencies are mandatory for building the interfaces.");
	}

	// Get current build stategy
	InterfacesBuildStrategy currentStrategy = getInterfacesBuildStrategy();

	// Early return if we don't need interfaces
	if (strategy == INTERFACES_NONE) {
		if (currentStrategy != INTERFACES_NONE) {
			destroyInterfaces();
		}

		return;
	}

	// Update the build strategy
	if (currentStrategy != strategy) {
		setInterfacesBuildStrategy(strategy);
	}

	// Reset interfaces
	resetInterfaces();

	// Update the interfaces
	updateInterfaces();
}

/*!
	Update the interfaces of the patch.

	\param forcedUpdated if set to true, bounding box information will be
	updated also if they are not marked as dirty
*/
void PatchKernel::updateInterfaces(bool forcedUpdated)
{
	// Early return if interfaces are not built
	InterfacesBuildStrategy currentStrategy = getInterfacesBuildStrategy();
	if (currentStrategy == INTERFACES_NONE) {
		return;
	}

	// Check if the interfaces are dirty
	bool interfacesDirty = areInterfacesDirty();
	if (!interfacesDirty && !forcedUpdated) {
		return;
	}

	// Interfaces need up-to-date adjacencies
	assert(getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);
	updateAdjacencies();

	// Update interfaces
	if (interfacesDirty) {
		// Enable advanced editing
		bool originalExpertStatus = isExpert();
		setExpert(true);

		// Prune stale interfaces
		pruneStaleInterfaces();

		// Update interfaces
		_updateInterfaces();

		// Interfaces are now updated
		unsetCellAlterationFlags(FLAG_INTERFACES_DIRTY);
		m_alteredInterfaces.clear();

		// Set original advanced editing status
		setExpert(originalExpertStatus);
	} else {
		initializeInterfaces(currentStrategy);
	}
}

/*!
	Destroy the interfaces.

	After deleting the interfaces, this function changes the build strategy
	to "None".
*/
void PatchKernel::destroyInterfaces()
{
	// Early return if no interfaces have been built
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Destroy the interfaces
	_resetInterfaces(true);

	// Clear list of cells with dirty adjacencies
	unsetCellAlterationFlags(FLAG_INTERFACES_DIRTY);

	// Clear list of altered interfaces
	m_alteredInterfaces.clear();

	// Set interface build strategy
	setInterfacesBuildStrategy(INTERFACES_NONE);
}

/*!
	Prune stale interfaces.

	The list of cells to process and the list of stale interfaces are filled
	during cell deletion.
*/
void PatchKernel::pruneStaleInterfaces()
{
	// Early return if no interfaces have been built
	if (getInterfacesBuildStrategy() == INTERFACES_NONE) {
		return;
	}

	// Remove dangling interfaces from cells
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_INTERFACES_DIRTY)) {
			continue;
		}

		long cellId = entry.first;
		Cell &cell = m_cells.at(cellId);
		if (cell.getInterfaceCount() == 0) {
			continue;
		}

		int nCellFaces = cell.getFaceCount();
		for (int face = nCellFaces - 1; face >= 0; --face) {
			long *faceInterfaces = cell.getInterfaces(face);
			int nFaceInterfaces = cell.getInterfaceCount(face);
			for (int i = nFaceInterfaces - 1; i >= 0; --i) {
				long interfaceId = faceInterfaces[i];
				if (!testInterfaceAlterationFlags(interfaceId, FLAG_DANGLING)) {
					continue;
				}

				// Delete the interface from the cell
				cell.deleteInterface(face, i);
			}
		}
	}

	// Delete dangling interfaces
	std::vector<long> danglingInterfaces;
	danglingInterfaces.reserve(m_alteredInterfaces.size());
	for (const auto &entry : m_alteredInterfaces) {
		AlterationFlags interfaceAlterationFlags = entry.second;
		if (!testAlterationFlags(interfaceAlterationFlags, FLAG_DANGLING)) {
			continue;
		}

		long interfaceId = entry.first;
		danglingInterfaces.push_back(interfaceId);
	}
	deleteInterfaces(danglingInterfaces);
}

/*!
	Internal function to update the interfaces of the patch.

	The function will process the cells whose interfaces have been marked as
	dirty.
*/
void PatchKernel::_updateInterfaces()
{
	//
	// Update interfaces
	//
	// Adjacencies and interfaces of a face are paired: the i-th face adjacency
	// corresponds to the i-th face interface. Moreover if we loop through the
	// adjacencies of a face, the adjacencies that have an interface are always
	// listed first. This meas that, to update the interfaces, we can count the
	// interfaces already associated to a face and loop only on the adjacencies
	// which have an index past the one of the last interface.
	//
	// On border faces of internal cells we need to build an interface, also
	// if there are no adjacencies.
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_INTERFACES_DIRTY)) {
			continue;
		}

		long cellId = entry.first;
		Cell &cell = m_cells.at(cellId);
		const int nCellFaces = cell.getFaceCount();
		for (int face = 0; face < nCellFaces; face++) {
			int nFaceInterfaces= cell.getInterfaceCount(face);

			bool isFaceBorder = cell.isFaceBorder(face);
			if (!isFaceBorder) {
				// Find the range of adjacencies that need an interface
				int updateEnd   = cell.getAdjacencyCount(face);
				int updateBegin = nFaceInterfaces;
				if (updateBegin == updateEnd) {
					continue;
				}

				// Build an interface for every adjacency
				const long *faceAdjacencies = cell.getAdjacencies(face);
				for (int k = updateBegin; k < updateEnd; ++k) {
					long neighId = faceAdjacencies[k];
					Cell *neigh  = &m_cells[neighId];

					int neighFace = findAdjoinNeighFace(cell, face, *neigh);

					buildCellInterface(&cell, face, neigh, neighFace);
				}
			} else if (nFaceInterfaces == 0) {
				// Internal borderes need an interface
				buildCellInterface(&cell, face, nullptr, -1);
			}
		}
	}
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
	// precise comparison, it's only necessary to define a repeatable order
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
			assert(ownerPairedAdjacencyIndex >= 0);
			intrOwner->setAdjacency(intrOwnerFace, ownerInterfaceIndex, intrNeighId);
			intrOwner->setAdjacency(intrOwnerFace, ownerPairedAdjacencyIndex, ownerPairedAdjacency);
		}

		// Fix adjacency order on the neighbour cell
		int neighInterfaceIndex   = intrNeigh->getInterfaceCount(intrNeighFace) - 1;
		long neighPairedAdjacency = intrNeigh->getAdjacency(intrNeighFace, neighInterfaceIndex);
		if (neighPairedAdjacency != intrOwnerId) {
			int neighPairedAdjacencyIndex = intrNeigh->findAdjacency(intrNeighFace, intrOwnerId);
			assert(neighPairedAdjacencyIndex >= 0);
			intrNeigh->setAdjacency(intrNeighFace, neighInterfaceIndex, intrOwnerId);
			intrNeigh->setAdjacency(intrNeighFace, neighPairedAdjacencyIndex, neighPairedAdjacency);
		}
	}

	return interfaceIterator;
}

/*!
	Finds the face of the supposed neighbour that adjoins the target face.

	\param cell is the cell
	\param cellFace is the target face of the cell
	\param neigh is the supposed neighbour of the cell
	\result The face of the neighbour which adjoins the target face. If the
	two cells are not neighbours, a negative number is returned.
 */
int PatchKernel::findAdjoinNeighFace(const Cell &cell, int cellFace, const Cell &neigh) const
{
	// Evaluate list of candidate faces
	//
	// The cells may be neighbours through multiple faces. Identify which face
	// matches the target one is quite expensive. Moreover, if the cells are
	// neighbours through a single face (which is the common case), the check
	// is not needed at all. Therefore, first we identify all the faces through
	// which the two cells are neighbours and then, if and only if there are
	// multiple candidates, we identify the face that matches the target one.
	long cellId = cell.getId();
	const int nNeighFaces = neigh.getFaceCount();

	int nCandidates = 0;
	BITPIT_CREATE_WORKSPACE(candidates, std::array<double BITPIT_COMMA 3>, nNeighFaces, ReferenceElementInfo::MAX_ELEM_FACES);
	for (int neighFace = 0; neighFace < nNeighFaces; neighFace++) {
		int nFaceAdjacencies = neigh.getAdjacencyCount(neighFace);
		const long *faceAdjacencies = neigh.getAdjacencies(neighFace);
		for (int k = 0; k < nFaceAdjacencies; ++k) {
			long geussId = faceAdjacencies[k];
			if (geussId == cellId) {
				(*candidates)[nCandidates] = neighFace;
				++nCandidates;
				break;
			}
		}
	}

	// Select the face of the neighbour which adjoin the given face
	if (nCandidates == 1) {
		return (*candidates)[0];
	} else {
		for (int i = 0; i < nCandidates; ++i) {
			int candidateFace = (*candidates)[i];
			if (isSameFace(cell, cellFace, neigh, candidateFace)) {
				return candidateFace;
			}
		}
	}

	return -1;
}

/*!
	Test the specified alteration flags for the given cell.

	\param id is the id of the cell
	\param flags are the flags that will be tested
	\result Return true if the flags are set, false otherwise.
*/
bool PatchKernel::testCellAlterationFlags(long id, AlterationFlags flags) const
{
	return testElementAlterationFlags(id, flags, m_alteredCells);
}

/*!
	Get the alteration flags of the specified cell.

	\param id is the id of the cell
	\result The alteration flags of the cell.
*/
PatchKernel::AlterationFlags PatchKernel::getCellAlterationFlags(long id) const
{
	return getElementAlterationFlags(id, m_alteredCells);
}

/*!
	Reset the alteration flags of the given cell.

	\param id is the id of the cell
	\param flags are the flags that will be set
*/
void PatchKernel::resetCellAlterationFlags(long id, AlterationFlags flags)
{
	resetElementAlterationFlags(id, flags, &m_alteredCells);
}

/*!
	Set the specified alteration flags for all the cells.

	\param flags are the flags that will be set
*/
void PatchKernel::setCellAlterationFlags(AlterationFlags flags)
{
	CellConstIterator endItr = cellConstEnd();
	for (CellConstIterator itr = cellConstBegin(); itr != endItr; ++itr) {
		setCellAlterationFlags(itr.getId(), flags);
	}
}

/*!
	Set the specified alteration flags for the given cell.

	\param id is the id of the cell
	\param flags are the flags that will be set
*/
void PatchKernel::setCellAlterationFlags(long id, AlterationFlags flags)
{
	setElementAlterationFlags(id, flags, &m_alteredCells);
}

/*!
	Unset the specified alteration flags for all the altered cells.

	\param flags are the flags that will be unset
*/
void PatchKernel::unsetCellAlterationFlags(AlterationFlags flags)
{
	unsetElementAlterationFlags(flags, &m_alteredCells);
}

/*!
	Unset the specified alteration flags for the given cell.

	\param id is the id of the cell
	\param flags are the flags that will be unset
*/
void PatchKernel::unsetCellAlterationFlags(long id, AlterationFlags flags)
{
	unsetElementAlterationFlags(id, flags, &m_alteredCells);
}

/*!
	Test the specified alteration flags for the given interface.

	\param id is the id of the interface
	\param flags are the flags that will be tested
	\result Return true if the flags are set, false otherwise.
*/
bool PatchKernel::testInterfaceAlterationFlags(long id, AlterationFlags flags) const
{
	return testElementAlterationFlags(id, flags, m_alteredInterfaces);
}

/*!
	Get the alteration flags of the specified interface.

	\param id is the id of the interface
	\result The alteration flags of the interface.
*/
PatchKernel::AlterationFlags PatchKernel::getInterfaceAlterationFlags(long id) const
{
	return getElementAlterationFlags(id, m_alteredInterfaces);
}

/*!
	Reset the alteration flags of the given interface.

	\param id is the id of the interface
	\param flags are the flags that will be set
*/
void PatchKernel::resetInterfaceAlterationFlags(long id, AlterationFlags flags)
{
	resetElementAlterationFlags(id, flags, &m_alteredInterfaces);
}

/*!
	Set the specified alteration flags for all the interfaces.

	\param flags are the flags that will be set
*/
void PatchKernel::setInterfaceAlterationFlags(AlterationFlags flags)
{
	InterfaceConstIterator endItr = interfaceConstEnd();
	for (InterfaceConstIterator itr = interfaceConstBegin(); itr != endItr; ++itr) {
		setInterfaceAlterationFlags(itr.getId(), flags);
	}
}

/*!
	Set the specified alteration flags for the given interface.

	\param id is the id of the interface
	\param flags are the flags that will be set
*/
void PatchKernel::setInterfaceAlterationFlags(long id, AlterationFlags flags)
{
	setElementAlterationFlags(id, flags, &m_alteredInterfaces);
}

/*!
	Unset the specified alteration flags for all the altered interfaces.

	\param flags are the flags that will be unset
*/
void PatchKernel::unsetInterfaceAlterationFlags(AlterationFlags flags)
{
	unsetElementAlterationFlags(flags, &m_alteredInterfaces);
}

/*!
	Unset the specified alteration flags for the given interface.

	\param id is the id of the interface
	\param flags are the flags that will be unset
*/
void PatchKernel::unsetInterfaceAlterationFlags(long id, AlterationFlags flags)
{
	unsetElementAlterationFlags(id, flags, &m_alteredInterfaces);
}

/*!
	Test the specified alteration flags for the given element.

	\param id is the id of the element
	\param flags are the flags that will be tested
	\param flagsStorage is the container of the elements' flags
	\result Return true if the flags are set, false otherwise.
*/
bool PatchKernel::testElementAlterationFlags(long id, AlterationFlags flags, const AlterationFlagsStorage &flagsStorage) const
{
	auto storedFlagsItr = flagsStorage.find(id);
	if (storedFlagsItr != flagsStorage.end()) {
		return testAlterationFlags(storedFlagsItr->second, flags);
	} else {
		return testAlterationFlags(FLAG_NONE, flags);
	}
}

/*!
	Test if the requested alteration flags are among the available flags.

	\param availableFlags are the available flags
	\param requestedFlags are the requested flags
	\result Return true if the flags are set, false otherwise.
*/
bool PatchKernel::testAlterationFlags(AlterationFlags availableFlags, AlterationFlags requestedFlags) const
{
	return ((availableFlags & requestedFlags) == requestedFlags);
}

/*!
	Get the alteration flags of the specified element.

	\param id is the id of the element
	\param flagsStorage is the container of the elements' flags
	\result The alteration flags of the element.
*/
PatchKernel::AlterationFlags PatchKernel::getElementAlterationFlags(long id, const AlterationFlagsStorage &flagsStorage) const
{
	auto storedFlagsItr = flagsStorage.find(id);
	if (storedFlagsItr != flagsStorage.end()) {
		return storedFlagsItr->second;
	} else {
		return FLAG_NONE;
	}
}

/*!
	Reset the alteration flags of the specified element.

	\param id is the id of the element
	\param flags are the flags that will be set
	\param flagsStorage is the container of the elements' flags
*/
void PatchKernel::resetElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const
{
	if (flags == FLAG_NONE) {
		flagsStorage->erase(id);
	} else {
		(*flagsStorage)[id] = flags;
	}
}

/*!
	Set the specified alteration flags of the specified element.

	\param id is the id of the element
	\param flags are the flags that will be set
	\param flagsStorage is the container of the elements' flags
*/
void PatchKernel::setElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const
{
	if (flags == FLAG_NONE) {
		return;
	}

	auto storedFlagsItr = flagsStorage->find(id);
	if (storedFlagsItr != flagsStorage->end()) {
		storedFlagsItr->second |= flags;
	} else {
		flagsStorage->insert({id, flags});
	}
}

/*!
	Unset the specified alteration flags for all the altered elements.

	\param flags are the flags that will be unset
	\param flagsStorage is the container of the elements' flags
*/
void PatchKernel::unsetElementAlterationFlags(AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const
{
	if (flags == FLAG_NONE) {
		return;
	}

	for (auto storedFlagsItr = flagsStorage->begin(); storedFlagsItr != flagsStorage->end();) {
		storedFlagsItr->second &= ~flags;
		if (storedFlagsItr->second == FLAG_NONE) {
			storedFlagsItr = flagsStorage->erase(storedFlagsItr);
		} else {
			++storedFlagsItr;
		}
	}
}

/*!
	Unset the specified alteration flags for the specified element.

	\param id is the id of the element
	\param flags are the flags that will be unset
	\param flagsStorage is the container of the elements' flags
*/
void PatchKernel::unsetElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const
{
	if (flags == FLAG_NONE) {
		return;
	}

	auto storedFlagsItr = flagsStorage->find(id);
	if (storedFlagsItr == flagsStorage->end()) {
		return;
	}

	storedFlagsItr->second &= ~flags;
	if (storedFlagsItr->second == FLAG_NONE) {
		flagsStorage->erase(storedFlagsItr);
	}
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
}

/*!
	Gets the previously stored patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void PatchKernel::getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const
{
	getBoundingBox(false, minPoint, maxPoint);
}

/*!
	Gets the previously stored patch bounding box.

	\param global if set to true, the bounding box will be evaluated globally
	across all the partitions
	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void PatchKernel::getBoundingBox(bool global, std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const
{
	minPoint = m_boxMinPoint;
	maxPoint = m_boxMaxPoint;

#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
		MPI_Comm communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, minPoint.data(), 3, MPI_DOUBLE, MPI_MIN, communicator);
		MPI_Allreduce(MPI_IN_PLACE, maxPoint.data(), 3, MPI_DOUBLE, MPI_MAX, communicator);
	}
#else
	BITPIT_UNUSED(global);
#endif
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

	\param global if set to true, the dirty status will be evaluated globally
	across all the partitions
	\result Returns true if the bounding box is dirty, false otherwise.
*/
bool PatchKernel::isBoundingBoxDirty(bool global) const
{
	bool isDirty = m_boxDirty;
#if BITPIT_ENABLE_MPI==1
	if (global && isPartitioned()) {
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

	\param forcedUpdated if set to true, bounding box information will be
	updated also if they are not marked as dirty
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

	for (size_t k = 0; k < point.size(); ++k) {
		double value = point[k];

		// Update maximum value
		if (utils::DoubleFloatingEqual()(value, m_boxMaxPoint[k], getTol())) {
			++m_boxMaxCounter[k];
		} else if (value > m_boxMaxPoint[k]) {
			m_boxMaxPoint[k]   = value;
			m_boxMaxCounter[k] = 1;
		}

		// Update minimum value
		if (utils::DoubleFloatingEqual()(value, m_boxMinPoint[k], getTol())) {
			++m_boxMinCounter[k];
		} else if (value < m_boxMinPoint[k]) {
			m_boxMinPoint[k]   = value;
			m_boxMinCounter[k] = 1;
		}
	}
}

/*!
	Update the bounding removing the specified point.

	The bounding box is not updated if it's set as frozen, or if it's in a
	dirty state.

	\param point is the point that will be removed from to the bounding box
*/
void PatchKernel::removePointFromBoundingBox(const std::array<double, 3> &point)
{
	if (isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		return;
	}

	double tolerance = getTol();
	for (size_t k = 0; k < point.size(); ++k) {
		double value = point[k];

		// Check if maximum value is still valid
		if (utils::DoubleFloatingEqual()(value, m_boxMaxPoint[k], tolerance)) {
			--m_boxMaxCounter[k];
			if (m_boxMaxCounter[k] == 0) {
				setBoundingBoxDirty(true);
				return;
			}
		} else if (value > m_boxMaxPoint[k]) {
			assert(false && "Bounding box is in inconsistent state.");
			setBoundingBoxDirty(true);
			return;
		}

		// Update minimum value
		if (utils::DoubleFloatingEqual()(value, m_boxMinPoint[k], tolerance)) {
			--m_boxMinCounter[k];
			if (m_boxMinCounter[k] == 0) {
				setBoundingBoxDirty(true);
				return;
			}
		} else if (value < m_boxMinPoint[k]) {
			assert(false && "Bounding box is in inconsistent state.");
			setBoundingBoxDirty(true);
			return;
		}
	}
}

/*!
	Get the coordinates of the specified element

	\param element is the element
	\result The coordinates of the element.
*/
ConstProxyVector<std::array<double, 3>> PatchKernel::getElementVertexCoordinates(const Element &element) const
{
	// Get the vertices ids
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	const int nElementVertices = elementVertexIds.size();

	// Build the proxy vector with the coordinates
	ConstProxyVector<std::array<double, 3>> vertexCoordinates(ConstProxyVector<std::array<double, 3>>::INTERNAL_STORAGE, nElementVertices);
	ConstProxyVector<std::array<double, 3>>::storage_pointer storage = vertexCoordinates.storedData();
	getVertexCoords(elementVertexIds.size(), elementVertexIds.data(), storage);

	return vertexCoordinates;
}

/*!
	Get vertex coordinates of the specified element

	\param element is the element
	\param[out] coordinates on output will contain the vertex coordinates
*/
void PatchKernel::getElementVertexCoordinates(const Element &element, std::unique_ptr<std::array<double, 3>[]> *coordinates) const
{
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	getVertexCoords(elementVertexIds.size(), elementVertexIds.data(), coordinates);
}

/*!
	Get vertex coordinates of the specified element

	\param element is the element
	\param[out] coordinates on output will contain the vertex coordinates, it
	is up to the caller to ensure that the storage has enough space for all
	the vertex coordinates
*/
void PatchKernel::getElementVertexCoordinates(const Element &element, std::array<double, 3> *coordinates) const
{
	ConstProxyVector<long> elementVertexIds = element.getVertexIds();
	getVertexCoords(elementVertexIds.size(), elementVertexIds.data(), coordinates);
}

/*!
    Group vertices on regular bins.

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
    double dx = std::max(1.0e-12, m_boxMaxPoint[0] - m_boxMinPoint[0]) / ((double) nBins);
    double dy = std::max(1.0e-12, m_boxMaxPoint[1] - m_boxMinPoint[1]) / ((double) nBins);
    double dz = std::max(1.0e-12, m_boxMaxPoint[2] - m_boxMinPoint[2]) / ((double) nBins);

    // Identify bins of vertices
    std::unordered_map<long, std::vector<long>> bins;
    for (const Vertex &vertex : vertices) {
        const std::array<double, 3> &coordinates = vertex.getCoords();

        int i = std::min(nBins - 1L, (long) ((coordinates[0] - m_boxMinPoint[0]) / dx));
        int j = std::min(nBins - 1L, (long) ((coordinates[1] - m_boxMinPoint[1]) / dy));
        int k = std::min(nBins - 1L, (long) ((coordinates[2] - m_boxMinPoint[2]) / dz));

        long binId = nBins * nBins * k + nBins * j + i;
        bins[binId].emplace_back(vertex.getId());
    }

    return bins;
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
*/
void PatchKernel::translate(const std::array<double, 3> &translation)
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
	Rotates the patch.

	\param[in] n0 is a first point on the rotation axis
	\param[in] n1 is a second point on the rotation axis
	\param[in] angle is the rotation angle, expressed in radiants and positive
	for counterclockwise rotations
*/
void PatchKernel::rotate(const std::array<double, 3> &n0, const std::array<double, 3> &n1, double angle)
{
	// Translate the patch
	for (auto &vertex : m_vertices) {
		vertex.rotate(n0, n1, angle);
	}

	// Update the bounding box
	if (!isBoundingBoxFrozen() && !isBoundingBoxDirty()) {
		// Save current bounding box
		std::array<std::array<double, 3>, 3> originalBox = {{m_boxMinPoint, m_boxMaxPoint}};

		// Clear bounding box
		clearBoundingBox();

		// Unset the dirty flag in order to be able to update the bounding box
		setBoundingBoxDirty(false);

		// Evaluate rotated bounding box
		for (int i = 0; i < 2; ++i) {
			double xCorner = originalBox[i][0];
			for (int j = 0; j < 2; ++j) {
				double yCorner = originalBox[j][1];
				for (int k = 0; k < 2; ++k) {
					double zCorner = originalBox[k][2];

					std::array<double, 3> corner = {{xCorner, yCorner, zCorner}};
					corner = CGElem::rotatePoint(corner, n0, n1, angle);
					addPointToBoundingBox(corner);
				}
			}
		}
	}
}

/*!
	Rotates the patch.

	\param[in] n0x is the x-component of a first point on the rotation axis
	\param[in] n0y is the y-component of a first point on the rotation axis
	\param[in] n0z is the z-component of a first point on the rotation axis
	\param[in] n1x is the x-component of a second point on the rotation axis
	\param[in] n1y is the y-component of a second point on the rotation axis
	\param[in] n1z is the z-component of a second point on the rotation axis
	\param[in] angle is the rotation angle, expressed in radiants.
	Counterclockwise rotations are considered positive.
*/
void PatchKernel::rotate(double n0x, double n0y, double n0z, double n1x,
						 double n1y, double n1z, double angle)
{
	rotate({{n0x, n0y, n0z}}, {{n1x, n1y, n1z}}, angle);
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor vector
*/
void PatchKernel::scale(const std::array<double, 3> &scaling)
{
	scale(scaling, m_boxMinPoint);
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor vector
	\param[in] center is the center of the scaling
*/
void PatchKernel::scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center)
{
	// Scale the patch
	for (auto &vertex : m_vertices) {
		vertex.scale(scaling, center);
	}

	// Update the bounding box
	if (!isBoundingBoxFrozen() || isBoundingBoxDirty()) {
		for (int k = 0; k < 3; ++k) {
			m_boxMinPoint[k] = center[k] + scaling[k] * (m_boxMinPoint[k] - center[k]);
			m_boxMaxPoint[k] = center[k] + scaling[k] * (m_boxMaxPoint[k] - center[k]);
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
	scale({{scaling, scaling, scaling}}, m_boxMinPoint);
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor
	\param[in] center is the center of the scaling
*/
void PatchKernel::scale(double scaling, const std::array<double, 3> &center)
{
	scale({{scaling, scaling, scaling}}, center);
}

/*!
	Scales the patch.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sz scaling factor along z direction
*/
void PatchKernel::scale(double sx, double sy, double sz)
{
	scale({{sx, sy, sz}}, m_boxMinPoint);
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sz scaling factor along z direction
	\param[in] center is the center of the scaling
*/
void PatchKernel::scale(double sx, double sy, double sz, const std::array<double, 3> &center)
{
	scale({{sx, sy, sz}}, center);
}

/*!
	Sets the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void PatchKernel::setTol(double tolerance)
{
	_setTol(tolerance);

	m_toleranceCustom = true;
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
	Gets the tolerance for geometrical checks.

	\result The tolerance for geometrical checks.
*/
double PatchKernel::getTol() const
{
	return m_tolerance;
}

/*!
	Resets the tolerance for geometrical checks.
*/
void PatchKernel::resetTol()
{
	_resetTol();

	m_toleranceCustom = false;
}

/*!
	Internal function to reset the tolerance for the geometrical checks.

	If a derived patch re-implements this function, it's up to the derived
	patch to initialize the tolerance of the newly created patches. Since the
	tolerance is initialized in the constructor, PatchKernel can only reset
	the tolerance using its base method. Derived classes need to explicitly
	initialize the tolerance calling the method they have re-implemented.
*/
void PatchKernel::_resetTol()
{
	const double DEFAULT_TOLERANCE = 1e-14;

	m_tolerance = DEFAULT_TOLERANCE;
}

/*!
	Checks if the tolerance for the geometrical checks has been customized
	by the user.

	\result True if the tolerance was customized by the user, false otherwise.
*/
bool PatchKernel::isTolCustomized() const
{
	return m_toleranceCustom;
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
	envelope.reserveVertices(envelope.getVertexCount() + countBorderVertices());
	envelope.reserveCells(envelope.getCellCount() + countBorderFaces());

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
	out << indent<< "Vertices --------------------------------"     << std::endl;
	out << indent<< "  # vertices        " << getVertexCount()      << std::endl;
	out << indent<< "  # orphan vertices " << countOrphanVertices() << std::endl;
	out << indent<< "  # border vertices " << countBorderVertices()   << std::endl;
        //out << indent<< "  # free vertices   " << countDoubleVertices()   << std::endl;

	// ====================================================================== //
	// FACE STATS                                                             //
	// ====================================================================== //
	out << indent<< "Faces -----------------------------------"     << std::endl;
	out << indent<< "  # faces           " << countFaces()          << std::endl;
	out << indent<< "  # border faces    " << countBorderFaces()    << std::endl;

	// ====================================================================== //
	// CELLS STATS                                                            //
	// ====================================================================== //
	out << indent<< "Cells -----------------------------------"     << std::endl;
	out << indent<< "  # cells           " << getCellCount()        << std::endl;
	out << indent<< "  # orphan cells    " << countOrphanCells()    << std::endl;
	out << indent<< "  # border cells    " << countBorderCells()    << std::endl;
        //out << indent<< "  # free vertices   " << countDoubleCells()   << std::endl;
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
		return CellConstRange(internalCellConstBegin(), internalCellConstEnd());
#endif
	} else {
		return CellConstRange(cellConstEnd(), cellConstEnd());
	}
}

/*!
	Replace the VTK streamer.

	\param original is the original VTK streamer
	\param updated is the updated VTK streamer
*/
void PatchKernel::replaceVTKStreamer(const VTKBaseStreamer *original, VTKBaseStreamer *updated)
{
	// Update the VTK streamer
	//
	// The pointer to VTK streamers are copied, if there are pointer to the
	// original object they have to be replace with a pointer to this object.
	// Geometry fields
	std::vector<std::string> streamedGeomFields;
	streamedGeomFields.reserve(m_vtk.getGeomDataCount());
	for (auto itr = m_vtk.getGeomDataBegin(); itr != m_vtk.getGeomDataEnd(); ++itr) {
		const VTKField &field = *itr;
		if (&field.getStreamer() != original) {
			continue;
		}

		streamedGeomFields.push_back(field.getName());
	}

	for (const std::string &name : streamedGeomFields) {
		const VTKField &field = *(m_vtk.findGeomData(name));
		VTKField updatedField(field);
		updatedField.setStreamer(*updated);

		m_vtk.setGeomData(std::move(updatedField));
	}

	std::vector<std::string> streamedDataFields;
	streamedDataFields.reserve(m_vtk.getDataCount());
	for (auto itr = m_vtk.getDataBegin(); itr != m_vtk.getDataEnd(); ++itr) {
		const VTKField &field = *itr;
		if (&field.getStreamer() != original) {
			continue;
		}

		streamedDataFields.push_back(field.getName());
	}

	for (const std::string &name : streamedDataFields) {
		const VTKField &field = *(m_vtk.findData(name));
		VTKField updatedField(field);
		updatedField.setStreamer(*updated);

		m_vtk.removeData(field.getName());
		m_vtk.addData(std::move(updatedField));
	}
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] name is the name of the data to be written. Either user
 *  data or patch data
 *  @param[in] format is the format that will be used for writing data. Only
 *  the "appended" format is supported. The "appended" format requires an
 *  unformatted binary stream
 */
void PatchKernel::flushData(std::fstream &stream, const std::string &name, VTKFormat format)
{
	if (name == "Points") {
		VertexConstIterator endItr = vertexConstEnd();
		for (VertexConstIterator itr = vertexConstBegin(); itr != endItr; ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap.rawAt(vertexRawId);
			if (vertexVTKId != Vertex::NULL_ID) {
				const Vertex &vertex = m_vertices.rawAt(vertexRawId);
				flushValue(stream, format,vertex.getCoords());
			}
		}
	} else if (name == "offsets") {
		long offset = 0;
		for (const Cell &cell : getVTKCellWriteRange()) {
			offset += cell.getVertexCount();
			flushValue(stream, format,offset);
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

			flushValue(stream, format,(int) VTKType);
		}
	} else if (name == "connectivity") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			const int nCellVertices = cellVertexIds.size();
			for (int k = 0; k < nCellVertices; ++k) {
				long vertexId = cellVertexIds[k];
				long vtkVertexId = m_vtkVertexMap.at(vertexId);
				flushValue(stream, format,vtkVertexId);
			}
		}
	} else if (name == "faces") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				flushValue(stream, format,(long) 0);
			} else {
				std::vector<long> faceStream = cell.getFaceStream();
				Cell::renumberFaceStream(m_vtkVertexMap, &faceStream);
				int faceStreamSize = faceStream.size();
				for (int k = 0; k < faceStreamSize; ++k) {
					flushValue(stream, format,faceStream[k]);
				}
			}
		}
	} else if (name == "faceoffsets") {
		long offset = 0;
		for (const Cell &cell : getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				offset += 1;
			} else {
				offset += cell.getFaceStreamSize();
			}

			flushValue(stream, format,offset);
		}
	} else if (name == "cellIndex") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			flushValue(stream, format,cell.getId());
		}
	} else if (name == "PID") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			flushValue(stream, format,cell.getPID());
		}
	} else if (name == "vertexIndex") {
		VertexConstIterator endItr = vertexConstEnd();
		for (VertexConstIterator itr = vertexConstBegin(); itr != endItr; ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap.rawAt(vertexRawId);
			if (vertexVTKId != Vertex::NULL_ID) {
				long vertexId = itr.getId();
				flushValue(stream, format,vertexId);
			}
		}
#if BITPIT_ENABLE_MPI==1
	} else if (name == "cellGlobalIndex") {
		PatchNumberingInfo numberingInfo(this);
		for (const Cell &cell : getVTKCellWriteRange()) {
			flushValue(stream, format,numberingInfo.getCellGlobalId(cell.getId()));
		}
	} else if (name == "cellRank") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			flushValue(stream, format,getCellOwner(cell.getId()));
		}
	} else if (name == "cellHaloLayer") {
		for (const Cell &cell : getVTKCellWriteRange()) {
			flushValue(stream, format,getCellHaloLayer(cell.getId()));
		}
	} else if (name == "vertexRank") {
		VertexConstIterator endItr = vertexConstEnd();
		for (VertexConstIterator itr = vertexConstBegin(); itr != endItr; ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap.rawAt(vertexRawId);
			if (vertexVTKId != Vertex::NULL_ID) {
				flushValue(stream, format,getVertexOwner(itr.getId()));
			}
		}
#endif
	}
}

/*!
 *  Renumbers vertices consecutively, starting from a given offset.
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

	// Reset index generator
	if (isVertexAutoIndexingEnabled()) {
		createVertexIndexGenerator(true);
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
	if (m_lastInternalCellId >= 0) {
		m_lastInternalCellId = map.at(m_lastInternalCellId);
	}

#if BITPIT_ENABLE_MPI==1
	if (m_firstGhostCellId >= 0) {
		m_firstGhostCellId = map.at(m_firstGhostCellId);
	}
#endif

	// Reset index generator
	if (isCellAutoIndexingEnabled()) {
		createCellIndexGenerator(true);
	}

#if BITPIT_ENABLE_MPI==1
	// Update partitioning information
	if (isPartitioned()) {
		updatePartitioningInfo(true);
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

	// Reset index generator
	if (isInterfaceAutoIndexingEnabled()) {
		createInterfaceIndexGenerator(true);
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
	const int KERNEL_DUMP_VERSION = 12;

	return (KERNEL_DUMP_VERSION + _getDumpVersion());
}

/*!
 *  Write the patch to the specified stream.
 *
 *  Dumping a patch that is not up-to-date is not supported. If the patch is
 *  not up-to-date, it will be automatically updated before dump it.
 *
 *  \param stream is the stream to write to
 *  \result Return true if the patch was successfully dumped, false otherwise.
 */
bool PatchKernel::dump(std::ostream &stream)
{
	// Update the patch
	update();

	// Dump the patch
	const PatchKernel *constPatch = this;
	return constPatch->dump(stream);
}

/*!
 *  Write the patch to the specified stream.
 *
 *  Dumping a patch that is not up-to-date is not supported. If the patch is
 *  not up-to-date and compilation of assertions is enabled, the function will
 *  assert, whereas if compilation of assertions is not enabled, the function
 *  is a no-op.
 *
 *  \param stream is the stream to write to
 *  \result Return true if the patch was successfully dumped, false otherwise.
 */
bool PatchKernel::dump(std::ostream &stream) const
{
	// Dumping a dirty patch is not supported.
	bool dirty = isDirty(true);
	assert(!dirty && "Dumping a patch that is not up-to-date is not supported.");
	if (dirty) {
		return false;
	}

	// Version
	utils::binary::write(stream, getDumpVersion());

	// Generic information
	utils::binary::write(stream, m_id);
	utils::binary::write(stream, m_dimension);
	utils::binary::write(stream, m_vtk.getName());
#if BITPIT_ENABLE_MPI==1
	utils::binary::write(stream, m_haloSize);
#else
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
	utils::binary::write(stream, (int) m_toleranceCustom);
	if (m_toleranceCustom) {
		utils::binary::write(stream, m_tolerance);
	}

	// Index generators
	bool hasVertexAutoIndexing = isVertexAutoIndexingEnabled();
	utils::binary::write(stream, hasVertexAutoIndexing);
	if (hasVertexAutoIndexing) {
		dumpVertexAutoIndexing(stream);
	}

	bool hasInterfaceAutoIndexing = isInterfaceAutoIndexingEnabled();
	utils::binary::write(stream, hasInterfaceAutoIndexing);
	if (hasInterfaceAutoIndexing) {
		dumpInterfaceAutoIndexing(stream);
	}

	bool hasCellAutoIndexing = isCellAutoIndexingEnabled();
	utils::binary::write(stream, hasCellAutoIndexing);
	if (hasCellAutoIndexing) {
		dumpCellAutoIndexing(stream);
	}

	// The patch has been dumped successfully
	return true;
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
		setId(id);
	}

	// Dimension
	int dimension;
	utils::binary::read(stream, dimension);
	setDimension(dimension);

	// Name
	std::string name;
	utils::binary::read(stream, name);
	m_vtk.setName(name);

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

	// Interfaces build strategy
	utils::binary::read(stream, m_interfacesBuildStrategy);

	// VTK data
	utils::binary::read(stream, m_vtkWriteTarget);

	// Specific restore
	_restore(stream);

#if BITPIT_ENABLE_MPI==1
	// Update partitioning information
	if (isPartitioned()) {
		updatePartitioningInfo(true);
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
	bool hasVertexAutoIndexing;
	utils::binary::read(stream, hasVertexAutoIndexing);
	if (hasVertexAutoIndexing) {
		restoreVertexAutoIndexing(stream);
	} else {
		setVertexAutoIndexing(false);
	}

	bool hasInterfaceAutoIndexing;
	utils::binary::read(stream, hasInterfaceAutoIndexing);
	if (hasInterfaceAutoIndexing) {
		restoreInterfaceAutoIndexing(stream);
	} else {
		setInterfaceAutoIndexing(false);
	}

	bool hasCellAutoIndexing;
	utils::binary::read(stream, hasCellAutoIndexing);
	if (hasCellAutoIndexing) {
		restoreCellAutoIndexing(stream);
	} else {
		setCellAutoIndexing(false);
	}
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
