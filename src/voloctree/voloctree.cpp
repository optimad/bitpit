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

#include <cassert>
#include <cmath>

#include "bitpit_common.hpp"

#include "logger.hpp"
#include "voloctree.hpp"

namespace bitpit {

/*!
	\class VolOctree
	\ingroup volumepatches

	\brief The VolOctree defines a Octree patch.

	VolOctree defines a Octree patch.
*/

/*!
	\struct VolOctree::OctantInfo
	\ingroup volumepatches

	OctantInfo allows to uniquely identify an octant of the underlying
	octree.
*/

/*!
	\struct VolOctree::OctantInfoHasher
	\ingroup volumepatches

	OctantInfoHasher allows to generate a hash for the OctantInfo structure.
*/

#if BITPIT_ENABLE_MPI==1
/*!
	Creates an uninitialized partitioned patch.

	Cells will be initialized the cells only on the process identified by the
	rank zero in the communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolOctree::VolOctree(MPI_Comm communicator, std::size_t haloSize)
	: VolumeKernel(communicator, haloSize, false)
#else
/*!
	Creates an uninitialized serial patch.
*/
VolOctree::VolOctree()
	: VolumeKernel(false)
#endif
{
	// Create the tree
#if BITPIT_ENABLE_MPI==1
	m_tree = std::unique_ptr<PabloUniform>(new PabloUniform(PabloUniform::DEFAULT_LOG_FILE, communicator));
#else
	m_tree = std::unique_ptr<PabloUniform>(new PabloUniform(PabloUniform::DEFAULT_LOG_FILE));
#endif

	// Initialize the tree
#if BITPIT_ENABLE_MPI==1
	initializeTree(nullptr, haloSize);
#else
	initializeTree(nullptr);
#endif

	// Initialize the patch
	initialize();

	// Reset
	//
	// The function that resets the patch is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__reset(false);
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	Cells will be initialized the cells only on the process identified by the
	rank zero in the communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolOctree::VolOctree(int dimension, const std::array<double, 3> &origin, double length, double dh, MPI_Comm communicator, std::size_t haloSize)
	: VolOctree(PatchManager::AUTOMATIC_ID, dimension, origin, length, dh, communicator, haloSize)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
*/
VolOctree::VolOctree(int dimension, const std::array<double, 3> &origin, double length, double dh)
	: VolOctree(PatchManager::AUTOMATIC_ID, dimension, origin, length, dh)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	Cells will be initialized the cells only on the process identified by the
	rank zero in the communicator.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolOctree::VolOctree(int id, int dimension, const std::array<double, 3> &origin, double length, double dh, MPI_Comm communicator, std::size_t haloSize)
	: VolumeKernel(id, dimension, communicator, haloSize, false)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
*/
VolOctree::VolOctree(int id, int dimension, const std::array<double, 3> &origin, double length, double dh)
	: VolumeKernel(id, dimension, false)
#endif
{
	// Create the tree
#if BITPIT_ENABLE_MPI==1
	m_tree = std::unique_ptr<PabloUniform>(
	    new PabloUniform(origin[0], origin[1], origin[2], length, dimension,
	                     PabloUniform::DEFAULT_LOG_FILE, communicator));
#else
	m_tree = std::unique_ptr<PabloUniform>(
	    new PabloUniform(origin[0], origin[1], origin[2], length, dimension,
	                     PabloUniform::DEFAULT_LOG_FILE));
#endif

	// Initialize the tree
#if BITPIT_ENABLE_MPI==1
	initializeTree(nullptr, haloSize);
#else
	initializeTree(nullptr);
#endif

	// Initialize the patch
	initialize();

	// Reset
	//
	// The function that resets the patch is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__reset(false);

	// Set the dimension
	//
	// The function that sets the dimension is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__setDimension(dimension);

	// Set the bounding
	setBoundingBox();

	// Initialize refinement markers
	if (m_tree->getNumOctants() > 0) {
		double initial_level = ceil(log2(std::max(1., length / dh)));
		m_tree->setMarker((uint32_t) 0, initial_level);
	}
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
VolOctree::VolOctree(std::istream &stream, MPI_Comm communicator, std::size_t haloSize)
	: VolumeKernel(communicator, haloSize, false)
#else
/*!
	Creates a patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
VolOctree::VolOctree(std::istream &stream)
	: VolumeKernel(false)
#endif
{
	// Initialize the tree
#if BITPIT_ENABLE_MPI==1
	initializeTree(nullptr, haloSize);
#else
	initializeTree(nullptr);
#endif

	// Initialize the patch
	initialize();

	// Restore the patch
	restore(stream);
}

/*!
	Creates a patch.
*/
#if BITPIT_ENABLE_MPI==1
/*!
	If the tree provides a valid communicator the patch will be considered
	partitioned, otherwise the patch will be serial.
*/
#endif
/*!
	\param tree is the tree that will be used
	\param adopter is a pointer to the tree adopter
*/
VolOctree::VolOctree(std::unique_ptr<PabloUniform> &&tree, std::unique_ptr<PabloUniform> *adopter)
	: VolOctree(PatchManager::AUTOMATIC_ID, std::move(tree), adopter)
{
}

/*!
	Creates a paritioned patch.

	Cells will be initialized the cells only on the process identified by the
	rank zero in the communicator.
*/
#if BITPIT_ENABLE_MPI==1
/*!
	If the tree provides a valid communicator the patch will be considered
	partitioned, otherwise the patch will be serial.
*/
#endif
/*!
	\param id is the id that will be assigned to the patch
	\param tree is the tree that will be used
	\param adopter is a pointer to the tree adopter
*/
VolOctree::VolOctree(int id, std::unique_ptr<PabloUniform> &&tree, std::unique_ptr<PabloUniform> *adopter)
#if BITPIT_ENABLE_MPI==1
	: VolumeKernel(id, tree->getDim(), tree->getComm(), tree->getNofGhostLayers(), false)
#else
	: VolumeKernel(id, tree->getDim(), false)
#endif
{
	// Associate the tree
	assert(tree);
	m_tree.swap(tree);

	// Initialize the tree
#if BITPIT_ENABLE_MPI==1
	initializeTree(adopter, m_tree->getNofGhostLayers());
#else
	initializeTree(adopter);
#endif

	// Initialize the patch
	initialize();

	// Reset
	//
	// The function that resets the patch is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__reset(false);

	// Set the dimension
	//
	// The function that sets the dimension is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__setDimension(m_tree->getDim());
}

/*!
	Copy constructor.

	\param other is another patch whose content is copied in this element
*/
VolOctree::VolOctree(const VolOctree &other)
	: VolumeKernel(other)
{
	m_octantLocalFacesOnVertex = other.m_octantLocalFacesOnVertex;
	m_octantLocalFacesOnEdge   = other.m_octantLocalFacesOnEdge;
	m_octantLocalEdgesOnVertex = other.m_octantLocalEdgesOnVertex;

	m_cellTypeInfo      = other.m_cellTypeInfo;
	m_interfaceTypeInfo = other.m_interfaceTypeInfo;

	m_cellToOctant = other.m_cellToOctant;
	m_cellToGhost  = other.m_cellToGhost;
	m_octantToCell = other.m_octantToCell;
	m_ghostToCell  = other.m_ghostToCell;

	if (other.m_tree) {
		m_tree = std::unique_ptr<PabloUniform>(new PabloUniform(*(other.m_tree)));
	}

	m_treeAdopter = other.m_treeAdopter;
}

/*!
	Copy assignment.

	\param other is another patch whose content is copied in this element
*/
VolOctree & VolOctree::operator=(const VolOctree &other)
{
	VolOctree temporary(other);
	std::swap(temporary, *this);

	return *this;
}

/*!
	Destructor.
*/
VolOctree::~VolOctree()
{
	if (m_treeAdopter) {
		m_treeAdopter->swap(m_tree);
	}
}

/*!
	Creates a clone of the pach.

	\result A clone of the pach.
*/
std::unique_ptr<PatchKernel> VolOctree::clone() const
{
	return std::unique_ptr<VolOctree>(new VolOctree(*this));
}

/*!
	Internal function to reset the patch.
*/
void VolOctree::reset()
{
	// Reset the patch kernel
	VolumeKernel::reset();

	// Reset the current patch
	__reset(true);
}

/*!
	Reset the patch.
*/
void VolOctree::__reset(bool resetTree)
{
	// Reset the tree
	if (resetTree) {
		m_tree->reset();
		m_tree->setOrigin(std::array<double, 3>{{0., 0., 0.}});
		m_tree->setL(1.);
	}

	// Reset cell-to-octants maps
	m_cellToOctant.clear();
	m_octantToCell.clear();

	m_cellToGhost.clear();
	m_ghostToCell.clear();
}

/*!
	Initialize the tree and prepare it to be used with the patch.

	\param adopter is a pointer to the tree adopter
*/
#if BITPIT_ENABLE_MPI==1
/*!
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void VolOctree::initializeTree(std::unique_ptr<PabloUniform> *adopter, std::size_t haloSize)
#else
void VolOctree::initializeTree(std::unique_ptr<PabloUniform> *adopter)
#endif
{
	// Initialize logger
#ifdef NDEBUG
	m_tree->getLog().setVerbosities(log::WARNING);
#endif

#if BITPIT_ENABLE_MPI==1
	// Initialize partitioning
    if (isCommunicatorSet()) {
		if (!m_tree->getParallel()) {
			// Initialize halo size
			if (haloSize != m_tree->getNofGhostLayers()) {
				initializeTreeHaloSize(haloSize);
			}

			// Initialize partitioning
			initializeTreePartitioning();
		} else {
			// Check if the current halo size is equal to the requested one
			if (haloSize != m_tree->getNofGhostLayers()) {
				throw std::runtime_error ("Unable to set the requested halo size.");
			}
		}
	}
#endif

	// Set tree adopter
	setTreeAdopter(adopter);
}

/*!
	Initialize the data structures of the patch.
*/
void VolOctree::initialize()
{
	log::cout() << ">> Initializing Octree mesh" << std::endl;

	// Reset the cell and interface type info
	m_cellTypeInfo      = nullptr;
	m_interfaceTypeInfo = nullptr;

	// This patch need to be spawn
	setSpawnStatus(SPAWN_NEEDED);

	// This patch supports adaption
	setAdaptionStatus(ADAPTION_CLEAN);

	// Initialize the tolerance
	//
	// Since the patch re-implements the function to reset the tolerance,
	// it has to explicitly initialize the tolerance.
	resetTol();

	// Set the bounding box as frozen
	setBoundingBoxFrozen(true);
}

/*!
	Sets the dimension of the patch.

	\param dimension the dimension of the patch
*/
void VolOctree::setDimension(int dimension)
{
	VolumeKernel::setDimension(dimension);

	__setDimension(dimension);
}

/*!
	Internal function to set the dimension of the patch.

	\param dimension the dimension of the patch
*/
void VolOctree::__setDimension(int dimension)
{
	if (m_tree->getDim() > 0 && dimension != m_tree->getDim()) {
		throw std::runtime_error ("The dimension does not match the dimension of the octree.");
	}

	// Initialize local edges/vertex/faces association
	std::vector<std::vector<int>>().swap(m_octantLocalFacesOnVertex);
	std::vector<std::vector<int>>().swap(m_octantLocalEdgesOnVertex);
	std::vector<std::vector<int>>().swap(m_octantLocalFacesOnEdge);

	if (dimension == 3) {
		m_octantLocalFacesOnVertex.reserve(8);
		m_octantLocalFacesOnVertex.push_back({{0, 2, 4}});
		m_octantLocalFacesOnVertex.push_back({{1, 2, 4}});
		m_octantLocalFacesOnVertex.push_back({{0, 3, 4}});
		m_octantLocalFacesOnVertex.push_back({{1, 3, 4}});
		m_octantLocalFacesOnVertex.push_back({{0, 2, 5}});
		m_octantLocalFacesOnVertex.push_back({{1, 2, 5}});
		m_octantLocalFacesOnVertex.push_back({{0, 3, 5}});
		m_octantLocalFacesOnVertex.push_back({{1, 3, 5}});

		m_octantLocalEdgesOnVertex.reserve(8);
		m_octantLocalEdgesOnVertex.push_back({{0, 2,  4}});
		m_octantLocalEdgesOnVertex.push_back({{1, 2,  5}});
		m_octantLocalEdgesOnVertex.push_back({{0, 3,  6}});
		m_octantLocalEdgesOnVertex.push_back({{1, 3,  7}});
		m_octantLocalEdgesOnVertex.push_back({{4, 8, 10}});
		m_octantLocalEdgesOnVertex.push_back({{5, 9, 10}});
		m_octantLocalEdgesOnVertex.push_back({{6, 8, 11}});
		m_octantLocalEdgesOnVertex.push_back({{7, 9, 11}});

		m_octantLocalFacesOnEdge.reserve(12);
		m_octantLocalFacesOnEdge.push_back({{0, 4}});
		m_octantLocalFacesOnEdge.push_back({{1, 4}});
		m_octantLocalFacesOnEdge.push_back({{2, 4}});
		m_octantLocalFacesOnEdge.push_back({{3, 4}});
		m_octantLocalFacesOnEdge.push_back({{0, 2}});
		m_octantLocalFacesOnEdge.push_back({{1, 2}});
		m_octantLocalFacesOnEdge.push_back({{0, 3}});
		m_octantLocalFacesOnEdge.push_back({{1, 3}});
		m_octantLocalFacesOnEdge.push_back({{0, 5}});
		m_octantLocalFacesOnEdge.push_back({{1, 5}});
		m_octantLocalFacesOnEdge.push_back({{2, 5}});
		m_octantLocalFacesOnEdge.push_back({{3, 5}});
	} else {
		m_octantLocalFacesOnVertex.reserve(4);
		m_octantLocalFacesOnVertex.push_back({{0, 2}});
		m_octantLocalFacesOnVertex.push_back({{1, 2}});
		m_octantLocalFacesOnVertex.push_back({{0, 3}});
		m_octantLocalFacesOnVertex.push_back({{1, 3}});
	}

	// Info of the cell type
	ElementType cellType;
	if (dimension == 3) {
		cellType = ElementType::VOXEL;
	} else {
		cellType = ElementType::PIXEL;
	}

	m_cellTypeInfo = &ReferenceElementInfo::getInfo(cellType);

	// Info on the interface type
	ElementType interfaceType;
	if (dimension == 3) {
		interfaceType = ElementType::PIXEL;
	} else {
		interfaceType = ElementType::LINE;
	}

	m_interfaceTypeInfo = &ReferenceElementInfo::getInfo(interfaceType);
}

/*!
	Set the bounding box
 */
void VolOctree::setBoundingBox()
{
	std::array<double, 3> minPoint;
	std::array<double, 3> maxPoint;

	// Get the bounding box from the tree
	m_tree->getBoundingBox(minPoint, maxPoint);

#if BITPIT_ENABLE_MPI==1
	// The tree is only evaluating the bounding box of the internal octants,
	// we need to consider also ghosts cells.
	CellConstIterator endItr = ghostCellConstEnd();
	for (CellConstIterator ghostCellItr = ghostCellConstBegin(); ghostCellItr != endItr; ++ghostCellItr) {
		ConstProxyVector<long> ghostVertexIds = ghostCellItr->getVertexIds();
		int nGhostCellVertices = ghostVertexIds.size();
		for (int i = 0; i < nGhostCellVertices; ++i) {
			const std::array<double, 3> coords = m_vertices[ghostVertexIds[i]].getCoords();
			for (int d = 0; d < 3; ++d) {
				minPoint[d] = std::min(coords[d], minPoint[d]);
				maxPoint[d] = std::max(coords[d], maxPoint[d]);
			}
		}
	}
#endif

	// Set the bounding box
	setBoundingBox(minPoint, maxPoint);
}

/*!
	Simulate the adaption of the specified cell.

	\param id is the id of the cell
	\param marker is the adaption marker of the simulated update
	\param[out] virtualCells are the virtual cells that would be outcome of the
	update
	\param[out] virtualVertices are the vertices of the virtual cells that
	would be outcome of the update
*/
void VolOctree::simulateCellUpdate(long id, adaption::Marker marker, std::vector<Cell> *virtualCells, PiercedVector<Vertex, long> *virtualVertices) const
{
	// Get virtual post-update octants
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	int8_t markerValue;
	switch (marker) {

	case adaption::MARKER_COARSEN:
		markerValue = -1;
		break;

	case adaption::MARKER_REFINE:
		markerValue = 1;
		break;

	default:
		markerValue = 0;

	}

	int nMaxVirtualOctants;
	if (getDimension() == 3) {
		nMaxVirtualOctants = 8;
	} else {
		nMaxVirtualOctants = 4;
	}

	std::vector<Octant> virtualOctants;
	virtualOctants.reserve(nMaxVirtualOctants);
	m_tree->expectedOctantAdapt(octant, markerValue, &virtualOctants);
	std::size_t nVirtualOctants = virtualOctants.size();

	// Create virtual post-update cells
	int nCellVertices = m_cellTypeInfo->nVertices;

	virtualVertices->clear();
	virtualVertices->reserve(nCellVertices * nVirtualOctants);

	virtualCells->clear();
	virtualCells->resize(nVirtualOctants);

	for (std::size_t k = 0; k < nVirtualOctants; ++k) {
		const Octant *virtualOctant = &(virtualOctants[k]);

		std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int i = 0; i < nCellVertices; ++i) {
			long vertexId = k * nCellVertices + i;
			connectStorage[i] = vertexId;
			virtualVertices->emplaceBack(vertexId, vertexId, m_tree->getNode(virtualOctant, i));
		}

		virtualCells->emplace_back(id, m_cellTypeInfo->type, std::move(connectStorage), true, false, false);
	}
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolOctree::evalCellVolume(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getVolume(octant);
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> VolOctree::evalCellCentroid(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getCenter(octant);
}

/*!
	Evaluates the bounding box of the specified cell.

	\param id is the id of the cell
	\param[out] minPoint is the minimum point of the bounding box
	\param[out] maxPoint is the maximum point of the bounding box
*/
void VolOctree::evalCellBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	*minPoint = m_tree->getNode(octant, 0);
	*maxPoint = m_tree->getNode(octant, m_cellTypeInfo->nVertices - 1);
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolOctree::evalCellSize(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getSize(octant);
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double VolOctree::evalInterfaceArea(long id) const
{
	const Interface &interface = getInterface(id);
	long owner = interface.getOwner();

	OctantInfo octantInfo = getCellOctant(owner);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getArea(octant);
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> VolOctree::evalInterfaceNormal(long id) const
{
	const Interface &interface = getInterface(id);
	int ownerFace = interface.getOwnerFace();
	long owner = interface.getOwner();

	OctantInfo octantInfo = getCellOctant(owner);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getNormal(octant, (uint8_t) ownerFace);
}

/*!
	Gets the octant of the cell with the specified id.

	\param id the id of the cell
	\result The octant info of the specified cell
*/
VolOctree::OctantInfo VolOctree::getCellOctant(long id) const
{
	// Search the cell among the internal octants
	auto octantItr = m_cellToOctant.find(id);
	if (octantItr != m_cellToOctant.end()) {
		return OctantInfo(octantItr->second, true);
	}

	// Search the cell among the ghosts
	return OctantInfo(m_cellToGhost.at(id), false);
}

/*!
	\brief Gets a reference to the octree associated with the patch.

	\result A reference to the octree associated to the patch.
*/
PabloUniform & VolOctree::getTree()
{
	return *m_tree;
}

/*!
	\brief Gets constant a reference to the octree associated with the patch.

	\result A constant reference to the octree associated to the patch.
*/
const PabloUniform & VolOctree::getTree() const
{
	return *m_tree;
}

/*!
	Sets the tree adopter, i.e., the pointer that will take the ownership of
	the tree when the patch is desctructed.

	\param adopter is a pointer to the tree adopter
*/
void VolOctree::setTreeAdopter(std::unique_ptr<PabloUniform> *adopter)
{
	m_treeAdopter = adopter;
}

/*!
	Gets the id of the specified octant.

	\param octantInfo the data of the octant
	\result The id of the specified octant
*/
long VolOctree::getOctantId(const OctantInfo &octantInfo) const
{
	std::unordered_map<uint32_t, long>::const_iterator octantItr;
	if (octantInfo.internal) {
		octantItr = m_octantToCell.find(octantInfo.id);
		if (octantItr == m_octantToCell.end()) {
			return Element::NULL_ID;
		}
	} else {
		octantItr = m_ghostToCell.find(octantInfo.id);
		if (octantItr == m_ghostToCell.end()) {
			return Element::NULL_ID;
		}
	}

	return octantItr->second;
}

/*!
	Gets a pointer to the specified octant.

	\param octantInfo the data of the octant
	\result The pointer of the specified octant
*/
Octant * VolOctree::getOctantPointer(const OctantInfo &octantInfo)
{
	Octant *octant;
	if (octantInfo.internal) {
		octant = m_tree->getOctant(octantInfo.id);
	} else {
		octant = m_tree->getGhostOctant(octantInfo.id);
	}

	return octant;
}

/*!
    Gets a constant pointer to the specified octant.

    \param octantInfo the data of the octant
    \result The constant pointer of the specified octant
*/
const Octant * VolOctree::getOctantPointer(const OctantInfo &octantInfo) const
{
    const Octant *octant;
    if (octantInfo.internal) {
        octant = m_tree->getOctant(octantInfo.id);
    } else {
        octant = m_tree->getGhostOctant(octantInfo.id);
    }

    return octant;
}

/*!
	Evaluates a unique hash for the octant.

	\param octantInfo the data of the octant
	\result A unique hash for the octant.
*/
VolOctree::OctantHash VolOctree::evaluateOctantHash(const OctantInfo &octantInfo)
{
	uint8_t level   = m_tree->getLevel(octantInfo.id);
	uint64_t morton = m_tree->getMorton(octantInfo.id);

	OctantHash octantHash;
	octantHash |= morton;
	octantHash <<= 8;
	octantHash |= level;

	return octantHash;
}

/*!
	Gets the refinement level of the cell with the specified id.

	\param id is the id of the cell
	\result The refinement level of the specified cell.
*/
int VolOctree::getCellLevel(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant* octant = getOctantPointer(octantInfo);

	return octant->getLevel();
}

/*!
	Get the local index of the family-splitting node for the specified cell.

	\param id is the id of the cell
	\result The local index of the family-splitting node for the specified cell.
*/
int VolOctree::getCellFamilySplitLocalVertex(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getFamilySplittingNode(octant);
}

/*!
	Generates the patch.

	\param trackSpawn if set to true the changes to the patch will be tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
std::vector<adaption::Info> VolOctree::_spawn(bool trackSpawn)
{
	std::vector<adaption::Info> updateInfo;

	// Perform initial import
	if (empty()) {
		m_tree->adapt();
		updateInfo = sync(trackSpawn);
	}

	return updateInfo;
}

/*!
	Prepares the patch for performing the adaption.

	NOTE: only cells are tracked.

	\param trackAdaption if set to true the function will return the changes
	that will be performed in the alter step
	\result If the adaption is tracked, returns a vector of adaption::Info that
	can be used to discover what changes will be performed in the alter step,
	otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> VolOctree::_adaptionPrepare(bool trackAdaption)
{
	BITPIT_UNUSED(trackAdaption);

	if (getSpawnStatus() == SPAWN_NEEDED) {
		throw std::runtime_error ("The initial import has not been performed.");
	}

	// Call pre-adapt routine
	m_tree->preadapt();

	// Track adaption changes
	adaption::InfoCollection adaptionData;
	if (trackAdaption) {
		// Current rank
		int currentRank = -1;
#if BITPIT_ENABLE_MPI==1
		currentRank = getRank();
#endif

		// Track internal and ghosts octants that will be coarsend/refined
		std::vector<uint32_t> treeIds;
		std::vector<int8_t> treeMarkers;
		std::vector<bool> treeGhosts;
		m_tree->getPreMapping(treeIds, treeMarkers, treeGhosts);
		std::size_t nUpdatedOctants = treeIds.size();

		adaption::Info *adaptionInfo = nullptr;
		uint64_t previousFatherMorton = PABLO::INVALID_MORTON;
		for (std::size_t n = 0; n < nUpdatedOctants; ++n) {
			OctantInfo octantInfo(treeIds[n], !treeGhosts[n]);
			Octant *octant = getOctantPointer(octantInfo);

			adaption::Type adaptionType;
			if (treeMarkers[n] > 0) {
				adaptionType = adaption::TYPE_REFINEMENT;
			} else {
				adaptionType = adaption::TYPE_COARSENING;
			}

			bool createAdaption = true;
			if (adaptionType == adaption::TYPE_COARSENING) {
				uint64_t fatherMorton = octant->computeFatherMorton();
				if (fatherMorton == previousFatherMorton) {
					createAdaption = false;
				} else {
					previousFatherMorton = fatherMorton;
				}
			}

			if (createAdaption) {
				std::size_t adaptionInfoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, currentRank);
				adaptionInfo = &(adaptionData[adaptionInfoId]);
			}

			assert(adaptionInfo);
			adaptionInfo->previous.emplace_back();
			long &cellId = adaptionInfo->previous.back();
			cellId = getOctantId(octantInfo);
		}

#if BITPIT_ENABLE_MPI==1
		// Ghost cells will be removed
		if (isPartitioned() && getGhostCellCount() > 0) {
			CellConstIterator beginItr = ghostCellConstBegin();
			CellConstIterator endItr   = ghostCellConstEnd();

			std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_DELETION, adaption::ENTITY_CELL, currentRank);
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
			adaptionInfo.previous.reserve(getGhostCellCount());
			for (CellConstIterator itr = beginItr; itr != endItr; ++itr) {
				adaptionInfo.previous.emplace_back();
				long &deletedGhostCellId = adaptionInfo.previous.back();
				deletedGhostCellId = itr.getId();
			}
		}
#endif
	}

	return adaptionData.dump();
}

/*!
	Alter the patch performing the adpation.

	\param trackAdaption if set to true the function will return the changes
	done to the patch during the adaption
	\result If the adaption is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the adaption, otherwise an
	empty vector will be returned.
*/
std::vector<adaption::Info> VolOctree::_adaptionAlter(bool trackAdaption)
{
	// Updating the tree
	log::cout() << ">> Adapting tree...";

	bool emtpyPatch = empty();
	bool buildMapping = !emtpyPatch;
	bool updated = m_tree->adapt(buildMapping);

	if (!updated && !emtpyPatch) {
		log::cout() << " Already updated" << std::endl;
		return std::vector<adaption::Info>();
	}
	log::cout() << " Done" << std::endl;

	// Sync the patch
	return sync(trackAdaption);
}

/*!
	Cleanup patch data structured after the adaption.
*/
void VolOctree::_adaptionCleanup()
{
	// Nothing to do
}

/*!
	Make the adaption markers set by the user consistent with the internal
	criteria defined by the patch.

	The patch will enforce the 2:1 balancing on the cells that have it enabled
	and will verify the feasability of coarsening. Adaption markers will be
	updated consequently.
*/
void VolOctree::settleAdaptionMarkers()
{
	m_tree->settleMarkers();
}

/*!
	Syncronizes the patch with the underlying octree.

	\param trackChanges if set to true the changes to the patch will be
	tracked
	\result Returns all the changes applied to the patch.
*/
std::vector<adaption::Info> VolOctree::sync(bool trackChanges)
{
	log::cout() << ">> Syncing patch..." << std::endl;

	// Detect if we are in import-from-scratch mode
	//
	// In import-from-scratch mode we start form an empty patch and we need
	// to import all the octants of the tree. If we are importing the tree
	// from scratch there are no cells to delete/renumber.
	bool importFromScratch = empty();

	// Last operation on the tree
	ParaTree::Operation lastTreeOperation = m_tree->getLastOperation();
	if (lastTreeOperation == ParaTree::OP_ADAPT_UNMAPPED && !importFromScratch) {
		throw std::runtime_error ("Unable to sync the patch after an unmapped adaption");
	}

	// Info on the tree
	long nOctants = m_tree->getNumOctants();
	long nPreviousOctants = m_octantToCell.size();

	log::cout() << ">> Number of octants : " << nOctants << std::endl;

	// Info on the tree
	long nGhostsOctants = m_tree->getNumGhosts();
	long nPreviousGhosts = m_ghostToCell.size();

	// Initialize tracking data
	adaption::InfoCollection adaptionData;

	// Current rank
	int currentRank = -1;
#if BITPIT_ENABLE_MPI==1
	currentRank = getRank();
#endif

	// Extract information for transforming the patch
	//
	// If there are no cells in the mesh we need to import all
	// octants.
	log::cout() << ">> Extract information for transforming the patch...";

	std::vector<bool> unmappedOctants(nPreviousOctants, true);
	std::vector<OctantInfo> addedOctants;
	std::vector<RenumberInfo> renumberedOctants;
	std::vector<DeleteInfo> deletedOctants;

	addedOctants.reserve(nOctants + nGhostsOctants);
	if (!importFromScratch) {
		renumberedOctants.reserve(nPreviousOctants + nPreviousGhosts);
		deletedOctants.reserve(nPreviousOctants + nPreviousGhosts);
	}

	uint32_t treeId = 0;
	std::vector<uint32_t> mapper_octantMap;
	std::vector<bool> mapper_ghostFlag;
	std::vector<int> mapper_octantRank;
	while (treeId < (uint32_t) nOctants) {
		// Octant mapping
		mapper_octantMap.clear();
		mapper_ghostFlag.clear();
		mapper_octantRank.clear();
		if (!importFromScratch) {
			m_tree->getMapping(treeId, mapper_octantMap, mapper_ghostFlag, mapper_octantRank);
		}

		// Adaption type
		adaption::Type adaptionType = adaption::TYPE_NONE;
		if (importFromScratch) {
			adaptionType = adaption::TYPE_CREATION;
		} else if (lastTreeOperation == ParaTree::OP_ADAPT_MAPPED) {
			bool isNewR = m_tree->getIsNewR(treeId);
			if (isNewR) {
				adaptionType = adaption::TYPE_REFINEMENT;
			} else {
				bool isNewC = m_tree->getIsNewC(treeId);
				if (isNewC) {
					adaptionType = adaption::TYPE_COARSENING;
				} else if (treeId != mapper_octantMap.front()) {
					adaptionType = adaption::TYPE_RENUMBERING;
				}
			}
#if BITPIT_ENABLE_MPI==1
		} else if (lastTreeOperation == ParaTree::OP_LOADBALANCE || lastTreeOperation == ParaTree::OP_LOADBALANCE_FIRST) {
			if (currentRank != mapper_octantRank.front()) {
				adaptionType = adaption::TYPE_PARTITION_RECV;
			} else if (treeId != mapper_octantMap.front()) {
				adaptionType = adaption::TYPE_RENUMBERING;
			}
#endif
		}

		// If the octant cell has not been modified we can skip to the next
		// octant.
		if (adaptionType == adaption::TYPE_NONE) {
			unmappedOctants[treeId] = false;
			++treeId;
			continue;
		}

		// Re-numbered cells just need to be added to the proper list.
		//
		// Renumbered cells are not tracked, because the re-numbering
		// only happens inside VolOctree.
		if (adaptionType == adaption::TYPE_RENUMBERING) {
			uint32_t previousTreeId = mapper_octantMap.front();
			OctantInfo previousOctantInfo(previousTreeId, !mapper_ghostFlag.front());
			long cellId = getOctantId(previousOctantInfo);
			renumberedOctants.emplace_back(cellId, treeId);
			unmappedOctants[previousTreeId] = false;

			// No more work needed, skip to the next octant
			++treeId;
			continue;
		}

		// Handle other kind of adaption
		//
		// New octants need to be imported into the patch,
		// whereas cells associated to previous octants
		// need to be removed.
		//
		// If the user want to track adaption, adaption
		// data needs to be filled.

		// Current tree ids that will be imported
		long nCurrentTreeIds;
		if (importFromScratch) {
			nCurrentTreeIds = nOctants - treeId;
		} else if (adaptionType == adaption::TYPE_REFINEMENT) {
			nCurrentTreeIds = pow(2, getDimension());
		} else {
			nCurrentTreeIds = 1;
		}

		const long lastCurrentTreeId = treeId + nCurrentTreeIds;
		for (int currentTreeId = treeId; currentTreeId < lastCurrentTreeId; ++currentTreeId) {
			addedOctants.emplace_back(currentTreeId, true);
		}

		// Cells that will be removed
		//
		// Mark the cells associated to previous local octants for deletion.
		if (!importFromScratch) {
			int nPreviousTreeIds = mapper_octantMap.size();
			for (int k = 0; k < nPreviousTreeIds; ++k) {
#if BITPIT_ENABLE_MPI==1
				// Only local cells can be deleted
				if (mapper_octantRank[k] != currentRank) {
					continue;
				}
#endif

				// Ghost octants will be processed later
				if (mapper_ghostFlag[k]) {
					continue;
				}

				// Mark previous octant for deletion
				uint32_t previousTreeId = mapper_octantMap[k];
				OctantInfo previousOctantInfo(previousTreeId, !mapper_ghostFlag[k]);
				long cellId = getOctantId(previousOctantInfo);
				deletedOctants.emplace_back(cellId, adaptionType);

				unmappedOctants[previousTreeId] = false;
			}
		}

		// Adaption tracking
		//
		// The adaption info associated to the octants that has been received
		// from external partitions will contain the current octants sorted by
		// their tree id (we are looping over the octants in that order), this
		// is the same order that will be used on the process that has sent
		// the octants. Since the order is the same, the two processes are able
		// to exchange cell data without any additional extra communication
		// (they already know the list of cells for which data is needed and
		// the order in which these data will be sent).
		if (trackChanges) {
			// Rank assocated to the adaption info
			int rank = currentRank;
#if BITPIT_ENABLE_MPI==1
			if (adaptionType == adaption::TYPE_PARTITION_RECV) {
				rank = mapper_octantRank[0];
			}
#endif

			// Get the adaption info
			std::size_t infoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, rank);
			adaption::Info &adaptionInfo = adaptionData[infoId];

			// Current status
			//
			// We don't know the id of the current status, because those
			// cells are not yet in the mesh. Store the trre id, and
			// make the translation later.
			//
			// WARNING: tree id are uint32_t wherase adaptionInfo stores
			//          id as long.
			adaptionInfo.current.reserve(nCurrentTreeIds);
			auto addedOctantsIter = addedOctants.cend() - nCurrentTreeIds;
			while (addedOctantsIter != addedOctants.cend()) {
				adaptionInfo.current.emplace_back();
				long &adaptionId = adaptionInfo.current.back();
				adaptionId = (*addedOctantsIter).id;

				addedOctantsIter++;
			}

			// Previous cells
			//
			// A coarsening can merge togheter cells of different processes.
			// However, since the coarsening is limited to one level, the
			// previous cells will always be internal or among the ghost of
			// the current process.
			int nPreviousCellIds = mapper_octantMap.size();
			adaptionInfo.previous.reserve(nPreviousCellIds);
			for (int k = 0; k < nPreviousCellIds; ++k) {
				long previousCellId;
#if BITPIT_ENABLE_MPI==1
				if (mapper_octantRank[k] != currentRank) {
					previousCellId = Cell::NULL_ID;
				} else
#endif
				{
					OctantInfo previousOctantInfo(mapper_octantMap[k], !mapper_ghostFlag[k]);
					previousCellId = getOctantId(previousOctantInfo);
				}

				adaptionInfo.previous.emplace_back();
				long &adaptionId = adaptionInfo.previous.back();
				adaptionId = previousCellId;
			}
		}

		// Incremente tree id
		treeId += nCurrentTreeIds;
	}

	log::cout() << " Done" << std::endl;

#if BITPIT_ENABLE_MPI==1
	// New ghost octants need to be added
	for (uint32_t treeId = 0; treeId < (uint32_t) nGhostsOctants; ++treeId) {
		addedOctants.emplace_back(treeId, false);
	}
#endif

	// Remove octants that are no more in the tree
	if (!importFromScratch) {
#if BITPIT_ENABLE_MPI==1
		// Cells that have been send to other processes need to be removed
		PabloUniform::LoadBalanceRanges loadBalanceRanges = m_tree->getLoadBalanceRanges();
		for (const auto &rankEntry : loadBalanceRanges.sendRanges) {
			int rank = rankEntry.first;
			if (rank == currentRank) {
				continue;
			}

			adaption::Type deletionType;
			if (loadBalanceRanges.sendAction == PabloUniform::LoadBalanceRanges::ACTION_DELETE) {
				deletionType = adaption::TYPE_DELETION;
			} else {
				deletionType = adaption::TYPE_PARTITION_SEND;
			}

			uint32_t beginTreeId = rankEntry.second[0];
			uint32_t endTreeId   = rankEntry.second[1];
			for (uint32_t treeId = beginTreeId; treeId < endTreeId; ++treeId) {
				OctantInfo octantInfo(treeId, true);
				long cellId = getOctantId(octantInfo);
				deletedOctants.emplace_back(cellId, deletionType, rank);
				unmappedOctants[treeId] = false;
			}
		}

		// Previous ghosts cells need to be removed
		if (nPreviousGhosts > 0) {
			for (uint32_t ghostTreeId = 0; ghostTreeId < nPreviousGhosts; ++ghostTreeId) {
				OctantInfo ghostOctantInfo(ghostTreeId, false);
				long ghostCellId = getOctantId(ghostOctantInfo);
				deletedOctants.emplace_back(ghostCellId, adaption::TYPE_DELETION);
			}
		}
#endif

		// Remove unmapped octants
		//
		// A coarsening that merges cells from different processes, can leave, on
		// the processes which own the ghost octants involved in the coarsening,
		// some octants that are not mapped.
		for (uint32_t previousTreeId = 0; previousTreeId < nPreviousOctants; ++previousTreeId) {
			if (unmappedOctants[previousTreeId]) {
				OctantInfo octantInfo = OctantInfo(previousTreeId, true);
				long cellId = getOctantId(octantInfo);
				deletedOctants.emplace_back(cellId, adaption::TYPE_DELETION);
			}
		}
	}

	// Enable advanced editing
	setExpert(true);

	// Renumber cells
	renumberCells(renumberedOctants);

	// Remove deleted cells
	StitchInfo stitchInfo;
	if (deletedOctants.size() > 0) {
		log::cout() << ">> Removing non-existing cells...";

		// Track changes
		//
		// The adaption info associated to the octants that has been sent
		// to external partitions will contain the current octants sorted by
		// their tree id (they were added to the deleted octants list in that
		// order), this is the same order that will be used on the process
		// that has received the octants. Since the order is the same, the two
		// processes are able to exchange cell data without any additional
		// extra communication (they already know the list of cells for which
		// data is needed and the order in which these data will be sent).
		if (trackChanges) {
			// Track deleted cells
			std::unordered_set<long> sendAdaptionInfo;
			std::unordered_set<long> removedInterfaces;
			for (const DeleteInfo &deleteInfo : deletedOctants) {
				// Cell id
				long cellId = deleteInfo.cellId;

				// Adaption tracking
				//
				// Only cells deleted from a real deletion or a partition send
				// needs to be tracked here, the other cells will be tracked
				// where the adaption that deleted the cell is tracked.
				adaption::Type adaptionType = deleteInfo.trigger;
				bool adaptionIsDeletion = (adaptionType == adaption::TYPE_DELETION);
				bool adaptionIsSend = (adaptionType == adaption::TYPE_PARTITION_SEND);

				if (adaptionIsDeletion || adaptionIsSend) {
					int rank = deleteInfo.rank;
					std::size_t adaptionInfoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, rank);
					adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
					adaptionInfo.previous.emplace_back();
					long &deletedId = adaptionInfo.previous.back();
					deletedId = cellId;

					// Keep track of adaption info for the send cells
					if (adaptionIsSend) {
						sendAdaptionInfo.insert(adaptionInfoId);
					}
				}

				// List of deleted interfaces
				const Cell &cell = m_cells.at(cellId);
				long nCellInterfaces = cell.getInterfaceCount();
				const long *interfaces = cell.getInterfaces();
				for (int k = 0; k < nCellInterfaces; ++k) {
					long interfaceId = interfaces[k];
					if (interfaceId >= 0) {
						removedInterfaces.insert(interfaceId);
					}
				}
			}

#if BITPIT_ENABLE_MPI==1
			// Sort sent cells
			//
			// We cannot use native functions to evaluate the position of the
			// cells because the octants associated to the cells no longer
			// exist on the octree. The cells are still there, therefore we
			// can evaluate the cell positions using generic patch functions.
			for (int adaptionInfoId : sendAdaptionInfo) {
				adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
				std::sort(adaptionInfo.previous.begin(), adaptionInfo.previous.end(), CellPositionLess(*this, false));
			}
#endif

			// Adaption info for the deleted interfaces
			std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_DELETION, adaption::ENTITY_INTERFACE, currentRank);
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
			adaptionInfo.previous.reserve(removedInterfaces.size());
			for (long interfaceId : removedInterfaces) {
				adaptionInfo.previous.emplace_back();
				long &deletedInterfaceId = adaptionInfo.previous.back();
				deletedInterfaceId = interfaceId;
			}
		}

		// Delete cells
		stitchInfo = deleteCells(deletedOctants);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Cells removed: " <<  deletedOctants.size() << std::endl;
	}

	std::vector<DeleteInfo>().swap(deletedOctants);

	// Import added cells
	std::vector<long> createdCells;
	if (addedOctants.size() > 0) {
		log::cout() << ">> Importing new octants...";

		createdCells = importCells(addedOctants, stitchInfo);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Octants imported: " <<  addedOctants.size() << std::endl;
	}

	StitchInfo().swap(stitchInfo);

	// Disable advanced editing
	setExpert(false);

	// Track mesh adaption
	if (trackChanges) {
		// Complete mesh adaption info for the cells
		for (auto &adaptionInfo : adaptionData.data()) {
			if (adaptionInfo.entity != adaption::ENTITY_CELL) {
				continue;
			}

			// Map ids of the added cells
			int nCurrentIds = adaptionInfo.current.size();
			for (int k = 0; k < nCurrentIds; ++k) {
				long cellId = m_octantToCell.at(adaptionInfo.current[k]);
				adaptionInfo.current[k] = cellId;
			}

#if BITPIT_ENABLE_MPI==1
			// Sort received cells
			//
			// To match the sorting done on the procesor that sent the cells,
			// we don't use the native functions to evaluate the position of
			// the cells.
			adaption::Type adaptionType = adaptionInfo.type;
			if (adaptionType == adaption::TYPE_PARTITION_RECV) {
				std::sort(adaptionInfo.current.begin(), adaptionInfo.current.end(), CellPositionLess(*this, false));
			}
#endif
		}

#if BITPIT_ENABLE_MPI==1
		// Track created ghosts cells
		if (nGhostsOctants > 0) {
			std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_CREATION, adaption::ENTITY_CELL, currentRank);
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];

			adaptionInfo.current.reserve(nGhostsOctants);
			auto cellIterator = m_cellToGhost.cbegin();
			while (cellIterator != m_cellToGhost.cend()) {
				adaptionInfo.current.emplace_back();
				long &adaptionId = adaptionInfo.current.back();
				adaptionId = cellIterator->first;

				cellIterator++;
			}
		}
#endif

		// Track created interfaces
		if (createdCells.size() > 0) {
			// List of unique interfaces that have been created
			std::unordered_set<long> createdInterfaces;
			for (const auto &cellId : createdCells) {
				const Cell &cell = m_cells.at(cellId);
				long nCellInterfaces = cell.getInterfaceCount();
				const long *interfaces = cell.getInterfaces();
				for (int k = 0; k < nCellInterfaces; ++k) {
					long interfaceId = interfaces[k];
					if (interfaceId >= 0) {
						createdInterfaces.insert(interfaceId);
					}
				}
			}

			// Adaption info
			std::size_t infoId = adaptionData.create(adaption::TYPE_CREATION, adaption::ENTITY_INTERFACE, currentRank);
			adaption::Info &adaptionInfo = adaptionData[infoId];
			adaptionInfo.current.reserve(createdInterfaces.size());
			for (long interfaceId : createdInterfaces) {
				adaptionInfo.current.emplace_back();
				long &createdInterfaceId = adaptionInfo.current.back();
				createdInterfaceId = interfaceId;
			}
		}
	}

	// Done
	return adaptionData.dump();
}


/*!
	Delete the specified cells.

	\param deletedOctants contains a list with the information about the
	deleted octants
	\result Returns the stitch information that can  the faces created
	after deleting the octants.
*/
VolOctree::StitchInfo VolOctree::deleteCells(const std::vector<DeleteInfo> &deletedOctants)
{
	// Info of the cells
	int nCellVertices = m_cellTypeInfo->nVertices;
	int nCellFaces    = m_cellTypeInfo->nFaces;

	// Delete the cells
	std::unordered_set<long> deadVertices;

	std::vector<long> deadCells;
	deadCells.reserve(deletedOctants.size());
	for (const DeleteInfo &deleteInfo : deletedOctants) {
		long cellId = deleteInfo.cellId;
		const Cell &cell = m_cells[cellId];

		// List vertices to remove
		//
		// For now, all cell vertices will be listed. Later, the vertex of
		// the dangling faces will be removed from the list.
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		for (int k = 0; k < nCellVertices; ++k) {
			long vertexId = cellVertexIds[k];
			deadVertices.insert(vertexId);
		}

		// Remove patch-tree associations
		const OctantInfo &octantInfo = getCellOctant(cellId);
		uint32_t treeId = octantInfo.id;

		if (cell.isInterior()) {
			m_cellToOctant.erase(cellId);

			auto octantToCellItr = m_octantToCell.find(treeId);
			if (octantToCellItr->second == cellId) {
				m_octantToCell.erase(octantToCellItr);
			}
		} else {
			m_cellToGhost.erase(cellId);

			auto ghostToCellItr = m_ghostToCell.find(treeId);
			if (ghostToCellItr->second == cellId) {
				m_ghostToCell.erase(ghostToCellItr);
			}
		}

		// Cell needs to be removed
		deadCells.push_back(cellId);
	}

	PatchKernel::deleteCells(deadCells);

	// Prune cell adjacencies and interfaces
	//
	// At this stage we cannot fully update adjacencies and interfaces, but
	// we need to remove stale adjacencies and interfaces.
	pruneStaleAdjacencies();

	pruneStaleInterfaces();

	// All the vertices belonging to the dangling cells has to be kept
	//
	// We need to keep all the vertices that lie on a dangling face and belong
	// to cells that were not deleted. These are the vertices of the dangling
	// faces, plus all the vertices that lie on the borders of the dangling
	// faces. Latter vertices arise from three-dimensional configurations where
	// a small cell has been deleted and, among its non-deleted neighburs, there
	// areboth larger cells and cells of the same size. In this case, a vertex
	// of the deleted small cell can lie on one of the edges of a bigger cell
	// and belong also to a smaller cell. If only the bigger cell is dangling,
	// that vertex will lie on a dangling face (it's on the edge of the cell,
	// hence on the border of the face), but it's not one of the vertices of
	// the dangling cells.
	//
	// To identify all vertices that need to be kept, we consider the vertices
	// of the dangling faces and the vertices of the smaller neighbouring faces
	// (i.e., the faces of the smaller neigbours of the dangling cells that lie
	// on the faces of the dangling cells).
	//
	// Morover we need to build a map between the patch numbering and the
	// octree numbering of the vertices of the dangling cells. This map will
	// be used when imprting the octants to stitch the imported octants to
	// the existing cells.
	StitchInfo stitchVertices;
	for (const auto &entry: m_alteredCells) {
		// Consider only dangling faces
		if (!testAlterationFlags(entry.second, FLAG_DANGLING)) {
			continue;
		}

		// Vertices of the cell
		long cellId = entry.first;
		const Cell &cell = m_cells[cellId];
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

		OctantInfo octantInfo = getCellOctant(cellId);
		Octant *octant = getOctantPointer(octantInfo);

		for (int k = 0; k < nCellVertices; ++k) {
			long vertexId = cellVertexIds[k];
			uint64_t vertexTreeKey = m_tree->computeNodePersistentKey(octant, k);
			stitchVertices.insert({vertexTreeKey, vertexId});
			deadVertices.erase(vertexId);
		}

		// Vertices of neighbours
		//
		// We need to select the neighbours smaller than the current cell and
		// keep the vertices of the faces that are shared with the current cell.
		int cellLevel = octant->getLevel();

		for (int face = 0; face < nCellFaces; ++face) {
			int nFaceAdjacencies = cell.getAdjacencyCount(face);
			const long *faceAdjacencies = cell.getAdjacencies(face);
			for (int i = 0; i < nFaceAdjacencies; ++i) {
				int neighId = faceAdjacencies[i];
				OctantInfo neighOctantInfo = getCellOctant(neighId);
				Octant *neighOctant = getOctantPointer(neighOctantInfo);
				int neighLevel = neighOctant->getLevel();
				if (neighLevel <= cellLevel) {
					continue;
				}

				const Cell &neigh = m_cells[neighId];
				int neighFace = findAdjoinNeighFace(cell, face, neigh);
				const int *localNeighFaceConnect = m_cellTypeInfo->faceConnectStorage[neighFace].data();
				ConstProxyVector<long> faceVertexIds = neigh.getFaceVertexIds(neighFace);
				std::size_t nFaceVertexIds = faceVertexIds.size();
				for (std::size_t k = 0; k < nFaceVertexIds; ++k) {
					long vertexId = faceVertexIds[k];
					uint64_t vertexTreeKey = m_tree->computeNodePersistentKey(neighOctant, localNeighFaceConnect[k]);
					stitchVertices.insert({vertexTreeKey, vertexId});
					deadVertices.erase(vertexId);
				}
			}
		}
	}

	// Delete the vertices
	PatchKernel::deleteVertices(deadVertices);

	// Done
	return stitchVertices;
}

/*!
	Renumber the specified cells.

	\param renumberedOctants contains the information about the renumbered
	octants
*/
void VolOctree::renumberCells(const std::vector<RenumberInfo> &renumberedOctants)
{
	// Remove old patch-to-tree and tree-to-patch associations
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		long cellId = renumberInfo.cellId;
		if (!m_cells[cellId].isInterior()) {
			continue;
		}

		const OctantInfo &previousOctantInfo = getCellOctant(cellId);
		uint32_t previousTreeId = previousOctantInfo.id;

		m_octantToCell.erase(previousTreeId);
	}

	// Create new patch-to-tree and tree-to-patch associations
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		long cellId = renumberInfo.cellId;
		if (!m_cells[cellId].isInterior()) {
			continue;
		}

		uint32_t treeId = renumberInfo.newTreeId;

		m_cellToOctant[cellId] = treeId;
		m_octantToCell[treeId] = cellId;
	}
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
	\param stitchInfo is the list of vertices that will be used to stitch the
	the new octants
	\param restoreStream is an optional stream from which the information about
	the restore will be read. If no restore stream is given the the cells will
	be created from scratch
*/
std::vector<long> VolOctree::importCells(const std::vector<OctantInfo> &octantInfoList,
                                         StitchInfo &stitchInfo, std::istream *restoreStream)
{
	// Create the new vertices
	int nCellVertices = m_cellTypeInfo->nVertices;
	for (const OctantInfo &octantInfo : octantInfoList) {
		Octant *octant = getOctantPointer(octantInfo);
		for (int k = 0; k < nCellVertices; ++k) {
			uint64_t vertexTreeKey = m_tree->computeNodePersistentKey(octant, k);
			if (stitchInfo.count(vertexTreeKey) == 0) {
				// Vertex coordinates
				std::array<double, 3> nodeCoords = m_tree->getNode(octant, k);

				// Create the vertex
				long vertexId;
				if (!restoreStream) {
#if BITPIT_ENABLE_MPI==1
					VertexIterator vertexIterator = addVertex(std::move(nodeCoords));
#else
					VertexIterator vertexIterator = addVertex(std::move(nodeCoords));
#endif
					vertexId = vertexIterator.getId();
				} else {
					utils::binary::read(*restoreStream, vertexId);

#if BITPIT_ENABLE_MPI==1
					int rank;
					utils::binary::read(*restoreStream, rank);

					restoreVertex(std::move(nodeCoords), rank, vertexId);
#else
					int dummtRank;
					utils::binary::read(*restoreStream, dummtRank);

					restoreVertex(std::move(nodeCoords), vertexId);
#endif
				}

				// Add the vertex to the stitching info
				stitchInfo[vertexTreeKey] = vertexId;
			}
		}
	}

	// Reserve space for the maps
	long nOctants = m_tree->getNumOctants();
	m_cellToOctant.reserve(nOctants);
	m_octantToCell.reserve(nOctants);

	long nGhostsOctants = m_tree->getNumGhosts();
	m_cellToGhost.reserve(nGhostsOctants);
	m_ghostToCell.reserve(nGhostsOctants);

	// Add the cells
	size_t octantInfoListSize = octantInfoList.size();
	m_cells.reserve(m_cells.size() + octantInfoListSize);

	std::vector<long> createdCells(octantInfoListSize);
	for (size_t i = 0; i < octantInfoListSize; ++i) {
		const OctantInfo &octantInfo = octantInfoList[i];

		// Octant connectivity
		Octant *octant = getOctantPointer(octantInfo);

		// Cell connectivity
		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int k = 0; k < nCellVertices; ++k) {
			uint64_t vertexTreeKey = m_tree->computeNodePersistentKey(octant, k);
			cellConnect[k] = stitchInfo.at(vertexTreeKey);
		}

#if BITPIT_ENABLE_MPI==1
		// Cell owner
		int rank;
		if (octantInfo.internal) {
			rank = getRank();
		} else {
			uint64_t globalTreeId = m_tree->getGhostGlobalIdx(octantInfo.id);
			rank = m_tree->getOwnerRank(globalTreeId);
		}
#endif

		// Add cell
		long cellId;
		if (!restoreStream) {
#if BITPIT_ENABLE_MPI==1
			CellIterator cellIterator = addCell(m_cellTypeInfo->type, std::move(cellConnect), rank);
#else
			CellIterator cellIterator = addCell(m_cellTypeInfo->type, std::move(cellConnect));
#endif
			cellId = cellIterator.getId();
		} else {
			utils::binary::read(*restoreStream, cellId);

#if BITPIT_ENABLE_MPI==1
			restoreCell(m_cellTypeInfo->type, std::move(cellConnect), rank, cellId);
#else
			restoreCell(m_cellTypeInfo->type, std::move(cellConnect), cellId);
#endif
		}

		// Create patch-tree associations
		if (octantInfo.internal) {
			m_cellToOctant.insert({{cellId, octantInfo.id}});
			m_octantToCell.insert({{octantInfo.id, cellId}});
		} else {
			m_cellToGhost.insert({{cellId, octantInfo.id}});
			m_ghostToCell.insert({{octantInfo.id, cellId}});
		}

		// Add the cell to the list of created cells
		createdCells[i] = cellId;
	}

	// Update first ghost and last internal information
	if (restoreStream) {
		updateLastInternalCellId();
#if BITPIT_ENABLE_MPI==1
		updateFirstGhostCellId();
#endif
	}

	// Update cell PIDs
	if (restoreStream) {
		for (Cell &cell : getCells()) {
			int PID;
			utils::binary::read(*restoreStream, PID);
			cell.setPID(PID);
		}
	}

	// Update adjacencies
	updateAdjacencies();

	// Update interfaces
	if (!restoreStream) {
		updateInterfaces();
	}

	// Done
	return createdCells;
}

/*!
	Internal function to update the adjacencies of the patch.

	In addition to the cells whose adjacencies are marked as dirty, also the
	adjacencies of their neighbours will be updated.
*/
void VolOctree::_updateAdjacencies()
{
	// Tree information
	int maxLevel = m_tree->getMaxDepth();

	// Face information
	int nCellFaces = m_cellTypeInfo->nFaces;
	const uint8_t *oppositeFace = m_tree->getOppface();

	// Count cells with dirty adjacencies
	long nDirtyAdjacenciesCells = 0;
	std::vector<std::size_t> nDirtyAdjacenciesCellsByLevel(maxLevel + 1, 0);
	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
			continue;
		}

		long cellId = entry.first;
		int cellLevel = getCellLevel(cellId);

		++nDirtyAdjacenciesCells;
		++nDirtyAdjacenciesCellsByLevel[cellLevel];
	}

	// Group altered cells by their tree level
	std::vector<std::vector<long>> hierarchicalCellIds(maxLevel + 1);
	for (int level = 0; level <= maxLevel; ++level) {
		hierarchicalCellIds[level].reserve(nDirtyAdjacenciesCellsByLevel[level]);
	}

	for (const auto &entry : m_alteredCells) {
		AlterationFlags cellAlterationFlags = entry.second;
		if (!testAlterationFlags(cellAlterationFlags, FLAG_ADJACENCIES_DIRTY)) {
			continue;
		}

		long cellId = entry.first;
		int cellLevel = getCellLevel(cellId);
		hierarchicalCellIds[cellLevel].push_back(cellId);
	}

	// Update the adjacencies
	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	for (int level = 0; level <= maxLevel; ++level) {
		for (long cellId : hierarchicalCellIds[level]) {
			Cell &cell = m_cells[cellId];
			OctantInfo octantInfo = getCellOctant(cellId);
			Octant *octant = getOctantPointer(octantInfo);
			for (int face = 0; face < nCellFaces; ++face) {
				// Check if the face needs to be processed
				//
				// If the face has no adjacencies, we need to process it to
				// figure out if it is a border or it has some neighbours.
				//
				// If the face has only one adjacency and the neighbour's level
				// is smaller or equal than the one of the current cell, all
				// face adjacencies have already been found.
				//
				// If the face has only one adjacency and the neighbour's level
				// is greater than the one of the current cell, there are still
				// adjacencies to be found. The adjacency associated with the
				// face is an adjacency that was not created during this update
				// (cells are processed in level increasing order, only cells
				// with a level smaller or equal to the one of the current cell
				// may already have been processed during this update). Since
				// the cell is bigger that its neighbour it cannot have only
				// one adjacency, there are other adjacencies to be found in
				// the cell marked as dirty.
				//
				// If the face has multiple adjacencies, we need to figure out
				// if all the adjacencies have already been found. Current
				// adjacencies were not created during this update (since there
				// are multiple adjacencies, those adjacencies are cells smaller
				// than the current one, however since cells are processed in
				// level increasing order, only cells with a level smaller or
				// equal to the one of the current cell may already have been
				// processed during this update). To check if all adjacencies
				// have been found, we check if the area coverd by the current
				// adjacencies is equal to the face area.
				int nFaceAdjacencies = cell.getAdjacencyCount(face);
				if (nFaceAdjacencies > 0) {
					if (nFaceAdjacencies == 1) {
						long cellLevel = octant->getLevel();
						long neighId = cell.getAdjacency(face, 0);
						long neighLevel = getCellLevel(neighId);
						if (neighLevel <= cellLevel) {
							continue;
						}
					} else {
						uint64_t neighsArea = 0;
						const long *faceAdjacencies = cell.getAdjacencies(face);
						for (int k = 0; k < nFaceAdjacencies; ++k) {
							long neighId = faceAdjacencies[k];
							OctantInfo neighOctantInfo = getCellOctant(neighId);
							Octant *neighOctant = getOctantPointer(neighOctantInfo);
							neighsArea += neighOctant->getLogicalArea();
						}

						uint64_t faceArea = octant->getLogicalArea();
						if (faceArea == neighsArea) {
							continue;
						}
					}
				}

				// Find cell neighbours
				neighTreeIds.clear();
				neighGhostFlags.clear();
				m_tree->findNeighbours(octant, face, 1, neighTreeIds, neighGhostFlags);

				// Set the adjacencies
				//
				// Some of the neighbours may already be among the adjacencies
				// of the cell, this is not a problem because the function that
				// pushs the adjacency will insert only unique adjacencies.
				int nNeighs = neighTreeIds.size();
				for (int k = 0; k < nNeighs; ++k) {
					OctantInfo neighOctantInfo(neighTreeIds[k], !neighGhostFlags[k]);
					long neighId = getOctantId(neighOctantInfo);

					// Set cell data
					cell.pushAdjacency(face, neighId);

					// Set neighbour data
					int neighFace = oppositeFace[face];
					Cell &neigh = m_cells[neighId];
					neigh.pushAdjacency(neighFace, cellId);
				}
			}
		}
	}
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
bool VolOctree::_markCellForRefinement(long id)
{
	return setMarker(id, 1);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
bool VolOctree::_markCellForCoarsening(long id)
{
	return setMarker(id, -1);
}

/*!
	Resets the adaption marker of the specified cell.

	\param id the cell to be refined
	\result Returns true if the marker was properly reset, false otherwise.
*/
bool VolOctree::_resetCellAdaptionMarker(long id)
{
	return setMarker(id, 0);
}

/*!
	Returns the adaption marker of the specified cell.

	\param id is the id of the cell
	\return The adaption marker of the cell.
*/
adaption::Marker VolOctree::_getCellAdaptionMarker(long id)
{
	OctantInfo octantInfo = getCellOctant(id);
	if (!octantInfo.internal) {
		return adaption::MARKER_UNDEFINED;
	}

	int8_t treeMarker = m_tree->getMarker(octantInfo.id);
	if (treeMarker > 0) {
		return adaption::MARKER_REFINE;
	} else if (treeMarker < 0) {
		return adaption::MARKER_COARSEN;
	} else {
		return adaption::MARKER_NONE;
	}
}

/*!
	Set the marker on a cell.

	\param id is the id of the cell
	\param value is the value of the marker
*/
bool VolOctree::setMarker(long id, int8_t value)
{
	OctantInfo octantInfo = getCellOctant(id);
	if (!octantInfo.internal) {
		return false;
	}

	m_tree->setMarker(octantInfo.id, value);

	return true;
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
	\result Returns true if the falg was properly set, false otherwise.
*/
bool VolOctree::_enableCellBalancing(long id, bool enabled)
{
	OctantInfo octantInfo = getCellOctant(id);
	if (!octantInfo.internal) {
		return false;
	}

	m_tree->setBalance(octantInfo.id, enabled);

	return true;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] point is the point to be checked
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolOctree::isPointInside(const std::array<double, 3> &point) const
{
	bool isGhost;

	return (m_tree->getPointOwner(point, isGhost) != nullptr);
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the idof the cell
	\param[in] point is the point to be checked
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolOctree::isPointInside(long id, const std::array<double, 3> &point) const
{
	const Cell &cell = m_cells[id];
	ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

    int lowerLeftVertex  = 0;
	int upperRightVertex = pow(2, getDimension()) - 1;

	const std::array<double, 3> &lowerLeft  = getVertexCoords(cellVertexIds[lowerLeftVertex]);
	const std::array<double, 3> &upperRight = getVertexCoords(cellVertexIds[upperRightVertex]);

	const double EPS = getTol();
    for (int d = 0; d < 3; ++d){
        if (point[d] < lowerLeft[d] - EPS || point[d] > upperRight[d] + EPS) {
            return false;
        }
    }

    return true;
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] point is the point to be checked
	\result Returns the id of the cell the contains the point. If the point
	is not inside the patch, the function returns the id of the null element.
*/
long VolOctree::locatePoint(const std::array<double, 3> &point) const
{
	bool isGhost;
	uint32_t treeId = m_tree->getPointOwnerIdx(point, isGhost);
	if (treeId == std::numeric_limits<uint32_t>::max()) {
		return Element::NULL_ID;
	}

	OctantInfo octantInfo(treeId, !isGhost);
	return getOctantId(octantInfo);
}

/*!
	Internal function to set the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void VolOctree::_setTol(double tolerance)
{
	m_tree->setTol(tolerance);

	VolumeKernel::_setTol(tolerance);
}

/*!
	Internal function to reset the tolerance for the geometrical checks.
*/
void VolOctree::_resetTol()
{
	m_tree->setTol();

	double tolerance = m_tree->getTol();
	VolumeKernel::_setTol(tolerance);
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int VolOctree::_getDumpVersion() const
{
	const int DUMP_VERSION = 4;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolOctree::_dump(std::ostream &stream) const
{
	// List all octants
	std::size_t nOctants       = m_tree->getNumOctants();
	std::size_t nGhostsOctants = m_tree->getNumGhosts();

	std::vector<OctantInfo> octantInfoList;
	octantInfoList.reserve(nOctants + nGhostsOctants);

	for (std::size_t n = 0; n < nOctants; ++n) {
		octantInfoList.emplace_back(n, true);
	}

	for (std::size_t n = 0; n < nGhostsOctants; ++n) {
		octantInfoList.emplace_back(n, false);
	}

	// Dump tree data
	m_tree->dump(stream);

	// Dump kernel of vertex's containers
	//
	// We want the kernel state to be kept after the restore.
	m_vertices.dumpKernel(stream);

	// Dump kernel of cell's containers
	//
	// We want the kernel state to be kept after the restore.
	m_cells.dumpKernel(stream);

	// Dump vertex information
	//
	// This indexes will be read when importing the cells to keep the
	// association between numeration of tree vertices and patch vertices.
	std::unordered_set<long> dumpedVertices;
	dumpedVertices.reserve(getVertexCount());
	for (const OctantInfo &octantInfo : octantInfoList) {
		const Cell &cell = m_cells.at(getOctantId(octantInfo));
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
		for (int k = 0; k < m_cellTypeInfo->nVertices; ++k) {
			long vertexId = cellVertexIds[k];
			if (dumpedVertices.count(vertexId) == 0) {
				utils::binary::write(stream, vertexId);
#if BITPIT_ENABLE_MPI==1
				utils::binary::write(stream, getVertexRank(vertexId));
#else
				int dummyRank = 0;
				utils::binary::write(stream, dummyRank);
#endif
				dumpedVertices.insert(vertexId);
			}
		}
	}

	// Dump cell information
	//
	// This indexes will be read when importing the cells to keep the
	// association between numeration of tree octants and patch cells.
	for (const OctantInfo &octantInfo : octantInfoList) {
		utils::binary::write(stream, getOctantId(octantInfo));
	}

	for (const Cell &cell : getCells()) {
		utils::binary::write(stream, cell.getPID());
	}

	// Dump interfaces
	dumpInterfaces(stream);
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void VolOctree::_restore(std::istream &stream)
{
	// Restore tree
	m_tree->restore(stream);

	// Restore kernel of vertex's containers
	m_vertices.restoreKernel(stream);

	// Restore kernel of cell's containers
	m_cells.restoreKernel(stream);

	// Activate expert mode
	setExpert(true);

	// Restore cells
	size_t nOctants       = m_tree->getNumOctants();
	size_t nGhostsOctants = m_tree->getNumGhosts();

	StitchInfo stitchInfo;

	std::vector<OctantInfo> octantInfoList;
	octantInfoList.reserve(nOctants + nGhostsOctants);
	for (std::size_t n = 0; n < nOctants; ++n) {
		octantInfoList.emplace_back(n, true);
	}
	for (std::size_t n = 0; n < nGhostsOctants; ++n) {
		octantInfoList.emplace_back(n, false);
	}

	importCells(octantInfoList, stitchInfo, &stream);

	// Restore interfaces
	restoreInterfaces(stream);

	// De-activate expert mode
	setExpert(false);

	//
	// Restore bounding box
	//
	// The bounding box is frozen, it is necessary to update it manually.
	setBoundingBox();
}

/*!
	Get the native index of a cell.

	\param id is the id of the cell
	\result The native index of a cell.
*/
long VolOctree::_getCellNativeIndex(long id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	return octantInfo.id;
}

/*!
	Gets the origin of the patch.

	The origin is the lower-left-back corner of the box that defines the patch
	domain.

	\return The origin of the patch.
*/
std::array<double, 3> VolOctree::getOrigin() const
{
	return m_tree->getOrigin();
}

/*!
	Sets the origin of the patch.

	The origin is the lower-left-back corner.

	\param origin is the new origin of the patch
*/
void VolOctree::setOrigin(const std::array<double, 3> &origin)
{
	std::array<double, 3> translation = origin - getOrigin();
	translate(translation);
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
 */
void VolOctree::translate(const std::array<double, 3> &translation)
{
	m_tree->setOrigin(m_tree->getOrigin() + translation);

	VolumeKernel::translate(translation);

	// The bounding box is frozen, it is not updated automatically
	setBoundingBox();
}

/*!
	Gets the reference length of the patch domain.

	\return The the reference length of the patch domain.
*/
double VolOctree::getLength() const
{
	return m_tree->getL();
}

/*!
	Sets the the reference length of the patch domain.

	\param length is the reference length of the patch domain
*/
void VolOctree::setLength(double length)
{
	// Set the length
	m_tree->setL(length);

	// Set the new bounding box
	setBoundingBox();

	// If needed, update the discretization
	if (m_vertices.size() > 0) {
		// Create tree connectivity
		m_tree->computeConnectivity();

		// Update vertex coordinates
		std::unordered_set<long> alreadyEvaluated;
		for (const Cell &cell : m_cells) {
			OctantInfo octantInfo = getCellOctant(cell.getId());
			Octant *octant = getOctantPointer(octantInfo);
			ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			for (int k = 0; k < m_cellTypeInfo->nVertices; ++k) {
				long vertexId = cellVertexIds[k];
				if (alreadyEvaluated.count(vertexId) == 0) {
					Vertex &vertex = m_vertices.at(vertexId);
					vertex.setCoords(m_tree->getNode(octant, k));
					alreadyEvaluated.insert(vertexId);
				}
			}
		}

		// Destroy tree connectivity
		m_tree->clearConnectivity();
	}
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor vector
	\param[in] center is the center of the scaling
 */
void VolOctree::scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center)
{
	bool uniformScaling = true;
	uniformScaling &= (std::abs(scaling[0] - scaling[1]) > 1e-14);
	uniformScaling &= (std::abs(scaling[0] - scaling[2]) > 1e-14);
	assert(uniformScaling);
	if (!uniformScaling) {
		log::cout() << "VolOctree patch only allows uniform scaling." << std::endl;
		return;
	}

	std::array<double, 3> origin = m_tree->getOrigin();
	for (int n = 0; n < 3; ++n) {
		origin[n] = center[n] + scaling[n] * (origin[n] - center[n]);
	}
	m_tree->setOrigin(origin);

	m_tree->setL(m_tree->getL() * scaling[0]);

	VolumeKernel::scale(scaling, center);

	// The bounding box is frozen, it is not updated automatically
	setBoundingBox();
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
void VolOctree::_findCellNeighs(long id, const std::vector<long> *blackList, std::vector<long> *neighs) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	m_tree->findAllCodimensionNeighbours(octant, neighTreeIds, neighGhostFlags);

	int nNeighs = neighTreeIds.size();
	for (int i = 0; i < nNeighs; ++i) {
		OctantInfo neighOctantInfo(neighTreeIds[i], !neighGhostFlags[i]);
		long neighId = getOctantId(neighOctantInfo);

		if (!blackList || utils::findInOrderedVector<long>(neighId, *blackList) == blackList->end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
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
void VolOctree::_findCellEdgeNeighs(long id, int edge, const std::vector<long> *blackList, std::vector<long> *neighs) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return;
	}

	// Get octant info
	const OctantInfo octantInfo = getCellOctant(id);

	// Get edge neighbours
	int codimension = getDimension() - 1;
	findOctantCodimensionNeighs(octantInfo, edge, codimension, blackList, neighs);

	// Add face neighbours
	//
	// Get all face neighbours and select the ones that contains the edge
	// for which the neighbours are requested. To correctly consider these
	// neighbours, the following logic can be used:
	//   - if a face/edge neighbour has the same level or a lower level than
	//     the current cell, then it certainly is also a vertex neighbour;
	//   - if a face/edge neighbour has a higher level than the current cell,
	//     it is necessary to check if the neighbour actually contains the
	//     edge.
	//
	std::vector<long> faceNeighs;
	const Octant *octant = getOctantPointer(octantInfo);
	int octantLevel = octant->getLevel();
	for (int face : m_octantLocalFacesOnEdge[edge]) {
		faceNeighs.clear();
		_findCellFaceNeighs(id, face, blackList, &faceNeighs);
		for (long neighId : faceNeighs) {
			const OctantInfo neighOctantInfo = getCellOctant(neighId);
			const Octant *neighOctant = getOctantPointer(neighOctantInfo);
			int neighOctantLevel = neighOctant->getLevel();
			if (neighOctantLevel <= octantLevel) {
				utils::addToOrderedVector<long>(neighId, *neighs);
			} else if (m_tree->isEdgeOnOctant(octant, edge, neighOctant)) {
				utils::addToOrderedVector<long>(neighId, *neighs);
			}
		}
	}
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

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
void VolOctree::_findCellVertexNeighs(long id, int vertex, const std::vector<long> *blackList, std::vector<long> *neighs) const
{
	// Get octant info
	const OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	// Get vertex neighbours
	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	m_tree->findAllNodeNeighbours(octant, vertex, neighTreeIds, neighGhostFlags);

	int nNeighs = neighTreeIds.size();
	for (int i = 0; i < nNeighs; ++i) {
		OctantInfo neighOctantInfo(neighTreeIds[i], !neighGhostFlags[i]);
		long neighId = getOctantId(neighOctantInfo);

		if (!blackList || utils::findInOrderedVector<long>(neighId, *blackList) == blackList->end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
}

/*!
	Finds the neighbours for the given co-dimension of the specified octant.

	Only the neighbours for the specified co-dimension are found, neighbours
	of higher co-dimensions are not inserted in the returned list.

	\param octantInfo the data of the octant
	\param codimension is the co-dimension
	\param index is the local index of the entity (vertex, edge or face)
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a pointer to a unique list of ordered cell ids
	or a null pointer if no cells should be excluded from the search
	\param[in,out] neighs is the vector were the neighbours will be stored.
	The vector is not cleared before adding the neighbours, it is extended
	by appending all the neighbours found by this function
*/
void VolOctree::findOctantCodimensionNeighs(const OctantInfo &octantInfo, int index, int codimension,
                                            const std::vector<long> *blackList, std::vector<long> *neighs) const
{
	int dimension = getDimension();
	if (codimension > dimension || codimension <= 0) {
		return;
	}

	const Octant *octant = getOctantPointer(octantInfo);

	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	m_tree->findNeighbours(octant, index, codimension, neighTreeIds, neighGhostFlags);

	int nNeighs = neighTreeIds.size();
	for (int i = 0; i < nNeighs; ++i) {
		OctantInfo neighOctantInfo(neighTreeIds[i], !neighGhostFlags[i]);
		long neighId = getOctantId(neighOctantInfo);

		if (!blackList || utils::findInOrderedVector<long>(neighId, *blackList) == blackList->end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
}

/*!
	Finds the face of the supposed neighbour that adjoins the target face.

	\param cellId is the id of the cell
	\param cellFace is the target face of the cell
	\param neighId is the id of a supposed neighbour of the cell
	\result The face of the neighbour which adjoins the target face. If the
	two cells are not neighbours, a negative number is returned.
 */
int VolOctree::findAdjoinNeighFace(const Cell &cell, int cellFace, const Cell &neigh) const
{
	long cellId = cell.getId();

	// Get the neighbour face
	int neighFace = m_tree->getOppface()[cellFace];

	// Check if the two cells are neighbours
	int nNeighFaceAdjacencies = neigh.getAdjacencyCount(neighFace);
	const long *neighFaceAdjacencies = neigh.getAdjacencies(neighFace);
	for (int k = 0; k < nNeighFaceAdjacencies; ++k) {
		if (neighFaceAdjacencies[k] == cellId) {
			return neighFace;
		}
	}

	return -1;
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
bool VolOctree::isSameFace(const Cell &cell_A, int face_A, const Cell &cell_B, int face_B) const
{
	// Check if the face matches
	if (m_tree->getOppface()[face_A] != face_B) {
		return false;
	}

	// Check if the two cells are neighbours
	long cellId_A = cell_A.getId();

	int nNeighFaceAdjacencies_B = cell_B.getAdjacencyCount(face_B);
	const long *neighFaceAdjacencies_B = cell_B.getAdjacencies(face_B);
	for (int k = 0; k < nNeighFaceAdjacencies_B; ++k) {
		if (neighFaceAdjacencies_B[k] == cellId_A) {
			return true;
		}
	}

	return false;
}

}
