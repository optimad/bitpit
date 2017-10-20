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

#include <cassert>
#include <cmath>

#include "bitpit_IO.hpp"

#include "voloctree.hpp"

namespace bitpit {

/*!
	\class VolOctree
	\ingroup volumepatches

	\brief The VolOctree defines a Octree patch.

	VolOctree defines a Octree patch.
*/

/*!
	Creates an uninitialized patch.
*/
VolOctree::VolOctree()
	: VolumeKernel(false)
{
	// Create the tree
#if BITPIT_ENABLE_MPI==1
	m_tree = std::unique_ptr<PabloUniform>(new PabloUniform(PabloUniform::DEFAULT_LOG_FILE, MPI_COMM_NULL));
#else
	m_tree = std::unique_ptr<PabloUniform>(new PabloUniform(PabloUniform::DEFAULT_LOG_FILE));
#endif

	// Initialize
	initialize();

	// Reset
	//
	// The function that resets the patch is virtual, but since is called
	// from the constructor of the patch kernel only the base function is
	// called.
	__reset(false);
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
*/
VolOctree::VolOctree(const int &dimension,
				     std::array<double, 3> origin, double length, double dh )
	: VolOctree(PatchManager::AUTOMATIC_ID, dimension, origin, length, dh)
{
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed cell size of the initial refinement
*/
VolOctree::VolOctree(const int &id, const int &dimension,
				 std::array<double, 3> origin, double length, double dh )
	: VolumeKernel(id, dimension, false)
{
	// Create the tree
#if BITPIT_ENABLE_MPI==1
	m_tree = std::unique_ptr<PabloUniform>(
	    new PabloUniform(origin[0], origin[1], origin[2], length, dimension,
	                     PabloUniform::DEFAULT_LOG_FILE, MPI_COMM_NULL));
#else
	m_tree = std::unique_ptr<PabloUniform>(
	    new PabloUniform(origin[0], origin[1], origin[2], length, dimension,
	                     PabloUniform::DEFAULT_LOG_FILE));
#endif

	// Initialize
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

	// Inizializzazione dell'octree
	double initial_level = ceil(log2(std::max(1., length / dh)));

	m_tree->setMarker((uint32_t) 0, initial_level);
}

/*!
	Creates a new patch.

	\param tree is the tree that will be used
	\param adopter is a pointer to the tree adopter
*/
VolOctree::VolOctree(std::unique_ptr<PabloUniform> &&tree, std::unique_ptr<PabloUniform> *adopter)
	: VolOctree(PatchManager::AUTOMATIC_ID, std::move(tree), adopter)
{
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param tree is the tree that will be used
	\param adopter is a pointer to the tree adopter
*/
VolOctree::VolOctree(const int &id, std::unique_ptr<PabloUniform> &&tree, std::unique_ptr<PabloUniform> *adopter)
	: VolumeKernel(id, tree->getDim(), false)
{
	// Associate the tree
	assert(tree);
	m_tree.swap(tree);

	// Initialize
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

#if BITPIT_ENABLE_MPI==1
	// Set the communicator
	PatchKernel::setCommunicator(m_tree->getComm());

	// Set the partitioned flag
	setPartitioned(!m_tree->getSerial());
#endif

	// Sync the patch with the tree
	sync(true, true, false);

	// Set the bounding
	setBoundingBox();

	// Set the aopter
	setTreeAdopter(adopter);
}

/*!
	Creates a new patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
VolOctree::VolOctree(std::istream &stream)
	: VolumeKernel(false)
{
	initialize();
	restore(stream);
}

/*!
	Copy constructor.

	\param stream is the stream to read from
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
	Initialize the data structures of the patch.
*/
void VolOctree::initialize()
{
	log::cout() << ">> Initializing Octree mesh" << std::endl;

	// Reset the cell and interface type info
	m_cellTypeInfo      = nullptr;
	m_interfaceTypeInfo = nullptr;

	// Reset the tree entruster
	m_treeAdopter = nullptr;

	// This patch need to be spawn
	setSpawnStatus(SPAWN_NEEDED);

	// This patch supports adaption
	setAdaptionStatus(ADAPTION_CLEAN);

#if BITPIT_ENABLE_MPI==1
	// This patch supports partitioning
	setPartitioningStatus(PARTITIONING_CLEAN);
#endif

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

	// The tree is only evaluating the bounding box of the internal octants,
	// we need to consider also ghosts cells.
	for (auto ghostItr = ghostBegin(); ghostItr != ghostEnd(); ++ghostItr) {
		ConstProxyVector<long> ghostVertexIds = ghostItr->getVertexIds();
		int nGhostVertices = ghostVertexIds.size();
		for (int i = 0; i < nGhostVertices; ++i) {
			const std::array<double, 3> coords = m_vertices[ghostVertexIds[i]].getCoords();
			for (int d = 0; d < 3; ++d) {
				minPoint[d] = std::min(coords[d], minPoint[d]);
				maxPoint[d] = std::max(coords[d], maxPoint[d]);
			}
		}
	}

	// Set the bounding box
	setBoundingBox(minPoint, maxPoint);
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolOctree::evalCellVolume(const long &id) const
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
std::array<double, 3> VolOctree::evalCellCentroid(const long &id) const
{
	OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);

	return m_tree->getCenter(octant);
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolOctree::evalCellSize(const long &id) const
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
double VolOctree::evalInterfaceArea(const long &id) const
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
std::array<double, 3> VolOctree::evalInterfaceNormal(const long &id) const
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
VolOctree::OctantInfo VolOctree::getCellOctant(const long &id) const
{
	OctantInfo octantInfo;
	octantInfo.internal = m_cells[id].isInterior();
	if (octantInfo.internal) {
		octantInfo.id = m_cellToOctant.at(id);
	} else {
		octantInfo.id = m_cellToGhost.at(id);
	}

	return octantInfo;
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
int VolOctree::getCellLevel(const long &id) const
{
	OctantInfo octantInfo = getCellOctant(id);

	const Octant* octant;
	if (octantInfo.internal) {
		octant = m_tree->getOctant(octantInfo.id);
	} else {
		octant = m_tree->getGhostOctant(octantInfo.id);
	}
	return m_tree->getLevel(octant);
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
	bool emtpyPatch = (getCellCount() == 0);
	ParaTree::Operation lastTreeOperation = m_tree->getLastOperation();
	if (lastTreeOperation == ParaTree::OP_INIT && emtpyPatch) {
		m_tree->adapt();
		updateInfo = sync(true, true, trackSpawn);
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

		// Track internal octants that will be coarsend/refined
		long nOctants = m_tree->getNumOctants();


		uint32_t treeId = 0;
		while (treeId < (uint32_t) nOctants) {
			int8_t marker = m_tree->getPreMarker(treeId);
			if (marker == 0) {
				treeId++;
				continue;
			}

			int nUpdatedOctants;
			adaption::Type adaptionType;
			if (marker > 0) {
				nUpdatedOctants = 1;
				adaptionType    = adaption::TYPE_REFINEMENT;
			} else {
				nUpdatedOctants = pow(2, getDimension());
				adaptionType    = adaption::TYPE_COARSENING;
			}

			std::size_t adaptionInfoId = adaptionData.create(adaptionType, adaption::ENTITY_CELL, currentRank);
			adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
			adaptionInfo.previous.reserve(nUpdatedOctants);
			for (int k = 0; k < nUpdatedOctants; ++k) {
				OctantInfo octantInfo(treeId, true);

				adaptionInfo.previous.emplace_back();
				long &cellId = adaptionInfo.previous.back();
				cellId = getOctantId(octantInfo);
				treeId++;
			}
		}

		// Ghost cells will be removed
		std::size_t adaptionInfoId = adaptionData.create(adaption::TYPE_DELETION, adaption::ENTITY_CELL, currentRank);
		adaption::Info &adaptionInfo = adaptionData[adaptionInfoId];
		adaptionInfo.previous.reserve(getGhostCount());
		for (auto itr = ghostBegin(); itr != ghostEnd(); ++itr) {
			adaptionInfo.previous.emplace_back();
			long &deletedGhostId = adaptionInfo.previous.back();
			deletedGhostId = itr.getId();
		}
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

	bool emtpyPatch   = (getCellCount() == 0);
	bool buildMapping = !emtpyPatch;
	bool updated = m_tree->adapt(buildMapping);

	if (!updated && !emtpyPatch) {
		log::cout() << " Already updated" << std::endl;
		return std::vector<adaption::Info>();
	}
	log::cout() << " Done" << std::endl;

	// Sync the patch
	return sync(true, true, trackAdaption);
}

/*!
	Cleanup patch data structured after the adaption.
*/
void VolOctree::_adaptionCleanup()
{
	// Nothing to do
}

/*!
	Syncronizes the patch with the underlying octree.

	\param updateOctantMaps if set to true the cell-to-octant maps will
	be updated, otherwise the function assumes that someone has already
	updated those maps
	\param generateInterfaces check if the interfaces of the imported cells
	will be generated
	\param trackChanges if set to true the changes to the patch will be
	tracked
	\result Returns all the changes applied to the patch.
*/
std::vector<adaption::Info> VolOctree::sync(bool updateOctantMaps, bool generateInterfaces, bool trackChanges)
{
	log::cout() << ">> Syncing patch..." << std::endl;

	// Detect if we are in import-from-scratch mode
	//
	// In import-from-scratch mode we start form an empty patch and we need
	// to import all the octants of the tree. If we are importing the tree
	// from scratch there are no cells to delete/renumber.
	bool importFromScratch = (getCellCount() == 0);

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
	while (treeId < (uint32_t) nOctants) {
		// Octant mapping
		std::vector<uint32_t> mapper_octantMap;
		std::vector<bool> mapper_ghostFlag;
		std::vector<int> mapper_octantRank;
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
		int nPreviousTreeIds;
		if (!importFromScratch) {
			nPreviousTreeIds = mapper_octantMap.size();
			for (int k = 0; k < nPreviousTreeIds; ++k) {
#if BITPIT_ENABLE_MPI==1
				// Only local cells can be deleted
				if (mapper_octantRank[k] != currentRank) {
					continue;
				}
#endif

				// Mark previous octant for deletion
				uint32_t previousTreeId = mapper_octantMap[k];
				OctantInfo previousOctantInfo(previousTreeId, !mapper_ghostFlag[k]);
				long cellId = getOctantId(previousOctantInfo);
				deletedOctants.emplace_back(cellId, adaptionType);

				unmappedOctants[previousTreeId] = false;
			}
		} else {
			nPreviousTreeIds = 0;
		}

		// Adaption tracking
		//
		// The adaption info associated to the octants that has been received
		// from external partitions will contain the current octants sorted by
		// their tree id (we are looping over the octants in that order), this
		// is the same order that will be used on the processor that has sent
		// the octants. Since the order is the same, the two processors are able
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
			// A coarsening can merge togheter cells of different processors.
			// However, since the coarsening is limited to one level, the
			// previous cells will always be internal or among the ghost of
			// the current processor.
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
		// Cells that have been send to other processors need to be removed
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
			clearGhostOwners(true);

			for (uint32_t ghostTreeId = 0; ghostTreeId < nPreviousGhosts; ++ghostTreeId) {
				OctantInfo ghostOctantInfo(ghostTreeId, false);
				long ghostId = getOctantId(ghostOctantInfo);
				deletedOctants.emplace_back(ghostId, adaption::TYPE_DELETION);
			}
		}
#endif

		// Remove unmapped octants
		//
		// A coarsening that merges cells from different processors, can leave, on
		// the processors which own the ghost octants involved in the coarsening,
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

	// Reset cell-to-octant and octant-to-cell map
	if (updateOctantMaps) {
		log::cout() << ">> Resetting cell-to-octant and octant-to-cell maps...";

		updateCellOctantMaps(deletedOctants, renumberedOctants, addedOctants);
		std::vector<RenumberInfo>().swap(renumberedOctants);

		log::cout() << " Done" << std::endl;
	}

	// Remove deleted octants
	StitchInfo stitchInfo;
	if (deletedOctants.size() > 0) {
		log::cout() << ">> Removing non-existing cells...";

		// Track changes
		//
		// The adaption info associated to the octants that has been sent
		// to external partitions will contain the current octants sorted by
		// their tree id (they were added to the deleted octants list in that
		// order), this is the same order that will be used on the processor
		// that has received the octants. Since the order is the same, the two
		// processors are able to exchange cell data without any additional
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
			for (const long &interfaceId : removedInterfaces) {
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

	// Import added octants
	std::vector<long> createdCells;
	if (addedOctants.size() > 0) {
		log::cout() << ">> Importing new octants...";

		createdCells = importCells(addedOctants, stitchInfo, generateInterfaces);

		log::cout() << " Done" << std::endl;
		log::cout() << ">> Octants imported: " <<  addedOctants.size() << std::endl;
	}

	StitchInfo().swap(stitchInfo);

	// Rebuild the ghost information
#if BITPIT_ENABLE_MPI==1
	buildGhostExchangeData();
#endif

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

		// Track created ghosts cells
#if BITPIT_ENABLE_MPI==1
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
			for (const long &interfaceId : createdInterfaces) {
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
	Update the cell-to-octants maps

	\param deletedOctants contains the information about the deleted octants
	\param renumberedOctants contains the information about the renumbered
	octants
	\param addedOctants contains the information about the added octants
*/
void VolOctree::updateCellOctantMaps(std::vector<DeleteInfo> &deletedOctants,
                                     std::vector<RenumberInfo> &renumberedOctants,
                                     std::vector<OctantInfo> &addedOctants)
{
	// Reset the ghost maps
	m_cellToGhost.clear();
	m_ghostToCell.clear();

	// Reserve space for the maps
	long nOctants = m_tree->getNumOctants();
	m_cellToOctant.reserve(nOctants);
	m_octantToCell.reserve(nOctants);

	long nGhostsOctants = m_tree->getNumGhosts();
	m_cellToGhost.reserve(nGhostsOctants);
	m_ghostToCell.reserve(nGhostsOctants);

	// Remove cell-to-tree associations for cells that will be deleted
	for (const DeleteInfo &deleteInfo : deletedOctants) {
		long cellId = deleteInfo.cellId;
		if (!m_cells[cellId].isInterior()) {
			continue;
		}

		const OctantInfo &octantInfo = getCellOctant(cellId);
		uint32_t treeId = octantInfo.id;

		m_cellToOctant.erase(cellId);
		m_octantToCell.erase(treeId);
	}

	// Remove old cell-to-tree associations for renumbered cells
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		long cellId = renumberInfo.cellId;
		if (!m_cells[cellId].isInterior()) {
			continue;
		}

		const OctantInfo &previousOctantInfo = getCellOctant(cellId);
		uint32_t previousTreeId = previousOctantInfo.id;

		m_octantToCell.erase(previousTreeId);
	}

	// Create cell-to-tree associations for renumbered cells
	for (const RenumberInfo &renumberInfo : renumberedOctants) {
		long cellId = renumberInfo.cellId;
		if (!m_cells[cellId].isInterior()) {
			continue;
		}

		uint32_t treeId = renumberInfo.newTreeId;

		m_cellToOctant[cellId] = treeId;
		m_octantToCell[treeId] = cellId;
	}

	// Create cell-to-tree associations for cells that wil be added
	for (OctantInfo &octantInfo : addedOctants) {
		long cellId = generateCellId();

		if (octantInfo.internal) {
			m_cellToOctant.insert({{cellId, octantInfo.id}});
			m_octantToCell.insert({{octantInfo.id, cellId}});
		} else {
			m_cellToGhost.insert({{cellId, octantInfo.id}});
			m_ghostToCell.insert({{octantInfo.id, cellId}});
		}
	}
}

/*!
	Imports a list of octants into the patch.

	\param octantInfoList is the list of octant to import
	\param stitchInfo is the list of vertices that will be used to stitch the
	the new octants
	\param generateInterfaces check if the interfaces of the imprted cells
	will be generated
*/
std::vector<long> VolOctree::importCells(std::vector<OctantInfo> &octantInfoList,
										 StitchInfo &stitchInfo, bool generateInterfaces)
{
	// Create the new vertices
	int nCellVertices = m_cellTypeInfo->nVertices;
	for (OctantInfo &octantInfo : octantInfoList) {
		Octant *octant = getOctantPointer(octantInfo);
		for (int k = 0; k < nCellVertices; ++k) {
			uint64_t vertexTreeMorton = m_tree->getNodeMorton(octant, k);
			if (stitchInfo.count(vertexTreeMorton) == 0) {
				// Vertex coordinates
				std::array<double, 3> nodeCoords = m_tree->getNode(octant, k);

				// Create the vertex
				VertexIterator vertexIterator = addVertex(std::move(nodeCoords));

				// Add the vertex to the stitching info
				stitchInfo[vertexTreeMorton] = vertexIterator->getId();
			}
		}
	}

	// Add the cells
	size_t octantInfoListSize = octantInfoList.size();
	m_cells.reserve(m_cells.size() + octantInfoListSize);

	std::vector<long> createdCells(octantInfoListSize);
	for (size_t k = 0; k < octantInfoListSize; ++k) {
		OctantInfo &octantInfo = octantInfoList[k];

		// Id that will be assigned to the cell
		long cellId = getOctantId(octantInfo);

		// Octant connectivity
		Octant *octant = getOctantPointer(octantInfo);

		// Cell connectivity
		std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
		for (int k = 0; k < nCellVertices; ++k) {
			uint64_t vertexTreeMorton = m_tree->getNodeMorton(octant, k);
			cellConnect[k] = stitchInfo.at(vertexTreeMorton);
		}

		// Add cell
		addCell(m_cellTypeInfo->type, octantInfo.internal, std::move(cellConnect), cellId);

		// If the cell is a ghost set its owner
#if BITPIT_ENABLE_MPI==1
		if (!octantInfo.internal) {
			uint64_t globalTreeId = m_tree->getGhostGlobalIdx(octantInfo.id);
			int rank = m_tree->getOwnerRank(globalTreeId);

			setGhostOwner(cellId, rank, false);
		}
#endif

		// Add the cell to the list of created cells
		createdCells[k] = cellId;
	}

	// Build adjacencies
	updateAdjacencies(createdCells, false);
	if (generateInterfaces) {
		updateInterfaces(createdCells, false);
	}

	// Done
	return createdCells;
}

/*!
	Remove a list of octants from the patch.

	\param deletedOctants contains a list with the information of the deleted
	octants
	\result Returns the vertices that will be used to stitch  faces created after deleting the octants.
*/
VolOctree::StitchInfo VolOctree::deleteCells(std::vector<DeleteInfo> &deletedOctants)
{
	// Info of the cells
	int nCellVertices = m_cellTypeInfo->nVertices;
	const std::vector<std::vector<int>> &cellLocalFaceConnect = m_cellTypeInfo->faceConnect;

	// Info on the faces
	const int &nInterfaceVertices = m_interfaceTypeInfo->nVertices;

	// List of cells ot delete
	std::unordered_set<long> deadCells;
	deadCells.reserve(deletedOctants.size());
	for (const DeleteInfo &deleteInfo : deletedOctants) {
		deadCells.insert(deleteInfo.cellId);
	}

	// Delete the cells
	std::unordered_set<long> deadVertices;
	std::unordered_set<long> deadInterfaces;
	std::unordered_set<long> danglingCells;
	for (long cellId : deadCells) {
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

		// List of interfaces to delete
		//
		// All the interfaces of the cell will be deleted, this means that the
		// neighbours that are not deleted will have a face not connected to
		// anything. Those faces are called dangling faces and a cell with
		// dangling faces is called dangling cell.
		int nCellInterfaces = cell.getInterfaceCount();
		const long *interfaces = cell.getInterfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];
			if (interfaceId < 0) {
				continue;
			}

			// Interfaces has to be considered just once
			if (deadInterfaces.count(interfaceId) > 0) {
				continue;
			}

			// Find if the face associated to the interface will be dangling
			Interface &interface = m_interfaces[interfaceId];

			int danglingSide = -1;
			if (!interface.isBorder()) {
				if (deadCells.count(interface.getOwner()) == 0) {
					danglingSide = 0;
				} else if (deadCells.count(interface.getNeigh()) == 0) {
					danglingSide = 1;
				}
			}

			// Handle dangling faces
			if (danglingSide >= 0) {
				// Info on the dangling face
				long danglingCellId;
				long danglingNeighId;
				int danglingCellFace;
				if (danglingSide == 0) {
					danglingCellId   = interface.getOwner();
					danglingNeighId  = interface.getNeigh();
					danglingCellFace = interface.getOwnerFace();
				} else {
					danglingCellId   = interface.getNeigh();
					danglingNeighId  = interface.getOwner();
					danglingCellFace = interface.getNeighFace();
				}

				Cell &danglingCell = m_cells[danglingCellId];
				danglingCells.insert(danglingCellId);

				// Since the dangling cell will not be deleted, we have to
				// updated its interface and adjacency data structures.
				int cellInterfaceIndex = danglingCell.findInterface(danglingCellFace, interfaceId);
				danglingCell.deleteInterface(danglingCellFace, cellInterfaceIndex);

				int cellAdjacencyIndex = danglingCell.findAdjacency(danglingCellFace, danglingNeighId);
				danglingCell.deleteAdjacency(danglingCellFace, cellAdjacencyIndex);
			}

			// Add the interface to the list of interfaces to delete
			deadInterfaces.insert(interfaceId);
		}

		// Delete cell
		deleteCell(cellId, false, true);
	}

	m_cells.flush();

	// Delete the interfaces
	for (long interfaceId : deadInterfaces) {
		deleteInterface(interfaceId, false, true);
	}

	m_interfaces.flush();

	// All the vertices belonging to the dangling cells has to be kept
	//
	// The vertices of the dangling faces need to be kept because there are
	// still cells using them. However it's not enough to consider only the
	// vertices on the dangling faces, we have to consider the vertices of
	// the whole cell. That's because we may need to keep vertices on the
	// edges of the cell, and those vertices may not be on interfaces of the
	// dangling faces.
	//
	// Morover we need to build a map between the patch numbering and the
	// octree numbering of the vertices of the dangling cells. This map will
	// be used when imprting the octants to stitch the imported octants to
	// the existing cells.
	StitchInfo stitchVertices;
	for (const long cellId : danglingCells) {
		// Vertices of the cell
		const Cell &cell = m_cells[cellId];
		ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

		OctantInfo octantInfo = getCellOctant(cellId);
		Octant *octant = getOctantPointer(octantInfo);

		for (int k = 0; k < nCellVertices; ++k) {
			long vertexId = cellVertexIds[k];
			uint64_t vertexTreeMorton = m_tree->getNodeMorton(octant, k);
			stitchVertices.insert({vertexTreeMorton, vertexId});
			deadVertices.erase(vertexId);
		}

		// Vertices of all other interfaces left of the cell
		int nCellInterfaces = cell.getInterfaceCount();
		const long *interfaces = cell.getInterfaces();
		for (int k = 0; k < nCellInterfaces; ++k) {
			long interfaceId = interfaces[k];
			if (interfaceId < 0) {
				continue;
			}

			const Interface &interface = m_interfaces[interfaceId];
			if (interface.isBorder()) {
				continue;
			}

			long ownerId  = interface.getOwner();
			int ownerFace = interface.getOwnerFace();

			const Cell &ownerCell = m_cells[ownerId];
			ConstProxyVector<long> ownerCellVertexIds = ownerCell.getVertexIds();

			OctantInfo ownerOctantInfo = getCellOctant(ownerId);
			Octant *ownerOctant = getOctantPointer(ownerOctantInfo);

			const std::vector<int> &localFaceConnect = cellLocalFaceConnect[ownerFace];
			for (int k = 0; k < nInterfaceVertices; ++k) {
				long vertexId = ownerCellVertexIds[localFaceConnect[k]];
				uint64_t vertexTreeMorton = m_tree->getNodeMorton(ownerOctant, localFaceConnect[k]);
				stitchVertices.insert({vertexTreeMorton, vertexId});
				deadVertices.erase(vertexId);
			}
		}
	}

	// Delete the vertices
	for (long vertexId : deadVertices) {
		deleteVertex(vertexId, true);
	}

	m_vertices.flush();

	// Done
	return stitchVertices;
}

/*!
	Build the adjacencies the cells.
*/
void VolOctree::updateAdjacencies(const std::vector<long> &cellIds, bool resetAdjacencies)
{
	// Face information
	int nCellFaces = 2 * getDimension();
	uint8_t oppositeFace[nCellFaces];
	m_tree->getOppface(oppositeFace);

	// Reset the adjacencies
	if (resetAdjacencies) {
		for (long cellId : cellIds) {
			m_cells[cellId].resetAdjacencies();
		}
	}

	// Sort the cells beased on their tree level
	int maxLevel = m_tree->getMaxDepth();
	size_t averageSize = cellIds.size() / (maxLevel + 1);
	std::vector<std::vector<long>> hierarchicalCellIds(maxLevel + 1);
	for (int level = 0; level <= maxLevel; ++level) {
		hierarchicalCellIds[level].reserve(averageSize);
	}

	for (long cellId : cellIds) {
		int cellLevel = getCellLevel(cellId);
		hierarchicalCellIds[cellLevel].push_back(cellId);
	}

	// Update the adjacencies
	FaceInfoSet processedFaces;
	processedFaces.reserve(cellIds.size() * getDimension());

	for (int level = 0; level <= maxLevel; ++level) {
		for (long cellId : hierarchicalCellIds[level]) {
			Cell &cell = m_cells[cellId];
			OctantInfo octantInfo = getCellOctant(cellId);
			for (int face = 0; face < nCellFaces; ++face) {
				FaceInfo currentFaceInfo(cellId, face);
				if (processedFaces.count(currentFaceInfo) > 0) {
					continue;
				}

				// Find cell neighbours
				std::vector<uint32_t> neighTreeIds;
				std::vector<bool> neighGhostFlags;
				if (octantInfo.internal) {
					m_tree->findNeighbours(octantInfo.id, face, 1, neighTreeIds, neighGhostFlags);
				} else {
					m_tree->findGhostNeighbours(octantInfo.id, face, 1, neighTreeIds, neighGhostFlags);
				}

				// Set the adjacencies
				//
				// Adjacencies will processed twice, once while processing the
				// current cell, and once while processing the neighbour cell.
				// However they will be set only once, because the function that
				// insert the adjacency in the cell will insert only unique
				// adjacencies.
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

					FaceInfo neighFaceInfo(neighId, neighFace);
					processedFaces.insert(neighFaceInfo);
				}
			}
		}
	}
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
bool VolOctree::_markCellForRefinement(const long &id)
{
	return set_marker(id, 1);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
bool VolOctree::_markCellForCoarsening(const long &id)
{
	return set_marker(id, -1);
}

/*!
	Set the marker on a cell.

	\param id is the id of the cell
	\param value is the value of the marker
*/
bool VolOctree::set_marker(const long &id, const int8_t &value)
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
*/
bool VolOctree::_enableCellBalancing(const long &id, bool enabled)
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
bool VolOctree::isPointInside(const std::array<double, 3> &point)
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
bool VolOctree::isPointInside(const long &id, const std::array<double, 3> &point)
{
	const Cell &cell = m_cells[id];
	ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

    int lowerLeftVertex  = 0;
	int upperRightVertex = pow(2, getDimension()) - 1;

	std::array<double, 3> lowerLeft  = getVertexCoords(cellVertexIds[lowerLeftVertex]);
	std::array<double, 3> upperRight = getVertexCoords(cellVertexIds[upperRightVertex]);

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
long VolOctree::locatePoint(const std::array<double, 3> &point)
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
	const int DUMP_VERSION = 2;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolOctree::_dump(std::ostream &stream) const
{
	// Dump tree data
	m_tree->dump(stream);

	size_t nOctants = m_octantToCell.size();
	for (size_t n = 0; n < nOctants; ++n) {
		utils::binary::write(stream, m_octantToCell.at(n));
	}

	size_t nGhosts = m_ghostToCell.size();
	for (size_t n = 0; n < nGhosts; ++n) {
		utils::binary::write(stream, m_ghostToCell.at(n));
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
	//
	// Restore tree data
	//
#if ENABLE_MPI==1
	m_tree->setComm(getCommunicator());
#endif
	m_tree->restore(stream);

	// Restore octant to cell map
	size_t nOctants = m_tree->getNumOctants();
	m_cellToOctant.reserve(nOctants);
	m_octantToCell.reserve(nOctants);
	for (size_t n = 0; n < nOctants; ++n) {
		long cellId;
		utils::binary::read(stream, cellId);

		m_cellToOctant.insert({cellId, n});
		m_octantToCell.insert({n, cellId});
	}

	// Restore ghost to cell map
	size_t nGhosts = m_tree->getNumGhosts();

	m_cellToGhost.reserve(nGhosts);
	m_ghostToCell.reserve(nGhosts);
	for (size_t n = 0; n < nGhosts; ++n) {
		long cellId;
		utils::binary::read(stream, cellId);

		m_cellToGhost.insert({cellId, n});
		m_ghostToCell.insert({n, cellId});
	}

	//
	// Sync the patch
	//
	sync(false, false, false);

	//
	// Restore the interfaces
	//
	setExpert(true);
	restoreInterfaces(stream);
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
void VolOctree::translate(std::array<double, 3> translation)
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
 */
void VolOctree::scale(std::array<double, 3> scaling)
{
	bool uniformScaling = true;
	uniformScaling &= (std::abs(scaling[0] - scaling[1]) > 1e-14);
	uniformScaling &= (std::abs(scaling[0] - scaling[2]) > 1e-14);
	assert(uniformScaling);
	if (!uniformScaling) {
		log::cout() << "octree patch only allows uniform scaling)" << std::endl;
		return;
	}

	m_tree->setL(m_tree->getL() * scaling[0]);

	VolumeKernel::scale(scaling);

	// The bounding box is frozen, it is not updated automatically
	setBoundingBox();
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
void VolOctree::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return;
	}

	// Get edge neighbours
	int codimension = getDimension() - 1;
	findCellCodimensionNeighs(id, edge, codimension, blackList, neighs);

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
	const OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);
	int octantLevel = m_tree->getLevel(octant);
	for (int face : m_octantLocalFacesOnEdge[edge]) {
		faceNeighs.clear();
		_findCellFaceNeighs(id, face, blackList, &faceNeighs);
		for (long neighId : faceNeighs) {
			const OctantInfo neighOctantInfo = getCellOctant(neighId);
			const Octant *neighOctant = getOctantPointer(neighOctantInfo);
			int neighOctantLevel = m_tree->getLevel(neighOctant);
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
	\param blackList is a list of cells that are excluded from the search
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void VolOctree::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	// Get vertex neighbours
	int codimension = getDimension();
	findCellCodimensionNeighs(id, vertex, codimension, blackList, neighs);

	// Add edge and face neighbours
	//
	// Get all face and edge neighbours and select the ones that contains the
	// vertex for which the neighbours are requested. On un-balanced trees
	// the vertex can be inside the face/edge of the neighbour (hanging nodes).
	// To correctly consider these neighbours, the following logic can be
	// used:
	//   - if a face/edge neighbour has the same level or a lower level than
	//     the current cell, then it certainly is also a vertex neighbour;
	//   - if a face/edge neighbour has a higher level than the current cell,
	//     it is necessary to check if the neighbour actually contains the
	//     vertex.
	//
	// NOTE: in three dimension the function "_findCellEdgeNeighs" will return
	// both edge and face neighbours.
	const OctantInfo octantInfo = getCellOctant(id);
	const Octant *octant = getOctantPointer(octantInfo);
	int octantLevel = m_tree->getLevel(octant);
	if (isThreeDimensional()) {
		std::vector<long> edgeNeighs;
		for (int edge : m_octantLocalEdgesOnVertex[vertex]) {
			edgeNeighs.clear();
			_findCellEdgeNeighs(id, edge, blackList, &edgeNeighs);
			for (long neighId : edgeNeighs) {
				const OctantInfo neighOctantInfo = getCellOctant(neighId);
				const Octant *neighOctant = getOctantPointer(neighOctantInfo);
				int neighOctantLevel = m_tree->getLevel(neighOctant);
				if (neighOctantLevel <= octantLevel) {
					utils::addToOrderedVector<long>(neighId, *neighs);
				} else if (m_tree->isNodeOnOctant(octant, vertex, neighOctant)) {
					utils::addToOrderedVector<long>(neighId, *neighs);
				}
			}
		}
	} else {
		std::vector<long> faceNeighs;
		for (int face : m_octantLocalFacesOnVertex[vertex]) {
			faceNeighs.clear();
			_findCellFaceNeighs(id, face, blackList, &faceNeighs);
			for (long neighId : faceNeighs) {
				const OctantInfo neighOctantInfo = getCellOctant(neighId);
				const Octant *neighOctant = getOctantPointer(neighOctantInfo);
				int neighOctantLevel = m_tree->getLevel(neighOctant);
				if (neighOctantLevel <= octantLevel) {
					utils::addToOrderedVector<long>(neighId, *neighs);
				} else if (m_tree->isNodeOnOctant(octant, vertex, neighOctant)) {
					utils::addToOrderedVector<long>(neighId, *neighs);
				}
			}
		}
	}
}

/*!
	Finds the neighbours for the given co-dimension of the specified cell.

	Only the neighbours for the specified co-dimension are found, neighbours
	of higher co-dimensions are not inserted in the returned list.

	\param id is the id of the cell
	\param codimension is the co-dimension
	\param index is the local index of the entity (vertex, edge or face)
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\result The neighbours for the given codimension of the specified cell.
*/
void VolOctree::findCellCodimensionNeighs(const long &id, const int &index,
                                          const int &codimension, const std::vector<long> &blackList,
                                          std::vector<long> *neighs) const
{
	int dimension = getDimension();
	if (codimension > dimension || codimension <= 0) {
		return;
	}

	OctantInfo octantInfo = getCellOctant(id);

	std::vector<uint32_t> neighTreeIds;
	std::vector<bool> neighGhostFlags;
	if (octantInfo.internal) {
		m_tree->findNeighbours(octantInfo.id, index, codimension, neighTreeIds, neighGhostFlags);
	} else {
		m_tree->findGhostNeighbours(octantInfo.id, index, codimension, neighTreeIds, neighGhostFlags);
	}

	int nNeighs = neighTreeIds.size();
	for (int i = 0; i < nNeighs; ++i) {
		OctantInfo neighOctantInfo(neighTreeIds[i], !neighGhostFlags[i]);
		long neighId = getOctantId(neighOctantInfo);

		if (utils::findInOrderedVector<long>(neighId, blackList) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
}

}
