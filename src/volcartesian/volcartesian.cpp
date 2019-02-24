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

#include <cmath>
#include <bitset>

#include "bitpit_common.hpp"

#include "volcartesian.hpp"

namespace bitpit {

/*!
	\class VolCartesian
	\ingroup volumepatches

	\brief The VolCartesian defines a Cartesian patch.

	VolCartesian defines a Cartesian patch.
*/

/*!
	Creates an uninitialized patch.
*/
VolCartesian::VolCartesian()
	: VolumeKernel(false)
{
	initialize();
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param lengths are the lengths of the domain
	\param nCells are the numbers of cells of the patch
*/
VolCartesian::VolCartesian(const int &dimension,
                               const std::array<double, 3> &origin,
                               const std::array<double, 3> &lengths,
                               const std::array<int, 3> &nCells)
	: VolCartesian(PatchManager::AUTOMATIC_ID, dimension, origin, lengths, nCells)
{
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param lengths are the lengths of the domain
	\param nCells are the number of cells along each direction
*/
VolCartesian::VolCartesian(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               const std::array<double, 3> &lengths,
                               const std::array<int, 3> &nCells)
	: VolumeKernel(id, dimension, false)
{
	initialize();

	setOrigin(origin);
	setLengths(lengths);

	setDiscretization(nCells);
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param nCells is the number of cells along each direction
*/
VolCartesian::VolCartesian(const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, int nCells)
	: VolCartesian(PatchManager::AUTOMATIC_ID, dimension, origin, length, nCells)
{
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param nCells is the number of cells along each direction
*/
VolCartesian::VolCartesian(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, int nCells)
	: VolumeKernel(id, dimension, false)
{
	initialize();

	setOrigin(origin);
	setLengths({{length, length, length}});

	setDiscretization({{nCells, nCells, nCells}});
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed mesh spacing
*/
VolCartesian::VolCartesian(const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, double dh)
	: VolCartesian(PatchManager::AUTOMATIC_ID, dimension, origin, length, dh)
{
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed mesh spacing
*/
VolCartesian::VolCartesian(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, double dh)
	: VolumeKernel(id, dimension, false)
{
	initialize();

	setOrigin(origin);
	setLengths({{length, length, length}});

	int nCells = (int) std::ceil(length / dh);
	setDiscretization({{nCells, nCells, nCells}});
}

/*!
	Creates a new patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
VolCartesian::VolCartesian(std::istream &stream)
	: VolumeKernel(false)
{
	initialize();
	restore(stream);
}

/*!
	Creates a clone of the pach.

	\result A clone of the pach.
*/
std::unique_ptr<PatchKernel> VolCartesian::clone()  const
{
	return std::unique_ptr<VolCartesian>(new VolCartesian(*this));
}

/*!
	Function to reset the patch.
*/
void VolCartesian::reset()
{
	// Switch to light memeory mode
	switchMemoryMode(MEMORY_LIGHT);

	// Reset geometry and discretization
	m_minCoords = {{0., 0., 0.}};
	m_maxCoords = {{1., 1., 1.}};

	m_nCells1D     = {{0, 0, 0}};
	m_nVertices1D  = {{0, 0, 0}};
	m_cellSpacings = {{0., 0., 0.}};

	for (int n = 0; n < 3; ++n) {
		std::vector<double>().swap(m_vertexCoords[n]);
		std::vector<double>().swap(m_cellCenters[n]);
	}

	m_nVertices   = 0;
	m_nCells      = 0;
	if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
		m_nInterfaces = 0;
	}
}

/*!
	Resest the interfaces of the patch.
*/
void VolCartesian::resetInterfaces()
{
	PatchKernel::resetInterfaces();

	m_nInterfaces = 0;
}

/*!
	Build interfaces among the cells.
*/
void VolCartesian::buildInterfaces()
{
	// Reset interfaces
	if (getInterfacesBuildStrategy() != INTERFACES_NONE) {
		clearInterfaces();
	}

	// Count the total number of interfaces
	m_nInterfaces = 0;
	for (int n = 0; n < getDimension(); n++) {
		std::array<int, 3> interfaceCount1D = getInterfaceCountDirection(n);

		int nDirectionInterfaces = 1;
		for (int n = 0; n < getDimension(); n++) {
			nDirectionInterfaces *= interfaceCount1D[n];
		}
		m_nInterfaces += nDirectionInterfaces;
	}

	// Build interfaces
	if (getMemoryMode() == MemoryMode::MEMORY_NORMAL) {
		// Enable advanced editing
		setExpert(true);

		// Update interfaces
		addInterfaces();

		// Disable advanced editing
		setExpert(false);
	} else {
		// Set interfaces build strategy
		setInterfacesBuildStrategy(INTERFACES_AUTOMATIC);
	}
}

/*!
	Update the interfaces.

	Although the function receives a list of cells to update, it is not
	possible to partially update the interfaces. The function will always
	update all the interfaces.

	\param[in] cellIds is the list of cell ids
*/
void VolCartesian::updateInterfaces(const std::vector<long> &cellIds)
{
	if (cellIds.size() != (std::size_t) m_nCells) {
		log::cout() << " It is not possible to partially update the interfaces.";
		log::cout() << " All interface will be updated.";
	}

	buildInterfaces();
}

/*!
	Initialize the data structures of the patch.
*/
void VolCartesian::initialize()
{
	// Normals
	int i = 0;
	for (int n = 0; n < 3; n++) {
		for (int k = -1; k <= 1; k += 2) {
			std::array<double, 3> normal = {{0.0, 0.0, 0.0}};
			normal[n] = k;

			m_normals[i++] = normal;
		}
	}

	// Deltas for the evaluation of the vertex neighbours
	m_vertexNeighDeltas = std::vector<std::array<int, 3>>(8);
	m_vertexNeighDeltas[0] = {{ 0,  0, 0}};
	m_vertexNeighDeltas[1] = {{-1,  0, 0}};
	m_vertexNeighDeltas[2] = {{ 0, -1, 0}};
	m_vertexNeighDeltas[3] = {{-1, -1, 0}};
	m_vertexNeighDeltas[4] = {{ 0,  0, -1}};
	m_vertexNeighDeltas[5] = {{-1,  0, -1}};
	m_vertexNeighDeltas[6] = {{ 0, -1, -1}};
	m_vertexNeighDeltas[7] = {{-1, -1, -1}};

	// Deltas for the evaluation of the edge neighbours
	m_edgeNeighDeltas = std::vector<std::array<int, 3>>(12);
	m_edgeNeighDeltas[ 0] = {{-1,  0, -1}};
	m_edgeNeighDeltas[ 1] = {{ 1,  0, -1}};
	m_edgeNeighDeltas[ 2] = {{ 0, -1, -1}};
	m_edgeNeighDeltas[ 3] = {{ 0,  1, -1}};
	m_edgeNeighDeltas[ 4] = {{-1, -1,  0}};
	m_edgeNeighDeltas[ 5] = {{ 1, -1,  0}};
	m_edgeNeighDeltas[ 6] = {{-1,  1,  0}};
	m_edgeNeighDeltas[ 7] = {{ 1,  1,  0}};
	m_edgeNeighDeltas[ 8] = {{-1,  0,  1}};
	m_edgeNeighDeltas[ 9] = {{ 1,  0,  1}};
	m_edgeNeighDeltas[10] = {{ 0, -1,  1}};
	m_edgeNeighDeltas[11] = {{ 0,  1,  1}};

	// Faces associated to the edges
	m_edgeFaces = std::vector<std::array<int, 2>>(12);
	m_edgeFaces[ 0] = {{ 0, 4}};
	m_edgeFaces[ 1] = {{ 1, 4}};
	m_edgeFaces[ 2] = {{ 2, 4}};
	m_edgeFaces[ 3] = {{ 3, 4}};
	m_edgeFaces[ 4] = {{ 0, 2}};
	m_edgeFaces[ 5] = {{ 1, 2}};
	m_edgeFaces[ 6] = {{ 0, 3}};
	m_edgeFaces[ 7] = {{ 1, 3}};
	m_edgeFaces[ 8] = {{ 0, 5}};
	m_edgeFaces[ 9] = {{ 1, 5}};
	m_edgeFaces[10] = {{ 2, 5}};
	m_edgeFaces[11] = {{ 3, 5}};

	// Set the bounding box as frozen
	setBoundingBoxFrozen(true);

	// This patch need to be spawn
	setSpawnStatus(SPAWN_NEEDED);

	// Set adjacencies build strategy
	setAdjacenciesBuildStrategy(ADJACENCIES_AUTOMATIC);

	// Set the light memory mode
	setMemoryMode(MemoryMode::MEMORY_LIGHT);

	// Reset the patch
	reset();
}

/*!
	Initializes cell volume
*/
void VolCartesian::initializeCellVolume()
{
	m_cellVolume = m_cellSpacings[Vertex::COORD_X] * m_cellSpacings[Vertex::COORD_Y];
	if (isThreeDimensional()) {
		m_cellVolume *= m_cellSpacings[Vertex::COORD_Z];
	}
}

/*!
	Initializes interface area
*/
void VolCartesian::initializeInterfaceArea()
{
	for (int n = 0; n < getDimension(); ++n) {
		m_interfaceArea[n] = m_cellVolume / m_cellSpacings[n];
	}

	if (!isThreeDimensional()) {
		m_interfaceArea[Vertex::COORD_Z] = 0.;
	}
}

/*!
	Discretizes the domain.

	\param nCells is the numbers of cells along each direction
*/
void VolCartesian::setDiscretization(const std::array<int, 3> &nCells)
{
	// Spacing
	for (int n = 0; n < getDimension(); ++n) {
		// Initialize cells
		m_nCells1D[n]     = nCells[n];
		m_cellSpacings[n] = (m_maxCoords[n] - m_minCoords[n]) / m_nCells1D[n];

		m_cellCenters[n].resize(m_nCells1D[n]);
		for (int i = 0; i < m_nCells1D[n]; i++) {
			m_cellCenters[n][i] = m_minCoords[n] + (0.5 + i) * m_cellSpacings[n];
		}

		log::cout() << "  - Cell count along direction " << n << " : " << m_nCells1D[n] << "\n";

		// Initialize vertices
		m_nVertices1D[n] = m_nCells1D[n] + 1;

		m_vertexCoords[n].resize(m_nVertices1D[n]);
		for (int i = 0; i < m_nVertices1D[n]; i++) {
			m_vertexCoords[n][i] = m_minCoords[n] + i * m_cellSpacings[n];
		}
	}

	if (!isThreeDimensional()) {
		m_nCells1D[Vertex::COORD_Z]     = 0;
		m_cellSpacings[Vertex::COORD_Z] = 0.;

		m_nVertices1D[Vertex::COORD_Z] = 0;
	}

	log::cout() << std::endl;

	// Count the total number of vertices
	m_nVertices = 1;
	for (int n = 0; n < getDimension(); ++n) {
		m_nVertices *= m_nVertices1D[n];
	}
	log::cout() << "  - Total vertex count: " << m_nVertices << "\n";

	// Count the total number of cells
	m_nCells = 1;
	for (int n = 0; n < getDimension(); ++n) {
		m_nCells    *= m_nCells1D[n];
	}
	log::cout() << "  - Total cell count: " << m_nCells << "\n";

	// Cell volume
	initializeCellVolume();

	// Interface area
	initializeInterfaceArea();

	// Update patch data structures
	if (getMemoryMode() == MemoryMode::MEMORY_NORMAL) {
		// Switch to light mode to reset patchdata structures
		switchMemoryMode(MemoryMode::MEMORY_LIGHT);

		// Swich back to normal mode to rebuild patchdata structures
		switchMemoryMode(MemoryMode::MEMORY_NORMAL);
	}
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long VolCartesian::getVertexCount() const
{
	return m_nVertices;
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long VolCartesian::getCellCount() const
{
	return m_nCells;
}

/*!
	Gets the element type for the cell with the specified id.

	\param id is the id of the requested cell
	\return The element type for the cell with the specified id.
*/
ElementType VolCartesian::getCellType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getCellType();
}

/*!
	Gets the element type for the cells in the patch.

	\return The element type for the cells in the patch.
*/
ElementType VolCartesian::getCellType() const
{
	if (isThreeDimensional()) {
		return ElementType::VOXEL;
	} else {
		return ElementType::PIXEL;
	}
}

/*!
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
long VolCartesian::getInterfaceCount() const
{
	return m_nInterfaces;
}

/*!
	Gets the element type for the interface with the specified id.

	\param id is the id of the requested interface
	\return The element type for the interface with the specified id.
*/
ElementType VolCartesian::getInterfaceType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getInterfaceType();
}

/*!
	Gets the element type for the interfaces in the patch.

	\return The element type for the interfaces in the patch.
*/
ElementType VolCartesian::getInterfaceType() const
{
	if (isThreeDimensional()) {
		return ElementType::PIXEL;
	} else {
		return ElementType::LINE;
	}
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolCartesian::evalCellVolume(const long &id) const
{
	BITPIT_UNUSED(id);

	return m_cellVolume;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolCartesian::evalCellSize(const long &id) const
{
	BITPIT_UNUSED(id);

	return pow(m_cellVolume, 1. / getDimension());
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double VolCartesian::evalInterfaceArea(const long &id) const
{
	const Interface &interface = getInterface(id);
	int ownerFace = interface.getOwnerFace();
	int direction = std::floor(ownerFace / 2.);

	return m_interfaceArea[direction];
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> VolCartesian::evalInterfaceNormal(const long &id) const
{
	const Interface &interface = getInterface(id);
	int ownerFace = interface.getOwnerFace();

	return m_normals[ownerFace];
}

/*!
	Get cell spacings of the patch.

	\result Cell spacings of the patch.
*/
std::array<double, 3> VolCartesian::getSpacing() const
{
	return m_cellSpacings;
}

/*!
	Switch to the specified memory mode.

	\param mode is the memory mode that will be set
*/
void VolCartesian::switchMemoryMode(MemoryMode mode)
{
	if (mode == getMemoryMode()) {
		return;
	}

	// Update the data structures
	switch (mode) {

	case MemoryMode::MEMORY_NORMAL:
		// Spawn the patch to activate normal memory mode
		spawn(false);

		break;

	case MemoryMode::MEMORY_LIGHT:
		// To put the patch in memory mode we need to reset the generic data
		// of the patch, therefore we can call the 'reset' implementation of
		// the kernel.
		VolumeKernel::reset();

		// Now the patch needs to be spawn
		setSpawnStatus(SPAWN_NEEDED);

		// Set the light memory mode
		setMemoryMode(mode);

		break;

	}
}

/*!
	Function to set the memory mode flag.

	This function just sets the flag to the specified value.

	\param mode is the memory mode that will be set
*/
void VolCartesian::setMemoryMode(MemoryMode mode)
{
	if (mode == m_memoryMode) {
		return;
	}

	m_memoryMode = mode;
}

/*!
	Get the current memory mode.

	\result The current memory mode.
*/
VolCartesian::MemoryMode VolCartesian::getMemoryMode()
{
	return m_memoryMode;
}

/*!
	Get cell spacing along the specificed direction.

	\param[in] direction is the direction along which the spacing is
	requested
	\result The cell spacing along the specificed direction.
*/
double VolCartesian::getSpacing(const int &direction) const
{
	return m_cellSpacings[direction];
}

/*!
	Generates the patch.

	\param trackSpawn if set to true the changes to the patch will be tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
std::vector<adaption::Info> VolCartesian::_spawn(bool trackSpawn)
{
	std::vector<adaption::Info> updateInfo;

	// If the patch is in 'normal' mode there is nothing to do.
	if (getMemoryMode() == MEMORY_NORMAL) {
		return updateInfo;
	}

	// Enable advanced editing
	setExpert(true);

	// Definition of the mesh
	addVertices();
	addCells();
	if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
		addInterfaces();
	}

	// Disable advanced editing
	setExpert(false);

	// Adaption info
	if (trackSpawn) {
		updateInfo.emplace_back();
		adaption::Info &adaptionCellInfo = updateInfo.back();
		adaptionCellInfo.type   = adaption::TYPE_CREATION;
		adaptionCellInfo.entity = adaption::ENTITY_CELL;
		adaptionCellInfo.current.reserve(m_cells.size());
		for (auto &cell : m_cells) {
			adaptionCellInfo.current.emplace_back();
			long &cellId = adaptionCellInfo.current.back();
			cellId = cell.getId();
		}

		updateInfo.emplace_back();
		adaption::Info &adaptionInterfaceInfo = updateInfo.back();
		adaptionInterfaceInfo.type   = adaption::TYPE_CREATION;
		adaptionInterfaceInfo.entity = adaption::ENTITY_INTERFACE;
		adaptionInterfaceInfo.current.reserve(m_interfaces.size());
		for (auto &interface : m_interfaces) {
			adaptionInterfaceInfo.current.emplace_back();
			long &interfaceId = adaptionInterfaceInfo.current.back();
			interfaceId = interface.getId();
		}
	} else {
		updateInfo.emplace_back();
	}

	// Updating the adaption brings the patch is in normal memory mode
	setMemoryMode(MemoryMode::MEMORY_NORMAL);

	// Done
	return updateInfo;
}

/*!
	Creates the vertices of the patch.
*/
void VolCartesian::addVertices()
{
	log::cout() << "  >> Creating vertices\n";

	log::cout() << "    - Vertex count: " << m_nVertices << "\n";

	m_vertices.reserve(m_nVertices);
	for (int k = 0; (isThreeDimensional()) ? (k < m_nVertices1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (int j = 0; j < m_nVertices1D[Vertex::COORD_Y]; j++) {
			for (int i = 0; i < m_nVertices1D[Vertex::COORD_X]; i++) {
				// Linear index of the vertex
				long id_vertex = getVertexLinearId(i, j, k);

				// Vertex coordinates
				std::array<double, 3> coords = evalVertexCoords(id_vertex);

				// Add vertex
				VolumeKernel::addVertex(std::move(coords), id_vertex);
			}
		}
	}
}

/*!
	Evaluates the coordinate of the specified vertex.

	\param id is the id of the vertex
	\result The coordinate of the specified vertex.
*/
std::array<double, 3> VolCartesian::evalVertexCoords(const long &id)
{
	std::array<int, 3> ijk = getVertexCartesianId(id);

	std::array<double, 3> coords;
	coords[Vertex::COORD_X] = m_vertexCoords[Vertex::COORD_X][ijk[0]];
	coords[Vertex::COORD_Y] = m_vertexCoords[Vertex::COORD_Y][ijk[1]];
	if (isThreeDimensional()) {
		coords[Vertex::COORD_Z] = m_vertexCoords[Vertex::COORD_Z][ijk[2]];
	} else {
		coords[Vertex::COORD_Z] = 0.0;
	}

	return coords;
}

/*!
	Creates the cells of the patch.
*/
void VolCartesian::addCells()
{
	log::cout() << "  >> Creating cells\n";

	// Info on the cells
	ElementType cellType = getCellType();

	// Create the cells
	log::cout() << "    - Cell count: " << m_nCells << "\n";

	m_cells.reserve(m_nCells);
	for (int k = 0; (isThreeDimensional()) ? (k < m_nCells1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; j++) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; i++) {
				long id_cell = getCellLinearId(i, j, k);
				CellIterator cellIterator = VolumeKernel::addCell(cellType, true, id_cell);
				Cell &cell = *cellIterator;

				// Connettività
				long *cellConnect = cell.getConnect();

				cellConnect[0] = getVertexLinearId(i,     j,     k);
				cellConnect[1] = getVertexLinearId(i + 1, j,     k);
				cellConnect[2] = getVertexLinearId(i,     j + 1, k);
				cellConnect[3] = getVertexLinearId(i + 1, j + 1, k);
				if (isThreeDimensional()) {
					cellConnect[4] = getVertexLinearId(i,     j,     k + 1);
					cellConnect[5] = getVertexLinearId(i + 1, j,     k + 1);
					cellConnect[6] = getVertexLinearId(i,     j + 1, k + 1);
					cellConnect[7] = getVertexLinearId(i + 1, j + 1, k + 1);
				}
			}
		}
	}
}

/*!
	Creates the interfaces of the patch.
*/
void VolCartesian::addInterfaces()
{
	log::cout() << "  >> Creating interfaces\n";

	log::cout() << "    - Interface count: " << m_nInterfaces << "\n";

	// Create interfaces
	m_interfaces.reserve(m_nInterfaces);
	for (int n = 0; n < getDimension(); n++) {
		addInterfacesDirection(n);
	}

	// Set interfaces build strategy
	setInterfacesBuildStrategy(INTERFACES_AUTOMATIC);
}

/*!
	Get the interface ount for the given direction.

	\param direction the method will count the interfaces normal to this
	                 direction
	\result The interface count for the given direction.
*/
std::array<int, 3> VolCartesian::getInterfaceCountDirection(const int &direction)
{
	std::array<int, 3> interfaceDirectionCount = m_nCells1D;
	interfaceDirectionCount[direction]++;

	return interfaceDirectionCount;
}

/*!
	Creates the interfaces normal to the given direction.

	\param direction the method will creat the interfaces normal to this
	                 direction
*/
void VolCartesian::addInterfacesDirection(const int &direction)
{
	log::cout() << "  >> Creating interfaces normal to direction " << direction << "\n";

	// Info on the interfaces
	ElementType interfaceType = getInterfaceType();

	const ReferenceElementInfo &interfaceTypeInfo = ReferenceElementInfo::getInfo(interfaceType);
	const int nInterfaceVertices = interfaceTypeInfo.nVertices;
	std::array<int, 3> interfaceCount1D = getInterfaceCountDirection(direction);

	// Counters
	std::array<int, 3> counters = {{0, 0, 0}};
	int &i = counters[Vertex::COORD_X];
	int &j = counters[Vertex::COORD_Y];
	int &k = counters[Vertex::COORD_Z];

	// Creation of the interfaces
	for (k = 0; (isThreeDimensional()) ? (k < interfaceCount1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (j = 0; j < interfaceCount1D[Vertex::COORD_Y]; j++) {
			for (i = 0; i < interfaceCount1D[Vertex::COORD_X]; i++) {
				InterfaceIterator interfaceIterator = VolumeKernel::addInterface(interfaceType);
				Interface &interface = *interfaceIterator;

				// Owner id
				std::array<int, 3> ownerIJK(counters);
				if (counters[direction] == (interfaceCount1D[direction] - 1)) {
					ownerIJK[direction] -= 1;
				}

				long ownerId = getCellLinearId(ownerIJK);

				// Neighbour id
				long neighId;
				if (counters[direction] != 0 && counters[direction] != interfaceCount1D[direction] - 1) {
					std::array<int, 3> neighIJK(counters);
					neighIJK[direction] -= 1;

					neighId = getCellLinearId(neighIJK);
				} else {
					neighId = Element::NULL_ID;
				}

				// Owner data
				Cell &owner = m_cells[ownerId];

				int ownerFace = 2 * direction;
				if (counters[direction] == interfaceCount1D[direction] - 1) {
					ownerFace++;
				}

				interface.setOwner(owner.getId(), ownerFace);
				owner.setInterface(ownerFace, 0, interface.getId());
				owner.setAdjacency(ownerFace, 0, neighId);

				// Neighbour data
				if (counters[direction] != 0 && counters[direction] != interfaceCount1D[direction] - 1) {
					Cell &neigh = m_cells[neighId];

					int neighFace = 2 * direction + 1;

					interface.setNeigh(neigh.getId(), neighFace);
					neigh.setInterface(neighFace, 0, interface.getId());
					neigh.setAdjacency(neighFace, 0, ownerId);
				} else {
					interface.unsetNeigh();
				}

				// Connectivity
				std::unique_ptr<long[]> connect = std::unique_ptr<long[]>(new long[nInterfaceVertices]);
				if (direction == Vertex::COORD_X) {
					connect[0] = getVertexLinearId(i, j,     k);
					connect[1] = getVertexLinearId(i, j + 1, k);
					if (isThreeDimensional()) {
						connect[2] = getVertexLinearId(i, j + 1, k + 1);
						connect[3] = getVertexLinearId(i, j,     k + 1);
					}
				} else if (direction == Vertex::COORD_Y) {
					connect[0] = getVertexLinearId(i,     j,     k);
					connect[1] = getVertexLinearId(i + 1, j,     k);
					if (isThreeDimensional()) {
						connect[2] = getVertexLinearId(i + 1, j, k + 1);
						connect[3] = getVertexLinearId(i,     j, k + 1);
					}
				} else if (direction == Vertex::COORD_Z) {
					connect[0] = getVertexLinearId(i,     j,     k);
					connect[1] = getVertexLinearId(i + 1, j,     k);
					if (isThreeDimensional()) {
						connect[2] = getVertexLinearId(i + 1, j + 1, k);
						connect[3] = getVertexLinearId(i,     j + 1, k);
					}
				}

				interface.setConnect(std::move(connect));
			}
		}
	}
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int VolCartesian::_getDumpVersion() const
{
	const int DUMP_VERSION = 2;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolCartesian::_dump(std::ostream &stream) const
{
	utils::binary::write(stream, m_minCoords[0]);
	utils::binary::write(stream, m_minCoords[1]);
	utils::binary::write(stream, m_minCoords[2]);

	utils::binary::write(stream, m_maxCoords[0] - m_minCoords[0]);
	utils::binary::write(stream, m_maxCoords[1] - m_minCoords[1]);
	utils::binary::write(stream, m_maxCoords[2] - m_minCoords[2]);

	utils::binary::write(stream, m_nCells1D[0]);
	utils::binary::write(stream, m_nCells1D[1]);
	utils::binary::write(stream, m_nCells1D[2]);

	utils::binary::write(stream, m_memoryMode);

	if (m_memoryMode == MEMORY_NORMAL) {
		for (const Cell &cell : getCells()) {
			utils::binary::write(stream, cell.getPID());
		}
	}

	if (m_memoryMode == MEMORY_NORMAL) {
		for (const Interface &interface : getInterfaces()) {
			utils::binary::write(stream, interface.getPID());
		}
	}
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void VolCartesian::_restore(std::istream &stream)
{
	// This patch need to be spawn
	setSpawnStatus(SPAWN_NEEDED);

	// Origin
	std::array<double, 3> origin;
	utils::binary::read(stream, origin[0]);
	utils::binary::read(stream, origin[1]);
	utils::binary::read(stream, origin[2]);

	setOrigin(origin);

	// Lengths
	std::array<double, 3> lengths;
	utils::binary::read(stream, lengths[0]);
	utils::binary::read(stream, lengths[1]);
	utils::binary::read(stream, lengths[2]);

	setLengths(lengths);

	// Discretization
	std::array<int, 3> nCells;
	utils::binary::read(stream, nCells[0]);
	utils::binary::read(stream, nCells[1]);
	utils::binary::read(stream, nCells[2]);

	setDiscretization(nCells);

	// Memory mode
	MemoryMode memoryMode;
	utils::binary::read(stream, memoryMode);
	switchMemoryMode(memoryMode);

	// Cell data
	if (m_memoryMode == MEMORY_NORMAL) {
		for (Cell &cell : getCells()) {
			int PID;
			utils::binary::read(stream, PID);
			cell.setPID(PID);
		}
	}

	// Interface data
	if (m_memoryMode == MEMORY_NORMAL) {
		for (Interface &interface : getInterfaces()) {
			int PID;
			utils::binary::read(stream, PID);
			interface.setPID(PID);
		}
	}
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the idof the cell
	\param[in] point is the point to be checked
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolCartesian::isPointInside(const long &id, const std::array<double, 3> &point)
{
	std::array<int, 3> cellIjk = getCellCartesianId(id);

	const double EPS = getTol();
    for (int d = 0; d < 3; ++d){
		long index = cellIjk[d];

		double spacing = m_cellSpacings[d];
		double cellMin = m_cellCenters[d][index] - 0.5 * spacing;
		double cellMax = m_cellCenters[d][index] + 0.5 * spacing;

        if (point[d]< cellMin - EPS || point[d] > cellMax + EPS) {
            return false;
        }
    }

    return true;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] point is the point to be checked
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolCartesian::isPointInside(const std::array<double, 3> &point)
{
	const double EPS = getTol();

	for (int n = 0; n < getDimension(); n++) {
		if (point[n] < m_minCoords[n] - EPS || point[n] > m_maxCoords[n] + EPS) {
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
	\result Returns the linear id of the cell the contains the point. If the
	point is not inside the patch, the function returns the id of the null
	element.
*/
long VolCartesian::locatePoint(const std::array<double, 3> &point)
{
	std::array<int, 3> pointIjk = locatePointCartesian(point);
	if (isCellCartesianIdValid(pointIjk)) {
		return getCellLinearId(pointIjk);
	} else {
		return Element::NULL_ID;
	}
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] point is the point to be checked
	\result Returns the set of cartesian id of the cell the contains the point.
	If the point is not inside the patch, the function returns negative ids.
*/
std::array<int, 3> VolCartesian::locatePointCartesian(const std::array<double, 3> &point)
{
	std::array<int, 3> ijk;
	if (!isPointInside(point)) {
		ijk[0] = -1;
		ijk[1] = -1;
		ijk[2] = -1;
	}

	ijk[0] = std::floor((point[Vertex::COORD_X] - m_minCoords[Vertex::COORD_X]) / m_cellSpacings[Vertex::COORD_X]);
	ijk[1] = std::floor((point[Vertex::COORD_Y] - m_minCoords[Vertex::COORD_Y]) / m_cellSpacings[Vertex::COORD_Y]);
	if (isThreeDimensional()) {
		ijk[2] = std::floor((point[Vertex::COORD_Z] - m_minCoords[Vertex::COORD_Z]) / m_cellSpacings[Vertex::COORD_Z]);
	} else {
		ijk[2] = -1;
	}

	return ijk;
}

/*!
	Locates the closest vertex of the given point.

	\param[in] point is the point
	\result The linear id of the closest vertex of the given point.
*/
long VolCartesian::locateClosestVertex(std::array<double,3> const &point) const
{
	return getVertexLinearId(locateClosestVertexCartesian(point));
}

/*!
	Locates the closest vertex of the given point.

	\param[in] point is the point
	\result The set of Cartesian id of the closest vertex of the given
	point.
*/
std::array<int, 3> VolCartesian::locateClosestVertexCartesian(std::array<double,3> const &point) const
{
	std::array<int, 3> ijk;
	for (int n = 0; n < getDimension(); ++n) {
		ijk[n] = std::min(m_nVertices1D[n] - 1, std::max(0, (int) round((point[n] - m_minCoords[n]) / m_cellSpacings[n])));
	}

	if (!isThreeDimensional()) {
		ijk[Vertex::COORD_Z] = -1;
	}

	return ijk;
}

/*!
	Locates the closest cell to point.

	If the point is not inside the patch, the function returns the id of the
	cell the point projects on.

	\param[in] point is the point to be checked
	\result Returns the linear id of the closest cell to the point.
*/
long VolCartesian::locateClosestCell(const std::array<double, 3> &point)
{
    std::array<int, 3> pointIjk = locateClosestCellCartesian(point);
    return getCellLinearId(pointIjk);
}

/*!
	Locates the cell the closest cell to the point.

	If the point is not inside the patch, the function returns the id of the
	cell the point projects on.

	\param[in] point is the point to be checked
	\result Returns the cartesian id of the closest cell to the point.
*/
std::array<int, 3> VolCartesian::locateClosestCellCartesian(const std::array<double, 3> &point)
{
	std::array<int,3> ijk({{0,0,0}});

    for( int i=0; i<getDimension(); ++i){
	    ijk[i] = std::floor((point[i] - m_minCoords[i]) / m_cellSpacings[i]);
        ijk[i] = std::max( ijk[i], 0 );
        ijk[i] = std::min( ijk[i], m_nCells1D[i]-1 );
    }

	return ijk;
}

/*!
	Converts the cell cartesian notation to a linear notation
*/
long VolCartesian::getCellLinearId(const int &i, const int &j, const int &k) const
{
	long id = i;
	id += m_nCells1D[Vertex::COORD_X] * j;
	if (getDimension() == 3) {
		id += m_nCells1D[Vertex::COORD_Y] * m_nCells1D[Vertex::COORD_X] * k;
	}

	return id;
}

/*!
	Converts the cell cartesian notation to a linear notation
*/
long VolCartesian::getCellLinearId(const std::array<int, 3> &ijk) const
{
	return getCellLinearId(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}

/*!
	Converts a cell linear index to a set of cartesian indices.

	No check on bounds is performed.

	\param[in] idx is the linear index of the cell
	\result Returns the set of cartesian indices of the cell.
*/
std::array<int, 3> VolCartesian::getCellCartesianId(long const &idx) const
{
	int offset_ij = m_nCells1D[0] * m_nCells1D[1];

	std::array<int, 3> id;
	id[2] = idx / offset_ij;
	id[1] = (idx - id[2] * offset_ij) / m_nCells1D[0];
	id[0] = idx - (idx / m_nCells1D[0]) * m_nCells1D[0];

	// Set to -1 the z-index for 2D patches
	unsigned int isThreeDimensionalFlag = isThreeDimensional();
	id[2] = id[2] * isThreeDimensionalFlag - !isThreeDimensionalFlag;

	return id;
}

/*!
	Checks if a cell cartesian index is valid.

	\param ijk is the set of cartesian indices of the cell
	\result Returns true if the index is valid, otherwise it returns false.
*/
bool VolCartesian::isCellCartesianIdValid(const std::array<int, 3> &ijk) const
{
	for (int k = 0; k < getDimension(); ++k) {
		if (ijk[k] < 0 || ijk[k] >= m_nCells1D[k]) {
			return false;
		}
	}

	return true;
}

/*!
	Converts the vertex cartesian notation to a linear notation
*/
long VolCartesian::getVertexLinearId(const int &i, const int &j, const int &k) const
{
	long id = i;
	id += m_nVertices1D[Vertex::COORD_X] * j;
	if (getDimension() == 3) {
		id += m_nVertices1D[Vertex::COORD_Y] * m_nVertices1D[Vertex::COORD_X] * k;
	}

	return id;
}

/*!
	Converts the vertex cartesian notation to a linear notation
*/
long VolCartesian::getVertexLinearId(const std::array<int, 3> &ijk) const
{
	return getVertexLinearId(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}

/*!
	Converts a vertex linear index to a set of cartesian indices.

	No check on bounds is performed.

	\param[in] idx is the linear index of the vertex
	\result Returns the set of cartesian indices of the vertex.
*/
std::array<int, 3> VolCartesian::getVertexCartesianId(long const &idx) const
{
	int offset_ij = m_nVertices1D[0] * m_nVertices1D[1];

	std::array<int, 3> id;
	id[2] = idx / offset_ij;
	id[1] = (idx - id[2] * offset_ij) / m_nVertices1D[0];
	id[0] = idx - (idx / m_nVertices1D[0]) * m_nVertices1D[0];

	// Set to -1 the z-index for 2D patches
	unsigned int isThreeDimensionalFlag = isThreeDimensional();
	id[2] = id[2] * isThreeDimensionalFlag - !isThreeDimensionalFlag;

	return id;
}

/*!
	Gets the cartesian indices of the specified local vertex.

	No check on bounds is performed.

	\param[in] cellIdx is the linear cell index
	\param[in] vertex is the local vertex
	\result Returns the set of cartesian indices of the vertex.
*/
std::array<int, 3> VolCartesian::getVertexCartesianId(long const &cellIdx, int const &vertex) const
{
	return getVertexCartesianId(getCellCartesianId(cellIdx), vertex);
}

/*!
	Gets the cartesian indices of the specified local vertex.

	No check on bounds is performed.

	\param[in] cellIjk is the Cartesian cell index
	\param[in] vertex is the local vertex
	\result Returns the set of cartesian indices of the vertex.
*/
std::array<int, 3> VolCartesian::getVertexCartesianId(const std::array<int, 3> &cellIjk, int const &vertex) const
{
	std::bitset<3> vertexBitset(vertex);
	std::array<int, 3> vertexIjk(cellIjk);
	for (int k = 0; k < 3; ++k) {
		vertexIjk[k] += vertexBitset[k];
	}

	return vertexIjk;
}

/*!
	Checks if a vertex cartesian index is valid.

	\param ijk is the set of cartesian indices of the vertex
	\result Returns true if the index is valid, otherwise it returns false.
*/
bool VolCartesian::isVertexCartesianIdValid(const std::array<int, 3> &ijk) const
{
	for (int k = 0; k < getDimension(); ++k) {
		if (ijk[k] < 0 || ijk[k] >= m_nVertices1D[k]) {
			return false;
		}
	}

	return true;
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours will be stored.
	The vector is not cleared before adding the neighbours, it is extended
	by appending all the neighbours found by this function
*/
void VolCartesian::_findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	int neighSide      = face % 2;
	int neighDirection = std::floor(face / 2);

	std::array<int, 3> neighIjk(getCellCartesianId(id));
	if (neighSide == 0) {
		neighIjk[neighDirection]--;
	} else {
		neighIjk[neighDirection]++;
	}

	if (isCellCartesianIdValid(neighIjk)) {
		long neighId = getCellLinearId(neighIjk);
		if (utils::findInOrderedVector<long>(neighId, blackList) == blackList.end()) {
			neighs->push_back(neighId);
		}
	}
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
void VolCartesian::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return;
	}

	// Diagonal neighbour
	std::array<int, 3> diagNeighIjk(getCellCartesianId(id) + m_edgeNeighDeltas[edge]);
	if (isCellCartesianIdValid(diagNeighIjk)) {
		long diagNeighId = getCellLinearId(diagNeighIjk);
		if (utils::findInOrderedVector<long>(diagNeighId, blackList) == blackList.end()) {
			utils::addToOrderedVector<long>(diagNeighId, *neighs);
		}
	}

	// Faces incident to the edge
	for (int face : m_edgeFaces[edge]) {
		_findCellFaceNeighs(id, face, blackList, neighs);
	}
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param blackList is a list of cells that are excluded from the search.
	The blacklist has to be a unique list of ordered cell ids.
	\param[in,out] neighs is the vector were the neighbours of the specified
	cell for the given vertex will be stored. The vector is not cleared before
	adding the neighbours, it is extended by appending all the neighbours
	found by this function
*/
void VolCartesian::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList, std::vector<long> *neighs) const
{
	std::array<int, 3> cellIjk   = getCellCartesianId(id);
	std::array<int, 3> vertexIjk = getVertexCartesianId(cellIjk, vertex);

	for (const auto &delta : m_vertexNeighDeltas) {
		// Get the Cartesian index
		std::array<int, 3> neighIjk(vertexIjk + delta);
		if (neighIjk == cellIjk) {
			continue;
		} else if (!isCellCartesianIdValid(neighIjk)) {
			continue;
		}

		// Get the linear neighbour index and, if it's not on the blacklist,
		// add it to the list of neighbours
		long neighId = getCellLinearId(neighIjk);
		if (utils::findInOrderedVector<long>(neighId, blackList) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, *neighs);
		}
	}
}

/*!
	Extract a cell subset.

	\param[in] ijkMin is the set of cartesian indices of the lower bound
	\param[in] ijkMax is the set of cartesian indices of the upper bound
	\result The linear indices of the cell subset.
*/
std::vector<long> VolCartesian::extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax)
{
	int nSubsetCells_x = ijkMax[0] - ijkMin[0] + 1;
	int nSubsetCells_y = ijkMax[1] - ijkMin[1] + 1;
	int nSubsetCells_z = ijkMax[2] - ijkMin[2] + 1;

	std::vector<long> ids;
	ids.resize(nSubsetCells_x * nSubsetCells_y * nSubsetCells_z);

	std::vector<long>::iterator it = ids.begin();
	for (int k = ijkMin[2]; k <= ijkMax[2]; ++k) {
		for (int j = ijkMin[1]; j <= ijkMax[1]; ++j) {
			for (int i = ijkMin[0]; i <= ijkMax[0]; ++i) {
				*it = getCellLinearId(i, j, k);
				++it;
			}
		}
	}

	return ids;
}

/*!
	Extract a cell subset.

	\param[in] idxMin is the linear index of the lower bound
	\param[in] idxMax is the linear index of the upper bound
	\result The linear indices of the cell subset.
*/
std::vector<long> VolCartesian::extractCellSubSet(int const &idxMin, int const &idxMax)
{
	return extractCellSubSet(getCellCartesianId(idxMin), getCellCartesianId(idxMax));
}

/*!
	Extract a cell subset.

	\param[in] pointMin is the lower bound
	\param[in] pointMax is the upper bound
	\result The linear indices of the cell subset.
*/
std::vector<long> VolCartesian::extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax)
{
	return extractCellSubSet(locatePointCartesian(pointMin), locatePointCartesian(pointMax));
}

/*!
	Extract a vertex subset.

	\param[in] ijkMin is the set of cartesian indices of the lower bound
	\param[in] ijkMax is the set of cartesian indices of the upper bound
	\result The linear indices of the vertex subset.
*/
std::vector<long> VolCartesian::extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax)
{
	int nSubsetVertices_x = ijkMax[0] - ijkMin[0] + 1;
	int nSubsetVertices_y = ijkMax[1] - ijkMin[1] + 1;
	int nSubsetVertices_z = ijkMax[2] - ijkMin[2] + 1;

	std::vector<long> ids;
	ids.resize(nSubsetVertices_x * nSubsetVertices_y * nSubsetVertices_z);

	std::vector<long>::iterator it = ids.begin();
	for (int k = ijkMin[2]; k <= ijkMax[2]; ++k) {
		for (int j = ijkMin[1]; j <= ijkMax[1]; ++j) {
			for (int i = ijkMin[0]; i <= ijkMax[0]; ++i) {
				*it = getVertexLinearId(i, j, k);
				++it;
			}
		}
	}

	return ids;
}

/*!
	Extract a vertex subset.

	\param[in] idxMin is the linear index of the lower bound
	\param[in] idxMax is the linear index of the upper bound
	\result The linear indices of the vertex subset.
*/
std::vector<long> VolCartesian::extractVertexSubSet(int const &idxMin, int const &idxMax)
{
	return extractVertexSubSet(getVertexCartesianId(idxMin), getVertexCartesianId(idxMax));
}

/*!
	Extract a vertex subset.

	\param[in] pointMin is the lower bound
	\param[in] pointMax is the upper bound
	\result The linear indices of the vertex subset.
*/
std::vector<long> VolCartesian::extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax)
{
	return extractVertexSubSet(locatePointCartesian(pointMin), locatePointCartesian(pointMax));
}

/*!
	Gets the origin of the patch.

	The origin is the lower-left-back corner of the box that defines the patch
	domain.

	\return The origin of the patch.
*/
std::array<double, 3> VolCartesian::getOrigin() const
{
	return m_minCoords;
}

/*!
	Sets the origin of the patch.

	The origin is the lower-left-back corner.

	\param origin is the new origin of the patch
*/
void VolCartesian::setOrigin(const std::array<double, 3> &origin)
{
	std::array<double, 3> translation = origin - getOrigin();
	translate(translation);
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
 */
void VolCartesian::translate(std::array<double, 3> translation)
{
	for (int n = 0; n < 3; ++n) {
		m_minCoords[n] += translation[n];
		m_maxCoords[n] += translation[n];
		for (int i = 1; i < m_nVertices1D[n]; ++i) {
			m_vertexCoords[n][i] += translation[n];
			m_cellCenters[n][i]  += translation[n];
		}
	}

	setBoundingBox(m_minCoords, m_maxCoords);

	VolumeKernel::translate(translation);
}

/*!
	Gets the lengths of the sides of the box that defines the patch domain.

	\return The lengths of the sides of the box that defines the patch domain.
*/
std::array<double, 3> VolCartesian::getLengths() const
{
	return (m_maxCoords - m_minCoords);
}

/*!
	Sets the lengths of the sides of the box that defines the patch domain.

	\param lengths is an array with the lengths of the sides of the box that
	defines the patch domain
*/
void VolCartesian::setLengths(const std::array<double, 3> &lengths)
{
	// Set the maximum coordinate
	m_maxCoords = m_minCoords + lengths;

	// Reset the bounding box
	setBoundingBox(m_minCoords, m_maxCoords);

	// If needed update the discretization
	if (m_nVertices > 0) {
		// Rebuild the discretization
		setDiscretization(m_nCells1D);

		// Rebuild vertex coordinates
		for (Vertex vertex : m_vertices) {
			long id = vertex.getId();
			vertex.setCoords(evalVertexCoords(id));
		}
	}
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor vector
 */
void VolCartesian::scale(std::array<double, 3> scaling)
{
	for (int n = 0; n < 3; ++n) {
		m_maxCoords[n] = m_minCoords[n] + scaling[n] * (m_maxCoords[n] - m_minCoords[n]);
		for (int i = 1; i < m_nVertices1D[n]; ++i) {
			m_vertexCoords[n][i] = m_minCoords[n] + scaling[n] * (m_vertexCoords[n][i] - m_minCoords[n]);
			m_cellCenters[n][i]  = m_minCoords[n] + scaling[n] * (m_cellCenters[n][i] - m_minCoords[n]);
		}

		m_cellSpacings[n] *= scaling[n];
	}

	initializeCellVolume();

	initializeInterfaceArea();

	setBoundingBox(m_minCoords, m_maxCoords);

	VolumeKernel::scale(scaling);
}

/*!
	Transform cell data to point data by calculating the mean of incident
	cells in each vertex.

	\param[in] cellData contains the data on cells
	\result The cell data converted to vertex data.
*/
std::vector<double> VolCartesian::convertToVertexData(const std::vector<double> &cellData) const
{
	int dimension = getDimension();

	int nContribs_x = 2;
	int nContribs_y = 2;
	int nContribs_z = dimension - 1;

	std::vector<int> nodeCounter(getVertexCount());
	std::vector<double> vertexData(getVertexCount());
	std::fill (vertexData.begin(), vertexData.end(), 0.);

	for (int k = 0; (isThreeDimensional()) ? (k < m_nCells1D[Vertex::COORD_Z]) : (k < 1); k++) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; ++j) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; ++i) {
				// Get cell contribution
				long cellId = getCellLinearId(i, j, k);
				double cellContrib = cellData[cellId];

				// Eval vertex data
				for (int n = 0; n < nContribs_z; n++) {
					int k_v = k + n;
					for (int m = 0; m < nContribs_y; ++m) {
						int j_v = j + m;
						for (int l = 0; l < nContribs_x; ++l) {
							int i_v = i + l;

							// Get vertex index
							long vertexId = getVertexLinearId(i_v, j_v, k_v);

							// Sum cell contribution
							vertexData[vertexId] += cellContrib;
							nodeCounter[vertexId]++;
						}
					}
				}
			}
		}
	}

	// Average values
	for (int i = 0; i < getVertexCount(); ++i) {
		vertexData[i] /= nodeCounter[i];
	}

	return vertexData;
}

/*!
	Transform node data to cell data by calculating the mean of incident
	vertices of each cell

	\param[in] vertexData contains the data on vertices
	\result The vertex data converted to cell data.
*/
std::vector<double> VolCartesian::convertToCellData(const std::vector<double> &vertexData) const
{
	int dimension = getDimension();

	int nContribs_x = 2;
	int nContribs_y = 2;
	int nContribs_z = dimension - 1;

	std::vector<double> cellData(getCellCount());
	std::fill (cellData.begin(), cellData.end(), 0.);

	for (int k = 0; (isThreeDimensional()) ? (k < m_nCells1D[Vertex::COORD_Z]) : (k < 1); k++) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; ++j) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; ++i) {
				// Eval cell data
				long cellId = getCellLinearId(i, j, k);
				for (int n = 0; n < nContribs_z; ++n) {
					int k_v = k + n;
					for (int m = 0; m < nContribs_y; ++m) {
						int j_v = j + m;
						for (int l = 0; l < nContribs_x; ++l) {
							int i_v = i + l;

							// Get vertex contribution
							long vertexId = getVertexLinearId(i_v, j_v, k_v);
							double vertexContrib = vertexData[vertexId];

							// Sum vertex contribution
							cellData[cellId] += vertexContrib;
						}
					}
				}
			}
		}
	}

	double weight = pow(0.5, dimension);
	for (int i = 0; i < getCellCount(); ++i) {
		cellData[i] *= weight;
	}

	return cellData;
}


/*!
	Calculates bi-/tri- linear interpolation stencil on cells for a
	given point.

	At boundaries stencil is reduced to assure positive weights. If the
	point is outside a null-stencil is returned

	\param[in] point are the point coordinates
	\param[out] stencil are the linear indices of the interpolation stencil
	\param[out] weights are the weights associated to stencil
	\return The number of cells used in the interpolation stencil. If the
	point is outside a null stencil size is returned.
*/
int VolCartesian::linearCellInterpolation(std::array<double,3> &point,
                                            std::vector<int> &stencil,
                                            std::vector<double> &weights)
{
	std::array<int, 3> ijk_point = locatePointCartesian(point);
	if (ijk_point[0] < 0) {
		stencil.clear();
		weights.clear();
		return 0;
	}

	int dimension = getDimension();

	int nContribs_x = 1;
	int nContribs_y = 1;
	int nContribs_z = 1;

	std::array< std::array<int,2>, 3> cStencil;
	std::array< std::array<double,2>, 3> cWeights;
	for (int d = 0; d < dimension; ++d) {
		// Find cell index
		int index_point = ijk_point[d];
		if (point[d] < m_cellCenters[d][index_point]) {
			index_point = index_point - 1;
		}

		int index_next = index_point + 1;

		if (index_point < 0) {
			cStencil[d][0] = 0.;
			cWeights[d][0] = 1.;
		} else if (index_next > m_nCells1D[d] - 1) {
			cStencil[d][0] = m_nCells1D[d] - 1;
			cWeights[d][0] = 1.;
		} else {
			if (d == 0) {
				nContribs_x = 2;
			} else if (d == 1) {
				nContribs_y = 2;
			} else {
				nContribs_z = 2;
			}

			cStencil[d][0] = index_point;
			cStencil[d][1] = index_next;

			cWeights[d][1] = (point[d] - m_cellCenters[d][index_point]) / m_cellSpacings[d]  ;
			cWeights[d][0] = 1.0 - cWeights[d][1];
		}
	}

	for (int d = dimension; d < 3; ++d) {
		cStencil[d][0] = 0.;
		cWeights[d][0] = 1.;
	}

	int stencilSize = nContribs_x * nContribs_y * nContribs_z;
	stencil.resize(stencilSize);
	weights.resize(stencilSize);

	std::vector<int>::iterator itrStencil = stencil.begin();
	std::vector<double>::iterator itrWeights = weights.begin();
	for (int k = 0; k < nContribs_z; ++k) {
		for (int j = 0; j < nContribs_y; ++j) {
			for (int i = 0; i < nContribs_x; ++i) {
				int &is = cStencil[0][i];
				int &js = cStencil[1][j];
				int &ks = cStencil[2][k];

				double &iw = cWeights[0][i];
				double &jw = cWeights[1][j];
				double &kw = cWeights[2][k];

				*itrStencil = getCellLinearId(is, js, ks);
				*itrWeights = iw *jw * kw;

				++itrStencil;
				++itrWeights;
			}
		}
	}

	return stencilSize;
}

/*!
	Calculates bi-/tri- linear interpolation stencil on vertex for a
	given point.

	At boundaries stencil is reduced to assure positive weights. If the
	point is outside a null-stencil is returned

	\param[in] point are the point coordinates
	\param[out] stencil are the linear indices of the interpolation stencil
	\param[out] weights are the weights associated to stencil
	\return The number of cells used in the interpolation stencil. If the
	point is outside a null stencil size is returned.
*/
int VolCartesian::linearVertexInterpolation(std::array<double,3> &point,
                                              std::vector<int> &stencil,
                                              std::vector<double> &weights)
{
	std::array<int, 3> ijk_point = locatePointCartesian(point);
	if (ijk_point[0] < 0) {
		stencil.clear();
		weights.clear();
		return 0;
	}

	int dimension = getDimension();

	int nContribs_x = 2;
	int nContribs_y = 2;
	int nContribs_z = dimension - 1;

	std::array< std::array<int,2>, 3> cStencil;
	std::array< std::array<double,2>, 3> cWeights;
	for (int d = 0; d < dimension; ++d) {
		int index_point = ijk_point[d];
		int index_next  = index_point +1;

		cStencil[d][0] = index_point;
		cStencil[d][1] = index_next;

		cWeights[d][1] = (point[d] - m_vertexCoords[d][index_point]) / m_cellSpacings[d];
		cWeights[d][0] = 1.0 - cWeights[d][1];
	}

	for (int d = dimension; d < 3; ++d) {
		cStencil[d][0] = 0;
		cWeights[d][0] = 1.;
	}

	int stencilSize = pow(2, dimension);
	stencil.resize(stencilSize);
	weights.resize(stencilSize);

	std::vector<int>::iterator itrStencil    = stencil.begin();
	std::vector<double>::iterator itrWeights = weights.begin();
	for (int k = 0; k < nContribs_z; ++k) {
		for (int j = 0; j < nContribs_y; ++j) {
			for (int i = 0; i < nContribs_x; ++i) {
				int &is = cStencil[0][i];
				int &js = cStencil[1][j];
				int &ks = cStencil[2][k];

				double &iw = cWeights[0][i];
				double &jw = cWeights[1][j];
				double &kw = cWeights[2][k];

				*itrStencil = getVertexLinearId(is, js, ks);
				*itrWeights = iw * jw * kw;

				++itrStencil;
				++itrWeights;
			}
		}
	}

	return stencilSize;
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> VolCartesian::evalCellCentroid(const long &id) const
{
	std::array<int, 3> ijk = getCellCartesianId(id);

	std::array<double, 3> centroid = {{0, 0, 0}};
	for (int n = 0; n < getDimension(); ++n) {
		centroid[n] = m_cellCenters[n][ijk[n]];
	}

	return centroid;
}

}
