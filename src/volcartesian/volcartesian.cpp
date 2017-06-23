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
	__reset();
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param lengths are the lengths of the domain
	\param nCells are the numbers of cells of the patch
*/
VolCartesian::VolCartesian(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               const std::array<double, 3> &lengths,
                               const std::array<int, 3> &nCells)
	: VolumeKernel(id, dimension, false)
{
	initialize();
	__reset();

	setOrigin(origin);
	setLengths(lengths);

	setDiscretization(nCells);
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param nCells1D is the number of cells along each direction
*/
VolCartesian::VolCartesian(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, int nCells)
	: VolumeKernel(id, dimension, false)
{
	initialize();
	__reset();

	setOrigin(origin);
	setLengths({{length, length, length}});

	setDiscretization({{nCells, nCells, nCells}});
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
	__reset();

	setOrigin(origin);
	setLengths({{length, length, length}});

	int nCells = (int) std::ceil(length / dh);
	setDiscretization({{nCells, nCells, nCells}});
}

/*!
	Creates a new patch restoring the patch saved in the specified stream.

	\param stream is the stream to read from
*/
VolCartesian::VolCartesian(std::istream stream)
	: VolumeKernel(false)
{
	initialize();
	restore(stream);
}

/*!
	Internal function to reset the patch.
*/
void VolCartesian::reset()
{
	// Reset the patch kernel
	VolumeKernel::reset();

	// Reset the current patch
	__reset();
}

/*!
	Reset the patch.
*/
void VolCartesian::__reset()
{
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
	m_nInterfaces = 0;

	// Set the light memory mode
	//
	// This will reset the patch data structures.
	setMemoryMode(MemoryMode::MEMORY_LIGHT);
}

/*!
	Initialize the data structures of the patch.
*/
void VolCartesian::initialize()
{
	// Set the light memory mode
	m_memoryMode = MemoryMode::MEMORY_LIGHT;

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
	log::cout() << "  - Total interface count: " << m_nInterfaces << "\n";

	// Cell volume
	initializeCellVolume();

	// Interface area
	initializeInterfaceArea();

	// Create cells, vertices and interfaces
	if (getMemoryMode() == MemoryMode::MEMORY_NORMAL) {
		update();
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
ElementInfo::Type VolCartesian::getCellType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getCellType();
}

/*!
	Gets the element type for the cells in the patch.

	\return The element type for the cells in the patch.
*/
ElementInfo::Type VolCartesian::getCellType() const
{
	if (isThreeDimensional()) {
		return ElementInfo::VOXEL;
	} else {
		return ElementInfo::PIXEL;
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
ElementInfo::Type VolCartesian::getInterfaceType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getInterfaceType();
}

/*!
	Gets the element type for the interfaces in the patch.

	\return The element type for the interfaces in the patch.
*/
ElementInfo::Type VolCartesian::getInterfaceType() const
{
	if (isThreeDimensional()) {
		return ElementInfo::PIXEL;
	} else {
		return ElementInfo::LINE;
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
	Set the memory mode.

	\param mode is the memory mode that will be set
*/
void VolCartesian::setMemoryMode(MemoryMode mode)
{
	setMemoryMode(mode, true);
}

/*!
	Set the memory mode.

	\param mode is the memory mode that will be set
	\param updatePatch if set to true, updates the data structures to switch
	into the specified mode, otherwise only the memory-mode flag is updated
	and the data structures needs to be updated externally
*/
void VolCartesian::setMemoryMode(MemoryMode mode, bool updatePatch)
{
	if (mode == m_memoryMode) {
		return;
	}

	// Set the mode
	m_memoryMode = mode;

	// Update the data structures
	if (updatePatch) {
		switch (mode) {

		case MemoryMode::MEMORY_NORMAL:
			update();
			break;

		case MemoryMode::MEMORY_LIGHT:
			VolumeKernel::reset();
			break;

		}
	}
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
	Updates the patch.

	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the update.
*/
std::vector<adaption::Info> VolCartesian::_updateAdaption(bool trackAdaption)
{
	log::cout() << ">> Updating cartesian mesh\n";

	// Enable advanced editing
	setExpert(true);

	// Reset the mesh
	//
	// We need to reset only the generic data of the patch, therefore we call
	// the 'reset' implementation of the kernel.
	VolumeKernel::reset();

	// Definition of the mesh
	addVertices();
	addCells();
	addInterfaces();

	// Disable advanced editing
	setExpert(false);

	// Adaption info
	std::vector<adaption::Info> adaptionData;
	if (trackAdaption) {
		adaptionData.emplace_back();
		adaption::Info &adaptionCellInfo = adaptionData.back();
		adaptionCellInfo.type   = adaption::TYPE_CREATION;
		adaptionCellInfo.entity = adaption::ENTITY_CELL;
		adaptionCellInfo.current.reserve(m_cells.size());
		for (auto &cell : m_cells) {
			adaptionCellInfo.current.emplace_back();
			long &cellId = adaptionCellInfo.current.back();
			cellId = cell.getId();
		}

		adaptionData.emplace_back();
		adaption::Info &adaptionInterfaceInfo = adaptionData.back();
		adaptionInterfaceInfo.type   = adaption::TYPE_CREATION;
		adaptionInterfaceInfo.entity = adaption::ENTITY_INTERFACE;
		adaptionInterfaceInfo.current.reserve(m_interfaces.size());
		for (auto &interface : m_interfaces) {
			adaptionInterfaceInfo.current.emplace_back();
			long &interfaceId = adaptionInterfaceInfo.current.back();
			interfaceId = interface.getId();
		}
	} else {
		adaptionData.emplace_back();
	}

	// Updating the adaption brings the patch is in normal memory mode
	setMemoryMode(MemoryMode::MEMORY_NORMAL, false);

	// Done
	return adaptionData;
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
	ElementInfo::Type cellType = getCellType();

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
				cell.setVertex(0, getVertexLinearId(i,     j,     k));
				cell.setVertex(1, getVertexLinearId(i + 1, j,     k));
				cell.setVertex(2, getVertexLinearId(i,     j + 1, k));
				cell.setVertex(3, getVertexLinearId(i + 1, j + 1, k));
				if (isThreeDimensional()) {
					cell.setVertex(4, getVertexLinearId(i,     j,     k + 1));
					cell.setVertex(5, getVertexLinearId(i + 1, j,     k + 1));
					cell.setVertex(6, getVertexLinearId(i,     j + 1, k + 1));
					cell.setVertex(7, getVertexLinearId(i + 1, j + 1, k + 1));
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

	m_interfaces.reserve(m_nInterfaces);
	for (int n = 0; n < getDimension(); n++) {
		addInterfacesDirection(n);
	}
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
	ElementInfo::Type interfaceType = getInterfaceType();

	const ElementInfo &interfaceTypeInfo = ElementInfo::getElementInfo(interfaceType);
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
				std::array<int, 3> ownerIJK;
				for (int n = 0; n < 3; n++) {
					ownerIJK[n] = counters[n];
				}
				if (counters[direction] == (interfaceCount1D[direction] - 1)) {
					ownerIJK[direction] -= 1;
				}

				long ownerId = getCellLinearId(ownerIJK);

				// Neighbour id
				long neighId;
				if (counters[direction] != 0 && counters[direction] != interfaceCount1D[direction] - 1) {
					std::array<int, 3> neighIJK;
					for (int n = 0; n < 3; n++) {
						neighIJK[n] = counters[n];
					}
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
	const int DUMP_VERSION = 1;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolCartesian::_dump(std::ostream &stream)
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
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void VolCartesian::_restore(std::istream &stream)
{
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
	setMemoryMode(memoryMode);
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
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> VolCartesian::_findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList) const
{
	int neighSide      = face % 2;
	int neighDirection = std::floor(face / 2);

	std::array<int, 3> neighIjk(getCellCartesianId(id));
	if (neighSide == 0) {
		neighIjk[neighDirection]--;
	} else {
		neighIjk[neighDirection]++;
	}

	std::vector<long> neighs;
	if (isCellCartesianIdValid(neighIjk)) {
		long neighId = getCellLinearId(neighIjk);
		if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
			neighs.push_back(neighId);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given edge.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param edge is an edge of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given edge.
*/
std::vector<long> VolCartesian::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return neighs;
	}

	// Diagonal neighbour
	std::array<int, 3> diagNeighIjk(getCellCartesianId(id) + m_edgeNeighDeltas[edge]);
	if (isCellCartesianIdValid(diagNeighIjk)) {
		long diagNeighId = getCellLinearId(diagNeighIjk);
		if (std::find(blackList.begin(), blackList.end(), diagNeighId) == blackList.end()) {
			utils::addToOrderedVector<long>(diagNeighId, neighs);
		}
	}

	// Faces incident to the edge
	for (const auto &face : m_edgeFaces[edge]) {
		std::vector<long> faceNeighIds = _findCellFaceNeighs(id, face, blackList);
		for (const auto &faceNeighId : faceNeighIds) {
			utils::addToOrderedVector<long>(faceNeighId, neighs);
		}
	}

	// Done
	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given local vertex.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> VolCartesian::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	std::array<int, 3> cellIjk   = getCellCartesianId(id);
	std::array<int, 3> vertexIjk = getVertexCartesianId(cellIjk, vertex);

	std::vector<long> neighs;
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
		if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
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

	std::vector<int> nodeCounter(getVertexCount());
	std::vector<double> vertexData(getVertexCount());
	std::fill (vertexData.begin(), vertexData.end(), 0.);

	for (int k = 0; k < m_nCells1D[Vertex::COORD_Z]; ++k) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; ++j) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; ++i) {
				// Cell index
				long cellId = getCellLinearId(i, j, k);

				for (int n = 0; n < dimension - 1; ++n) {
					for (int m = 0; m < 2; ++m) {
						for (int l = 0; l < 2; ++l) {
							// Vertex index
							int i_v = i + l;
							int j_v = j + m;
							int k_v = j + n;
							long vertexId = getVertexLinearId(i_v, j_v, k_v);

							vertexData[vertexId] += cellData[cellId];
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

	std::vector<double> cellData(getCellCount());
	std::fill (cellData.begin(), cellData.end(), 0.);

	for (int k = 0; k < m_nCells1D[Vertex::COORD_Z]; ++k) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; ++j) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; ++i) {
    				// Cell index
				long cellId = getCellLinearId(i, j, k);

				for (int n = 0; n < dimension - 1; ++n) {
					for (int m = 0; m < 2; ++m) {
						for (int l = 0; l < 2; ++l) {
							// Vertex index
							int i_v = i + l;
							int j_v = j + m;
							int k_v = j + n;
							long vertexId = getVertexLinearId(i_v, j_v, k_v);

							cellData[cellId] += vertexData[vertexId];
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

	std::array<int, 3> ijk_next;
	std::array<int, 3> nS = {{1,1,1}};
	std::array< std::array<int,2>, 3> cStencil;
	std::array< std::array<double,2>, 3> cWeights;
	for (int d = 0; d < dimension; ++d) {
		// Find cell index
		if (point[d] < m_cellCenters[d][ijk_point[d]]) {
			ijk_point[d] = ijk_point[d] - 1;
		}

		ijk_next[d] = ijk_point[d] + 1;

		if (ijk_point[d] < 0) {
			cStencil[d][0] = 0.;
			cWeights[d][0] = 1.;
		} else if (ijk_next[d] > m_nCells1D[d] - 1) {
			cStencil[d][0] = m_nCells1D[d] - 1;
			cWeights[d][0] = 1.;
		} else {
			nS[d] = 2;

			cStencil[d][0] = ijk_point[d];
			cStencil[d][1] = ijk_next[d];

			cWeights[d][1] = (point[d] - m_cellCenters[d][ijk_point[d]]) / m_cellSpacings[d]  ;
			cWeights[d][0] = 1.0 - cWeights[d][1];
		}
	}

	for (int d = dimension; d < 3; ++d) {
		cStencil[d][0] = 0.;
		cWeights[d][0] = 1.;
	}

	int stencilSize = nS[0] * nS[1] * nS[2];
	stencil.resize(stencilSize);
	weights.resize(stencilSize);

	std::vector<int>::iterator itrStencil = stencil.begin();
	std::vector<double>::iterator itrWeights = weights.begin();
	for (int k = 0; k < nS[2]; ++k) {
		for (int j = 0; j < nS[1]; ++j) {
			for (int i = 0; i < nS[0]; ++i) {
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

	std::array<int, 3> ijk_next;
	std::array< std::array<int,2>, 3> cStencil;
	std::array< std::array<double,2>, 3> cWeights;
	for (int d = 0; d < dimension; ++d) {
		ijk_next[d] = ijk_point[d] +1;

		cStencil[d][0] = ijk_point[d];
		cStencil[d][1] = ijk_next[d];

		cWeights[d][1] = (point[d] - m_vertexCoords[d][ijk_point[d]]) / m_cellSpacings[d];
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
	for (int k = 0; k < dimension - 1; ++k) {
		for (int j = 0; j < 2; ++j) {
			for (int i = 0; i < 2; ++i) {
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
