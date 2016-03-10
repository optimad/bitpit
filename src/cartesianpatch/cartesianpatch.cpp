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

#include <cmath>
#include <bitset>

#include "bitpit_common.hpp"

#include "cartesianpatch.hpp"

namespace bitpit {

/*!
	\ingroup cartesianpatch
	@{
*/

/*!
	\class CartesianPatch

	\brief The CartesianPatch defines a Cartesian patch.

	CartesianPatch defines a Cartesian patch.
*/

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param lengths are the lengths of the domain
	\param nCells are the numbers of cells of the patch
*/
CartesianPatch::CartesianPatch(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               const std::array<double, 3> &lengths,
                               const std::array<int, 3> &nCells)
	: VolumeKernel(id, dimension, false)
{
	initialize(origin, lengths, nCells);
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param nCells1D is the number of cells along each direction
*/
CartesianPatch::CartesianPatch(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, int nCells1D)
	: VolumeKernel(id, dimension, false)
{
	// Number of cells
	std::array<int, 3> nCells;
	for (int n = 0; n < dimension; n++) {
		nCells[n] = nCells1D;
	}

	if (!isThreeDimensional()) {
		nCells[Vertex::COORD_Z] = 0;
	}

	// Domain lengths
	std::array<double, 3> lengths;
	for (int n = 0; n < dimension; n++) {
		lengths[n] = length;
	}

	if (!isThreeDimensional()) {
		lengths[Vertex::COORD_Z] = 0;
	}

	// Patch initialization
	initialize(origin, lengths, nCells);
}

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param origin is the origin of the domain
	\param length is the length of the domain
	\param dh is the maximum allowed mesh spacing
*/
CartesianPatch::CartesianPatch(const int &id, const int &dimension,
                               const std::array<double, 3> &origin,
                               double length, double dh)
	: VolumeKernel(id, dimension, false)
{
	// Number of cells
	std::array<int, 3> nCells;
	for (int n = 0; n < dimension; n++) {
		nCells[n] = (int) std::ceil(length / dh);
	}

	if (!isThreeDimensional()) {
		nCells[Vertex::COORD_Z] = 0;
	}

	// Domain lengths
	std::array<double, 3> lengths;
	for (int n = 0; n < dimension; n++) {
		lengths[n] = length;
	}

	if (!isThreeDimensional()) {
		lengths[Vertex::COORD_Z] = 0;
	}

	// Patch initialization
	initialize(origin, lengths, nCells);
}

/*!
	Initializes the patch

	\param origin is the origin of the domain
	\param lengths are the lengths of the domain
	\param nCells are the numbers of cells of the patch
*/
void CartesianPatch::initialize(const std::array<double, 3> &origin,
                                const std::array<double, 3> &lengths,
                                const std::array<int, 3> &nCells)
{
	log::cout() << ">> Initializing cartesian patch\n";

	// Info sulle celle
	for (int n = 0; n < getDimension(); ++n) {
		// Initialize cells
		m_nCells1D[n]     = nCells[n];
		m_minCoords[n]    = origin[n];
		m_maxCoords[n]    = origin[n] + lengths[n];
		m_cellSpacings[n] = lengths[n] / m_nCells1D[n];

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
		m_minCoords[Vertex::COORD_Z]    = 0.;
		m_maxCoords[Vertex::COORD_Z]    = 0.;
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

	// Interface area
	initializeInterfaceArea();

	// Cell volume
	initializeCellVolume();

	// Normals
	int i = 0;
	for (int n = 0; n < getDimension(); n++) {
		for (int k = -1; k <= 1; k += 2) {
			std::array<double, 3> normal = {{0.0, 0.0, 0.0}};
			normal[n] = k;

			m_normals[i++] = normal;
		}
	}

	// Deltas for the evaluation of the vertex neighbours
	m_vertexNeighDeltas = std::vector<std::array<int, 3>>(std::pow(2, getDimension()));
	m_vertexNeighDeltas[0] = {{ 0,  0, 0}};
	m_vertexNeighDeltas[1] = {{-1,  0, 0}};
	m_vertexNeighDeltas[2] = {{ 0, -1, 0}};
	m_vertexNeighDeltas[3] = {{-1, -1, 0}};
	if (isThreeDimensional()) {
		m_vertexNeighDeltas[4] = {{ 0,  0, -1}};
		m_vertexNeighDeltas[5] = {{-1,  0, -1}};
		m_vertexNeighDeltas[6] = {{ 0, -1, -1}};
		m_vertexNeighDeltas[7] = {{-1, -1, -1}};
	}

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
}

/*!
	Initializes cell volume
*/
void CartesianPatch::initializeCellVolume()
{
	m_cellVolume = m_cellSpacings[Vertex::COORD_X] * m_cellSpacings[Vertex::COORD_Y];
	if (isThreeDimensional()) {
		m_cellVolume *= m_cellSpacings[Vertex::COORD_Z];
	}
}

/*!
	Initializes interface area
*/
void CartesianPatch::initializeInterfaceArea()
{
	for (int n = 0; n < getDimension(); ++n) {
		m_interfaceArea[n] = m_cellVolume / m_cellSpacings[n];
	}

	if (!isThreeDimensional()) {
		m_interfaceArea[Vertex::COORD_Z] = 0.;
	}
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long CartesianPatch::getVertexCount() const
{
	return m_nVertices;
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long CartesianPatch::getCellCount() const
{
	return m_nCells;
}

/*!
	Gets the element type for the cell with the specified id.

	\param id is the id of the requested cell
	\return The element type for the cell with the specified id.
*/
ElementInfo::Type CartesianPatch::getCellType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getCellType();
}

/*!
	Gets the element type for the cells in the patch.

	\return The element type for the cells in the patch.
*/
ElementInfo::Type CartesianPatch::getCellType() const
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
long CartesianPatch::getInterfaceCount() const
{
	return m_nInterfaces;
}

/*!
	Gets the element type for the interface with the specified id.

	\param id is the id of the requested interface
	\return The element type for the interface with the specified id.
*/
ElementInfo::Type CartesianPatch::getInterfaceType(const long &id) const
{
	BITPIT_UNUSED(id);

	return getInterfaceType();
}

/*!
	Gets the element type for the interfaces in the patch.

	\return The element type for the interfaces in the patch.
*/
ElementInfo::Type CartesianPatch::getInterfaceType() const
{
	if (isThreeDimensional()) {
		return ElementInfo::PIXEL;
	} else {
		return ElementInfo::LINE;
	}
}

/*!
	Destroys the patch.
*/
CartesianPatch::~CartesianPatch()
{
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double CartesianPatch::evalCellVolume(const long &id)
{
	BITPIT_UNUSED(id);

	return m_cellVolume;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double CartesianPatch::evalCellSize(const long &id)
{
	BITPIT_UNUSED(id);

	return pow(m_cellVolume, 1. / getDimension());
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double CartesianPatch::evalInterfaceArea(const long &id)
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
std::array<double, 3> CartesianPatch::evalInterfaceNormal(const long &id)
{
	const Interface &interface = getInterface(id);
	int ownerFace = interface.getOwnerFace();

	return m_normals[ownerFace];
}

/*!
	Get cell spacings of the patch.

	\result Cell spacings of the patch.
*/
std::array<double, 3> CartesianPatch::getSpacing() const
{
	return m_cellSpacings;
}

/*!
	Get cell spacing along the specificed direction.

	\param[in] direction is the direction along which the spacing is
	requested
	\result The cell spacing along the specificed direction.
*/
double CartesianPatch::getSpacing(const int &direction) const
{
	return m_cellSpacings[direction];
}

/*!
	Updates the patch.

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> CartesianPatch::_updateAdaption(bool trackAdaption)
{
	log::cout() << ">> Updating cartesian mesh\n";

	// Enable advanced editing
	setExpert(true);

	// Reset the mesh
	reset();

	// Definition of the mesh
	addVertices();
	addCells();
	addInterfaces();

	// Disable advanced editing
	setExpert(false);

	// Adaption info
	std::vector<Adaption::Info> adaptionData;
	if (trackAdaption) {
		adaptionData.emplace_back();
		Adaption::Info &adaptionCellInfo = adaptionData.back();
		adaptionCellInfo.type   = Adaption::TYPE_CREATION;
		adaptionCellInfo.entity = Adaption::ENTITY_CELL;
		adaptionCellInfo.current.reserve(m_cells.size());
		for (auto &cell : m_cells) {
			adaptionCellInfo.current.emplace_back();
			unsigned long &cellId = adaptionCellInfo.current.back();
			cellId = cell.getId();
		}

		adaptionData.emplace_back();
		Adaption::Info &adaptionInterfaceInfo = adaptionData.back();
		adaptionInterfaceInfo.type   = Adaption::TYPE_CREATION;
		adaptionInterfaceInfo.entity = Adaption::ENTITY_INTERFACE;
		adaptionInterfaceInfo.current.reserve(m_interfaces.size());
		for (auto &interface : m_interfaces) {
			adaptionInterfaceInfo.current.emplace_back();
			unsigned long &interfaceId = adaptionInterfaceInfo.current.back();
			interfaceId = interface.getId();
		}
	} else {
		adaptionData.emplace_back();
	}

	// Done
	return adaptionData;
}

/*!
	Creates the vertices of the patch.
*/
void CartesianPatch::addVertices()
{
	log::cout() << "  >> Creating vertices\n";

	log::cout() << "    - Vertex count: " << m_nVertices << "\n";

	m_vertices.reserve(m_nVertices);
	for (int k = 0; (isThreeDimensional()) ? (k < m_nVertices1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (int j = 0; j < m_nVertices1D[Vertex::COORD_Y]; j++) {
			for (int i = 0; i < m_nVertices1D[Vertex::COORD_X]; i++) {
				// Vertex coordinates
				std::array<double, 3> coords;
				coords[Vertex::COORD_X] = m_vertexCoords[Vertex::COORD_X][i];
				coords[Vertex::COORD_Y] = m_vertexCoords[Vertex::COORD_Y][j];
				if (isThreeDimensional()) {
					coords[Vertex::COORD_Z] = m_vertexCoords[Vertex::COORD_Z][k];
				} else {
					coords[Vertex::COORD_Z] = 0.0;
				}

				// Add vertex
				long id_vertex = getVertexLinearId(i, j, k);
				VolumeKernel::addVertex(std::move(coords), id_vertex);
			}
		}
	}
}

/*!
	Creates the cells of the patch.
*/
void CartesianPatch::addCells()
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
void CartesianPatch::addInterfaces()
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
std::array<int, 3> CartesianPatch::getInterfaceCountDirection(const int &direction)
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
void CartesianPatch::addInterfacesDirection(const int &direction)
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
				if (counters[direction] == 0) {
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
	Marks a cell for refinement.

	This is a void function since mesh refinement is not implemented
	for Cartesian meshes.

	\param id is the id of the cell that needs to be refined
*/
bool CartesianPatch::_markCellForRefinement(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Marks a cell for coarsening.

	This is a void function since mesh coarsening is not implemented
	for Cartesian meshes.

	\param id the cell to be refined
*/
bool CartesianPatch::_markCellForCoarsening(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Enables cell balancing.

	This is a void function since mesh coarsening is not implemented
	for Cartesian meshes.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
bool CartesianPatch::_enableCellBalancing(const long &id, bool enabled)
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(enabled);

	return false;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] point is the point to be checked
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool CartesianPatch::isPointInside(const std::array<double, 3> &point)
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
long CartesianPatch::locatePoint(const std::array<double, 3> &point)
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
std::array<int, 3> CartesianPatch::locatePointCartesian(const std::array<double, 3> &point)
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
long CartesianPatch::locateClosestVertex(std::array<double,3> const &point) const
{
	return getVertexLinearId(locateClosestVertexCartesian(point));
}

/*!
	Locates the closest vertex of the given point.

	\param[in] point is the point
	\result The set of Cartesian id of the closest vertex of the given
	point.
*/
std::array<int, 3> CartesianPatch::locateClosestVertexCartesian(std::array<double,3> const &point) const
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
	Converts the cell cartesian notation to a linear notation
*/
long CartesianPatch::getCellLinearId(const int &i, const int &j, const int &k) const
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
long CartesianPatch::getCellLinearId(const std::array<int, 3> &ijk) const
{
	return getCellLinearId(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}

/*!
	Converts a cell linear index to a set of cartesian indices.

	No check on bounds is performed.

	\param[in] idx is the linear index of the cell
	\result Returns the set of cartesian indices of the cell.
*/
std::array<int, 3> CartesianPatch::getCellCartesianId(long const &idx) const
{
	int ijPlane = m_nCells1D[0] * m_nCells1D[1];

	std::array<int, 3> id;
	id[0] = idx % m_nCells1D[0];
	if (isThreeDimensional()) {
		id[2] = idx / ijPlane;
		id[1] = (idx - id[2] * ijPlane) / m_nCells1D[0];
	} else {
		id[2] = -1;
		id[1] = idx / m_nCells1D[0];
	}

	return id;
}

/*!
	Checks if a cell cartesian index is valid.

	\param ijk is the set of cartesian indices of the cell
	\result Returns true if the index is valid, otherwise it returns false.
*/
bool CartesianPatch::isCellCartesianIdValid(const std::array<int, 3> &ijk) const
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
long CartesianPatch::getVertexLinearId(const int &i, const int &j, const int &k) const
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
long CartesianPatch::getVertexLinearId(const std::array<int, 3> &ijk) const
{
	return getVertexLinearId(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}

/*!
	Converts a vertex linear index to a set of cartesian indices.

	No check on bounds is performed.

	\param[in] idx is the linear index of the vertex
	\result Returns the set of cartesian indices of the vertex.
*/
std::array<int, 3> CartesianPatch::getVertexCartesianId(long const &idx) const
{
	int ijPlane = m_nVertices1D[0] * m_nVertices1D[1];

	std::array<int, 3> id;
	id[0] = idx % m_nVertices1D[0];
	if (isThreeDimensional()) {
		id[2] = idx / ijPlane;
		id[1] = (idx - id[2] * ijPlane) / m_nVertices1D[0];
	} else {
		id[2] = -1;
		id[1] = idx / m_nVertices1D[0];
	}

	return id;
}

/*!
	Gets the cartesian indices of the specified local vertex.

	No check on bounds is performed.

	\param[in] cellIdx is the linear cell index
	\param[in] vertex is the local vertex
	\result Returns the set of cartesian indices of the vertex.
*/
std::array<int, 3> CartesianPatch::getVertexCartesianId(long const &cellIdx, int const &vertex) const
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
std::array<int, 3> CartesianPatch::getVertexCartesianId(const std::array<int, 3> &cellIjk, int const &vertex) const
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
bool CartesianPatch::isVertexCartesianIdValid(const std::array<int, 3> &ijk) const
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
std::vector<long> CartesianPatch::_findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList) const
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
std::vector<long> CartesianPatch::_findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
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
std::vector<long> CartesianPatch::_findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
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
	Evalautes patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void CartesianPatch::evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	minPoint = m_minCoords;
	maxPoint = m_maxCoords;
}

/*!
	Extract a cell subset.

	\param[in] ijkMin is the set of cartesian indices of the lower bound
	\param[in] ijkMax is the set of cartesian indices of the upper bound
	\result The linear indices of the cell subset.
*/
std::vector<long> CartesianPatch::extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax)
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
std::vector<long> CartesianPatch::extractCellSubSet(int const &idxMin, int const &idxMax)
{
	return extractCellSubSet(getCellCartesianId(idxMin), getCellCartesianId(idxMax));
}

/*!
	Extract a cell subset.

	\param[in] pointMin is the lower bound
	\param[in] pointMax is the upper bound
	\result The linear indices of the cell subset.
*/
std::vector<long> CartesianPatch::extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax)
{
	return extractCellSubSet(locatePointCartesian(pointMin), locatePointCartesian(pointMax));
}

/*!
	Extract a vertex subset.

	\param[in] ijkMin is the set of cartesian indices of the lower bound
	\param[in] ijkMax is the set of cartesian indices of the upper bound
	\result The linear indices of the vertex subset.
*/
std::vector<long> CartesianPatch::extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax)
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
std::vector<long> CartesianPatch::extractVertexSubSet(int const &idxMin, int const &idxMax)
{
	return extractVertexSubSet(getVertexCartesianId(idxMin), getVertexCartesianId(idxMax));
}

/*!
	Extract a vertex subset.

	\param[in] pointMin is the lower bound
	\param[in] pointMax is the upper bound
	\result The linear indices of the vertex subset.
*/
std::vector<long> CartesianPatch::extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax)
{
	return extractVertexSubSet(locatePointCartesian(pointMin), locatePointCartesian(pointMax));
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
 */
void CartesianPatch::translate(std::array<double, 3> translation)
{
	for (int n = 1; n < 3; ++n) {
		m_minCoords[n] += translation[n];
		m_maxCoords[n] += translation[n];
		for (int i = 1; i < m_nVertices1D[n]; ++i) {
			m_vertexCoords[n][i] += translation[n];
			m_cellCenters[n][i]  += translation[n];
		}
	}

	VolumeKernel::translate(translation);
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor vector
 */
void CartesianPatch::scale(std::array<double, 3> scaling)
{
	for (int n = 1; n < 3; ++n) {
		m_maxCoords[n] = m_minCoords[n] + scaling[n] * (m_maxCoords[n] - m_minCoords[n]);
		for (int i = 1; i < m_nVertices1D[n]; ++i) {
			m_vertexCoords[n][i] = m_minCoords[n] + scaling[n] * (m_vertexCoords[n][i] - m_minCoords[n]);
			m_cellCenters[n][i]  = m_minCoords[n] + scaling[n] * (m_cellCenters[n][i] - m_minCoords[n]);
		}

		m_cellSpacings[n] *= scaling[n];
	}

	initializeInterfaceArea();

	initializeCellVolume();

	VolumeKernel::scale(scaling);
}

/*!
	Transform cell data to point data by calculating the mean of incident
	cells in each vertex.

	\param[in] cellData contains the data on cells
	\result The cell data converted to vertex data.
*/
std::vector<double> CartesianPatch::convertToVertexData(const std::vector<double> &cellData) const
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
std::vector<double> CartesianPatch::convertToCellData(const std::vector<double> &vertexData) const
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
int CartesianPatch::linearCellInterpolation(std::array<double,3> &point,
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
int CartesianPatch::linearVertexInterpolation(std::array<double,3> &point,
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

	for (int d = 0; d < dimension; ++d) {
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
std::array<double, 3> CartesianPatch::evalCellCentroid(const long &id)
{
	std::array<int, 3> ijk = getCellCartesianId(id);

	std::array<double, 3> centroid;
	for (int n = 0; n < 3; ++n) {
		centroid[n] = m_cellCenters[n][ijk[n]];
	}

	return centroid;
}

/*!
	@}
*/

}
