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

#include <math.h>

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
	\param length is the length of the domain
	\param dh is the maximum allowed mesh spacing
*/
CartesianPatch::CartesianPatch(const int &id, const int &dimension,
                               std::array<double, 3> origin, double length, double dh)
	: Patch(id, dimension)
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
void CartesianPatch::initialize(std::array<double, 3> origin, std::array<double, 3> lengths,
                                std::array<int, 3> nCells)
{
	std::cout << ">> Initializing cartesian patch\n";

	// Info sulle celle
	for (int n = 0; n < getDimension(); ++n) {
		// Initialize cells
		m_nCells1D[n]     = nCells[n];
		m_minCoords[n]    = origin[n] - 0.5 * lengths[n];
		m_maxCoords[n]    = origin[n] + 0.5 * lengths[n];
		m_cellSpacings[n] = lengths[n] / m_nCells1D[n];

		m_cellCenters[n].resize(m_nCells1D[n]);
		for (int i = 0; i < m_nCells1D[n]; i++) {
			m_cellCenters[n][i] = m_minCoords[n] + (0.5 + i) * m_cellSpacings[n];
		}

		std::cout << "  - Cell count along direction " << n << " : " << m_nCells1D[n] << "\n";

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

	// Interface area
	initializeInterfaceArea();

	// Cell volume
	initializeCellVolume();

	// Normals
	int i = 0;
	for (int n = 0; n < getDimension(); n++) {
		for (int k = -1; k <= 1; k += 2) {
			std::array<double, 3> normal = {0.0, 0.0, 0.0};
			normal[n] = k;

			m_normals[i++] = normal;
		}
	}
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
const std::vector<Adaption::Info> CartesianPatch::_update(bool trackAdaption)
{
	if (!isDirty()) {
		return std::vector<Adaption::Info>();
	}

	std::cout << ">> Updating cartesian mesh\n";

	// Reset the mesh
	reset();

	// Definition of the mesh
	createVertices();
	createCells();
	createInterfaces();

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
			cellId = cell.get_id();
		}

		adaptionData.emplace_back();
		Adaption::Info &adaptionInterfaceInfo = adaptionData.back();
		adaptionInterfaceInfo.type   = Adaption::TYPE_CREATION;
		adaptionInterfaceInfo.entity = Adaption::ENTITY_INTERFACE;
		adaptionInterfaceInfo.current.reserve(m_interfaces.size());
		for (auto &interface : m_interfaces) {
			adaptionInterfaceInfo.current.emplace_back();
			unsigned long &interfaceId = adaptionInterfaceInfo.current.back();
			interfaceId = interface.get_id();
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
void CartesianPatch::createVertices()
{
	std::cout << "  >> Creating vertices\n";

	long nTotalVertices = 1;
	for (int n = 0; n < getDimension(); n++) {
		nTotalVertices *= m_nVertices1D[n];
	}

	std::cout << "    - Vertex count: " << nTotalVertices << "\n";

	m_vertices.reserve(nTotalVertices);
	for (int k = 0; (isThreeDimensional()) ? (k < m_nVertices1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (int j = 0; j < m_nVertices1D[Vertex::COORD_Y]; j++) {
			for (int i = 0; i < m_nVertices1D[Vertex::COORD_X]; i++) {
				long id_vertex = getVertexLinearId(i, j, k);
				Patch::createVertex(id_vertex);
				Vertex &vertex = m_vertices[id_vertex];

				// Coordinate
				std::array<double, 3> coords;
				coords[Vertex::COORD_X] = m_vertexCoords[Vertex::COORD_X][i];
				coords[Vertex::COORD_Y] = m_vertexCoords[Vertex::COORD_Y][j];
				if (isThreeDimensional()) {
					coords[Vertex::COORD_Z] = m_vertexCoords[Vertex::COORD_Z][k];
				} else {
					coords[Vertex::COORD_Z] = 0.0;
				}

				vertex.setCoords(coords);
			}
		}
	}
}

/*!
	Creates the cells of the patch.
*/
void CartesianPatch::createCells()
{
	std::cout << "  >> Creating cells\n";

	// Info on the cells
	ElementInfo::Type cellType;
	if (isThreeDimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(cellType);
	const int &nCellVertices = cellTypeInfo.nVertices;

	// Count the cells
	long nTotalCells = 1;
	for (int n = 0; n < getDimension(); n++) {
		nTotalCells *= m_nCells1D[n];
	}

	std::cout << "    - Cell count: " << nTotalCells << "\n";

	m_cells.reserve(nTotalCells);

	// Create the cells
	for (int k = 0; (isThreeDimensional()) ? (k < m_nCells1D[Vertex::COORD_Z]) : (k <= 0); k++) {
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; j++) {
			for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; i++) {
				long id_cell = getCellLinearId(i, j, k);
				Patch::createCell(id_cell);
				Cell &cell = m_cells[id_cell];

				// Initialize the cell
				cell.initialize(cellType, 1);

				// Interior flag
				cell.setInterior(true);

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
void CartesianPatch::createInterfaces()
{
	std::cout << "  >> Creating interfaces\n";

	// Count the interfaces
	long nTotalInterfaces = 0;
	for (int n = 0; n < getDimension(); n++) {
		std::array<int, 3> interfaceCount1D = getInterfaceCountDirection(n);

		int nDirectionInterfaces = 1;
		for (int n = 0; n < getDimension(); n++) {
			nDirectionInterfaces *= interfaceCount1D[n];
		}
		nTotalInterfaces += nDirectionInterfaces;
	}

	std::cout << "    - Interface count: " << nTotalInterfaces << "\n";

	// Create the interfaces
	m_interfaces.reserve(nTotalInterfaces);
	for (int n = 0; n < getDimension(); n++) {
		createInterfacesDirection(n);
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
void CartesianPatch::createInterfacesDirection(const int &direction)
{
	std::cout << "  >> Creating interfaces normal to direction " << direction << "\n";

	// Info on the interfaces
	ElementInfo::Type interfaceType;
	if (isThreeDimensional()) {
		interfaceType = ElementInfo::PIXEL;
	} else {
		interfaceType = ElementInfo::LINE;
	}

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
				long id_interface = Patch::createInterface();
				Interface &interface = m_interfaces[id_interface];

				// Interface type
				if (isThreeDimensional()) {
					interface.setType(ElementInfo::PIXEL);
				} else {
					interface.setType(ElementInfo::LINE);
				}

				// Owner
				std::array<int, 3> ownerIJK;
				for (int n = 0; n < 3; n++) {
					ownerIJK[n] = counters[n];
				}
				if (counters[direction] > 0) {
					ownerIJK[direction] -= 1;
				}
				Cell &owner = m_cells[getCellLinearId(ownerIJK)];

				int ownerFace = 2 * direction;
				if (counters[direction] == 0) {
					ownerFace++;
				}

				interface.setOwner(owner.get_id(), ownerFace);
				owner.setInterface(ownerFace, 0, interface.get_id());

				// Neighbour
				if (counters[direction] != 0 && counters[direction] != interfaceCount1D[direction] - 1) {
					std::array<int, 3> neighIJK;
					for (int n = 0; n < 3; n++) {
						neighIJK[n] = counters[n];
					}

					Cell &neigh = m_cells[getCellLinearId(neighIJK)];

					int neighFace = 2 * direction + 1;

					interface.setNeigh(neigh.get_id(), neighFace);
					neigh.setInterface(neighFace, 0, interface.get_id());
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
	if (!isPointInside(point)) {
		return Element::NULL_ELEMENT_ID;
	}

	return getCellLinearId(locatePointCartesian(point));
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

	ijk[0] = std::floor(point[Vertex::COORD_X] - m_minCoords[Vertex::COORD_X]) / m_cellSpacings[Vertex::COORD_X];
	ijk[1] = std::floor(point[Vertex::COORD_Y] - m_minCoords[Vertex::COORD_Y]) / m_cellSpacings[Vertex::COORD_Y];
	if (isThreeDimensional()) {
		ijk[2] = std::floor(point[Vertex::COORD_Z] - m_minCoords[Vertex::COORD_Z]) / m_cellSpacings[Vertex::COORD_Z];
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

	if (isThreeDimensional()) {
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
	id[2] = idx / ijPlane;
	id[1] = (idx - id[2] * ijPlane) / m_nCells1D[0];

	return id;
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
	id[2] = idx / ijPlane;
	id[1] = (idx - id[2] * ijPlane) / m_nVertices1D[0];

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

	Patch::translate(translation);
}

/*!
	Scales the patch.

	\param[in] scaling is the scaling factor vector
 */
void CartesianPatch::scale(std::array<double, 3> scaling)
{
	for (int n = 1; n < 3; ++n) {
		m_minCoords[n] *= scaling[n];
		m_maxCoords[n] *= scaling[n];
		for (int i = 1; i < m_nVertices1D[n]; ++i) {
			m_vertexCoords[n][i] *= scaling[n];
			m_cellCenters[n][i]  *= scaling[n];
		}

		m_cellSpacings[n] *= scaling[n];
	}

	initializeInterfaceArea();

	initializeCellVolume();

	Patch::scale(scaling);
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
	@}
*/

}
