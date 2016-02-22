//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "patch_cartesian.hpp"

#include <bitpit_common.hpp>
#include <math.h>

namespace pman {

/*!
	\ingroup PatchMan
	@{
*/

/*!
	\class PatchCartesian

	\brief The PatchCartesian defines a Cartesian patch.

	PatchCartesian defines a Cartesian patch.
*/

const int PatchCartesian::SPACE_MAX_DIM = 3;

/*!
	Creates a new patch.
*/
PatchCartesian::PatchCartesian(const int &id, const int &dimension,
					 std::array<double, 3> origin, double length, double dh)
	: Patch(id, dimension)
{
	std::cout << ">> Initializing cartesian mesh\n";

	// Info sulle celle
	m_cellSize.resize(dimension);
	m_minCoord.resize(dimension);
	m_nCells1D.resize(dimension);
	for (int n = 0; n < dimension; n++) {
		// Dimensioni della cella
		m_cellSize[n] = dh;

		// Numero di celle
		m_nCells1D[n] = (int) ceil(length / m_cellSize[n]);

		std::cout << "  - Cell count along direction " << n << " : " << m_nCells1D[n] << "\n";

		// Minima coordinata del dominio
		m_minCoord[n] = origin[n] - 0.5 * (m_nCells1D[n] * m_cellSize[n]);
	}

	m_cell_volume = m_cellSize[Vertex::COORD_X] * m_cellSize[Vertex::COORD_Y];
	if (is_three_dimensional()) {
		m_cell_volume *= m_cellSize[Vertex::COORD_Z];
	}

	// Info sui vertici
	m_nVertices1D.resize(dimension);
	for (int n = 0; n < dimension; n++) {
		if (!is_three_dimensional() && n == Vertex::COORD_Z) {
			m_nVertices1D[n] = 0;
		} else {
			m_nVertices1D[n] = m_nCells1D[n] + 1;
		}
	}

	// Info sulle interfacce
	m_x_nInterfaces1D.resize(dimension);
	for (int n = 0; n < dimension; n++) {
		m_x_nInterfaces1D[n] = m_nCells1D[n];
		if (n == Vertex::COORD_X) {
			m_x_nInterfaces1D[n]++;
		}
	}

	m_y_nInterfaces1D.resize(dimension);
	for (int n = 0; n < dimension; n++) {
		m_y_nInterfaces1D[n] = m_nCells1D[n];
		if (n == Vertex::COORD_Y) {
			m_y_nInterfaces1D[n]++;
		}
	}


	if (is_three_dimensional()) {
		m_z_nInterfaces1D.resize(dimension);
		for (int n = 0; n < dimension; n++) {
			m_z_nInterfaces1D[n] = m_nCells1D[n];
			if (n == Vertex::COORD_Z) {
				m_z_nInterfaces1D[n]++;
			}
		}
	}

	m_x_interface_area = m_cellSize[Vertex::COORD_Y];
	m_y_interface_area = m_cellSize[Vertex::COORD_X];
	if (is_three_dimensional()) {
		m_x_interface_area *= m_cellSize[Vertex::COORD_Z];
		m_y_interface_area *= m_cellSize[Vertex::COORD_Z];
		m_z_interface_area  = m_cellSize[Vertex::COORD_X] * m_cellSize[Vertex::COORD_Y];
	}

	for (int i = 0; i < dimension; i++) {
		for (int n = -1; n <= 1; n += 2) {
			std::array<double, 3> normal = {0.0, 0.0, 0.0};
			normal[i] = n;

			m_normals.push_back(normal);
		}
	}
}

/*!
	Destroys the patch.
*/
PatchCartesian::~PatchCartesian()
{
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double PatchCartesian::eval_cell_volume(const long &id)
{
	BITPIT_UNUSED(id);

	return m_cell_volume;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double PatchCartesian::eval_cell_size(const long &id)
{
	BITPIT_UNUSED(id);

	double cellSize = 0;
	for(int i = 0; i < get_dimension(); ++i) {
		cellSize *= m_cellSize[i];
	}
	cellSize = pow(cellSize, 1. / get_dimension());

	return cellSize;
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double PatchCartesian::eval_interface_area(const long &id)
{
	long offset_x = 1;
	for (int i = 0; i < get_dimension(); ++i) {
		offset_x *= m_x_nInterfaces1D[i];
	}

	if (id < offset_x) {
		return m_x_interface_area;
	}

	long offset_y = 1;
	for (int i = 0; i < get_dimension(); ++i) {
		offset_y *= m_y_nInterfaces1D[i];
	}

	if (id < (offset_x + offset_y)) {
		return m_y_interface_area;
	}

	return m_z_interface_area;
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> PatchCartesian::eval_interface_normal(const long &id)
{
	const Interface &interface = get_interface(id);
	int ownerFace = interface.get_owner_face();

	return m_normals[ownerFace];
}

/*!
	Updates the patch.

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> PatchCartesian::_update(bool trackAdaption)
{
	if (!is_dirty()) {
		return std::vector<Adaption::Info>();
	}

	std::cout << ">> Updating cartesian mesh\n";

	// Reset the mesh
	reset();

	// Definition of the mesh
	create_vertices();
	create_cells();
	create_interfaces();

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
void PatchCartesian::create_vertices()
{
	std::cout << "  >> Creating vertices\n";

	// Definition of the vertices
	m_x = std::vector<double>(m_nVertices1D[Vertex::COORD_X]);
	for (int i = 0; i < m_nVertices1D[Vertex::COORD_X]; i++) {
		m_x[i] = m_minCoord[Vertex::COORD_X] + i * m_cellSize[Vertex::COORD_X];
	}

	m_y = std::vector<double>(m_nVertices1D[Vertex::COORD_Y]);
	for (int j = 0; j < m_nVertices1D[Vertex::COORD_Y]; j++) {
		m_y[j] = m_minCoord[Vertex::COORD_Y] + j * m_cellSize[Vertex::COORD_Y];
	}

	if (is_three_dimensional()) {
		m_z = std::vector<double>(m_nVertices1D[Vertex::COORD_Z]);
		for (int k = 0; k < m_nVertices1D[Vertex::COORD_Z]; k++) {
			m_z[k] = m_minCoord[Vertex::COORD_Z] + k * m_cellSize[Vertex::COORD_Z];
		}
	} else {
		m_z = std::vector<double>(0);
	}

	long nTotalVertices = 1;
	for (int n = 0; n < get_dimension(); n++) {
		nTotalVertices *= m_nVertices1D[n];
	}

	std::cout << "    - Vertex count: " << nTotalVertices << "\n";

	m_vertices.reserve(nTotalVertices);

	for (int i = 0; i < m_nVertices1D[Vertex::COORD_X]; i++) {
		for (int j = 0; j < m_nVertices1D[Vertex::COORD_Y]; j++) {
			for (int k = 0; (is_three_dimensional()) ? (k < m_nVertices1D[Vertex::COORD_Z]) : (k <= 0); k++) {
				long id_vertex = vertex_ijk_to_id(i, j, k);
				Patch::create_vertex(id_vertex);
				Vertex &vertex = m_vertices[id_vertex];

				// Coordinate
				std::array<double, 3> coords;
				coords[Vertex::COORD_X] = m_x[i];
				coords[Vertex::COORD_Y] = m_y[j];
				if (is_three_dimensional()) {
					coords[Vertex::COORD_Z] = m_z[k];
				} else {
					coords[Vertex::COORD_Z] = 0.0;
				}

				vertex.set_coords(coords);
			}
		}
	}
}

/*!
	Creates the cells of the patch.
*/
void PatchCartesian::create_cells()
{
	std::cout << "  >> Creating cells\n";

	// Info on the cells
	ElementInfo::Type cellType;
	if (is_three_dimensional()) {
		cellType = ElementInfo::VOXEL;
	} else {
		cellType = ElementInfo::PIXEL;
	}

	const ElementInfo &cellTypeInfo = ElementInfo::get_element_info(cellType);
	const int &nCellVertices = cellTypeInfo.nVertices;

	// Count the cells
	long nTotalCells = 1;
	for (int n = 0; n < get_dimension(); n++) {
		nTotalCells *= m_nCells1D[n];
	}

	std::cout << "    - Cell count: " << nTotalCells << "\n";

	m_cells.reserve(nTotalCells);

	// Create the cells
	std::array<double, 3> centroid = {0.0, 0.0, 0.0};
	for (int i = 0; i < m_nCells1D[Vertex::COORD_X]; i++) {
		centroid[Vertex::COORD_X] = 0.5 * (m_x[i] + m_x[i+1]);
		for (int j = 0; j < m_nCells1D[Vertex::COORD_Y]; j++) {
			centroid[Vertex::COORD_Y] = 0.5 * (m_y[j] + m_y[j+1]);
			for (int k = 0; (is_three_dimensional()) ? (k < m_nCells1D[Vertex::COORD_Z]) : (k <= 0); k++) {
				long id_cell = cell_ijk_to_id(i, j, k);
				Patch::create_cell(id_cell);
				Cell &cell = m_cells[id_cell];

				// Initialize the cell
				cell.initialize(cellType, 1);

				// Interior flag
				cell.set_interior(true);

				// Connettività
				cell.set_vertex(0, vertex_ijk_to_id(i,     j,     k));
				cell.set_vertex(1, vertex_ijk_to_id(i + 1, j,     k));
				cell.set_vertex(2, vertex_ijk_to_id(i,     j + 1, k));
				cell.set_vertex(3, vertex_ijk_to_id(i + 1, j + 1, k));
				if (is_three_dimensional()) {
					cell.set_vertex(4, vertex_ijk_to_id(i,     j,     k + 1));
					cell.set_vertex(5, vertex_ijk_to_id(i + 1, j,     k + 1));
					cell.set_vertex(6, vertex_ijk_to_id(i,     j + 1, k + 1));
					cell.set_vertex(7, vertex_ijk_to_id(i + 1, j + 1, k + 1));
				}
			}
		}
	}
}

/*!
	Creates the interfaces of the patch.
*/
void PatchCartesian::create_interfaces()
{
	std::cout << "  >> Creating interfaces\n";

	// Count the interfaces
	long nTotalInterfaces = 0;
	nTotalInterfaces += count_interfaces_direction(Vertex::COORD_X);
	nTotalInterfaces += count_interfaces_direction(Vertex::COORD_Y);
	if (is_three_dimensional()) {
		nTotalInterfaces += count_interfaces_direction(Vertex::COORD_Z);
	}

	std::cout << "    - Interface count: " << nTotalInterfaces << "\n";

	// Create the interfaces
	m_interfaces.reserve(nTotalInterfaces);

	create_interfaces_direction(Vertex::COORD_X);
	create_interfaces_direction(Vertex::COORD_Y);
	if (is_three_dimensional()) {
		create_interfaces_direction(Vertex::COORD_Z);
	}
}

/*!
	Counts the interfaces normal to the given direction.

	\param direction the method will count the interfaces normal to this
	                 direction
*/
int PatchCartesian::count_interfaces_direction(const Vertex::Coordinate &direction)
{
	std::vector<int> *interfaceCount1D;
	switch (direction)  {

	case Vertex::COORD_X:
		interfaceCount1D = &m_x_nInterfaces1D;
		break;

	case Vertex::COORD_Y:
		interfaceCount1D = &m_y_nInterfaces1D;
		break;

	case Vertex::COORD_Z:
		interfaceCount1D = &m_z_nInterfaces1D;
		break;

	}

	int nInterfaces = 1;
	for (int n = 0; n < get_dimension(); n++) {
		nInterfaces *= (*interfaceCount1D)[n];
	}

	return nInterfaces;
}

/*!
	Creates the interfaces normal to the given direction.

	\param direction the method will creat the interfaces normal to this
	                 direction
*/
void PatchCartesian::create_interfaces_direction(const Vertex::Coordinate &direction)
{
	std::cout << "  >> Creating interfaces normal to direction " << direction << "\n";

	// Info on the interfaces
	ElementInfo::Type interfaceType;
	if (is_three_dimensional()) {
		interfaceType = ElementInfo::PIXEL;
	} else {
		interfaceType = ElementInfo::LINE;
	}

	const ElementInfo &interfaceTypeInfo = ElementInfo::get_element_info(interfaceType);
	const int nInterfaceVertices = interfaceTypeInfo.nVertices;

	double *area;
	std::vector<int> *interfaceCount1D;
	switch (direction)  {

	case Vertex::COORD_X:
		area = &m_x_interface_area;
		interfaceCount1D = &m_x_nInterfaces1D;
		break;

	case Vertex::COORD_Y:
		area = &m_y_interface_area;
		interfaceCount1D = &m_y_nInterfaces1D;
		break;

	case Vertex::COORD_Z:
		area = &m_z_interface_area;
		interfaceCount1D = &m_z_nInterfaces1D;
		break;

	}

	// Counters
	int counters[SPACE_MAX_DIM] = {0, 0, 0};
	int &i = counters[Vertex::COORD_X];
	int &j = counters[Vertex::COORD_Y];
	int &k = counters[Vertex::COORD_Z];

	// Creation of the interfaces
	for (i = 0; i < (*interfaceCount1D)[Vertex::COORD_X]; i++) {
		for (j = 0; j < (*interfaceCount1D)[Vertex::COORD_Y]; j++) {
			for (k = 0; (is_three_dimensional()) ? (k < (*interfaceCount1D)[Vertex::COORD_Z]) : (k <= 0); k++) {
				long id_interface = interface_nijk_to_id(direction, i, j, k);
				Patch::create_interface(id_interface);
				Interface &interface = m_interfaces[id_interface];

				// Interface type
				if (is_three_dimensional()) {
					interface.set_type(ElementInfo::PIXEL);
				} else {
					interface.set_type(ElementInfo::LINE);
				}

				// Owner
				int ownerIJK[SPACE_MAX_DIM];
				for (int n = 0; n < SPACE_MAX_DIM; n++) {
					ownerIJK[n] = counters[n];
				}
				if (counters[direction] > 0) {
					ownerIJK[direction] -= 1;
				}
				Cell &owner = m_cells[cell_ijk_to_id(ownerIJK)];

				int ownerFace = 2 * direction;
				if (counters[direction] == 0) {
					ownerFace++;
				}

				interface.set_owner(owner.get_id(), ownerFace);
				owner.set_interface(ownerFace, 0, interface.get_id());

				// Neighbour
				if (counters[direction] != 0 && counters[direction] != (*interfaceCount1D)[direction] - 1) {
					int neighIJK[SPACE_MAX_DIM];
					for (int n = 0; n < SPACE_MAX_DIM; n++) {
						neighIJK[n] = counters[n];
					}

					Cell &neigh = m_cells[cell_ijk_to_id(neighIJK)];

					int neighFace = 2 * direction + 1;

					interface.set_neigh(neigh.get_id(), neighFace);
					neigh.set_interface(neighFace, 0, interface.get_id());
				} else {
					interface.unset_neigh();
				}

				// Connectivity
				std::unique_ptr<long[]> connect = std::unique_ptr<long[]>(new long[nInterfaceVertices]);
				if (direction == Vertex::COORD_X) {
					connect[0] = vertex_ijk_to_id(i, j,     k);
					connect[1] = vertex_ijk_to_id(i, j + 1, k);
					if (is_three_dimensional()) {
						connect[2] = vertex_ijk_to_id(i, j + 1, k + 1);
						connect[3] = vertex_ijk_to_id(i, j,     k + 1);
					}
				} else if (direction == Vertex::COORD_Y) {
					connect[0] = vertex_ijk_to_id(i,     j,     k);
					connect[1] = vertex_ijk_to_id(i + 1, j,     k);
					if (is_three_dimensional()) {
						connect[2] = vertex_ijk_to_id(i + 1, j, k + 1);
						connect[3] = vertex_ijk_to_id(i,     j, k + 1);
					}
				} else if (direction == Vertex::COORD_Z) {
					connect[0] = vertex_ijk_to_id(i,     j,     k);
					connect[1] = vertex_ijk_to_id(i + 1, j,     k);
					if (is_three_dimensional()) {
						connect[2] = vertex_ijk_to_id(i + 1, j + 1, k);
						connect[3] = vertex_ijk_to_id(i,     j + 1, k);
					}
				}

				interface.set_connect(std::move(connect));
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
bool PatchCartesian::_mark_cell_for_refinement(const long &id)
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
bool PatchCartesian::_mark_cell_for_coarsening(const long &id)
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
bool PatchCartesian::_enable_cell_balancing(const long &id, bool enabled)
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(enabled);

	return false;
}

/*!
	Converts the cell (i, j, k) notation to a single index notation
*/
long PatchCartesian::cell_ijk_to_id(const int &i, const int &j, const int &k) const
{
	long id = i;
	id += m_nCells1D[Vertex::COORD_X] * j;
	if (get_dimension() == 3) {
		id += m_nCells1D[Vertex::COORD_Y] * m_nCells1D[Vertex::COORD_X] * k;
	}

	return id;
}

/*!
	Converts the cell (i, j, k) notation to a single index notation
*/
long PatchCartesian::cell_ijk_to_id(const int ijk[]) const
{
	return cell_ijk_to_id(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}


/*!
	Converts the vertex (i, j, k) notation to a single index notation
*/
long PatchCartesian::vertex_ijk_to_id(const int &i, const int &j, const int &k) const
{
	long id = i;
	id += m_nVertices1D[Vertex::COORD_X] * j;
	if (get_dimension() == 3) {
		id += m_nVertices1D[Vertex::COORD_Y] * m_nVertices1D[Vertex::COORD_X] * k;
	}

	return id;
}

/*!
	Converts the vertex (i, j, k) notation to a single index notation
*/
long PatchCartesian::vertex_ijk_to_id(const int ijk[]) const
{
	return vertex_ijk_to_id(ijk[Vertex::COORD_X], ijk[Vertex::COORD_Y], ijk[Vertex::COORD_Z]);
}

/*!
	Converts the interface (normal, i, j, k) notation to a single index notation
*/
long PatchCartesian::interface_nijk_to_id(const int &normal, const int &i, const int &j, const int &k) const
{
	const int *nInterfaces;
	switch (normal) {

	case Vertex::COORD_X:
		nInterfaces = m_x_nInterfaces1D.data();
		break;

	case Vertex::COORD_Y:
		nInterfaces = m_y_nInterfaces1D.data();
		break;

	case Vertex::COORD_Z:
		nInterfaces = m_z_nInterfaces1D.data();
		break;

	}

	long offset = 0;
	if (normal <=  Vertex::COORD_Y) {
		long offset_y = 1;
		for (int i = 0; i < get_dimension(); ++i) {
			offset_y *= m_y_nInterfaces1D[i];
		}
		offset += offset_y;
	}
	if (normal <=  Vertex::COORD_X) {
		long offset_x = 1;
		for (int i = 0; i < get_dimension(); ++i) {
			offset_x *= m_x_nInterfaces1D[i];
		}
		offset += offset_x;
	}

	long id = offset + i;
	id += nInterfaces[Vertex::COORD_X] * j;
	if (is_three_dimensional()) {
		id += nInterfaces[Vertex::COORD_Y] * nInterfaces[Vertex::COORD_X] * k;
	}

	return id;
}

/*!
	@}
*/

}
