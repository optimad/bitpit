//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <sstream>

#include "patch.hpp"

namespace pman {

/*!
    \ingroup    PatchMan
    @{
	\class      Patch
	\brief      The Patch class provides an interface for defining patches.

	Patch is the base class for defining patches like .
*/

/*!
	Creates a new patch.
*/
Patch::Patch(const int &id, const int &dimension)
	: m_dirty(true)
{
	set_id(id) ;
	set_dimension(dimension);

	std::ostringstream convert;
	convert << get_id();
	set_name(convert.str());
}

/*!
	Destroys the patch.
*/
Patch::~Patch()
{
	reset();
}

/*!
	Gets a pointer to the the opposite normal.

	\param normal is a pointer to the normal
	\result A pointer to the opposite normal.
 */
std::array<double, 3> & Patch::get_opposite_normal(std::array<double, 3> &normal)
{
	return _get_opposite_normal(normal);
}

/*!
	Updates the mesh

	\result Returns true if the mesh was updated, false otherwise.
*/
bool Patch::update()
{
	std::vector<uint32_t> mapper;

	return update(mapper);
}

/*!
	Updates the mesh

	\result Returns true if the mesh was updated, false otherwise.
*/
bool Patch::update(std::vector<uint32_t> &cellMapping)
{
	// Update the mesh
	bool updated = _update(cellMapping);
	set_dirty(false);
	if (!updated) {
		return updated;
	}

	// Create the output data structures
	initialize_output();

	return updated;
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void Patch::mark_cell_for_refinement(const long &id)
{
	bool updated = _mark_cell_for_refinement(id);

	set_dirty(updated);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void Patch::mark_cell_for_coarsening(const long &id)
{
	bool updated = _mark_cell_for_coarsening(id);

	set_dirty(updated);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void Patch::enable_cell_balancing(const long &id, bool enabled)
{
	bool updated = _enable_cell_balancing(id, enabled);

	set_dirty(updated);
}

/*!
	Resest the patch.
*/
void Patch::reset()
{
	reset_vertices();
	reset_cells();
	reset_interfaces();
	reset_interfaces();
}

/*!
	Resest the vertices of the patch.
*/
void Patch::reset_vertices()
{
	m_vertices.clear();
	PiercedVector<Node>().swap(m_vertices);

	for (auto &cell : m_cells) {
		cell.unset_connect();
	}
}

/*!
	Resest the cells of the patch.
*/
void Patch::reset_cells()
{
	m_cells.clear();
	PiercedVector<Cell>().swap(m_cells);

	for (auto &interface : m_interfaces) {
		interface.unset_neigh();
		interface.unset_owner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void Patch::reset_interfaces()
{
	m_interfaces.clear();
	PiercedVector<Interface>().swap(m_interfaces);

	for (auto &cell : m_cells) {
		cell.unset_interfaces();
	}
}

/*!
	Resest the output manager of the patch.
*/
void Patch::reset_output()
{
	m_output_manager->Delete();
}

/*!
	Initializes output dataset.
*/
void Patch::initialize_output()
{
	long nVertices = m_vertices.size();
	long nCells = m_cells.size();

	// Create object
	m_output_manager = vtkSmartPointer<OutputManager>::New();
	m_output_manager->initialize(nCells, nVertices);

	// Vertices
	for (auto &vertex : m_vertices) {
		m_output_manager->InsertNextVertex(vertex);
	}

	// Cells
	for (auto &cell : m_cells) {
		m_output_manager->InsertNextCell(cell);
	}

	// Fields
	m_output_manager->resetFields();

	// Finalize
	m_output_manager->finalize();
}

/*!
	Writes the mesh to filename specified in input.

	\param filename the filename where the mesh will be written to
*/
void Patch::write_mesh(std::string filename)
{
	m_output_manager->write(filename);
}

/*!
	Writes the mesh a filename with the same name of the mesh
*/
void Patch::write_mesh()
{
	m_output_manager->write(get_name());
}

/*!
	Writes a field defined on the patch.

	\param name is the name of the field
	\param type is the type of the field, a field can be defined either
	on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::write_field(std::string name, int type, std::vector<double> values)
{
	write_field(get_name(), name, type, values);
}

/*!
	Writes a field defined on the patch using the specified filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param type is the type of the field, a field can be defined either
	on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::write_field(std::string filename, std::string name, int type, std::vector<double> values)
{
	int index = m_output_manager->addField(type, name.c_str());
	m_output_manager->addFieldValues(type, index, values.data());

	m_output_manager->write(filename);

	m_output_manager->resetFields();
}

/*!
	Writes a field defined on the cells of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_cell_field(std::string name, std::vector<double> values)
{
	write_cell_field(get_name(), name, values);
}

/*!
	Writes a field defined on the cells of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_cell_field(std::string filename, std::string name, std::vector<double> values)
{
	int index = m_output_manager->addField(OutputManager::DATA_ON_CELL, name.c_str());
	m_output_manager->addFieldValues(OutputManager::DATA_ON_CELL, index, values.data());

	m_output_manager->write(filename);

	m_output_manager->resetFields();
}

/*!
	Writes a field defined on the vertices of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_vertex_field(std::string name, std::vector<double> values)
{
	write_vertex_field(get_name(), name, values);
}

/*!
	Writes a field defined on the vertices of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::write_vertex_field(std::string filename, std::string name, std::vector<double> values)
{
	int index = m_output_manager->addField(OutputManager::DATA_ON_CELL, name.c_str());
	m_output_manager->addFieldValues(OutputManager::DATA_ON_CELL, index, values.data());

	m_output_manager->write(filename);

	m_output_manager->resetFields();
}

/*!
	Gets the output dataset.

	\return A pointer to the output datase
*/
OutputManager & Patch::get_output_manager()
{
	return *m_output_manager;

}

/*!
	Flags the mesh for update.

	\param dirty if true, then mesh is informed that something in the mesh
	definition has changed and thus the current data structures are not
	valid anymore.
*/
void Patch::set_dirty(bool dirty)
{
	if (m_dirty == dirty) {
		return;
	}

	m_dirty = dirty;
}

/*!
	Returns true if the the mesh needs to update its data strucutres.

	\return This method returns true to indicate the mesh needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool Patch::is_dirty() const
{
	return m_dirty;
}

/*!
	Sets the ID of the patch.

	\param id the ID of the patch
*/
void Patch::set_id(int id)
{
	m_id = id;
}

/*!
	Gets the ID of the patch.

	\return The ID of the patch
*/
int Patch::get_id() const
{
	return m_id;
}

/*!
	Sets the dimension of the patch.

	\param dimension the dimension of the patch
*/
void Patch::set_dimension(int dimension)
{
	m_dimension = dimension;
}

/*!
	Gets the dimension of the patch.

	\return The dimension of the patch
*/
int Patch::get_dimension() const
{
	return m_dimension;
}

/*!
	Returns true if the mesh is a three-dimensional mesh.

	\return This method returns true to indicate the mesh is three-dimensional
*/
bool Patch::is_three_dimensional() const
{
	return (m_dimension == 3);
}

/*!
	Sets the name of the patch.

	\param id the name of the patch
*/
void Patch::set_name(std::string name)
{
	m_name = name;
}

/*!
	Gets the name of the patch.

	\return The name of the patch
*/
std::string Patch::get_name() const
{
	return m_name;
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long Patch::get_vertex_count() const
{
	return m_vertices.size();
}

/*!
	Gets the nodes owned by the patch.

	\return The nodes owned by the patch.
*/
PiercedVector<Node> & Patch::vertices()
{
	return m_vertices;
}

/*!
	Gets the vertex with the specified id.

	\param id is the id of the requested vertex
	\return The vertex with the specified id.
*/
Node & Patch::get_vertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long Patch::get_cell_count() const
{
	return m_cells.size();
}

/*!
	Gets the cells owned by the patch.

	\return The cells owned by the patch.
*/
PiercedVector<Cell> & Patch::cells()
{
	return m_cells;
}

/*!
	Gets the cell with the specified id.

	\param id is the id of the requested cell
	\return The cell with the specified id.
*/
Cell & Patch::get_cell(const long &id)
{
	return m_cells[id];
}

/*!
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
long Patch::get_interface_count() const
{
	return m_interfaces.size();
}

/*!
	Gets the interfaces owned by the patch.

	\return The interfaces owned by the patch.
*/
PiercedVector<Interface> & Patch::interfaces()
{
	return m_interfaces;
}

/*!
	Gets the interface with the specified id.

	\param id is the id of the requested interface
	\return The interface with the specified id.
*/
Interface & Patch::get_interface(const long &id)
{
	return m_interfaces[id];
}

/*!
	Sorts the internal storage for cells, vertices and interfaces in
	ascending id order.
*/
void Patch::sort()
{
	m_vertices.sort();
	m_cells.sort();
	m_interfaces.sort();
}

/*!
	Requests the patch to compact the data structures and reduce its capacity
	to fit its size.

	The request is non-binding, and after the function call the patch can
	still occupy more memory than it actually needs.
*/
void Patch::squeeze()
{
	m_vertices.squeeze();
	m_cells.squeeze();
	m_interfaces.squeeze();
}

}

/*! @} */
