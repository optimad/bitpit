//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <sstream>
#include <unordered_map>

#include "patch.hpp"

namespace pman {

// Explicit instantiation of the PiercedVectors
template class PiercedVector<Cell>;
template class PiercedVector<Interface>;
template class PiercedVector<Vertex>;

/*!
	\class Patch

	\brief The Patch class provides an interface for defining patches.

	Patch is the base class for defining patches like .
*/

/*!
	Creates a new patch.
*/
Patch::Patch(const int &id, const int &dimension)
	: m_dirty(true), m_dirty_output(true), m_output_manager(nullptr)
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

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> Patch::update(bool trackAdaption)
{
	const std::vector<Adaption::Info> adaptionInfo = _update(trackAdaption);

	m_cells.flush();
	m_interfaces.flush();
	m_vertices.flush();

	set_dirty(false);

	return adaptionInfo;
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
	PiercedVector<Vertex>().swap(m_vertices);

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
	if (m_output_manager == nullptr) {
		return;
	}

	get_output_manager().Delete();
}

/*!
	Initializes output dataset.
*/
void Patch::update_output_manager()
{
	long nVertices = m_vertices.size();
	long nCells = m_cells.size();

	// Create object
	m_output_manager = vtkSmartPointer<OutputManager>::New();
	m_output_manager->initialize(nCells, nVertices);

	// Vertices
	std::unordered_map<long, vtkIdType> vertexMap;
	for (auto &vertex : m_vertices) {
		vertexMap[vertex.get_id()] = m_output_manager->InsertNextVertex(vertex);
	}

	// Cells
	for (auto &cell : m_cells) {
		m_output_manager->InsertNextCell(cell, vertexMap);
	}

	// Fields
	m_output_manager->resetFields();

	// Finalize
	m_output_manager->finalize();

	// The output is not dirty anymore
	m_dirty_output = false;
}

/*!
	Writes the mesh to filename specified in input.

	\param filename the filename where the mesh will be written to
*/
void Patch::write_mesh(std::string filename)
{
	get_output_manager().write(filename);
}

/*!
	Writes the mesh a filename with the same name of the mesh
*/
void Patch::write_mesh()
{
	get_output_manager().write(get_name());
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
	OutputManager &outputManager = get_output_manager();

	int index = outputManager.addField(type, name.c_str());
	outputManager.addFieldValues(type, index, values.data());

	outputManager.write(filename);

	outputManager.resetFields();
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
	OutputManager &outputManager = get_output_manager();

	int index = outputManager.addField(OutputManager::DATA_ON_CELL, name.c_str());
	outputManager.addFieldValues(OutputManager::DATA_ON_CELL, index, values.data());

	outputManager.write(filename);

	outputManager.resetFields();
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
	OutputManager &outputManager = get_output_manager();

	int index = outputManager.addField(OutputManager::DATA_ON_CELL, name.c_str());
	outputManager.addFieldValues(OutputManager::DATA_ON_CELL, index, values.data());

	outputManager.write(filename);

	outputManager.resetFields();
}

/*!
	Gets the output dataset.

	\return A pointer to the output datase
*/
OutputManager & Patch::get_output_manager()
{
	if (is_output_dirty()) {
		update_output_manager();
	}

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
	if (m_dirty) {
		m_dirty_output = true;
	}
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
	Returns true if the the output needs to update its data strucutres.

	\return This method returns true to indicate the output needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool Patch::is_output_dirty() const
{
	return m_dirty_output;
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
PiercedVector<Vertex> & Patch::vertices()
{
	return m_vertices;
}

/*!
	Gets the vertex with the specified id.

	\param id is the id of the requested vertex
	\return The vertex with the specified id.
*/
Vertex & Patch::get_vertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Creates a new vertex with the specified id.

	\param id is the id of the new vertex
*/
long Patch::create_vertex(const long &id)
{
	m_vertices.emplace(id);

	return id;
}

/*!
	Creates a new vertex.
*/
long Patch::create_vertex()
{
	long id;
	if (m_unusedVertexIds.empty()) {
		id = m_vertices.size();
	} else {
		id = m_unusedVertexIds.front();
		m_unusedVertexIds.pop_front();
	}

	return create_vertex(id);
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
*/
void Patch::delete_vertex(const long &id, bool delayed)
{
	m_vertices.erase(id, delayed);
	m_unusedVertexIds.push_back(id);
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
	Creates a new cell with the specified id.

	\param id is the id of the new cell
	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::create_cell(const long &id, bool internal, ElementInfo::Type type)
{
	if (internal) {
		m_cells.emplace(id, type, this);
	} else {
		m_cells.emplace_back(id, type, this);
	}

	return id;
}

/*!
	Creates a new cell.

	\param internal is true if the cell is an internal cell, false otherwise
*/
long Patch::create_cell(bool internal, ElementInfo::Type type)
{
	long id;
	if (m_unusedCellIds.empty()) {
		id = m_cells.size();
	} else {
		id = m_unusedCellIds.front();
		m_unusedCellIds.pop_front();
	}

	return create_cell(id, internal, type);
}

/*!
	Deletes a cell.

	\param id is the id of the cell
*/
void Patch::delete_cell(const long &id, bool delayed)
{
	m_cells.erase(id, delayed);
	m_unusedCellIds.push_back(id);
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
	Creates a new interface with the specified id.

	\param id is the id of the new interface
*/
long Patch::create_interface(const long &id, ElementInfo::Type type)
{
	m_interfaces.emplace(id, type, this);

	return id;
}

/*!
	Creates a new interface.
*/
long Patch::create_interface(ElementInfo::Type type)
{
	long id;
	if (m_unusedInterfaceIds.empty()) {
		id = m_interfaces.size();
	} else {
		id = m_unusedInterfaceIds.front();
		m_unusedInterfaceIds.pop_front();
	}

	return create_interface(id, type);
}

/*!
	Deletes an interface.

	\param id is the id of the interface
*/
void Patch::delete_interface(const long &id, bool delayed)
{
	m_interfaces.erase(id, delayed);
	m_unusedInterfaceIds.push_back(id);
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
