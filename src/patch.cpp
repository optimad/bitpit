//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\class Patch

	\brief The Patch class provides an interface for defining patches.

	Patch is the base class for defining patches like .
*/

#include <sstream> 

#include "patch.hpp"

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
	Updates the mesh
*/
void Patch::update()
{
	update(std::vector<uint32_t>());
}

/*!
	Updates the mesh
*/
void Patch::update(const std::vector<uint32_t> &cellMapping)
{
	// Update the mesh
	_update(cellMapping);
	set_dirty(false);

	// Create the output data structures
	output_initialize();
}

/*!
	Marks a cell for refinement.

	\param cell the cell to be refined
*/
void Patch::mark_for_refinement(Cell &cell)
{
	_mark_for_refinement(cell);

	set_dirty(true);
}

/*!
	Marks a cell for refinement.

	\param cellId index to the cell to be refined
*/
void Patch::mark_for_refinement(const int &cellId)
{
	mark_for_refinement(m_cells[cellId]);
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
	std::vector<Node>().swap(m_vertices);

	for (unsigned int n = 0; n < m_cells.size(); n++) {
		m_cells[n].unset_connect();
	}
}

/*!
	Resest the cells of the patch.
*/
void Patch::reset_cells()
{
	m_cells.clear();
	std::vector<Cell>().swap(m_cells);

	for (unsigned int n = 0; n < m_interfaces.size(); n++) {
		m_interfaces[n].unset_neigh();
		m_interfaces[n].unset_owner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void Patch::reset_interfaces()
{
	m_interfaces.clear();
	std::vector<Interface>().swap(m_interfaces);

	for (unsigned int n = 0; n < m_cells.size(); n++) {
		m_cells[n].unset_interfaces();
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
void Patch::output_initialize()
{
	int nVertices = m_vertices.size();
	int nCells = m_cells.size();

	// Create object
	m_output_manager = vtkSmartPointer<OutputManager>::New();
	m_output_manager->initialize(nCells, nVertices);

	// Vertices
	for (int i = 0; i < nVertices; i++) {
		m_output_manager->InsertNextVertex(m_vertices[i], get_dimension());
	}

	// Cells
	for (int i = 0; i < nCells; i++) {
		m_output_manager->InsertNextCell(m_cells[i]);
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
void Patch::output_write(std::string filename)
{
	m_output_manager->write(filename);
}

/*!
	Writes the mesh a filename with the same name of the mesh
*/
void Patch::output_write()
{
	m_output_manager->write(get_name());
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
int Patch::get_vertex_count() const
{
	return m_vertices.size();
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
int Patch::get_cell_count() const
{
	return m_cells.size();
}

/*!
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
int Patch::get_interface_count() const
{
	return m_interfaces.size();
}
