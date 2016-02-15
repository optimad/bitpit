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

#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>

#include "bitpit_SA.hpp"

#include "patch.hpp"
#include "utils.hpp"

namespace bitpit {


/*!
	\ingroup patch
	@{
*/

/*!
	\class IndexGenerator

	\brief The IndexGenerator class allows to generate unique ids.
*/

/*!
	Creates a new generator.
*/
IndexGenerator::IndexGenerator()
	: m_next(0), m_depleted(false)
{

}

/*!
	Generates a unique index.

	If the trash is empty a new index is generated, otherwise an index taken
	from the trash is recycled.

	\return A new unique index.
*/
long IndexGenerator::generateId()
{
	long id;
	if (m_trash.empty()) {
		assert(!m_depleted);
		if (m_next == std::numeric_limits<long>::max()) {
			m_depleted = true;
		}

		id = m_next++;
	} else {
		id = m_trash.front();
		m_trash.pop_front();
	}

	return id;
}

/*!
	Gets the last assigned id.

	\return The last assigned index.
*/
long IndexGenerator::getLastId()
{
	return (m_next - 1);
}

/*!
	Trashes an index.

	A trashed index is an index no more used that can be recycled.

	\param id is the index that will be trashed
*/
void IndexGenerator::trashId(const long &id)
{
	m_trash.push_back(id);
}

/*!
	Reset the generator.
*/
void IndexGenerator::reset()
{
	m_next = 0;
	m_trash.clear();
}

/*!
	\ingroup patch
	@{
*/

/*!
	\class Patch

	\brief The Patch class provides an interface for defining patches.

	Patch is the base class for defining patches.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/

/*!
	Creates a new patch.
*/
Patch::Patch(const int &id, const int &dimension, bool expert)
	: m_nVertices(0), m_nInternals(0), m_nGhosts(0), m_nInterfaces(0),
	  m_last_internal_id(Element::NULL_ID),
	  m_first_ghost_id(Element::NULL_ID),
	  m_dirty(true), m_expert(expert), m_hasCustomTolerance(false),
	  m_rank(0), m_nProcessors(1)
#if ENABLE_MPI==1
	  , m_communicator(MPI_COMM_NULL)
#endif
{
	set_id(id) ;
	setDimension(dimension);

	std::ostringstream convert;
	convert << get_id();
	setName(convert.str());

	VTKUnstructuredGrid::setCodex(VTKFormat::APPENDED);
}

/*!
	Destroys the patch.
*/
Patch::~Patch()
{
	reset();
}

/*!
	Updates the patch

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> Patch::update(bool trackAdaption)
{
	const std::vector<Adaption::Info> adaptionInfo = _update(trackAdaption);

	m_cells.flush();
	m_interfaces.flush();
	m_vertices.flush();

	updateBoundingBox();

	setDirty(false);

	return adaptionInfo;
}

/*!
	Marks a cell for refinement.

	\param id is the id of the cell that needs to be refined
*/
void Patch::markCellForRefinement(const long &id)
{
	bool updated = _markCellForRefinement(id);

	setDirty(updated);
}

/*!
	Marks a cell for coarsening.

	\param id is the id of the cell that needs to be coarsened
*/
void Patch::markCellForCoarsening(const long &id)
{
	bool updated = _markCellForCoarsening(id);

	setDirty(updated);
}

/*!
	Enables cell balancing.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
void Patch::enableCellBalancing(const long &id, bool enabled)
{
	bool updated = _enableCellBalancing(id, enabled);

	setDirty(updated);
}

/*!
	Resest the patch.
*/
void Patch::reset()
{
	resetVertices();
	resetCells();
	resetInterfaces();
}

/*!
	Resest the vertices of the patch.
*/
void Patch::resetVertices()
{
	m_vertices.clear();
	PiercedVector<Vertex>().swap(m_vertices);
	m_vertexIdGenerator.reset();
	m_nVertices = 0;

	for (auto &cell : m_cells) {
		cell.unsetConnect();
	}
}

/*!
	Resest the cells of the patch.
*/
void Patch::resetCells()
{
	m_cells.clear();
	PiercedVector<Cell>().swap(m_cells);
	m_cellIdGenerator.reset();
	m_nInternals = 0;
	m_nGhosts = 0;

	for (auto &interface : m_interfaces) {
		interface.unsetNeigh();
		interface.unsetOwner();
	}
}

/*!
	Resest the interfaces of the patch.
*/
void Patch::resetInterfaces()
{
	m_interfaces.clear();
	PiercedVector<Interface>().swap(m_interfaces);
	m_interfaceIdGenerator.reset();
	m_nInterfaces = 0;

	for (auto &cell : m_cells) {
		cell.resetInterfaces();
	}
}

/*!
    Reserve memory for vertex storage.

    If the reserve size is smaller than the number of vertices currently stored
    within the patch no action will be taken.

    If instead, the reserve size is greater than the current number of vertices,
    reserve might cause re-location of the internal container into memory,
    potentially invalidating pointers and iterators to vertex entities.

    \param[in] nVertices size of memory reserve (in terms of number of
    vertices).
*/
bool Patch::reserveVertices(size_t nVertices)
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.reserve(nVertices);

	return true;
}

/*!
	Reserve memory for cell storage.

	If the reserve size is smaller than the number of cells currently stored
	within the patch no action will be taken.

	If instead, the reserve size is greater than the current number of cells,
	reserve might cause re-location of the internal container into memory,
	potentially invalidating pointers and iterators to cell entities.

	\param[in] nCells is size of memory reserve (in terms of number of cells).
*/
bool Patch::reserveCells(size_t nCells)
{
	if (!isExpert()) {
		return false;
	}

	m_cells.reserve(nCells);

	return true;
}

/*!
	Reserve memory for interface storage.

	If the reserve size is smaller than the number of interfaces currently
	stored within the patch no action will be taken.

	If instead, the reserve size is greater than the current number of
	interfaces, reserve might cause re-location of the internal container
	into memory, potentially invalidating pointers and iterators to cell
	entities.

	\param[in] nInterfaces is size of memory reserve (in terms of number of
	interfaces).
*/
bool Patch::reserveInterfaces(size_t nInterfaces)
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.reserve(nInterfaces);

	return true;
}

/*!
	Writes the patch to filename specified in input.

	\param filename the filename where the patch will be written to
*/
void Patch::write(std::string filename)
{
	std::string oldFilename = VTKUnstructuredGrid::getName();

	VTKUnstructuredGrid::setName(filename);
	write();
	VTKUnstructuredGrid::setName(oldFilename);
}

/*!
	Writes the patch a filename with the same name of the patch
*/
void Patch::write()
{
	// Add basic patch data
	std::vector<long> cellIndex;
	cellIndex.reserve(getCellCount());
	for (const Cell &cell : m_cells) {
		cellIndex.emplace_back(cell.get_id());
	}
	VTKUnstructuredGrid::addData<long>("cellIndex", VTKFieldType::SCALAR, VTKLocation::CELL, cellIndex);

	std::vector<long> vertexIndex;
	vertexIndex.reserve(getVertexCount());
	for (const Vertex &vertex : m_vertices) {
		vertexIndex.emplace_back(vertex.get_id());
	}
	VTKUnstructuredGrid::addData<long>("vertexIndex", VTKFieldType::SCALAR, VTKLocation::POINT, vertexIndex);

#if ENABLE_MPI==1
	std::vector<long> rank(getCellCount(), m_rank);
	VTKUnstructuredGrid::addData<long>("rank", VTKFieldType::SCALAR, VTKLocation::CELL, rank);
#endif

	// Write mesh
	VTKUnstructuredGrid::write();
}

/*!
	Writes a field defined on the patch.

	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::writeField(std::string name, VTKLocation location, std::vector<double> &values)
{
	writeField(getName(), name, location, values);
}

/*!
	Writes a field defined on the patch using the specified filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param location is the location of the field, a field can be defined
	either on the vertices of on the cells
	\param values is a vector with the values of the field
*/
void Patch::writeField(std::string filename, std::string name, VTKLocation location, std::vector<double> &values)
{
	VTKUnstructuredGrid::addData<double>(name, VTKFieldType::SCALAR, location, values);

	write(filename);

	VTKUnstructuredGrid::removeData(name);
}

/*!
	Writes a field defined on the cells of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeCellField(std::string name, std::vector<double> &values)
{
	writeCellField(getName(), name, values);
}

/*!
	Writes a field defined on the cells of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeCellField(std::string filename, std::string name, std::vector<double> &values)
{
	writeField(filename, name, VTKLocation::CELL, values);
}

/*!
	Writes a field defined on the vertices of the patch.

	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeVertexField(std::string name, std::vector<double> &values)
{
	writeVertexField(getName(), name, values);
}

/*!
	Writes a field defined on the vertices of the patch using the specified
	filename.

	\param filename is the name of the file to write
	\param name is the name of the field
	\param values is a vector with the values of the field
*/
void Patch::writeVertexField(std::string filename, std::string name, std::vector<double> &values)
{
	writeField(filename, name, VTKLocation::POINT, values);
}

/*!
	Flags the patch for update.

	\param dirty if true, then patch is informed that something in the
	patch definition has changed and thus the current data structures
	are not valid anymore.
*/
void Patch::setDirty(bool dirty)
{
	if (m_dirty == dirty) {
		return;
	}

	m_dirty = dirty;
}

/*!
	Returns true if the the patch needs to update its data strucutres.

	\return This method returns true to indicate the patch needs to update
	its data strucutres. Otherwise, it returns false.
*/
bool Patch::isDirty() const
{
	return m_dirty;
}

/*!
	Enables or disables expert mode.

	When expert mode is enabled, it will be possible to change the
	patch using low level functions (e.g., it will be possible to
	add individual cells, add vertices, delete cells, ...).

	\param expert if true, the expert mode will be enabled
*/
void Patch::setExpert(bool expert)
{
	if (isExpert() == expert) {
		return;
	}

	m_expert = expert;
}

/*!
	Checks if the expert mode is enabled.

	When expert mode is enabled, it will be possible to change the
	patch using low level functions (e.g., it will be possible to
	add individual cells, add vertices, delete cells, ...).

	\return This method returns true when the expert is enabled,
	otherwise it returns false.
*/
bool Patch::isExpert() const
{
	return m_expert;
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
void Patch::setDimension(int dimension)
{
	m_dimension = dimension;
}

/*!
	Gets the dimension of the patch.

	\return The dimension of the patch
*/
int Patch::getDimension() const
{
	return m_dimension;
}

/*!
	Returns true if the patch is a three-dimensional patch.

	\return This method returns true to indicate the patch is
	three-dimensional
*/
bool Patch::isThreeDimensional() const
{
	return (m_dimension == 3);
}

/*!
	Gets the number of vertices in the patch.

	\return The number of vertices in the patch
*/
long Patch::getVertexCount() const
{
	return m_nVertices;
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
	Gets a reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A reference to the vertex with the specified id.
*/
Vertex & Patch::getVertex(const long &id)
{
	return m_vertices[id];
}

/*!
	Gets a constant reference to the vertex with the specified id.

	\param id is the id of the requested vertex
	\return A constant reference to the vertex with the specified id.
*/
const Vertex & Patch::getVertex(const long &id) const
{
	return m_vertices[id];
}

/*!
	Returns an iterator pointing to the specified vertex.

	\result An iterator to the specified vertex.
*/
Patch::VertexIterator Patch::getVertexIterator(const long &id)
{
	return VertexIterator(m_vertices.raw_begin() + m_vertices.raw_index(id));
}

/*!
	Returns iterator pointing to the first vertex.

	\result An iterator to the first vertex.
*/
Patch::VertexIterator Patch::vertexBegin()
{
	return m_vertices.begin();
}

/*!
	Returns iterator pointing to last vertex.

	\result An iterator to the last vertex.
*/
Patch::VertexIterator Patch::vertexEnd()
{
	return m_vertices.end();
}

/*!
	Generates a new unique id for the vertices.

	\result A new unique id for the vertices.
*/
long Patch::generateVertexId()
{
	if (!isExpert()) {
		return Vertex::NULL_ID;
	}

	return m_vertexIdGenerator.generateId();
}

/*!
	Creates a new vertex with the specified id.

	\param id is the id of the new vertex
	\return An iterator pointing to the newly created vertex.
*/
Patch::VertexIterator Patch::createVertex(long id)
{
	if (id < 0) {
		id = generateVertexId();
	}

	PiercedVector<Vertex>::iterator iterator = m_vertices.reclaim(id);
	m_nVertices++;

	return iterator;
}

/*!
	Adds a new vertex with the specified id.

	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
Patch::VertexIterator Patch::addVertex(const long &id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	return createVertex(id);
}

/*!
	Adds a new vertex with the specified coordinates.

	\param coords are the coordinates of the vertex
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
Patch::VertexIterator Patch::addVertex(const std::array<double, 3> &coords, const long &id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = createVertex(id);
	Vertex &vertex = (*iterator);
	vertex.setCoords(coords);

	return iterator;
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
Patch::VertexIterator Patch::addVertex(const Vertex &source, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	VertexIterator iterator = createVertex(id);
	Vertex &vertex = (*iterator);
	id = vertex.get_id();
	vertex = source;
	vertex.set_id(id);

	return iterator;
}

/*!
	Adds the specified vertex to the patch.

	\param source is the vertex that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added vertex.
*/
Patch::VertexIterator Patch::addVertex(Vertex &&source, long id)
{
	if (!isExpert()) {
		return vertexEnd();
	}

	if (id < 0) {
		id = source.get_id();
	}

	VertexIterator iterator = createVertex(std::max(source.get_id(), id));
	Vertex &vertex = (*iterator);
	id = vertex.get_id();
	vertex = std::move(source);
	vertex.set_id(id);

	return iterator;
}

/*!
	Deletes a vertex.

	\param id is the id of the vertex
	\param delayed is true a delayed delete will be performed
*/
bool Patch::deleteVertex(const long &id, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.erase(id, delayed);
	m_vertexIdGenerator.trashId(id);
	m_nVertices--;

	return true;
}

/*!
	Deletes a list of vertices.

	\param ids are the ids of the vertices to be deleted
	\param delayed is true a delayed delete will be performed
*/
bool Patch::deleteVertices(const std::vector<long> &ids, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteVertex(*i, delayed);
	}

	return true;
}

/*!
	Counts free vertices within the patch.

	A free vertex is a vertex on a free face.

	\return The number of free vertices.
*/
long Patch::countFreeVertices() const
{
	std::unordered_set<long> freeVertices;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				continue;
			}

			std::vector<int> faceLocalConnect = cell.getFaceLocalConnect(i);
			for (unsigned int j = 0; j < faceLocalConnect.size(); ++j) {
				freeVertices.insert(cell.getVertex(faceLocalConnect[j]));
			}
		}
	}

        return freeVertices.size();
}

/*!
	Count orphan vertices in the patch.

	An orphan vertex is a vertex not linked by any cells.

	\result The number of orphan vertices.
*/
long Patch::countOrphanVertices() const
{
	std::unordered_set<long> usedVertices;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int i = 0; i < nCellVertices; ++i) {
			usedVertices.insert(cell.getVertex(i));
		}
	}

	return (m_nVertices - usedVertices.size());
}

/*!
	Find orphan vertices in the patch.

	An orphan vertex is a vertex not linked by any cells.

	\result The list of orphan vertice.
*/
std::vector<long> Patch::findOrphanVertices()
{
	// Add all the vertices to the list
	std::unordered_set<long> vertexSet;
	for (const Vertex &vertex : m_vertices) {
		vertexSet.insert(vertex.get_id());
	}

	// Remove used vertices
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int i = 0; i < nCellVertices; ++i) {
			vertexSet.erase(cell.getVertex(i));
		}
	}

	// Build a list
	std::vector<long> vertexList;
	vertexList.reserve(vertexSet.size());
	for (const long &id : vertexSet) {
		vertexList.emplace_back();
		long &lastId = vertexList.back();
		lastId = id;
	}

	return vertexList;
}

/*!
	Remove orphan vertices
*/
bool Patch::deleteOrphanVertices()
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> list = findOrphanVertices();
	deleteVertices(list);

	return true;
}

/*!
	Find and collapse coincident vertices. Cell connectivity is
	automatically updated.

	\param[in] nBins (default = 128) is the number of bins used by
	bin-sorting algorithm to sort tasselation vertices
	\result The list of the of the collapsed vertices.
*/
std::vector<long> Patch::collapseCoincidentVertices(int nBins)
{
	std::vector<long> collapsedVertices;
	if (!isExpert()) {
		return collapsedVertices;
	}

	std::vector<std::vector<std::array<long, 2> > > bins;

	// ====================================================================== //
	// INITIALIZE LOCAL VARIABLES                                             //
	// ====================================================================== //

	// Random number generator
	srand(1223145611);

	// Resize variables
	bins.resize(nBins * nBins * nBins);

	// ====================================================================== //
	// SORT VERTICES ON BINS                                                  //
	// ====================================================================== //

	// Sort vertices
	std::unordered_map<long, long> bin_index = binSortVertex(nBins);

	// Sort cells
	std::array<long, 2> binEntry;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; ++j) {
			binEntry[0] = cell.get_id();
			binEntry[1] = j;

			long vertexId = cell.getVertex(j);
			long binId = bin_index[vertexId];
			bins[binId].push_back(binEntry);
		}
	}

	// Free memory
	bin_index.clear();

	// ====================================================================== //
	// COLLAPSE DOUBLE VERTICES                                               //
	// ====================================================================== //
	long collapsedVertexId;
	std::vector<bool> flag(m_nVertices, false);
	for (auto &bin : bins) {
		int nBinCells = bin.size();
		if (nBinCells > 0) {
			// Randomize vertex insertion
			std::vector<int> list;
			utils::extractWithoutReplacement(nBinCells, nBinCells - 1, list);

			// Vertex insertion
			KdTree<3, Vertex, long> kd(nBinCells);
			for (int j = 0; j < nBinCells; ++j) {
				long cellId   = bin[list[j]][0];
				long k        = bin[list[j]][1];
				long vertexId = m_cells[cellId].getVertex(k);

				if (kd.exist(&m_vertices[vertexId], collapsedVertexId) >= 0) {
					m_cells[cellId].setVertex(k, collapsedVertexId);
					if (!flag[vertexId]) {
						flag[vertexId] = true;
						collapsedVertices.push_back(vertexId);
					}
				} else {
					flag[vertexId] = true;
					kd.insert(&m_vertices[vertexId], vertexId);
				}
			}
		}
	}

	return collapsedVertices;
}

/*!
	Remove coincident vertices from the patch.

	\param[in] nBins (default = 128) is the number of bins used by bin
	sorting algotrithm to sort patch vertices.
*/
bool Patch::deleteCoincidentVertex(int nBins)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long> verticesToDelete = collapseCoincidentVertices(nBins);
	deleteVertices(verticesToDelete);

	return true;
}

/*!
	Gets the coordinates of the specified vertex.

	\param id is the id of the vertex
	\result The coordinates of the specified vertex.
*/
const std::array<double, 3> & Patch::getVertexCoords(const long &id) const
{
	return getVertex(id).getCoords();
}

/*!
	Gets the number of cells in the patch.

	\return The number of cells in the patch
*/
long Patch::getCellCount() const
{
	return m_cells.size();
}

/*!
	Gets the number of internal cells in the patch.

	\return The number of internal cells in the patch
*/
long Patch::getInternalCount() const
{
	return m_nInternals;
}

/*!
	Gets the number of ghost cells in the patch.

	\return The number of ghost cells in the patch
*/
long Patch::getGhostCount() const
{
	return m_nGhosts;
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
	Gets a reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A reference to the cell with the specified id.
*/
Cell & Patch::getCell(const long &id)
{
	return m_cells[id];
}

/*!
	Gets a constant reference to the cell with the specified id.

	\param id is the id of the requested cell
	\return A constant reference to the cell with the specified id.
*/
const Cell & Patch::getCell(const long &id) const
{
	return m_cells[id];
}

/*!
	Gets a reference to the last internal cell.

	\return A reference to the last internal cell.
*/
Cell & Patch::getLastInternal()
{
	return m_cells[m_last_internal_id];
}

/*!
	Gets a constant reference to the last internal cell.

	\return A constant reference to the last internal cell.
*/
const Cell & Patch::getLastInternal() const
{
	return m_cells[m_last_internal_id];
}

/*!
	Gets a reference to the first ghost cell.

	\return A reference to the first ghost cell.
*/
Cell & Patch::getFirstGhost()
{
	return m_cells[m_first_ghost_id];
}

/*!
	Gets a constant reference to the first ghost cell.

	\return A constant reference to the first ghost cell.
*/
const Cell & Patch::getFirstGhost() const
{
	return m_cells[m_first_ghost_id];
}

/*!
	Returns an iterator pointing to the specified cell.

	\result An iterator to the specified cell.
*/
Patch::CellIterator Patch::getCellIterator(const long &id)
{
	return CellIterator(m_cells.raw_begin() + m_cells.raw_index(id));
}

/*!
	Returns iterator pointing to the first cell.

	\result An iterator to the first cell.
*/
Patch::CellIterator Patch::cellBegin()
{
	return m_cells.begin();
}

/*!
	Returns iterator pointing to last cell.

	\result An iterator to the last cell.
*/
Patch::CellIterator Patch::cellEnd()
{
	return m_cells.end();
}

/*!
	Returns iterator pointing to the first internal cell.

	\result An iterator to the first internal cell.
*/
Patch::CellIterator Patch::internalBegin()
{
	return m_cells.begin();
}

/*!
	Returns iterator pointing to the end of the list of internal cells.

	\result An iterator to the end of the list of internal cells.
*/
Patch::CellIterator Patch::internalEnd()
{
	return ++CellIterator(m_cells.raw_begin() + m_cells.raw_index(m_last_internal_id));
}

/*!
    Returns iterator to the first ghost cells within the cell list.

    \result An iterator to the first ghost cell.
*/
Patch::CellIterator Patch::ghostBegin()
{
    return CellIterator(m_cells.raw_begin() + m_cells.raw_index(m_first_ghost_id));
}

/*!
	Returns iterator to the end of the list of ghost cells.

	\result An iterator to the end of the list of ghost cell.
*/
Patch::CellIterator Patch::ghostEnd()
{
	return m_cells.end();
}

/*!
	Generates a new unique id for the cells.

	\result A new unique id for the cells.
*/
long Patch::generateCellId()
{
	if (!isExpert()) {
		return Element::NULL_ID;
	}

	return m_cellIdGenerator.generateId();
}

/*!
	Creates a new cell with the specified id.

	\param type is the type of the cell
	\param id is the id of the new cell
	\param interior is true if the cell is an interior cell, false otherwise
	\return An iterator pointing to the newly created cell.
*/
Patch::CellIterator Patch::createCell(ElementInfo::Type type, bool interior, long id)
{
	if (id < 0) {
		id = generateCellId();
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(type);
	if (cellTypeInfo.dimension > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator;
	if (interior) {
		// Create an internal cell
		//
		// If there are ghosts cells, the internal cell should be inserted
		// before the first ghost cell.
		if (m_first_ghost_id < 0) {
			iterator = m_cells.reclaim(id);
		} else {
			iterator = m_cells.reclaim_before(m_first_ghost_id, id);
		}
		m_nInternals++;

		// Update the id of the last internal cell
		if (m_last_internal_id < 0) {
			m_last_internal_id = id;
		} else if (m_cells.raw_index(m_last_internal_id) < m_cells.raw_index(id)) {
			m_last_internal_id = id;
		}
	} else {
		// Create a ghost cell
		//
		// If there are internal cells, the ghost cell should be inserted
		// after the last internal cell.
		if (m_last_internal_id < 0) {
			iterator = m_cells.reclaim(id);
		} else {
			iterator = m_cells.reclaim_after(m_last_internal_id, id);
		}
		m_nGhosts++;

		// Update the id of the first ghost cell
		if (m_first_ghost_id < 0) {
			m_first_ghost_id = id;
		} else if (m_cells.raw_index(m_first_ghost_id) > m_cells.raw_index(id)) {
			m_first_ghost_id = id;
		}
	}

	return iterator;
}

/*!
	Adds a new cell with the specified id.

	\param type is the type of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(ElementInfo::Type type, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	return createCell(type, true, id);
}

/*!
	Adds a new cell with the specified id and type.

	\param type is the type of the cell
	\param interior is true if the cell is the interior of the patch,
	false otherwise
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(ElementInfo::Type type, bool interior, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = createCell(type, interior, id);
	Cell &cell = (*iterator);
	cell.initialize(type, interior);

	return iterator;
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connect is the connectivity of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(ElementInfo::Type type, bool interior,
                                   std::unique_ptr<long[]> &connect, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = addCell(type, interior, id);
	Cell &cell = (*iterator);
	cell.setConnect(std::move(connect));

	return iterator;
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param interior defines if the cell is in the interior of the patch
	or if it's a ghost cell
	\param connect is the connectivity of the cell
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(ElementInfo::Type type, bool interior,
								   const std::vector<long> &connect, const long &id)
{
	if (!isExpert()) {
		return cellEnd();
	}


	// Add the cell
	CellIterator iterator = addCell(type, interior, id);

	// Set the connectivity
	Cell &cell = (*iterator);
	int nCellVertices = cell.getVertexCount();
	std::unique_ptr<long[]> cellConnect = std::unique_ptr<long[]>(new long[nCellVertices]);
	std::copy(connect.data(), connect.data() + nCellVertices, cellConnect.get());
	cell.setConnect(std::move(cellConnect));

	return iterator;
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(const Cell &source, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	CellIterator iterator = createCell(source.getType(), source.isInterior(), id);
	Cell &cell = (*iterator);
	id = cell.get_id();
	cell = source;
	cell.set_id(id);

	return iterator;
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added cell.
*/
Patch::CellIterator Patch::addCell(Cell &&source, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (id < 0) {
		id = source.get_id();
	}

	CellIterator iterator = createCell(source.getType(), source.isInterior(), id);
	Cell &cell = (*iterator);
	id = cell.get_id();
	cell = std::move(source);
	cell.set_id(id);

	return iterator;
}

/*!
	Deletes a cell.

	\param id is the id of the cell
	\param updateNeighs if true the neighbour data will be updated after
	removing the cell
	\param delayed is true a delayed delete will be performed
*/
bool Patch::deleteCell(const long &id, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	// Update neighbours
	if (updateNeighs) {
		const Cell &cell = m_cells[id];
		int nCellFaces = m_cells[id].getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			int nFaceAdjacencies = cell.getAdjacencyCount(i);
			for (int k = 0; k < nFaceAdjacencies; ++k) {


				// Update adjacency of the neighbours
				long neighId = cell.getAdjacency(i,k);
                                if (neighId >= 0) {
                                    Cell &neigh = m_cells[neighId];

/*TODO AGGIORNARE */            int neighFace;
/*TODO AGGIORNARE */            if (interface.getOwner() == neighId) {
/*TODO AGGIORNARE */            	neighFace = interface.getOwnerFace();
/*TODO AGGIORNARE */            } else {
/*TODO AGGIORNARE */            	neighFace = interface.getNeighFace();
/*TODO AGGIORNARE */            }

/*TODO AGGIORNARE */            int adjacenyId = 0;
/*TODO AGGIORNARE */            while (neigh.getAdjacency(neighId, adjacenyId) != id) {
/*TODO AGGIORNARE */                    ++adjacenyId;
/*TODO AGGIORNARE */            }
                                    FindNeigh(neighId, id, neighFace, adjacencyId);
                                    if (neighFace >= 0) neigh.deleteAdjacency(neighFace, adjacenyId);

                                }

				// Update interface
                                long interfaceId = cell.getInterface(i,k);
                                if (interfaceId >= 0) {
                                    Interface &interface = m_interfaces[interfaceId];
                                    if (interface.getOwner() == id) {
                                            interface.unsetOwner();
                                    } else {
                                            interface.unsetNeigh();
                                    }
                                }
			}
		}
	}

	// Delete cell
	bool isInternal = m_cells.at(id).isInterior();
	m_cells.erase(id, delayed);
	m_cellIdGenerator.trashId(id);
	if (isInternal) {
		m_nInternals--;
		m_last_internal_id = m_cells.get_size_marker(m_nInternals - 1, Element::NULL_ID);
	} else {
		m_nGhosts--;
		m_first_ghost_id = m_cells.get_size_marker(m_nInternals, Element::NULL_ID);
	}

	return true;
}

/*!
	Deletes a list of cells.

	\param ids are the ids of the cells to be deleted
	\param updateNeighs if true the neighbour data will be updated after
	removing the cell
	\param delayed is true a delayed delete will be performed
 */
bool Patch::deleteCells(const std::vector<long> &ids, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteCell(*i, updateNeighs, delayed);
	}

	return true;
}

/*!
	Sets the internal flag of a cell.

	\param[in] id is the index of the cell
	\param[in] isInternal is the internal flag that will be set
*/
bool Patch::setCellInternal(const long &id, bool isInternal)
{
	if (!isExpert()) {
		return false;
	}

	if (m_cells[id].isInterior() == isInternal) {
		return true;
	} else if (isInternal) {
		moveGhost2Internal(id);
	} else {
		moveInternal2Ghost(id);
	}

	return true;
}

/*!
	Converts an internal cell to a ghost cell.

	\param[in] id is the index of the cell
*/
Patch::CellIterator Patch::moveInternal2Ghost(const long &id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// If we are moving the last internal cell we can just update the
	// last internal and first ghost markers. Otherwise the cell needs
	// to be moved.
	PiercedVector<Cell>::iterator iterator;
	if (id == m_last_internal_id) {
		// Cell iterator
		iterator = CellIterator(m_cells.raw_begin() + m_cells.raw_index(id));

		// Update markers
		m_first_ghost_id   = id;
		m_last_internal_id = m_cells.get_size_marker(m_nInternals - 1, Element::NULL_ID);
	} else {
		// Move the cell
		iterator = m_cells.move_after(m_last_internal_id, id);

		// Update markers
		if (m_first_ghost_id < 0) {
			m_first_ghost_id = id;
		} else if (m_cells.raw_index(m_first_ghost_id) > m_cells.raw_index(id)) {
			m_first_ghost_id = id;
		}
	}
	iterator->setInterior(false);

	// Update counters
	--m_nInternals;
	++m_nGhosts;

	// Return the iterator to the new position
	return(iterator);
}

/*!
	Converts a ghost cell to an internal cell.

	\param[in] id is the index of the cell
*/
Patch::CellIterator Patch::moveGhost2Internal(const long &id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// If we are moving the first ghost cell we can just update the
	// last internal and first ghost markers. Otherwise the cell needs
	// to be moved.
	PiercedVector<Cell>::iterator iterator;
	if (id == m_first_ghost_id) {
		// Cell iterator
		iterator = CellIterator(m_cells.raw_begin() + m_cells.raw_index(id));

		// Update markers
		m_last_internal_id = id;
		m_first_ghost_id   = m_cells.get_size_marker(m_nInternals, Element::NULL_ID);
	} else {
		// Move cell
		iterator = m_cells.move_before(m_first_ghost_id, id);

		// Update the id of the last internal cell
		if (m_last_internal_id < 0) {
			m_last_internal_id = id;
		} else if (m_cells.raw_index(m_last_internal_id) < m_cells.raw_index(id)) {
			m_last_internal_id = id;
		}
	}
	iterator->setInterior(true);

	// Update counters
	++m_nInternals;
	--m_nGhosts;

	// Return the iterator to the new position
	return(iterator);
}

/*!
	Counts free cells within the patch.

	A cell is free if contains at least one free face.

	\return The number of free cells.
*/
long Patch::countFreeCells() const
{
	double nFreeCells = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nFreeCells;
				break;
			}
		}
	}

	return nFreeCells;
}

/*!
	Counts orphan cells within the patch.

	A cell is orphan if not adjacent to any cell in the patch (neither
	along an edge, nor at vertex)

	\return The number of orphan cells.
*/
long Patch::countOrphanCells() const
{
	// Compute vertex valence
	std::unordered_map<long, short> vertexValence;
	for (const Cell &cell : m_cells) {
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; j++) {
			vertexValence[cell.getVertex(j)] += 1;
		}
	}

	// Loop over cells
	long nOrphanCells = 0;
	for (const Cell &cell : m_cells) {
		long isIsolated = true;
		int nCellVertices = cell.getVertexCount();
		for (int j = 0; j < nCellVertices; j++) {
			long vertexId = cell.getVertex(j);
			if (vertexValence[vertexId] > 1) {
				isIsolated = false;
				break;
			}
		}

		if (isIsolated) {
			++nOrphanCells;
		}
        }

	return nOrphanCells;
}

/*!
	Extracts the neighbours of all the faces of the specified cell.

	\param id is the id of the cell
	\result The neighbours of all the faces of the specified cell.
*/
std::vector<long> Patch::findCellFaceNeighs(const long &id) const
{
	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getFaceCount(); ++i) {
		std::vector<long> faceNeighs = findCellFaceNeighs(id, i);
		for (auto &neighId : faceNeighs) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts all the neighbours of the specified cell

	\param id is the id of the cell
	\result All the neighbours of the specified cell.
*/
std::vector<long> Patch::findCellNeighs(const long &id) const
{
	return findCellVertexNeighs(id);
}

/*!
	Extracts all the neighbours of the specified cell for the given
	codimension.

	\param id is the id of the cell
	\param codimension the codimension for which the neighbours
	are requested. For a three-dimensional cell a codimension
	equal 1 will extract the face neighbours, a codimension equal
	2 will extract the edge negihbours and a codimension equal
	3 will extract the vertex neighbours. For a two-dimensional
	cell a codimension qual 1 will extract the face neighbours,
	and a codimension equal 2 will extract the vertex neighbours.
	\param complete controls if the list of neighbours should contain
	only the neighbours for the specified codimension, or should contain
	also the neighbours for lower codimensions.
	\result The neighbours for the specified codimension.
*/
std::vector<long> Patch::findCellNeighs(const long &id, int codimension, bool complete) const
{
	assert(codimension >= 1 && codimension <= getDimension());

	if (codimension == 1) {
		return findCellFaceNeighs(id);
	} else if (codimension == getDimension()) {
		return findCellVertexNeighs(id, complete);
	} else if (codimension == 2) {
		return findCellEdgeNeighs(id, complete);
	} else {
		return std::vector<long>();
	}
}

/*!
	Extracts the neighbours of the specified cell for the given face.

	\param id is the id of the cell
	\param face is a face of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given face.
*/
std::vector<long> Patch::findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList) const
{
	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getAdjacencyCount(face); ++i) {
		long neighId = cell.getAdjacency(face, i);
		if (neighId < 0) {
			continue;
		}

		if (std::find(blackList.begin(), blackList.end(), neighId) == blackList.end()) {
			utils::addToOrderedVector<long>(neighId, neighs);
		}
	}

	return neighs;
}

/*!
        Stores the local index of the face shared by cell_idx and neigh_idx
        into face_loc_idx.
        If cell cell_idx and neigh_idx do not share any face, -1 is stored in
        face_loc_idx.

        \param[in] cell_idx cell index
        \param[in] neigh_idx neighbour index
        \param[in,out] face_loc_idx on output stores the local index (on cell
        cell_idx) shared by cell_idx and neigh_idx. If cells cell_idx and neigh_idx
        do not share any face, -1 is stored in face_loc_idx.
        \param[in,out] intf_loc_idx on output stores the index of the adjacency (on face
        face_loc_idx of cell cell_idx). If cells cell_idx and neigh_idx do not share
        any face, -1 is stored into intf_loc_idx.
*/
void Patch::findFaceNeighCell(const long &cell_idx, const long &neigh_idx, int &face_loc_idx, int &intf_loc_idx)
{

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //
    bool                        loop_continue = true;
    int                         n_faces, n_adj;
    int                         j, k;
    Cell                       &cell_ = cells[cell_idx];

    // ====================================================================== //
    // LOOP OVER ADJACENCIES                                                  //
    // ====================================================================== //
    n_faces = cell_->getFaceCount();
    j = 0;
    while ( loop_continue && (j < n_faces) ) {
        n_adj = cell_->getAdjacencyCount(j);
        k = 0;
        while ( loop_continue && (k < n_adj) ) {
            loop_continue = ( cell_->getAdjacency( j, k ) != neigh_idx );
            ++k;
        } //next k
        ++j;
    } //next j

    if ( loop_continue) { face_loc_idx = intf_loc_idx = -1; }
    else                { face_loc_idx = --j; intf_loc_idx = --k; }

    return;
}

/*!
	Extracts the neighbours of all the edges of the specified cell.

	This function can be only used with three-dimensional cells.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified edge, or should
	contain also neighbours that share an entire face
	\result The neighbours of all the edges of the specified cell.
*/
std::vector<long> Patch::findCellEdgeNeighs(const long &id, bool complete) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return std::vector<long>();
	}

	std::vector<long> blackList;
	if (!complete) {
		blackList = findCellFaceNeighs(id);
	}

	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getEdgeCount(); ++i) {
		for (auto &neigh : findCellEdgeNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
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
std::vector<long> Patch::findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList) const
{
	assert(isThreeDimensional());
	if (!isThreeDimensional()) {
		return std::vector<long>();
	}

	const Cell &cell = getCell(id);
	std::vector<int> vertices = cell.getEdgeLocalConnect(edge);

	return findCellVertexNeighs(id, vertices, blackList);
}

/*!
	Extracts the neighbours of all the vertices of the specified cell.

	\param id is the id of the cell
	\param complete controls if the list of neighbours should contain
	only the neighbours that share just the specified vertex, or should
	contain also neighbours that share an entire face or an entire edge
	\result The neighbours of all the vertices of the specified cell.
*/
std::vector<long> Patch::findCellVertexNeighs(const long &id, bool complete) const
{
	std::vector<long> blackList;
	if (!complete) {
		if (isThreeDimensional()) {
			blackList = findCellEdgeNeighs(id);
		} else {
			blackList = findCellFaceNeighs(id);
		}
	}

	std::vector<long> neighs;
	const Cell &cell = getCell(id);
	for (int i = 0; i < cell.getVertexCount(); ++i) {
		for (auto &neigh : findCellVertexNeighs(id, i, blackList)) {
			utils::addToOrderedVector<long>(neigh, neighs);
		}
	}

	return neighs;
}

/*!
	Extracts the neighbours of the specified cell for the given vertex.

	Cells that has only a vertex in common are considered neighbours only
	if there are other cells "connecting" them.

	                  .-----.                   .-----.
	                  |     |                   |     |
	                V | A1  |                 V | A2  |
	            .-----+-----.             .-----+-----.
	            |     |                   |     |     |
	            | B1  |                   | B2  | C2  |
	            .-----.                   .-----.-----.

	For example, A1 and B1 are not neighbours (although they share the
	vertex V), whereas A2 and B2 are neighbours.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertex.
*/
std::vector<long> Patch::findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList) const
{
	std::vector<int> vertexList(1);
	vertexList[0] = vertex;

	return findCellVertexNeighs(id, vertexList, blackList);
}

/*!
	Extracts the neighbours of the specified cell for the given vertices.

	Cells that has only a vertex in common are considered neighbours only
	if there are other cells "connecting" them.

	                  .-----.                   .-----.
	                  |     |                   |     |
	                V | A1  |                 V | A2  |
	            .-----+-----.             .-----+-----.
	            |     |                   |     |     |
	            | B1  |                   | B2  | C2  |
	            .-----.                   .-----.-----.

	For example, A1 and B1 are not neighbours (although they share the
	vertex V), whereas A2 and B2 are neighbours.

	\param id is the id of the cell
	\param vertices is the list of vertices of the cell
	\param blackList is a list of cells that are excluded from the search
	\result The neighbours of the specified cell for the given vertices.
*/
std::vector<long> Patch::findCellVertexNeighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList) const
{
	const Cell &cell = getCell(id);

	int nVerticesToFound = vertices.size();
	std::vector<long> verticesToFind(nVerticesToFound);
	for (int k = 0; k < nVerticesToFound; ++k) {
		verticesToFind[k] = cell.getVertex(vertices[k]);
	}

	std::vector<long> neighs;
	std::unordered_set<long> scanQueue;
	std::unordered_set<long> alreadyScan;
	scanQueue.insert(cell.get_id());
	while (!scanQueue.empty()) {
		// Pop a cell to process
		long scanId = *(scanQueue.begin());
		const Cell &scanCell = getCell(scanId);

		scanQueue.erase(scanId);
		alreadyScan.insert(scanId);

		// Info on the cell
		const ElementInfo &cellTypeInfo = scanCell.getInfo();
		const std::vector<std::vector<int>> &cellLocalFaceConnect = cellTypeInfo.faceConnect;
		const long *scanCellConnect = scanCell.getConnect();

		// Find the faces that share the vertices
		std::vector<long> faceList;
		for (int i = 0; i < scanCell.getFaceCount(); ++i) {
			// Info on the face
			ElementInfo::Type faceType = scanCell.getFaceType(i);
			const ElementInfo &faceTypeInfo = ElementInfo::getElementInfo(faceType);
			const std::vector<int> &faceLocalConnect = cellLocalFaceConnect[i];

			// Check if the face shares all the vertices
			int nCommonVertices = 0;
			for (int k = 0; k < faceTypeInfo.nVertices; ++k) {
				long faceVertexId = scanCellConnect[faceLocalConnect[k]];
				for (int n = 0; n < nVerticesToFound; ++n) {
					if (faceVertexId == verticesToFind[n]) {
						nCommonVertices++;
						break;
					}
				}
			}

			if (nCommonVertices == nVerticesToFound) {
				faceList.push_back(i);
			}
		}

		// If there are no faces that share the vertices go to the next face
		if (faceList.empty()) {
			continue;
		}

		// Add the current cell to the neighoburs
		if (scanId != id && std::find(blackList.begin(), blackList.end(), scanId) == blackList.end()) {
			utils::addToOrderedVector<long>(scanId, neighs);
		}

		// Add the neighbours of the faces to the scan list
		for (const long &face : faceList) {
			int nFaceNeighs = scanCell.getAdjacencyCount(face);
			for (int k = 0; k < nFaceNeighs; ++k) {
				long neighId = scanCell.getAdjacency(face, k);
				if (neighId >= 0 && alreadyScan.count(neighId) == 0) {
					scanQueue.insert(neighId);
				}
			}
		}
	}

	return neighs;
}

/*!
	Finds the one-ring of the specified vertex of the cell.

	\param id is the id of the cell
	\param vertex is a vertex of the cell
	\result The one-ring of the specified vertex of the cell.
*/
std::vector<long> Patch::findCellVertexOneRing(const long &id, const int &vertex) const
{
	std::vector<long> oneRing = findCellVertexNeighs(id, vertex);
	utils::addToOrderedVector<long>(id, oneRing);

	return oneRing;
}

/*!
	Gets the number of interfaces in the patch.

	\return The number of interfaces in the patch
*/
long Patch::getInterfaceCount() const
{
	return m_nInterfaces;
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
	Gets a reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A reference to the interface with the specified id.
*/
Interface & Patch::getInterface(const long &id)
{
	return m_interfaces[id];
}

/*!
	Gets a constant reference to the interface with the specified id.

	\param id is the id of the requested interface
	\return A constant reference to the interface with the specified id.
*/
const Interface & Patch::getInterface(const long &id) const
{
	return m_interfaces[id];
}

/*!
	Returns an iterator pointing to the specified interface.

	\result An iterator to the specified interface.
*/
Patch::InterfaceIterator Patch::getInterfaceIterator(const long &id)
{
	return InterfaceIterator(m_interfaces.raw_begin() + m_interfaces.raw_index(id));
}

/*!
	Returns iterator pointing to the first interface.

	\result An iterator to the first interface.
*/
Patch::InterfaceIterator Patch::interfaceBegin()
{
	return m_interfaces.begin();
}

/*!
	Returns iterator pointing to last interface.

	\result An iterator to the last interface.
*/
Patch::InterfaceIterator Patch::interfaceEnd()
{
	return m_interfaces.end();
}

/*!
 * Generates a new unique id for the interfaces.
 *
 * \result A new unique id for the interfaces.
 */
long Patch::generateInterfaceId()
{
	if (!isExpert()) {
		return Element::NULL_ID;
	}

	return m_interfaceIdGenerator.generateId();
}

/*!
	Creates a new interface with the specified id.

	\param type is the type of the interface
	\param id is the id of the new interface
	\return An iterator pointing to the newly created interface.
*/
Patch::InterfaceIterator Patch::createInterface(ElementInfo::Type type, long id)
{
	if (id < 0) {
		id = generateInterfaceId();
	}

	const ElementInfo &cellTypeInfo = ElementInfo::getElementInfo(type);
	if (cellTypeInfo.dimension > (getDimension() - 1)) {
		return interfaceEnd();
	}

	PiercedVector<Interface>::iterator iterator = m_interfaces.reclaim(id);
	m_nInterfaces++;

	return iterator;
}

/*!
	Adds a new interface with the specified id.

	\param type is the type of the interface
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
Patch::InterfaceIterator Patch::addInterface(ElementInfo::Type type, const long &id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = createInterface(type, id);
	Interface &interface = (*iterator);
	interface.initialize(type);

	return iterator;
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id of the new cell. If a negative id value is
	specified, ad new unique id will be generated
	\return An iterator pointing to the added interface.
*/
Patch::InterfaceIterator Patch::addInterface(const Interface &source, long id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	InterfaceIterator iterator = createInterface(source.getType(), id);
	Interface &interface = (*iterator);
	id = interface.get_id();
	interface = source;
	interface.set_id(id);

	return iterator;
}

/*!
	Adds the specified interface to the patch.

	\param source is the interface that will be added
	\param id is the id that will be assigned to the interface
	\return An iterator pointing to the added interface.

*/
Patch::InterfaceIterator Patch::addInterface(Interface &&source, long id)
{
	if (!isExpert()) {
		return interfaceEnd();
	}

	if (id < 0) {
		id = source.get_id();
	}

	InterfaceIterator iterator = createInterface(source.getType(), id);
	Interface &interface = (*iterator);
	id = interface.get_id();
	interface = std::move(source);
	interface.set_id(id);

	return iterator;
}

/*!
	Deletes an interface.

	\param id is the id of the interface
	\param updateNeighs if true the neighbour data will be updated after
	removing the interface
	\param delayed is true a delayed delete will be performed
*/
bool Patch::deleteInterface(const long &id, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	// Update neighbours
	if (updateNeighs) {
		Interface &interface = m_interfaces[id];

		// Update owner
		long ownerId = interface.getOwner();
		Cell &owner = m_cells[ownerId];
		int ownerFace = interface.getOwnerFace();

		int ownerInterfaceId = 0;
		while (owner.getInterface(ownerFace, ownerInterfaceId) != id) {
			++ownerInterfaceId;
		}
		owner.deleteInterface(ownerFace, ownerInterfaceId);

		// Update neighbour
		long neighId = interface.getNeigh();
		if (neighId >= 0) {
			Cell &neigh = m_cells[neighId];
			int neighFace = interface.getNeighFace();

			int neighInterfaceId = 0;
			while (neigh.getInterface(neighFace, neighInterfaceId) != id) {
				++neighInterfaceId;
			}
			neigh.deleteInterface(neighFace, neighInterfaceId);
		}
	}

	// Delete interface
	m_interfaces.erase(id, delayed);
	m_interfaceIdGenerator.trashId(id);
	m_nInterfaces--;

	return true;
}

/*!
	Deletes a list of interfaces.

	\param ids are the ids of the interfaces to be deleted
	\param updateNeighs if true the neighbour data will be updated after
	removing the interface
	\param delayed is true a delayed delete will be performed
*/
bool Patch::deleteInterfaces(const std::vector<long> &ids, bool updateNeighs, bool delayed)
{
	if (!isExpert()) {
		return false;
	}

	std::vector<long>::const_iterator end = ids.cend();
	for (std::vector<long>::const_iterator i = ids.cbegin(); i != end; ++i) {
		deleteInterface(*i, updateNeighs, delayed);
	}

	return true;
}

/*!
	Counts free interfaces within the patch.

	An interface is free if belongs to just one cell.

	\result The number of free interfaces.
*/
long Patch::countFreeInterfaces() const
{
	long nFreeInterfaces = 0;
	for (const Interface &interface : m_interfaces) {
		if (interface.getNeigh() < 0) {
			++nFreeInterfaces;
		}
        }

	return nFreeInterfaces;
}

/*!
	Counts orphan interfaces within the patch.

	An interface is orphan if not linked to any cell in the patch.

	\return The number of orphan interfaces.
*/
long Patch::countOrphanInterfaces() const
{
	long nOrphanInterfaces = 0;
	for (const Interface &interface : m_interfaces) {
		if (interface.getOwner() < 0 && interface.getNeigh() < 0) {
			++nOrphanInterfaces;
		}
        }

	return nOrphanInterfaces;
}

/*!
	Count faces within the patch.

	\result The total number of faces in the patch.
*/
long Patch::countFaces() const
{
	double nFaces = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				nFaces += 1. / (cell.getAdjacencyCount(i) + 1);
			} else {
				nFaces += 1.;
			}
		}
	}

	return ((long) round(nFaces));
}

/*!
	Counts free faces within the patch.

	A face is free if a cell has no adjacent along that faces.

	\result The number of free faces.
*/
long Patch::countFreeFaces() const
{
	double nFreeFaces = 0;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (cell.isFaceBorder(i)) {
				++nFreeFaces;
			}
		}
	}

	return nFreeFaces;
}

/*!
	Sorts internal vertex storage in ascending id order.
*/
bool Patch::sortVertices()
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.sort();

	return true;
}

/*!
	Sorts internal cell storage in ascending id order.
*/
bool Patch::sortCells()
{
	if (!isExpert()) {
		return false;
	}

	m_cells.sort();

	return true;
}

/*!
	Sorts internal interface storage in ascending id order.
*/
bool Patch::sortInterfaces()
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.sort();

	return true;
}

/*!
	Sorts internal storage for cells, vertices and interfaces in
	ascending id order.
*/
bool Patch::sort()
{
	bool status = sortVertices();
	status |= sortCells();
	status |= sortInterfaces();

	return status;
}

/*!
	Requests the patch to compact the vertex data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the vertex
	data structure can still occupy more memory than it actually needs.
*/
bool Patch::squeezeVertices()
{
	if (!isExpert()) {
		return false;
	}

	m_vertices.squeeze();

	return true;
}

/*!
	Requests the patch to compact the cell data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the cell
	data structure can still occupy more memory than it actually needs.
*/
bool Patch::squeezeCells()
{
	if (!isExpert()) {
		return false;
	}

	m_cells.squeeze();

	return true;
}

/*!
	Requests the patch to compact the interface data structure and reduce
	its capacity to fit its size.

	The request is non-binding, and after the function call the interface
	data structure can still occupy more memory than it actually needs.
*/
bool Patch::squeezeInterfaces()
{
	if (!isExpert()) {
		return false;
	}

	m_interfaces.squeeze();

	return true;
}

/*!
	Requests the patch to compact the data structures and reduce its
	capacity to fit its size.

	The request is non-binding, and after the function call the patch
	can still occupy more memory than it actually needs.
*/
bool Patch::squeeze()
{
	bool status = squeezeVertices();
	status |= squeezeCells();
	status |= squeezeInterfaces();

	return status;
}

/*!
	Evaluates the centroid of the specified cell.

	\param id is the id of the cell
	\result The centroid of the specified cell.
*/
std::array<double, 3> Patch::evalCellCentroid(const long &id)
{
	Cell &cell = getCell(id);

	return evalElementCentroid(cell);
}

/*!
	Evaluates the centroid of the specified interface.

	\param id is the id of the interface
	\result The centroid of the specified interface.
*/
std::array<double, 3> Patch::evalInterfaceCentroid(const long &id)
{
	Interface &interface = getInterface(id);

	return evalElementCentroid(interface);
}

/*!
	Evaluates the centroid of the specified element.

	\param element is the element
	\result The centroid of the specified element.
*/
std::array<double, 3> Patch::evalElementCentroid(const Element &element)
{
	const int nDimensions = 3;

	const long *elementConnect = element.getConnect();
	const ElementInfo &elementInfo = element.getInfo();

	std::array<double, nDimensions> centroid = {{0., 0., 0.}};
	for (int i = 0; i < elementInfo.nVertices; ++i) {
		Vertex &vertex = getVertex(elementConnect[i]);
		const std::array<double, nDimensions> &vertexCoords = vertex.getCoords();
		for (int k = 0; k < nDimensions; ++k) {
			centroid[k] += vertexCoords[k];
		}
	}

	for (int k = 0; k < nDimensions; ++k) {
		centroid[k] /= elementInfo.nVertices;
	}

	return centroid;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool Patch::isPointInside(const double &x, const double &y, const double &z)
{
	return isPointInside({{x, y, z}});
}

/*!
	Locates the cell the contains the point.

	If the point is not inside the patch, the function returns the id of the
	null element.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns the id of the cell the contains the point. If the point
	is not inside the patch, the function returns the id of the null element.
*/
long Patch::locatePoint(const double &x, const double &y, const double &z)
{
	return locatePoint({{x, y, z}});
}

/*!
 * Check whether the i-th face on cell "cell_1" is the same as the j-th face
 * on cell "cell_2".
 * 
 * \param[in] cell_1 global ID of the 1st cell
 * \param[in] i      local index of face to be checked on cell_1
 * \param[in] cell_1 global ID of the 2nd cell
 * \param[in] i      local index of face to be checked on cell_2
 * 
 * \result returns true if face (cell_1, i) and face (cell_2, j) are the same.
*/
bool Patch::isSameFace(
    const long                  &cell_1,
    const int                   &i,
    const long                  &cell_2,
    const int                   &j
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                            check = false;
std::vector<int>                face_loc_connect_A, face_loc_connect_B;
Cell                            *cell_1_ = &m_cells[cell_1], *cell_2_ = &m_cells[cell_2];

// Counters
int                             k;

// ========================================================================== //
// CHECK FOR COINCIDENT FACES                                                 //
// ========================================================================== //
face_loc_connect_A = cell_1_->getFaceLocalConnect(i);
face_loc_connect_B = cell_2_->getFaceLocalConnect(j);
if (face_loc_connect_A.size() == face_loc_connect_B.size()) {
    for (k = 0; k < face_loc_connect_A.size(); ++k) {
        face_loc_connect_A[k] = cell_1_->getVertex(face_loc_connect_A[k]);
    } //next i
    for (k = 0; k < face_loc_connect_B.size(); ++k) {
        face_loc_connect_B[k] = cell_2_->getVertex(face_loc_connect_B[k]);
    } //next i
    std::sort(face_loc_connect_A.begin(), face_loc_connect_A.end());
    std::sort(face_loc_connect_B.begin(), face_loc_connect_B.end());
    check = (face_loc_connect_A == face_loc_connect_B);
}

return(check);
    
};

/*!
	Updates the stored patch bounding box.
*/
void Patch::updateBoundingBox()
{
	evalBoundingBox(m_minPoint, m_maxPoint);
}

/*!
	Gets the previously stored patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void Patch::getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	minPoint = m_minPoint;
	maxPoint = m_maxPoint;
}

/*!
	Evalautes patch bounding box.

	\param[out] minPoint on output stores the minimum point of the patch
	\param[out] maxPoint on output stores the maximum point of the patch
*/
void Patch::evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint)
{
	// Initialize bounding box
	for (int k = 0; k < 3; ++k) {
		minPoint[k] =   std::numeric_limits<double>::max();
		maxPoint[k] = - std::numeric_limits<double>::max();
	}

	// Compute bounding box limits
	PiercedVector<Vertex>::iterator e = m_vertices.end();
	for (PiercedVector<Vertex>::iterator i = m_vertices.begin(); i != e; ++i) {
		for (int k = 0; k < 3; ++k) {
			double value = i->getCoords()[k];

			minPoint[k] = std::min(minPoint[k], value);
			maxPoint[k] = std::max(maxPoint[k], value);
		}
	}
}

/*!
	Sort patch vertices on regular bins.

	\param[in] nBins (default = 128) is the number of bins (on each space
	direction)
	\result Returns the bin index associated to each vertex.
*/
std::unordered_map<long, long> Patch::binSortVertex(int nBins)
{
	// ====================================================================== //
	// VARIABLES DECLARATION                                                  //
	// ====================================================================== //

	// Local variables
	double                              dx, dy, dz;

	// Counters
	long                                i, j, k;
	PiercedVector<Vertex>::iterator     V, E = m_vertices.end();

	// ====================================================================== //
	// ASSOCIATE EACH VERTEX WITH A BIN                                       //
	// ====================================================================== //

	// Bin's spacing
	dx = max(1.0e-12, m_maxPoint[0] - m_minPoint[0]) / ((double) nBins);
	dy = max(1.0e-12, m_maxPoint[1] - m_minPoint[1]) / ((double) nBins);
	dz = max(1.0e-12, m_maxPoint[2] - m_minPoint[2]) / ((double) nBins);

	// Loop over vertices
	std::unordered_map<long, long> bin_index;
	for (V = m_vertices.begin(); V != E; ++V) {
		i = std::min(nBins - 1L, long((V->getCoords()[0] - m_minPoint[0]) / dx));
		j = std::min(nBins - 1L, long((V->getCoords()[1] - m_minPoint[1]) / dy));
		k = std::min(nBins - 1L, long((V->getCoords()[2] - m_minPoint[2]) / dz));
		bin_index[V->get_id()] = nBins * nBins * k + nBins * j + i;
	}

	return bin_index;
}

/*!
	Translates the patch.

	\param[in] translation is the translation vector
*/
void Patch::translate(std::array<double, 3> translation)
{
	for (auto &vertex : m_vertices) {
		vertex.translate(translation);
	}
}

/*!
	Translates the patch.

	\param[in] sx translation along x direction
	\param[in] sy translation along y direction
	\param[in] sz translation along z direction
*/
void Patch::translate(double sx, double sy, double sz)
{
	translate({{sx, sy, sz}});
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor vector
*/
void Patch::scale(std::array<double, 3> scaling)
{
	for (auto &vertex : m_vertices) {
		vertex.scale(scaling, m_minPoint);
	}
}

/*!
	Scales the patch.

	The patch is scaled about the lower-left point of the bounding box.

	\param[in] scaling is the scaling factor
*/
void Patch::scale(double scaling)
{
	scale({{scaling, scaling, scaling}});
}

/*!
	Scales the patch.

	\param[in] sx scaling factor along x direction
	\param[in] sy scaling factor along y direction
	\param[in] sz scaling factor along z direction
*/
void Patch::scale(double sx, double sy, double sz)
{
	scale({{sx, sy, sz}});
}

/*!
	Sets the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void Patch::setTol(double tolerance)
{
	_setTol(tolerance);

	m_hasCustomTolerance = true;
}

/*!
	Internal function to set the tolerance for the geometrical checks.

	\param tolerance is the tolerance that will be used for the geometrical
	checks
*/
void Patch::_setTol(double tolerance)
{
	m_tolerance = tolerance;
}

/*!
	Gets the tolerance for the geometrical checks.

	\result The tolerance fot the geometrical checks.
*/
double Patch::getTol() const
{
	return m_tolerance;
}

/*!
	Resets the tolerance for the geometrical checks.
*/
void Patch::resetTol()
{
	_resetTol();

	m_hasCustomTolerance = false;
}

/*!
	Internal function to reset the tolerance for the geometrical checks.
*/
void Patch::_resetTol()
{
	m_tolerance = 1;
	for (int k = 0; k < 3; ++k) {
		m_tolerance = std::max(m_maxPoint[0] - m_minPoint[0], m_tolerance);
	}
	m_tolerance *= 1e-14;
}

/*!
	Checks if the tolerance for the geometrical checks has been customized
	by the user.

	\result True if the tolerance was customized by the user, false otherwise.
*/
bool Patch::isTolCustomized() const
{
	return m_hasCustomTolerance;
}

/*!
	Extracts the external envelope and appends it to the given patch.

	The external envelope is composed by all the free faces of the patch.

	\param[in,out] envelope is the patch to which the external envelope
	will be appended
*/
void Patch::extractEnvelope(Patch &envelope) const
{
	// ====================================================================== //
	// RESIZE DATA STRUCTURES                                                 //
	// ====================================================================== //
	envelope.reserveVertices(envelope.getVertexCount() + countFreeVertices());
	envelope.reserveCells(envelope.getCellCount() + countFreeFaces());

	// ====================================================================== //
	// LOOP OVER CELLS                                                        //
	// ====================================================================== //
	std::unordered_map<long, long> vertexMap;
	for (const Cell &cell : m_cells) {
		int nCellFaces = cell.getFaceCount();
		for (int i = 0; i < nCellFaces; ++i) {
			if (!cell.isFaceBorder(i)) {
				continue;
			}

			// Add face vertices to the envelope and get face
			// connectivity in the envelope
			std::vector<int> faceLocalConnect = cell.getFaceLocalConnect(i);
			int nFaceVertices = faceLocalConnect.size();

			std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);
			for (int j = 0; j < nFaceVertices; ++j) {
				long vertexId = cell.getVertex(faceLocalConnect[j]);

				// If the vertex is not yet in the envelope
				// add it.
				if (vertexMap.count(vertexId) == 0) {
					const Vertex &vertex = getVertex(vertexId);
					VertexIterator envelopeVertex = envelope.addVertex(vertex);
					vertexMap[vertexId] = envelopeVertex->get_id();
				}

				// Update face ace connectivity in the envelope
				faceEnvelopeConnect[j] = vertexMap.at(vertexId);
			}

			// Add face to envelope
			ElementInfo::Type faceType = cell.getFaceType(i);
			envelope.addCell(faceType, true, faceEnvelopeConnect);
		}
	}
}

/*!
	Display patch statistics.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void Patch::displayTopologyStats(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');

	// ====================================================================== //
	// VERTEX STATS                                                           //
	// ====================================================================== //
	out << indent<< "Vertices --------------------------------"     << endl;
	out << indent<< "  # vertices        " << getVertexCount()      << endl;
	out << indent<< "  # orphan vertices " << countOrphanVertices() << endl;
	out << indent<< "  # free vertices   " << countFreeVertices()   << endl;
        //out << indent<< "  # free vertices   " << countDoubleVertices()   << endl;

	// ====================================================================== //
	// FACE STATS                                                             //
	// ====================================================================== //
	out << indent<< "Faces -----------------------------------"     << endl;
	out << indent<< "  # faces           " << countFaces()          << endl;
	out << indent<< "  # free faces      " << countFreeFaces()      << endl;

	// ====================================================================== //
	// CELLS STATS                                                            //
	// ====================================================================== //
	out << indent<< "Cells -----------------------------------"     << endl;
	out << indent<< "  # cells           " << getCellCount()        << endl;
	out << indent<< "  # orphan cells    " << countOrphanCells()    << endl;
	out << indent<< "  # free cells      " << countFreeCells()      << endl;
        //out << indent<< "  # free vertices   " << countDoubleCells()   << endl;
}

/*!
	Display all the vertices currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void Patch::displayVertices(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Vertex &vertex : m_vertices) {
		out << indent << "vertex: " << std::endl;
		vertex.display(out, padding + 2);
	}
}

/*!
	Display all the cells currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void Patch::displayCells(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Cell &cell : m_cells) {
		out << indent << "cell: " << std::endl;
		cell.display(out, padding + 2);
	}
}

/*!
	Display all the interfaces currently stored within the patch.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void Patch::displayInterfaces(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');
	for (const Interface &interface : m_interfaces) {
		out << indent << "interface: " << std::endl;
		interface.display(out, padding + 2);
	}
}

/*!
 *  Interface method for obtaining field meta Data
 *
 *  @param[in] name is the name of the field to be written
 *  @return Returns a VTKFieldMetaData struct containing the metadata
 *  of the requested custom data.
 */
const VTKFieldMetaData Patch::getMetaData(std::string name)
{
	if (name == "Points") {
		return VTKFieldMetaData(3 * m_vertices.size(), typeid(double));
	} else if (name == "offsets") {
		return VTKFieldMetaData(m_cells.size(), typeid(int));
	} else if (name == "types") {
		return VTKFieldMetaData(m_cells.size(), typeid(VTKElementType));
	} else if (name == "connectivity") {
		long connectSize = 0;
		for (Cell &cell : m_cells) {
			connectSize += cell.getInfo().nVertices;
		}

		return VTKFieldMetaData(connectSize, typeid(long));
	}

	// This code should never be reached
	assert(false);
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] format is the format which must be used. Supported options
 *  are "ascii" or "appended". For "appended" type an unformatted binary
 *  stream must be used
 *  @param[in] name is the name of the data to be written. Either user
 *  data or patch data
 */
void Patch::flushData(std::fstream &stream, VTKFormat format, std::string name)
{
	assert(format == VTKFormat::APPENDED);

	static std::unordered_map<long, long> vertexMap;

	if (name == "Points") {
		long vertexId = 0;
		for (Vertex &vertex : m_vertices) {
			vertexMap[vertex.get_id()] = vertexId++;

			genericIO::flushBINARY(stream, vertex.getCoords());
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (Cell &cell : m_cells) {
			offset += cell.getInfo().nVertices;
			genericIO::flushBINARY(stream, offset);
		}
	} else if (name == "types") {
		for (Cell &cell : m_cells) {
			VTKElementType VTKType;
			switch (cell.getType())  {

			case ElementInfo::VERTEX:
				VTKType = VTKElementType::VERTEX;
				break;

			case ElementInfo::LINE:
				VTKType = VTKElementType::LINE;
				break;

			case ElementInfo::TRIANGLE:
				VTKType = VTKElementType::TRIANGLE;
				break;

			case ElementInfo::PIXEL:
				VTKType = VTKElementType::PIXEL;
				break;

			case ElementInfo::QUAD:
				VTKType = VTKElementType::QUAD;
				break;

			case ElementInfo::TETRA:
				VTKType = VTKElementType::TETRA;
				break;

			case ElementInfo::VOXEL:
				VTKType = VTKElementType::VOXEL;
				break;

			case ElementInfo::HEXAHEDRON:
				VTKType = VTKElementType::HEXAHEDRON;
				break;

			case ElementInfo::WEDGE:
				VTKType = VTKElementType::WEDGE;
				break;

			case ElementInfo::PYRAMID:
				VTKType = VTKElementType::PYRAMID;
				break;

			default:
				VTKType = VTKElementType::UNDEFINED;
				break;

			}

			genericIO::flushBINARY(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (Cell &cell : m_cells) {
			for (int i = 0; i < cell.getInfo().nVertices; ++i) {
				genericIO::flushBINARY(stream, vertexMap.at(cell.getVertex(i)));
			}
		}

		vertexMap.clear();
		std::unordered_map<long, long>().swap(vertexMap);
	}
}

/*!
	@}
*/

}
