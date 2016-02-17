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

#ifndef __BITPIT_PATCH_HPP__
#define __BITPIT_PATCH_HPP__

#include <cstddef>
#include <memory>
#if ENABLE_MPI==1
#	include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_IO.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class IndexGenerator {

public:
	IndexGenerator();

	long generateId();
	long getLastId();
	void trashId(const long &id);
	void reset();

private:
	long m_next;
	bool m_depleted;
	std::deque<long> m_trash;

};

class Patch : public VTKUnstructuredGrid {

public:
	typedef PiercedVector<Vertex>::iterator VertexIterator;
	typedef PiercedVector<Cell>::iterator CellIterator;
	typedef PiercedVector<Interface>::iterator InterfaceIterator;

	Patch(const int &id, const int &dimension, bool epxert);

	virtual ~Patch();

	void reset();
	void resetVertices();
	void resetCells();
	void resetInterfaces();

	bool reserveVertices(size_t nVertices);
	bool reserveCells(size_t nCells);
	bool reserveInterfaces(size_t nInterfaces);

	const std::vector<Adaption::Info> update(bool trackAdaption = true);

	void markCellForRefinement(const long &id);
	void markCellForCoarsening(const long &id);
	void enableCellBalancing(const long &id, bool enabled);

	bool isDirty() const;
	bool isExpert() const;

	int get_id() const;
	int getDimension() const;
	bool isThreeDimensional() const;

	long getVertexCount() const;
	PiercedVector<Vertex> &vertices();
	Vertex &getVertex(const long &id);
	const Vertex & getVertex(const long &id) const;
	const std::array<double, 3> & getVertexCoords(const long &id) const;
	long generateVertexId();
	VertexIterator addVertex(const long &id = Vertex::NULL_ID);
	VertexIterator addVertex(const std::array<double, 3> &coords, const long &id = Vertex::NULL_ID);
	VertexIterator addVertex(Vertex source);
	VertexIterator addVertex(Vertex &&source, long id = Vertex::NULL_ID);
	long countFreeVertices() const;
	long countOrphanVertices() const;
	std::vector<long> findOrphanVertices();
	bool deleteOrphanVertices();
	std::vector<long> collapseCoincidentVertices(int nBins = 128);
	bool deleteCoincidentVertex(int nBins = 128);

	VertexIterator getVertexIterator(const long &id);
	VertexIterator vertexBegin();
	VertexIterator vertexEnd();

	long getCellCount() const;
	long getInternalCount() const;
	long getGhostCount() const;
	PiercedVector<Cell> &cells();
	Cell &getCell(const long &id);
	const Cell &getCell(const long &id) const;
	Cell &getLastInternal();
	const Cell &getLastInternal() const;
	Cell &getFirstGhost();
	const Cell &getFirstGhost() const;
	long generateCellId();
	CellIterator addCell(const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, std::unique_ptr<long[]> &connect, const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, const std::vector<long> &connect, const long &id = Element::NULL_ID);
	CellIterator addCell(Cell source);
	CellIterator addCell(Cell &&source, long id = Element::NULL_ID);
	bool deleteCell(const long &id, bool updateNeighs = true, bool delayed = false);
	bool deleteCells(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
	bool setCellInternal(const long &id, bool isInternal);
	CellIterator moveGhost2Internal(const long &id);
	CellIterator moveInternal2Ghost(const long &id);
	virtual double evalCellVolume(const long &id) = 0;
	virtual double evalCellSize(const long &id) = 0;
	long countFreeCells() const;
	long countOrphanCells() const;
	virtual std::array<double, 3> evalCellCentroid(const long &id);
	std::vector<long> findCellNeighs(const long &id) const;
	std::vector<long> findCellNeighs(const long &id, int codimension, bool complete = true) const;
	std::vector<long> findCellFaceNeighs(const long &id) const;
	std::vector<long> findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> findCellEdgeNeighs(const long &id, bool complete = true) const;
	std::vector<long> findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> findCellVertexNeighs(const long &id, bool complete = true) const;
	std::vector<long> findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> findCellVertexNeighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> findCellVertexOneRing(const long &id, const int &vertex) const;

	CellIterator getCellIterator(const long &id);
	CellIterator cellBegin();
	CellIterator cellEnd();
	CellIterator internalBegin();
	CellIterator internalEnd();
	CellIterator ghostBegin();
	CellIterator ghostEnd();

	long getInterfaceCount() const;
	PiercedVector<Interface> &interfaces();
	Interface &getInterface(const long &id);
	const Interface &getInterface(const long &id) const;
	long generateInterfaceId();
	InterfaceIterator addInterface(const long &id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementInfo::Type type, const long &id = Element::NULL_ID);
	InterfaceIterator addInterface(Interface source);
	InterfaceIterator addInterface(Interface &&source, long id = Element::NULL_ID);
	bool deleteInterface(const long &id, bool updateNeighs = true, bool delayed = false);
	bool deleteInterfaces(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
	long countFreeInterfaces() const;
	long countOrphanInterfaces() const;
	virtual double evalInterfaceArea(const long &id) = 0;
	virtual std::array<double, 3> evalInterfaceCentroid(const long &id);
	virtual std::array<double, 3> evalInterfaceNormal(const long &id) = 0;

	InterfaceIterator getInterfaceIterator(const long &id);
	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

	long countFaces() const;
	long countFreeFaces() const;

	bool sort();
	bool sortVertices();
	bool sortCells();
	bool sortInterfaces();

	bool squeeze();
	bool squeezeVertices();
	bool squeezeCells();
	bool squeezeInterfaces();

	bool isPointInside(const double &x, const double &y, const double &z);
	virtual bool isPointInside(const std::array<double, 3> &point) = 0;
	long locatePoint(const double &x, const double &y, const double &z);
	virtual long locatePoint(const std::array<double, 3> &point) = 0;

	void updateBoundingBox();
	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);
	std::unordered_map<long, long> binSortVertex(int nBins = 128);

	virtual void translate(std::array<double, 3> translation);
	void translate(double sx, double sy, double sz);
	virtual void scale(std::array<double, 3> scaling);
	void scale(double scaling);
	void scale(double sx, double sy, double sz);

	void setTol(double tolerance);
	double getTol() const;
	void resetTol();
	bool isTolCustomized() const;

	void extractEnvelope(Patch &envelope) const;

	void displayStats(std::ostream &out, unsigned int padding = 0) const;
	void displayVertices(std::ostream &out, unsigned int padding = 0) const;
	void displayCells(std::ostream &out, unsigned int padding = 0) const;
	void displayInterfaces(std::ostream &out, unsigned int padding = 0) const;

	void write();
	void write(std::string name);
	void writeField(std::string name, VTKLocation location, std::vector<double> &values);
	void writeField(std::string filename, std::string name, VTKLocation location, std::vector<double> &values);
	void writeCellField(std::string name, std::vector<double> &values);
	void writeCellField(std::string filename, std::string name, std::vector<double> &values);
	void writeVertexField(std::string name, std::vector<double> &values);
	void writeVertexField(std::string filename, std::string name, std::vector<double> &values);

	const VTKFieldMetaData getMetaData(std::string name);
	void flushData(std::fstream &stream, VTKFormat format, std::string name);

#if ENABLE_MPI==1
	void setCommunicator(MPI::Intracomm *communicator);
	void unsetCommunicator();
	MPI::Comm & getCommunicator() const;
	int getRank() const;
	int getProcessorCount() const;
#endif

protected:
	PiercedVector<Vertex> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	IndexGenerator m_vertexIdGenerator;
	IndexGenerator m_interfaceIdGenerator;
	IndexGenerator m_cellIdGenerator;

	long m_nVertices;
	long m_nInternals;
	long m_nGhosts;
	long m_nInterfaces;

	long m_last_internal_id;
	long m_first_ghost_id;

	std::array<double, 3> m_minPoint;
	std::array<double, 3> m_maxPoint;

	bool deleteVertex(const long &id, bool delayed = false);
	bool deleteVertices(const std::vector<long> &ids, bool delayed = false);

	virtual void evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);

	virtual const std::vector<Adaption::Info> _update(bool trackAdaption) = 0;
	virtual bool _markCellForRefinement(const long &id) = 0;
	virtual bool _markCellForCoarsening(const long &id) = 0;
	virtual bool _enableCellBalancing(const long &id, bool enabled) = 0;
	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	void setDirty(bool dirty);
	void setExpert(bool expert);

private:
	bool m_dirty;
	bool m_expert;

	int m_id;
	int m_dimension;

	bool m_hasCustomTolerance;
	double m_tolerance;

	int m_rank;
	int m_nProcessors;
#if ENABLE_MPI==1
	MPI::Comm *m_communicator;
#endif

	VertexIterator createVertex(long id = Vertex::NULL_ID);
	InterfaceIterator createInterface(long id = Element::NULL_ID);
	CellIterator createCell(bool interior, long id = Element::NULL_ID);

	void set_id(int id);
	void setDimension(int dimension);

	std::array<double, 3> evalElementCentroid(const Element &element);
};

}

#endif
