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
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_IO.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

typedef PiercedVector<Vertex>::iterator VertexIterator;
typedef PiercedVector<Cell>::iterator CellIterator;
typedef PiercedVector<Interface>::iterator InterfaceIterator;

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
	Patch(const int &id, const int &dimension);

	virtual ~Patch();

	void reset();
	void resetVertices();
	void resetCells();
	void resetInterfaces();

	const std::vector<Adaption::Info> update(bool trackAdaption = true);

	void markCellForRefinement(const long &id);
	void markCellForCoarsening(const long &id);
	void enableCellBalancing(const long &id, bool enabled);

	bool isDirty() const;

	int get_id() const;
	int getDimension() const;
	bool isThreeDimensional() const;

	std::string getName() const;
	void setName(std::string name);

	long getVertexCount() const;
	PiercedVector<Vertex> &vertices();
	Vertex &getVertex(const long &id);
	const Vertex & getVertex(const long &id) const;
	const std::array<double, 3> & getVertexCoords(const long &id) const;

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
	virtual double evalCellVolume(const long &id) = 0;
	virtual double evalCellSize(const long &id) = 0;
	virtual std::array<double, 3> evalCellCentroid(const long &id);
	std::vector<long> extractCellNeighs(const long &id) const;
	std::vector<long> extractCellNeighs(const long &id, int codimension, bool complete = true) const;
	std::vector<long> extractCellFaceNeighs(const long &id) const;
	std::vector<long> extractCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extractCellEdgeNeighs(const long &id, bool complete = true) const;
	std::vector<long> extractCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extractCellVertexNeighs(const long &id, bool complete = true) const;
	std::vector<long> extractCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extractCellVertexNeighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList = std::vector<long>()) const;

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
	virtual double evalInterfaceArea(const long &id) = 0;
	virtual std::array<double, 3> evalInterfaceCentroid(const long &id);
	virtual std::array<double, 3> evalInterfaceNormal(const long &id) = 0;

	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

	void sort();
	void sortVertices();
	void sortCells();
	void sortInterfaces();

	void squeeze();
	void squeezeVertices();
	void squeezeCells();
	void squeezeInterfaces();

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

	void write();
	void write(std::string name);
	void writeField(std::string name, VTKLocation location, const std::vector<double> &values);
	void writeField(std::string filename, std::string name, VTKLocation location, const std::vector<double> &values);
	void writeCellField(std::string name, const std::vector<double> &values);
	void writeCellField(std::string filename, std::string name, const std::vector<double> &values);
	void writeVertexField(std::string name, const std::vector<double> &values);
	void writeVertexField(std::string filename, std::string name, const std::vector<double> &values);

	const VTKFieldMetaData getMetaData(std::string name);
	void flushData(std::fstream &stream, VTKFormat format, std::string name);

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

	long addVertex(const long &id = Vertex::NULL_VERTEX_ID);
	long addVertex(Vertex source);
	long addVertex(Vertex &&source, long id = Vertex::NULL_VERTEX_ID);
	void deleteVertex(const long &id, bool delayed = false);

	long addInterface(const long &id = Element::NULL_ELEMENT_ID);
	long addInterface(ElementInfo::Type type, const long &id = Element::NULL_ELEMENT_ID);
	long addInterface(Interface source);
	long addInterface(Interface &&source, long id = Element::NULL_ELEMENT_ID);
	void deleteInterface(const long &id, bool delayed = false);

	long addCell(const long &id = Element::NULL_ELEMENT_ID);
	long addCell(ElementInfo::Type type, bool interior, const long &id = Element::NULL_ELEMENT_ID);
	long addCell(Cell source);
	long addCell(Cell &&source, long id = Element::NULL_ELEMENT_ID);
	void deleteCell(const long &id, bool delayed = false);
	void setCellInternal(const long &id, bool isInternal);
	CellIterator moveGhost2Internal(const long &id);
	CellIterator moveInternal2Ghost(const long &id);

	virtual void evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);

	virtual const std::vector<Adaption::Info> _update(bool trackAdaption) = 0;
	virtual bool _markCellForRefinement(const long &id) = 0;
	virtual bool _markCellForCoarsening(const long &id) = 0;
	virtual bool _enableCellBalancing(const long &id, bool enabled) = 0;
	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	void setDirty(bool dirty);

private:
	bool m_dirty;

	int m_id;
	int m_dimension;
	std::string m_name;

	bool m_hasCustomTolerance;
	double m_tolerance;

	Vertex & createVertex(long id = Vertex::NULL_VERTEX_ID);
	Interface & createInterface(long id = Element::NULL_ELEMENT_ID);
	Cell & createCell(bool interior, long id = Element::NULL_ELEMENT_ID);

	std::unordered_map<std::string, const std::vector<double> *> m_dataFields;
	std::unordered_map<std::string, VTKLocation> m_dataLocations;
	std::unordered_map<std::string, VTKFieldType> m_dataType;

	void set_id(int id);
	void setDimension(int dimension);

	std::array<double, 3> evalElementCentroid(const Element &element);
};

}

#endif
