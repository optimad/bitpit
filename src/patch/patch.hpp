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

#include "bitpit_IO.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class Patch : public bitpit::VTKUnstructuredGrid {

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
	bitpit::PiercedVector<Vertex> &vertices();
	Vertex &getVertex(const long &id);
	const Vertex & getVertex(const long &id) const;
	const std::array<double, 3> & getVertexCoords(const long &id) const;

	long getCellCount() const;
	bitpit::PiercedVector<Cell> &cells();
	Cell &getCell(const long &id);
	const Cell &getCell(const long &id) const;
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

	long getInterfaceCount() const;
	bitpit::PiercedVector<Interface> &interfaces();
	Interface &getInterface(const long &id);
	const Interface &getInterface(const long &id) const;
	virtual double evalInterfaceArea(const long &id) = 0;
	virtual std::array<double, 3> evalInterfaceCentroid(const long &id);
	virtual std::array<double, 3> evalInterfaceNormal(const long &id) = 0;

	void sort();
	void squeeze();

	bool isPointInside(const double &x, const double &y, const double &z);
	virtual bool isPointInside(const std::array<double, 3> &point) = 0;
	long locatePoint(const double &x, const double &y, const double &z);
	virtual long locatePoint(const std::array<double, 3> &point) = 0;

	void updateBoundingBox();
	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);

	void writeMesh();
	void writeMesh(std::string name);
	void writeField(std::string name, bitpit::VTKLocation location, const std::vector<double> &values);
	void writeField(std::string filename, std::string name, bitpit::VTKLocation location, const std::vector<double> &values);
	void writeCellField(std::string name, const std::vector<double> &values);
	void writeCellField(std::string filename, std::string name, const std::vector<double> &values);
	void writeVertexField(std::string name, const std::vector<double> &values);
	void writeVertexField(std::string filename, std::string name, const std::vector<double> &values);

	const bitpit::VTKFieldMetaData getMetaData(std::string name);
	void flushData(std::fstream &stream, bitpit::VTKFormat format, std::string name);

protected:
	bitpit::PiercedVector<Vertex> m_vertices;
	bitpit::PiercedVector<Cell> m_cells;
	bitpit::PiercedVector<Interface> m_interfaces;

	std::deque<long> m_unusedVertexIds;
	std::deque<long> m_unusedInterfaceIds;
	std::deque<long> m_unusedCellIds;

	std::array<double, 3> m_minPoint;
	std::array<double, 3> m_maxPoint;

	long createVertex();
	long createVertex(const long &id);
	void deleteVertex(const long &id, bool delayed = false);

	long createInterface(ElementInfo::Type type = ElementInfo::UNDEFINED);
	long createInterface(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);
	void deleteInterface(const long &id, bool delayed = false);

	long createCell(bool internal = true, ElementInfo::Type type = ElementInfo::UNDEFINED);
	long createCell(const long &id, bool internal = true, ElementInfo::Type type = ElementInfo::UNDEFINED);
	void deleteCell(const long &id, bool delayed = false);

	virtual void evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);

	virtual const std::vector<Adaption::Info> _update(bool trackAdaption) = 0;
	virtual bool _markCellForRefinement(const long &id) = 0;
	virtual bool _markCellForCoarsening(const long &id) = 0;
	virtual bool _enableCellBalancing(const long &id, bool enabled) = 0;

	void setDirty(bool dirty);

private:
	bool m_dirty;

	int m_id;
	int m_dimension;
	std::string m_name;

	std::unordered_map<std::string, const std::vector<double> *> m_dataFields;
	std::unordered_map<std::string, bitpit::VTKLocation> m_dataLocations;
	std::unordered_map<std::string, bitpit::VTKFieldType> m_dataType;

	void set_id(int id);
	void setDimension(int dimension);

	std::array<double, 3> evalElementCentroid(const Element &element);
};

}

#endif
