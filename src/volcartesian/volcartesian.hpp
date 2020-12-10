/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#ifndef __BITPIT_VOLCARTESIAN_HPP__
#define __BITPIT_VOLCARTESIAN_HPP__

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "bitpit_patchkernel.hpp"

namespace bitpit {

class VolCartesian : public VolumeKernel {

public:
	using VolumeKernel::isPointInside;
	using PatchKernel::locatePoint;
	using PatchKernel::getCellType;
	using PatchKernel::getInterfaceType;

	enum MemoryMode {
		MEMORY_NORMAL,
		MEMORY_LIGHT
	};

	VolCartesian();
	VolCartesian(int dimension, const std::array<double, 3> &origin,
			   const std::array<double, 3> &lengths, const std::array<int, 3> &nCells);
	VolCartesian(int id, int dimension, const std::array<double, 3> &origin,
			   const std::array<double, 3> &lengths, const std::array<int, 3> &nCells);
	VolCartesian(int dimension, const std::array<double, 3> &origin,
			   double length, int nCells1D);
	VolCartesian(int id, int dimension, const std::array<double, 3> &origin,
			   double length, int nCells1D);
	VolCartesian(int dimension, const std::array<double, 3> &origin,
			   double length, double dh);
	VolCartesian(int id, int dimension, const std::array<double, 3> &origin,
			   double length, double dh);
	VolCartesian(std::istream &stream);

	std::unique_ptr<PatchKernel> clone() const override;

	void reset() override;

	void resetInterfaces() override;

	long getVertexCount() const override;
	int getVertexCount(int direction) const;

	long getCellCount() const override;
	int getCellCount(int direction) const;

	ElementType getCellType() const;
	ElementType getCellType(long id) const override;

	long getInterfaceCount() const override;
	ElementType getInterfaceType() const;
	ElementType getInterfaceType(long id) const override;

	std::array<double, 3> evalVertexCoords(long id) const;
	std::array<double, 3> evalVertexCoords(const std::array<int, 3> &ijk) const;
	const std::vector<double> & getVertexCoords(int direction) const;

	double evalCellVolume(long id) const override;
	double evalCellVolume(const std::array<int, 3> &ijk) const;
	double evalCellSize(long id) const override;
	double evalCellSize(const std::array<int, 3> &ijk) const;
	std::array<double, 3> evalCellCentroid(long id) const override;
	std::array<double, 3> evalCellCentroid(const std::array<int, 3> &ijk) const;
	const std::vector<double> & getCellCentroids(int direction) const;

	double evalInterfaceArea(long id) const override;
	std::array<double, 3> evalInterfaceNormal(long id) const override;

	std::array<double, 3> getSpacing() const;
	double getSpacing(int direction) const;

	void switchMemoryMode(MemoryMode mode);
	MemoryMode getMemoryMode();

	bool isPointInside(const std::array<double, 3> &point) const override;
	bool isPointInside(long id, const std::array<double, 3> &point) const override;
	long locatePoint(const std::array<double, 3> &point) const override;
	std::array<int, 3> locatePointCartesian(const std::array<double, 3> &point) const;
	long locateClosestVertex(std::array<double,3> const &point) const;
	std::array<int, 3> locateClosestVertexCartesian(std::array<double,3> const &point) const;
	long locateClosestCell(std::array<double,3> const &point);
	std::array<int, 3> locateClosestCellCartesian(std::array<double,3> const &point) const;

	std::vector<long> extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax);
	std::vector<long> extractCellSubSet(int idMin, int idMax);
	std::vector<long> extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax);
	std::vector<long> extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax);
	std::vector<long> extractVertexSubSet(int idMin, int idMax);
	std::vector<long> extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax);

	std::array<double, 3> getOrigin() const;
	void setOrigin(const std::array<double, 3> &origin);
	void translate(const std::array<double, 3> &translation) override;
	std::array<double, 3> getLengths() const;
	void setLengths(const std::array<double, 3> &lengths);
	void scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center) override;

	std::vector<double> convertToVertexData(const std::vector<double> &cellData) const;
	std::vector<double> convertToCellData(const std::vector<double> &nodeData) const;

	int linearCellInterpolation(std::array<double,3> &point, std::vector<int> &stencil, std::vector<double> &weights);
	int linearVertexInterpolation(std::array<double,3> &point, std::vector<int> &stencil, std::vector<double> &weights);

	long getCellLinearId(int i, int j, int k) const;
	long getCellLinearId(const std::array<int, 3> &ijk) const;
	std::array<int, 3> getCellCartesianId(long id) const;
	bool isCellCartesianIdValid(const std::array<int, 3> &ijk) const;
	long getVertexLinearId(int i, int j, int k) const;
	long getVertexLinearId(const std::array<int, 3> &ijk) const;
	std::array<int, 3> getVertexCartesianId(long id) const;
	std::array<int, 3> getVertexCartesianId(long cellId, int vertex) const;
	std::array<int, 3> getVertexCartesianId(const std::array<int, 3> &cellIjk, int vertex) const;
	bool isVertexCartesianIdValid(const std::array<int, 3> &ijk) const;

protected:
	VolCartesian(const VolCartesian &other) = default;

	std::vector<adaption::Info> _spawn(bool trackSpawn) override;

	void _updateInterfaces() override;

	int _getDumpVersion() const override;
	void _dump(std::ostream &stream) const override;
	void _restore(std::istream &stream) override;

	void _findCellFaceNeighs(long id, int face, const std::vector<long> &blackList, std::vector<long> *neighs) const override;
	void _findCellEdgeNeighs(long id, int edge, const std::vector<long> &blackList, std::vector<long> *neighs) const override;
	void _findCellVertexNeighs(long id, int vertex, const std::vector<long> &blackList, std::vector<long> *neighs) const override;

private:
	MemoryMode m_memoryMode;

	std::array<double, 3> m_cellSpacings;
	std::array<double, 3> m_minCoords;
	std::array<double, 3> m_maxCoords;

	std::array<double, 3> m_directionOrdering;

	std::array<std::vector<double>, 3> m_vertexCoords;
	std::array<std::vector<double>, 3> m_cellCenters;

	std::array<int, 3> m_nCells1D;
	std::array<int, 3> m_nVertices1D;

	long m_nVertices;
	long m_nCells;
	long m_nInterfaces;

	double m_cellVolume;
	std::array<double, 3> m_interfaceArea;
	std::array<std::array<double, 3>, 6> m_normals;

	std::vector<std::array<int, 3>> m_vertexNeighDeltas;
	std::vector<std::array<int, 3>> m_edgeNeighDeltas;
	std::vector<std::array<int, 2>> m_edgeFaces;

	void initialize();
	void initializeInterfaceArea();
	void initializeCellVolume();

	void setDiscretization(const std::array<int, 3> &nCells);

	void setMemoryMode(MemoryMode mode);

	void addVertices();

	void addCells();

	void addInterfaces();
	std::array<int, 3> getInterfaceCountDirection(int direction);
	void addInterfacesDirection(int direction);

	double evalCellVolume() const;
	double evalCellSize() const;

};

}

#endif
