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

#ifndef __BITPIT_CARTESIANPATCH_HPP__
#define __BITPIT_CARTESIANPATCH_HPP__

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "bitpit_patch.hpp"

namespace bitpit {

class CartesianPatch : public Patch {

public:
	CartesianPatch(const int &id, const int &dimension, const std::array<double, 3> &origin,
			   const std::array<double, 3> &lengths, const std::array<int, 3> &nCells);
	CartesianPatch(const int &id, const int &dimension, const std::array<double, 3> &origin,
			   double length, int nCells1D);
	CartesianPatch(const int &id, const int &dimension, const std::array<double, 3> &origin,
			   double length, double dh);

	~CartesianPatch();

	double evalCellVolume(const long &id);
	double evalCellSize(const long &id);
	std::array<double, 3> evalCellCentroid(const long &id);

	double evalInterfaceArea(const long &id);
	std::array<double, 3> evalInterfaceNormal(const long &id);

	std::array<double, 3> getSpacing() const;
	double getSpacing(const int &direction) const;

	bool isPointInside(const std::array<double, 3> &point);
	long locatePoint(const std::array<double, 3> &point);
	std::array<int, 3> locatePointCartesian(const std::array<double, 3> &point);
	long locateClosestVertex(std::array<double,3> const &point) const;
	std::array<int, 3> locateClosestVertexCartesian(std::array<double,3> const &point) const;

	std::vector<long> extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax);
	std::vector<long> extractCellSubSet(int const &idxMin, int const &idxMax);
	std::vector<long> extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax);
	std::vector<long> extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax);
	std::vector<long> extractVertexSubSet(int const &idxMin, int const &idxMax);
	std::vector<long> extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax);

	void translate(std::array<double, 3> translation);
	void scale(std::array<double, 3> scaling);

	std::vector<double> convertToVertexData(const std::vector<double> &cellData) const;
	std::vector<double> convertToCellData(const std::vector<double> &nodeData) const;

	int linearCellInterpolation(std::array<double,3> &point,
		std::vector<int> &stencil, std::vector<double> &weights);
	int linearVertexInterpolation(std::array<double,3> &point,
		std::vector<int> &stencil, std::vector<double> &weights);

protected:
	const std::vector<Adaption::Info> _update(bool trackAdaption);
	bool _markCellForRefinement(const long &id);
	bool _markCellForCoarsening(const long &id);
	bool _enableCellBalancing(const long &id, bool enabled);

	void evalBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);

private:
	std::array<double, 3> m_cellSpacings;
	std::array<double, 3> m_minCoords;
	std::array<double, 3> m_maxCoords;

	std::array<std::vector<double>, 3> m_vertexCoords;
	std::array<std::vector<double>, 3> m_cellCenters;

	std::array<int, 3> m_nCells1D;
	std::array<int, 3> m_nVertices1D;

	double m_cellVolume;
	std::array<double, 3> m_interfaceArea;
	std::array<std::array<double, 3>, 6> m_normals;

	void initialize(const std::array<double, 3> &origin, const std::array<double, 3> &lengths,
	                const std::array<int, 3> &nCells);

	void initializeInterfaceArea();
	void initializeCellVolume();

	void createVertices();

	void createCells();

	void createInterfaces();
	std::array<int, 3> getInterfaceCountDirection(const int &direction);
	void createInterfacesDirection(const int &direction);

	long getCellLinearId(const int &i, const int &j, const int &k) const;
	long getCellLinearId(const std::array<int, 3> &ijk) const;
	std::array<int, 3> getCellCartesianId(long const &idx) const;
	long getVertexLinearId(const int &i, const int &j, const int &k) const;
	long getVertexLinearId(const std::array<int, 3> &ijk) const;
	std::array<int, 3> getVertexCartesianId(long const &idx) const;

};

}

#endif
