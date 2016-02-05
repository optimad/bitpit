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
	CartesianPatch(const int &id, const int &dimension, std::array<double, 3> origin,
			   double length, double dh);

	~CartesianPatch();

	double evalCellVolume(const long &id);
	double evalCellSize(const long &id);

	double evalInterfaceArea(const long &id);
	std::array<double, 3> evalInterfaceNormal(const long &id);

protected:
	const std::vector<Adaption::Info> _update(bool trackAdaption);
	bool _markCellForRefinement(const long &id);
	bool _markCellForCoarsening(const long &id);
	bool _enableCellBalancing(const long &id, bool enabled);

private:
	std::vector<double> m_cellSize;
	std::vector<double> m_minCoord;

	std::array<std::vector<double>, 3> m_vertexCoords;

	std::vector<int> m_nCells1D;
	std::array<int, 3> m_nVertices1D;
	std::vector<int> m_nInterfacesX1D;
	std::vector<int> m_nInterfacesY1D;
	std::vector<int> m_nInterfacesZ1D;

	double m_cellVolume;
	std::array<double, 3> m_interfaceArea;

	std::vector<std::array<double, 3> > m_normals;

	void createVertices();

	void createCells();

	void createInterfaces();
	int countInterfacesDirection(const Vertex::Coordinate &direction);
	void createInterfacesDirection(const Vertex::Coordinate &direction);

	long getCellLinearId(const int &i, const int &j, const int &k) const;
	long getCellLinearId(const std::array<int, 3> &ijk) const;
	long getVertexLinearId(const int &i, const int &j, const int &k) const;
	long getVertexLinearId(const std::array<int, 3> &ijk) const;
	long getInterfaceLinearId(const int &normal, const int &i, const int &j, const int &k) const;

};

}

#endif
