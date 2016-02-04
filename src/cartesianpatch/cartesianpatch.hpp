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

/*! \file */

#include <cstddef>
#include <memory>
#include <vector>

#include "patch.hpp"
#include "utils.hpp"

namespace bitpit {

/*!
	\ingroup cartesianpatch
	@{
*/

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
	static const int SPACE_MAX_DIM;

	std::vector<double> m_cellSize;
	std::vector<double> m_minCoord;

	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_z;

	std::vector<int> m_nCells1D;
	std::vector<int> m_nVertices1D;
	std::vector<int> m_x_nInterfaces1D;
	std::vector<int> m_y_nInterfaces1D;
	std::vector<int> m_z_nInterfaces1D;

	double m_cell_volume;
	double m_x_interface_area;
	double m_y_interface_area;
	double m_z_interface_area;

	std::vector<std::array<double, 3> > m_normals;

	void createVertices();

	void createCells();

	void createInterfaces();
	int countInterfacesDirection(const Vertex::Coordinate &direction);
	void createInterfaces_direction(const Vertex::Coordinate &direction);

	long cell_cartesianToLinear(const int &i, const int &j, const int &k) const;
	long cell_cartesianToLinear(const int ijk[]) const;
	long vertex_cartesianToLinear(const int &i, const int &j, const int &k) const;
	long vertex_cartesianToLinear(const int ijk[]) const;
	long interface_cartesianToLinear(const int &normal, const int &i, const int &j, const int &k) const;

};

/*!
	@}
*/

}

#endif
