//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PATCH_CARTESIAN_HPP__
#define __PATCHMAN_PATCH_CARTESIAN_HPP__

/*! \file */

#include "patch.hpp"
#include "utils.hpp"

#include <cstddef>
#include <memory>
#include <vector>

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
