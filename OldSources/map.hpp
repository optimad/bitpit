/*
 * map.hpp
 *
 *  Created on: Mar 18, 2014
 *      Author: edoardo
 */

#ifndef MAP_HPP_
#define MAP_HPP_

#include "global.hpp"
#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <map>
#include <iostream>
#include <vector>
#include <string.h>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //

// Class_Map allows only translation and scaling
// Default : unit cube with origin in (0,0,0)
class Class_Map{

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	double X0, Y0, Z0;
	double L;

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Map();
	Class_Map(double & X, double & Y, double & Z, double & LL);

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	double mapX(uint32_t const & X);
	double mapY(uint32_t const & Y);
	double mapZ(uint32_t const & Z);
	uint32_t mapX(double const & X);
	uint32_t mapY(double const & Y);
	uint32_t mapZ(double const & Z);
	double mapSize(uint32_t const & size);
	double mapArea(uint32_t const & Area);
	double mapVolume(uint32_t const & Volume);
	void mapCenter(double* & center,
					vector<double> & mapcenter);
	void mapNodes(uint32_t (*nodes)[DIM],
					vector<vector<double> > & mapnodes);
	void mapNodesIntersection(uint32_t (*nodes)[DIM],
					vector<vector<double> > & mapnodes);
	void mapNormals(vector<int8_t> normal,
					vector<double> & mapnormal);
	// ------------------------------------------------------------------------------- //

};
// end class_map

#endif /* MAP_HPP_ */
