/*
 * map.hpp
 *
 *  Created on: Mar 18, 2014
 *      Author: edoardo
 */

#ifndef MAP_HPP_
#define MAP_HPP_

#include "global.hpp"

inline double mapX(uint32_t const & X){
	return (X0 + L/double(max_length) * double(X));
}
inline double mapY(uint32_t const & Y){
	return (Y0 + L/double(max_length) * double(Y));
}
inline double mapZ(uint32_t const & Z){
	return (Z0 + L/double(max_length) * double(Z));
}
inline double mapSize(uint32_t const & size){
	return ((L/double(max_length))*double(size));
}
inline double mapArea(uint32_t const & Area){
	return ((pow(L,2.0)/pow(double(max_length),2.0))*double(Area));
}
inline double mapVolume(uint32_t const & Volume){
	return ((pow(L,3.0)/pow(double(max_length),3.0))*double(Volume));
}
inline void mapCenter(double* & center, vector<double> & mapcenter){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	if (DIM == 3)
		orig.push_back(Z0);
	orig.shrink_to_fit();
	mapcenter.resize(DIM);
	for (int i=0; i<DIM; i++){
		mapcenter[i] = orig[i] + L/double(max_length) * center[i];
	}
	mapcenter.shrink_to_fit();
}
inline void mapNodes(uint32_t (*nodes)[DIM], vector<vector<double> > & mapnodes){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	if (DIM == 3)
		orig.push_back(Z0);
	orig.shrink_to_fit();
	mapnodes.resize(nnodes);
	for (int i=0; i<nnodes; i++){
		mapnodes[i].resize(DIM);
		for (int j=0; j<DIM; j++){
			mapnodes[i][j] = orig[j] + L/double(max_length) * double(nodes[i][j]);
		}
		mapnodes[i].shrink_to_fit();
	}
	mapnodes.shrink_to_fit();
}

#endif /* MAP_HPP_ */
