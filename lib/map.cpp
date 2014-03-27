
/*
 * map.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: edoardo
 */

#include "map.hpp"

Class_Map::Class_Map() {
	X0 = Y0 = Z0 = 0.0;
	L = 1.0;
}
Class_Map::Class_Map(double & X, double & Y, double & Z, double & LL) {
	X0 = X;
	Y0 = Y;
	Z0 = Z;
	L = LL;
}

double Class_Map::mapX(uint32_t const & X){
	return (X0 + L/double(max_length) * double(X));
}
double Class_Map::mapY(uint32_t const & Y){
	return (Y0 + L/double(max_length) * double(Y));
}
double Class_Map::mapZ(uint32_t const & Z){
	return (Z0 + L/double(max_length) * double(Z));
}
double Class_Map::mapSize(uint32_t const & size){
	return ((L/double(max_length))*double(size));
}
double Class_Map::mapArea(uint32_t const & Area){
	return ((pow(L,2.0)/pow(double(max_length),2.0))*double(Area));
}
double Class_Map::mapVolume(uint32_t const & Volume){
	return ((pow(L,3.0)/pow(double(max_length),3.0))*double(Volume));
}

void Class_Map::mapCenter(double* & center, vector<double> & mapcenter){
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

void Class_Map::mapNodes(uint32_t (*nodes)[DIM], vector<vector<double> > & mapnodes){
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

void Class_Map::mapNormals(vector<int8_t> normal_, vector<double> & mapnormal){
	mapnormal = vector<double>(normal_.begin(), normal_.end());
	mapnormal.shrink_to_fit();
}

// ------------------------------------------------------------------------------- //

