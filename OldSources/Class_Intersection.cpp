/*
 * Class_Intersection.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: edoardo
 */

#include "Class_Intersection.hpp"

Class_Intersection::Class_Intersection(){
	owners[0] = 0;
	owners[1] = 0;
	iface = 0;
	isnew = false;
	isghost = false;
	finer = 0;
	octree = NULL;
}


Class_Intersection::Class_Intersection(Class_Local_Tree & tree){
	owners[0] = 0;
	owners[1] = 0;
	iface = 0;
	isnew = false;
	isghost = false;
	finer = 0;
	octree = &tree;
}

Class_Intersection::~Class_Intersection() {
	octree = NULL;
}

Class_Intersection::Class_Intersection(const Class_Intersection & intersection){
	owners[0] = intersection.owners[0];
	owners[1] = intersection.owners[1];
	iface = intersection.iface;
	isnew = intersection.isnew;
	isghost = intersection.isghost;
	finer = intersection.finer;
	octree = intersection.octree;
}

Class_Intersection& Class_Intersection::operator =(const Class_Intersection & intersection){
	owners[0] = intersection.owners[0];
	owners[1] = intersection.owners[1];
	iface = intersection.iface;
	isnew = intersection.isnew;
	isghost = intersection.isghost;
	finer = intersection.finer;
	octree = intersection.octree;
	return *this;
}

bool Class_Intersection::operator ==(const Class_Intersection & intersection){
	bool check = true;
	check = check && (owners[0] == intersection.owners[0]);
	check = check && (owners[1] == intersection.owners[1]);
	check = check && (iface == intersection.iface);
	check = check && (isnew == intersection.isnew);
	check = check && (isghost == intersection.isghost);
	check = check && (finer == intersection.finer);
	check = check && (octree == intersection.octree);
	return check;
}



uint32_t Class_Intersection::getOut() {
	return owners[0];
}

uint32_t Class_Intersection::getIn() {
	return owners[1];
}

void Class_Intersection::getNormal(int8_t normal[DIM]) {
	for (int i=0; i<DIM; i++){
		normal[i] = normals[iface][i];
	}
}

uint32_t Class_Intersection::getSize() {
	Class_Octant oct = octree->extractOctant(owners[finer]);
	return oct.getSize();
}

uint64_t Class_Intersection::getArea() {
	Class_Octant oct = octree->extractOctant(owners[finer]);
	return oct.getArea();
}

void Class_Intersection::getNodes(u32vector2D & nodes) {
	uint8_t		i;
	uint32_t	dh;
	Class_Octant oct = octree->extractOctant(owners[finer]);
	uint32_t (*nodes_)[DIM] = oct.getNodes();

	dh = oct.getSize();
	nodes.clear();
	nodes.resize(nnodesperface);

	for (i = 0; i < nnodesperface; i++){
		nodes[i].resize(DIM);
		nodes[i][0] = nodes_[facenode[iface][i]][0];
		nodes[i][1] = nodes_[facenode[iface][i]][1];
		if (DIM ==3){
			nodes[i][2] = nodes_[facenode[iface][i]][2];
		}
		nodes[i].shrink_to_fit();
	}
	nodes.shrink_to_fit();
	delete [] nodes_;
	nodes_ = NULL;
}

uint32_t (*Class_Intersection::getNodes())[DIM] {
	uint8_t		i;
	uint32_t	dh;
	Class_Octant oct = octree->extractOctant(owners[finer]);
	uint32_t (*nodes_)[DIM] = oct.getNodes();

	dh = oct.getSize();
	uint32_t (*nodes)[DIM] = new uint32_t[nnodesperface][DIM];

	for (i = 0; i < nnodesperface; i++){
		nodes[i][0] = nodes_[facenode[iface][i]][0];
		nodes[i][1] = nodes_[facenode[iface][i]][1];
		if (DIM ==3){
			nodes[i][2] = nodes_[facenode[iface][i]][2];
		}
	}
	delete [] nodes_;
	nodes_ = NULL;
	return nodes;
}


void Class_Intersection::getCenter(vector<double> & center) {
	uint8_t		i;
	double	dh;
	int8_t normal[DIM];
	Class_Octant oct = octree->extractOctant(owners[finer]);
	uint32_t (*nodes_)[DIM] = oct.getNodes();

	dh = double(oct.getSize())/2.0;
	getNormal(normal);
	center.resize(DIM);
	center[0] = double(nodes_[facenode[iface][0]][0]) + dh*double(!normal[0]);
	center[1] = double(nodes_[facenode[iface][0]][1]) + dh*double(!normal[1]);
#if DIM == 3
	center[1] = double(nodes_[facenode[iface][0]][2]) + dh*double(!normal[2]);
#endif
	delete [] nodes_;
	nodes_ = NULL;
}

double*	Class_Intersection::getCenter() {
	uint8_t		i;
	double	dh;

	int8_t normal[DIM];
	Class_Octant oct = octree->extractOctant(owners[finer]);
	uint32_t (*nodes_)[DIM] = oct.getNodes();
	dh = double(oct.getSize())/2.0;
	getNormal(normal);

	double *center = new double[DIM];
	center[0] = double(nodes_[facenode[iface][0]][0]) + dh*double(!normal[0]);
	center[1] = double(nodes_[facenode[iface][0]][1]) + dh*double(!normal[1]);
#if DIM == 3
	center[1] = double(nodes_[facenode[iface][0]][2]) + dh*double(!normal[2]);
#endif
	delete [] nodes_;
	nodes_ = NULL;
	return center;
}





