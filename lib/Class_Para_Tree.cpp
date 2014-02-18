/*
 * Class_Para_Tree.cpp
 *
 *  Created on: 12/feb/2014
 *      Author: Marco Cisternino
 */

#include "Class_Para_Tree.hpp"

Class_Para_Tree::Class_Para_Tree() {
	serial = true;
	error_flag = 0;
	max_depth = 0;
	global_num_octants = octree.getnumoctants();
	error_flag = MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	partition_range_position = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];
	for(int p = 0; p < nproc; ++p){
		partition_range_globalidx[p] = 0;
		partition_range_position[p] = 0;
	}

}
Class_Para_Tree::~Class_Para_Tree() {
	delete [] partition_range_position;
	delete [] partition_range_globalidx;
}

void Class_Para_Tree::loadBalance(){
	serial = false;
	int divisionResult = 0;
	int remind = 0;
	divisionResult = global_num_octants/nproc;
	remind = global_num_octants%nproc;

}

void Class_Para_Tree::update() {
	if(serial)
	{}
}

void Class_Para_Tree::refine() {
	octree.refine();
	update();
}
