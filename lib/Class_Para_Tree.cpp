/*
 * Class_Para_Tree.cpp
 *
 *  Created on: 12/feb/2014
 *      Author: Marco Cisternino
 */

#include "Class_Para_Tree.hpp"

Class_Para_Tree::Class_Para_Tree() {
	// TODO Auto-generated constructor stub
	error_flag = 0;
	max_depth = 0;
	global_num_octants = 0;
	error_flag = MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	partition_range_position = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];

}
Class_Para_Tree::~Class_Para_Tree() {
	// TODO Auto-generated destructor stub
	delete [] partition_range_position;
	delete [] partition_range_globalidx;
}

