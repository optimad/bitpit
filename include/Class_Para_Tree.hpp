/*
 * Class_Para_Tree.hpp
 *
 *  Created on: 12/feb/2014
 *      Author: Marco Cisternino
 */

#ifndef CLASS_PARA_TREE_H_
#define CLASS_PARA_TREE_H_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "Class_Octant.hpp"
#include "Class_Local_Tree.hpp"
#include <cstdint>
#include "mpi.h"


// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

class Class_Para_Tree {
	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	//undistributed members
	uint64_t* partition_range_position; //global array containing position of the last existing octant in each processor
	uint64_t* partition_range_globalidx; //global array containing global index of the last existing octant in each processor
	uint64_t global_num_octants; // global number of octants in the parallel octree
	int nproc;
	uint8_t max_depth; // global max existing level in the parallel octree

	//distributed members
	int rank;
	Class_Local_Tree octree; // local tree in each processor

	//auxiliary members
	int error_flag; // MPI error flag
	bool serial; // 1 if the octree is the same on each processor, 0 if the octree is distributed



	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Para_Tree();
	~Class_Para_Tree();
	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	void loadBalance();
	void refine();
	void update();
};

#endif /* CLASS_PARA_TREE_H_ */
