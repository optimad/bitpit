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
#include "Class_Comm_Buffer.hpp"
#include "Class_Array.hpp"
#include <cstdint>
#include <iterator>
#include <set>
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
	uint64_t* partition_last_desc; 				//global array containing position of the last existing octant in each processor
	uint64_t* partition_range_globalidx;	 	//global array containing global index of the last existing octant in each processor
	uint64_t global_num_octants;   				// global number of octants in the parallel octree
	int nproc;
	uint8_t max_depth;							// global max existing level in the parallel octree

	//distributed members
	int rank;
	Class_Local_Tree octree;					// local tree in each processor

	//auxiliary members
	int error_flag;								// MPI error flag
	bool serial;								// 1 if the octree is the same on each processor, 0 if the octree is distributed



	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Para_Tree();
	~Class_Para_Tree();
	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	void loadBalance();							//assign the octants to the processes following a computed partition
	void adapt();  								//call refine and coarse on the local tree
	void updateAdapt();							//update Class_Para_Tree members after a refine and/or coarse
	void updateLoadBalance();					//update Class_Para_Tree members after a load balance
	void computePartition(uint32_t* partition); // compute octant partition giving the same number of octant to each process and redistributing the reminder
	void buildGhosts();
	int findOwner(const uint64_t & morton);		// given the morton of an octant it find the process owning that octant
	void setPboundGhosts(); 			 		// set pbound and build ghosts after static load balance
};

#endif /* CLASS_PARA_TREE_H_ */
