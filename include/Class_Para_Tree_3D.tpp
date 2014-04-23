/*
 * Class_Para_Tree_3D.tpp
 *
 *  Created on: 23/apr/2014
 *      Author: Marco Cisternino
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

template<>
class Class_Para_Tree<3>{
	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant<3> > 		OctantsType;
	typedef vector<uint32_t>			u32vector;
	typedef vector<double>				dvector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;
	typedef vector<vector<double>	>	dvector2D;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	//undistributed members
	uint64_t* partition_first_desc; 			//global array containing position of the first possible octant in each processor
	uint64_t* partition_last_desc; 				//global array containing position of the last possible octant in each processor
	uint64_t* partition_range_globalidx;	 	//global array containing global index of the last existing octant in each processor
	uint64_t global_num_octants;   				// global number of octants in the parallel octree
	map<int,vector<uint32_t> > bordersPerProc;	//local indices of border octants per process
	int nproc;
	uint8_t max_depth;							// global max existing level in the parallel octree

	//distributed members
	int rank;
	Class_Local_Tree<3> octree;					// local tree in each processor

	//auxiliary members
	int error_flag;								// MPI error flag
	bool serial;								// 1 if the octree is the same on each processor, 0 if the octree is distributed

	//map member
	Class_Map<3> trans;

	// connectivity
	dvector2D					nodes;				// Local vector of nodes (x,y,z) ordered with Morton Number
	u32vector2D					connectivity;		// Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	// The nodes are stored as index of vector nodes
	dvector2D					ghostsnodes;		// Local vector of ghosts nodes (x,y,z) ordered with Morton Number
	u32vector2D					ghostsconnectivity;	// Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	// The nodes are stored as index of vector nodes

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:


};


