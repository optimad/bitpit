/*
 * Class_Local_Tree.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#ifndef CLASS_LOCAL_TREE_HPP_
#define CLASS_LOCAL_TREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "Class_Octant.hpp"
#include <math.h>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <map>
#include <iostream>


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //


class Class_Local_Tree{

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	vector<Class_Octant>		octants;			// Local vector of octants ordered with Morton Number
	vector<Class_Octant>		ghosts;				// Local vector of ghost octants ordered with Morton Number
	Class_Octant 				first_desc;			// First (Morton order) most refined octant possible in local partition
	Class_Octant 				last_desc;			// Last (Morton order) most refined octant possible in local partition
	uint32_t 					size_ghosts;		// Size of vector of ghost octants
	uint8_t						local_max_depth;	// Reached max depth in local tree
public:
	vector<vector<uint32_t>	>	nodes;				// Local vector of nodes (x,y,z) ordered with Morton Number
	vector<vector<uint64_t>	>	connectivity;		// Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
													//The nodes are stored as index of vector nodes

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //

public:
	Class_Local_Tree();
	~Class_Local_Tree();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

public:
	Class_Octant  getFirstDesc() const;
	Class_Octant  getLastDesc() const;
	uint32_t  	  getSizeGhost() const;
	uint64_t  	  getNumOctants() const;
	uint8_t       getLocalMaxDepth() const;						// Get max depth reached in local tree
	uint8_t       getMarker(int64_t idx);						// Get refinement/coarsening marker for idx-th octant
	bool          getBalance(int64_t idx);						// Get if balancing-blocked idx-th octant

	void		  setMarker(int64_t idx, int8_t marker);		// Set refinement/coarsening marker for idx-th octant
	void		  setBalance(int64_t idx, bool balance);		// Set if balancing-blocked idx-th octant

	//-------------------------------------------------------------------------------- //
	// Debug methods ----------------------------------------------------------------- //
	void addOctantToTree(Class_Octant octant);

private:

	//-------------------------------------------------------------------------------- //
	// Other Get/Set methods --------------------------------------------------------- //

public:

private:

	//-------------------------------------------------------------------------------- //
	// Other methods ----------------------------------------------------------------- //

public:
	Class_Octant&	extractOctant(uint64_t idx);
	void			refine();							// Refine local tree: refine one time octants with marker >0
	void			computeConnectivity();				// Computes nodes vector and connectivity of octants of local tree
	void			clearConnectivity();				// Clear nodes vector and connectivity of octants of local tree

private:

	// ------------------------------------------------------------------------------- //


};//end Class_Local_Tree;




#endif /* CLASS_LOCAL_TREE_HPP_ */
