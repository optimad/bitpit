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
	// FRIENDSHIPS ------------------------------------------------------------------- //
	friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant> 		OctantsType;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;


	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	OctantsType					octants;			// Local vector of octants ordered with Morton Number
	OctantsType					ghosts;				// Local vector of ghost octants ordered with Morton Number
	Class_Octant 				first_desc;			// First (Morton order) most refined octant possible in local partition
	Class_Octant 				last_desc;			// Last (Morton order) most refined octant possible in local partition
	uint32_t 					size_ghosts;		// Size of vector of ghost octants
	uint8_t						local_max_depth;	// Reached max depth in local tree
public:
	u32vector2D					nodes;				// Local vector of nodes (x,y,z) ordered with Morton Number
	u32vector2D					connectivity;		// Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
													// The nodes are stored as index of vector nodes
	u32vector2D					ghostsnodes;		// Local vector of ghosts nodes (x,y,z) ordered with Morton Number
	u32vector2D					ghostsconnectivity;	// Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
													// The nodes are stored as index of vector nodes

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //

public:
	Class_Local_Tree();
	~Class_Local_Tree();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

public:
	const Class_Octant &  getFirstDesc() const;
	const Class_Octant &  getLastDesc() const;
	uint32_t  	  getSizeGhost() const;
	uint32_t  	  getNumOctants() const;
	uint8_t       getLocalMaxDepth() const;						// Get max depth reached in local tree
	uint8_t       getMarker(int32_t idx);						// Get refinement/coarsening marker for idx-th octant
	bool          getBalance(int32_t idx);						// Get if balancing-blocked idx-th octant

	void		  setMarker(int32_t idx, int8_t marker);		// Set refinement/coarsening marker for idx-th octant
	void		  setBalance(int32_t idx, bool balance);		// Set if balancing-blocked idx-th octant
	void  		  setFirstDesc();
	void  		  setLastDesc();

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
	const Class_Octant&	extractOctant(uint32_t idx) const ;
	bool			refine();									// Refine local tree: refine one time octants with marker >0
	bool			coarse();									// Coarse local tree: coarse one time family of octants with marker <0
																// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
	void       		updateLocalMaxDepth();						// Update max depth reached in local tree
	void			computeConnectivity();						// Computes nodes vector and connectivity of octants of local tree
	void			clearConnectivity();						// Clear nodes vector and connectivity of octants of local tree
	void			updateConnectivity();						// Updates nodes vector and connectivity of octants of local tree
	void			computeghostsConnectivity();				// Computes ghosts nodes vector and connectivity of ghosts octants of local tree
	void			clearghostsConnectivity();					// Clear ghosts nodes vector and connectivity of ghosts octants of local tree
	void			updateghostsConnectivity();					// Update ghosts nodes vector and connectivity of ghosts octants of local tree

	uint32_t*		findNeighbours(uint32_t idx,				// Finds n=sizeneigh neighbours of idx-th octant through iface in vector octants.
								  uint8_t iface,				// Returns a pointer to an array of size sizeneigh with the index of neighbours
								  uint8_t & sizeneigh,			// in their structure (octants or ghosts) and sets isghost = true if the
								  bool isghost);				// neighbours are ghost in the local tree

private:

	// ------------------------------------------------------------------------------- //


};//end Class_Local_Tree;




#endif /* CLASS_LOCAL_TREE_HPP_ */
