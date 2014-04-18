/*
 * Class_Local_Tree_3D.tpp
 *
 *  Created on: 18/apr/2014
 *      Author: Marco Cisternino
 */
// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

template<>
class Class_Local_Tree<3>{
	// ------------------------------------------------------------------------------- //
	// FRIENDSHIPS ------------------------------------------------------------------- //
	friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector< Class_Octant<3> > 	OctantsType;
	typedef vector<uint32_t>			u32vector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;


	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	OctantsType					octants;			// Local vector of octants ordered with Morton Number
	OctantsType					ghosts;				// Local vector of ghost octants ordered with Morton Number
	u32vector 					pborders;			// Local vector of pborder octants ordered with Morton Number
	Class_Octant<3> 			first_desc;			// First (Morton order) most refined octant possible in local partition
	Class_Octant<3> 			last_desc;			// Last (Morton order) most refined octant possible in local partition
	uint32_t 					size_ghosts;		// Size of vector of ghost octants
	uint8_t						local_max_depth;	// Reached max depth in local tree
public:
	u32vector2D					nodes;				// Local vector of nodes (x,y,z) ordered with Morton Number
	u32vector2D					connectivity;		// Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
													// The nodes are stored as index of vector nodes
	u32vector2D					ghostsnodes;		// Local vector of ghosts nodes (x,y,z) ordered with Morton Number
	u32vector2D					ghostsconnectivity;	// Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
													// The nodes are stored as index of vector nodes

//	// ------------------------------------------------------------------------------- //
//	// CONSTRUCTORS ------------------------------------------------------------------ //
//
//public:
//	Class_Local_Tree();
//	~Class_Local_Tree();
//
//	// ------------------------------------------------------------------------------- //
//	// METHODS ----------------------------------------------------------------------- //
//
//	// Basic Get/Set methods --------------------------------------------------------- //
//
//public:
//	const Class_Octant &  getFirstDesc() const;
//	const Class_Octant &  getLastDesc() const;
//	uint32_t  	  getSizeGhost() const;
//	uint32_t  	  getNumOctants() const;
//	uint8_t       getLocalMaxDepth() const;						// Get max depth reached in local tree
//	uint8_t       getMarker(int32_t idx);						// Get refinement/coarsening marker for idx-th octant
//	bool          getBalance(int32_t idx);						// Get if balancing-blocked idx-th octant
//
//	void		  setMarker(int32_t idx, int8_t marker);		// Set refinement/coarsening marker for idx-th octant
//	void		  setBalance(int32_t idx, bool balance);		// Set if balancing-blocked idx-th octant
//	void  		  setFirstDesc();
//	void  		  setLastDesc();
//
//	//-------------------------------------------------------------------------------- //
//	// Debug methods ----------------------------------------------------------------- //
//	void addOctantToTree(Class_Octant octant);
//
//private:
//
//	//-------------------------------------------------------------------------------- //
//	// Other Get/Set methods --------------------------------------------------------- //
//
//public:
//
//private:
//
//	//-------------------------------------------------------------------------------- //
//	// Other methods ----------------------------------------------------------------- //
//
//public:
//	const Class_Octant&	extractOctant(uint32_t idx) const ;
//	bool			refine();									// Refine local tree: refine one time octants with marker >0
//	bool			coarse();									// Coarse local tree: coarse one time family of octants with marker <0
//																// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
//	bool			refine(u32vector & mapidx);					// Refine local tree: refine one time octants with marker >0
//																// mapidx[i] = index in old octants vector of the i-th octant (index of father if octant is new after)
//	bool			coarse(u32vector & mapidx);					// Coarse local tree: coarse one time family of octants with marker <0
//																// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
//																// mapidx[i] = index in old octants vector of the i-th octant (index of father if octant is new after)
//	void			checkCoarse(uint64_t lastDescPre,			// Delete overlapping octants after coarse local tree. Check first and last descendants
//								uint64_t firstDescPost);		// of process before and after the local process
//	void       		updateLocalMaxDepth();						// Update max depth reached in local tree
//	void			computeConnectivity();						// Computes nodes vector and connectivity of octants of local tree
//	void			clearConnectivity();						// Clear nodes vector and connectivity of octants of local tree
//	void			updateConnectivity();						// Updates nodes vector and connectivity of octants of local tree
//	void			computeghostsConnectivity();				// Computes ghosts nodes vector and connectivity of ghosts octants of local tree
//	void			clearghostsConnectivity();					// Clear ghosts nodes vector and connectivity of ghosts octants of local tree
//	void			updateghostsConnectivity();					// Update ghosts nodes vector and connectivity of ghosts octants of local tree
//
//	void			findNeighbours(uint32_t idx,				// Finds neighbours of idx-th octant through iface in vector octants.
//								uint8_t iface,					// Returns a vector (empty if iface is a bound face) with the index of neighbours
//								u32vector & neighbours,			// in their structure (octants or ghosts) and sets isghost[i] = true if the
//								vector<bool> & isghost);		// i-th neighbour is ghost in the local tree
//
//	void			findGhostNeighbours(uint32_t const idx,		// Finds neighbours of idx-th ghost octant through iface in vector octants.
//								uint8_t iface,					// Returns a vector (empty if iface is not the pbound face for ghost) with the index of neighbours
//								u32vector & neighbours);		// in the structure octants
//
////	void			findNeighbours(Class_Octant const & oct,	// Finds neighbours of octant oct through iface in vector octants
////								uint8_t iface,					// It may be slow due to the research of oct in octants.
////								u32vector & neighbours,			// Returns a vector (empty if iface is a bound face) with the index of neighbours
////								vector<bool> & isghost);		// in their structure (octants or ghosts) and sets isghost[i] = true if the
////																// i-th neighbour is ghost in the local tree
//
//	bool			localBalance(bool doInterior);				// 2:1 balancing on level a local tree already adapted (balance only the octants with info[14] = false) (refinement wins!)
//																// Return true if balanced done with some markers modification
//																// Seto doInterior = false if the interior octants are already balanced
//
//	void			findEdgeNeighbours(uint32_t idx,			// Finds neighbours of idx-th octant through iedge in vector octants.
//									uint8_t iedge,				// Returns a vector (empty if iedge is a bound edge) with the index of neighbours
//									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
//									vector<bool> & isghost);	// i-th neighbour is ghost in the local tree
//
//	void			findNodeNeighbours(uint32_t idx,			// Finds neighbours of idx-th octant through inode in vector octants.
//									uint8_t inode,				// Returns a vector (empty if inode is a bound node) with the index of neighbours
//									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
//									vector<bool> & isghost);	// i-th neighbour is ghost in the local tree
//
//private:
//
//	// ------------------------------------------------------------------------------- //
//

};//end Class_Local_Tree;


