#ifndef CLASSLOCALTREE_HPP_
#define CLASSLOCALTREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classOctant.hpp"
#include "classIntersection.hpp"
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <map>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <bitset>
#include <array>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //

typedef vector<classOctant>		 		octvector;
typedef vector<classIntersection>	 	intervector;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 * classLocalTree.hpp
 *
 *  \ingroup    PABLO
 *  @{
 *  \date		15/dec/2015
 *	\authors	Edoardo Lombardi
 *	\authors	Marco Cisternino
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Local octree portion for each process
 *
 *	Local tree consists mainly of two vectors with:
 *	- actual octants stored on current process;
 *	- ghost octants neighbours of the first ones.
 *
 *	The octants (and ghosts) are ordered following the Z-curve defined by the Morton index.
 *
 *	Optionally in local tree three vectors of intersections are stored:
 *	- intersections located on the physical domain of the octree;
 *	- intersections of process bord (i.e. between octants and ghosts);
 *	- intersections completely located in the domain of the process (i.e. between actual octants).
 *
 *	Class classLocalTree is built with a dimensional parameter int dim and it accepts only two values: dim=2 and dim=3, for 2D and 3D respectively.
 */
class classLocalTree{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class classParaTree;

	// =================================================================================== //
	// TYPEDEFS
	// =================================================================================== //
public:
	typedef vector<classOctant>		 		octvector;
	typedef vector<classIntersection>	 	intervector;
	typedef vector<uint8_t>					u8vector;
	typedef vector<uint32_t>				u32vector;
	typedef vector<uint64_t>				u64vector;
	typedef vector<u32array3>				u32arr3vector;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //

private:
	octvector					octants;			/**< Local vector of octants ordered with Morton Number */
	octvector					ghosts;				/**< Local vector of ghost octants ordered with Morton Number */
	intervector					intersections;		/**< Local vector of intersections */
	u64vector 					globalidx_ghosts;	/**< Global index of the ghost octants (size = size_ghosts) */
	classOctant 				first_desc;			/**< First (Morton order) most refined octant possible in local partition */
	classOctant		 			last_desc;			/**< Last (Morton order) most refined octant possible in local partition */
	uint32_t 					size_ghosts;		/**< Size of vector of ghost octants */
	uint8_t						local_max_depth;	/**< Reached max depth in local tree */
	uint8_t 					balance_codim;		/**<Maximum codimension of the entity for 2:1 balancing (1 = 2:1 balance through faces (default);
	 	 	 	 	 	 	 	 	 	 	 	 	 	 2 = 2:1 balance through edges and faces;
	 	 	 	 	 	 	 	 	 	 	 	 	 	 3 = 2:1 balance through nodes, edges and faces)*/
	u32vector 					last_ghost_bros;	/**<Index of ghost brothers in case of broken family coarsened*/

	// connectivity
	//	u32vector2D					nodes;				/**<Local vector of nodes (x,y,z) ordered with Morton Number*/
	u32vector2D					connectivity;		/**<Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	//	u32vector2D					ghostsnodes;		/**<Local vector of ghosts nodes (x,y,z) ordered with Morton Number*/
	u32vector2D					ghostsconnectivity;	/**<Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	u32arr3vector				nodes;				/**<Local vector of nodes (x,y,z) ordered with Morton Number*/
	u32arr3vector				ghostsnodes;		/**<Local vector of ghosts nodes (x,y,z) ordered with Morton Number*/

	uint8_t						dim;				/**<Space dimension. Only 2D or 3D space accepted*/
	classGlobal					global;				/**<Global variables*/

	// =================================================================================== //
	// CONSTRUCTORS
	// =================================================================================== //

public:
	/*!Dimensional and default constructor.
	 */
	classLocalTree(int8_t maxlevel, uint8_t dim_);

	/*!Default destructor.
	 */
	~classLocalTree();

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
private:
	const classOctant &  getFirstDesc() const;
	const classOctant &  getLastDesc() const;
	uint32_t getSizeGhost() const;
	uint32_t getNumOctants() const;
	uint8_t getLocalMaxDepth() const;
	int8_t getMarker(int32_t idx);
	uint8_t getLevel(int32_t idx);
	uint8_t getGhostLevel(int32_t idx);
	bool getBalance(int32_t idx);
	/*! Get the codimension for 2:1 balancing
	 * \return Maximum codimension of the entity through which the 2:1 balance is performed.
	 */
	uint8_t getBalanceCodim() const;
	void setMarker(int32_t idx, int8_t marker);
	void setBalance(int32_t idx, bool balance);
	/*! Set the codimension for 2:1 balancing
	 * \param[in] Maximum codimension of the entity through which the 2:1 balance is performed.
	 */
	void setBalanceCodim(uint8_t b21codim);
	void setFirstDesc();
	void setLastDesc();

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

	// =================================================================================== //
	// OTHER METHODS
	// =================================================================================== //

	classOctant& extractOctant(uint32_t idx);
	const classOctant&	extractOctant(uint32_t idx) const;
	classOctant& extractGhostOctant(uint32_t idx);
	const classOctant& extractGhostOctant(uint32_t idx) const;


	bool refine(u32vector & mapidx);
	bool coarse(u32vector & mapidx);
	bool globalRefine(u32vector & mapidx);
	bool globalCoarse(u32vector & mapidx);
	void checkCoarse(uint64_t lastDescPre, uint64_t firstDescPost, u32vector & mapidx);
	void updateLocalMaxDepth();


	void findNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours,
			vector<bool> & isghost);
	void findNeighbours(classOctant* oct, uint8_t iface, u32vector & neighbours,
			vector<bool> & isghost);
	void findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours);
	void findEdgeNeighbours(uint32_t idx, uint8_t iedge,
			u32vector & neighbours, vector<bool> & isghost);
	void findEdgeNeighbours(classOctant* oct, uint8_t iedge,
			u32vector & neighbours, vector<bool> & isghost);
	void findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge,
			u32vector & neighbours);
	void findNodeNeighbours(classOctant* oct, uint8_t inode,
			u32vector & neighbours, vector<bool> & isghost);
	void findNodeNeighbours(uint32_t idx, uint8_t inode,
			u32vector & neighbours, vector<bool> & isghost);
	void findGhostNodeNeighbours(uint32_t idx, uint8_t inode,
			u32vector & neighbours);


	void preBalance21(bool internal);
	void preBalance21(u32vector& newmodified);
	bool localBalance(bool doInterior);
	bool localBalanceAll(bool doInterior);

	void computeIntersections();

	uint32_t findMorton(uint64_t Morton);
	uint32_t findGhostMorton(uint64_t Morton);

	void computeConnectivity();
	void clearConnectivity();
	void updateConnectivity();
	void computeGhostsConnectivity();
	void clearGhostsConnectivity();
	void updateGhostsConnectivity();

	// =================================================================================== //


};



#endif /* CLASSLOCALTREE_HPP_ */
