#ifndef CLASSLOCALTREE_HPP_
#define CLASSLOCALTREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ClassGlobal.hpp"
#include "ClassOctant.hpp"
#include "ClassIntersection.hpp"

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //

typedef std::vector<ClassOctant>		 	octvector;
typedef std::vector<ClassIntersection>	 	intervector;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
/*!
 * ClassLocalTree.hpp
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
 *	Class ClassLocalTree is built with a dimensional parameter int dim and it accepts only two values: dim=2 and dim=3, for 2D and 3D respectively.
 */
class ClassLocalTree{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class ClassParaTree;

	// =================================================================================== //
	// TYPEDEFS
	// =================================================================================== //
public:
	typedef std::vector<ClassOctant>		 	octvector;
	typedef std::vector<ClassIntersection>	 	intervector;
	typedef std::vector<uint8_t>				u8vector;
	typedef std::vector<uint32_t>				u32vector;
	typedef std::vector<uint64_t>				u64vector;
	typedef std::vector<u32array3>				u32arr3vector;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //

private:
	octvector				m_octants;				/**< Local vector of octants ordered with Morton Number */
	octvector				m_ghosts;				/**< Local vector of ghost octants ordered with Morton Number */
	intervector				m_intersections;		/**< Local vector of intersections */
	u64vector 				m_globalIdxGhosts;		/**< Global index of the ghost octants (size = size_ghosts) */
	ClassOctant 			m_firstDesc;			/**< First (Morton order) most refined octant possible in local partition */
	ClassOctant		 		m_lastDesc;				/**< Last (Morton order) most refined octant possible in local partition */
	uint32_t 				m_sizeGhosts;			/**< Size of vector of ghost octants */
	uint8_t					m_localMaxDepth;		/**< Reached max depth in local tree */
	uint8_t 				m_balanceCodim;			/**<Maximum codimension of the entity for 2:1 balancing (1 = 2:1 balance through faces (default);
	 	 	 	 	 	 	 	 	 	 	 	 	 	 2 = 2:1 balance through edges and faces;
	 	 	 	 	 	 	 	 	 	 	 	 	 	 3 = 2:1 balance through nodes, edges and faces)*/
	u32vector 				m_lastGhostBros;		/**<Index of ghost brothers in case of broken family coarsened*/
	u32vector2D				m_connectivity;			/**<Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	u32vector2D				m_ghostsConnectivity;	/**<Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	u32arr3vector			m_nodes;				/**<Local vector of nodes (x,y,z) ordered with Morton Number*/
	u32arr3vector			m_ghostsNodes;			/**<Local vector of ghosts nodes (x,y,z) ordered with Morton Number*/

	uint8_t					m_dim;					/**<Space dimension. Only 2D or 3D space accepted*/
	ClassGlobal				m_global;				/**<Global variables*/

	// =================================================================================== //
	// CONSTRUCTORS
	// =================================================================================== //
	ClassLocalTree(int8_t maxlevel, uint8_t dim);
	~ClassLocalTree();

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
private:
	const ClassOctant&	getFirstDesc() const;
	const ClassOctant&	getLastDesc() const;
	uint32_t 		getSizeGhost() const;
	uint32_t 		getNumOctants() const;
	uint8_t 		getLocalMaxDepth() const;
	int8_t 			getMarker(int32_t idx);
	uint8_t 		getLevel(int32_t idx);
	uint64_t 		computeMorton(int32_t idx);
	uint8_t 		getGhostLevel(int32_t idx);
	uint64_t 		computeGhostMorton(int32_t idx);
	bool 			getBalance(int32_t idx);
	uint8_t 		getBalanceCodim() const;
	void 			setMarker(int32_t idx, int8_t marker);
	void 			setBalance(int32_t idx, bool balance);
	void 			setBalanceCodim(uint8_t b21codim);
	void 			setFirstDesc();
	void 			setLastDesc();

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

	// =================================================================================== //
	// OTHER METHODS
	// =================================================================================== //

	ClassOctant& 		extractOctant(uint32_t idx);
	const ClassOctant&	extractOctant(uint32_t idx) const;
	ClassOctant& 		extractGhostOctant(uint32_t idx);
	const ClassOctant&	extractGhostOctant(uint32_t idx) const;


	bool 		refine(u32vector & mapidx);
	bool 		coarse(u32vector & mapidx);
	bool 		globalRefine(u32vector & mapidx);
	bool 		globalCoarse(u32vector & mapidx);
	void 		checkCoarse(uint64_t lastDescPre, uint64_t firstDescPost, u32vector & mapidx);
	void 		updateLocalMaxDepth();


	void 		findNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours,
								std::vector<bool> & isghost);
	void 		findNeighbours(ClassOctant* oct, uint8_t iface, u32vector & neighbours,
								std::vector<bool> & isghost);
	void 		findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours);
	void 		findEdgeNeighbours(uint32_t idx, uint8_t iedge,
								u32vector & neighbours, std::vector<bool> & isghost);
	void 		findEdgeNeighbours(ClassOctant* oct, uint8_t iedge,
								u32vector & neighbours, std::vector<bool> & isghost);
	void 		findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge,
								u32vector & neighbours);
	void 		findNodeNeighbours(ClassOctant* oct, uint8_t inode,
								u32vector & neighbours, std::vector<bool> & isghost);
	void 		findNodeNeighbours(uint32_t idx, uint8_t inode,
								u32vector & neighbours, std::vector<bool> & isghost);
	void 		findGhostNodeNeighbours(uint32_t idx, uint8_t inode,
								u32vector & neighbours);


	void 		preBalance21(bool internal);
	void 		preBalance21(u32vector& newmodified);
	bool 		localBalance(bool doInterior);
	bool 		localBalanceAll(bool doInterior);

	void 		computeIntersections();

	uint32_t 	findMorton(uint64_t Morton);
	uint32_t 	findGhostMorton(uint64_t Morton);

	void 		computeConnectivity();
	void 		clearConnectivity();
	void 		updateConnectivity();
	void 		computeGhostsConnectivity();
	void 		clearGhostsConnectivity();
	void 		updateGhostsConnectivity();

	// =================================================================================== //

};

#endif /* CLASSLOCALTREE_HPP_ */
