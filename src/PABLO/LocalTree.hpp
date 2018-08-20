/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BITPIT_PABLO_LOCALTREE_HPP__
#define __BITPIT_PABLO_LOCALTREE_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "tree_constants.hpp"
#include "Octant.hpp"
#include "Intersection.hpp"

namespace bitpit {

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //

typedef std::vector<Octant>		 		octvector;
typedef std::vector<Intersection>	 	intervector;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
/*!
 * LocalTree.hpp
 *
 *	\ingroup	PABLO
 *	\date		15/dec/2015
 *	\authors	Edoardo Lombardi
 *	\authors	Marco Cisternino
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
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
 *	Class LocalTree is built with a dimensional parameter int dim and it accepts only two values: dim=2 and dim=3, for 2D and 3D respectively.
 */
class LocalTree{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class ParaTree;

	// =================================================================================== //
	// TYPEDEFS
	// =================================================================================== //
public:

	/*!Vector of Octants.
	 */
	typedef std::vector<Octant>				 	octvector;

	/*!Vector of Intersections.
	 */
	typedef std::vector<Intersection>	 		intervector;

	/*!Vector of boolean values.
	 */
	typedef std::vector<bool>					bvector;

	/*!Vector of unsigned int (8-bit).
	 */
	typedef std::vector<uint8_t>				u8vector;

	/*!Vector of usigned int (32-bit).
	 */
	typedef std::vector<uint32_t>				u32vector;

	/*!Vector of unsigned int (64-bit).
	 */
	typedef std::vector<uint64_t>				u64vector;

	/*!Vector of three-dimensional arrays of unsigned int (32-bit).
	 */
	typedef std::vector<u32array3>				u32arr3vector;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //

private:
	octvector				m_octants;				/**< Local vector of octants ordered with Morton Number */
	octvector				m_ghosts;				/**< Local vector of ghost octants ordered with Morton Number */
	intervector				m_intersections;		/**< Local vector of intersections */
	u64vector 				m_globalIdxGhosts;		/**< Global index of the ghost octants (size = size_ghosts) */
	uint64_t 				m_firstDescMorton;		/**< Morton number of first (Morton order) most refined octant possible in local partition */
	uint64_t			 	m_lastDescMorton;		/**< Morton number of last (Morton order) most refined octant possible in local partition */
	uint32_t 				m_sizeGhosts;			/**< Size of vector of ghost octants */
	uint32_t 				m_sizeOctants;			/**< Size of vector of local octants */
	uint8_t					m_localMaxDepth;		/**< Reached max depth in local tree */
	uint8_t 				m_balanceCodim;			/**<Maximum codimension of the entity for 2:1 balancing (1 = 2:1 balance through faces (default);
	 	 	 	 	 	 	 	 	 	 	 	 	 	 2 = 2:1 balance through edges and faces;
	 	 	 	 	 	 	 	 	 	 	 	 	 	 3 = 2:1 balance through nodes, edges and faces)*/
	u32vector 				m_lastGhostBros;		/**<Index of ghost brothers in case of broken family coarsened (tail of local octants)*/
	u32vector 				m_firstGhostBros;		/**<Index of ghost brothers in case of broken family coarsened (head of local octants)*/
	u32vector2D				m_connectivity;			/**<Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	u32vector2D				m_ghostsConnectivity;	/**<Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 The nodes are stored as index of vector nodes*/
	u32arr3vector			m_nodes;				/**<Local vector of nodes (x,y,z) ordered with Morton Number*/

	uint8_t					m_dim;					/**<Space dimension. Only 2D or 3D space accepted*/
	TreeConstants			m_global;				/**<Global variables*/
	bvector 				m_periodic;				/**<Boolvector: i-th element is true if the i-th boundary face is a periodic interface.*/


	// =================================================================================== //
	// CONSTRUCTORS
	// =================================================================================== //
	LocalTree();
	LocalTree(uint8_t dim);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
private:
	uint64_t		getFirstDescMorton() const;
	uint64_t		getLastDescMorton() const;
	uint32_t 		getNumGhosts() const;
	uint32_t 		getNumOctants() const;
	uint8_t 		getLocalMaxDepth() const;
	int8_t 			getMarker(int32_t idx) const;
	uint8_t 		getLevel(int32_t idx) const;
	uint64_t 		computeMorton(int32_t idx) const;
	uint64_t 		computeNodeMorton(int32_t idx, uint8_t inode) const;
	uint8_t 		getGhostLevel(int32_t idx) const;
	uint64_t 		computeGhostMorton(int32_t idx) const;
	uint64_t 		computeGhostNodeMorton(int32_t idx, uint8_t inode) const;
	bool 			getBalance(int32_t idx) const;
	uint8_t 		getBalanceCodim() const;
	void 			setMarker(int32_t idx, int8_t marker);
	void 			setBalance(int32_t idx, bool balance);
	void 			setBalanceCodim(uint8_t b21codim);
	void 			setFirstDescMorton();
	void 			setLastDescMorton();
	void 			setPeriodic(bvector & periodic);

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

	// =================================================================================== //
	// OTHER METHODS
	// =================================================================================== //

	void	initialize();
	void	initialize(uint8_t dim);
	void	reset(bool createRoot);

	Octant& 		extractOctant(uint32_t idx);
	const Octant&	extractOctant(uint32_t idx) const;
	Octant& 		extractGhostOctant(uint32_t idx);
	const Octant&	extractGhostOctant(uint32_t idx) const;


	bool 		refine(u32vector & mapidx);
	bool 		coarse(u32vector & mapidx);
	bool 		globalRefine(u32vector & mapidx);
	bool 		globalCoarse(u32vector & mapidx);
	void 		checkCoarse(uint64_t partLastDesc, u32vector & mapidx);
	void 		updateLocalMaxDepth();

    
    void        findNeighbours(const Octant* oct, bool haveIidx, uint32_t idx, uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findEdgeNeighbours(const Octant* oct, bool haveIidx, uint32_t idx, uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findNodeNeighbours(const Octant* oct, bool haveIidx, uint32_t idx, uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;

    void        findNeighbours(uint32_t idx, bool amIghost, uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findEdgeNeighbours(uint32_t idx, bool amIghost, uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findNodeNeighbours(uint32_t idx, bool amIghost, uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findEdgeNeighbours(const Octant* oct, uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;
    void        findNodeNeighbours(const Octant* oct, uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal) const;

	void 		findPeriodicNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost) const;
	void 		findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours, bvector & isghost) const;
	void 		findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours, bvector & isghost) const;
	void 		findGhostNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours, bvector & isghost) const;
	void 		findGhostPeriodicNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours) const;


    void        findNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours, bvector & isghost) const;
    void        findEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours, bvector & isghost) const;
    void        findNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours, bvector & isghost) const;
    void        findNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost) const;
    void        findEdgeNeighbours(const Octant* oct, uint8_t iedge, u32vector & neighbours, bvector & isghost) const;
    void        findNodeNeighbours(const Octant* oct, uint8_t inode, u32vector & neighbours, bvector & isghost) const;
	void 		findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours) const;
	void 		findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours) const;
	void 		findGhostNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours) const;


	void 		preBalance21(bool internal);
	void 		preBalance21(u32vector& newmodified);
	bool 		localBalance(bool doNew, bool doInterior);

	void 		computeIntersections();

	uint32_t 	findMorton(uint64_t Morton) const;
	uint32_t 	findGhostMorton(uint64_t Morton) const;

	void 		computeConnectivity();
	void 		clearConnectivity();
	void 		updateConnectivity();

	// =================================================================================== //

};

}

#endif /* __BITPIT_PABLO_LOCALTREE_HPP__ */
