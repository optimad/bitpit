#ifndef CLASSOCTANT_HPP_
#define CLASSOCTANT_HPP_

// INCLUDES                                                                            //
#include "inlinedFunct.hpp"
#include <vector>
#include <bitset>
#include <array>

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //
typedef std::vector<uint8_t>				u8vector;
typedef std::vector<uint32_t>				u32vector;
typedef std::vector<uint64_t>				u64vector;
typedef std::vector<double>					dvector;
typedef std::array<double, 3>				darray3;
typedef std::array<int8_t, 3>				i8array3;
typedef std::array<uint32_t, 3>				u32array3;
typedef std::vector<std::vector<uint32_t> >	u32vector2D;
typedef std::vector<std::vector<uint64_t> >	u64vector2D;
typedef std::vector<u32array3>				u32arr3vector;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			10/dec/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Octant class definition
 *
 *	Octants are the grid elements of PABLO. In the logical domain octants are
 *	squares or cubes, depending on dimension, with size function of their level.
 *	Each octant has nodes and faces ordered with Z-order.
 *
 *	The main feature of each octant are:
 *	- x,y,z        : coordinates of the node 0 of the octant;
 *	- Morton index : classical Morton index defined anly by the coordinates
 *	(info about level used additionally for equality operator);
 *	- marker       : refinement marker can assume negative, positive or zero values, wich mean
 *	a coarsening, refinement and none adaptation respectively;
 *	- level        : octant level in the octree, zero for the first upper level.
 *	- balance      : flag to fix if the octant has to 2:1 balanced with respect
 *	to its face neighbours.
 *
 */
class ClassOctant{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class ClassLocalTree;
	friend class ClassParaTree;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	uint32_t		  	m_x;			/**< Coordinate x */
	uint32_t		  	m_y;			/**< Coordinate y */
	uint32_t		  	m_z;			/**< Coordinate z (2D case = 0)*/
	uint8_t   			m_level;		/**< Refinement level (0=root) */
	int8_t    			m_marker;		/**< Set for Refinement(m>0) or Coarsening(m<0) |m|-times */
	std::bitset<17> 	m_info;			/**< -Info[0..#faces]: true if 0..#faces face is a boundary face [bound] \n
										-Info[6..#faces+5]: true if 0..#faces face is a process boundary face [pbound] \n
										-Info[12/13]: true if octant is new after refinement/coarsening \n
										-Info[14]   : true if balancing is not required for this octant \n
										-Info[15]   : Aux
										-Info[16]   : true if octant is a scary ghost */
	uint8_t				m_dim;			/**< Dimension of octant (2D/3D) */

	//TODO add bitset for edge & node

	// =================================================================================== //
	// STATIC MEMBERS
	// =================================================================================== //

	static constexpr int 	sm_CoeffFaceCenter[6][3] = {{0,1,1},{2,1,1},{1,0,1},
														{1,2,1},{1,1,0},{1,1,2}}; /**< Static member for internal use. */

	static constexpr int 	sm_CoeffEdgeCenter[12][3] =	{{0,1,0},{2,1,0},{1,0,0},
														{1,2,0},{0,0,1},{2,0,1},
														{0,2,1},{2,2,1},{0,1,2},
														{2,1,2},{1,0,2},{1,2,2}};  /**< Static member for internal use. */

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:
	ClassOctant();
	ClassOctant(const ClassOctant &octant);
private:
	ClassOctant(uint8_t dim);
	ClassOctant(uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z = 0);
	ClassOctant(bool bound, uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z = 0);
	bool operator ==(const ClassOctant & oct2);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
	uint32_t	getDim() const;
	u32array3	getCoordinates() const;
	uint32_t	getX() const;
	uint32_t	getY() const;
	uint32_t	getZ() const;
	u32array3	getCoord();
	uint8_t		getLevel() const;
	int8_t		getMarker() const;
	bool		getBound(uint8_t face) const;
	bool		getBound() const;
	void		setBound(uint8_t face);
	bool		getPbound(uint8_t face) const;
	bool		getPbound() const;
	bool		getIsNewR() const;
	bool		getIsNewC() const;
	bool		getIsGhost() const;
	bool		getNotBalance() const;
	bool		getBalance() const;
	void		setMarker(int8_t marker);
	void		setBalance(bool balance);
	void		setLevel(uint8_t level);
	void 		setPbound(uint8_t face, bool flag);

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //
	uint32_t	getSize(int8_t maxlevel) const;
	uint64_t	getArea(int8_t maxlevel) const;
	uint64_t	getVolume(int8_t maxlevel) const;
	darray3		getCenter(int8_t maxlevel) const;
	darray3		getFaceCenter(uint8_t iface, int8_t maxlevel) const;
	darray3		getEdgeCenter(uint8_t iedge, int8_t maxlevel) const;
	void		getNodes(u32arr3vector & nodes, int8_t maxlevel) const;
	void		getNode(u32array3 & node, uint8_t inode, int8_t maxlevel) const;
	u32array3		getNode(uint8_t inode, int8_t maxlevel) const;
	void		getNormal(uint8_t & iface, i8array3 & normal, int8_t (&normals)[6][3], int8_t maxlevel) const;
	uint64_t	computeMorton() const;
	uint64_t	computeMorton();

	// =================================================================================== //
	// OTHER METHODS												    			   //
	// =================================================================================== //
	ClassOctant	buildLastDesc(int8_t & maxlevel);
	ClassOctant	buildFather(int8_t & maxlevel);
	std::vector< ClassOctant >	buildChildren(int8_t & maxlevel);
	std::vector<uint64_t> 		computeHalfSizeMorton(uint8_t iface, uint32_t & sizehf,
			int8_t & maxlevel);
	std::vector<uint64_t>		computeMinSizeMorton(uint8_t iface, const uint8_t & maxdepth,
			uint32_t & sizem, int8_t & maxlevel);
	std::vector<uint64_t> 		computeVirtualMorton(uint8_t iface, const uint8_t & maxdepth,
			uint32_t & sizeneigh, int8_t & maxlevel);
	std::vector<uint64_t> 		computeEdgeHalfSizeMorton(uint8_t iedge, uint32_t & sizehf,
			int8_t & maxlevel, uint8_t (&edgeface)[12][2]);
	std::vector<uint64_t> 		computeEdgeMinSizeMorton(uint8_t iedge, const uint8_t & maxdepth,
			uint32_t & sizem, int8_t & maxlevel, uint8_t (&edgeface)[12][2]);
	std::vector<uint64_t>		computeEdgeVirtualMorton(uint8_t iedge, const uint8_t & maxdepth,
			uint32_t & sizeneigh, uint8_t balance_codim, int8_t & maxlevel, uint8_t (&edgeface)[12][2]);
	uint64_t 		computeNodeMinSizeMorton(uint8_t inode, const uint8_t & maxdepth,
			uint32_t & sizehf, int8_t & maxlevel, uint8_t (&nodeface)[8][3]);
	uint64_t 		computeNodeVirtualMorton(uint8_t inode, const uint8_t & maxdepth,
			uint32_t & sizeneigh, int8_t & maxlevel, uint8_t (&nodeface)[8][3]);

};


/*  @} */

#endif /* CLASSOCTANT_HPP_ */
