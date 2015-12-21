#ifndef CLASSOCTANT_HPP_
#define CLASSOCTANT_HPP_

// INCLUDES                                                                            //
#include "classGlobal.hpp"
#include "inlinedFunct.hpp"
#include "Class_Array.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include <bitset>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //
typedef vector<uint8_t>				u8vector;
typedef vector<uint32_t>			u32vector;
typedef vector<uint64_t>			u64vector;
typedef vector<double>				dvector;
typedef array<double, 3>			darray3;
typedef array<uint32_t, 3>			u32array3;
typedef vector<vector<uint32_t>	>	u32vector2D;
typedef vector<vector<uint64_t>	>	u64vector2D;
typedef vector<u32array3>			u32arr3vector;
typedef classGlobal					CG;

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
 *	- Morton index : classical Morton index defined anly by the coordinates (info about level used additionally for equality operator);
 *	- marker       : refinement marker can assume negative, positive or zero values, wich mean
 *	a coarsening, refinement and none adaptation respectively;
 *	- level        : octant level in the octree, zero for the first upper level.
 *	- balance      : flag to fix if the octant has to 2:1 balanced with respect to its face neighbours.
 *
 */
class classOctant{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class classLocalTree;
	friend class classParaTree;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	uint32_t  	x;				/**< Coordinate x */
	uint32_t  	y;				/**< Coordinate y */
	uint32_t  	z;				/**< Coordinate z (2D case = 0)*/
	uint8_t   	level;			/**< Refinement level (0=root) */
	int8_t    	marker;			/**< Set for Refinement(m>0) or Coarsening(m<0) |m|-times */
	bitset<17>  info;			/**< -Info[0..#faces]: true if 0..#faces face is a boundary face [bound] \n
								-Info[6..#faces+5]: true if 0..#faces face is a process boundary face [pbound] \n
								-Info[12/13]: true if octant is new after refinement/coarsening \n
								-Info[14]   : true if balancing is not required for this octant \n
								-Info[15]   : Aux
								-Info[16]   : true if octant is a scary ghost */
	uint8_t		dim;			/**< Dimension of octant (2D/3D) */

	//TODO add bitset for edge & node

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:
	classOctant();

	classOctant(uint8_t dim);

	classOctant(uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z = 0);

	classOctant(bool bound, uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z = 0);

	classOctant(const classOctant &octant);

	/*! Check if two octants are equal (no check on info)
	 */
	bool operator ==(const classOctant & oct2);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
public:
	/*! Get the dimension of an octant, 2 or 3 if 2D or 3D octant.
	 * \return Dimension of octant.
	 */
	uint32_t	getDim() const;

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate X of node 0.
	 */
	uint32_t	getX() const;

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate Y of node 0.
	 */
	uint32_t	getY() const;

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate Z of node 0.
	 */
	uint32_t	getZ() const;

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinates of node 0.
	 */
	u32array3	getCoord();

	/*! Get the level of an octant.
	 * \return Level of octant.
	 */
	uint8_t		getLevel() const;

	/*! Get the refinement marker of an octant.
	 * \return Marker of octant.
	 */
	int8_t		getMarker() const;

	/*! Get the bound flag on an octant face.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a boundary face.
	 */
	bool		getBound(uint8_t face) const;

private:
	void		setBound(uint8_t face);

public:
	/*! Get the pbound flag on an octant face.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a process boundary face.
	 */
	bool		getPbound(uint8_t face) const;

	/*! Get if the octant is new after a refinement.
	 * \return true if the the octant is new after a refinement.
	 */
	bool		getIsNewR() const;

	/*! Get if the octant is new after a coarsening.
	 * \return true if the the octant is new after a coarsening.
	 */
	bool		getIsNewC() const;

	/*! Get if the octant is a scary ghost octant.
	 * \return true if the octant is a ghost octant.
	 */
	bool		getIsGhost() const;

	/*! Get if the octant is a balancing-blocked octant.
	 * \return false if the octant has to be balanced.
	 */
	bool		getNotBalance() const;

	/*! Get if the octant has to be balanced.
	 * \return true if the octant has to be balanced.
	 */
	bool		getBalance() const;

	/*! Set the refinement marker of an octant.
	 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
	 */
	void		setMarker(int8_t marker);

	/*! Set the balancing condition of an octant.
	 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
	 */
	void		setBalance(bool balance);

private:
	void		setLevel(uint8_t level);

	void 		setPbound(uint8_t face, bool flag);

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

public:

	/*! Get the size of an octant in logical domain, i.e. the side length.
	 * \return Size of octant.
	 */
	uint32_t	getSize(int8_t maxlevel) const;

	/*! Get the area of an octant in logical domain .
	 * \return Area of octant.
	 */
	uint64_t	getArea(int8_t maxlevel) const;

	/*! Get the volume of an octant in logical domain.
	 * \return Volume of octant.
	 */
	uint64_t	getVolume(int8_t maxlevel) const;


	/*! Get the coordinates of the center of an octant in logical domain.
	 * \return Array[3] with the coordinates of the center of octant.
	 */
	vector<double>	getCenter(int8_t maxlevel) const;

	/*! Get the coordinates of the center of a face of an octant in logical domain.
	 * \return Array[3] with the coordinates of the center of the octant face.
	 */
	vector<double>	getFaceCenter(uint8_t iface, int8_t maxlevel) const;

	/*! Get the coordinates of the center of a edge of an octant in logical domain.
	 * \return Array[3] with the coordinates of the center of the octant edge.
	 */
	vector<double>	getEdgeCenter(uint8_t iedge, int8_t maxlevel) const;

	/*! Get the coordinates of the nodes of an octant in logical domain.
	 * \param[out] nodes Vector[nnodes][3] with the coordinates of the nodes of octant.
	 */
	void		getNodes(u32vector2D & nodes, int8_t maxlevel) const;

	/*! Get the coordinates of a nodes of an octant in logical domain.
	 * \param[in] inode Local index of the node
	 * \param[out] node dim-vector with the logical coordinates of the node of the octant.
	 */
	void		getNode(u32vector & node, uint8_t inode, int8_t maxlevel) const;

	/*! Get the coordinates of a nodes of an octant in logical domain.
	 * \param[in] inode Local index of the node
	 * \param[out] node dim-vector with the logical coordinates of the node of the octant.
	 */
	u32vector		getNode(uint8_t inode, int8_t maxlevel) const;

	/*! Get the normal of a face of an octant in logical domain.
	 * \param[in] iface Index of the face for normal computing.
	 * \param[out] normal Vector[3] with components (with z=0) of the normal of face.
	 */
	void		getNormal(uint8_t & iface, vector<int8_t> & normal, int8_t (&normals)[6][3], int8_t maxlevel) const;

	/** Compute the Morton index of the octant (without level).
	 * \return morton Morton index of the octant.
	 */
	uint64_t	computeMorton() const;

	/** Compute the Morton index of the octant (without level).
	 * \return morton Morton index of the octant.
	 */
	uint64_t	computeMorton();

	// =================================================================================== //
	// OTHER METHODS												    			   //
	// =================================================================================== //

private:

	classOctant	buildLastDesc(int8_t & maxlevel);

	classOctant	buildFather(int8_t & maxlevel);

	/** Builds children of octant.
	 *   \return Ordered (by Z-index) vector of children[nchildren] (info update)
	 */
	vector< classOctant >	buildChildren(int8_t & maxlevel);

	vector<uint64_t> 		computeHalfSizeMorton(uint8_t iface, uint32_t & sizehf, int8_t & maxlevel);

	vector<uint64_t>		computeMinSizeMorton(uint8_t iface, const uint8_t & maxdepth,
			uint32_t & sizem, int8_t & maxlevel);

	vector<uint64_t> 		computeVirtualMorton(uint8_t iface, const uint8_t & maxdepth,
			uint32_t & sizeneigh, int8_t & maxlevel);

	vector<uint64_t> 		computeEdgeHalfSizeMorton(uint8_t iedge, uint32_t & sizehf,
			int8_t & maxlevel, uint8_t (&edgeface)[12][2]);

	vector<uint64_t> 		computeEdgeMinSizeMorton(uint8_t iedge, const uint8_t & maxdepth,
			uint32_t & sizem, int8_t & maxlevel, uint8_t (&edgeface)[12][2]);

	vector<uint64_t>		computeEdgeVirtualMorton(uint8_t iedge, const uint8_t & maxdepth,
			uint32_t & sizeneigh, uint8_t balance_codim, int8_t & maxlevel, uint8_t (&edgeface)[12][2]);

	uint64_t 		computeNodeHalfSizeMorton(uint8_t inode, uint32_t & sizehf,
			int8_t & maxlevel, uint8_t (&nodeface)[8][3]);

	uint64_t 		computeNodeMinSizeMorton(uint8_t inode, const uint8_t & maxdepth,
			uint32_t & sizehf, int8_t & maxlevel, uint8_t (&nodeface)[8][3]);

	uint64_t 		computeNodeVirtualMorton(uint8_t inode, const uint8_t & maxdepth,
			uint32_t & sizeneigh, int8_t & maxlevel, uint8_t (&nodeface)[8][3]);

};


/*  @} */

#endif /* CLASSOCTANT_HPP_ */
