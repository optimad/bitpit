#ifndef CLASSMAP_HPP_
#define CLASSMAP_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classGlobal.hpp"
#include <vector>
#include <iostream>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //
typedef vector<double>				dvector;
typedef vector<dvector>				dvector2D;
typedef vector<uint32_t>			u32vector;
typedef vector<u32vector>			u32vector2D;
typedef vector<uint64_t>			u64vector;
typedef vector<u64vector>			u64vector2D;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			17/dec/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Transformation Mapper
 *
 *	Definition of the transformation from the logical domain to the physical domain.
 *	It contains a default (temporary) implementation of a scaling and translation mapper
 *	of logical octree.
 *	classMap has to be implemented and customized by the user for different applications.
 */
class classMap{

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
public:
	double 		X0;						/**<Coordinate X of the origin of the octree in the physical domain*/
	double 		Y0;						/**<Coordinate Y of the origin of the octree in the physical domain*/
	double 		Z0;						/**<Coordinate Z of the origin of the octree in the physical domain*/
	double 		L;						/**<Side length of octree in the physical domain*/
	uint8_t		dim;					/**<Space Dimension*/
	uint8_t		nnodes;					/**<Number of nodes*/
	uint8_t		nnodesperface;			/**<Number of nodes for each face*/
	uint32_t	max_length;				/**< Length of the logical domain */

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:
	/*!Default constructor. Origin of octree in physical domain in (0,0,0)
	 * side length 1 and 2D space.
	 */
	classMap(int8_t maxlevel, uint8_t dim_);

	/*!Customized constructor with origin of octree in physical
	 * domain side length provided by the user.
	 * \param[in] X Coordinate X of the origin.
	 * \param[in] Y Coordinate Y of the origin.
	 * \param[in] Z Coordinate Z of the origin.
	 * \param[in] LL Side length of domain.
	 * \param[in] dim Space dimension 2D/3D (default=2).
	 */
	classMap(double & X, double & Y, double & Z, double & LL, int8_t maxlevel, uint8_t dim_);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	/*! Transformation of coordinate X.
	 * \param[in] X Coordinate X from logical domain.
	 * \return Coordinate X in physical domain.
	 */
	double mapX(uint32_t const & X);

	/*! Transformation of coordinate Y.
	 * \param[in] Y Coordinate Y from logical domain.
	 * \return Coordinate Y in physical domain.
	 */
	double mapY(uint32_t const & Y);

	/*! Transformation of coordinate Z.
	 * \param[in] Z Coordinate Z from logical domain.
	 * \return Coordinate Z in physical domain.
	 */
	double mapZ(uint32_t const & Z);

	/*! Transformation of coordinate X.
	 * \param[in] X Coordinate X from physical domain.
	 * \return Coordinate X in logical domain.
	 */
	uint32_t mapX(double const & X);

	/*! Transformation of coordinate Y.
	 * \param[in] Y Coordinate Y from physical domain.
	 * \return Coordinate Y in logical domain.
	 */
	uint32_t mapY(double const & Y);

	/*! Transformation of coordinate Z.
	 * \param[in] Z Coordinate Z from physical domain.
	 * \return Coordinate Z in logical domain.
	 */
	uint32_t mapZ(double const & Z);

	/*! Transformation of size of an octant.
	 * \param[in] size Size of octant from logical domain.
	 * \return Size of octant in physical domain.
	 */
	double mapSize(uint32_t const & size);

	/*! Transformation of area of an octant.
	 * \param[in] area Area of octant from logical domain.
	 * \return Area of octant in physical domain.
	 */
	double mapArea(uint64_t const & area);

	/*! Transformation of volume of an octant.
	 * \param[in] volume Volume of octant from logical domain.
	 * \return Coordinate Volume of octant in physical domain.
	 */
	double mapVolume(uint64_t const & volume);

	/*! Transformation of coordinates of center of an octant.
	 * \param[in] center Pointer to coordinates of center from logical domain.
	 * \param[out] mapcenter Coordinates of center in physical domain.
	 */
	void mapCenter(double* & center,
			dvector & mapcenter);

	/*! Transformation of coordinates of center of an octant.
	 * \param[in] center Vector of coordinates of center from logical domain.
	 * \param[out] mapcenter Coordinates of center in physical domain.
	 */
	void mapCenter(dvector & center,
			dvector & mapcenter);

	/*! Transformation of coordinates of nodes of an octant.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodes(uint32_t (*nodes)[3],
			dvector2D & mapnodes);

	/*! Transformation of coordinates of nodes of an octant.
	 * \param[in] nodes Vector of coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodes(u32vector2D nodes,
			dvector2D & mapnodes);

	/*! Transformation of coordinates of a node of an octant.
	 * \param[in] node Coordinates of  the node from logical domain.
	 * \param[out] mapnodes Coordinates of the node in physical domain.
	 */
	void mapNode(u32vector & node,
			dvector & mapnode);

	/*! Transformation of coordinates of nodes of an intersection.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodesIntersection(uint32_t (*nodes)[3],
			dvector2D & mapnodes);

	/*! Transformation of coordinates of nodes of an intersection.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodesIntersection(u32vector2D nodes,
			dvector2D & mapnodes);

	/*! Transformation of components of normal of an intersection.
	 * \param[in] nodes Pointer to components of normal from logical domain.
	 * \param[out] mapnodes components of normal in physical domain.
	 */
	void mapNormals(vector<int8_t> normal,
			dvector & mapnormal);

};

/* @} */

#endif /* CLASSMAP_HPP_ */
