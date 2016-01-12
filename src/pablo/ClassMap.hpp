#ifndef CLASSMAP_HPP_
#define CLASSMAP_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include <vector>
#include <iostream>
#include <array>

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //
typedef std::vector<double>			dvector;
typedef std::vector<dvector>		dvector2D;
typedef std::vector<uint32_t>		u32vector;
typedef std::vector<u32vector>		u32vector2D;
typedef std::vector<uint64_t>		u64vector;
typedef std::vector<u64vector>		u64vector2D;
typedef std::array<double, 3>		darray3;
typedef std::array<int8_t, 3>		i8array3;
typedef std::array<uint32_t, 3>		u32array3;
typedef std::vector<u32array3>		u32arr3vector;
typedef std::vector<darray3>		darr3vector;

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
 *	ClassMap has to be implemented and customized by the user for different applications.
 */
class ClassMap{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class ClassParaTree;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
//	double 		X0;						/**<Coordinate X of the origin of the octree in the physical domain*/
//	double 		Y0;						/**<Coordinate Y of the origin of the octree in the physical domain*/
//	double 		Z0;						/**<Coordinate Z of the origin of the octree in the physical domain*/
	darray3 	m_origin;				/**<Coordinate X,Y,Z of the origin of the octree in the physical domain*/
	double 		m_L;					/**<Side length of octree in the physical domain*/
	uint8_t		m_dim;					/**<Space Dimension*/
	uint8_t		m_nnodes;				/**<Number of nodes*/
	uint8_t		m_nnodesPerFace;		/**<Number of nodes for each face*/
	uint32_t	m_maxLength;			/**< Length of the logical domain */

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:
	/*!Default constructor. Origin of octree in physical domain in (0,0,0)
	 * side length 1 and 2D space.
	 */
	ClassMap(int8_t maxlevel, uint8_t dim_);

	/*!Customized constructor with origin of octree in physical
	 * domain side length provided by the user.
	 * \param[in] X Coordinate X of the origin.
	 * \param[in] Y Coordinate Y of the origin.
	 * \param[in] Z Coordinate Z of the origin.
	 * \param[in] LL Side length of domain.
	 * \param[in] dim Space dimension 2D/3D (default=2).
	 */
	ClassMap(double & X, double & Y, double & Z, double & LL, int8_t maxlevel, uint8_t dim_);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	darray3 mapCoordinates(u32array3 const & X);

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

	u32array3 mapCoordinates(darray3 const & X);

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
			darray3 & mapcenter);

	/*! Transformation of coordinates of center of an octant.
	 * \param[in] center Array of coordinates of center from logical domain.
	 * \param[out] mapcenter Coordinates of center in physical domain.
	 */
	void mapCenter(darray3 & center,
			darray3 & mapcenter);

	/*! Transformation of coordinates of nodes of an octant.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodes(uint32_t (*nodes)[3],
			darr3vector & mapnodes);

	/*! Transformation of coordinates of nodes of an octant.
	 * \param[in] nodes Vector of coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodes(u32arr3vector nodes,
			darr3vector & mapnodes);

	/*! Transformation of coordinates of a node of an octant.
	 * \param[in] node Coordinates of  the node from logical domain.
	 * \param[out] mapnodes Coordinates of the node in physical domain.
	 */
	void mapNode(u32array3 & node,
			darray3 & mapnode);

	/*! Transformation of coordinates of nodes of an intersection.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodesIntersection(uint32_t (*nodes)[3],
			darr3vector & mapnodes);

	/*! Transformation of coordinates of nodes of an intersection.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodesIntersection(u32arr3vector nodes,
			darr3vector & mapnodes);

	/*! Transformation of components of normal of an intersection.
	 * \param[in] nodes Pointer to components of normal from logical domain.
	 * \param[out] mapnodes components of normal in physical domain.
	 */
	void mapNormals(i8array3 normal,
			darray3 & mapnormal);

};

/* @} */

#endif /* CLASSMAP_HPP_ */
