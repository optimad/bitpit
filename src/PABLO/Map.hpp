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

#ifndef __BITPIT_PABLO_MAP_HPP__
#define __BITPIT_PABLO_MAP_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include <vector>
#include <iostream>
#include <array>

namespace bitpit {

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
 *	\ingroup		PABLO
 *	\date			17/dec/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *
 *	\brief Transformation Mapper
 *
 *	Definition of the transformation from the logical domain to the physical reference domain.
 *	It contains a default implementation of a scaling and translation mapper
 *	of logical octree in the reference domain with origin in (0,0,0) and size 1.
 *	Map has to be implemented and customized by the user for different applications as a derived
 *	class of ParaTree (see PabloUniform for a basic example).
 */
class Map{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //

	friend class ParaTree;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	darray3 	m_origin;				/**<Coordinate X,Y,Z of the origin of the octree in the physical domain*/
	double 		m_L;					/**<Side length of octree in the physical domain*/
	uint8_t		m_dim;					/**<Space Dimension*/
	uint8_t		m_nnodes;				/**<Number of nodes*/
	uint8_t		m_nnodesPerFace;		/**<Number of nodes for each face*/
	uint32_t	m_maxLength;			/**< Length of the logical domain */
	double		m_maxLength_1;			/**< 1/Length of the logical domain */

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
	Map();
	Map(uint8_t dim);
//	Map(double & X, double & Y, double & Z, double & LL, uint8_t dim);

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	void initialize();
	void initialize(uint8_t dim);

	darray3 mapCoordinates(u32array3 const & X) const;
	double mapX(uint32_t const & X) const;
	double mapY(uint32_t const & Y) const;
	double mapZ(uint32_t const & Z) const;
	u32array3 mapCoordinates(darray3 const & X) const;
	uint32_t mapX(double const & X) const;
	uint32_t mapY(double const & Y) const;
	uint32_t mapZ(double const & Z) const;
	double mapSize(uint32_t const & size) const;
	double mapArea(uint64_t const & area) const;
	double mapVolume(uint64_t const & volume) const;
	void mapCenter(double* & center, darray3 & mapcenter) const;
	void mapCenter(darray3 & center, darray3 & mapcenter) const;
	void mapNodes(uint32_t (*nodes)[3], darr3vector & mapnodes) const;
	void mapNodes(u32arr3vector nodes, darr3vector & mapnodes) const;
	void mapNode(u32array3 & node, darray3 & mapnode) const;
	void mapNodesIntersection(uint32_t (*nodes)[3], darr3vector & mapnodes) const;
	void mapNodesIntersection(u32arr3vector nodes, darr3vector & mapnodes) const;
	void mapNormals(i8array3 normal, darray3 & mapnormal) const;

};

}

#endif /* __BITPIT_PABLO_MAP_HPP__ */
