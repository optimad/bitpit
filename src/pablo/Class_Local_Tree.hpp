#ifndef CLASS_LOCAL_TREE_HPP_
#define CLASS_LOCAL_TREE_HPP_


// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Octant.hpp"
#include "classOctant.hpp"
#include "Class_Intersection.hpp"
#include <math.h>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <map>
#include <iostream>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 * Class_Local_Tree.hpp
 *
 *  \ingroup    PABLO
 *  @{
 *  \date		23/apr/2014
 *	\authors	Edoardo Lombardi
 *	\authors	Marco Cisternino
 *	\version	0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
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
 *	Class Local_Tree is a templated class in dimensional parameter int dim and it accepts only two values: dim=2 and dim=3, obviously for 2D and 3D respectively.
 */
template<int dim>
class Class_Local_Tree{};

#include "Class_Local_Tree_3D.tpp"
#include "Class_Local_Tree_2D.tpp"

/*  @} */

#endif /* CLASS_LOCAL_TREE_HPP_ */
