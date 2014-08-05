#ifndef CLASS_INTERSECTION_HPP_
#define CLASS_INTERSECTION_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "logFunct.hpp"
#include <vector>
#include <string.h>
#include "inlinedFunct.hpp"

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *
 *	\brief Intersections class definition
 *
 *	Intersections are the faces (edges in 2D) or portion of faces shared by two octants. An intersection is defined
 *	by :
 *	- the owner octants, i.e. the octants sharing the intersection, identified by a couple (array[2]) of indices;
 *	- the index of the face, that contains the intersection, of the first owner;
 *	- an identifier of the octant in the couple with higher level of refinement (0/1) [if same level identifier =0];
 *	- a flag stating if an owner is ghost;
 *	- a flag to communicate if the intersection is new after a mesh refinement.
 *
 */
template<int dim>
class Class_Intersection{};

#include "Class_Intersection_3D.tpp"
#include "Class_Intersection_2D.tpp"


#endif /* CLASS_INTERSECTION_HPP_ */
