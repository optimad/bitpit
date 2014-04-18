/*
 * octant.hpp
 *
 *  Created on: Feb 10, 2014
 *      Author: Edoardo Lombardi
 */

#ifndef OCTANT_HPP_
#define OCTANT_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "global.hpp"
#include "logFunct.hpp"
#include <vector>
#include <string.h>
#include "inlinedFunct.hpp"
//#include "Class_Local_Tree.hpp"

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
template<int dim>
class Class_Octant{};

#include "Class_Octant_3D.tpp"
#include "Class_Octant_2D.tpp"

#endif /* OCTREE_HPP_ */
