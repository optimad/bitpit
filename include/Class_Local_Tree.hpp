/*
 * Class_Local_Tree.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#ifndef CLASS_LOCAL_TREE_HPP_
#define CLASS_LOCAL_TREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "global.hpp"
#include "Class_Octant.hpp"
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

template<int dim>
class Class_Local_Tree{};

#include "Class_Local_Tree_3D.tpp"
//#include "Class_Local_Tree_2D.tpp"

#endif /* CLASS_LOCAL_TREE_HPP_ */
