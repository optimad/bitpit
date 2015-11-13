#ifndef CLASS_GLOBAL_HPP_
#define CLASS_GLOBAL_HPP_

#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>

/*!
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Global variables used in PABLO
 *
 *	Global variables are used in PABLO everywhere and they are public, i.e. each
 *	global variable can be used as constant by external codes.
 *
 *	Class Class_Global is a templated class in dimensional parameter int dim and it accepts only two values: dim=2 and dim=3, obviously for 2D and 3D respectively.
 *
 */
template <int dim>
class Class_Global;

#include "Class_Global_2D.tpp"
#include "Class_Global_3D.tpp"

extern const Class_Global<2> global2D;
extern const Class_Global<3> global3D;


#endif /* CLASS_GLOBAL_HPP_ */
