/*
 * Class_Global.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: Marco
 */

#ifndef CLASS_GLOBAL_HPP_
#define CLASS_GLOBAL_HPP_

#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>

template <int dim>
class Class_Global;

#include "Class_Global_2D.tpp"
#include "Class_Global_3D.tpp"

extern const Class_Global<2> global2D;
extern const Class_Global<3> global3D;


#endif /* CLASS_GLOBAL_HPP_ */
