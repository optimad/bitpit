/*
 * Class_Global.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */
#include "Class_Global.hpp"

//Explicit Instantiation
template class Class_Global<2>;
template class Class_Global<3>;

const Class_Global<2> global2D;
const Class_Global<3> global3D;
