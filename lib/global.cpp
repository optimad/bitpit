/*
 * global.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */
#include "global.hpp"

//Explicit Instantiation
template class global<2>;
template class global<3>;

const global<2> global2D;
const global<3> global3D;
