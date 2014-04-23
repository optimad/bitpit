/*
 * global.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */
#include "global.hpp"

//Explicit Instantiation
template class Global<2>;
template class Global<3>;

const Global<2> global2D;
const Global<3> global3D;
