/*
 * global.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>

// =================================================================================== //
// PARAMETERS DEFINITION                                                               //
// =================================================================================== //

extern const uint32_t max_length;
extern const uint8_t  nchildren;
extern const uint8_t  nface;
extern const uint8_t  nnodes;
extern const uint8_t  octantBytes;
extern const uint8_t  markerBytes;
extern const uint8_t  levelBytes;
extern const uint8_t  boolBytes;
extern const uint8_t  oppface[2*DIM];
#if DIM == 3
extern const uint8_t  edgeface[12][2];
extern const uint8_t  nodeface[8][3];
extern const uint8_t  normals[6][3];
#else
#endif



#endif /* GLOBAL_HPP_ */
