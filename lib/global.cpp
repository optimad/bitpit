/*
 * global.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //

#include "global.hpp"

const uint32_t	max_length	=	uint32_t(pow(2.0,MAX_LEVEL));
const uint8_t	nchildren	=	uint8_t(pow(2.0,DIM));
const uint8_t	nface		=	uint8_t (2.0*DIM);
const uint8_t	nnodes		=	uint8_t(pow(2.0,DIM));
const uint8_t   octantBytes =   uint8_t(sizeof(uint32_t)*DIM + sizeof(uint8_t) + sizeof(int8_t) + (2*nface+4)*sizeof(bool));
const uint8_t   markerBytes =   uint8_t(sizeof(int8_t));
const uint8_t   levelBytes	=   uint8_t(sizeof(uint8_t));
const uint8_t 	oppface[2*DIM] = {1,0,3,2,5,4};
