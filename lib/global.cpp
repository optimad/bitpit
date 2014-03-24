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
const uint8_t   boolBytes	=   uint8_t(sizeof(bool));
const uint8_t 	oppface[2*DIM] = {1,0,3,2,5,4};
#if DIM == 3
const uint8_t   edgeface[12][2] = {{0,4},{1,4},{2,4},{3,4},{0,2},{1,2},{0,3},{1,3},{0,5},{1,5},{2,5},{3,5}};
const uint8_t   nodeface[8][3] = {{0,2,4},{1,2,4},{0,3,4},{1,3,4},{0,2,5},{1,2,5},{0,3,5},{1,3,5}};
const uint8_t   normals[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
#else
#endif
