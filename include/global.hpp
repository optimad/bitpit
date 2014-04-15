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

template <int dim>
class global{};

template<>
class global<2>
{
public:
	global() : max_length (uint32_t(pow(2.0,MAX_LEVEL_2D))),
	nchildren(4),
	nfaces(4),
	nnodes(4),
	octantBytes(uint8_t(sizeof(uint32_t)*2 + sizeof(uint8_t) + sizeof(int8_t) + (12)*sizeof(bool))),
	markerBytes(sizeof(int8_t)),
	levelBytes(sizeof(uint8_t)),
	boolBytes(sizeof(bool)),
	oppface{1,0,3,2},
	nodeface{{0,2},{1,2},{0,3},{1,3}},
	facenode{{0,2},{1,3},{0,1},{2,3}},
	normals{{-1,0},{1,0},{0,-1},{0,1}}
	{};
	uint8_t puppa() const {return nfaces;};
	const uint32_t max_length;
	const uint8_t  nchildren;
	const uint8_t  nfaces;
	const uint8_t  nnodes;
	const uint8_t  octantBytes;
	const uint8_t  markerBytes;
	const uint8_t  levelBytes;
	const uint8_t  boolBytes;
	const uint8_t  oppface[4];
	const uint8_t  nodeface[4][2];
	const uint8_t  facenode[4][2];
	const int8_t   normals[4][2];
};

template<>
class global<3>
{
public:
	global() : max_length (uint32_t(pow(2.0,MAX_LEVEL_3D))),
	nchildren(8),
	nfaces(6),
	nnodes(8),
	nedges(12),
	octantBytes(uint8_t(sizeof(uint32_t)*3 + sizeof(uint8_t) + sizeof(int8_t) + (16)*sizeof(bool))),
	markerBytes(sizeof(int8_t)),
	levelBytes(sizeof(uint8_t)),
	boolBytes(sizeof(bool)),
	oppface{1,0,3,2,5,4},
	nodeface{{0,2,4},{1,2,4},{0,3,4},{1,3,4},{0,2,5},{1,2,5},{0,3,5},{1,3,5}},
	facenode{{0,2,4,6},{1,3,5,7},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7}},
	normals{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}
	{};

	const uint32_t max_length;
	const uint8_t  nchildren;
	const uint8_t  nfaces;
	const uint8_t  nedges;
	const uint8_t  nnodes;
	const uint8_t  octantBytes;
	const uint8_t  markerBytes;
	const uint8_t  levelBytes;
	const uint8_t  boolBytes;
	const uint8_t  oppface[6];
	const uint8_t  nodeface[8][3];
	const uint8_t  facenode[6][4];
	const int8_t   normals[6][3];
};




extern const global<2> global2D;
extern const global<3> global3D;


#endif /* GLOBAL_HPP_ */
