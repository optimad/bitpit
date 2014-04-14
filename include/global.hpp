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
	global() : max_length (uint32_t(pow(2.0,MAX_LEVEL_2D))){};
	void puppa(){};

	const uint32_t max_length;
//	const uint8_t  nchildren;
//	const uint8_t  nfaces;
//	const uint8_t  nedges;
//	const uint8_t  nnodes;
//	const uint8_t  octantBytes;
//	const uint8_t  markerBytes;
//	const uint8_t  levelBytes;
//	const uint8_t  boolBytes;
//	const uint8_t  oppface[4];
//	const uint8_t  nodeface[nnodes][2];
//	const uint8_t  facenode[nfaces][2];
//	const int8_t   normals[nfaces][2];
};

//template<>
//class global<3>
//{
//	global();
//
//	const uint32_t max_length;
//	const uint8_t  nchildren;
//	const uint8_t  nfaces;
//	const uint8_t  nedges;
//	const uint8_t  nnodes;
//	const uint8_t  octantBytes;
//	const uint8_t  markerBytes;
//	const uint8_t  levelBytes;
//	const uint8_t  boolBytes;
//	const uint8_t  oppface[6];
//	const uint8_t  nodeface[nnodes][3];
//	const uint8_t  facenode[nfaces][4];
//	const int8_t   normals[nfaces][2];
//};

#include "global.tpp"

//extern const global<2> global2D;
//extern const global<3> global3D;


#endif /* GLOBAL_HPP_ */
