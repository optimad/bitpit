/*
 * Class_Global_2D.tpp
 *
 *  Created on: 18/apr/2014
 *      Author: Marco Cisternino
 */
#ifndef CLASS_GLOBAL_2D_TPP
#define CLASS_GLOBAL_2D_TPP

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //
template<>
class Class_Global<2>
{
public:
	Class_Global() : max_length (uint32_t(pow(2.0,MAX_LEVEL_2D))),
	nchildren(4),
	nfaces(4),
	nnodes(4),
	nnodesperface(2),
	octantBytes(uint8_t(sizeof(uint32_t)*2 + sizeof(uint8_t) + sizeof(int8_t) + (12)*sizeof(bool))),
	globalIndexBytes(uint8_t(sizeof(uint64_t))),
	markerBytes(sizeof(int8_t)),
	levelBytes(sizeof(uint8_t)),
	boolBytes(sizeof(bool)),
	oppface{1,0,3,2},
	nodeface{{0,2},{1,2},{0,3},{1,3}},
	facenode{{0,2},{1,3},{0,1},{2,3}},
	normals{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0}}
	{};
	uint8_t puppa() const {return nfaces;};
	const uint32_t max_length;
	const uint8_t  nchildren;
	const uint8_t  nfaces;
	const uint8_t  nnodes;
	const uint8_t  nnodesperface;
	const uint8_t  octantBytes;
	const uint8_t  globalIndexBytes;
	const uint8_t  markerBytes;
	const uint8_t  levelBytes;
	const uint8_t  boolBytes;
	const uint8_t  oppface[4];
	const uint8_t  nodeface[4][2];
	const uint8_t  facenode[4][2];
	const int8_t   normals[4][3];
};

#endif
