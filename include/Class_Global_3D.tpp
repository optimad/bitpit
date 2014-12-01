#ifndef CLASS_GLOBAL_3D_TPP
#define CLASS_GLOBAL_3D_TPP

/*!
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Global variables used in PABLO - 3D specialization
 *
 *	Global variables are used in PABLO everywhere and they are public, i.e. each
 *	global variable can be used as constant by external codes.
 *
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //
template<>
class Class_Global<3>
{
public:
	Class_Global() : max_length (uint32_t(pow(2.0,MAX_LEVEL_3D))),
	nchildren(8),
	nfaces(6),
	nnodes(8),
	nedges(12),
	nnodesperface(4),
	octantBytes(uint8_t(sizeof(uint32_t)*3 + sizeof(uint8_t) + sizeof(int8_t) + (16)*sizeof(bool))),
	globalIndexBytes(uint8_t(sizeof(uint64_t))),
	markerBytes(sizeof(int8_t)),
	levelBytes(sizeof(uint8_t)),
	boolBytes(sizeof(bool)),
	oppface{1,0,3,2,5,4},
	nodeface{{0,2,4},{1,2,4},{0,3,4},{1,3,4},{0,2,5},{1,2,5},{0,3,5},{1,3,5}},
	facenode{{0,2,4,6},{1,3,5,7},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7}},
	edgeface{{0,4},{1,4},{2,4},{3,4},{0,2},{1,2},{0,3},{1,3},{0,5},{1,5},{2,5},{3,5}},
	normals{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}
	{};

	const uint32_t max_length;		/**< Length of the logical domain */
	const uint8_t  nchildren;		/**< Number of children of an octant */
	const uint8_t  nfaces;			/**< Number of faces of an octant */
	const uint8_t  nedges;			/**< Number of edges of an octant */
	const uint8_t  nnodes;			/**< Number of nodes of an octant */
	const uint8_t  nnodesperface;	/**< Number of nodes per face of an octant */
	const uint8_t  octantBytes;		/**< Bytes occupation of an octant */
	const uint8_t  globalIndexBytes;/**< Bytes occupation of the index of an octant */
	const uint8_t  markerBytes;		/**< Bytes occupation of the refinement marker of an octant */
	const uint8_t  levelBytes;		/**< Bytes occupation of the level of an octant */
	const uint8_t  boolBytes;		/**< Bytes occupation of a boolean */
	const uint8_t  oppface[6];		/**< oppface[i] = Index of the face of an octant neighbour through the i-th face of the current octant */
	const uint8_t  nodeface[8][3];	/**< nodeface[i][0:1] = local indices of faces sharing the i-th node of an octant */
	const uint8_t  facenode[6][4];	/**< facenode[i][0:1] = local indices of nodes of the i-th face of an octant */
	const uint8_t  edgeface[12][2];	/**< edgeface[i][0:1] = local indices of faces sharing the i-th edge of an octant */
	const int8_t   normals[6][3];	/**< Components (x,y,z) of the normals per face */

};

#endif



