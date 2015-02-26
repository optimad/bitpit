#ifndef CLASS_GLOBAL_2D_TPP
#define CLASS_GLOBAL_2D_TPP

/*!
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Global variables used in PABLO - 2D specialization
 *
 *	Global variables are used in PABLO everywhere and they are public, i.e. each
 *	global variable can be used as constant by external codes.
 *
 */


// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //
template<>
class Class_Global<2>
{
public:
	Class_Global() :
		max_length (uint32_t(pow(2.0,MAX_LEVEL_2D))),
		nchildren(4),
		nfaces(4),
		nnodes(4),
		nnodesperface(2),
		octantBytes(uint8_t(sizeof(uint32_t)*2 + sizeof(uint8_t) + sizeof(int8_t) + (12)*sizeof(bool))),
		globalIndexBytes(uint8_t(sizeof(uint64_t))),
		markerBytes(sizeof(int8_t)),
		levelBytes(sizeof(uint8_t)),
		boolBytes(sizeof(bool))
{
		oppface[0] = 1;
		oppface[1] = 0;
		oppface[2] = 3;
		oppface[3] = 2;

		nodeface[0][0] = 0;
		nodeface[0][1] = 2;
		nodeface[1][0] = 1;
		nodeface[1][1] = 2;
		nodeface[2][0] = 0;
		nodeface[2][1] = 3;
		nodeface[3][0] = 1;
		nodeface[3][1] = 3;

		facenode[0][0] = 0;
		facenode[0][1] = 2;
		facenode[1][0] = 1;
		facenode[1][1] = 3;
		facenode[2][0] = 0;
		facenode[2][1] = 1;
		facenode[3][0] = 2;
		facenode[3][1] = 3;

		normals[0][0] = -1;
		normals[0][1] =  0;
		normals[0][2] =  0;
		normals[1][0] =  1;
		normals[1][1] =  0;
		normals[1][2] =  0;
		normals[2][0] =  0;
		normals[2][1] = -1;
		normals[2][2] =  0;
		normals[3][0] =  0;
		normals[3][1] =  1;
		normals[3][2] =  0;
};


	uint32_t max_length;			/**< Length of the logical domain */
	uint8_t  nchildren;				/**< Number of children of an octant */
	uint8_t  nfaces;				/**< Number of faces of an octant */
	uint8_t  nnodes;				/**< Number of nodes of an octant */
	uint8_t  nnodesperface;			/**< Number of nodes per face of an octant */
	uint8_t  octantBytes;			/**< Bytes occupation of an octant */
	uint8_t  globalIndexBytes;		/**< Bytes occupation of the index of an octant */
	uint8_t  markerBytes;			/**< Bytes occupation of the refinement marker of an octant */
	uint8_t  levelBytes;			/**< Bytes occupation of the level of an octant */
	uint8_t  boolBytes;				/**< Bytes occupation of a boolean */
	uint8_t  oppface[4];			/**< oppface[i] = Index of the face of an octant neighbour through the i-th face of the current octant */
	uint8_t  nodeface[4][2];		/**< nodeface[i][0:1] = local indices of faces sharing the i-th node of an octant */
	uint8_t  facenode[4][2];		/**< facenode[i][0:1] = local indices of nodes of the i-th face of an octant */
	int8_t   normals[4][3];			/**< Components (x,y,z) of the normals per face (z=0 in 2D) */
};

#endif
