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
	boolBytes(sizeof(bool))
{
		oppface[0] = 1;
		oppface[1] = 0;
		oppface[2] = 3;
		oppface[3] = 2;
		oppface[4] = 5;
		oppface[5] = 4;

		nodeface[0][0] = 0;
		nodeface[0][1] = 2;
		nodeface[0][2] = 4;
		nodeface[1][0] = 1;
		nodeface[1][1] = 2;
		nodeface[1][2] = 4;
		nodeface[2][0] = 0;
		nodeface[2][1] = 3;
		nodeface[2][2] = 4;
		nodeface[3][0] = 1;
		nodeface[3][1] = 3;
		nodeface[3][2] = 4;
		nodeface[4][0] = 0;
		nodeface[4][1] = 2;
		nodeface[4][2] = 5;
		nodeface[5][0] = 1;
		nodeface[5][1] = 2;
		nodeface[5][2] = 5;
		nodeface[6][0] = 0;
		nodeface[6][1] = 3;
		nodeface[6][2] = 5;
		nodeface[7][0] = 1;
		nodeface[7][1] = 3;
		nodeface[7][2] = 5;

		facenode[0][0] = 0;
		facenode[0][1] = 2;
		facenode[0][2] = 4;
		facenode[0][3] = 6;
		facenode[1][0] = 1;
		facenode[1][1] = 3;
		facenode[1][2] = 5;
		facenode[1][3] = 7;
		facenode[2][0] = 0;
		facenode[2][1] = 1;
		facenode[2][2] = 4;
		facenode[2][3] = 5;
		facenode[3][0] = 2;
		facenode[3][1] = 3;
		facenode[3][2] = 6;
		facenode[3][3] = 7;
		facenode[4][0] = 0;
		facenode[4][1] = 1;
		facenode[4][2] = 2;
		facenode[4][3] = 3;
		facenode[5][0] = 4;
		facenode[5][1] = 5;
		facenode[5][2] = 6;
		facenode[5][3] = 7;

		edgeface[0][0] = 0;
		edgeface[0][1] = 4;
		edgeface[1][0] = 1;
		edgeface[1][1] = 4;
		edgeface[2][0] = 2;
		edgeface[2][1] = 4;
		edgeface[3][0] = 3;
		edgeface[3][1] = 4;
		edgeface[4][0] = 0;
		edgeface[4][1] = 2;
		edgeface[5][0] = 1;
		edgeface[5][1] = 2;
		edgeface[6][0] = 0;
		edgeface[6][1] = 3;
		edgeface[7][0] = 1;
		edgeface[7][1] = 3;
		edgeface[8][0] = 0;
		edgeface[8][1] = 5;
		edgeface[9][0] = 1;
		edgeface[9][1] = 5;
		edgeface[10][0] = 2;
		edgeface[10][1] = 5;
		edgeface[11][0] = 3;
		edgeface[11][1] = 5;

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
		normals[4][0] =  0;
		normals[4][1] =  0;
		normals[4][2] = -1;
		normals[5][0] =  0;
		normals[5][1] =  0;
		normals[5][2] =  1;

		edgecoeffs[0][0] = -1;
		edgecoeffs[0][1] =  0;
		edgecoeffs[0][2] = -1;
		edgecoeffs[1][0] =  1;
		edgecoeffs[1][1] =  0;
		edgecoeffs[1][2] = -1;
		edgecoeffs[2][0] =  0;
		edgecoeffs[2][1] = -1;
		edgecoeffs[2][2] = -1;
		edgecoeffs[3][0] =  0;
		edgecoeffs[3][1] =  1;
		edgecoeffs[3][2] = -1;
		edgecoeffs[4][0] = -1;
		edgecoeffs[4][1] = -1;
		edgecoeffs[4][2] =  0;
		edgecoeffs[5][0] =  1;
		edgecoeffs[5][1] = -1;
		edgecoeffs[5][2] =  0;
		edgecoeffs[6][0] = -1;
		edgecoeffs[6][1] =  1;
		edgecoeffs[6][2] =  0;
		edgecoeffs[7][0] =  1;
		edgecoeffs[7][1] =  1;
		edgecoeffs[7][2] =  0;
		edgecoeffs[8][0] = -1;
		edgecoeffs[8][1] =  0;
		edgecoeffs[8][2] =  1;
		edgecoeffs[9][0] =  1;
		edgecoeffs[9][1] =  0;
		edgecoeffs[9][2] =  1;
		edgecoeffs[10][0] =  0;
		edgecoeffs[10][1] = -1;
		edgecoeffs[10][2] =  1;
		edgecoeffs[11][0] =  0;
		edgecoeffs[11][1] =  1;
		edgecoeffs[11][2] =  1;
};


	uint32_t max_length;		/**< Length of the logical domain */
	uint8_t  nchildren;			/**< Number of children of an octant */
	uint8_t  nfaces;			/**< Number of faces of an octant */
	uint8_t  nedges;			/**< Number of edges of an octant */
	uint8_t  nnodes;			/**< Number of nodes of an octant */
	uint8_t  nnodesperface;		/**< Number of nodes per face of an octant */
	uint8_t  octantBytes;		/**< Bytes occupation of an octant */
	uint8_t  globalIndexBytes;	/**< Bytes occupation of the index of an octant */
	uint8_t  markerBytes;		/**< Bytes occupation of the refinement marker of an octant */
	uint8_t  levelBytes;		/**< Bytes occupation of the level of an octant */
	uint8_t  boolBytes;			/**< Bytes occupation of a boolean */
	uint8_t  oppface[6];		/**< oppface[i] = Index of the face of an octant neighbour through the i-th face of the current octant */
	uint8_t  nodeface[8][3];	/**< nodeface[i][0:1] = local indices of faces sharing the i-th node of an octant */
	uint8_t  facenode[6][4];	/**< facenode[i][0:1] = local indices of nodes of the i-th face of an octant */
	uint8_t  edgeface[12][2];	/**< edgeface[i][0:1] = local indices of faces sharing the i-th edge of an octant */
	int8_t   normals[6][3];		/**< Components (x,y,z) of the normals per face */
	int8_t   edgecoeffs[12][3];	/**< Components (x,y,z) of the "normals" per edge */

};

#endif



