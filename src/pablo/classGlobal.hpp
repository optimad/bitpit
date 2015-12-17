#ifndef CLASSGLOBAL_HPP_
#define CLASSGLOBAL_HPP_

#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>

/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Global variables used in PABLO
 *
 *	Global variables are used in PABLO everywhere and they are public, i.e. each
 *	global variable can be used as constant by external codes.
 *
 *	Class classGlobal is a class with static members initialized during the construction of a paratree object.
 *
 */
class classGlobal{

public:
static	uint32_t max_length;			/**< Length of the logical domain */
static	uint8_t  nchildren;				/**< Number of children of an octant */
static	uint8_t  nfaces;				/**< Number of faces of an octant */
static  uint8_t  nedges;				/**< Number of edges of an octant */
static	uint8_t  nnodes;				/**< Number of nodes of an octant */
static	uint8_t  nnodesperface;			/**< Number of nodes per face of an octant */
static	uint8_t  octantBytes;			/**< Bytes occupation of an octant */
static	uint8_t  globalIndexBytes;		/**< Bytes occupation of the index of an octant */
static	uint8_t  markerBytes;			/**< Bytes occupation of the refinement marker of an octant */
static	uint8_t  levelBytes;			/**< Bytes occupation of the level of an octant */
static	uint8_t  boolBytes;				/**< Bytes occupation of a boolean */
static	uint8_t  oppface[4];			/**< oppface[i] = Index of the face of an octant neighbour through the i-th face of the current octant */
static	uint8_t  nodeface[8][3];		/**< nodeface[i][0:1] = local indices of faces sharing the i-th node of an octant */
static	uint8_t  facenode[6][3];		/**< facenode[i][0:1] = local indices of nodes of the i-th face of an octant */
static  uint8_t  edgeface[12][2];		/**< edgeface[i][0:1] = local indices of faces sharing the i-th edge of an octant */
static	int8_t   normals[6][3];			/**< Components (x,y,z) of the normals per face (z=0 in 2D) */
static  int8_t   edgecoeffs[12][3];		/**< Components (x,y,z) of the "normals" per edge */
static	int8_t   nodecoeffs[8][3];		/**< Components (x,y,z) of the "normals" per node */
static  int8_t  MAX_LEVEL;				/**< Maximum allowed refinement level of octree */

// =================================================================================== //

/*! Initialization of static global variables.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
static void setGlobal(uint8_t dim);

// =================================================================================== //

};


/*  @} */
 
#endif /* CLASSGLOBAL_HPP_ */
