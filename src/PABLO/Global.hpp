#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include <stdint.h>

/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Global variables used in PABLO
 *
 *	Global variables are used in PABLO everywhere and they are public, i.e. each
 *	global variable can be used asant by external codes.
 *
 *	Class Global is a class with static members initialized during the construction
 *	of a paratree object.
 *
 */
class Global{

	// =================================================================================== //
	// FRIENDSHIPS
	// =================================================================================== //
	friend class ParaTree;
	friend class LocalTree;

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	uint32_t m_maxLength;			/**< Length of the logical domain */
	uint8_t  m_nchildren;			/**< Number of children of an octant */
	uint8_t  m_nfaces;				/**< Number of faces of an octant */
	uint8_t  m_nedges;				/**< Number of edges of an octant */
	uint8_t  m_nnodes;				/**< Number of nodes of an octant */
	uint8_t  m_nnodesPerFace;		/**< Number of nodes per face of an octant */
	uint8_t  m_octantBytes;			/**< Bytes occupation of an octant */
	uint8_t  m_globalIndexBytes;	/**< Bytes occupation of the index of an octant */
	uint8_t  m_markerBytes;			/**< Bytes occupation of the refinement marker of an octant */
	uint8_t  m_levelBytes;			/**< Bytes occupation of the level of an octant */
	uint8_t  m_boolBytes;			/**< Bytes occupation of a boolean */
	uint8_t  m_oppFace[6];			/**< oppface[i] = Index of the face of an octant neighbour through the i-th face of the current octant */
	uint8_t  m_nodeFace[8][3];		/**< nodeface[i][0:1] = local indices of faces sharing the i-th node of an octant */
	uint8_t  m_faceNode[6][4];		/**< facenode[i][0:1] = local indices of nodes of the i-th face of an octant */
	uint8_t  m_edgeFace[12][2];		/**< edgeface[i][0:1] = local indices of faces sharing the i-th edge of an octant */
	int8_t   m_normals[6][3];		/**< Components (x,y,z) of the normals per face (z=0 in 2D) */
	int8_t   m_edgeCoeffs[12][3];	/**< Components (x,y,z) of the "normals" per edge */
	int8_t   m_nodeCoeffs[8][3];	/**< Components (x,y,z) of the "normals" per node */
	int8_t   m_maxLevel;			/**< Maximum allowed refinement level of octree */

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //
	uint8_t 	getBoolBytes();
	void 		getEdgecoeffs(int8_t edgecoeffs[12][3]);
	void 		getEdgeface(uint8_t edgeface[12][2]);
	void 		getFacenode(uint8_t facenode[6][4]);
	uint8_t 	getGlobalIndexBytes();
	uint8_t 	getLevelBytes();
	uint8_t 	getMarkerBytes();
	uint32_t 	getMaxLength();
	int8_t 		getMaxLevel();
	uint8_t 	getNchildren();
	uint8_t 	getNedges();
	uint8_t 	getNfaces();
	uint8_t 	getNnodes();
	uint8_t 	getNnodesperface();
	void 		getNodecoeffs(int8_t nodecoeffs[8][3]);
	void 		getNodeface(uint8_t nodeface[8][3]);
	void 		getNormals(int8_t normals[6][3]);
	uint8_t 	getOctantBytes();
	void 		getOppface(uint8_t oppface[6]);

	void 		setGlobal(int8_t maxlevel, uint8_t dim);

};

/*  @} */

#endif /* GLOBAL_HPP_ */
