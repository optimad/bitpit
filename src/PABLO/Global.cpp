/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "Global.hpp"

namespace bitpit {

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

const int8_t Global::m_maxLevel = 20;
const uint32_t Global::m_maxLength = uint32_t(1 << Global::m_maxLevel);


/*!Get the bytes occupation of a boolean.
 * \return Bytes occupation of a boolean.
 */
uint8_t
Global::getBoolBytes() const {
	return m_boolBytes;
}

/*! Get the components of the "normals" per edge (bisector of edge)
 * \return Components (x,y,z) of the "normals" per edge.
 */
void
Global::getEdgecoeffs(int8_t edgecoeffs_[12][3]) const {
	for (int i=0; i<12; i++){
		for (int j=0; j<3; j++){
		edgecoeffs_[i][j] = m_edgeCoeffs[i][j];
		}
	}
}

/*! Get the connectivity edge-face.
 * \param[out] edgeface_ edgeface[i][0:1] = local indices of faces sharing
 * the i-th edge of an octant.
 */
void
Global::getEdgeface(uint8_t edgeface_[12][2]) const {
	for (int i=0; i<12; i++){
		for (int j=0; j<2; j++){
			edgeface_[i][j] = m_edgeFace[i][j];
		}
	}
}

/*! Get the connectivity edge-node.
 * \param[out] edgeNode edgeNode[i][0:1] = local indices of nodes sharing
 * the i-th edge of an octant.
 */
void
Global::getEdgenode(uint8_t edgeNode[12][2]) const {
    for(int i = 0; i < 12; ++i){
        for(int j = 0; j < 2; ++j){
            edgeNode[i][j] = m_edgeNode[i][j];
        }
    }
};

/*! Get the connectivity edge-node.
 * \param[in] edge Edge index
 * \param[out] edgeNode edgeNode[0:1] = local indices of nodes sharing
 * the edge edge of an octant.
 */
void
Global::getEdgenode(uint8_t edge, uint8_t edgeNode[2]) const {
    for(int j = 0; j < 2; ++j){
        edgeNode[j] = m_edgeNode[edge][j];
    }
};

/*! Get the connectivity face-node.
 * \param[out] facenode_ facenode[i][0:1] = local indices of nodes
 * of the i-th face of an octant.
 */
void
Global::getFacenode(uint8_t facenode_[6][4]) const {
	for (int i=0; i<6; i++){
		for (int j=0; j<4; j++){
			facenode_[i][j] = m_faceNode[i][j];
		}
	}
};

/*! Get the connectivity node-edge.
 * \param[out] nodeEdge_ nodeedge[i][0:2] = local indices of edges
 * sharing the i-th node of an octant.
 */
void
Global::getNodeedge(uint8_t nodeedge_[8][3]) const {
    for (int i=0; i<8; i++){
        for (int j=0; j<3; j++){
            nodeedge_[i][j] = m_nodeEdge[i][j];
        }
    }
}

/*!Get the bytes occupation of the global index of an octant.
 * \return Bytes occupation of a global index.
 */
uint8_t
Global::getGlobalIndexBytes() const {
	return m_globalIndexBytes;
}

/*!Get the bytes occupation of the level of an octant.
 * \return Bytes occupation of level.
 */
uint8_t
Global::getLevelBytes() const {
	return m_levelBytes;
}

/*!Get the bytes occupation of the marker of an octant.
 * \return Bytes occupation of marker.
 */
uint8_t
Global::getMarkerBytes() const {
	return m_markerBytes;
}

/*!Get the length of the logical domain.
 * \return Length of the logical domain.
 */
uint32_t
Global::getMaxLength()  {
	return m_maxLength;
}

/*!Get the maximum allowed refinement level of octree.
 * \return Maximum allowed refinement level of octree.
 */
int8_t
Global::getMaxLevel()  {
	return m_maxLevel;
}

/*!Get the number of children of an octant.
 * \return Number of children of an octant.
 */
uint8_t
Global::getNchildren() const {
	return m_nchildren;
}

/*!Get the number of edges of an octant.
 * \return Number of edges of an octant.
 */
uint8_t
Global::getNedges() const {
	return m_nedges;
}

/*!Get the number of faces of an octant.
 * \return Number of faces of an octant.
 */
uint8_t
Global::getNfaces() const {
	return m_nfaces;
}

/*!Get the number of nodes of an octant.
 * \return Number of nodes of an octant.
 */
uint8_t
Global::getNnodes() const {
	return m_nnodes;
}

/*!Get the number of nodes of a face an octant.
 * \return Number of nodes of a face of an octant.
 */
uint8_t
Global::getNnodesperface() const {
	return m_nnodesPerFace;
}

/*! Get the components of the "normals" per node (bisector of node)
 * \return Components (x,y,z) of the "normals" per node.
 */
void
Global::getNodecoeffs(int8_t nodecoeffs_[8][3]) const {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodecoeffs_[i][j] = m_nodeCoeffs[i][j];
		}
	}
}

/*! Get the connectivity node-face.
 * \param[out] nodeface_ nodeface[i][0:1] = local indices of faces
 * sharing the i-th node of an octant.
 */
void
Global::getNodeface(uint8_t nodeface_[8][3]) const {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodeface_[i][j] = m_nodeFace[i][j];
		}
	}
}

/*! Get the components of the normals of faces.
 * \return Components (x,y,z) of the normals per face (z=0 in 2D).
 */
void
Global::getNormals(int8_t normals_[6][3]) const {
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			normals_[i][j] = m_normals[i][j];
		}
	}
}

/*!Get the bytes occupation of an octant.
 * \return Bytes occupation of an octant.
 */
uint8_t
Global::getOctantBytes() const {
	return m_octantBytes;
}
/*! Get the index of the opposite face of each face.
* \return oppface[i] = index of the face of an octant neighbour
* through the i-th face of the current octant.
*/
void
Global::getOppface(uint8_t oppface_[6]) const {
	for (int j=0; j<6; j++){
		oppface_[j] = m_oppFace[j];
	}
}

/*! Get the indeces of the parallel edges for each edge.
* \param[out] parallelEdges[i][j] = matrix containing parallelEdges j for each edge i
*/
void
Global::getParallelEdges(uint8_t parallelEdges[12][3]) const {
    for(int i = 0; i < 12; ++i){
        for(int j = 0; j < 3; ++j){
            parallelEdges[i][j] = m_parallelEdges[i][j];
        }
    }
}

/*! Get the indeces of the parallel edges for an edge.
* \param[in] edge Edge index
* \param[out] parallelEdges matrix containing parallelEdges for each edge
*/
void
Global::getParallelEdges(uint8_t edge, uint8_t parallelEdges[]) const {
    for(int j = 0; j < 3; ++j){
        parallelEdges[j] = m_parallelEdges[edge][j];
    }
}

/*! Initialization of static global variables.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
void
Global::initialize(uint8_t dim) {

	if (dim>3) dim = 3;
	if (dim<2) dim = 2;

	m_nchildren 			= 1<<dim;
	m_nfaces 				= 2*dim;
	m_nnodes				= 1<<dim;
	m_nedges 				= (dim-2)*12;
	m_nnodesPerFace 		= 1<<(dim-1);
	m_octantBytes 			= uint8_t(sizeof(uint32_t)*3 + sizeof(uint8_t) + sizeof(int8_t) + (17)*sizeof(bool));
	m_globalIndexBytes 		= uint8_t(sizeof(uint64_t));
	m_markerBytes 			= sizeof(int8_t);
	m_levelBytes 			= sizeof(uint8_t);
	m_boolBytes 			= sizeof(bool);

	m_oppFace[0] = 1;
	m_oppFace[1] = 0;
	m_oppFace[2] = 3;
	m_oppFace[3] = 2;
	m_oppFace[4] = 5;
	m_oppFace[5] = 4;

	m_nodeFace[0][0] = 0;
	m_nodeFace[0][1] = 2;
	m_nodeFace[0][2] = 4;
	m_nodeFace[1][0] = 1;
	m_nodeFace[1][1] = 2;
	m_nodeFace[1][2] = 4;
	m_nodeFace[2][0] = 0;
	m_nodeFace[2][1] = 3;
	m_nodeFace[2][2] = 4;
	m_nodeFace[3][0] = 1;
	m_nodeFace[3][1] = 3;
	m_nodeFace[3][2] = 4;
	m_nodeFace[4][0] = 0;
	m_nodeFace[4][1] = 2;
	m_nodeFace[4][2] = 5;
	m_nodeFace[5][0] = 1;
	m_nodeFace[5][1] = 2;
	m_nodeFace[5][2] = 5;
	m_nodeFace[6][0] = 0;
	m_nodeFace[6][1] = 3;
	m_nodeFace[6][2] = 5;
	m_nodeFace[7][0] = 1;
	m_nodeFace[7][1] = 3;
	m_nodeFace[7][2] = 5;

	m_nodeEdge[0][0] = 0;
	m_nodeEdge[0][1] = 2;
	m_nodeEdge[0][2] = 4;
	m_nodeEdge[1][0] = 1;
	m_nodeEdge[1][1] = 2;
	m_nodeEdge[1][2] = 5;
	m_nodeEdge[2][0] = 0;
	m_nodeEdge[2][1] = 3;
	m_nodeEdge[2][2] = 6;
	m_nodeEdge[3][0] = 1;
	m_nodeEdge[3][1] = 3;
	m_nodeEdge[3][2] = 7;
	m_nodeEdge[4][0] = 4;
	m_nodeEdge[4][1] = 8;
	m_nodeEdge[4][2] = 10;
	m_nodeEdge[5][0] = 5;
	m_nodeEdge[5][1] = 9;
	m_nodeEdge[5][2] = 10;
	m_nodeEdge[6][0] = 6;
	m_nodeEdge[6][1] = 8;
	m_nodeEdge[6][2] = 11;
	m_nodeEdge[7][0] = 7;
	m_nodeEdge[7][1] = 9;
	m_nodeEdge[7][2] = 11;

	m_edgeNode[0][0] = 0;
	m_edgeNode[0][1] = 2;
	m_edgeNode[1][0] = 1;
	m_edgeNode[1][1] = 3;
	m_edgeNode[2][0] = 0;
	m_edgeNode[2][1] = 1;
	m_edgeNode[3][0] = 2;
	m_edgeNode[3][1] = 3;
	m_edgeNode[4][0] = 0;
	m_edgeNode[4][1] = 4;
	m_edgeNode[5][0] = 1;
	m_edgeNode[5][1] = 5;
	m_edgeNode[6][0] = 2;
	m_edgeNode[6][1] = 6;
	m_edgeNode[7][0] = 3;
	m_edgeNode[7][1] = 7;
	m_edgeNode[8][0] = 4;
	m_edgeNode[8][1] = 6;
	m_edgeNode[9][0] = 5;
	m_edgeNode[9][1] = 7;
	m_edgeNode[10][0] = 4;
	m_edgeNode[10][1] = 5;
	m_edgeNode[11][0] = 6;
	m_edgeNode[11][1] = 7;

	m_faceNode[0][0] = 0;
	m_faceNode[0][1] = 2;
	m_faceNode[0][2] = 4;
	m_faceNode[0][3] = 6;
	m_faceNode[1][0] = 1;
	m_faceNode[1][1] = 3;
	m_faceNode[1][2] = 5;
	m_faceNode[1][3] = 7;
	m_faceNode[2][0] = 0;
	m_faceNode[2][1] = 1;
	m_faceNode[2][2] = 4;
	m_faceNode[2][3] = 5;
	m_faceNode[3][0] = 2;
	m_faceNode[3][1] = 3;
	m_faceNode[3][2] = 6;
	m_faceNode[3][3] = 7;
	m_faceNode[4][0] = 0;
	m_faceNode[4][1] = 1;
	m_faceNode[4][2] = 2;
	m_faceNode[4][3] = 3;
	m_faceNode[5][0] = 4;
	m_faceNode[5][1] = 5;
	m_faceNode[5][2] = 6;
	m_faceNode[5][3] = 7;

	m_edgeFace[0][0] = 0;
	m_edgeFace[0][1] = 4;
	m_edgeFace[1][0] = 1;
	m_edgeFace[1][1] = 4;
	m_edgeFace[2][0] = 2;
	m_edgeFace[2][1] = 4;
	m_edgeFace[3][0] = 3;
	m_edgeFace[3][1] = 4;
	m_edgeFace[4][0] = 0;
	m_edgeFace[4][1] = 2;
	m_edgeFace[5][0] = 1;
	m_edgeFace[5][1] = 2;
	m_edgeFace[6][0] = 0;
	m_edgeFace[6][1] = 3;
	m_edgeFace[7][0] = 1;
	m_edgeFace[7][1] = 3;
	m_edgeFace[8][0] = 0;
	m_edgeFace[8][1] = 5;
	m_edgeFace[9][0] = 1;
	m_edgeFace[9][1] = 5;
	m_edgeFace[10][0] = 2;
	m_edgeFace[10][1] = 5;
	m_edgeFace[11][0] = 3;
	m_edgeFace[11][1] = 5;

	m_normals[0][0] = -1;
	m_normals[0][1] =  0;
	m_normals[0][2] =  0;
	m_normals[1][0] =  1;
	m_normals[1][1] =  0;
	m_normals[1][2] =  0;
	m_normals[2][0] =  0;
	m_normals[2][1] = -1;
	m_normals[2][2] =  0;
	m_normals[3][0] =  0;
	m_normals[3][1] =  1;
	m_normals[3][2] =  0;
	m_normals[4][0] =  0;
	m_normals[4][1] =  0;
	m_normals[4][2] = -1;
	m_normals[5][0] =  0;
	m_normals[5][1] =  0;
	m_normals[5][2] =  1;

	m_edgeCoeffs[0][0] = -1;
	m_edgeCoeffs[0][1] =  0;
	m_edgeCoeffs[0][2] = -1;
	m_edgeCoeffs[1][0] =  1;
	m_edgeCoeffs[1][1] =  0;
	m_edgeCoeffs[1][2] = -1;
	m_edgeCoeffs[2][0] =  0;
	m_edgeCoeffs[2][1] = -1;
	m_edgeCoeffs[2][2] = -1;
	m_edgeCoeffs[3][0] =  0;
	m_edgeCoeffs[3][1] =  1;
	m_edgeCoeffs[3][2] = -1;
	m_edgeCoeffs[4][0] = -1;
	m_edgeCoeffs[4][1] = -1;
	m_edgeCoeffs[4][2] =  0;
	m_edgeCoeffs[5][0] =  1;
	m_edgeCoeffs[5][1] = -1;
	m_edgeCoeffs[5][2] =  0;
	m_edgeCoeffs[6][0] = -1;
	m_edgeCoeffs[6][1] =  1;
	m_edgeCoeffs[6][2] =  0;
	m_edgeCoeffs[7][0] =  1;
	m_edgeCoeffs[7][1] =  1;
	m_edgeCoeffs[7][2] =  0;
	m_edgeCoeffs[8][0] = -1;
	m_edgeCoeffs[8][1] =  0;
	m_edgeCoeffs[8][2] =  1;
	m_edgeCoeffs[9][0] =  1;
	m_edgeCoeffs[9][1] =  0;
	m_edgeCoeffs[9][2] =  1;
	m_edgeCoeffs[10][0] =  0;
	m_edgeCoeffs[10][1] = -1;
	m_edgeCoeffs[10][2] =  1;
	m_edgeCoeffs[11][0] =  0;
	m_edgeCoeffs[11][1] =  1;
	m_edgeCoeffs[11][2] =  1;

	m_nodeCoeffs[0][0] = -1;
	m_nodeCoeffs[0][1] = -1;
	m_nodeCoeffs[0][2] = -1;
	m_nodeCoeffs[1][0] =  1;
	m_nodeCoeffs[1][1] = -1;
	m_nodeCoeffs[1][2] = -1;
	m_nodeCoeffs[2][0] = -1;
	m_nodeCoeffs[2][1] =  1;
	m_nodeCoeffs[2][2] = -1;
	m_nodeCoeffs[3][0] =  1;
	m_nodeCoeffs[3][1] =  1;
	m_nodeCoeffs[3][2] = -1;
	m_nodeCoeffs[4][0] = -1;
	m_nodeCoeffs[4][1] = -1;
	m_nodeCoeffs[4][2] =  1;
	m_nodeCoeffs[5][0] =  1;
	m_nodeCoeffs[5][1] = -1;
	m_nodeCoeffs[5][2] =  1;
	m_nodeCoeffs[6][0] = -1;
	m_nodeCoeffs[6][1] =  1;
	m_nodeCoeffs[6][2] =  1;
	m_nodeCoeffs[7][0] =  1;
	m_nodeCoeffs[7][1] =  1;
	m_nodeCoeffs[7][2] =  1;

	m_parallelEdges[0][0] = 1;
	m_parallelEdges[0][1] = 8;
	m_parallelEdges[0][2] = 9;
	m_parallelEdges[1][0] = 0;
	m_parallelEdges[1][1] = 8;
	m_parallelEdges[1][2] = 9;
	m_parallelEdges[2][0] = 3;
	m_parallelEdges[2][1] = 10;
	m_parallelEdges[2][2] = 11;
	m_parallelEdges[3][0] = 2;
	m_parallelEdges[3][1] = 10;
	m_parallelEdges[3][2] = 11;
	m_parallelEdges[4][0] = 5;
	m_parallelEdges[4][1] = 6;
	m_parallelEdges[4][2] = 7;
	m_parallelEdges[5][0] = 4;
	m_parallelEdges[5][1] = 6;
	m_parallelEdges[5][2] = 7;
	m_parallelEdges[6][0] = 4;
	m_parallelEdges[6][1] = 5;
	m_parallelEdges[6][2] = 7;
	m_parallelEdges[7][0] = 4;
	m_parallelEdges[7][1] = 5;
	m_parallelEdges[7][2] = 6;
	m_parallelEdges[8][0] = 0;
	m_parallelEdges[8][1] = 1;
	m_parallelEdges[8][2] = 9;
	m_parallelEdges[9][0] = 0;
	m_parallelEdges[9][1] = 1;
	m_parallelEdges[9][2] = 8;
	m_parallelEdges[10][0] = 2;
	m_parallelEdges[10][1] = 3;
	m_parallelEdges[10][2] = 11;
	m_parallelEdges[11][0] = 2;
	m_parallelEdges[11][1] = 3;
	m_parallelEdges[11][2] = 10;

}

// =================================================================================== //

}
