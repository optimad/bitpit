/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

uint8_t
Global::getBoolBytes()  {
	return m_boolBytes;
}

void
Global::getEdgecoeffs(int8_t edgecoeffs_[12][3])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<3; j++){
		edgecoeffs_[i][j] = m_edgeCoeffs[i][j];
		}
	}
}

void
Global::getEdgeface(uint8_t edgeface_[12][2])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<2; j++){
			edgeface_[i][j] = m_edgeFace[i][j];
		}
	}
}

void
Global::getFacenode(uint8_t facenode_[6][4])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<4; j++){
			facenode_[i][j] = m_faceNode[i][j];
		}
	}
}

uint8_t
Global::getGlobalIndexBytes()  {
	return m_globalIndexBytes;
}

uint8_t
Global::getLevelBytes()  {
	return m_levelBytes;
}

uint8_t
Global::getMarkerBytes()  {
	return m_markerBytes;
}

uint32_t
Global::getMaxLength()  {
	return m_maxLength;
}

int8_t
Global::getMaxLevel()  {
	return m_maxLevel;
}

uint8_t
Global::getNchildren()  {
	return m_nchildren;
}

uint8_t
Global::getNedges()  {
	return m_nedges;
}

uint8_t
Global::getNfaces()  {
	return m_nfaces;
}

uint8_t
Global::getNnodes()  {
	return m_nnodes;
}

uint8_t
Global::getNnodesperface()  {
	return m_nnodesPerFace;
}

void
Global::getNodecoeffs(int8_t nodecoeffs_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodecoeffs_[i][j] = m_nodeCoeffs[i][j];
		}
	}
}

void
Global::getNodeface(uint8_t nodeface_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodeface_[i][j] = m_nodeFace[i][j];
		}
	}
}

void
Global::getNormals(int8_t normals_[6][3])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			normals_[i][j] = m_normals[i][j];
		}
	}
}

uint8_t
Global::getOctantBytes()  {
	return m_octantBytes;
}

void
Global::getOppface(uint8_t oppface_[6])  {
	for (int j=0; j<6; j++){
		oppface_[j] = m_oppFace[j];
	}
}

/*! Initialization of static global variables.
 * \param[in] maxlevel Space dimension; only value equal to 2 or 3 are accepted.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
void
Global::setGlobal(int8_t maxlevel, uint8_t dim){

	if (dim>3) dim = 3;
	if (dim<2) dim = 2;

	m_maxLevel = maxlevel;

	m_maxLength 			= uint32_t(1<<m_maxLevel);
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

}

// =================================================================================== //

}
