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

#include "tree_constants.hpp"
#include "Octant.hpp"

namespace bitpit {

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

/*! Get an instance of global variables.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
const TreeConstants &
TreeConstants::instance(uint8_t dim) {
	assert(dim == 2 || dim == 3);
	return instances()[dim];
}

/*! Get the instances of global variables.
 */
const TreeConstants::Instances &
TreeConstants::instances() {
	static TreeConstants::Instances instances = {{
		TreeConstants(0),
		TreeConstants(0),
		TreeConstants(2),
		TreeConstants(3)
	}};

	return instances;
}

/*! Default constructor.
 */
TreeConstants::TreeConstants() : TreeConstants(0) {
}

/*! Constructor.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
TreeConstants::TreeConstants(uint8_t dim) {
	initialize(dim);
}

/*! Initialization of static global variables.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
void
TreeConstants::initialize(uint8_t dim) {

	if (dim>3) dim = 3;
	if (dim<2) dim = 2;

	nChildren 			= uint8_t(1)<<dim;
	nFaces 				= 2*dim;
	nNodes				= uint8_t(1)<<dim;
	nEdges 				= (dim-2)*12;
	nNodesPerFace 		= uint8_t(1)<<(dim-1);

	oppositeFace[0] = 1;
	oppositeFace[1] = 0;
	oppositeFace[2] = 3;
	oppositeFace[3] = 2;
	oppositeFace[4] = 5;
	oppositeFace[5] = 4;

	nodeFace[0][0] = 0;
	nodeFace[0][1] = 2;
	nodeFace[0][2] = 4;
	nodeFace[1][0] = 1;
	nodeFace[1][1] = 2;
	nodeFace[1][2] = 4;
	nodeFace[2][0] = 0;
	nodeFace[2][1] = 3;
	nodeFace[2][2] = 4;
	nodeFace[3][0] = 1;
	nodeFace[3][1] = 3;
	nodeFace[3][2] = 4;
	nodeFace[4][0] = 0;
	nodeFace[4][1] = 2;
	nodeFace[4][2] = 5;
	nodeFace[5][0] = 1;
	nodeFace[5][1] = 2;
	nodeFace[5][2] = 5;
	nodeFace[6][0] = 0;
	nodeFace[6][1] = 3;
	nodeFace[6][2] = 5;
	nodeFace[7][0] = 1;
	nodeFace[7][1] = 3;
	nodeFace[7][2] = 5;

	nodeEdge[0][0] = 0;
	nodeEdge[0][1] = 2;
	nodeEdge[0][2] = 4;
	nodeEdge[1][0] = 1;
	nodeEdge[1][1] = 2;
	nodeEdge[1][2] = 5;
	nodeEdge[2][0] = 0;
	nodeEdge[2][1] = 3;
	nodeEdge[2][2] = 6;
	nodeEdge[3][0] = 1;
	nodeEdge[3][1] = 3;
	nodeEdge[3][2] = 7;
	nodeEdge[4][0] = 4;
	nodeEdge[4][1] = 8;
	nodeEdge[4][2] = 10;
	nodeEdge[5][0] = 5;
	nodeEdge[5][1] = 9;
	nodeEdge[5][2] = 10;
	nodeEdge[6][0] = 6;
	nodeEdge[6][1] = 8;
	nodeEdge[6][2] = 11;
	nodeEdge[7][0] = 7;
	nodeEdge[7][1] = 9;
	nodeEdge[7][2] = 11;

	edgeNode[0][0] = 0;
	edgeNode[0][1] = 2;
	edgeNode[1][0] = 1;
	edgeNode[1][1] = 3;
	edgeNode[2][0] = 0;
	edgeNode[2][1] = 1;
	edgeNode[3][0] = 2;
	edgeNode[3][1] = 3;
	edgeNode[4][0] = 0;
	edgeNode[4][1] = 4;
	edgeNode[5][0] = 1;
	edgeNode[5][1] = 5;
	edgeNode[6][0] = 2;
	edgeNode[6][1] = 6;
	edgeNode[7][0] = 3;
	edgeNode[7][1] = 7;
	edgeNode[8][0] = 4;
	edgeNode[8][1] = 6;
	edgeNode[9][0] = 5;
	edgeNode[9][1] = 7;
	edgeNode[10][0] = 4;
	edgeNode[10][1] = 5;
	edgeNode[11][0] = 6;
	edgeNode[11][1] = 7;

	faceNode[0][0] = 0;
	faceNode[0][1] = 2;
	faceNode[0][2] = 4;
	faceNode[0][3] = 6;
	faceNode[1][0] = 1;
	faceNode[1][1] = 3;
	faceNode[1][2] = 5;
	faceNode[1][3] = 7;
	faceNode[2][0] = 0;
	faceNode[2][1] = 1;
	faceNode[2][2] = 4;
	faceNode[2][3] = 5;
	faceNode[3][0] = 2;
	faceNode[3][1] = 3;
	faceNode[3][2] = 6;
	faceNode[3][3] = 7;
	faceNode[4][0] = 0;
	faceNode[4][1] = 1;
	faceNode[4][2] = 2;
	faceNode[4][3] = 3;
	faceNode[5][0] = 4;
	faceNode[5][1] = 5;
	faceNode[5][2] = 6;
	faceNode[5][3] = 7;

	edgeFace[0][0] = 0;
	edgeFace[0][1] = 4;
	edgeFace[1][0] = 1;
	edgeFace[1][1] = 4;
	edgeFace[2][0] = 2;
	edgeFace[2][1] = 4;
	edgeFace[3][0] = 3;
	edgeFace[3][1] = 4;
	edgeFace[4][0] = 0;
	edgeFace[4][1] = 2;
	edgeFace[5][0] = 1;
	edgeFace[5][1] = 2;
	edgeFace[6][0] = 0;
	edgeFace[6][1] = 3;
	edgeFace[7][0] = 1;
	edgeFace[7][1] = 3;
	edgeFace[8][0] = 0;
	edgeFace[8][1] = 5;
	edgeFace[9][0] = 1;
	edgeFace[9][1] = 5;
	edgeFace[10][0] = 2;
	edgeFace[10][1] = 5;
	edgeFace[11][0] = 3;
	edgeFace[11][1] = 5;

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

	edgeCoeffs[0][0] = -1;
	edgeCoeffs[0][1] =  0;
	edgeCoeffs[0][2] = -1;
	edgeCoeffs[1][0] =  1;
	edgeCoeffs[1][1] =  0;
	edgeCoeffs[1][2] = -1;
	edgeCoeffs[2][0] =  0;
	edgeCoeffs[2][1] = -1;
	edgeCoeffs[2][2] = -1;
	edgeCoeffs[3][0] =  0;
	edgeCoeffs[3][1] =  1;
	edgeCoeffs[3][2] = -1;
	edgeCoeffs[4][0] = -1;
	edgeCoeffs[4][1] = -1;
	edgeCoeffs[4][2] =  0;
	edgeCoeffs[5][0] =  1;
	edgeCoeffs[5][1] = -1;
	edgeCoeffs[5][2] =  0;
	edgeCoeffs[6][0] = -1;
	edgeCoeffs[6][1] =  1;
	edgeCoeffs[6][2] =  0;
	edgeCoeffs[7][0] =  1;
	edgeCoeffs[7][1] =  1;
	edgeCoeffs[7][2] =  0;
	edgeCoeffs[8][0] = -1;
	edgeCoeffs[8][1] =  0;
	edgeCoeffs[8][2] =  1;
	edgeCoeffs[9][0] =  1;
	edgeCoeffs[9][1] =  0;
	edgeCoeffs[9][2] =  1;
	edgeCoeffs[10][0] =  0;
	edgeCoeffs[10][1] = -1;
	edgeCoeffs[10][2] =  1;
	edgeCoeffs[11][0] =  0;
	edgeCoeffs[11][1] =  1;
	edgeCoeffs[11][2] =  1;

	nodeCoeffs[0][0] = -1;
	nodeCoeffs[0][1] = -1;
	nodeCoeffs[0][2] = -1;
	nodeCoeffs[1][0] =  1;
	nodeCoeffs[1][1] = -1;
	nodeCoeffs[1][2] = -1;
	nodeCoeffs[2][0] = -1;
	nodeCoeffs[2][1] =  1;
	nodeCoeffs[2][2] = -1;
	nodeCoeffs[3][0] =  1;
	nodeCoeffs[3][1] =  1;
	nodeCoeffs[3][2] = -1;
	nodeCoeffs[4][0] = -1;
	nodeCoeffs[4][1] = -1;
	nodeCoeffs[4][2] =  1;
	nodeCoeffs[5][0] =  1;
	nodeCoeffs[5][1] = -1;
	nodeCoeffs[5][2] =  1;
	nodeCoeffs[6][0] = -1;
	nodeCoeffs[6][1] =  1;
	nodeCoeffs[6][2] =  1;
	nodeCoeffs[7][0] =  1;
	nodeCoeffs[7][1] =  1;
	nodeCoeffs[7][2] =  1;

	parallelEdges[0][0] = 1;
	parallelEdges[0][1] = 8;
	parallelEdges[0][2] = 9;
	parallelEdges[1][0] = 0;
	parallelEdges[1][1] = 8;
	parallelEdges[1][2] = 9;
	parallelEdges[2][0] = 3;
	parallelEdges[2][1] = 10;
	parallelEdges[2][2] = 11;
	parallelEdges[3][0] = 2;
	parallelEdges[3][1] = 10;
	parallelEdges[3][2] = 11;
	parallelEdges[4][0] = 5;
	parallelEdges[4][1] = 6;
	parallelEdges[4][2] = 7;
	parallelEdges[5][0] = 4;
	parallelEdges[5][1] = 6;
	parallelEdges[5][2] = 7;
	parallelEdges[6][0] = 4;
	parallelEdges[6][1] = 5;
	parallelEdges[6][2] = 7;
	parallelEdges[7][0] = 4;
	parallelEdges[7][1] = 5;
	parallelEdges[7][2] = 6;
	parallelEdges[8][0] = 0;
	parallelEdges[8][1] = 1;
	parallelEdges[8][2] = 9;
	parallelEdges[9][0] = 0;
	parallelEdges[9][1] = 1;
	parallelEdges[9][2] = 8;
	parallelEdges[10][0] = 2;
	parallelEdges[10][1] = 3;
	parallelEdges[10][2] = 11;
	parallelEdges[11][0] = 2;
	parallelEdges[11][1] = 3;
	parallelEdges[11][2] = 10;

	for (int level = 0; level < MAX_LEVEL; ++level) {
		lengths[level] = uint32_t(1) << (MAX_LEVEL - level);
		areas[level]   = uint64_t(1) << ((dim - 1) * (MAX_LEVEL - level));
		volumes[level] = uint64_t(1) << (dim * (MAX_LEVEL - level));
	}

}

// =================================================================================== //

}
