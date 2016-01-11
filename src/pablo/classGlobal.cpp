#include "classGlobal.hpp"

// =================================================================================== //

uint8_t
classGlobal::getBoolBytes()  {
	return boolBytes;
}

void
classGlobal::getEdgecoeffs(int8_t edgecoeffs_[12][3])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<3; j++){
		edgecoeffs_[i][j] = edgecoeffs[i][j];
		}
	}
}

void
classGlobal::getEdgeface(uint8_t edgeface_[12][2])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<2; j++){
			edgeface_[i][j] = edgeface[i][j];
		}
	}
}

void
classGlobal::getFacenode(uint8_t facenode_[6][4])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<4; j++){
			facenode_[i][j] = facenode[i][j];
		}
	}
}

uint8_t
classGlobal::getGlobalIndexBytes()  {
	return globalIndexBytes;
}

uint8_t
classGlobal::getLevelBytes()  {
	return levelBytes;
}

uint8_t
classGlobal::getMarkerBytes()  {
	return markerBytes;
}

uint32_t
classGlobal::getMaxLength()  {
	return max_length;
}

int8_t
classGlobal::getMaxLevel()  {
	return MAX_LEVEL;
}

uint8_t
classGlobal::getNchildren()  {
	return nchildren;
}

uint8_t
classGlobal::getNedges()  {
	return nedges;
}

uint8_t
classGlobal::getNfaces()  {
	return nfaces;
}

uint8_t
classGlobal::getNnodes()  {
	return nnodes;
}

uint8_t
classGlobal::getNnodesperface()  {
	return nnodesperface;
}

void
classGlobal::getNodecoeffs(int8_t nodecoeffs_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodecoeffs_[i][j] = nodecoeffs[i][j];
		}
	}
}

void
classGlobal::getNodeface(uint8_t nodeface_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodeface_[i][j] = nodeface[i][j];
		}
	}
}

void
classGlobal::getNormals(int8_t normals_[6][3])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			normals_[i][j] = normals[i][j];
		}
	}
}

uint8_t
classGlobal::getOctantBytes()  {
	return octantBytes;
}

void
classGlobal::getOppface(uint8_t oppface_[6])  {
	for (int j=0; j<6; j++){
		oppface_[j] = oppface[j];
	}
}

/*! Initialization of static global variables.
 * \param[in] maxlevel Space dimension; only value equal to 2 or 3 are accepted.
 * \param[in] dim Space dimension; only value equal to 2 or 3 are accepted.
 */
void
classGlobal::setGlobal(int8_t maxlevel, uint8_t dim){

	if (dim>3) dim = 3;
	if (dim<2) dim = 2;

	MAX_LEVEL = maxlevel;

	max_length 			= uint32_t(1<<MAX_LEVEL);
	nchildren 			= 1<<dim;
	nfaces 				= 2*dim;
	nnodes				= 1<<dim;
	nedges 				= (dim-2)*12;
	nnodesperface 		= 1<<(dim-1);
	octantBytes 		= uint8_t(sizeof(uint32_t)*3 + sizeof(uint8_t) + sizeof(int8_t) + (17)*sizeof(bool));
	globalIndexBytes 	= uint8_t(sizeof(uint64_t));
	markerBytes 		= sizeof(int8_t);
	levelBytes 			= sizeof(uint8_t);
	boolBytes 			= sizeof(bool);

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

	nodecoeffs[0][0] = -1;
	nodecoeffs[0][1] = -1;
	nodecoeffs[0][2] = -1;
	nodecoeffs[1][0] =  1;
	nodecoeffs[1][1] = -1;
	nodecoeffs[1][2] = -1;
	nodecoeffs[2][0] = -1;
	nodecoeffs[2][1] =  1;
	nodecoeffs[2][2] = -1;
	nodecoeffs[3][0] =  1;
	nodecoeffs[3][1] =  1;
	nodecoeffs[3][2] = -1;
	nodecoeffs[4][0] = -1;
	nodecoeffs[4][1] = -1;
	nodecoeffs[4][2] =  1;
	nodecoeffs[5][0] =  1;
	nodecoeffs[5][1] = -1;
	nodecoeffs[5][2] =  1;
	nodecoeffs[6][0] = -1;
	nodecoeffs[6][1] =  1;
	nodecoeffs[6][2] =  1;
	nodecoeffs[7][0] =  1;
	nodecoeffs[7][1] =  1;
	nodecoeffs[7][2] =  1;

}




// =================================================================================== //



