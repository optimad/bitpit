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

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Global.hpp"
#include "Octant.hpp"
#include "bitpit_common.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>


/*!
 * Input stream operator for class Octant. Stream cell data from memory
 * input stream to container.
 *
 * \param[in] buffer is the input stream from memory
 * \param[in] octant is the octant object
 * \result Returns the same input stream received in input.
 */
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::Octant &octant)
{
    uint8_t dimensions;
    buffer >> dimensions;

    uint8_t level;
    buffer >> level;

    octant.initialize(dimensions, level, true);

    buffer >> octant.m_x;
    buffer >> octant.m_y;
    buffer >> octant.m_z;

    buffer >> octant.m_marker;

    buffer >> octant.m_ghost;

    for(int i = 0; i < bitpit::Octant::INFO_ITEM_COUNT; ++i){
        bool value;
        buffer >> value;
        octant.m_info[i] = value;
    }

    return buffer;
}

/*!
 * Output stream operator for class Cell. Stream octant data from container
 * to output stream.
 *
 * \param[in] buffer is the output stream from memory
 * \param[in] octant is the octant object
 * \result Returns the same output stream received in input.
 */
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const bitpit::Octant &octant)
{
    buffer << octant.m_dim;
    buffer << octant.m_level;

    buffer << octant.m_x;
    buffer << octant.m_y;
    buffer << octant.m_z;

    buffer << octant.m_marker;

    buffer << octant.m_ghost;

    for(int i = 0; i < bitpit::Octant::INFO_ITEM_COUNT; ++i){
        buffer << (bool) octant.m_info[i];
    }

    return buffer;
}

namespace bitpit {

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

// =================================================================================== //
// STATIC AND CONSTANT
// =================================================================================== //

constexpr int Octant::sm_CoeffNode[8][3];
constexpr int Octant::sm_CoeffFaceCenter[6][3];
constexpr int Octant::sm_CoeffEdgeCenter[12][3];

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*! Create a dummy octant.
 */
Octant::Octant(){
	initialize();
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D zero-level octant with origin in (0,0,0).
 * \param[in] dim Dimension of octant (2/3 for 2D/3D octant).
 */
Octant::Octant(uint8_t dim){
	initialize(dim, 0, true);
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D octant with user defined origin and level.
 * \param[in] dim Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level Refinement level of octant (0 for root octant).
 * \param[in] x X-coordinates of the origin of the octant.
 * \param[in] y Y-coordinates of the origin of the octant.
 * \param[in] z Z-Coordinates of the origin of the octant (The default value is 0).
 */
Octant::Octant(uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z){
	initialize(dim, level, true);

	// Set the coordinates
	m_x = x;
	m_y = y;
	m_z = (m_dim-2) * z;
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D octant with user defined origin, level and boundary conditions.
 * \param[in] bound Boundary condition for the faces of the octant (the same for each face).
 * \param[in] dim Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level Refinement level of octant (0 for root octant).
 * \param[in] x X-coordinates of the origin of the octant.
 * \param[in] y Y-coordinates of the origin of the octant.
 * \param[in] z Z-Coordinates of the origin of the octant (The default value is 0).
 */
Octant::Octant(bool bound, uint8_t dim, uint8_t level, int32_t x, int32_t y, int32_t z){
	initialize(dim, level, bound);

	// Set the coordinates
	m_x = x;
	m_y = y;
	m_z = (m_dim-2) * z;
};

/*! Copy constructor of an octant.
 */
Octant::Octant(const Octant &octant){
	m_dim = octant.m_dim;
	m_x = octant.m_x;
	m_y = octant.m_y;
	m_z = octant.m_z;
	m_level = octant.m_level;
	m_marker = octant.m_marker;
	m_info = octant.m_info;
	m_ghost = octant.m_ghost;
};

/*! Check if two octants are equal (no check on info)
 */
bool Octant::operator ==(const Octant & oct2){
	bool check = true;
	check = check && (m_dim == oct2.m_dim);
	check = check && (m_x == oct2.m_x);
	check = check && (m_y == oct2.m_y);
	check = check && (m_z == oct2.m_z);
	check = check && (m_level == oct2.m_level);
	return check;
}

// =================================================================================== //
// METHODS
// =================================================================================== //

/*! Initialize a dummy octant.
 */
void
Octant::initialize() {
	initialize(0, 0, false);
}

/*! Initialize the octant.
 * \param[in] dim Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level Refinement level of octant (0 for root octant).
 * \param[in] bound Boundary condition for the faces of the octant (the same for each face).
 */
void
Octant::initialize(uint8_t dim, uint8_t level, bool bound) {
	m_dim   = dim;
	m_level = level;

	// Reset the marker
	m_marker = 0;

	// Set the coordinates
	m_x = 0;
	m_y = 0;
	m_z = 0;

	// Initialize octant info
	m_info.reset();
	m_info[OctantInfo::INFO_BALANCED] = true;
	m_ghost = -1;

	// If this is the root octant we need to set the boundary condition bound
	// for faces
	if (m_dim >= 2 && m_level == 0) {
		uint8_t nf = m_dim*2;
		for (uint8_t i=0; i<nf; i++){
			m_info[i] = bound;
		}
	}
};

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //

/*! Get the dimension of an octant, 2 or 3 if 2D or 3D octant.
 * \return Dimension of octant.
 */
uint32_t
Octant::getDim() const{return m_dim;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinates of node 0.
 */
u32array3
Octant::getCoordinates() const{
	u32array3 xx;
	xx[0] = m_x;
	xx[1] = m_y;
	xx[2] = m_z;
	return xx;
};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate X of node 0.
 */
uint32_t
Octant::getX() const{return m_x;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Y of node 0.
 */
uint32_t
Octant::getY() const{return m_y;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Z of node 0.
 */
uint32_t
Octant::getZ() const{return m_z;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinates of node 0.
 */
u32array3
Octant::getCoord() const {
	u32array3 coord;
	coord[0] = m_x;
	coord[1] = m_y;
	coord[2] = m_z;
	return coord;
};

/*! Get the level of an octant.
 * \return Level of octant.
 */
uint8_t
Octant::getLevel() const{return m_level;};

/*! Get the refinement marker of an octant.
 * \return Marker of octant.
 */
int8_t
Octant::getMarker() const{return m_marker;};

/*! Get the bound flag on an octant face.
 * \param[in] face local index of the face.
 * \return true if the face face is a boundary face.
 */
bool
Octant::getBound(uint8_t face) const{
	return m_info[face];
};

/*! Get the bound flag on an octant.
 * \return true if the octant is a boundary octant.
 */
bool
Octant::getBound() const{
	return m_info[OctantInfo::INFO_BOUNDFACE0]||m_info[OctantInfo::INFO_BOUNDFACE1]||
	        m_info[OctantInfo::INFO_BOUNDFACE2]||m_info[OctantInfo::INFO_BOUNDFACE3]||
	        ((m_dim-2)&&(m_info[OctantInfo::INFO_BOUNDFACE4]||m_info[OctantInfo::INFO_BOUNDFACE5]));
};

/*! Set the boundary flag to true on an octant face.
 * \param[in] face local index of the boundary face.
 */
void
Octant::setBound(uint8_t face) {
	m_info[face] = true;
};

/*! Get the pbound flag on an octant face.
 * \param[in] face local index of the face.
 * \return true if the face-th face is a process boundary face.
 */
bool
Octant::getPbound(uint8_t face) const{
	return m_info[6+face];
};

/*! Get the pbound flag on an octant face.
 * \return true if the octant has at least one face that is a process boundary face.
 */
bool
Octant::getPbound() const{
	return m_info[OctantInfo::INFO_PBOUNDFACE0]||m_info[OctantInfo::INFO_PBOUNDFACE1]||
	        m_info[OctantInfo::INFO_PBOUNDFACE2]||m_info[OctantInfo::INFO_PBOUNDFACE3]||
	        ((m_dim-2)*(m_info[OctantInfo::INFO_PBOUNDFACE4]||m_info[OctantInfo::INFO_PBOUNDFACE5]));
};

/*! Get if the octant is new after a refinement.
 * \return true if the the octant is new after a refinement.
 */
bool
Octant::getIsNewR() const{return m_info[OctantInfo::INFO_NEW4REFINEMENT];};

/*! Get if the octant is new after a coarsening.
 * \return true if the the octant is new after a coarsening.
 */
bool
Octant::getIsNewC() const{return m_info[OctantInfo::INFO_NEW4COARSENING];};

/*! Get if the octant is a scary ghost octant.
 * \return true if the octant is a ghost octant.
 */
bool
Octant::getIsGhost() const{return (m_ghost >= 0);};

/*! Get the layer number of the ghost halo an octant belong to.
 * \return the layer in the ghost halo. 0 is for internal (non-ghost) octant.
 */
int
Octant::getGhostLayer() const{return m_ghost;};

/*! Get if the octant has to be balanced.
 * \return true if the octant has to be balanced.
 */
bool
Octant::getBalance() const{return (m_info[OctantInfo::INFO_BALANCED]);};

/*! Set the refinement marker of an octant.
 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
 */
void
Octant::setMarker(int8_t marker){
	if (marker != m_marker)
		m_info[OctantInfo::INFO_AUX] = true;
	this->m_marker = marker;
};

/*! Set the balancing condition of an octant.
 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
 */
void
Octant::setBalance(bool balance){
	if (balance != m_info[OctantInfo::INFO_BALANCED])
		m_info[OctantInfo::INFO_AUX] = true;
	m_info[OctantInfo::INFO_BALANCED] = balance;
};

/*! Set the level of an octant.
 * \param[in] level New level of the octant.
 */
void
Octant::setLevel(uint8_t level){
	this->m_level = level;
};

/*! Set the process boundary condition of an octant.
 * \param[in] face Index of the target face to set.
 * \param[in] flag Condition true/false if the target face is /is not a process boundary face.
 */
void
Octant::setPbound(uint8_t face, bool flag){
	m_info[6+face] = flag;
};

/*! Set the ghost specifier of an octant.
 * \param[in] ghostLayer Number defining the layer of the ghost halo. 0 is for internal octants, i.e. non-ghost.
 * ghostLayer > 0 means that the octant is ghost in the ghostLayer layer of the ghost halo.
 */
void
Octant::setGhostLayer(int ghostLayer){
    m_ghost = ghostLayer;
};


// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

/*! Get the size of an octant in logical domain, i.e. the side length.
 * \return Size of octant.
 */
uint32_t
Octant::getSize() const{
	uint32_t size = uint32_t(1) << (Global::getMaxLevel() - m_level);
	return size;
};

/*! Get the area of an octant in logical domain .
 * \return Area of octant.
 */
uint64_t
Octant::getArea() const{
	uint64_t area = uint64_t(1) << ((m_dim-1) * (Global::getMaxLevel() - m_level));
	return area;
};

/*! Get the volume of an octant in logical domain.
 * \return Volume of octant.
 */
uint64_t
Octant::getVolume() const{
	uint64_t volume = uint64_t(1) << (m_dim * (Global::getMaxLevel() - m_level));
	return volume;
};

// =================================================================================== //

/*! Get the coordinates of the center of an octant in logical domain.
 * \return Array[3] with the coordinates of the center of octant.
 */
darray3
Octant::getCenter() const{
	double	dh;
	darray3 center;

	dh = double(getSize())*0.5;
	center[0] = (double)m_x + dh;
	center[1] = (double)m_y + dh;
	center[2] = (double)m_z + double(m_dim-2)*dh;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a face of an octant in logical domain.
 * \param[in] iface Local index of the face
 * \return Array[3] with the coordinates of the center of the octant face.
 */
darray3
Octant::getFaceCenter(uint8_t iface) const{
	double	dh_2;
	darray3 center;

	assert(iface < m_dim*2);

	dh_2 = double(getSize())*0.5;
	center[0] = (double)m_x + (double)sm_CoeffFaceCenter[iface][0] * dh_2;
	center[1] = (double)m_y + (double)sm_CoeffFaceCenter[iface][1] * dh_2;
	center[2] = (double)m_z + double(m_dim-2) * (double)sm_CoeffFaceCenter[iface][2] * dh_2;

	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a edge of an octant in logical domain.
 * \param[in] iedge Local index of the edge
 * \return Array[3] with the coordinates of the center of the octant edge.
 */
darray3
Octant::getEdgeCenter(uint8_t iedge) const{
	double	dh_2;
	darray3 center;

	dh_2 = double(getSize())*0.5;
	center[0] = (double)m_x + (double)sm_CoeffEdgeCenter[iedge][0] * dh_2;
	center[1] = (double)m_y + (double)sm_CoeffEdgeCenter[iedge][1] * dh_2;
	center[2] = (double)m_z + double(m_dim-2) * (double)sm_CoeffEdgeCenter[iedge][2] * dh_2;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \param[out] nodes Vector of arrays [nnodes][3] with the coordinates of the nodes of octant.
 */
void
Octant::getNodes(u32arr3vector & nodes) const{
	uint8_t		i;
	uint32_t	dh;
	uint8_t nn = uint8_t(1)<<m_dim;

	dh = getSize();
	nodes.resize(nn);

	for (i = 0; i < nn; i++){
		nodes[i][0] = m_x + uint32_t(sm_CoeffNode[i][0])*dh;
		nodes[i][1] = m_y + uint32_t(sm_CoeffNode[i][1])*dh;
		nodes[i][2] = m_z + uint32_t(sm_CoeffNode[i][2])*dh;
	}
};

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \return Vector of arrays [nnodes][3] with the coordinates of the nodes of octant.
 */
u32arr3vector
Octant::getNodes() const{
	uint8_t		i;
	uint32_t	dh;
	uint8_t nn = uint8_t(1)<<m_dim;
	u32arr3vector nodes;

	dh = getSize();
	nodes.resize(nn);

	for (i = 0; i < nn; i++){
		nodes[i][0] = m_x + sm_CoeffNode[i][0]*dh;
		nodes[i][1] = m_y + sm_CoeffNode[i][1]*dh;
		nodes[i][2] = m_z + sm_CoeffNode[i][2]*dh;
	}

	return nodes;
};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[out] node Array[3] with the logical coordinates of the node of the octant.
 */
void		Octant::getNode(u32array3 & node, uint8_t inode) const{
	uint32_t	dh;

	dh = getSize();
	node[0] = m_x + sm_CoeffNode[inode][0]*dh;
	node[1] = m_y + sm_CoeffNode[inode][1]*dh;
	node[2] = m_z + sm_CoeffNode[inode][2]*dh;

};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \return Array[3] with the logical coordinates of the node of the octant.
 */
u32array3		Octant::getNode(uint8_t inode) const{
	u32array3 	node;
	uint32_t	dh;

	dh = getSize();
	node[0] = m_x + sm_CoeffNode[inode][0]*dh;
	node[1] = m_y + sm_CoeffNode[inode][1]*dh;
	node[2] = m_z + sm_CoeffNode[inode][2]*dh;
	return node;
};

/*! Get the normal of a face of an octant in logical domain.
 * \param[in] iface Index of the face for normal computing.
 * \param[out] normal Array[3] with components (with z=0) of the normal of face.
 * \param[in] normals Global structure with components of face normals of a reference octant.
 */
void		Octant::getNormal(uint8_t iface, i8array3 & normal, const int8_t (&normals)[6][3]) const{
	uint8_t		i;
	for (i = 0; i < 3; i++){
		normal[i] = normals[iface][i];
	}
};


/** Compute the Morton index of the octant (without level).
 * \return morton Morton index of the octant.
 */
uint64_t	Octant::computeMorton() const{
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->m_x,this->m_y,this->m_z);
	return morton;
};

/** Compute the Morton index of the given node (without level).
 * \return morton Morton index of the node.
 */
uint64_t	Octant::computeNodeMorton(uint8_t inode) const{

	u32array3 node = getNode(inode);

	return keyXYZ(node[0], node[1], node[2], Global::getMaxLevel());
};

/** Get the size of the buffer required to communicate the octant.
 * \return Returns the buffer size (in bytes).
 */
unsigned int Octant::getBinarySize()
{
    unsigned int binarySize = 0;
    binarySize += sizeof(uint8_t); // dimensions
    binarySize += sizeof(uint8_t); // level
    binarySize += 3 * sizeof(uint32_t); // 3 coordinates
    binarySize += sizeof(int8_t); // marker
    binarySize += sizeof(int); // ghost layer
    binarySize += INFO_ITEM_COUNT * sizeof(bool); // info

    return binarySize;
}

// =================================================================================== //
// OTHER METHODS
// =================================================================================== //

/** Build the last descendant octant of this octant.
 * \return Last descendant octant.
 */
Octant	Octant::buildLastDesc() const {
	u32array3 delta = { {0,0,0} };
	for (int i=0; i<m_dim; i++){
		delta[i] = (uint32_t(1) << (Global::getMaxLevel() - m_level)) - 1;
	}
	Octant last_desc(m_dim, Global::getMaxLevel(), (m_x+delta[0]), (m_y+delta[1]), (m_z+delta[2]));
	return last_desc;
};

// =================================================================================== //

/** Build the father octant of this octant.
 * \return Father octant.
 */
Octant	Octant::buildFather() const {
	uint32_t delta[3];
	uint32_t xx[3];
	xx[0] = m_x;
	xx[1] = m_y;
	xx[2] = m_z;
	delta[2] = 0;
	for (int i=0; i<m_dim; i++){
		delta[i] = xx[i]%(uint32_t(1) << (Global::getMaxLevel() - max(0,(m_level-1))));
	}
	Octant father(m_dim, max(0,m_level-1), m_x-delta[0], m_y-delta[1], m_z-delta[2]);
	return father;
};

// =================================================================================== //

/** Builds children of octant.
 *   \return Ordered (by Z-index) vector of children[nchildren] (info update)
 */
vector< Octant >	Octant::buildChildren() const {
	uint8_t xf,yf,zf;
	int nchildren = 1<<m_dim;

	if (this->m_level < Global::getMaxLevel()){
		vector< Octant > children(nchildren, Octant(m_dim));
		for (int i=0; i<nchildren; i++){
			switch (i) {
			case 0 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=5;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[0] = oct;
			}
			break;
			case 1 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_x += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=5;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[1] = oct;
			}
			break;
			case 2 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_y += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=5;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[2] = oct;
			}
			break;
			case 3 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_x += dh;
				oct.m_y += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=5;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[3] = oct;
			}
			break;
			case 4 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_z += dh;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=4;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[4] = oct;
			}
			break;
			case 5 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_x += dh;
				oct.m_z += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=4;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[5] = oct;
			}
			break;
			case 6 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_y += dh;
				oct.m_z += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=4;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[6] = oct;
			}
			break;
			case 7 :
			{
				Octant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[OctantInfo::INFO_NEW4REFINEMENT]=true;
				uint32_t dh = oct.getSize();
				oct.m_x += dh;
				oct.m_y += dh;
				oct.m_z += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=4;
				oct.m_info[xf] = oct.m_info[xf+6] = false;
				oct.m_info[yf] = oct.m_info[yf+6] = false;
				oct.m_info[zf] = oct.m_info[zf+6] = false;
				children[7] = oct;
			}
			break;
			}
		}
		return children;
	}
	else{
		vector< Octant > children(0, Octant(m_dim));
		return children;
	}
};

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant
 * throught face iface (sizehf=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[out] sizehf Number of possible neighbours.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> Octant::computeHalfSizeMorton(uint8_t iface, uint32_t & sizehf) const {
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i,cx,cy,cz;
	int nchildren = 1<<m_dim;

	nneigh = (m_level < Global::getMaxLevel()) ? nchildren/2 : 1;
	dh = (m_level < Global::getMaxLevel()) ? getSize()/2 : getSize();
	dh2 = getSize();

	if (m_info[iface]){
		sizehf = 0;
		vector<uint64_t> Morton(0);
		return Morton;
	}
	else{
		vector<uint64_t> Morton(nneigh);
		switch (iface) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i==1)||(i==3);
				cz = (m_dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->m_x-dh,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i==1)||(i==3);
				cz = (m_dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (m_dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y-dh,this->m_z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (m_dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2,this->m_z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2);
			}
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}


};

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant
 * throught face iface (sizem=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizem Number of possible neighbours.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> Octant::computeMinSizeMorton(uint8_t iface, const uint8_t & maxdepth, uint32_t & sizem) const {
	uint32_t dh,dh2;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;

	nneigh = (m_level < Global::getMaxLevel()) ? uint32_t(1)<<((m_dim-1)*(maxdepth-m_level)) : 1;
	dh = (m_level < Global::getMaxLevel()) ? uint32_t(1)<<(Global::getMaxLevel() - maxdepth) : getSize();
	dh2 = getSize();
	nline = uint32_t(1)<<(maxdepth-m_level);

	if (m_info[iface]){
		sizem = 0;
		vector<uint64_t> Morton(0);
		return Morton;
	}
	else{
		vector<uint64_t> Morton(nneigh);
		switch (iface) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cz = (m_dim-2)*(i%nline);
				cy = (m_dim==2)*(i%nline) + (m_dim-2)*(i/nline);
				Morton[i] = mortonEncode_magicbits(this->m_x-dh,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cz = (m_dim-2)*(i%nline);
				cy = (m_dim==2)*(i%nline) + (m_dim-2)*(i/nline);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cz = (m_dim-2)*(i%nline);
				cx = (m_dim==2)*(i%nline) + (m_dim-2)*(i/nline);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y-dh,this->m_z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cz = (m_dim-2)*(i%nline);
				cx = (m_dim==2)*(i%nline) + (m_dim-2)*(i/nline);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2,this->m_z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2);
			}
		}
		break;
		}
		sizem = nneigh;
		sort(Morton.begin(), Morton.end());
		return Morton;
	}

};

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught iface. Checks if balanced or not and uses half-size or min-size method
 * (sizeneigh=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizeneigh Number of possible neighbours.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> Octant::computeVirtualMorton(uint8_t iface, const uint8_t & maxdepth, uint32_t & sizeneigh) const {
	vector<uint64_t> Morton;
	if (!getBalance()){
		return computeMinSizeMorton(iface,
				maxdepth,
				sizeneigh);
	}
	else{
		return computeHalfSizeMorton(iface,
				sizeneigh);
	}
};

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * edge iedge (sizehf=0 if boundary octant)
 * \param[in] iedge Local index of the edge target.
 * \param[out] sizehf Number of possible neighbours.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> Octant::computeEdgeHalfSizeMorton(uint8_t iedge, uint32_t & sizehf, uint8_t (&edgeface)[12][2]) const {
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i;
	int32_t cx,cy,cz;
	uint8_t iface1, iface2;

	nneigh = (m_level < Global::getMaxLevel()) ? 2 : 1;
	dh = (m_level < Global::getMaxLevel()) ? getSize()/2 : getSize();
	dh2 = getSize();
	iface1 = edgeface[iedge][0];
	iface2 = edgeface[iedge][1];

	if (m_info[iface1] || m_info[iface2]){
		sizehf = 0;
		vector<uint64_t> Morton(0);
		return Morton;
	}
	else{
		vector<uint64_t> Morton(nneigh);
		switch (iedge) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = (i==1);
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = (i==1);
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = -1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = 1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = -1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = -1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 6 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = 1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 7 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = 1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 8 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = (i==1);
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 9 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = (i==1);
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 10 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = -1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 11 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = 1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh2*cz);
			}
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}


};

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * edge iedge (sizem=0 if boundary octant)
 * \param[in] iedge Local index of the edge target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizem Number of possible neighbours.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> 		Octant::computeEdgeMinSizeMorton(uint8_t iedge, const uint8_t & maxdepth, uint32_t & sizem, uint8_t (&edgeface)[12][2]) const {
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i;
	int32_t cx,cy,cz;
	uint8_t iface1, iface2;


	nneigh = (m_level < Global::getMaxLevel()) ? uint32_t(1)<<(maxdepth-m_level) : 1;
	dh = (m_level < Global::getMaxLevel()) ? uint32_t(1)<<(Global::getMaxLevel() - maxdepth) : getSize();
	dh2 = getSize();
	iface1 = edgeface[iedge][0];
	iface2 = edgeface[iedge][1];

	if (m_info[iface1] || m_info[iface2]){
		sizem = 0;
		vector<uint64_t> Morton(0);
		return Morton;
	}
	else{
		vector<uint64_t> Morton(nneigh);
		switch (iedge) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = i;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = i;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = -1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = 1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = -1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = -1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 6 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = 1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 7 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = 1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
			}
		}
		break;
		case 8 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = i;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 9 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = i;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 10 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = -1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
			}
		}
		break;
		case 11 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = 1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh2*cz);
			}
		}
		break;
		}
		sizem = nneigh;
		sort(Morton.begin(),Morton.end());
		return Morton;
	}
};

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught iedge. Checks if balanced or not and uses half-size or min-size method
 * (sizeneigh=0 if boundary octant).
 * \param[in] iedge Local index of the edge target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizeneigh Number of possible neighbours.
 * \param[in] balance_codim Maximum codimension of the entity for 2:1 balancing.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t>		Octant::computeEdgeVirtualMorton(uint8_t iedge, const uint8_t & maxdepth, uint32_t & sizeneigh, uint8_t balance_codim, uint8_t (&edgeface)[12][2]) const {

	if(getBalance() && balance_codim > 1){
		return computeEdgeHalfSizeMorton(iedge,
				sizeneigh, edgeface);
	}
	else{
		return computeEdgeMinSizeMorton(iedge,
				maxdepth, sizeneigh, edgeface);
	}
};

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * node inode (sizem=0 if boundary octant)
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizem Number of possible neighbours (1 or 0).
 * \param[in] nodeface Local node-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
uint64_t 		Octant::computeNodeMinSizeMorton(uint8_t inode, const uint8_t & maxdepth,uint32_t & sizem, uint8_t (&nodeface)[8][3]) const {

	uint32_t dh,dh2;
	uint32_t nneigh;
	int8_t cx,cy,cz;

	for (int i=0; i<m_dim; i++){
		uint8_t iface = nodeface[inode][i];
		if (m_info[iface]) {
			sizem = 0;
			return this->computeMorton();
		}
	}

	nneigh = 1;
	dh = (m_level < Global::getMaxLevel()) ? uint32_t(1)<<(Global::getMaxLevel() - maxdepth) : getSize();
	dh2 = getSize();

	uint64_t Morton;
	switch (inode) {
	case 0 :
	{
		cx = -1;
		cy = -1;
		cz = -1*(m_dim-2);
		Morton = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh*cz);
	}
	break;
	case 1 :
	{
		cx = 1;
		cy = -1;
		cz = -1*(m_dim-2);
		Morton = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh*cz);
	}
	break;
	case 2 :
	{
		cx = -1;
		cy = 1;
		cz = -1*(m_dim-2);
		Morton = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
	}
	break;
	case 3 :
	{
		cx = 1;
		cy = 1;
		cz = -1*(m_dim-2);
		Morton = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh2*cy,this->m_z+dh*cz);
	}
	break;
	case 4 :
	{
		cx = -1;
		cy = -1;
		cz = 1;
		Morton = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
	}
	break;
	case 5 :
	{
		cx = 1;
		cy = -1;
		cz = 1;
		Morton = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh*cy,this->m_z+dh2*cz);
	}
	break;
	case 6 :
	{
		cx = -1;
		cy = 1;
		cz = 1;
		Morton = mortonEncode_magicbits(this->m_x+dh*cx,this->m_y+dh2*cy,this->m_z+dh2*cz);
	}
	break;
	case 7 :
	{
		cx = 1;
		cy = 1;
		cz = 1;
		Morton = mortonEncode_magicbits(this->m_x+dh2*cx,this->m_y+dh2*cy,this->m_z+dh2*cz);
	}
	break;
	default:
		BITPIT_UNREACHABLE("The maximum number of nodes is 8.");
	}
	sizem = nneigh;
	return Morton;

};

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught inode. Uses min-size method (sizeneigh=0 if boundary octant).
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizeneigh Number of possible neighbours (1).
 * \param[in] nodeface Local node-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
uint64_t 		Octant::computeNodeVirtualMorton(uint8_t inode, const uint8_t & maxdepth, uint32_t & sizeneigh, uint8_t (&nodeface)[8][3]) const {

	return computeNodeMinSizeMorton(inode, maxdepth,
			sizeneigh, nodeface);

 };

/*! Computes Morton index (without level) of same size
 * periodic neighbour of octant throught face iface.
 * \param[in] iface Local index of the face target.
 * \return Periodic neighbour morton number.
 */
uint64_t Octant::computePeriodicMorton(uint8_t iface) const {
	uint64_t Morton;
	uint32_t dh;
	dh = getSize();
	uint32_t maxLength = uint32_t(1)<<Global::getMaxLevel();

	if (!m_info[iface]){
		return this->computeMorton();
	}
	else{
		switch (iface) {
		case 0 :
		{
			Morton = mortonEncode_magicbits(maxLength-dh,this->m_y,this->m_z);
		}
		break;
		case 1 :
		{
			Morton = mortonEncode_magicbits(0,this->m_y,this->m_z);
		}
		break;
		case 2 :
		{
			Morton = mortonEncode_magicbits(this->m_x,maxLength-dh,this->m_z);
		}
		break;
		case 3 :
		{
			Morton = mortonEncode_magicbits(this->m_x,0,this->m_z);
		}
		break;
		case 4 :
		{
			Morton = mortonEncode_magicbits(this->m_x,this->m_y,maxLength-dh);
		}
		break;
		case 5 :
		{
			Morton = mortonEncode_magicbits(this->m_x,this->m_y,0);
		}
		break;
		default:
			BITPIT_UNREACHABLE("The maximum number of faces is 6.");
		}
		return Morton;
	}
};

/** Build a same size periodic octant of this octant throught face iface.
 * \return Periodic octant of the same size (note: it is a stand-alone octant,
 * may be not living in octree).
 */
Octant Octant::computePeriodicOctant(uint8_t iface) const {
	Octant degOct(this->m_dim, this->m_level, this->m_x, this->m_y, this->m_z);
	uint32_t maxLength = uint32_t(1)<<Global::getMaxLevel();
	uint32_t dh = this->getSize();

	if (!m_info[iface]){
		return this->computeMorton();
	}
	else{
		switch (iface) {
		case 0 :
		{
			degOct.m_x = maxLength-dh;
		}
		break;
		case 1 :
		{
			degOct.m_x = 0;
		}
		break;
		case 2 :
		{
			degOct.m_y = maxLength-dh;
		}
		break;
		case 3 :
		{
			degOct.m_y = 0;
		}
		break;
		case 4 :
		{
			degOct.m_z = maxLength-dh;
		}
		break;
		case 5 :
		{
			degOct.m_y = 0;
		}
		break;
		}
		degOct.m_level = this->m_level;
		degOct.m_info = false;
		return degOct;
	}


};

/*! Get the coordinates of the octant shifted throught face iface and
 * near the opposite periodic boundary (i.e. the coordinates considering this octant
 * as a ghost periodic octant).
 * \param[in] iface Local index of the face target.
 * \return Coordinates of octant considered as periodic ghost out of the logical domain.
 */
array<int64_t,3> Octant::getPeriodicCoord(uint8_t iface) const {
	array<int64_t,3> coord;
	coord[0] = this->m_x;
	coord[1] = this->m_y;
	coord[2] = this->m_z;
	int64_t dh = this->getSize();
	int64_t maxLength = int64_t(1)<<Global::getMaxLevel();

	switch (iface) {
	case 0 :
	{
		coord[0] = maxLength;
	}
	break;
	case 1 :
	{
		coord[0]  = -dh;
	}
	break;
	case 2 :
	{
		coord[1]  = maxLength;
	}
	break;
	case 3 :
	{
		coord[1] = -dh;
	}
	break;
	case 4 :
	{
		coord[2] = maxLength;
	}
	break;
	case 5 :
	{
		coord[2] = -dh;
	}
	break;
	}
	return coord;

};



}
