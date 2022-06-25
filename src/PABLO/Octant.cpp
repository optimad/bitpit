/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "Octant.hpp"
#include "bitpit_common.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace bitpit {

/*!
 * Input stream operator for class Octant. Stream cell data from memory
 * input stream to container.
 *
 * \param[in] buffer is the input stream from memory
 * \param[in] octant is the octant object
 * \result Returns the same input stream received in input.
 */
IBinaryStream& operator>>(IBinaryStream &buffer, Octant &octant)
{
    uint8_t dimensions;
    buffer >> dimensions;

    uint8_t level;
    buffer >> level;

    octant.initialize(dimensions, level, true);

    buffer >> octant.m_morton;

    buffer >> octant.m_marker;

    buffer >> octant.m_ghost;

    for(int i = 0; i < Octant::INFO_ITEM_COUNT; ++i){
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
OBinaryStream& operator<<(OBinaryStream  &buffer, const Octant &octant)
{
    buffer << octant.m_dim;
    buffer << octant.m_level;

    buffer << octant.m_morton;

    buffer << octant.m_marker;

    buffer << octant.m_ghost;

    for(int i = 0; i < Octant::INFO_ITEM_COUNT; ++i){
        buffer << (bool) octant.m_info[i];
    }

    return buffer;
}

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

const TreeConstants::Instances Octant::sm_treeConstants = TreeConstants::instances();

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
 * It builds a 2D or 3D octant with user defined Morton index and level.
 * \param[in] dim Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level Refinement level of octant (0 for root octant).
 * \param[in] morton Morton number associated with the octant
 */
Octant::Octant(uint8_t dim, uint8_t level, uint64_t morton){
	initialize(dim, level, true);

	// Set the morton
	m_morton = morton;

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

	// Set the morton
	m_morton = PABLO::computeMorton(m_dim, x, y, z);

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

	// Set the morton
	m_morton = PABLO::computeMorton(m_dim, x, y, z);

};

/*! Check if two octants are equal (no check on info)
 */
bool Octant::operator ==(const Octant & oct2){
	bool check = true;
	check = check && (m_dim == oct2.m_dim);
	check = check && (m_morton == oct2.m_morton);
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

	// Reset the morton
	m_morton = 0;

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
Octant::getLogicalCoordinates() const{
	u32array3 coords;
	coords[0] = getLogicalX();
	coords[1] = getLogicalY();
	coords[2] = getLogicalZ();
	return coords;
};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate X of node 0.
 */
uint32_t
Octant::getLogicalX() const{
	return PABLO::computeCoordinate(m_dim, m_morton, 0);
};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Y of node 0.
 */
uint32_t
Octant::getLogicalY() const{
	return PABLO::computeCoordinate(m_dim, m_morton, 1);
};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Z of node 0.
 */
uint32_t
Octant::getLogicalZ() const{
	return PABLO::computeCoordinate(m_dim, m_morton, 2);
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
 * \return true if the face is a boundary.
 */
bool
Octant::getBound(uint8_t face) const{
	return (m_info[OctantInfo::INFO_BOUNDFACE0 + face]);
};

/*! Get the bound flag on an octant edge.
 * \param[in] edge local index of the edge.
 * \return true if the edge is on a boundary.
 */
bool
Octant::getEdgeBound(uint8_t edge) const{
	for (int face : sm_treeConstants[m_dim].edgeFace[edge]) {
		if (getBound(face)) {
			return true;
		}
	}

	return false;
};

/*! Get the bound flag on an octant node.
 * \param[in] node local index of the node.
 * \return true if the node is on a boundary.
 */
bool
Octant::getNodeBound(uint8_t node) const{
	for (int face : sm_treeConstants[m_dim].nodeFace[node]) {
		if (getBound(face)) {
			return true;
		}
	}

	return false;
};

/*! Get the bound flag on an octant.
 * \return true if the octant is a boundary octant.
 */
bool
Octant::getBound() const{
	for (int i = 0; i < sm_treeConstants[m_dim].nFaces; ++i) {
		if (getBound(i)) {
			return true;
		}
	}

	return false;
};

/*! Set the boundary flag to true on an octant face.
 * \param[in] face local index of the boundary face.
 */
void
Octant::setBound(uint8_t face) {
	m_info[INFO_BOUNDFACE0 + face] = true;
};

/*! Get the pbound flag on an octant face.
 * \param[in] face local index of the face.
 * \return true if the face-th face is a process boundary face.
 */
bool
Octant::getPbound(uint8_t face) const{
	return m_info[INFO_PBOUNDFACE0 + face];
};

/*! Get the pbound flag on an octant face.
 * \return true if the octant has at least one face that is a process boundary face.
 */
bool
Octant::getPbound() const{
	for (int i = 0; i < sm_treeConstants[m_dim].nFaces; ++i) {
		if (getPbound(i)) {
			return true;
		}
	}

	return false;
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
	m_info[INFO_PBOUNDFACE0 + face] = flag;
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
Octant::getLogicalSize() const{
	return sm_treeConstants[m_dim].lengths[m_level];
};

/*! Get the area of an octant in logical domain .
 * \return Area of octant.
 */
uint64_t
Octant::getLogicalArea() const{
	return sm_treeConstants[m_dim].areas[m_level];
};

/*! Get the volume of an octant in logical domain.
 * \return Volume of octant.
 */
uint64_t
Octant::getLogicalVolume() const{
	return sm_treeConstants[m_dim].volumes[m_level];
};

// =================================================================================== //

/*! Get the coordinates of the center of an octant in logical domain.
 * \return Array[3] with the coordinates of the center of octant.
 */
darray3
Octant::getLogicalCenter() const{
	double	dh;
	darray3 center;

	dh = double(getLogicalSize())*0.5;
	center[0] = getLogicalX() + dh;
	center[1] = getLogicalY() + dh;
	center[2] = getLogicalZ() + double(m_dim-2)*dh;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a face of an octant in logical domain.
 * \param[in] iface Local index of the face
 * \return Array[3] with the coordinates of the center of the octant face.
 */
darray3
Octant::getLogicalFaceCenter(uint8_t iface) const{
	double	dh_2;
	darray3 center;

	assert(iface < m_dim*2);

	dh_2 = double(getLogicalSize())*0.5;
	center[0] = getLogicalX() + (double)sm_treeConstants[m_dim].faceDisplacements[iface][0] * dh_2;
	center[1] = getLogicalY() + (double)sm_treeConstants[m_dim].faceDisplacements[iface][1] * dh_2;
	center[2] = getLogicalZ() + double(m_dim-2) * (double)sm_treeConstants[m_dim].faceDisplacements[iface][2] * dh_2;

	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a edge of an octant in logical domain.
 * \param[in] iedge Local index of the edge
 * \return Array[3] with the coordinates of the center of the octant edge.
 */
darray3
Octant::getLogicalEdgeCenter(uint8_t iedge) const{
	double	dh_2;
	darray3 center;

	dh_2 = double(getLogicalSize())*0.5;
	center[0] = getLogicalX() + (double)sm_treeConstants[m_dim].edgeDisplacements[iedge][0] * dh_2;
	center[1] = getLogicalY() + (double)sm_treeConstants[m_dim].edgeDisplacements[iedge][1] * dh_2;
	center[2] = getLogicalZ() + double(m_dim-2) * (double)sm_treeConstants[m_dim].edgeDisplacements[iedge][2] * dh_2;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \param[out] nodes Vector of arrays [nnodes][3] with the coordinates of the nodes of octant.
 */
void
Octant::getLogicalNodes(u32arr3vector & nodes) const{
	uint8_t		i;
	uint32_t	dh;
	uint8_t nn = uint8_t(1)<<m_dim;

	dh = getLogicalSize();
	nodes.resize(nn);

	u32array3 coords = getLogicalCoordinates();

	for (i = 0; i < nn; i++){
		nodes[i][0] = coords[0] + uint32_t(sm_treeConstants[m_dim].nodeCoordinates[i][0])*dh;
		nodes[i][1] = coords[1] + uint32_t(sm_treeConstants[m_dim].nodeCoordinates[i][1])*dh;
		nodes[i][2] = coords[2] + uint32_t(sm_treeConstants[m_dim].nodeCoordinates[i][2])*dh;
	}
};

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \return Vector of arrays [nnodes][3] with the coordinates of the nodes of octant.
 */
u32arr3vector
Octant::getLogicalNodes() const{
	u32arr3vector nodes;
	getLogicalNodes(nodes);
	return nodes;
};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[out] node Array[3] with the logical coordinates of the node of the octant.
 */
void		Octant::getLogicalNode(u32array3 & node, uint8_t inode) const{
	uint32_t	dh;

	dh = getLogicalSize();

	node = getLogicalCoordinates();

	node[0] += sm_treeConstants[m_dim].nodeCoordinates[inode][0]*dh;
	node[1] += sm_treeConstants[m_dim].nodeCoordinates[inode][1]*dh;
	node[2] += sm_treeConstants[m_dim].nodeCoordinates[inode][2]*dh;

};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \return Array[3] with the logical coordinates of the node of the octant.
 */
u32array3		Octant::getLogicalNode(uint8_t inode) const{
	u32array3 node;
	getLogicalNode(node, inode);
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


/** Get the Morton index of the octant (without level).
 * \return morton Morton index of the octant.
 */
uint64_t	Octant::getMorton() const{
	return m_morton;
};

/** Compute the Morton index of the last descendant octant of this octant.
 * \return Morton index of the last descendant octant.
 */
uint64_t	Octant::computeLastDescMorton() const {
	u32array3 lastDescCoordinates = computeLastDescCoordinates();
	return PABLO::computeMorton(m_dim, lastDescCoordinates[0], lastDescCoordinates[1], lastDescCoordinates[2]);
};

/** Compute the coordinates (i.e. the coordinates of the node 0) of the last
 * descendant octant of this octant.
 * \return The coordinates (i.e. the coordinates of the node 0) of the last
 * descendant octant of this octant.
 */
u32array3	Octant::computeLastDescCoordinates() const {
	u32array3 lastDescCoords = getLogicalCoordinates();
	for (int i=0; i<m_dim; i++){
		lastDescCoords[i] += (uint32_t(1) << (sm_treeConstants[m_dim].maxLevel - m_level)) - 1;
	}

	return lastDescCoords;
};

/** Compute the Morton index of the father of this octant.
 * \return Morton index of the father of this octant.
 */
uint64_t	Octant::computeFatherMorton() const {
	u32array3 fatherCoordinates = computeFatherCoordinates();
	return PABLO::computeMorton(m_dim, fatherCoordinates[0], fatherCoordinates[1], fatherCoordinates[2]);
};

/** Compute the coordinates (i.e. the coordinates of the node 0) of the father
 * of this octant.
 * \return The coordinates (i.e. the coordinates of the node 0) of the father
 * of this octant.
 */
u32array3	Octant::computeFatherCoordinates() const {
	u32array3 fatherCoordinates = getLogicalCoordinates();
	for (int i=0; i<m_dim; i++){
		fatherCoordinates[i] -= fatherCoordinates[i]%(uint32_t(1) << (sm_treeConstants[m_dim].maxLevel - max(0,(m_level-1))));
	}
	return fatherCoordinates;
};

/** Compute the persistent XYZ key of the given node (without level).
 * The XYZ key combines the 32bit coordinates and generates a unique 64bit value
 * (in two dimensions this is always possible, in three dimensions this is possible
 * because with a 64bit wide Morton number not all 32 bits of the coordinates are
 * used).
 * \param[in] inode Local index of the node
 * \return persistent XYZ key of the node.
 */
uint64_t	Octant::computeNodePersistentKey(uint8_t inode) const{

	u32array3 node = getLogicalNode(inode);

	return computeNodePersistentKey(node);
};

/** Compute the persistent XYZ key of the given node (without level).
 * The XYZ key combines the 32bit coordinates and generates a unique 64bit value
 * (in two dimensions this is always possible, in three dimensions this is possible
 * because with a 64bit wide Morton number not all 32 bits of the coordinates are
 * used).
 * \param[in] node Logical coordinates of the node
 * \return persistent XYZ key of the node.
 */
uint64_t	Octant::computeNodePersistentKey(const u32array3 &node) const{

	return PABLO::computeXYZKey(m_dim, node[0], node[1], node[2]);
};

/** Get the size of the buffer required to communicate the octant.
 * \return Returns the buffer size (in bytes).
 */
unsigned int Octant::getBinarySize()
{
    unsigned int binarySize = 0;
    binarySize += sizeof(uint8_t); // dimensions
    binarySize += sizeof(uint8_t); // level
    binarySize += sizeof(uint64_t); // morton
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
	uint64_t lastDescMorton = computeLastDescMorton();
	Octant lastDesc(m_dim, sm_treeConstants[m_dim].maxLevel, lastDescMorton);
	return lastDesc;
};

// =================================================================================== //

/** Build the father octant of this octant.
 * \return Father octant.
 */
Octant	Octant::buildFather() const {
	uint64_t fatherMorton = computeFatherMorton();
	Octant father(m_dim, max(0,m_level-1), fatherMorton);
	return father;
};

// =================================================================================== //

/** Count children of octant.
 *   \return The number of children of the octant.
 */
uint8_t	Octant::countChildren() const {
	if (this->m_level < sm_treeConstants[m_dim].maxLevel){
		return sm_treeConstants[m_dim].nChildren;
	} else {
		return 0;
	}
}

// =================================================================================== //

/** Builds children of octant.
 *   \return The children of the octant ordered by Z-index (info update)
 */
vector< Octant >	Octant::buildChildren() const {

	uint8_t nchildren = countChildren();
	std::vector< Octant > children(nchildren);
	buildChildren(children.data());

	return children;
}

// =================================================================================== //

/** Builds children of octant.
 *   \param[out] children On output will containt the children of the octant
 *   ordered by Z-index, it's up to the caller allocate enough space for all
 *   the possible children.
 */
void	Octant::buildChildren(Octant *children) const {

	u32array3 coords = getLogicalCoordinates();

	int nChildren = countChildren();
	for (int i=0; i<nChildren; ++i){
		// Octant information
		uint8_t xf;
		uint8_t yf;
		uint8_t zf;
		uint8_t dx;
		uint8_t dy;
		uint8_t dz;
		switch (i) {

		case 0 :
			dx = 0;
			dy = 0;
			dz = 0;

			xf = 1;
			yf = 3;
			zf = 5;

			break;

		case 1 :
			dx = 1;
			dy = 0;
			dz = 0;

			xf = 0;
			yf = 3;
			zf = 5;

			break;

		case 2 :
			dx = 0;
			dy = 1;
			dz = 0;

			xf = 1;
			yf = 2;
			zf = 5;

			break;

		case 3 :
			dx = 1;
			dy = 1;
			dz = 0;

			xf = 0;
			yf = 2;
			zf = 5;

			break;

		case 4 :
			dx = 0;
			dy = 0;
			dz = 1;

			xf = 1;
			yf = 3;
			zf = 4;

			break;

		case 5 :
			dx = 1;
			dy = 0;
			dz = 1;

			xf = 0;
			yf = 3;
			zf = 4;

			break;

		case 6 :
			dx = 0;
			dy = 1;
			dz = 1;

			xf = 1;
			yf = 2;
			zf = 4;

			break;

		case 7 :
			dx = 1;
			dy = 1;
			dz = 1;

			xf = 0;
			yf = 2;
			zf = 4;

			break;

		default:
			BITPIT_UNREACHABLE("The maximum number of children is 8.");

		}

		// Create octant
		children[i] = Octant(*this);
		Octant &oct = children[i];

		oct.setMarker(std::max(0, oct.m_marker - 1));
		oct.setLevel(oct.m_level + 1);

		uint32_t dh = oct.getLogicalSize();
		oct.m_morton = PABLO::computeMorton(m_dim, coords[0] + dh * dx, coords[1] + dh * dy, coords[2] + dh * dz);

		oct.m_info[OctantInfo::INFO_NEW4REFINEMENT] = true;

		oct.m_info[INFO_BOUNDFACE0 + xf] = false;
		oct.m_info[INFO_BOUNDFACE0 + yf] = false;
		oct.m_info[INFO_BOUNDFACE0 + zf] = false;

		oct.m_info[INFO_PBOUNDFACE0 + xf] = false;
		oct.m_info[INFO_PBOUNDFACE0 + yf] = false;
		oct.m_info[INFO_PBOUNDFACE0 + zf] = false;
	}
};

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant
 * throught face iface (sizehf=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeHalfSizeMortons(uint8_t iface, uint32_t *nMortons, std::vector<uint64_t> *mortons) const {

	if (m_info[iface]) {
		*nMortons = 0;
		mortons->clear();
		return;
	}

	int nchildren = 1<<m_dim;
	*nMortons = (m_level < sm_treeConstants[m_dim].maxLevel) ? nchildren/2 : 1;
	mortons->resize(*nMortons);

	u32array3 coords = getLogicalCoordinates();

	uint32_t dh  = (m_level < sm_treeConstants[m_dim].maxLevel) ? getLogicalSize()/2 : getLogicalSize();
	uint32_t dh2 = getLogicalSize();
	switch (iface) {
	case 0 :
	{
		uint32_t x = coords[0] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + ((i == 1) || (i == 3)) * dh;
			uint32_t z = coords[2] + ((m_dim  ==  3) && ((i == 2) || (i == 3))) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 1 :
	{
		uint32_t x = coords[0] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + ((i == 1) || (i == 3)) * dh;
			uint32_t z = coords[2] + ((m_dim  ==  3) && ((i == 2) || (i == 3))) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 2 :
	{
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((i == 1) || (i == 3)) * dh;
			uint32_t z = coords[2] + ((m_dim  ==  3) && ((i == 2) || (i == 3))) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 3 :
	{
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((i == 1) || (i == 3)) * dh;
			uint32_t z = coords[2] + ((m_dim  ==  3) && ((i == 2) || (i == 3))) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 4 :
	{
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((i == 1) || (i == 3)) * dh;
			uint32_t y = coords[1] + ((i == 2) || (i == 3)) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 5 :
	{
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((i == 1) || (i == 3)) * dh;
			uint32_t y = coords[1] + ((i == 2) || (i == 3)) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	}
}

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant
 * throught face iface (sizem=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeMinSizeMortons(uint8_t iface, uint8_t maxdepth, uint32_t *nMortons, std::vector<uint64_t> *mortons) const {

	if (m_info[iface]) {
		*nMortons = 0;
		mortons->clear();
		return;
	}

	*nMortons = (m_level < sm_treeConstants[m_dim].maxLevel) ? uint32_t(1)<<((m_dim-1)*(maxdepth-m_level)) : 1;
	mortons->resize(*nMortons);

	u32array3 coords = getLogicalCoordinates();

	uint32_t dh    = (m_level < sm_treeConstants[m_dim].maxLevel) ? uint32_t(1)<<(sm_treeConstants[m_dim].maxLevel - maxdepth) : getLogicalSize();
	uint32_t dh2   = getLogicalSize();
	uint32_t nline = uint32_t(1)<<(maxdepth-m_level);
	switch (iface) {
	case 0 :
	{
		uint32_t x = coords[0] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + ((m_dim == 2) * (i % nline) + (m_dim - 2) * (i / nline)) * dh;
			uint32_t z = coords[2] + (m_dim - 2) * (i % nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 1 :
	{
		uint32_t x = coords[0] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + ((m_dim == 2) * (i % nline) + (m_dim - 2) * (i / nline)) * dh;
			uint32_t z = coords[2] + (m_dim - 2) * (i % nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 2 :
	{
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((m_dim == 2) * (i % nline) + (m_dim - 2) * (i / nline)) * dh;
			uint32_t z = coords[2] + (m_dim - 2) * (i % nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 3 :
	{
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + ((m_dim == 2) * (i%nline) + (m_dim - 2) * (i / nline)) * dh;
			uint32_t z = coords[2] + (m_dim - 2) * (i%nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 4 :
	{
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i / nline) * dh;
			uint32_t y = coords[1] + (i % nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 5 :
	{
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i / nline) * dh;
			uint32_t y = coords[1] + (i % nline) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	}

	std::sort(mortons->begin(), mortons->end());
}

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught iface. Checks if balanced or not and uses half-size or min-size method
 * (sizeneigh=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeFaceVirtualMortons(uint8_t iface, uint8_t maxdepth, uint32_t *nMortons, std::vector<uint64_t> *mortons) const {
	if (!getBalance()){
		computeMinSizeMortons(iface, maxdepth, nMortons, mortons);
	} else{
		computeHalfSizeMortons(iface, nMortons, mortons);
	}
}

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * edge iedge (sizehf=0 if boundary octant)
 * \param[in] iedge Local index of the edge target.
 * \param[in] edgeface Local edge-face connectivity.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeEdgeHalfSizeMortons(uint8_t iedge, const uint8_t (&edgeface)[12][2], uint32_t *nMortons, std::vector<uint64_t> *mortons) const {

	uint32_t iface1 = edgeface[iedge][0];
	uint32_t iface2 = edgeface[iedge][1];
	if (m_info[iface1] || m_info[iface2]) {
		*nMortons = 0;
		mortons->clear();
		return;
	}

	*nMortons = (m_level < sm_treeConstants[m_dim].maxLevel) ? 2 : 1;
	mortons->resize(*nMortons);

	u32array3 coords = getLogicalCoordinates();

	uint32_t dh = (m_level < sm_treeConstants[m_dim].maxLevel) ? getLogicalSize()/2 : getLogicalSize();
	uint32_t dh2 = getLogicalSize();

	switch (iedge) {
	case 0 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 1 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 2 :
	{
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 3 :
	{
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 4 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 5 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 6 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 7 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 8 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 9 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 10 :
	{
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i == 1) * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 11 :
	{
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + (i == 1) *  dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	}
}

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * edge iedge (sizem=0 if boundary octant)
 * \param[in] iedge Local index of the edge target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[in] edgeface Local edge-face connectivity.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeEdgeMinSizeMortons(uint8_t iedge, uint8_t maxdepth, const uint8_t (&edgeface)[12][2], uint32_t *nMortons, std::vector<uint64_t> *mortons) const {

	uint8_t iface1 = edgeface[iedge][0];
	uint8_t iface2 = edgeface[iedge][1];
	if (m_info[iface1] || m_info[iface2]) {
		*nMortons = 0;
		mortons->clear();
		return;
	}

	*nMortons = (m_level < sm_treeConstants[m_dim].maxLevel) ? uint32_t(1)<<(maxdepth-m_level) : 1;
	mortons->resize(*nMortons);

	u32array3 coords = getLogicalCoordinates();

	uint32_t dh = (m_level < sm_treeConstants[m_dim].maxLevel) ? uint32_t(1)<<(sm_treeConstants[m_dim].maxLevel - maxdepth) : getLogicalSize();
	uint32_t dh2 = getLogicalSize();
	switch (iedge) {
	case 0 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + i * dh;
			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 1 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 2 :
	{
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 3 :
	{
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 4 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 5 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] - dh;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 6 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 7 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t z = coords[2] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 8 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 9 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t y = coords[1] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 10 :
	{
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	case 11 :
	{
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] + dh2;
		for (uint32_t i = 0; i < *nMortons; ++i) {
			uint32_t x = coords[0] + i * dh;

			(*mortons)[i] = PABLO::computeMorton(m_dim, x, y, z);
		}
	}
	break;
	}

	std::sort(mortons->begin(), mortons->end());
}

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught iedge. Checks if balanced or not and uses half-size or min-size method
 * (sizeneigh=0 if boundary octant).
 * \param[in] iedge Local index of the edge target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[in] balance_codim Maximum codimension of the entity for 2:1 balancing.
 * \param[in] edgeface Local edge-face connectivity.
 * \param[out] nMortons number of morton numbers.
 * \param[out] mortons are the requestd morton numbers.
 */
void Octant::computeEdgeVirtualMortons(uint8_t iedge, uint8_t maxdepth, const uint8_t balance_codim, const uint8_t (&edgeface)[12][2], uint32_t *nMortons, std::vector<uint64_t> *mortons) const {

	if (getBalance() && balance_codim > 1) {
		computeEdgeHalfSizeMortons(iedge, edgeface, nMortons, mortons);
	} else{
		computeEdgeMinSizeMortons(iedge, maxdepth, edgeface, nMortons, mortons);
	}
}

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * node inode (sizem=0 if boundary octant)
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[in] nodeface Local node-face connectivity.
 * \param[out] hasMorton true if the request morton exists
 * \param[out] morton the requested morton number.
 */
void Octant::computeNodeMinSizeMorton(uint8_t inode, uint8_t maxdepth, const uint8_t (&nodeface)[8][3], bool *hasMorton, uint64_t *morton) const {

	for (int i=0; i<m_dim; i++) {
		uint8_t iface = nodeface[inode][i];
		if (m_info[iface]) {
			*hasMorton = false;
			*morton = this->getMorton();
			return;
		}
	}

	*hasMorton = true;

	u32array3 coords = getLogicalCoordinates();

	uint32_t dh  = (m_level < sm_treeConstants[m_dim].maxLevel) ? uint32_t(1)<<(sm_treeConstants[m_dim].maxLevel - maxdepth) : getLogicalSize();
	uint32_t dh2 = getLogicalSize();
	switch (inode) {
	case 0 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] - (m_dim - 2) * dh;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 1 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] - (m_dim - 2) * dh;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 2 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] - (m_dim - 2) * dh;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 3 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] - (m_dim - 2) * dh;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 4 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] + dh2;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 5 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] - dh;
		uint32_t z = coords[2] + dh2;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 6 :
	{
		uint32_t x = coords[0] - dh;
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] + dh2;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	case 7 :
	{
		uint32_t x = coords[0] + dh2;
		uint32_t y = coords[1] + dh2;
		uint32_t z = coords[2] + dh2;

		*morton = PABLO::computeMorton(m_dim, x, y, z);
	}
	break;
	default:
		BITPIT_UNREACHABLE("The maximum number of nodes is 8.");
	}
}

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught inode. Uses min-size method (sizeneigh=0 if boundary octant).
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[in] nodeface Local node-face connectivity.
 * \param[out] hasMorton true if the request morton exists
 * \param[out] morton the requested morton number.
 */
void Octant::computeNodeVirtualMorton(uint8_t inode, uint8_t maxdepth, const uint8_t (&nodeface)[8][3], bool *hasMorton, uint64_t *morton) const {

	computeNodeMinSizeMorton(inode, maxdepth, nodeface, hasMorton, morton);

}

/*! Computes Morton index (without level) of same size
 * periodic neighbour of octant throught face iface.
 * \param[in] iface Local index of the face target.
 * \return Periodic neighbour morton number.
 */
uint64_t Octant::computePeriodicMorton(uint8_t iface) const {
	uint64_t Morton;
	uint32_t dh;
	dh = getLogicalSize();
	uint32_t maxLength = uint32_t(1)<<sm_treeConstants[m_dim].maxLevel;

	u32array3 coords = getLogicalCoordinates();

	if (!m_info[iface]){
		return this->getMorton();
	}
	else{
		switch (iface) {
		case 0 :
		{
			Morton = PABLO::computeMorton(m_dim,maxLength-dh,coords[1],coords[2]);
		}
		break;
		case 1 :
		{
			Morton = PABLO::computeMorton(m_dim,0,coords[1],coords[2]);
		}
		break;
		case 2 :
		{
			Morton = PABLO::computeMorton(m_dim,coords[0],maxLength-dh,coords[2]);
		}
		break;
		case 3 :
		{
			Morton = PABLO::computeMorton(m_dim,coords[0],0,coords[2]);
		}
		break;
		case 4 :
		{
			Morton = PABLO::computeMorton(m_dim,coords[0],coords[1],maxLength-dh);
		}
		break;
		case 5 :
		{
			Morton = PABLO::computeMorton(m_dim,coords[0],coords[1],0);
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
	uint32_t maxLength = uint32_t(1)<<sm_treeConstants[m_dim].maxLevel;
	uint32_t dh = this->getLogicalSize();

	u32array3 periodicCoords = getLogicalCoordinates();
	if (!m_info[iface]){
		return *this;
	}
	else{
		switch (iface) {
		case 0 :
		{
			periodicCoords[0] = maxLength-dh;
		}
		break;
		case 1 :
		{
			periodicCoords[0] = 0;
		}
		break;
		case 2 :
		{
			periodicCoords[1] = maxLength-dh;
		}
		break;
		case 3 :
		{
			periodicCoords[1] = 0;
		}
		break;
		case 4 :
		{
			periodicCoords[2] = maxLength-dh;
		}
		break;
		case 5 :
		{
			periodicCoords[2] = 0;
		}
		break;
		}

		Octant degOct(m_dim, m_level, periodicCoords[0], periodicCoords[1], periodicCoords[2]);
		degOct.m_info = false;
		return degOct;
	}

};

/** Build a same size periodic octant of this octant throught node inode.
 * \return Periodic octant of the same size (note: it is a stand-alone octant,
 * may be not living in octree).
 * \param[in] inode Local index of the node target.
 */
Octant Octant::computeNodePeriodicOctant(uint8_t inode) const {
    uint32_t maxLength = sm_treeConstants[m_dim].lengths[0];
    uint32_t dh = this->getLogicalSize();

    uint8_t iface1 = sm_treeConstants[m_dim].nodeFace[inode][0];
    uint8_t iface2 = sm_treeConstants[m_dim].nodeFace[inode][1];
    uint8_t iface3 = sm_treeConstants[m_dim].nodeFace[inode][m_dim-1];

    u32array3 periodicCoords = getLogicalCoordinates();
    if (!m_info[iface1] && !m_info[iface2] && !m_info[iface3]){
        return *this;
    }
    else{
        switch (iface1) {
        case 0 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE0]){
                periodicCoords[0] = maxLength-dh;
            }
            else{
                periodicCoords[0] -= dh;
            }
        }
        break;
        case 1 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE1]){
                periodicCoords[0] = 0;
            }
            else{
                periodicCoords[0] += dh;
            }
        }
        break;
        }

        switch (iface2) {
        case 2 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE2]){
                periodicCoords[1] = maxLength-dh;
            }
            else{
                periodicCoords[1] -= dh;
            }
        }
        break;
        case 3 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE3]){
                periodicCoords[1] = 0;
            }
            else{
                periodicCoords[1] += dh;
            }
        }
        break;
        }

        switch (iface3) {
        case 4 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE4]){
                periodicCoords[2] = maxLength-dh;
            }
            else{
                periodicCoords[2] -= dh;
            }
        }
        break;
        case 5 :
        {
            if (m_info[OctantInfo::INFO_BOUNDFACE5]){
                periodicCoords[2] = 0;
            }
            else{
                periodicCoords[2] += dh;
            }
        }
        break;
        }

        Octant degOct(m_dim, m_level, periodicCoords[0], periodicCoords[1], periodicCoords[2]);
        degOct.m_info = false;
        return degOct;
    }

};

/** Build a same size periodic octant of this octant throught edge inode.
 * \return Periodic octant of the same size (note: it is a stand-alone octant,
 * may be not living in octree).
 * \param[in] iedge Local index of the edge target.
 */
Octant Octant::computeEdgePeriodicOctant(uint8_t iedge) const {
    uint32_t maxLength = uint32_t(1)<<sm_treeConstants[m_dim].maxLevel;
    uint32_t dh = this->getLogicalSize();

    std::array<uint8_t,2> iface;
    iface[0] = sm_treeConstants[m_dim].edgeFace[iedge][0];
    iface[1] = sm_treeConstants[m_dim].edgeFace[iedge][1];

    int8_t          cxyz[3] = {0,0,0};
    for (int idim=0; idim<m_dim; idim++){
        cxyz[idim] = sm_treeConstants[m_dim].edgeCoeffs[iedge][idim];
    }

    u32array3 periodicCoords = getLogicalCoordinates();
    if (!m_info[iface[0]] && !m_info[iface[1]]){
        return *this;
    }
    else{
        for (std::size_t i=0; i<2; i++){
            switch (iface[i]) {
            case 0 :
                if (cxyz[0] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE0]){
                        periodicCoords[0] = maxLength-dh;
                    }
                    else{
                        periodicCoords[0] -= dh;
                    }
                }
                break;
            case 1 :
            {
                if (cxyz[0] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE1]){
                        periodicCoords[0] = 0;
                    }
                    else{
                        periodicCoords[0] += dh;
                    }
                }
            }
            break;
            case 2 :
            {
                if (cxyz[1] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE2]){
                        periodicCoords[1] = maxLength-dh;
                    }
                    else{
                        periodicCoords[1] -= dh;
                    }
                }
            }
            break;
            case 3 :
            {
                if (cxyz[1] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE3]){
                        periodicCoords[1] = 0;
                    }
                    else{
                        periodicCoords[1] += dh;
                    }
                }
            }
            break;
            case 4 :
            {
                if (cxyz[2] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE4]){
                        periodicCoords[2] = maxLength-dh;
                    }
                    else{
                        periodicCoords[2] -= dh;
                    }
                }
            }
            break;
            case 5 :
            {
                if (cxyz[2] != 0){
                    if (m_info[OctantInfo::INFO_BOUNDFACE5]){
                        periodicCoords[2] = 0;
                    }
                    else{
                        periodicCoords[2] += dh;
                    }
                }
            }
            break;
            }
        } // end loop on i

        Octant degOct(m_dim, m_level, periodicCoords[0], periodicCoords[1], periodicCoords[2]);
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
	coord[0] = getLogicalX();
	coord[1] = getLogicalY();
	coord[2] = getLogicalZ();
	int64_t dh = this->getLogicalSize();
	int64_t maxLength = int64_t(1)<<sm_treeConstants[m_dim].maxLevel;

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

/*! Get the coordinates of the octant shifted through node inode and
 * near the opposite periodic boundary (i.e. the coordinates considering this octant
 * as a ghost periodic octant).
 * \param[in] inode Local index of the node target.
 * \return Coordinates of octant considered as periodic ghost out of the logical domain.
 */
array<int64_t,3> Octant::getNodePeriodicCoord(uint8_t inode) const {
    array<int64_t,3> coord;
    coord[0] = getLogicalX();
    coord[1] = getLogicalY();
    coord[2] = getLogicalZ();
    int64_t dh = this->getLogicalSize();
    int64_t maxLength = int64_t(1)<<sm_treeConstants[m_dim].maxLevel;

    uint8_t iface1 = sm_treeConstants[m_dim].nodeFace[inode][0];
    uint8_t iface2 = sm_treeConstants[m_dim].nodeFace[inode][1];
    uint8_t iface3 = sm_treeConstants[m_dim].nodeFace[inode][m_dim-1];

    switch (iface1) {
    case 0 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE0]){
            coord[0] = maxLength;
        }
    }
    break;
    case 1 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE1]){
            coord[0]  = -dh;
        }
    }
    break;
    }

    switch (iface2) {
    case 2 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE2]){
            coord[1]  = maxLength;
        }
    }
    break;
    case 3 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE3]){
            coord[1] = -dh;
        }
    }
    break;
    }

    switch (iface3) {
    case 4 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE4]){
            coord[2] = maxLength;
        }
    }
    break;
    case 5 :
    {
        if (m_info[OctantInfo::INFO_BOUNDFACE5]){
            coord[2] = -dh;
        }
    }
    break;
    }

    return coord;

};

/*! Get the coordinates of the octant shifted through edge iedge and
 * near the opposite periodic boundary (i.e. the coordinates considering this octant
 * as a ghost periodic octant).
 * \param[in] iedge Local index of the edge target.
 * \return Coordinates of octant considered as periodic ghost out of the logical domain.
 */
array<int64_t,3> Octant::getEdgePeriodicCoord(uint8_t iedge) const {
    array<int64_t,3> coord;
    coord[0] = getLogicalX();
    coord[1] = getLogicalY();
    coord[2] = getLogicalZ();
    int64_t dh = this->getLogicalSize();
    int64_t maxLength = int64_t(1)<<sm_treeConstants[m_dim].maxLevel;

    std::array<uint8_t,2> iface;
    iface[0] = sm_treeConstants[m_dim].edgeFace[iedge][0];
    iface[1] = sm_treeConstants[m_dim].edgeFace[iedge][1];

    int8_t          cxyz[3] = {0,0,0};
    for (int idim=0; idim<m_dim; idim++){
        cxyz[idim] = sm_treeConstants[m_dim].edgeCoeffs[iedge][idim];
    }

    for (std::size_t i=0; i<2; i++){
        switch (iface[i]) {
        case 0 :
        {
            if (cxyz[0] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE0]){
                    coord[0] = maxLength;
                }
            }
        }
        break;
        case 1 :
        {
            if (cxyz[0] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE1]){
                    coord[0]  = -dh;
                }
            }
        }
        break;
        case 2 :
        {
            if (cxyz[1] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE2]){
                    coord[1]  = maxLength;
                }
            }
        }
        break;
        case 3 :
        {
            if (cxyz[1] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE3]){
                    coord[1] = -dh;
                }
            }
        }
        break;
        case 4 :
        {
            if (cxyz[2] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE4]){
                    coord[2] = maxLength;
                }
            }
        }
        break;
        case 5 :
        {
            if (cxyz[2] != 0){
                if (m_info[OctantInfo::INFO_BOUNDFACE5]){
                    coord[2] = -dh;
                }
            }
        }
        break;
        }
    } // end loop on i
    return coord;

};

/** Get the local index of the node corresponding to the splitting node of the octant family; i.e. the index of the local node
 * coincident with the center point of the father.
 * \return Local index of octant node corresponding to the splitting family node.
 */
uint8_t Octant::getFamilySplittingNode() const {

	u32array3 xx = getLogicalCoordinates();

	//Delta to build father (use the boolean negation to identify the splitting node)
	bool delta[3];
	delta[2] = 0;
	for (int i=0; i<m_dim; i++){
		delta[i] = ((xx[i]%(uint32_t(1) << (sm_treeConstants[m_dim].maxLevel - max(0,(m_level-1))))) == 0);
	}
	return sm_treeConstants[m_dim].nodeFromCoordinates[delta[0]][delta[1]][delta[2]];
};



}
