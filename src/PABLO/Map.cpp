// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Global.hpp"
#include "Operators.hpp"
#include "Map.hpp"
#include <cmath>

namespace bitpit {

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*!Default constructor.
 */
Map::Map(){
    initialize();
};

/*!Creates a new transformation of coordinates.
 * Origin of octree in reference domain in (0,0,0) and side length 1.
 * \param[in] dim The space dimension of the m_octree.
 */
Map::Map(uint8_t dim){
    initialize(dim);
};

// =================================================================================== //
// METHODS
// =================================================================================== //

/*! Initialize a dummy transformation of coordinates.
 */
void Map::initialize(){
    initialize(0);
}

/*! Initialize the transformation of coordinates.
 *  Origin of octree in reference domain in (0,0,0) and side length 1.
 * \param[in] dim The space dimension of the m_octree.
 */
void Map::initialize(uint8_t dim){
    m_dim = dim;

    m_origin[0] = 0.;
    m_origin[1] = 0.;
    m_origin[2] = 0.;

    m_L = 1.0;

    if (m_dim > 0) {
        m_nnodes        = uint8_t(1)<<m_dim;
        m_nnodesPerFace = uint8_t(1)<<(m_dim-1);

        m_maxLength = Global::getMaxLength();
        m_maxArea   = uipow<uint64_t>(m_maxLength, m_dim - 1);
        m_maxVolume = uipow<uint64_t>(m_maxLength, m_dim);

        m_maxLength_1 = 1. / double(m_maxLength);
        m_maxArea_1   = 1. / double(m_maxArea);
        m_maxVolume_1 = 1. / double(m_maxVolume);
    } else {
        m_nnodes        = 0;
        m_nnodesPerFace = 0;

        m_maxLength = 0;
        m_maxArea   = 0;
        m_maxVolume = 0;

        m_maxLength_1 = std::numeric_limits<double>::infinity();
        m_maxArea_1   = std::numeric_limits<double>::infinity();
        m_maxVolume_1 = std::numeric_limits<double>::infinity();
    }

};

/*! Transformation of coordinates X,Y,Z (logical->physical).
 * \param[in] X Coordinates from logical domain.
 * \return Coordinates in physical domain.
 */
darray3 Map::mapCoordinates(u32array3 const & X) const {
	darray3 coords;
	for (int i=0; i<3; ++i){
		coords[i] = m_maxLength_1 * double(X[i]);
	}
	return coords;
};

/*! Transformation of coordinate X (logical->physical).
 * \param[in] X Coordinate X from logical domain.
 * \return Coordinate X in physical domain.
 */
double Map::mapX(uint32_t const & X) const {
	return m_maxLength_1 * double(X);
};

/*! Transformation of coordinate Y (logical->physical).
 * \param[in] Y Coordinate Y from logical domain.
 * \return Coordinate Y in physical domain.
 */
double Map::mapY(uint32_t const & Y) const {
	return m_maxLength_1 * double(Y);
};

/*! Transformation of coordinate Z (logical->physical).
 * \param[in] Z Coordinate Z from logical domain.
 * \return Coordinate Z in physical domain.
 */
double Map::mapZ(uint32_t const & Z) const {
	return m_maxLength_1 * double(Z);
};

/*! Transformation of coordinates X,Y,Z (physical->logical).
 * \param[in] X Coordinates from physical domain.
 * \return Coordinates in logical domain.
 */
u32array3 Map::mapCoordinates(darray3 const & X) const {
	u32array3 coords;
	for (int i=0; i<3; ++i){
		coords[i] = (uint32_t)(double(m_maxLength) * X[i]);
	}
	return coords;
};

/*! Transformation of coordinate X (physical->logical).
 * \param[in] X Coordinate X from physical domain.
 * \return Coordinate X in logical domain.
 */
uint32_t Map::mapX(double const & X) const {
	return (uint32_t)(double(m_maxLength) * X);
};

/*! Transformation of coordinate Y (physical->logical).
 * \param[in] Y Coordinate Y from physical domain.
 * \return Coordinate Y in logical domain.
 */
uint32_t Map::mapY(double const & Y) const {
	return (uint32_t)(double(m_maxLength) * Y);
};

/*! Transformation of coordinate Z (physical->logical).
 * \param[in] Z Coordinate Z from physical domain.
 * \return Coordinate Z in logical domain.
 */
uint32_t Map::mapZ(double const & Z) const {
	return (uint32_t)(double(m_maxLength) * Z);
};

/*! Transformation of size of an octant (logical->physical).
 * \param[in] size Size of octant from logical domain.
 * \return Size of octant in physical domain.
 */
double Map::mapSize(uint32_t const & size) const {
	return m_maxLength_1 *double(size);
};

/*! Transformation of area of an octant (logical->physical).
 * \param[in] area Area of octant from logical domain.
 * \return Area of octant in physical domain.
 */
double Map::mapArea(uint64_t const & area) const {
	return m_maxArea_1*double(area);
};

/*! Transformation of volume of an octant (logical->physical).
 * \param[in] volume Volume of octant from logical domain.
 * \return Coordinate Volume of octant in physical domain.
 */
double Map::mapVolume(uint64_t const & volume) const {
	return m_maxVolume_1*double(volume);
};

/*! Transformation of coordinates of center of an octant (logical->physical).
 * \param[in] center Pointer to coordinates of center from logical domain.
 * \param[out] mapcenter Coordinates of center in physical domain.
 */
void Map::mapCenter(double* & center, darray3 & mapcenter) const {
	for (int i=0; i<3; i++){
		mapcenter[i] = m_maxLength_1 * center[i];
	}
};

/*! Transformation of coordinates of center of an octant (logical->physical).
 * \param[in] center Array of coordinates of center from logical domain.
 * \param[out] mapcenter Coordinates of center in physical domain.
 */
void Map::mapCenter(darray3 & center, darray3 & mapcenter) const {
	for (int i=0; i<3; i++){
		mapcenter[i] = m_maxLength_1 * center[i];
	}
};

/*! Transformation of coordinates of nodes of an octant (logical->physical).
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void Map::mapNodes(uint32_t (*nodes)[3], darr3vector & mapnodes) const {
	mapnodes.resize(m_nnodes);
	for (int i=0; i<m_nnodes; i++){
		for (int j=0; j<3; j++){
			mapnodes[i][j] = m_maxLength_1 * double(nodes[i][j]);
		}
	}
};

/*! Transformation of coordinates of nodes of an octant (logical->physical).
 * \param[in] nodes Vector of coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void Map::mapNodes(u32arr3vector nodes, darr3vector & mapnodes) const {
	mapnodes.resize(m_nnodes);
	for (int i=0; i<m_nnodes; i++){
		for (int j=0; j<3; j++){
			mapnodes[i][j] = m_maxLength_1 * double(nodes[i][j]);
		}
	}
};

/*! Transformation of coordinates of a node of an octant (logical->physical).
 * \param[in] node Coordinates of  the node from logical domain.
 * \param[out] mapnode Coordinates of the node in physical domain.
 */
void Map::mapNode(u32array3 & node, darray3 & mapnode) const {
	for (int j=0; j<3; j++){
		mapnode[j] = m_maxLength_1 * double(node[j]);
	}
};

/*! Transformation of coordinates of nodes of an intersection (logical->physical).
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void Map::mapNodesIntersection(uint32_t (*nodes)[3], darr3vector & mapnodes) const {
	mapnodes.resize(m_nnodesPerFace);
	for (int i=0; i<m_nnodesPerFace; i++){
		for (int j=0; j<3; j++){
			mapnodes[i][j] = m_maxLength_1 * double(nodes[i][j]);
		}
	}
};

/*! Transformation of coordinates of nodes of an intersection (logical->physical).
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void Map::mapNodesIntersection(u32arr3vector nodes, darr3vector & mapnodes) const {
	mapnodes.resize(m_nnodesPerFace);
	for (int i=0; i<m_nnodesPerFace; i++){
		for (int j=0; j<3; j++){
			mapnodes[i][j] = m_maxLength_1 * double(nodes[i][j]);
		}
	}
};

/*! Transformation of components of normal of an intersection (logical->physical).
 * \param[in] normal Pointer to components of normal from logical domain.
 * \param[out] mapnormal components of normal in physical domain.
 */
void Map::mapNormals(i8array3 normal, darray3 & mapnormal) const {
	mapnormal[0] = double(normal[0]);
	mapnormal[1] = double(normal[1]);
	mapnormal[2] = double(normal[2]);
};

}
