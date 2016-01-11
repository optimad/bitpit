// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classMap.hpp"

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*!Default constructor. Origin of octree in physical domain in (0,0,0)
 * and side length 1.
 */
classMap::classMap(int8_t maxlevel, uint8_t dim_){
	maxlevel = int8_t(max(0,min(int(maxlevel),21)));
	X0 = Y0 = Z0 = 0.0;
	L = 1.0;
	dim = dim_;
	nnodes = 1<<dim;
	nnodesperface = 1<<(dim-1);
	max_length = uint32_t(1<<maxlevel);

};

/*!Customized constructor with origin of octree in physical
 * domain side length provided by the user.
 * \param[in] X Coordinate X of the origin.
 * \param[in] Y Coordinate Y of the origin.
 * \param[in] Z Coordinate Z of the origin.
 * \param[in] LL Side length of domain.
 */
classMap::classMap(double & X, double & Y, double & Z, double & LL, int8_t maxlevel, uint8_t dim_){
	maxlevel = int8_t(max(0,min(int(maxlevel),21)));
	X0 = X;
	Y0 = Y;
	Z0 = Z;
	L = LL;
	dim = dim_;
	nnodes = 1<<dim;
	nnodesperface = 1<<(dim-1);
	max_length = uint32_t(1<<maxlevel);
};

// =================================================================================== //
// METHODS
// =================================================================================== //

/*! Transformation of coordinate X.
 * \param[in] X Coordinate X from logical domain.
 * \return Coordinate X in physical domain.
 */
double classMap::mapX(uint32_t const & X){
	return (X0 + L/double(max_length) * double(X));
};

/*! Transformation of coordinate Y.
 * \param[in] Y Coordinate Y from logical domain.
 * \return Coordinate Y in physical domain.
 */
double classMap::mapY(uint32_t const & Y){
	return (Y0 + L/double(max_length) * double(Y));
};

/*! Transformation of coordinate Z.
 * \param[in] Z Coordinate Z from logical domain.
 * \return Coordinate Z in physical domain.
 */
double classMap::mapZ(uint32_t const & Z){
	return (Z0 + L/double(max_length) * double(Z));
};

/*! Transformation of coordinate X.
 * \param[in] X Coordinate X from physical domain.
 * \return Coordinate X in logical domain.
 */
uint32_t classMap::mapX(double const & X){
	return (uint32_t)(double(max_length)/L * (X - X0));
};

/*! Transformation of coordinate Y.
 * \param[in] Y Coordinate Y from physical domain.
 * \return Coordinate Y in logical domain.
 */
uint32_t classMap::mapY(double const & Y){
	return (uint32_t)(double(max_length)/L * (Y - Y0));
};

/*! Transformation of coordinate Z.
 * \param[in] Z Coordinate Z from physical domain.
 * \return Coordinate Z in logical domain.
 */
uint32_t classMap::mapZ(double const & Z){
	return (uint32_t)(double(max_length)/L * (Z - Z0));
};

/*! Transformation of size of an octant.
 * \param[in] size Size of octant from logical domain.
 * \return Size of octant in physical domain.
 */
double classMap::mapSize(uint32_t const & size){
	return ((L/double(max_length))*double(size));
};

/*! Transformation of area of an octant.
 * \param[in] area Area of octant from logical domain.
 * \return Area of octant in physical domain.
 */
double classMap::mapArea(uint64_t const & Area){
	return ((pow(L,dim-1)/pow(double(max_length),dim-1))*double(Area));
};

/*! Transformation of volume of an octant.
 * \param[in] volume Volume of octant from logical domain.
 * \return Coordinate Volume of octant in physical domain.
 */
double classMap::mapVolume(uint64_t const & Volume){
	return ((pow(L,dim)/pow(double(max_length),dim))*double(Volume));
};

/*! Transformation of coordinates of center of an octant.
 * \param[in] center Pointer to coordinates of center from logical domain.
 * \param[out] mapcenter Coordinates of center in physical domain.
 */
void classMap::mapCenter(double* & center,
		dvector & mapcenter){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapcenter.resize(3);
	mapcenter = orig;
	for (int i=0; i<dim; i++){
		mapcenter[i] = mapcenter[i] + L/double(max_length) * center[i];
	}
	dvector(mapcenter).swap(mapcenter);
};

/*! Transformation of coordinates of center of an octant.
 * \param[in] center Vector of coordinates of center from logical domain.
 * \param[out] mapcenter Coordinates of center in physical domain.
 */
void classMap::mapCenter(dvector & center,
		dvector & mapcenter){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapcenter.resize(3);
	mapcenter = orig;
	for (int i=0; i<dim; i++){
		mapcenter[i] = mapcenter[i] + L/double(max_length) * center[i];
	}
	dvector(mapcenter).swap(mapcenter);
};

/*! Transformation of coordinates of nodes of an octant.
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void classMap::mapNodes(uint32_t (*nodes)[3],
		dvector2D & mapnodes){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapnodes.resize(nnodes);
	for (int i=0; i<nnodes; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(max_length) * double(nodes[i][j]);
		}
		dvector(mapnodes[i]).swap(mapnodes[i]);
	}
	dvector2D(mapnodes).swap(mapnodes);
};

/*! Transformation of coordinates of nodes of an octant.
 * \param[in] nodes Vector of coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void classMap::mapNodes(u32vector2D nodes,
		dvector2D & mapnodes){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapnodes.resize(nnodes);
	for (int i=0; i<nnodes; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(max_length) * double(nodes[i][j]);
		}
		dvector(mapnodes[i]).swap(mapnodes[i]);
	}
	dvector2D(mapnodes).swap(mapnodes);
};

/*! Transformation of coordinates of a node of an octant.
 * \param[in] node Coordinates of  the node from logical domain.
 * \param[out] mapnodes Coordinates of the node in physical domain.
 */
void classMap::mapNode(u32vector & node,
		dvector & mapnode){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapnode.resize(3);
	for (int j=0; j<3; j++){
		mapnode[j] = orig[j] + L/double(max_length) * double(node[j]);
	}
	dvector(mapnode).swap(mapnode);
};

/*! Transformation of coordinates of nodes of an intersection.
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void classMap::mapNodesIntersection(uint32_t (*nodes)[3],
		dvector2D & mapnodes){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapnodes.resize(nnodesperface);
	for (int i=0; i<nnodesperface; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(max_length) * double(nodes[i][j]);
		}
		dvector(mapnodes[i]).swap(mapnodes[i]);
	}
	dvector2D(mapnodes).swap(mapnodes);
};

/*! Transformation of coordinates of nodes of an intersection.
 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
 * \param[out] mapnodes Coordinates of nodes in physical domain.
 */
void classMap::mapNodesIntersection(u32vector2D nodes,
		dvector2D & mapnodes){
	dvector orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	dvector(orig).swap(orig);
	mapnodes.resize(nnodesperface);
	for (int i=0; i<nnodesperface; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(max_length) * double(nodes[i][j]);
		}
		dvector(mapnodes[i]).swap(mapnodes[i]);
	}
	dvector2D(mapnodes).swap(mapnodes);
};

/*! Transformation of components of normal of an intersection.
 * \param[in] nodes Pointer to components of normal from logical domain.
 * \param[out] mapnodes components of normal in physical domain.
 */
void classMap::mapNormals(vector<int8_t> normal,
		dvector & mapnormal){
	mapnormal = dvector(normal.begin(), normal.end());
	dvector(mapnormal).swap(mapnormal);
};
