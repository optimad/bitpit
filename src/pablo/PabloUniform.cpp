// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "PabloUniform.hpp"

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
#if ENABLE_MPI==1
PabloUniform::PabloUniform(uint8_t dim, int8_t maxlevel, std::string logfile, MPI_Comm comm):ParaTree(dim,maxlevel,logfile,comm){
	m_origin = { {0,0,0} };
	m_L = 1.0;
};

PabloUniform::PabloUniform(double X, double Y, double Z, double L, uint8_t dim, int8_t maxlevel, std::string logfile, MPI_Comm comm):ParaTree(dim,maxlevel,logfile,comm){
	m_origin[0] = X;
	m_origin[1] = Y;
	m_origin[2] = Z;
	m_L = L;
};
#else
PabloUniform::PabloUniform(uint8_t dim, int8_t maxlevel, std::string logfile):ParaTree(dim,maxlevel,logfile){
	m_origin = { {0,0,0} };
	m_L = 1.0;
};

PabloUniform::PabloUniform(double X, double Y, double Z, double L, uint8_t dim, int8_t maxlevel, std::string logfile):ParaTree(dim,maxlevel,logfile){
	m_origin[0] = X;
	m_origin[1] = Y;
	m_origin[2] = Z;
	m_L = L;
};
#endif

// =================================================================================== //
// METHODS
// =================================================================================== //

// =================================================================================== //
// BASIC GET/SET METHODS															   //
// =================================================================================== //
darray3
PabloUniform::getOrigin(){
	return m_origin;
};

double
PabloUniform::getX0(){
	return m_origin[0];
};

double
PabloUniform::getY0(){
	return m_origin[1];
};

double
PabloUniform::getZ0(){
	return m_origin[2];
};

double
PabloUniform::getL(){
	return m_L;
};

/*! Set the length of the domain.
 * \param[in] Length of the octree.
 */
void
PabloUniform::setL(double L){
	m_L = L;
};

/*! Set the origin of the domain.
 * \param[in] Oriin of the octree.
 */
void
PabloUniform::setOrigin(darray3 origin){
	m_origin = origin;
};

/*! Get the size of an octant corresponding to a target level.
 * \param[in] idx Input level.
 * \return Size of an octant of input level.
 */
double
PabloUniform::levelToSize(uint8_t & level) {
	double size = ParaTree::levelToSize(level);
	return m_L *size;
}

// =================================================================================== //
// INDEX BASED METHODS																   //
// =================================================================================== //
darray3
PabloUniform::getCoordinates(uint32_t idx){
	darray3 coords, coords_;
	coords_ = ParaTree::getCoordinates(idx);
	for (int i=0; i<3; i++){
		coords[i] = m_origin[i] + m_L * coords_[i];
	}
	return coords;
};

double
PabloUniform::getX(uint32_t idx){
	double X, X_;
	X_ = ParaTree::getX(idx);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getY(uint32_t idx){
	double X, X_;
	X_ = ParaTree::getY(idx);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getZ(uint32_t idx){
	double X, X_;
	X_ = ParaTree::getZ(idx);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getSize(uint32_t idx){
	return m_L * ParaTree::getSize(idx);
};

double
PabloUniform::getArea(uint32_t idx){
	return m_L * m_L * ParaTree::getArea(idx);
};

double
PabloUniform::getVolume(uint32_t idx){
	return m_L * m_L * m_L* ParaTree::getArea(idx);
};

void
PabloUniform::getCenter(uint32_t idx, darray3& center){
	darray3 center_ = ParaTree::getCenter(idx);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
};

darray3
PabloUniform::getCenter(uint32_t idx){
	darray3 center, center_ = ParaTree::getCenter(idx);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
	return center;
};

void
PabloUniform::getFaceCenter(uint32_t idx, uint8_t iface, darray3& center){
	darray3 center_ = ParaTree::getFaceCenter(idx, iface);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
};

darray3
PabloUniform::getFaceCenter(uint32_t idx, uint8_t iface){
	darray3 center, center_ = ParaTree::getFaceCenter(idx, iface);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
	return center;
};

darray3
PabloUniform::getNode(uint32_t idx, uint8_t inode){
	darray3 node, node_ = ParaTree::getNode(idx, inode);
	for (int i=0; i<3; i++){
		node[i] = m_origin[i] + m_L * node_[i];
	}
	return node;
};

void
PabloUniform::getNode(uint32_t idx, uint8_t inode, darray3& node){
	darray3 node_ = ParaTree::getNode(idx, inode);
	for (int i=0; i<3; i++){
		node[i] = m_origin[i] + m_L * node_[i];
	}
};

void
PabloUniform::getNodes(uint32_t idx, darr3vector & nodes){
	darray3vector nodes_ = ParaTree::getNodes(idx);
	nodes.resize(ParaTree::getNnodes());
	for (int j=0; j<ParaTree::getNnodes(); j++){
		for (int i=0; i<3; i++){
			nodes[j][i] = m_origin[i] + m_L * nodes_[j][i];
		}
	}
};

darr3vector
PabloUniform::getNodes(uint32_t idx){
	darray3vector nodes, nodes_ = ParaTree::getNodes(idx);
	nodes.resize(ParaTree::getNnodes());
	for (int j=0; j<ParaTree::getNnodes(); j++){
		for (int i=0; i<3; i++){
			nodes[j][i] = m_origin[i] + m_L * nodes_[j][i];
		}
	}
	return nodes;
};

void
PabloUniform::getNormal(uint32_t idx, uint8_t & iface, darray3 & normal) {
	ParaTree::getNormal(idx, iface, normal);
}

darray3
PabloUniform::getNormal(uint32_t idx, uint8_t & iface){
	return ParaTree::getNormal(idx, iface);
}

// =================================================================================== //
// POINTER BASED METHODS															   //
// =================================================================================== //
darray3
PabloUniform::getCoordinates(Octant* oct){
	darray3 coords, coords_;
	coords_ = ParaTree::getCoordinates(oct);
	for (int i=0; i<3; i++){
		coords[i] = m_origin[i] + m_L * coords_[i];
	}
	return coords;
};

double
PabloUniform::getX(Octant* oct){
	double X, X_;
	X_ = ParaTree::getX(oct);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getY(Octant* oct){
	double X, X_;
	X_ = ParaTree::getY(oct);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getZ(Octant* oct){
	double X, X_;
	X_ = ParaTree::getZ(oct);
	X = m_origin[0] + m_L * X_;
	return X;
};

double
PabloUniform::getSize(Octant* oct){
	return m_L * ParaTree::getSize(oct);
};

double
PabloUniform::getArea(Octant* oct){
	return m_L * m_L * ParaTree::getArea(oct);
};

double
PabloUniform::getVolume(Octant* oct){
	return m_L * m_L * m_L* ParaTree::getArea(oct);
};

void
PabloUniform::getCenter(Octant* oct, darray3& center){
	darray3 center_ = ParaTree::getCenter(oct);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
};

darray3
PabloUniform::getCenter(Octant* oct){
	darray3 center, center_ = ParaTree::getCenter(oct);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
	return center;
};

void
PabloUniform::getFaceCenter(Octant* oct, uint8_t iface, darray3& center){
	darray3 center_ = ParaTree::getFaceCenter(oct, iface);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
};

darray3
PabloUniform::getFaceCenter(Octant* oct, uint8_t iface){
	darray3 center, center_ = ParaTree::getFaceCenter(oct, iface);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center_[i];
	}
	return center;
};

darray3
PabloUniform::getNode(Octant* oct, uint8_t inode){
	darray3 node, node_ = ParaTree::getNode(oct, inode);
	for (int i=0; i<3; i++){
		node[i] = m_origin[i] + m_L * node_[i];
	}
	return node;
};

void
PabloUniform::getNode(Octant* oct, uint8_t inode, darray3& node){
	darray3 node_ = ParaTree::getNode(oct, inode);
	for (int i=0; i<3; i++){
		node[i] = m_origin[i] + m_L * node_[i];
	}
};

void
PabloUniform::getNodes(Octant* oct, darr3vector & nodes){
	darray3vector nodes_ = ParaTree::getNodes(oct);
	nodes.resize(ParaTree::getNnodes());
	for (int j=0; j<ParaTree::getNnodes(); j++){
		for (int i=0; i<3; i++){
			nodes[j][i] = m_origin[i] + m_L * nodes_[j][i];
		}
	}
};

darr3vector
PabloUniform::getNodes(Octant* oct){
	darray3vector nodes, nodes_ = ParaTree::getNodes(oct);
	nodes.resize(ParaTree::getNnodes());
	for (int j=0; j<ParaTree::getNnodes(); j++){
		for (int i=0; i<3; i++){
			nodes[j][i] = m_origin[i] + m_L * nodes_[j][i];
		}
	}
	return nodes;
};

void
PabloUniform::getNormal(Octant* oct, uint8_t & iface, darray3 & normal) {
	ParaTree::getNormal(oct, iface, normal);
}

darray3
PabloUniform::getNormal(Octant* oct, uint8_t & iface){
	return ParaTree::getNormal(oct, iface);
}

// =================================================================================== //
// LOCAL TREE GET/SET METHODS														   //
// =================================================================================== //
double
PabloUniform::getLocalMaxSize(){
	return m_L * ParaTree::getLocalMaxSize();
};

double
PabloUniform::getLocalMinSize(){
	return m_L * ParaTree::getLocalMinSize();
};


/*! Get the coordinates of the extreme points of a bounding box containing the local tree
 *  \param[out] P0 Array with coordinates of the first point (lowest coordinates);
 *  \param[out] P1 Array with coordinates of the last point (highest coordinates).
 */
void
PabloUniform::getBoundingBox(darray3 & P0, darray3 & P1){
	darray3		cnode, cnode0, cnode1;
	uint32_t 	nocts = ParaTree::getNumOctants();
	uint32_t	id = 0;
	uint8_t 	nnodes = ParaTree::getNnodes();

	P0 = ParaTree::getNode(id, 0);
	P1 = ParaTree::getNode(nocts-1, nnodes-1);

	for (id=0; id<nocts; id++){
		cnode0 = ParaTree::getNode(id, 0);
		cnode1 = ParaTree::getNode(id, nnodes-1);
		for (int i=0; i<3; i++){
			P0[i] = min(P0[i], (double)cnode0[i]);
			P1[i] = max(P1[i], (double)cnode1[i]);
		}
	}
	for (int i=0; i<3; i++){
		P0[i] = m_origin[i] + m_L * P0[i];
		P1[i] = m_origin[i] + m_L * P1[i];
	}
};


// =================================================================================== //
// INTERSECTION GET/SET METHODS														   //
// =================================================================================== //
double
PabloUniform::getSize(Intersection* inter){
	return m_L * ParaTree::getSize(inter);
};

double
PabloUniform::getArea(Intersection* inter){
	return m_L * m_L * ParaTree::getArea(inter);
};

darray3
PabloUniform::getCenter(Intersection* inter){
	darray3 center = ParaTree::getCenter(inter);
	for (int i=0; i<3; i++){
		center[i] = m_origin[i] + m_L * center[i];
	}
	return center;
}

// =================================================================================== //
// OTHER OCTANT BASED METHODS												    	   //
// =================================================================================== //
Octant* PabloUniform::getPointOwner(darray3 & point){
	for (int i=0; i<3; i++){
		point[i] = (point[i] - m_origin[i])/m_L;
	}
	return ParaTree::getPointOwner(point);
};

uint32_t PabloUniform::getPointOwnerIdx(darray3 & point){
	for (int i=0; i<3; i++){
		point[i] = (point[i] - m_origin[i])/m_L;
	}
	return ParaTree::getPointOwnerIdx(point);
};

