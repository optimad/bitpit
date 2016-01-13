// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ClassOctant.hpp"
#include <algorithm>

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

constexpr int ClassOctant::sm_CoeffFaceCenter[6][3];
constexpr int ClassOctant::sm_CoeffEdgeCenter[12][3];

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*! Default constructor of an octant.
 * It builds a 2D zero-level octant with origin in (0,0,0).
 */
ClassOctant::ClassOctant(){
	m_dim = 2;
	m_x = m_y = m_z = 0;
	m_level = 0;
	m_marker = 0;
	//default constructor of bitset is zero value -> set boundary condition true for faces
	for (uint8_t i=0; i<4; i++){
		m_info[i] = true;
	}
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D zero-level octant with origin in (0,0,0).
 * \param[in] dim_ Dimension of octant (2/3 for 2D/3D octant).
 */
ClassOctant::ClassOctant(uint8_t dim_){
	m_dim = dim_;
	m_x = m_y = m_z = 0;
	m_level = 0;
	m_marker = 0;
	uint8_t nf = m_dim*2;
	//default constructor of bitset is zero value -> set boundary condition true for faces
	for (uint8_t i=0; i<nf; i++){
		m_info[i] = true;
	}
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D octant with user defined origin and level.
 * \param[in] dim_ Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level_ Refinement level of octant (0 for root octant).
 * \param[in] x_,y_,z_ Coordinates of the origin of the octant (default values for z=0).
 */
ClassOctant::ClassOctant(uint8_t dim_, uint8_t level_, int32_t x_, int32_t y_, int32_t z_){
	m_dim = dim_;
	m_x = x_;
	m_y = y_;
	m_z = (m_dim-2)*z_;
	m_level = level_;
	m_marker = 0;
	//default constructor of bitset is zero value -> set boundary condition true for faces
	if (m_level==0){
		uint8_t nf = m_dim*2;
		for (uint8_t i=0; i<nf; i++){
			m_info[i] = true;
		}
	}
};

/*! Custom constructor of an octant.
 * It builds a 2D or 3D octant with user defined origin, level and boundary conditions.
 * \param[in] bound Boundary condition for the faces of the octant (the same for each face).
 * \param[in] dim_ Dimension of octant (2/3 for 2D/3D octant).
 * \param[in] level_ Refinement level of octant (0 for root octant).
 * \param[in] x_,y_,z_ Coordinates of the origin of the octant (default values for z=0).
 */
ClassOctant::ClassOctant(bool bound, uint8_t dim_, uint8_t level_, int32_t x_, int32_t y_, int32_t z_){
	m_dim = dim_;
	m_x = x_;
	m_y = y_;
	m_z = (m_dim-2)*z_;
	m_level = level_;
	m_marker = 0;
	//default constructor of bitset is zero value -> set boundary condition bound for faces
	if (m_level==0){
		uint8_t nf = m_dim*2;
		for (uint8_t i=0; i<nf; i++){
			m_info[i] = bound;
		}
	}
};

/*! Copy constructor of an octant.
 */
ClassOctant::ClassOctant(const ClassOctant &octant){
	m_dim = octant.m_dim;
	m_x = octant.m_x;
	m_y = octant.m_y;
	m_z = octant.m_z;
	m_level = octant.m_level;
	m_marker = octant.m_marker;
	m_info = octant.m_info;
};

/*! Check if two octants are equal (no check on info)
 */
bool ClassOctant::operator ==(const ClassOctant & oct2){
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

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //

/*! Get the dimension of an octant, 2 or 3 if 2D or 3D octant.
 * \return Dimension of octant.
 */
uint32_t
ClassOctant::getDim() const{return m_dim;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinates of node 0.
 */
u32array3
ClassOctant::getCoordinates() const{
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
ClassOctant::getX() const{return m_x;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Y of node 0.
 */
uint32_t
ClassOctant::getY() const{return m_y;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Z of node 0.
 */
uint32_t
ClassOctant::getZ() const{return m_z;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinates of node 0.
 */
u32array3
ClassOctant::getCoord(){
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
ClassOctant::getLevel() const{return m_level;};

/*! Get the refinement marker of an octant.
 * \return Marker of octant.
 */
int8_t
ClassOctant::getMarker() const{return m_marker;};

/*! Get the bound flag on an octant face.
 * \param[in] iface local index of the face.
 * \return true if the iface face is a boundary face.
 */
bool
ClassOctant::getBound(uint8_t face) const{
	return m_info[face];
};

/*! Get the bound flag on an octant.
 * \return true if the octant is a boundary octant.
 */
bool
ClassOctant::getBound() const{
	return m_info[0]||m_info[1]||m_info[2]||m_info[3]||((m_dim-2)*(m_info[4]||m_info[5]));
};

/*! Set the boundary flag to true on an octant face.
 * \param[in] face local index of the boundary face.
 */
void
ClassOctant::setBound(uint8_t face) {
	m_info[face] = true;
};

/*! Get the pbound flag on an octant face.
 * \param[in] face local index of the face.
 * \return true if the face-th face is a process boundary face.
 */
bool
ClassOctant::getPbound(uint8_t face) const{
	return m_info[6+face];
};

/*! Get the pbound flag on an octant face.
 * \param[in] iface local index of the face.
 * \return true if the iface face is a process boundary face.
 */
bool
ClassOctant::getPbound() const{
	return m_info[6]||m_info[7]||m_info[8]||m_info[9]||((m_dim-2)*(m_info[10]||m_info[11]));
};

/*! Get if the octant is new after a refinement.
 * \return true if the the octant is new after a refinement.
 */
bool
ClassOctant::getIsNewR() const{return m_info[12];};

/*! Get if the octant is new after a coarsening.
 * \return true if the the octant is new after a coarsening.
 */
bool
ClassOctant::getIsNewC() const{return m_info[13];};

/*! Get if the octant is a scary ghost octant.
 * \return true if the octant is a ghost octant.
 */
bool
ClassOctant::getIsGhost() const{return m_info[16];};

/*! Get if the octant is a balancing-blocked octant.
 * \return false if the octant has to be balanced.
 */
bool
ClassOctant::getNotBalance() const{return m_info[14];};

/*! Get if the octant has to be balanced.
 * \return true if the octant has to be balanced.
 */
bool
ClassOctant::getBalance() const{return (!m_info[14]);};

/*! Set the refinement marker of an octant.
 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
 */
void
ClassOctant::setMarker(int8_t marker){
	this->m_marker = marker;
};

/*! Set the balancing condition of an octant.
 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
 */
void
ClassOctant::setBalance(bool balance){
	m_info[14] = balance;
};

/*! Set the level of an octant.
 * \param[in] level New level of the octant.
 */
void
ClassOctant::setLevel(uint8_t level){
	this->m_level = level;
};

/*! Set the process boundary condition of an octant.
 * \param[in] balance Is the octant a process boundary octant?
 */
void
ClassOctant::setPbound(uint8_t face, bool flag){
	m_info[6+face] = flag;
};

// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

/*! Get the size of an octant in logical domain, i.e. the side length.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Size of octant.
 */
uint32_t
ClassOctant::getSize(int8_t maxlevel) const{
	uint32_t size = uint32_t(1<<(maxlevel-m_level));
	return size;
};

/*! Get the area of an octant in logical domain .
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Area of octant.
 */
uint64_t
ClassOctant::getArea(int8_t maxlevel) const{
	uint64_t area = uint64_t(pow(double(getSize(maxlevel)),double(m_dim-1)));
	return area;
};

/*! Get the volume of an octant in logical domain.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Volume of octant.
 */
uint64_t
ClassOctant::getVolume(int8_t maxlevel) const{
	uint64_t volume = uint64_t(pow(double(getSize(maxlevel)),double(m_dim)));
	return volume;
};

// =================================================================================== //

/*! Get the coordinates of the center of an octant in logical domain.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Array[3] with the coordinates of the center of octant.
 */
darray3
ClassOctant::getCenter(int8_t maxlevel) const{
	double	dh;
	darray3 center;

	dh = double(getSize(maxlevel))/2.0;
	center[0] = (double)m_x + dh;
	center[1] = (double)m_y + dh;
	center[2] = (double)m_z + dh;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a face of an octant in logical domain.
 * \param[in] iface Local index of the face
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Array[3] with the coordinates of the center of the octant face.
 */
darray3
ClassOctant::getFaceCenter(uint8_t iface, int8_t maxlevel) const{
	double	dh_2;
	darray3 center;

	dh_2 = double(getSize(maxlevel))/2.0;
	uint8_t nf = m_dim*2;
	if (iface < nf){
		center[0] = (double)m_x + (double)sm_CoeffFaceCenter[iface][0] * dh_2;
		center[1] = (double)m_y + (double)sm_CoeffFaceCenter[iface][1] * dh_2;
		center[2] = (double)m_z + (double)sm_CoeffFaceCenter[iface][2] * dh_2;
	}
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a edge of an octant in logical domain.
 * \param[in] iedge Local index of the edge
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Array[3] with the coordinates of the center of the octant edge.
 */
darray3
ClassOctant::getEdgeCenter(uint8_t iedge, int8_t maxlevel) const{
	double	dh_2;
	darray3 center;

	dh_2 = double(getSize(maxlevel))/2.0;
	center[0] = (double)m_x + (double)sm_CoeffEdgeCenter[iedge][0] * dh_2;
	center[1] = (double)m_y + (double)sm_CoeffEdgeCenter[iedge][1] * dh_2;
	center[2] = (double)m_z + (double)sm_CoeffEdgeCenter[iedge][2] * dh_2;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \param[out] nodes Vector of arrays [nnodes][3] with the coordinates of the nodes of octant.
 * \param[in] maxlevel Maximum refinement level of the octree.
 */
void
ClassOctant::getNodes(u32arr3vector & nodes, int8_t maxlevel) const{
	uint8_t		i, cx, cy, cz;
	uint32_t	dh;
	uint8_t nn = 1<<m_dim;

	dh = getSize(maxlevel);
	nodes.resize(nn);

	for (i = 0; i < nn; i++){
		cx = uint8_t(i%2);
		cy = uint8_t((i-4*(i/4))/2);
		cz = uint8_t(i/4);
		nodes[i][0] = m_x + cx*dh;
		nodes[i][1] = m_y + cy*dh;
		nodes[i][2] = m_z + cz*dh;
	}
};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[out] node Array[3] with the logical coordinates of the node of the octant.
 * \param[in] maxlevel Maximum refinement level of the octree.
 */
void		ClassOctant::getNode(u32array3 & node, uint8_t inode, int8_t maxlevel) const{
	uint8_t		cx, cy, cz;
	uint32_t	dh;

	dh = getSize(maxlevel);
	cx = inode%2;
	cy = (inode-4*(inode/4))/2;
	cz = inode/4;
	node[0] = m_x + cx*dh;
	node[1] = m_y + cy*dh;
	node[2] = m_z + cz*dh;

};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Array[3] with the logical coordinates of the node of the octant.
 */
u32array3		ClassOctant::getNode(uint8_t inode, int8_t maxlevel) const{
	u32array3 	node;
	uint8_t		cx, cy, cz;
	uint32_t	dh;

	dh = getSize(maxlevel);
	cx = inode%2;
	cy = (inode-4*(inode/4))/2;
	cz = inode/4;
	node[0] = m_x + cx*dh;
	node[1] = m_y + cy*dh;
	node[2] = m_z + cz*dh;
	return node;
};

/*! Get the normal of a face of an octant in logical domain.
 * \param[in] iface Index of the face for normal computing.
 * \param[out] normal Array[3] with components (with z=0) of the normal of face.
 * \param[in] maxlevel Maximum refinement level of the octree.
 */
void		ClassOctant::getNormal(uint8_t & iface, i8array3 & normal, int8_t (&normals)[6][3], int8_t maxlevel) const{
	uint8_t		i;
	for (i = 0; i < 3; i++){
		normal[i] = normals[iface][i];
	}
};


/** Compute the Morton index of the octant (without level).
 * \return morton Morton index of the octant.
 */
uint64_t	ClassOctant::computeMorton() const{
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->m_x,this->m_y,this->m_z);
	return morton;
};


/** Compute the Morton index of the octant (without level).
 * \return Morton index of the octant.
 */
uint64_t	ClassOctant::computeMorton(){
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->m_x,this->m_y,this->m_z);
	return morton;
};

// =================================================================================== //
// OTHER METHODS
// =================================================================================== //

/** Build the last descendant octant of this octant.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Last descendant octant.
 */
ClassOctant	ClassOctant::buildLastDesc(int8_t & maxlevel){
	uint32_t delta[3] = {0,0,0};
	for (int i=0; i<m_dim; i++){
		delta[i] = (uint32_t)(1 << (maxlevel - m_level)) - 1;
	}

	ClassOctant last_desc(m_dim, maxlevel, m_x+delta[0], m_y+delta[1], m_z+delta[2]);
	return last_desc;
};

// =================================================================================== //

/** Build the father octant of this octant.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Father octant.
 */
ClassOctant	ClassOctant::buildFather(int8_t & maxlevel){
	uint32_t delta[3];
	uint32_t xx[3];
	xx[0] = m_x;
	xx[1] = m_y;
	xx[2] = m_z;
	delta[2] = 0;
	for (int i=0; i<m_dim; i++){
		delta[i] = xx[i]%(uint32_t(1 << (maxlevel - max(0,(m_level-1)))));
	}
	ClassOctant father(m_dim, max(0,m_level-1), m_x-delta[0], m_y-delta[1], m_z-delta[2]);
	return father;
};

// =================================================================================== //

/** Builds children of octant.
 * \param[in] maxlevel Maximum refinement level of the octree.
 *   \return Ordered (by Z-index) vector of children[nchildren] (info update)
 */
vector< ClassOctant >	ClassOctant::buildChildren(int8_t & maxlevel){
	uint8_t xf,yf,zf;
	int nchildren = 1<<m_dim;

	if (this->m_level < maxlevel){
		vector< ClassOctant > children(nchildren);
		for (int i=0; i<nchildren; i++){
			switch (i) {
			case 0 :
			{
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
				ClassOctant oct(*this);
				oct.setMarker(max(0,oct.m_marker-1));
				oct.setLevel(oct.m_level+1);
				oct.m_info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
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
		vector< ClassOctant > children(0);
		return children;
	}
};

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant
 * throught face iface (sizehf=0 if boundary octant).
 * \param[in] iface Local index of the face target.
 * \param[out] sizehf Number of possible neighbours.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> ClassOctant::computeHalfSizeMorton(uint8_t iface, uint32_t & sizehf, int8_t & maxlevel){
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i,cx,cy,cz;
	int nchildren = 1<<m_dim;

	nneigh = (m_level < maxlevel) ? nchildren/2 : 1;
	dh = (m_level < maxlevel) ? getSize(maxlevel)/2 : getSize(maxlevel);
	dh2 = getSize(maxlevel);

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
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> ClassOctant::computeMinSizeMorton(uint8_t iface, const uint8_t & maxdepth, uint32_t & sizem, int8_t & maxlevel){
	uint32_t dh,dh2;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;

	nneigh = (m_level < maxlevel) ? uint32_t(1<<((m_dim-1)*(maxdepth-m_level))) : 1;
	dh = (m_level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	nline = uint32_t(1<<(maxdepth-m_level));

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
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> ClassOctant::computeVirtualMorton(uint8_t iface, const uint8_t & maxdepth, uint32_t & sizeneigh, int8_t & maxlevel){
	vector<uint64_t> Morton;
	if (getNotBalance()){
		return computeMinSizeMorton(iface,
				maxdepth,
				sizeneigh, maxlevel);
	}
	else{
		return computeHalfSizeMorton(iface,
				sizeneigh, maxlevel);
	}
};

/*! Computes Morton index (without level) of "n=sizehf" half-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * edge iedge (sizehf=0 if boundary octant)
 * \param[in] iedge Local index of the edge target.
 * \param[out] sizehf Number of possible neighbours.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> ClassOctant::computeEdgeHalfSizeMorton(uint8_t iedge, uint32_t & sizehf, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i,cx,cy,cz;
	uint8_t iface1, iface2;

	nneigh = (m_level < maxlevel) ? 2 : 1;
	dh = (m_level < maxlevel) ? getSize(maxlevel)/2 : getSize(maxlevel);
	dh2 = getSize(maxlevel);
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
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t> 		ClassOctant::computeEdgeMinSizeMorton(uint8_t iedge, const uint8_t & maxdepth, uint32_t & sizem, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){
	uint32_t dh,dh2;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;
	uint8_t iface1, iface2;


	nneigh = (m_level < maxlevel) ? uint32_t(1<<(maxdepth-m_level)) : 1;
	dh = (m_level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	nline = uint32_t(1<<(maxdepth-m_level));
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
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] edgeface Local edge-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
vector<uint64_t>		ClassOctant::computeEdgeVirtualMorton(uint8_t iedge, const uint8_t & maxdepth, uint32_t & sizeneigh, uint8_t balance_codim, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){

	if(!getNotBalance() && balance_codim > 1){
		return computeEdgeHalfSizeMorton(iedge,
				sizeneigh, maxlevel, edgeface);
	}
	else{
		return computeEdgeMinSizeMorton(iedge,
				maxdepth, sizeneigh, maxlevel, edgeface);
	}
};

/*! Computes Morton index (without level) of "n=sizem" min-size
 * (or same size if level=maxlevel) possible neighbours of octant throught
 * node inode (sizem=0 if boundary octant)
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizem Number of possible neighbours (1 or 0).
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] nodeface Local node-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
uint64_t 		ClassOctant::computeNodeMinSizeMorton(uint8_t inode, const uint8_t & maxdepth,uint32_t & sizem, int8_t & maxlevel, uint8_t (&nodeface)[8][3]){

	uint32_t dh,dh2;
	uint32_t nneigh;
	int8_t cx,cy,cz;
	uint8_t iface[3];

	nneigh = 1;
	dh = (m_level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	for (int i=0; i<m_dim; i++){
		iface[i] = nodeface[inode][i];
	}

	if (m_info[iface[0]] || m_info[iface[1]] || m_info[iface[m_dim-1]]){
		sizem = 0;
		return this->computeMorton();
	}
	else{
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
		}
		sizem = nneigh;
		return Morton;
	}

};

/*! Computes Morton index (without level) of possible (virtual) neighbours of octant
 * throught inode. Uses min-size method (sizeneigh=0 if boundary octant).
 * \param[in] inode Local index of the node target.
 * \param[in] maxdepth Maximum refinement level currently reached in the octree.
 * \param[out] sizeneigh Number of possible neighbours (1).
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] nodeface Local node-face connectivity.
 * \return Vector of neighbours morton numbers.
 */
uint64_t 		ClassOctant::computeNodeVirtualMorton(uint8_t inode, const uint8_t & maxdepth, uint32_t & sizeneigh, int8_t & maxlevel, uint8_t (&nodeface)[8][3]){

	return computeNodeMinSizeMorton(inode, maxdepth,
			sizeneigh, maxlevel, nodeface);

 };



