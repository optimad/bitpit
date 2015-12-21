// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classOctant.hpp"

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

classOctant::classOctant(){
	dim = 2;
	x = y = z = 0;
	level = 0;
	marker = 0;
	for (uint8_t i=0; i<4; i++){
		info[i] = true;
	}
	info[14] = false;
};

classOctant::classOctant(uint8_t dim_){
	dim = dim_;
	x = y = z = 0;
	level = 0;
	marker = 0;
	uint8_t nf = dim*2;
	for (uint8_t i=0; i<nf; i++){
		info[i] = true;
	}
	info[14] = false;
};

classOctant::classOctant(uint8_t dim_, uint8_t level_, int32_t x_, int32_t y_, int32_t z_){
	dim = dim_;
	x = x_;
	y = y_;
	z = (dim-2)*z_;
	level = level_;
	marker = 0;
	if (level==0){
		uint8_t nf = dim*2;
		for (uint8_t i=0; i<nf; i++){
			info[i] = true;
		}
	}
	info[14] = false;
};

classOctant::classOctant(bool bound, uint8_t dim_, uint8_t level_, int32_t x_, int32_t y_, int32_t z_){
	dim = dim_;
	x = x_;
	y = y_;
	z = (dim-2)*z_;
	level = level_;
	marker = 0;
	if (level==0){
		uint8_t nf = dim*2;
		for (uint8_t i=0; i<nf; i++){
			info[i] = bound;
		}
	}
	info[14] = false;
};

classOctant::classOctant(const classOctant &octant){
	dim = octant.dim;
	x = octant.x;
	y = octant.y;
	z = octant.z;
	level = octant.level;
	marker = octant.marker;
	info = octant.info;
};

/*! Check if two octants are equal (no check on info)
 */
bool classOctant::operator ==(const classOctant & oct2){
	bool check = true;
	check = check && (dim == oct2.dim);
	check = check && (x == oct2.x);
	check = check && (y == oct2.y);
	check = check && (z == oct2.z);
	check = check && (level == oct2.level);
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
uint32_t	classOctant::getDim() const{return dim;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate X of node 0.
 */
uint32_t	classOctant::getX() const{return x;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Y of node 0.
 */
uint32_t	classOctant::getY() const{return y;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinate Z of node 0.
 */
uint32_t	classOctant::getZ() const{return z;};

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \return Coordinates of node 0.
 */
u32array3	classOctant::getCoord(){

	u32array3 coord;

	coord[0] = x;
	coord[1] = y;
	coord[2] = z;

	return coord;
};

/*! Get the level of an octant.
 * \return Level of octant.
 */
uint8_t		classOctant::getLevel() const{return level;};

/*! Get the refinement marker of an octant.
 * \return Marker of octant.
 */
int8_t		classOctant::getMarker() const{return marker;};

/*! Get the bound flag on an octant face.
 * \param[in] iface local index of the face.
 * \return true if the iface face is a boundary face.
 */
bool		classOctant::getBound(uint8_t face) const{
	return info[face];
};

//private:
void		classOctant::setBound(uint8_t face) {
	info[face] = true;
};

//public:
/*! Get the pbound flag on an octant face.
 * \param[in] iface local index of the face.
 * \return true if the iface face is a process boundary face.
 */
bool		classOctant::getPbound(uint8_t face) const{
	return info[6+face];
};

/*! Get if the octant is new after a refinement.
 * \return true if the the octant is new after a refinement.
 */
bool		classOctant::getIsNewR() const{return info[12];};

/*! Get if the octant is new after a coarsening.
 * \return true if the the octant is new after a coarsening.
 */
bool		classOctant::getIsNewC() const{return info[13];};

/*! Get if the octant is a scary ghost octant.
 * \return true if the octant is a ghost octant.
 */
bool		classOctant::getIsGhost() const{return info[16];};

/*! Get if the octant is a balancing-blocked octant.
 * \return false if the octant has to be balanced.
 */
bool		classOctant::getNotBalance() const{return info[14];};

/*! Get if the octant has to be balanced.
 * \return true if the octant has to be balanced.
 */
bool		classOctant::getBalance() const{return (!info[14]);};

/*! Set the refinement marker of an octant.
 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
 */
void		classOctant::setMarker(int8_t marker){
	this->marker = marker;
};

/*! Set the balancing condition of an octant.
 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
 */
void		classOctant::setBalance(bool balance){
	info[14] = balance;
};

//private:
void		classOctant::setLevel(uint8_t level){
	this->level = level;
};

void 		classOctant::setPbound(uint8_t face, bool flag){
	info[6+face] = flag;
};

// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

/*! Get the size of an octant in logical domain, i.e. the side length.
 * \return Size of octant.
 */
uint32_t	classOctant::getSize(int8_t maxlevel) const{
	uint32_t size = uint32_t(1<<(maxlevel-level));
	return size;
};

/*! Get the area of an octant in logical domain .
 * \return Area of octant.
 */
uint64_t	classOctant::getArea(int8_t maxlevel) const{
	uint64_t area = uint64_t(pow(double(getSize(maxlevel)),double(dim-1)));
	return area;
};

/*! Get the volume of an octant in logical domain.
 * \return Volume of octant.
 */
uint64_t	classOctant::getVolume(int8_t maxlevel) const{
	uint64_t volume = uint64_t(pow(double(getSize(maxlevel)),double(dim)));
	return volume;
};

// =================================================================================== //

/*! Get the coordinates of the center of an octant in logical domain.
 * \return Array[3] with the coordinates of the center of octant.
 */
vector<double>	classOctant::getCenter(int8_t maxlevel) const{
	double	dh;

	dh = double(getSize(maxlevel))/2.0;
	vector<double> center(3);

	center[0] = (double)x + dh;
	center[1] = (double)y + dh;
	center[2] = (double)z + dh;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a face of an octant in logical domain.
 * \return Array[3] with the coordinates of the center of the octant face.
 */
vector<double>	classOctant::getFaceCenter(uint8_t iface, int8_t maxlevel) const{
	double	dh_2;

	int A[6][3] = { {0,1,1} , {2,1,1} , {1,0,1} , {1,2,1} , {1,1,0} , {1,1,2} };

	dh_2 = double(getSize(maxlevel))/2.0;
	vector<double> center(3);
	uint8_t nf = dim*2;
	if (iface < nf){
		center[0] = (double)x + (double)A[iface][0] * dh_2;
		center[1] = (double)y + (double)A[iface][1] * dh_2;
		center[2] = (double)z + (double)A[iface][2] * dh_2;
	}
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the center of a edge of an octant in logical domain.
 * \return Array[3] with the coordinates of the center of the octant edge.
 */
vector<double>	classOctant::getEdgeCenter(uint8_t iedge, int8_t maxlevel) const{
	double	dh_2;

	int A[12][3] = { {0,1,0},{2,1,0},{1,0,0},{1,2,0},{0,0,1},{2,0,1},{0,2,1},{2,2,1},{0,1,2},{2,1,2},{1,0,2},{1,2,2} };//{ {0,1,1} , {2,1,1} , {1,0,1} , {1,2,1} , {1,1,0} , {1,1,2} };

	dh_2 = double(getSize(maxlevel))/2.0;
	vector<double> center(3);

	center[0] = (double)x + (double)A[iedge][0] * dh_2;
	center[1] = (double)y + (double)A[iedge][1] * dh_2;
	center[2] = (double)z + (double)A[iedge][2] * dh_2;
	return center;
};

// =================================================================================== //

/*! Get the coordinates of the nodes of an octant in logical domain.
 * \param[out] nodes Vector[nnodes][3] with the coordinates of the nodes of octant.
 */
void		classOctant::getNodes(u32vector2D & nodes, int8_t maxlevel) const{
	uint8_t		i, cx, cy, cz;
	uint32_t	dh;

	dh = getSize(maxlevel);
	nodes.clear();
	nodes.resize(1<<dim);

	uint8_t nn = 1<<dim;
	for (i = 0; i < nn; i++){
		nodes[i].resize(3);
		cx = uint8_t(i%2);
		cy = uint8_t((i-4*(i/4))/2);
		cz = uint8_t(i/4);
		nodes[i][0] = x + cx*dh;
		nodes[i][1] = y + cy*dh;
		nodes[i][2] = z + cz*dh;
		u32vector(nodes[i]).swap(nodes[i]);
	}
	u32vector2D(nodes).swap(nodes);
};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[out] node dim-vector with the logical coordinates of the node of the octant.
 */
void		classOctant::getNode(u32vector & node, uint8_t inode, int8_t maxlevel) const{
	uint8_t		cx, cy, cz;
	uint32_t	dh;

	dh = getSize(maxlevel);
	node.clear();
	node.resize(3);
	cx = inode%2;
	cy = (inode-4*(inode/4))/2;
	cz = inode/4;
	node[0] = x + cx*dh;
	node[1] = y + cy*dh;
	node[2] = z + cz*dh;

};

/*! Get the coordinates of a nodes of an octant in logical domain.
 * \param[in] inode Local index of the node
 * \param[out] node dim-vector with the logical coordinates of the node of the octant.
 */
u32vector		classOctant::getNode(uint8_t inode, int8_t maxlevel) const{
	u32vector node;
	uint8_t		cx, cy, cz;
	uint32_t	dh;

	dh = getSize(maxlevel);
	node.clear();
	node.resize(3);
	cx = inode%2;
	cy = (inode-4*(inode/4))/2;
	cz = inode/4;
	node[0] = x + cx*dh;
	node[1] = y + cy*dh;
	node[2] = z + cz*dh;
	return node;
};

//TODO getNormal in localtree
/*! Get the normal of a face of an octant in logical domain.
 * \param[in] iface Index of the face for normal computing.
 * \param[out] normal Vector[3] with components (with z=0) of the normal of face.
 */
void		classOctant::getNormal(uint8_t & iface,
		vector<int8_t> & normal, int8_t (&normals)[6][3], int8_t maxlevel) const{
	uint8_t		i;

	normal.clear();
	normal.resize(3);
	for (i = 0; i < 3; i++){
		normal[i] = normals[iface][i];
	}
	vector<int8_t>(normal).swap(normal);
};


/** Compute the Morton index of the octant (without level).
 * \return morton Morton index of the octant.
 */
uint64_t	classOctant::computeMorton() const{
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->x,this->y,this->z);
	return morton;
};


/** Compute the Morton index of the octant (without level).
 * \return morton Morton index of the octant.
 */
uint64_t	classOctant::computeMorton(){
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->x,this->y,this->z);
	return morton;
};

// =================================================================================== //
// OTHER METHODS
// =================================================================================== //

classOctant	classOctant::buildLastDesc(int8_t & maxlevel){
	uint32_t delta[3] = {0,0,0};
	for (int i=0; i<dim; i++){
		delta[i] = (uint32_t)(1 << (maxlevel - level)) - 1;
	}

	classOctant last_desc(dim, maxlevel, x+delta[0], y+delta[1], z+delta[2]);
	return last_desc;
};

// =================================================================================== //

classOctant	classOctant::buildFather(int8_t & maxlevel){
	uint32_t delta[3] = {0,0,0};
	uint32_t xx[3] = {x, y, z};
	for (int i=0; i<dim; i++){
		delta[i] = xx[i]%(uint32_t(1 << (maxlevel - (level-1))));
	}
	classOctant father(dim, level-1, x-delta[0], y-delta[1], z-delta[2]);
	return father;
};

// =================================================================================== //

/** Builds children of octant.
 *   \return Ordered (by Z-index) vector of children[nchildren] (info update)
 */
vector< classOctant >	classOctant::buildChildren(int8_t & maxlevel){
	uint8_t xf,yf,zf;
	int nchildren = 1<<dim;
//	int nfaces = 2*dim;

	if (this->level < maxlevel){
		vector< classOctant > children(nchildren);
		for (int i=0; i<nchildren; i++){
			switch (i) {
			case 0 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=5;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[0] = oct;
			}
			break;
			case 1 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.x += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=5;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[1] = oct;
			}
			break;
			case 2 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.y += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=5;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[2] = oct;
			}
			break;
			case 3 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.x += dh;
				oct.y += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=5;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[3] = oct;
			}
			break;
			case 4 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.z += dh;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=4;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[4] = oct;
			}
			break;
			case 5 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.x += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=4;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[5] = oct;
			}
			break;
			case 6 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.y += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=4;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[6] = oct;
			}
			break;
			case 7 :
			{
				classOctant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize(maxlevel);
				oct.x += dh;
				oct.y += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=4;
				oct.info[xf] = oct.info[xf+6] = false;
				oct.info[yf] = oct.info[yf+6] = false;
				oct.info[zf] = oct.info[zf+6] = false;
				children[7] = oct;
			}
			break;
			}
		}
		return children;
	}
	else{
		vector< classOctant > children(0);
		//writeLog("Max level reached ---> No Children Built");
		return children;
	}
};

vector<uint64_t> 		classOctant::computeHalfSizeMorton(uint8_t iface, 			// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
		uint32_t & sizehf, int8_t & maxlevel){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i,cx,cy,cz;
	int nchildren = 1<<dim;

	nneigh = (level < maxlevel) ? nchildren/2 : 1;
	dh = (level < maxlevel) ? getSize(maxlevel)/2 : getSize(maxlevel);
	dh2 = getSize(maxlevel);

	if (info[iface]){
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
				cz = (dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i==1)||(i==3);
				cz = (dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (dim == 3) && ((i==2)||(i==3));
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2);
			}
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}


};

vector<uint64_t>		classOctant::computeMinSizeMorton(uint8_t iface, 			// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
		const uint8_t & maxdepth,	// possible neighbours of octant throught face iface (sizem=0 if boundary octant)
		uint32_t & sizem, int8_t & maxlevel){
	uint32_t dh,dh2;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;

	//		nneigh = (level < CG::MAX_LEVEL) ? uint32_t(pow(2.0,double((dim-1)*(maxdepth-level)))) : 1;
	//		dh = (level < CG::MAX_LEVEL) ? uint32_t(pow(2.0,double(CG::MAX_LEVEL - maxdepth))) : getSize(maxlevel);
	//		dh2 = getSize(maxlevel);
	//		nline = uint32_t(pow(2.0,double((maxdepth-level))));
	nneigh = (level < maxlevel) ? uint32_t(1<<((dim-1)*(maxdepth-level))) : 1;
	dh = (level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	nline = uint32_t(1<<(maxdepth-level));

	if (info[iface]){
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
				cz = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cz = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cz = (i/nline);
				cx = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cz = (i/nline);
				cx = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i/nline);
				cx = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i/nline);
				cx = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2);
			}
		}
		break;
		}
		sizem = nneigh;
		//		sort(Morton,Morton+nneigh);
		sort(Morton.begin(), Morton.end());
		return Morton;
	}

};

vector<uint64_t> 		classOctant::computeVirtualMorton(uint8_t iface, 			// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
		const uint8_t & maxdepth,	// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
		uint32_t & sizeneigh, int8_t & maxlevel){
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

vector<uint64_t> 		classOctant::computeEdgeHalfSizeMorton(uint8_t iedge, 		// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
		uint32_t & sizehf, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
	uint32_t dh,dh2;
	uint32_t nneigh;
	uint32_t i,cx,cy,cz;
	uint8_t iface1, iface2;

	nneigh = (level < maxlevel) ? 2 : 1;
	dh = (level < maxlevel) ? getSize(maxlevel)/2 : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	iface1 = edgeface[iedge][0];
	iface2 = edgeface[iedge][1];

	if (info[iface1] || info[iface2]){
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
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = (i==1);
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = -1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = 1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = -1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = -1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 6 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = 1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 7 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = 1;
				cz = (i==1);
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 8 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = (i==1);
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 9 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = (i==1);
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 10 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = -1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 11 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1);
				cy = 1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh2*cz);
			}
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}


};

vector<uint64_t> 		classOctant::computeEdgeMinSizeMorton(uint8_t iedge, 		// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
		const uint8_t & maxdepth,	// possible neighbours of octant throught edge iedge (sizem=0 if boundary octant)
		uint32_t & sizem, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){
	uint32_t dh,dh2;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;
	uint8_t iface1, iface2;


	nneigh = (level < maxlevel) ? uint32_t(1<<(maxdepth-level)) : 1;
	dh = (level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	nline = uint32_t(1<<(maxdepth-level));
	iface1 = edgeface[iedge][0];
	iface2 = edgeface[iedge][1];

	if (info[iface1] || info[iface2]){
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
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = i;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = -1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = 1;
				cz = -1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = -1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = -1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 6 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = 1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 7 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = 1;
				cz = i;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh*cz);
			}
		}
		break;
		case 8 :
		{
			for (i=0; i<nneigh; i++){
				cx = -1;
				cy = i;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 9 :
		{
			for (i=0; i<nneigh; i++){
				cx = 1;
				cy = i;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 10 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = -1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
			}
		}
		break;
		case 11 :
		{
			for (i=0; i<nneigh; i++){
				cx = i;
				cy = 1;
				cz = 1;
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh2*cz);
			}
		}
		break;
		}
		sizem = nneigh;
		sort(Morton.begin(),Morton.end());
		return Morton;
	}
};

vector<uint64_t>		classOctant::computeEdgeVirtualMorton(uint8_t iedge, 		// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
		const uint8_t & maxdepth,	// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
		uint32_t & sizeneigh,
		uint8_t balance_codim, int8_t & maxlevel, uint8_t (&edgeface)[12][2]){

	//	return computeEdgeMinSizeMorton(iedge,
	//			maxdepth, sizeneigh, maxlevel, edgeface);

	if(!getNotBalance() && balance_codim > 1){
		return computeEdgeHalfSizeMorton(iedge,
				sizeneigh, maxlevel, edgeface);
	}
	else{
		return computeEdgeMinSizeMorton(iedge,
				maxdepth, sizeneigh, maxlevel, edgeface);
	}
};

uint64_t 		classOctant::computeNodeHalfSizeMorton(uint8_t inode, 		// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
		uint32_t & sizehf, int8_t & maxlevel, uint8_t (&nodeface)[8][3]){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
	uint32_t dh,dh2;
	uint32_t nneigh;
	int8_t cx,cy,cz;
	uint8_t iface[3];
	nneigh = 1;
	dh = (level < maxlevel) ? getSize(maxlevel)/2 : getSize(maxlevel);
	dh2 = getSize(maxlevel);

	for (int i=0; i<dim; i++){
		iface[i] = nodeface[inode][i];
	}

	if (info[iface[0]] || info[iface[1]] || info[iface[2]]){
		sizehf = 0;
		return this->computeMorton();
	}
	else{
		uint64_t Morton;
		switch (inode) {
		case 0 :
		{
			cx = -1;
			cy = -1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
		}
		break;
		case 1 :
		{
			cx = 1;
			cy = -1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
		}
		break;
		case 2 :
		{
			cx = -1;
			cy = 1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
		}
		break;
		case 3 :
		{
			cx = 1;
			cy = 1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh*cz);
		}
		break;
		case 4 :
		{
			cx = -1;
			cy = -1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
		}
		break;
		case 5 :
		{
			cx = 1;
			cy = -1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh2*cz);
		}
		break;
		case 6 :
		{
			cx = -1;
			cy = 1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh2*cz);
		}
		break;
		case 7 :
		{
			cx = 1;
			cy = 1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh2*cz);
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}
};

// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
// possible neighbours of octant throught face iface (sizem=0 if boundary octant)
uint64_t 		classOctant::computeNodeMinSizeMorton(uint8_t inode,
		const uint8_t & maxdepth,uint32_t & sizehf,
		int8_t & maxlevel, uint8_t (&nodeface)[8][3]){

	uint32_t dh,dh2;
	uint32_t nneigh;
	int8_t cx,cy,cz;
	uint8_t iface[3];

	nneigh = 1;
	dh = (level < maxlevel) ? uint32_t(1<<(maxlevel - maxdepth)) : getSize(maxlevel);
	dh2 = getSize(maxlevel);
	for (int i=0; i<dim; i++){
		iface[i] = nodeface[inode][i];
	}

	if (info[iface[0]] || info[iface[1]] || info[iface[dim-1]]){
		sizehf = 0;
		return this->computeMorton();
	}
	else{
		uint64_t Morton;
		switch (inode) {
		case 0 :
		{
			cx = -1;
			cy = -1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
		}
		break;
		case 1 :
		{
			cx = 1;
			cy = -1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
		}
		break;
		case 2 :
		{
			cx = -1;
			cy = 1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
		}
		break;
		case 3 :
		{
			cx = 1;
			cy = 1;
			cz = -1*(dim-2);
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh*cz);
		}
		break;
		case 4 :
		{
			cx = -1;
			cy = -1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2*cz);
		}
		break;
		case 5 :
		{
			cx = 1;
			cy = -1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh2*cz);
		}
		break;
		case 6 :
		{
			cx = -1;
			cy = 1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh2*cz);
		}
		break;
		case 7 :
		{
			cx = 1;
			cy = 1;
			cz = 1;
			Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh2*cy,this->z+dh2*cz);
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}

};

// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
uint64_t 		classOctant::computeNodeVirtualMorton(uint8_t inode,
		const uint8_t & maxdepth, uint32_t & sizeneigh, int8_t & maxlevel,
		uint8_t (&nodeface)[8][3]){

	//		if (getNotBalance()){
	return computeNodeMinSizeMorton(inode, maxdepth,
			sizeneigh, maxlevel, nodeface);
	//		}
	//		else{
	//			return computeNodeHalfSizeMorton(inode,
	//					sizeneigh);
	//		}
};



