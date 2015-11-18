/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Octant class definition - 3D specialization
 *
 *	Octants are the grid elements of PABLO. In the logical domain octants are, in 3-dimensional case,
 *	cubes with size function of their level. Each octant has 8 nodes, 12 edge and 6 faces ordered with Z-order as shown below.
 *
 *	\image html PabloOctant3D.png
 *
 *
 *	The main feature of each octant are:
 *	- x,y,z        : coordinates of the node 0 of the octant;
 *	- Morton index : classical Morton index defined anly by the coordinates (info about level used additionally for equality operator);
 *	- marker       : refinement marker can assume negative, positive or zero values, wich mean
 *	a coarsening, refinement and none adaptation respectively;
 *	- level        : octant level in the octree, zero for the first upper level.
 *	- balance      : flag to fix if the octant has to 2:1 balanced with respect to its face neighbours.
 *
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //
template<>
class Class_Octant<3>{
	// ------------------------------------------------------------------------------- //
	// FRIENDSHIPS ------------------------------------------------------------------- //
	template<int dim> friend class Class_Local_Tree;
	template<int dim> friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant<3> > 	OctantsType;
	typedef vector<double>				dvector;
	typedef vector<uint32_t>			u32vector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	uint32_t  	x;				/**<  Coordinate x */
	uint32_t  	y;				/**<  Coordinate y */
	uint32_t  	z;				/**<  Coordinate z */
	uint8_t   	level;			/**< Refinement level (0=root) */
	int8_t    	marker;			/**< Set for Refinement(m>0) or Coarsening(m<0) |m|-times */
	bitset<16>  info;			/**< -Info[0..5] : true if 0..3 face is a boundary face [bound] \n
								-Info[6..11]: true if 0..3 face is a process boundary face [pbound] \n
								-Info[12/13]: true if octant is new after refinement/coarsening \n
								-Info[14]   : true if balancing is not required for this octant \n
								-Info[15]   : Aux (before : true if octant is a scary ghost) */


	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS AND OPERATORS----------------------------------------------------- //

public:
	Class_Octant(){
		x = y = z = 0;
		level = 0;
		marker = 0;
		for (int i=0; i<global3D.nfaces; i++){
			info[i] = true;
		}
	};

	Class_Octant(int8_t level, int32_t x, int32_t y, int32_t z){
		this->x = x;
		this->y = y;
		this->z = z;
		this->level = level;
		marker = 0;
		if (level==0){
			for (int i=0; i<global3D.nfaces; i++){
				info[i] = true;
			}
		}

	};

	Class_Octant(int8_t level, int32_t x, int32_t y, int32_t z, bool bound){
		this->x = x;
		this->y = y;
		this->z = z;
		this->level = level;
		marker = 0;
		if (level==0){
			for (int i=0; i<global3D.nfaces; i++){
				info[i] = bound;
			}
		}

	};


	Class_Octant(const Class_Octant<3> &octant){
		x = octant.x;
		y = octant.y;
		z = octant.z;
		level = octant.level;
		marker = octant.marker;
		info = octant.info;
	};

	/*! Check if two octants are equal (no check on info)
	 */
	bool operator ==(const Class_Octant<3> & oct2){
		bool check = true;
		check = check && (x == oct2.x);
		check = check && (y == oct2.y);
		check = check && (z == oct2.z);
		check = check && (level == oct2.level);
		return check;
	}

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

public:
	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate X of node 0.
	 */
	uint32_t	getX() const{return x;};

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate Y of node 0.
	 */
	uint32_t	getY() const{return y;};

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \return Coordinate Z of node 0.
	 */
	uint32_t	getZ() const{return z;};

	/*! Get the level of an octant.
	 * \return Level of octant.
	 */
	uint8_t		getLevel() const{return level;};

	/*! Get the refinement marker of an octant.
	 * \return Marker of octant.
	 */
	int8_t		getMarker() const{return marker;};

	/*! Get the bound flag on an octant face.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a boundary face.
	 */
	bool		getBound(uint8_t face) const{
		return info[face];
	};

private:
	void		setBound(uint8_t face) {
		info[face] = true;
	};

public:
	/*! Get the pbound flag on an octant face.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a process boundary face.
	 */
	bool		getPbound(uint8_t face) const{
		return info[global3D.nfaces+face];
	};

	/*! Get if the octant is new after a refinement.
	 * \return true if the the octant is new after a refinement.
	 */
	bool		getIsNewR() const{return info[12];};

	/*! Get if the octant is new after a coarsening.
	 * \return true if the the octant is new after a coarsening.
	 */
	bool		getIsNewC() const{return info[13];};

	/*! Get if the octant is a balancing-blocked octant.
	 * \return false if the octant has to be balanced.
	 */
	bool		getNotBalance() const{return info[14];};

	/*! Get if the octant has to be balanced.
	 * \return true if the octant has to be balanced.
	 */
	bool		getBalance() const{return (!info[10]);};

	/*! Set the refinement marker of an octant.
	 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
	 */
	void		setMarker(int8_t marker){
		this->marker = marker;
	};

	/*! Set the balancing condition of an octant.
	 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
	 */
	void		setBalance(bool balance){
		info[14] = balance;
	};

private:
	void		setLevel(uint8_t level){
		this->level = level;
	};

	void 		setPbound(uint8_t face, bool flag){
		info[global3D.nfaces+face] = flag;
	};

	//-------------------------------------------------------------------------------- //
	// Other Get/Set methods --------------------------------------------------------- //

public:
	/*! Get the size of an octant in logical domain, i.e. the side length.
	 * \return Size of octant.
	 */
	uint32_t	getSize() const{
		uint32_t size = uint32_t(pow(double(2),double(MAX_LEVEL_3D-level)));
		return size;
	};

	/*! Get the area of an octant in logical domain .
	 * \return Area of octant.
	 */
	uint64_t	getArea() const{
		uint64_t area = uint64_t(pow(double(getSize()),2.0));
		return area;
	};

	/*! Get the volume of an octant in logical domain.
	 * \return Volume of octant.
	 */
	uint64_t	getVolume() const{
		uint64_t volume = uint64_t(pow(double(getSize()),3.0));
		return volume;
	};

	// ------------------------------------------------------------------------------- //

	/*! Get the coordinates of the center of an octant in logical domain.
	 * \return Vector[3] with the coordinates of the center of octant.
	 */
	dvector	getCenter(){
		double	dh;

		dh = double(getSize())/2.0;
		vector<double> center(3);

		center[0] = (double)x + dh;
		center[1] = (double)y + dh;
		center[2] = (double)z + dh;
		return center;
	};

	// ------------------------------------------------------------------------------- //

	/*! Get the coordinates of the center of a face of an octant in logical domain.
	 * \return Vector[3] with the coordinates of the center of the octant face.
	 */
	dvector	getFaceCenter(uint8_t iface){
		double	dh_2;

		int A[6][3] = { {0,1,1} , {2,1,1} , {1,0,1} , {1,2,1} , {1,1,0} , {1,1,2} };

		dh_2 = double(getSize())/2.0;
		vector<double> center(3);

		if (iface < global3D.nfaces){
			center[0] = (double)x + (double)A[iface][0] * dh_2;
			center[1] = (double)y + (double)A[iface][1] * dh_2;
			center[2] = (double)z + (double)A[iface][2] * dh_2;
		}
		return center;
	};

	// ------------------------------------------------------------------------------- //

	/*! Get the coordinates of the center of a edge of an octant in logical domain.
	 * \return Vector[3] with the coordinates of the center of the octant edge.
	 */
	dvector	getEdgeCenter(uint8_t iedge){
		double	dh_2;

		int A[12][3] = { {0,1,0},{2,1,0},{1,0,0},{1,2,0},{0,0,1},{2,0,1},{0,2,1},{2,2,1},{0,1,2},{2,1,2},{1,0,2},{1,2,2} };//{ {0,1,1} , {2,1,1} , {1,0,1} , {1,2,1} , {1,1,0} , {1,1,2} };

		dh_2 = double(getSize())/2.0;
		vector<double> center(3);

		if (iedge < global3D.nedges){
			center[0] = (double)x + (double)A[iedge][0] * dh_2;
			center[1] = (double)y + (double)A[iedge][1] * dh_2;
			center[2] = (double)z + (double)A[iedge][2] * dh_2;
		}
		return center;
	};

	// ------------------------------------------------------------------------------- //

	/*! Get the coordinates of the nodes of an octant in logical domain.
	 * \param[out] nodes Vector[4][3] with the coordinates of the nodes of octant.
	 */
	void		getNodes(u32vector2D & nodes){
		uint8_t		i, cx, cy, cz;
		uint32_t	dh;

		dh = getSize();
		nodes.clear();
		nodes.resize(global3D.nnodes);

		for (i = 0; i < global3D.nnodes; i++){
			nodes[i].resize(3);
			cx = uint8_t(i%2);
			cy = uint8_t((i-4*(i/4))/2);
			cz = uint8_t(i/4);
			nodes[i][0] = x + cx*dh;
			nodes[i][1] = y + cy*dh;
			nodes[i][2] = z + cz*dh;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#else
			nodes[i].shrink_to_fit();
#endif
		}
#if defined(__INTEL_COMPILER) || defined(__ICC)
#else
		nodes.shrink_to_fit();
#endif
	};

	/*! Get the coordinates of a nodes of an octant in logical domain.
	 * \param[in] inode Local index of the node
	 * \param[out] node dim-vector with the logical coordinates of the node of the octant.
	 */
	void		getNode(u32vector & node, uint8_t inode){
		uint8_t		cx, cy, cz;
		uint32_t	dh;

		dh = getSize();
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
	u32vector		getNode(uint8_t inode){
		u32vector node;
		uint8_t		cx, cy, cz;
		uint32_t	dh;

		dh = getSize();
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

	// ------------------------------------------------------------------------------- //

	/*! Get the normal of a face of an octant in logical domain.
	 * \param[in] iface Index of the face for normal computing.
	 * \param[out] normal Vector[3] with components (with z=0) of the normal of face.
	 */
	void		getNormal(uint8_t & iface,
			vector<int8_t> & normal){
		uint8_t		i;

		normal.clear();
		normal.resize(3);
		for (i = 0; i < 3; i++){
			normal[i] = global3D.normals[iface][i];
		}
#if defined(__INTEL_COMPILER) || defined(__ICC)
#else
		normal.shrink_to_fit();
#endif
	};

	// ------------------------------------------------------------------------------- //

	/** Compute the Morton index of the octant (without level).
	 * \return morton Morton index of the octant.
	 */
	uint64_t	computeMorton() const{							uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->x,this->y,this->z);
	return morton;
	};

	// ------------------------------------------------------------------------------- //

	/** Compute the Morton index of the octant (without level).
	 * \return morton Morton index of the octant.
	 */
	uint64_t	computeMorton(){
		uint64_t morton = 0;
		morton = mortonEncode_magicbits(this->x,this->y,this->z);
		return morton;
	};

	// =================================================================================== //
	// Other methods													    			   //
	// =================================================================================== //

private:
	Class_Octant<3>	buildLastDesc(){								// Build last descendant of octant and return the last descendant octant (no info update)
		uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - level)) - 1;
		Class_Octant<3> last_desc(MAX_LEVEL_3D,x+delta,y+delta,z+delta);
		return last_desc;
	};

	// =================================================================================== //

	Class_Octant<3>	buildFather(){									// Build father of octant and return the father octant (no info update)
		uint32_t deltax = x%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - (level-1)))));
		uint32_t deltay = y%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - (level-1)))));
		uint32_t deltaz = z%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - (level-1)))));
		Class_Octant<3> father(level-1, x-deltax, y-deltay, z-deltaz);
		return father;
	};

	// =================================================================================== //

	//	 ------------------------------------------------------------------------------- //

	/** Builds children of octant.
	 *   \return Ordered (by Z-index) vector of children[nchildren] (info update)
	 */
	vector< Class_Octant<3> >	buildChildren(){
		uint8_t xf,yf,zf;

		if (this->level < MAX_LEVEL_3D){
			vector< Class_Octant<3> > children(global3D.nchildren);
			for (int i=0; i<global3D.nchildren; i++){
				switch (i) {
				case 0 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					// Update interior face bound and pbound
					xf=1; yf=3; zf=5;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[0] = oct;
				}
				break;
				case 1 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.x += dh;
					// Update interior face bound and pbound
					xf=0; yf=3; zf=5;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[1] = oct;
				}
				break;
				case 2 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.y += dh;
					// Update interior face bound and pbound
					xf=1; yf=2; zf=5;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[2] = oct;
				}
				break;
				case 3 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.x += dh;
					oct.y += dh;
					// Update interior face bound and pbound
					xf=0; yf=2; zf=5;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[3] = oct;
				}
				break;
				case 4 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.z += dh;
					// Update interior face bound and pbound
					xf=1; yf=3; zf=4;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[4] = oct;
				}
				break;
				case 5 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.x += dh;
					oct.z += dh;
					// Update interior face bound and pbound
					xf=0; yf=3; zf=4;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[5] = oct;
				}
				break;
				case 6 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.y += dh;
					oct.z += dh;
					// Update interior face bound and pbound
					xf=1; yf=2; zf=4;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[6] = oct;
				}
				break;
				case 7 :
				{
					Class_Octant<3> oct(*this);
					oct.setMarker(max(0,oct.marker-1));
					oct.setLevel(oct.level+1);
					oct.info[12]=true;
					uint32_t dh = oct.getSize();
					oct.x += dh;
					oct.y += dh;
					oct.z += dh;
					// Update interior face bound and pbound
					xf=0; yf=2; zf=4;
					oct.info[xf] = oct.info[xf+global3D.nfaces] = false;
					oct.info[yf] = oct.info[yf+global3D.nfaces] = false;
					oct.info[zf] = oct.info[zf+global3D.nfaces] = false;
					children[7] = oct;
				}
				break;
				}
			}
			return children;
		}
		else{
			vector< Class_Octant<3> > children(0);
			//writeLog("Max level reached ---> No Children Built");
			return children;
		}
	};

	// ------------------------------------------------------------------------------- //

	vector<uint64_t> 		computeHalfSizeMorton(uint8_t iface, 			// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
			uint32_t & sizehf){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
		uint32_t dh,dh2;
		uint32_t nneigh;
		uint32_t i,cx,cy,cz;

		nneigh = (level < MAX_LEVEL_3D) ? global3D.nchildren/2 : 1;
		dh = (level < MAX_LEVEL_3D) ? getSize()/2 : getSize();
		dh2 = getSize();

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
					cz = (i==2)||(i==3);
					Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
				}
			}
			break;
			case 1 :
			{
				for (i=0; i<nneigh; i++){
					cy = (i==1)||(i==3);
					cz = (i==2)||(i==3);
					Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
				}
			}
			break;
			case 2 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i==1)||(i==3);
					cz = (i==2)||(i==3);
					Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
				}
			}
			break;
			case 3 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i==1)||(i==3);
					cz = (i==2)||(i==3);
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

	// ------------------------------------------------------------------------------- //

	vector<uint64_t>		computeMinSizeMorton(uint8_t iface, 			// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
			const uint8_t & maxdepth,	// possible neighbours of octant throught face iface (sizem=0 if boundary octant)
			uint32_t & sizem){
		uint32_t dh,dh2;
		uint32_t nneigh, nline;
		uint32_t i,cx,cy,cz;

		nneigh = (level < MAX_LEVEL_3D) ? uint32_t(pow(2.0,double(2*(maxdepth-level)))) : 1;
		dh = (level < MAX_LEVEL_3D) ? uint32_t(pow(2.0,double(MAX_LEVEL_3D - maxdepth))) : getSize();
		dh2 = getSize();
		nline = uint32_t(pow(2.0,double((maxdepth-level))));

		if (info[iface]){
			//		uint64_t* Morton = new uint64_t[0];
			sizem = 0;
			//		return Morton;
			vector<uint64_t> Morton(0);
			return Morton;
		}
		else{
			vector<uint64_t> Morton(nneigh);
			switch (iface) {
			case 0 :
			{
				for (i=0; i<nneigh; i++){
					cy = (i/nline);
					cz = (i%nline);
					Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
				}
			}
			break;
			case 1 :
			{
				for (i=0; i<nneigh; i++){
					cy = (i/nline);
					cz = (i%nline);
					Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
				}
			}
			break;
			case 2 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i/nline);
					cz = (i%nline);
					Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
				}
			}
			break;
			case 3 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i/nline);
					cz = (i%nline);
					Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2,this->z+dh*cz);
				}
			}
			break;
			case 4 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i/nline);
					cy = (i%nline);
					Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z-dh);
				}
			}
			break;
			case 5 :
			{
				for (i=0; i<nneigh; i++){
					cx = (i/nline);
					cy = (i%nline);
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

	// ------------------------------------------------------------------------------- //

	vector<uint64_t> 		computeVirtualMorton(uint8_t iface, 			// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
			const uint8_t & maxdepth,	// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
			uint32_t & sizeneigh){
		vector<uint64_t> Morton;
		if (getNotBalance()){
			return computeMinSizeMorton(iface,
					maxdepth,
					sizeneigh);
		}
		else{
			return computeHalfSizeMorton(iface,
					sizeneigh);
		}
	};

	// ------------------------------------------------------------------------------- //

	vector<uint64_t> 		computeEdgeHalfSizeMorton(uint8_t iedge, 		// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
			uint32_t & sizehf){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
		uint32_t dh,dh2;
		uint32_t nneigh;
		uint32_t i,cx,cy,cz;
		uint8_t iface1, iface2;

		nneigh = (level < MAX_LEVEL_3D) ? 2 : 1;
		dh = (level < MAX_LEVEL_3D) ? getSize()/2 : getSize();
		dh2 = getSize();
		iface1 = global3D.edgeface[iedge][0];
		iface2 = global3D.edgeface[iedge][1];

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

	// ------------------------------------------------------------------------------- //

	vector<uint64_t> 		computeEdgeMinSizeMorton(uint8_t iedge, 		// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
			const uint8_t & maxdepth,	// possible neighbours of octant throught edge iedge (sizem=0 if boundary octant)
			uint32_t & sizem){
		uint32_t dh,dh2;
		uint32_t nneigh, nline;
		uint32_t i,cx,cy,cz;
		uint8_t iface1, iface2;


		nneigh = (level < MAX_LEVEL_3D) ? uint32_t(pow(2.0,double((maxdepth-level)))) : 1;
		dh = (level < MAX_LEVEL_3D) ? uint32_t(pow(2.0,double(MAX_LEVEL_3D - maxdepth))) : getSize();
		dh2 = getSize();
		nline = uint32_t(pow(2.0,double((maxdepth-level))));
		iface1 = global3D.edgeface[iedge][0];
		iface2 = global3D.edgeface[iedge][1];

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
			//sort(Morton,Morton+nneigh);
			sort(Morton.begin(),Morton.end());
			return Morton;
		}
	};

	// ------------------------------------------------------------------------------- //

	vector<uint64_t>		computeEdgeVirtualMorton(uint8_t iedge, 		// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
			const uint8_t & maxdepth,	// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
			uint32_t & sizeneigh,
			uint8_t balance_codim){
//		if (getNotBalance()){
		return computeEdgeMinSizeMorton(iedge,
					maxdepth,
					sizeneigh);
//		}
//		else{
//			return computeEdgeHalfSizeMorton(iedge,
//					sizeneigh);
//		}
		if(!getNotBalance() && balance_codim > 1){
			return computeEdgeHalfSizeMorton(iedge,
					sizeneigh);
		}
		else{
			return computeEdgeMinSizeMorton(iedge,
					maxdepth,
					sizeneigh);
		}
	};

	// ------------------------------------------------------------------------------- //


	uint64_t 		computeNodeHalfSizeMorton(uint8_t inode, 		// Computes Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
			uint32_t & sizehf){		// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)
		uint32_t dh,dh2;
		uint32_t nneigh;
		int8_t cx,cy,cz;
		uint8_t iface1, iface2, iface3;
		nneigh = 1;
		dh = (level < MAX_LEVEL_3D) ? getSize()/2 : getSize();
		dh2 = getSize();
		iface1 = global3D.nodeface[inode][0];
		iface2 = global3D.nodeface[inode][1];
		iface3 = global3D.nodeface[inode][2];

		if (info[iface1] || info[iface2] || info[iface3]){
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
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
			break;
			case 1 :
			{
				cx = 1;
				cy = -1;
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
			break;
			case 2 :
			{
				cx = -1;
				cy = 1;
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
			break;
			case 3 :
			{
				cx = 1;
				cy = 1;
				cz = -1;
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

	// ------------------------------------------------------------------------------- //

	uint64_t 		computeNodeMinSizeMorton(uint8_t inode, 		// Computes Morton index (without level) of "n=sizem" min-size (or same size if level=maxlevel)
			const uint8_t & maxdepth,	// possible neighbours of octant throught face iface (sizem=0 if boundary octant)
			uint32_t & sizehf){
		uint32_t dh,dh2;
		uint32_t nneigh;
		int8_t cx,cy,cz;
		uint8_t iface1, iface2, iface3;

		nneigh = 1;
		dh = (level < MAX_LEVEL_3D) ? uint32_t(pow(2.0,double(MAX_LEVEL_3D - maxdepth))) : getSize();
		dh2 = getSize();
		iface1 = global3D.nodeface[inode][0];
		iface2 = global3D.nodeface[inode][1];
		iface3 = global3D.nodeface[inode][2];

		if (info[iface1] || info[iface2] || info[iface3]){
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
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh*cz);
			}
			break;
			case 1 :
			{
				cx = 1;
				cy = -1;
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh2*cx,this->y+dh*cy,this->z+dh*cz);
			}
			break;
			case 2 :
			{
				cx = -1;
				cy = 1;
				cz = -1;
				Morton = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2*cy,this->z+dh*cz);
			}
			break;
			case 3 :
			{
				cx = 1;
				cy = 1;
				cz = -1;
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

	// ------------------------------------------------------------------------------- //

	uint64_t 		computeNodeVirtualMorton(uint8_t inode, 		// Computes Morton index (without level) of possible (virtual) neighbours of octant throught iface
			const uint8_t & maxdepth,	// Checks if balanced or not and uses half-size or min-size method (sizeneigh=0 if boundary octant)
			uint32_t & sizeneigh){
//		if (getNotBalance()){
			return computeNodeMinSizeMorton(inode,
					maxdepth,
					sizeneigh);
//		}
//		else{
//			return computeNodeHalfSizeMorton(inode,
//					sizeneigh);
//		}
	};

	// ------------------------------------------------------------------------------- //


};//end Class_Octant<3>

/*  @} */


