// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classParaTree.hpp"

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //
/*! Default Constructor of Para_Tree.
 * It builds one octant with node 0 in the Origin (0,0,0)
 * and side of length 1
 * \param[in] dim_ The space dimension of the octree. PABLO.log is the default value
 * \param[in] logfile The file name for the log of this object. PABLO.log is the default value
 * \param[in] logfile The file name for the log of this object. PABLO.log is the default value
 *
 */
#if NOMPI==0
classParaTree::classParaTree(uint8_t dim_, int8_t maxlevel, string logfile, MPI_Comm comm_) : dim(uint8_t(min(max(2,int(dim_)),3))),log(logfile,comm_),comm(comm_),trans(maxlevel,dim_),octree(maxlevel,dim_){
#else
classParaTree::classParaTree(uint8_t dim_, int8_t maxlevel, string logfile ) : dim(uint8_t(min(max(2,int(dim_)),3))),log(logfile),trans(maxlevel, dim_),octree(maxlevel,dim_){
#endif
//	dim = dim_;
	global.setGlobal(maxlevel, dim);
	serial = true;
	error_flag = 0;
	max_depth = 0;
	global_num_octants = octree.getNumOctants();
#if NOMPI==0
	error_flag = MPI_Comm_size(comm,&nproc);
	error_flag = MPI_Comm_rank(comm,&rank);
#else
	rank = 0;
	nproc = 1;
#endif
	partition_first_desc = new uint64_t[nproc];
	partition_last_desc = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];
	uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
	uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
	for(int p = 0; p < nproc; ++p){
		partition_range_globalidx[p] = 0;
		partition_last_desc[p] = lastDescMorton;
		partition_last_desc[p] = firstDescMorton;
	}
	// Write info log
	log.writeLog("---------------------------------------------");
	log.writeLog("- PABLO PArallel Balanced Linear Octree -");
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
	log.writeLog("---------------------------------------------");
	log.writeLog(" Number of proc		:	" + to_string(static_cast<unsigned long long>(nproc)));
	log.writeLog(" Dimension		:	" + to_string(static_cast<unsigned long long>(dim)));
	log.writeLog(" Max allowed level	:	" + to_string(static_cast<unsigned long long>(global.MAX_LEVEL)));
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
#if NOMPI==0
	MPI_Barrier(comm);
#endif
};

// =============================================================================== //

/*! Constructor of Para_Tree with input parameters.
 * It builds one octant with :
 * \param[in] X Coordinate X of node 0,
 * \param[in] Y Coordinate Y of node 0,
 * \param[in] Z Coordinate Z of node 0,
 * \param[in] L Side length of the octant.
 * \param[in] logfile The file name for the log of this object. PABLO.log is the default value
 */
#if NOMPI==0
classParaTree::classParaTree(double X, double Y, double Z, double L, uint8_t dim_, int8_t maxlevel, string logfile, MPI_Comm comm_):dim(uint8_t(min(max(2,int(dim_)),3))),trans(X,Y,Z,L,maxlevel,dim_),log(logfile,comm_),comm(comm_),octree(maxlevel,dim_){
#else
classParaTree::classParaTree(double X, double Y, double Z, double L, uint8_t dim_, int8_t maxlevel, string logfile):dim(uint8_t(min(max(2,int(dim_)),3))),trans(X,Y,Z,L,maxlevel,dim_),log(logfile),octree(maxlevel,dim_){
#endif
	global.setGlobal(maxlevel, dim);
	serial = true;
	error_flag = 0;
	max_depth = 0;
	global_num_octants = octree.getNumOctants();
#if NOMPI==0
	error_flag = MPI_Comm_size(comm,&nproc);
	error_flag = MPI_Comm_rank(comm,&rank);
#else
	rank = 0;
	nproc = 1;
#endif
	partition_first_desc = new uint64_t[nproc];
	partition_last_desc = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];
	uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
	uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
	for(int p = 0; p < nproc; ++p){
		partition_range_globalidx[p] = 0;
		partition_last_desc[p] = lastDescMorton;
		partition_last_desc[p] = firstDescMorton;
	}
	// Write info log
	log.writeLog("---------------------------------------------");
	log.writeLog("- PABLO PArallel Balanced Linear Octree -");
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
	log.writeLog("---------------------------------------------");
	log.writeLog(" Number of proc		:	" + to_string(static_cast<unsigned long long>(nproc)));
	log.writeLog(" Dimension		:	" + to_string(static_cast<unsigned long long>(dim)));
	log.writeLog(" Max allowed level	:	" + to_string(static_cast<unsigned long long>(global.MAX_LEVEL)));
	log.writeLog(" Domain Origin		:	" + to_string(static_cast<unsigned long long>(X)));
	log.writeLog("				" + to_string(static_cast<unsigned long long>(Y)));
	log.writeLog("				" + to_string(static_cast<unsigned long long>(Z)));
	log.writeLog(" Domain Size		:	" + to_string(static_cast<unsigned long long>(L)));
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
#if NOMPI==0
	MPI_Barrier(comm);
#endif
};

// =============================================================================== //

/*! Constructor of Para_Tree for restart a simulation with input parameters.
 * For each process it builds a vector of octants. The input parameters are :
 * \param[in] X Physical Coordinate X of node 0,
 * \param[in] Y Physical Coordinate Y of node 0,
 * \param[in] Z Physical Coordinate Z of node 0,
 * \param[in] L Physical Side length of the domain,
 * \param[in] XY Coordinates of octants (node 0) in logical domain,
 * \param[in] levels Level of each octant.
 * \param[in] logfile The file name for the log of this object. PABLO.log is the default value
 */
#if NOMPI==0
classParaTree::classParaTree(double X, double Y, double Z, double L, u32vector2D & XYZ, u8vector & levels, uint8_t dim_, int8_t maxlevel, string logfile, MPI_Comm comm_):dim(uint8_t(min(max(2,int(dim_)),3))),trans(X,Y,Z,L,maxlevel,dim_),log(logfile,comm_),comm(comm_),octree(maxlevel,dim_){
#else
classParaTree::classParaTree(double X, double Y, double Z, double L, u32vector2D & XYZ, u8vector & levels, uint8_t dim_, int8_t maxlevel, string logfile):dim(uint8_t(min(max(2,int(dim_)),3))),trans(X,Y,Z,L,maxlevel,dim_),log(logfile),octree(maxlevel,dim_){
#endif
	uint8_t lev, iface;
	uint32_t x0, y0, z0;
	uint32_t NumOctants = XYZ.size();
	dim = dim_;
	global.setGlobal(maxlevel, dim);
	octree.octants.resize(NumOctants);
	for (uint32_t i=0; i<NumOctants; i++){
		lev = uint8_t(levels[i]);
		x0 = uint32_t(XYZ[i][0]);
        y0 = uint32_t(XYZ[i][1]);
        z0 = uint32_t(XYZ[i][2]);
		classOctant oct(false, dim, lev, x0, y0, z0);
		oct.setBalance(false);
		if (x0 == 0){
			iface = 0;
			oct.setBound(iface);
		}
		else if (x0 == global.max_length - oct.getSize(global.MAX_LEVEL)){
			iface = 1;
			oct.setBound(iface);
		}
        if (y0 == 0){
            iface = 2;
            oct.setBound(iface);
        }
        else if (y0 == global.max_length - oct.getSize(global.MAX_LEVEL)){
            iface = 3;
            oct.setBound(iface);
        }
        if (z0 == 0){
            iface = 4;
            oct.setBound(iface);
        }
        else if (z0 == global.max_length - oct.getSize(global.MAX_LEVEL)){
            iface = 5;
            oct.setBound(iface);
        }
		octree.octants[i] = oct;
	}

	//ATTENTO if nompi deve aver l'else
#if NOMPI==0
	error_flag = MPI_Comm_size(comm,&nproc);
	error_flag = MPI_Comm_rank(comm,&rank);
	serial = true;
	if (nproc > 1 ) serial = false;
#else
	serial = true;
	nproc = 1;
	rank = 0;
#endif
	partition_first_desc = new uint64_t[nproc];
	partition_last_desc = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];

	setFirstDesc();
	setLastDesc();
	octree.updateLocalMaxDepth();
	updateAdapt();
#if NOMPI==0
	setPboundGhosts();
#endif
	// Write info log
	log.writeLog("---------------------------------------------");
	log.writeLog("- PABLO PArallel Balanced Linear Octree -");
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
	log.writeLog("---------------------------------------------");
	log.writeLog("- PABLO restart -");
	log.writeLog("---------------------------------------------");
	log.writeLog(" Number of proc		:	" + to_string(static_cast<unsigned long long>(nproc)));
	log.writeLog(" Dimension		:	" + to_string(static_cast<unsigned long long>(dim)));
	log.writeLog(" Max allowed level	:	" + to_string(static_cast<unsigned long long>(global.MAX_LEVEL)));
	log.writeLog(" Domain Origin		:	" + to_string(static_cast<unsigned long long>(X)));
	log.writeLog("				" + to_string(static_cast<unsigned long long>(Y)));
	log.writeLog("				" + to_string(static_cast<unsigned long long>(Z)));
	log.writeLog(" Domain Size		:	" + to_string(static_cast<unsigned long long>(L)));
	log.writeLog(" Number of octants	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
	log.writeLog("---------------------------------------------");
	log.writeLog(" ");
#if NOMPI==0
	MPI_Barrier(comm);
#endif
};

// =============================================================================== //

classParaTree::~classParaTree(){
	log.writeLog("---------------------------------------------");
	log.writeLog("--------------- R.I.P. PABLO ----------------");
	log.writeLog("---------------------------------------------");
	log.writeLog("---------------------------------------------");
};


// =================================================================================== //
// METHODS
// =================================================================================== //

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //

int
classParaTree::getRank(){
	return rank;
};

int
classParaTree::getNproc(){
	return nproc;
};

int
classParaTree::getMaxLevel(){
	return global.MAX_LEVEL;
};

uint32_t
classParaTree::getMaxLength()  {
	return global.max_length;
}


uint8_t
classParaTree::getNnodes()  {
	return global.nnodes;
}

uint8_t
classParaTree::getNfaces()  {
	return global.nfaces;
}

uint8_t
classParaTree::getNedges()  {
	return global.nedges;
}

uint8_t
classParaTree::getNchildren()  {
	return global.nchildren;
}


uint8_t
classParaTree::getNnodesperface()  {
	return global.nnodesperface;
}

void
classParaTree::getNormals(int8_t normals_[6][3])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			normals_[i][j] = global.normals[i][j];
		}
	}
}

void
classParaTree::getOppface(uint8_t oppface_[4])  {
	for (int j=0; j<4; j++){
		oppface_[j] = global.oppface[j];
	}
}

void
classParaTree::getFacenode(uint8_t facenode_[6][3])  {
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			facenode_[i][j] = global.facenode[i][j];
		}
	}
}

void
classParaTree::getNodeface(uint8_t nodeface_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodeface_[i][j] = global.nodeface[i][j];
		}
	}
}

void
classParaTree::getEdgeface(uint8_t edgeface_[12][2])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<2; j++){
			edgeface_[i][j] = global.edgeface[i][j];
		}
	}
}

void
classParaTree::getNodecoeffs(int8_t nodecoeffs_[8][3])  {
	for (int i=0; i<8; i++){
		for (int j=0; j<3; j++){
			nodecoeffs_[i][j] = global.nodecoeffs[i][j];
		}
	}
}

void
classParaTree::getEdgecoeffs(int8_t edgecoeffs_[12][3])  {
	for (int i=0; i<12; i++){
		for (int j=0; j<3; j++){
		edgecoeffs_[i][j] = global.edgecoeffs[i][j];
		}
	}
}

void
classParaTree::setMaxLevel(int8_t maxlevel){
	global.MAX_LEVEL = maxlevel;
};

// =================================================================================== //
// INDEX BASED METHODS
// =================================================================================== //

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate X of node 0.
 */
double
classParaTree::getX(uint32_t idx) {
	return trans.mapX(octree.octants[idx].getX());
}

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate Y of node 0.
 */
double
classParaTree::getY(uint32_t idx) {
	return trans.mapY(octree.octants[idx].getY());
}

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate Z of node 0.
 */
double
classParaTree::getZ(uint32_t idx) {
	return trans.mapZ(octree.octants[idx].getZ());
}

/*! Get the size of an octant, i.e. the side length.
 * \param[in] idx Local index of target octant.
 * \return Area of octant.
 */
double
classParaTree::getSize(uint32_t idx) {
	return trans.mapSize(octree.octants[idx].getSize(global.MAX_LEVEL));
}

/*! Get the area of an octant (for 2D case the same value of getSize).
 * \param[in] idx Local index of target octant.
 * \return Area of octant.
 */
double
classParaTree::getArea(uint32_t idx) {
	return trans.mapSize(octree.octants[idx].getArea(global.MAX_LEVEL));
}

/*! Get the volume of an octant.
 * \param[in] idx Local index of target octant.
 * \return Volume of octant.
 */
double
classParaTree::getVolume(uint32_t idx) {
	return trans.mapArea(octree.octants[idx].getVolume(global.MAX_LEVEL));
}

/*! Get the coordinates of the center of an octant.
 * \param[in] idx Local index of target octant.
 * \param[out] center Coordinates of the center of octant.
 */
void
classParaTree::getCenter(uint32_t idx,
		dvector& center) {
	dvector center_ = octree.octants[idx].getCenter(global.MAX_LEVEL);
	trans.mapCenter(center_, center);
}

/*! Get the coordinates of the center of an octant.
 * \param[in] idx Local index of target octant.
 * \return center Coordinates of the center of octant.
 */
dvector
classParaTree::getCenter(uint32_t idx) {
	dvector center;
	dvector center_ = octree.octants[idx].getCenter(global.MAX_LEVEL);
	trans.mapCenter(center_, center);
	return center;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \return center Coordinates of the center of the iface-th face af octant.
 */
dvector
classParaTree::getFaceCenter(uint32_t idx, uint8_t iface) {
	dvector center;
	dvector center_ = octree.octants[idx].getFaceCenter(iface, global.MAX_LEVEL);
	trans.mapCenter(center_, center);
	return center;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \param[out] center Coordinates of the center of the iface-th face af octant.
 */
void
classParaTree::getFaceCenter(uint32_t idx, uint8_t iface, dvector& center) {
	dvector center_ = octree.octants[idx].getFaceCenter(iface, global.MAX_LEVEL);
	trans.mapCenter(center_, center);
}

/*! Get the coordinates of single node of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] inode Index of the target node.
 * \return center Coordinates of the center of the iface-th face af octant.
 */
dvector
classParaTree::getNode(uint32_t idx, uint8_t inode) {
	dvector node;
	u32vector node_ = octree.octants[idx].getNode(inode, global.MAX_LEVEL);
	trans.mapNode(node_, node);
	return node;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \param[out] center Coordinates of the center of the iface-th face af octant.
 */
void
classParaTree::getNode(uint32_t idx, uint8_t inode, dvector& node) {
	u32vector node_ = octree.octants[idx].getNode(inode, global.MAX_LEVEL);
	trans.mapNode(node_, node);
}

/*! Get the coordinates of the nodes of an octant.
 * \param[in] idx Local index of target octant.
 * \param[out] nodes Coordinates of the nodes of octant.
 */
void
classParaTree::getNodes(uint32_t idx,
		dvector2D & nodes) {
	u32vector2D nodes_;
	octree.octants[idx].getNodes(nodes_, global.MAX_LEVEL);
	trans.mapNodes(nodes_, nodes);
}

/*! Get the coordinates of the nodes of an octant.
 * \param[in] idx Local index of target octant.
 * \return nodes Coordinates of the nodes of octant.
 */
dvector2D
classParaTree::getNodes(uint32_t idx){
	dvector2D nodes;
	u32vector2D nodes_;
	octree.octants[idx].getNodes(nodes_, global.MAX_LEVEL);
	trans.mapNodes(nodes_, nodes);
	return nodes;
}

/*! Get the normal of a face of an octant.
 * \param[in] Local index of target octant.
 * \param[in] iface Index of the face for normal computing.
 * \param[out] normal Coordinates of the normal of face.
 */
void
classParaTree::getNormal(uint32_t idx,
		uint8_t & iface,
		dvector & normal) {
	vector<int8_t> normal_;
	octree.octants[idx].getNormal(iface, normal_, global.normals, global.MAX_LEVEL);
	trans.mapNormals(normal_, normal);
}

/*! Get the normal of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the face for normal computing.
 * \return normal Coordinates of the normal of face.
 */
dvector
classParaTree::getNormal(uint32_t idx,
		uint8_t & iface){
	dvector normal;
	vector<int8_t> normal_;
	octree.octants[idx].getNormal(iface, normal_, global.normals, global.MAX_LEVEL);
	trans.mapNormals(normal_, normal);
	return normal;
}

/*! Get the refinement marker of an octant.
 * \param[in] idx Local index of target octant.
 * \return Marker of octant.
 */
int8_t
classParaTree::getMarker(uint32_t idx){
	return octree.getMarker(idx);
};

/*! Get the level of an octant.
 * \param[in] idx Local index of target octant.
 * \return Level of octant.
 */
uint8_t
classParaTree::getLevel(uint32_t idx){
	return octree.getLevel(idx);
};

/*! Get the balancing condition of an octant.
 * \param[in] idx Local index of target octant.
 * \return Has octant to be balanced?
 */
bool
classParaTree::getBalance(uint32_t idx){
	return !octree.getBalance(idx);
};

#if NOMPI==0
/*! Get the nature of an octant.
 * \param[in] idx Local index of target octant.
 * \return Is octant ghost?
 */
bool
classParaTree::getIsGhost(uint32_t idx){
	return (findOwner(octree.octants[idx].computeMorton()) != rank);
};
#endif

/*! Get if the octant is new after refinement.
 * \param[in] idx Local index of target octant.
 * \return Is octant new?
 */
bool
classParaTree::getIsNewR(uint32_t idx){
	return octree.octants[idx].getIsNewR();
};

/*! Get if the octant is new after coarsening.
 * \param[in] idx Local index of target octant.
 * \return Is octant new?
 */
bool
classParaTree::getIsNewC(uint32_t idx){
	return octree.octants[idx].getIsNewC();
};

/*! Get the global index of an octant.
 * \param[in] idx Local index of target octant.
 * \return Global index of octant.
 */
uint64_t
classParaTree::getGlobalIdx(uint32_t idx){
	if (rank){
		return partition_range_globalidx[rank-1] + uint64_t(idx + 1);
	}
	else{
		return uint64_t(idx);
	};
	return global_num_octants;
};

/*! Get the global index of a ghost octant.
 * \param[in] idx Local index of target ghost octant.
 * \return Global index of ghost octant.
 */
uint64_t
classParaTree::getGhostGlobalIdx(uint32_t idx){
	if (idx<octree.size_ghosts){
		return octree.globalidx_ghosts[idx];
	};
	return uint64_t(octree.size_ghosts);
};


/*! Set the refinement marker of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
 */
void
classParaTree::setMarker(uint32_t idx, int8_t marker){
	octree.setMarker(idx, marker);
};

/*! Set the balancing condition of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
 */
void
classParaTree::setBalance(uint32_t idx, bool balance){
	octree.setBalance(idx, !balance);
};

// =================================================================================== //
// POINTER BASED METHODS
// =================================================================================== //

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate X of node 0.
 */
double
classParaTree::getX(classOctant* oct) {
	return trans.mapX(oct->getX());
}

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate Y of node 0.
 */
double
classParaTree::getY(classOctant* oct) {
	return trans.mapY(oct->getY());
}

/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
 * \param[in] idx Local index of target octant.
 * \return Coordinate Z of node 0.
 */
double
classParaTree::getZ(classOctant* oct) {
	return trans.mapZ(oct->getZ());
}

/*! Get the size of an octant, i.e. the side length.
 * \param[in] idx Local index of target octant.
 * \return Area of octant.
 */
double
classParaTree::getSize(classOctant* oct) {
	return trans.mapSize(oct->getSize(global.MAX_LEVEL));
}

/*! Get the area of an octant (for 2D case the same value of getSize).
 * \param[in] idx Local index of target octant.
 * \return Area of octant.
 */
double
classParaTree::getArea(classOctant* oct) {
	return trans.mapSize(oct->getArea(global.MAX_LEVEL));
}

/*! Get the volume of an octant.
 * \param[in] idx Local index of target octant.
 * \return Volume of octant.
 */
double
classParaTree::getVolume(classOctant* oct) {
	return trans.mapArea(oct->getVolume(global.MAX_LEVEL));
}

/*! Get the coordinates of the center of an octant.
 * \param[in] idx Local index of target octant.
 * \param[out] center Coordinates of the center of octant.
 */
void
classParaTree::getCenter(classOctant* oct,
		dvector& center) {
	dvector center_ = oct->getCenter(global.MAX_LEVEL);
	trans.mapCenter(center_, center);
}

/*! Get the coordinates of the center of an octant.
 * \param[in] idx Local index of target octant.
 * \return center Coordinates of the center of octant.
 */
dvector
classParaTree::getCenter(classOctant* oct) {
	dvector center;
	dvector center_ = oct->getCenter(global.MAX_LEVEL);
	trans.mapCenter(center_, center);
	return center;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \return center Coordinates of the center of the iface-th face af octant.
 */
dvector
classParaTree::getFaceCenter(classOctant* oct, uint8_t iface) {
	dvector center;
	dvector center_ = oct->getFaceCenter(iface, global.MAX_LEVEL);
	trans.mapCenter(center_, center);
	return center;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \param[out] center Coordinates of the center of the iface-th face af octant.
 */
void
classParaTree::getFaceCenter(classOctant* oct, uint8_t iface, dvector& center) {
	dvector center_ = oct->getFaceCenter(iface, global.MAX_LEVEL);
	trans.mapCenter(center_, center);
}

/*! Get the coordinates of single node of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] inode Index of the target node.
 * \return center Coordinates of the center of the iface-th face af octant.
 */
dvector
classParaTree::getNode(classOctant* oct, uint8_t inode) {
	dvector node;
	u32vector node_ = oct->getNode(inode, global.MAX_LEVEL);
	trans.mapNode(node_, node);
	return node;
}

/*! Get the coordinates of the center of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the target face.
 * \param[out] center Coordinates of the center of the iface-th face af octant.
 */
void
classParaTree::getNode(classOctant* oct, uint8_t inode, dvector& node) {
	u32vector node_ = oct->getNode(inode, global.MAX_LEVEL);
	trans.mapNode(node_, node);
}

/*! Get the coordinates of the nodes of an octant.
 * \param[in] idx Local index of target octant.
 * \param[out] nodes Coordinates of the nodes of octant.
 */
void
classParaTree::getNodes(classOctant* oct,
		dvector2D & nodes) {
	u32vector2D nodes_;
	oct->getNodes(nodes_, global.MAX_LEVEL);
	trans.mapNodes(nodes_, nodes);
}

/*! Get the coordinates of the nodes of an octant.
 * \param[in] idx Local index of target octant.
 * \return nodes Coordinates of the nodes of octant.
 */
dvector2D
classParaTree::getNodes(classOctant* oct){
	dvector2D nodes;
	u32vector2D nodes_;
	oct->getNodes(nodes_, global.MAX_LEVEL);
	trans.mapNodes(nodes_, nodes);
	return nodes;
}

/*! Get the normal of a face of an octant.
 * \param[in] Local index of target octant.
 * \param[in] iface Index of the face for normal computing.
 * \param[out] normal Coordinates of the normal of face.
 */
void
classParaTree::getNormal(classOctant* oct,
		uint8_t & iface,
		dvector & normal) {
	vector<int8_t> normal_;
	oct->getNormal(iface, normal_, global.normals, global.MAX_LEVEL);
	trans.mapNormals(normal_, normal);
}

/*! Get the normal of a face of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] iface Index of the face for normal computing.
 * \return normal Coordinates of the normal of face.
 */
dvector
classParaTree::getNormal(classOctant* oct,
		uint8_t & iface){
	dvector normal;
	vector<int8_t> normal_;
	oct->getNormal(iface, normal_, global.normals, global.MAX_LEVEL);
	trans.mapNormals(normal_, normal);
	return normal;
}

/*! Get the refinement marker of an octant.
 * \param[in] idx Local index of target octant.
 * \return Marker of octant.
 */
int8_t
classParaTree::getMarker(classOctant* oct){
	return oct->getMarker();
};

/*! Get the level of an octant.
 * \param[in] idx Local index of target octant.
 * \return Level of octant.
 */
uint8_t
classParaTree::getLevel(classOctant* oct){
	return oct->getLevel();
};

/*! Get the balancing condition of an octant.
 * \param[in] idx Local index of target octant.
 * \return Has octant to be balanced?
 */
bool
classParaTree::getBalance(classOctant* oct){
	return !oct->getBalance();
};

/*! Get if the octant is new after refinement.
 * \param[in] idx Local index of target octant.
 * \return Is octant new?
 */
bool
classParaTree::getIsNewR(classOctant* oct){
	return oct->getIsNewR();
};

/*! Get if the octant is new after coarsening.
 * \param[in] idx Local index of target octant.
 * \return Is octant new?
 */
bool
classParaTree::getIsNewC(classOctant* oct){
	return oct->getIsNewC();
};

/*! Set the refinement marker of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
 */
void
classParaTree::setMarker(classOctant* oct, int8_t marker){
	oct->setMarker(marker);
};

/*! Set the balancing condition of an octant.
 * \param[in] idx Local index of target octant.
 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
 */
void
classParaTree::setBalance(classOctant* oct, bool balance){
	oct->setBalance(!balance);
};

// =================================================================================== //
// LOCAL TREE GET/SET METHODS
// =================================================================================== //

/*! Get the status label of the octree.
 * 	\return Status.
 */
uint64_t
classParaTree::getStatus(){
	return status;
}

/*! Get the local number of octants.
 * \return Local number of octants.
 */
uint32_t
classParaTree::getNumOctants() const{
	return octree.getNumOctants();
};

/*! Get the local number of ghost octants.
 * \return Local number of ghost octants.
 */
uint32_t
classParaTree::getNumGhosts() const{
	return octree.getSizeGhost();
};

/** Get the local number of nodes.
 */
uint32_t
classParaTree::getNumNodes() const{
	return octree.nodes.size();
}

/*! Get the local depth of octree.
 * \return Local depth of octree.
 */
uint8_t
classParaTree::getLocalMaxDepth() const{
	return octree.getLocalMaxDepth();
};

/*! Get the codimension for 2:1 balancing
 * \return Maximum codimension of the entity through which the 2:1 balance is performed.
 */
uint8_t
classParaTree::getBalanceCodimension() const{
	return octree.getBalanceCodim();
};

/*! Get the coordinates of the extreme points of a bounding box containing the local tree
 *  \param[out] P0 Vector with coordinates of the first point (lowest coordinates);
 *  \param[out] P1 Vector with coordinates of the last point (highest coordinates).
 */
void
classParaTree::getBoundingBox(dvector & P0, dvector & P1){
	dvector	cnode;
	uint32_t 	nocts = getNumOctants();
	uint32_t id = 0;
	P0 = getNode(id, 0);
	id = nocts-1;
	P1 = getNode(id, global.nnodes-1);
	for (uint32_t idx=0; idx<global.nnodes; idx++){
		for (uint8_t inode=0; inode<nocts; inode++){
			cnode = getNode(idx, inode);
			for (uint8_t i=0; i<3; i++){
				P0[i] = min(P0[i], cnode[i]);
				P1[i] = max(P1[i], cnode[i]);
			}
		}
	}
};

/*! Get the coordinates of the extreme points of a bounding box containing the local tree
 *  \param[out] P0 Array with coordinates of the first point (lowest coordinates);
 *  \param[out] P1 Array with coordinates of the last point (highest coordinates).
 */
void
classParaTree::getBoundingBox(darray3 & P0, darray3 & P1){
	dvector	cnode, cnode0, cnode1;
	uint32_t 	nocts = getNumOctants();
	uint32_t	id = 0;
	cnode0 = getNode(id, 0);
	id = nocts-1;
	cnode1 = getNode(id, global.nnodes-1);
	for (uint8_t i=0; i<3; i++){
		P0[i] = cnode0[i];
		P1[i] = cnode1[i];
	}
	for (uint32_t idx=0; idx<global.nnodes; idx++){
		for (uint8_t inode=0; inode<nocts; inode++){
			cnode = getNode(idx, inode);
			for (uint8_t i=0; i<3; i++){
				P0[i] = min(P0[i], cnode[i]);
				P1[i] = max(P1[i], cnode[i]);
			}
		}
	}
};

const
classOctant & classParaTree::getFirstDesc() const{
	return octree.getFirstDesc();
};

const
classOctant & classParaTree::getLastDesc() const{
	return octree.getLastDesc();
};

uint64_t
classParaTree::getLastDescMorton(uint32_t idx) {
	return octree.octants[idx].buildLastDesc(global.MAX_LEVEL).computeMorton();
};

/*! Set the codimension for 2:1 balancing
 * \param[in] Maximum codimension of the entity through which the 2:1 balance is performed (1 = 2:1 balance through edges (default); 2 = 2:1 balance through nodes and edges).
 */
void
classParaTree::setBalanceCodimension(uint8_t b21codim){
	octree.setBalanceCodim(b21codim);
};

// =================================================================================== //
// INTERSECTION GET/SET METHODS
// =================================================================================== //

/*! Get the local number of intersections.
 * \return Local number of intersections.
 */
uint32_t
classParaTree::getNumIntersections() {
	return octree.intersections.size();
}

/*! Get a pointer to target intersection.
 * \param[in] idx Local index of intersection.
 * \return Pointer to target intersection.
 */
classIntersection*
classParaTree::getIntersection(uint32_t idx) {
	if (idx < octree.intersections.size()){
		return &octree.intersections[idx];
	}
	return NULL;
}

/*! Get the level of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \return Level of intersection.
 */
uint8_t
classParaTree::getLevel(classIntersection* inter) {
	if(inter->finer && inter->isghost)
		return octree.extractGhostOctant(inter->owners[inter->finer]).getLevel();
	else
		return octree.extractOctant(inter->owners[inter->finer]).getLevel();
}

/*! Get the finer owner octant of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \return The finer octant of the owners of intersection (false/true = 0/1).
 */
bool
classParaTree::getFiner(classIntersection* inter) {
	return inter->finer;
}

/*! Get if an intersection is a boundary domain intersection.
 * \param[in] inter Pointer to target intersection.
 * \return Boundary or not boundary?.
 */
bool
classParaTree::getBound(classIntersection* inter) {
	return inter->getBound();
}

/*! Get if an intersection is an intersection between an internal and a ghost element.
 * \param[in] inter Pointer to target intersection.
 * \return Ghost or not ghost?.
 */
bool
classParaTree::getIsGhost(classIntersection* inter) {
	return inter->getIsGhost();
}

/*! Get if an intersection is a boundary intersection for a process.
 * \param[in] inter Pointer to target intersection.
 * \return Process boundary or not boundary?.
 */
bool
classParaTree::getPbound(classIntersection* inter) {
	return inter->getPbound();
}

/*! Get the face index of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \return Face index of the first octant owner of intersection (owners[0]).
 */
uint8_t
classParaTree::getFace(classIntersection* inter) {
	return inter->iface;
}

/*! Get the owner octants of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \return A couple of octants owners of intersection.
 */
u32vector
classParaTree::getOwners(classIntersection* inter) {
	u32vector owners(2);
	owners[0] = inter->owners[0];
	owners[1] = inter->owners[1];
	return owners;
}

/*! Get the owner octant of an intersection with inner normal.
 * \param[in] inter Pointer to target intersection.
 * \return Index of the octant owner with inner normal.
 */
uint32_t
classParaTree::getIn(classIntersection* inter) {
	return inter->getIn();
}

/*! Get the owner octant of an intersection with outer normal.
 * \param[in] inter Pointer to target intersection.
 * \return Index of the octant owner with outer normal.
 */
uint32_t
classParaTree::getOut(classIntersection* inter) {
	return inter->getOut();
}

/*! Get the size of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \return Size of intersection.
 */
double
classParaTree::getSize(classIntersection* inter) {
	uint32_t Size;
	if(inter->finer && inter->isghost)
		Size = octree.extractGhostOctant(inter->owners[inter->finer]).getSize(global.MAX_LEVEL);
	else
		Size = octree.extractOctant(inter->owners[inter->finer]).getSize(global.MAX_LEVEL);
	return trans.mapSize(Size);
}

/*! Get the area of an intersection (for 2D case the same value of getSize).
 * \param[in] inter Pointer to target intersection.
 * \return Area of intersection.
 */
double
classParaTree::getArea(classIntersection* inter) {
	uint32_t Area;
	if(inter->finer && inter->isghost)
		Area = octree.extractGhostOctant(inter->owners[1]).getArea(global.MAX_LEVEL);
	else
		Area = octree.extractOctant(inter->owners[inter->finer]).getArea(global.MAX_LEVEL);
	return trans.mapSize(Area);
}

/*! Get the coordinates of the center of an intersection.
 * \param[in] inter Pointer to target intersection.
 * \param[out] center Coordinates of the center of intersection.
 */
vector<double>
classParaTree::getCenter(classIntersection* inter){
	vector<double> center;
	classOctant oct;
	if(inter->finer && inter->isghost)
		oct = octree.extractGhostOctant(inter->owners[inter->finer]);
	else
		oct = octree.extractOctant(inter->owners[inter->finer]);
	vector<double>  center_ = oct.getCenter(global.MAX_LEVEL);
	int sign = ( int(2*((inter->iface)%2)) - 1);
	double deplace = double (sign * int(oct.getSize(global.MAX_LEVEL))) / 2;
	center_[inter->iface/2] = uint32_t(int(center_[inter->iface/2]) + deplace);
	trans.mapCenter(center_, center);
	return center;
}

/*! Get the coordinates of the nodes of an intersection.
 * \param[in] oct Pointer to target intersection.
 * \return nodes Coordinates of the nodes of intersection.
 */
dvector2D
classParaTree::getNodes(classIntersection* inter){
	dvector2D nodes;
	classOctant oct;
	if(inter->finer && inter->isghost)
		oct = octree.extractGhostOctant(inter->owners[inter->finer]);
	else
		oct = octree.extractOctant(inter->owners[inter->finer]);
	uint8_t iface = inter->iface;
	u32vector2D nodes_all;
	oct.getNodes(nodes_all, global.MAX_LEVEL);
	u32vector2D nodes_(global.nnodesperface, u32vector(3));
	for (int i=0; i<global.nnodesperface; i++){
		for (int j=0; j<3; j++){
			nodes_[i][j] = nodes_all[global.facenode[iface][i]][j];
		}
	}
	trans.mapNodesIntersection(nodes_, nodes);
	return nodes;
}

/*! Get the normal of an intersection.
 * \param[in] oct Pointer to target intersection.
 * \param[out] normal Coordinates of the normal of intersection.
 */
dvector
classParaTree::getNormal(classIntersection* inter){
	dvector normal;
	classOctant oct;
	if(inter->finer && inter->isghost)
		oct = octree.extractGhostOctant(inter->owners[inter->finer]);
	else
		oct = octree.extractOctant(inter->owners[inter->finer]);
	uint8_t iface = inter->iface;
	vector<int8_t> normal_;
	oct.getNormal(iface, normal_, global.normals, global.MAX_LEVEL);
	trans.mapNormals(normal_, normal);
	return normal;
}

// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

/** Get an octant as pointer to the target octant.
 * \param[in] idx Local index of target octant.
 * \return Pointer to target octant.
 */
classOctant*
classParaTree::getOctant(uint32_t idx) {
	if (idx < octree.getNumOctants()){
		return &octree.octants[idx] ;
	}
	return NULL;
};

/** Get a ghost octant as pointer to the target octant.
 * \param[in] idx Local index (in ghosts structure) of target ghost octant.
 * \return Pointer to target ghost octant.
 */
classOctant*
classParaTree::getGhostOctant(uint32_t idx) {
	if (idx < octree.getSizeGhost()){
		return &octree.ghosts[idx] ;
	}
	return NULL;
};

/*! Get the global index of an octant.
 * \param[in] idx Local index of target octant.
 * \return Global index of octant.
 */
uint64_t
classParaTree::getGlobalIdx(classOctant* oct){
#if NOMPI==0
	if (getIsGhost(oct)){
		uint32_t idx = octree.findGhostMorton(oct->computeMorton());
		return octree.globalidx_ghosts[idx];
	}
#endif
	uint32_t idx = octree.findMorton(oct->computeMorton());
	if (rank){
		return partition_range_globalidx[rank-1] + uint64_t(idx + 1);
	}
	return uint64_t(idx);
};

/*! Get the local index of an octant.
 * \param[in] oct Pointer to target octant.
 * \return Local index of octant.
 */
uint32_t
classParaTree::getIdx(classOctant* oct){
#if NOMPI==0
	if (getIsGhost(oct)){
		return octree.findGhostMorton(oct->computeMorton());
	}
#endif
	return octree.findMorton(oct->computeMorton());
};

/*! Get the local index of an octant.
 * \param[in] oct Target octant.
 * \return Local index of octant.
 */
uint32_t
classParaTree::getIdx(classOctant oct){
#if NOMPI==0
	if (getIsGhost(oct)){
		return octree.findGhostMorton(oct.computeMorton());
	}
	else{
#endif
		return octree.findMorton(oct.computeMorton());
#if NOMPI==0
	};
#endif
	return octree.getNumOctants();
};

#if NOMPI==0
/*! Get the nature of an octant.
 * \param[in] oct Pointer to target octant.
 * \return Is octant ghost?
 */
bool
classParaTree::getIsGhost(classOctant* oct){
	if (serial)
		return false;
	return (findOwner(oct->computeMorton()) != rank);
};

/*! Get the nature of an octant.
 * \param[in] oct Target octant.
 * \return Is octant ghost?
 */
bool
classParaTree::getIsGhost(classOctant oct){
	if (serial)
		return false;
	return (findOwner(oct.computeMorton()) != rank);
};
#endif

// =================================================================================== //
// PRIVATE GET/SET METHODS
// =================================================================================== //

void
classParaTree::setFirstDesc(){
	octree.setFirstDesc();
};

void
classParaTree::setLastDesc(){
	octree.setLastDesc();
};

// =================================================================================== //
// OTHER METHODS												    			   //
// =================================================================================== //

// =================================================================================== //
// OTHER OCTANT BASED METHODS												    			   //
// =================================================================================== //

/** Finds neighbours of octant through iface in vector octants.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] idx Index of current octant
 * \param[in] iface Index of face/edge/node passed through for neighbours finding
 * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
 * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
 * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
void
classParaTree::findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim,
		u32vector & neighbours, vector<bool> & isghost){

	bool	Fedge = ((codim>1) && (dim==3));
	bool	Fnode = (codim == dim);

	if (codim == 1){
		octree.findNeighbours(idx, iface, neighbours, isghost);
	}
	else if (Fedge){
		octree.findEdgeNeighbours(idx, iface, neighbours, isghost);
	}
	else if (Fnode){
		octree.findNodeNeighbours(idx, iface, neighbours, isghost);
	}
	else {
		neighbours.clear();
		isghost.clear();
	}
};

/** Finds neighbours of octant through iface in vector octants.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] oct Pointer to current octant
 * \param[in] iface Index of face/edge/node passed through for neighbours finding
 * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
 * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
 * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
void
classParaTree::findNeighbours(classOctant* oct, uint8_t iface, uint8_t codim,
		u32vector & neighbours, vector<bool> & isghost){

	bool	Fedge = ((codim>1) && (dim==3));
	bool	Fnode = (codim == dim);

	if (codim == 1){
		octree.findNeighbours(oct, iface, neighbours, isghost);
	}
	else if (Fedge){
		octree.findEdgeNeighbours(oct, iface, neighbours, isghost);
	}
	else if (Fnode){
		octree.findNodeNeighbours(oct, iface, neighbours, isghost);
	}
	else {
		neighbours.clear();
		isghost.clear();
	}

};

/** Finds neighbours of ghost octant through iface in vector octants.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure ( only local octants ).
 * \param[in] idx Index of current octant
 * \param[in] iface Index of face/edge/node passed through for neighbours finding
 * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
 * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
 */
void
classParaTree::findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours){

	bool	Fedge = ((codim>1) && (dim==3));
	bool	Fnode = (codim == dim);

	if (codim == 1){
		octree.findGhostNeighbours(idx, iface, neighbours);
	}
	else if (Fedge){
		octree.findGhostEdgeNeighbours(idx, iface, neighbours);
	}
	else if (Fnode){
		octree.findGhostNodeNeighbours(idx, iface, neighbours);
	}
	else {
		neighbours.clear();
	}
};

/** Get the octant owner of an input point.
 * \param[in] point Coordinates of target point.
 * \return Pointer to octant owner of target point (=NULL if point is outside of the domain).
 */
classOctant*
classParaTree::getPointOwner(dvector & point){
	uint32_t noctants = octree.octants.size();
	uint32_t idxtry = noctants/2;
	uint32_t x, y, z;
	uint64_t morton, mortontry;
	int powner = 0;

	x = trans.mapX(point[0]);
	y = trans.mapX(point[1]);
	z = trans.mapX(point[2]);
	if ((x > global.max_length) || (y > global.max_length) || (z > global.max_length))
		return NULL;

	if (x == global.max_length) x = x - 1;
	if (y == global.max_length) y = y - 1;
	if (z == global.max_length) z = z - 1;
	morton = mortonEncode_magicbits(x,y,z);

	powner = 0;
	if (!serial) powner = findOwner(morton);

	if ((powner!=rank) && (!serial))
		return NULL;

	int32_t jump = idxtry;
	while(abs(jump) > 0){
		mortontry = octree.octants[idxtry].computeMorton();
		jump = ((mortontry<morton)-(mortontry>morton))*abs(jump)/2;
		idxtry += jump;
		if (idxtry > noctants-1){
			if (jump > 0){
				idxtry = noctants - 1;
				jump = 0;
			}
			else if (jump < 0){
				idxtry = 0;
				jump = 0;
			}
		}
	}
	if(octree.octants[idxtry].computeMorton() == morton){
		return &octree.octants[idxtry];
	}
	else{
		// Step until the mortontry lower than morton (one idx of distance)
		{
			while(octree.octants[idxtry].computeMorton() < morton){
				idxtry++;
				if(idxtry > noctants-1){
					idxtry = noctants-1;
					break;
				}
			}
			while(octree.octants[idxtry].computeMorton() > morton){
				idxtry--;
				if(idxtry > noctants-1){
					idxtry = 0;
					break;
				}
			}
		}
		return &octree.octants[idxtry];
	}

};

/** Get the octant owner of an input point.
 * \param[in] point Coordinates of target point.
 * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
 */
uint32_t
classParaTree::getPointOwnerIdx(dvector & point){
	uint32_t noctants = octree.octants.size();
	uint32_t idxtry = noctants/2;
	uint32_t x, y, z;
	uint64_t morton, mortontry;
	int powner = 0;

	x = trans.mapX(point[0]);
	y = trans.mapY(point[1]);
	z = trans.mapZ(point[2]);

	if ((x > global.max_length) || (y > global.max_length) || (z > global.max_length)
			|| (point[0] < trans.X0) || (point[1] < trans.Y0) || (point[2] < trans.Z0)){
		return -1;
	}

	if (x == global.max_length) x = x - 1;
	if (y == global.max_length) y = y - 1;
	if (z == global.max_length) z = z - 1;
	morton = mortonEncode_magicbits(x,y,z);


	powner = 0;
	if(!serial) powner = findOwner(morton);

	if ((powner!=rank) && (!serial))
		return -1;

	int32_t jump = idxtry;
	while(abs(jump) > 0){

		mortontry = octree.octants[idxtry].computeMorton();
		jump = ((mortontry<morton)-(mortontry>morton))*abs(jump)/2;
		idxtry += jump;
		if (idxtry > noctants-1){
			if (jump > 0){
				idxtry = noctants - 1;
				jump = 0;
			}
			else if (jump < 0){
				idxtry = 0;
				jump = 0;
			}
		}
	}
	if(octree.octants[idxtry].computeMorton() == morton){
		return idxtry;
	}
	else{
		// Step until the mortontry lower than morton (one idx of distance)
		{
			while(octree.octants[idxtry].computeMorton() < morton){
				idxtry++;
				if(idxtry > noctants-1){
					idxtry = noctants-1;
					break;
				}
			}
			while(octree.octants[idxtry].computeMorton() > morton){
				idxtry--;
				if(idxtry > noctants-1){
					idxtry = 0;
					break;
				}
			}
		}
		return idxtry;
	}
};

/** Get mapping info of an octant after an adapting with tracking changes.
 * \param[in] idx Index of new octant.
 * \param[out] mapper Mapper from new octants to old octants. I.e. mapper[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
 * if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
 * if the i-th octant is new after coarsening the j-th old octant was a child of the new octant (mapper size = 4).
 * \param[out] isghost Info on ghostness of old octants.
 * I.e. isghost[i] = true/false -> the mapper[i] = j-th old octant was a local/ghost octant.
 */
void
classParaTree::getMapping(uint32_t & idx, u32vector & mapper, vector<bool> & isghost){

	uint32_t	i, nocts = getNumOctants();
	uint32_t	nghbro = octree.last_ghost_bros.size();;

	mapper.clear();
	isghost.clear();

	mapper.push_back(mapidx[idx]);
	isghost.push_back(false);
	if (getIsNewC(idx)){
		if (idx < nocts-1 || !nghbro){
			for (i=1; i<global.nchildren; i++){
				mapper.push_back(mapidx[idx]+i);
				isghost.push_back(false);
			}
		}
		else if (idx == nocts-1 && nghbro){
			for (i=1; i<global.nchildren-nghbro; i++){
				mapper.push_back(mapidx[idx]+i);
				isghost.push_back(false);
			}
			for (i=0; i<nghbro; i++){
				mapper.push_back(octree.last_ghost_bros[i]);
				isghost.push_back(true);
			}
		}
	}

};

// =================================================================================== //
// OTHER PARATREE BASED METHODS												    			   //
// =================================================================================== //

/** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
 * \param[in] mapper_flag True/False if you want/don't want to track the changes in structure octant by a mapper.
 * \return Boolean if adapt has done something.
 */
bool
classParaTree::adapt(bool mapper_flag){

	bool done = false;

//	if (mapper_flag){
		done = private_adapt_mapidx(mapper_flag);
		status += done;
		return done;
//	}
//	else{
//		done = private_adapt();
//		status += done;
//		return done;
//	}

};

/** Adapt the octree mesh refining all the octants by one level.
 * Optionally track the changes in structure octant by a mapper.
 * \param[in] mapper_flag True/false for tracking/not tracking the changes in structure octant .
 */
bool
classParaTree::adaptGlobalRefine(bool mapper_flag) {
	//TODO recoding for adapting with abs(marker) > 1
	bool globalDone = false, localDone = false;
	uint32_t nocts = octree.getNumOctants();
	vector<classOctant>::iterator iter, iterend = octree.octants.end();

	for (iter = octree.octants.begin(); iter != iterend; iter++){
		iter->info[12] = false;
		iter->info[13] = false;
		iter->info[15] = false;
	}

	mapidx.clear();
	if (mapper_flag){
		// mapidx init
		mapidx.resize(nocts);
		u32vector(mapidx).swap(mapidx);

		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
		}
	}
#if NOMPI==0
	if(serial){
#endif
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Global Refine)");
		log.writeLog(" ");

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));

		// Refine
		if (mapper_flag){
			while(octree.globalRefine(mapidx));
		}
		else{
			while(octree.globalRefine(mapidx));
		}

		if (octree.getNumOctants() > nocts)
			localDone = true;
		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));
		nocts = octree.getNumOctants();
		updateAdapt();

#if NOMPI==0
		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
#endif
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
#if NOMPI==0
	}
	else{
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Global Refine)");
		log.writeLog(" ");

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));

		// Refine
		if (mapper_flag){
			while(octree.globalRefine(mapidx));
		}
		else{
			while(octree.globalRefine(mapidx));
		}

		if (octree.getNumOctants() > nocts)
			localDone = true;
		updateAdapt();
		setPboundGhosts();
		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
		nocts = octree.getNumOctants();

		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
	}
	return globalDone;
#else
	return localDone;
#endif
}

/** Adapt the octree mesh coarsening all the octants by one level.
 * Optionally track the changes in structure octant by a mapper.
 * \param[in] mapper_flag True/false for tracking/not tracking the changes in structure octant .
 */
bool
classParaTree::adaptGlobalCoarse(bool mapper_flag) {
	//TODO recoding for adapting with abs(marker) > 1
	bool globalDone = false, localDone = false;
	uint32_t nocts = octree.getNumOctants();
	vector<classOctant>::iterator iter, iterend = octree.octants.end();

	for (iter = octree.octants.begin(); iter != iterend; iter++){
		iter->info[12] = false;
		iter->info[13] = false;
		iter->info[15] = false;
	}

	mapidx.clear();
	if (mapper_flag){
		// mapidx init
		mapidx.resize(nocts);
		u32vector(mapidx).swap(mapidx);

		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
		}
	}
#if NOMPI==0
	if(serial){
#endif
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Global Coarse)");
		log.writeLog(" ");

		// 2:1 Balance
		balance21(true);

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));

		// Coarse
		if (mapper_flag){
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse(mapidx);
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
		}
		else{
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse();
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
		}

		if (octree.getNumOctants() < nocts){
			localDone = true;
		}
		nocts = octree.getNumOctants();

		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(nocts)));
#if NOMPI==0
		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
#endif
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
#if NOMPI==0
	}
	else{
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Global Coarse)");
		log.writeLog(" ");

		// 2:1 Balance
		balance21(true);

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));

		// Coarse
		if (mapper_flag){
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse(mapidx);
			setPboundGhosts();
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
		}
		else{
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse();
			setPboundGhosts();
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
		}
		setPboundGhosts();
		if (octree.getNumOctants() < nocts){
			localDone = true;
		}
		nocts = octree.getNumOctants();

		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
	}
	return globalDone;
#else
	return localDone;
#endif
}

/** It finds the process owning the element definded by the Morton number passed as argument
 * The Morton number can be computed using the method classOctant#computeMorton().
 * \param[in] morton is the Morton number of the element you want find the owner of
 * \return it returns the rank of the process owning the element
 */
int
classParaTree::findOwner(const uint64_t & morton) {
	int p = -1;
	int length = nproc;
	int beg = 0;
	int end = nproc -1;
	int seed = nproc/2;
	while(beg != end){
		if(morton <= partition_last_desc[seed]){
			end = seed;
			if(morton > partition_last_desc[seed-1])
				beg = seed;
		}
		else{
			beg = seed;
			if(morton <= partition_last_desc[seed+1])
				beg = seed + 1;
		}
		length = end - beg;
		seed = beg + length/2;
	}
	p = beg;
	return p;
}

/** Compute the connectivity of octants and store the coordinates of nodes.
 */
void
classParaTree::computeConnectivity() {
	octree.computeConnectivity();
}

/** Clear the connectivity of octants.
 */
void
classParaTree::clearConnectivity() {
	octree.clearConnectivity();
}

/** Update the connectivity of octants.
 */
void
classParaTree::updateConnectivity() {
	octree.updateConnectivity();
}

/** Get the connectivity of the octants
 * \return connectivity Matrix of noctants*4 with the connectivity of each octant (4 indices of nodes).
 */
const u32vector2D &
classParaTree::getConnectivity(){
	return octree.connectivity;
}

/** Get the local connectivity of an octant
 * \param[in] idx Local index of octant
 * \return connectivity Connectivity of the octant (4 indices of nodes).
 */
const u32vector &
classParaTree::getConnectivity(uint32_t idx){
	return octree.connectivity[idx];
}

/** Get the local connectivity of an octant
 * \param[in] oct Pointer to an octant
 * \return connectivity Connectivity of the octant (4 indices of nodes).
 */
const u32vector &
classParaTree::getConnectivity(classOctant* oct){
	return octree.connectivity[getIdx(oct)];
}

/** Get the logical coordinates of the nodes
 * \return nodes Matrix of nnodes*3 with the coordinates of the nodes.
 */
const u32arr3vector &
classParaTree::getNodes(){
	return octree.nodes;
}

/** Get the logical coordinates of a node
 * \param[in] inode Local index of node
 * \return nodes Vector with the coordinates of the node.
 */
const u32array3 &
classParaTree::getNodeLogicalCoordinates(uint32_t inode){
	return octree.nodes[inode];
}

/** Get the physical coordinates of a node
 * \param[in] inode Local index of node
 * \return nodes Vector with the coordinates of the node.
 */
dvector
classParaTree::getNodeCoordinates(uint32_t inode){
	vector<double> coords(3,0);
	coords[0] = trans.mapX(octree.nodes[inode][0]);
	coords[1] = trans.mapY(octree.nodes[inode][1]);
	coords[2] = trans.mapZ(octree.nodes[inode][2]);
	return coords;
}

/** Compute the connectivity of ghost octants and store the coordinates of nodes.
 */
void
classParaTree::computeGhostsConnectivity() {
	octree.computeGhostsConnectivity();
}

/** Clear the connectivity of ghost octants.
 */
void
classParaTree::clearGhostsConnectivity() {
	octree.clearGhostsConnectivity();
}

/** Update the connectivity of ghost octants.
 */
void
classParaTree::updateGhostsConnectivity() {
	octree.updateGhostsConnectivity();
}

/** Get the connectivity of the ghost octants
 * \return connectivity Matrix of nghostoctants*4 with the connectivity of each octant (4 indices of nodes).
 */
const u32vector2D &
classParaTree::getGhostConnectivity(){
	return octree.ghostsconnectivity;
}

/** Get the local connectivity of a ghost octant
 * \param[in] idx Local index of ghost octant
 * \return connectivity Connectivity of the ghost octant (4 indices of nodes).
 */
const u32vector &
classParaTree::getGhostConnectivity(uint32_t idx){
	return octree.ghostsconnectivity[idx];
}

/** Get the local connectivity of a ghost octant
 * \param[in] oct Pointer to a ghost octant
 * \return connectivity Connectivity of the ghost octant (4 indices of nodes).
 */
const u32vector &
classParaTree::getGhostConnectivity(classOctant* oct){
	return octree.ghostsconnectivity[getIdx(oct)];
}

/** Get the logical coordinates of the ghost nodes
 * \return nodes Matrix of nghostnodes*3 with the coordinates of the nodes.
 */
const u32arr3vector &
classParaTree::getGhostNodes(){
	return octree.ghostsnodes;
}

/** Get the logical coordinates of a ghost node
 * \param[in] inode Local index of node
 * \return nodes Vector with the coordinates of the node.
 */
const u32array3 &
classParaTree::getGhostNodeLogicalCoordinates(uint32_t inode){
	return octree.ghostsnodes[inode];
}

/** Get the physical coordinates of a ghost node
 * \param[in] inode Local index of node
 * \return nodes Vector with the coordinates of the node.
 */
dvector
classParaTree::getGhostNodeCoordinates(uint32_t inode){
	vector<double> coords(3,0);
	coords[0] = trans.mapX(octree.ghostsnodes[inode][0]);
	coords[1] = trans.mapY(octree.ghostsnodes[inode][1]);
	coords[2] = trans.mapZ(octree.ghostsnodes[inode][2]);
	return coords;
}

#if NOMPI==0
/** Distribute Load-Balancing the octants of the whole tree over
 * the processes of the job following the Morton order.
 * Until loadBalance is not called for the first time the mesh is serial.
 */
void
classParaTree::loadBalance(){

	//Write info on log
	log.writeLog("---------------------------------------------");
	log.writeLog(" LOAD BALANCE ");

	uint32_t* partition = new uint32_t [nproc];
	computePartition(partition);
	if(serial)
	{
		log.writeLog(" ");
		log.writeLog(" Initial Serial distribution : ");
		for(int ii=0; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]+1)));
		}

		uint32_t stride = 0;
		for(int i = 0; i < rank; ++i)
			stride += partition[i];
		classLocalTree::octvector octantsCopy = octree.octants;
		classLocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
		classLocalTree::octvector::const_iterator last = first + partition[rank];
		octree.octants.assign(first, last);

		octvector(octree.octants).swap(octree.octants);

		first = octantsCopy.end();
		last = octantsCopy.end();

		//Update and ghosts here
		updateLoadBalance();
		setPboundGhosts();

	}
	else
	{
		log.writeLog(" ");
		log.writeLog(" Initial Parallel partition : ");
		log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[0]+1)));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]-partition_range_globalidx[ii-1])));
		}

		//empty ghosts
		octree.ghosts.clear();
		octree.size_ghosts = 0;
		//compute new partition range globalidx
		uint64_t* newPartitionRangeGlobalidx = new uint64_t[nproc];
		for(int p = 0; p < nproc; ++p){
			newPartitionRangeGlobalidx[p] = 0;
			for(int pp = 0; pp <= p; ++pp)
				newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
			--newPartitionRangeGlobalidx[p];
		}

		//find resident octants local offset lastHead(lh) and firstTail(ft)
		int32_t lh,ft;
		if(rank == 0)
			lh = -1;
		else{
			lh = (int32_t)(newPartitionRangeGlobalidx[rank-1] + 1 - partition_range_globalidx[rank-1] - 1 - 1);
		}
		if(lh < 0)
			lh = - 1;
		else if(lh > (int32_t)(octree.octants.size() - 1))
			lh = octree.octants.size() - 1;

		if(rank == nproc - 1)
			ft = octree.octants.size();
		else if(rank == 0)
			ft = (int32_t)(newPartitionRangeGlobalidx[rank] + 1);
		else{
			ft = (int32_t)(newPartitionRangeGlobalidx[rank] - partition_range_globalidx[rank -1]);
		}
		if(ft > (int32_t)(octree.octants.size() - 1))
			ft = octree.octants.size();
		else if(ft < 0)
			ft = 0;

		//compute size Head and size Tail
		uint32_t headSize = (uint32_t)(lh + 1);
		uint32_t tailSize = (uint32_t)(octree.octants.size() - ft);
		uint32_t headOffset = headSize;
		uint32_t tailOffset = tailSize;

		//build send buffers
		map<int,Class_Comm_Buffer> sendBuffers;

		//Compute first predecessor and first successor to send buffers to
		int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
		int64_t globalLastHead = (int64_t) lh;
		int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
		int firstPredecessor = -1;
		int firstSuccessor = nproc;
		if(rank != 0){
			firstOctantGlobalIdx = (int64_t)(partition_range_globalidx[rank-1] + 1);
			globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
			globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
			for(int pre = rank - 1; pre >=0; --pre){
				if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
					firstPredecessor = pre;
			}
			for(int post = rank + 1; post < nproc; ++post){
				if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
					firstSuccessor = post;
			}
		}
		else if(rank == 0){
			firstSuccessor = 1;
		}
		MPI_Barrier(comm); //da spostare prima della prima comunicazione

		uint32_t x,y,z;
		uint8_t l;
		int8_t m;
		bool info[17];
		int intBuffer = 0;
		int contatore = 0;
		//build send buffers from Head
		uint32_t nofElementsFromSuccessiveToPrevious = 0;
		if(headSize != 0){
			for(int p = firstPredecessor; p >= 0; --p){
				if(headSize < partition[p]){

					intBuffer = (newPartitionRangeGlobalidx[p] - partition[p] );
					intBuffer = abs(intBuffer);
					nofElementsFromSuccessiveToPrevious = globalLastHead - intBuffer;
					if(nofElementsFromSuccessiveToPrevious > headSize || contatore == 1)
						nofElementsFromSuccessiveToPrevious  = headSize;

					//						int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					//for(uint32_t i = 0; i <= (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); ++i){
					for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
						//PACK octants from 0 to lh in sendBuffer[p]
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int ii = 0; ii < 17; ++ii)
							info[ii] = octant.info[ii];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}

					if(nofElementsFromSuccessiveToPrevious == headSize)
						break;

					lh -= nofElementsFromSuccessiveToPrevious;
					globalLastHead -= nofElementsFromSuccessiveToPrevious;
					headSize = lh + 1;
					++contatore;

				}
				else{
					nofElementsFromSuccessiveToPrevious = globalLastHead - (newPartitionRangeGlobalidx[p] - partition[p]);
					int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
						//pack octants from lh - partition[p] to lh
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int i = 0; i < 17; ++i)
							info[i] = octant.info[i];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					lh -= nofElementsFromSuccessiveToPrevious;
					globalLastHead -= nofElementsFromSuccessiveToPrevious;
					headSize = lh + 1;
					if(headSize == 0)
						break;
				}
			}

		}
		uint32_t nofElementsFromPreviousToSuccessive = 0;
		contatore = 0;
		//build send buffers from Tail
		if(tailSize != 0){
			for(int p = firstSuccessor; p < nproc; ++p){
				if(tailSize < partition[p]){

					nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
					if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
						nofElementsFromPreviousToSuccessive = tailSize;

					int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					uint32_t octantsSize = (uint32_t)octree.octants.size();
					for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
						//PACK octants from ft to octantsSize-1
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int ii = 0; ii < 17; ++ii)
							info[ii] = octant.info[ii];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}

					if(nofElementsFromPreviousToSuccessive == tailSize)
						break;
					ft += nofElementsFromPreviousToSuccessive;
					globalFirstTail += nofElementsFromPreviousToSuccessive;
					tailSize -= nofElementsFromPreviousToSuccessive;
					++contatore;
				}
				else{
					nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
					//						int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					//int buffSize = partition[p] * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive  - 1;
					int pos = 0;
					for(uint32_t i = ft; i <= endOctants; ++i ){
						//PACK octants from ft to ft + partition[p] -1
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int i = 0; i < 17; ++i)
							info[i] = octant.info[i];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					//ft += partition[p];
					//tailSize -= partition[p];
					ft += nofElementsFromPreviousToSuccessive;
					globalFirstTail += nofElementsFromPreviousToSuccessive;
					tailSize -= nofElementsFromPreviousToSuccessive;
					if(tailSize == 0)
						break;
				}
			}
		}

		//Build receiver sources
		vector<Class_Array> recvs(nproc);
		recvs[rank] = Class_Array((uint32_t)sendBuffers.size()+1,-1);
		recvs[rank].array[0] = rank;
		int counter = 1;
		map<int,Class_Comm_Buffer>::iterator sitend = sendBuffers.end();
		for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
			recvs[rank].array[counter] = sit->first;
			++counter;
		}
		int* nofRecvsPerProc = new int[nproc];
		error_flag = MPI_Allgather(&recvs[rank].arraySize,1,MPI_INT,nofRecvsPerProc,1,MPI_INT,comm);
		int globalRecvsBuffSize = 0;
		int* displays = new int[nproc];
		for(int pp = 0; pp < nproc; ++pp){
			displays[pp] = 0;
			globalRecvsBuffSize += nofRecvsPerProc[pp];
			for(int ppp = 0; ppp < pp; ++ppp){
				displays[pp] += nofRecvsPerProc[ppp];
			}
		}
		int* globalRecvsBuff = new int[globalRecvsBuffSize];
		error_flag = MPI_Allgatherv(recvs[rank].array,recvs[rank].arraySize,MPI_INT,globalRecvsBuff,nofRecvsPerProc,displays,MPI_INT,comm);

		vector<set<int> > sendersPerProc(nproc);
		for(int pin = 0; pin < nproc; ++pin){
			for(int k = displays[pin]+1; k < displays[pin] + nofRecvsPerProc[pin]; ++k){
				sendersPerProc[globalRecvsBuff[k]].insert(globalRecvsBuff[displays[pin]]);
			}
		}

		//Communicate Octants (size)
		MPI_Request* req = new MPI_Request[sendBuffers.size()+sendersPerProc[rank].size()];
		MPI_Status* stats = new MPI_Status[sendBuffers.size()+sendersPerProc[rank].size()];
		int nReq = 0;
		map<int,int> recvBufferSizePerProc;
		set<int>::iterator senditend = sendersPerProc[rank].end();
		for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
			recvBufferSizePerProc[*sendit] = 0;
			error_flag = MPI_Irecv(&recvBufferSizePerProc[*sendit],1,MPI_UINT32_T,*sendit,rank,comm,&req[nReq]);
			++nReq;
		}
		map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
		for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//COMMUNICATE THE BUFFERS TO THE RECEIVERS
		//recvBuffers structure is declared and each buffer is initialized to the right size
		//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
		uint32_t nofNewHead = 0;
		uint32_t nofNewTail = 0;
		map<int,Class_Comm_Buffer> recvBuffers;
		map<int,int>::iterator ritend = recvBufferSizePerProc.end();
		for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
			recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a',comm);
			uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)global.octantBytes / (double)(CHAR_BIT/8)));
			if(rit->first < rank)
				nofNewHead += nofNewPerProc;
			else if(rit->first > rank)
				nofNewTail += nofNewPerProc;
		}
		nReq = 0;
		for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
			//nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
			error_flag = MPI_Irecv(recvBuffers[*sendit].commBuffer,recvBuffers[*sendit].commBufferSize,MPI_PACKED,*sendit,rank,comm,&req[nReq]);
			++nReq;
		}
		for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//MOVE RESIDENT TO BEGIN IN OCTANTS
		uint32_t resEnd = octree.getNumOctants() - tailOffset;
		uint32_t nofResidents = resEnd - headOffset;
		int octCounter = 0;
		for(uint32_t i = headOffset; i < resEnd; ++i){
			octree.octants[octCounter] = octree.octants[i];
			++octCounter;
		}
		uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
		octree.octants.resize(newCounter);
		//MOVE RESIDENTS IN RIGHT POSITION
		uint32_t resCounter = nofNewHead + nofResidents - 1;
		for(uint32_t k = 0; k < nofResidents ; ++k){
			octree.octants[resCounter - k] = octree.octants[nofResidents - k - 1];
		}

		//UNPACK BUFFERS AND BUILD NEW OCTANTS
		newCounter = 0;
		bool jumpResident = false;
		map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
		for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
			uint32_t nofNewPerProc = (uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)global.octantBytes / (double)(CHAR_BIT/8)));
			int pos = 0;
			if(rbit->first > rank && !jumpResident){
				newCounter += nofResidents ;
				jumpResident = true;
			}
			for(int i = nofNewPerProc - 1; i >= 0; --i){
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&z,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
				octree.octants[newCounter] = classOctant(dim,l,x,y,z);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
				octree.octants[newCounter].setMarker(m);
				for(int j = 0; j < 17; ++j){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,comm);
					octree.octants[newCounter].info[j] = info[j];
				}
				++newCounter;
			}
		}
		octvector(octree.octants).swap(octree.octants);

		delete [] newPartitionRangeGlobalidx; newPartitionRangeGlobalidx = NULL;
		delete [] nofRecvsPerProc; nofRecvsPerProc = NULL;
		delete [] displays; displays = NULL;
		delete [] req; req = NULL;
		delete [] stats; stats = NULL;
		delete [] globalRecvsBuff; globalRecvsBuff = NULL;
		//Update and ghosts here
		updateLoadBalance();
		setPboundGhosts();

	}
	delete [] partition;
	partition = NULL;

	//Write info of final partition on log
	log.writeLog(" ");
	log.writeLog(" Final Parallel partition : ");
	log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[0]+1)));
	for(int ii=1; ii<nproc; ii++){
		log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]-partition_range_globalidx[ii-1])));
	}
	log.writeLog(" ");
	log.writeLog("---------------------------------------------");

}

/** Distribute Load-Balanced the octants of the whole tree over
 * the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
 * The families of octants of a desired level are retained compact on the same process.
 * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
 */
void
classParaTree::loadBalance(uint8_t & level){

	//Write info on log
	log.writeLog("---------------------------------------------");
	log.writeLog(" LOAD BALANCE ");

	uint32_t* partition = new uint32_t [nproc];
	computePartition(partition, level);
	if(serial)
	{
		log.writeLog(" ");
		log.writeLog(" Initial Serial distribution : ");
		for(int ii=0; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]+1)));
		}

		uint32_t stride = 0;
		for(int i = 0; i < rank; ++i)
			stride += partition[i];
		classLocalTree::octvector octantsCopy = octree.octants;
		classLocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
		classLocalTree::octvector::const_iterator last = first + partition[rank];
		octree.octants.assign(first, last);
		octvector(octree.octants).swap(octree.octants);

		first = octantsCopy.end();
		last = octantsCopy.end();

		//Update and ghosts here
		updateLoadBalance();
		setPboundGhosts();

	}
	else
	{
		log.writeLog(" ");
		log.writeLog(" Initial Parallel partition : ");
		log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[0]+1)));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]-partition_range_globalidx[ii-1])));
		}

		//empty ghosts
		octree.ghosts.clear();
		octree.size_ghosts = 0;
		//compute new partition range globalidx
		uint64_t* newPartitionRangeGlobalidx = new uint64_t[nproc];
		for(int p = 0; p < nproc; ++p){
			newPartitionRangeGlobalidx[p] = 0;
			for(int pp = 0; pp <= p; ++pp)
				newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
			--newPartitionRangeGlobalidx[p];
		}

		//find resident octants local offset lastHead(lh) and firstTail(ft)
		int32_t lh,ft;
		if(rank == 0)
			lh = -1;
		else{
			lh = (int32_t)(newPartitionRangeGlobalidx[rank-1] + 1 - partition_range_globalidx[rank-1] - 1 - 1);
		}
		if(lh < 0)
			lh = - 1;
		else if(lh > (int32_t)(octree.octants.size() - 1))
			lh = octree.octants.size() - 1;

		if(rank == nproc - 1)
			ft = octree.octants.size();
		else if(rank == 0)
			ft = (int32_t)(newPartitionRangeGlobalidx[rank] + 1);
		else{
			ft = (int32_t)(newPartitionRangeGlobalidx[rank] - partition_range_globalidx[rank -1]);
		}
		if(ft > (int32_t)(octree.octants.size() - 1))
			ft = octree.octants.size();
		else if(ft < 0)
			ft = 0;

		//compute size Head and size Tail
		uint32_t headSize = (uint32_t)(lh + 1);
		uint32_t tailSize = (uint32_t)(octree.octants.size() - ft);
		uint32_t headOffset = headSize;
		uint32_t tailOffset = tailSize;

		//build send buffers
		map<int,Class_Comm_Buffer> sendBuffers;

		//Compute first predecessor and first successor to send buffers to
		int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
		int64_t globalLastHead = (int64_t) lh;
		int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
		int firstPredecessor = -1;
		int firstSuccessor = nproc;
		if(rank != 0){
			firstOctantGlobalIdx = (int64_t)(partition_range_globalidx[rank-1] + 1);
			globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
			globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
			for(int pre = rank - 1; pre >=0; --pre){
				if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
					firstPredecessor = pre;
			}
			for(int post = rank + 1; post < nproc; ++post){
				if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
					firstSuccessor = post;
			}
		}
		else if(rank == 0){
			firstSuccessor = 1;
		}
		MPI_Barrier(comm); //da spostare prima della prima comunicazione

		uint32_t x,y,z;
		uint8_t l;
		int8_t m;
		bool info[17];
		int intBuffer = 0;
		int contatore = 0;
		//build send buffers from Head
		uint32_t nofElementsFromSuccessiveToPrevious = 0;
		if(headSize != 0){
			for(int p = firstPredecessor; p >= 0; --p){
				if(headSize < partition[p]){
					intBuffer = (newPartitionRangeGlobalidx[p] - partition[p] );
					intBuffer = abs(intBuffer);
					nofElementsFromSuccessiveToPrevious = globalLastHead - intBuffer;
					if(nofElementsFromSuccessiveToPrevious > headSize || contatore == 1)
						nofElementsFromSuccessiveToPrevious  = headSize;

					int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
						//PACK octants from 0 to lh in sendBuffer[p]
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int ii = 0; ii < 17; ++ii)
							info[ii] = octant.info[ii];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					if(nofElementsFromSuccessiveToPrevious == headSize)
						break;

					lh -= nofElementsFromSuccessiveToPrevious;
					globalLastHead -= nofElementsFromSuccessiveToPrevious;
					headSize = lh + 1;
					++contatore;
				}
				else{
					nofElementsFromSuccessiveToPrevious = globalLastHead - (newPartitionRangeGlobalidx[p] - partition[p]);
					int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
						//pack octants from lh - partition[p] to lh
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int i = 0; i < 17; ++i)
							info[i] = octant.info[i];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					lh -= nofElementsFromSuccessiveToPrevious;
					globalLastHead -= nofElementsFromSuccessiveToPrevious;
					headSize = lh + 1;
					if(headSize == 0)
						break;
				}
			}

		}
		uint32_t nofElementsFromPreviousToSuccessive = 0;
		contatore = 0;
		//build send buffers from Tail
		if(tailSize != 0){
			for(int p = firstSuccessor; p < nproc; ++p){
				if(tailSize < partition[p]){
					nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
					if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
						nofElementsFromPreviousToSuccessive = tailSize;

					int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					uint32_t octantsSize = (uint32_t)octree.octants.size();
					for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
						//PACK octants from ft to octantsSize-1
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int ii = 0; ii < 17; ++ii)
							info[ii] = octant.info[ii];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					if(nofElementsFromPreviousToSuccessive == tailSize)
						break;
					ft += nofElementsFromPreviousToSuccessive;
					globalFirstTail += nofElementsFromPreviousToSuccessive;
					tailSize -= nofElementsFromPreviousToSuccessive;
					++contatore;
				}
				else{
					nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
					int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
					uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
					int pos = 0;
					for(uint32_t i = ft; i <= endOctants; ++i ){
						//PACK octants from ft to ft + partition[p] -1
						const classOctant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						for(int ii = 0; ii < 17; ++ii)
							info[ii] = octant.info[ii];
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						for(int j = 0; j < 17; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
						}
					}
					ft += nofElementsFromPreviousToSuccessive;
					globalFirstTail += nofElementsFromPreviousToSuccessive;
					tailSize -= nofElementsFromPreviousToSuccessive;
					if(tailSize == 0)
						break;
				}
			}
		}

		//Build receiver sources
		vector<Class_Array> recvs(nproc);
		recvs[rank] = Class_Array((uint32_t)sendBuffers.size()+1,-1);
		recvs[rank].array[0] = rank;
		int counter = 1;
		map<int,Class_Comm_Buffer>::iterator sitend = sendBuffers.end();
		for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
			recvs[rank].array[counter] = sit->first;
			++counter;
		}
		int* nofRecvsPerProc = new int[nproc];
		error_flag = MPI_Allgather(&recvs[rank].arraySize,1,MPI_INT,nofRecvsPerProc,1,MPI_INT,comm);
		int globalRecvsBuffSize = 0;
		int* displays = new int[nproc];
		for(int pp = 0; pp < nproc; ++pp){
			displays[pp] = 0;
			globalRecvsBuffSize += nofRecvsPerProc[pp];
			for(int ppp = 0; ppp < pp; ++ppp){
				displays[pp] += nofRecvsPerProc[ppp];
			}
		}
		int* globalRecvsBuff = new int[globalRecvsBuffSize];
		error_flag = MPI_Allgatherv(recvs[rank].array,recvs[rank].arraySize,MPI_INT,globalRecvsBuff,nofRecvsPerProc,displays,MPI_INT,comm);

		vector<set<int> > sendersPerProc(nproc);
		for(int pin = 0; pin < nproc; ++pin){
			for(int k = displays[pin]+1; k < displays[pin] + nofRecvsPerProc[pin]; ++k){
				sendersPerProc[globalRecvsBuff[k]].insert(globalRecvsBuff[displays[pin]]);
			}
		}

		//Communicate Octants (size)
		MPI_Request* req = new MPI_Request[sendBuffers.size()+sendersPerProc[rank].size()];
		MPI_Status* stats = new MPI_Status[sendBuffers.size()+sendersPerProc[rank].size()];
		int nReq = 0;
		map<int,int> recvBufferSizePerProc;
		set<int>::iterator senditend = sendersPerProc[rank].end();
		for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
			recvBufferSizePerProc[*sendit] = 0;
			error_flag = MPI_Irecv(&recvBufferSizePerProc[*sendit],1,MPI_UINT32_T,*sendit,rank,comm,&req[nReq]);
			++nReq;
		}
		map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
		for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//COMMUNICATE THE BUFFERS TO THE RECEIVERS
		//recvBuffers structure is declared and each buffer is initialized to the right size
		//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
		uint32_t nofNewHead = 0;
		uint32_t nofNewTail = 0;
		map<int,Class_Comm_Buffer> recvBuffers;
		map<int,int>::iterator ritend = recvBufferSizePerProc.end();
		for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
			recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a',comm);
			uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)global.octantBytes / (double)(CHAR_BIT/8)));
			if(rit->first < rank)
				nofNewHead += nofNewPerProc;
			else if(rit->first > rank)
				nofNewTail += nofNewPerProc;
		}
		nReq = 0;
		for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
			//nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
			error_flag = MPI_Irecv(recvBuffers[*sendit].commBuffer,recvBuffers[*sendit].commBufferSize,MPI_PACKED,*sendit,rank,comm,&req[nReq]);
			++nReq;
		}
		for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//MOVE RESIDENT TO BEGIN IN OCTANTS
		uint32_t resEnd = octree.getNumOctants() - tailOffset;
		uint32_t nofResidents = resEnd - headOffset;
		int octCounter = 0;
		for(uint32_t i = headOffset; i < resEnd; ++i){
			octree.octants[octCounter] = octree.octants[i];
			++octCounter;
		}
		uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
		octree.octants.resize(newCounter);
		//MOVE RESIDENTS IN RIGHT POSITION
		uint32_t resCounter = nofNewHead + nofResidents - 1;
		for(uint32_t k = 0; k < nofResidents ; ++k){
			octree.octants[resCounter - k] = octree.octants[nofResidents - k - 1];
		}

		//UNPACK BUFFERS AND BUILD NEW OCTANTS
		newCounter = 0;
		bool jumpResident = false;
		map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
		for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
			uint32_t nofNewPerProc = (uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)global.octantBytes / (double)(CHAR_BIT/8)));
			int pos = 0;
			if(rbit->first > rank && !jumpResident){
				newCounter += nofResidents ;
				jumpResident = true;
			}
			for(int i = nofNewPerProc - 1; i >= 0; --i){
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&z,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
				octree.octants[newCounter] = classOctant(dim,l,x,y,z);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
				octree.octants[newCounter].setMarker(m);
				for(int j = 0; j < 17; ++j){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,comm);
					octree.octants[newCounter].info[j] = info[j];
				}
				++newCounter;
			}
		}
		octvector(octree.octants).swap(octree.octants);

		delete [] newPartitionRangeGlobalidx; newPartitionRangeGlobalidx = NULL;
		delete [] nofRecvsPerProc; nofRecvsPerProc = NULL;
		delete [] displays; displays = NULL;
		delete [] req; req = NULL;
		delete [] stats; stats = NULL;
		delete [] globalRecvsBuff; globalRecvsBuff = NULL;
		//Update and ghosts here
		updateLoadBalance();
		setPboundGhosts();

	}
	delete [] partition;
	partition = NULL;

	//Write info of final partition on log
	log.writeLog(" ");
	log.writeLog(" Final Parallel partition : ");
	log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[0]+1)));
	for(int ii=1; ii<nproc; ii++){
		log.writeLog(" Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(partition_range_globalidx[ii]-partition_range_globalidx[ii-1])));
	}
	log.writeLog(" ");
	log.writeLog("---------------------------------------------");

}
#endif

// =================================================================================== //
// OTHER INTERSECTION BASED METHODS												    			   //
// =================================================================================== //

/** Compute the intersection of octants (intersections of bord, of inner domain and with ghost octants).
 */
void
classParaTree::computeIntersections(){
	octree.computeIntersections();
}

// =================================================================================== //
// OTHER PRIVATE METHODS												    			   //
// =================================================================================== //

classOctant&
classParaTree::extractOctant(uint32_t idx) {
	return octree.extractOctant(idx) ;
};

/** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
 */
//bool
//classParaTree::private_adapt() {
//	bool globalDone = false, localDone = false, cDone = false;
//	uint32_t nocts = octree.getNumOctants();
//	vector<classOctant >::iterator iter, iterend = octree.octants.end();
//
//	mapidx.clear();
//
//	for (iter = octree.octants.begin(); iter != iterend; iter++){
//		iter->info[12] = false;
//		iter->info[13] = false;
//		iter->info[15] = false;
//	}
//#if NOMPI==0
//	if(serial){
//#endif
//		log.writeLog("---------------------------------------------");
//		log.writeLog(" ADAPT (Refine/Coarse)");
//		log.writeLog(" ");
//
//		// 2:1 Balance
//		balance21(true);
//
//		log.writeLog(" ");
//		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));
//
//		// Refine
//		while(octree.refine(mapidx));
//
//		if (octree.getNumOctants() > nocts)
//			localDone = true;
//		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));
//		nocts = octree.getNumOctants();
//		updateAdapt();
//
//		// Coarse
//		while(octree.coarse(mapidx));
//		updateAfterCoarse();
//		//			balance21(false);
//		//			while(octree.refine());
//		//			updateAdapt();
//		if (octree.getNumOctants() < nocts){
//			localDone = true;
//		}
//		nocts = octree.getNumOctants();
//
//		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(nocts)));
//#if NOMPI==0
//		MPI_Barrier(comm);
//		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
//#endif
//		log.writeLog(" ");
//		log.writeLog("---------------------------------------------");
//#if NOMPI==0
//	}
//	else{
//		log.writeLog("---------------------------------------------");
//		log.writeLog(" ADAPT (Refine/Coarse)");
//		log.writeLog(" ");
//
//		// 2:1 Balance
//		balance21(true);
//
//		log.writeLog(" ");
//		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
//
//		// Refine
//		while(octree.refine(mapidx));
//		if (octree.getNumOctants() > nocts)
//			localDone = true;
//		//cout << rank << " refine done " << localDone << endl;
//		updateAdapt();
//		setPboundGhosts();
//		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
//		nocts = octree.getNumOctants();
//
//		// Coarse
//		while(octree.coarse(mapidx));
//		//		cout << rank << " out coarse " << endl;
//		updateAfterCoarse();
//		//		cout << rank << " out updatecoarse " << endl;
//		setPboundGhosts();
//		//		cout << rank << " out setghosts " << endl;
//		//			balance21(false);
//		//			while(octree.refine());
//		//			updateAdapt();
//		//			setPboundGhosts();
//		if (octree.getNumOctants() < nocts){
//			localDone = true;
//			//cout << rank << " coarse done " << localDone << endl;
//		}
//		nocts = octree.getNumOctants();
//
//		MPI_Barrier(comm);
//		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
//		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
//		log.writeLog(" ");
//		log.writeLog("---------------------------------------------");
//	}
//	status += globalDone;
//	return globalDone;
//#else
//	status += localDone;
//	return localDone;
//#endif
//}

/** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
 * Track the changes in structure octant by a mapper.
 */
bool
classParaTree::private_adapt_mapidx(bool mapflag) {
	//TODO recoding for adapting with abs(marker) > 1

	bool globalDone = false, localDone = false;
	bool refine = true, coarse = true, globalCoarse = true;
	uint32_t nocts = octree.getNumOctants();
	u32vector mapidx_temp, mapidx_temp2;
	vector<classOctant >::iterator iter, iterend = octree.octants.end();

	for (iter = octree.octants.begin(); iter != iterend; iter++){
		iter->info[12] = false;
		iter->info[13] = false;
		iter->info[15] = false;
	}

	// mapidx init
	u32vector().swap(mapidx);
	u32vector().swap(mapidx_temp);
	u32vector().swap(mapidx_temp2);
	if (mapflag) {
		mapidx.resize(nocts);
		mapidx_temp.resize(nocts);
		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
			mapidx_temp[i] = i;
		}
	}

#if NOMPI==0
	if(serial){
#endif
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Refine/Coarse)");
		log.writeLog(" ");

		// 2:1 Balance
		balance21(true);

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));

		// Refine
		//		while(octree.refine(mapidx));
		while (refine) {
			refine = octree.refine(mapidx_temp);
			mapidx_temp2.resize(octree.getNumOctants());
			for (uint32_t i=0; i<octree.getNumOctants(); i++){
				mapidx_temp2[mapidx_temp[i]] = mapidx[mapidx_temp[i]];
			}
			mapidx.clear();
			mapidx = mapidx_temp2;
		}
		if (octree.getNumOctants() > nocts)
			localDone = true;
		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(octree.getNumOctants())));
		nocts = octree.getNumOctants();
		updateAdapt();

		// Coarse
		while(octree.coarse(mapidx));
		updateAfterCoarse(mapidx);
		//			balance21(false);
		//			while(octree.refine(mapidx));
		//			updateAdapt();
		if (octree.getNumOctants() < nocts){
			localDone = true;
		}
		nocts = octree.getNumOctants();

		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(nocts)));
#if NOMPI==0
		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
#endif
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
#if NOMPI==0
	}
	else{
		log.writeLog("---------------------------------------------");
		log.writeLog(" ADAPT (Refine/Coarse)");
		log.writeLog(" ");

		// 2:1 Balance
		balance21(true);

		log.writeLog(" ");
		log.writeLog(" Initial Number of octants	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));

		// Refine
//		while(octree.refine(mapidx));
		while (refine) {
			refine = octree.refine(mapidx_temp);
			if (mapflag){
				mapidx_temp2.resize(octree.getNumOctants());
				for (uint32_t i=0; i<octree.getNumOctants(); i++){
					mapidx_temp2[mapidx_temp[i]] = mapidx[mapidx_temp[i]];
				}
				mapidx.clear();
				mapidx = mapidx_temp2;
			}
		}
		if (octree.getNumOctants() > nocts)
			localDone = true;
		updateAdapt();
		setPboundGhosts();
		log.writeLog(" Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
		nocts = octree.getNumOctants();


		// Coarse
//		while(octree.coarse(mapidx));
		while (globalCoarse) {
			coarse = octree.coarse(mapidx_temp);
			if (mapflag){
				mapidx_temp2.resize(octree.getNumOctants());
				for (uint32_t i=0; i<octree.getNumOctants(); i++){
					mapidx_temp2[mapidx_temp[i]] = mapidx[mapidx_temp[i]];
				}
				mapidx.clear();
				mapidx = mapidx_temp2;
			}
			updateAfterCoarse(mapidx);
			setPboundGhosts();
			globalCoarse = false;
			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&coarse,&globalCoarse,1,MPI::BOOL,MPI_LOR,comm);
		}
//		updateAfterCoarse(mapidx);
//		setPboundGhosts();
		//			balance21(false);
		//			while(octree.refine(mapidx));
		//			updateAdapt();
		//			setPboundGhosts();
		if (octree.getNumOctants() < nocts){
			localDone = true;
		}
		nocts = octree.getNumOctants();

		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
		log.writeLog(" Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(global_num_octants)));
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");
	}
	return globalDone;
#else
	return localDone;
#endif
}

void
classParaTree::updateAdapt(){
#if NOMPI==0
	if(serial)
	{
#endif
		max_depth = octree.local_max_depth;
		global_num_octants = octree.getNumOctants();
		for(int p = 0; p < nproc; ++p){
			partition_range_globalidx[p] = global_num_octants - 1;
		}
#if NOMPI==0
	}
	else
	{
		//update max_depth
		error_flag = MPI_Allreduce(&octree.local_max_depth,&max_depth,1,MPI_UINT8_T,MPI_MAX,comm);
		//update global_num_octants
		uint64_t local_num_octants = octree.getNumOctants();
		error_flag = MPI_Allreduce(&local_num_octants,&global_num_octants,1,MPI_UINT64_T,MPI_SUM,comm);
		//update partition_range_globalidx
		uint64_t* rbuff = new uint64_t[nproc];
		error_flag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,rbuff,1,MPI_UINT64_T,comm);
		for(int p = 0; p < nproc; ++p){
			partition_range_globalidx[p] = 0;
			for(int pp = 0; pp <=p; ++pp)
				partition_range_globalidx[p] += rbuff[pp];
			--partition_range_globalidx[p];
		}
		//update partition_range_position
		uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
		error_flag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,partition_last_desc,1,MPI_UINT64_T,comm);
		uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
		error_flag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,partition_first_desc,1,MPI_UINT64_T,comm);
		delete [] rbuff; rbuff = NULL;
	}
#endif
}

#if NOMPI==0
void
classParaTree::computePartition(uint32_t* partition){

	uint32_t division_result = 0;
	uint32_t remind = 0;

	division_result = uint32_t(global_num_octants/(uint64_t)nproc);
	remind = (uint32_t)(global_num_octants%(uint64_t)nproc);

	for(uint32_t i = 0; i < (uint32_t)nproc; ++i)
		if(i<remind)
			partition[i] = division_result + 1;
		else
			partition[i] = division_result;

}

/*! compute octant partition giving the same weight to
 * each process and redistributing the reminder
 */
void
classParaTree::computePartition(uint32_t* partition, dvector* weight){
	if(serial){

		double division_result = 0;
		double global_weight = 0.0;
		for (int i=0; i<weight->size(); i++){
			global_weight += (*weight)[i];
		}
		division_result = global_weight/(double)nproc;

		//Estimate resulting weight distribution starting from proc 0 (sending tail)
		//Estimate sending weight by each proc in initial conf (sending tail)
		uint32_t i = 0, tot = 0;
		int iproc = 0;
		while (iproc < nproc-1){
			double partial_weight = 0.0;
			partition[iproc] = 0;
			while(partial_weight < division_result){
				partial_weight += (*weight)[i];
				tot++;
				partition[iproc]++;
				i++;
			}
			iproc++;
		}
		partition[nproc-1] = weight->size() - tot;
	}
	else{

		double division_result = 0;
		double remind = 0;
		dvector local_weight(nproc,0.0);
		dvector temp_local_weight(nproc,0.0);
		dvector2D sending_weight(nproc, dvector(nproc,0.0));
		double* rbuff = new double[nproc];
		double global_weight = 0.0;
		for (int i=0; i<weight->size(); i++){
			local_weight[rank] += (*weight)[i];
		}
		error_flag = MPI_Allgather(&local_weight[rank],1,MPI_DOUBLE,rbuff,1,MPI_DOUBLE,comm);
		for (int i=0; i<nproc; i++){
			local_weight[i] = rbuff[i];
			global_weight += rbuff[i];
		}
		delete [] rbuff; rbuff = NULL;
		division_result = global_weight/(double)nproc;

		//Estimate resulting weight distribution starting from proc 0 (sending tail)

		temp_local_weight = local_weight;
		//Estimate sending weight by each proc in initial conf (sending tail)

		for (int iter = 0; iter < 1; iter++){

			vector<double> delta(nproc);
			for (int i=0; i<nproc; i++){
				delta[i] = temp_local_weight[i] - division_result;
			}

			for (int i=0; i<nproc-1; i++){

				double post_weight = 0.0;
				for (int j=i+1; j<nproc; j++){
					post_weight += temp_local_weight[j];
				}
				if (temp_local_weight[i] > division_result){

					delta[i] = temp_local_weight[i] - division_result;
					if (post_weight < division_result*(nproc-i-1)){

						double post_delta =  division_result*(nproc-i-1) - post_weight;
						double delta_sending = min(local_weight[i], min(delta[i], post_delta));
						int jproc = i+1;
						double sending = 0;
						while (delta_sending > 0 && jproc<nproc){
							sending = min(division_result, delta_sending);
							sending = min(sending, (division_result-temp_local_weight[jproc]));
							sending = max(sending, 0.0);
							sending_weight[i][jproc] += sending;
							temp_local_weight[jproc] += sending;
							temp_local_weight[i] -= sending;
							delta_sending -= sending;
							delta[i] -= delta_sending;
							jproc++;
						}
					} //post
				}//weight>
			}//iproc

			for (int i = nproc-1; i>0; i--){

				double pre_weight = 0.0;
				for (int j=i-1; j>=0; j--){
					pre_weight += temp_local_weight[j];
				}
				if (temp_local_weight[i] > division_result){

					delta[i] = temp_local_weight[i] - division_result;
					if (pre_weight < division_result*(i)){

						double pre_delta =  division_result*(i) - pre_weight;
						double delta_sending = min(local_weight[i], min(temp_local_weight[i], min(delta[i], pre_delta)));
						int jproc = i-1;
						double sending = 0;
						while (delta_sending > 0 && jproc >=0){
							sending = min(division_result, delta_sending);
							sending = min(sending, (division_result-temp_local_weight[jproc]));
							sending = max(sending, 0.0);
							sending_weight[i][jproc] += sending;
							temp_local_weight[jproc] += sending;
							temp_local_weight[i] -= sending;
							delta_sending -= sending;
							delta[i] -= delta_sending;
							jproc--;
						}
					}//pre
				}//weight>
			}//iproc
		}//iter

		//Update partition locally
		//to send
		u32vector sending_cell(nproc,0);
		int i = getNumOctants();;
		for (int jproc=nproc-1; jproc>rank; jproc--){
			double pack_weight = 0.0;
			while(pack_weight < sending_weight[rank][jproc] && i > 0){
				i--;
				pack_weight += (*weight)[i];
				sending_cell[jproc]++;
			}
		}
		partition[rank] = i;
		i = 0;
		for (int jproc=0; jproc<rank; jproc++){
			double pack_weight = 0.0;
			while(pack_weight < sending_weight[rank][jproc] && i <  getNumOctants()-1){
				i++;
				pack_weight += (*weight)[i];
				sending_cell[jproc]++;
			}
		}
		partition[rank] -= i;

		//to receive
		u32vector rec_cell(nproc,0);
		MPI_Request* req = new MPI_Request[nproc*10];
		MPI_Status* stats = new MPI_Status[nproc*10];
		int nReq = 0;
		for (int iproc=0; iproc<nproc; iproc++){
			error_flag = MPI_Irecv(&rec_cell[iproc],1,MPI_UINT32_T,iproc,rank,comm,&req[nReq]);
			++nReq;
		}
		for (int iproc=0; iproc<nproc; iproc++){
			error_flag =  MPI_Isend(&sending_cell[iproc],1,MPI_UINT32_T,iproc,iproc,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		delete [] req; req = NULL;
		delete [] stats; stats = NULL;

		i = 0;
		for (int jproc=0; jproc<nproc; jproc++){
			i+= rec_cell[jproc];
		}
		partition[rank] += i;
		uint32_t part = partition[rank];
		error_flag = MPI_Allgather(&part,1,MPI_UINT32_T,partition,1,MPI_UINT32_T,comm);
	}
};

void
classParaTree::computePartition(uint32_t* partition, uint8_t & level_) {

	uint8_t level = uint8_t(min(int(max(int(max_depth) - int(level_), int(1))) , int(global.MAX_LEVEL)));
	uint32_t* partition_temp = new uint32_t[nproc];
	uint8_t* boundary_proc = new uint8_t[nproc-1];
	uint8_t dimcomm, indcomm;
	uint8_t* glbdimcomm = new uint8_t[nproc];
	uint8_t* glbindcomm = new uint8_t[nproc];

	uint32_t division_result = 0;
	uint32_t remind = 0;
	uint32_t Dh = uint32_t(pow(double(2),double(global.MAX_LEVEL-level)));
	uint32_t istart, nocts, rest, forw, backw;
	uint32_t i = 0, iproc, j;
	uint64_t sum;
	int32_t* pointercomm;
	int32_t* deplace = new int32_t[nproc-1];
	division_result = uint32_t(global_num_octants/(uint64_t)nproc);
	remind = (uint32_t)(global_num_octants%(uint64_t)nproc);
	for(uint32_t i = 0; i < (uint32_t)nproc; ++i)
		if(i<remind)
			partition_temp[i] = division_result + 1;
		else
			partition_temp[i] = division_result;

	j = 0;
	sum = 0;
	for (iproc=0; iproc<(uint32_t)(nproc-1); iproc++){
		sum += partition_temp[iproc];
		while(sum > partition_range_globalidx[j]){
			j++;
		}
		boundary_proc[iproc] = j;
	}
	nocts = octree.octants.size();
	sum = 0;
	dimcomm = 0;
	indcomm = 0;
	for (iproc=0; iproc<(uint32_t)(nproc-1); iproc++){
		deplace[iproc] = 1;
		sum += partition_temp[iproc];
		if (boundary_proc[iproc] == rank){
			if (dimcomm == 0){
				indcomm = iproc;
			}
			dimcomm++;
			if (rank!=0)
				istart = sum - partition_range_globalidx[rank-1] - 1;
			else
				istart = sum;

			i = istart;
			rest = octree.octants[i].getX()%Dh + octree.octants[i].getY()%Dh;
			while(rest!=0){
				if (i==nocts){
					i = istart + nocts;
					break;
				}
				i++;
				rest = octree.octants[i].getX()%Dh + octree.octants[i].getY()%Dh;
			}
			forw = i - istart;
			i = istart;
			rest = octree.octants[i].getX()%Dh + octree.octants[i].getY()%Dh;
			while(rest!=0){
				if (i==0){
					i = istart - nocts;
					break;
				}
				i--;
				rest = octree.octants[i].getX()%Dh + octree.octants[i].getY()%Dh;
			}
			backw = istart - i;
			if (forw<backw)
				deplace[iproc] = forw;
			else
				deplace[iproc] = -(int32_t)backw;
		}
	}

	error_flag = MPI_Allgather(&dimcomm,1,MPI_UINT8_T,glbdimcomm,1,MPI_UINT8_T,comm);
	error_flag = MPI_Allgather(&indcomm,1,MPI_UINT8_T,glbindcomm,1,MPI_UINT8_T,comm);
	for (iproc=0; iproc<(uint32_t)(nproc); iproc++){
		pointercomm = &deplace[glbindcomm[iproc]];
		error_flag = MPI_Bcast(pointercomm, glbdimcomm[iproc], MPI_INT32_T, iproc, comm);
	}

	for (iproc=0; iproc<(uint32_t)(nproc); iproc++){
		if (iproc < (uint32_t)(nproc-1))
			partition[iproc] = partition_temp[iproc] + deplace[iproc];
		else
			partition[iproc] = partition_temp[iproc];
		if (iproc !=0)
			partition[iproc] = partition[iproc] - deplace[iproc-1];
	}

	delete [] partition_temp; partition_temp = NULL;
	delete [] boundary_proc; boundary_proc = NULL;
	delete [] glbdimcomm; glbdimcomm = NULL;
	delete [] glbindcomm; glbindcomm = NULL;
	delete [] deplace; deplace = NULL;
}

void
classParaTree::updateLoadBalance() {
	octree.updateLocalMaxDepth();
	uint64_t* rbuff = new uint64_t[nproc];
	uint64_t local_num_octants = octree.getNumOctants();
	error_flag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,rbuff,1,MPI_UINT64_T,comm);
	for(int p = 0; p < nproc; ++p){
		partition_range_globalidx[p] = 0;
		for(int pp = 0; pp <=p; ++pp)
			partition_range_globalidx[p] += rbuff[pp];
		--partition_range_globalidx[p];
	}
	//update first last descendant
	octree.setFirstDesc();
	octree.setLastDesc();
	//update partition_range_position
	uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
	error_flag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,partition_last_desc,1,MPI_UINT64_T,comm);
	uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
	error_flag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,partition_first_desc,1,MPI_UINT64_T,comm);
	serial = false;
	delete [] rbuff; rbuff = NULL;
}

void
classParaTree::setPboundGhosts() {
	//BUILD BORDER OCTANT INDECES VECTOR (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//find local octants to be sent as ghost to the right processes
	//it visits the local octants building virtual neighbors on each octant face
	//find the owner of these virtual neighbor and build a map (process,border octants)
	//this map contains the local octants as ghosts for neighbor processes

	// NO PBORDERS !
	classLocalTree::octvector::iterator end = octree.octants.end();
	classLocalTree::octvector::iterator begin = octree.octants.begin();
	bordersPerProc.clear();
	int count = 0;
	for(classLocalTree::octvector::iterator it = begin; it != end; ++it){
		set<int> procs;
		//if (rank!=0) cout << rank << " border face " << endl;
        //if (rank!=0 && it->getZ()==0) cout << "oct " << count << " info " << it->info << endl;
        //if (rank!=0 && it->getZ()==0) cout << it->getX() << " " << it->getY() << "  " << it->getZ() << endl;
		//Virtual Face Neighbors
		for(uint8_t i = 0; i < global.nfaces; ++i){
			if(it->getBound(i) == false){
				//if (rank!=0 && it->getZ()==0) cout << "oct " << count << " face " << int(i) << " bound " << it->getBound(i) << endl;
				uint32_t virtualNeighborsSize = 0;
				vector<uint64_t> virtualNeighbors = it->computeVirtualMorton(i,max_depth,virtualNeighborsSize,global.MAX_LEVEL);
				uint32_t maxDelta = virtualNeighborsSize/2;
				for(uint32_t j = 0; j <= maxDelta; ++j){
					//if (rank!=0 && it->getZ()==0) cout << "oct " << count << " j " << j << " face " << int(i) << " searching for pBegin " << virtualNeighbors[j] << endl;
					int pBegin = findOwner(virtualNeighbors[j]);
					//if (rank!=0) cout << j << " pBegin " << pBegin << endl;
					//if (rank!=0) cout << j << " searching for pEnd " << virtualNeighbors[virtualNeighborsSize - 1 - j] << endl;
					int pEnd = findOwner(virtualNeighbors[virtualNeighborsSize - 1 - j]);
					//if (rank!=0) cout << j << " pEnd " << pEnd << endl;
					procs.insert(pBegin);
					procs.insert(pEnd);
					if(pBegin != rank || pEnd != rank){
						it->setPbound(i,true);
					}
					else{
						it->setPbound(i,false);
					}
				}
			}
		}
		// if (rank!=0) cout << rank << " border edge " << endl;
		//Virtual Edge Neighbors
		for(uint8_t e = 0; e < global.nedges; ++e){
			uint32_t virtualEdgeNeighborSize = 0;
			vector<uint64_t> virtualEdgeNeighbors = it->computeEdgeVirtualMorton(e,max_depth,virtualEdgeNeighborSize,octree.balance_codim, global.MAX_LEVEL, global.edgeface);
			uint32_t maxDelta = virtualEdgeNeighborSize/2;
			if(virtualEdgeNeighborSize){
				for(uint32_t ee = 0; ee <= maxDelta; ++ee){
					int pBegin = findOwner(virtualEdgeNeighbors[ee]);
					int pEnd = findOwner(virtualEdgeNeighbors[virtualEdgeNeighborSize - 1- ee]);
					procs.insert(pBegin);
					procs.insert(pEnd);
				}
			}
		}
		//if (rank!=0) cout << rank << " border corner " << endl;
		//Virtual Corner Neighbors
		for(uint8_t c = 0; c < global.nnodes; ++c){
			if(!it->getBound(global.nodeface[c][0]) && !it->getBound(global.nodeface[c][1])){
				uint32_t virtualCornerNeighborSize = 0;
				uint64_t virtualCornerNeighbor = it ->computeNodeVirtualMorton(c,max_depth,virtualCornerNeighborSize, global.MAX_LEVEL, global.nodeface);
				if(virtualCornerNeighborSize){
					int proc = findOwner(virtualCornerNeighbor);
					procs.insert(proc);
				}
			}
		}

		set<int>::iterator pitend = procs.end();
		for(set<int>::iterator pit = procs.begin(); pit != pitend; ++pit){
			int p = *pit;
			if(p != rank){
				//TODO better reserve to avoid if
				bordersPerProc[p].push_back(distance(begin,it));
				vector<uint32_t> & bordersSingleProc = bordersPerProc[p];
				if(bordersSingleProc.capacity() - bordersSingleProc.size() < 2)
					bordersSingleProc.reserve(2*bordersSingleProc.size());
			}
		}
		count++;
	}

	//cout << rank << " border end " << endl;

	MPI_Barrier(comm);

	//	cout << rank << " pack " << endl;
	//PACK (mpi) BORDER OCTANTS IN CHAR BUFFERS WITH SIZE (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//it visits every element in bordersPerProc (one for every neighbor proc)
	//for every element it visits the border octants it contains and pack them in a new structure, sendBuffers
	//this map has an entry Class_Comm_Buffer for every proc containing the size in bytes of the buffer and the octants
	//to be sent to that proc packed in a char* buffer
	uint64_t global_index;
	uint32_t x,y,z;
	uint8_t l;
	int8_t m;
	bool info[17];
	map<int,Class_Comm_Buffer> sendBuffers;
	map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
	uint32_t pbordersOversize = 0;
	for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
		pbordersOversize += bit->second.size();
		//int buffSize = bit->second.size() * (int)ceil((double)(global.octantBytes + global.globalIndexBytes) / (double)(CHAR_BIT/8));
		int buffSize = bit->second.size() * (int)ceil((double)(global.octantBytes + global.globalIndexBytes) / (double)(CHAR_BIT/8));
		int key = bit->first;
		const vector<uint32_t> & value = bit->second;
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
		int pos = 0;
		int nofBorders = value.size();
		for(int i = 0; i < nofBorders; ++i){
			//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
			const classOctant & octant = octree.octants[value[i]];
			x = octant.getX();
			y = octant.getY();
			z = octant.getZ();
			l = octant.getLevel();
			m = octant.getMarker();
			global_index = getGlobalIdx(value[i]);
			for(int i = 0; i < 17; ++i)
				info[i] = octant.info[i];
			error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			for(int j = 0; j < 17; ++j){
				MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			}
			error_flag = MPI_Pack(&global_index,1,MPI_UINT64_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
		}
	}

	//COMMUNICATE THE SIZE OF BUFFER TO THE RECEIVERS
	//the size of every borders buffer is communicated to the right process in order to build the receive buffer
	//and stored in the recvBufferSizePerProc structure
	MPI_Request* req = new MPI_Request[sendBuffers.size()*20];
	MPI_Status* stats = new MPI_Status[sendBuffers.size()*20];
	int nReq = 0;
	map<int,int> recvBufferSizePerProc;
	map<int,Class_Comm_Buffer>::iterator sitend = sendBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		recvBufferSizePerProc[sit->first] = 0;
		error_flag = MPI_Irecv(&recvBufferSizePerProc[sit->first],1,MPI_UINT32_T,sit->first,rank,comm,&req[nReq]);
		++nReq;
	}
	map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,comm,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//COMMUNICATE THE BUFFERS TO THE RECEIVERS
	//recvBuffers structure is declared and each buffer is initialized to the right size
	//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
	//at the same time every process compute the size in bytes of all the ghost octants
	uint32_t nofBytesOverProc = 0;
	map<int,Class_Comm_Buffer> recvBuffers;
	map<int,int>::iterator ritend = recvBufferSizePerProc.end();
	for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
		recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a',comm);
	}
	nReq = 0;
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
		error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,comm,&req[nReq]);
		++nReq;
	}
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//COMPUTE GHOSTS SIZE IN BYTES
	//number of ghosts in every process is obtained through the size in bytes of the single octant
	//and ghost vector in local tree is resized
	//uint32_t nofGhosts = nofBytesOverProc / (uint32_t)(global.octantBytes + global.globalIndexBytes);
	uint32_t nofGhosts = nofBytesOverProc / (uint32_t)(global.octantBytes + global.globalIndexBytes);
	octree.size_ghosts = nofGhosts;
	octree.ghosts.clear();
	octree.ghosts.resize(nofGhosts);
	octree.globalidx_ghosts.resize(nofGhosts);

	//UNPACK BUFFERS AND BUILD GHOSTS CONTAINER OF CLASS_LOCAL_TREE
	//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked octant by octant.
	//every ghost octant is built and put in the ghost vector
	uint32_t ghostCounter = 0;
	map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
		int pos = 0;
		//			int nofGhostsPerProc = int(rrit->second.commBufferSize / (uint32_t) (global.octantBytes + global.globalIndexBytes));
		int nofGhostsPerProc = int(rrit->second.commBufferSize / (uint32_t) (global.octantBytes + global.globalIndexBytes));
		for(int i = 0; i < nofGhostsPerProc; ++i){
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&z,1,MPI_UINT32_T,comm);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
			octree.ghosts[ghostCounter] = classOctant(dim,l,x,y,z);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
			octree.ghosts[ghostCounter].setMarker(m);
			for(int j = 0; j < 17; ++j){
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,comm);
				octree.ghosts[ghostCounter].info[j] = info[j];
			}
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&global_index,1,MPI_UINT64_T,comm);
			octree.globalidx_ghosts[ghostCounter] = global_index;
			++ghostCounter;
		}
	}
	recvBuffers.clear();
	sendBuffers.clear();
	recvBufferSizePerProc.clear();

	delete [] req; req = NULL;
	delete [] stats; stats = NULL;

}

void
classParaTree::commMarker() {
	//PACK (mpi) LEVEL AND MARKER OF BORDER OCTANTS IN CHAR BUFFERS WITH SIZE (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//it visits every element in bordersPerProc (one for every neighbor proc)
	//for every element it visits the border octants it contains and pack its marker in a new structure, sendBuffers
	//this map has an entry Class_Comm_Buffer for every proc containing the size in bytes of the buffer and the octants marker
	//to be sent to that proc packed in a char* buffer
	int8_t marker;
	bool mod;
	map<int,Class_Comm_Buffer> sendBuffers;
	map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
	uint32_t pbordersOversize = 0;
	for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
		pbordersOversize += bit->second.size();
		int buffSize = bit->second.size() * (int)ceil((double)(global.markerBytes + global.boolBytes) / (double)(CHAR_BIT/8));
		int key = bit->first;
		const vector<uint32_t> & value = bit->second;
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
		int pos = 0;
		int nofBorders = value.size();
		for(int i = 0; i < nofBorders; ++i){
			//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
			const classOctant & octant = octree.octants[value[i]];
			marker = octant.getMarker();
			mod	= octant.info[15];
			error_flag = MPI_Pack(&marker,1,MPI_INT8_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
			error_flag = MPI_Pack(&mod,1,MPI::BOOL,sendBuffers[key].commBuffer,buffSize,&pos,comm);
		}
	}

	//COMMUNICATE THE SIZE OF BUFFER TO THE RECEIVERS
	//the size of every borders buffer is communicated to the right process in order to build the receive buffer
	//and stored in the recvBufferSizePerProc structure
	MPI_Request* req = new MPI_Request[sendBuffers.size()*2];
	MPI_Status* stats = new MPI_Status[sendBuffers.size()*2];
	int nReq = 0;
	map<int,int> recvBufferSizePerProc;
	map<int,Class_Comm_Buffer>::iterator sitend = sendBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		recvBufferSizePerProc[sit->first] = 0;
		error_flag = MPI_Irecv(&recvBufferSizePerProc[sit->first],1,MPI_UINT32_T,sit->first,rank,comm,&req[nReq]);
		++nReq;
	}
	map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,comm,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//COMMUNICATE THE BUFFERS TO THE RECEIVERS
	//recvBuffers structure is declared and each buffer is initialized to the right size
	//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
	//at the same time every process compute the size in bytes of all the level and marker of ghost octants
	uint32_t nofBytesOverProc = 0;
	map<int,Class_Comm_Buffer> recvBuffers;
	map<int,int>::iterator ritend = recvBufferSizePerProc.end();
	for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
		recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a',comm);
	}
	nReq = 0;
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
		error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,comm,&req[nReq]);
		++nReq;
	}
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//UNPACK BUFFERS AND BUILD GHOSTS CONTAINER OF CLASS_LOCAL_TREE
	//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked octant by octant.
	//every ghost octant is built and put in the ghost vector
	uint32_t ghostCounter = 0;
	map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
		int pos = 0;
		int nofGhostsPerProc = int(rrit->second.commBufferSize / ((uint32_t) (global.markerBytes + global.boolBytes)));
		for(int i = 0; i < nofGhostsPerProc; ++i){
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&marker,1,MPI_INT8_T,comm);
			octree.ghosts[ghostCounter].setMarker(marker);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&mod,1,MPI::BOOL,comm);
			octree.ghosts[ghostCounter].info[15] = mod;
			++ghostCounter;
		}
	}
	recvBuffers.clear();
	sendBuffers.clear();
	recvBufferSizePerProc.clear();
	delete [] req; req = NULL;
	delete [] stats; stats = NULL;

}
#endif

void
classParaTree::updateAfterCoarse(){
	mapidx.clear();
#if NOMPI==0
	if(serial){
#endif
		updateAdapt();
#if NOMPI==0
	}
	else{
		//Only if parallel
		updateAdapt();
		uint64_t lastDescMortonPre, firstDescMortonPost;
		lastDescMortonPre = (rank!=0) * partition_last_desc[rank-1];
		firstDescMortonPost = (rank<nproc-1)*partition_first_desc[rank+1] + (rank==nproc-1)*partition_last_desc[rank];
		octree.checkCoarse(lastDescMortonPre, firstDescMortonPost, mapidx);
		updateAdapt();
	}
#endif
}

void
classParaTree::updateAfterCoarse(u32vector & mapidx){
#if NOMPI==0
	if(serial){
#endif
		updateAdapt();
#if NOMPI==0
	}
	else{
		//Only if parallel
		updateAdapt();
		uint64_t lastDescMortonPre, firstDescMortonPost;
		lastDescMortonPre = (rank!=0) * partition_last_desc[rank-1];
		firstDescMortonPost = (rank<nproc-1)*partition_first_desc[rank+1] + (rank==nproc-1)*partition_last_desc[rank];
		octree.checkCoarse(lastDescMortonPre, firstDescMortonPost, mapidx);
		updateAdapt();
	}

#endif
}

void
classParaTree::balance21(bool const first){
#if NOMPI==0
	bool globalDone = true, localDone = false;
	int  iteration  = 0;

	commMarker();
	octree.preBalance21(true);

	if (first){
		log.writeLog("---------------------------------------------");
		log.writeLog(" 2:1 BALANCE (balancing Marker before Adapt)");
		log.writeLog(" ");
		log.writeLog(" Iterative procedure	");
		log.writeLog(" ");
		log.writeLog(" Iteration	:	" + to_string(static_cast<unsigned long long>(iteration)));

		commMarker();

		localDone = octree.localBalance(true);
		commMarker();
		octree.preBalance21(false);
		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);

		while(globalDone){
			iteration++;
			log.writeLog(" Iteration	:	" + to_string(static_cast<unsigned long long>(iteration)));
			commMarker();
			localDone = octree.localBalance(false);
			commMarker();
			octree.preBalance21(false);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
		}

		commMarker();
		log.writeLog(" Iteration	:	Finalizing ");
		log.writeLog(" ");
		//localDone = octree.localBalance(false);
		//commMarker();
		//octree.preBalance21(true);
		//commMarker();

		log.writeLog(" 2:1 Balancing reached ");
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");

	}
	else{

		commMarker();
		MPI_Barrier(comm);
		localDone = octree.localBalanceAll(true);
		commMarker();
		octree.preBalance21(false);
		MPI_Barrier(comm);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);

		while(globalDone){
			iteration++;
			commMarker();
			localDone = octree.localBalanceAll(false);
			commMarker();
			octree.preBalance21(false);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
		}

		commMarker();
		//			localDone = octree.localBalance(false);
		//			commMarker();
		//			octree.preBalance21(false);
		//			commMarker();

	}
#else
	bool localDone = false;
	int  iteration  = 0;

	octree.preBalance21(true);

	if (first){
		log.writeLog("---------------------------------------------");
		log.writeLog(" 2:1 BALANCE (balancing Marker before Adapt)");
		log.writeLog(" ");
		log.writeLog(" Iterative procedure	");
		log.writeLog(" ");
		log.writeLog(" Iteration	:	" + to_string(static_cast<unsigned long long>(iteration)));

		localDone = octree.localBalance(true);
		octree.preBalance21(false);

		while(localDone){
			iteration++;
			log.writeLog(" Iteration	:	" + to_string(static_cast<unsigned long long>(iteration)));
			localDone = octree.localBalance(false);
			octree.preBalance21(false);
		}

		log.writeLog(" Iteration	:	Finalizing ");
		log.writeLog(" ");
		//			localDone = octree.localBalance(false);
		//			octree.preBalance21(false);

		log.writeLog(" 2:1 Balancing reached ");
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");

	}
	else{

		localDone = octree.localBalanceAll(true);
		octree.preBalance21(false);

		while(localDone){
			iteration++;
			localDone = octree.localBalanceAll(false);
			octree.preBalance21(false);
		}

		//			localDone = octree.localBalance(false);
		//			octree.preBalance21(false);

	}

#endif /* NOMPI */
}

// =================================================================================== //
// TESTING OUTPUT METHODS												    			   //
// =================================================================================== //

/** Write the physical octree mesh in .vtu format in a user-defined file.
 * If the connectivity is not stored, the method temporary computes it.
 * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
 * \param[in] filename Seriously?....
 */
void
classParaTree::write(string filename) {

	bool clear = false;
	if (octree.connectivity.size() == 0) {
		octree.computeConnectivity();
		clear = true;
	}

	stringstream name;
	name << "s" << std::setfill('0') << std::setw(4) << nproc << "-p" << std::setfill('0') << std::setw(4) << rank << "-" << filename << ".vtu";

	ofstream out(name.str().c_str());
	if(!out.is_open()){
		stringstream ss;
		ss << filename << "*.vtu cannot be opened and it won't be written.";
		log.writeLog(ss.str());
		return;
	}
	int nofNodes = octree.nodes.size();
	int nofGhostNodes = octree.ghostsnodes.size();
	int nofOctants = octree.connectivity.size();
	int nofGhosts = octree.ghostsconnectivity.size();
	int nofAll = nofGhosts + nofOctants;
	out << "<?xml version=\"1.0\"?>" << endl
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
			<< "  <UnstructuredGrid>" << endl
			<< "    <Piece NumberOfCells=\"" << octree.connectivity.size() + octree.ghostsconnectivity.size() << "\" NumberOfPoints=\"" << octree.nodes.size() + octree.ghostsnodes.size() << "\">" << endl;
	out << "      <Points>" << endl
			<< "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
			<< "          " << std::fixed;
	for(int i = 0; i < nofNodes; i++)
	{
		for(int j = 0; j < 3; ++j){
			if (j==0) out << std::setprecision(6) << trans.mapX(octree.nodes[i][j]) << " ";
			if (j==1) out << std::setprecision(6) << trans.mapY(octree.nodes[i][j]) << " ";
			if (j==2) out << std::setprecision(6) << trans.mapZ(octree.nodes[i][j]) << " ";
		}
		if((i+1)%4==0 && i!=nofNodes-1)
			out << endl << "          ";
	}
	for(int i = 0; i < nofGhostNodes; i++)
	{
		for(int j = 0; j < 3; ++j){
			if (j==0) out << std::setprecision(6) << trans.mapX(octree.ghostsnodes[i][j]) << " ";
			if (j==1) out << std::setprecision(6) << trans.mapY(octree.ghostsnodes[i][j]) << " ";
			if (j==2) out << std::setprecision(6) << trans.mapZ(octree.ghostsnodes[i][j]) << " ";
		}
		if((i+1)%4==0 && i!=nofNodes-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "      </Points>" << endl
			<< "      <Cells>" << endl
			<< "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofOctants; i++)
	{
		for(int j = 0; j < global.nnodes; j++)
		{
			int jj = j;
			if (dim==2){
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
				}
			}
			out << octree.connectivity[i][jj] << " ";
		}
		if((i+1)%3==0 && i!=nofOctants-1)
			out << endl << "          ";
	}
	for(int i = 0; i < nofGhosts; i++)
	{
		for(int j = 0; j < global.nnodes; j++)
		{
			int jj = j;
			if (dim==2){
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
				}
			}
			out << octree.ghostsconnectivity[i][jj] + nofNodes << " ";
		}
		if((i+1)%3==0 && i!=nofGhosts-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofAll; i++)
	{
		out << (i+1)*global.nnodes << " ";
		if((i+1)%12==0 && i!=nofAll-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofAll; i++)
	{
		int type;
		type = 5 + (dim*2);
		out << type << " ";
		if((i+1)%12==0 && i!=nofAll-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "      </Cells>" << endl
			<< "    </Piece>" << endl
			<< "  </UnstructuredGrid>" << endl
			<< "</VTKFile>" << endl;


	if(rank == 0){
		name.str("");
		name << "s" << std::setfill('0') << std::setw(4) << nproc << "-" << filename << ".pvtu";
		ofstream pout(name.str().c_str());
		if(!pout.is_open()){
			stringstream ss;
			ss << filename << "*.pvtu cannot be opened and it won't be written.";
			log.writeLog(ss.str());
			return;
		}

		pout << "<?xml version=\"1.0\"?>" << endl
				<< "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
				<< "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
				<< "    <PPointData>" << endl
				<< "    </PPointData>" << endl
				<< "    <PCellData Scalars=\"\">" << endl;
		pout << "    </PCellData>" << endl
				<< "    <PPoints>" << endl
				<< "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
				<< "    </PPoints>" << endl;
		for(int i = 0; i < nproc; i++)
			pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << nproc << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
		pout << "  </PUnstructuredGrid>" << endl
				<< "</VTKFile>";

		pout.close();

	}
#if NOMPI==0
	MPI_Barrier(comm);
#endif

}

/** Write the physical octree mesh in .vtu format with data for test in a user-defined file.
 * If the connectivity is not stored, the method temporary computes it.
 * The method doesn't write the ghosts on file.
 * \param[in] filename Seriously?....
 */
void
classParaTree::writeTest(string filename, vector<double> data) {

	bool clear = false;
	if (octree.connectivity.size() == 0) {
		octree.computeConnectivity();
		clear = true;
	}

	stringstream name;
	name << "s" << std::setfill('0') << std::setw(4) << nproc << "-p" << std::setfill('0') << std::setw(4) << rank << "-" << filename << ".vtu";

	ofstream out(name.str().c_str());
	if(!out.is_open()){
		stringstream ss;
		ss << filename << "*.vtu cannot be opened and it won't be written.";
		log.writeLog(ss.str());
		return;
	}
	int nofNodes = octree.nodes.size();
	int nofOctants = octree.connectivity.size();
	int nofAll = nofOctants;
	out << "<?xml version=\"1.0\"?>" << endl
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
			<< "  <UnstructuredGrid>" << endl
			<< "    <Piece NumberOfCells=\"" << octree.connectivity.size() << "\" NumberOfPoints=\"" << octree.nodes.size() << "\">" << endl;
	out << "      <CellData Scalars=\"Data\">" << endl;
	out << "      <DataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          " << std::fixed;
	int ndata = octree.connectivity.size();
	for(int i = 0; i < ndata; i++)
	{
		out << std::setprecision(6) << data[i] << " ";
		if((i+1)%4==0 && i!=ndata-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "      </CellData>" << endl
			<< "      <Points>" << endl
			<< "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
			<< "          " << std::fixed;
	for(int i = 0; i < nofNodes; i++)
	{
		for(int j = 0; j < 3; ++j){
			if (j==0) out << std::setprecision(6) << trans.mapX(octree.nodes[i][j]) << " ";
			if (j==1) out << std::setprecision(6) << trans.mapY(octree.nodes[i][j]) << " ";
			if (j==2) out << std::setprecision(6) << trans.mapZ(octree.nodes[i][j]) << " ";
		}
		if((i+1)%4==0 && i!=nofNodes-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "      </Points>" << endl
			<< "      <Cells>" << endl
			<< "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofOctants; i++)
	{
		for(int j = 0; j < global.nnodes; j++)
		{
			int jj = j;
			if (dim==2){
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
				}
			}
			out << octree.connectivity[i][jj] << " ";
		}
		if((i+1)%3==0 && i!=nofOctants-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofAll; i++)
	{
		out << (i+1)*global.nnodes << " ";
		if((i+1)%12==0 && i!=nofAll-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
			<< "          ";
	for(int i = 0; i < nofAll; i++)
	{
		int type;
		type = 5 + (dim*2);
		out << type << " ";
		if((i+1)%12==0 && i!=nofAll-1)
			out << endl << "          ";
	}
	out << endl << "        </DataArray>" << endl
			<< "      </Cells>" << endl
			<< "    </Piece>" << endl
			<< "  </UnstructuredGrid>" << endl
			<< "</VTKFile>" << endl;


	if(rank == 0){
		name.str("");
		name << "s" << std::setfill('0') << std::setw(4) << nproc << "-" << filename << ".pvtu";
		ofstream pout(name.str().c_str());
		if(!pout.is_open()){
			stringstream ss;
			ss << filename << "*.pvtu cannot be opened and it won't be written.";
			log.writeLog(ss.str());
			return;
		}

		pout << "<?xml version=\"1.0\"?>" << endl
				<< "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
				<< "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
				<< "    <PPointData>" << endl
				<< "    </PPointData>" << endl
				<< "    <PCellData Scalars=\"Data\">" << endl
				<< "      <PDataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\"/>" << endl
				<< "    </PCellData>" << endl
				<< "    <PPoints>" << endl
				<< "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
				<< "    </PPoints>" << endl;
		for(int i = 0; i < nproc; i++)
			pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << nproc << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
		pout << "  </PUnstructuredGrid>" << endl
				<< "</VTKFile>";

		pout.close();

	}
#if NOMPI==0
	MPI_Barrier(comm);
#endif

}

// =============================================================================== //




