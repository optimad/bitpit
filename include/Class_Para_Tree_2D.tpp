/*!
 *	\date			23/apr/2014
 *	\authors		Marco Cisternino
 *	\authors		Edoardo Lombardi
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Parallel Octree Manager Class - 2D specialization
 *
 *	Para Tree is the user interface class. One user should (read can...) work only
 *	with this Class and its methods.
 *	The sizes are intended in physical domain. The transformation from the logical
 *	domain to the physical domain is defined by Class_Map<2> trans.
 *
 *	The partition of the octree is performed by following the Z-curve defined by the Morton
 *	index of the octants. By default it is a balanced partition over the number of octants for each
 *	process.
 *
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

template<>
class Class_Para_Tree<2>{
	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant<2> > 	OctantsType;
	typedef vector<uint32_t>			u32vector;
	typedef vector<double>				dvector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;
	typedef vector<vector<double>	>	dvector2D;
	typedef vector<int>					ivector;
	typedef vector<vector<int>	>		ivector2D;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	//undistributed members
	uint64_t* partition_first_desc; 			/**<Global array containing position of the first possible octant in each processor*/
	uint64_t* partition_last_desc; 				/**<Global array containing position of the last possible octant in each processor*/
	uint64_t* partition_range_globalidx;	 	/**<Global array containing global index of the last existing octant in each processor*/
	uint64_t global_num_octants;   				/**<Global number of octants in the parallel octree*/
	map<int,vector<uint32_t> > bordersPerProc;	/**<Local indices of border octants per process*/
	int nproc;									/**<Number of processes of the job*/
	uint8_t max_depth;							/**<Global max existing level in the parallel octree*/

	//distributed members
	int rank;									/**<Local rank of process*/
	Class_Local_Tree<2> octree;					/**<Local tree in each processor*/

	//auxiliary members
	int error_flag;								/**<MPI error flag*/
	bool serial;								/**<True if the octree is the same on each processor, False if the octree is distributed*/

	//map member
	Class_Map<2> trans;							/**<Transformation map from logical to physical domain*/

	//log member
	Class_Log log;								/**<Log object*/

#if NOMPI==0
	MPI_Comm comm;								/**<MPI communicator*/
#endif

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	/*! Default Constructor of Para_Tree.
	 * It builds one octant with node 0 in the Origin (0,0,0)
	 * and side of length 1
	 * \param[in] logfile The file name for the log of this object. PABLO.log is the default value*/
#if NOMPI==0
	Class_Para_Tree(string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD) : log(logfile,comm_),comm(comm_){
#else
	Class_Para_Tree(string logfile="PABLO.log") : log(logfile){
#endif
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
		log.writeLog(" Number of proc		:	" + to_string(nproc));
		log.writeLog(" Dimension		:	" + to_string(2));
		log.writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL_2D));
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
	Class_Para_Tree(double X, double Y, double Z, double L,string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD):trans(X,Y,Z,L),log(logfile,comm_),comm(comm_){
#else
	Class_Para_Tree(double X, double Y, double Z, double L,string logfile="PABLO.log"):trans(X,Y,Z,L),log(logfile){
#endif
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
		log.writeLog(" Number of proc		:	" + to_string(nproc));
		log.writeLog(" Dimension		:	" + to_string(2));
		log.writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL_2D));
		log.writeLog(" Domain Origin		:	" + to_string(X));
		log.writeLog("				" + to_string(Y));
		log.writeLog("				" + to_string(Z));
		log.writeLog(" Domain Size		:	" + to_string(L));
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
	Class_Para_Tree(double X, double Y, double Z, double L, ivector2D & XY, ivector & levels,string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD):trans(X,Y,Z,L),log(logfile,comm_),comm(comm_){
#else
	Class_Para_Tree(double X, double Y, double Z, double L, ivector2D & XY, ivector & levels,string logfile="PABLO.log"):trans(X,Y,Z,L),log(logfile){
#endif
		uint8_t lev, iface;
		uint32_t x0, y0;
		uint32_t NumOctants = XY.size();
		octree.octants.resize(NumOctants);
		for (uint32_t i=0; i<NumOctants; i++){
			lev = uint8_t(levels[i]);
			x0 = uint32_t(XY[i][0]);
			y0 = uint32_t(XY[i][1]);
			Class_Octant<2> oct(lev, x0, y0);
			if (x0 == 0){
				iface = 0;
				oct.setBound(iface);
			}
			else if (x0 == global2D.max_length - oct.getSize()){
				iface = 1;
				oct.setBound(iface);
			}
			if (y0 == 0){
				iface = 2;
				oct.setBound(iface);
			}
			else if (y0 == global2D.max_length - oct.getSize()){
				iface = 3;
				oct.setBound(iface);
			}
			octree.octants[i] = oct;
		}

		setFirstDesc();
		setLastDesc();

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
		log.writeLog(" Number of proc		:	" + to_string(nproc));
		log.writeLog(" Dimension		:	" + to_string(2));
		log.writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL_2D));
		log.writeLog(" Domain Origin		:	" + to_string(X));
		log.writeLog("				" + to_string(Y));
		log.writeLog("				" + to_string(Z));
		log.writeLog(" Domain Size		:	" + to_string(L));
		log.writeLog(" Number of octants	:	" + to_string(global_num_octants));
		log.writeLog("---------------------------------------------");
		log.writeLog(" ");
#if NOMPI==0
		MPI_Barrier(comm);
#endif
	};

	// =============================================================================== //

	~Class_Para_Tree(){
		log.writeLog("---------------------------------------------");
		log.writeLog("--------------- R.I.P. PABLO ----------------");
		log.writeLog("---------------------------------------------");
		log.writeLog("---------------------------------------------");
	};

	// =============================================================================== //
	// GET/SET METHODS --------------------------------------------------------------- //

	// Octant get/set Methods
	/*! Get the global class object
	 * \return A reference to the global class object.
	 */
	const Class_Global<2>& getGlobal(){
		return global2D;
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Pointer to target octant.
	 * \return Coordinate X of node 0.
	 */
	double getX(Class_Octant<2>* oct) {
		return trans.mapX(oct->getX());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Pointer to target octant.
	 * \return Coordinate Y of node 0.
	 */
	double getY(Class_Octant<2>* oct) {
		return trans.mapY(oct->getY());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Pointer to target octant.
	 * \return Coordinate Z of node 0.
	 */
	double getZ(Class_Octant<2>* oct) {
		return trans.mapZ(oct->getZ());
	}

	/*! Get the size of an octant, i.e. the side length.
	 * \param[in] oct Pointer to target octant.
	 * \return Size of octant.
	 */
	double getSize(Class_Octant<2>* oct) {
		return trans.mapSize(oct->getSize());
	}

	/*! Get the area of an octant (for 2D case the same value of getSize).
	 * \param[in] oct Pointer to target octant.
	 * \return Size of octant.
	 */
	double getArea(Class_Octant<2>* oct) {
		return trans.mapSize(oct->getArea());
	}

	/*! Get the volume of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Volume of octant.
	 */
	double getVolume(Class_Octant<2>* oct) {
		return trans.mapArea(oct->getVolume());
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[out] center Coordinates of the center of octant.
	 */
	void getCenter(Class_Octant<2>* oct,
			vector<double>& center) {
		vector<double> center_ = oct->getCenter();
		trans.mapCenter(center_, center);
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return center Coordinates of the center of octant.
	 */
	vector<double> getCenter(Class_Octant<2>* oct) {
		vector<double> center;
		vector<double> center_ = oct->getCenter();
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the center of a face of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface Index of the target face.
	 * \return center Coordinates of the center of the iface-th face af octant.
	 */
	vector<double> getFaceCenter(Class_Octant<2>* oct, uint8_t iface) {
		vector<double> center;
		vector<double> center_ = oct->getFaceCenter(iface);
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the center of a face of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface Index of the target face.
	 * \param[out] center Coordinates of the center of the iface-th face af octant.
	 */
	void getFaceCenter(Class_Octant<2>* oct, uint8_t iface, vector<double>& center) {
		vector<double> center_ = oct->getFaceCenter(iface);
		trans.mapCenter(center_, center);
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[out] nodes Coordinates of the nodes of octant.
	 */
	void getNodes(Class_Octant<2>* oct,
			dvector2D & nodes) {
		u32vector2D nodes_;
		oct->getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return nodes Coordinates of the nodes of octant.
	 */
	dvector2D getNodes(Class_Octant<2>* oct){
		dvector2D nodes;
		u32vector2D nodes_;
		oct->getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
		return nodes;
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \param[out] normal Coordinates of the normal of face.
	 */
	void getNormal(Class_Octant<2>* oct,
			uint8_t & iface,
			dvector & normal) {
		vector<int8_t> normal_;
		oct->getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \return normal Coordinates of the normal of face.
	 */
	dvector getNormal(Class_Octant<2>* oct,
			uint8_t & iface){
		dvector normal;
		vector<int8_t> normal_;
		oct->getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
		return normal;
	}

	/*! Get the refinement marker of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Marker of octant.
	 */
	uint8_t getMarker(Class_Octant<2>* oct){								// Get refinement/coarsening marker for idx-th octant
		return oct->getMarker();
	};

	/*! Get the level of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Level of octant.
	 */
	uint8_t getLevel(Class_Octant<2>* oct){								// Get refinement/coarsening marker for idx-th octant
		return oct->getLevel();
	};

	/*! Get the bound flag on an octant face.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a boundary face.
	 */
	bool getBound(Class_Octant<2>* oct, uint8_t iface){								// Get refinement/coarsening marker for idx-th octant
		return oct->getBound(iface);
	};

	/*! Get the pbound flag on an octant face.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] iface local index of the face.
	 * \return true if the iface face is a process boundary face.
	 */
	bool getPbound(Class_Octant<2>* oct, uint8_t iface){								// Get refinement/coarsening marker for idx-th octant
		return oct->getPbound(iface);
	};

	/*! Get the union of every bound flags on faces
	 * \param[in] oct Pointer to target octant.
	 * \return true if the octant has at least a boundary face.
	 */
	bool getBound(Class_Octant<2>* oct){
		int temp = 0;
		for(int i = 0; i < global2D.nfaces; ++i)
			temp += oct->getBound(i);
		return temp != 0;
	};

	/*! Get the union of every pbound flags on faces
	 * \param[in] oct Pointer to target octant.
	 * \return true if the octant has at least a process boundary face.
	 */
	bool getPbound(Class_Octant<2>* oct){								// Get refinement/coarsening marker for idx-th octant
		int temp = 0;
		for(int i = 0; i < global2D.nfaces; ++i)
			temp += oct->getPbound(i);
		return temp != 0;
	};

	/*! Get the balancing condition of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Has octant to be balanced?
	 */
	bool getBalance(Class_Octant<2>* oct){								// Get if balancing-blocked idx-th octant
		return !oct->getNotBalance();
	};

#if NOMPI==0
	/*! Get the nature of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Is octant ghost?
	 */
	bool getIsGhost(Class_Octant<2>* oct){
		if (serial)
			return false;
		return (findOwner(oct->computeMorton()) != rank);
	};
#endif

	/*! Get if the octant is new after refinement.
	 * \param[in] oct Pointer to target octant.
	 * \return Is octant new?
	 */
	bool getIsNewR(Class_Octant<2>* oct){
		return oct->getIsNewR();
	};

	/*! Get if the octant is new after coarsening.
	 * \param[in] oct Pointer to target octant.
	 * \return Is octant new?
	 */
	bool getIsNewC(Class_Octant<2>* oct){
		return oct->getIsNewC();
	};

	/*! Get the global index of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \return Global index of octant.
	 */
	uint64_t getGlobalIdx(Class_Octant<2>* oct){
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
	uint32_t getIdx(Class_Octant<2>* oct){
#if NOMPI==0
		if (getIsGhost(oct)){
			return octree.findGhostMorton(oct->computeMorton());
		}
#endif
		return octree.findMorton(oct->computeMorton());
	};

	/*! Set the refinement marker of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
	 */
	void setMarker(Class_Octant<2>* oct, int8_t marker){					// Set refinement/coarsening marker for idx-th octant
		oct->setMarker(marker);
	};

	/*! Set the balancing condition of an octant.
	 * \param[in] oct Pointer to target octant.
	 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
	 */
	void setBalance(Class_Octant<2>* oct, bool balance){					// Set if balancing-blocked idx-th octant
		oct->setBalance(!balance);
	};


private:
	// ------------------------------------------------------------------------------- //
	//No pointer Octants get/set Methods

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Target octant.
	 * \return Coordinate X of node 0.
	 */
	double getX(Class_Octant<2> oct) {
		return trans.mapX(oct.getX());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Target octant.
	 * \return Coordinate Y of node 0.
	 */
	double getY(Class_Octant<2> oct) {
		return trans.mapY(oct.getY());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] oct Target octant.
	 * \return Coordinate Z of node 0.
	 */
	double getZ(Class_Octant<2> oct) {
		return trans.mapZ(oct.getZ());
	}

	/*! Get the size of an octant, i.e. the side length.
	 * \param[in] oct Target octant.
	 * \return Size of octant.
	 */
	double getSize(Class_Octant<2> oct) {
		return trans.mapSize(oct.getSize());
	}

	/*! Get the area of an octant (for 2D case the same value of getSize).
	 * \param[in] oct Target octant.
	 * \return Area of octant.
	 */
	double getArea(Class_Octant<2> oct) {
		return trans.mapSize(oct.getArea());
	}

	/*! Get the volume of an octant.
	 * \param[in] oct Target octant.
	 * \return Volume of octant.
	 */
	double getVolume(Class_Octant<2> oct) {
		return trans.mapArea(oct.getVolume());
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] oct Target octant.
	 * \param[out] center Coordinates of the center of octant.
	 */
	void getCenter(Class_Octant<2> oct,
			vector<double>& center) {
		vector<double> center_ = oct.getCenter();
		trans.mapCenter(center_, center);
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] oct Target octant.
	 * \return center Coordinates of the center of octant.
	 */
	vector<double> getCenter(Class_Octant<2> oct) {
		vector<double> center;
		vector<double> center_ = oct.getCenter();
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] oct Target octant.
	 * \param[out] nodes Coordinates of the nodes of octant.
	 */
	void getNodes(Class_Octant<2> oct,
			dvector2D & nodes) {
		u32vector2D nodes_;
		oct.getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] oct Target octant.
	 * \return nodes Coordinates of the nodes of octant.
	 */
	dvector2D getNodes(Class_Octant<2> oct){
		dvector2D nodes;
		u32vector2D nodes_;
		oct.getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
		return nodes;
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] oct Target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \param[out] normal Coordinates of the normal of face.
	 */
	void getNormal(Class_Octant<2> oct,
			uint8_t & iface,
			dvector & normal) {
		vector<int8_t> normal_;
		oct.getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] oct Target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \return normal Coordinates of the normal of face.
	 */
	dvector getNormal(Class_Octant<2> oct,
			uint8_t & iface){
		dvector normal;
		vector<int8_t> normal_;
		oct.getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
		return normal;
	}

	/*! Get the refinement marker of an octant.
	 * \param[in] oct Target octant.
	 * \return Marker of octant.
	 */
	uint8_t getMarker(Class_Octant<2> oct){								// Get refinement/coarsening marker for idx-th octant
		return oct.getMarker();
	};

	/*! Get the level of an octant.
	 * \param[in] oct Target octant.
	 * \return Level of octant.
	 */
	uint8_t getLevel(Class_Octant<2> oct){								// Get refinement/coarsening marker for idx-th octant
		return oct.getLevel();
	};

	/*! Get the balancing condition of an octant.
	 * \param[in] oct Target octant.
	 * \return Has octant to be balanced?
	 */
	bool getBalance(Class_Octant<2> oct){								// Get if balancing-blocked idx-th octant
		return !oct.getNotBalance();
	};

#if NOMPI==0
	/*! Get the nature of an octant.
	 * \param[in] oct Target octant.
	 * \return Is octant ghost?
	 */
	bool getIsGhost(Class_Octant<2> oct){
		return (findOwner(oct.computeMorton()) != rank);
	};
#endif

	/*! Get the global index of an octant.
	 * \param[in] oct Target octant.
	 * \return Global index of octant.
	 */
	uint64_t getGlobalIdx(Class_Octant<2> oct){
#if NOMPI==0
		if (getIsGhost(oct)){
			uint32_t idx = octree.findGhostMorton(oct.computeMorton());
			return octree.globalidx_ghosts[idx];
		}
		else{
#endif
			uint32_t idx = octree.findMorton(oct.computeMorton());
			if (rank){
				return partition_range_globalidx[rank-1] + uint64_t(idx + 1);
			}
			else{
				return uint64_t(idx);
			};
#if NOMPI==0
		};
#endif
		return global_num_octants;
	};

	/*! Get the local index of an octant.
	 * \param[in] oct Target octant.
	 * \return Local index of octant.
	 */
	uint32_t getIdx(Class_Octant<2> oct){
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

	/*! Set the refinement marker of an octant.
	 * \param[in] oct Target octant.
	 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
	 */
	void setMarker(Class_Octant<2> oct, int8_t marker){					// Set refinement/coarsening marker for idx-th octant
		oct.setMarker(marker);
	};

	/*! Set the balancing condition of an octant.
	 * \param[in] oct Target octant.
	 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
	 */
	void setBalance(Class_Octant<2> oct, bool balance){					// Set if balancing-blocked idx-th octant
		oct.setBalance(!balance);
	};

	// ------------------------------------------------------------------------------- //
	// Index get/set Methods

public:
	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] idx Local index of target octant.
	 * \return Coordinate X of node 0.
	 */
	double getX(uint32_t idx) {
		return trans.mapX(octree.octants[idx].getX());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] idx Local index of target octant.
	 * \return Coordinate Y of node 0.
	 */
	double getY(uint32_t idx) {
		return trans.mapY(octree.octants[idx].getY());
	}

	/*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
	 * \param[in] idx Local index of target octant.
	 * \return Coordinate Z of node 0.
	 */
	double getZ(uint32_t idx) {
		return trans.mapZ(octree.octants[idx].getZ());
	}

	/*! Get the size of an octant, i.e. the side length.
	 * \param[in] idx Local index of target octant.
	 * \return Area of octant.
	 */
	double getSize(uint32_t idx) {
		return trans.mapSize(octree.octants[idx].getSize());
	}

	/*! Get the area of an octant (for 2D case the same value of getSize).
	 * \param[in] idx Local index of target octant.
	 * \return Area of octant.
	 */
	double getArea(uint32_t idx) {
		return trans.mapSize(octree.octants[idx].getArea());
	}

	/*! Get the volume of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Volume of octant.
	 */
	double getVolume(uint32_t idx) {
		return trans.mapArea(octree.octants[idx].getVolume());
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[out] center Coordinates of the center of octant.
	 */
	void getCenter(uint32_t idx,
			vector<double>& center) {
		vector<double> center_ = octree.octants[idx].getCenter();
		trans.mapCenter(center_, center);
	}

	/*! Get the coordinates of the center of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return center Coordinates of the center of octant.
	 */
	vector<double> getCenter(uint32_t idx) {
		vector<double> center;
		vector<double> center_ = octree.octants[idx].getCenter();
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the center of a face of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] iface Index of the target face.
	 * \return center Coordinates of the center of the iface-th face af octant.
	 */
	vector<double> getFaceCenter(uint32_t idx, uint8_t iface) {
		vector<double> center;
		vector<double> center_ = octree.octants[idx].getFaceCenter(iface);
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the center of a face of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] iface Index of the target face.
	 * \param[out] center Coordinates of the center of the iface-th face af octant.
	 */
	void getFaceCenter(uint32_t idx, uint8_t iface, vector<double>& center) {
		vector<double> center_ = octree.octants[idx].getFaceCenter(iface);
		trans.mapCenter(center_, center);
	}

	/*! Get the coordinates of single node of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] inode Index of the target node.
	 * \return center Coordinates of the center of the iface-th face af octant.
	 */
	vector<double> getNode(uint32_t idx, uint8_t inode) {
		vector<double> node;
		u32vector node_ = octree.octants[idx].getNode(inode);
		trans.mapNode(node_, node);
		return node;
	}

	/*! Get the coordinates of the center of a face of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] iface Index of the target face.
	 * \param[out] center Coordinates of the center of the iface-th face af octant.
	 */
	void getNode(uint32_t idx, uint8_t inode, vector<double>& node) {
		u32vector node_ = octree.octants[idx].getNode(inode);
		trans.mapNode(node_, node);
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[out] nodes Coordinates of the nodes of octant.
	 */
	void getNodes(uint32_t idx,
			dvector2D & nodes) {
		u32vector2D nodes_;
		octree.octants[idx].getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
	}

	/*! Get the coordinates of the nodes of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return nodes Coordinates of the nodes of octant.
	 */
	dvector2D getNodes(uint32_t idx){
		dvector2D nodes;
		u32vector2D nodes_;
		octree.octants[idx].getNodes(nodes_);
		trans.mapNodes(nodes_, nodes);
		return nodes;
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] Local index of target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \param[out] normal Coordinates of the normal of face.
	 */
	void getNormal(uint32_t idx,
			uint8_t & iface,
			dvector & normal) {
		vector<int8_t> normal_;
		octree.octants[idx].getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
	}

	/*! Get the normal of a face of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] iface Index of the face for normal computing.
	 * \return normal Coordinates of the normal of face.
	 */
	dvector getNormal(uint32_t idx,
			uint8_t & iface){
		dvector normal;
		vector<int8_t> normal_;
		octree.octants[idx].getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
		return normal;
	}

	/*! Get the refinement marker of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Marker of octant.
	 */
	uint8_t getMarker(uint32_t idx){							// Get refinement/coarsening marker for idx-th octant
		return octree.getMarker(idx);
	};

	/*! Get the level of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Level of octant.
	 */
	uint8_t getLevel(uint32_t idx){								// Get refinement/coarsening marker for idx-th octant
		return octree.getLevel(idx);
	};

	/*! Get the balancing condition of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Has octant to be balanced?
	 */
	bool getBalance(uint32_t idx){								// Get if balancing-blocked idx-th octant
		return !octree.getBalance(idx);
	};

#if NOMPI==0
	/*! Get the nature of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Is octant ghost?
	 */
	bool getIsGhost(uint32_t idx){
		return (findOwner(octree.octants[idx].computeMorton()) != rank);
	};
#endif

	/*! Get if the octant is new after refinement.
	 * \param[in] idx Local index of target octant.
	 * \return Is octant new?
	 */
	bool getIsNewR(uint32_t idx){
		return octree.octants[idx].getIsNewR();
	};

	/*! Get if the octant is new after coarsening.
	 * \param[in] idx Local index of target octant.
	 * \return Is octant new?
	 */
	bool getIsNewC(uint32_t idx){
		return octree.octants[idx].getIsNewC();
	};

	/*! Get the global index of an octant.
	 * \param[in] idx Local index of target octant.
	 * \return Global index of octant.
	 */
	uint64_t getGlobalIdx(uint32_t idx){
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
	uint64_t getGhostGlobalIdx(uint32_t idx){
		if (idx<octree.size_ghosts){
			return octree.globalidx_ghosts[idx];
		};
		return uint64_t(octree.size_ghosts);
	};

	/*! Set the refinement marker of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
	 */
	void setMarker(uint32_t idx, int8_t marker){					// Set refinement/coarsening marker for idx-th octant
		octree.setMarker(idx, marker);
	};

	/*! Set the balancing condition of an octant.
	 * \param[in] idx Local index of target octant.
	 * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
	 */
	void setBalance(uint32_t idx, bool balance){					// Set if balancing-blocked idx-th octant
		octree.setBalance(idx, !balance);
	};

	// ------------------------------------------------------------------------------- //
	// Local Tree get/set Methods

	/*! Get the local number of octants.
	 * \return Local number of octants.
	 */
	uint32_t getNumOctants() const{
		return octree.getNumOctants();
	};

	/*! Get the local number of ghost octants.
	 * \return Local number of ghost octants.
	 */
	uint32_t getNumGhosts() const{
		return octree.getSizeGhost();
	};

	/*! Get the local depth of octree.
	 * \return Local depth of octree.
	 */
	uint8_t getLocalMaxDepth() const{							// Get max depth reached in local tree
		return octree.getLocalMaxDepth();
	};

	/*! Get the codimension for 2:1 balancing
	 * \return Maximum codimension of the entity through which the 2:1 balance is performed.
	 */
	uint8_t getBalanceCodimension() const{
		return octree.getBalanceCodim();
	};

	/*! Set the codimension for 2:1 balancing
	 * \param[in] Maximum codimension of the entity through which the 2:1 balance is performed (1 = 2:1 balance through edges (default); 2 = 2:1 balance through nodes and edges).
	 */
	void setBalanceCodimension(uint8_t b21codim){
		octree.setBalanceCodim(b21codim);
	};


	// --------------------------------
private:

	const Class_Octant<2> &  getFirstDesc() const{
		return octree.getFirstDesc();
	};

	const Class_Octant<2> &  getLastDesc() const{
		return octree.getLastDesc();
	};

	void setFirstDesc(){
		octree.setFirstDesc();
	};

	void setLastDesc(){
		octree.setLastDesc();
	};

	Class_Octant<2>& extractOctant(uint32_t idx) {
		return octree.extractOctant(idx) ;
	};

	// --------------------------------

public:

	/** Get an octant as pointer to the target octant.
	 * \param[in] idx Local index of target octant.
	 * \return Pointer to target octant.
	 */
	Class_Octant<2>* getOctant(uint32_t idx) {
		if (idx < octree.getNumOctants()){
			return &octree.octants[idx] ;
		}
		return NULL;
	};

	/** Get a ghost octant as pointer to the target octant.
	 * \param[in] idx Local index (in ghosts structure) of target ghost octant.
	 * \return Pointer to target ghost octant.
	 */
	Class_Octant<2>* getGhostOctant(uint32_t idx) {
		if (idx < octree.getSizeGhost()){
			return &octree.ghosts[idx] ;
		}
		return NULL;
	};

	/** Finds neighbours of octant through iface in vector octants.
	 * Returns a vector (empty if iface is a bound face) with the index of neighbours
	 * in their structure (octants or ghosts) and sets isghost[i] = true if the
	 * i-th neighbour is ghost in the local tree.
	 * \param[in] idx Index of current octant
	 * \param[in] iface Index of face/edge/node passed through for neighbours finding
	 * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
	 * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
	 * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
	void findNeighbours(uint32_t idx,
			uint8_t iface,
			uint8_t codim,
			u32vector & neighbours,
			vector<bool> & isghost){

		if (codim == 1){
			octree.findNeighbours(idx, iface, neighbours, isghost);
		}
		else if (codim == 2){
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
	void findNeighbours(Class_Octant<2>* oct,
			uint8_t iface,
			uint8_t codim,
			u32vector & neighbours,
			vector<bool> & isghost){

		if (codim == 1){
			octree.findNeighbours(oct, iface, neighbours, isghost);
		}
		else if (codim == 2){
			octree.findNodeNeighbours(oct, iface, neighbours, isghost);
		}
		else {
			neighbours.clear();
			isghost.clear();
		}

	};

private:
	/** Finds neighbours of octant through iface in vector octants.
	 * Returns a vector (empty if iface is a bound face) with the index of neighbours
	 * in their structure (octants or ghosts) and sets isghost[i] = true if the
	 * i-th neighbour is ghost in the local tree.
	 * \param[in] oct Current octant
	 * \param[in] iface Index of face/edge/node passed through for neighbours finding
	 * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
	 * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
	 * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
	void findNeighbours(Class_Octant<2> oct,
			uint8_t iface,
			uint8_t codim,
			u32vector & neighbours,
			vector<bool> & isghost){

		if (codim == 1){
			octree.findNeighbours(&oct, iface, neighbours, isghost);
		}
		else if (codim == 2){
			octree.findNodeNeighbours(&oct, iface, neighbours, isghost);
		}
		else {
			neighbours.clear();
			isghost.clear();
		}
	};

public:
	//-------------------------------------------------------------------------------- //
	// Intersections get Methods

	/*! Get the local number of intersections.
	 * \return Local number of intersections.
	 */
	uint32_t getNumIntersections() {
		return octree.intersections.size();
	}

	/*! Get a pointer to target intersection.
	 * \param[in] idx Local index of intersection.
	 * \return Pointer to target intersection.
	 */
	Class_Intersection<2>* getIntersection(uint32_t idx) {
		if (idx < octree.intersections.size()){
			return &octree.intersections[idx];
		}
		return NULL;
	}

	/*! Get the level of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return Level of intersection.
	 */
	uint8_t getLevel(Class_Intersection<2>* inter) {
		if(inter->finer && inter->isghost)
			return octree.extractGhostOctant(inter->owners[inter->finer]).getLevel();
		else
			return octree.extractOctant(inter->owners[inter->finer]).getLevel();
	}

	/*! Get the finer owner octant of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return The finer octant of the owners of intersection (false/true = 0/1).
	 */
	bool getFiner(Class_Intersection<2>* inter) {
		return inter->finer;
	}

	/*! Get if an intersection is a boundary domain intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return Boundary or not boundary?.
	 */
	bool getBound(Class_Intersection<2>* inter) {
		return inter->getBound();
	}

	/*! Get if an intersection is an intersection between an internal and a ghost element.
	 * \param[in] inter Pointer to target intersection.
	 * \return Ghost or not ghost?.
	 */
	bool getIsGhost(Class_Intersection<2>* inter) {
		return inter->getIsGhost();
	}

	/*! Get if an intersection is a boundary intersection for a process.
	 * \param[in] inter Pointer to target intersection.
	 * \return Process boundary or not boundary?.
	 */
	bool getPbound(Class_Intersection<2>* inter) {
		return inter->getPbound();
	}

	/*! Get the face index of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return Face index of the first octant owner of intersection (owners[0]).
	 */
	uint8_t getFace(Class_Intersection<2>* inter) {
		return inter->iface;
	}

	/*! Get the owner octants of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return A couple of octants owners of intersection.
	 */
	u32vector getOwners(Class_Intersection<2>* inter) {
		u32vector owners(2);
		owners[0] = inter->owners[0];
		owners[1] = inter->owners[1];
		return owners;
	}

	/*! Get the size of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \return Size of intersection.
	 */
	double getSize(Class_Intersection<2>* inter) {
		uint32_t Size;
		if(inter->finer && inter->isghost)
			Size = octree.extractGhostOctant(inter->owners[inter->finer]).getSize();
		else
			Size = octree.extractOctant(inter->owners[inter->finer]).getSize();
		return trans.mapSize(Size);
	}

	/*! Get the area of an intersection (for 2D case the same value of getSize).
	 * \param[in] inter Pointer to target intersection.
	 * \return Area of intersection.
	 */
	double getArea(Class_Intersection<2>* inter) {
		uint32_t Area;
		if(inter->finer && inter->isghost)
			Area = octree.extractGhostOctant(inter->owners[1]).getArea();
		else
			Area = octree.extractOctant(inter->owners[inter->finer]).getArea();
		return trans.mapSize(Area);
	}

	/*! Get the coordinates of the center of an intersection.
	 * \param[in] inter Pointer to target intersection.
	 * \param[out] center Coordinates of the center of intersection.
	 */
	vector<double> getCenter(Class_Intersection<2>* inter){
		vector<double> center;
		Class_Octant<2> oct;
		if(inter->finer && inter->isghost)
			oct = octree.extractGhostOctant(inter->owners[inter->finer]);
		else
			oct = octree.extractOctant(inter->owners[inter->finer]);
		vector<double>  center_ = oct.getCenter();
		int sign = ( int(2*((inter->iface)%2)) - 1);
		double deplace = double (sign * int(oct.getSize())) / 2;
		center_[inter->iface/2] = uint32_t(int(center_[inter->iface/2]) + deplace);
		trans.mapCenter(center_, center);
		return center;
	}

	/*! Get the coordinates of the nodes of an intersection.
	 * \param[in] oct Pointer to target intersection.
	 * \return nodes Coordinates of the nodes of intersection.
	 */
	dvector2D getNodes(Class_Intersection<2>* inter){
		dvector2D nodes;
		Class_Octant<2> oct;
		if(inter->finer && inter->isghost)
			oct = octree.extractGhostOctant(inter->owners[inter->finer]);
		else
			oct = octree.extractOctant(inter->owners[inter->finer]);
		uint8_t iface = inter->iface;
		u32vector2D nodes_all;
		oct.getNodes(nodes_all);
		u32vector2D nodes_(global2D.nnodesperface, u32vector(3));
		for (int i=0; i<global2D.nnodesperface; i++){
			for (int j=0; j<3; j++){
				nodes_[i][j] = nodes_all[global2D.facenode[iface][i]][j];
			}
		}
		trans.mapNodesIntersection(nodes_, nodes);
		return nodes;
	}

	/*! Get the normal of an intersection.
	 * \param[in] oct Pointer to target intersection.
	 * \param[out] normal Coordinates of the normal of intersection.
	 */
	dvector getNormal(Class_Intersection<2>* inter){
		dvector normal;
		Class_Octant<2> oct;
		if(inter->finer && inter->isghost)
			oct = octree.extractGhostOctant(inter->owners[inter->finer]);
		else
			oct = octree.extractOctant(inter->owners[inter->finer]);
		uint8_t iface = inter->iface;
		vector<int8_t> normal_;
		oct.getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
		return normal;
	}

	//-------------------------------------------------------------------------------- //
	// No Pointer Intersections get Methods

private:
	double getSize(Class_Intersection<2> inter) {
		uint32_t Size;
		if(inter.finer && inter.isghost)
			Size = octree.extractGhostOctant(inter.owners[inter.finer]).getSize();
		else
			Size = octree.extractOctant(inter.owners[inter.finer]).getSize();
		return trans.mapSize(Size);
	}

	double getArea(Class_Intersection<2> inter) {
		uint32_t Area;
		if(inter.finer && inter.isghost)
			Area = octree.extractGhostOctant(inter.owners[inter.finer]).getArea();
		else
			Area = octree.extractOctant(inter.owners[inter.finer]).getArea();
		return trans.mapSize(Area);
	}

	void getCenter(Class_Intersection<2> inter,dvector & center){
		Class_Octant<2> oct;
		if(inter.finer && inter.isghost)
			oct = octree.extractGhostOctant(inter.owners[inter.finer]);
		else
			oct = octree.extractOctant(inter.owners[inter.finer]);
		vector<double>center_ = oct.getCenter();
		int sign = ( int(2*((inter.iface)%2)) - 1);
		double deplace = double (sign * int(oct.getSize())) / 2;
		center_[inter.iface/2] = uint32_t(int(center_[inter.iface/2]) + deplace);
		trans.mapCenter(center_, center);
	}

	void getNodes(Class_Intersection<2> inter,
			dvector2D & nodes) {
		Class_Octant<2> oct;
		if(inter.finer && inter.isghost)
			oct = octree.extractGhostOctant(inter.owners[inter.finer]);
		else
			oct = octree.extractOctant(inter.owners[inter.finer]);
		uint8_t iface = inter.iface;
		u32vector2D nodes_all;
		oct.getNodes(nodes_all);
		u32vector2D nodes_(global2D.nnodesperface, u32vector(3));
		for (int i=0; i<global2D.nnodesperface; i++){
			for (int j=0; j<3; j++){
				nodes_[i][j] = nodes_all[global2D.facenode[iface][i]][j];
			}
		}
		trans.mapNodesIntersection(nodes_, nodes);
	}

	void getNormal(Class_Intersection<2> inter,
			dvector & normal) {
		Class_Octant<2> oct;
		if(inter.finer && inter.isghost)
			oct = octree.extractGhostOctant(inter.owners[inter.finer]);
		else
			oct = octree.extractOctant(inter.owners[inter.finer]);
		uint8_t iface = inter.iface;
		vector<int8_t> normal_;
		oct.getNormal(iface, normal_);
		trans.mapNormals(normal_, normal);
	}

	// =============================================================================== //

public:
	/** Compute the intersection of octants (intersections of bord, of inner domain and with ghost octants).
	 */
	void computeIntersections(){
		octree.computeIntersections();
	}

	// =============================================================================== //

	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point.
	 * \return Pointer to octant owner of target point (=NULL if point outside of the domain).
	 */
	Class_Octant<2>* getPointOwner(dvector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = trans.mapX(point[0]);
		y = trans.mapY(point[1]);
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (x > global2D.max_length) || (y > global2D.max_length))
			return NULL;

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
	}

	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point.
	 * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
	 */
	uint32_t getPointOwnerIdx(dvector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner = 0;

		x = trans.mapX(point[0]);
		y = trans.mapY(point[1]);

		if ((x > global3D.max_length) || (y > global3D.max_length)
				|| (point[0] < trans.X0) || (point[1] < trans.Y0))
			return -1;


#if NOMPI==0
		if(!serial) powner = findOwner(morton);
#else
		powner = 0;
#endif
		//if ((powner!=rank) || (x > global2D.max_length) || (y > global2D.max_length))
		if ((powner!=rank) && (!serial))
			return -1;

		if (x >= global2D.max_length) x = x - 1;
		if (y >= global2D.max_length) y = y - 1;
		morton = mortonEncode_magicbits(x,y);


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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return idxtry;
		}
	}

private:
	Class_Octant<2> getPointOwner2(dvector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = trans.mapX(point[0]);
		y = trans.mapY(point[1]);

		if ((x > global3D.max_length) || (y > global3D.max_length)
				|| (point[0] < trans.X0) || (point[1] < trans.Y0))
			return -1;


		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		if (!serial) powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) && (!serial)){
			Class_Octant<2> oct0;
			return oct0;
		}

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
			return octree.octants[idxtry];
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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return octree.octants[idxtry];
		}
	}

public:
	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point in logical domain.
	 * \return Pointer to octant owner of target point (=NULL if point outside of the domain).
	 */
	Class_Octant<2>* getPointOwner(u32vector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = point[0];
		y = point[1];
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (x > global2D.max_length) || (y > global2D.max_length))
			return NULL;

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return &octree.octants[idxtry];
		}
	}

	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point in logical domain.
	 * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
	 */
	uint32_t getPointOwnerIdx(u32vector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = point[0];
		y = point[1];
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (x > global2D.max_length) || (y > global2D.max_length))
			return -1;

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return idxtry;
		}
	}

	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point in logical domain.
	 * \return Pointer to octant owner of target point (=NULL if point outside of the domain).
	 */
	Class_Octant<2>* getLogicalPointOwner(dvector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = uint32_t(point[0]);
		y = uint32_t(point[1]);
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (point[0] < 0) || (point[0] > double(global2D.max_length)) || (point[1] < 0) || (point[1] > double(global2D.max_length)))
			return NULL;

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return &octree.octants[idxtry];
		}
	}

	/** Get the octant owner of an input point.
	 * \param[in] point Coordinates of target point in logical domain.
	 * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
	 */
	uint32_t getLogicalPointOwnerIdx(dvector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = uint32_t(point[0]);
		y = uint32_t(point[1]);
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (point[0] < 0) || (point[0] > double(global2D.max_length)) || (point[1] < 0) || (point[1] > double(global2D.max_length)))
			return -1;

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return idxtry;
		}
	}

private:
	Class_Octant<2> getPointOwner2(u32vector & point){
		uint32_t noctants = octree.octants.size();
		uint32_t idxtry = noctants/2;
		uint32_t x, y;
		uint64_t morton, mortontry;
		int powner;

		x = point[0];
		y = point[1];
		morton = mortonEncode_magicbits(x,y);

#if NOMPI==0
		powner = findOwner(morton);
#else
		powner = 0;
#endif
		if ((powner!=rank) || (x > global2D.max_length) || (y > global2D.max_length)){
			Class_Octant<2> oct0;
			return oct0;
		}

		if (x == global2D.max_length) x = x - 1;
		if (y == global2D.max_length) y = y - 1;

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
			return octree.octants[idxtry];
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
						idxtry = noctants-1;
						break;
					}
				}
			}
			return octree.octants[idxtry];
		}
	}

	// =============================================================================== //
	// PARATREE METHODS ----------------------------------------------------------------------- //

#if NOMPI==0
	void computePartition(uint32_t* partition) {
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

	// =============================================================================== //

	void computePartition(uint32_t* partition, uint8_t & level_) {
		uint8_t level = uint8_t(min(int(max(int(max_depth) - int(level_), int(1))) , MAX_LEVEL_2D));
		uint32_t* partition_temp = new uint32_t[nproc];
		uint8_t* boundary_proc = new uint8_t[nproc-1];
		uint8_t dimcomm, indcomm;
		uint8_t* glbdimcomm = new uint8_t[nproc];
		uint8_t* glbindcomm = new uint8_t[nproc];

		uint32_t division_result = 0;
		uint32_t remind = 0;
		uint32_t Dh = uint32_t(pow(double(2),double(MAX_LEVEL_2D-level)));
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

	// =============================================================================== //

	void updateLoadBalance() {
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

	// =============================================================================== //

	int findOwner(const uint64_t & morton) {
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

	// =============================================================================== //

	void setPboundGhosts() {
		//BUILD BORDER OCTANT INDECES VECTOR (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
		//find local octants to be sent as ghost to the right processes
		//it visits the local octants building virtual neighbors on each octant face
		//find the owner of these virtual neighbor and build a map (process,border octants)
		//this map contains the local octants as ghosts for neighbor processes

		// NO PBORDERS !
		Class_Local_Tree<2>::OctantsType::iterator end = octree.octants.end();
		Class_Local_Tree<2>::OctantsType::iterator begin = octree.octants.begin();
		bordersPerProc.clear();
		for(Class_Local_Tree<2>::OctantsType::iterator it = begin; it != end; ++it){
			set<int> procs;
			//Virtual Face Neighbors
			for(uint8_t i = 0; i < global2D.nfaces; ++i){
				if(it->getBound(i) == false){
					uint32_t virtualNeighborsSize = 0;
					vector<uint64_t> virtualNeighbors = it->computeVirtualMorton(i,max_depth,virtualNeighborsSize);
					uint32_t maxDelta = virtualNeighborsSize/2;
					for(uint32_t j = 0; j <= maxDelta; ++j){
						int pBegin = findOwner(virtualNeighbors[j]);
						int pEnd = findOwner(virtualNeighbors[virtualNeighborsSize - 1 - j]);
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
			//Virtual Corner Neighbors
			for(uint8_t c = 0; c < global2D.nnodes; ++c){
				if(!it->getBound(global2D.nodeface[c][0]) && !it->getBound(global2D.nodeface[c][1])){
					uint32_t virtualCornerNeighborSize = 0;
					uint64_t virtualCornerNeighbor = it ->computeNodeVirtualMorton(c,max_depth,virtualCornerNeighborSize);
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
		}
		MPI_Barrier(comm);

		//PACK (mpi) BORDER OCTANTS IN CHAR BUFFERS WITH SIZE (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
		//it visits every element in bordersPerProc (one for every neighbor proc)
		//for every element it visits the border octants it contains and pack them in a new structure, sendBuffers
		//this map has an entry Class_Comm_Buffer for every proc containing the size in bytes of the buffer and the octants
		//to be sent to that proc packed in a char* buffer
		uint64_t global_index;
		uint32_t x,y;
		uint8_t l;
		int8_t m;
		bool info[12];
		map<int,Class_Comm_Buffer> sendBuffers;
		map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
		uint32_t pbordersOversize = 0;
		for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
			pbordersOversize += bit->second.size();
			int buffSize = bit->second.size() * (int)ceil((double)(global2D.octantBytes + global2D.globalIndexBytes) / (double)(CHAR_BIT/8));
			int key = bit->first;
			const vector<uint32_t> & value = bit->second;
			sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
			int pos = 0;
			int nofBorders = value.size();
			for(int i = 0; i < nofBorders; ++i){
				//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
				const Class_Octant<2> & octant = octree.octants[value[i]];
				x = octant.getX();
				y = octant.getY();
				l = octant.getLevel();
				m = octant.getMarker();
				global_index = getGlobalIdx(value[i]);
				for(int i = 0; i < 12; ++i)
					info[i] = octant.info[i];
				error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
				error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
				error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
				error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
				for(int j = 0; j < 12; ++j){
					MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[key].commBuffer,buffSize,&pos,comm);
				}
				error_flag = MPI_Pack(&global_index,1,MPI_UINT64_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
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
		uint32_t nofGhosts = nofBytesOverProc / (uint32_t)(global2D.octantBytes + global2D.globalIndexBytes);
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
			int nofGhostsPerProc = int(rrit->second.commBufferSize / (uint32_t) (global2D.octantBytes + global2D.globalIndexBytes));
			for(int i = 0; i < nofGhostsPerProc; ++i){
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
				octree.ghosts[ghostCounter] = Class_Octant<2>(l,x,y);
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
				octree.ghosts[ghostCounter].setMarker(m);
				for(int j = 0; j < 12; ++j){
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

	// =============================================================================== //

public:
	/** Distribute Load-Balancing the octants of the whole tree over
	 * the processes of the job following the Morton order.
	 * Until loadBalance is not called for the first time the mesh is serial.
	 */
	void loadBalance(){

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
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]+1));
			}

			uint32_t stride = 0;
			for(int i = 0; i < rank; ++i)
				stride += partition[i];
			Class_Local_Tree<2>::OctantsType octantsCopy = octree.octants;
			Class_Local_Tree<2>::OctantsType::const_iterator first = octantsCopy.begin() + stride;
			Class_Local_Tree<2>::OctantsType::const_iterator last = first + partition[rank];
			octree.octants.assign(first, last);
			octree.octants.shrink_to_fit();
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
			log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
			for(int ii=1; ii<nproc; ii++){
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
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

			uint32_t x,y;
			uint8_t l;
			int8_t m;
			bool info[12];
			//build send buffers from Head
			if(headSize != 0){
				for(int p = firstPredecessor; p >= 0; --p){
					if(headSize <=partition[p]){
						int buffSize = headSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						for(uint32_t i = 0; i <= (uint32_t)lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						for(uint32_t i = (uint32_t)(lh - partition[p] + 1); i <= (uint32_t)lh; ++i){
							//pack octants from lh - partition[p] to lh
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						lh -= partition[p];
						headSize = lh + 1;
					}
				}

			}
			//build send buffers from Tail
			if(tailSize != 0){
				for(int p = firstSuccessor; p < nproc; ++p){
					if(tailSize <= partition[p]){
						int buffSize = tailSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						uint32_t octantsSize = (uint32_t)octree.octants.size();
						for(uint32_t i = ft; i < octantsSize; ++i){
							//PACK octants from ft to octantsSize-1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						uint32_t endOctants = ft + partition[p] - 1;
						int pos = 0;
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						ft += partition[p];
						tailSize -= partition[p];
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
				uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8)));
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
				uint32_t nofNewPerProc = (uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8)));
				int pos = 0;
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
					octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 12; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					++newCounter;
				}
			}
			octree.octants.shrink_to_fit();

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
		log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
		}
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");

	}

	// =============================================================================== //

	/** Distribute Load-Balanced the octants of the whole tree over
	 * the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
	 * The families of octants of a desired level are retained compact on the same process.
	 * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
	 */
	void loadBalance(uint8_t & level){

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
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]+1));
			}

			uint32_t stride = 0;
			for(int i = 0; i < rank; ++i)
				stride += partition[i];
			Class_Local_Tree<2>::OctantsType octantsCopy = octree.octants;
			Class_Local_Tree<2>::OctantsType::const_iterator first = octantsCopy.begin() + stride;
			Class_Local_Tree<2>::OctantsType::const_iterator last = first + partition[rank];
			octree.octants.assign(first, last);
			octree.octants.shrink_to_fit();
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
			log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
			for(int ii=1; ii<nproc; ii++){
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
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

			uint32_t x,y;
			uint8_t l;
			int8_t m;
			bool info[12];
			//build send buffers from Head
			if(headSize != 0){
				for(int p = firstPredecessor; p >= 0; --p){
					if(headSize <=partition[p]){
						int buffSize = headSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						for(uint32_t i = 0; i <= (uint32_t)lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						for(uint32_t i = (uint32_t)(lh - partition[p] + 1); i <= (uint32_t)lh; ++i){
							//pack octants from lh - partition[p] to lh
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						lh -= partition[p];
						headSize = lh + 1;
					}
				}

			}
			//build send buffers from Tail
			if(tailSize != 0){
				for(int p = firstSuccessor; p < nproc; ++p){
					if(tailSize <= partition[p]){
						int buffSize = tailSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						int pos = 0;
						uint32_t octantsSize = (uint32_t)octree.octants.size();
						for(uint32_t i = ft; i < octantsSize; ++i){
							//PACK octants from ft to octantsSize-1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						uint32_t endOctants = ft + partition[p] - 1;
						int pos = 0;
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,comm);
							}
						}
						ft += partition[p];
						tailSize -= partition[p];
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
				uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8)));
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
				uint32_t nofNewPerProc = (uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8)));
				int pos = 0;
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,comm);
					octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 12; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					++newCounter;
				}
			}
			octree.octants.shrink_to_fit();

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
		log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
		}
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");

	}

	// =============================================================================== //

	/** Distribute Load-Balancing the octants of the whole tree and data provided by the user
	 * over the processes of the job following the Morton order.
	 * Until loadBalance is not called for the first time the mesh is serial.
	 */
	template<class Impl>
	void loadBalance(Class_Data_LB_Interface<Impl> & userData){
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
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]+1));
			}

			uint32_t stride = 0;
			for(int i = 0; i < rank; ++i)
				stride += partition[i];
			Class_Local_Tree<2>::OctantsType octantsCopy = octree.octants;
			Class_Local_Tree<2>::OctantsType::const_iterator first = octantsCopy.begin() + stride;
			Class_Local_Tree<2>::OctantsType::const_iterator last = first + partition[rank];
			octree.octants.assign(first, last);
			octree.octants.shrink_to_fit();
			first = octantsCopy.end();
			last = octantsCopy.end();

			userData.assign(stride,partition[rank]);

			//Update and build ghosts here
			updateLoadBalance();
			setPboundGhosts();
		}
		else
		{
			log.writeLog(" ");
			log.writeLog(" Initial Parallel partition : ");
			log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
			for(int ii=1; ii<nproc; ii++){
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
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
			else if(lh > octree.octants.size() - 1)
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

			uint32_t x,y;
			uint8_t l;
			int8_t m;
			bool info[12];
			//build send buffers from Head
			if(headSize != 0){
				for(int p = firstPredecessor; p >= 0; --p){
					if(headSize <=partition[p]){
						int buffSize = headSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * headSize;
						}
						else{
							for(uint32_t i = 0; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&headSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = 0; i <= lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);

							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * partition[p];
						}
						else{
							for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
							//pack octants from lh - partition[p] to lh
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						lh -= partition[p];
						headSize = lh + 1;
					}
				}

			}
			//build send buffers from Tail
			if(tailSize != 0){
				for(int p = firstSuccessor; p < nproc; ++p){
					if(tailSize <= partition[p]){
						uint32_t octantsSize = (uint32_t)octree.octants.size();
						int buffSize = tailSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * tailSize;
						}
						else{
							for(uint32_t i = ft; i <= octantsSize; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&tailSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = ft; i < octantsSize; ++i){
							//PACK octants from ft to octantsSize-1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						break;
					}
					else{
						uint32_t endOctants = ft + partition[p] - 1;
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * partition[p];
						}
						else{
							for(uint32_t i = ft; i <= endOctants; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						ft += partition[p];
						tailSize -= partition[p];
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
			//int globalRecvsBuff[globalRecvsBuffSize];
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
			}

			nReq = 0;
			for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
				error_flag = MPI_Irecv(recvBuffers[*sendit].commBuffer,recvBuffers[*sendit].commBufferSize,MPI_PACKED,*sendit,rank,comm,&req[nReq]);
				++nReq;
			}
			for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
				error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
				++nReq;
			}
			MPI_Waitall(nReq,req,stats);

			//Unpack number of octants per sender
			map<int,uint32_t> nofNewOverProcs;
			map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
			for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
				uint32_t nofNewPerProc;
				MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&nofNewPerProc,1,MPI_UINT32_T,comm);
				nofNewOverProcs[rbit->first] = nofNewPerProc;
				if(rbit->first < rank)
					nofNewHead += nofNewPerProc;
				else if(rbit->first > rank)
					nofNewTail += nofNewPerProc;
			}

			//MOVE RESIDENT TO BEGIN IN OCTANTS
			uint32_t resEnd = octree.getNumOctants() - tailOffset;
			uint32_t nofResidents = resEnd - headOffset;
			uint32_t octCounter = 0;
			for(uint32_t i = headOffset; i < resEnd; ++i){
				octree.octants[octCounter] = octree.octants[i];
				//TODO move data - DONE
				userData.move(i,octCounter);
				++octCounter;
			}
			uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
			octree.octants.resize(newCounter);
			userData.resize(newCounter);
			//MOVE RESIDENTS IN RIGHT POSITION
			uint32_t resCounter = nofNewHead + nofResidents - 1;
			for(uint32_t k = 0; k < nofResidents ; ++k){
				octree.octants[resCounter - k] = octree.octants[nofResidents - k - 1];
				//TODO move data - DON
				userData.move(nofResidents - k - 1,resCounter - k);
			}

			//UNPACK BUFFERS AND BUILD NEW OCTANTS
			newCounter = 0;
			bool jumpResident = false;

			for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
				//TODO change new octants counting, probably you have to communicate the number of news per proc
				uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&x,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&y,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&l,1,MPI_UINT8_T,comm);
					octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 12; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					//TODO Unpack data
					userData.scatter(rbit->second,newCounter);
					++newCounter;
				}
			}
			octree.octants.shrink_to_fit();
			userData.shrink();

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
		log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
		}
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");


	}

	// =============================================================================== //

	/** Distribute Load-Balanced the octants of the whole tree and data provided by the user
	 * over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
	 * The families of octants of a desired level are retained compact on the same process.
	 * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
	 */
	template<class Impl>
	void loadBalance(Class_Data_LB_Interface<Impl> & userData, uint8_t & level){

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
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]+1));
			}

			uint32_t stride = 0;
			for(int i = 0; i < rank; ++i)
				stride += partition[i];
			Class_Local_Tree<2>::OctantsType octantsCopy = octree.octants;
			Class_Local_Tree<2>::OctantsType::const_iterator first = octantsCopy.begin() + stride;
			Class_Local_Tree<2>::OctantsType::const_iterator last = first + partition[rank];
			octree.octants.assign(first, last);
			octree.octants.shrink_to_fit();
			first = octantsCopy.end();
			last = octantsCopy.end();

			userData.assign(stride,partition[rank]);

			//Update and build ghosts here
			updateLoadBalance();
			setPboundGhosts();

		}
		else
		{
			log.writeLog(" ");
			log.writeLog(" Initial Parallel partition : ");
			log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
			for(int ii=1; ii<nproc; ii++){
				log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
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
			else if(lh > octree.octants.size() - 1)
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

			uint32_t x,y;
			uint8_t l;
			int8_t m;
			bool info[12];
			//build send buffers from Head
			if(headSize != 0){
				for(int p = firstPredecessor; p >= 0; --p){
					if(headSize <=partition[p]){
						int buffSize = headSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * headSize;
						}
						else{
							for(uint32_t i = 0; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&headSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = 0; i <= lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);

							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * partition[p];
						}
						else{
							for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
							//pack octants from lh - partition[p] to lh
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						lh -= partition[p];
						headSize = lh + 1;
					}
				}

			}
			//build send buffers from Tail
			if(tailSize != 0){
				for(int p = firstSuccessor; p < nproc; ++p){
					if(tailSize <= partition[p]){
						uint32_t octantsSize = (uint32_t)octree.octants.size();
						int buffSize = tailSize * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * tailSize;
						}
						else{
							for(uint32_t i = ft; i <= octantsSize; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&tailSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = ft; i < octantsSize; ++i){
							//PACK octants from ft to octantsSize-1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						break;
					}
					else{
						uint32_t endOctants = ft + partition[p] - 1;
						int buffSize = partition[p] * (int)ceil((double)global2D.octantBytes / (double)(CHAR_BIT/8));
						//TODO loop over head octants and add data size to buffer size - DONE
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * partition[p];
						}
						else{
							for(uint32_t i = ft; i <= endOctants; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							const Class_Octant<2> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int i = 0; i < 12; ++i)
								info[i] = octant.info[i];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 12; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							//TODO call gather to pack user data - DONE
							userData.gather(sendBuffers[p],i);
						}
						ft += partition[p];
						tailSize -= partition[p];
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
			}

			nReq = 0;
			for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
				error_flag = MPI_Irecv(recvBuffers[*sendit].commBuffer,recvBuffers[*sendit].commBufferSize,MPI_PACKED,*sendit,rank,comm,&req[nReq]);
				++nReq;
			}
			for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
				error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
				++nReq;
			}
			MPI_Waitall(nReq,req,stats);

			//Unpack number of octants per sender
			map<int,uint32_t> nofNewOverProcs;
			map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
			for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
				uint32_t nofNewPerProc;
				MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&nofNewPerProc,1,MPI_UINT32_T,comm);
				nofNewOverProcs[rbit->first] = nofNewPerProc;
				if(rbit->first < rank)
					nofNewHead += nofNewPerProc;
				else if(rbit->first > rank)
					nofNewTail += nofNewPerProc;
			}

			//MOVE RESIDENT TO BEGIN IN OCTANTS
			uint32_t resEnd = octree.getNumOctants() - tailOffset;
			uint32_t nofResidents = resEnd - headOffset;
			uint32_t octCounter = 0;
			for(uint32_t i = headOffset; i < resEnd; ++i){
				octree.octants[octCounter] = octree.octants[i];
				//TODO move data - DONE
				userData.move(i,octCounter);
				++octCounter;
			}
			uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
			octree.octants.resize(newCounter);
			userData.resize(newCounter);
			//MOVE RESIDENTS IN RIGHT POSITION
			uint32_t resCounter = nofNewHead + nofResidents - 1;
			for(uint32_t k = 0; k < nofResidents ; ++k){
				octree.octants[resCounter - k] = octree.octants[nofResidents - k - 1];
				//TODO move data - DON
				userData.move(nofResidents - k - 1,resCounter - k);
			}

			//UNPACK BUFFERS AND BUILD NEW OCTANTS
			newCounter = 0;
			bool jumpResident = false;

			for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
				//TODO change new octants counting, probably you have to communicate the number of news per proc
				uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&x,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&y,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&l,1,MPI_UINT8_T,comm);
					octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 12; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					//TODO Unpack data
					userData.scatter(rbit->second,newCounter);
					++newCounter;
				}
			}
			octree.octants.shrink_to_fit();
			userData.shrink();

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
		log.writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
		for(int ii=1; ii<nproc; ii++){
			log.writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
		}
		log.writeLog(" ");
		log.writeLog("---------------------------------------------");


	}
#endif /* NOMPI */
	// =============================================================================== //

private:
	void updateAdapt() {
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

	// =============================================================================== //

	void updateAfterCoarse(){
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
			octree.checkCoarse(lastDescMortonPre, firstDescMortonPost);
			updateAdapt();
		}
#endif
	}

	// =============================================================================== //

	void updateAfterCoarse(u32vector & mapidx){
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

	//TODO Update after coarse with intersections

	// =============================================================================== //

#if NOMPI==0
	void commMarker() {
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
			int buffSize = bit->second.size() * (int)ceil((double)(global2D.markerBytes + global2D.boolBytes) / (double)(CHAR_BIT/8));
			int key = bit->first;
			const vector<uint32_t> & value = bit->second;
			sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
			int pos = 0;
			int nofBorders = value.size();
			for(int i = 0; i < nofBorders; ++i){
				//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
				const Class_Octant<2> & octant = octree.octants[value[i]];
				marker = octant.getMarker();
				mod	= octant.info[11];
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
			int nofGhostsPerProc = int(rrit->second.commBufferSize / ((uint32_t) (global2D.markerBytes + global2D.boolBytes)));
			for(int i = 0; i < nofGhostsPerProc; ++i){
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&marker,1,MPI_INT8_T,comm);
				octree.ghosts[ghostCounter].setMarker(marker);
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&mod,1,MPI::BOOL,comm);
				octree.ghosts[ghostCounter].info[11] = mod;
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

	//==============================================================

	void balance21(bool const first){
#if NOMPI==0
		bool globalDone = true, localDone = false;
		int  iteration  = 0;

		if (first){
			log.writeLog("---------------------------------------------");
			log.writeLog(" 2:1 BALANCE (balancing Marker before Adapt)");
			log.writeLog(" ");
			log.writeLog(" Iterative procedure	");
			log.writeLog(" ");
			log.writeLog(" Iteration	:	" + to_string(iteration));

			commMarker();

			localDone = octree.localBalance(true);
			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);

			while(globalDone){
				iteration++;
				log.writeLog(" Iteration	:	" + to_string(iteration));
				commMarker();
				localDone = octree.localBalance(false);
				error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			}

			commMarker();
			log.writeLog(" Iteration	:	Finalizing ");
			log.writeLog(" ");
			localDone = octree.localBalance(false);
			commMarker();

			log.writeLog(" 2:1 Balancing reached ");
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");

		}
		else{

			commMarker();
			MPI_Barrier(comm);
			localDone = octree.localBalanceAll(true);
			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);

			while(globalDone){
				iteration++;
				commMarker();
				localDone = octree.localBalanceAll(false);
				error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			}

			commMarker();
			localDone = octree.localBalance(false);
			commMarker();

		}
#else
		bool localDone = false;
		int  iteration  = 0;


		if (first){
			log.writeLog("---------------------------------------------");
			log.writeLog(" 2:1 BALANCE (balancing Marker before Adapt)");
			log.writeLog(" ");
			log.writeLog(" Iterative procedure	");
			log.writeLog(" ");
			log.writeLog(" Iteration	:	" + to_string(iteration));


			localDone = octree.localBalance(true);

			while(localDone){
				iteration++;
				log.writeLog(" Iteration	:	" + to_string(iteration));
				localDone = octree.localBalance(false);
			}

			log.writeLog(" Iteration	:	Finalizing ");
			log.writeLog(" ");
			localDone = octree.localBalance(false);

			log.writeLog(" 2:1 Balancing reached ");
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");

		}
		else{

			localDone = octree.localBalanceAll(true);

			while(localDone){
				iteration++;
				localDone = octree.localBalanceAll(false);
			}

			localDone = octree.localBalance(false);

		}

#endif /* NOMPI */
	}

	// =============================================================================== //

public:
	/** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
	 */
	bool adapt() {
		bool globalDone = false, localDone = false, cDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
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
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Refine
			while(octree.refine());

			if (octree.getNumOctants() > nocts)
				localDone = true;
			log.writeLog(" Number of octants after Refine	:	" + to_string(octree.getNumOctants()));
			nocts = octree.getNumOctants();
			updateAdapt();

			// Coarse
			while(octree.coarse());
			updateAfterCoarse();
			balance21(false);
			while(octree.refine());
			updateAdapt();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			log.writeLog(" Number of octants after Coarse	:	" + to_string(nocts));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Refine
			while(octree.refine());
			if (octree.getNumOctants() > nocts)
				localDone = true;
			updateAdapt();
			//setPboundGhosts();
			log.writeLog(" Number of octants after Refine	:	" + to_string(global_num_octants));
			nocts = octree.getNumOctants();

			// Coarse
			while(octree.coarse());
			updateAfterCoarse();
			setPboundGhosts();
			balance21(false);
			while(octree.refine());
			updateAdapt();
			setPboundGhosts();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			log.writeLog(" Number of octants after Coarse	:	" + to_string(global_num_octants));
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");
		}
		return globalDone;
#else
		return localDone;
#endif
	}

	// =============================================================================== //

	/** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
	 * Track the changes in structure octant by a mapper.
	 * \param[out] mapidx Mapper from new octants to old octants.
	 * mapidx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
	 * if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
	 * if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 */
	bool adapt(u32vector & mapidx) {

		bool globalDone = false, localDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
		}

		// mapidx init
		mapidx.clear();
		mapidx.resize(nocts);
		mapidx.shrink_to_fit();
		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
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
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Refine
			while(octree.refine(mapidx));

			if (octree.getNumOctants() > nocts)
				localDone = true;
			log.writeLog(" Number of octants after Refine	:	" + to_string(octree.getNumOctants()));
			nocts = octree.getNumOctants();
			updateAdapt();

			// Coarse
			while(octree.coarse(mapidx));
			updateAfterCoarse(mapidx);
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			log.writeLog(" Number of octants after Coarse	:	" + to_string(nocts));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Refine
			while(octree.refine(mapidx));
			if (octree.getNumOctants() > nocts)
				localDone = true;
			updateAdapt();
			log.writeLog(" Number of octants after Refine	:	" + to_string(global_num_octants));
			nocts = octree.getNumOctants();

			// Coarse
			while(octree.coarse(mapidx));
			updateAfterCoarse(mapidx);
			setPboundGhosts();
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
			setPboundGhosts();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			log.writeLog(" Number of octants after Coarse	:	" + to_string(global_num_octants));
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");
		}
		return globalDone;
#else
		return localDone;
#endif
	}

	// =============================================================================== //

	/** Adapt the octree mesh refining all the octants by one level.
	 */
	bool adaptGlobalRefine() {
		bool globalDone = false, localDone = false, cDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
		}
#if NOMPI==0
		if(serial){
#endif
			log.writeLog("---------------------------------------------");
			log.writeLog(" ADAPT (GlobalRefine)");
			log.writeLog(" ");

			log.writeLog(" ");
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Refine
			while(octree.globalRefine());

			if (octree.getNumOctants() > nocts)
				localDone = true;
			log.writeLog(" Number of octants after Refine	:	" + to_string(octree.getNumOctants()));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Refine
			while(octree.globalRefine());
			if (octree.getNumOctants() > nocts)
				localDone = true;
			updateAdapt();
			setPboundGhosts();
			log.writeLog(" Number of octants after Refine	:	" + to_string(global_num_octants));
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

	// =============================================================================== //

	/** Adapt the octree mesh refining all the octants by one level.
	 * Track the changes in structure octant by a mapper.
	 * \param[out] mapidx Mapper from new octants to old octants.
	 * mapidx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
	 * if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
	 * if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 */
	bool adaptGlobalRefine(u32vector & mapidx) {

		bool globalDone = false, localDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
		}

		// mapidx init
		mapidx.clear();
		mapidx.resize(nocts);
		mapidx.shrink_to_fit();
		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
		}
#if NOMPI==0
		if(serial){
#endif
			log.writeLog("---------------------------------------------");
			log.writeLog(" ADAPT (Global Refine)");
			log.writeLog(" ");

			log.writeLog(" ");
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Refine
			while(octree.globalRefine(mapidx));

			if (octree.getNumOctants() > nocts)
				localDone = true;
			log.writeLog(" Number of octants after Refine	:	" + to_string(octree.getNumOctants()));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Refine
			while(octree.globalRefine(mapidx));
			if (octree.getNumOctants() > nocts)
				localDone = true;
			updateAdapt();
			setPboundGhosts();
			log.writeLog(" Number of octants after Refine	:	" + to_string(global_num_octants));
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

	// =============================================================================== //

	/** Adapt the octree mesh coarsening all the octants by one level.
	 */
	bool adaptGlobalCoarse() {
		bool globalDone = false, localDone = false, cDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
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
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Coarse
			while(octree.globalCoarse());
			updateAfterCoarse();
			balance21(false);
			while(octree.refine());
			updateAdapt();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			log.writeLog(" Number of octants after Coarse	:	" + to_string(nocts));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Coarse
			while(octree.globalCoarse());
			updateAfterCoarse();
			setPboundGhosts();
			balance21(false);
			while(octree.refine());
			updateAdapt();
			setPboundGhosts();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			log.writeLog(" Number of octants after Coarse	:	" + to_string(global_num_octants));
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");
		}
		return globalDone;
#else
		return localDone;
#endif
	}

	// =============================================================================== //

	/** Adapt the octree mesh coarsening all the octants by one level.
	 * Track the changes in structure octant by a mapper.
	 * \param[out] mapidx Mapper from new octants to old octants.
	 * mapidx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
	 * if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
	 * if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 */
	bool adaptGlobalCoarse(u32vector & mapidx) {

		bool globalDone = false, localDone = false;
		uint32_t nocts = octree.getNumOctants();
		vector<Class_Octant<2> >::iterator iter, iterend = octree.octants.end();

		for (iter = octree.octants.begin(); iter != iterend; iter++){
			iter->info[8] = false;
			iter->info[9] = false;
			iter->info[11] = false;
		}

		// mapidx init
		mapidx.clear();
		mapidx.resize(nocts);
		mapidx.shrink_to_fit();
		for (uint32_t i=0; i<nocts; i++){
			mapidx[i] = i;
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
			log.writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));

			// Coarse
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse(mapidx);
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			log.writeLog(" Number of octants after Coarse	:	" + to_string(nocts));
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
			log.writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));

			// Coarse
			while(octree.globalCoarse(mapidx));
			updateAfterCoarse(mapidx);
			setPboundGhosts();
			balance21(false);
			while(octree.refine(mapidx));
			updateAdapt();
			setPboundGhosts();
			if (octree.getNumOctants() < nocts){
				localDone = true;
			}
			nocts = octree.getNumOctants();

			MPI_Barrier(comm);
			error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,comm);
			log.writeLog(" Number of octants after Coarse	:	" + to_string(global_num_octants));
			log.writeLog(" ");
			log.writeLog("---------------------------------------------");
		}
		return globalDone;
#else
		return localDone;
#endif
	}

#if NOMPI==0
	// =============================================================================== //

	/** Communicate data provided by the user between the processes.
	 */
	template<class Impl>
	void communicate(Class_Data_Comm_Interface<Impl> & userData){
		//BUILD SEND BUFFERS
		map<int,Class_Comm_Buffer> sendBuffers;
		size_t fixedDataSize = userData.fixedSize();
		map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
		map<int,vector<uint32_t> >::iterator bitbegin = bordersPerProc.begin();
		for(map<int,vector<uint32_t> >::iterator bit = bitbegin; bit != bitend; ++bit){
			const int & key = bit->first;
			const vector<uint32_t> & pborders = bit->second;
			size_t buffSize = 0;
			size_t nofPbordersPerProc = pborders.size();
			if(fixedDataSize != 0){
				buffSize = fixedDataSize*nofPbordersPerProc;
			}
			else{
				for(size_t i = 0; i < nofPbordersPerProc; ++i){
					buffSize += userData.size(pborders[i]);
				}
			}
			//enlarge buffer to store number of pborders from this proc
			buffSize += sizeof(int);
			//build buffer for this proc
			sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
			//store number of pborders from this proc at the begining
			MPI_Pack(&nofPbordersPerProc,1,MPI_INT,sendBuffers[key].commBuffer,sendBuffers[key].commBufferSize,&sendBuffers[key].pos,comm);

			//WRITE SEND BUFFERS
			for(size_t j = 0; j < nofPbordersPerProc; ++j){
				userData.gather(sendBuffers[key],pborders[j]);
			}
		}

		//Communicate Buffers Size
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

		//Communicate Buffers
		map<int,Class_Comm_Buffer> recvBuffers;
		map<int,int>::iterator ritend = recvBufferSizePerProc.end();
		for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
			recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a',comm);
		}
		nReq = 0;
		for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
			error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,comm,&req[nReq]);
			++nReq;
		}
		for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//READ RECEIVE BUFFERS
		int ghostOffset = 0;
		map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
		map<int,Class_Comm_Buffer>::iterator rbitbegin = recvBuffers.begin();
		for(map<int,Class_Comm_Buffer>::iterator rbit = rbitbegin; rbit != rbitend; ++rbit){
			int nofGhostFromThisProc = 0;
			MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&nofGhostFromThisProc,1,MPI_INT,comm);
			for(int k = 0; k < nofGhostFromThisProc; ++k){
				userData.scatter(rbit->second, k+ghostOffset);
			}
			ghostOffset += nofGhostFromThisProc;
		}
		delete [] req; req = NULL;
		delete [] stats; stats = NULL;

	};
#endif /* NOMPI */
	// =============================================================================== //

	/** Compute the connectivity of octants and store the coordinates of nodes.
	 */
	void computeConnectivity() {
		octree.computeConnectivity();
	}

	// =================================================================================== //

	/** Clear the connectivity of octants.
	 */
	void clearConnectivity() {
		octree.clearConnectivity();
	}

	// =================================================================================== //

	/** Update the connectivity of octants.
	 */
	void updateConnectivity() {
		octree.updateConnectivity();
	}

	// =================================================================================== //

	/** Compute the connectivity of ghost octants and store the coordinates of nodes.
	 */
	void computeGhostsConnectivity() {
		octree.computeghostsConnectivity();
	}

	// =================================================================================== //

	/** Clear the connectivity of ghost octants.
	 */
	void clearGhostsConnectivity() {
		octree.clearghostsConnectivity();
	}

	// =================================================================================== //

	/** Update the connectivity of ghost octants.
	 */
	void updateGhostsConnectivity() {
		octree.updateghostsConnectivity();
	}

	// =================================================================================== //

	/** Get the local number of nodes.
	 */
	uint32_t getNumNodes() {
		return octree.nodes.size();
	}

	// =============================================================================== //

	/** Get the connectivity of the octants
	 * \return connectivity Matrix of noctants*4 with the connectivity of each octant (4 indices of nodes).
	 */
	const u32vector2D & getConnectivity(){
		return octree.connectivity;
	}

	// =============================================================================== //

	/** Get the connectivity of the ghost octants
	 * \return connectivity Matrix of nghostoctants*4 with the connectivity of each octant (4 indices of nodes).
	 */
	const u32vector2D & getGhostConnectivity(){
		return octree.ghostsconnectivity;
	}

	// =============================================================================== //
	/** Get the local connectivity of an octant
	 * \param[in] idx Local index of octant
	 * \return connectivity Connectivity of the octant (4 indices of nodes).
	 */
	u32vector getOctantConnectivity(uint32_t idx){
		return octree.connectivity[idx];
	}

	// =============================================================================== //
	/** Get the local connectivity of an octant
	 * \param[in] oct Pointer to an octant
	 * \return connectivity Connectivity of the octant (4 indices of nodes).
	 */
	u32vector getOctantConnectivity(Class_Octant<2>* oct){
		return octree.connectivity[getIdx(oct)];
	}

	// =============================================================================== //

	/** Get the local connectivity of a ghost octant
	 * \param[in] idx Local index of ghost octant
	 * \return connectivity Connectivity of the ghost octant (4 indices of nodes).
	 */
	u32vector getGhostOctantConnectivity(uint32_t idx){
		return octree.ghostsconnectivity[idx];
	}

	// =============================================================================== //

	/** Get the local connectivity of a ghost octant
	 * \param[in] oct Pointer to a ghost octant
	 * \return connectivity Connectivity of the ghost octant (4 indices of nodes).
	 */
	u32vector getGhostOctantConnectivity(Class_Octant<2>* oct){
		return octree.ghostsconnectivity[getIdx(oct)];
	}

	// =============================================================================== //

	/** Get the logical coordinates of the nodes
	 * \return nodes Matrix of nnodes*3 with the coordinates of the nodes.
	 */
	const u32vector2D & getNodes(){
		return octree.nodes;
	}

	// =============================================================================== //

	/** Get the logical coordinates of a node
	 * \param[in] inode Local index of node
	 * \return nodes Vector with the coordinates of the node.
	 */
	u32vector getNodeLogicalCoordinates(uint32_t inode){
		return octree.nodes[inode];
	}

	// =============================================================================== //

	/** Get the logical coordinates of the ghost nodes
	 * \return nodes Matrix of nghostnodes*3 with the coordinates of the nodes.
	 */
	const u32vector2D & getGhostNodes(){
		return octree.ghostsnodes;
	}

	// =============================================================================== //

	/** Get the physical coordinates of a node
	 * \param[in] inode Local index of node
	 * \return nodes Vector with the coordinates of the node.
	 */
	dvector getNodeCoordinates(uint32_t inode){
		vector<double> coords(3,0);
		coords[0] = trans.mapX(octree.nodes[inode][0]);
		coords[1] = trans.mapY(octree.nodes[inode][1]);
		coords[2] = trans.mapZ(octree.nodes[inode][2]);
		return coords;
	}

	// =============================================================================== //

	/** Get the logical coordinates of a ghost node
	 * \param[in] inode Local index of node
	 * \return nodes Vector with the coordinates of the node.
	 */
	u32vector getGhostNodeLogicalCoordinates(uint32_t inode){
		return octree.ghostsnodes[inode];
	}

	// =============================================================================== //

	/** Get the physical coordinates of a ghost node
	 * \param[in] inode Local index of node
	 * \return nodes Vector with the coordinates of the node.
	 */
	dvector getGhostNodeCoordinates(uint32_t inode){
		vector<double> coords(3,0);
		coords[0] = trans.mapX(octree.ghostsnodes[inode][0]);
		coords[1] = trans.mapY(octree.ghostsnodes[inode][1]);
		coords[2] = trans.mapZ(octree.ghostsnodes[inode][2]);
		return coords;
	}

	// =============================================================================== //

	/** Map the elements of the actual octree mesh to the elements of another one.
	 * If the connectivity is not stored, the method temporary computes it.
	 * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
	 * \param[in] ptree Second octree. The map goes from the firs one to the second one.
	 * \param[out] Map between octrees. Each i-th pair gives the first index and the last one of the elements
	 * of the second octree which lie in the i-th element of the first octree. If the indices are equal, then
	 * the element of the second octree is of the same level or lower (bigger size).  Each second i-th pair gives
	 * the indices of the two processes owning the two quadrants defined by the first pair.
	 */
	vector<pair<pair<uint32_t, uint32_t>, pair<int, int> > > mapPablos(Class_Para_Tree<2> & ptree){
		//TODO DO IT WITH ITERATORS
		vector<pair<pair<uint32_t, uint32_t>, pair<int, int> > > mapper;
		uint64_t morton2 = 0, morton1 = 0, mortonlastdesc = 0, mortonfirstdesc = 0;
		uint32_t idx1 = 0, idx2 = 0;
		uint32_t nocts = octree.getNumOctants();
		uint32_t nocts2 = ptree.octree.getNumOctants();
		int owner;
		mapper.resize(nocts);
#if NOMPI==0
		if (ptree.serial){
#endif
			for (uint32_t i=0; i<nocts; i++){
				mapper[i].first.first = idx1;
				mapper[i].first.second = idx2;
				mapper[i].second.first = rank;
				mapper[i].second.second = rank;
				mortonfirstdesc = octree.octants[i].computeMorton();
				mortonlastdesc = octree.octants[i].buildLastDesc().computeMorton();
				while(morton1 <= mortonfirstdesc && idx1 < nocts2){
					mapper[i].first.first = idx1;
					idx1++;
					if (idx1 < nocts2)
						morton1 = ptree.getOctant(idx1)->computeMorton();
				}
				if(idx1 > 0){
					idx1--;
					morton1 = ptree.getOctant(idx1)->computeMorton();
				}
				while(morton2 <= mortonlastdesc && idx2 < nocts2){
					mapper[i].first.second = idx2;
					idx2++;
					if (idx2 < nocts2)
						morton2 = ptree.getOctant(idx2)->computeMorton();
				}
				if (idx2 > 0){
					idx2--;
					morton2 = ptree.getOctant(idx2)->computeMorton();
				}
			}
#if NOMPI==0
		}
		else{
			map<int,vector<uint64_t> > FirstMortonperproc, SecondMortonperproc;
			map<int,vector<uint64_t> > FirstMortonReceived, SecondMortonReceived;
			map<int,vector<uint32_t> > FirstIndexperproc, SecondIndexperproc;
			map<int,vector<uint32_t> > FirstLocalIndex, SecondLocalIndex;
			idx1 = 0;
			morton1 = 0;
			idx2 = 0;
			morton2 = 0;
			for (uint32_t i=0; i<nocts; i++){
				mortonfirstdesc = octree.octants[i].computeMorton();
				owner = ptree.findOwner(mortonfirstdesc);
				if (rank == owner){
					mapper[i].second.first = rank;
					while(morton1 <= mortonfirstdesc && idx1 < nocts2){
						mapper[i].first.first = idx1;
						idx1++;
						if (idx1 < nocts2)
							morton1 = ptree.getOctant(idx1)->computeMorton();
					}
					if(idx1 > 0){
						idx1--;
						morton1 = ptree.getOctant(idx1)->computeMorton();
					}
					mortonlastdesc = octree.octants[i].buildLastDesc().computeMorton();
					owner = ptree.findOwner(mortonlastdesc);
					if (rank == owner){
						mapper[i].second.second = rank;
						mapper[i].first.second = idx2;
						while(morton2 <= mortonlastdesc && idx2 < nocts2){
							mapper[i].first.second = idx2;
							idx2++;
							if (idx2 < nocts2)
								morton2 = ptree.getOctant(idx2)->computeMorton();
						}
						if(idx2 > 0){
							idx2--;
							morton2 = ptree.getOctant(idx2)->computeMorton();
						}
					}
					else{
						mapper[i].second.second = owner;
						SecondMortonperproc[owner].push_back(mortonfirstdesc);
						SecondLocalIndex[owner].push_back(i);
					}
				}
				else{
					mapper[i].second.first = owner;
					FirstMortonperproc[owner].push_back(mortonfirstdesc);
					FirstLocalIndex[owner].push_back(i);
					mortonlastdesc = octree.octants[i].buildLastDesc().computeMorton();
					owner = ptree.findOwner(mortonlastdesc);
					if (rank == owner){
						mapper[i].second.second = rank;
						mapper[i].first.second = idx2;
						while(morton2 <= mortonlastdesc && idx2 < nocts2){
							mapper[i].first.second = idx2;
							idx2++;
							if (idx2 < nocts2)
								morton2 = ptree.getOctant(idx2)->computeMorton();
						}
						if(idx2 > 0){
							idx2--;
							morton2 = ptree.getOctant(idx2)->computeMorton();
						}
					}
					else{
						mapper[i].second.second = owner;
						SecondMortonperproc[owner].push_back(mortonfirstdesc);
						SecondLocalIndex[owner].push_back(i);
					}
				}
			}

			MPI_Barrier(comm);


			for(int iproc=0; iproc<nproc; iproc++){
				FirstMortonperproc[iproc].push_back(-1);
				SecondMortonperproc[iproc].push_back(-1);
			}

			{

				//COMM FIRST MORTON PER PROC
				map<int,Class_Comm_Buffer> sendBuffers;
				map<int,vector<uint64_t> >::iterator bitend = FirstMortonperproc.end();
				for(map<int,vector<uint64_t> >::iterator bit = FirstMortonperproc.begin(); bit != bitend; ++bit){
					int buffSize = bit->second.size() * (int)ceil((double)(sizeof(uint64_t)) / (double)(CHAR_BIT/8));
					int key = bit->first;
					vector<uint64_t> & value = bit->second;
					sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					int nofMortons = value.size();
					for(int i = 0; i < nofMortons; ++i){
						//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
						uint64_t Morton = value[i];
						error_flag = MPI_Pack(&Morton,1,MPI_UINT64_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
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
				//at the same time every process compute the size in bytes
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

				//UNPACK BUFFERS AND BUILD CONTAINER OF RECEIVED MORTON
				//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked.
				//every Morton is built and put in the MorontReceived vector
				uint32_t Mortoncounter = 0;
				uint64_t Morton = 0;
				map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
				for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
					int pos = 0;
					int nofMortonPerProc = int(rrit->second.commBufferSize / (uint32_t) (sizeof(uint64_t)));
					for(int i = 0; i < nofMortonPerProc-1; ++i){
						error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&Morton,1,MPI_UINT64_T,comm);
						FirstMortonReceived[rrit->first].push_back(Morton);
						++Mortoncounter;
					}
				}

				recvBuffers.clear();
				sendBuffers.clear();
				recvBufferSizePerProc.clear();
				delete [] req; req = NULL;
				delete [] stats; stats = NULL;

			}

			{
				//COMM SECOND MORTON PER PROC
				map<int,Class_Comm_Buffer> sendBuffers;
				map<int,vector<uint64_t> >::iterator bitend = SecondMortonperproc.end();
				uint32_t pbordersOversize = 0;
				for(map<int,vector<uint64_t> >::iterator bit = SecondMortonperproc.begin(); bit != bitend; ++bit){
					pbordersOversize += bit->second.size();
					int buffSize = bit->second.size() * (int)ceil((double)(sizeof(uint64_t)) / (double)(CHAR_BIT/8));
					int key = bit->first;
					vector<uint64_t> & value = bit->second;
					sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					int nofMortons = value.size();
					for(int i = 0; i < nofMortons; ++i){
						//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
						uint64_t Morton = value[i];
						error_flag = MPI_Pack(&Morton,1,MPI_UINT64_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
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
				//at the same time every process compute the size in bytes
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

				//UNPACK BUFFERS AND BUILD CONTAINER OF RECEIVED MORTON
				//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked.
				//every Morton is built and put in the MorontReceived vector
				uint32_t Mortoncounter = 0;
				uint64_t Morton = 0;
				map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
				for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
					int pos = 0;
					int nofMortonPerProc = int(rrit->second.commBufferSize / (uint32_t) (sizeof(uint64_t)));
					for(int i = 0; i < nofMortonPerProc-1; ++i){
						error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&Morton,1,MPI_UINT64_T,comm);
						SecondMortonReceived[rrit->first].push_back(Morton);
						++Mortoncounter;
					}
				}
				recvBuffers.clear();
				sendBuffers.clear();
				recvBufferSizePerProc.clear();
				delete [] req; req = NULL;
				delete [] stats; stats = NULL;
			}

			//FIND FIRST INDEX FOR FIRST MORTONS IN EACH PROCESS
			for (int iproc=0; iproc<nproc; iproc++){
				vector<Class_Octant<2> >::iterator oend = octree.octants.end();
				vector<Class_Octant<2> >::iterator obegin = octree.octants.begin();
				vector<Class_Octant<2> >::iterator it = obegin;
				int nmortons = FirstMortonReceived[iproc].size();
				FirstIndexperproc[iproc].resize(nmortons);
				for (int idx=0; idx<nmortons; idx++){
					FirstIndexperproc[iproc][idx] = octree.getNumOctants()-1;
					uint32_t idx1 = 0;
					mortonfirstdesc = FirstMortonReceived[iproc][idx];
					morton1 = it->computeMorton();
					while(morton1 <= mortonfirstdesc && it != oend){
						idx1++;
						FirstIndexperproc[iproc][idx] = idx1;
						it++;
						if (it != oend)
							morton1 = it->computeMorton();
					}
					if(idx1 > 0){
						idx1--;
						it--;
						morton1 = ptree.getOctant(idx1)->computeMorton();
					}
				}
			}

			//FIND SECOND INDEX FOR SECOND MORTONS IN EACH PROCESS
			for (int iproc=0; iproc<nproc; iproc++){
				vector<Class_Octant<2> >::iterator oend = octree.octants.end();
				vector<Class_Octant<2> >::iterator obegin = octree.octants.begin();
				vector<Class_Octant<2> >::iterator it = obegin;
				int nmortons = SecondMortonReceived[iproc].size();
				SecondIndexperproc[iproc].resize(nmortons);
				for (int idx=0; idx<nmortons; idx++){
					SecondIndexperproc[iproc][idx] = octree.getNumOctants()-1;
					uint32_t idx2 = 0;
					mortonlastdesc = SecondMortonReceived[iproc][idx];
					morton2 = it->computeMorton();
					while(morton2 <= mortonlastdesc && it != oend){
						SecondIndexperproc[iproc][idx] = idx2;
						idx2++;
						it++;
						if (it != oend)
							morton2 = it->computeMorton();
					}
					if(idx2 > 0){
						idx2--;
						it--;
						morton2 = ptree.getOctant(idx2)->computeMorton();
					}
				}
			}

			for(int iproc=0; iproc<nproc; iproc++){
				FirstIndexperproc[iproc].push_back(-1);
				SecondIndexperproc[iproc].push_back(-1);
			}


			{
				//COMM BACK FIRST INDEX PER PROC
				map<int,Class_Comm_Buffer> sendBuffers;
				map<int,vector<uint32_t> >::iterator bitend = FirstIndexperproc.end();
				for(map<int,vector<uint32_t> >::iterator bit = FirstIndexperproc.begin(); bit != bitend; ++bit){
					int buffSize = bit->second.size() * (int)ceil((double)(sizeof(uint32_t)) / (double)(CHAR_BIT/8));
					int key = bit->first;
					vector<uint32_t> & value = bit->second;
					sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					int nofIndices = value.size();
					for(int i = 0; i < nofIndices; ++i){
						//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
						uint32_t Index = value[i];
						error_flag = MPI_Pack(&Index,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
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
				//at the same time every process compute the size in bytes
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

				//UNPACK BUFFERS AND BUILD CONTAINER OF RECEIVED MORTON
				//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked.
				//every Index is built and put in the mapper
				uint32_t Indexcounter = 0;
				uint32_t Index = 0;
				map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
				for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
					int pos = 0;
					int nofIndexPerProc = int(rrit->second.commBufferSize / (uint32_t) (sizeof(uint32_t)));
					for(int i = 0; i < nofIndexPerProc-1; ++i){
						error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&Index,1,MPI_UINT32_T,comm);
						mapper[FirstLocalIndex[rrit->first][i]].first.first = Index;
						++Indexcounter;
					}

				}
				recvBuffers.clear();
				sendBuffers.clear();
				recvBufferSizePerProc.clear();
				delete [] req; req = NULL;
				delete [] stats; stats = NULL;
			}

			{
				//COMM BACK SECOND INDEX PER PROC
				map<int,Class_Comm_Buffer> sendBuffers;
				map<int,vector<uint32_t> >::iterator bitend = SecondIndexperproc.end();
				uint32_t pbordersOversize = 0;
				for(map<int,vector<uint32_t> >::iterator bit = SecondIndexperproc.begin(); bit != bitend; ++bit){
					pbordersOversize += bit->second.size();
					int buffSize = bit->second.size() * (int)ceil((double)(sizeof(uint32_t)) / (double)(CHAR_BIT/8));
					int key = bit->first;
					vector<uint32_t> & value = bit->second;
					sendBuffers[key] = Class_Comm_Buffer(buffSize,'a',comm);
					int pos = 0;
					int nofIndices = value.size();
					for(int i = 0; i < nofIndices; ++i){
						//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
						uint64_t Index = value[i];
						error_flag = MPI_Pack(&Index,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,comm);
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
				//at the same time every process compute the size in bytes
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

				//UNPACK BUFFERS AND BUILD CONTAINER OF RECEIVED MORTON
				//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked.
				//every Index is built and put in the mapper
				uint32_t Indexcounter = 0;
				uint32_t Index = 0;
				map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
				for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
					int pos = 0;
					int nofIndexPerProc = int(rrit->second.commBufferSize / (uint32_t) (sizeof(uint32_t)));
					for(int i = 0; i < nofIndexPerProc-1; ++i){
						error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&Index,1,MPI_UINT32_T,comm);
						mapper[SecondLocalIndex[rrit->first][i]].first.second = Index;
						++Indexcounter;
					}
				}
				recvBuffers.clear();
				sendBuffers.clear();
				recvBufferSizePerProc.clear();
				delete [] req; req = NULL;
				delete [] stats; stats = NULL;
			}
		}
#endif /* NOMPI */
		//TODO PARALLEL VERSION - (BUGS!!!)
		return mapper;
	}

	// =============================================================================== //

	/** Write the logical octree mesh in .vtu format in a user-defined file.
	 * If the connectivity is not stored, the method temporary computes it.
	 * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
	 * \param[in] filename Seriously?....
	 */
	void writeLogical(string filename) {

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
			for(int j = 0; j < 3; ++j)
				out << std::setprecision(6) << octree.nodes[i][j] << " ";
			if((i+1)%4==0 && i!=nofNodes-1)
				out << endl << "          ";
		}
		for(int i = 0; i < nofGhostNodes; i++)
		{
			for(int j = 0; j < 3; ++j)
				out << std::setprecision(6) << octree.ghostsnodes[i][j] << " ";
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
			for(int j = 0; j < global2D.nnodes; j++)
			{
				int jj;
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
				}
				out << octree.connectivity[i][jj] << " ";
			}
			if((i+1)%3==0 && i!=nofOctants-1)
				out << endl << "          ";
		}
		for(int i = 0; i < nofGhosts; i++)
		{
			for(int j = 0; j < global2D.nnodes; j++)
			{
				int jj;
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
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
			out << (i+1)*global2D.nnodes << " ";
			if((i+1)%12==0 && i!=nofAll-1)
				out << endl << "          ";
		}
		out << endl << "        </DataArray>" << endl
				<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
				<< "          ";
		for(int i = 0; i < nofAll; i++)
		{
			int type;
			type = 9;
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

		if (clear){
			octree.clearConnectivity();
		}


	}

	// ----------------------------------------------------------------------------------- //

	/** Write the physical octree mesh in .vtu format in a user-defined file.
	 * If the connectivity is not stored, the method temporary computes it.
	 * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
	 * \param[in] filename Seriously?....
	 */
	void write(string filename) {

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
			for(int j = 0; j < global2D.nnodes; j++)
			{
				int jj;
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
				}
				out << octree.connectivity[i][jj] << " ";
			}
			if((i+1)%3==0 && i!=nofOctants-1)
				out << endl << "          ";
		}
		for(int i = 0; i < nofGhosts; i++)
		{
			for(int j = 0; j < global2D.nnodes; j++)
			{
				int jj;
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
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
			out << (i+1)*global2D.nnodes << " ";
			if((i+1)%12==0 && i!=nofAll-1)
				out << endl << "          ";
		}
		out << endl << "        </DataArray>" << endl
				<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
				<< "          ";
		for(int i = 0; i < nofAll; i++)
		{
			int type;
			type = 9;
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

	// =============================================================================== //

	/** Write the physical octree mesh in .vtu format with data for test in a user-defined file.
	 * If the connectivity is not stored, the method temporary computes it.
	 * The method doesn't write the ghosts on file.
	 * \param[in] filename Seriously?....
	 */
	void writeTest(string filename, vector<double> data) {

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
			for(int j = 0; j < global2D.nnodes; j++)
			{
				int jj;
				if (j<2){
					jj = j;
				}
				else if(j==2){
					jj = 3;
				}
				else if(j==3){
					jj = 2;
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
			out << (i+1)*global2D.nnodes << " ";
			if((i+1)%12==0 && i!=nofAll-1)
				out << endl << "          ";
		}
		out << endl << "        </DataArray>" << endl
				<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
				<< "          ";
		for(int i = 0; i < nofAll; i++)
		{
			int type;
			type = 9;
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

};



