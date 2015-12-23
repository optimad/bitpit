#ifndef CLASSPARATREE_HPP_
#define CLASSPARATREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#if NOMPI==0
#include <mpi.h>
#endif
#include "classGlobal.hpp"
#include "classOctant.hpp"
#include "classLocalTree.hpp"
#include "Class_Comm_Buffer.hpp"
#include "classMap.hpp"
#include "Class_Data_Comm_Interface.hpp"
#include "Class_Data_LB_Interface.hpp"
#include "Class_Log.hpp"
#include <cstdint>
#include <iterator>
#include <set>
#include <algorithm>
#include <string>
#include <functional>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <bitset>
#include <array>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			17/dec/2015
 *	\authors		Marco Cisternino
 *	\authors		Edoardo Lombardi
 *
 *	\brief Para Tree is the user interface class
 *
 *	The user should (read can...) work only
 *	with this Class and its methods.
 *	The sizes are intended in physical domain. The transformation from the logical
 *	domain to the physical domain is defined by classMap trans.
 *
 *	The partition of the octree is performed by following the Z-curve defined by the Morton
 *	index of the octants. By default it is a balanced partition over the number of octants for each
 *	process.
 *
 *	Class classParaTree has a dimensional parameter int dim and it accepts only two values: dim=2 (Class_Para_Tree<2>)and dim=3 (Class_Para_Tree<3>), obviously for 2D and 3D respectively.
 */
class classParaTree{

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	//undistributed members
	uint64_t* 			partition_first_desc; 			/**<Global array containing position of the first possible octant in each processor*/
	uint64_t*			partition_last_desc; 			/**<Global array containing position of the last possible octant in each processor*/
	uint64_t* 			partition_range_globalidx;	 	/**<Global array containing global index of the last existing octant in each processor*/
	uint64_t 			global_num_octants;   			/**<Global number of octants in the parallel octree*/
	map<int,u32vector> 	bordersPerProc;					/**<Local indices of border octants per process*/
	int 				nproc;							/**<Number of processes of the job*/
	uint8_t 			max_depth;						/**<Global max existing level in the parallel octree*/
	classGlobal			global;							/**<Global variables*/

	//distributed members
	int 				rank;							/**<Local rank of process*/
	classLocalTree 		octree;							/**<Local tree in each processor*/

	//distributed adpapting memebrs
	u32vector 			mapidx;							/**<Local mapper for adapting. Mapper from new octants to old octants.
														mapidx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
														if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
														if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 */

	//auxiliary members
	int 				error_flag;						/**<MPI error flag*/
	bool 				serial;							/**<True if the octree is the same on each processor, False if the octree is distributed*/

	//map members
	classMap 			trans;							/**<Transformation map from logical to physical domain*/
	uint8_t				dim;							/**<Space dimension of the octree object (2D/3D).*/

	//info member
	uint64_t			status;							/**<Label of actual status of octree (incremental after an adpat
														with at least one modifyed element).*/

	//log member
	Class_Log 			log;							/**<Log object*/

#if NOMPI==0
	MPI_Comm 			comm;							/**<MPI communicator*/
#endif

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:

#if NOMPI==0
	classParaTree(uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);// : log(logfile,comm_),comm(comm_);
	classParaTree(double X, double Y, double Z, double L, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);//:dim(2),trans(X,Y,Z,L),log(logfile,comm_),comm(comm_);
	classParaTree(double X, double Y, double Z, double L, u32vector2D & XYZ, u8vector & levels, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);//:trans(X,Y,Z,L),log(logfile,comm_),comm(comm_);
#else
	classParaTree(uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");// : log(logfile);
	classParaTree(double X, double Y, double Z, double L, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");//:dim(2),trans(X,Y,Z,L),log(logfile);
	classParaTree(double X, double Y, double Z, double L, u32vector2D & XYZ, u8vector & levels, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");//:trans(X,Y,Z,L),log(logfile);
#endif

	~classParaTree();

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //

	int 		getRank();
	int 		getMaxLevel();
	uint32_t 	getMaxLength();
	uint8_t 	getNnodes();
	uint8_t 	getNfaces();
	uint8_t 	getNedges();
	uint8_t 	getNchildren();
	uint8_t 	getNnodesperface();
	void 		getNormals(int8_t normals[6][3]);
	void 		getOppface(uint8_t oppface[4]);
	void 		getFacenode(uint8_t facenode[6][3]);
	void 		getNodeface(uint8_t nodeface[8][3]);
	void 		getEdgeface(uint8_t edgeface[12][2]);
	void 		getNodecoeffs(int8_t nodecoeffs[8][3]);
	void 		getEdgecoeffs(int8_t edgecoeffs[12][3]);
	void 		setMaxLevel(int8_t maxlevel);

	// =================================================================================== //
	// INDEX BASED METHODS
	// =================================================================================== //

	double 		getX(uint32_t idx);
	double 		getY(uint32_t idx);
	double 		getZ(uint32_t idx);
	double 		getSize(uint32_t idx);
	double 		getArea(uint32_t idx);
	double 		getVolume(uint32_t idx);
	void 		getCenter(uint32_t idx, dvector& center);
	dvector 	getCenter(uint32_t idx);
	dvector 	getFaceCenter(uint32_t idx, uint8_t iface);
	void 		getFaceCenter(uint32_t idx, uint8_t iface, dvector& center);
	dvector 	getNode(uint32_t idx, uint8_t inode);
	void 		getNode(uint32_t idx, uint8_t inode, dvector& node);
	void 		getNodes(uint32_t idx, dvector2D & nodes);
	dvector2D 	getNodes(uint32_t idx);
	void 		getNormal(uint32_t idx, uint8_t & iface, dvector & normal);
	dvector 	getNormal(uint32_t idx, uint8_t & iface);
	int8_t 		getMarker(uint32_t idx);
	uint8_t 	getLevel(uint32_t idx);
	bool 		getBalance(uint32_t idx);
#if NOMPI==0
	bool 		getIsGhost(uint32_t idx);
#endif
	bool 		getIsNewR(uint32_t idx);
	bool 		getIsNewC(uint32_t idx);
	uint64_t 	getGlobalIdx(uint32_t idx);
	uint64_t 	getGhostGlobalIdx(uint32_t idx);
	void 		setMarker(uint32_t idx, int8_t marker);
	void 		setBalance(uint32_t idx, bool balance);

	// =================================================================================== //
	// POINTER BASED METHODS
	// =================================================================================== //

	double 		getX(classOctant* oct);
	double 		getY(classOctant* oct);
	double 		getZ(classOctant* oct);
	double 		getSize(classOctant* oct);
	double 		getArea(classOctant* oct);
	double 		getVolume(classOctant* oct);
	void 		getCenter(classOctant* oct, dvector& center);
	dvector 	getCenter(classOctant* oct);
	dvector 	getFaceCenter(classOctant* oct, uint8_t iface);
	void 		getFaceCenter(classOctant* oct, uint8_t iface, dvector& center);
	dvector 	getNode(classOctant* oct, uint8_t inode);
	void 		getNode(classOctant* oct, uint8_t inode, dvector& node);
	void 		getNodes(classOctant* oct, dvector2D & nodes);
	dvector2D 	getNodes(classOctant* oct);
	void 		getNormal(classOctant* oct, uint8_t & iface, dvector & normal);
	dvector 	getNormal(classOctant* oct, uint8_t & iface);
	int8_t 		getMarker(classOctant* oct);
	uint8_t 	getLevel(classOctant* oct);
	bool 		getBalance(classOctant* oct);
	bool 		getIsNewR(classOctant* oct);
	bool 		getIsNewC(classOctant* oct);
	void 		setMarker(classOctant* oct, int8_t marker);
	void 		setBalance(classOctant* oct, bool balance);

	// =================================================================================== //
	// LOCAL TREE GET/SET METHODS
	// =================================================================================== //

	uint64_t 	getStatus();
	uint32_t 	getNumOctants() const;
	uint32_t 	getNumGhosts() const;
	uint32_t 	getNumNodes() const;
	uint8_t 	getLocalMaxDepth() const;
	uint8_t 	getBalanceCodimension() const;
	void 		getBoundingBox(dvector & P0, dvector & P1);
	void 		getBoundingBox(darray3 & P0, darray3 & P1);
	void 		setBalanceCodimension(uint8_t b21codim);
	const classOctant & getFirstDesc() const;
	const classOctant & getLastDesc() const;
	uint64_t 	getLastDescMorton(uint32_t idx);

	// =================================================================================== //
	// INTERSECTION GET/SET METHODS
	// =================================================================================== //

	uint32_t 	getNumIntersections();
	classIntersection* getIntersection(uint32_t idx);
	uint8_t 	getLevel(classIntersection* inter);
	bool 		getFiner(classIntersection* inter);
	bool 		getBound(classIntersection* inter);
	bool 		getIsGhost(classIntersection* inter);
	bool 		getPbound(classIntersection* inter);
	uint8_t 	getFace(classIntersection* inter);
	u32vector 	getOwners(classIntersection* inter);
	uint32_t 	getIn(classIntersection* inter);
	uint32_t 	getOut(classIntersection* inter);
	double 		getSize(classIntersection* inter);
	double 		getArea(classIntersection* inter);
	dvector 	getCenter(classIntersection* inter);
	dvector2D 	getNodes(classIntersection* inter);
	dvector 	getNormal(classIntersection* inter);

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

	classOctant* getOctant(uint32_t idx);
	classOctant* getGhostOctant(uint32_t idx);
	uint64_t 	getGlobalIdx(classOctant* oct);
	uint32_t 	getIdx(classOctant* oct);
	uint32_t 	getIdx(classOctant oct);
#if NOMPI==0
	bool 		getIsGhost(classOctant* oct);
	bool 		getIsGhost(classOctant oct);
#endif

	// =================================================================================== //
	// PRIVATE GET/SET METHODS
	// =================================================================================== //
private:
	void 		setFirstDesc();
	void 		setLastDesc();

	// =================================================================================== //
	// OTHER METHODS												    			   //
	// =================================================================================== //

	// =================================================================================== //
	// OTHER OCTANT BASED METHODS												    			   //
	// =================================================================================== //
public:
	void 		findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, vector<bool> & isghost);
	void 		findNeighbours(classOctant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, vector<bool> & isghost);
	void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours);
	classOctant* getPointOwner(dvector & point);
	uint32_t 	getPointOwnerIdx(dvector & point);
	void 		getMapping(uint32_t & idx, u32vector & mapper, vector<bool> & isghost);

	// =================================================================================== //
	// OTHER PARATREE BASED METHODS												    			   //
	// =================================================================================== //

	int 		findOwner(const uint64_t & morton);
	bool 		adapt(bool mapper_flag = false);
	bool 		adaptGlobalRefine(bool mapper_flag = false);
	bool 		adaptGlobalCoarse(bool mapper_flag = false);
	void 		computeConnectivity();
	void 		clearConnectivity();
	void 		updateConnectivity();
	const u32vector2D & getConnectivity();
	const u32vector & getConnectivity(uint32_t idx);
	const u32vector & getConnectivity(classOctant* oct);
	const u32arr3vector & getNodes();
	const u32array3 & getNodeLogicalCoordinates(uint32_t inode);
	dvector 	getNodeCoordinates(uint32_t inode);
	void 		computeGhostsConnectivity();
	void 		clearGhostsConnectivity();
	void 		updateGhostsConnectivity();
	const u32vector2D & getGhostConnectivity();
	const u32vector & getGhostConnectivity(uint32_t idx);
	const u32vector & getGhostConnectivity(classOctant* oct);
	const u32arr3vector & getGhostNodes();
	const u32array3 & getGhostNodeLogicalCoordinates(uint32_t inode);
	dvector 	getGhostNodeCoordinates(uint32_t inode);
	//TODO MapPablos
#if NOMPI==0
	void 		loadBalance();
	void 		loadBalance(uint8_t & level);
#endif

	// =================================================================================== //
	// OTHER INTERSECTION BASED METHODS												    			   //
	// =================================================================================== //

	void 		computeIntersections();

	// =================================================================================== //
	// OTHER PRIVATE METHODS												    			   //
	// =================================================================================== //
private:
	classOctant& extractOctant(uint32_t idx);
	bool 		private_adapt();
	bool 		private_adapt_mapidx(bool mapflag);
	void 		updateAdapt();
#if NOMPI==0
	void 		computePartition(uint32_t* partition);
	void 		computePartition(uint32_t* partition, dvector* weight);
	void 		computePartition(uint32_t* partition, uint8_t & level_);
	void 		updateLoadBalance();
	void 		setPboundGhosts();
	void 		commMarker();
#endif
	void 		updateAfterCoarse();
	void 		updateAfterCoarse(u32vector & mapidx);
	void 		balance21(bool const first);

	// =================================================================================== //
	// TESTING OUTPUT METHODS												    			   //
	// =================================================================================== //
public:
	void 		write(string filename);
	void 		writeTest(string filename, vector<double> data);

	// =================================================================================== //
	// TEMPLATE METHODS												    			       //
	// =================================================================================== //

#if NOMPI==0

	/** Communicate data provided by the user between the processes.
	 */
	template<class Impl>
	void
	communicate(Class_Data_Comm_Interface<Impl> & userData){
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

	/** Distribute Load-Balancing the octants of the whole tree and data provided by the user
	 * over the processes of the job following the Morton order.
	 * Until loadBalance is not called for the first time the mesh is serial.
	 */
	template<class Impl>
	void
	loadBalance(Class_Data_LB_Interface<Impl> & userData, dvector* weight = NULL){
		//Write info on log
		log.writeLog("---------------------------------------------");
		log.writeLog(" LOAD BALANCE ");

		uint32_t* partition = new uint32_t [nproc];
		if (weight == NULL)
			computePartition(partition);
		else
			computePartition(partition, weight);

		weight = NULL;

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

			userData.assign(stride,partition[rank]);

			//Update and build ghosts here
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
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
						}
						else{
							for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);

							}
							userData.gather(sendBuffers[p],i);
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
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
						}
						else{
							for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
							//pack octants from lh - partition[p] to lh
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
						}
						lh -= nofElementsFromSuccessiveToPrevious;
						globalLastHead -= nofElementsFromSuccessiveToPrevious;
						headSize = lh + 1;
						if(headSize == 0)
							break;
					}
				}

			}
//			cout << "first" << endl;
			uint32_t nofElementsFromPreviousToSuccessive = 0;
			contatore = 0;
			//build send buffers from Tail
			if(tailSize != 0){
				for(int p = firstSuccessor; p < nproc; ++p){
					if(tailSize < partition[p]){
						nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
						if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
							nofElementsFromPreviousToSuccessive = tailSize;

						uint32_t octantsSize = (uint32_t)octree.octants.size();
						int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
						}
						else{
							for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
							//PACK octants from ft to octantsSize-1
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
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
						uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
						int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
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
						MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
						}
						ft += nofElementsFromPreviousToSuccessive;
						globalFirstTail += nofElementsFromPreviousToSuccessive;
						tailSize -= nofElementsFromPreviousToSuccessive;
						if(tailSize == 0)
							break;
					}
				}
			}
//			cout << "second" << endl;

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
				userData.move(nofResidents - k - 1,resCounter - k);
			}

			//UNPACK BUFFERS AND BUILD NEW OCTANTS
			newCounter = 0;
			bool jumpResident = false;

			for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
				uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&x,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&y,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&z,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&l,1,MPI_UINT8_T,comm);
					//octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					octree.octants[newCounter] = classOctant(dim,l,x,y,z);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 17; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					userData.scatter(rbit->second,newCounter);
					++newCounter;
				}
			}
			octvector(octree.octants).swap(octree.octants);
//			cout << "third" << endl;

			userData.shrink();

			delete [] newPartitionRangeGlobalidx; newPartitionRangeGlobalidx = NULL;
			delete [] nofRecvsPerProc; nofRecvsPerProc = NULL;
			delete [] displays; displays = NULL;
			delete [] req; req = NULL;
			delete [] stats; stats = NULL;
			delete [] globalRecvsBuff; globalRecvsBuff = NULL;

			//Update and ghosts here
	//		cout << "in update" << endl;
			updateLoadBalance();
	//		cout << "in setpbound" << endl;
			setPboundGhosts();
			uint32_t nofGhosts = getNumGhosts();
			userData.resizeGhost(nofGhosts);
	//		cout << "fourth" << endl;

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

	/** Distribute Load-Balanced the octants of the whole tree and data provided by the user
	 * over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
	 * The families of octants of a desired level are retained compact on the same process.
	 * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
	 */
	template<class Impl>
	void
	loadBalance(Class_Data_LB_Interface<Impl> & userData, uint8_t & level){

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

			userData.assign(stride,partition[rank]);

			//Update and build ghosts here
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
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
						}
						else{
							for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);

							}
							userData.gather(sendBuffers[p],i);
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
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
						}
						else{
							for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
							//pack octants from lh - partition[p] to lh
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
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

						uint32_t octantsSize = (uint32_t)octree.octants.size();
						int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
						}
						else{
							for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
								buffSize += userData.size(i);
							}
						}
						//add room for int, number of octants in this buffer
						buffSize += sizeof(int);
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a',comm);
						//store the number of octants at the beginning of the buffer
						MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						//USE BUFFER POS
						for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
							//PACK octants from ft to octantsSize-1
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
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
						uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
						int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)global.octantBytes / (double)(CHAR_BIT/8));
						//compute size of data in buffers
						if(userData.fixedSize()){
							buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
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
						MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,comm);
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							//const Class_Octant<2> & octant = octree.octants[i];
							const classOctant & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							for(int j = 0; j < 17; ++j)
								info[j] = octant.info[j];
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							for(int j = 0; j < 17; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,comm);
							}
							userData.gather(sendBuffers[p],i);
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
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&z,1,MPI_UINT32_T,comm);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&l,1,MPI_UINT8_T,comm);
					//octree.octants[newCounter] = Class_Octant<2>(l,x,y);
					octree.octants[newCounter] = classOctant(dim,l,x,y,z);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&m,1,MPI_INT8_T,comm);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 17; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&info[j],1,MPI::BOOL,comm);
						octree.octants[newCounter].info[j] = info[j];
					}
					//TODO Unpack data
					userData.scatter(rbit->second,newCounter);
					++newCounter;
				}
			}
			octvector(octree.octants).swap(octree.octants);

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
			uint32_t nofGhosts = getNumGhosts();
			userData.resizeGhost(nofGhosts);

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



	// =============================================================================== //


};


/*  @}  */

#endif /* CLASSPARATREE_HPP_ */
