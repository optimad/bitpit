/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#ifndef __BITPIT_PARA_TREE_HPP__
#define __BITPIT_PARA_TREE_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#include "CommBuffer.hpp"
#include "DataLBInterface.hpp"
#include "DataCommInterface.hpp"
#endif
#include "Global.hpp"
#include "Array.hpp"
#include "Octant.hpp"
#include "LocalTree.hpp"
#include "Map.hpp"
#include "bitpit_IO.hpp"
#include <map>
#include <unordered_map>
#include <set>
#include <bitset>
#include <algorithm>

namespace bitpit {

// =================================================================================== //
// TYPEDEFS																			   //
// =================================================================================== //
typedef std::vector<bool>				bvector;
typedef std::vector<int>				ivector;
typedef std::bitset<72>					octantID;
typedef std::vector<Octant*>			ptroctvector;
typedef ptroctvector::iterator			octantIterator;

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
 *	The user should (read can...) work only with this class and its methods.
 *	The sizes are intended in reference physical domain with limits [0,1]. The transformation from the logical
 *	domain to the physical domain is defined by an internal mapping.
 *
 *	The partition of the octree is performed by following the Z-curve defined by the Morton
 *	index of the octants. By default it is a balanced partition over the number of octants for each
 *	process.
 *
 *	Class ParaTree has a dimensional parameter int dim and it accepts only two
 *	 values: dim=2 and dim=3, for 2D and 3D respectively.
 */
class ParaTree{

	// =================================================================================== //
	// MEMBERS																			   //
	// =================================================================================== //
private:
	//undistributed members
	uint64_t* 				m_partitionFirstDesc; 			/**<Global array containing position of the first possible octant in each processor*/
	uint64_t*				m_partitionLastDesc; 			/**<Global array containing position of the last possible octant in each processor*/
	uint64_t* 				m_partitionRangeGlobalIdx;	 	/**<Global array containing global index of the last existing octant in each processor*/
	uint64_t* 				m_partitionRangeGlobalIdx0;	 	/**<Global array containing global index of the last existing octant in each processor before the last loadBalance (after an adapt is set equal to the actual.)*/
	uint64_t 				m_globalNumOctants;   			/**<Global number of octants in the parallel octree*/
	int 					m_nproc;						/**<Number of processes of the job*/
	uint8_t 				m_maxDepth;						/**<Global max existing level in the parallel octree*/
	Global					m_global;						/**<Global variables*/

	//distributed members
	int 					m_rank;							/**<Local m_rank of process*/
	LocalTree 				m_octree;						/**<Local tree in each processor*/
	std::map<int,u32vector> m_bordersPerProc;				/**<Local indices of border octants per process*/
	ptroctvector 			m_internals;					/**<Local pointers to internal octants*/
	ptroctvector 			m_pborders;						/**<Local pointers to border of process octants*/

	//distributed adpapting memebrs
	u32vector 				m_mapIdx;						/**<Local mapper for adapting. Mapper from new octants to old octants.
															m_mapIdx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
															if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
															if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 */
	//elements sent during last loadbalance operation
	std::unordered_map<int,std::array<uint32_t,4> >	   m_sentIdx;											  /**<Local mapper for sent elements. Each element refers to the receiver rank and collect the */

	//auxiliary members
	int 					m_errorFlag;					/**<MPI error flag*/
	bool 					m_serial;						/**<True if the octree is the same on each processor, False if the octree is distributed*/
	double					m_tol;							/**<Tolerance for geometric operations.*/

	//map members
	Map 					m_trans;						/**<Transformation map from m_logical to physical domain*/
	uint8_t					m_dim;							/**<Space dimension of the octree object (2D/3D).*/

	//boundary conditions members
	bvector 				m_periodic;						/**<Boolvector: i-th element is true if the i-th boundary face is a periodic interface.*/

	//info member
	uint64_t				m_status;						/**<Label of actual m_status of octree (incremental after an adpat
															with at least one modifyed element).*/
	std::string				m_lastOp;						/**<Last adapting operation type (adapt or loadbalance).*/

	//log member
	Logger* 				m_log;							/**<Log object pointer*/

	//communicator
#if BITPIT_ENABLE_MPI==1
	//TODO Duplicate communicator
	MPI_Comm 				m_comm;							/**<MPI communicator*/
#endif

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS														   //
	// =================================================================================== //
public:
#if BITPIT_ENABLE_MPI==1
	ParaTree(uint8_t dim = 2, int8_t maxlevel = 20, std::string logfile = "PABLO", MPI_Comm comm = MPI_COMM_WORLD);
	ParaTree(u32vector2D & XYZ, u8vector & levels, uint8_t dim = 2, int8_t maxlevel = 20,  std::string logfile = "PABLO", MPI_Comm comm = MPI_COMM_WORLD);
#else
	ParaTree(uint8_t dim = 2, int8_t maxlevel = 20,  std::string logfile = "PABLO");
	ParaTree(u32vector2D & XYZ, u8vector & levels, uint8_t dim = 2, int8_t maxlevel = 20,  std::string logfile = "PABLO");
#endif
	~ParaTree();

	// =================================================================================== //
	// METHODS																			   //
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS															   //
	// =================================================================================== //
	uint8_t 	getDim();
	uint64_t 	getGlobalNumOctants();
	bool		getSerial();
	bool		getParallel();
	int 		getRank();
	int 		getNproc();
	Logger& 	getLog();
#if BITPIT_ENABLE_MPI==1
	void		setComm(MPI_Comm communicator);
	MPI_Comm	getComm() const;
	bool		isCommSet() const;
#endif
	uint64_t*	getPartitionRangeGlobalIdx();
	darray3		getOrigin();
	double		getX0();
	double		getY0();
	double		getZ0();
	double		getL();
	int 		getMaxLevel();
	uint32_t 	getMaxLength();
	uint8_t 	getNnodes();
	uint8_t 	getNfaces();
	uint8_t 	getNedges();
	uint8_t 	getNchildren();
	uint8_t 	getNnodesperface();
	void 		getNormals(int8_t normals[6][3]);
	void 		getOppface(uint8_t oppface[6]);
	void 		getFacenode(uint8_t facenode[6][4]);
	void 		getNodeface(uint8_t nodeface[8][3]);
	void 		getEdgeface(uint8_t edgeface[12][2]);
	void 		getNodecoeffs(int8_t nodecoeffs[8][3]);
	void 		getEdgecoeffs(int8_t edgecoeffs[12][3]);
	int8_t		(*getNormals())[3];
	uint8_t		(*getOppface());
	uint8_t 	(*getFacenode())[4];
	uint8_t 	(*getNodeface())[3];
	uint8_t 	(*getEdgeface())[2];
	int8_t 		(*getNodecoeffs())[3];
	int8_t 		(*getEdgecoeffs())[3];
	bvector		getPeriodic();
	double		getTol();
	bool		getPeriodic(uint8_t i);
	void 		setMaxLevel(int8_t maxlevel);
	void		setPeriodic(uint8_t i);
	void		setTol(double tol = 1.0e-14);

	// =================================================================================== //
	// INDEX BASED METHODS																   //
	// =================================================================================== //
	darray3 	getCoordinates(uint32_t idx);
	double 		getX(uint32_t idx);
	double 		getY(uint32_t idx);
	double 		getZ(uint32_t idx);
	double 		getSize(uint32_t idx);
	double 		getArea(uint32_t idx);
	double 		getVolume(uint32_t idx);
	void 		getCenter(uint32_t idx, darray3& center);
	darray3 	getCenter(uint32_t idx);
	darray3 	getFaceCenter(uint32_t idx, uint8_t iface);
	void 		getFaceCenter(uint32_t idx, uint8_t iface, darray3& center);
	darray3 	getNode(uint32_t idx, uint8_t inode);
	void 		getNode(uint32_t idx, uint8_t inode, darray3& node);
	void 		getNodes(uint32_t idx, darr3vector & nodes);
	darr3vector getNodes(uint32_t idx);
	void 		getNormal(uint32_t idx, uint8_t & iface, darray3 & normal);
	darray3 	getNormal(uint32_t idx, uint8_t & iface);
	int8_t 		getMarker(uint32_t idx);
	uint8_t 	getLevel(uint32_t idx);
	uint64_t 	getMorton(uint32_t idx);
	bool 		getBalance(uint32_t idx);
	bool		getBound(uint32_t idx, uint8_t iface);
	bool		getBound(uint32_t idx);
	bool		getPbound(uint32_t idx, uint8_t iface);
	bool		getPbound(uint32_t idx);
	bool 		getIsNewR(uint32_t idx);
	bool 		getIsNewC(uint32_t idx);
	uint64_t 	getGlobalIdx(uint32_t idx);
	uint64_t 	getGhostGlobalIdx(uint32_t idx);
    uint32_t    getLocalIdx(uint64_t gidx);
    uint32_t    getGhostLocalIdx(uint64_t gidx);
	octantID	getPersistentIdx(uint32_t idx);
	void 		setMarker(uint32_t idx, int8_t marker);
	void 		setBalance(uint32_t idx, bool balance);

	// =================================================================================== //
	// POINTER BASED METHODS															   //
	// =================================================================================== //
	darray3 	getCoordinates(Octant* oct);
	double 		getX(Octant* oct);
	double 		getY(Octant* oct);
	double 		getZ(Octant* oct);
	double 		getSize(Octant* oct);
	double 		getArea(Octant* oct);
	double 		getVolume(Octant* oct);
	void 		getCenter(Octant* oct, darray3& center);
	darray3 	getCenter(Octant* oct);
	darray3 	getFaceCenter(Octant* oct, uint8_t iface);
	void 		getFaceCenter(Octant* oct, uint8_t iface, darray3& center);
	darray3 	getNode(Octant* oct, uint8_t inode);
	void 		getNode(Octant* oct, uint8_t inode, darray3& node);
	void 		getNodes(Octant* oct, darr3vector & nodes);
	darr3vector getNodes(Octant* oct);
	void 		getNormal(Octant* oct, uint8_t & iface, darray3 & normal);
	darray3 	getNormal(Octant* oct, uint8_t & iface);
	int8_t 		getMarker(Octant* oct);
	uint8_t 	getLevel(Octant* oct);
	uint64_t 	getMorton(Octant* oct);
	bool 		getBalance(Octant* oct);
	bool		getBound(Octant* oct, uint8_t iface);
	bool		getBound(Octant* oct);
	bool		getPbound(Octant* oct, uint8_t iface);
	bool		getPbound(Octant* oct);
	bool 		getIsNewR(Octant* oct);
	bool 		getIsNewC(Octant* oct);
	uint32_t 	getIdx(Octant* oct);
	uint64_t 	getGlobalIdx(Octant* oct);
	octantID	getPersistentIdx(Octant* oct);
	void 		setMarker(Octant* oct, int8_t marker);
	void 		setBalance(Octant* oct, bool balance);

	// =================================================================================== //
	// LOCAL TREE GET/SET METHODS														   //
	// =================================================================================== //
	uint64_t 	getStatus();
	uint32_t 	getNumOctants() const;
	uint32_t 	getNumGhosts() const;
	uint32_t 	getNumNodes() const;
	uint8_t 	getLocalMaxDepth() const;
	double	 	getLocalMaxSize();
	double	 	getLocalMinSize();
	uint8_t 	getBalanceCodimension() const;
	const Octant & getFirstDesc() const;
	const Octant & getLastDesc() const;
	uint64_t 	getLastDescMorton(uint32_t idx);
	octantIterator	getInternalOctantsBegin();
	octantIterator	getInternalOctantsEnd();
	octantIterator	getPboundOctantsBegin();
	octantIterator	getPboundOctantsEnd();
	void 		setBalanceCodimension(uint8_t b21codim);

	// =================================================================================== //
	// INTERSECTION GET/SET METHODS														   //
	// =================================================================================== //
	uint32_t 	getNumIntersections();
	Intersection* getIntersection(uint32_t idx);
	uint8_t 	getLevel(Intersection* inter);
	bool 		getFiner(Intersection* inter);
	bool 		getBound(Intersection* inter);
	bool 		getIsGhost(Intersection* inter);
	bool 		getPbound(Intersection* inter);
	uint8_t 	getFace(Intersection* inter);
	u32vector 	getOwners(Intersection* inter);
	uint32_t 	getIn(Intersection* inter);
	uint32_t 	getOut(Intersection* inter);
	bool		getOutIsGhost(Intersection* inter);
	double 		getSize(Intersection* inter);
	double 		getArea(Intersection* inter);
	darray3 	getCenter(Intersection* inter);
	darr3vector getNodes(Intersection* inter);
	darray3 	getNormal(Intersection* inter);

	// =================================================================================== //
	// OTHER GET/SET METHODS															   //
	// =================================================================================== //
	Octant*	getOctant(uint32_t idx);
	Octant*	getGhostOctant(uint32_t idx);
	uint32_t 		getIdx(Octant oct);
	bool 			getIsGhost(Octant* oct);
	bool 			getIsGhost(Octant oct);
	const std::unordered_map<int,std::array<uint32_t,4> > & getSentIdx();

	// =================================================================================== //
	// PRIVATE GET/SET METHODS															   //
	// =================================================================================== //
private:
	void 		setFirstDesc();
	void 		setLastDesc();

	// =================================================================================== //
	// OTHER METHODS												    			       //
	// =================================================================================== //

	// =================================================================================== //
	// OTHER OCTANT BASED METHODS												    	   //
	// =================================================================================== //
public:
	void 		findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
	void 		findNeighbours(Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const ;
	void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours) const;
	Octant* 	getPointOwner(dvector point);
	uint32_t 	getPointOwnerIdx(dvector point);
	Octant* 	getPointOwner(darray3 point);
	uint32_t 	getPointOwnerIdx(darray3 point);
	void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost);
	void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost, ivector & rank);

	// =================================================================================== //
	// OTHER PARATREE BASED METHODS												    	   //
	// =================================================================================== //
	uint8_t		getMaxDepth() const;
	int 		findOwner(const uint64_t & morton);
        int             getOwnerRank(const uint64_t & globalIdx);
	bool 		adapt(bool mapper_flag = false);
	bool 		adaptGlobalRefine(bool mapper_flag = false);
	bool 		adaptGlobalCoarse(bool mapper_flag = false);
	void 		computeConnectivity();
	void 		clearConnectivity();
	void 		updateConnectivity();
	const u32vector2D & getConnectivity();
	const u32vector & getConnectivity(uint32_t idx);
	const u32vector & getConnectivity(Octant* oct);
	const u32arr3vector & getNodes();
	const u32array3 & getNodeLogicalCoordinates(uint32_t inode);
	darray3 	getNodeCoordinates(uint32_t inode);
	const u32vector2D & getGhostConnectivity();
	const u32vector & getGhostConnectivity(uint32_t idx);
	const u32vector & getGhostConnectivity(Octant* oct);
#if BITPIT_ENABLE_MPI==1
	void 		loadBalance(dvector* weight = NULL);
	void 		loadBalance(uint8_t & level, dvector* weight = NULL);
private:
	void 		privateLoadBalance(uint32_t* partition);
#endif
public:
	double		levelToSize(uint8_t & level);

	// =================================================================================== //
	// OTHER INTERSECTION BASED METHODS										     		   //
	// =================================================================================== //
	void 		computeIntersections();

	// =================================================================================== //
	// OTHER PRIVATE METHODS												    		   //
	// =================================================================================== //
private:
	Octant& extractOctant(uint32_t idx);
	bool 		private_adapt_mapidx(bool mapflag);
	void 		updateAdapt();
#if BITPIT_ENABLE_MPI==1
	void 		computePartition(uint32_t* partition);
	void 		computePartition(uint32_t* partition, dvector* weight);
	void 		computePartition(uint32_t* partition, uint8_t & level_, dvector* weight);
	void 		updateLoadBalance();
	void 		setPboundGhosts();
	void 		commMarker();
#endif
	void 		updateAfterCoarse();
	void 		updateAfterCoarse(u32vector & mapidx);
	void 		balance21(bool const first);
	void		createPartitionInfo(bool deletePrevious);
	void		deletePartitionInfo();

	// =================================================================================== //
	// TESTING OUTPUT METHODS												    			   //
	// =================================================================================== //
public:
	void 		write(std::string filename);
	void 		writeTest(std::string filename, dvector data);

	// =================================================================================== //
	// TEMPLATE METHODS												    			       //
	// =================================================================================== //
#if BITPIT_ENABLE_MPI==1

	/** Communicate data provided by the user between the processes.
	 */
	template<class Impl>
	void
	communicate(DataCommInterface<Impl> & userData){
		//BUILD SEND BUFFERS
		std::map<int,CommBuffer> sendBuffers;
		size_t fixedDataSize = userData.fixedSize();
		std::map<int,u32vector >::iterator bitend = m_bordersPerProc.end();
		std::map<int,u32vector >::iterator bitbegin = m_bordersPerProc.begin();
		for(std::map<int,u32vector >::iterator bit = bitbegin; bit != bitend; ++bit){
			const int & key = bit->first;
			const u32vector & pborders = bit->second;
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
			sendBuffers[key] = CommBuffer(buffSize,'a',m_comm);
			//store number of pborders from this proc at the begining
			MPI_Pack(&nofPbordersPerProc,1,MPI_INT,sendBuffers[key].m_commBuffer,sendBuffers[key].m_commBufferSize,&sendBuffers[key].m_pos,m_comm);

			//WRITE SEND BUFFERS
			for(size_t j = 0; j < nofPbordersPerProc; ++j){
				userData.gather(sendBuffers[key],pborders[j]);
			}
		}

		//Communicate Buffers Size
		MPI_Request* req = new MPI_Request[sendBuffers.size()*2];
		MPI_Status* stats = new MPI_Status[sendBuffers.size()*2];
		int nReq = 0;
		std::map<int,int> recvBufferSizePerProc;
		std::map<int,CommBuffer>::iterator sitend = sendBuffers.end();
		for(std::map<int,CommBuffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
			recvBufferSizePerProc[sit->first] = 0;
			m_errorFlag = MPI_Irecv(&recvBufferSizePerProc[sit->first],1,MPI_UINT32_T,sit->first,m_rank,m_comm,&req[nReq]);
			++nReq;
		}
		std::map<int,CommBuffer>::reverse_iterator rsitend = sendBuffers.rend();
		for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			m_errorFlag =  MPI_Isend(&rsit->second.m_commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,m_comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//Communicate Buffers
		std::map<int,CommBuffer> recvBuffers;
		std::map<int,int>::iterator ritend = recvBufferSizePerProc.end();
		for(std::map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
			recvBuffers[rit->first] = CommBuffer(rit->second,'a',m_comm);
		}
		nReq = 0;
		for(std::map<int,CommBuffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
			m_errorFlag = MPI_Irecv(recvBuffers[sit->first].m_commBuffer,recvBuffers[sit->first].m_commBufferSize,MPI_PACKED,sit->first,m_rank,m_comm,&req[nReq]);
			++nReq;
		}
		for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
			m_errorFlag =  MPI_Isend(rsit->second.m_commBuffer,rsit->second.m_commBufferSize,MPI_PACKED,rsit->first,rsit->first,m_comm,&req[nReq]);
			++nReq;
		}
		MPI_Waitall(nReq,req,stats);

		//READ RECEIVE BUFFERS
		int ghostOffset = 0;
		std::map<int,CommBuffer>::iterator rbitend = recvBuffers.end();
		std::map<int,CommBuffer>::iterator rbitbegin = recvBuffers.begin();
		for(std::map<int,CommBuffer>::iterator rbit = rbitbegin; rbit != rbitend; ++rbit){
			int nofGhostFromThisProc = 0;
			MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&nofGhostFromThisProc,1,MPI_INT,m_comm);
			for(int k = 0; k < nofGhostFromThisProc; ++k){
				userData.scatter(rbit->second, k+ghostOffset);
			}
			ghostOffset += nofGhostFromThisProc;
		}
		delete [] req; req = NULL;
		delete [] stats; stats = NULL;

	};

	/** Distribute Load-Balancing the octants (with user defined weights) of the whole tree and data provided by the user
	 * over the processes of the job following the Morton order.
	 * Until loadBalance is not called for the first time the mesh is serial.
	 * Even distribute data provided by the user between the processes.
	 * \param[in] userData User interface to distribute the data during loadBalance.
	 * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
	 */
	template<class Impl>
	void
	loadBalance(DataLBInterface<Impl> & userData, dvector* weight = NULL){
		//Write info on m_log
		(*m_log) << "---------------------------------------------" << endl;
		(*m_log) << " LOAD BALANCE " << endl;

		m_sentIdx.clear();
        std::array<uint32_t,4> limits = {{0,0,0,0}};
		if (m_nproc>1){

			uint32_t* partition = new uint32_t [m_nproc];
			if (weight == NULL)
				computePartition(partition);
			else
				computePartition(partition, weight);

			weight = NULL;

			if(m_serial)
			{
				(*m_log) << " " << endl;
				(*m_log) << " Initial Serial distribution : " << endl;
				for(int ii=0; ii<m_nproc; ii++){
					(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]+1)) << endl;
				}

				uint32_t stride = 0;
				for(int i = 0; i < m_rank; ++i)
					stride += partition[i];
				LocalTree::octvector octantsCopy = m_octree.m_octants;
				LocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
				LocalTree::octvector::const_iterator last = first + partition[m_rank];

                limits[1] = stride;
                limits[2] = limits[1] + partition[m_rank];
                limits[3] = m_octree.m_octants.size();
                std::pair<int,std::array<uint32_t,4> > procLimits(m_rank,limits);
                m_sentIdx.insert(procLimits);

				m_octree.m_octants.assign(first, last);
				octvector(m_octree.m_octants).swap(m_octree.m_octants);

				first = octantsCopy.end();
				last = octantsCopy.end();

				userData.assign(stride,partition[m_rank]);

				//Update and build ghosts here
				updateLoadBalance();
				setPboundGhosts();
			}
			else
			{
				(*m_log) << " " << endl;
				(*m_log) << " Initial Parallel partition : " << endl;
				(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
				for(int ii=1; ii<m_nproc; ii++){
					(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
				}

				//empty ghosts
				m_octree.m_ghosts.clear();
				m_octree.m_sizeGhosts = 0;
				//compute new partition range globalidx
				uint64_t* newPartitionRangeGlobalidx = new uint64_t[m_nproc];
				for(int p = 0; p < m_nproc; ++p){
					newPartitionRangeGlobalidx[p] = 0;
					for(int pp = 0; pp <= p; ++pp)
						newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
					--newPartitionRangeGlobalidx[p];
				}

				//find resident octants local offset lastHead(lh) and firstTail(ft)
				int32_t lh,ft;
				if(m_rank == 0)
					lh = -1;
				else{
					lh = (int32_t)(newPartitionRangeGlobalidx[m_rank-1] + 1 - m_partitionRangeGlobalIdx[m_rank-1] - 1 - 1);
				}
				if(lh < 0)
					lh = - 1;
				else if(lh > m_octree.m_octants.size() - 1)
					lh = m_octree.m_octants.size() - 1;

				if(m_rank == m_nproc - 1)
					ft = m_octree.m_octants.size();
				else if(m_rank == 0)
					ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] + 1);
				else{
					ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1]);
				}
				if(ft > (int32_t)(m_octree.m_octants.size() - 1))
					ft = m_octree.m_octants.size();
				else if(ft < 0)
					ft = 0;

				//compute size Head and size Tail
				uint32_t headSize = (uint32_t)(lh + 1);
				uint32_t tailSize = (uint32_t)(m_octree.m_octants.size() - ft);
				uint32_t headOffset = headSize;
				uint32_t tailOffset = tailSize;

				//build send buffers
				std::map<int,CommBuffer> sendBuffers;

				//Compute first predecessor and first successor to send buffers to
				int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
				int64_t globalLastHead = (int64_t) lh;
				int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
				int firstPredecessor = -1;
				int firstSuccessor = m_nproc;
				if(m_rank != 0){
					firstOctantGlobalIdx = (int64_t)(m_partitionRangeGlobalIdx[m_rank-1] + 1);
					globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
					globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
					for(int pre = m_rank - 1; pre >=0; --pre){
						if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
							firstPredecessor = pre;
					}
					for(int post = m_rank + 1; post < m_nproc; ++post){
						if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
							firstSuccessor = post;
					}
				}
				else if(m_rank == 0){
					firstSuccessor = 1;
				}
				MPI_Barrier(m_comm); //da spostare prima della prima comunicazione

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

							int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
							limits[1] = (uint32_t)lh + 1;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
								//PACK octants from 0 to lh in sendBuffer[p]
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);

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
							int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
							limits[1] = (uint32_t)lh + 1;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
								//pack octants from lh - partition[p] to lh
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
					for(int p = firstSuccessor; p < m_nproc; ++p){
						if(tailSize < partition[p]){
							nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
							if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
								nofElementsFromPreviousToSuccessive = tailSize;

							int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = ft;
							limits[1] = ft + nofElementsFromPreviousToSuccessive;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
								//PACK octants from ft to octantsSize-1
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
							int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);

							limits[0] = ft;
							limits[1] = ft + nofElementsFromPreviousToSuccessive;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = ft; i <= endOctants; ++i ){
								//PACK octants from ft to ft + partition[p] -1
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
				std::vector<Array> recvs(m_nproc);
				recvs[m_rank] = Array((uint32_t)sendBuffers.size()+1,-1);
				recvs[m_rank].m_array[0] = m_rank;
				int counter = 1;
				std::map<int,CommBuffer>::iterator sitend = sendBuffers.end();
				for(std::map<int,CommBuffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
					recvs[m_rank].m_array[counter] = sit->first;
					++counter;
				}
				int* nofRecvsPerProc = new int[m_nproc];
				m_errorFlag = MPI_Allgather(&recvs[m_rank].m_arraySize,1,MPI_INT,nofRecvsPerProc,1,MPI_INT,m_comm);
				int globalRecvsBuffSize = 0;
				int* displays = new int[m_nproc];
				for(int pp = 0; pp < m_nproc; ++pp){
					displays[pp] = 0;
					globalRecvsBuffSize += nofRecvsPerProc[pp];
					for(int ppp = 0; ppp < pp; ++ppp){
						displays[pp] += nofRecvsPerProc[ppp];
					}
				}
				//int globalRecvsBuff[globalRecvsBuffSize];
				int* globalRecvsBuff = new int[globalRecvsBuffSize];
				m_errorFlag = MPI_Allgatherv(recvs[m_rank].m_array,recvs[m_rank].m_arraySize,MPI_INT,globalRecvsBuff,nofRecvsPerProc,displays,MPI_INT,m_comm);

				std::vector<std::set<int> > sendersPerProc(m_nproc);
				for(int pin = 0; pin < m_nproc; ++pin){
					for(int k = displays[pin]+1; k < displays[pin] + nofRecvsPerProc[pin]; ++k){
						sendersPerProc[globalRecvsBuff[k]].insert(globalRecvsBuff[displays[pin]]);
					}
				}

				//Communicate Octants (size)
				MPI_Request* req = new MPI_Request[sendBuffers.size()+sendersPerProc[m_rank].size()];
				MPI_Status* stats = new MPI_Status[sendBuffers.size()+sendersPerProc[m_rank].size()];
				int nReq = 0;
				std::map<int,int> recvBufferSizePerProc;
				std::set<int>::iterator senditend = sendersPerProc[m_rank].end();
				for(std::set<int>::iterator sendit = sendersPerProc[m_rank].begin(); sendit != senditend; ++sendit){
					recvBufferSizePerProc[*sendit] = 0;
					m_errorFlag = MPI_Irecv(&recvBufferSizePerProc[*sendit],1,MPI_UINT32_T,*sendit,m_rank,m_comm,&req[nReq]);
					++nReq;
				}
				std::map<int,CommBuffer>::reverse_iterator rsitend = sendBuffers.rend();
				for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
					m_errorFlag =  MPI_Isend(&rsit->second.m_commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,m_comm,&req[nReq]);
					++nReq;
				}
				MPI_Waitall(nReq,req,stats);

				//COMMUNICATE THE BUFFERS TO THE RECEIVERS
				//recvBuffers structure is declared and each buffer is initialized to the right size
				//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
				uint32_t nofNewHead = 0;
				uint32_t nofNewTail = 0;
				std::map<int,CommBuffer> recvBuffers;

				std::map<int,int>::iterator ritend = recvBufferSizePerProc.end();
				for(std::map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
					recvBuffers[rit->first] = CommBuffer(rit->second,'a',m_comm);
				}

				nReq = 0;
				for(std::set<int>::iterator sendit = sendersPerProc[m_rank].begin(); sendit != senditend; ++sendit){
					m_errorFlag = MPI_Irecv(recvBuffers[*sendit].m_commBuffer,recvBuffers[*sendit].m_commBufferSize,MPI_PACKED,*sendit,m_rank,m_comm,&req[nReq]);
					++nReq;
				}
				for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
					m_errorFlag =  MPI_Isend(rsit->second.m_commBuffer,rsit->second.m_commBufferSize,MPI_PACKED,rsit->first,rsit->first,m_comm,&req[nReq]);
					++nReq;
				}
				MPI_Waitall(nReq,req,stats);

				//Unpack number of octants per sender
				std::map<int,uint32_t> nofNewOverProcs;
				std::map<int,CommBuffer>::iterator rbitend = recvBuffers.end();
				for(std::map<int,CommBuffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
					uint32_t nofNewPerProc;
					MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&nofNewPerProc,1,MPI_UINT32_T,m_comm);
					nofNewOverProcs[rbit->first] = nofNewPerProc;
					if(rbit->first < m_rank)
						nofNewHead += nofNewPerProc;
					else if(rbit->first > m_rank)
						nofNewTail += nofNewPerProc;
				}

				//MOVE RESIDENT TO BEGIN IN OCTANTS
				uint32_t resEnd = m_octree.getNumOctants() - tailOffset;
				uint32_t nofResidents = resEnd - headOffset;
				uint32_t octCounter = 0;
				for(uint32_t i = headOffset; i < resEnd; ++i){
					m_octree.m_octants[octCounter] = m_octree.m_octants[i];
					userData.move(i,octCounter);
					++octCounter;
				}
				uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
				m_octree.m_octants.resize(newCounter, Octant(m_dim, m_global.m_maxLevel));
				userData.resize(newCounter);
				//MOVE RESIDENTS IN RIGHT POSITION
				uint32_t resCounter = nofNewHead + nofResidents - 1;
				for(uint32_t k = 0; k < nofResidents ; ++k){
					m_octree.m_octants[resCounter - k] = m_octree.m_octants[nofResidents - k - 1];
					userData.move(nofResidents - k - 1,resCounter - k);
				}

				//UNPACK BUFFERS AND BUILD NEW OCTANTS
				newCounter = 0;
				bool jumpResident = false;

				for(std::map<int,CommBuffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
					uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];
					if(rbit->first > m_rank && !jumpResident){
						newCounter += nofResidents ;
						jumpResident = true;
					}
					for(int i = nofNewPerProc - 1; i >= 0; --i){
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&x,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&y,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&z,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&l,1,MPI_UINT8_T,m_comm);
						//m_octree.m_octants[newCounter] = Class_Octant<2>(l,x,y);
						m_octree.m_octants[newCounter] = Octant(m_dim,l,x,y,z,m_global.m_maxLevel);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&m,1,MPI_INT8_T,m_comm);
						m_octree.m_octants[newCounter].setMarker(m);
						for(int j = 0; j < 17; ++j){
							m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&info[j],1,MPI_C_BOOL,m_comm);
							m_octree.m_octants[newCounter].m_info[j] = info[j];
						}
						userData.scatter(rbit->second,newCounter);
						++newCounter;
					}
				}
				octvector(m_octree.m_octants).swap(m_octree.m_octants);

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

			//Write info of final partition on m_log
			(*m_log) << " " << endl;
			(*m_log) << " Final Parallel partition : " << endl;
			(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
			for(int ii=1; ii<m_nproc; ii++){
				(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
			}
			(*m_log) << " " << endl;
			(*m_log) << "---------------------------------------------" << endl;

		}
		else{
			(*m_log) << " " << endl;
			(*m_log) << " Serial partition : " << endl;
			(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
			(*m_log) << " " << endl;
			(*m_log) << "---------------------------------------------" << endl;
		}
	}

	/** Distribute Load-Balanced the octants (with user defined weights) of the whole tree and data provided by the user
	 * over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
	 * The families of octants of a desired level are retained compact on the same process.
	 * Even distribute data provided by the user between the processes.
	 * \param[in] userData User interface to distribute the data during loadBalance.
	 * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
	 * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
	 */
	template<class Impl>
	void
	loadBalance(DataLBInterface<Impl> & userData, uint8_t & level, dvector* weight = NULL){

		//Write info on m_log
		(*m_log) << "---------------------------------------------" << endl;
		(*m_log) << " LOAD BALANCE " << endl;

		m_sentIdx.clear();
        std::array<uint32_t,4> limits = {{0,0,0,0}};
		if (m_nproc>1){

			uint32_t* partition = new uint32_t [m_nproc];
			computePartition(partition, level, weight);

			if(m_serial)
			{
				(*m_log) << " " << endl;
				(*m_log) << " Initial Serial distribution : " << endl;
				for(int ii=0; ii<m_nproc; ii++){
					(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]+1)) << endl;
				}

				uint32_t stride = 0;
				for(int i = 0; i < m_rank; ++i)
					stride += partition[i];
				LocalTree::octvector octantsCopy = m_octree.m_octants;
				LocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
				LocalTree::octvector::const_iterator last = first + partition[m_rank];

                limits[1] = stride;
                limits[2] = limits[1] + partition[m_rank];
                limits[3] = m_octree.m_octants.size();
                std::pair<int,std::array<uint32_t,4> > procLimits(m_rank,limits);
                m_sentIdx.insert(procLimits);

				m_octree.m_octants.assign(first, last);
				octvector(m_octree.m_octants).swap(m_octree.m_octants);


				first = octantsCopy.end();
				last = octantsCopy.end();

				userData.assign(stride,partition[m_rank]);

				//Update and build ghosts here
				updateLoadBalance();
				setPboundGhosts();

			}
			else
			{
				(*m_log) << " " << endl;
				(*m_log) << " Initial Parallel partition : " << endl;
				(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
				for(int ii=1; ii<m_nproc; ii++){
					(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
				}

				//empty ghosts
				m_octree.m_ghosts.clear();
				m_octree.m_sizeGhosts = 0;
				//compute new partition range globalidx
				uint64_t* newPartitionRangeGlobalidx = new uint64_t[m_nproc];
				for(int p = 0; p < m_nproc; ++p){
					newPartitionRangeGlobalidx[p] = 0;
					for(int pp = 0; pp <= p; ++pp)
						newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
					--newPartitionRangeGlobalidx[p];
				}

				//find resident octants local offset lastHead(lh) and firstTail(ft)
				int32_t lh,ft;
				if(m_rank == 0)
					lh = -1;
				else{
					lh = (int32_t)(newPartitionRangeGlobalidx[m_rank-1] + 1 - m_partitionRangeGlobalIdx[m_rank-1] - 1 - 1);
				}
				if(lh < 0)
					lh = - 1;
				else if(lh > m_octree.m_octants.size() - 1)
					lh = m_octree.m_octants.size() - 1;

				if(m_rank == m_nproc - 1)
					ft = m_octree.m_octants.size();
				else if(m_rank == 0)
					ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] + 1);
				else{
					ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1]);
				}
				if(ft > (int32_t)(m_octree.m_octants.size() - 1))
					ft = m_octree.m_octants.size();
				else if(ft < 0)
					ft = 0;

				//compute size Head and size Tail
				uint32_t headSize = (uint32_t)(lh + 1);
				uint32_t tailSize = (uint32_t)(m_octree.m_octants.size() - ft);
				uint32_t headOffset = headSize;
				uint32_t tailOffset = tailSize;

				//build send buffers
				std::map<int,CommBuffer> sendBuffers;

				//Compute first predecessor and first successor to send buffers to
				int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
				int64_t globalLastHead = (int64_t) lh;
				int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
				int firstPredecessor = -1;
				int firstSuccessor = m_nproc;
				if(m_rank != 0){
					firstOctantGlobalIdx = (int64_t)(m_partitionRangeGlobalIdx[m_rank-1] + 1);
					globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
					globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
					for(int pre = m_rank - 1; pre >=0; --pre){
						if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
							firstPredecessor = pre;
					}
					for(int post = m_rank + 1; post < m_nproc; ++post){
						if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
							firstSuccessor = post;
					}
				}
				else if(m_rank == 0){
					firstSuccessor = 1;
				}
				MPI_Barrier(m_comm); //da spostare prima della prima comunicazione

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

							int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
							limits[1] = (uint32_t)lh + 1;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
								//PACK octants from 0 to lh in sendBuffer[p]
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);

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
							int buffSize = nofElementsFromSuccessiveToPrevious * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromSuccessiveToPrevious,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
							limits[1] = (uint32_t)lh + 1;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
								//pack octants from lh - partition[p] to lh
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
					for(int p = firstSuccessor; p < m_nproc; ++p){
						if(tailSize < partition[p]){
							nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
							if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
								nofElementsFromPreviousToSuccessive = tailSize;

							int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);
							//USE BUFFER POS

							limits[0] = ft;
							limits[1] = ft + nofElementsFromPreviousToSuccessive;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
								//PACK octants from ft to octantsSize-1
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
							int buffSize = nofElementsFromPreviousToSuccessive * (int)ceil((double)m_global.m_octantBytes / (double)(CHAR_BIT/8));
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
							sendBuffers[p] = CommBuffer(buffSize,'a',m_comm);
							//store the number of octants at the beginning of the buffer
							MPI_Pack(&nofElementsFromPreviousToSuccessive,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,sendBuffers[p].m_commBufferSize,&sendBuffers[p].m_pos,m_comm);

							limits[0] = ft;
							limits[1] = endOctants + 1;
							std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);
							m_sentIdx.insert(procLimits);

							for(uint32_t i = ft; i <= endOctants; ++i ){
								//PACK octants from ft to ft + partition[p] -1
								//const Class_Octant<2> & octant = m_octree.m_octants[i];
								const Octant & octant = m_octree.m_octants[i];
								x = octant.getX();
								y = octant.getY();
								z = octant.getZ();
								l = octant.getLevel();
								m = octant.getMarker();
								for(int j = 0; j < 17; ++j)
									info[j] = octant.m_info[j];
								m_errorFlag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								m_errorFlag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
								for(int j = 0; j < 17; ++j){
									MPI_Pack(&info[j],1,MPI_C_BOOL,sendBuffers[p].m_commBuffer,buffSize,&sendBuffers[p].m_pos,m_comm);
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
				std::vector<Array> recvs(m_nproc);
				recvs[m_rank] = Array((uint32_t)sendBuffers.size()+1,-1);
				recvs[m_rank].m_array[0] = m_rank;
				int counter = 1;
				std::map<int,CommBuffer>::iterator sitend = sendBuffers.end();
				for(std::map<int,CommBuffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
					recvs[m_rank].m_array[counter] = sit->first;
					++counter;
				}
				int* nofRecvsPerProc = new int[m_nproc];
				m_errorFlag = MPI_Allgather(&recvs[m_rank].m_arraySize,1,MPI_INT,nofRecvsPerProc,1,MPI_INT,m_comm);
				int globalRecvsBuffSize = 0;
				int* displays = new int[m_nproc];
				for(int pp = 0; pp < m_nproc; ++pp){
					displays[pp] = 0;
					globalRecvsBuffSize += nofRecvsPerProc[pp];
					for(int ppp = 0; ppp < pp; ++ppp){
						displays[pp] += nofRecvsPerProc[ppp];
					}
				}
				int* globalRecvsBuff = new int[globalRecvsBuffSize];
				m_errorFlag = MPI_Allgatherv(recvs[m_rank].m_array,recvs[m_rank].m_arraySize,MPI_INT,globalRecvsBuff,nofRecvsPerProc,displays,MPI_INT,m_comm);

				std::vector<std::set<int> > sendersPerProc(m_nproc);
				for(int pin = 0; pin < m_nproc; ++pin){
					for(int k = displays[pin]+1; k < displays[pin] + nofRecvsPerProc[pin]; ++k){
						sendersPerProc[globalRecvsBuff[k]].insert(globalRecvsBuff[displays[pin]]);
					}
				}

				//Communicate Octants (size)
				MPI_Request* req = new MPI_Request[sendBuffers.size()+sendersPerProc[m_rank].size()];
				MPI_Status* stats = new MPI_Status[sendBuffers.size()+sendersPerProc[m_rank].size()];
				int nReq = 0;
				std::map<int,int> recvBufferSizePerProc;
				std::set<int>::iterator senditend = sendersPerProc[m_rank].end();
				for(std::set<int>::iterator sendit = sendersPerProc[m_rank].begin(); sendit != senditend; ++sendit){
					recvBufferSizePerProc[*sendit] = 0;
					m_errorFlag = MPI_Irecv(&recvBufferSizePerProc[*sendit],1,MPI_UINT32_T,*sendit,m_rank,m_comm,&req[nReq]);
					++nReq;
				}
				std::map<int,CommBuffer>::reverse_iterator rsitend = sendBuffers.rend();
				for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
					m_errorFlag =  MPI_Isend(&rsit->second.m_commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,m_comm,&req[nReq]);
					++nReq;
				}
				MPI_Waitall(nReq,req,stats);

				//COMMUNICATE THE BUFFERS TO THE RECEIVERS
				//recvBuffers structure is declared and each buffer is initialized to the right size
				//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
				uint32_t nofNewHead = 0;
				uint32_t nofNewTail = 0;
				std::map<int,CommBuffer> recvBuffers;

				std::map<int,int>::iterator ritend = recvBufferSizePerProc.end();
				for(std::map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
					recvBuffers[rit->first] = CommBuffer(rit->second,'a',m_comm);
				}

				nReq = 0;
				for(std::set<int>::iterator sendit = sendersPerProc[m_rank].begin(); sendit != senditend; ++sendit){
					m_errorFlag = MPI_Irecv(recvBuffers[*sendit].m_commBuffer,recvBuffers[*sendit].m_commBufferSize,MPI_PACKED,*sendit,m_rank,m_comm,&req[nReq]);
					++nReq;
				}
				for(std::map<int,CommBuffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
					m_errorFlag =  MPI_Isend(rsit->second.m_commBuffer,rsit->second.m_commBufferSize,MPI_PACKED,rsit->first,rsit->first,m_comm,&req[nReq]);
					++nReq;
				}
				MPI_Waitall(nReq,req,stats);

				//Unpack number of octants per sender
				std::map<int,uint32_t> nofNewOverProcs;
				std::map<int,CommBuffer>::iterator rbitend = recvBuffers.end();
				for(std::map<int,CommBuffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
					uint32_t nofNewPerProc;
					MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&nofNewPerProc,1,MPI_UINT32_T,m_comm);
					nofNewOverProcs[rbit->first] = nofNewPerProc;
					if(rbit->first < m_rank)
						nofNewHead += nofNewPerProc;
					else if(rbit->first > m_rank)
						nofNewTail += nofNewPerProc;
				}

				//MOVE RESIDENT TO BEGIN IN OCTANTS
				uint32_t resEnd = m_octree.getNumOctants() - tailOffset;
				uint32_t nofResidents = resEnd - headOffset;
				uint32_t octCounter = 0;
				for(uint32_t i = headOffset; i < resEnd; ++i){
					m_octree.m_octants[octCounter] = m_octree.m_octants[i];
					//TODO move data - DONE
					userData.move(i,octCounter);
					++octCounter;
				}
				uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
				m_octree.m_octants.resize(newCounter, Octant(m_dim, m_global.m_maxLevel));
				userData.resize(newCounter);
				//MOVE RESIDENTS IN RIGHT POSITION
				uint32_t resCounter = nofNewHead + nofResidents - 1;
				for(uint32_t k = 0; k < nofResidents ; ++k){
					m_octree.m_octants[resCounter - k] = m_octree.m_octants[nofResidents - k - 1];
					//TODO move data - DON
					userData.move(nofResidents - k - 1,resCounter - k);
				}

				//UNPACK BUFFERS AND BUILD NEW OCTANTS
				newCounter = 0;
				bool jumpResident = false;

				for(std::map<int,CommBuffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
					//TODO change new octants counting, probably you have to communicate the number of news per proc
					uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];
					if(rbit->first > m_rank && !jumpResident){
						newCounter += nofResidents ;
						jumpResident = true;
					}
					for(int i = nofNewPerProc - 1; i >= 0; --i){
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&x,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&y,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&z,1,MPI_UINT32_T,m_comm);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&l,1,MPI_UINT8_T,m_comm);
						//m_octree.m_octants[newCounter] = Class_Octant<2>(l,x,y);
						m_octree.m_octants[newCounter] = Octant(m_dim,l,x,y,z,m_global.m_maxLevel);
						m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&m,1,MPI_INT8_T,m_comm);
						m_octree.m_octants[newCounter].setMarker(m);
						for(int j = 0; j < 17; ++j){
							m_errorFlag = MPI_Unpack(rbit->second.m_commBuffer,rbit->second.m_commBufferSize,&rbit->second.m_pos,&info[j],1,MPI_C_BOOL,m_comm);
							m_octree.m_octants[newCounter].m_info[j] = info[j];
						}
						//TODO Unpack data
						userData.scatter(rbit->second,newCounter);
						++newCounter;
					}
				}
				octvector(m_octree.m_octants).swap(m_octree.m_octants);

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

			//Write info of final partition on m_log
			(*m_log) << " " << endl;
			(*m_log) << " Final Parallel partition : " << endl;
			(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
			for(int ii=1; ii<m_nproc; ii++){
				(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
			}
			(*m_log) << " " << endl;
			(*m_log) << "---------------------------------------------" << endl;

		}
		else{
			(*m_log) << " " << endl;
			(*m_log) << " Serial partition : " << endl;
			(*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
			(*m_log) << " " << endl;
			(*m_log) << "---------------------------------------------" << endl;
		}
	}

#endif

	// =============================================================================== //


};

/*  @}  */

}

#endif /* __BITPIT_PARA_TREE_HPP__ */
