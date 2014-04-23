/*
 * Class_Para_Tree_3D.tpp
 *
 *  Created on: 23/apr/2014
 *      Author: Marco Cisternino
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

template<>
class Class_Para_Tree<3>{
	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant<3> > 		OctantsType;
	typedef vector<uint32_t>			u32vector;
	typedef vector<double>				dvector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;
	typedef vector<vector<double>	>	dvector2D;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	//undistributed members
	uint64_t* partition_first_desc; 			//global array containing position of the first possible octant in each processor
	uint64_t* partition_last_desc; 				//global array containing position of the last possible octant in each processor
	uint64_t* partition_range_globalidx;	 	//global array containing global index of the last existing octant in each processor
	uint64_t global_num_octants;   				// global number of octants in the parallel octree
	map<int,vector<uint32_t> > bordersPerProc;	//local indices of border octants per process
	int nproc;
	uint8_t max_depth;							// global max existing level in the parallel octree

	//distributed members
	int rank;
	Class_Local_Tree<3> octree;					// local tree in each processor

	//auxiliary members
	int error_flag;								// MPI error flag
	bool serial;								// 1 if the octree is the same on each processor, 0 if the octree is distributed

	//map member
	Class_Map<3> trans;

	// connectivity
	dvector2D					nodes;				// Local vector of nodes (x,y,z) ordered with Morton Number
	u32vector2D					connectivity;		// Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	// The nodes are stored as index of vector nodes
	dvector2D					ghostsnodes;		// Local vector of ghosts nodes (x,y,z) ordered with Morton Number
	u32vector2D					ghostsconnectivity;	// Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	// The nodes are stored as index of vector nodes

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Para_Tree(){
		serial = true;
		error_flag = 0;
		max_depth = 0;
		global_num_octants = octree.getNumOctants();
		error_flag = MPI_Comm_size(MPI_COMM_WORLD,&nproc);
		error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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
		if(rank==0){
			int sysError = system("rm PABLO.log");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		writeLog("---------------------------------------------");
		writeLog("- PABLO PArallel Balanced Linear Octree -");
		writeLog("---------------------------------------------");
		writeLog(" ");
		writeLog("---------------------------------------------");
		writeLog(" Number of proc		:	" + to_string(nproc));
		writeLog(" Dimension		:	" + to_string(3));
		writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL_3D));
		writeLog("---------------------------------------------");
		writeLog(" ");

	};
	Class_Para_Tree(double & X, double & Y, double & Z, double & L){
		serial = true;
		error_flag = 0;
		max_depth = 0;
		global_num_octants = octree.getNumOctants();
		error_flag = MPI_Comm_size(MPI_COMM_WORLD,&nproc);
		error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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
		if(rank==0){
			int sysError = system("rm PABLO.log");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		writeLog("---------------------------------------------");
		writeLog("- PABLO PArallel Balanced Linear Octree -");
		writeLog("---------------------------------------------");
		writeLog(" ");
		writeLog("---------------------------------------------");
		writeLog(" Number of proc		:	" + to_string(nproc));
		writeLog(" Dimension		:	" + to_string(3));
		writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL_3D));
		writeLog(" Domain Origin		:	" + to_string(X));
		writeLog("				" + to_string(Y));
		writeLog("				" + to_string(Z));
		writeLog(" Domain Size		:	" + to_string(L));
		writeLog("---------------------------------------------");
		writeLog(" ");

	};
	~Class_Para_Tree(){
		writeLog("---------------------------------------------");
		writeLog("--------------- R.I.P. PABLO ----------------");
		writeLog("---------------------------------------------");
		writeLog("---------------------------------------------");
	};
	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	void computePartition(uint32_t* partition){ 	// compute octant partition giving the same number of octant to each process and redistributing the reminder

	};
	//	void computePartition(uint32_t* partition,  // compute octant partition giving almost the same number of octant to each process
	//						uint8_t & level);   // with complete families contained in octants of n "level" over the leaf in each process
	//	void updateLoadBalance();					//update Class_Para_Tree members after a load balance
	//	void setPboundGhosts(); 			 		// set pbound and build ghosts after static load balance
	void loadBalance(){								//assign the octants to the processes following a computed partition

		//Write info on log
		writeLog("---------------------------------------------");
		writeLog(" LOAD BALANCE ");

		uint32_t* partition = new uint32_t [nproc];
		computePartition(partition);
		if(serial)
		{
			writeLog(" ");
			writeLog(" Initial Serial distribution : ");
			for(int ii=0; ii<nproc; ii++){
				writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]+1));
			}

			uint32_t stride = 0;
			for(int i = 0; i < rank; ++i)
				stride += partition[i];
			Class_Local_Tree<3>::OctantsType::const_iterator first = octree.octants.begin() + stride;
			Class_Local_Tree<3>::OctantsType::const_iterator last = first + partition[rank];
			octree.octants.assign(first, last);
			octree.octants.shrink_to_fit();
			first = octree.octants.end();
			last = octree.octants.end();

			//Update and ghosts here
			updateLoadBalance();
			setPboundGhosts();

		}
		else
		{
			writeLog(" ");
			writeLog(" Initial Parallel partition : ");
			writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
			for(int ii=1; ii<nproc; ii++){
				writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
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
			MPI_Barrier(MPI_COMM_WORLD); //da spostare prima della prima comunicazione

			uint32_t x,y,z;
			uint8_t l;
			int8_t m;
			bool info[16];
			//build send buffers from Head
			if(headSize != 0){
				for(int p = firstPredecessor; p >= 0; --p){
					if(headSize <=partition[p]){
						int buffSize = headSize * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
						int pos = 0;
						for(uint32_t i = 0; i <= lh; ++i){
							//PACK octants from 0 to lh in sendBuffer[p]
							const Class_Octant<3> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							memcpy(info,octant.info,16);
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							for(int j = 0; j < 16; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
						int pos = 0;
						for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
							//pack octants from lh - partition[p] to lh
							const Class_Octant<3> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							memcpy(info,octant.info,16);
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							for(int j = 0; j < 16; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
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
						int buffSize = tailSize * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
						int pos = 0;
						uint32_t octantsSize = (uint32_t)octree.octants.size();
						for(uint32_t i = ft; i < octantsSize; ++i){
							//PACK octants from ft to octantsSize-1
							const Class_Octant<3> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							memcpy(info,octant.info,16);
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							for(int j = 0; j < 16; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							}
						}
						break;
					}
					else{
						int buffSize = partition[p] * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
						sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
						uint32_t endOctants = ft + partition[p] - 1;
						int pos = 0;
						for(uint32_t i = ft; i <= endOctants; ++i ){
							//PACK octants from ft to ft + partition[p] -1
							const Class_Octant<3> & octant = octree.octants[i];
							x = octant.getX();
							y = octant.getY();
							z = octant.getZ();
							l = octant.getLevel();
							m = octant.getMarker();
							memcpy(info,octant.info,16);
							error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
							for(int j = 0; j < 16; ++j){
								MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
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
			int nofRecvsPerProc[nproc];
			error_flag = MPI_Allgather(&recvs[rank].arraySize,1,MPI_INT,nofRecvsPerProc,1,MPI_INT,MPI_COMM_WORLD);
			int globalRecvsBuffSize = 0;
			int displays[nproc];
			for(int pp = 0; pp < nproc; ++pp){
				displays[pp] = 0;
				globalRecvsBuffSize += nofRecvsPerProc[pp];
				for(int ppp = 0; ppp < pp; ++ppp){
					displays[pp] += nofRecvsPerProc[ppp];
				}
			}
			int globalRecvsBuff[globalRecvsBuffSize];
			error_flag = MPI_Allgatherv(recvs[rank].array,recvs[rank].arraySize,MPI_INT,globalRecvsBuff,nofRecvsPerProc,displays,MPI_INT,MPI_COMM_WORLD);

			vector<set<int> > sendersPerProc(nproc);
			for(int pin = 0; pin < nproc; ++pin){
				for(int k = displays[pin]+1; k < displays[pin] + nofRecvsPerProc[pin]; ++k){
					sendersPerProc[globalRecvsBuff[k]].insert(globalRecvsBuff[displays[pin]]);
				}
			}

			//Communicate Octants (size)
			MPI_Request req[sendBuffers.size()+sendersPerProc[rank].size()];
			MPI_Status stats[sendBuffers.size()+sendersPerProc[rank].size()];
			int nReq = 0;
			map<int,int> recvBufferSizePerProc;
			set<int>::iterator senditend = sendersPerProc[rank].end();
			for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
				recvBufferSizePerProc[*sendit] = 0;
				error_flag = MPI_Irecv(&recvBufferSizePerProc[*sendit],1,MPI_UINT32_T,*sendit,rank,MPI_COMM_WORLD,&req[nReq]);
				++nReq;
			}
			map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
			for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
				error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
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
				recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a');
				uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)octantBytes / (double)(CHAR_BIT/8)));
				if(rit->first < rank)
					nofNewHead += nofNewPerProc;
				else if(rit->first > rank)
					nofNewTail += nofNewPerProc;
			}
			nReq = 0;
			for(set<int>::iterator sendit = sendersPerProc[rank].begin(); sendit != senditend; ++sendit){
				//nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
				error_flag = MPI_Irecv(recvBuffers[*sendit].commBuffer,recvBuffers[*sendit].commBufferSize,MPI_PACKED,*sendit,rank,MPI_COMM_WORLD,&req[nReq]);
				++nReq;
			}
			for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
				error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
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
				uint32_t nofNewPerProc = (uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)octantBytes / (double)(CHAR_BIT/8)));
				int pos = 0;
				if(rbit->first > rank && !jumpResident){
					newCounter += nofResidents ;
					jumpResident = true;
				}
				for(int i = nofNewPerProc - 1; i >= 0; --i){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,MPI_COMM_WORLD);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,MPI_COMM_WORLD);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&z,1,MPI_UINT32_T,MPI_COMM_WORLD);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,MPI_COMM_WORLD);
					octree.octants[newCounter] = Class_Octant(l,x,y,z);
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,MPI_COMM_WORLD);
					octree.octants[newCounter].setMarker(m);
					for(int j = 0; j < 16; ++j){
						error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,MPI_COMM_WORLD);
						octree.octants[newCounter].info[j] = info[j];
					}
					++newCounter;
				}
			}
			octree.octants.shrink_to_fit();
			octree.pborders.clear();

			delete [] newPartitionRangeGlobalidx;
			newPartitionRangeGlobalidx = NULL;

			//Update and ghosts here
			updateLoadBalance();
			setPboundGhosts();

		}
		delete [] partition;
		partition = NULL;

		//Write info of final partition on log
		writeLog(" ");
		writeLog(" Final Parallel partition : ");
		writeLog(" Octants for proc	"+ to_string(0)+"	:	" + to_string(partition_range_globalidx[0]+1));
		for(int ii=1; ii<nproc; ii++){
			writeLog(" Octants for proc	"+ to_string(ii)+"	:	" + to_string(partition_range_globalidx[ii]-partition_range_globalidx[ii-1]));
		}
		writeLog(" ");
		writeLog("---------------------------------------------");

	};
	//void loadBalance(uint8_t & level);			//assign the octants to the processes following a computed partition with complete families contained in octants of n "level" over the leaf in each process
	//	template<class UserDataComm>
	//	void loadBalance(UserDataComm & userData);
	//	bool adapt();  								//call refine and coarse on the local tree
	//	bool adapt(u32vector & mapidx);  			//call refine and coarse on the local tree
	//												// mapidx[i] = index in old octants vector of the i-th octant (index of father or first child if octant is new after refine or coarse)
	//	bool adapt(u32vector & mapidx,
	//				u32vector & mapinters_int,
	//				u32vector & mapinters_ghost,
	//				u32vector & mapinters_bord);
	//	void updateAdapt();							//update Class_Para_Tree members after a refine and/or coarse
	//	void updateAfterCoarse();					//update Class_Para_Tree members and delete overlapping octants after a coarse
	//	void updateAfterCoarse(u32vector & mapidx);	//update Class_Para_Tree members and delete overlapping octants after a coarse
	//	int findOwner(const uint64_t & morton);		// given the morton of an octant it finds the process owning that octant
	//	void commMarker();							// communicates marker of ghosts
	//	void balance21();							// 2:1 balancing of parallel octree
	//	template<class UserDataComm>
	//	void communicate(UserDataComm & userData);
	//
	//	void computeConnectivity();						// Computes nodes vector and connectivity of octants of local tree
	//	void clearConnectivity();						// Clear nodes vector and connectivity of octants of local tree
	//	void updateConnectivity();						// Updates nodes vector and connectivity of octants of local tree
	//	void computeghostsConnectivity();				// Computes ghosts nodes vector and connectivity of ghosts octants of local tree
	//	void clearghostsConnectivity();					// Clear ghosts nodes vector and connectivity of ghosts octants of local tree
	//	void updateghostsConnectivity();					// Update ghosts nodes vector and connectivity of ghosts octants of local tree
	//
	//	// --------------------------------------------------------------------------------------------- //
	//	// Basic Get Methods --------------------------------------------------------------------------- //
	//
	//public:
	//	double			getX(Class_Octant* const oct);
	//	double			getY(Class_Octant* const oct);
	//	double			getZ(Class_Octant* const oct);
	//	double			getSize(Class_Octant* const oct);		// Get the size of octant if mapped in hypercube
	//	double			getArea(Class_Octant* const oct);		// Get the face area of octant
	//	double			getVolume(Class_Octant* const oct);		// Get the volume of octant
	//	void			getCenter(Class_Octant* oct, 			// Get a vector of DIM with the coordinates of the center of octant
	//					dvector & center);
	//	void			getNodes(Class_Octant* oct, 			// Get a vector of vector (size [nnodes][DIM]) with the nodes of octant
	//					dvector2D & nodes);
	//	void			getNormal(Class_Octant* oct, 			// Get a vector of vector (size [DIM]) with the normal of the iface
	//					uint8_t & iface,
	//					dvector & normal);
	//
	//	Class_Octant*	getPointOwner(dvector & point);			// Get the pointer to the octant owner of an input point
	//															// (vector<double> with x,y,z). If the point is out of process
	//															// return NULL.
	//
	//	double			getSize(Class_Intersection* const inter);		// Get the size of intersection if mapped in hypercube
	//	double			getArea(Class_Intersection* const inter);		// Get the face area of intersection
	//	void			getCenter(Class_Intersection* const inter, 		// Get a vector of DIM with the coordinates of the center of intersection
	//					dvector & center);
	//	void			getNodes(Class_Intersection* const inter, 		// Get a vector of vector (size [nnodes][DIM]) with the nodes of intersection
	//					dvector2D & nodes);
	//	void			getNormal(Class_Intersection* const inter, 		// Get a vector of vector (size [DIM]) with the normal of the intersection
	//					dvector & normal);


};


