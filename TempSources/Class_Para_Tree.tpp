/*
 * Class_Para_Tree.tpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */


template<class UserDataComm>
void Class_Para_Tree::communicate(UserDataComm & userData){

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
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a');
		//store number of pborders from this proc at the begining
		MPI_Pack(&nofPbordersPerProc,1,MPI_INT,sendBuffers[key].commBuffer,sendBuffers[key].commBufferSize,&sendBuffers[key].pos,MPI_COMM_WORLD);

		//WRITE SEND BUFFERS
		for(size_t j = 0; j < nofPbordersPerProc; ++j){
			userData.gather(sendBuffers[key],pborders[j]);
		}
	}

	//Communicate Buffers Size
	MPI_Request req[sendBuffers.size()*2];
	MPI_Status stats[sendBuffers.size()*2];
	int nReq = 0;
	map<int,int> recvBufferSizePerProc;
	map<int,Class_Comm_Buffer>::iterator sitend = sendBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		recvBufferSizePerProc[sit->first] = 0;
		error_flag = MPI_Irecv(&recvBufferSizePerProc[sit->first],1,MPI_UINT32_T,sit->first,rank,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	map<int,Class_Comm_Buffer>::reverse_iterator rsitend = sendBuffers.rend();
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(&rsit->second.commBufferSize,1,MPI_UINT32_T,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//Communicate Buffers
	//uint32_t nofBytesOverProc = 0;
	map<int,Class_Comm_Buffer> recvBuffers;
	map<int,int>::iterator ritend = recvBufferSizePerProc.end();
	for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
		recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a');
	}
	nReq = 0;
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		//nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
		error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//READ RECEIVE BUFFERS
	int ghostOffset = 0;
	map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
	map<int,Class_Comm_Buffer>::iterator rbitbegin = recvBuffers.begin();
	for(map<int,Class_Comm_Buffer>::iterator rbit = rbitbegin; rbit != rbitend; ++rbit){
		int nofGhostFromThisProc = 0;
		MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&nofGhostFromThisProc,1,MPI_INT,MPI_COMM_WORLD);
		for(int k = 0; k < nofGhostFromThisProc; ++k){
			userData.scatter(rbit->second, k+ghostOffset);
		}
		ghostOffset += nofGhostFromThisProc;
	}


};

template<class UserDataComm>
void Class_Para_Tree::loadBalance(UserDataComm & userData){

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
		Class_Local_Tree::OctantsType::const_iterator first = octree.octants.begin() + stride;
		Class_Local_Tree::OctantsType::const_iterator last = first + partition[rank];
		typename UserDataComm::Data::iterator firstData = userData.data.begin() + stride;
		typename UserDataComm::Data::iterator lastData = firstData + partition[rank];
		octree.octants.assign(first, last);
		userData.data.assign(firstData,lastData);
		octree.octants.shrink_to_fit();
		userData.data.shrink_to_fit();
		first = octree.octants.end();
		last = octree.octants.end();

		//Update and build ghosts here
		updateLoadBalance();
		setPboundGhosts();

//		userData.ghostData.resize(octree.size_ghosts,0.0);
//		communicate(userData);

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
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
					//store the number of octants at the beginning of the buffer
					MPI_Pack(&headSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
					//USE BUFFER POS
					//int pos = 0;
					for(uint32_t i = 0; i <= lh; ++i){
						//PACK octants from 0 to lh in sendBuffer[p]
						const Class_Octant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						memcpy(info,octant.info,16);
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						for(int j = 0; j < 16; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);

						}
						//TODO call gather to pack user data - DONE
						userData.gather(sendBuffers[p],i);
					}
					break;
				}
				else{
					int buffSize = partition[p] * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
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
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
					//store the number of octants at the beginning of the buffer
					MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
					//USE BUFFER POS
					//int pos = 0;
					for(uint32_t i = lh - partition[p] + 1; i <= lh; ++i){
						//pack octants from lh - partition[p] to lh
						const Class_Octant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						memcpy(info,octant.info,16);
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						for(int j = 0; j < 16; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
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
					int buffSize = tailSize * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
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
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
					//store the number of octants at the beginning of the buffer
					MPI_Pack(&tailSize,1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
					//USE BUFFER POS
					//int pos = 0;
					for(uint32_t i = ft; i < octantsSize; ++i){
						//PACK octants from ft to octantsSize-1
						const Class_Octant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						memcpy(info,octant.info,16);
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						for(int j = 0; j < 16; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						}
						//TODO call gather to pack user data - DONE
						userData.gather(sendBuffers[p],i);
					}
					break;
				}
				else{
					uint32_t endOctants = ft + partition[p] - 1;
					int buffSize = partition[p] * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));
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
					sendBuffers[p] = Class_Comm_Buffer(buffSize,'a');
					//store the number of octants at the beginning of the buffer
					MPI_Pack(&partition[p],1,MPI_UINT32_T,sendBuffers[p].commBuffer,sendBuffers[p].commBufferSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
					//int pos = 0;
					for(uint32_t i = ft; i <= endOctants; ++i ){
						//PACK octants from ft to ft + partition[p] -1
						const Class_Octant & octant = octree.octants[i];
						x = octant.getX();
						y = octant.getY();
						z = octant.getZ();
						l = octant.getLevel();
						m = octant.getMarker();
						memcpy(info,octant.info,16);
						error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
						for(int j = 0; j < 16; ++j){
							MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[p].commBuffer,buffSize,&sendBuffers[p].pos,MPI_COMM_WORLD);
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
//			uint32_t nofNewPerProc = (uint32_t)(rit->second / (uint32_t)ceil((double)octantBytes / (double)(CHAR_BIT/8)));
//			if(rit->first < rank)
//				nofNewHead += nofNewPerProc;
//			else if(rit->first > rank)
//				nofNewTail += nofNewPerProc;
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

		//Unpack number of octants per sender
		map<int,uint32_t> nofNewOverProcs;
		map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
		for(map<int,Class_Comm_Buffer>::iterator rbit = recvBuffers.begin(); rbit != rbitend; ++rbit){
			uint32_t nofNewPerProc;
			MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&nofNewPerProc,1,MPI_UINT32_T,MPI_COMM_WORLD);
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
		userData.data.resize(newCounter);
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
			uint32_t nofNewPerProc = nofNewOverProcs[rbit->first];//(uint32_t)(rbit->second.commBufferSize / (uint32_t)ceil((double)octantBytes / (double)(CHAR_BIT/8)));
			//int pos = 0;
			if(rbit->first > rank && !jumpResident){
				newCounter += nofResidents ;
				jumpResident = true;
			}
			for(int i = nofNewPerProc - 1; i >= 0; --i){
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&x,1,MPI_UINT32_T,MPI_COMM_WORLD);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&y,1,MPI_UINT32_T,MPI_COMM_WORLD);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&z,1,MPI_UINT32_T,MPI_COMM_WORLD);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&l,1,MPI_UINT8_T,MPI_COMM_WORLD);
				octree.octants[newCounter] = Class_Octant(l,x,y,z);
				error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&m,1,MPI_INT8_T,MPI_COMM_WORLD);
				octree.octants[newCounter].setMarker(m);
				for(int j = 0; j < 16; ++j){
					error_flag = MPI_Unpack(rbit->second.commBuffer,rbit->second.commBufferSize,&rbit->second.pos,&info[j],1,MPI::BOOL,MPI_COMM_WORLD);
					octree.octants[newCounter].info[j] = info[j];
				}
				//TODO Unpack data
				userData.scatter(rbit->second,newCounter);
				++newCounter;
			}
		}
		octree.octants.shrink_to_fit();
		userData.data.shrink_to_fit();
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


}




