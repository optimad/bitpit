/*
 * Class_Para_Tree.cpp
 *
 *  Created on: 12/feb/2014
 *      Author: Marco Cisternino
 */

#include "Class_Para_Tree.hpp"

Class_Para_Tree::Class_Para_Tree() {
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
		system("rm PABLO.log");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	writeLog("---------------------------------------------");
	writeLog("- PABLO PArallel Balanced Linear Octree -");
	writeLog("---------------------------------------------");
	writeLog(" ");
	writeLog("---------------------------------------------");
	writeLog(" Number of proc		:	" + to_string(nproc));
	writeLog(" Dimension		:	" + to_string(DIM));
	writeLog(" Max allowed level	:	" + to_string(MAX_LEVEL));
	writeLog("---------------------------------------------");
	writeLog(" ");

}

Class_Para_Tree::~Class_Para_Tree() {
	writeLog("---------------------------------------------");
	writeLog("--------------- R.I.P. PABLO ----------------");
	writeLog("---------------------------------------------");
	writeLog("---------------------------------------------");
	//delete [] partition_last_desc;
	//delete [] partition_range_globalidx;
}


void Class_Para_Tree::computePartition(uint32_t* partition) {
	uint32_t division_result = 0;
	uint32_t remind = 0;
	division_result = uint32_t(global_num_octants/(uint64_t)nproc);
	remind = (uint32_t)(global_num_octants%(uint64_t)nproc);
	for(int i = 0; i < nproc; ++i)
		if(i<remind)
			partition[i] = division_result + 1;
		else
			partition[i] = division_result;

}

void Class_Para_Tree::updateLoadBalance() {
	octree.updateLocalMaxDepth();
	//update partition_range_globalidx
	uint64_t rbuff [nproc];
	uint64_t local_num_octants = octree.getNumOctants();
	error_flag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,&rbuff,1,MPI_UINT64_T,MPI_COMM_WORLD);
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
	error_flag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,partition_last_desc,1,MPI_UINT64_T,MPI_COMM_WORLD);
	uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
	error_flag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,partition_first_desc,1,MPI_UINT64_T,MPI_COMM_WORLD);
	serial = false;
}

void Class_Para_Tree::loadBalance(){

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
						const Class_Octant & octant = octree.octants[i];
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
						const Class_Octant & octant = octree.octants[i];
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
						const Class_Octant & octant = octree.octants[i];
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
					for(uint32_t i = ft; i < endOctants; ++i ){
						//PACK octants from ft to ft + partition[p] -1
						const Class_Octant & octant = octree.octants[i];
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


}

void Class_Para_Tree::updateAdapt() {
	if(serial)
	{
		max_depth = octree.local_max_depth;
		global_num_octants = octree.getNumOctants();
		for(int p = 0; p < nproc; ++p){
			partition_range_globalidx[p] = global_num_octants - 1;
			//partition_range_position[p] = octree.last_desc.computeMorton();
		}
	}
	else
	{
		//update max_depth
		error_flag = MPI_Allreduce(&octree.local_max_depth,&max_depth,1,MPI_UINT8_T,MPI_MAX,MPI_COMM_WORLD);
		//update global_num_octants
		uint64_t local_num_octants = octree.getNumOctants();
		error_flag = MPI_Allreduce(&local_num_octants,&global_num_octants,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
		//update partition_range_globalidx
		uint64_t rbuff [nproc];
		error_flag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,rbuff,1,MPI_UINT64_T,MPI_COMM_WORLD);
		for(int p = 0; p < nproc; ++p){
			partition_range_globalidx[p] = 0;
			for(int pp = 0; pp <=p; ++pp)
				partition_range_globalidx[p] += rbuff[pp];
			--partition_range_globalidx[p];
		}
		//update partition_range_position
		uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
		error_flag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,partition_last_desc,1,MPI_UINT64_T,MPI_COMM_WORLD);
		uint64_t firstDescMorton = octree.getFirstDesc().computeMorton();
		error_flag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,partition_first_desc,1,MPI_UINT64_T,MPI_COMM_WORLD);
	}
}

void Class_Para_Tree::updateAfterCoarse(){
	if(serial){
		updateAdapt();
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
}

void Class_Para_Tree::adapt() {
	if(serial){
		writeLog("---------------------------------------------");
		writeLog(" ADAPT (Refine/Coarse)");
		writeLog(" ");
		writeLog(" Initial Number of octants	:	" + to_string(octree.getNumOctants()));
		while(octree.refine());
		writeLog(" Number of octants after Refine	:	" + to_string(octree.getNumOctants()));
		while(octree.coarse());
		writeLog(" Number of octants after Coarse	:	" + to_string(octree.getNumOctants()));
		updateAdapt();
		writeLog(" ");
		writeLog("---------------------------------------------");
	}
	else{
		writeLog("---------------------------------------------");
		writeLog(" ADAPT (Refine/Coarse)");
		writeLog(" ");
		writeLog(" Initial Number of octants	:	" + to_string(global_num_octants));
		updateAdapt();			// Togliere se non necessario
		setPboundGhosts();		// Togliere se non necessario
		while(octree.refine());
		updateAdapt();
		setPboundGhosts();
		writeLog(" Number of octants after Refine	:	" + to_string(global_num_octants));
		while(octree.coarse());
		writeLog("coarse done ");
		updateAfterCoarse();
		writeLog(" Number of octants after Coarse	:	" + to_string(global_num_octants));
		writeLog(" ");
		setPboundGhosts();
		writeLog("---------------------------------------------");
	}
}

void Class_Para_Tree::buildGhosts() {

}

int Class_Para_Tree::findOwner(const uint64_t & morton) {
	int p = -1;
	int length = nproc;
	int beg = 0;
	int end = nproc -1;
	int seed = nproc/2;
	while(beg != end){
		if(morton <= partition_last_desc[seed]){
			end = seed;
//			length = seed + 1;
			if(morton > partition_last_desc[seed-1])
				beg = seed;
		}
		else{
			beg = seed;
			if(morton <= partition_last_desc[seed+1])
				beg = seed + 1;
	//	length = end - seed -1;
		}
		length = end - beg;
		seed = beg + length/2;
	}
	p = beg;
	return p;
}

void Class_Para_Tree::setPboundGhosts() {
	//BUILD BORDER OCTANT INDECES VECTOR (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//find local octants to be sent as ghost to the right processes
	//it visits the local octants building virtual neighbors on each octant face
	//find the owner of these virtual neighbor and build a map (process,border octants)
	//this map contains the local octants as ghosts for neighbor processes
	if(octree.pborders.size() == 0){
		Class_Local_Tree::OctantsType::iterator end = octree.octants.end();
		Class_Local_Tree::OctantsType::iterator begin = octree.octants.begin();
		bordersPerProc.clear();
		for(Class_Local_Tree::OctantsType::iterator it = begin; it != end; ++it){
			set<int> procs;
			for(uint8_t i = 0; i < nface; ++i){
				if(it->getBound(i) == false){
					uint32_t virtualNeighborsSize = 0;
					uint64_t* virtualNeighbors = it->computeVirtualMorton(i,max_depth,virtualNeighborsSize);
					uint32_t maxDelta = virtualNeighborsSize/2;
					for(int j = 0; j <= maxDelta; ++j){
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
//						if(pBegin == pEnd || pBegin == pEnd - 1)
//							break;
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
	}
	else{
		Class_Local_Tree::u32vector::iterator end = octree.pborders.end();
		Class_Local_Tree::u32vector::iterator begin = octree.pborders.begin();
		bordersPerProc.clear();
		for(Class_Local_Tree::u32vector::iterator it = begin; it != end; ++it){
			Class_Octant & oct = octree.octants[*it];
			set<int> procs;
			for(uint8_t i = 0; i < nface; ++i){
				if(oct.getBound(i) == false){
					uint32_t virtualNeighborsSize = 0;
					uint64_t* virtualNeighbors = oct.computeVirtualMorton(i,max_depth,virtualNeighborsSize);
					uint32_t maxDelta = virtualNeighborsSize/2;
					for(int j = 0; j <= maxDelta; ++j){
						int pBegin = findOwner(virtualNeighbors[j]);
						int pEnd = findOwner(virtualNeighbors[virtualNeighborsSize - 1 - j]);
						procs.insert(pBegin);
						procs.insert(pEnd);
						if(pBegin == pEnd || pBegin == pEnd - 1)
							break;
					}
				}
			}
			set<int>::iterator pitend = procs.end();
			for(set<int>::iterator pit = procs.begin(); pit != pitend; ++pit){
				int p = *pit;
				if(p != rank){
					//TODO better reserve to avoid if
					bordersPerProc[p].push_back(*it);
					vector<uint32_t> & bordersSingleProc = bordersPerProc[p];
					if(bordersSingleProc.capacity() - bordersSingleProc.size() < 2)
						bordersSingleProc.reserve(2*bordersSingleProc.size());
				}
			}
		}

	}
	//PACK (mpi) BORDER OCTANTS IN CHAR BUFFERS WITH SIZE (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//it visits every element in bordersPerProc (one for every neighbor proc)
	//for every element it visits the border octants it contains and pack them in a new structure, sendBuffers
	//this map has an entry Class_Comm_Buffer for every proc containing the size in bytes of the buffer and the octants
	//to be sent to that proc packed in a char* buffer
	uint32_t x,y,z;
	uint8_t l;
	int8_t m;
	bool info[16];
	map<int,Class_Comm_Buffer> sendBuffers;
	map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
	uint32_t pbordersOversize = 0;
	for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
		pbordersOversize += bit->second.size();
		int buffSize = bit->second.size() * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));// + (int)ceil((double)sizeof(int)/(double)(CHAR_BIT/8));
		int key = bit->first;
		const vector<uint32_t> & value = bit->second;
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a');
		int pos = 0;
		int nofBorders = value.size();
		for(int i = 0; i < nofBorders; ++i){
			//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
			const Class_Octant & octant = octree.octants[value[i]];
			x = octant.getX();
			y = octant.getY();
			z = octant.getZ();
			l = octant.getLevel();
			m = octant.getMarker();
			memcpy(info,octant.info,16);
			error_flag = MPI_Pack(&x,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			error_flag = MPI_Pack(&y,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			error_flag = MPI_Pack(&z,1,MPI_UINT32_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			error_flag = MPI_Pack(&l,1,MPI_UINT8_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			error_flag = MPI_Pack(&m,1,MPI_INT8_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			for(int j = 0; j < 16; ++j){
				MPI_Pack(&info[j],1,MPI::BOOL,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			}
		}
	}

	//Build pborders
	octree.pborders.clear();
//	octree.pborders.reserve(pbordersOversize);
//	for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
//		set_union(bit->second.begin(),bit->second.end(),octree.pborders.begin(),octree.pborders.end(),octree.pborders.begin());
//	}

	//COMMUNICATE THE SIZE OF BUFFER TO THE RECEIVERS
	//the size of every borders buffer is communicated to the right process in order to build the receive buffer
	//and stored in the recvBufferSizePerProc structure
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

	//COMMUNICATE THE BUFFERS TO THE RECEIVERS
	//recvBuffers structure is declared and each buffer is initialized to the right size
	//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
	//at the same time every process compute the size in bytes of all the ghost octants
	uint32_t nofBytesOverProc = 0;
	map<int,Class_Comm_Buffer> recvBuffers;
	map<int,int>::iterator ritend = recvBufferSizePerProc.end();
	for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
		recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a');
	}
	nReq = 0;
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
		error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

	//COMPUTE GHOSTS SIZE IN BYTES
	//number of ghosts in every process is obtained through the size in bytes of the single octant
	//and ghost vector in local tree is resized
	uint32_t nofGhosts = nofBytesOverProc / (uint32_t)octantBytes;
	octree.size_ghosts = nofGhosts;
	octree.ghosts.resize(nofGhosts);

	//UNPACK BUFFERS AND BUILD GHOSTS CONTAINER OF CLASS_LOCAL_TREE
	//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked octant by octant.
	//every ghost octant is built and put in the ghost vector
	uint32_t ghostCounter = 0;
	map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
		int pos = 0;
		int nofGhostsPerProc = int(rrit->second.commBufferSize / (uint32_t) octantBytes);
		for(int i = 0; i < nofGhostsPerProc; ++i){
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&x,1,MPI_UINT32_T,MPI_COMM_WORLD);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&y,1,MPI_UINT32_T,MPI_COMM_WORLD);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&z,1,MPI_UINT32_T,MPI_COMM_WORLD);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&l,1,MPI_UINT8_T,MPI_COMM_WORLD);
			octree.ghosts[ghostCounter] = Class_Octant(l,x,y,z);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&m,1,MPI_INT8_T,MPI_COMM_WORLD);
			octree.ghosts[ghostCounter].setMarker(m);
			for(int j = 0; j < 16; ++j){
				error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&info[j],1,MPI::BOOL,MPI_COMM_WORLD);
				octree.ghosts[ghostCounter].info[j] = info[j];
			}
//			octree.ghosts[ghostCounter].info[15] = true;
			++ghostCounter;
		}
	}
	recvBuffers.clear();
	sendBuffers.clear();
	recvBufferSizePerProc.clear();
}

//==============================================================

void Class_Para_Tree::commMarker() {
	// borderPerProcs has to be built

	//PACK (mpi) LEVEL AND MARKER OF BORDER OCTANTS IN CHAR BUFFERS WITH SIZE (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
	//it visits every element in bordersPerProc (one for every neighbor proc)
	//for every element it visits the border octants it contains and pack its marker in a new structure, sendBuffers
	//this map has an entry Class_Comm_Buffer for every proc containing the size in bytes of the buffer and the octants marker
	//to be sent to that proc packed in a char* buffer
	int8_t marker;
	bool info[16], mod;
	map<int,Class_Comm_Buffer> sendBuffers;
	map<int,vector<uint32_t> >::iterator bitend = bordersPerProc.end();
	uint32_t pbordersOversize = 0;
	for(map<int,vector<uint32_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
		pbordersOversize += bit->second.size();
//		int buffSize = bit->second.size() * (int)ceil((double)(markerBytes) / (double)(CHAR_BIT/8));
		int buffSize = bit->second.size() * (int)ceil((double)(markerBytes + boolBytes) / (double)(CHAR_BIT/8));
		int key = bit->first;
		const vector<uint32_t> & value = bit->second;
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a');
		int pos = 0;
		int nofBorders = value.size();
		for(int i = 0; i < nofBorders; ++i){
			//the use of auxiliary variable can be avoided passing to MPI_Pack the members of octant but octant in that case cannot be const
			const Class_Octant & octant = octree.octants[value[i]];
			marker = octant.getMarker();
			mod	= octant.info[15];
			error_flag = MPI_Pack(&marker,1,MPI_INT8_T,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
			error_flag = MPI_Pack(&mod,1,MPI::BOOL,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
		}
	}

	//COMMUNICATE THE SIZE OF BUFFER TO THE RECEIVERS
	//the size of every borders buffer is communicated to the right process in order to build the receive buffer
	//and stored in the recvBufferSizePerProc structure
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

	//COMMUNICATE THE BUFFERS TO THE RECEIVERS
	//recvBuffers structure is declared and each buffer is initialized to the right size
	//then, sendBuffers are communicated by senders and stored in recvBuffers in the receivers
	//at the same time every process compute the size in bytes of all the level and marker of ghost octants
	uint32_t nofBytesOverProc = 0;
	map<int,Class_Comm_Buffer> recvBuffers;
	map<int,int>::iterator ritend = recvBufferSizePerProc.end();
	for(map<int,int>::iterator rit = recvBufferSizePerProc.begin(); rit != ritend; ++rit){
		recvBuffers[rit->first] = Class_Comm_Buffer(rit->second,'a');
	}
	nReq = 0;
	for(map<int,Class_Comm_Buffer>::iterator sit = sendBuffers.begin(); sit != sitend; ++sit){
		nofBytesOverProc += recvBuffers[sit->first].commBufferSize;
		error_flag = MPI_Irecv(recvBuffers[sit->first].commBuffer,recvBuffers[sit->first].commBufferSize,MPI_PACKED,sit->first,rank,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	for(map<int,Class_Comm_Buffer>::reverse_iterator rsit = sendBuffers.rbegin(); rsit != rsitend; ++rsit){
		error_flag =  MPI_Isend(rsit->second.commBuffer,rsit->second.commBufferSize,MPI_PACKED,rsit->first,rsit->first,MPI_COMM_WORLD,&req[nReq]);
		++nReq;
	}
	MPI_Waitall(nReq,req,stats);

/*
	//COMPUTE GHOSTS SIZE IN BYTES
	//number of ghosts in every process is obtained through the size in bytes of the single octant
	//and ghost vector in local tree is resized
	uint32_t nofGhosts = nofBytesOverProc / (uint32_t)(levelBytes+markerBytes);
	octree.size_ghosts = nofGhosts;
	octree.ghosts.resize(nofGhosts);
*/

	//UNPACK BUFFERS AND BUILD GHOSTS CONTAINER OF CLASS_LOCAL_TREE
	//every entry in recvBuffers is visited, each buffers from neighbor processes is unpacked octant by octant.
	//every ghost octant is built and put in the ghost vector
	uint32_t ghostCounter = 0;
	map<int,Class_Comm_Buffer>::iterator rritend = recvBuffers.end();
	for(map<int,Class_Comm_Buffer>::iterator rrit = recvBuffers.begin(); rrit != rritend; ++rrit){
		int pos = 0;
		int nofGhostsPerProc = int(rrit->second.commBufferSize / ((uint32_t) (markerBytes + boolBytes)));
		for(int i = 0; i < nofGhostsPerProc; ++i){
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&marker,1,MPI_INT8_T,MPI_COMM_WORLD);
			octree.ghosts[ghostCounter].setMarker(marker);
			error_flag = MPI_Unpack(rrit->second.commBuffer,rrit->second.commBufferSize,&pos,&mod,1,MPI::BOOL,MPI_COMM_WORLD);
			octree.ghosts[ghostCounter].info[15] = mod;
			++ghostCounter;
		}
	}
	recvBuffers.clear();
	sendBuffers.clear();
	recvBufferSizePerProc.clear();
}

//==============================================================

void Class_Para_Tree::balance21(){
	bool globalDone = true, localDone = false;
	int  iteration  = 0;


	writeLog("---------------------------------------------");
	writeLog(" 2:1 BALANCE (balancing Marker before Adapt)");
	writeLog(" ");
	writeLog(" Iterative procedure	");
	writeLog(" ");
	writeLog(" Iteration	:	" + to_string(iteration));


	commMarker();
	localDone = octree.localBalance(true);
	MPI_Barrier(MPI_COMM_WORLD);
	error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,MPI_COMM_WORLD);

	while(globalDone){
		iteration++;
		writeLog(" Iteration	:	" + to_string(iteration));
		commMarker();
		localDone = octree.localBalance(false);
		error_flag = MPI_Allreduce(&localDone,&globalDone,1,MPI::BOOL,MPI_LOR,MPI_COMM_WORLD);
	}

	commMarker();
	writeLog(" Iteration	:	Finalizing ");
	writeLog(" ");
	localDone = octree.localBalance(false);
	commMarker();

	writeLog(" 2:1 Balancing reached ");
	writeLog(" ");
	writeLog("---------------------------------------------");

}

//==============================================================
