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
			size_t buffSize = 0;
			for(size_t i = 0; i < nofPbordersPerProc; ++i){
				buffSize += userData.size(pborders[i]);
			}
		}
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a');

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
	map<int,Class_Comm_Buffer>::iterator rbitend = recvBuffers.end();
	map<int,Class_Comm_Buffer>::iterator rbitbegin = recvBuffers.begin();
	for(map<int,Class_Comm_Buffer>::iterator rbit = rbitbegin; rbit != rbitend; ++rbit){
		int ghostSize = octree.ghosts.size();
		for(int k = 0; k < ghostSize; ++k){
			userData.scatter(rbit->second, k);
		}
	}


};
