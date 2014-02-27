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
	partition_last_desc = new uint64_t[nproc];
	partition_range_globalidx = new uint64_t[nproc];
	uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
	for(int p = 0; p < nproc; ++p){
		partition_range_globalidx[p] = 0;
		partition_last_desc[p] = lastDescMorton;
	}

}

Class_Para_Tree::~Class_Para_Tree() {
	//delete [] partition_last_desc;
	//delete [] partition_range_globalidx;
}


void Class_Para_Tree::computePartition(uint64_t* partition) {
	int division_result = 0;
	int remind = 0;
	division_result = global_num_octants/nproc;
	remind = global_num_octants%nproc;
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
	serial = false;
}

void Class_Para_Tree::loadBalance(){
	uint64_t* partition = new uint64_t [nproc];
	computePartition(partition);
	uint64_t stride = 0;
	if(serial)
	{
		for(int i = 0; i < rank; ++i)
			stride += partition[i];
		Class_Local_Tree::OctantsType::const_iterator first = octree.octants.begin() + stride;
		Class_Local_Tree::OctantsType::const_iterator last = first + partition[rank];
		octree.octants.assign(first, last);
		octree.octants.shrink_to_fit();
		first = octree.octants.end();
		last = octree.octants.end();
	}
	else
	{
	}
	delete [] partition;
}

void Class_Para_Tree::updateRefine() {
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
		int rbuff [nproc];
		error_flag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,&rbuff,1,MPI_UINT64_T,MPI_COMM_WORLD);
		for(int p = 0; p < nproc; ++p){
			partition_range_globalidx[p] = 0;
			for(int pp = 0; pp <=p; ++pp)
				partition_range_globalidx[p] += rbuff[pp];
			--partition_range_globalidx[p];
		}
		//update partition_range_position
		uint64_t lastDescMorton = octree.getLastDesc().computeMorton();
		error_flag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,&partition_last_desc,1,MPI_UINT64_T,MPI_COMM_WORLD);
	}
}

void Class_Para_Tree::adapt() {
	octree.refine();
	updateRefine();
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
	Class_Local_Tree::OctantsType::iterator end = octree.octants.end();
	Class_Local_Tree::OctantsType::iterator begin = octree.octants.begin();
	map<int,vector<uint64_t> > bordersPerProc;
	for(Class_Local_Tree::OctantsType::iterator it = begin; it != end; ++it){
		set<int> procs;
		for(uint8_t i = 0; i < nface; ++i){
			if(it->getBound(i) == false){
				uint8_t virtualNeighborsSize = 0;
				uint64_t* virtualNeighbors = it->computeVirtualMorton(i,max_depth,virtualNeighborsSize);
				uint8_t maxDelta = virtualNeighborsSize/2;
				for(int j = 0; j <= maxDelta; ++j){
					int pBegin = findOwner(virtualNeighbors[j]);
					int pEnd = findOwner(virtualNeighbors[virtualNeighborsSize - 1 - j]);
					procs.insert(pBegin);
					procs.insert(pEnd);
					if(pBegin != rank || pEnd != rank)
						it->setPbound(i,true);
					if(pBegin == pEnd || pBegin == pEnd - 1)
						break;
				}
			}
		}
		set<int>::iterator pitend = procs.end();
		for(set<int>::iterator pit = procs.begin(); pit != pitend; ++pit){
			int p = *pit;
			if(p != rank){
				bordersPerProc[p].push_back(distance(begin,it));
				vector<uint64_t> & bordersSingleProc = bordersPerProc[p];
				if(bordersSingleProc.capacity() - bordersSingleProc.size() < 2)
					bordersSingleProc.reserve(2*bordersSingleProc.size());
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
	map<int,vector<uint64_t> >::iterator bitend = bordersPerProc.end();
	for(map<int,vector<uint64_t> >::iterator bit = bordersPerProc.begin(); bit != bitend; ++bit){
		int buffSize = bit->second.size() * (int)ceil((double)octantBytes / (double)(CHAR_BIT/8));// + (int)ceil((double)sizeof(int)/(double)(CHAR_BIT/8));
		int key = bit->first;
		const vector<uint64_t> & value = bit->second;
		sendBuffers[key] = Class_Comm_Buffer(buffSize,'a');
		int pos = 0;
		int nofBorders = value.size();
		//MPI_Pack(&nofBorders,1,MPI_INT,sendBuffers[key].commBuffer,buffSize,&pos,MPI_COMM_WORLD);
		for(int i = 0; i < nofBorders; ++i){
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

	cout << "Communicate sizes" << endl;

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

	cout << "Communicate buffers" << endl;

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
	cout << "rank: " << rank << " nofGhosts: " << nofGhosts << endl;
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
			++ghostCounter;
		}
	}


	cout << "size Ghosts " << octree.size_ghosts <<  endl;

	recvBuffers.clear();
	sendBuffers.clear();
	recvBufferSizePerProc.clear();

	cout << "Exiting" << endl;
}
