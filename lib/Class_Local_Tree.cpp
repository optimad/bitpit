/*
 * Class_Local_Tree.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#include "Class_Local_Tree.hpp"

Class_Local_Tree::Class_Local_Tree() {
	Class_Octant oct0;
	Class_Octant octf(MAX_LEVEL,0,0,0);
	Class_Octant octl(MAX_LEVEL,max_length-1,max_length-1,max_length-1);
	octants.resize(1);
	octants[0] = oct0;
	first_desc = octf;
	last_desc = octl;
	size_ghosts = 0;
	local_max_depth = 0;

}

Class_Local_Tree::~Class_Local_Tree() {
}

//-------------------------------------------------------------------------------- //
// Basic Get/Set methods --------------------------------------------------------- //

uint64_t Class_Local_Tree::getNumOctants() const {
	return octants.size();
}

uint8_t Class_Local_Tree::getLocalMaxDepth() const {
	return local_max_depth;
}

uint8_t Class_Local_Tree::getMarker(int64_t idx) {
	return octants[idx].getMarker();
}

bool Class_Local_Tree::getBalance(int64_t idx) {
	return octants[idx].getNotBalance();
}

const Class_Octant & Class_Local_Tree::getFirstDesc() const {
	return first_desc;
}

const Class_Octant & Class_Local_Tree::getLastDesc() const {
	return last_desc;
}

uint32_t Class_Local_Tree::getSizeGhost() const {
	return size_ghosts;
}


void Class_Local_Tree::setMarker(int64_t idx, int8_t marker) {
	octants[idx].setMarker(marker);
}

void Class_Local_Tree::setBalance(int64_t idx, bool balance) {
	octants[idx].setBalance(balance);
}

void Class_Local_Tree::setFirstDesc() {
	OctantsType::const_iterator firstOctant = octants.begin();
	first_desc = Class_Octant(MAX_LEVEL,firstOctant->x,firstOctant->y,firstOctant->z);
}

void Class_Local_Tree::setLastDesc() {
	OctantsType::const_iterator lastOctant = octants.end() - 1;
	uint32_t x,y,z,delta;
	delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - lastOctant->level)) - 1;
	x = lastOctant->x + delta;
	y = lastOctant->y + delta;
	z = lastOctant->z + delta;
	last_desc = Class_Octant(MAX_LEVEL,x,y,z);
}

//-------------------------------------------------------------------------------- //
// Debug methods ----------------------------------------------------------------- //

void Class_Local_Tree::addOctantToTree(Class_Octant octant){
	octants.push_back(octant);
	octants.shrink_to_fit();
}

const Class_Octant& Class_Local_Tree::extractOctant(uint64_t idx) const {
	return octants[idx];
}

//-------------------------------------------------------------------------------- //
// Other methods ----------------------------------------------------------------- //

void Class_Local_Tree::refine() {

	// Local variables
	vector<uint64_t> last_child_index;
	Class_Octant* children;
	uint64_t idx, ich, nocts;
	uint64_t offset = 0, blockidx;
	uint8_t nchm1 = nchildren-1;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker() && octants[idx].getLevel() < MAX_LEVEL){
			last_child_index.push_back(idx+nchm1+offset);
			offset += nchm1;
		}
		else{
			octants[idx].info[12] = false;
		}
	}
	if (offset >0){
		octants.resize(octants.size()+offset);
		blockidx = last_child_index[0]-nchm1;
		idx = octants.size();
		while (idx>blockidx){
			idx--;
			if(octants[idx-offset].getMarker() && octants[idx-offset].getLevel() < MAX_LEVEL){
				children = octants[idx-offset].buildChildren();
				for (ich=0; ich<nchildren; ich++){
					octants[idx-ich] = (children[nchm1-ich]);
				}
				offset -= nchm1;
				idx -= nchm1;
				//Update local max depth
				if (children[0].getLevel() > local_max_depth){
					local_max_depth = children[0].getLevel();
				}
				delete []children;
			}
			else {
				octants[idx] = octants[idx-offset];
			}
		}
	}
	octants.shrink_to_fit();
}

//-------------------------------------------------------------------------------- //

void Class_Local_Tree::computeConnectivity() {
	map<uint64_t, vector<uint64_t> > mapnodes;
	map<uint64_t, vector<uint64_t> >::iterator iter, iterend;
	uint64_t i, k, morton, counter;
	uint64_t noctants = getNumOctants();
	uint32_t (*octnodes)[DIM];
	uint8_t j;

	if (nodes.size() == 0){
		connectivity.resize(noctants);
		for (i = 0; i < noctants; i++){
			octnodes = octants[i].getNodes();
			for (j = 0; j < nnodes; j++){
#if DIM == 3
				morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
#else
#endif
				if (mapnodes[morton].size()==0){
					for (k = 0; k < DIM; k++){
						mapnodes[morton].push_back(octnodes[j][k]);
					}
				}
				mapnodes[morton].push_back(i);
			}
			delete []octnodes;
		}
		iter	= mapnodes.begin();
		iterend	= mapnodes.end();
		counter = 0;
		while (iter != iterend){
			vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+DIM);
			nodes.push_back(nodecasting);
			for(vector<uint64_t>::iterator iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
				connectivity[(*iter2)].push_back(counter);
			}
			mapnodes.erase(iter++);
			counter++;
		}
		nodes.shrink_to_fit();
		connectivity.shrink_to_fit();
	}
	iter = mapnodes.end();
}

void Class_Local_Tree::clearConnectivity() {
	u32vector2D().swap(nodes);
	u64vector2D().swap(connectivity);
}

//-------------------------------------------------------------------------------- //

void Class_Local_Tree::computeghostsConnectivity() {
	map<uint64_t, vector<uint64_t> > mapnodes;
	map<uint64_t, vector<uint64_t> >::iterator iter, iterend;
	uint64_t i, k, morton, counter;
	uint64_t noctants = size_ghosts;
	uint32_t (*octnodes)[DIM];
	uint8_t j;

	if (nodes.size() == 0){
		connectivity.resize(noctants);
		for (i = 0; i < noctants; i++){
			octnodes = ghosts[i].getNodes();
			for (j = 0; j < nnodes; j++){
#if DIM == 3
				morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
#else
#endif
				if (mapnodes[morton].size()==0){
					for (k = 0; k < DIM; k++){
						mapnodes[morton].push_back(octnodes[j][k]);
					}
				}
				mapnodes[morton].push_back(i);
			}
			delete []octnodes;
		}
		iter	= mapnodes.begin();
		iterend	= mapnodes.end();
		counter = 0;
		while (iter != iterend){
			vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+DIM);
			ghostsnodes.push_back(nodecasting);
			for(vector<uint64_t>::iterator iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
				ghostsconnectivity[(*iter2)].push_back(counter);
			}
			mapnodes.erase(iter++);
			counter++;
		}
		ghostsnodes.shrink_to_fit();
		ghostsconnectivity.shrink_to_fit();
	}
	iter = mapnodes.end();
}

void Class_Local_Tree::clearghostsConnectivity() {
	u32vector2D().swap(ghostsnodes);
	u64vector2D().swap(ghostsconnectivity);
}


//-------------------------------------------------------------------------------- //

void Class_Local_Tree::updateLocalMaxDepth() {
	uint64_t noctants = getNumOctants();
	uint64_t i;

	local_max_depth = 0;
	for(i = 0; i < noctants; i++){
		if(octants[i].getLevel() > local_max_depth){
			local_max_depth = octants[i].getLevel();
		}
	}
}

uint64_t* Class_Local_Tree::findNeighbours(uint64_t idx, uint8_t iface,
		uint8_t& sizeneigh, bool isghost) {

//	uint64_t* NeighIdx = new uint64_t[];
	uint64_t  noctants = getNumOctants();
	uint64_t  Morton, Mortontry, idxtry, idxtry_old, idxtry_old_;
	Class_Octant* oct = octants[idx];
	uint32_t size = oct->getSize();

	// Check if octants face is a process boundary
	if (oct->info[nface+iface] == false){

		switch (iface) {
		case 0 :
		{
			//Build Morton number of virtual neigh of same size
			Morton = mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
			//TODO make a method of this...
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int8_t direction = -1;
			int8_t diff = 100;
			idxtry = uint64_t((idx + (1+direction)*noctants ) / 2);
			idxtry_old = uint64_t((1+direction)*noctants );
			while(diff > 1){
				Mortontry = octants[idxtry].computeMorton();
				if (Mortontry > Morton){
					idxtry_old_ = idxtry;
					idxtry += (idxtry + idxtry_old)/2;
					idxtry_old = idxtry_old_;
					diff = abs(idxtry - idxtry_old);
				}
			}
			if(octants[idxtry].computeMorton() == Morton)




		}
		break;
		}


	}
	else{
		// If octants face is a process boundary search in ghosts

	}


}
//-------------------------------------------------------------------------------- //
