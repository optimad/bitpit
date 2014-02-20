/*
 * Class_Local_Tree.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#include "Class_Local_Tree.hpp"

Class_Local_Tree::Class_Local_Tree() {
	// TODO Auto-generated constructor stub
	Class_Octant oct0;
	Class_Octant octf(0,0,0,MAX_LEVEL);
	Class_Octant octl(max_length-1,max_length-1,max_length-1,MAX_LEVEL);
	octants.resize(1);
	octants[0] = oct0;
	first_desc = octf;
	last_desc = octl;
	size_ghosts = 0;
	local_max_depth = 0;

}

Class_Local_Tree::~Class_Local_Tree() {
	// TODO Auto-generated destructor stub
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
	return octants[idx].getBalance();
}

void Class_Local_Tree::setMarker(int64_t idx, int8_t marker) {
	octants[idx].setMarker(marker);
}

void Class_Local_Tree::setBalance(int64_t idx, bool balance) {
	octants[idx].setBalance(balance);
}

//-------------------------------------------------------------------------------- //
// Debug methods ----------------------------------------------------------------- //

void Class_Local_Tree::addOctantToTree(Class_Octant octant){
	octants.push_back(octant);
	octants.shrink_to_fit();
}

Class_Octant& Class_Local_Tree::extractOctant(uint64_t idx) {
	return octants[idx];
}

Class_Octant Class_Local_Tree::getFirstDesc() const {
	OctantsType::const_iterator firstOctant = octants.begin();
	return Class_Octant(MAX_LEVEL,firstOctant->x,firstOctant->y,firstOctant->z);
}

Class_Octant Class_Local_Tree::getLastDesc() const {
	OctantsType::const_iterator lastOctant = octants.end() - 1;
	uint32_t x,y,z,delta;
	delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - lastOctant->level)) - 1;
	x = lastOctant->x + delta;
	y = lastOctant->y + delta;
	z = lastOctant->z + delta;
	return Class_Octant(MAX_LEVEL,x,y,z);
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

//-------------------------------------------------------------------------------- //
