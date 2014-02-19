/*
 * Class_Local_Tree.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#include "Class_Local_Tree.hpp"

Class_Local_Tree::Class_Local_Tree() {
	// TODO Auto-generated constructor stub
	Class_Octant oct0(0,0,0,0);
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
}

//-------------------------------------------------------------------------------- //
// Other methods ----------------------------------------------------------------- //

void Class_Local_Tree::refine() {

	// Local variables
	vector<uint64_t> last_child_index;
	vector<Class_Octant> children;
	uint64_t idx, ich, nocts;
	uint64_t offset = 0, blockidx;
	uint8_t nchm1 = nchildren-1;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker()){
			last_child_index.push_back(idx+nchm1+offset);
			offset += nchm1;
		}
	}
	octants.resize(octants.size()+offset);
	blockidx = last_child_index[0]-nchm1;
	idx = octants.size();
	while (idx>blockidx){
		idx--;
		if(octants[idx-offset].getMarker()){
			octants[idx-offset].buildChildren(children);
			for (ich=0; ich<nchildren; ich++){
				octants[idx-ich]=children[nchm1-ich];
			}
			offset -= nchm1;
			idx -= nchm1;
			//Update local max depth
			if (children[0].getLevel() > local_max_depth){
				local_max_depth = children[0].getLevel();
			}
		}
		else {
			octants[idx] = octants[idx-offset];
		}
	}
}
