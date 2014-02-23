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

// =================================================================================== //

uint64_t* Class_Local_Tree::findNeighbours(uint64_t idx, uint8_t iface,
		uint8_t& sizeneigh, bool isghost) {

	uint64_t  noctants = getNumOctants();
	uint64_t  Morton, Mortontry, idxtry, idxtry_old, idxtry_old_;
	Class_Octant* oct = &octants[idx];
	uint32_t size = oct->getSize();

	//Alternative to switch case
	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > nface){
		writeLog("Face index out of range in find neighbours !!!");
		isghost = false;
		sizeneigh = 0;
		uint64_t* NeighIdx = new uint64_t[sizeneigh];
		return NeighIdx;
	}

	// Check if octants face is a process boundary
	if (oct->info[nface+iface] == false){

		// Check if octants face is a boundary
		if (oct->info[iface] == false){

			//Build Morton number of virtual neigh of same size
			Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int64_t jump = (oct->computeMorton() > Morton) ? int64_t(idx/2+1) : int64_t((noctants -idx)/2+1);
			idxtry = uint64_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
			//idxtry_old = uint64_t((1+direction)*noctants );
			while(abs(jump) > 0){
				Mortontry = octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
				idxtry += jump;
			}
			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
				//Found neighbour of same size
				sizeneigh = 1;
				uint64_t* NeighIdx = new uint64_t[1];
				NeighIdx[0] = idxtry;
				return NeighIdx;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(octants[idxtry].computeMorton() < Morton){
						idxtry++;
					}
					while(octants[idxtry].computeMorton() > Morton){
						idxtry--;
					}
				}
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					sizeneigh = 1;
					uint64_t* NeighIdx = new uint64_t[1];
					NeighIdx[0] = idxtry;
					return NeighIdx;
				}
				// Compute Last discendent of virtual octant of same size
				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
				Class_Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				vector<uint64_t> bufferidx;
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dh;
				int32_t eqcoord;
				while(Mortontry < Mortonlast){
					Dh = cx*(int32_t(oct->x) - int32_t(octants[idxtry].x));
					Dh += cy*(int32_t(oct->y) - int32_t(octants[idxtry].y));
					Dh += cz*(int32_t(oct->z) - int32_t(octants[idxtry].z));
					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
						bufferidx.push_back(idxtry);
					}
					idxtry++;
					Mortontry = octants[idxtry].computeMorton();
				}
				sizeneigh = bufferidx.size();
				uint64_t* NeighIdx = new uint64_t[sizeneigh];
				for (int i = 0; i < sizeneigh; i++){
					NeighIdx[i] = bufferidx[i];
				}
				return NeighIdx;
			}
		}
		else{
			// Boundary Face
			sizeneigh = 0;
			isghost = false;
			uint64_t* NeighIdx = new uint64_t[sizeneigh];
			return NeighIdx;
		}
	}
	else{
		// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH IN GHOSTS
		uint64_t idxghost = uint64_t(size_ghosts/2);
		Class_Octant* octghost = &ghosts[idxghost];

		//Build Morton number of virtual neigh of same size
		Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int64_t jump = (octghost->computeMorton() > Morton) ? int64_t(idxghost/2+1) : int64_t((size_ghosts -idxghost)/2+1);
		idxtry = uint64_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
		while(abs(jump) > 0){
			Mortontry = ghosts[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
			idxtry += jump;
		}
		if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
			//Found neighbour of same size
			sizeneigh = 1;
			uint64_t* NeighIdx = new uint64_t[1];
			NeighIdx[0] = idxtry;
			return NeighIdx;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(ghosts[idxtry].computeMorton() < Morton){
					idxtry++;
				}
				while(ghosts[idxtry].computeMorton() > Morton){
					idxtry--;
				}
			}
			if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
				//Found neighbour of same size
				sizeneigh = 1;
				uint64_t* NeighIdx = new uint64_t[1];
				NeighIdx[0] = idxtry;
				return NeighIdx;
			}
			// Compute Last discendent of virtual octant of same size
			uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
			Class_Octant last_desc = samesizeoct.buildLastDesc();
			uint64_t Mortonlast = last_desc.computeMorton();
			vector<uint64_t> bufferidx;
			Mortontry = ghosts[idxtry].computeMorton();
			int32_t Dh;
			int32_t eqcoord;
			while(Mortontry < Mortonlast){
				Dh = cx*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
				Dh += cy*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
				Dh += cz*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
				if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
					bufferidx.push_back(idxtry);
				}
				idxtry++;
				Mortontry = ghosts[idxtry].computeMorton();
			}
			sizeneigh = bufferidx.size();
			uint64_t* NeighIdx = new uint64_t[sizeneigh];
			for (int i = 0; i < sizeneigh; i++){
				NeighIdx[i] = bufferidx[i];
			}
			return NeighIdx;
		}
	}
}
//-------------------------------------------------------------------------------- //
