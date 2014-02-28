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

// =================================================================================== //
// Basic Get/Set methods ============================================================= //

uint32_t Class_Local_Tree::getNumOctants() const {
	return octants.size();
}

uint8_t Class_Local_Tree::getLocalMaxDepth() const {
	return local_max_depth;
}

uint8_t Class_Local_Tree::getMarker(int32_t idx) {
	return octants[idx].getMarker();
}

bool Class_Local_Tree::getBalance(int32_t idx) {
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


void Class_Local_Tree::setMarker(int32_t idx, int8_t marker) {
	octants[idx].setMarker(marker);
}

void Class_Local_Tree::setBalance(int32_t idx, bool balance) {
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

// =================================================================================== //
// Debug methods ===================================================================== //

void Class_Local_Tree::addOctantToTree(Class_Octant octant){
	octants.push_back(octant);
	octants.shrink_to_fit();
}


// =================================================================================== //
// Other methods ===================================================================== //

const Class_Octant& Class_Local_Tree::extractOctant(uint32_t idx) const {
	return octants[idx];
}

// =================================================================================== //

bool Class_Local_Tree::refine() {

	// Local variables
	vector<uint32_t> last_child_index;
	Class_Octant* children;
	uint32_t idx, nocts;
	uint32_t offset = 0, blockidx;
	uint8_t nchm1 = nchildren-1, ich;
	bool dorefine = false;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker() && octants[idx].getLevel() < MAX_LEVEL){
			last_child_index.push_back(idx+nchm1+offset);
			offset += nchm1;
		}
		else{
			octants[idx].info[12] = false;
			octants[idx].marker = 0;
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
				if (children[0].getMarker() > 0){
					dorefine = true;
				}
				delete []children;
			}
			else {
				octants[idx] = octants[idx-offset];
			}
		}
	}
	octants.shrink_to_fit();
	return dorefine;
}

// =================================================================================== //

bool Class_Local_Tree::coarse() {

/*
	{

		//TODO DA TOGLIERE DOPO DEBUG CERTO DEL NUOVO ALGORITMO
		// Local variables
		vector<uint32_t> first_child_index;
		vector<uint32_t> first_child_index_for_ghosts;
		vector<uint8_t>  nbrothers_for_ghosts;
		Class_Octant father;
		uint64_t idx, idx2, ich, nocts, nghosts;
		uint64_t offset = 0, blockidx;
		int8_t markerfather, nbrothers;
		uint8_t nchm1 = nchildren-1, nmarker, iface;
		uint32_t nidx = 0;
		bool docoarse = false;

		nocts   = octants.size();
		nghosts = ghosts.size();


		// Check and coarse in ghost
		// If refined the father goes to the lower processor ...
		for (idx=0; idx<nghosts; idx++){
			if(ghosts[idx].getMarker() < 0 && ghosts[idx].getLevel() > 0){
				nmarker = 0;
				father = ghosts[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+nchildren; idx2++){
					if (idx2<nghosts){
						if(ghosts[idx2].getMarker() < 0 && ghosts[idx2].buildFather() == father){
							nmarker++;
						}
					}
				}
				if (nmarker != nchildren){
					bool first_child = false;
					nbrothers = 0;
					for (idx2=0; idx2<nchildren; idx2++){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nmarker++;
							if (!first_child){
								first_child_index_for_ghosts.push_back(idx2);
								first_child = true;
								nbrothers++;
							}
							if (nmarker == nchildren){
								nbrothers_for_ghosts.push_back(nbrothers);
								nidx += nbrothers + 1;
							}
						}
					}
					first_child = false;
					nbrothers = 0;
					for (idx2=nocts-nchildren; idx2<nocts; idx2++){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nmarker++;
							if (!first_child){
								first_child_index_for_ghosts.push_back(idx2);
								first_child = true;
								nbrothers++;
							}
							if (nmarker == nchildren){
								nbrothers_for_ghosts.push_back(nbrothers);
								nidx += nbrothers;
							}
						}
					}
				}
			}
		}
		if (nidx != 0){
			uint32_t nblock = nocts - nidx;
			nidx = 0;
			for (idx=0; idx<nocts; idx++){
				for (idx=first_child_index_for_ghosts[0]; idx<nblock; idx++){
					if (idx+offset == first_child_index_for_ghosts[nidx]){
						markerfather = -MAX_LEVEL;
						father = octants[idx+offset].buildFather();
						for(idx2=0; idx2<nbrothers_for_ghosts[nidx]+1; idx2++){
							if (markerfather < octants[idx+offset+idx2].getMarker()-1){
								markerfather = octants[idx+offset+idx2].getMarker()-1;
							}
							for (iface=0; iface<nface; iface++){
								father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
								father.info[iface+nface] = (father.info[iface+nface] || octants[idx+offset+idx2].info[iface+nface]);
							}
							father.info[13] = true;
							father.setMarker(markerfather);
							if (markerfather < 0){
								docoarse = true;
							}
							if(idx+offset<nchildren){
								offset += nbrothers_for_ghosts[nidx];
								octants[idx] = octants[idx+offset];
								nidx++;
							}
							else{
								octants[idx] = father;
								offset += nbrothers_for_ghosts[nidx];
								nidx++;
							}
						}
					}
					else{
						octants[idx] = octants[idx+offset];
					}
				}
			}
		}
		octants.resize(nocts-offset);
		octants.shrink_to_fit();

		// Check and coarse internal octants
		offset = 0;
		for (idx=0; idx<nocts; idx++){
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nmarker = 0;
				father = octants[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+nchildren; idx2++){
					if (idx2<nocts){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nmarker++;
						}
					}
				}
				if (nmarker == nchildren){
					nidx++;
					first_child_index.push_back(idx);
					idx = idx2-1;
				}
				else{
					octants[idx].setMarker(0);
				}
			}
			else{
				octants[idx].info[13] = false;
			}
		}
		if (nidx!=0){
			uint32_t nblock = nocts - nidx*nchm1;
			nidx = 0;
			for (idx=first_child_index[0]; idx<nblock; idx++){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -MAX_LEVEL;
					father = octants[idx+offset].buildFather();
					for(idx2=0; idx2<nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()-1){
							markerfather = octants[idx+offset+idx2].getMarker()-1;
						}
						for (iface=0; iface<nface; iface++){
							father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
							father.info[iface+nface] = (father.info[iface+nface] || octants[idx+offset+idx2].info[iface+nface]);
						}
					}
					father.info[13] = true;
					father.setMarker(markerfather);
					if (markerfather < 0){
						docoarse = true;
					}
					octants[idx] = father;
					offset += nchm1;
					nidx++;
				}
				else{
					octants[idx] = octants[idx+offset];
				}
			}
		}
		octants.resize(nocts-offset);
		octants.shrink_to_fit();
		return docoarse;

	}
*/


	// Local variables
	vector<uint32_t> first_child_index;
	Class_Octant father;
	uint32_t idx, idx2, ich, nocts, nghosts;
	uint32_t offset;
	uint32_t idx1_gh, idx2_gh;
	uint32_t nidx;
	int8_t markerfather, marker;
	uint8_t nbro, nstart, nend;
	uint8_t nchm1 = nchildren-1, iface;
	uint8_t oppface[nface];
	bool docoarse = false;

	//------------------------------------------ //
	// Initialization

	nbro = nstart = nend = 0;
	nidx = offset = 0;

	idx1_gh = idx2_gh = 0;
	// Set matrix of opposite faces
	for (iface=0; iface<DIM; iface++){
		oppface[2*iface]   = 2*iface+1;
		oppface[2*iface+1] = 2*iface;
	}

	nocts   = octants.size();
	nghosts = ghosts.size();

	// Init first and last desc (even if already calculated)
	setFirstDesc();
	setLastDesc();

	//------------------------------------------ //

	// Set index for start and end check for ghosts
	if (ghosts.size()){
		while(ghosts[idx1_gh].computeMorton() < first_desc.computeMorton()){
			idx1_gh++;
		}
		idx1_gh--;
		while(ghosts[idx2_gh].computeMorton() < last_desc.computeMorton()){
			idx2_gh++;
		}

		// Start on ghosts
		if ((ghosts[idx1_gh].getMarker() < 0) & (octants[0].getMarker() < 0)){
			father = ghosts[idx1_gh].buildFather();
			nbro = 0;
			idx = idx1_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 & ghosts[idx].buildFather() == father){
				nbro++;
				marker = ghosts[idx].getMarker();
				idx--;
			}
			nstart = 0;
			idx = 0;
			marker = octants[idx].getMarker();
			while(marker<0 & octants[idx].buildFather() == father){
				nbro++;
				marker = octants[idx].getMarker();
				nstart++;
				idx++;
			}
			if (nbro == nchildren){
				offset = nstart;
				// For update pbound of neighbours only check
				// the odd faces of new father (placed nstart-times
				// in the first nstart positions of octants)
				// If there is father after coarse will be the first
				// element of local octants (lowest Morton)
				for (int i=0; i<nstart; i++){
					octants[i] = father;
				}
				uint32_t	 sizeneigh;
				u32vector    neigh;
				vector<bool> isghost;
				for (iface=0; iface<DIM; iface++){
					uint8_t oddface = ((iface*2)+1);
					findNeighbours(nstart-1, oddface, neigh, isghost);
					sizeneigh = neigh.size();
					for(int i=0; i<sizeneigh; i++){
						if (!isghost[i])
							octants[neigh[i]].setPbound(oppface[oddface], true);
					}
				}
			}
			else{
				nstart == 0;
			}
		}
	}


	// Check and coarse internal octants
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
			nbro = 0;
			father = octants[idx].buildFather();
			// Check if family is to be refined
			for (idx2=idx; idx2<idx+nchildren; idx2++){
				if (idx2<nocts){
					if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
						nbro++;
					}
				}
			}
			if (nbro == nchildren){
				nidx++;
				first_child_index.push_back(idx);
				idx = idx2-1;
			}
			else{
				octants[idx].setMarker(0);
			}
		}
		else{
			octants[idx].info[13] = false;
		}
	}
	//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
	if (nidx!=0){
		uint32_t nblock = nocts - nidx*nchm1 - nstart;
		nidx = 0;
		for (idx=0; idx<nblock; idx++){
			if (idx+offset == first_child_index[nidx]){
				markerfather = -MAX_LEVEL;
				father = octants[idx+offset].buildFather();
				for(idx2=0; idx2<nchildren; idx2++){
					if (markerfather < octants[idx+offset+idx2].getMarker()+1){
						markerfather = octants[idx+offset+idx2].getMarker()+1;
					}
					for (iface=0; iface<nface; iface++){
						father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
						father.info[iface+nface] = (father.info[iface+nface] || octants[idx+offset+idx2].info[iface+nface]);
					}
				}
				father.info[13] = true;
				father.setMarker(markerfather);
				if (markerfather < 0){
					docoarse = true;
				}
				octants[idx] = father;
				offset += nchm1;
				nidx++;
			}
			else{
				octants[idx] = octants[idx+offset];
			}
		}
	}
	octants.resize(nocts-offset);
	octants.shrink_to_fit();
	nocts = octants.size();


	// End on ghosts
	if (ghosts.size()){
		if ((ghosts[idx2_gh].getMarker() < 0) & (octants[nocts-1].getMarker() < 0)){
			father = ghosts[idx2_gh].buildFather();
			markerfather = -MAX_LEVEL;
			nbro = 0;
			idx = idx2_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 & ghosts[idx].buildFather() == father){
				nbro++;
				marker = ghosts[idx].getMarker();
				if (markerfather < octants[idx+offset+idx2].getMarker()+1){
					markerfather = octants[idx+offset+idx2].getMarker()+1;
				}
				idx++;
			}
			nend = 0;
			idx = nocts-1;
			marker = octants[idx].getMarker();
			while(marker<0 & octants[idx].buildFather() == father){
				nbro++;
				marker = octants[idx].getMarker();
				if (markerfather < octants[idx+offset+idx2].getMarker()+1){
					markerfather = octants[idx+offset+idx2].getMarker()+1;
				}
				nend++;
				idx--;
			}
			if (nbro == nchildren){
				offset = nend;
			}
			else{
				nend == 0;
			}
		}
		if (nend != 0){
			for (idx=0; idx < nend; idx++){
				for (iface=0; iface<nface; iface++){
					father.info[iface] = (father.info[iface] || octants[nocts-idx].info[iface]);
					father.info[iface+nface] = (father.info[iface+nface] || octants[nocts-idx].info[iface+nface]);
				}
			}
			father.info[13] = true;
			father.setMarker(markerfather);
			if (markerfather < 0){
				docoarse = true;
			}
			octants.resize(nocts-offset);
			octants.push_back(father);
		}
	}
	return docoarse;
}

// =================================================================================== //

void Class_Local_Tree::computeConnectivity() {
	map<uint64_t, vector<uint32_t> > mapnodes;
	map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
	uint32_t i, k, counter;
	uint64_t morton;
	uint32_t noctants = getNumOctants();
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
			for(vector<uint32_t>::iterator iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
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
	u32vector2D().swap(connectivity);
}

void Class_Local_Tree::updateConnectivity() {
	clearConnectivity();
	computeConnectivity();
}

// =================================================================================== //

void Class_Local_Tree::computeghostsConnectivity() {
	map<uint64_t, vector<uint32_t> > mapnodes;
	map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
	uint32_t i, k, counter;
	uint64_t morton;
	uint32_t noctants = size_ghosts;
	uint32_t (*octnodes)[DIM];
	uint8_t j;

	if (ghostsnodes.size() == 0){
		ghostsconnectivity.resize(noctants);
		for (i = 0; i < noctants; i++){
			octnodes = ghosts[i].getNodes();
			for (j = 0; j < nnodes; j++){
#if DIM == 3
				morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
#else
				puppa
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
			for(vector<uint32_t>::iterator iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
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
	u32vector2D().swap(ghostsconnectivity);
}

void Class_Local_Tree::updateghostsConnectivity() {
	clearghostsConnectivity();
	computeghostsConnectivity();
}

// =================================================================================== //

void Class_Local_Tree::updateLocalMaxDepth() {
	uint32_t noctants = getNumOctants();
	uint32_t i;

	local_max_depth = 0;
	for(i = 0; i < noctants; i++){
		if(octants[i].getLevel() > local_max_depth){
			local_max_depth = octants[i].getLevel();
		}
	}
}

// =================================================================================== //

void Class_Local_Tree::findNeighbours(uint32_t idx, uint8_t iface,
		u32vector & neighbours, vector<bool> & isghost) {

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry, idxtry_old, idxtry_old_;
	Class_Octant* oct = &octants[idx];
	uint32_t size = oct->getSize();

	//Alternative to switch case
	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > nface){
		writeLog("Face index out of range in find neighbours !!!");
		return;
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
			int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
			idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
			//idxtry_old = uint64_t((1+direction)*noctants );
			while(abs(jump) > 0){
				Mortontry = octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
				idxtry += jump;
			}
			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
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
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
				Class_Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				vector<uint32_t> bufferidx;
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dh;
				int32_t eqcoord;
				while(Mortontry < Mortonlast){
					Dh = cx*(int32_t(oct->x) - int32_t(octants[idxtry].x));
					Dh += cy*(int32_t(oct->y) - int32_t(octants[idxtry].y));
					Dh += cz*(int32_t(oct->z) - int32_t(octants[idxtry].z));
					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
					}
					idxtry++;
					Mortontry = octants[idxtry].computeMorton();
				}
				return;
			}
		}
		else{
			// Boundary Face
			return;
		}
	}
	//--------------------------------------------------------------- //
	//--------------------------------------------------------------- //
	else{
		// Check if octants face is a boundary
		if (oct->info[iface] == false){
			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS

			// Search in ghosts

			uint32_t idxghost = uint32_t(size_ghosts/2);
			Class_Octant* octghost = &ghosts[idxghost];

			//Build Morton number of virtual neigh of same size
			Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
				idxtry += jump;
			}
			if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
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
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
				Class_Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				vector<uint32_t> bufferidx;
				Mortontry = ghosts[idxtry].computeMorton();
				int32_t Dh;
				int32_t eqcoord;
				while(Mortontry < Mortonlast){
					Dh = cx*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
					Dh += cy*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
					Dh += cz*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
					if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
						neighbours.push_back(idxtry);
						isghost.push_back(true);
					}
					idxtry++;
					Mortontry = ghosts[idxtry].computeMorton();
				}
			}
			uint32_t lengthneigh = 0;
			uint32_t sizeneigh = neighbours.size();
			for (idx=0; idx<sizeneigh; idx++){
				lengthneigh += ghosts[neighbours[idx]].getSize();
			}
			if (lengthneigh < oct->getSize()){
				// Search in octants

				// Check if octants face is a boundary
				if (oct->info[iface] == false){

					//Build Morton number of virtual neigh of same size
					Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
					Morton = samesizeoct.computeMorton();
					// Search morton in octants
					// If a even face morton is lower than morton of oct, if odd higher
					// ---> can i search only before or after idx in octants
					int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
					idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
					//idxtry_old = uint64_t((1+direction)*noctants );
					while(abs(jump) > 0){
						Mortontry = octants[idxtry].computeMorton();
						jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
						idxtry += jump;
					}
					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(false);
						neighbours.push_back(idxtry);
						writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
						return;
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
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
							return;
						}
						// Compute Last discendent of virtual octant of same size
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
						Class_Octant last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						int32_t Dh;
						int32_t eqcoord;
						while(Mortontry < Mortonlast){
							Dh = cx*(int32_t(oct->x) - int32_t(octants[idxtry].x));
							Dh += cy*(int32_t(oct->y) - int32_t(octants[idxtry].y));
							Dh += cz*(int32_t(oct->z) - int32_t(octants[idxtry].z));
							if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
							idxtry++;
							Mortontry = octants[idxtry].computeMorton();
						}
					}
				}
			}
			return;
		}
		else{
			// Boundary Face
			return;
		}

	}
}

// =================================================================================== //

void Class_Local_Tree::findNeighbours(Class_Octant const & oct, uint8_t iface,
		u32vector & neighbours, vector<bool> & isghost) {

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry, idxtry_old, idxtry_old_;
	uint32_t size = oct.getSize();

	//Alternative to switch case
	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > nface){
		writeLog("Face index out of range in find neighbours !!!");
		return;
	}

	// Check if octants face is a process boundary
	if (oct.info[nface+iface] == false){

		// Check if octants face is a boundary
		if (oct.info[iface] == false){

			// Find octant in octants
			uint32_t idx;
			{
			vector<Class_Octant>::iterator it = find(octants.begin(), octants.end(), oct);
			idx = distance(octants.begin(), it);
			it = octants.end();
			}

			//Build Morton number of virtual neigh of same size
			Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// --. can i search only before or after idx in octants
			int32_t jump = (oct.computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
			idxtry = uint32_t(idx +((oct.computeMorton()<Morton)-(oct.computeMorton()>Morton))*jump);
			//idxtry_old = uint64_t((1+direction)*noctants );
			while(abs(jump) > 0){
				Mortontry = octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
				idxtry += jump;
			}
			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
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
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
				Class_Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				vector<uint32_t> bufferidx;
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dh;
				int32_t eqcoord;
				while(Mortontry < Mortonlast){
					Dh = cx*(int32_t(oct.x) - int32_t(octants[idxtry].x));
					Dh += cy*(int32_t(oct.y) - int32_t(octants[idxtry].y));
					Dh += cz*(int32_t(oct.z) - int32_t(octants[idxtry].z));
					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
					}
					idxtry++;
					Mortontry = octants[idxtry].computeMorton();
				}
				return;
			}
		}
		else{
			// Boundary Face
			return;
		}
	}
	//--------------------------------------------------------------- //
	//--------------------------------------------------------------- //
	else{
		// Check if octants face is a boundary
		if (oct.info[iface] == false){

			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS

			// Find octant in octants
			uint32_t idx;
			{
			vector<Class_Octant>::iterator it = find(octants.begin(), octants.end(), oct);
			idx = distance(octants.begin(), it);
			it = octants.end();
			}

			// Search in ghosts

			uint32_t idxghost = uint32_t(size_ghosts/2);
			Class_Octant* octghost = &ghosts[idxghost];

			//Build Morton number of virtual neigh of same size
			Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct.x-size,oct.y,oct.z);
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// --. can i search only before or after idx in octants
			int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
				idxtry += jump;
			}
			if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct.level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
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
				if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct.level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
				Class_Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				vector<uint32_t> bufferidx;
				Mortontry = ghosts[idxtry].computeMorton();
				int32_t Dh;
				int32_t eqcoord;
				while(Mortontry < Mortonlast){
					Dh = cx*(int32_t(oct.x) - int32_t(ghosts[idxtry].x));
					Dh += cy*(int32_t(oct.y) - int32_t(ghosts[idxtry].y));
					Dh += cz*(int32_t(oct.z) - int32_t(ghosts[idxtry].z));
					if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
						neighbours.push_back(idxtry);
						isghost.push_back(true);
					}
					idxtry++;
					Mortontry = ghosts[idxtry].computeMorton();
				}
			}
			uint32_t areaneigh = 0;
			uint32_t sizeneigh = neighbours.size();
			for (idx=0; idx<sizeneigh; idx++){
				areaneigh += ghosts[neighbours[idx]].getArea();
			}
			if (areaneigh < oct.getArea()){
				// Search in octants

				// Check if octants face is a boundary
				if (oct.info[iface] == false){

					//Build Morton number of virtual neigh of same size
					Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
					Morton = samesizeoct.computeMorton();
					// Search morton in octants
					// If a even face morton is lower than morton of oct, if odd higher
					// --. can i search only before or after idx in octants
					int32_t jump = (oct.computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
					idxtry = uint32_t(idx +((oct.computeMorton()<Morton)-(oct.computeMorton()>Morton))*jump);
					//idxtry_old = uint64_t((1+direction)*noctants );
					while(abs(jump) > 0){
						Mortontry = octants[idxtry].computeMorton();
						jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
						idxtry += jump;
					}
					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
						//Found neighbour of same size
						isghost.push_back(false);
						neighbours.push_back(idxtry);
						writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
						return;
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
						if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
							return;
						}
						// Compute Last discendent of virtual octant of same size
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
						Class_Octant last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						int32_t Dh;
						int32_t eqcoord;
						while(Mortontry < Mortonlast){
							Dh = cx*(int32_t(oct.x) - int32_t(octants[idxtry].x));
							Dh += cy*(int32_t(oct.y) - int32_t(octants[idxtry].y));
							Dh += cz*(int32_t(oct.z) - int32_t(octants[idxtry].z));
							if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
							idxtry++;
							Mortontry = octants[idxtry].computeMorton();
						}
					}
				}
			}
			return;
		}
		else{
			// Boundary Face
			return;
		}

	}
}

// =================================================================================== //

void Class_Local_Tree::localBalance(){

	//TODO DA FINIRE!!


	// Local variables
	uint32_t noctants = getNumOctants();
	uint32_t sizeneigh;
	vector<uint32_t> neigh;
	vector<uint32_t> modified, newmodified, pborder;
	uint8_t i, idx, iface;
	vector<bool> isghost;

	// First loop on the octants
	for(idx=0 ; idx<noctants, idx++){
		if (!octants[idx].getNotBalance()){
			for (iface=0; iface<nface; iface++){
				findNeighbours(idx, iface, neigh, isghost);
				for(i=0; i<sizeneigh; i++){
					if (!isghost[i]){
						{
							if(octants[neigh[i]].getMarker() > octants[idx.getMarker() + 1]){
								octants[idx].setMarker(octants[neigh[i]].getMarker()-1);
								modified.push_back(idx);
							}
							else if(octants[neigh[i]].getMarker() < octants[idx.getMarker() - 1]){
								octants[neigh[i]].setMarker(octants[idx].getMarker()-1);
								modified.push_back(neigh[i]);
							}
						};
					}
					else{
						if(ghosts[neigh[i]].getMarker() > octants[idx.getMarker() + 1]){
							octants[idx].setMarker(ghosts[neigh[i]].getMarker()-1);
							modified.push_back(idx);
							pborder.push_back(idx);
						}
						else if(ghosts[neigh[i]].getMarker() < octants[idx.getMarker() - 1]){
							ghosts[neigh[i]].setMarker(octants[idx].getMarker()-1);
						}
					}
				}
			}
		}
	}







}



// =================================================================================== //
