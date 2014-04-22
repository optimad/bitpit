/*
 * Class_Local_Tree.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#include "Class_Local_Tree.hpp"

// =================================================================================== //
// EXPLICIT INSTANTIATION                                                              //
// =================================================================================== //

template class Class_Local_Tree<2>;
template class Class_Local_Tree<3>;


//Class_Local_Tree::Class_Local_Tree() {
//	Class_Octant oct0;
//	Class_Octant octf(MAX_LEVEL,0,0,0);
//	Class_Octant octl(MAX_LEVEL,max_length-1,max_length-1,max_length-1);
//	octants.resize(1);
//	octants[0] = oct0;
//	first_desc = octf;
//	last_desc = octl;
//	size_ghosts = 0;
//	local_max_depth = 0;
//
//}
//
//Class_Local_Tree::~Class_Local_Tree() {
//}
//
//// =================================================================================== //
//// Basic Get/Set methods ============================================================= //
//
//uint32_t Class_Local_Tree::getNumOctants() const {
//	return octants.size();
//}
//
//uint8_t Class_Local_Tree::getLocalMaxDepth() const {
//	return local_max_depth;
//}
//
//uint8_t Class_Local_Tree::getMarker(int32_t idx) {
//	return octants[idx].getMarker();
//}
//
//bool Class_Local_Tree::getBalance(int32_t idx) {
//	return octants[idx].getNotBalance();
//}
//
//const Class_Octant & Class_Local_Tree::getFirstDesc() const {
//	return first_desc;
//}
//
//const Class_Octant & Class_Local_Tree::getLastDesc() const {
//	return last_desc;
//}
//
//uint32_t Class_Local_Tree::getSizeGhost() const {
//	return size_ghosts;
//}
//
//
//void Class_Local_Tree::setMarker(int32_t idx, int8_t marker) {
//	octants[idx].setMarker(marker);
//}
//
//void Class_Local_Tree::setBalance(int32_t idx, bool balance) {
//	octants[idx].setBalance(balance);
//}
//
//void Class_Local_Tree::setFirstDesc() {
//	OctantsType::const_iterator firstOctant = octants.begin();
//	first_desc = Class_Octant(MAX_LEVEL,firstOctant->x,firstOctant->y,firstOctant->z);
//}
//
//void Class_Local_Tree::setLastDesc() {
//	OctantsType::const_iterator lastOctant = octants.end() - 1;
//	uint32_t x,y,z,delta;
//	delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - lastOctant->level)) - 1;
//	x = lastOctant->x + delta;
//	y = lastOctant->y + delta;
//	z = lastOctant->z + delta;
//	last_desc = Class_Octant(MAX_LEVEL,x,y,z);
//}
//
//// =================================================================================== //
//// Debug methods ===================================================================== //
//
//void Class_Local_Tree::addOctantToTree(Class_Octant octant){
//	octants.push_back(octant);
//	octants.shrink_to_fit();
//}
//
//
//// =================================================================================== //
//// Other methods ===================================================================== //
//
//const Class_Octant& Class_Local_Tree::extractOctant(uint32_t idx) const {
//	return octants[idx];
//}
//
//// =================================================================================== //
//
//bool Class_Local_Tree::refine() {
//
//	// Local variables
//	vector<uint32_t> last_child_index;
//	Class_Octant* children;
//	uint32_t idx, nocts;
//	uint32_t offset = 0, blockidx;
//	uint8_t nchm1 = nchildren-1, ich, iface;
//	bool dorefine = false;
//
//	nocts = octants.size();
//	for (idx=0; idx<nocts; idx++){
//		if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL){
//			last_child_index.push_back(idx+nchm1+offset);
//			offset += nchm1;
//		}
//		else{
////			octants[idx].info[12] = false;
//			if (octants[idx].marker > 0)
//				octants[idx].marker = 0;
//		}
//	}
//	if (offset > 0){
//		octants.resize(octants.size()+offset);
//		blockidx = last_child_index[0]-nchm1;
//		idx = octants.size();
//		while (idx>blockidx){
//			idx--;
//			//TODO Sostituire questo if con il controllo su last_index_child
//			if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL){
//				children = octants[idx-offset].buildChildren();
//				for (ich=0; ich<nchildren; ich++){
//					octants[idx-ich] = (children[nchm1-ich]);
//				}
//				offset -= nchm1;
//				idx -= nchm1;
//				//Update local max depth
//				if (children[0].getLevel() > local_max_depth){
//					local_max_depth = children[0].getLevel();
//				}
//				if (children[0].getMarker() > 0){
//					//More Refinement to do
//					dorefine = true;
//				}
//				delete []children;
//			}
//			else {
//				octants[idx] = octants[idx-offset];
//			}
//		}
//	}
//	octants.shrink_to_fit();
//
//	//Update pborders (adesso inefficiente, loop di nuovo su tutti gli elementi)
//	//Si può trovare la maniera di inserirlo nel loop precedente
//	pborders.clear();
//	nocts = octants.size();
//	pborders.reserve(int(pow(double(nocts),2.0/3.0)*double(nfaces)));
//	for(idx=0; idx<nocts; idx++){
//		for(iface=0; iface<nfaces; iface++){
//			if (octants[idx].info[iface+nfaces]){
//				pborders.push_back(idx);
//				break;
//			}
//		}
//	}
//	pborders.shrink_to_fit();
//
//	setFirstDesc();
//	setLastDesc();
//
//	return dorefine;
//}
//
//// =================================================================================== //
//
//bool Class_Local_Tree::coarse() {
//
///*
//	{
//
//		//TODO DA TOGLIERE DOPO DEBUG CERTO DEL NUOVO ALGORITMO
//		// Local variables
//		vector<uint32_t> first_child_index;
//		vector<uint32_t> first_child_index_for_ghosts;
//		vector<uint8_t>  nbrothers_for_ghosts;
//		Class_Octant father;
//		uint64_t idx, idx2, ich, nocts, nghosts;
//		uint64_t offset = 0, blockidx;
//		int8_t markerfather, nbrothers;
//		uint8_t nchm1 = nchildren-1, nmarker, iface;
//		uint32_t nidx = 0;
//		bool docoarse = false;
//
//		nocts   = octants.size();
//		nghosts = ghosts.size();
//
//
//		// Check and coarse in ghost
//		// If refined the father goes to the lower processor ...
//		for (idx=0; idx<nghosts; idx++){
//			if(ghosts[idx].getMarker() < 0 && ghosts[idx].getLevel() > 0){
//				nmarker = 0;
//				father = ghosts[idx].buildFather();
//				// Check if family is to be refined
//				for (idx2=idx; idx2<idx+nchildren; idx2++){
//					if (idx2<nghosts){
//						if(ghosts[idx2].getMarker() < 0 && ghosts[idx2].buildFather() == father){
//							nmarker++;
//						}
//					}
//				}
//				if (nmarker != nchildren){
//					bool first_child = false;
//					nbrothers = 0;
//					for (idx2=0; idx2<nchildren; idx2++){
//						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
//							nmarker++;
//							if (!first_child){
//								first_child_index_for_ghosts.push_back(idx2);
//								first_child = true;
//								nbrothers++;
//							}
//							if (nmarker == nchildren){
//								nbrothers_for_ghosts.push_back(nbrothers);
//								nidx += nbrothers + 1;
//							}
//						}
//					}
//					first_child = false;
//					nbrothers = 0;
//					for (idx2=nocts-nchildren; idx2<nocts; idx2++){
//						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
//							nmarker++;
//							if (!first_child){
//								first_child_index_for_ghosts.push_back(idx2);
//								first_child = true;
//								nbrothers++;
//							}
//							if (nmarker == nchildren){
//								nbrothers_for_ghosts.push_back(nbrothers);
//								nidx += nbrothers;
//							}
//						}
//					}
//				}
//			}
//		}
//		if (nidx != 0){
//			uint32_t nblock = nocts - nidx;
//			nidx = 0;
//			for (idx=0; idx<nocts; idx++){
//				for (idx=first_child_index_for_ghosts[0]; idx<nblock; idx++){
//					if (idx+offset == first_child_index_for_ghosts[nidx]){
//						markerfather = -MAX_LEVEL;
//						father = octants[idx+offset].buildFather();
//						for(idx2=0; idx2<nbrothers_for_ghosts[nidx]+1; idx2++){
//							if (markerfather < octants[idx+offset+idx2].getMarker()-1){
//								markerfather = octants[idx+offset+idx2].getMarker()-1;
//							}
//							for (iface=0; iface<nface; iface++){
//								father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
//								father.info[iface+nface] = (father.info[iface+nface] || octants[idx+offset+idx2].info[iface+nface]);
//							}
//							father.info[13] = true;
//							father.setMarker(markerfather);
//							if (markerfather < 0){
//								docoarse = true;
//							}
//							if(idx+offset<nchildren){
//								offset += nbrothers_for_ghosts[nidx];
//								octants[idx] = octants[idx+offset];
//								nidx++;
//							}
//							else{
//								octants[idx] = father;
//								offset += nbrothers_for_ghosts[nidx];
//								nidx++;
//							}
//						}
//					}
//					else{
//						octants[idx] = octants[idx+offset];
//					}
//				}
//			}
//		}
//		octants.resize(nocts-offset);
//		octants.shrink_to_fit();
//
//		// Check and coarse internal octants
//		offset = 0;
//		for (idx=0; idx<nocts; idx++){
//			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
//				nmarker = 0;
//				father = octants[idx].buildFather();
//				// Check if family is to be refined
//				for (idx2=idx; idx2<idx+nchildren; idx2++){
//					if (idx2<nocts){
//						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
//							nmarker++;
//						}
//					}
//				}
//				if (nmarker == nchildren){
//					nidx++;
//					first_child_index.push_back(idx);
//					idx = idx2-1;
//				}
//				else{
//					octants[idx].setMarker(0);
//				}
//			}
//			else{
//				octants[idx].info[13] = false;
//			}
//		}
//		if (nidx!=0){
//			uint32_t nblock = nocts - nidx*nchm1;
//			nidx = 0;
//			for (idx=first_child_index[0]; idx<nblock; idx++){
//				if (idx+offset == first_child_index[nidx]){
//					markerfather = -MAX_LEVEL;
//					father = octants[idx+offset].buildFather();
//					for(idx2=0; idx2<nchildren; idx2++){
//						if (markerfather < octants[idx+offset+idx2].getMarker()-1){
//							markerfather = octants[idx+offset+idx2].getMarker()-1;
//						}
//						for (iface=0; iface<nface; iface++){
//							father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
//							father.info[iface+nface] = (father.info[iface+nface] || octants[idx+offset+idx2].info[iface+nface]);
//						}
//					}
//					father.info[13] = true;
//					father.setMarker(markerfather);
//					if (markerfather < 0){
//						docoarse = true;
//					}
//					octants[idx] = father;
//					offset += nchm1;
//					nidx++;
//				}
//				else{
//					octants[idx] = octants[idx+offset];
//				}
//			}
//		}
//		octants.resize(nocts-offset);
//		octants.shrink_to_fit();
//		return docoarse;
//
//	}
//*/
//
//
//	// Local variables
//	vector<uint32_t> first_child_index;
//	Class_Octant father;
//	uint32_t ich, nocts, nghosts;
//	int32_t idx, idx2;
//	uint32_t offset;
//	int32_t idx1_gh, idx2_gh;
//	uint32_t nidx;
//	int8_t markerfather, marker;
//	uint8_t nbro, nstart, nend;
//	uint8_t nchm1 = nchildren-1, iface;
//	bool docoarse = false;
//
//	//------------------------------------------ //
//	// Initialization
//
//	nbro = nstart = nend = 0;
//	nidx = offset = 0;
//
//	idx1_gh = idx2_gh = 0;
//
//	nocts   = octants.size();
//	size_ghosts = ghosts.size();
//
//
//	// Init first and last desc (even if already calculated)
//	setFirstDesc();
//	setLastDesc();
//
//	//------------------------------------------ //
//
//	// Set index for start and end check for ghosts
//	if (ghosts.size()){
//		while(idx1_gh < size_ghosts && ghosts[idx1_gh].computeMorton() < first_desc.computeMorton()){
//			idx1_gh++;
//		}
//		idx1_gh = max(0, idx1_gh-1);
//		while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() < last_desc.computeMorton()){
//			idx2_gh++;
//		}
//		idx2_gh = min(int(size_ghosts-1), idx2_gh);
//
//		// Start on ghosts
//		if ((ghosts[idx1_gh].getMarker() < 0) & (octants[0].getMarker() < 0)){
//			father = ghosts[idx1_gh].buildFather();
//			nbro = 0;
//			idx = idx1_gh;
//			marker = ghosts[idx].getMarker();
//			while(marker < 0 & ghosts[idx].buildFather() == father){
//				nbro++;
//				marker = ghosts[idx].getMarker();
//				idx--;
//				if (idx<0){
//					break;
//				}
//			}
//			nstart = 0;
//			idx = 0;
//			marker = octants[idx].getMarker();
//			while(marker<0 & octants[idx].buildFather() == father){
//				nbro++;
//				marker = octants[idx].getMarker();
//				nstart++;
//				idx++;
//				if (idx==nocts){
//					break;
//				}
//			}
//			if (nbro == nchildren){
////				offset = nstart;
//				// For update pbound of neighbours only check
//				// the odd faces of new father (placed nstart-times
//				// in the first nstart positions of octants)
//				// If there is father after coarse will be the first
//				// element of local octants (lowest Morton)
///*
//				for (int i=0; i<nstart; i++){
//					octants[i] = father;
//				}
//*/
//				uint32_t	 sizeneigh;
//				u32vector    neigh;
//				vector<bool> isghost;
//				for (iface=0; iface<DIM; iface++){
//					uint8_t oddface = ((iface*2)+1);
//					findNeighbours(nstart-1, oddface, neigh, isghost);
//					sizeneigh = neigh.size();
//					for(int i=0; i<sizeneigh; i++){
//						if (!isghost[i])
//							octants[neigh[i]].setPbound(oppface[oddface], true);
//					}
//				}
//			}
//			else{
//				nstart = 0;
//			}
//		}
//	}
//
//	// Check and coarse internal octants
//	for (idx=0; idx<nocts; idx++){
//		if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
//			nbro = 0;
//			father = octants[idx].buildFather();
//			// Check if family is to be refined
//			for (idx2=idx; idx2<idx+nchildren; idx2++){
//				if (idx2<nocts-1){
//					if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
//						nbro++;
//					}
//				}
//			}
//			if (nbro == nchildren){
//				nidx++;
//				first_child_index.push_back(idx);
//				idx = idx2-1;
//			}
//			else{
//				if (idx < (nocts>nchildren)*(nocts-nchildren)){
//					octants[idx].setMarker(0);
//				}
//			}
//		}
//		else{
////			octants[idx].info[13] = false;
//		}
//	}
//	//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
//	if (nidx!=0){
//		uint32_t nblock = nocts - nidx*nchm1 - nstart;
//		nidx = 0;
//		for (idx=0; idx<nblock; idx++){
//			if (idx+offset == first_child_index[nidx]){
//				markerfather = -MAX_LEVEL;
//				father = octants[idx+offset].buildFather();
//				for(idx2=0; idx2<nchildren; idx2++){
//					if (markerfather < octants[idx+offset+idx2].getMarker()+1){
//						markerfather = octants[idx+offset+idx2].getMarker()+1;
//					}
//					for (iface=0; iface<nfaces; iface++){
//						father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
//						father.info[iface+nfaces] = (father.info[iface+nfaces] || octants[idx+offset+idx2].info[iface+nfaces]);
//					}
//				}
//				father.info[13] = true;
//				father.setMarker(markerfather);
//				if (markerfather < 0){
//					docoarse = true;
//				}
//				octants[idx] = father;
//				offset += nchm1;
//				nidx++;
//			}
//			else{
//				octants[idx] = octants[idx+offset];
//			}
//		}
//	}
//	octants.resize(nocts-offset);
//	octants.shrink_to_fit();
//	nocts = octants.size();
//
//
//	// End on ghosts
//	if (ghosts.size() && nocts > 0){
//		if ((ghosts[idx2_gh].getMarker() < 0) & (octants[nocts-1].getMarker() < 0)){
//			father = ghosts[idx2_gh].buildFather();
//			markerfather = -MAX_LEVEL;
//			nbro = 0;
//			idx = idx2_gh;
//			marker = ghosts[idx].getMarker();
//			while(marker < 0 & ghosts[idx].buildFather() == father){
//				nbro++;
//				marker = ghosts[idx].getMarker();
//				if (markerfather < ghosts[idx].getMarker()+1){
//					markerfather = ghosts[idx].getMarker()+1;
//				}
//				idx++;
//				if(idx == size_ghosts){
//					break;
//				}
//			}
//			nend = 0;
//			idx = nocts-1;
//			marker = octants[idx].getMarker();
//			while(marker < 0 & octants[idx].buildFather() == father & idx >= 0){
//				nbro++;
//				marker = octants[idx].getMarker();
//				if (markerfather < octants[idx+offset+idx2].getMarker()+1){
//					markerfather = octants[idx+offset+idx2].getMarker()+1;
//				}
//				nend++;
//				idx--;
//				if (idx<0){
//					break;
//				}
//			}
//			if (nbro == nchildren){
//				offset = nend;
//			}
//			else{
//				nend = 0;
//				for(int ii=nocts-nchildren; ii<nocts; ii++){
//					octants[ii].setMarker(0);
//				}
//			}
//		}
//		if (nend != 0){
//			for (idx=0; idx < nend; idx++){
//				for (iface=0; iface<nfaces; iface++){
//					father.info[iface] = (father.info[iface] || octants[nocts-idx].info[iface]);
//					father.info[iface+nfaces] = (father.info[iface+nfaces] || octants[nocts-idx].info[iface+nfaces]);
//				}
//			}
//			father.info[13] = true;
//			father.setMarker(markerfather);
//			if (markerfather < 0){
//				docoarse = true;
//			}
//			octants.resize(nocts-offset);
//			octants.push_back(father);
//			octants.shrink_to_fit();
//			nocts = octants.size();
//		}
//
//	}
//
//	//Update pborders (adesso inefficiente, loop di nuovo su tutti gli elementi)
//	//Si può trovare la maniera di inserirlo nel loop precedente
//	pborders.clear();
//	nocts = octants.size();
//	pborders.reserve(int(pow(double(nocts),2.0/3.0)*double(nfaces)));
//	for(idx=0; idx<nocts; idx++){
//		for(iface=0; iface<nfaces; iface++){
//			if (octants[idx].info[iface+nfaces]){
//				pborders.push_back(idx);
//				break;
//			}
//		}
//	}
//	pborders.shrink_to_fit();
//
//	// Set final first and last desc
//	if(nocts>0){
//		setFirstDesc();
//		setLastDesc();
//	}
//	return docoarse;
//}
//
//// =================================================================================== //
//
//bool Class_Local_Tree::refine(u32vector & mapidx) {
//
//	// Local variables
//	vector<uint32_t> last_child_index;
//	Class_Octant* children;
//	uint32_t idx, nocts;
//	uint32_t offset = 0, blockidx;
//	uint8_t nchm1 = nchildren-1, ich, iface;
//	bool dorefine = false;
//
//	nocts = octants.size();
//	for (idx=0; idx<nocts; idx++){
//		if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL){
//			last_child_index.push_back(idx+nchm1+offset);
//			offset += nchm1;
//		}
//		else{
////			octants[idx].info[12] = false;
//			if (octants[idx].marker > 0)
//				octants[idx].marker = 0;
//		}
//	}
//	if (offset > 0){
//		mapidx.resize(octants.size()+offset);
//		mapidx.shrink_to_fit();
//
//		octants.resize(octants.size()+offset);
//		blockidx = last_child_index[0]-nchm1;
//		idx = octants.size();
//		//while (idx>blockidx){
//		while (idx>0){
//			idx--;
//			//TODO Sostituire questo if con il controllo su last_index_child
//			if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL){
//				children = octants[idx-offset].buildChildren();
//				for (ich=0; ich<nchildren; ich++){
//					octants[idx-ich] = (children[nchm1-ich]);
//					mapidx[idx-ich]  = mapidx[idx-offset];
//				}
//				offset -= nchm1;
//				idx -= nchm1;
//				//Update local max depth
//				if (children[0].getLevel() > local_max_depth){
//					local_max_depth = children[0].getLevel();
//				}
//				if (children[0].getMarker() > 0){
//					//More Refinement to do
//					dorefine = true;
//				}
//				delete []children;
//			}
//			else {
//				octants[idx] = octants[idx-offset];
//				mapidx[idx]  = mapidx[idx-offset];
//			}
//		}
//	}
//	octants.shrink_to_fit();
//
//	//Update pborders (adesso inefficiente, loop di nuovo su tutti gli elementi)
//	//Si può trovare la maniera di inserirlo nel loop precedente
//	pborders.clear();
//	nocts = octants.size();
//	pborders.reserve(int(pow(double(nocts),2.0/3.0)*double(nfaces)));
//	for(idx=0; idx<nocts; idx++){
//		for(iface=0; iface<nfaces; iface++){
//			if (octants[idx].info[iface+nfaces]){
//				pborders.push_back(idx);
//				break;
//			}
//		}
//	}
//	pborders.shrink_to_fit();
//
//	setFirstDesc();
//	setLastDesc();
//
//	return dorefine;
//}
//
//// =================================================================================== //
//
//bool Class_Local_Tree::coarse(u32vector & mapidx) {
//
//	// Local variables
//	vector<uint32_t> first_child_index;
//	Class_Octant father;
//	uint32_t ich, nocts, nghosts, nocts0;
//	int32_t idx, idx2;
//	uint32_t offset;
//	int32_t idx1_gh, idx2_gh;
//	uint32_t nidx;
//	int8_t markerfather, marker;
//	uint8_t nbro, nstart, nend;
//	uint8_t nchm1 = nchildren-1, iface;
//	bool docoarse = false;
//
//	//------------------------------------------ //
//	// Initialization
//
//	nbro = nstart = nend = 0;
//	nidx = offset = 0;
//
//	idx1_gh = idx2_gh = 0;
//
//	nocts = nocts0 = octants.size();
//	size_ghosts = ghosts.size();
//
//
//	// Init first and last desc (even if already calculated)
//	setFirstDesc();
//	setLastDesc();
//
//	//------------------------------------------ //
//
//	// Set index for start and end check for ghosts
//	if (ghosts.size()){
//		while(idx1_gh < size_ghosts && ghosts[idx1_gh].computeMorton() < first_desc.computeMorton()){
//			idx1_gh++;
//		}
//		idx1_gh = max(0, idx1_gh-1);
//		while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() < last_desc.computeMorton()){
//			idx2_gh++;
//		}
//		idx2_gh = min(int(size_ghosts-1), idx2_gh);
//
//		// Start on ghosts
//		if ((ghosts[idx1_gh].getMarker() < 0) & (octants[0].getMarker() < 0)){
//			father = ghosts[idx1_gh].buildFather();
//			nbro = 0;
//			idx = idx1_gh;
//			marker = ghosts[idx].getMarker();
//			while(marker < 0 & ghosts[idx].buildFather() == father){
//				nbro++;
//				marker = ghosts[idx].getMarker();
//				idx--;
//				if (idx<0){
//					break;
//				}
//			}
//			nstart = 0;
//			idx = 0;
//			marker = octants[idx].getMarker();
//			while(marker<0 & octants[idx].buildFather() == father){
//				nbro++;
//				marker = octants[idx].getMarker();
//				nstart++;
//				idx++;
//				if (idx==nocts){
//					break;
//				}
//			}
//			if (nbro == nchildren){
////				offset = nstart;
//				// For update pbound of neighbours only check
//				// the odd faces of new father (placed nstart-times
//				// in the first nstart positions of octants)
//				// If there is father after coarse will be the first
//				// element of local octants (lowest Morton)
///*
//				for (int i=0; i<nstart; i++){
//					octants[i] = father;
//				}
//*/
//				uint32_t	 sizeneigh;
//				u32vector    neigh;
//				vector<bool> isghost;
//				for (iface=0; iface<DIM; iface++){
//					uint8_t oddface = ((iface*2)+1);
//					findNeighbours(nstart-1, oddface, neigh, isghost);
//					sizeneigh = neigh.size();
//					for(int i=0; i<sizeneigh; i++){
//						if (!isghost[i])
//							octants[neigh[i]].setPbound(oppface[oddface], true);
//					}
//				}
//			}
//			else{
//				nstart = 0;
//			}
//		}
//	}
//
//	// Check and coarse internal octants
//	for (idx=0; idx<nocts; idx++){
//		if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
//			nbro = 0;
//			father = octants[idx].buildFather();
//			// Check if family is to be refined
//			for (idx2=idx; idx2<idx+nchildren; idx2++){
//				if (idx2<nocts-1){
//					if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
//						nbro++;
//					}
//				}
//			}
//			if (nbro == nchildren){
//				nidx++;
//				first_child_index.push_back(idx);
//				idx = idx2-1;
//			}
//			else{
//				if (idx < (nocts>nchildren)*(nocts-nchildren)){
//					octants[idx].setMarker(0);
//				}
//			}
//		}
//		else{
////			octants[idx].info[13] = false;
//		}
//	}
//	//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
//	if (nidx!=0){
//		uint32_t nblock = nocts - nidx*nchm1 - nstart;
//		nidx = 0;
//		//for (idx=0; idx<nblock; idx++){
//		for (idx=0; idx<nocts-offset; idx++){
//			if (idx+offset == first_child_index[nidx]){
//				markerfather = -MAX_LEVEL;
//				father = octants[idx+offset].buildFather();
//				for(idx2=0; idx2<nchildren; idx2++){
//					if (markerfather < octants[idx+offset+idx2].getMarker()+1){
//						markerfather = octants[idx+offset+idx2].getMarker()+1;
//					}
//					for (iface=0; iface<nfaces; iface++){
//						father.info[iface] = (father.info[iface] || octants[idx+offset+idx2].info[iface]);
//						father.info[iface+nfaces] = (father.info[iface+nfaces] || octants[idx+offset+idx2].info[iface+nfaces]);
//					}
//				}
//				father.info[13] = true;
//				father.setMarker(markerfather);
//				if (markerfather < 0){
//					docoarse = true;
//				}
//				octants[idx] = father;
//				mapidx[idx] = mapidx[idx+offset];
//				offset += nchm1;
//				nidx++;
//			}
//			else{
//				octants[idx] = octants[idx+offset];
//				mapidx[idx] = mapidx[idx+offset];
//			}
//		}
//	}
//	octants.resize(nocts-offset);
//	octants.shrink_to_fit();
//	nocts = octants.size();
//	mapidx.resize(nocts-offset);
//	mapidx.shrink_to_fit();
//
//
//	// End on ghosts
//	if (ghosts.size() && nocts > 0){
//		if ((ghosts[idx2_gh].getMarker() < 0) & (octants[nocts-1].getMarker() < 0)){
//			father = ghosts[idx2_gh].buildFather();
//			markerfather = -MAX_LEVEL;
//			nbro = 0;
//			idx = idx2_gh;
//			marker = ghosts[idx].getMarker();
//			while(marker < 0 & ghosts[idx].buildFather() == father){
//				nbro++;
//				marker = ghosts[idx].getMarker();
//				if (markerfather < ghosts[idx].getMarker()+1){
//					markerfather = ghosts[idx].getMarker()+1;
//				}
//				idx++;
//				if(idx == size_ghosts){
//					break;
//				}
//			}
//			nend = 0;
//			idx = nocts-1;
//			marker = octants[idx].getMarker();
//			while(marker < 0 & octants[idx].buildFather() == father & idx >= 0){
//				nbro++;
//				marker = octants[idx].getMarker();
//				if (markerfather < octants[idx+offset+idx2].getMarker()+1){
//					markerfather = octants[idx+offset+idx2].getMarker()+1;
//				}
//				nend++;
//				idx--;
//				if (idx<0){
//					break;
//				}
//			}
//			if (nbro == nchildren){
//				offset = nend;
//			}
//			else{
//				nend = 0;
//				for(int ii=nocts-nchildren; ii<nocts; ii++){
//					octants[ii].setMarker(0);
//				}
//			}
//		}
//		if (nend != 0){
//			for (idx=0; idx < nend; idx++){
//				for (iface=0; iface<nfaces; iface++){
//					father.info[iface] = (father.info[iface] || octants[nocts-idx].info[iface]);
//					father.info[iface+nfaces] = (father.info[iface+nfaces] || octants[nocts-idx].info[iface+nfaces]);
//				}
//			}
//			father.info[13] = true;
//			father.setMarker(markerfather);
//			if (markerfather < 0){
//				docoarse = true;
//			}
//			octants.resize(nocts-offset);
//			octants.push_back(father);
//			octants.shrink_to_fit();
//			nocts = octants.size();
//			mapidx.resize(nocts-offset);
//			mapidx.push_back(nocts0-nend);
//			mapidx.shrink_to_fit();
//		}
//
//	}
//
//	//Update pborders (adesso inefficiente, loop di nuovo su tutti gli elementi)
//	//Si può trovare la maniera di inserirlo nel loop precedente
//	pborders.clear();
//	nocts = octants.size();
//	pborders.reserve(int(pow(double(nocts),2.0/3.0)*double(nfaces)));
//	for(idx=0; idx<nocts; idx++){
//		for(iface=0; iface<nfaces; iface++){
//			if (octants[idx].info[iface+nfaces]){
//				pborders.push_back(idx);
//				break;
//			}
//		}
//	}
//	pborders.shrink_to_fit();
//
//	// Set final first and last desc
//	if(nocts>0){
//		setFirstDesc();
//		setLastDesc();
//	}
//	return docoarse;
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::checkCoarse(uint64_t lastDescPre,uint64_t firstDescPost){
//	int32_t idx;
//	uint32_t nocts;
//	uint64_t Morton;
//	uint8_t toDelete = 0;
//
//	nocts = getNumOctants();
//	idx = 0;
//	Morton = octants[idx].computeMorton();
//	while(Morton <= lastDescPre & idx < nocts & Morton != 0){
//		// To delete, the father is in proc before me
//		toDelete++;
//		idx++;
//		Morton = octants[idx].computeMorton();
//	}
//	for(idx=0; idx<nocts-toDelete; idx++){
//		octants[idx] = octants[idx+toDelete];
//	}
//	octants.resize(nocts-toDelete);
//	octants.shrink_to_fit();
//
//	toDelete = 0;
//	Morton = last_desc.computeMorton();
//	if(int(firstDescPost  - Morton) > 1){
//		// To insert, the father is not yet here!!
//		idx = nocts - 1;
//		Class_Octant father = octants[idx].buildFather();
//		while(octants[idx].buildFather() == father & idx >= 0){
//			toDelete++;
//			idx--;
//		}
//		father.info[nfaces+1] = father.info[nfaces+3] = true;
//		if(nfaces == 6)
//			father.info[nfaces+5] = true;
//		octants.resize(nocts-toDelete);
//		octants.push_back(father);
//		octants.shrink_to_fit();
//	}
//
//	//Update pborders (adesso inefficiente, loop di nuovo su tutti gli elementi)
//	//Si può trovare la maniera di inserirlo nel loop precedente
//	uint8_t iface;
//	pborders.clear();
//	nocts = octants.size();
//	pborders.reserve(int(pow(double(nocts),2.0/3.0)*double(nfaces)));
//	for(idx=0; idx<nocts; idx++){
//		for(iface=0; iface<nfaces; iface++){
//			if (octants[idx].info[iface+nfaces]){
//				pborders.push_back(idx);
//				break;
//			}
//		}
//	}
//	pborders.shrink_to_fit();
//
//	setFirstDesc();
//	setLastDesc();
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::computeConnectivity() {
//	map<uint64_t, vector<uint32_t> > mapnodes;
//	map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
//	uint32_t i, k, counter;
//	uint64_t morton;
//	uint32_t noctants = getNumOctants();
//	uint32_t (*octnodes)[DIM];
//	uint8_t j;
//
//	//TODO Reserve for vector for 2D and 3D
//
//	if (nodes.size() == 0){
//		connectivity.resize(noctants);
//		for (i = 0; i < noctants; i++){
//			octnodes = octants[i].getNodes();
//			for (j = 0; j < nnodes; j++){
//#if DIM == 3
//				morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
//#else
//#endif
//				if (mapnodes[morton].size()==0){
//					mapnodes[morton].reserve(12);
//					for (k = 0; k < DIM; k++){
//						mapnodes[morton].push_back(octnodes[j][k]);
//					}
//				}
//				mapnodes[morton].push_back(i);
//			}
//			delete []octnodes;
//		}
//		iter	= mapnodes.begin();
//		iterend	= mapnodes.end();
//		counter = 0;
//		uint32_t numnodes = mapnodes.size();
//		nodes.resize(numnodes);
//		vector<uint32_t>::iterator iter2;
//		while (iter != iterend){
//			vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+DIM);
//			nodes[counter] = nodecasting;
//			nodes[counter].shrink_to_fit();
//			for(iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
//				if (connectivity[(*iter2)].size()==0){
//					connectivity[(*iter2)].reserve(8);
//				}
//				connectivity[(*iter2)].push_back(counter);
//			}
//			mapnodes.erase(iter++);
//			counter++;
//		}
//		nodes.shrink_to_fit();
//		//Lento. Solo per risparmiare memoria
//		for (int ii=0; ii<noctants; ii++){
//			connectivity[ii].shrink_to_fit();
//		}
//		connectivity.shrink_to_fit();
//	}
//	map<uint64_t, vector<uint32_t> >().swap(mapnodes);
//	iter = mapnodes.end();
//}
//
//void Class_Local_Tree::clearConnectivity() {
//	u32vector2D().swap(nodes);
//	u32vector2D().swap(connectivity);
//}
//
//void Class_Local_Tree::updateConnectivity() {
//	clearConnectivity();
//	computeConnectivity();
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::computeghostsConnectivity() {
//	map<uint64_t, vector<uint32_t> > mapnodes;
//	map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
//	uint32_t i, k, counter;
//	uint64_t morton;
//	uint32_t noctants = size_ghosts;
//	uint32_t (*octnodes)[DIM];
//	uint8_t j;
//
//	if (ghostsnodes.size() == 0){
//		ghostsconnectivity.resize(noctants);
//		for (i = 0; i < noctants; i++){
//			octnodes = ghosts[i].getNodes();
//			for (j = 0; j < nnodes; j++){
//#if DIM == 3
//				morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
//#else
//				puppa
//#endif
//				if (mapnodes[morton].size()==0){
//					for (k = 0; k < DIM; k++){
//						mapnodes[morton].push_back(octnodes[j][k]);
//					}
//				}
//				mapnodes[morton].push_back(i);
//			}
//			delete []octnodes;
//		}
//		iter	= mapnodes.begin();
//		iterend	= mapnodes.end();
//		uint32_t numnodes = mapnodes.size();
//		ghostsnodes.resize(numnodes);
//		counter = 0;
//		vector<uint32_t>::iterator iter2;
//		while (iter != iterend){
//			vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+DIM);
//			ghostsnodes[counter] = nodecasting;
//			ghostsnodes[counter].shrink_to_fit();
//			for(iter2 = iter->second.begin()+DIM; iter2 != iter->second.end(); iter2++){
//				if (ghostsconnectivity[(*iter2)].size()==0){
//					ghostsconnectivity[(*iter2)].reserve(8);
//				}
//				ghostsconnectivity[(*iter2)].push_back(counter);
//			}
//			mapnodes.erase(iter++);
//			counter++;
//		}
//		ghostsnodes.shrink_to_fit();
//		//Lento. Solo per risparmiare memoria
//		for (int ii=0; ii<noctants; ii++){
//			ghostsconnectivity[ii].shrink_to_fit();
//		}
//		ghostsconnectivity.shrink_to_fit();
//	}
//	iter = mapnodes.end();
//}
//
//void Class_Local_Tree::clearghostsConnectivity() {
//	u32vector2D().swap(ghostsnodes);
//	u32vector2D().swap(ghostsconnectivity);
//}
//
//void Class_Local_Tree::updateghostsConnectivity() {
//	clearghostsConnectivity();
//	computeghostsConnectivity();
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::updateLocalMaxDepth() {
//	uint32_t noctants = getNumOctants();
//	uint32_t i;
//
//	local_max_depth = 0;
//	for(i = 0; i < noctants; i++){
//		if(octants[i].getLevel() > local_max_depth){
//			local_max_depth = octants[i].getLevel();
//		}
//	}
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::findNeighbours(uint32_t const idx, uint8_t iface,
//		u32vector & neighbours, vector<bool> & isghost) {
//
//	uint64_t  Morton, Mortontry;
//	uint32_t  noctants = getNumOctants();
//	uint32_t idxtry;
//	Class_Octant* oct = &octants[idx];
//	uint32_t size = oct->getSize();
//
//	//Alternative to switch case
//	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
//	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
//	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));
//
//	isghost.clear();
//	neighbours.clear();
//
//	// Default if iface is nface<iface<0
//	if (iface < 0 || iface > nfaces){
//		writeLog("Face index out of range in find neighbours !!!");
//		return;
//	}
//
//	// Check if octants face is a process boundary
//	if (oct->info[nfaces+iface] == false){
//
//		// Check if octants face is a boundary
//		if (oct->info[iface] == false){
//
//			//Build Morton number of virtual neigh of same size
//			Class_Octant samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
//			Morton = samesizeoct.computeMorton();
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
//			idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
//			Mortontry = oct->computeMorton();
//			while(abs(jump) > 0){
//				Mortontry = octants[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//				if (idxtry > noctants-1){
//					if (jump > 0){
//						idxtry = noctants - 1;
//						jump = 0;
//					}
//					else if (jump < 0){
//						idxtry = 0;
//						jump = 0;
//					}
//				}
//			}
//			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//				//Found neighbour of same size
//				isghost.push_back(false);
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(octants[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > noctants-1){
//							idxtry = noctants-1;
//							return;
//						}
//					}
//					while(octants[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > noctants-1){
//							idxtry = noctants-1;
//							return;
//						}
//					}
//				}
//				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//					//Found neighbour of same size
//					isghost.push_back(false);
//					neighbours.push_back(idxtry);
//					return;
//				}
//				// Compute Last discendent of virtual octant of same size
//				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//				Class_Octant last_desc = samesizeoct.buildLastDesc();
//				uint64_t Mortonlast = last_desc.computeMorton();
//				vector<uint32_t> bufferidx;
//				Mortontry = octants[idxtry].computeMorton();
//				int32_t Dh;
//				int32_t eqcoord;
//				while(Mortontry < Mortonlast & idxtry < noctants){
//					Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
//					Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
//					Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
//					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
//						neighbours.push_back(idxtry);
//						isghost.push_back(false);
//					}
//					idxtry++;
//					if(idxtry>noctants-1){
//						break;
//					}
//					Mortontry = octants[idxtry].computeMorton();
//				}
//				return;
//			}
//		}
//		else{
//			// Boundary Face
//			return;
//		}
//	}
//	//--------------------------------------------------------------- //
//	//--------------------------------------------------------------- //
//	else{
//		cout << "idx   " << int(idx) << "  iface  " << int(iface) << " is pbound" << endl;
//		// Check if octants face is a boundary
//		if (oct->info[iface] == false){
//			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS
//
//			if (ghosts.size()>0){
//				// Search in ghosts
//
//				uint32_t idxghost = uint32_t(size_ghosts/2);
//				Class_Octant* octghost = &ghosts[idxghost];
//
//				//Build Morton number of virtual neigh of same size
//				Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
//				// Search morton in octants
//				// If a even face morton is lower than morton of oct, if odd higher
//				// ---> can i search only before or after idx in octants
//				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
//				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
//				while(abs(jump) > 0){
//					Mortontry = ghosts[idxtry].computeMorton();
//					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//					idxtry += jump;
//					if (idxtry > ghosts.size()-1){
//						if (jump > 0){
//							idxtry = ghosts.size() - 1;
//							jump = 0;
//						}
//						else if (jump < 0){
//							idxtry = 0;
//							jump = 0;
//						}
//					}
//				}
//				if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//					//Found neighbour of same size
//					isghost.push_back(true);
//					neighbours.push_back(idxtry);
//					return;
//				}
//				else{
//					// Step until the mortontry lower than morton (one idx of distance)
//					{
//						while(ghosts[idxtry].computeMorton() < Morton){
//							idxtry++;
//							if(idxtry > ghosts.size()-1){
//								idxtry = ghosts.size();
//								break;
//							}
//						}
//						while(ghosts[idxtry].computeMorton() > Morton){
//							idxtry--;
//							if(idxtry > ghosts.size()-1){
//								idxtry = ghosts.size();
//								break;
//							}
//						}
//					}
//					if(idxtry < size_ghosts){
//						if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//							//Found neighbour of same size
//							isghost.push_back(true);
//							neighbours.push_back(idxtry);
//							return;
//						}
//						// Compute Last discendent of virtual octant of same size
//						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//						Class_Octant last_desc = samesizeoct.buildLastDesc();
//						uint64_t Mortonlast = last_desc.computeMorton();
//						vector<uint32_t> bufferidx;
//						Mortontry = ghosts[idxtry].computeMorton();
//						int32_t Dh;
//						int32_t eqcoord;
//						while(Mortontry < Mortonlast & idxtry < size_ghosts){
//							Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
//							Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
//							Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
//							if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
//								neighbours.push_back(idxtry);
//								isghost.push_back(true);
//							}
//							idxtry++;
//							if(idxtry>size_ghosts-1){
//								break;
//							}
//							Mortontry = ghosts[idxtry].computeMorton();
//						}
//					}
//				}
//
//				uint32_t lengthneigh = 0;
//				uint32_t sizeneigh = neighbours.size();
//				for (idxtry=0; idxtry<sizeneigh; idxtry++){
//					lengthneigh += ghosts[neighbours[idxtry]].getSize();
//				}
//				if (lengthneigh < oct->getSize()){
//					// Search in octants
//
//					// Check if octants face is a boundary
//					if (oct->info[iface] == false){
//
//						//Build Morton number of virtual neigh of same size
//						Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//						Morton = samesizeoct.computeMorton();
//						// Search morton in octants
//						// If a even face morton is lower than morton of oct, if odd higher
//						// ---> can i search only before or after idx in octants
//						int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
//						idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
//						while(abs(jump) > 0){
//							Mortontry = octants[idxtry].computeMorton();
//							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//							idxtry += jump;
//						}
//						if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//							//Found neighbour of same size
//							isghost.push_back(false);
//							neighbours.push_back(idxtry);
//							cout << " idx  " << idx << "   idxtry  "<<idxtry << "/" << noctants << "   face   " << int(iface) << endl;
//							writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
//							return;
//						}
//						else{
//							// Step until the mortontry lower than morton (one idx of distance)
//							{
//								while(octants[idxtry].computeMorton() < Morton){
//									idxtry++;
//									if(idxtry > noctants-1){
//										idxtry = noctants;
//										break;
//									}
//								}
//								while(octants[idxtry].computeMorton() > Morton){
//									idxtry--;
//									if(idxtry > noctants-1){
//										idxtry = noctants;
//										break;
//									}
//								}
//							}
//							if (idxtry < noctants){
//								if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//									//Found neighbour of same size
//									isghost.push_back(false);
//									neighbours.push_back(idxtry);
//									cout << " idx  " << idx << "   idxtry  "<<idxtry << "/" << noctants << "   face   " << int(iface) << endl;
//									writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
//									return;
//								}
//								// Compute Last discendent of virtual octant of same size
//								uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//								Class_Octant last_desc = samesizeoct.buildLastDesc();
//								uint64_t Mortonlast = last_desc.computeMorton();
//								vector<uint32_t> bufferidx;
//								Mortontry = octants[idxtry].computeMorton();
//								int32_t Dh;
//								int32_t eqcoord;
//								while(Mortontry < Mortonlast & idxtry < noctants-1){
//									Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
//									Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
//									Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
//									if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
//										neighbours.push_back(idxtry);
//										isghost.push_back(false);
//									}
//									idxtry++;
//									Mortontry = octants[idxtry].computeMorton();
//								}
//							}
//						}
//					}
//				}
//				return;
//			}
//		}
//		else{
//			// Boundary Face
//			return;
//		}
//
//	}
//}
//
//// =================================================================================== //
//
////void Class_Local_Tree::findNeighbours(Class_Octant const & oct, uint8_t iface,
////		u32vector & neighbours, vector<bool> & isghost) {
////
////	//TODO DA RIVEDERE E RIFARE !!!
////	uint64_t  Morton, Mortontry;
////	uint32_t  noctants = getNumOctants();
////	uint32_t idxtry, idxtry_old, idxtry_old_;
////	uint32_t size = oct.getSize();
////
////	//Alternative to switch case
////	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
////	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
////	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));
////
////	isghost.clear();
////	neighbours.clear();
////
////	// Default if iface is nface<iface<0
////	if (iface < 0 || iface > nface){
////		writeLog("Face index out of range in find neighbours !!!");
////		return;
////	}
////
////	// Check if octants face is a process boundary
////	if (oct.info[nface+iface] == false){
////
////		// Check if octants face is a boundary
////		if (oct.info[iface] == false){
////
////			// Find octant in octants
////			uint32_t idx;
////			{
////			vector<Class_Octant>::iterator it = find(octants.begin(), octants.end(), oct);
////			idx = distance(octants.begin(), it);
////			it = octants.end();
////			}
////
////			//Build Morton number of virtual neigh of same size
////			Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
////			Morton = samesizeoct.computeMorton();
////			// Search morton in octants
////			// If a even face morton is lower than morton of oct, if odd higher
////			// --. can i search only before or after idx in octants
////			int32_t jump = (oct.computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
////			idxtry = uint32_t(idx +((oct.computeMorton()<Morton)-(oct.computeMorton()>Morton))*jump);
////			//idxtry_old = uint64_t((1+direction)*noctants );
////			while(abs(jump) > 0){
////				Mortontry = octants[idxtry].computeMorton();
////				jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
////				idxtry += jump;
////				if (idxtry > noctants-1){
////					if (jump > 0){
////						idxtry = noctants - 1;
////						jump = 0;
////					}
////					else if (jump < 0){
////						idxtry = 0;
////						jump = 0;
////					}
////				}
////			}
////			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
////				//Found neighbour of same size
////				isghost.push_back(false);
////				neighbours.push_back(idxtry);
////				return;
////			}
////			else{
////				// Step until the mortontry lower than morton (one idx of distance)
////				{
////					while(octants[idxtry].computeMorton() < Morton){
////						idxtry++;
////					}
////					while(octants[idxtry].computeMorton() > Morton){
////						idxtry--;
////					}
////				}
////				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
////					//Found neighbour of same size
////					isghost.push_back(false);
////					neighbours.push_back(idxtry);
////					return;
////				}
////				// Compute Last discendent of virtual octant of same size
////				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
////				Class_Octant last_desc = samesizeoct.buildLastDesc();
////				uint64_t Mortonlast = last_desc.computeMorton();
////				vector<uint32_t> bufferidx;
////				Mortontry = octants[idxtry].computeMorton();
////				int32_t Dh;
////				int32_t eqcoord;
////				while(Mortontry < Mortonlast & idxtry < noctants-1){
////					Dh = int32_t(cx)*(int32_t(oct.x) - int32_t(octants[idxtry].x));
////					Dh += int32_t(cy)*(int32_t(oct.y) - int32_t(octants[idxtry].y));
////					Dh += int32_t(cz)*(int32_t(oct.z) - int32_t(octants[idxtry].z));
////					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
////						neighbours.push_back(idxtry);
////						isghost.push_back(false);
////					}
////					idxtry++;
////					Mortontry = octants[idxtry].computeMorton();
////				}
////				return;
////			}
////		}
////		else{
////			// Boundary Face
////			return;
////		}
////	}
////	//--------------------------------------------------------------- //
////	//--------------------------------------------------------------- //
////	else{
////		// Check if octants face is a boundary
////		if (oct.info[iface] == false){
////
////			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS
////
////			// Find octant in octants
////			uint32_t idx;
////			{
////			vector<Class_Octant>::iterator it = find(octants.begin(), octants.end(), oct);
////			idx = distance(octants.begin(), it);
////			it = octants.end();
////			}
////
////			if (ghosts.size()>0){
////				// Search in ghosts
////
////				uint32_t idxghost = uint32_t(size_ghosts/2);
////				Class_Octant* octghost = &ghosts[idxghost];
////
////				//Build Morton number of virtual neigh of same size
////				Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
////				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct.x-size,oct.y,oct.z);
////				// Search morton in octants
////				// If a even face morton is lower than morton of oct, if odd higher
////				// --. can i search only before or after idx in octants
////				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
////				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
////				while(abs(jump) > 0){
////					Mortontry = ghosts[idxtry].computeMorton();
////					jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
////					idxtry += jump;
////					if (idxtry > ghosts.size()-1){
////						if (jump > 0){
////							idxtry = ghosts.size() - 1;
////							jump = 0;
////						}
////						else if (jump < 0){
////							idxtry = 0;
////							jump = 0;
////						}
////					}
////				}
////				if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct.level){
////					//Found neighbour of same size
////					isghost.push_back(true);
////					neighbours.push_back(idxtry);
////					return;
////				}
////				else{
////					// Step until the mortontry lower than morton (one idx of distance)
////					{
////						while(ghosts[idxtry].computeMorton() < Morton){
////							idxtry++;
////						}
////						while(ghosts[idxtry].computeMorton() > Morton){
////							idxtry--;
////						}
////					}
////					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct.level){
////						//Found neighbour of same size
////						isghost.push_back(true);
////						neighbours.push_back(idxtry);
////						return;
////					}
////					// Compute Last discendent of virtual octant of same size
////					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
////					Class_Octant last_desc = samesizeoct.buildLastDesc();
////					uint64_t Mortonlast = last_desc.computeMorton();
////					vector<uint32_t> bufferidx;
////					Mortontry = ghosts[idxtry].computeMorton();
////					int32_t Dh;
////					int32_t eqcoord;
////					while(Mortontry < Mortonlast & idxtry < size_ghosts-1){
////						Dh = int32_t(cx)*(int32_t(oct.x) - int32_t(ghosts[idxtry].x));
////						Dh += int32_t(cy)*(int32_t(oct.y) - int32_t(ghosts[idxtry].y));
////						Dh += int32_t(cz)*(int32_t(oct.z) - int32_t(ghosts[idxtry].z));
////						if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
////							neighbours.push_back(idxtry);
////							isghost.push_back(true);
////						}
////						idxtry++;
////						Mortontry = ghosts[idxtry].computeMorton();
////					}
////				}
////				uint32_t areaneigh = 0;
////				uint32_t sizeneigh = neighbours.size();
////				for (idx=0; idx<sizeneigh; idx++){
////					areaneigh += ghosts[neighbours[idx]].getArea();
////				}
////				if (areaneigh < oct.getArea()){
////					// Search in octants
////
////					// Check if octants face is a boundary
////					if (oct.info[iface] == false){
////
////						//Build Morton number of virtual neigh of same size
////						Class_Octant samesizeoct(oct.level, oct.x+cx*size, oct.y+cy*size, oct.z+cz*size);
////						Morton = samesizeoct.computeMorton();
////						// Search morton in octants
////						// If a even face morton is lower than morton of oct, if odd higher
////						// --. can i search only before or after idx in octants
////						int32_t jump = (oct.computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
////						idxtry = uint32_t(idx +((oct.computeMorton()<Morton)-(oct.computeMorton()>Morton))*jump);
////						//idxtry_old = uint64_t((1+direction)*noctants );
////						while(abs(jump) > 0){
////							Mortontry = octants[idxtry].computeMorton();
////							jump = ((Mortontry<Morton)-(Mortontry>Morton))*jump/2;
////							idxtry += jump;
////						}
////						if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
////							//Found neighbour of same size
////							isghost.push_back(false);
////							neighbours.push_back(idxtry);
////							writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
////							return;
////						}
////						else{
////							// Step until the mortontry lower than morton (one idx of distance)
////							{
////								while(octants[idxtry].computeMorton() < Morton){
////									idxtry++;
////								}
////								while(octants[idxtry].computeMorton() > Morton){
////									idxtry--;
////								}
////							}
////							if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct.level){
////								//Found neighbour of same size
////								isghost.push_back(false);
////								neighbours.push_back(idxtry);
////								writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
////								return;
////							}
////							// Compute Last discendent of virtual octant of same size
////							uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
////							Class_Octant last_desc = samesizeoct.buildLastDesc();
////							uint64_t Mortonlast = last_desc.computeMorton();
////							vector<uint32_t> bufferidx;
////							Mortontry = octants[idxtry].computeMorton();
////							int32_t Dh;
////							int32_t eqcoord;
////							while(Mortontry < Mortonlast & idxtry < noctants-1){
////								Dh = int32_t(cx)*(int32_t(oct.x) - int32_t(octants[idxtry].x));
////								Dh += int32_t(cy)*(int32_t(oct.y) - int32_t(octants[idxtry].y));
////								Dh += int32_t(cz)*(int32_t(oct.z) - int32_t(octants[idxtry].z));
////								if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
////									neighbours.push_back(idxtry);
////									isghost.push_back(false);
////								}
////								idxtry++;
////								Mortontry = octants[idxtry].computeMorton();
////							}
////						}
////					}
////				}
////				return;
////			}
////		}
////		else{
////			// Boundary Face
////			return;
////		}
////
////	}
////}
//
//// =================================================================================== //
//
//void Class_Local_Tree::findGhostNeighbours(uint32_t idx, uint8_t iface,
//		u32vector & neighbours) {
//
//	uint64_t  Morton, Mortontry;
//	uint32_t  noctants = getNumOctants();
//	uint32_t idxtry;
//	Class_Octant* oct = &ghosts[idx];
//	uint32_t size = oct->getSize();
//
//	//Alternative to switch case
//	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
//	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
//	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));
//
//	neighbours.clear();
//
//	// Default if iface is nface<iface<0
//	if (iface < 0 || iface > nfaces){
//		writeLog("Face index out of range in find neighbours !!!");
//		return;
//	}
//
//	// Check if octants face is a process boundary
//	if (oct->info[nfaces+iface] == true){
//
//			//Build Morton number of virtual neigh of same size
//			Class_Octant samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
//			Morton = samesizeoct.computeMorton();
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump;
//			idxtry = uint32_t(getNumOctants()/2);
//			Mortontry = octants[idxtry].computeMorton();
//			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//			while(abs(jump) > 0){
//				Mortontry = octants[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//				if (idxtry > noctants-1){
//					if (jump > 0){
//						idxtry = noctants - 1;
//						jump = 0;
//					}
//					else if (jump < 0){
//						idxtry = 0;
//						jump = 0;
//					}
//				}
//			}
//			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//				//Found neighbour of same size
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(octants[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > noctants-1){
//							idxtry = noctants-1;
//							return;
//						}
//					}
//					while(octants[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > noctants-1){
//							idxtry = noctants-1;
//							return;
//						}
//					}
//				}
//				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//					//Found neighbour of same size
//					neighbours.push_back(idxtry);
//					return;
//				}
//				// Compute Last discendent of virtual octant of same size
//				uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//				Class_Octant last_desc = samesizeoct.buildLastDesc();
//				uint64_t Mortonlast = last_desc.computeMorton();
//				vector<uint32_t> bufferidx;
//				Mortontry = octants[idxtry].computeMorton();
//				int32_t Dh;
//				int32_t eqcoord;
//				while(Mortontry < Mortonlast & idxtry < noctants){
//					Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
//					Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
//					Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
//					if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
//						neighbours.push_back(idxtry);
//					}
//					idxtry++;
//					if(idxtry>noctants-1){
//						break;
//					}
//					Mortontry = octants[idxtry].computeMorton();
//				}
//				return;
//			}
//	}
//	//--------------------------------------------------------------- //
//	//-----Not Pbound face------------- ----------------------------- //
//	else{
//		return;
//	}
//}
//
//// =================================================================================== //
//
//bool Class_Local_Tree::localBalance(bool doInterior){
//
//	// Local variables
//	uint32_t 			noctants = getNumOctants();
//	uint32_t			sizeneigh, modsize;
//	u32vector		 	neigh;
//	u32vector		 	modified, newmodified;
//	uint32_t 			i, idx, imod;
//	uint8_t				iface;
//	int8_t				targetmarker;
//	vector<bool> 		isghost;
//	bool				Bdone = false;
//
//	OctantsType::iterator 	obegin, oend, it;
//	u32vector::iterator 	ibegin, iend, iit;
//
//	//If interior octants have to be balanced
//	if(doInterior){
//		// First loop on the octants
//		/*		for(idx=0 ; idx<noctants; idx++){
//			if (!octants[idx].getNotBalance()){
//				for (iface=1; iface<nface; iface+=2){
//					if(!octants[idx].getPbound(iface)){
//						findNeighbours(idx, iface, neigh, isghost);
//						sizeneigh = neigh.size();
//						for(i=0; i<sizeneigh; i++){
//							if (!isghost[i]){
//								{
//									if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (octants[idx].getLevel() + octants[idx].getMarker() + 1) ){
//										octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
//										modified.push_back(idx);
//										Bdone = true;
//									}
//									else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (octants[idx].getLevel() + octants[idx].getMarker() - 1)){
//										octants[neigh[i]].setMarker(octants[idx].getLevel()+octants[idx].getMarker()-octants[neigh[i]].getLevel()-1);
//										modified.push_back(neigh[i]);
//										Bdone = true;
//									}
//								};
//							}
//
//						else{
//							if(ghosts.size()>0){
//								if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker())> (octants[idx].getLevel() + octants[idx].getMarker() + 1)){
//									octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-octants[idx].getLevel()-1);
//									modified.push_back(idx);
//									Bdone = true;
//								}
//							}
//						}
//
//						}
//					}
//				}
//			}
//		}*/
//		obegin = octants.begin();
//		oend = octants.end();
//		idx = 0;
//		for (it=obegin; it!=oend; it++){
//			it->info[15] = false;
//			if (!it->getNotBalance() && it->getMarker() != 0){
//				targetmarker = min(MAX_LEVEL, (octants[idx].getLevel() + octants[idx].getMarker()));
//				for (iface=0; iface<nfaces; iface++){
//					if(!it->getPbound(iface)){
//						findNeighbours(idx, iface, neigh, isghost);
//						sizeneigh = neigh.size();
//						for(i=0; i<sizeneigh; i++){
//							if (!isghost[i]){
//								{
//									if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
//										octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
//										octants[idx].info[15] = true;
//										modified.push_back(idx);
//										Bdone = true;
//									}
//									else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
//										octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
//										octants[neigh[i]].info[15] = true;
//										modified.push_back(neigh[i]);
//										Bdone = true;
//									}
//								};
//							}
//						}
//					}
//				}
//			}
//			idx++;
//		}
//		// Loop on ghost octants (influence over interior borders)
//		obegin = ghosts.begin();
//		oend = ghosts.end();
//		idx = 0;
//		for (it=obegin; it!=oend; it++){
//			if (!it->getNotBalance() && it->getMarker() != 0){
//				targetmarker = min(MAX_LEVEL, (it->getLevel()+it->getMarker()));
//				for (iface=0; iface<nfaces; iface++){
//					if(it->getPbound(iface) == true){
//						neigh.clear();
//						findGhostNeighbours(idx, iface, neigh);
//						sizeneigh = neigh.size();
//						for(i=0; i<sizeneigh; i++){
//							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
//								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
//								octants[neigh[i]].info[15] = true;
//								modified.push_back(neigh[i]);
//								Bdone = true;
//							}
//						}
//					}
//				}
//			}
//			idx++;
//		}
//
//		// While loop for iterative balancing
//		u32vector().swap(newmodified);
//		modsize = modified.size();
//		while(modsize!=0){
//			ibegin = modified.begin();
//			iend = modified.end();
//			for (iit=ibegin; iit!=iend; iit++){
//				idx = *iit;
//				if (!octants[idx].getNotBalance()){
//					targetmarker = min(MAX_LEVEL, (octants[idx].getLevel()+octants[idx].getMarker()));
//					for (iface=0; iface<nfaces; iface++){
//						if(!octants[idx].getPbound(iface)){
//							findNeighbours(idx, iface, neigh, isghost);
//							sizeneigh = neigh.size();
//							for(i=0; i<sizeneigh; i++){
//								if (!isghost[i]){
//									{
//										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
//											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
//											octants[idx].info[15] = true;
//											newmodified.push_back(idx);
//											Bdone = true;
//										}
//										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
//											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
//											octants[neigh[i]].info[15] = true;
//											newmodified.push_back(neigh[i]);
//											Bdone = true;
//										}
//									};
//								}
//							}
//						}
//					}
//				}
//			}
//			u32vector().swap(modified);
//			swap(modified,newmodified);
//			modsize = modified.size();
//		}// end while
//
//	}
//	else{
//
//		// Loop on ghost octants (influence over interior borders)
//		/*	for (idx=0; idx<size_ghosts; idx++){
//		if (!ghosts[idx].getNotBalance()){
//			for (iface=0; iface<nface; iface++){
//				if(ghosts[idx].getPbound(iface) == true){
//					neigh.clear();
//					findGhostNeighbours(idx, iface, neigh);
//					sizeneigh = neigh.size();
//					for(i=0; i<sizeneigh; i++){
//						if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (ghosts[idx].getLevel() + ghosts[idx].getMarker() - 1)){
//							octants[neigh[i]].setMarker(ghosts[idx].getLevel()+ghosts[idx].getMarker()-octants[neigh[i]].getLevel()-1);
//							modified.push_back(neigh[i]);
//							Bdone = true;
//						}
//					}
//				}
//			}
//		}
//	}*/
//		obegin = ghosts.begin();
//		oend = ghosts.end();
//		idx = 0;
//		for (it=obegin; it!=oend; it++){
//			if (!it->getNotBalance() && it->info[15]){
//				targetmarker = min(MAX_LEVEL, (it->getLevel()+it->getMarker()));
//				for (iface=0; iface<nfaces; iface++){
//					if(it->getPbound(iface) == true){
//						neigh.clear();
//						findGhostNeighbours(idx, iface, neigh);
//						sizeneigh = neigh.size();
//						for(i=0; i<sizeneigh; i++){
//							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
//								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
//								octants[neigh[i]].info[15] = true;
//								modified.push_back(neigh[i]);
//								Bdone = true;
//							}
//						}
//					}
//				}
//			}
//			idx++;
//		}
//
//		// While loop for iterative balancing
//		u32vector().swap(newmodified);
//		modsize = modified.size();
//		while(modsize!=0){
//			ibegin = modified.begin();
//			iend = modified.end();
//			for (iit=ibegin; iit!=iend; iit++){
//				idx = *iit;
//				if (!octants[idx].getNotBalance()){
//					targetmarker = min(MAX_LEVEL, (octants[idx].getLevel()+octants[idx].getMarker()));
//					for (iface=0; iface<nfaces; iface++){
//						if(!octants[idx].getPbound(iface)){
//							findNeighbours(idx, iface, neigh, isghost);
//							sizeneigh = neigh.size();
//							for(i=0; i<sizeneigh; i++){
//								if (!isghost[i]){
//									{
//										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
//											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
//											octants[idx].info[15] = true;
//											newmodified.push_back(idx);
//											Bdone = true;
//										}
//										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
//											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
//											octants[neigh[i]].info[15] = true;
//											newmodified.push_back(neigh[i]);
//											Bdone = true;
//										}
//									};
//								}
//							}
//						}
//					}
//				}
//			}
//			u32vector().swap(modified);
//			swap(modified,newmodified);
//			modsize = modified.size();
//		}// end while
//		obegin = oend = octants.end();
//		ibegin = iend = modified.end();
//	}
//	return Bdone;
//	// Pay attention : info[15] may be true after local balance for some octants
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::findEdgeNeighbours(uint32_t const idx, uint8_t iedge,
//		u32vector & neighbours, vector<bool> & isghost) {
//
//	uint64_t  Morton, Mortontry;
//	uint32_t  noctants = getNumOctants();
//	uint32_t idxtry;
//	Class_Octant* oct = &octants[idx];
//	uint32_t size = oct->getSize();
//	uint8_t iface1, iface2;
//	int32_t Dhx, Dhy, Dhz;
//	int32_t Dhxref, Dhyref, Dhzref;
//
//	//Alternative to switch case
//	int8_t Cx[12] = {-1,1,0,0,-1,1,-1,1,-1,1,0,0};
//	int8_t Cy[12] = {0,0,-1,1,-1,-1,1,1,0,0,-1,1};
//	int8_t Cz[12] = {-1,-1,-1,-1,0,0,0,0,1,1,1,1};
//	int8_t cx = Cx[iedge];
//	int8_t cy = Cy[iedge];
//	int8_t cz = Cz[iedge];
//
//	isghost.clear();
//	neighbours.clear();
//
//	// Default if iedge is nface<iedge<0
//	if (iedge < 0 || iedge > nfaces*2){
//		writeLog("Edge index out of range in find neighbours !!!");
//		return;
//	}
//
//	// Check if octants edge is a process boundary
//	iface1 = edgeface[iedge][0];
//	iface2 = edgeface[iedge][1];
//
//	// Check if octants edge is a boundary
//	if (oct->info[iface1] == false && oct->info[iface2] == false){
//
//		//Build Morton number of virtual neigh of same size
//		Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
//
//		//SEARCH IN GHOSTS
//
//		if (ghosts.size()>0){
//			// Search in ghosts
//
//			uint32_t idxghost = uint32_t(size_ghosts/2);
//			Class_Octant* octghost = &ghosts[idxghost];
//
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
//			idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
//			while(abs(jump) > 0){
//				Mortontry = ghosts[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//				if (idxtry > ghosts.size()-1){
//					if (jump > 0){
//						idxtry = ghosts.size() - 1;
//						jump = 0;
//					}
//					else if (jump < 0){
//						idxtry = 0;
//						jump = 0;
//					}
//				}
//			}
//			if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//				//Found neighbour of same size
//				isghost.push_back(true);
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(ghosts[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > ghosts.size()-1){
//							idxtry = ghosts.size();
//							break;
//						}
//					}
//					while(ghosts[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > ghosts.size()-1){
//							idxtry = ghosts.size() - 1;
//							break;
//						}
//						if(idxtry < 0){
//							idxtry = 0;
//							break;
//						}
//					}
//				}
//				if(idxtry < size_ghosts){
//					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//						//Found neighbour of same size
//						isghost.push_back(true);
//						neighbours.push_back(idxtry);
//						return;
//					}
//					// Compute Last discendent of virtual octant of same size
//					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//					Class_Octant last_desc = samesizeoct.buildLastDesc();
//					uint64_t Mortonlast = last_desc.computeMorton();
//					vector<uint32_t> bufferidx;
//					Mortontry = ghosts[idxtry].computeMorton();
//					while(Mortontry < Mortonlast & idxtry < ghosts.size()){
//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
//						Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
//						Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
//						Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
//						if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
//							neighbours.push_back(idxtry);
//							isghost.push_back(true);
//						}
//						idxtry++;
//						if(idxtry>size_ghosts-1){
//							break;
//						}
//						Mortontry = ghosts[idxtry].computeMorton();
//					}
//				}
//			}
//		}
//		uint32_t lengthneigh = 0;
//		uint32_t sizeneigh = neighbours.size();
//		for (idxtry=0; idxtry<sizeneigh; idxtry++){
//			lengthneigh += ghosts[neighbours[idxtry]].getSize();
//		}
//		if (lengthneigh < oct->getSize()){
//
//			// Search in octants
//
//			//Build Morton number of virtual neigh of same size
//			//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//			//Morton = samesizeoct.computeMorton();
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
//			idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
//			if (idxtry > noctants-1)
//				idxtry = noctants-1;
//			while(abs(jump) > 0){
//				Mortontry = octants[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//				if (idxtry > octants.size()-1){
//					if (jump > 0){
//						idxtry = octants.size() - 1;
//						jump = 0;
//					}
//					else if (jump < 0){
//						idxtry = 0;
//						jump = 0;
//					}
//				}
//			}
//			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//				//Found neighbour of same size
//				isghost.push_back(false);
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(octants[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > noctants-1){
//							idxtry = noctants;
//							break;
//						}
//					}
//					while(octants[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > noctants-1){
//							idxtry = 0;
//							break;
//						}
//					}
//				}
//				if (idxtry < noctants){
//					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//						//Found neighbour of same size
//						isghost.push_back(false);
//						neighbours.push_back(idxtry);
//						return;
//					}
//					// Compute Last discendent of virtual octant of same size
//					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//					Class_Octant last_desc = samesizeoct.buildLastDesc();
//					uint64_t Mortonlast = last_desc.computeMorton();
//					vector<uint32_t> bufferidx;
//					Mortontry = octants[idxtry].computeMorton();
//					while(Mortontry < Mortonlast & idxtry < noctants-1){
//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
//						Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
//						Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
//						Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
//						if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
//							neighbours.push_back(idxtry);
//							isghost.push_back(false);
//						}
//						idxtry++;
//						if(idxtry>noctants-1){
//							break;
//						}
//						Mortontry = octants[idxtry].computeMorton();
//					}
//				}
//			}
//			return;
//		}
//	}
//	else{
//		// Boundary Face
//		return;
//	}
//}
//
//// =================================================================================== //
//
//void Class_Local_Tree::findNodeNeighbours(uint32_t const idx, uint8_t inode,
//		u32vector & neighbours, vector<bool> & isghost) {
//
//	uint64_t  Morton, Mortontry;
//	uint32_t  noctants = getNumOctants();
//	uint32_t idxtry;
//	Class_Octant* oct = &octants[idx];
//	uint32_t size = oct->getSize();
//	uint8_t iface1, iface2, iface3;
//	int32_t Dhx, Dhy, Dhz;
//	int32_t Dhxref, Dhyref, Dhzref;
//
//	//Alternative to switch case
//	int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
//	int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
//	int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
//	int8_t cx = Cx[inode];
//	int8_t cy = Cy[inode];
//	int8_t cz = Cz[inode];
//
//	isghost.clear();
//	neighbours.clear();
//
//	// Default if inode is nnodes<inode<0
//	if (inode < 0 || inode > nnodes){
//		writeLog("Node index out of range in find neighbours !!!");
//		return;
//	}
//
//	// Check if octants node is a boundary
//	iface1 = nodeface[inode][0];
//	iface2 = nodeface[inode][1];
//	iface3 = nodeface[inode][2];
//
//	// Check if octants node is a boundary
//	if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){
//
//		//Build Morton number of virtual neigh of same size
//		Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
//
//		//SEARCH IN GHOSTS
//
//		if (ghosts.size()>0){
//			// Search in ghosts
//
//			uint32_t idxghost = uint32_t(size_ghosts/2);
//			Class_Octant* octghost = &ghosts[idxghost];
//
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
//			idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
//			if (idxtry > size_ghosts-1)
//				idxtry = size_ghosts-1;
//			while(abs(jump) > 0){
//				Mortontry = ghosts[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//				if (idxtry > ghosts.size()-1){
//					if (jump > 0){
//						idxtry = ghosts.size() - 1;
//						jump = 0;
//					}
//					else if (jump < 0){
//						idxtry = 0;
//						jump = 0;
//					}
//				}
//			}
//			if(octants[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//				//Found neighbour of same size
//				isghost.push_back(true);
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(ghosts[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > ghosts.size()-1){
//							idxtry = ghosts.size();
//							break;
//						}
//					}
//					while(ghosts[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > ghosts.size()-1){
//							idxtry = ghosts.size();
//							break;
//						}
//					}
//				}
//				if(idxtry < size_ghosts){
//					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
//						//Found neighbour of same size
//						isghost.push_back(true);
//						neighbours.push_back(idxtry);
//						return;
//					}
//					// Compute Last discendent of virtual octant of same size
//					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//					Class_Octant last_desc = samesizeoct.buildLastDesc();
//					uint64_t Mortonlast = last_desc.computeMorton();
//					vector<uint32_t> bufferidx;
//					Mortontry = ghosts[idxtry].computeMorton();
//					while(Mortontry < Mortonlast & idxtry < size_ghosts){
//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
//						Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
//						Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
//						Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
//						if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
//							neighbours.push_back(idxtry);
//							isghost.push_back(true);
//						}
//						idxtry++;
//						if(idxtry>size_ghosts-1){
//							break;
//						}
//						Mortontry = ghosts[idxtry].computeMorton();
//					}
//				}
//			}
//		}
//		uint32_t lengthneigh = 0;
//		uint32_t sizeneigh = neighbours.size();
//		for (idxtry=0; idxtry<sizeneigh; idxtry++){
//			lengthneigh += ghosts[neighbours[idxtry]].getSize();
//		}
//		if (lengthneigh < oct->getSize()){
//
//			// Search in octants
//
//			//Build Morton number of virtual neigh of same size
//			//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
//			//Morton = samesizeoct.computeMorton();
//			// Search morton in octants
//			// If a even face morton is lower than morton of oct, if odd higher
//			// ---> can i search only before or after idx in octants
//			int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
//			idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
//			if (idxtry > noctants-1)
//				idxtry = noctants-1;
//			while(abs(jump) > 0){
//				Mortontry = octants[idxtry].computeMorton();
//				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
//				idxtry += jump;
//			}
//			if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//				//Found neighbour of same size
//				isghost.push_back(false);
//				neighbours.push_back(idxtry);
//				return;
//			}
//			else{
//				// Step until the mortontry lower than morton (one idx of distance)
//				{
//					while(octants[idxtry].computeMorton() < Morton){
//						idxtry++;
//						if(idxtry > noctants-1){
//							idxtry = noctants;
//							break;
//						}
//					}
//					while(octants[idxtry].computeMorton() > Morton){
//						idxtry--;
//						if(idxtry > noctants-1){
//							idxtry = noctants;
//							break;
//						}
//					}
//				}
//				if (idxtry < noctants){
//					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
//						//Found neighbour of same size
//						isghost.push_back(false);
//						neighbours.push_back(idxtry);
//						return;
//					}
//					// Compute Last discendent of virtual octant of same size
//					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - samesizeoct.level)) - 1;
//					Class_Octant last_desc = samesizeoct.buildLastDesc();
//					uint64_t Mortonlast = last_desc.computeMorton();
//					vector<uint32_t> bufferidx;
//					Mortontry = octants[idxtry].computeMorton();
//					while(Mortontry < Mortonlast & idxtry < noctants-1){
//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
//						Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
//						Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
//						Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
//						if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
//							neighbours.push_back(idxtry);
//							isghost.push_back(false);
//						}
//						idxtry++;
//						Mortontry = octants[idxtry].computeMorton();
//					}
//				}
//			}
//			return;
//		}
//	}
//	else{
//		// Boundary Node
//		return;
//	}
//}
//
//// =================================================================================== //
