// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classLocalTree.hpp"

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*!Dimensional and default constructor.
 * \param[in] dim_ Space dimension of octree.
 */
classLocalTree::classLocalTree(int8_t maxlevel, uint8_t dim_){
	dim = dim_;
	global.setGlobal(maxlevel, dim);
	classOctant oct0(dim);
	classOctant octf(dim,global.MAX_LEVEL,0,0,0);
	classOctant octl(dim,global.MAX_LEVEL,global.max_length-1,global.max_length-1,(dim-2)*(global.max_length-1));
	octants.resize(1);
	octants[0] = oct0;
	first_desc = octf;
	last_desc = octl;
	size_ghosts = 0;
	local_max_depth = 0;
	balance_codim = 1;
};

/*!Default destructor.
 */
classLocalTree::~classLocalTree(){};

// =================================================================================== //
// METHODS
// =================================================================================== //

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //

const classOctant &  classLocalTree::getFirstDesc() const{
	return first_desc;
};
const classOctant &  classLocalTree::getLastDesc() const{
	return last_desc;
};
uint32_t classLocalTree::getSizeGhost() const{
	return size_ghosts;
};
uint32_t classLocalTree::getNumOctants() const{
	return octants.size();
};
uint8_t classLocalTree::getLocalMaxDepth() const{							// Get max depth reached in local tree
	return local_max_depth;
};
int8_t classLocalTree::getMarker(int32_t idx){								// Get refinement/coarsening marker for idx-th octant
	return octants[idx].getMarker();
};
uint8_t classLocalTree::getLevel(int32_t idx){								// Get refinement/coarsening marker for idx-th octant
	return octants[idx].getLevel();
};
uint8_t classLocalTree::getGhostLevel(int32_t idx){								// Get refinement/coarsening marker for idx-th ghost octant
	return ghosts[idx].getLevel();
};
bool classLocalTree::getBalance(int32_t idx){								// Get if balancing-blocked idx-th octant
	return octants[idx].getNotBalance();
};

/*! Get the codimension for 2:1 balancing
 * \return Maximum codimension of the entity through which the 2:1 balance is performed.
 */
uint8_t classLocalTree::getBalanceCodim() const{
	return balance_codim;
};

void classLocalTree::setMarker(int32_t idx, int8_t marker){					// Set refinement/coarsening marker for idx-th octant
	octants[idx].setMarker(marker);
};
void classLocalTree::setBalance(int32_t idx, bool balance){					// Set if balancing-blocked idx-th octant
	octants[idx].setBalance(balance);
};

/*! Set the codimension for 2:1 balancing
 * \param[in] Maximum codimension of the entity through which the 2:1 balance is performed.
 */
void classLocalTree::setBalanceCodim(uint8_t b21codim){
	balance_codim = b21codim;
};

void classLocalTree::setFirstDesc(){
	octvector::const_iterator firstOctant = octants.begin();
	first_desc = classOctant(dim, global.MAX_LEVEL, firstOctant->x, firstOctant->y, firstOctant->z);
};

void classLocalTree::setLastDesc(){
	octvector::const_iterator lastOctant = octants.end() - 1;
	uint32_t x,y,z,delta;
	//delta = (uint32_t)pow(2.0,(double)((uint8_t)global.MAX_LEVEL - lastOctant->level)) - 1;
	delta = (uint32_t)(1<<((uint8_t)global.MAX_LEVEL - lastOctant->level)) - 1;
	x = lastOctant->x + delta;
	y = lastOctant->y + delta;
	z = lastOctant->z + (dim-2)*delta;
	last_desc = classOctant(dim, global.MAX_LEVEL,x,y,z);
};


// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

// =================================================================================== //
// OTHER METHODS
// =================================================================================== //

classOctant& classLocalTree::extractOctant(uint32_t idx){
	return octants[idx];
};

const classOctant&	classLocalTree::extractOctant(uint32_t idx) const{
	return octants[idx];
};

classOctant& classLocalTree::extractGhostOctant(uint32_t idx) {
	return ghosts[idx];
};

const classOctant& classLocalTree::extractGhostOctant(uint32_t idx) const{
	return ghosts[idx];
};

// =================================================================================== //

/*! Refine local tree: refine one time octants with marker >0
 * \param[in] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
 * \return	true if refinement done
 */
bool classLocalTree::refine(u32vector* mapidx){

	u32vector		last_child_index;
	octvector 		children;
	uint32_t 		idx, nocts, ilastch;
	uint32_t 		offset = 0, blockidx;
	uint8_t 		nchm1 = global.nchildren-1, ich;
	bool 			dorefine = false;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < global.MAX_LEVEL){
			last_child_index.push_back(idx+nchm1+offset);
			offset += nchm1;
		}
		else{
			if (octants[idx].marker > 0){
				octants[idx].marker = 0;
				octants[idx].info[15] = false;
			}
		}
	}
	if (offset > 0){
		if(mapidx != NULL){
			mapidx->resize(octants.size()+offset);
			u32vector((*mapidx)).swap((*mapidx));
		}
		octants.resize(octants.size()+offset);
		blockidx = last_child_index[0]-nchm1;
		idx = octants.size();
		ilastch = last_child_index.size()-1;
		while (idx>blockidx){
			idx--;
			if(idx == last_child_index[ilastch]){
				children = octants[idx-offset].buildChildren(global.MAX_LEVEL);
				for (ich=0; ich<global.nchildren; ich++){
					octants[idx-ich] = (children[nchm1-ich]);
					if(mapidx != NULL) mapidx[idx-ich]  = mapidx[idx-offset];
				}
				offset -= nchm1;
				idx -= nchm1;
				//Update local max depth
				if (children[0].getLevel() > local_max_depth){
					local_max_depth = children[0].getLevel();
				}
				if (children[0].getMarker() > 0){
					//More Refinement to do
					dorefine = true;
				}
				//delete []children;
				if (ilastch != 0){
					ilastch--;
				}
			}
			else {
				octants[idx] = octants[idx-offset];
				if(mapidx != NULL) mapidx[idx]  = mapidx[idx-offset];
			}
		}
	}
	octvector(octants).swap(octants);
	nocts = octants.size();
	if(mapidx != NULL) {
		mapidx->resize(nocts);
		u32vector((*mapidx)).swap((*mapidx));
	}

	setFirstDesc();
	setLastDesc();

	return dorefine;

};

// =================================================================================== //
/*! Coarse local tree: coarse one time family of octants with marker <0
 * (if at least one octant of family has marker>=0 set marker=0 for the entire family)
 * \param[in] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
 * \return	true is coarsening done
 */
bool classLocalTree::coarse(u32vector* mapidx){

	u32vector		first_child_index;
	classOctant		father;
	uint32_t 		nocts, nocts0;
	uint32_t 		idx, idx2;
	uint32_t 		offset;
	uint32_t 		idx2_gh;
	uint32_t 		nidx;
	int8_t 			markerfather, marker;
	uint8_t 		nbro, nend;
	uint8_t 		nchm1 = global.nchildren-1;
	bool 			docoarse = false;
	bool 			wstop = false;

	//------------------------------------------ //
	// Initialization

	nbro = nend = 0;
	nidx = offset = 0;

	idx2_gh = 0;

	nocts = nocts0 = octants.size();
	size_ghosts = ghosts.size();


	// Init first and last desc (even if already calculated)
	setFirstDesc();
	setLastDesc();

	//------------------------------------------ //

	// Set index for start and end check for ghosts
	if (ghosts.size()){
		while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() < last_desc.computeMorton()){
			idx2_gh++;
		}
		idx2_gh = min((size_ghosts-1), idx2_gh);
	}

	// Check and coarse internal octants
	for (idx=0; idx<nocts; idx++){
		if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
			nbro = 0;
			father = octants[idx].buildFather(global.MAX_LEVEL);
			// Check if family is to be refined
			for (idx2=idx; idx2<idx+global.nchildren; idx2++){
				if (idx2<nocts){
					if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather(global.MAX_LEVEL) == father){
						nbro++;
					}
				}
			}
			if (nbro == global.nchildren){
				nidx++;
				first_child_index.push_back(idx);
				idx = idx2-1;
			}
		}
	}
	uint32_t nblock = nocts;
	uint32_t nfchild = first_child_index.size();
	if (nidx!=0){
		nblock = nocts - nidx*nchm1;
		nidx = 0;
		for (idx=0; idx<nblock; idx++){
			if (nidx < nfchild){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -global.MAX_LEVEL;
					father = octants[idx+offset].buildFather(global.MAX_LEVEL);
					for (uint32_t iii=0; iii<17; iii++){
						father.info[iii] = false;
					}
					for(idx2=0; idx2<global.nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()+1){
							markerfather = octants[idx+offset+idx2].getMarker()+1;
						}
						for (uint32_t iii=0; iii<17; iii++){
							father.info[iii] = father.info[iii] || octants[idx+offset+idx2].info[iii];
						}
					}
					father.info[13] = true;
					father.info[15] = true;
					father.setMarker(markerfather);
					if (markerfather < 0){
						docoarse = true;
					}
					octants[idx] = father;
					if(mapidx != NULL) mapidx[idx] = mapidx[idx+offset];
					offset += nchm1;
					nidx++;
				}
				else{
					octants[idx] = octants[idx+offset];
					if(mapidx != NULL) mapidx[idx] = mapidx[idx+offset];
				}
			}
			else{
				octants[idx] = octants[idx+offset];
				if(mapidx != NULL) mapidx[idx] = mapidx[idx+offset];
			}
		}
	}
	octants.resize(nblock);
	octvector(octants).swap(octants);
	nocts = octants.size();
	if(mapidx != NULL){
		mapidx->resize(nocts);
		u32vector((*mapidx)).swap((*mapidx));
	}

	// End on ghosts
	if (ghosts.size() && nocts > 0){
		if (ghosts[idx2_gh].buildFather(global.MAX_LEVEL) == octants[nocts-1].buildFather(global.MAX_LEVEL)){
			father = ghosts[idx2_gh].buildFather(global.MAX_LEVEL);
			for (uint32_t iii=0; iii<17; iii++){
				father.info[iii] = false;
			}
			markerfather = ghosts[idx2_gh].getMarker()+1;
			nbro = 0;
			idx = idx2_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 && ghosts[idx].buildFather(global.MAX_LEVEL) == father){
				nbro++;
				if (markerfather < ghosts[idx].getMarker()+1){
					markerfather = ghosts[idx].getMarker()+1;
				}
				for (uint32_t iii=0; iii<global.nfaces; iii++){
					father.info[iii] = father.info[iii] || ghosts[idx].info[iii];
				}
				father.info[14] = father.info[14] || ghosts[idx].info[14];
				idx++;
				if(idx == size_ghosts){
					break;
				}
				marker = ghosts[idx].getMarker();
			}
			nend = 0;
			idx = nocts-1;
			marker = octants[idx].getMarker();
			while(marker < 0 && octants[idx].buildFather(global.MAX_LEVEL) == father && idx >= 0){
				nbro++;
				nend++;
				if (markerfather < octants[idx].getMarker()+1){
					markerfather = octants[idx].getMarker()+1;
				}
				idx--;
				marker = octants[idx].getMarker();
				if (wstop){
					break;
				}
				if (idx==0){
					wstop = true;
				}
			}
			if (nbro == global.nchildren){
				offset = nend;
			}
			else{
				nend = 0;
			}
		}
		if (nend != 0){
			for (idx=0; idx < nend; idx++){
				for (uint32_t iii=0; iii<16; iii++){
					father.info[iii] = father.info[iii] || octants[nocts-idx-1].info[iii];
				}
			}
			father.info[13] = true;
			father.info[15] = true;
			if (markerfather < 0){
				docoarse = true;
			}
			father.setMarker(markerfather);
			octants.resize(nocts-offset);
			octants.push_back(father);
			octvector(octants).swap(octants);
			nocts = octants.size();
			if(mapidx != NULL){
				mapidx->resize(nocts);
				u32vector((*mapidx)).swap((*mapidx));
			}
		}

	}

	// Set final first and last desc
	if(nocts>0){
		setFirstDesc();
		setLastDesc();
	}
	return docoarse;

};

// =================================================================================== //

/*! Refine local tree: refine one time all the octants
 * \param[in] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
 * \return	true if refinement done
 */
bool classLocalTree::globalRefine(u32vector* mapidx){

	uint32_t 	idx, nocts;
	bool 		dorefine = false;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		octants[idx].setMarker(1);
	}

	dorefine = refine(mapidx);

	return dorefine;

};

// =================================================================================== //

/*! Refine local tree: corse one time all the octants
 * \param[in] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
 * \return	true if refinement done
 */
bool classLocalTree::globalCoarse(u32vector* mapidx){

	uint32_t 	idx, nocts;
	bool 		dorefine = false;

	nocts = octants.size();
	for (idx=0; idx<nocts; idx++){
		octants[idx].setMarker(-1);
	}

	dorefine = coarse(mapidx);

	return dorefine;

};

// =================================================================================== //
/*! Delete overlapping octants after coarse local tree. Check first and last descendants
 * of process before and after the local process
 */
void classLocalTree::checkCoarse(uint64_t lastDescPre,
		uint64_t firstDescPost,
		u32vector* mapidx){

	uint32_t		idx;
	uint32_t 		nocts;
	uint64_t 		Morton;
	uint8_t 		toDelete = 0;

	nocts = getNumOctants();
	idx = 0;
	Morton = octants[idx].computeMorton();
	while(Morton <= lastDescPre && idx < nocts && Morton != 0){
		// To delete, the father is in proc before me
		toDelete++;
		idx++;
		Morton = octants[idx].computeMorton();
	}
	for(idx=0; idx<nocts-toDelete; idx++){
		octants[idx] = octants[idx+toDelete];
		if (mapidx != NULL) mapidx[idx] = mapidx[idx+toDelete];
	}
	octants.resize(nocts-toDelete);
	octvector(octants).swap(octants);
	if (mapidx != NULL){
		mapidx->resize(nocts-toDelete);
		u32vector((*mapidx)).swap((*mapidx));
	}
	nocts = getNumOctants();

	setFirstDesc();
	setLastDesc();

};

// =================================================================================== //
/*! Update max depth reached in local tree
 */
void classLocalTree::updateLocalMaxDepth(){

	uint32_t noctants = getNumOctants();
	uint32_t i;

	local_max_depth = 0;
	for(i = 0; i < noctants; i++){
		if(octants[i].getLevel() > local_max_depth){
			local_max_depth = octants[i].getLevel();
		}
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iface in vector octants.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findNeighbours(uint32_t idx, uint8_t iface,
		u32vector & neighbours, vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t 		noctants = getNumOctants();
	uint32_t 		idxtry;
	classOctant* 	oct = &octants[idx];
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);

	//	int8_t 			cx = global.normals[iface][0];
	//	int8_t 			cy = global.normals[iface][1];
	//	int8_t 			cz = global.normals[iface][2];
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.normals[iface][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > global.nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->info[6+iface] == false){

		// Check if octants face is a boundary
		if (oct->info[iface] == false){

			//Build Morton number of virtual neigh of same size
			classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = oct->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
			idxtry = uint32_t(idx + ((Mortontry < Morton) - (Mortontry > Morton))*jump);
			if (idxtry > noctants-1) idxtry = noctants-1;
			while(abs(jump) > 0){
				Mortontry = octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > noctants-1){
					if (jump > 0){
						idxtry = noctants - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			Mortontry = octants[idxtry].computeMorton();
			if(Mortontry == Morton && octants[idxtry].level == oct->level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(Mortontry < Morton){
						idxtry++;
						if(idxtry > noctants-1){
							idxtry = noctants-1;
							Mortontry = octants[idxtry].computeMorton();
							break;
						}
						Mortontry = octants[idxtry].computeMorton();
					}
					while(Mortontry > Morton){
						idxtry--;
						if(idxtry > noctants-1){
							idxtry = 0;
							Mortontry = octants[idxtry].computeMorton();
							break;
						}
						Mortontry = octants[idxtry].computeMorton();
					}
				}
				if(Mortontry == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				//				int32_t Dx, Dy, Dz;
				//				int32_t Dxstar, Dystar, Dzstar;
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = octants[idxtry].getCoord();
				u32array3 coord1 = {1,1,1};
				u32array3 coordtry1 = {1,1,1};
				uint8_t level = oct->level;
				uint8_t leveltry = octants[idxtry].getLevel();
				while(Mortontry < Mortonlast && idxtry < noctants){
					//					Dx = int32_t(abs(cxyz[0]))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
					//					Dy = int32_t(abs(cxyz[1]))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
					//					Dz = int32_t(abs(cxyz[2]))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
					//					Dxstar = int32_t((cxyz[0]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[0]+1)/2)*size;
					//					Dystar = int32_t((cxyz[1]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[1]+1)/2)*size;
					//					Dzstar = int32_t((cxyz[2]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[2]+1)/2)*size;
					for (int idim=0; idim<dim; idim++){
						Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
					}

					//					uint32_t x0 = oct->x;
					//					uint32_t x1 = x0 + size;
					//					uint32_t y0 = oct->y;
					//					uint32_t y1 = y0 + size;
					//					uint32_t z0 = oct->z;
					//					uint32_t z1 = z0 + size;
					//					uint32_t x0try = octants[idxtry].x;
					//					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
					//					uint32_t y0try = octants[idxtry].y;
					//					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
					//					uint32_t z0try = octants[idxtry].z;
					//					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
					leveltry = octants[idxtry].getLevel();


					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
						if (leveltry > level){
							//							if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
							if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							//							if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
							if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
					coordtry = octants[idxtry].getCoord();
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
			if (ghosts.size()>0){

				// Search in ghosts
				uint32_t idxghost = uint32_t(size_ghosts/2);
				classOctant* octghost = &ghosts[idxghost];

				//Build Morton number of virtual neigh of same size
				//classOctant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				Mortontry = octghost->computeMorton();
				int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
				if (idxtry > ghosts.size()-1) idxtry = ghosts.size()-1;
				while(abs(jump) > 0){
					Mortontry = ghosts[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
					idxtry += jump;
					if (idxtry > ghosts.size()-1){
						if (jump > 0){
							idxtry = ghosts.size() - 1;
							Mortontry = ghosts[idxtry].computeMorton();
							jump = 0;
						}
						else if (jump < 0){
							idxtry = 0;
							Mortontry = ghosts[idxtry].computeMorton();
							jump = 0;
						}
					}
				}
				Mortontry = ghosts[idxtry].computeMorton();
				if(Mortontry == Morton && ghosts[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				else{
					// Step until the mortontry lower than morton (one idx of distance)
					{
						while(Mortontry < Morton){
							idxtry++;
							if(idxtry > ghosts.size()-1){
								idxtry = ghosts.size()-1;
								Mortontry = ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
						}
						while(ghosts[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > ghosts.size()-1){
								idxtry = 0;
								Mortontry = ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
						}
					}
					if(idxtry < size_ghosts){
						if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
							//Found neighbour of same size
							isghost.push_back(true);
							neighbours.push_back(idxtry);
							return;
						}
						// Compute Last discendent of virtual octant of same size
						classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
						uint64_t Mortonlast = last_desc.computeMorton();
						Mortontry = ghosts[idxtry].computeMorton();
						//						int32_t Dx, Dy, Dz;
						//						int32_t Dxstar, Dystar, Dzstar;
						int32_t Dx[3] = {0,0,0};
						int32_t Dxstar[3] = {0,0,0};
						u32array3 coord = oct->getCoord();
						u32array3 coordtry = ghosts[idxtry].getCoord();
						u32array3 coord1 = {1,1,1};
						u32array3 coordtry1 = {1,1,1};
						uint8_t level = oct->level;
						uint8_t leveltry = octants[idxtry].getLevel();
						while(Mortontry < Mortonlast && idxtry < size_ghosts){
							//							Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(ghosts[idxtry].x));
							//							Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(ghosts[idxtry].y));
							//							Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(ghosts[idxtry].z));
							//							Dxstar = int32_t((cx-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
							//							Dystar = int32_t((cy-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
							//							Dzstar = int32_t((cz-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;
							for (int idim=0; idim<dim; idim++){
								Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
								Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
								coord1[idim] 	= coord[idim] + size;
								coordtry1[idim] = coordtry[idim] + ghosts[idxtry].getSize(global.MAX_LEVEL);
							}

							//							uint32_t x0 = oct->x;
							//							uint32_t x1 = x0 + size;
							//							uint32_t y0 = oct->y;
							//							uint32_t y1 = y0 + size;
							//							uint32_t z0 = oct->z;
							//							uint32_t z1 = z0 + size;
							//							uint32_t x0try = ghosts[idxtry].x;
							//							uint32_t x1try = x0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint32_t y0try = ghosts[idxtry].y;
							//							uint32_t y1try = y0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint32_t z0try = ghosts[idxtry].z;
							//							uint32_t z1try = z0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint8_t level = oct->level;
							leveltry = ghosts[idxtry].getLevel();

							if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
								if (leveltry > level){
									//									if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
									if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
								if (leveltry < level){
									//									if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
									if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
							}

							idxtry++;
							if(idxtry>size_ghosts-1){
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
							coordtry = ghosts[idxtry].getCoord();
						}
					}
				}

				uint32_t lengthneigh = 0;
				uint32_t sizeneigh = neighbours.size();
				for (idxtry=0; idxtry<sizeneigh; idxtry++){
					lengthneigh += ghosts[neighbours[idxtry]].getArea(global.MAX_LEVEL);
				}
				if (lengthneigh < oct->getArea(global.MAX_LEVEL)){
					// Search in octants

					// Check if octants face is a boundary
					if (oct->info[iface] == false){

						//Build Morton number of virtual neigh of same size
						classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
						Morton = samesizeoct.computeMorton();
						// Search morton in octants
						// If a even face morton is lower than morton of oct, if odd higher
						// ---> can i search only before or after idx in octants
						Mortontry = oct->computeMorton();
						int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
						idxtry = uint32_t(idx + ((Mortontry < Morton) - (Mortontry > Morton))*jump);
						if (idxtry > noctants-1) idxtry = noctants-1;
						while(abs(jump) > 0){
							Mortontry = octants[idxtry].computeMorton();
							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
							idxtry += jump;
							if (idxtry > noctants-1){
								if (jump > 0){
									idxtry = noctants - 1;
									jump = 0;
								}
								else if (jump < 0){
									idxtry = 0;
									jump = 0;
								}
							}
						}
						Mortontry = octants[idxtry].computeMorton();
						if(Mortontry == Morton && octants[idxtry].level == oct->level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							return;
						}
						else{
							// Step until the mortontry lower than morton (one idx of distance)
							{
								while(Mortontry < Morton){
									idxtry++;
									if(idxtry > noctants-1){
										idxtry = noctants-1;
										Mortontry = octants[idxtry].computeMorton();
										break;
									}
									Mortontry = octants[idxtry].computeMorton();
								}
								while(Mortontry > Morton){
									idxtry--;
									if(idxtry > noctants-1){
										idxtry = 0;
										Mortontry = octants[idxtry].computeMorton();
										break;
									}
									Mortontry = octants[idxtry].computeMorton();
								}
							}
							if(Mortontry == Morton && octants[idxtry].level == oct->level){
								//Found neighbour of same size
								isghost.push_back(false);
								neighbours.push_back(idxtry);
								return;
							}
							// Compute Last discendent of virtual octant of same size
							classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
							uint64_t Mortonlast = last_desc.computeMorton();
							Mortontry = octants[idxtry].computeMorton();
							//				int32_t Dx, Dy, Dz;
							//				int32_t Dxstar, Dystar, Dzstar;
							int32_t Dx[3] = {0,0,0};
							int32_t Dxstar[3] = {0,0,0};
							u32array3 coord = oct->getCoord();
							u32array3 coordtry = octants[idxtry].getCoord();
							u32array3 coord1 = {1,1,1};
							u32array3 coordtry1 = {1,1,1};
							uint8_t level = oct->level;
							uint8_t leveltry = octants[idxtry].getLevel();
							while(Mortontry < Mortonlast && idxtry < noctants){
								//					Dx = int32_t(abs(cxyz[0]))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
								//					Dy = int32_t(abs(cxyz[1]))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
								//					Dz = int32_t(abs(cxyz[2]))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
								//					Dxstar = int32_t((cxyz[0]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[0]+1)/2)*size;
								//					Dystar = int32_t((cxyz[1]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[1]+1)/2)*size;
								//					Dzstar = int32_t((cxyz[2]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[2]+1)/2)*size;
								for (int idim=0; idim<dim; idim++){
									Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
									Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
									coord1[idim] 	= coord[idim] + size;
									coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
								}

								//					uint32_t x0 = oct->x;
								//					uint32_t x1 = x0 + size;
								//					uint32_t y0 = oct->y;
								//					uint32_t y1 = y0 + size;
								//					uint32_t z0 = oct->z;
								//					uint32_t z1 = z0 + size;
								//					uint32_t x0try = octants[idxtry].x;
								//					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
								//					uint32_t y0try = octants[idxtry].y;
								//					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
								//					uint32_t z0try = octants[idxtry].z;
								//					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
								leveltry = octants[idxtry].getLevel();


								if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
									if (leveltry > level){
										//							if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
										if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
									if (leveltry < level){
										//							if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
										if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
								}

								idxtry++;
								if(idxtry>noctants-1){
									break;
								}
								Mortontry = octants[idxtry].computeMorton();
								coordtry = octants[idxtry].getCoord();
							}
							return;
						}
					}
				}
				return;
			}
		}
		else{
			// Boundary Face
			return;
		}
	}
};

// =================================================================================== //

/*! Finds neighbours of octant through iface in vector octants.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findNeighbours(classOctant* oct,
		uint8_t iface,
		u32vector & neighbours,
		vector<bool> & isghost){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	uint32_t size = oct->getSize(global.MAX_LEVEL);

	//	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
	//	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
	//	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.normals[iface][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > global.nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->info[6+iface] == false){

		// Check if octants face is a boundary
		if (oct->info[iface] == false){

			//Build Morton number of virtual neigh of same size
			//			classOctant samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
			classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = int32_t((noctants)/2+1);
			idxtry = uint32_t(jump);
			Mortontry = oct->computeMorton();
			while(abs(jump) > 0){
				Mortontry = octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > noctants-1){
					if (jump > 0){
						idxtry = noctants - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			Mortontry = octants[idxtry].computeMorton();
			if(Mortontry == Morton && octants[idxtry].level == oct->level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(Mortontry < Morton){
						idxtry++;
						if(idxtry > noctants-1){
							idxtry = noctants-1;
							Mortontry = octants[idxtry].computeMorton();
							break;
						}
						Mortontry = octants[idxtry].computeMorton();
					}
					while(Mortontry > Morton){
						idxtry--;
						if(idxtry > noctants-1){
							idxtry = 0;
							Mortontry = octants[idxtry].computeMorton();
							break;
						}
						Mortontry = octants[idxtry].computeMorton();
					}
				}
				if(Mortontry == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				//				int32_t Dx, Dy, Dz;
				//				int32_t Dxstar, Dystar, Dzstar;
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = octants[idxtry].getCoord();
				u32array3 coord1 = {1,1,1};
				u32array3 coordtry1 = {1,1,1};
				uint8_t level = oct->level;
				uint8_t leveltry = octants[idxtry].getLevel();
				while(Mortontry < Mortonlast && idxtry < noctants){
					//					Dx = int32_t(abs(cxyz[0]))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
					//					Dy = int32_t(abs(cxyz[1]))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
					//					Dz = int32_t(abs(cxyz[2]))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
					//					Dxstar = int32_t((cxyz[0]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[0]+1)/2)*size;
					//					Dystar = int32_t((cxyz[1]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[1]+1)/2)*size;
					//					Dzstar = int32_t((cxyz[2]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[2]+1)/2)*size;
					for (int idim=0; idim<dim; idim++){
						Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
					}

					//					uint32_t x0 = oct->x;
					//					uint32_t x1 = x0 + size;
					//					uint32_t y0 = oct->y;
					//					uint32_t y1 = y0 + size;
					//					uint32_t z0 = oct->z;
					//					uint32_t z1 = z0 + size;
					//					uint32_t x0try = octants[idxtry].x;
					//					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
					//					uint32_t y0try = octants[idxtry].y;
					//					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
					//					uint32_t z0try = octants[idxtry].z;
					//					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
					leveltry = octants[idxtry].getLevel();


					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
						if (leveltry > level){
							//							if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
							if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							//							if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
							if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
					coordtry = octants[idxtry].getCoord();
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
			if (ghosts.size()>0){

				// Search in ghosts
				uint32_t idxghost = uint32_t(size_ghosts/2);
				classOctant* octghost = &ghosts[idxghost];

				//Build Morton number of virtual neigh of same size
				//classOctant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				Mortontry = octghost->computeMorton();
				int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
				if (idxtry > ghosts.size()-1) idxtry = ghosts.size()-1;
				while(abs(jump) > 0){
					Mortontry = ghosts[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
					idxtry += jump;
					if (idxtry > ghosts.size()-1){
						if (jump > 0){
							idxtry = ghosts.size() - 1;
							Mortontry = ghosts[idxtry].computeMorton();
							jump = 0;
						}
						else if (jump < 0){
							idxtry = 0;
							Mortontry = ghosts[idxtry].computeMorton();
							jump = 0;
						}
					}
				}
				Mortontry = ghosts[idxtry].computeMorton();
				if(Mortontry == Morton && ghosts[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				else{
					// Step until the mortontry lower than morton (one idx of distance)
					{
						while(Mortontry < Morton){
							idxtry++;
							if(idxtry > ghosts.size()-1){
								idxtry = ghosts.size()-1;
								Mortontry = ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
						}
						while(ghosts[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > ghosts.size()-1){
								idxtry = 0;
								Mortontry = ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
						}
					}
					if(idxtry < size_ghosts){
						if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
							//Found neighbour of same size
							isghost.push_back(true);
							neighbours.push_back(idxtry);
							return;
						}
						// Compute Last discendent of virtual octant of same size
						classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
						uint64_t Mortonlast = last_desc.computeMorton();
						Mortontry = ghosts[idxtry].computeMorton();
						//						int32_t Dx, Dy, Dz;
						//						int32_t Dxstar, Dystar, Dzstar;
						int32_t Dx[3] = {0,0,0};
						int32_t Dxstar[3] = {0,0,0};
						u32array3 coord = oct->getCoord();
						u32array3 coordtry = ghosts[idxtry].getCoord();
						u32array3 coord1 = {1,1,1};
						u32array3 coordtry1 = {1,1,1};
						uint8_t level = oct->level;
						uint8_t leveltry = octants[idxtry].getLevel();
						while(Mortontry < Mortonlast && idxtry < size_ghosts){
							//							Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(ghosts[idxtry].x));
							//							Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(ghosts[idxtry].y));
							//							Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(ghosts[idxtry].z));
							//							Dxstar = int32_t((cx-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
							//							Dystar = int32_t((cy-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
							//							Dzstar = int32_t((cz-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;
							for (int idim=0; idim<dim; idim++){
								Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
								Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
								coord1[idim] 	= coord[idim] + size;
								coordtry1[idim] = coordtry[idim] + ghosts[idxtry].getSize(global.MAX_LEVEL);
							}

							//							uint32_t x0 = oct->x;
							//							uint32_t x1 = x0 + size;
							//							uint32_t y0 = oct->y;
							//							uint32_t y1 = y0 + size;
							//							uint32_t z0 = oct->z;
							//							uint32_t z1 = z0 + size;
							//							uint32_t x0try = ghosts[idxtry].x;
							//							uint32_t x1try = x0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint32_t y0try = ghosts[idxtry].y;
							//							uint32_t y1try = y0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint32_t z0try = ghosts[idxtry].z;
							//							uint32_t z1try = z0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
							//							uint8_t level = oct->level;
							leveltry = ghosts[idxtry].getLevel();

							if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
								if (leveltry > level){
									//									if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
									if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
								if (leveltry < level){
									//									if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
									if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
							}

							idxtry++;
							if(idxtry>size_ghosts-1){
								break;
							}
							Mortontry = ghosts[idxtry].computeMorton();
							coordtry = ghosts[idxtry].getCoord();
						}
					}
				}

				uint32_t lengthneigh = 0;
				uint32_t sizeneigh = neighbours.size();
				for (idxtry=0; idxtry<sizeneigh; idxtry++){
					lengthneigh += ghosts[neighbours[idxtry]].getArea(global.MAX_LEVEL);
				}
				if (lengthneigh < oct->getArea(global.MAX_LEVEL)){
					// Search in octants

					// Check if octants face is a boundary
					if (oct->info[iface] == false){

						//Build Morton number of virtual neigh of same size
						//classOctant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
						classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
						Morton = samesizeoct.computeMorton();
						// Search morton in octants
						// If a even face morton is lower than morton of oct, if odd higher
						// ---> can i search only before or after idx in octants
						int32_t jump = (int32_t((noctants)/2+1));
						idxtry = uint32_t(jump);
						while(abs(jump) > 0){
							Mortontry = octants[idxtry].computeMorton();
							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
							idxtry += jump;
							if (idxtry > noctants-1){
								if (jump > 0){
									idxtry = noctants - 1;
									jump = 0;
								}
								else if (jump < 0){
									idxtry = 0;
									jump = 0;
								}
							}
						}
						Mortontry = octants[idxtry].computeMorton();
						if(Mortontry == Morton && octants[idxtry].level == oct->level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							return;
						}
						else{
							// Step until the mortontry lower than morton (one idx of distance)
							{
								while(Mortontry < Morton){
									idxtry++;
									if(idxtry > noctants-1){
										idxtry = noctants-1;
										Mortontry = octants[idxtry].computeMorton();
										break;
									}
									Mortontry = octants[idxtry].computeMorton();
								}
								while(Mortontry > Morton){
									idxtry--;
									if(idxtry > noctants-1){
										idxtry = 0;
										Mortontry = octants[idxtry].computeMorton();
										break;
									}
									Mortontry = octants[idxtry].computeMorton();
								}
							}
							if(Mortontry == Morton && octants[idxtry].level == oct->level){
								//Found neighbour of same size
								isghost.push_back(false);
								neighbours.push_back(idxtry);
								return;
							}
							// Compute Last discendent of virtual octant of same size
							classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
							uint64_t Mortonlast = last_desc.computeMorton();
							Mortontry = octants[idxtry].computeMorton();
							//				int32_t Dx, Dy, Dz;
							//				int32_t Dxstar, Dystar, Dzstar;
							int32_t Dx[3] = {0,0,0};
							int32_t Dxstar[3] = {0,0,0};
							u32array3 coord = oct->getCoord();
							u32array3 coordtry = octants[idxtry].getCoord();
							u32array3 coord1 = {1,1,1};
							u32array3 coordtry1 = {1,1,1};
							uint8_t level = oct->level;
							uint8_t leveltry = octants[idxtry].getLevel();
							while(Mortontry < Mortonlast && idxtry < noctants){
								//					Dx = int32_t(abs(cxyz[0]))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
								//					Dy = int32_t(abs(cxyz[1]))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
								//					Dz = int32_t(abs(cxyz[2]))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
								//					Dxstar = int32_t((cxyz[0]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[0]+1)/2)*size;
								//					Dystar = int32_t((cxyz[1]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[1]+1)/2)*size;
								//					Dzstar = int32_t((cxyz[2]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[2]+1)/2)*size;
								for (int idim=0; idim<dim; idim++){
									Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
									Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
									coord1[idim] 	= coord[idim] + size;
									coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
								}

								//					uint32_t x0 = oct->x;
								//					uint32_t x1 = x0 + size;
								//					uint32_t y0 = oct->y;
								//					uint32_t y1 = y0 + size;
								//					uint32_t z0 = oct->z;
								//					uint32_t z1 = z0 + size;
								//					uint32_t x0try = octants[idxtry].x;
								//					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
								//					uint32_t y0try = octants[idxtry].y;
								//					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
								//					uint32_t z0try = octants[idxtry].z;
								//					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
								leveltry = octants[idxtry].getLevel();


								if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
									if (leveltry > level){
										//							if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
										if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
									if (leveltry < level){
										//							if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
										if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
								}

								idxtry++;
								if(idxtry>noctants-1){
									break;
								}
								Mortontry = octants[idxtry].computeMorton();
								coordtry = octants[idxtry].getCoord();
							}
							return;
						}
					}
				}
				return;
			}
		}
		else{
			// Boundary Face
			return;
		}
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost octant through iface in vector octants.
 * Returns a vector (empty if iface is not the pbound face for ghost) with the index of neighbours
 * in the structure octants
 */
void classLocalTree::findGhostNeighbours(uint32_t const idx,
		uint8_t iface,
		u32vector & neighbours){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	classOctant* oct = &ghosts[idx];
	uint32_t size = oct->getSize(global.MAX_LEVEL);

	//	int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
	//	int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
	//	int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.normals[iface][idim];
	}

	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface < 0 || iface > global.nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->info[6+iface] == true){

		//Build Morton number of virtual neigh of same size
		//classOctant samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
		classOctant samesizeoct(dim, oct->level, int32_t(oct->x)+int32_t(cxyz[0]*size), int32_t(oct->y)+int32_t(cxyz[1]*size), int32_t(oct->z)+int32_t(cxyz[2]*size));
		Morton = samesizeoct.computeMorton();
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = getNumOctants()/2;
		idxtry = uint32_t(getNumOctants()/2);
		Mortontry = octants[idxtry].computeMorton();
		jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > noctants-1){
				if (jump > 0){
					idxtry = noctants - 1;
					Mortontry = octants[idxtry].computeMorton();
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					Mortontry = octants[idxtry].computeMorton();
					jump = 0;
				}
			}
		}
		Mortontry = octants[idxtry].computeMorton();
		if(Mortontry == Morton && octants[idxtry].level == oct->level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(Mortontry < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						Mortontry = octants[idxtry].computeMorton();
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
				while(Mortontry > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						Mortontry = octants[idxtry].computeMorton();
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
			if(Mortontry == Morton && octants[idxtry].level == oct->level){
				//Found neighbour of same size
				neighbours.push_back(idxtry);
				return;
			}
			// Compute Last discendent of virtual octant of same size
			classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
			uint64_t Mortonlast = last_desc.computeMorton();
			Mortontry = octants[idxtry].computeMorton();
			//			int32_t Dx, Dy, Dz;
			//			int32_t Dxstar, Dystar, Dzstar;
			int32_t Dx[3] = {0,0,0};
			int32_t Dxstar[3] = {0,0,0};
			u32array3 coord = oct->getCoord();
			u32array3 coordtry = octants[idxtry].getCoord();
			u32array3 coord1 = {1,1,1};
			u32array3 coordtry1 = {1,1,1};
			uint8_t level = oct->level;
			uint8_t leveltry = octants[idxtry].getLevel();
			while(Mortontry < Mortonlast && idxtry < noctants){
				//				Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
				//				Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
				//				Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
				//				Dxstar = int32_t((cx-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
				//				Dystar = int32_t((cy-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
				//				Dzstar = int32_t((cz-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;
				for (int idim=0; idim<dim; idim++){
					Dx[idim] 		= int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]);
					Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
					coord1[idim] 	= coord[idim] + size;
					coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
				}

				//				uint32_t x0 = oct->x;
				//				uint32_t x1 = x0 + size;
				//				uint32_t y0 = oct->y;
				//				uint32_t y1 = y0 + size;
				//				uint32_t z0 = oct->z;
				//				uint32_t z1 = z0 + size;
				//				uint32_t x0try = octants[idxtry].x;
				//				uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
				//				uint32_t y0try = octants[idxtry].y;
				//				uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
				//				uint32_t z0try = octants[idxtry].z;
				//				uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
				//				uint8_t level = oct->level;
				uint8_t leveltry = octants[idxtry].getLevel();

				if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
					if (leveltry > level){
						//							if((abs(cx)*((y0try>=y0)*(y0try<y1))*((z0try>=z0)*(z0try<z1))) + (abs(cy)*((x0try>=x0)*(x0try<x1))*((z0try>=z0)*(z0try<z1))) + (abs(cz)*((x0try>=x0)*(x0try<x1))*((y0try>=y0)*(y0try<y1)))){
						if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
							neighbours.push_back(idxtry);
						}
					}
					if (leveltry < level){
						//							if((abs(cx)*((y0>=y0try)*(y0<y1try))*((z0>=z0try)*(z0<z1try))) + (abs(cy)*((x0>=x0try)*(x0<x1try))*((z0>=z0try)*(z0<z1try))) + (abs(cz)*((x0>=x0try)*(x0<x1try))*((y0>=y0try)*(y0<y1try)))){
						if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
							neighbours.push_back(idxtry);
						}
					}
				}

				idxtry++;
				if(idxtry>noctants-1){
					break;
				}
				Mortontry = octants[idxtry].computeMorton();
				coordtry = octants[idxtry].getCoord();
			}
			return;
		}
	}
	//--------------------------------------------------------------- //
	//-----Not Pbound face------------- ----------------------------- //
	else{
		return;
	}

};

// =================================================================================== //

void classLocalTree::preBalance21(bool internal){

	classOctant 	father, lastdesc;
	uint64_t 		mortonld;
	uint32_t 		nocts;
	uint32_t 		idx, idx2, idx0, last_idx;
	uint32_t 		idx1_gh, idx2_gh;
	int8_t 			markerfather, marker;
	uint8_t 		nbro;
	uint8_t 		nchm1 = global.nchildren-1;
	bool 			Bdone = false;

	//------------------------------------------ //
	// Initialization

	nbro = 0;
	idx2_gh = idx0 = 0;
	idx1_gh=0;

	nocts   = octants.size();
	size_ghosts = ghosts.size();
	last_idx=nocts-1;

	//Clean index of ghost brothers in case of coarsening a broken family
	last_ghost_bros.clear();


	// Set index for start and end check for ghosts
	if (ghosts.size()){
		while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() <= last_desc.computeMorton()){
			idx2_gh++;
		}
		idx2_gh = min((size_ghosts-1), idx2_gh);

		while(idx1_gh < size_ghosts && ghosts[idx1_gh].computeMorton() <= octants[0].computeMorton()){
			idx1_gh++;
		}
		idx1_gh-=1;
		if (idx1_gh==-1) idx1_gh=0;
	}

	// End on ghosts
	if (ghosts.size() && nocts > 0){
		if (ghosts[idx1_gh].buildFather(global.MAX_LEVEL)==octants[0].buildFather(global.MAX_LEVEL)){
			father = ghosts[idx1_gh].buildFather(global.MAX_LEVEL);
			nbro = 0;
			idx = idx1_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 && ghosts[idx].buildFather(global.MAX_LEVEL) == father){
				nbro++;
				if (idx==0)
					break;
				idx--;
				marker = ghosts[idx].getMarker();
			}
			idx = 0;
			while(idx<nocts && octants[idx].buildFather(global.MAX_LEVEL) == father){
				if(octants[idx].getMarker()<0)
					nbro++;
				idx++;
				if(idx==nocts)
					break;
			}
			if (nbro != global.nchildren && idx!=nocts-1){
				for(uint32_t ii=0; ii<idx; ii++){
					if (octants[ii].getMarker()<0){
						octants[ii].setMarker(0);
						octants[ii].info[15]=true;
						Bdone=true;
					}
				}
			}
		}

		if (ghosts[idx2_gh].buildFather(global.MAX_LEVEL)==octants[nocts-1].buildFather(global.MAX_LEVEL)){
			father = ghosts[idx2_gh].buildFather(global.MAX_LEVEL);
			nbro = 0;
			idx = idx2_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 && ghosts[idx].buildFather(global.MAX_LEVEL) == father){

				//Add ghost index to structure for mapper in case of coarsening a broken family
				last_ghost_bros.push_back(idx);

				nbro++;
				idx++;
				if(idx == size_ghosts){
					break;
				}
				marker = ghosts[idx].getMarker();
			}
			idx = nocts-1;
			while(octants[idx].buildFather(global.MAX_LEVEL) == father ){
				if (octants[idx].getMarker()<0)
					nbro++;
				if (idx==0)
					break;
				idx--;
			}
			last_idx=idx;
			if (nbro != global.nchildren && idx!=nocts-1){
				for(uint32_t ii=idx+1; ii<nocts; ii++){
					if (octants[ii].getMarker()<0){
						octants[ii].setMarker(0);
						octants[ii].info[15]=true;
						Bdone=true;
					}
				}
				//Clean ghost index to structure for mapper in case of coarsening a broken family
				last_ghost_bros.clear();
			}
		}
	}

	// Check first internal octants
	if (internal){
		father = octants[0].buildFather(global.MAX_LEVEL);
		lastdesc = father.buildLastDesc(global.MAX_LEVEL);
		mortonld = lastdesc.computeMorton();
		nbro = 0;
		for (idx=0; idx<global.nchildren; idx++){
			// Check if family is complete or to be checked in the internal loop (some brother refined)
			if (octants[idx].computeMorton() <= mortonld){
				nbro++;
			}
		}
		if (nbro != global.nchildren)
			idx0 = nbro;

		// Check and coarse internal octants
		for (idx=idx0; idx<nocts; idx++){
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nbro = 0;
				father = octants[idx].buildFather(global.MAX_LEVEL);
				// Check if family is to be coarsened
				for (idx2=idx; idx2<idx+global.nchildren; idx2++){
					if (idx2<nocts){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather(global.MAX_LEVEL) == father){
							nbro++;
						}
					}
				}
				if (nbro == global.nchildren){
					idx = idx2-1;
				}
				else{
					if (idx<=last_idx){
						octants[idx].setMarker(0);
						octants[idx].info[15]=true;
						Bdone=true;
					}
				}
			}
		}
	}
};

// =================================================================================== //

void classLocalTree::preBalance21(u32vector& newmodified){

	classOctant 		father, lastdesc;
	uint64_t 			mortonld;
	uint32_t 			nocts;
	uint32_t 			idx, idx2, idx0, last_idx;
	uint32_t 			idx1_gh, idx2_gh;
	int8_t 				markerfather, marker;
	uint8_t 			nbro;
	uint8_t 			nchm1 = global.nchildren-1;
	bool 				Bdone = false;

	//------------------------------------------ //
	// Initialization

	nbro = 0;
	idx2_gh = idx0 = 0;
	idx1_gh=0;

	nocts   = octants.size();
	size_ghosts = ghosts.size();
	last_idx=nocts-1;

	//Clean index of ghost brothers in case of coarsening a broken family
	last_ghost_bros.clear();

	// Set index for start and end check for ghosts
	if (ghosts.size()){
		while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() <= last_desc.computeMorton()){
			idx2_gh++;
		}
		idx2_gh = min((size_ghosts-1), idx2_gh);

		while(idx1_gh < size_ghosts && ghosts[idx1_gh].computeMorton() <= octants[0].computeMorton()){
			idx1_gh++;
		}
		idx1_gh-=1;
		if (idx1_gh==-1) idx1_gh=0;
	}

	// End on ghosts
	if (ghosts.size() && nocts > 0){
		if (ghosts[idx1_gh].buildFather(global.MAX_LEVEL)==octants[0].buildFather(global.MAX_LEVEL)){
			father = ghosts[idx1_gh].buildFather(global.MAX_LEVEL);
			nbro = 0;
			idx = idx1_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 && ghosts[idx].buildFather(global.MAX_LEVEL) == father){

				//Add ghost index to structure for mapper in case of coarsening a broken family
				last_ghost_bros.push_back(idx);

				nbro++;
				if (idx==0)
					break;
				idx--;
				marker = ghosts[idx].getMarker();
			}
			idx = 0;
			while(idx<nocts && octants[idx].buildFather(global.MAX_LEVEL) == father){
				if (octants[idx].getMarker()<0)
					nbro++;
				idx++;
				if(idx==nocts)
					break;
			}
			if (nbro != global.nchildren && idx!=nocts-1){
				for(uint32_t ii=0; ii<idx; ii++){
					if (octants[ii].getMarker()<0){
						octants[ii].setMarker(0);
						octants[ii].info[15]=true;
						Bdone=true;
						newmodified.push_back(ii);
					}
				}
				//Clean index of ghost brothers in case of coarsening a broken family
				last_ghost_bros.clear();
			}
		}

		if (ghosts[idx2_gh].buildFather(global.MAX_LEVEL)==octants[nocts-1].buildFather(global.MAX_LEVEL)){
			father = ghosts[idx2_gh].buildFather(global.MAX_LEVEL);
			nbro = 0;
			idx = idx2_gh;
			marker = ghosts[idx].getMarker();
			while(marker < 0 && ghosts[idx].buildFather(global.MAX_LEVEL) == father){
				nbro++;
				idx++;
				if(idx == size_ghosts){
					break;
				}
				marker = ghosts[idx].getMarker();
			}
			idx = nocts-1;
			while(octants[idx].buildFather(global.MAX_LEVEL) == father){
				if (octants[idx].getMarker()<0)
					nbro++;
				idx--;
				if (idx==0)
					break;
			}
			last_idx=idx;
			if (nbro != global.nchildren && idx!=nocts-1){
				for(uint32_t ii=idx+1; ii<nocts; ii++){
					if (octants[ii].getMarker()<0){
						octants[ii].setMarker(0);
						octants[ii].info[15]=true;
						Bdone=true;
						newmodified.push_back(ii);
					}
				}
			}
		}
	}

	// Check first internal octants
	father = octants[0].buildFather(global.MAX_LEVEL);
	lastdesc = father.buildLastDesc(global.MAX_LEVEL);
	mortonld = lastdesc.computeMorton();
	nbro = 0;
	for (idx=0; idx<global.nchildren; idx++){
		// Check if family is complete or to be checked in the internal loop (some brother refined)
		if (octants[idx].computeMorton() <= mortonld){
			nbro++;
		}
	}
	if (nbro != global.nchildren)
		idx0 = nbro;

	// Check and coarse internal octants
	for (idx=idx0; idx<nocts; idx++){
		if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
			nbro = 0;
			father = octants[idx].buildFather(global.MAX_LEVEL);
			// Check if family is to be coarsened
			for (idx2=idx; idx2<idx+global.nchildren; idx2++){
				if (idx2<nocts){
					if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather(global.MAX_LEVEL) == father){
						nbro++;
					}
				}
			}
			if (nbro == global.nchildren){
				idx = idx2-1;
			}
			else{
				if (idx<=last_idx){
					octants[idx].setMarker(0);
					octants[idx].info[15]=true;
					Bdone=true;
					newmodified.push_back(idx);
				}
			}
		}
	}
};

// =================================================================================== //

/*! 2:1 balancing on level a local tree already adapted (balance only the octants with info[14] = false) (refinement wins!)
 * Return true if balanced done with some markers modification
 * Seto doInterior = false if the interior octants are already balanced
 */
bool classLocalTree::localBalance(bool doInterior){

	uint32_t			sizeneigh, modsize;
	u32vector		 	neigh;
	u32vector		 	modified, newmodified;
	uint32_t 			i, idx;
	uint8_t				iface, iedge, inode;
	int8_t				targetmarker;
	vector<bool> 		isghost;
	bool				Bdone = false;
	bool				Bedge = ((balance_codim>1) && (dim==3));
	bool				Bnode = (balance_codim==dim);

	octvector::iterator 	obegin, oend, it;
	u32vector::iterator 	ibegin, iend, iit;

	//If interior octants have to be balanced
	if(doInterior){
		// First loop on the octants
		obegin = octants.begin();
		oend = octants.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && it->getMarker() != 0){
				targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel() + octants[idx].getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(!it->getBound(iface)){
						findNeighbours(idx, iface, neigh, isghost);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if (!isghost[i]){
								{
									if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
									else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
										octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
										octants[neigh[i]].info[15] = true;
										modified.push_back(neigh[i]);
										Bdone = true;
									}
								};
							}
							else{
								{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								};
							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(!it->getBound(global.edgeface[iedge][0]) && !it->getBound(global.edgeface[iedge][1])){
							findEdgeNeighbours(idx, iedge, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
											octants[idx].info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(!it->getBound(global.nodeface[inode][0]) && !it->getBound(global.nodeface[inode][1]) && !it->getBound(global.nodeface[inode][dim-1])){
							findNodeNeighbours(idx, inode, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
											octants[idx].info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
					}
				}

			}
			idx++;
		}
		// Loop on ghost octants (influence over interior borders)
		obegin = ghosts.begin();
		oend = ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && it->getMarker() != 0){
				targetmarker = min(global.MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
								octants[neigh[i]].info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(it->getPbound(global.edgeface[iedge][0]) == true || it->getPbound(global.edgeface[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(it->getPbound(global.nodeface[inode][0]) == true || it->getPbound(global.nodeface[inode][1]) == true || it->getPbound(global.nodeface[inode][dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!octants[idx].getNotBalance()){
					targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel()+octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<global.nfaces; iface++){
						if(!octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
											octants[idx].info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<global.nedges; iedge++){
							//if(!octants[idx].getPbound(global.edgeface[iedge][0]) || !octants[idx].getPbound(global.edgeface[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<global.nnodes; inode++){
							//if(!octants[idx].getPbound(global.nodeface[inode][0]) || !octants[idx].getPbound(global.nodeface[inode][1]) || !octants[idx].getPbound(global.nodeface[inode][dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while

	}
	else{

		// Loop on ghost octants (influence over interior borders)
		obegin = ghosts.begin();
		oend = ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && it->info[15]){
				targetmarker = min(global.MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
								octants[neigh[i]].info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(it->getPbound(global.edgeface[iedge][0]) == true || it->getPbound(global.edgeface[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(it->getPbound(global.nodeface[inode][0]) == true || it->getPbound(global.nodeface[inode][1]) == true || it->getPbound(global.nodeface[inode][dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!octants[idx].getNotBalance()){
					targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel()+octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<global.nfaces; iface++){
						if(!octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
											octants[idx].info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<global.nedges; iedge++){
							//if(!octants[idx].getPbound(global.edgeface[iedge][0]) || !octants[idx].getPbound(global.edgeface[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<global.nnodes; inode++){
							//if(!octants[idx].getPbound(global.nodeface[inode][0]) || !octants[idx].getPbound(global.nodeface[inode][1]) || !octants[idx].getPbound(global.nodeface[inode][dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while
		obegin = oend = octants.end();
		ibegin = iend = modified.end();
	}
	return Bdone;
	// Pay attention : info[15] may be true after local balance for some octants
};

// =================================================================================== //

/*! 2:1 balancing on level a local tree already adapted (balance only the octants with info[14] = false) (refinement wins!)
 * Return true if balanced done with some markers modification
 * Seto doInterior = false if the interior octants are already balanced
 */
bool classLocalTree::localBalanceAll(bool doInterior){
	// Local variables
	uint32_t			sizeneigh, modsize;
	u32vector		 	neigh;
	u32vector		 	modified, newmodified;
	uint32_t 			i, idx;
	uint8_t				iface, iedge, inode;
	int8_t				targetmarker;
	vector<bool> 		isghost;
	bool				Bdone = false;
	bool				Bedge = ((balance_codim>1) && (dim==3));
	bool				Bnode = (balance_codim==dim);

	octvector::iterator 	obegin, oend, it;
	u32vector::iterator 	ibegin, iend, iit;


	//If interior octants have to be balanced
	if(doInterior){
		// First loop on the octants
		obegin = octants.begin();
		oend = octants.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if ((!it->getNotBalance()) && ((it->info[15]) || (it->getMarker()!=0) || ((it->getIsNewC()) || (it->getIsNewR())))){
				targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel() + octants[idx].getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(!it->getBound(iface)){
						findNeighbours(idx, iface, neigh, isghost);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if (!isghost[i]){
								{
									if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
									else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
										octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
										octants[neigh[i]].info[15] = true;
										modified.push_back(neigh[i]);
										Bdone = true;
									}
								};
							}
							else{
								{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								};

							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(!it->getBound(global.edgeface[iedge][0]) && !it->getBound(global.edgeface[iedge][1])){
							findEdgeNeighbours(idx, iedge, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
											octants[idx].info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(!it->getBound(global.nodeface[inode][0]) && !it->getBound(global.nodeface[inode][1]) && !it->getBound(global.nodeface[inode][dim-1])){
							findNodeNeighbours(idx, inode, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
											octants[idx].info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-1-octants[idx].getLevel());
										octants[idx].info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
					}
				}

			}
			idx++;
		}
		// Loop on ghost octants (influence over interior borders)
		obegin = ghosts.begin();
		oend = ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && (it->info[15] || (it->getIsNewC() || it->getIsNewR()))){
				targetmarker = min(global.MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
								octants[neigh[i]].info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(it->getPbound(global.edgeface[iedge][0]) == true || it->getPbound(global.edgeface[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(it->getPbound(global.nodeface[inode][0]) == true || it->getPbound(global.nodeface[inode][1]) == true || it->getPbound(global.nodeface[inode][dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!octants[idx].getNotBalance()){
					targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel()+octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<global.nfaces; iface++){
						if(!octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
											octants[idx].info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<global.nedges; iedge++){
							//if(!octants[idx].getPbound(global.edgeface[iedge][0]) || !octants[idx].getPbound(global.edgeface[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<global.nnodes; inode++){
							//if(!octants[idx].getPbound(global.nodeface[inode][0]) || !octants[idx].getPbound(global.nodeface[inode][1]) || !octants[idx].getPbound(global.nodeface[inode][dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while

	}
	else{

		// Loop on ghost octants (influence over interior borders)
		obegin = ghosts.begin();
		oend = ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && (it->info[15] || (it->getIsNewC() || it->getIsNewR()))){
				targetmarker = min(global.MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<global.nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
								octants[neigh[i]].info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<global.nedges; iedge++){
						//if(it->getPbound(global.edgeface[iedge][0]) == true || it->getPbound(global.edgeface[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<global.nnodes; inode++){
						//if(it->getPbound(global.nodeface[inode][0]) == true || it->getPbound(global.nodeface[inode][1]) == true || it->getPbound(global.nodeface[inode][dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
									octants[neigh[i]].info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!octants[idx].getNotBalance()){
					targetmarker = min(global.MAX_LEVEL, int8_t(octants[idx].getLevel()+octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<global.nfaces; iface++){
						if(!octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
											octants[idx].info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
											octants[neigh[i]].info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<global.nedges; iedge++){
							//if(!octants[idx].getPbound(global.edgeface[iedge][0]) || !octants[idx].getPbound(global.edgeface[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<global.nnodes; inode++){
							//if(!octants[idx].getPbound(global.nodeface[inode][0]) || !octants[idx].getPbound(global.nodeface[inode][1]) || !octants[idx].getPbound(global.nodeface[inode][dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-octants[idx].getLevel()-1);
												octants[idx].info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												octants[neigh[i]].setMarker(targetmarker-octants[neigh[i]].getLevel()-1);
												octants[neigh[i]].info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while
		obegin = oend = octants.end();
		ibegin = iend = modified.end();
	}
	return Bdone;
	// Pay attention : info[15] may be true after local balance for some octants
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iedge in vector octants.
 * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findEdgeNeighbours(uint32_t idx,
		uint8_t iedge,
		u32vector & neighbours,
		vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	classOctant* 	oct = &octants[idx];
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = global.edgecoeffs[iedge][0];
	int8_t cy = global.edgecoeffs[iedge][1];
	int8_t cz = global.edgecoeffs[iedge][2];

	isghost.clear();
	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge < 0 || iedge > global.nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = global.edgeface[iedge][0];
	iface2 = global.edgeface[iedge][1];

	// Check if octants edge is a boundary
	if (oct->info[iface1] == false && oct->info[iface2] == false){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

		//SEARCH IN GHOSTS

		if (ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(size_ghosts/2);
			classOctant* octghost = &ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = int32_t(idxghost/2);
			idxtry = uint32_t(((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
			if (idxtry > ghosts.size()-1) idxtry = ghosts.size()-1;
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > ghosts.size()-1){
					if (jump > 0){
						idxtry = ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
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
						if(idxtry > ghosts.size()-1){
							idxtry = ghosts.size()-1;
							break;
						}
					}
					while(ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < size_ghosts){
					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = ghosts[idxtry].computeMorton();
					while(Mortontry < Mortonlast && idxtry < ghosts.size()){
						Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(ghosts[idxtry].x));
						Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(ghosts[idxtry].y));
						Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(ghosts[idxtry].z));
						Dxstar = int32_t((cx-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
						Dystar = int32_t((cy-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
						Dzstar = int32_t((cz-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;

						uint32_t x0 = oct->x;
						uint32_t x1 = x0 + size;
						uint32_t y0 = oct->y;
						uint32_t y1 = y0 + size;
						uint32_t z0 = oct->z;
						uint32_t z1 = z0 + size;
						uint32_t x0try = ghosts[idxtry].x;
						uint32_t x1try = x0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint32_t y0try = ghosts[idxtry].y;
						uint32_t y1try = y0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint32_t z0try = ghosts[idxtry].z;
						uint32_t z1try = z0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint8_t level = oct->level;
						uint8_t leveltry = ghosts[idxtry].getLevel();

						if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
							if (leveltry > level){
								if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
							if (leveltry < level){
								if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
						}

						idxtry++;
						if(idxtry>size_ghosts-1){
							break;
						}
						Mortontry = ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		Mortontry = oct->computeMorton();
		int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
		idxtry = uint32_t(idx +((Mortontry<Morton)-(Mortontry>Morton))*jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > octants.size()-1){
				if (jump > 0){
					idxtry = octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
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
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
					Dxstar = int32_t((cx-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = octants[idxtry].x;
					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t y0try = octants[idxtry].y;
					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t z0try = octants[idxtry].z;
					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint8_t level = oct->level;
					uint8_t leveltry = octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iedge in vector octants.
 * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findEdgeNeighbours(classOctant* oct,
		uint8_t iedge,
		u32vector & neighbours,
		vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = global.edgecoeffs[iedge][0];
	int8_t cy = global.edgecoeffs[iedge][1];
	int8_t cz = global.edgecoeffs[iedge][2];

	isghost.clear();
	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge < 0 || iedge > global.nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = global.edgeface[iedge][0];
	iface2 = global.edgeface[iedge][1];

	// Check if octants edge is a boundary
	if (oct->info[iface1] == false && oct->info[iface2] == false){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

		//SEARCH IN GHOSTS

		if (ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(size_ghosts/2);
			classOctant* octghost = &ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = int32_t(idxghost/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > ghosts.size()-1){
					if (jump > 0){
						idxtry = ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
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
						if(idxtry > ghosts.size()-1){
							idxtry = ghosts.size()-1;
							break;
						}
					}
					while(ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < size_ghosts){
					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = ghosts[idxtry].computeMorton();
					while(Mortontry < Mortonlast && idxtry < ghosts.size()){
						Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(ghosts[idxtry].x));
						Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(ghosts[idxtry].y));
						Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(ghosts[idxtry].z));
						Dxstar = int32_t((cx-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
						Dystar = int32_t((cy-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
						Dzstar = int32_t((cz-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;

						uint32_t x0 = oct->x;
						uint32_t x1 = x0 + size;
						uint32_t y0 = oct->y;
						uint32_t y1 = y0 + size;
						uint32_t z0 = oct->z;
						uint32_t z1 = z0 + size;
						uint32_t x0try = ghosts[idxtry].x;
						uint32_t x1try = x0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint32_t y0try = ghosts[idxtry].y;
						uint32_t y1try = y0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint32_t z0try = ghosts[idxtry].z;
						uint32_t z1try = z0try + ghosts[idxtry].getSize(global.MAX_LEVEL);
						uint8_t level = oct->level;
						uint8_t leveltry = ghosts[idxtry].getLevel();

						if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
							if (leveltry > level){
								if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
							if (leveltry < level){
								if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
						}

						idxtry++;
						if(idxtry>size_ghosts-1){
							break;
						}
						Mortontry = ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = int32_t((noctants)/2+1);
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > octants.size()-1){
				if (jump > 0){
					idxtry = octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
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
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
					Dxstar = int32_t((cx-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = octants[idxtry].x;
					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t y0try = octants[idxtry].y;
					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t z0try = octants[idxtry].z;
					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint8_t level = oct->level;
					uint8_t leveltry = octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}

};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost through iedge in vector octants.
 * Returns a vector (empty if iedge is not a pbound edge) with the index of neighbours
 * in the structure octants.
 */
void classLocalTree::findGhostEdgeNeighbours(uint32_t idx,
		uint8_t iedge,
		u32vector & neighbours){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	classOctant* 	oct = &ghosts[idx];
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = global.edgecoeffs[iedge][0];
	int8_t cy = global.edgecoeffs[iedge][1];
	int8_t cz = global.edgecoeffs[iedge][2];

	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge < 0 || iedge > global.nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = global.edgeface[iedge][0];
	iface2 = global.edgeface[iedge][1];

	// Check if octants edge is a pboundary edge
	if (oct->info[iface1+6] == true || oct->info[iface2+6] == true){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
		Morton = samesizeoct.computeMorton();

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = getNumOctants()/2;
		idxtry = uint32_t(getNumOctants()/2);
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > octants.size()-1){
				if (jump > 0){
					idxtry = octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->x) + int32_t(octants[idxtry].x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->y) + int32_t(octants[idxtry].y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->z) + int32_t(octants[idxtry].z));
					Dxstar = int32_t((cx-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = octants[idxtry].x;
					uint32_t x1try = x0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t y0try = octants[idxtry].y;
					uint32_t y1try = y0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint32_t z0try = octants[idxtry].z;
					uint32_t z1try = z0try + octants[idxtry].getSize(global.MAX_LEVEL);
					uint8_t level = oct->level;
					uint8_t leveltry = octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
							}
						}
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through inode in vector octants.
 * Returns a vector (empty if inode is a bound node) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findNodeNeighbours(classOctant* oct,
		uint8_t inode,
		u32vector & neighbours,
		vector<bool> & isghost){

	uint64_t  	Morton, Mortontry;
	uint32_t  	noctants = getNumOctants();
	uint32_t 	idxtry;
	uint32_t 	size = oct->getSize(global.MAX_LEVEL);
	uint8_t 	iface1, iface2, iface3;
	//	int32_t 	Dhx, Dhy, Dhz;
	//	int32_t 	Dhxref, Dhyref, Dhzref;

	//Alternative to switch case
	//	int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
	//	int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
	//	int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
	//	int8_t cx = Cx[inode];
	//	int8_t cy = Cy[inode];
	//	int8_t cz = Cz[inode];
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.nodecoeffs[inode][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode < 0 || inode > global.nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = global.nodeface[inode][0];
	iface2 = global.nodeface[inode][1];
	iface3 = global.nodeface[inode][dim-1];

	// Check if octants node is a boundary
	if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cxyz[0]*size, oct->y+cxyz[1]*size, oct->z+cxyz[2]*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

		//SEARCH IN GHOSTS

		if (ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(size_ghosts/2);
			classOctant* octghost = &ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			if (idxtry > size_ghosts-1)
				idxtry = size_ghosts-1;
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > ghosts.size()-1){
					if (jump > 0){
						idxtry = ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
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
						if(idxtry > ghosts.size()-1){
							idxtry = ghosts.size()-1;
							break;
						}
					}
					while(ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < size_ghosts){
					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = ghosts[idxtry].computeMorton();
					int32_t Dx[3] = {0,0,0};
					int32_t Dxstar[3] = {0,0,0};
					u32array3 coord = oct->getCoord();
					u32array3 coordtry = ghosts[idxtry].getCoord();
					u32array3 coord1 = {1,1,1};
					u32array3 coordtry1 = {1,1,1};
					while(Mortontry < Mortonlast && idxtry < size_ghosts){
						//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
						//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
						//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
						//						Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cx>0)*size;
						//						Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cy>0)*size;
						//						Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cz>0)*size;
						for (int idim=0; idim<dim; idim++){
							Dx[idim] 		= abs(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
							Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
							coord1[idim] 	= coord[idim] + size;
							coordtry1[idim] = coordtry[idim] + ghosts[idxtry].getSize(global.MAX_LEVEL);
						}
						if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
							neighbours.push_back(idxtry);
							isghost.push_back(true);
							return;
						}
						idxtry++;
						if(idxtry>size_ghosts-1){
							break;
						}
						Mortontry = ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = (int32_t((noctants)/2+1));
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > octants.size()-1){
				if (jump > 0){
					idxtry = octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
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
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = octants[idxtry].getCoord();
				u32array3 coord1 = {1,1,1};
				u32array3 coordtry1 = {1,1,1};
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					//					Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
					//					Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
					//					Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
					//					Dhxref = int32_t(cx<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cx>0)*size;
					//					Dhyref = int32_t(cy<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cy>0)*size;
					//					Dhzref = int32_t(cz<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cz>0)*size;
					for (int idim=0; idim<dim; idim++){
						Dx[idim] 		= abs(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
						return;
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through inode in vector octants.
 * Returns a vector (empty if inode is a bound node) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree
 */
void classLocalTree::findNodeNeighbours(uint32_t idx,
		uint8_t inode,
		u32vector & neighbours,
		vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	classOctant* 	oct = &octants[idx];
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);
	uint8_t 		iface1, iface2, iface3;
	//	int32_t Dhx, Dhy, Dhz;
	//	int32_t Dhxref, Dhyref, Dhzref;

	//Alternative to switch case
	//	int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
	//	int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
	//	int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
	//	int8_t cx = Cx[inode];
	//	int8_t cy = Cy[inode];
	//	int8_t cz = Cz[inode];
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.nodecoeffs[inode][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode < 0 || inode > global.nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = global.nodeface[inode][0];
	iface2 = global.nodeface[inode][1];
	iface3 = global.nodeface[inode][dim-1];

	// Check if octants node is a boundary
	if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cxyz[0]*size, oct->y+cxyz[1]*size, oct->z+cxyz[2]*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

		//SEARCH IN GHOSTS

		if (ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(size_ghosts/2);
			classOctant* octghost = &ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			if (idxtry > size_ghosts-1)
				idxtry = size_ghosts-1;
			while(abs(jump) > 0){
				Mortontry = ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > ghosts.size()-1){
					if (jump > 0){
						idxtry = ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
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
						if(idxtry > ghosts.size()-1){
							idxtry = ghosts.size()-1;
							break;
						}
					}
					while(ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < size_ghosts){
					if(ghosts[idxtry].computeMorton() == Morton && ghosts[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = ghosts[idxtry].computeMorton();
					int32_t Dx[3] = {0,0,0};
					int32_t Dxstar[3] = {0,0,0};
					u32array3 coord = oct->getCoord();
					u32array3 coordtry = ghosts[idxtry].getCoord();
					u32array3 coord1 = {1,1,1};
					u32array3 coordtry1 = {1,1,1};
					while(Mortontry < Mortonlast && idxtry < size_ghosts){
						//						Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
						//						Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
						//						Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
						//						Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cx>0)*size;
						//						Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cy>0)*size;
						//						Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize(global.MAX_LEVEL) + int32_t(cz>0)*size;
						for (int idim=0; idim<dim; idim++){
							Dx[idim] 		= abs(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
							Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(ghosts[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
							coord1[idim] 	= coord[idim] + size;
							coordtry1[idim] = coordtry[idim] + ghosts[idxtry].getSize(global.MAX_LEVEL);
						}
						if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
							neighbours.push_back(idxtry);
							isghost.push_back(true);
							return;
						}
						idxtry++;
						if(idxtry>size_ghosts-1){
							break;
						}
						Mortontry = ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = (int32_t((noctants)/2+1));
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > octants.size()-1){
				if (jump > 0){
					idxtry = octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
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
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = octants[idxtry].getCoord();
				u32array3 coord1 = {1,1,1};
				u32array3 coordtry1 = {1,1,1};
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					//					Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
					//					Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
					//					Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
					//					Dhxref = int32_t(cx<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cx>0)*size;
					//					Dhyref = int32_t(cy<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cy>0)*size;
					//					Dhzref = int32_t(cz<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cz>0)*size;
					for (int idim=0; idim<dim; idim++){
						Dx[idim] 		= abs(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
						return;
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost through inode in vector octants.
 * Returns a vector (empty if inode is not a pbound node) with the index of neighbours
 * in the structure octants.
 */
void classLocalTree::findGhostNodeNeighbours(uint32_t idx,
		uint8_t inode,
		u32vector & neighbours){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	classOctant* 	oct = &ghosts[idx];
	uint32_t 		size = oct->getSize(global.MAX_LEVEL);
	uint8_t 		iface1, iface2, iface3;
	//int32_t Dhx, Dhy, Dhz;
	//int32_t Dhxref, Dhyref, Dhzref;

	//Alternative to switch case
	//int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
	//int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
	//int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
	//int8_t cx = Cx[inode];
	//int8_t cy = Cy[inode];
	//int8_t cz = Cz[inode];
	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<dim; idim++){
		cxyz[idim] = global.nodecoeffs[inode][idim];
	}

	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode < 0 || inode > global.nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = global.nodeface[inode][0];
	iface2 = global.nodeface[inode][1];
	iface3 = global.nodeface[inode][dim-1];

	//		// Check if octants node is a boundary
	//		if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){
	// Check if octants node is a pboundary node
	if (oct->info[iface1+6] == true || oct->info[iface2+6] == true || oct->info[iface3+6] == true){

		//Build Morton number of virtual neigh of same size
		classOctant samesizeoct(dim, oct->level, oct->x+cxyz[0]*size, oct->y+cxyz[1]*size, oct->z+cxyz[2]*size);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
		int32_t jump = noctants/2;
		idxtry = jump;
		while(abs(jump) > 0){
			Mortontry = octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
		}
		if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
					//Found neighbour of same size
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				classOctant last_desc = samesizeoct.buildLastDesc(global.MAX_LEVEL);
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = octants[idxtry].getCoord();
				u32array3 coord1 = {1,1,1};
				u32array3 coordtry1 = {1,1,1};
				while(Mortontry < Mortonlast && idxtry <= noctants-1){
					//				Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
					//				Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
					//				Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
					//				Dhxref = int32_t(cx<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cx>0)*size;
					//				Dhyref = int32_t(cy<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cy>0)*size;
					//				Dhzref = int32_t(cz<0)*octants[idxtry].getSize(global.MAX_LEVEL) + int32_t(cz>0)*size;
					for (int idim=0; idim<dim; idim++){
						Dx[idim] 		= abs(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(octants[idxtry].getSize(global.MAX_LEVEL)) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + octants[idxtry].getSize(global.MAX_LEVEL);
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[dim-1] == Dxstar[dim-1]){
						neighbours.push_back(idxtry);
					}
					idxtry++;
					Mortontry = octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

void classLocalTree::computeIntersections() {

		octvector::iterator 	it, obegin, oend;
		classIntersection 		intersection;
		u32vector 				neighbours;
		vector<bool>			isghost;
		uint32_t 				counter, idx;
		uint32_t 				i, nsize;
		uint8_t 				iface, iface2;

		intersections.clear();
		intersections.reserve(2*3*octants.size());

		counter = idx = 0;

		// Loop on ghosts
		obegin = ghosts.begin();
		oend = ghosts.end();
		for (it = obegin; it != oend; it++){
			for (iface = 0; iface < dim; iface++){
				iface2 = iface*2;
				findGhostNeighbours(idx, iface2, neighbours);
				nsize = neighbours.size();
				for (i = 0; i < nsize; i++){
					intersection.dim = dim;
					intersection.finer = getGhostLevel(idx) >= getLevel((int)neighbours[i]);
					intersection.owners[0]  = neighbours[i];
					intersection.owners[1] = idx;
					intersection.iface = global.oppface[iface2] - (getGhostLevel(idx) >= getLevel((int)neighbours[i]));
					intersection.isnew = false;
					intersection.isghost = true;
					intersection.bound = false;
					intersection.pbound = true;
					intersections.push_back(intersection);
					counter++;
				}
			}
			idx++;
		}

		// Loop on octants
		idx=0;
		obegin = octants.begin();
		oend = octants.end();
		for (it = obegin; it != oend; it++){
			for (iface = 0; iface < dim; iface++){
				iface2 = iface*2;
				findNeighbours(idx, iface2, neighbours, isghost);
				nsize = neighbours.size();
				if (nsize) {
					for (i = 0; i < nsize; i++){
						if (isghost[i]){
							intersection.dim = dim;
							intersection.owners[0] = idx;
							intersection.owners[1] = neighbours[i];
							intersection.finer = (nsize>1);
							intersection.iface = iface2 + (nsize>1);
							intersection.isnew = false;
							intersection.isghost = true;
							intersection.bound = false;
							intersection.pbound = true;
							intersections.push_back(intersection);
							counter++;
						}
						else{
							intersection.dim = dim;
							intersection.owners[0] = idx;
							intersection.owners[1] = neighbours[i];
							intersection.finer = (nsize>1);
							intersection.iface = iface2 + (nsize>1);
							intersection.isnew = false;
							intersection.isghost = false;
							intersection.bound = false;
							intersection.pbound = false;
							intersections.push_back(intersection);
							counter++;
						}
					}
				}
				else{
					intersection.dim = dim;
					intersection.owners[0] = idx;
					intersection.owners[1] = idx;
					intersection.finer = 0;
					intersection.iface = iface2;
					intersection.isnew = false;
					intersection.isghost = false;
					intersection.bound = true;
					intersection.pbound = false;
					intersections.push_back(intersection);
					counter++;
				}
				if (it->info[iface2+1]){
					intersection.dim = dim;
					intersection.owners[0] = idx;
					intersection.owners[1] = idx;
					intersection.finer = 0;
					intersection.iface = iface2+1;
					intersection.isnew = false;
					intersection.isghost = false;
					intersection.bound = true;
					intersection.pbound = false;
					intersections.push_back(intersection);
					counter++;
				}
			}
			idx++;
		}
		intervector(intersections).swap(intersections);
	}

// =================================================================================== //
/*! Find an input Morton in octants and return the local idx
 * Return nocts if target Morton not found.
*/
uint32_t classLocalTree::findMorton(uint64_t Morton){

	uint32_t 		nocts = octants.size();
	uint32_t 		idx = nocts/2;
	uint64_t 		Mortontry = octants[idx].computeMorton();
	int32_t 		jump = nocts/2;

	while(abs(jump)>0){
		if (Mortontry == Morton){
			return idx;
		}
		Mortontry = octants[idx].computeMorton();
		jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
		idx += jump;
		if (idx > nocts){
			return nocts-1;
		}
	}
	if (Mortontry<Morton){
		for (uint32_t idx2=idx; idx2<nocts; idx2++){
			Mortontry = octants[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	else{
		for(uint32_t idx2=0; idx2<idx+1; idx2++){
			Mortontry = octants[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	return nocts-1;
};

// =================================================================================== //
/*! Find an input Morton in ghosts and return the local idx
* Return nghosts if target Morton not found.
*/
uint32_t classLocalTree::findGhostMorton(uint64_t Morton){
	uint32_t 		nocts = ghosts.size();
	uint32_t 		idx = nocts/2;
	uint64_t 		Mortontry = ghosts[idx].computeMorton();
	int32_t 		jump = nocts/2;

	while(abs(jump)>0){
		if (Mortontry == Morton){
			return idx;
		}
		Mortontry = ghosts[idx].computeMorton();
		jump = (Mortontry<Morton)*jump/4;
		idx += jump;
		if (idx > nocts){
			return nocts;
		}
	}
	if (Mortontry<Morton){
		for (uint32_t idx2=idx; idx2<nocts; idx2++){
			Mortontry = ghosts[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	else{
		for(uint32_t idx2=0; idx2<idx; idx2++){
			Mortontry = ghosts[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	return nocts;
};

// =================================================================================== //

/** Compute the connectivity of octants and store the coordinates of nodes.
 */
void classLocalTree::computeConnectivity(){

	map<uint64_t, vector<uint32_t> > 			mapnodes;
	map<uint64_t, vector<uint32_t> >::iterator 	iter, iterend;
	uint32_t 									i, k, counter;
	uint64_t 									morton;
	uint32_t 									noctants = getNumOctants();
	u32vector2D 								octnodes;
	uint8_t 									j;

	clearConnectivity();
	octnodes.reserve(global.nnodes);

	if (nodes.size() == 0){
		connectivity.resize(noctants);
		for (i = 0; i < noctants; i++){
			octants[i].getNodes(octnodes, global.MAX_LEVEL);
			for (j = 0; j < global.nnodes; j++){
				morton = keyXYZ(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
				if (mapnodes[morton].size()==0){
					mapnodes[morton].reserve(16);
					for (k = 0; k < 3; k++){
						mapnodes[morton].push_back(octnodes[j][k]);
					}
				}
				mapnodes[morton].push_back(i);
			}
			u32vector2D().swap(octnodes);
		}
		iter	= mapnodes.begin();
		iterend	= mapnodes.end();
		counter = 0;
		uint32_t numnodes = mapnodes.size();
		nodes.resize(numnodes);
		while (iter != iterend){
			u32vector nodecasting(iter->second.begin(), iter->second.begin()+3);
			for (k=0; k<3; k++){
				nodes[counter][k] = nodecasting[k];
			}
			u32array3(nodes[counter]).swap(nodes[counter]);

			for(u32vector::iterator iter2 = iter->second.begin()+3; iter2 != iter->second.end(); iter2++){
				if (connectivity[(*iter2)].size()==0){
					connectivity[(*iter2)].reserve(8);
				}
				connectivity[(*iter2)].push_back(counter);
			}
			mapnodes.erase(iter++);
			counter++;
		}
		u32arr3vector(nodes).swap(nodes);

		u32vector2D(connectivity).swap(connectivity);

	}

	map<uint64_t, vector<uint32_t> >().swap(mapnodes);
	iter = mapnodes.end();

};

// =================================================================================== //
/*! Clear nodes vector and connectivity of octants of local tree
*/
void classLocalTree::clearConnectivity(){
	u32arr3vector().swap(nodes);
	u32vector2D().swap(connectivity);
};

// =================================================================================== //
/*! Updates nodes vector and connectivity of octants of local tree
*/
void classLocalTree::updateConnectivity(){
	clearConnectivity();
	computeConnectivity();
};

// =================================================================================== //
/*! Computes ghosts nodes vector and connectivity of ghosts octants of local tree
*/
void classLocalTree::computeGhostsConnectivity(){

	map<uint64_t, vector<uint32_t> > 			mapnodes;
	map<uint64_t, vector<uint32_t> >::iterator 	iter, iterend;
	uint32_t 									i, k, counter;
	uint64_t 									morton;
	uint32_t 									noctants = size_ghosts;
	u32vector2D									octnodes;
	uint8_t 									j;

	octnodes.reserve(global.nnodes);
	if (ghostsnodes.size() == 0){
		ghostsconnectivity.resize(noctants);
		for (i = 0; i < noctants; i++){
			ghosts[i].getNodes(octnodes, global.MAX_LEVEL);
			for (j = 0; j < global.nnodes; j++){
				morton = keyXYZ(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
				if (mapnodes[morton].size()==0){
					mapnodes[morton].reserve(16);
					for (k = 0; k < 3; k++){
						mapnodes[morton].push_back(octnodes[j][k]);
					}
				}
				mapnodes[morton].push_back(i);
			}
			u32vector2D().swap(octnodes);
		}
		iter	= mapnodes.begin();
		iterend	= mapnodes.end();
		uint32_t numnodes = mapnodes.size();
		ghostsnodes.resize(numnodes);
		counter = 0;
		while (iter != iterend){
			vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+3);
			for (k=0; k<3; k++){
				ghostsnodes[counter][k] = nodecasting[k];
			}
			u32array3(ghostsnodes[counter]).swap(ghostsnodes[counter]);
			for(vector<uint32_t>::iterator iter2 = iter->second.begin()+3; iter2 != iter->second.end(); iter2++){
				if (ghostsconnectivity[(*iter2)].size()==0){
					ghostsconnectivity[(*iter2)].reserve(8);
				}
				ghostsconnectivity[(*iter2)].push_back(counter);
			}
			mapnodes.erase(iter++);
			counter++;
		}
		u32arr3vector(ghostsnodes).swap(ghostsnodes);

		u32vector2D(ghostsconnectivity).swap(ghostsconnectivity);

	}
	iter = mapnodes.end();

};

// =================================================================================== //
/*! Clear ghosts nodes vector and connectivity of ghosts octants of local tree
*/
void classLocalTree::clearGhostsConnectivity(){
	u32arr3vector().swap(ghostsnodes);
	u32vector2D().swap(ghostsconnectivity);
};

// =================================================================================== //
/*! Update ghosts nodes vector and connectivity of ghosts octants of local tree
*/
void classLocalTree::updateGhostsConnectivity(){
	clearGhostsConnectivity();
	computeGhostsConnectivity();
};

// =================================================================================== //


