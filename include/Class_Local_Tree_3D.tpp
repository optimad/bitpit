/*!
 * Class_Local_Tree_3D.tpp
 *
 *  \date		23/apr/2014
 *	\authors	Edoardo Lombardi
 *	\authors	Marco Cisternino
 *	\version	0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Local octree portion for each process - 3D specialization
 *
 *	Local tree consists mainly of two vectors with:
 *	- actual octants stored on current process;
 *	- ghost octants neighbours of the first ones.
 *
 *	The octants (and ghosts) are ordered following the Z-curve defined by the Morton index.
 *
 *	Optionally in local tree three vectors of intersections are stored:
 *	- intersections located on the bord of the physical domain of the octree;
 *	- intersections of process bord (i.e. between octants and ghosts);
 *	- intersections completely located in the domain of the process (i.e. between actual octants).
 *
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //
template<>
class Class_Local_Tree<3>{
	// ------------------------------------------------------------------------------- //
	// FRIENDSHIPS ------------------------------------------------------------------- //
	template<int dim> friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<Class_Octant<3> > 		OctantsType;
	typedef vector<Class_Intersection<3> > 	IntersectionsType;
	typedef vector<uint32_t>				u32vector;
	typedef vector<uint64_t>				u64vector;
	typedef vector<vector<uint32_t>	>		u32vector2D;
	typedef vector<vector<uint64_t>	>		u64vector2D;


	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	OctantsType					octants;			/**< Local vector of octants ordered with Morton Number */
	OctantsType					ghosts;				/**< Local vector of ghost octants ordered with Morton Number */
	IntersectionsType			intersections;		/**< Local vector of intersections */
	u64vector 					globalidx_ghosts;	/**< Global index of the ghost octants (size = size_ghosts) */
	Class_Octant<3> 			first_desc;			/**< First (Morton order) most refined octant possible in local partition */
	Class_Octant<3> 			last_desc;			/**< Last (Morton order) most refined octant possible in local partition */
	uint32_t 					size_ghosts;		/**< Size of vector of ghost octants */
	uint8_t						local_max_depth;	/**< Reached max depth in local tree */

	// connectivity
	u32vector2D					nodes;				/**<Local vector of nodes (x,y,z) ordered with Morton Number*/
	u32vector2D					connectivity;		/**<Local vector of connectivity (node1, node2, ...) ordered with Morton-order.
	 *The nodes are stored as index of vector nodes*/
	u32vector2D					ghostsnodes;		/**<Local vector of ghosts nodes (x,y,z) ordered with Morton Number*/
	u32vector2D					ghostsconnectivity;	/**<Local vector of ghosts connectivity (node1, node2, ...) ordered with Morton-order.
	 *The nodes are stored as index of vector nodes*/

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //

public:
	Class_Local_Tree(){
		Class_Octant<3> oct0;
		Class_Octant<3> octf(MAX_LEVEL_3D,0,0,0);
		Class_Octant<3> octl(MAX_LEVEL_3D,global3D.max_length-1,global3D.max_length-1,global3D.max_length-1);
		octants.resize(1);
		octants[0] = oct0;
		first_desc = octf;
		last_desc = octl;
		size_ghosts = 0;
		local_max_depth = 0;

	};
	~Class_Local_Tree(){};

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

private:
	const Class_Octant<3> &  getFirstDesc() const{
		return first_desc;
	};
	const Class_Octant<3> &  getLastDesc() const{
		return last_desc;
	};
	uint32_t getSizeGhost() const{
		return size_ghosts;
	};
	uint32_t getNumOctants() const{
		return octants.size();
	};
	uint8_t getLocalMaxDepth() const{							// Get max depth reached in local tree
		return local_max_depth;
	};
	uint8_t getMarker(int32_t idx){								// Get refinement/coarsening marker for idx-th octant
		return octants[idx].getMarker();
	};
	uint8_t getLevel(int32_t idx){								// Get refinement/coarsening marker for idx-th octant
		return octants[idx].getLevel();
	};
	bool getBalance(int32_t idx){								// Get if balancing-blocked idx-th octant
		return octants[idx].getNotBalance();
	};

	void setMarker(int32_t idx, int8_t marker){					// Set refinement/coarsening marker for idx-th octant
		octants[idx].setMarker(marker);
	};
	void setBalance(int32_t idx, bool balance){					// Set if balancing-blocked idx-th octant
		octants[idx].setBalance(balance);
	};
	void setFirstDesc(){
		OctantsType::const_iterator firstOctant = octants.begin();
		first_desc = Class_Octant<3>(MAX_LEVEL_3D,firstOctant->x,firstOctant->y,firstOctant->z);
	};
	void setLastDesc(){
		OctantsType::const_iterator lastOctant = octants.end() - 1;
		uint32_t x,y,z,delta;
		delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - lastOctant->level)) - 1;
		x = lastOctant->x + delta;
		y = lastOctant->y + delta;
		z = lastOctant->z + delta;
		last_desc = Class_Octant<3>(MAX_LEVEL_3D,x,y,z);

	};

	//-------------------------------------------------------------------------------- //
	// Other methods ----------------------------------------------------------------- //

	Class_Octant<3>& extractOctant(uint32_t idx){
		return octants[idx];
	};

	const Class_Octant<3>&	extractOctant(uint32_t idx) const{
		return octants[idx];
	};

	// =================================================================================== //

	bool refine(){									// Refine local tree: refine one time octants with marker >0

		// Local variables
		vector<uint32_t> last_child_index;
		Class_Octant<3>* children;
		uint32_t idx, nocts, ilastch;
		uint32_t offset = 0, blockidx;
		uint8_t nchm1 = global3D.nchildren-1, ich, iface;
		bool dorefine = false;

		nocts = octants.size();
		for (idx=0; idx<nocts; idx++){
			if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL_3D){
				last_child_index.push_back(idx+nchm1+offset);
				offset += nchm1;
			}
			else{
				//			octants[idx].info[12] = false;
				if (octants[idx].marker > 0)
					octants[idx].marker = 0;
				octants[idx].info[15] = true;
			}
		}
		if (offset > 0){
			octants.resize(octants.size()+offset);
			blockidx = last_child_index[0]-nchm1;
			idx = octants.size();
			ilastch = last_child_index.size()-1;
			while (idx>blockidx){
				idx--;
				//TODO Sostituire questo if con il controllo su last_index_child
				if(idx == last_child_index[ilastch]){
					//				if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL_3D){
					children = octants[idx-offset].buildChildren();
					for (ich=0; ich<global3D.nchildren; ich++){
						octants[idx-ich] = (children[nchm1-ich]);
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
					delete []children;
					if (ilastch != 0){
						ilastch--;
					}
				}
				else {
					octants[idx] = octants[idx-offset];
				}
			}
		}
		octants.shrink_to_fit();
		nocts = octants.size();

		setFirstDesc();
		setLastDesc();

		return dorefine;

	};

	// =================================================================================== //

	bool coarse(){												// Coarse local tree: coarse one time family of octants with marker <0
		// Local variables										// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
		vector<uint32_t> first_child_index;
		Class_Octant<3> father;
		uint32_t ich, nocts, nghosts;
		int32_t idx, idx2;
		uint32_t offset;
		int32_t idx1_gh, idx2_gh;
		uint32_t nidx;
		int8_t markerfather, marker;
		uint8_t nbro, nstart, nend;
		uint8_t nchm1 = global3D.nchildren-1, iface;
		bool docoarse = false;

		//------------------------------------------ //
		// Initialization

		nbro = nstart = nend = 0;
		nidx = offset = 0;

		idx1_gh = idx2_gh = 0;

		nocts   = octants.size();
		size_ghosts = ghosts.size();


		// Init first and last desc (even if already calculated)
		setFirstDesc();
		setLastDesc();

		//------------------------------------------ //

		// Set index for start and end check for ghosts
		if (ghosts.size()){
			while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() <= last_desc.computeMorton()){
				idx2_gh++;
			}
			idx2_gh = min(int(size_ghosts-1), idx2_gh);
		}

		// Check and coarse internal octants
		for (idx=0; idx<nocts; idx++){
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nbro = 0;
				father = octants[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+global3D.nchildren; idx2++){
					if (idx2<nocts){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nbro++;
						}
					}
				}
				if (nbro == global3D.nchildren){
					nidx++;
					first_child_index.push_back(idx);
					idx = idx2-1;
				}
				else{
					if (idx < (nocts>global3D.nchildren)*(nocts-global3D.nchildren)){
						octants[idx].setMarker(0);
						octants[idx].info[11] = true;
					}
				}
			}
			//			else{
			//	//			octants[idx].info[13] = false;
			//			}
		}
		//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
		if (nidx!=0){
			uint32_t nblock = nocts - nidx*nchm1 - nstart;
			nidx = 0;
			for (idx=0; idx<nblock; idx++){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -MAX_LEVEL_3D;
					father = octants[idx+offset].buildFather();
					for(idx2=0; idx2<global3D.nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()+1){
							markerfather = octants[idx+offset+idx2].getMarker()+1;
						}
						for (int iii=0; iii<16; iii++){
							father.info[iii] = father.info[iii] || octants[idx+offset+idx2].info[iii];
						}
					}
					father.info[13] = true;
					if (markerfather < 0){
						docoarse = true;
					}
					father.setMarker(markerfather);
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
		if (ghosts.size() && nocts > 0){
			if ((ghosts[idx2_gh].getMarker() < 0) && (octants[nocts-1].getMarker() < 0)){
				father = ghosts[idx2_gh].buildFather();
				markerfather = ghosts[idx2_gh].getMarker()+1;//-MAX_LEVEL_3D;
				nbro = 0;
				idx = idx2_gh;
				marker = ghosts[idx].getMarker();
				while(marker < 0 && ghosts[idx].buildFather() == father){
					nbro++;
					//TODO CAMBIATO IDX DA CAMBIARE ANCHE NELLE ALTRE COARSE!!!
					if (markerfather < ghosts[idx].getMarker()+1){
						markerfather = ghosts[idx].getMarker()+1;
					}
					idx++;
					if(idx == size_ghosts){
						break;
					}
					marker = ghosts[idx].getMarker();
				}
				nend = 0;
				idx = nocts-1;
				marker = octants[idx].getMarker();
				while(marker < 0 && octants[idx].buildFather() == father && idx >= 0){
					nbro++;
					nend++;
					if (markerfather < octants[idx].getMarker()+1){
						markerfather = octants[idx].getMarker()+1;
					}
					idx--;
					marker = octants[idx].getMarker();
					if (idx<0){
						break;
					}
				}
				if (nbro == global3D.nchildren){
					offset = nend;
				}
				else{
					nend = 0;
					for(int ii=nocts-global3D.nchildren; ii<nocts; ii++){
						octants[ii].setMarker(0);
						octants[ii].info[15] = true;
					}
				}
			}
			if (nend != 0){
				for (int iii=0; iii<16; iii++){
					father.info[iii] = false;
				}
				for (idx=0; idx < nend; idx++){
					for (int iii=0; iii<16; iii++){
						father.info[iii] = father.info[iii] || octants[nocts-idx-1].info[iii];
					}
				}
				father.info[13] = true;
				if (markerfather < 0){
					docoarse = true;
				}
				father.setMarker(markerfather);
				octants.resize(nocts-offset);
				octants.push_back(father);
				octants.shrink_to_fit();
				nocts = octants.size();
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

	bool refine(u32vector & mapidx){							// Refine local tree: refine one time octants with marker >0
		// mapidx[i] = index in old octants vector of the i-th octant (index of father if octant is new after)
		// Local variables
		vector<uint32_t> last_child_index;
		Class_Octant<3>* children;
		uint32_t idx, nocts, ilastch;
		uint32_t offset = 0, blockidx;
		uint8_t nchm1 = global3D.nchildren-1, ich, iface;
		bool dorefine = false;

		nocts = octants.size();
		for (idx=0; idx<nocts; idx++){
			if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL_3D){
				last_child_index.push_back(idx+nchm1+offset);
				offset += nchm1;
			}
			else{
				//			octants[idx].info[12] = false;
				if (octants[idx].marker > 0){
					octants[idx].marker = 0;
					octants[idx].info[15] = false;
				}
			}
		}
		if (offset > 0){
			mapidx.resize(octants.size()+offset);
			mapidx.shrink_to_fit();

			octants.resize(octants.size()+offset);
			blockidx = last_child_index[0]-nchm1;
			idx = octants.size();
			ilastch = last_child_index.size()-1;
			while (idx>blockidx){
				//			while (idx>0){
				idx--;
				//				if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL_3D){
				if(idx == last_child_index[ilastch]){
					children = octants[idx-offset].buildChildren();
					for (ich=0; ich<global3D.nchildren; ich++){
						octants[idx-ich] = (children[nchm1-ich]);
						mapidx[idx-ich]  = mapidx[idx-offset];
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
					delete []children;
					if (ilastch != 0){
						ilastch--;
					}
				}
				else {
					octants[idx] = octants[idx-offset];
					mapidx[idx]  = mapidx[idx-offset];
				}
			}
		}
		octants.shrink_to_fit();
		nocts = octants.size();

		setFirstDesc();
		setLastDesc();

		return dorefine;

	};

	// =================================================================================== //

	bool coarse(u32vector & mapidx){							// Coarse local tree: coarse one time family of octants with marker <0
		// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
		// mapidx[i] = index in old octants vector of the i-th octant (index of father if octant is new after)
		// Local variables
		vector<uint32_t> first_child_index;
		Class_Octant<3> father;
		uint32_t ich, nocts, nghosts, nocts0;
		int32_t idx, idx2;
		uint32_t offset;
		int32_t idx1_gh, idx2_gh;
		uint32_t nidx;
		int8_t markerfather, marker;
		uint8_t nbro, nstart, nend;
		uint8_t nchm1 = global3D.nchildren-1, iface;
		bool docoarse = false;

		//------------------------------------------ //
		// Initialization

		nbro = nstart = nend = 0;
		nidx = offset = 0;

		idx1_gh = idx2_gh = 0;

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
			idx2_gh = min(int(size_ghosts-1), idx2_gh);
		}

		// Check and coarse internal octants
		for (idx=0; idx<nocts; idx++){
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nbro = 0;
				father = octants[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+global3D.nchildren; idx2++){
					if (idx2<nocts){
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nbro++;
						}
					}
				}
				if (nbro == global3D.nchildren){
					nidx++;
					first_child_index.push_back(idx);
					idx = idx2-1;
				}
				else{
					if (idx < (nocts>global3D.nchildren)*(nocts-global3D.nchildren)){
						octants[idx].setMarker(0);
						octants[idx].info[15] = true;
					}
				}
			}
			//			else{
			//	//			octants[idx].info[13] = false;
			//			}
		}
		//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
		if (nidx!=0){
			uint32_t nblock = nocts - nidx*nchm1 - nstart;
			nidx = 0;
			//for (idx=0; idx<nblock; idx++){
			for (idx=0; idx<nocts-offset; idx++){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -MAX_LEVEL_3D;
					father = octants[idx+offset].buildFather();
					for (int iii=0; iii<16; iii++){
						father.info[iii] = false;
					}
					for(idx2=0; idx2<global3D.nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()+1){
							markerfather = octants[idx+offset+idx2].getMarker()+1;
						}
						for (int iii=0; iii<16; iii++){
							father.info[iii] = father.info[iii] || octants[idx+offset+idx2].info[iii];
						}
					}
					father.info[13] = true;
					father.setMarker(markerfather);
					if (markerfather < 0){
						docoarse = true;
					}
					octants[idx] = father;
					mapidx[idx] = mapidx[idx+offset];
					offset += nchm1;
					nidx++;
				}
				else{
					octants[idx] = octants[idx+offset];
					mapidx[idx] = mapidx[idx+offset];
				}
			}
		}
		octants.resize(nocts-offset);
		octants.shrink_to_fit();
		nocts = octants.size();
		mapidx.resize(nocts);
		mapidx.shrink_to_fit();


		// End on ghosts
		if (ghosts.size() && nocts > 0){
			if ((ghosts[idx2_gh].getMarker() < 0) && (octants[nocts-1].getMarker() < 0)){
				father = ghosts[idx2_gh].buildFather();
				markerfather = ghosts[idx2_gh].getMarker()+1;//-MAX_LEVEL_3D;
				nbro = 0;
				idx = idx2_gh;
				marker = ghosts[idx].getMarker();
				while(marker < 0 && ghosts[idx].buildFather() == father){
					nbro++;
					//TODO CAMBIATO IDX DA CAMBIARE ANCHE NELLE ALTRE COARSE!!!
					if (markerfather < ghosts[idx].getMarker()+1){
						markerfather = ghosts[idx].getMarker()+1;
					}
					idx++;
					if(idx == size_ghosts){
						break;
					}
					marker = ghosts[idx].getMarker();
				}
				nend = 0;
				idx = nocts-1;
				marker = octants[idx].getMarker();
				while(marker < 0 && octants[idx].buildFather() == father && idx >= 0){
					nbro++;
					nend++;
					if (markerfather < octants[idx].getMarker()+1){
						markerfather = octants[idx].getMarker()+1;
					}
					idx--;
					marker = octants[idx].getMarker();
					if (idx<0){
						break;
					}
				}
				if (nbro == global3D.nchildren){
					offset = nend;
				}
				else{
					nend = 0;
					for(int ii=nocts-global3D.nchildren; ii<nocts; ii++){
						octants[ii].setMarker(0);
						octants[ii].info[15] = true;
					}
				}
			}
			if (nend != 0){
				for (int iii=0; iii<16; iii++){
					father.info[iii] = false;
				}
				for (idx=0; idx < nend; idx++){
					for (int iii=0; iii<16; iii++){
						father.info[iii] = father.info[iii] || octants[nocts-idx-1].info[iii];
					}
				}
				father.info[13] = true;
				if (markerfather < 0){
					docoarse = true;
				}
				father.setMarker(markerfather);
				octants.resize(nocts-offset);
				octants.push_back(father);
				octants.shrink_to_fit();
				nocts = octants.size();
				mapidx.resize(nocts-offset);
				mapidx.push_back(nocts0-nend);
				mapidx.shrink_to_fit();
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

	// Global refine of octree (one level every element)
	bool globalRefine(){

		// Local variables
		vector<uint32_t> last_child_index;
		Class_Octant<3>* children;
		uint32_t idx, nocts, ilastch;
		uint32_t offset = 0, blockidx;
		uint8_t nchm1 = global3D.nchildren-1, ich, iface;
		bool dorefine = false;

		nocts = octants.size();
		for (idx=0; idx<nocts; idx++){
			octants[idx].setMarker(1);
			if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL_3D){
				last_child_index.push_back(idx+nchm1+offset);
				offset += nchm1;
			}
			else{
				if (octants[idx].marker > 0)
					octants[idx].marker = 0;
				octants[idx].info[15] = true;
			}
		}
		if (offset > 0){
			octants.resize(octants.size()+offset);
			blockidx = last_child_index[0]-nchm1;
			idx = octants.size();
			ilastch = last_child_index.size()-1;
			while (idx>blockidx){
				idx--;
				//TODO Sostituire questo if con il controllo su last_index_child
				if(idx == last_child_index[ilastch]){
					//				if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL_3D){
					children = octants[idx-offset].buildChildren();
					for (ich=0; ich<global3D.nchildren; ich++){
						octants[idx-ich] = (children[nchm1-ich]);
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
					delete []children;
					if (ilastch != 0){
						ilastch--;
					}
				}
				else {
					octants[idx] = octants[idx-offset];
				}
			}
		}
		octants.shrink_to_fit();
		nocts = octants.size();

		setFirstDesc();
		setLastDesc();

		return dorefine;

	};

	// =================================================================================== //

	// Global coarse of octree (every marker set =-1)
	bool globalCoarse(){
		//Local Variables
		vector<uint32_t> first_child_index;
		Class_Octant<3> father;
		uint32_t ich, nocts, nghosts;
		int32_t idx, idx2;
		uint32_t offset;
		int32_t idx1_gh, idx2_gh;
		uint32_t nidx;
		int8_t markerfather, marker;
		uint8_t nbro, nstart, nend;
		uint8_t nchm1 = global3D.nchildren-1, iface;
		bool docoarse = false;

		//------------------------------------------ //
		// Initialization

		nbro = nstart = nend = 0;
		nidx = offset = 0;

		idx1_gh = idx2_gh = 0;

		nocts   = octants.size();
		size_ghosts = ghosts.size();


		// Init first and last desc (even if already calculated)
		setFirstDesc();
		setLastDesc();

		//------------------------------------------ //

		// Set index for start and end check for ghosts
		if (ghosts.size()){
			while(idx2_gh < size_ghosts && ghosts[idx2_gh].computeMorton() <= last_desc.computeMorton()){
				idx2_gh++;
			}
			idx2_gh = min(int(size_ghosts-1), idx2_gh);
		}

		// Check and coarse internal octants
		for (idx=0; idx<nocts; idx++){
			octants[idx].setMarker(-1);
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nbro = 0;
				father = octants[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+global3D.nchildren; idx2++){
					if (idx2<nocts){
						octants[idx2].setMarker(-1);
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nbro++;
						}
					}
				}
				if (nbro == global3D.nchildren){
					nidx++;
					first_child_index.push_back(idx);
					idx = idx2-1;
				}
				else{
					if (idx < (nocts>global3D.nchildren)*(nocts-global3D.nchildren)){
						octants[idx].setMarker(0);
						octants[idx].info[11] = true;
					}
				}
			}
			//			else{
			//	//			octants[idx].info[13] = false;
			//			}
		}
		//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
		if (nidx!=0){
			uint32_t nblock = nocts - nidx*nchm1 - nstart;
			nidx = 0;
			for (idx=0; idx<nblock; idx++){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -MAX_LEVEL_3D;
					father = octants[idx+offset].buildFather();
					for(idx2=0; idx2<global3D.nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()+1){
							markerfather = octants[idx+offset+idx2].getMarker()+1;
						}
						for (int iii=0; iii<16; iii++){
							father.info[iii] = father.info[iii] || octants[idx+offset+idx2].info[iii];
						}
					}
					father.info[13] = true;
					if (markerfather < 0){
						docoarse = true;
					}
					father.setMarker(markerfather);
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
		if (ghosts.size() && nocts > 0){
			ghosts[idx2_gh].setMarker(-1);
			if ((ghosts[idx2_gh].getMarker() < 0) && (octants[nocts-1].getMarker() < 0)){
				father = ghosts[idx2_gh].buildFather();
				markerfather = ghosts[idx2_gh].getMarker()+1;//-MAX_LEVEL_3D;
				nbro = 0;
				idx = idx2_gh;
				ghosts[idx].setMarker(-1);
				marker = ghosts[idx].getMarker();
				while(marker < 0 && ghosts[idx].buildFather() == father){
					nbro++;
					//TODO CAMBIATO IDX DA CAMBIARE ANCHE NELLE ALTRE COARSE!!!
					if (markerfather < ghosts[idx].getMarker()+1){
						markerfather = ghosts[idx].getMarker()+1;
					}
					idx++;
					if(idx == size_ghosts){
						break;
					}
					ghosts[idx].setMarker(-1);
					marker = ghosts[idx].getMarker();
				}
				nend = 0;
				idx = nocts-1;
				octants[idx].setMarker(-1);
				marker = octants[idx].getMarker();
				while(marker < 0 && octants[idx].buildFather() == father && idx >= 0){
					nbro++;
					nend++;
					if (markerfather < octants[idx].getMarker()+1){
						markerfather = octants[idx].getMarker()+1;
					}
					idx--;
					octants[idx].setMarker(-1);
					marker = octants[idx].getMarker();
					if (idx<0){
						break;
					}
				}
				if (nbro == global3D.nchildren){
					offset = nend;
				}
				else{
					nend = 0;
					for(int ii=nocts-global3D.nchildren; ii<nocts; ii++){
						octants[ii].setMarker(0);
						octants[ii].info[15] = true;
					}
				}
			}
			if (nend != 0){
				for (int iii=0; iii<16; iii++){
					father.info[iii] = false;
				}
				for (idx=0; idx < nend; idx++){
					for (int iii=0; iii<16; iii++){
						father.info[iii] = father.info[iii] || octants[nocts-idx-1].info[iii];
					}
				}
				father.info[13] = true;
				if (markerfather < 0){
					docoarse = true;
				}
				father.setMarker(markerfather);
				octants.resize(nocts-offset);
				octants.push_back(father);
				octants.shrink_to_fit();
				nocts = octants.size();
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

	// One level global refine with mapidx
	bool globalRefine(u32vector & mapidx){

		vector<uint32_t> last_child_index;
		Class_Octant<3>* children;
		uint32_t idx, nocts, ilastch;
		uint32_t offset = 0, blockidx;
		uint8_t nchm1 = global3D.nchildren-1, ich, iface;
		bool dorefine = false;

		nocts = octants.size();
		for (idx=0; idx<nocts; idx++){
			octants[idx].setMarker(1);
			if(octants[idx].getMarker() > 0 && octants[idx].getLevel() < MAX_LEVEL_3D){
				last_child_index.push_back(idx+nchm1+offset);
				offset += nchm1;
			}
			else{
				//			octants[idx].info[12] = false;
				if (octants[idx].marker > 0){
					octants[idx].marker = 0;
					octants[idx].info[15] = false;
				}
			}
		}
		if (offset > 0){
			mapidx.resize(octants.size()+offset);
			mapidx.shrink_to_fit();

			octants.resize(octants.size()+offset);
			blockidx = last_child_index[0]-nchm1;
			idx = octants.size();
			ilastch = last_child_index.size()-1;
			while (idx>blockidx){
				//			while (idx>0){
				idx--;
				//				if(octants[idx-offset].getMarker() > 0 && octants[idx-offset].getLevel() < MAX_LEVEL_3D){
				if(idx == last_child_index[ilastch]){
					children = octants[idx-offset].buildChildren();
					for (ich=0; ich<global3D.nchildren; ich++){
						octants[idx-ich] = (children[nchm1-ich]);
						mapidx[idx-ich]  = mapidx[idx-offset];
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
					delete []children;
					if (ilastch != 0){
						ilastch--;
					}
				}
				else {
					octants[idx] = octants[idx-offset];
					mapidx[idx]  = mapidx[idx-offset];
				}
			}
		}
		octants.shrink_to_fit();
		nocts = octants.size();

		setFirstDesc();
		setLastDesc();

		return dorefine;

	};

	// =================================================================================== //

	bool globalCoarse(u32vector & mapidx){							// Coarse local tree: coarse one time family of octants with marker <0
		// (if at least one octant of family has marker>=0 set marker=0 for the entire family)
		// mapidx[i] = index in old octants vector of the i-th octant (index of father if octant is new after)
		// Local variables
		vector<uint32_t> first_child_index;
		Class_Octant<3> father;
		uint32_t ich, nocts, nghosts, nocts0;
		int32_t idx, idx2;
		uint32_t offset;
		int32_t idx1_gh, idx2_gh;
		uint32_t nidx;
		int8_t markerfather, marker;
		uint8_t nbro, nstart, nend;
		uint8_t nchm1 = global3D.nchildren-1, iface;
		bool docoarse = false;

		//------------------------------------------ //
		// Initialization

		nbro = nstart = nend = 0;
		nidx = offset = 0;

		idx1_gh = idx2_gh = 0;

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
			idx2_gh = min(int(size_ghosts-1), idx2_gh);
		}

		// Check and coarse internal octants
		for (idx=0; idx<nocts; idx++){
			octants[idx].setMarker(-1);
			if(octants[idx].getMarker() < 0 && octants[idx].getLevel() > 0){
				nbro = 0;
				father = octants[idx].buildFather();
				// Check if family is to be refined
				for (idx2=idx; idx2<idx+global3D.nchildren; idx2++){
					if (idx2<nocts){
						octants[idx2].setMarker(-1);
						if(octants[idx2].getMarker() < 0 && octants[idx2].buildFather() == father){
							nbro++;
						}
					}
				}
				if (nbro == global3D.nchildren){
					nidx++;
					first_child_index.push_back(idx);
					idx = idx2-1;
				}
				else{
					if (idx < (nocts>global3D.nchildren)*(nocts-global3D.nchildren)){
						octants[idx].setMarker(0);
						octants[idx].info[15] = true;
					}
				}
			}
			//			else{
			//	//			octants[idx].info[13] = false;
			//			}
		}
		//TODO Da mettere dentro il primo ciclo per renderlo meno costoso
		if (nidx!=0){
			uint32_t nblock = nocts - nidx*nchm1 - nstart;
			nidx = 0;
			//for (idx=0; idx<nblock; idx++){
			for (idx=0; idx<nocts-offset; idx++){
				if (idx+offset == first_child_index[nidx]){
					markerfather = -MAX_LEVEL_3D;
					father = octants[idx+offset].buildFather();
					for (int iii=0; iii<16; iii++){
						father.info[iii] = false;
					}
					for(idx2=0; idx2<global3D.nchildren; idx2++){
						if (markerfather < octants[idx+offset+idx2].getMarker()+1){
							markerfather = octants[idx+offset+idx2].getMarker()+1;
						}
						for (int iii=0; iii<16; iii++){
							father.info[iii] = father.info[iii] || octants[idx+offset+idx2].info[iii];
						}
					}
					father.info[13] = true;
					father.setMarker(markerfather);
					if (markerfather < 0){
						docoarse = true;
					}
					octants[idx] = father;
					mapidx[idx] = mapidx[idx+offset];
					offset += nchm1;
					nidx++;
				}
				else{
					octants[idx] = octants[idx+offset];
					mapidx[idx] = mapidx[idx+offset];
				}
			}
		}
		octants.resize(nocts-offset);
		octants.shrink_to_fit();
		nocts = octants.size();
		mapidx.resize(nocts);
		mapidx.shrink_to_fit();


		// End on ghosts
		if (ghosts.size() && nocts > 0){
			ghosts[idx2_gh].setMarker(-1);
			if ((ghosts[idx2_gh].getMarker() < 0) && (octants[nocts-1].getMarker() < 0)){
				father = ghosts[idx2_gh].buildFather();
				markerfather = ghosts[idx2_gh].getMarker()+1;//-MAX_LEVEL_3D;
				nbro = 0;
				idx = idx2_gh;
				ghosts[idx].setMarker(-1);
				marker = ghosts[idx].getMarker();
				while(marker < 0 && ghosts[idx].buildFather() == father){
					nbro++;
					//TODO CAMBIATO IDX DA CAMBIARE ANCHE NELLE ALTRE COARSE!!!
					if (markerfather < ghosts[idx].getMarker()+1){
						markerfather = ghosts[idx].getMarker()+1;
					}
					idx++;
					if(idx == size_ghosts){
						break;
					}
					ghosts[idx].setMarker(-1);
					marker = ghosts[idx].getMarker();
				}
				nend = 0;
				idx = nocts-1;
				octants[idx].setMarker(-1);
				marker = octants[idx].getMarker();
				while(marker < 0 && octants[idx].buildFather() == father && idx >= 0){
					nbro++;
					nend++;
					if (markerfather < octants[idx].getMarker()+1){
						markerfather = octants[idx].getMarker()+1;
					}
					idx--;
					octants[idx].setMarker(-1);
					marker = octants[idx].getMarker();
					if (idx<0){
						break;
					}
				}
				if (nbro == global3D.nchildren){
					offset = nend;
				}
				else{
					nend = 0;
					for(int ii=nocts-global3D.nchildren; ii<nocts; ii++){
						octants[ii].setMarker(0);
						octants[ii].info[15] = true;
					}
				}
			}
			if (nend != 0){
				for (int iii=0; iii<16; iii++){
					father.info[iii] = false;
				}
				for (idx=0; idx < nend; idx++){
					for (int iii=0; iii<16; iii++){
						father.info[iii] = father.info[iii] || octants[nocts-idx-1].info[iii];
					}
				}
				father.info[13] = true;
				if (markerfather < 0){
					docoarse = true;
				}
				father.setMarker(markerfather);
				octants.resize(nocts-offset);
				octants.push_back(father);
				octants.shrink_to_fit();
				nocts = octants.size();
				mapidx.resize(nocts-offset);
				mapidx.push_back(nocts0-nend);
				mapidx.shrink_to_fit();
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

	void checkCoarse(uint64_t lastDescPre,						// Delete overlapping octants after coarse local tree. Check first and last descendants
								uint64_t firstDescPost){		// of process before and after the local process
		int32_t idx;
		uint32_t nocts;
		uint64_t Morton;
		uint8_t toDelete = 0;

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
		}
		octants.resize(nocts-toDelete);
		octants.shrink_to_fit();
		nocts = getNumOctants();

//		toDelete = 0;
//		Morton = last_desc.computeMorton();
//		if(int(firstDescPost  - Morton) > 1){
//			// To insert, the father is not yet here!!
//			idx = nocts - 1;
//			Class_Octant<3> father = octants[idx].buildFather();
//			while(octants[idx].buildFather() == father && idx >= 0){
//				toDelete++;
//				idx--;
//			}
//			father.info[global3D.nfaces+1] = father.info[global3D.nfaces+3] = true;
////			if(global3D.nfaces == 6)
//			father.info[global3D.nfaces+5] = true;
//			octants.resize(nocts-toDelete);
//			octants.push_back(father);
//			octants.shrink_to_fit();
//		}

		setFirstDesc();
		setLastDesc();

	};

	// =================================================================================== //

	void checkCoarse(uint64_t lastDescPre,						// Delete overlapping octants after coarse local tree. Check first and last descendants
					uint64_t firstDescPost,
					u32vector & mapidx){					// of process before and after the local process
		int32_t idx;
		uint32_t nocts;
		uint64_t Morton;
		uint8_t toDelete = 0;

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
		}
		octants.resize(nocts-toDelete);
		octants.shrink_to_fit();
		mapidx.resize(nocts-toDelete);
		mapidx.shrink_to_fit();
		nocts = getNumOctants();

		setFirstDesc();
		setLastDesc();

	};

	// =================================================================================== //

	void updateLocalMaxDepth(){						// Update max depth reached in local tree
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

	void findNeighbours(uint32_t idx,							// Finds neighbours of idx-th octant through iface in vector octants.
								uint8_t iface,					// Returns a vector (empty if iface is a bound face) with the index of neighbours
								u32vector & neighbours,			// in their structure (octants or ghosts) and sets isghost[i] = true if the
								vector<bool> & isghost){		// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		Class_Octant<3>* oct = &octants[idx];
		uint32_t size = oct->getSize();

		//Alternative to switch case
		int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
		int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
		int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

		isghost.clear();
		neighbours.clear();

		// Default if iface is nface<iface<0
		if (iface < 0 || iface > global3D.nfaces){
			writeLog("Face index out of range in find neighbours !!!");
			return;
		}

		// Check if octants face is a process boundary
		if (oct->info[global3D.nfaces+iface] == false){

			// Check if octants face is a boundary
			if (oct->info[iface] == false){

				//Build Morton number of virtual neigh of same size
				Class_Octant<3> samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
				Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
								return;
							}
						}
						while(octants[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > noctants-1){
								idxtry = 0;
								return;
							}
						}
					}
					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(false);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
					Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					vector<uint32_t> bufferidx;
					Mortontry = octants[idxtry].computeMorton();
					int32_t Dh;
					int32_t eqcoord;
					while(Mortontry < Mortonlast && idxtry < noctants){
						Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
						Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
						Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
						if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
							neighbours.push_back(idxtry);
							isghost.push_back(false);
						}
						idxtry++;
						if(idxtry>noctants-1){
							break;
						}
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

				if (ghosts.size()>0){
					// Search in ghosts

					uint32_t idxghost = uint32_t(size_ghosts/2);
					Class_Octant<3>* octghost = &ghosts[idxghost];

					//Build Morton number of virtual neigh of same size
					Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
					Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
					// Search morton in octants
					// If a even face morton is lower than morton of oct, if odd higher
					// ---> can i search only before or after idx in octants
					int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
					idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
							uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
							Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
							uint64_t Mortonlast = last_desc.computeMorton();
							vector<uint32_t> bufferidx;
							Mortontry = ghosts[idxtry].computeMorton();
							int32_t Dh;
							int32_t eqcoord;
							while(Mortontry < Mortonlast && idxtry < size_ghosts){
								Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
								Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
								Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
								if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
								idxtry++;
								if(idxtry>size_ghosts-1){
									break;
								}
								Mortontry = ghosts[idxtry].computeMorton();
							}
						}
					}

					uint32_t lengthneigh = 0;
					uint32_t sizeneigh = neighbours.size();
					for (idxtry=0; idxtry<sizeneigh; idxtry++){
						lengthneigh += ghosts[neighbours[idxtry]].getSize();
					}
					if (lengthneigh < oct->getSize()){
						// Search in octants

						// Check if octants face is a boundary
						if (oct->info[iface] == false){

							//Build Morton number of virtual neigh of same size
							Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
							Morton = samesizeoct.computeMorton();
							// Search morton in octants
							// If a even face morton is lower than morton of oct, if odd higher
							// ---> can i search only before or after idx in octants
							int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
							idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
										writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
										return;
									}
									// Compute Last discendent of virtual octant of same size
									uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
									Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
									uint64_t Mortonlast = last_desc.computeMorton();
									vector<uint32_t> bufferidx;
									Mortontry = octants[idxtry].computeMorton();
									int32_t Dh;
									int32_t eqcoord;
									while(Mortontry < Mortonlast && idxtry < noctants-1){
										Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
										Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
										Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
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
//
	// =================================================================================== //

	void findNeighbours(Class_Octant<3>* oct,					// Finds neighbours of octant through iface in vector octants.
			uint8_t iface,					// Returns a vector (empty if iface is a bound face) with the index of neighbours
			u32vector & neighbours,			// in their structure (octants or ghosts) and sets isghost[i] = true if the
			vector<bool> & isghost){		// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		uint32_t size = oct->getSize();
		//Find octant in octants
		OctantsType::iterator itoct = find(octants.begin(), octants.end(), (*oct));
		if (itoct == octants.end()){
			return;
		}
		uint32_t idx = distance(octants.begin(), itoct);

		// TODO Create a global matrix
		//Alternative to switch case
		int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
		int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
		int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

		isghost.clear();
		neighbours.clear();

		// Default if iface is nface<iface<0
		if (iface < 0 || iface > global3D.nfaces){
			writeLog("Face index out of range in find neighbours !!!");
			return;
		}

		// Check if octants face is a process boundary
		if (oct->info[global3D.nfaces+iface] == false){

			// Check if octants face is a boundary
			if (oct->info[iface] == false){

				//Build Morton number of virtual neigh of same size
				Class_Octant<3> samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
				Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
						//Found neighbour of same size
						isghost.push_back(false);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
					Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					vector<uint32_t> bufferidx;
					Mortontry = octants[idxtry].computeMorton();
					int32_t Dh;
					int32_t eqcoord;
					while(Mortontry < Mortonlast && idxtry < noctants){
						Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
						Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
						Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
						if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
							neighbours.push_back(idxtry);
							isghost.push_back(false);
						}
						idxtry++;
						if(idxtry>noctants-1){
							break;
						}
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

				if (ghosts.size()>0){
					// Search in ghosts

					uint32_t idxghost = uint32_t(size_ghosts/2);
					Class_Octant<3>* octghost = &ghosts[idxghost];

					//Build Morton number of virtual neigh of same size
					Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
					Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);
					// Search morton in octants
					// If a even face morton is lower than morton of oct, if odd higher
					// ---> can i search only before or after idx in octants
					int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
					idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
							uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
							Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
							uint64_t Mortonlast = last_desc.computeMorton();
							vector<uint32_t> bufferidx;
							Mortontry = ghosts[idxtry].computeMorton();
							int32_t Dh;
							int32_t eqcoord;
							while(Mortontry < Mortonlast & idxtry < size_ghosts){
								Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
								Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
								Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
								if ((abs(Dh) == ((1-(iface%2))*ghosts[idxtry].getSize() + (iface%2)*size))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
								idxtry++;
								if(idxtry>size_ghosts-1){
									break;
								}
								Mortontry = ghosts[idxtry].computeMorton();
							}
						}
					}

					uint32_t lengthneigh = 0;
					uint32_t sizeneigh = neighbours.size();
					for (idxtry=0; idxtry<sizeneigh; idxtry++){
						lengthneigh += ghosts[neighbours[idxtry]].getSize();
					}
					if (lengthneigh < oct->getSize()){
						// Search in octants

						// Check if octants face is a boundary
						if (oct->info[iface] == false){

							//Build Morton number of virtual neigh of same size
							Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
							Morton = samesizeoct.computeMorton();
							// Search morton in octants
							// If a even face morton is lower than morton of oct, if odd higher
							// ---> can i search only before or after idx in octants
							int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
							idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
										writeLog("Face marked pbound but only a non-ghost neighbour found!!!");
										return;
									}
									// Compute Last discendent of virtual octant of same size
									uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
									Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
									uint64_t Mortonlast = last_desc.computeMorton();
									vector<uint32_t> bufferidx;
									Mortontry = octants[idxtry].computeMorton();
									int32_t Dh;
									int32_t eqcoord;
									while(Mortontry < Mortonlast & idxtry < noctants-1){
										Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
										Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
										Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
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

	void findGhostNeighbours(uint32_t const idx,		// Finds neighbours of idx-th ghost octant through iface in vector octants.
								uint8_t iface,					// Returns a vector (empty if iface is not the pbound face for ghost) with the index of neighbours
								u32vector & neighbours){		// in the structure octants

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		Class_Octant<3>* oct = &ghosts[idx];
		uint32_t size = oct->getSize();

		//Alternative to switch case
		int8_t cx = int8_t((iface<2)*(int8_t(2*iface-1)));
		int8_t cy = int8_t((iface<4)*(int8_t(iface/2))*(int8_t(2*iface-5)));
		int8_t cz = int8_t((int8_t(iface/4))*(int8_t(2*iface-9)));

		neighbours.clear();

		// Default if iface is nface<iface<0
		if (iface < 0 || iface > global3D.nfaces){
			writeLog("Face index out of range in find neighbours !!!");
			return;
		}

		// Check if octants face is a process boundary
		if (oct->info[global3D.nfaces+iface] == true){

				//Build Morton number of virtual neigh of same size
				Class_Octant<3> samesizeoct(oct->level, int32_t(oct->x)+int32_t(cx*size), int32_t(oct->y)+int32_t(cy*size), int32_t(oct->z)+int32_t(cz*size));
				Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump;
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
								return;
							}
						}
						while(octants[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > noctants-1){
								idxtry = 0;
								return;
							}
						}
					}
					if(octants[idxtry].computeMorton() == Morton && octants[idxtry].level == oct->level){
						//Found neighbour of same size
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
					Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					vector<uint32_t> bufferidx;
					Mortontry = octants[idxtry].computeMorton();
					int32_t Dh;
					int32_t eqcoord;
					while(Mortontry < Mortonlast && idxtry < noctants){
						Dh = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
						Dh += int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
						Dh += int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
						if ((abs(Dh) == ((1-(iface%2))*octants[idxtry].getSize() + (iface%2)*size))){
							neighbours.push_back(idxtry);
						}
						idxtry++;
						if(idxtry>noctants-1){
							break;
						}
						Mortontry = octants[idxtry].computeMorton();
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

	bool localBalance(bool doInterior){				// 2:1 balancing on level a local tree already adapted (balance only the octants with info[14] = false) (refinement wins!)
																// Return true if balanced done with some markers modification
																// Seto doInterior = false if the interior octants are already balanced
		// Local variables
		uint32_t 			noctants = getNumOctants();
		uint32_t			sizeneigh, modsize;
		u32vector		 	neigh;
		u32vector		 	modified, newmodified;
		uint32_t 			i, idx, imod;
		uint8_t				iface;
		int8_t				targetmarker;
		vector<bool> 		isghost;
		bool				Bdone = false;

		OctantsType::iterator 	obegin, oend, it;
		u32vector::iterator 	ibegin, iend, iit;

		//If interior octants have to be balanced
		if(doInterior){
			// First loop on the octants
			/*		for(idx=0 ; idx<noctants; idx++){
				if (!octants[idx].getNotBalance()){
					for (iface=1; iface<nface; iface+=2){
						if(!octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) > (octants[idx].getLevel() + octants[idx].getMarker() + 1) ){
											octants[idx].setMarker(octants[neigh[i]].getLevel()+octants[neigh[i]].getMarker()-1-octants[idx].getLevel());
											modified.push_back(idx);
											Bdone = true;
										}
										else if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (octants[idx].getLevel() + octants[idx].getMarker() - 1)){
											octants[neigh[i]].setMarker(octants[idx].getLevel()+octants[idx].getMarker()-octants[neigh[i]].getLevel()-1);
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}

							else{
								if(ghosts.size()>0){
									if((ghosts[neigh[i]].getLevel() + ghosts[neigh[i]].getMarker())> (octants[idx].getLevel() + octants[idx].getMarker() + 1)){
										octants[idx].setMarker(ghosts[neigh[i]].getLevel()+ghosts[neigh[i]].getMarker()-octants[idx].getLevel()-1);
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}

							}
						}
					}
				}
			}*/
			obegin = octants.begin();
			oend = octants.end();
			idx = 0;
			for (it=obegin; it!=oend; it++){
				it->info[15] = false;
				if (!it->getNotBalance() && it->getMarker() != 0){
					targetmarker = min(MAX_LEVEL_3D, (octants[idx].getLevel() + octants[idx].getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
						if(!it->getPbound(iface)){
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
											octants[idx].info[11] = true;
											modified.push_back(idx);
											Bdone = true;
										}
									};
								}
							}
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
					targetmarker = min(MAX_LEVEL_3D, (it->getLevel()+it->getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
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
						targetmarker = min(MAX_LEVEL_3D, (octants[idx].getLevel()+octants[idx].getMarker()));
						for (iface=0; iface<global3D.nfaces; iface++){
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
					}
				}
				u32vector().swap(modified);
				swap(modified,newmodified);
				modsize = modified.size();
				u32vector().swap(newmodified);
			}// end while

		}
		else{

			// Loop on ghost octants (influence over interior borders)
			/*	for (idx=0; idx<size_ghosts; idx++){
			if (!ghosts[idx].getNotBalance()){
				for (iface=0; iface<nface; iface++){
					if(ghosts[idx].getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((octants[neigh[i]].getLevel() + octants[neigh[i]].getMarker()) < (ghosts[idx].getLevel() + ghosts[idx].getMarker() - 1)){
								octants[neigh[i]].setMarker(ghosts[idx].getLevel()+ghosts[idx].getMarker()-octants[neigh[i]].getLevel()-1);
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
				}
			}
		}*/
			obegin = ghosts.begin();
			oend = ghosts.end();
			idx = 0;
			for (it=obegin; it!=oend; it++){
				if (!it->getNotBalance() && it->info[15]){
					targetmarker = min(MAX_LEVEL_3D, (it->getLevel()+it->getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
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
						targetmarker = min(MAX_LEVEL_3D, (octants[idx].getLevel()+octants[idx].getMarker()));
						for (iface=0; iface<global3D.nfaces; iface++){
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
					}
				}
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

	bool localBalanceAll(bool doInterior){				// 2:1 balancing on level a local tree already adapted (balance only the octants with info[14] = false) (refinement wins!)
																// Return true if balanced done with some markers modification
																// Seto doInterior = false if the interior octants are already balanced
		// Local variables
		uint32_t 			noctants = getNumOctants();
		uint32_t			sizeneigh, modsize;
		u32vector		 	neigh;
		u32vector		 	modified, newmodified;
		uint32_t 			i, idx, imod;
		int32_t				idx1_gh = 0, idx2_gh = 0;
		uint8_t				iface;
		int8_t				targetmarker;
		vector<bool> 		isghost;
		bool				Bdone = false;

		OctantsType::iterator 	obegin, oend, it;
		u32vector::iterator 	ibegin, iend, iit;


		//If interior octants have to be balanced
		if(doInterior){
			// First loop on the octants
			obegin = octants.begin();
			oend = octants.end();
			idx = 0;
			for (it=obegin; it!=oend; it++){
				if ((!it->getNotBalance()) && ((it->info[15]) || (it->getMarker()!=0) || ((it->getIsNewC()) || (it->getIsNewR())))){
					targetmarker = min(MAX_LEVEL_3D, int(octants[idx].getLevel()) + int(octants[idx].getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
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
				}
				idx++;
			}
			// Loop on ghost octants (influence over interior borders)
			obegin = ghosts.begin();
			oend = ghosts.end();
			idx = 0;
			for (it=obegin; it!=oend; it++){
				if (!it->getNotBalance() && (it->info[15] || (it->getIsNewC() || it->getIsNewR()))){
					targetmarker = min(MAX_LEVEL_3D, (it->getLevel()+it->getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
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
						targetmarker = min(MAX_LEVEL_3D, (octants[idx].getLevel()+octants[idx].getMarker()));
						for (iface=0; iface<global3D.nfaces; iface++){
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
					}
				}
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
					targetmarker = min(MAX_LEVEL_3D, (it->getLevel()+it->getMarker()));
					for (iface=0; iface<global3D.nfaces; iface++){
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
						targetmarker = min(MAX_LEVEL_3D, (octants[idx].getLevel()+octants[idx].getMarker()));
						for (iface=0; iface<global3D.nfaces; iface++){
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
					}
				}
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

	void			findEdgeNeighbours(uint32_t idx,			// Finds neighbours of idx-th octant through iedge in vector octants.
									uint8_t iedge,				// Returns a vector (empty if iedge is a bound edge) with the index of neighbours
									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
									vector<bool> & isghost){	// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		Class_Octant<3>* oct = &octants[idx];
		uint32_t size = oct->getSize();
		uint8_t iface1, iface2;
		int32_t Dhx, Dhy, Dhz;
		int32_t Dhxref, Dhyref, Dhzref;

		//Alternative to switch case
		int8_t Cx[12] = {-1,1,0,0,-1,1,-1,1,-1,1,0,0};
		int8_t Cy[12] = {0,0,-1,1,-1,-1,1,1,0,0,-1,1};
		int8_t Cz[12] = {-1,-1,-1,-1,0,0,0,0,1,1,1,1};
		int8_t cx = Cx[iedge];
		int8_t cy = Cy[iedge];
		int8_t cz = Cz[iedge];

		isghost.clear();
		neighbours.clear();

		// Default if iedge is nface<iedge<0
		if (iedge < 0 || iedge > global3D.nfaces*2){
			writeLog("Edge index out of range in find neighbours !!!");
			return;
		}

		// Check if octants edge is a process boundary
		iface1 = global3D.edgeface[iedge][0];
		iface2 = global3D.edgeface[iedge][1];

		// Check if octants edge is a boundary
		if (oct->info[iface1] == false && oct->info[iface2] == false){

			//Build Morton number of virtual neigh of same size
			Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

			//SEARCH IN GHOSTS

			if (ghosts.size()>0){
				// Search in ghosts

				uint32_t idxghost = uint32_t(size_ghosts/2);
				Class_Octant<3>* octghost = &ghosts[idxghost];

				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
							if(idxtry < 0){
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = ghosts[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < ghosts.size()){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
							Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(true);
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
			uint32_t lengthneigh = 0;
//			uint32_t sizeneigh = neighbours.size();
//			for (idxtry=0; idxtry<sizeneigh; idxtry++){
//				lengthneigh += ghosts[neighbours[idxtry]].getSize();
//			}
			if (lengthneigh < oct->getSize()){

				// Search in octants

				//Build Morton number of virtual neigh of same size
				//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				//Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < noctants-1){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
							Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
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
		}
		else{
			// Boundary Face
			return;
		}

	};

	// =================================================================================== //

	void			findEdgeNeighbours(Class_Octant<3>* oct,			// Finds neighbours of idx-th octant through iedge in vector octants.
									uint8_t iedge,				// Returns a vector (empty if iedge is a bound edge) with the index of neighbours
									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
									vector<bool> & isghost){	// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		//Class_Octant<3>* oct = &octants[idx];
		uint32_t size = oct->getSize();
		uint8_t iface1, iface2;
		int32_t Dhx, Dhy, Dhz;
		int32_t Dhxref, Dhyref, Dhzref;

		//Alternative to switch case
		int8_t Cx[12] = {-1,1,0,0,-1,1,-1,1,-1,1,0,0};
		int8_t Cy[12] = {0,0,-1,1,-1,-1,1,1,0,0,-1,1};
		int8_t Cz[12] = {-1,-1,-1,-1,0,0,0,0,1,1,1,1};
		int8_t cx = Cx[iedge];
		int8_t cy = Cy[iedge];
		int8_t cz = Cz[iedge];

		isghost.clear();
		neighbours.clear();

		// Default if iedge is nface<iedge<0
		if (iedge < 0 || iedge > global3D.nfaces*2){
			writeLog("Edge index out of range in find neighbours !!!");
			return;
		}

		//Find octant in octants
		OctantsType::iterator itoct = find(octants.begin(), octants.end(), (*oct));
		if (itoct == octants.end()){
			return;
		}
		uint32_t idx = distance(octants.begin(), itoct);

		// Check if octants edge is a process boundary
		iface1 = global3D.edgeface[iedge][0];
		iface2 = global3D.edgeface[iedge][1];

		// Check if octants edge is a boundary
		if (oct->info[iface1] == false && oct->info[iface2] == false){

			//Build Morton number of virtual neigh of same size
			Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

			//SEARCH IN GHOSTS

			if (ghosts.size()>0){
				// Search in ghosts

				uint32_t idxghost = uint32_t(size_ghosts/2);
				Class_Octant<3>* octghost = &ghosts[idxghost];

				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
							if(idxtry < 0){
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = ghosts[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < ghosts.size()){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
							Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(true);
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
			uint32_t lengthneigh = 0;
//			uint32_t sizeneigh = neighbours.size();
//			for (idxtry=0; idxtry<sizeneigh; idxtry++){
//				lengthneigh += ghosts[neighbours[idxtry]].getSize();
//			}
			if (lengthneigh < oct->getSize()){

				// Search in octants

				//Build Morton number of virtual neigh of same size
				//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				//Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < noctants-1){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
							Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
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
		}
		else{
			// Boundary Face
			return;
		}

	};

	// =================================================================================== //

	void			findNodeNeighbours(uint32_t idx,			// Finds neighbours of idx-th octant through inode in vector octants.
									uint8_t inode,				// Returns a vector (empty if inode is a bound node) with the index of neighbours
									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
									vector<bool> & isghost){	// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
		Class_Octant<3>* oct = &octants[idx];
		uint32_t size = oct->getSize();
		uint8_t iface1, iface2, iface3;
		int32_t Dhx, Dhy, Dhz;
		int32_t Dhxref, Dhyref, Dhzref;

		//Alternative to switch case
		int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
		int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
		int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
		int8_t cx = Cx[inode];
		int8_t cy = Cy[inode];
		int8_t cz = Cz[inode];

		isghost.clear();
		neighbours.clear();

		// Default if inode is nnodes<inode<0
		if (inode < 0 || inode > global3D.nnodes){
			writeLog("Node index out of range in find neighbours !!!");
			return;
		}

		// Check if octants node is a boundary
		iface1 = global3D.nodeface[inode][0];
		iface2 = global3D.nodeface[inode][1];
		iface3 = global3D.nodeface[inode][2];

		// Check if octants node is a boundary
		if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){

			//Build Morton number of virtual neigh of same size
			Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

			//SEARCH IN GHOSTS

			if (ghosts.size()>0){
				// Search in ghosts

				uint32_t idxghost = uint32_t(size_ghosts/2);
				Class_Octant<3>* octghost = &ghosts[idxghost];

				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = ghosts[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < size_ghosts){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
							Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(true);
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
			uint32_t lengthneigh = 0;
//			uint32_t sizeneigh = neighbours.size();
//			for (idxtry=0; idxtry<sizeneigh; idxtry++){
//				lengthneigh += ghosts[neighbours[idxtry]].getSize();
//			}
			if (lengthneigh < oct->getSize()){

				// Search in octants

				//Build Morton number of virtual neigh of same size
				//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				//Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
				if (idxtry > noctants-1)
					idxtry = noctants-1;
				while(abs(jump) > 0){
					Mortontry = octants[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < noctants-1){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
							Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
							idxtry++;
							Mortontry = octants[idxtry].computeMorton();
						}
					}
				}
				return;
			}
		}
		else{
			// Boundary Node
			return;
		}

	};

	// =================================================================================== //

	void			findNodeNeighbours(Class_Octant<3>* oct,			// Finds neighbours of idx-th octant through inode in vector octants.
									uint8_t inode,				// Returns a vector (empty if inode is a bound node) with the index of neighbours
									u32vector & neighbours,		// in their structure (octants or ghosts) and sets isghost[i] = true if the
									vector<bool> & isghost){	// i-th neighbour is ghost in the local tree

		uint64_t  Morton, Mortontry;
		uint32_t  noctants = getNumOctants();
		uint32_t idxtry;
//		Class_Octant<3>* oct = &octants[idx];
		uint32_t size = oct->getSize();
		uint8_t iface1, iface2, iface3;
		int32_t Dhx, Dhy, Dhz;
		int32_t Dhxref, Dhyref, Dhzref;

		//Alternative to switch case
		int8_t Cx[8] = {-1,1,-1,1,-1,1,-1,1};
		int8_t Cy[8] = {-1,-1,1,1,-1,-1,1,1};
		int8_t Cz[8] = {-1,-1,-1,-1,1,1,1,1};
		int8_t cx = Cx[inode];
		int8_t cy = Cy[inode];
		int8_t cz = Cz[inode];

		isghost.clear();
		neighbours.clear();

		// Default if inode is nnodes<inode<0
		if (inode < 0 || inode > global3D.nnodes){
			writeLog("Node index out of range in find neighbours !!!");
			return;
		}

		//Find octant in octants
		OctantsType::iterator itoct = find(octants.begin(), octants.end(), (*oct));
		if (itoct == octants.end()){
			return;
		}
		uint32_t idx = distance(octants.begin(), itoct);

		// Check if octants node is a boundary
		iface1 = global3D.nodeface[inode][0];
		iface2 = global3D.nodeface[inode][1];
		iface3 = global3D.nodeface[inode][2];

		// Check if octants node is a boundary
		if (oct->info[iface1] == false && oct->info[iface2] == false && oct->info[iface3] == false){

			//Build Morton number of virtual neigh of same size
			Class_Octant<3> samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
			Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->x-size,oct->y,oct->z);

			//SEARCH IN GHOSTS

			if (ghosts.size()>0){
				// Search in ghosts

				uint32_t idxghost = uint32_t(size_ghosts/2);
				Class_Octant<3>* octghost = &ghosts[idxghost];

				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (octghost->computeMorton() > Morton) ? int32_t(idxghost/2+1) : int32_t((size_ghosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = ghosts[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < size_ghosts){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(ghosts[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(ghosts[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(ghosts[idxtry].z));
							Dhxref = int32_t(cx<0)*ghosts[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*ghosts[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*ghosts[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(true);
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
			uint32_t lengthneigh = 0;
//			uint32_t sizeneigh = neighbours.size();
//			for (idxtry=0; idxtry<sizeneigh; idxtry++){
//				lengthneigh += ghosts[neighbours[idxtry]].getSize();
//			}
			if (lengthneigh < oct->getSize()){

				// Search in octants

				//Build Morton number of virtual neigh of same size
				//Class_Octant samesizeoct(oct->level, oct->x+cx*size, oct->y+cy*size, oct->z+cz*size);
				//Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				int32_t jump = (oct->computeMorton() > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
				idxtry = uint32_t(idx +((oct->computeMorton()<Morton)-(oct->computeMorton()>Morton))*jump);
				if (idxtry > noctants-1)
					idxtry = noctants-1;
				while(abs(jump) > 0){
					Mortontry = octants[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
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
						uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL_3D - samesizeoct.level)) - 1;
						Class_Octant<3> last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						vector<uint32_t> bufferidx;
						Mortontry = octants[idxtry].computeMorton();
						while(Mortontry < Mortonlast && idxtry < noctants-1){
							Dhx = int32_t(cx)*(int32_t(oct->x) - int32_t(octants[idxtry].x));
							Dhy = int32_t(cy)*(int32_t(oct->y) - int32_t(octants[idxtry].y));
							Dhz = int32_t(cz)*(int32_t(oct->z) - int32_t(octants[idxtry].z));
							Dhxref = int32_t(cx<0)*octants[idxtry].getSize() + int32_t(cx>0)*size;
							Dhyref = int32_t(cy<0)*octants[idxtry].getSize() + int32_t(cy>0)*size;
							Dhzref = int32_t(cz<0)*octants[idxtry].getSize() + int32_t(cz>0)*size;
							if ((abs(Dhx) == Dhxref) && (abs(Dhy) == Dhyref) && (abs(Dhz) == Dhzref)){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
							idxtry++;
							Mortontry = octants[idxtry].computeMorton();
						}
					}
				}
				return;
			}
		}
		else{
			// Boundary Node
			return;
		}

	};


	// =================================================================================== //

	void computeIntersections() {

		OctantsType::iterator it, obegin, oend;
		Class_Intersection<3> intersection;
		u32vector neighbours;
		vector<bool> isghost;
		uint32_t counter, idx;
		uint32_t i, j, k, nsize;
		uint8_t iface, iface2;

		intersections.clear();
		intersections.reserve(2*3*octants.size());

		counter = idx = 0;

		// Loop on ghosts
		obegin = ghosts.begin();
		oend = ghosts.end();
		for (it = obegin; it != oend; it++){
			for (iface = 0; iface < 3; iface++){
				iface2 = iface*2;
				findGhostNeighbours(idx, iface2, neighbours);
				nsize = neighbours.size();
				for (i = 0; i < nsize; i++){
					intersection.finer = (nsize==1);
					intersection.owners[0]  = neighbours[i];
					intersection.owners[1] = idx;
					intersection.iface = global3D.oppface[iface2] - (nsize==1);
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
			for (iface = 0; iface < 3; iface++){
				iface2 = iface*2;
				findNeighbours(idx, iface2, neighbours, isghost);
				nsize = neighbours.size();
				if (nsize) {
					for (i = 0; i < nsize; i++){
						if (isghost[i]){
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
		intersections.shrink_to_fit();
	}

	// =================================================================================== //

	//TODO Update intersections killed
//	void updateIntersections(u32vector & mapidx,
//							u32vector & mapinters_int,
//							u32vector & mapinters_ghost,
//							u32vector & mapinters_bord) {
//
//		map<uint32_t, uint32_t> invmapidx;
//		vector<uint32_t> newocts;
//		OctantsType::iterator it, obegin, oend;
//		Class_Intersection<3> intersection;
//		u32vector neighbours;
//		vector<bool> isghost;
//		uint32_t counter_g, idx;
//		uint32_t i, j, nsize, msize, isize, osize, offset;
//		uint8_t iface, iface2;
//
//
//		if (intersections_int.size()==0){
//			computeIntersections();
//			mapinters_int.clear();
//			mapinters_ghost.clear();
//			mapinters_bord.clear();
//		}
//		else{
//			msize = mapidx.size();
//			for (i=0; i<msize; i++){
//				invmapidx[mapidx[i]] = i;
//			}
//
//			//Internal Intersections
//			isize = intersections_int.size();
//			mapinters_int.clear();
//			mapinters_int.resize(isize);
//			offset = 0;
//			newocts.clear();
//			for (i=0; i<isize-offset; i++){
//				if (octants[invmapidx[intersections_int[i].owners[0]]].info[12] ||
//						octants[invmapidx[intersections_int[i].owners[0]]].info[13]){
//					offset++;
//					newocts.push_back(invmapidx[intersections_int[i].owners[0]]);
//				}
//				intersections_int[i] = intersections_int[i+offset];
//				intersections_int[i].owners[0] = invmapidx[intersections_int[i+offset].owners[0]];
//				intersections_int[i].owners[1] = invmapidx[intersections_int[i+offset].owners[1]];
//				mapinters_int[i] = i+offset;
//			}
//			intersections_int.resize(isize-offset);
//			isize = isize-offset;
//			mapinters_int.resize(isize);
//			osize = newocts.size();
//			for (j=0; j<osize; j++){
//				idx = newocts[j];
//				for (iface = 0; iface < 3; iface++){
//					iface2 = iface*2;
//					findNeighbours(idx, iface2, neighbours, isghost);
//					nsize = neighbours.size();
//					if (nsize) {
//						for (i = 0; i < nsize; i++){
//							if (isghost[i]){
//								intersection.owners[0] = idx;
//								intersection.owners[1] = neighbours[i];
//								intersection.finer = (i>=1);
//								intersection.iface = iface2;
//								intersection.isnew = true;
//								intersection.isghost = false;
//								intersections_int.push_back(intersection);
//							}
//						}
//					}
//				}
//			}
//
//			//Border Intersections
//			isize = intersections_bord.size();
//			mapinters_bord.clear();
//			mapinters_bord.resize(isize);
//			offset = 0;
//			newocts.clear();
//			for (i=0; i<isize-offset; i++){
//				if (octants[invmapidx[intersections_bord[i].owners[0]]].info[12] ||
//						octants[invmapidx[intersections_bord[i].owners[0]]].info[13]){
//					offset++;
//					newocts.push_back(invmapidx[intersections_bord[i].owners[0]]);
//				}
//				intersections_bord[i] = intersections_bord[i+offset];
//				intersections_bord[i].owners[0] = invmapidx[intersections_bord[i+offset].owners[0]];
//				intersections_bord[i].owners[1] = invmapidx[intersections_bord[i+offset].owners[1]];
//				mapinters_bord[i] = i+offset;
//			}
//			intersections_bord.resize(isize-offset);
//			isize = isize-offset;
//			mapinters_bord.resize(isize);
//			osize = newocts.size();
//			for (j=0; j<osize; j++){
//				idx = newocts[j];
//				for (iface = 0; iface < 3; iface++){
//					if (octants[idx].info[iface]){
//						intersection.owners[0] = idx;
//						intersection.owners[1] = idx;
//						intersection.finer = 0;
//						intersection.iface = iface;
//						intersection.isnew = true;
//						intersection.isghost = false;
//						intersections_bord.push_back(intersection);
//					}
//				}
//			}
//
//			// TODO GHOSTS INTERSECTIONS. NOW ALL ARE SET AS NEW
//			//Ghost Intersections
//			mapinters_ghost.clear();
//			intersections_ghost.clear();
//			intersections_ghost.reserve(2*3*ghosts.size());
//
//			counter_g = idx = 0;
//
//			// Loop on ghosts
//			obegin = ghosts.begin();
//			oend = ghosts.end();
//			for (it = obegin; it != oend; it++){
//				for (iface = 0; iface < 2*3; iface++){
//					findGhostNeighbours(idx, iface, neighbours);
//					nsize = neighbours.size();
//					for (i = 0; i < nsize; i++){
//						intersection.finer = (i<1);
//						intersection.owners[0]  = neighbours[i];
//						intersection.owners[1] = idx;
//						intersection.iface = global3D.oppface[iface];
//						if (octants[idx].info[12] ||
//								octants[idx].info[13] ||
//								it->info[12] || it->info[13]){
//							intersection.isnew = true;
//						}
//						else{
//							intersection.isnew = false;
//						}
//						intersection.isghost = true;
//						intersections_ghost.push_back(intersection);
//						counter_g++;
//					}
//				}
//				idx++;
//			}
//
//			intersections_int.shrink_to_fit();
//			intersections_ghost.shrink_to_fit();
//			intersections_bord.shrink_to_fit();
//		}
//
//	}

	// =================================================================================== //

	uint32_t findMorton(uint64_t Morton){				// Find an input Morton in octants and return the local idx
		uint32_t nocts = octants.size();
		uint32_t idx = nocts/2;
		uint64_t Mortontry = octants[idx].computeMorton();
		int32_t jump = nocts/2;
		while(abs(jump)>0){
			if (Mortontry == Morton){
				return idx;
			}
			Mortontry = octants[idx].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idx += jump;
			if (idx > nocts){
				return nocts-1;   // return nocts if not found the Morton
			}
		}
		if (Mortontry<Morton){
			for (int idx2=idx; idx2<nocts; idx2++){
				Mortontry = octants[idx2].computeMorton();
				if (Mortontry == Morton){
					return idx2;
				}
			}
		}
		else{
			for(int idx2=0; idx2<idx+1; idx2++){
				Mortontry = octants[idx2].computeMorton();
				if (Mortontry == Morton){
					return idx2;
				}
			}
		}
		return nocts-1;
	};

	// =================================================================================== //

	uint32_t findGhostMorton(uint64_t Morton){			// Find an input Morton in ghosts and return the local idx
		uint32_t nocts = ghosts.size();
		uint32_t idx = nocts/2;
		uint64_t Mortontry = ghosts[idx].computeMorton();
		int32_t jump = nocts/2;
		while(abs(jump)>0){
			if (Mortontry == Morton){
				return idx;
			}
			Mortontry = ghosts[idx].computeMorton();
			jump = (Mortontry<Morton)*jump/4;
			idx += jump;
			if (idx > nocts){
				return nocts;   // return nocts if not found the Morton
			}
		}
		if (Mortontry<Morton){
			for (int idx2=idx; idx2<nocts; idx2++){
				Mortontry = ghosts[idx2].computeMorton();
				if (Mortontry == Morton){
					return idx2;
				}
			}
		}
		else{
			for(int idx2=0; idx2<idx; idx2++){
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
	void computeConnectivity(){						// Computes nodes vector and connectivity of octants of local tree
		map<uint64_t, vector<uint32_t> > mapnodes;
		map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
		uint32_t i, k, counter;
		uint64_t morton;
		uint32_t noctants = getNumOctants();
		u32vector2D octnodes;
//		uint32_t (*octnodes)[3];
		uint8_t j;

		clearConnectivity();
		octnodes.reserve(global3D.nnodes);

		if (nodes.size() == 0){
			connectivity.resize(noctants);
			for (i = 0; i < noctants; i++){
				octants[i].getNodes(octnodes);
				for (j = 0; j < global3D.nnodes; j++){
					morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
					if (mapnodes[morton].size()==0){
						mapnodes[morton].reserve(12);
						for (k = 0; k < 3; k++){
							mapnodes[morton].push_back(octnodes[j][k]);
						}
					}
					mapnodes[morton].push_back(double(i));
				}
				u32vector2D().swap(octnodes);
			}
			iter	= mapnodes.begin();
			iterend	= mapnodes.end();
			counter = 0;
			uint32_t numnodes = mapnodes.size();
			nodes.resize(numnodes);
			while (iter != iterend){
				vector<uint32_t> nodecasting(iter->second.begin(), iter->second.begin()+3);
				nodes[counter] = nodecasting;
				nodes[counter].shrink_to_fit();
				for(vector<uint32_t>::iterator iter2 = iter->second.begin()+3; iter2 != iter->second.end(); iter2++){
					if (connectivity[(*iter2)].size()==0){
						connectivity[(*iter2)].reserve(8);
					}
					connectivity[(*iter2)].push_back(counter);
				}
				mapnodes.erase(iter++);
				counter++;
			}
			nodes.shrink_to_fit();
			//Lento. Solo per risparmiare memoria
			for (int ii=0; ii<noctants; ii++){
				connectivity[ii].shrink_to_fit();
			}
			connectivity.shrink_to_fit();
		}
		map<uint64_t, vector<uint32_t> >().swap(mapnodes);
		iter = mapnodes.end();

	};

	// =================================================================================== //

	void clearConnectivity(){						// Clear nodes vector and connectivity of octants of local tree
		u32vector2D().swap(nodes);
		u32vector2D().swap(connectivity);
	};

	// =================================================================================== //

	void updateConnectivity(){						// Updates nodes vector and connectivity of octants of local tree
		clearConnectivity();
		computeConnectivity();
	};

	// =================================================================================== //

	void computeGhostsConnectivity(){				// Computes ghosts nodes vector and connectivity of ghosts octants of local tree
		map<uint64_t, vector<uint32_t> > mapnodes;
		map<uint64_t, vector<uint32_t> >::iterator iter, iterend;
		uint32_t i, k, counter;
		uint64_t morton;
		uint32_t noctants = size_ghosts;
		u32vector2D octnodes;
		uint8_t j;

		octnodes.reserve(global3D.nnodes);
		if (ghostsnodes.size() == 0){
			ghostsconnectivity.resize(noctants);
			for (i = 0; i < noctants; i++){
				ghosts[i].getNodes(octnodes);
				for (j = 0; j < global3D.nnodes; j++){
					morton = mortonEncode_magicbits(octnodes[j][0], octnodes[j][1], octnodes[j][2]);
					if (mapnodes[morton].size()==0){
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
				ghostsnodes[counter] = nodecasting;
				ghostsnodes[counter].shrink_to_fit();
				for(vector<uint32_t>::iterator iter2 = iter->second.begin()+3; iter2 != iter->second.end(); iter2++){
					if (ghostsconnectivity[(*iter2)].size()==0){
						ghostsconnectivity[(*iter2)].reserve(8);
					}
					ghostsconnectivity[(*iter2)].push_back(counter);
				}
				mapnodes.erase(iter++);
				counter++;
			}
			ghostsnodes.shrink_to_fit();
			//Lento. Solo per risparmiare memoria
			for (int ii=0; ii<noctants; ii++){
				ghostsconnectivity[ii].shrink_to_fit();
			}
			ghostsconnectivity.shrink_to_fit();
		}
		iter = mapnodes.end();

	};

	// =================================================================================== //

	void clearGhostsConnectivity(){					// Clear ghosts nodes vector and connectivity of ghosts octants of local tree
		u32vector2D().swap(ghostsnodes);
		u32vector2D().swap(ghostsconnectivity);
	};

	// =================================================================================== //

	void updateGhostsConnectivity(){					// Update ghosts nodes vector and connectivity of ghosts octants of local tree
		clearGhostsConnectivity();
		computeGhostsConnectivity();
	};

	// =================================================================================== //

};//end Class_Local_Tree<3> specialization;


