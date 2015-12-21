/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Intersection class definition - 3D specialization
 *
 *	The intersection is the face or portion of face shared by two octants. An intersection is defined
 *	by :
 *	- the owner octants, i.e. the octants sharing the intersection, identified by a couple (array[2]) of indices;
 *	- the index of the face, that contains the intersection, of the first owner;
 *	- an identifier of the octant with higher level of refinement (0/1) [if same level identifier =0];
 *	- a flag stating if an owner is ghost;
 *	- a flag to communicate if the intersection is new after a mesh refinement.
 *
 */

// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

template<>
class Class_Intersection<3> {
	// ------------------------------------------------------------------------------- //
	// FRIENDSHIPS ------------------------------------------------------------------- //
	template<int dim> friend class Class_Local_Tree;
	template<int dim> friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// TYPEDEFS ----------------------------------------------------------------------- //
public:
	typedef vector<uint32_t>			u32vector;
	typedef vector<vector<uint32_t>	>	u32vector2D;
	typedef vector<vector<uint64_t>	>	u64vector2D;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	uint32_t 	owners[2];			/**< Owner octants of the intersection (first is the internal octant) */
	uint8_t   	iface;				/**< Index of the face of the first owner */
	bool		finer;				/**< 0/1 finer octant (if same level =0) */
	//TODO proposta: eliminare isghost, contiene la stessa info in pbound
	bool		isghost;			/**< The intersection has a member ghost */
	bool		isnew;				/**< The intersection is new after a mesh adapting? */
	bool		bound;				/**< The intersection is a boundary intersection of the whole domain */
	bool		pbound;				/**< The intersection is a boundary intersection of a process domain */

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS AND OPERATORS----------------------------------------------------- //

public:
	Class_Intersection(){
		owners[0] = 0;
		owners[1] = 0;
		iface = 0;
		isnew = false;
		isghost = false;
		finer = 0;
		bound = pbound = false;
	};
	~Class_Intersection(){};
	Class_Intersection(const Class_Intersection<3> & intersection){
		owners[0] = intersection.owners[0];
		owners[1] = intersection.owners[1];
		iface = intersection.iface;
		isnew = intersection.isnew;
		isghost = intersection.isghost;
		finer = intersection.finer;
		bound = intersection.bound;
		pbound = intersection.pbound;
	};
	Class_Intersection<3>& operator =(const Class_Intersection<3> & intersection){
		owners[0] = intersection.owners[0];
		owners[1] = intersection.owners[1];
		iface = intersection.iface;
		isnew = intersection.isnew;
		isghost = intersection.isghost;
		finer = intersection.finer;
		bound = intersection.bound;
		pbound = intersection.pbound;
		return *this;
	};
	bool operator ==(const Class_Intersection<3> & intersection){
		bool check = true;
		check = check && (owners[0] == intersection.owners[0]);
		check = check && (owners[1] == intersection.owners[1]);
		check = check && (iface == intersection.iface);
		check = check && (isnew == intersection.isnew);
		check = check && (isghost == intersection.isghost);
		check = check && (finer == intersection.finer);
		check = check && (bound == intersection.bound);
		check = check && (pbound == intersection.pbound);
		return check;

	};

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

	uint32_t getOut(){						// Get the owner with exiting normal
		return owners[0];
	};
	uint32_t getIn(){						// Get the owner with entering normal
		return owners[1];
	};
	void getNormal(int8_t normal[3]){		// Get the normal of the intersection
		for (int i=0; i<3; i++){
			normal[i] = global3D.normals[iface][i];
		}
	};
	bool getBound(){
		return bound;
	}

	bool getIsGhost(){
		return isghost;
	}

	bool getPbound(){
		return pbound;
	}
}; // end of Class_Intersection_3D.tpp

/*  @}  */
