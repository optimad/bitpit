/*
 * Class_Intersection_3D.tpp
 *
 *  Created on: 22/apr/2014
 *      Author: Marco Cisternino
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
