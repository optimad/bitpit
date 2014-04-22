/*
 * Class_Intersection.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: edoardo
 */

#ifndef CLASS_INTERSECTION_HPP_
#define CLASS_INTERSECTION_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Class_Local_Tree.hpp"
#include "Class_Octant.hpp"
#include "preprocessor_defines.dat"
#include "global.hpp"
#include "logFunct.hpp"
#include <vector>
#include <string.h>
#include "inlinedFunct.hpp"

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
template<int dim>
class Class_Intersection{};

#include "Class_Intersection_3D.tpp"
//#include "Class_Intersection_2D.tpp"


//
//class Class_Local_Tree;
//
//class Class_Intersection {
//	// ------------------------------------------------------------------------------- //
//	// FRIENDSHIPS ------------------------------------------------------------------- //
//	friend class Class_Local_Tree;
//	friend class Class_Para_Tree;
//
//	// ------------------------------------------------------------------------------- //
//	// TYPEDEFS ----------------------------------------------------------------------- //
//public:
//	typedef vector<uint32_t>			u32vector;
//	typedef vector<vector<uint32_t>	>	u32vector2D;
//	typedef vector<vector<uint64_t>	>	u64vector2D;
//
//	// ------------------------------------------------------------------------------- //
//	// MEMBERS ----------------------------------------------------------------------- //
//
//	uint32_t 	owners[2];			// Owner octants of the intersection (first is the internal octant)
//	uint8_t   	iface;				// Index of the face of the first owner
//	bool		finer;				// 0/1 finer octant (if same level =0)
//	bool		isghost;			// The intersection has a member ghost
//	bool		isnew;				// The intersection is new after a mesh adapting?
//
//	Class_Local_Tree *octree;
//
//	// ------------------------------------------------------------------------------- //
//	// CONSTRUCTORS AND OPERATORS----------------------------------------------------- //
//
//public:
//	Class_Intersection(Class_Local_Tree & tree);
//	Class_Intersection();
//	~Class_Intersection();
//	Class_Intersection(const Class_Intersection & intersection);
//	Class_Intersection& operator =(const Class_Intersection & intersection);
//	bool operator ==(const Class_Intersection & intersection);
//
//	// ------------------------------------------------------------------------------- //
//	// METHODS ----------------------------------------------------------------------- //
//
//	// Basic Get/Set methods --------------------------------------------------------- //
//
//	uint32_t getOut();						// Get the owner with exiting normal
//	uint32_t getIn();						// Get the owner with entering normal
//	void getNormal(int8_t normal[DIM]);		// Get the normal of the intersection
//	uint32_t getSize();
//	uint64_t getArea();
//	void getNodes(u32vector2D & nodes);
//	uint32_t (*getNodes())[DIM];
//	void getCenter(vector<double> & center);
//	double* getCenter();
//
//};


#endif /* CLASS_INTERSECTION_HPP_ */
