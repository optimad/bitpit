#ifndef CLASS_ARRAY_HPP_
#define CLASS_ARRAY_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include <cstdint>
#include "mpi.h"

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 *
 *	\date		03/mar/2014
 *  \authors	Marco Cisternino
 *  \authors	Edoardo Lombardi
 *  \version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *  \brief Customized array definition
 *
 *  Class_Array contains a pointer to an array of integer values. Implemented here for fast using in PABLO.
 *
 */
class Class_Array {

	template<int dim> friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
	uint32_t arraySize;		/**< Size of array */
	int* array;				/**< Pointer to array of integers */

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Array();
	Class_Array(uint32_t size, int value);
	Class_Array(const Class_Array& other);
	~Class_Array();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	Class_Array& operator=(const Class_Array& rhs);
};

#endif /* CLASS_ARRAY_HPP_ */
