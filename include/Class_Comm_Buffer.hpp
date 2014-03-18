/*
 * ClassCommBuffer.hpp
 *
 *  Created on: 25/feb/2014
 *      Author: Marco Cisternino
 */

#ifndef CLASSCOMMBUFFER_HPP_
#define CLASSCOMMBUFFER_HPP_

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

class Class_Comm_Buffer {

	friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
	uint32_t commBufferSize;
	char* commBuffer;

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Comm_Buffer();
	Class_Comm_Buffer(uint32_t size, char value);
	Class_Comm_Buffer(const Class_Comm_Buffer& other);
	~Class_Comm_Buffer();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	//TODO routines write and read to write and read POD types in buffer
	Class_Comm_Buffer& operator=(const Class_Comm_Buffer& rhs);
};

#endif /* CLASSCOMMBUFFER_HPP_ */
