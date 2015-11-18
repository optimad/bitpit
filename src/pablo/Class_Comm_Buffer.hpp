#if NOMPI==0
#ifndef CLASSCOMMBUFFER_HPP_
#define CLASSCOMMBUFFER_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "mpi.h"
#include "preprocessor_defines.dat"
#include <cstdint>
#include <typeinfo>
#include <algorithm>
#include "mpi_datatype_conversion.hpp"

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //
/*!
 * \ingroup PABLO
 * @{
 *
 *	\date			09/sep/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Bundle char container for communications
 *
 *	This calls is intended to provide the user with a basic container for data MPI communications.
 *
 *	The user should use this container implementing his communications interface specializations.
 *
 *	More precisely, he has to call read/write methods to read/write every MPI-compatible POD datum in the buffer. By this way, data communications are data independent.
 */
class Class_Comm_Buffer {

	template<int dim> friend class Class_Para_Tree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
	uint32_t commBufferSize;
	char* commBuffer;
	int pos;
	MPI_Comm comm;

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Comm_Buffer();
	Class_Comm_Buffer(MPI_Comm comm_);
	Class_Comm_Buffer(uint32_t size, char value, MPI_Comm comm_);
	Class_Comm_Buffer(const Class_Comm_Buffer& other);
	~Class_Comm_Buffer();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	//TODO routines write and read to write and read POD types in buffer
	Class_Comm_Buffer& operator=(const Class_Comm_Buffer& rhs);

	/*! This method writes a MPI-compatible POD datum of type T in commBuffer
	 * \param[in] val The values that has to be written in the buffer.
	 */
	template<class T>
	void write(T& val);

	/*! This method reads from commBuffer the user MPI-compatible POD datum of type T.
	 * \param[in] val The values that has to be read from the buffer.
	 */
	template<class T>
	void read(T& val);
};

/* @} */

#include "Class_Comm_Buffer.tpp"

#endif /* CLASSCOMMBUFFER_HPP_ */
#endif /* NOMPI */
