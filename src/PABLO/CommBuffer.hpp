/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#if ENABLE_MPI==1
#ifndef COMMBUFFER_HPP_
#define COMMBUFFER_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //

#include "mpi.h"
#include "mpi_datatype_conversion.hpp"

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
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Bundle char container for communications
 *
 *	This calls is intended to provide the user with a basic container
 *	for data MPI communications.
 *
 *	The user should use this container implementing his communications
 *	interface specializations.
 *
 *	More precisely, he has to call read/write methods to read/write
 *	every MPI-compatible POD datum in the buffer.
 *	By this way, data communications are data independent.
 */
class CommBuffer {

	friend class ParaTree;

	// =================================================================================== //
	// MEMBERS																			   //
	// =================================================================================== //
	uint32_t m_commBufferSize;
	char* m_commBuffer;
	int m_pos;
	MPI_Comm m_comm;

	// =================================================================================== //
	// CONSTRUCTORS 																	   //
	// =================================================================================== //
public:
	CommBuffer();
	CommBuffer(MPI_Comm comm_);
	CommBuffer(uint32_t size, char value, MPI_Comm comm_);
	CommBuffer(const CommBuffer& other);
	~CommBuffer();

	// =================================================================================== //
	// METHODS                                                                		       //
	// =================================================================================== //
	//TODO routines write and read to write and read POD types in buffer
	CommBuffer& operator=(const CommBuffer& rhs);

	// =================================================================================== //
	// TEMPLATE METHODS                                                                    //
	// =================================================================================== //

	/*! This method writes a MPI-compatible POD datum of type T in commBuffer
	 * \param[in] val The values that has to be written in the buffer.
	*/
	template<class T>
	void write(T& val) {
		MPI_Datatype datatype = convert<T>();
		int error = MPI_Pack(&val,1,datatype,m_commBuffer,m_commBufferSize,&m_pos,m_comm);
	};

	/*! This method reads from commBuffer the user MPI-compatible POD datum of type T.
	 * \param[in] val The values that has to be read from the buffer.
	 */
	template<class T>
	void read(T& val) {
		MPI_Datatype datatype = convert<T>();
		int error = MPI_Unpack(m_commBuffer,m_commBufferSize,&m_pos,&val,1,datatype,m_comm);
	};

};

/* @} */

#endif /* COMMBUFFER_HPP_ */
#endif /* NOMPI */
