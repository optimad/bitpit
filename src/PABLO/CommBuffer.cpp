/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#if BITPIT_ENABLE_MPI==1

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "CommBuffer.hpp"
#include <cstdint>
#include <algorithm>

namespace bitpit {

    // =================================================================================== //
    // NAME SPACES                                                                         //
    // =================================================================================== //
    using namespace std;


    /*! Default constructor. It constructs an empty buffer with MPI_COMM_WORLD as MPI-communicator
     */
    CommBuffer::CommBuffer(){

        m_commBufferSize = 0;
        m_commBuffer = NULL;
        m_pos = 0;
        m_comm = MPI_COMM_WORLD;
    }

    /*! Constructor. It constructs an empty buffer with comm_ as MPI-communicator
     * \param[in] comm_ The user MPI-communicator
     */
    CommBuffer::CommBuffer(MPI_Comm comm_) : m_comm(comm_){

        m_commBufferSize = 0;
        m_commBuffer = NULL;
        m_pos = 0;
    }

    /*! Constructor. It constructs an value initilized buffer of size size with comm_ as MPI-communicator
     * \param[in] size The size of the buffer
     * \param[in] value The initilization value for the buffer
     * \param[in] comm_ The user MPI-communicator
     */
    CommBuffer::CommBuffer(uint32_t size, char value, MPI_Comm comm_) : m_comm(comm_){

        m_commBufferSize = size;
        m_commBuffer = new char [size];
        for(uint32_t i = 0; i < size; ++i)
            m_commBuffer[i] = value;
        m_pos = 0;
    }

    /*! Copy constructor. It constructs a buffer by copy
     * \param[in] other The buffer to be copied in
     */
    CommBuffer::CommBuffer(const CommBuffer& other) {

        m_commBufferSize = other.m_commBufferSize;
        //	if(commBuffer != NULL){
        //		delete [] commBuffer;
        //		commBuffer = NULL;
        //	}
        m_commBuffer = new char [m_commBufferSize];
        for(uint32_t i = 0; i < m_commBufferSize; ++i)
            m_commBuffer[i] = other.m_commBuffer[i];
        m_pos = other.m_pos;
        m_comm = other.m_comm;
    }

    /*! Destructor. It dinamically deletes the buffer
     * \param[in] comm_ The user MPI-communicator
     */
    CommBuffer::~CommBuffer() {
		delete [] m_commBuffer;
		m_commBuffer = NULL;
    }

    /*! Assignement operator. It assigns a buffer all the memebers of another buffer
     * \param[in] rhs The other buffer whose member have to be copied in
     */
    CommBuffer& CommBuffer::operator =(const CommBuffer& rhs) {
        if(this != &rhs)
            {
                char* new_array = new char[rhs.m_commBufferSize];
                copy(rhs.m_commBuffer,rhs.m_commBuffer+rhs.m_commBufferSize,new_array);

                delete [] m_commBuffer;

                m_commBuffer = new_array;
                m_commBufferSize = rhs.m_commBufferSize;
                m_pos = rhs.m_pos;
                m_comm = rhs.m_comm;
            }
        return *this;

    }

}

#endif
