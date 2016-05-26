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

#ifndef __BITPIT_COMMUNICATIONS_BUFFERS_HPP__
#define __BITPIT_COMMUNICATIONS_BUFFERS_HPP__

#include <vector>

#include "bitpit_containers.hpp"

// Declarations of the stream operators
namespace bitpit{
    class SendBuffer;
    class RecvBuffer;
};

template<typename T>
bitpit::SendBuffer & operator<<(bitpit::SendBuffer &buffer, const T &value);

template<typename T>
bitpit::RecvBuffer & operator>>(bitpit::RecvBuffer &buffer, T &value);

// Declaration of the communications buffers
namespace bitpit {

class DataCommunicator;

typedef OBinaryStream RawSendBuffer;
typedef IBinaryStream RawRecvBuffer;

/*!
    \ingroup communications

    @brief Buffer to be used for data communications.

    @tparam RawBufferType is the type of raw buffer that will be used
*/
template<typename RawBufferType>
class CommunicationBuffer
{

public:
    CommunicationBuffer(size_t capacity = 0, bool doubleBuffer = false);

    void setCapacity(size_t capacity);
    size_t capacity() const;

    bool seekg (size_t pos);
    std::ifstream::pos_type tellg(void) const;

    bool isDouble() const;

protected:
    RawBufferType & getFront();
    RawBufferType & getBack();

    void swap();

private:
    std::vector<RawBufferType> m_buffers;
    RawBufferType *m_front;
    RawBufferType *m_back;

};

class SendBuffer : public CommunicationBuffer<RawSendBuffer>
{
    friend DataCommunicator;

    template<typename T>
    friend SendBuffer & (::operator<<) (SendBuffer &buffer, const T &value);

public:
    SendBuffer(size_t capacity = 0, bool doubleBuffer = false);

};

class RecvBuffer : public CommunicationBuffer<RawRecvBuffer>
{
    friend DataCommunicator;

    template<typename T>
    friend RecvBuffer & (::operator>>) (RecvBuffer &buffer, T &value);

public:
    RecvBuffer(size_t capacity = 0, bool doubleBuffer = false);

};

}

#include "communications_buffers.tpp"

# endif
