/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#ifndef __BITPIT_COMMUNICATIONS_BUFFERS_HPP__
#define __BITPIT_COMMUNICATIONS_BUFFERS_HPP__

#include <vector>

#include "bitpit_containers.hpp"

// Declarations of the stream operators
namespace bitpit{
    class SendBuffer;
    class RecvBuffer;
};

// Declaration of the communications buffers
namespace bitpit {

class DataCommunicator;

typedef OBinaryStream RawSendBuffer;
typedef IBinaryStream RawRecvBuffer;

template<typename T>
SendBuffer & operator<<(SendBuffer &buffer, const T &value);

template<typename T>
RecvBuffer & operator>>(RecvBuffer &buffer, T &value);

/*!
    \ingroup communications

    @brief Buffer to be used for data communications.

    @tparam RawBufferType is the type of raw buffer that will be used
*/
template<typename RawBufferType>
class CommunicationBuffer
{

public:
    CommunicationBuffer(size_t size = 0, bool doubleBuffer = false);

    size_t getSize() const;
    void setSize(size_t size);

    bool seekg (size_t pos);
    std::ifstream::pos_type tellg(void) const;

    bool isDouble() const;

protected:
    RawBufferType & getFront();
    RawBufferType & getBack();

    void swap();

    std::vector<RawBufferType> & getBuffers();
    const std::vector<RawBufferType> & getBuffers() const;

private:
    std::vector<RawBufferType> m_buffers;
    RawBufferType *m_front;
    RawBufferType *m_back;

};

class SendBuffer : public CommunicationBuffer<RawSendBuffer>
{
    friend DataCommunicator;

    template<typename T>
    friend SendBuffer & (operator<<) (SendBuffer &buffer, const T &value);

public:
    SendBuffer(size_t size = 0, bool doubleBuffer = false);

    void squeeze();

    void write(const char *data, std::size_t size);

};

class RecvBuffer : public CommunicationBuffer<RawRecvBuffer>
{
    friend DataCommunicator;

    template<typename T>
    friend RecvBuffer & (operator>>) (RecvBuffer &buffer, T &value);

public:
    RecvBuffer(size_t size = 0, bool doubleBuffer = false);

    void read(char *data, std::size_t size);

};

}

#include "communications_buffers.tpp"

# endif
