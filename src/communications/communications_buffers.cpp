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

#include "communications_buffers.hpp"

namespace bitpit {

/*!
    \class SendBuffer
    \ingroup communications

    @brief Buffer to be used for send communications.
*/

/*!
    Creates a new send buffer
*/
SendBuffer::SendBuffer(size_t size, bool doubleBuffer)
    : CommunicationBuffer(size, doubleBuffer)
{
}

/*!
    Requests the buffer to reduce its size to fit the data currently
    contained in the stream.
*/
void SendBuffer::squeeze()
{
    for (RawSendBuffer & buffer : getBuffers()) {
        buffer.squeeze();
    }
}

/*!
* Write data into the buffer.
*
* \param[in] data is the memory that contain the data
* \param[in] size is the size (in bytes) of the data to be written into
* the stream
*/
void SendBuffer::write(const char *data, std::size_t size)
{
    getFront().write(data, size);
}

/*!
    \class RecvBuffer
    \ingroup communications

    @brief Buffer to be used for receive communications.
*/

/*!
    Creates a new receive buffer
*/
RecvBuffer::RecvBuffer(size_t size, bool doubleBuffer)
    : CommunicationBuffer(size, doubleBuffer)
{
}

/*!
* Read data from the buffer.
*
* \param[out] data is the memory location that will contain the data
* \param[in] size is the size (in bytes) of the data to be read from the
* stream
*/
void RecvBuffer::read(char *data, std::size_t size)
{
    getFront().read(data, size);
}

}
