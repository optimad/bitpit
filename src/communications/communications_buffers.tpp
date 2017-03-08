/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#ifndef __BITPIT_COMMUNICATIONS_BUFFERS_TPP__
#define __BITPIT_COMMUNICATIONS_BUFFERS_TPP__

namespace bitpit {

/*!
    Create a new communication buffer
 */
template<typename RawBufferType>
CommunicationBuffer<RawBufferType>::CommunicationBuffer(size_t capacity, bool doubleBuffer)
{
    // Create the buffers
    m_buffers.emplace_back(capacity);
    if (doubleBuffer) {
        m_buffers.emplace_back(capacity);
    }

    // Assign the front and back buffers
    m_front = &m_buffers[doubleBuffer ? 1 : 0];
    m_back  = &m_buffers[0];
}

/*!
    Set the capacity of the buffer

    \param capacity is the capacity of the buffer
 */
template<typename RawBufferType>
void CommunicationBuffer<RawBufferType>::setCapacity(size_t capacity)
{
    for (RawBufferType & buffer : m_buffers) {
        buffer.setCapacity(capacity);
    }
}

/*!
    Get the capacity of the buffer

    \result The capacity of the buffer.
 */
template<typename RawBufferType>
size_t CommunicationBuffer<RawBufferType>::capacity() const
{
    return m_front->capacity();
}

/*!
    Set cursor position within the front buffer

    \param pos is the new cursor position
    \result Return true if new position is valid, false otherwise
*/
template<typename RawBufferType>
bool CommunicationBuffer<RawBufferType>::seekg(size_t pos)
{
    return m_front->seekg(pos);
}

/*!
    Get cursor position within the front buffer

    \result Returns the cursor position.
*/
template<typename RawBufferType>
std::ifstream::pos_type CommunicationBuffer<RawBufferType>::tellg() const
{
    return m_front->tellg();
}

/*!
    Checks if the buffer is a double buffer.

    \result Returns true if the buffer is a double buffer, false otherwise.
 */
template<typename RawBufferType>
bool CommunicationBuffer<RawBufferType>::isDouble() const
{
    return (m_buffers.size() == 2);
}

/*!
    Get the front buffer

    \result The front buffer.
 */
template<typename RawBufferType>
RawBufferType & CommunicationBuffer<RawBufferType>::getFront()
{
    return *m_front;
}

/*!
    Get the back buffer

    \result The back buffer.
 */
template<typename RawBufferType>
RawBufferType & CommunicationBuffer<RawBufferType>::getBack()
{
    return *m_back;
}

/*!
    Swap front and back buffers.
 */
template<typename RawBufferType>
void CommunicationBuffer<RawBufferType>::swap()
{
    std::swap(m_front, m_back);
}

/*!
    Get a reference to the raw buffers
 */
template<typename RawBufferType>
std::vector<RawBufferType> & CommunicationBuffer<RawBufferType>::getBuffers()
{
    return m_buffers;
}

/*!
    Get a constant reference to the raw buffers
 */
template<typename RawBufferType>
const std::vector<RawBufferType> & CommunicationBuffer<RawBufferType>::getBuffers() const
{
    return m_buffers;
}

}

/*!
    Output stream operator for the send buffer.

    \param[in] buffer is the send buffer
    \param[in] value is the value to be streamed
    \result A reference of the buffer.
*/
template<typename T>
bitpit::SendBuffer & operator<<(bitpit::SendBuffer &buffer, const T &value)
{
    buffer.getFront() << value;

    return buffer;
}

/*!
    Input stream operator for the receive buffer.

    \param[in] buffer is the send buffer
    \param[in] value is the value to be streamed
    \result A reference of the buffer.
*/
template<typename T>
bitpit::RecvBuffer & operator>>(bitpit::RecvBuffer &buffer, T &value)
{
    buffer.getFront() >> value;

    return buffer;
}

# endif
