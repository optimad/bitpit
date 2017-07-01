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

#include "binary_stream.hpp"

/*!
* Stream a string from the binary stream.
*
* \param[in] stream input stream
* \param[in] value is the string to be streamed
*/
template<>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &stream, std::string &value)
{
    int size = 0;
    stream.read(size);

    if (size > 0) {
        std::vector<char> buffer(size);
        stream.read(buffer);
        value.assign(buffer.data(), size);
    }

    return stream;
}

/*!
* Stream a string to the binary stream.
*
* \param[in] stream is the output stream
* \param[in] value is the string to be streamed
*/
template<>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &stream, const std::string &value)
{
    int size = value.size();
    stream.write(size);

    if (size > 0) {
        stream.write(value.c_str(), size);
    }

    return stream;
}


namespace bitpit {

/*!
* \class IBinaryStream
* \ingroup BinaryStream
*
* \brief Output binary stream.
*/

/*!
* Constructor.
*
* Initialize an empty binary stream.
*/
IBinaryStream::IBinaryStream()
    : m_pos(0)
{
}

/*!
* Constructor.
*
* Initialize an empty binary stream with the assigned capacity.
*
* \param[in] capacity is the stream capacity
*/
IBinaryStream::IBinaryStream(std::size_t capacity)
    : m_pos(0)
{
    setCapacity(capacity);
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] capacity is the buffer capacity
*/
IBinaryStream::IBinaryStream(const char *buffer, std::size_t capacity)
    : m_pos(0)
{
    setCapacity(capacity);
    m_buffer.assign(buffer, buffer + capacity);
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
*/
IBinaryStream::IBinaryStream(const std::vector<char> &buffer)
    : m_pos(0)
{
    setCapacity(buffer.size());
    m_buffer.assign(buffer.begin(), buffer.end());
}

/*!
* Open a binary stream initializing it with the data contained in the
* specified buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] capacity is the buffer capacity
*/
void IBinaryStream::open(const char *buffer, std::size_t capacity)
{
    m_pos = 0;

    setCapacity(capacity);
    m_buffer.assign(buffer, buffer + capacity);
}

/*!
* Returns true if end of file condition is met.
*
* \result Returns true if end of file condition is met, false otherwise.
*/
bool IBinaryStream::eof() const
{
    return (m_pos >= m_buffer.size());
}

/*!
* Returns the cursor position within the stream.
*
* \result Returns the cursor position within the stream.
*/
std::ifstream::pos_type IBinaryStream::tellg() const
{
    return m_pos;
}

/*!
* Set the cursor position within the stream.
*
* \param[in] pos is the cursor new position
* \result Returns true if new position is valid, false otherwise.
*/
bool IBinaryStream::seekg(std::size_t pos)
{
    if (pos >= m_buffer.size()) {
        return false;
    }

    m_pos = pos;

    return true;
}

/*!
* Get the internal data.
*
* \return The internal data.
*/
const std::vector<char> & IBinaryStream::data()
{
    return m_buffer;
}

/*!
* Get the internal data.
*
* \return The internal data.
*/
char * IBinaryStream::rawData()
{
    return m_buffer.data();
}

/*!
* Set the cursor position within the stream.
*
* \param[in] offset is the offset value, relative to the way parameter
* \param[in] way is the offset direction
* \result Returns rue if new position is valid, false otherwise.
*/
bool IBinaryStream::seekg(std::streamoff offset, std::ios_base::seekdir way)
{
    if ((way == std::ios_base::beg) && (offset < (long) m_buffer.size())) {
        m_pos = offset;
    } else if ((way == std::ios_base::cur) && (m_pos + offset < m_buffer.size())) {
        m_pos += offset;
    } else if ((way == std::ios_base::end) && ((long) m_buffer.size() - offset >= 0)) {
        m_pos = m_buffer.size() - offset;
    } else {
        return false;
    }

    return true;
}

/*!
* Set the capacity of the stream.
*
* \param[in] capacity is the new capacity (in bytes) of the stream
*/
void IBinaryStream::setCapacity(std::size_t capacity)
{
    m_buffer.resize(capacity);
}

/*!
* Get the capacity of the stream.
*
* \return The capacity of the stream
*/
std::size_t IBinaryStream::capacity() const
{
    return m_buffer.size();
}

/*!
* Read data from the stream.
*
* \param[out] data is the memory location that will contain the data
* \param[in] size is the size (in bytes) of the data to be read from the
* stream
*/
void IBinaryStream::read(char *data, std::size_t size)
{
    if (eof() || (m_pos + size) > m_buffer.size()) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>(data), m_buffer.data() + m_pos, size);
    m_pos += size;
}

/*!
* \class OBinaryStream
* \ingroup BinaryStream
*
* \brief Output binary stream.
*/

/*!
* Constructor.
*
* Initialize an empty binary stream.
*/
OBinaryStream::OBinaryStream()
    : m_pos(0)
{
}

/*!
* Constructor.
*
* Initialize a binary stream with the specified capacity.
*
* \param[in] capacity is the capacity of the stream
*/
OBinaryStream::OBinaryStream(std::size_t capacity)
    : m_pos(0)
{
    open(capacity);
}

/*!
* Open a stream with the specified capacity.
*
* \param[in] capacity is the capacity (in bytes) of the stream
*/
void OBinaryStream::open(std::size_t capacity)
{
    setCapacity(capacity);
}

/*!
* Returns true if end of file condition is met.
*
* \result Returns true if end of file condition is met, false otherwise.
*/
bool OBinaryStream::eof() const
{
    return (m_pos >= m_buffer.size());
}

/*!
* Returns the cursor position within the stream.
*
* \result Returns the cursor position within the stream.
*/
std::ofstream::pos_type OBinaryStream::tellg() const
{
    return m_pos;
}

/*!
* Set the cursor position within the stream.
*
* \param[in] pos is the cursor new position
* \result Returns true if new position is valid, false otherwise.
*/
bool OBinaryStream::seekg(std::size_t pos)
{
    if (pos >= m_buffer.size()) {
        return false;
    }

    m_pos = pos;

    return true;
}

/*!
* Set the cursor position within the stream.
*
* \param[in] offset is the offset value, relative to the way parameter
* \param[in] way is the offset direction
* \result Returns rue if new position is valid, false otherwise.
*/
bool OBinaryStream::seekg(std::streamoff offset, std::ios_base::seekdir way)
{
    if ((way == std::ios_base::beg) && (offset < (long) m_buffer.size())) {
        m_pos = offset;
    } else if ((way == std::ios_base::cur) && (m_pos + offset < m_buffer.size())) {
        m_pos += offset;
    } else if ((way == std::ios_base::end) && ((long) m_buffer.size() - offset >= 0)) {
        m_pos = m_buffer.size() - offset;
    } else {
        return false;
    }

    return true;
}

/*!
* Requests the stream to reduce its capacity to fit the data currently
* contained in the stream.
*/
void OBinaryStream::squeeze()
{
    setCapacity(m_pos);
}

/*!
* Get the internal data.
*
* \return The internal data.
*/
const std::vector<char> & OBinaryStream::data()
{
    return m_buffer;
}

/*!
* Get the internal data.
*
* \return The internal data.
*/
char * OBinaryStream::rawData()
{
    return m_buffer.data();
}

/*!
* Set the capacity of the stream.
*
* \param[in] capacity is the new capacity (in bytes) of the stream
*/
void OBinaryStream::setCapacity(std::size_t capacity)
{
    m_buffer.resize(capacity);
}

/*!
* Capacity of the stream.
*
* \return The capacity (in bytes) of the stream.
*/
std::size_t OBinaryStream::capacity() const
{
    return m_buffer.size();
}

/*!
* Write data into the stream.
*
* \param[in] data is the memory that contain the data
* \param[in] size is the size (in bytes) of the data to be written into
* the stream
*/
void OBinaryStream::write(const char *data, std::size_t size)
{
    if (m_buffer.size() - m_pos < size) {
        // We should take into account the current position in the buffer only
        // when the buffer is not empty. That's because when the buffer is
        // empty, the current position is set to 0 but there are no elements
        // in the buffer.
        std::size_t bufferSize = size;
        if (!m_buffer.empty()) {
            bufferSize += (m_pos + 1);
        }

        setCapacity(bufferSize);
    }

    for (std::size_t i = 0;  i < size; ++i) {
        m_buffer[m_pos] = data[i];
        ++m_pos;
    }
}

}
