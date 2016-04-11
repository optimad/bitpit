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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include "binary_stream.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

namespace bitpit{
// ========================================================================== //
// IMPLEMENTATIONS OF METHODS FOR CLASS IBinaryStream                         //
// ========================================================================== //

/*!
 * @ingroup BinaryStream
 * @{
 */

/*!
 * @class   IBinaryStream
 * @brief   creates input binary stream
 */


// Constructor(s) =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        Default constructor. Initialize an empty object of class IBinaryStream

*/
IBinaryStream::IBinaryStream(
    void
) {
    current_pos = 0;
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #1. Initialize an empty object of class IBinaryStream
        with assigned size

        \param[in] size buffer size
*/
IBinaryStream::IBinaryStream(
    size_t                      size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #2. Initialize a object of class IBinaryStream
        pointing to a memory location with specified size

        \param[in] buf_ pointer to memory location
        \param[in] size buffer size

*/
IBinaryStream::IBinaryStream(
    const char                  *buf_,
    size_t                       size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    buffer.assign(buf_, buf_ + size);
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #3. Initialize a object of class IBinaryStream
        from a std::vector of char

        \param[in] vec input vector

*/
IBinaryStream::IBinaryStream(
    const std::vector<char>          &vec
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(vec.size());
    buffer.assign(vec.begin(), vec.end());
}

// Destructor(s) ============================================================ //
// default

// Assignament operator(s) ================================================== //
// disabled

// Public methods =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        Resize the buffer stream

        \param[in] size is the new size (in bytes) of the stream

*/
void IBinaryStream::resize(
    size_t                       size
) {
    buffer.resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Size the buffer stream

        \return The size of the buffer stream

*/
size_t IBinaryStream::size(
    void
) {
    return buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Open stream from memory

        \param[in] mem pointer to memory location
        \param[in] size size (in bytes) of memory location to be streamed

*/
void IBinaryStream::open(
    const char                  *mem,
    size_t                       size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    buffer.assign(mem, mem + size);
}

// -------------------------------------------------------------------------- //
/*!
        Close stream from memory

*/
void IBinaryStream::close(
    void
)
{
    buffer.clear();
}

// -------------------------------------------------------------------------- //
/*!
        Returns true if end of file condition is met.

        \result boolean flag (true) if end of file is reached, (false) otherwise

*/
bool IBinaryStream::eof(
    void
) const
{
    return current_pos >= buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Returns cursor position within the current buffer

        \result cursor position

*/
ifstream::pos_type IBinaryStream::tellg(
    void
) {
    return current_pos;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] pos cursor new position

        \result (true) if new position is valid, (false) otherwise

*/
bool IBinaryStream::seekg (
    size_t                       pos
) {
    if(pos<buffer.size())
        current_pos = pos;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] offset new offset position
        \param[in] way offset direction

        \result (true) if new position is valid, (false) otherwise

*/
bool IBinaryStream::seekg (
    std::streamoff               offset,
    std::ios_base::seekdir       way
) {
    if ( ( way == ios_base::beg ) && ( offset < (long) buffer.size() ) )
        current_pos = offset;
    else if ( ( way == ios_base::cur ) && ( current_pos + offset < buffer.size() ) )
        current_pos += offset;
    else if ( ( way == ios_base::end ) && ( (long) buffer.size() - offset >= 0 ) )
        current_pos = buffer.size() - offset;
    else
        return false;

    return true;
}

// Private method(s) ======================================================== //

// -------------------------------------------------------------------------- //
/*!
        Read data from memory location pointed by p and store into stream buffer

        \param[in] p pointer to memory location where data are stored
        \param[in] size size (in bytes) of data chunk to be read

*/
void IBinaryStream::read(
    char                        *p,
    size_t                       size
) {
    if ( eof() || (current_pos + size) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>( p ), &buffer[current_pos], size);
    current_pos += size;
}

// -------------------------------------------------------------------------- //
/*!
        Read data from std::vector<char> and store content into stream buffer

        \param[in] vec vector of char

*/
void IBinaryStream::read(
    std::vector<char>           &vec
) {
    if ( eof() || ( current_pos + vec.size() ) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>(&vec[0]), &buffer[current_pos], vec.size());
    current_pos += vec.size();
}

/*!
 * @}
 */

// ========================================================================== //
// IMPLEMENTATIONS OF METHODS FOR CLASS OBinaryStream                           //
// ========================================================================== //

/*!
 * @ingroup BinaryStream
 * @{
 */

/*!
 * @class   OBinaryStream
 * @brief   creates output binary stream
 */

// Constructor(s) =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        \ingroup BinaryStream
        Default constructor. Initialize an empty object
*/
OBinaryStream::OBinaryStream(
    void
) {
    current_pos = 0;
}

// -------------------------------------------------------------------------- //
/*!
        Default constructor. Initialize an empty object with buffer of specified
        size.

        \param[in] size buffer size
*/
OBinaryStream::OBinaryStream(
    size_t                       size
) {
    current_pos = 0;
    open(size);
}

// Destructor(s) ============================================================ //
// default

// Assignament operator(s) ================================================== //
// disabled

// Public method(s) ========================================================= //

// -------------------------------------------------------------------------- //
/*!
        Resize the buffer stream

        \param[in] size is the new size (in bytes) of the stream

*/
void OBinaryStream::resize(
    size_t                       size
) {
    buffer.resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Size the buffer stream

        \return The size of the buffer stream

*/
size_t OBinaryStream::size(
    void
) {
    return buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Open output stream

        \param[in] size stream size

*/
void OBinaryStream::open(
    size_t                       size
) {
    buffer.reserve(size);
    buffer.resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Close output stream

*/
void OBinaryStream::close(
    void
) {
    buffer.clear();
}

// -------------------------------------------------------------------------- //
/*!
        Returns true if end of file condition is met.

        \result boolean flag (true) if end of file is reached, (false) otherwise

*/
bool OBinaryStream::eof(
    void
) const
{
    return current_pos >= buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Returns cursor position within the current buffer

        \result cursor position

*/
ifstream::pos_type OBinaryStream::tellg(
    void
) {
    return current_pos;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] pos cursor new position

        \result (true) if new position is valid, (false) otherwise

*/
bool OBinaryStream::seekg (
    size_t                       pos
) {
    if(pos < buffer.size())
        current_pos = pos;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] offset new offset position
        \param[in] way offset direction

        \result (true) if new position is valid, (false) otherwise

*/
bool OBinaryStream::seekg (
    std::streamoff               offset,
    std::ios_base::seekdir       way
) {
    if ( ( way == ios_base::beg ) && ( offset < (long) buffer.size() ) )
        current_pos = offset;
    else if ( ( way == ios_base::cur ) && ( current_pos + offset < buffer.size() ) )
        current_pos += offset;
    else if ( ( way == ios_base::end ) && ( (long) buffer.size() - offset >= 0 ) )
        current_pos = buffer.size() - offset;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Write char array to internal buffer

        \param[in] p pointer to char array
        \param[in] size size of memory chunck to be writting to the internal
        buffer

*/
void OBinaryStream::write(
    const char                  *p,
    size_t                       size
) {
    if ( buffer.size() - current_pos < size ) {
        buffer.resize( size + current_pos );
    }
    for (size_t i = 0;  i < size; ++i) {
        buffer[current_pos] = p[i];
        ++current_pos;
    } // next i
}

// -------------------------------------------------------------------------- //
/*!
        Write vector of char to internal buffer

        \param[in] vec vector of char to be written in the internal buffer
        buffer

*/

void OBinaryStream::write(
    const vector<char>          &vec
) {
    if ( buffer.size() - current_pos < vec.size() ) {
        buffer.resize( vec.size() + current_pos );
    }
    for(size_t i = 0; i < vec.size(); ++i) {
        buffer[current_pos] = vec[i];
        ++current_pos;
    } //next i
}

}

// ========================================================================== //
// OPERATORS
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Explicit template specialization for stream operator for class IBinaryStream

        \param[in] istm input stream
        \param[in] val std::string to be streamed

*/
template<>
bitpit::IBinaryStream& operator>>(
        bitpit::IBinaryStream             &istm,
        std::string                  &val)
{
    int                 size = 0;

    istm.read(size);

    if(size<=0)         return istm;

    std::vector<char> vec((size_t)size);
    istm.read(vec);
    val.assign(&vec[0], (size_t)size);

    return istm;
}

// -------------------------------------------------------------------------- //
/*!
        Stream std::string to internal buffer.

        \param[in] ostm output stream
        \param[in] val string

*/
template<>
bitpit::OBinaryStream& operator << (
    bitpit::OBinaryStream                 & ostm,
    const std::string                & val
) {
    int size = val.size();

    ostm.write(size);

    if(val.size()<=0)
        return ostm;

    ostm.write(val.c_str(), val.size());

    return ostm;
}

// -------------------------------------------------------------------------- //
/*!
        Stream char array to internal buffer.

        \param[in] ostm output stream
        \param[in] val pointer to char array

*/
bitpit::OBinaryStream& operator<<(
    bitpit::OBinaryStream                 &ostm,
    const char                  *val
) {
    int size = strlen(val);

    ostm.write(size);

    if(size<=0)
        return ostm;

    ostm.write(val, size);

    return ostm;
}


/*!
 * @}
 */
