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

namespace bitpit {

/*!
* Read the specified value from the stream.
*
* \param[in] stream is the input stream
* \param[in] value is the value to be streamed
* \result Returns the updated input stream.
*/
template<typename T>
IBinaryStream& operator>>(IBinaryStream &stream, T &value)
{
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));

    return stream;
}

/*!
* Read the specified vector from the stream.
*
* Input vector will be resized to match the size of the vector stored in the
* stream.
*
* \param[in] stream is the input stream
* \param[in] vector is the vector to be streamed
* \result Returns the updated input stream.
*/
template<typename T>
IBinaryStream& operator>>(IBinaryStream &stream, std::vector<T> &vector)
{
    std::size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));
    vector.resize(size);

    stream.read(reinterpret_cast<char *>(vector.data()), size * sizeof(T));

    return stream;
}

/*!
* Write the specified value into the stream.
*
* \param[in] stream is the output stream
* \param[in] value is the value to be streamed
* \result Returns the updated output stream.
*/
template<typename T>
OBinaryStream& operator<<(OBinaryStream &stream, const T &value)
{
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));

    return stream;
}

/*!
* Write the specified vector into the stream.
*
* Along with vector data, also the size of the vector is stored into the
* stream.
*
* \param[in] vector is the vector to be streamed
* \result Returns the updated output stream.
*/
template<typename T>
OBinaryStream& operator<<(OBinaryStream &stream, const std::vector<T> &vector)
{
    std::size_t size = vector.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

    stream.write(reinterpret_cast<const char *>(vector.data()), size * sizeof(T));

    return stream;
}

}
