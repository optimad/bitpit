/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

/*!
* Write the specified value into the stream.
*
* \param[in] stream is the input stream
* \param[in] value is the value to be streamed
* \result Returns the updated input stream.
*/
template<typename T>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &stream, T &value)
{
    stream.read(value);

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
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &stream, const T &value)
{
    stream.write(value);

    return stream;
}

namespace bitpit {

/*!
* Read the specified value from the stream.
*
* \param[in] value is the value that will be read
*/
template<typename T>
void IBinaryStream::read(T &value)
{
    read(reinterpret_cast<char *>(&value), sizeof(T));
}

/*!
* Write the specified value into the stream.
*
* \param[in] value is the value that will be written
*/
template<typename T>
void OBinaryStream::write(const T &value)
{
    write(reinterpret_cast<const char *>(&value), sizeof(T));
}

}
