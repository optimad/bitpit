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

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS IBinaryStream                             //
// ========================================================================== //

namespace bitpit{

// -------------------------------------------------------------------------- //
/*!
        Read data from memory location and store in the stream buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void IBinaryStream::read(
    T                           &t
) {
    if ( eof() || ( current_pos + sizeof(T) ) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>( &t ), &buffer[current_pos], sizeof(T));
    current_pos += sizeof(T);
}

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS OBinaryStream                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Write data from memory to internal buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void OBinaryStream::write(
    const T                     &t
) {
    std::vector<char> vec(sizeof(T));
    memcpy( reinterpret_cast<void*>( &vec[0] ), reinterpret_cast<const void*>( &t ), sizeof(T) );
    write(vec);
}

}

// ========================================================================== //
// OPERATORS                                                                  //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Input stream operator for class IBinaryStream.

        \param[in] istm input stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
bitpit::IBinaryStream& operator>> (
    bitpit::IBinaryStream       &istm,
    T                           &val
) {
    istm.read(val);
    return istm;
}

// -------------------------------------------------------------------------- //
/*!
        Output stream operator for class IBinaryStream.

        \param[in] ostm output stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
bitpit::OBinaryStream& operator<< (
    bitpit::OBinaryStream       &ostm,
    const T                     &val
) {
    ostm.write(val);

    return ostm;
}
