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

#include "GenericIO.hpp"

namespace bitpit {

namespace genericIO {

/*!
 * Writes uint_8 data as formatted integer 
 * @param[in] str file stream to be copied; file needs already to be opened
 * @param[in] data data to be written
 */
template<>
void flushASCII(std::fstream &str, const uint8_t &data)
{
    std::ios::fmtflags streamFlags(str.flags());

    str << std::setprecision(8) << std::scientific;
    str << unsigned(data) << " ";

    str.flags(streamFlags);
}

/*!
 *
 * Copies entire file into a char array.
 *
 * NOTE: only valid characters will be copied.
 *
 * NOTE: On Windows OS, read-from-file \r\n characters take 2 position in
 * the stream, but are accounted as a unique character (endline) after C++
 * translation with read function.
 *
 * @param[in] str file stream to be copied; file needs already to be opened
 * @param[out] buffer char array containing entire array (only valid characters)
 * @param[out] length number of elements which have been copied
 */
void copyUntilEOFInString(std::fstream &str, char *&buffer, int &length)
{
    std::fstream::pos_type position_insert, position_eof;

    // Ask for actual position and end-of-file positions
    position_insert = str.tellg();
    str.seekg(0, std::ios::end);
    position_eof = str.tellg();

    // Evaluate the length on file. This is a prediction of the maximum
    // number of charactes readable in the current file slot.
    length = position_eof - position_insert;

    // Instantiate a char vector for reading purposes
    std::vector<char> trybuf(length);

    // Clear and get again the stream to the initial insert position.
    str.clear();
    str.seekg(position_insert);

    // Read data
    str.read(trybuf.data(), length);

    // Resize the trybuf with the effective valid character read.
    //
    // This number can be less than length, especially while reading OS-Windows
    // files.
    trybuf.resize(str.gcount());
    length = trybuf.size();

    // Allocate your real buffer and copy effective data in it.
    buffer = new char[trybuf.size()];
    std::copy(trybuf.begin(), trybuf.end(), buffer);

    // Clear the stream and get the read stream input to the initial position.
    str.clear();
    str.seekg(position_insert);
}

}

}
