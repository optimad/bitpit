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

#include "binaryUtils.hpp"

namespace bitpit {

namespace utils {

namespace binary {

/*!
    Write the given container to the specified stream in binary format.

    \param stream is the stream to write to
    \param container is the container to write
*/
template<>
void write(std::ostream &stream, const std::vector<bool> &container)
{
    std::copy(container.begin(), container.end(), std::ostream_iterator<bool>(stream));
}

/*!
    Write the given string to the specified stream in binary format.

    \param stream is the stream to write to
    \param string is the string to write
*/
void write(std::ostream &stream, const std::string &string)
{
    size_t size = string.length() + 1;
    write(stream, size);
    write(stream, string.c_str(), size);
}

/*!
    Read the given container to the specified stream in binary format.

    \param stream is the stream to write to
    \param container is the container to write
*/
template<>
void read(std::istream &stream, std::vector<bool> &container)
{
    std::copy(std::istream_iterator<bool>(stream), {}, std::back_inserter(container));
}

/*!
    Read the given string to the specified stream in binary format.

    \param stream is the stream to read from
    \param[out] string on output it will contain the read string
*/
void read(std::istream &stream, std::string &string)
{
    size_t size;
    read(stream, size);
    char * cstring = new char [size];
    read(stream, cstring, size);
    string = cstring;
}

}

}

}
