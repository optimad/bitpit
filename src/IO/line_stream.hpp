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

#ifndef __BITPIT_LINE_STREAM_HPP__
#define __BITPIT_LINE_STREAM_HPP__

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace bitpit {

class LineBuffer : public std::streambuf {

public:
    int readLine(std::ifstream &fileHandle);
    void copyLine(std::string *line) const;

protected:
    LineBuffer();

private:
    static const int CHUNK_SIZE = 1024;

    std::vector<char> m_buffer;

};

class LineStream : protected LineBuffer, public std::istream {

public:
    LineStream();

    int readLine(std::ifstream &fileHandle);

    using LineBuffer::copyLine;

};

}

#endif
