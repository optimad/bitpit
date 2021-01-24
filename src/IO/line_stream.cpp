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

#include <cstring>

#include "line_stream.hpp"

namespace bitpit {

/*!
    \class LineBuffer
    \brief LineBuffer defines a buffer for reading lines.
*/

/*!
    Constructor.
*/
LineBuffer::LineBuffer()
    : m_buffer(CHUNK_SIZE)
{
}

/*!
    Read a line from the file.

    \param fileHandle is the file handle
*/
int LineBuffer::readLine(std::ifstream &fileHandle)
{
    int lineLength = 0;
    while (true) {
        if (m_buffer.size() - lineLength < CHUNK_SIZE) {
            m_buffer.resize(m_buffer.size() + CHUNK_SIZE);
        }

        fileHandle.getline(m_buffer.data() + lineLength, CHUNK_SIZE);
        lineLength = strlen(m_buffer.data());

        // If we haven't find the end of the line, we keep reading
        if (fileHandle.good() || fileHandle.eof() || fileHandle.bad()) {
            break;
        } else if (fileHandle.fail()) {
            fileHandle.clear();
        }
    }

    setg(m_buffer.data(), m_buffer.data(), m_buffer.data() + lineLength);

    if ((fileHandle.eof() || fileHandle.bad()) && lineLength == 0) {
        return -1;
    } else {
        return lineLength;
    }
}

/*!
    Get a copy of the current line.

    \param line on output will contain the current line
*/
void LineBuffer::copyLine(std::string *line) const
{
    line->assign(eback(), egptr() - eback());
}

/*!
    \class LineStream
    \brief LineStream defines a stream for reading lines.
*/

/*!
    Constructor.
*/
LineStream::LineStream()
    : LineBuffer(), std::istream(this)
{
}

/*!
    Read a line from the file.

    \param fileHandle is the file handle
*/
int LineStream::readLine(std::ifstream &fileHandle)
{
    std::istream::clear();

    return static_cast<LineBuffer *>(std::istream::rdbuf())->readLine(fileHandle);
}

}
