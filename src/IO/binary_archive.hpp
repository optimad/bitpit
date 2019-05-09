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
#ifndef __BITPIT_BINARY_ARCHIVE_HPP__
#define __BITPIT_BINARY_ARCHIVE_HPP__

#include <array>
#include <iostream>
#include <limits>
#include <fstream>

namespace bitpit {

class BinaryArchive : protected std::fstream
{

public:
    static const int HEADER_SIZE       = 1024;
    static const int VERSION_UNDEFINED = - std::numeric_limits<int>::max();

    static const std::string EXTENSION_DEFAULT;

    using std::fstream::close;

    BinaryArchive();
    ~BinaryArchive();

    int getVersion() const;
    std::string getHeader() const;

protected:
    int m_version;
    std::string m_header;

    void open(const std::string &name, const std::string &extension,
              ios_base::openmode mode, int block = -1);

};

class IBinaryArchive : public BinaryArchive
{

public:
    IBinaryArchive(const std::string &name, int block = -1);
    IBinaryArchive(const std::string &name, const std::string &extension, int block = -1);

    void open(const std::string &name, int block = -1);
    void open(const std::string &name, const std::string &extension, int block = -1);

    std::istream & getStream();

    bool checkVersion(int version);

    using BinaryArchive::operator>>;
    using BinaryArchive::read;

};

class OBinaryArchive : public BinaryArchive
{

public:
    OBinaryArchive(const std::string &name, int version,
                   int block = -1);
    OBinaryArchive(const std::string &name, int version, const std::string &header,
                   int block = -1);
    OBinaryArchive(const std::string &name, const std::string &extension, int version,
                   int block = -1);
    OBinaryArchive(const std::string &name, const std::string &extension, int version,
                   const std::string &header, int block = -1);

    void open(const std::string &name, const int version,
              int block = -1);
    void open(const std::string &name, int version, const std::string &header,
              int block = -1);
    void open(const std::string &name, const std::string &extension, int version,
              int block = -1);
    void open(const std::string &name, const std::string &extension, int version,
              const std::string &header, int block = -1);

    std::ostream & getStream();

    using BinaryArchive::operator<<;
    using BinaryArchive::write;

};

}

#endif
