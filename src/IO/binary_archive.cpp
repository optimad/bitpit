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

#include "bitpit_common.hpp"

#include "binary_archive.hpp"
#include "FileHandler.hpp"

namespace bitpit {

/*!
    \class BinaryArchive
    \ingroup Binary

    \brief Base class for binary archives.

    The BinaryArchive is the blas class from which input and output binary
    archives are derived from.
*/

const std::string BinaryArchive::EXTENSION_DEFAULT = "dat";

/*!
    Default constructor
*/
BinaryArchive::BinaryArchive()
    : m_version(VERSION_UNDEFINED)
{
}

/*!
    Destructor
*/
BinaryArchive::~BinaryArchive()
{
    close();
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param mode specifies the opening mode
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void BinaryArchive::open(const std::string &name, const std::string &extension,
                         ios_base::openmode mode, int block)
{
    if (is_open()) {
        close();
    }

    FileHandler fileHandler("", name, extension);
    if (block >= 0) {
        fileHandler.setParallel(true);
        fileHandler.setBlock(block);
    }

    std::fstream::open(fileHandler.getPath().c_str(), std::ios::binary | mode);
    if (fail() && bad()) {
        return;
    }
}

/*!
    Gets the version associated to the archive.

    \result The version associated to the archive.
*/
int BinaryArchive::getVersion() const
{
    return m_version;
}

/*!
    Gets the header associated to the archive.

    \result The header associated to the archive.
*/
std::string BinaryArchive::getHeader() const
{
    return m_header;
}

/*!
    \ingroup binary
    \class IBinaryArchive
    \brief Input binary archive.

    The IBinaryArchive class can read a binary archive, the binary archive
    has a fixed ASCII header.
*/

/*!
    Creates a new input archive.

    \param name is the name of the file
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
IBinaryArchive::IBinaryArchive(const std::string &name, int block)
{
    open(name, EXTENSION_DEFAULT, block);
}

/*!
    Creates a new input archive.

    \param name is the name of the file
    \param extension is the extension of the file
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
IBinaryArchive::IBinaryArchive(const std::string &name, const std::string &extension, int block)
{
    open(name, extension, block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void IBinaryArchive::open(const std::string &name, int block)
{
    open(name, EXTENSION_DEFAULT, block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void IBinaryArchive::open(const std::string &name, const std::string &extension, int block)
{
    // Reset the header
    m_header.clear();

    // Open the stream
    BinaryArchive::open(name, extension, std::ios::in | std::ios_base::binary, block);

    // Read the header
    utils::binary::read(*this, m_header);

    // Read the version
    utils::binary::read(*this, m_version);
}

/*!
    Get a reference to the input stream associated to the archive.

    \result A reference to the input stream associated to the archive.
*/
std::istream & IBinaryArchive::getStream()
{
    return *this;
}

/*!
    Checks if the version of the archive matches the specified version.

    \param version is the requested version
    \result Retuns true if the version of the archive matches the specified
    version, otherwise it returns false.
*/
bool IBinaryArchive::checkVersion(int version)
{
    return (version == m_version);
}

/*!
    \ingroup binary
    \class OBinaryArchive
    \brief Output binary archive.

    The IBinaryArchive class can write a binary archive, the binary archive
    has a fixed ASCII header.
*/

/*!
    Creates a new output archive.

    \param name is the name of the file
    \param version is the version of the archive
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
OBinaryArchive::OBinaryArchive(const std::string &name, int version, int block)
{
    open(name, EXTENSION_DEFAULT, version, "", block);
}

/*!
    Creates a new output archive.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param header is the header of the archive, the length of the header is
    limited to HEADER_SIZE characters.
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
OBinaryArchive::OBinaryArchive(const std::string &name, int version,
                               const std::string &header, int block)
{
    open(name, EXTENSION_DEFAULT, version, header, block);
}

/*!
    Creates a new output archive.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
OBinaryArchive::OBinaryArchive(const std::string &name, const std::string &extension,
                               int version, int block)
{
    open(name, extension, version, "", block);
}


/*!
    Creates a new output archive.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param header is the header of the archive, the length of the header is
    limited to HEADER_SIZE characters.
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
OBinaryArchive::OBinaryArchive(const std::string &name, const std::string &extension,
                               int version, const std::string &header, int block)
{
    open(name, extension, version, header, block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void OBinaryArchive::open(const std::string &name, int version, int block)
{
    open(name, EXTENSION_DEFAULT, version, "", block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param header is the header of the archive, the length of the header is
    limited to HEADER_SIZE characters.
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void OBinaryArchive::open(const std::string &name, int version,
                          const std::string &header, int block)
{
    open(name, EXTENSION_DEFAULT, version, header, block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void OBinaryArchive::open(const std::string &name, const std::string &extension,
                          int version, int block)
{
    open(name, extension, version, "", block);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param name is the name of the file
    \param extension is the extension of the file
    \param version is the version of the archive
    \param header is the header of the archive, the length of the header is
    limited to HEADER_SIZE characters.
    \param block is the parallel block thefile belongs to, a negative value
    mean that the file is serial
*/
void OBinaryArchive::open(const std::string &name, const std::string &extension,
                          int version, const std::string &header, int block)
{
    // Open the stream
    BinaryArchive::open(name, extension, std::ios::out | std::ios_base::binary, block);

    // Write the header
    std::string archiveHeader(header);
    archiveHeader.resize(HEADER_SIZE, ' ');
    utils::binary::write(*this, archiveHeader);

    // Write the version
    utils::binary::write(*this, version);
}

/*!
    Get a reference to the input stream associated to the archive.

    \result A reference to the input stream associated to the archive.
*/
std::ostream & OBinaryArchive::getStream()
{
    return *this;
}

}
