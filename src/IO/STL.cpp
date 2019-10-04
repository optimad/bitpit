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

#include <cassert>

#include "bitpit_common.hpp"
#include "bitpit_operators.hpp"

#include "STL.hpp"
#include "logger.hpp"

namespace bitpit {

/*!
    \class STLBase
    \brief Base class for the STL writer and the STL reader.
*/

const std::size_t STLBase::BINARY_HEADER_SIZE  = 80 * sizeof(STLBase::BINARY_UINT8);
const std::size_t STLBase::BINARY_MINIMUM_SIZE = STLBase::BINARY_HEADER_SIZE + sizeof(STLBase::BINARY_UINT32);

const std::string STLBase::ASCII_SOLID_BEGIN  = "solid";
const std::string STLBase::ASCII_SOLID_END    = "endsolid";
const std::string STLBase::ASCII_FACET_BEGIN  = "facet";
const std::string STLBase::ASCII_FACET_END    = "endfacet";
const std::string STLBase::ASCII_FILE_BEGIN   = STLBase::ASCII_SOLID_BEGIN + " ";
const std::string STLBase::ASCII_FILE_END     = STLBase::ASCII_SOLID_END;
const std::size_t STLBase::ASCII_MINIMUM_SIZE = STLBase::ASCII_FILE_BEGIN.length() + STLBase::ASCII_FILE_END.length();

/*!
    Constructor.

    \param filename is the name of the STL file
*/
STLBase::STLBase(const std::string &filename)
{
    setFilename(filename);
    setFormat(FormatUnknown);
}

/*!
    Constructor.

    \param filename is the name of the STL file
    \param format is the format of the STL file
*/
STLBase::STLBase(const std::string &filename, Format format)
{
    setFilename(filename);
    setFormat(format);
}

/*!
    Get name of the STL file.

    \result The name of the STL file.
*/
const std::string & STLBase::getFilename() const
{
    return m_filename;
}

/*!
    Set the name of the STL file.

    \param filename is the name that will be set
*/
void STLBase::setFilename(const std::string &filename)
{
    m_filename = filename;
}

/*!
    Get format of the STL file.

    \result The format of the STL file.
*/
STLBase::Format STLBase::getFormat() const
{
    return m_format;
}

/*!
    Set the format of the STL file.

    \param format is the format that will be set
*/
void STLBase::setFormat(STLBase::Format format)
{
    m_format = format;
}

/*!
    \class STLReader
    \brief Class for reading ASCII and binary STL files.
*/

/*!
    Constructor.

    \param filename is the name of the STL file
    \param format is the format of the STL file
*/
STLReader::STLReader(const std::string &filename, Format format)
    : STLBase(filename)
{
    if (format == FormatUnknown) {
        format = detectFormat(filename);
    }

    if (format == FormatUnknown) {
        throw std::runtime_error("Invalid STL format.");
    }

    setFormat(format);
}

/*!
    Detects if the specified STL file is in binary format.

    \param filename is the name of the STL file
*/
STLReader::Format STLReader::detectFormat(const std::string &filename)
{
    std::ifstream fileStream;

    // Get the file size
    fileStream.open(filename, std::ifstream::ate | std::ifstream::binary);
    if (!fileStream.good()) {
        fileStream.close();
        throw std::runtime_error("Invalid STL.");
    }

    std::size_t fileSize = fileStream.tellg();
    fileStream.close();
    fileStream.clear();

    //
    // ASCII check
    //

    // Check if the size is compatible with an ASCII STL file.
    //
    // An ASCII contains at least the "solid " and "endsolid" markers, therefore
    // the minimum size of an empty ASCII file is 14 bytes.
    if (fileSize < ASCII_MINIMUM_SIZE) {
        return FormatUnknown;
    }

    // If a files starts with "solid" and ends with "endsolid" is an ASCII file.
    //
    // Binary files should never start with "solid ", but that's not mandatory.
    // We need to check both the beginning and the and of the file to be sure
    // of the file format.
    char c;
    std::size_t bufferPos;

    bufferPos = 0;
    std::string beginString(ASCII_FILE_BEGIN.size(), ' ');
    fileStream.open(filename, std::ifstream::binary);
    while (fileStream.get(c)) {
        if (bufferPos == 0 && (std::isblank(c) || std::isspace(c))) {
            continue;
        }

        beginString.at(bufferPos) = tolower(c);
        ++bufferPos;
        if (bufferPos == ASCII_FILE_BEGIN.size()) {
            break;
        }
    }
    fileStream.close();
    fileStream.clear();

    bool maybeASCII = (beginString.compare(ASCII_FILE_BEGIN) == 0);
    if (maybeASCII) {
        // Open the file
        fileStream.open(filename, std::ifstream::ate | std::ifstream::binary);

        // Move the cursor at the beginning of the last non-empty line
        bool empty = true;
        fileStream.seekg(-1, std::ios_base::cur);
        while (fileStream.get(c)) {
            if (c == '\n') {
                if (!empty) {
                    break;
                }
            }

            if (empty) {
                empty = (!std::isblank(c) && !std::isspace(c));
            }

            fileStream.seekg(-2, std::ios_base::cur);
        }

        // Search the end-line keyword
        bufferPos = 0;
        std::string endString(ASCII_SOLID_END.size(), ' ');
        while (fileStream.get(c)) {
            if (bufferPos == 0 && (std::isblank(c) || std::isspace(c))) {
                continue;
            }

            endString.at(bufferPos) = tolower(c);
            ++bufferPos;
            if (bufferPos == ASCII_SOLID_END.size()) {
                break;
            }
        }

        // Close the file
        fileStream.close();
        fileStream.clear();

        // Check if the end-solid keyword was found
        bool isASCII = (endString.compare(ASCII_SOLID_END) == 0);
        if (isASCII) {
            return FormatASCII;
        }
    }

    //
    // Binary check
    //

    // Check if the size is compatible with a binary STL file.
    //
    // An empty binary file contains the header and the number of facets,
    // therefore the minimum size of an empty binary file is 84 bytes.
    if (fileSize < BINARY_MINIMUM_SIZE) {
        return FormatUnknown;
    }

    // Read the number of facets
    std::uint32_t nFacets;

    fileStream.open(filename, std::ifstream::binary);
    fileStream.seekg(BINARY_HEADER_SIZE);
    fileStream.read(reinterpret_cast<char *>(&nFacets), sizeof(BINARY_UINT32));
    fileStream.close();
    fileStream.clear();

    // Check that the size of the file is compatiblewith the number of facets
    //
    // Each facet has three facet and each facet contains:
    //  - Normal: 3 float_32
    //  - Vertices' coordinates: 3x float_32
    //  - Attribute byte count: 1 unit_16
    const std::size_t BINARY_FACET_SIZE = 3 * sizeof(BINARY_REAL32) +
                                          3 * 3 * sizeof(BINARY_REAL32) +
                                          sizeof(BINARY_UINT16);

    std::size_t expectedFileSize = BINARY_HEADER_SIZE + sizeof(BINARY_UINT32) + (nFacets * BINARY_FACET_SIZE);
    if (fileSize == expectedFileSize) {
        return FormatBinary;
    }

    return FormatUnknown;
}

/*!
    Inspect the STL file.

    \param info on output it will contain the information gathered during
    the inspection
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: error opening the STL file
        - error = -2: the file is already open
*/
int STLReader::inspect(InspectionInfo *info)
{
    // Begin read
    int beginError = readBegin();
    if (beginError != 0) {
        return beginError;
    }

    // Get format
    Format format = getFormat();

    // Inspect file
    int inspectionError;
    if (format == FormatASCII) {
        inspectionError = inspectASCII(info);
    } else {
        inspectionError = inspectBinary(info);
    }

    // End read
    int endError = readEnd();

    // Done
    if (inspectionError != 0) {
        return inspectionError;
    } else if (endError != 0) {
        return endError;
    } else {
        return 0;
    }
}

/*!
    Inspect an ASCII STL file.

    \param[out] info on output it will contain the information gathered during
    the inspection
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to scan STL file
*/
int STLReader::inspectASCII(InspectionInfo *info)
{
    // Initialize info
    info->nSolids = 0;
    info->solidErrors.clear();
    info->solidNames.clear();
    info->solidFacetCount.clear();
    info->solidVertexCount.clear();

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Set cursor at file begin
    m_fileHandle.clear();
    long start_pos = m_fileHandle.tellg();
    m_fileHandle.seekg(0);

    // Scan file
    std::string line;
    std::string word;
    std::stringstream sline;
    std::stringstream sname;

    int inspectionError = 0;
    while (getline(m_fileHandle, line)) {
        line = utils::string::trim(line);
        sline.clear(),
        sline.str(line);

        // Get keyword
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            // Initialize inspection info
            int solidIndex = info->nSolids;

            ++info->nSolids;
            info->solidErrors.emplace_back();
            info->solidNames.emplace_back("");
            info->solidFacetCount.emplace_back(0);
            info->solidVertexCount.emplace_back(0);

            // Get solid name
            sname.str("");
            while (sline >> word) {
                sname << word << " ";
            }
            std::string name = sname.str();
            utils::string::trim(name);
            info->solidNames[solidIndex] = name;

            // Get solid info
            inspectionError = inspectSolidASCII(info->solidFacetCount.data() + solidIndex, info->solidErrors.data() + solidIndex);
            if (inspectionError != 0) {
                break;
            }

            info->solidVertexCount[solidIndex] = 3 * info->solidFacetCount[solidIndex];
        }
    }

    // Restore cursor position
    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    return inspectionError;
}

/*!
    Inspect solid data of an ASCII STL file.

    \param[out] nFacets on output it will contain the number of facets of
    the solid
    \param[out] errors bitmask for format errors encountered in the solid
    data section. The meaning of the error codes is the following:
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to inspect the STL file
*/
int STLReader::inspectSolidASCII(std::size_t *nFacets, std::array<bool, 6> *errors)
{
    // Initialize solid information
    (*nFacets) = 0;
    errors->fill(false);

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Inspect the solid
    std::string line;
    std::string word;
    std::stringstream sline;
    long last_valid_pos;

    int inspectError = 0;
    while (true) {
        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        getline(m_fileHandle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }

        // Exit conditions
        if (m_fileHandle.eof()) {
            break;
        } else if (word.compare(ASCII_SOLID_END) == 0) {
            break;
        } else if (word.compare(ASCII_SOLID_BEGIN) == 0) {
            break;
        }

        // Look for keyword "facet"
        if (word.compare(ASCII_FACET_BEGIN) == 0) {
            ++(*nFacets);

            m_fileHandle.seekg(last_valid_pos);
            inspectError = inspectFacetASCII(errors);
            if (inspectError != 0) {
                break;
            }
        }
    }

    // Check block temination
    if (word.compare(ASCII_SOLID_END) != 0) {
        (*errors)[0] = true;
        m_fileHandle.clear();
        m_fileHandle.seekg(last_valid_pos);
    }

    return inspectError;
}

/*!
    Inspect facet data of an ASCII STL file.

    \param[out] errors bitmask for format errors encountered in the facet data section.
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to scan/check STL file
*/
int STLReader::inspectFacetASCII(std::array<bool, 6> *errors)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Check facet data
    std::string line;
    std::string word;
    std::stringstream sline;
    long last_valid_pos;

    last_valid_pos = m_fileHandle.tellg();
    getline(m_fileHandle, line);
    line = utils::string::trim(line);
    sline.clear(),
    sline.str(line);
    if ((!(sline >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
        word = "begin";
    }

    std::size_t nV = 0;
    bool normal_found = false;
    while ((!m_fileHandle.eof())
           && ((word.compare(ASCII_FACET_END) != 0)
           &&  (word.compare(ASCII_FACET_BEGIN)    != 0)
           &&  (word.compare(ASCII_SOLID_END) != 0)
           &&  (word.compare(ASCII_SOLID_BEGIN)    != 0))) {

        // Check facet normal or facet vertices
        if (word.compare("begin") == 0) {
            if ((sline >> word) && (word.compare("normal") == 0)) {
                normal_found = true;

                int nxyz = 0;
                while (sline >> word) {
                    nxyz++;
                }

                if (nxyz != 3) {
                    (*errors)[3] = true;
                }
            }
        }
        else if (word.compare("vertex") == 0) {
            nV++;

            int nxyz = 0;
            while (sline >> word) {
                nxyz++;
            }

            if (nxyz != 3) {
                (*errors)[5] = true;
            }
        }

        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        getline(m_fileHandle, line);
        line = utils::string::trim(line);
        sline.clear(),
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }
    }

    // Check if facket section is properly closed
    if (word.compare(ASCII_FACET_END) != 0) {
        (*errors)[1] = true;
        m_fileHandle.clear(),
        m_fileHandle.seekg(last_valid_pos);
    }

    // Check if normal is valid
    if (!normal_found) {
        (*errors)[2] = true;
    }

    // Check if number of vertices is valid
    if (nV != 3) {
        (*errors)[4] = true;
    }

    return 0;
}

/*!
    Inspect a binary STL file.

    \param[out] info on output it will contain the information gathered during
    the inspection
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to scan STL file
*/
int STLReader::inspectBinary(InspectionInfo *info)
{
    // Initialize info
    info->nSolids = 0;
    info->solidErrors.clear();
    info->solidNames.clear();
    info->solidFacetCount.clear();
    info->solidVertexCount.clear();

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Set cursor at file begin
    m_fileHandle.clear();
    long start_pos = m_fileHandle.tellg();
    m_fileHandle.seekg(0);

    // Inspect header
    for (std::size_t i = 0; i < BINARY_HEADER_SIZE / sizeof(BINARY_UINT8); ++i) {
        BINARY_UINT8 headerCharacter;
        m_fileHandle.read(reinterpret_cast<char*>(&headerCharacter), sizeof(BINARY_UINT8));
    }

    if (m_fileHandle.eof()) {
        info->solidErrors[0][0] = true;
        return 0;
    }

    // Binary files does not contain information about the solid name
    info->solidNames.push_back("");

    // Read number of facets
    BINARY_UINT32 nFacets;
    m_fileHandle.read(reinterpret_cast<char *>(&nFacets), sizeof(BINARY_UINT32));
    info->solidFacetCount.push_back((std::size_t) nFacets);
    info->solidVertexCount.push_back((std::size_t) (3 * nFacets));

    if (m_fileHandle.eof()) {
        info->solidErrors[0][0] = true;
        return 0;
    }

    // Check facet data
    std::size_t n = 0;
    while ((!m_fileHandle.eof()) && (n < nFacets)) {
        // Read normal
        for (int j = 0; j < 3; ++j) {
            BINARY_REAL32 N_ij;
            m_fileHandle.read(reinterpret_cast<char *>(&N_ij), sizeof(BINARY_REAL32));
        }

        // Read vertex coordinates
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                BINARY_REAL32 V_ijk;
                m_fileHandle.read(reinterpret_cast<char *>(&V_ijk), sizeof(BINARY_REAL32));
            }
        }

        // Attribute byte count
        BINARY_UINT16 attributeByteCount;
        m_fileHandle.read(reinterpret_cast<char *>(&attributeByteCount), sizeof(BINARY_UINT16));

        // Update facet counter
        n++;
    }

    // Check number of facets
    if (n < nFacets) {
        info->solidErrors[0][1] = true;
    }

    // Reset cursor position
    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    return 0;
}

/*!
    Display inspection information.

    \param info are the information that will be displayed
    \param[in,out] out output stream
*/
void STLReader::displayInspectionInfo(const InspectionInfo &info, std::ostream &out) const
{
    out << "Inspection info" << std::endl;

    out << "  Filename : " << getFilename() << std::endl;

    Format format = getFormat();
    if (format == FormatBinary) {
        out << "  Format   : binary" << std::endl;
    } else {
        out << "  Format   : ASCII" << std::endl;
    }

    out << std::endl;
    if (info.nSolids > 0) {
        out << "  Solid count    : " << info.nSolids << std::endl;

        for (int i = 0; i < info.nSolids; ++i) {
            out << std::endl;
            out << "  Solid index    : " << i << std::endl;
            out << "  Solid name     : " << info.solidNames[i] << std::endl;
            out << "  Solid facets   : " << info.solidFacetCount[i] << std::endl;
            out << "  Solid vertices : " << info.solidVertexCount[i] << std::endl;

            if (info.solidErrors[i][0]) {
                out << "    **ERROR** Unterminated solid block." << std::endl;
            }
            if (info.solidErrors[i][1]) {
                out << "    **ERROR** Unterminated facet block." << std::endl;
            }
            if (info.solidErrors[i][2]) {
                out << "    **ERROR** Normal data are missing." << std::endl;
            }
            if (info.solidErrors[i][3]) {
                out << "    **ERROR** Wrong number of components for normal data." << std::endl;
            }
            if (info.solidErrors[i][4]) {
                out << "    **ERROR** Wrong number of vertices in facet block." << std::endl;
            }
            if (info.solidErrors[i][5]) {
                out << "    **ERROR** Wrong number of coordinates for vertice." << std::endl;
            }
        }
    } else {
        out << "  STL file contains no solids." << std::endl;
    }
}

/*!
    Begin reading the file.

    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: error opening the STL file
        - error = -2: the file is already open
*/
int STLReader::readBegin()
{
    if (m_fileHandle.is_open()) {
        return -2;
    }

    m_fileHandle.open(getFilename(), std::ifstream::in | std::ifstream::binary);
    if (!m_fileHandle.good()) {
        return -1;
    }

    return 0;
}

/*!
    Close the current stream to STL file.
*/
int STLReader::readEnd()
{
    m_fileHandle.close();

    return 0;
}

/*!
    Read solid data from the STL file.

    This routine assumes that the file stream is already open.

    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the STL file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the STL file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit std::vector to
    each facet acquired from the STL file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the STL file. New connectivity entries
    are appended at the end of T.
*/
int STLReader::readSolid(std::string *name, std::size_t *nV, std::size_t *nT,
                         std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                         std::vector<std::array<std::size_t, 3>> *T)
{
    Format format = getFormat();
    if (format == FormatASCII) {
        return readSolidASCII("", false, name, nV, nT, V, N, T);
    } else {
        return readSolidBinary(name, nV, nT, V, N, T);
    }
}

/*!
    Read solid data from the STL file.

    This routine assumes that the file stream is already open.

    \param[in] solid is the name of the solid that will to be read, if the
    name is empty, the first solid found will be read. Loading solids with
    a specific name is only supported for ASCII files
    \param[out] name on output will contain the name of the solid that has
    been actually read. Solid names are only supported for ASCII files
    \param[in,out] nV on input stores the current number of vertices hosted
    in V, on output stores the input values incremented by the number of
    vertices read from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T, on output stores the input value incremented by the
    number of facets read from the STL file.
    \param[in,out] V is the list of vertex coordinates, vertices read from the
    STL file are appended at the end of V
    \param[in,out] N is the list of facet normals, normals read from the STL
    file are appended at the end of V
    \param[in,out] T is the facet->vertex connectivity, entries read from the
    STL file are appended at the end of T
    \result Returns a negative number if an error occured, zero otherwise.
*/
int STLReader::readSolid(const std::string &solid, std::string *name, std::size_t *nV, std::size_t *nT,
                         std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                         std::vector<std::array<std::size_t, 3>> *T)
{
    Format format = getFormat();
    if (format == FormatASCII) {
        return readSolidASCII(solid, true, name, nV, nT, V, N, T);
    } else {
        std::string trimmedSolid = solid;
        utils::string::trim(trimmedSolid);
        if (!trimmedSolid.empty()) {
            log::cout() << "WARNING: loading solids with a specific name is only supported for ASCII files." << std::endl;
            log::cout() << "         The reader will read the next solid." << std::endl;
        }

        return readSolidBinary(name, nV, nT, V, N, T);
    }
}

/*!
    Read solid data from an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param[in] solid is the name of the solid that will to be read, if the
    name is empty, the first solid found will be read. Loading solids with
    a specific name is only supported for ASCII files
    \param[in] wrapAround controls if, when the end of file is reached, the
    search for the solid to read will begin again from the beginning of the
    file
    \param[out] name on output will contain the name of the solid that has
    been actually read. Solid names are only supported for ASCII files
    \param[in,out] nV on input stores the current number of vertices hosted
    in V, on output stores the input values incremented by the number of
    vertices read from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T, on output stores the input value incremented by the
    number of facets read from the STL file.
    \param[in,out] V is the list of vertex coordinates, vertices read from the
    STL file are appended at the end of V
    \param[in,out] N is the list of facet normals, normals read from the STL
    file are appended at the end of V
    \param[in,out] T is the facet->vertex connectivity, entries read from the
    STL file are appended at the end of T
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
        - error = -2: failed to read a facet
*/
int STLReader::readSolidASCII(const std::string &solid, bool wrapAround, std::string *name,
                              std::size_t *nV, std::size_t *nT, std::vector<std::array<double, 3>> *V,
                              std::vector<std::array<double, 3>> *N, std::vector<std::array<std::size_t, 3>> *T)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Get solid key
    std::string solidKey = solid;
    utils::string::trim(solidKey);
    solidKey = ASCII_SOLID_BEGIN + " " + solidKey;
    utils::string::trim(solidKey);

    // Scan file until solid is found
    std::string line;
    std::string word;
    std::stringstream sline;

    long start_pos   = m_fileHandle.tellg();
    long current_pos = start_pos + 1;

    bool check = false;
    while (!check && (start_pos != current_pos)) {
        // Get current line
        getline(m_fileHandle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Check end of file
        if (m_fileHandle.eof()) {
            if (wrapAround) {
                m_fileHandle.clear();
                m_fileHandle.seekg(0);
            } else {
                check = false;
                break;
            }
        }
        current_pos = m_fileHandle.tellg();

        // Look for keyword "solid"
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            if (solidKey.compare(ASCII_SOLID_BEGIN) == 0 || line.compare(solidKey) == 0) {
                *name = line.erase(0, ASCII_SOLID_BEGIN.size());
                *name = utils::string::trim(*name);

                start_pos = current_pos;
                check = true;
            }
        }
    }

    if (!check) {
        return 0;
    }

    // Read number of facets
    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    std::size_t nSolidFacets;
    std::array<bool, 6> solidErrors;
    inspectSolidASCII(&nSolidFacets, &solidErrors);

    // Read solid data
    std::array<double, 3> dummyDoubleArray;
    dummyDoubleArray.fill(0.);

    std::array<std::size_t, 3> dummyIntArray;
    dummyIntArray.fill(-1);

    V->resize(*nV + 3 * nSolidFacets, dummyDoubleArray);
    N->resize(*nT + nSolidFacets, dummyDoubleArray);
    T->resize(*nT + nSolidFacets, dummyIntArray);

    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);
    word = "begin";
    while ((!m_fileHandle.eof())
        && (word.compare(ASCII_SOLID_END) != 0)
        && (word.compare(ASCII_SOLID_BEGIN) != 0)) {

        // Get current line
        current_pos = m_fileHandle.tellg();
        getline(m_fileHandle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Look for keyword "facet"
        if ((sline >> word) && (word.compare(ASCII_FACET_BEGIN) == 0)) {
            m_fileHandle.seekg(current_pos);
            int readFacetError = readFacetASCII(nV, nT, V, N, T);
            if (readFacetError != 0) {
                return -2;
            }
        }
    }

    if (word.compare(ASCII_SOLID_END) != 0) {
        m_fileHandle.clear();
        m_fileHandle.seekg(current_pos);
    }

    return 0;
}

/*!
    Read facet data from ASCII STL file.

    This routine assumes that the file stream is already open.

    \param[in,out] nV on input stores the current number of vertices hosted
    in V, on output stores the input values incremented by the number of
    vertices read from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T, on output stores the input value incremented by the
    number of facets read from the STL file.
    \param[in,out] V is the list of vertex coordinates, vertices read from the
    STL file are appended at the end of V
    \param[in,out] N is the list of facet normals, normals read from the STL
    file are appended at the end of V
    \param[in,out] T is the facet->vertex connectivity, entries read from the
    STL file are appended at the end of T
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
        - error = -2: one facet has more than 3 vertices
*/
int STLReader::readFacetASCII(std::size_t *nV, std::size_t *nT, std::vector<std::array<double, 3>> *V,
                              std::vector<std::array<double, 3>> *N, std::vector<std::array<std::size_t, 3>> *T)
{
    // Read facet data
    std::string line;
    std::string word;
    std::stringstream sline;
    long last_valid_pos;

    last_valid_pos = m_fileHandle.tellg();
    getline(m_fileHandle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if ((!(sline >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
        word = "begin";
    }

    int nFacetVertices = 0;
    while ((!m_fileHandle.eof())
           && ((word.compare(ASCII_FACET_END) != 0)
           &&  (word.compare(ASCII_FACET_BEGIN)    != 0)
           &&  (word.compare(ASCII_SOLID_BEGIN)    != 0)
           &&  (word.compare(ASCII_SOLID_END) != 0))) {

        // Read facet normal or facet vertices
        if (word.compare("begin") == 0) {
            if ((sline >> word) && (word.compare("normal") == 0)) {
                sline >> (*N)[*nT];
            }
        }
        else if (word.compare("vertex") == 0) {
            if (nFacetVertices >= 3) {
                return 2;
            }

            sline >> (*V)[*nV + nFacetVertices];
            nFacetVertices++;
        }

        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        getline(m_fileHandle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }
    }

    // Restor cursor position
    if (word.compare(ASCII_FACET_END) != 0) {
        m_fileHandle.clear();
        m_fileHandle.seekg(last_valid_pos);
    }

    // Update facet-vertex connectivity
    for (int i = 0; i < nFacetVertices; ++i) {
        (*T)[*nT][i] = *nV + i;
    }

    // Update facet/vertex counters
    *nV += nFacetVertices;
    (*nT)++;

    return 0;
}

/*!
    Read solid data from a binary STL file.

    This routine assumes that the file stream is already open.

    \param[out] name since binary files have no information about solid names,
    on output it will contain an empty string
    \param[in,out] nV on input stores the current number of vertices hosted
    in V, on output stores the input values incremented by the number of
    vertices read from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T, on output stores the input value incremented by the
    number of facets read from the STL file.
    \param[in,out] V is the list of vertex coordinates, vertices read from the
    STL file are appended at the end of V
    \param[in,out] N is the list of facet normals, normals read from the STL
    file are appended at the end of V
    \param[in,out] T is the facet->vertex connectivity, entries read from the
    STL file are appended at the end of T
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
*/
int STLReader::readSolidBinary(std::string *name, std::size_t *nV, std::size_t *nT,
                               std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                               std::vector<std::array<std::size_t, 3>> *T)
{
    BITPIT_UNUSED(name);

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Set cursor at file begin
    m_fileHandle.clear();
    long start_pos = m_fileHandle.tellg();
    m_fileHandle.seekg(0);

    // Skip header
    for (std::size_t i = 0; i < BINARY_HEADER_SIZE / sizeof(BINARY_UINT8); ++i) {
        BINARY_UINT8 headerCharacter;
        m_fileHandle.read(reinterpret_cast<char *>(&headerCharacter), sizeof(BINARY_UINT8));
    }

    // Read number of facets
    BINARY_UINT32 nSolidFacets;
    m_fileHandle.read(reinterpret_cast<char *>(&nSolidFacets), sizeof(BINARY_UINT32));

    // Read data
    std::array<double, 3> dummyDoubleArray;
    dummyDoubleArray.fill(0.);

    std::array<std::size_t, 3> dummyIntArray;
    dummyIntArray.fill(-1);

    V->resize(*nV + 3 * nSolidFacets, dummyDoubleArray);
    N->resize(*nT + nSolidFacets, dummyDoubleArray);
    T->resize(*nT + nSolidFacets, dummyIntArray);

    for (std::size_t i = 0; i < nSolidFacets; ++i) {
        // Read normal
        for (int j = 0; j < 3; ++j) {
            BINARY_REAL32 N_ij;
            m_fileHandle.read(reinterpret_cast<char *>(&N_ij), sizeof(BINARY_REAL32));
            (*N)[i][j] = (double) N_ij;
        }

        // Read vertex coordinates
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                BINARY_REAL32 V_ijk;
                m_fileHandle.read(reinterpret_cast<char *>(&V_ijk), sizeof(BINARY_REAL32));
                (*V)[*nV][k] = (double) V_ijk;
            }

            (*nV)++;
        }

        // Facet-vertex connectivity
        (*T)[i][0] = *nV - 3;
        (*T)[i][1] = *nV - 2;
        (*T)[i][2] = *nV - 1;

        // Attribute byte count
        BINARY_UINT16 attributeByteCount;
        m_fileHandle.read(reinterpret_cast<char *>(&attributeByteCount), sizeof(BINARY_UINT16));
    }

    // Restore cursor position
    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    return 0;
}

/*!
    \class STLWriter
    \brief Class for writing ASCII and binary STL files.
*/

/*!
    Constructor.

    \param filename is the name of the STL file
    \param format is the format of the STL file
*/
STLWriter::STLWriter(const std::string &filename, Format format)
    : STLBase(filename, format)
{
    if (format == FormatUnknown) {
        throw std::runtime_error("Invalid STL format.");
    }
}

/*!
    Begin writing the file.

    \param writeMode is the write mode that will be used
    \result Returns a negative number if an error occured, zero otherwise.
*/
int STLWriter::writeBegin(WriteMode writeMode)
{
    if (m_fileHandle.is_open()) {
        return -2;
    }

    Format format = getFormat();

    std::ios_base::openmode openMode;
    if (format == FormatBinary) {
        if (writeMode != WriteOverwrite) {
            throw std::runtime_error("Specified write mode is not supported for binary files.");
        }

        openMode = std::ifstream::in | std::ifstream::binary;
    } else {
        if (writeMode == WriteOverwrite) {
            openMode = std::ifstream::in;
        } else if (writeMode == WriteAppend) {
            openMode = std::ifstream::app;
        } else {
            throw std::runtime_error("Specified write mode is not supported for ASCII files.");
        }

    }

    m_fileHandle.open(getFilename(), openMode);
    if (!m_fileHandle.good()) {
        return -1;
    }

    return 0;
}

/*!
    Close the current stream to STL file.
*/
int STLWriter::writeEnd()
{
    m_fileHandle.close();

    return 0;
}

/*!
    Write the specified solid data to the STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nV are the number of vertices of the solid
    \param nT are the number of facets of the solid
    \param V is the list of vertex coordinates
    \param N is the list of facet normals
    \param T is the facet->vertex connectivity
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
        - error = -2: input variable are not self-consistent
*/
int STLWriter::writeSolid(const std::string &name, std::size_t nV, std::size_t nT,
                          const std::vector<std::array<double,3>> &V, const std::vector<std::array<double,3>> &N,
                          const std::vector<std::array<std::size_t,3>> &T)
{
    Format format = getFormat();
    if (format == FormatASCII) {
        return writeSolidASCII(name, nV, nT, V, N, T);
    } else {
        return writeSolidBinary(name, nV, nT, V, N, T);
    }
}

/*!
    Write the specified solid data to an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nV are the number of vertices of the solid
    \param nT are the number of facets of the solid
    \param V is the list of vertex coordinates
    \param N is the list of facet normals
    \param T is the facet->vertex connectivity
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
        - error = -2: input variable are not self-consistent
*/
int STLWriter::writeSolidASCII(const std::string &name, std::size_t nV, std::size_t nT,
                               const std::vector<std::array<double, 3>> &V, const std::vector<std::array<double, 3>> &N,
                               const std::vector<std::array<std::size_t, 3>> &T)
{
    // Check input variable
    if (V.size() < nV) {
        return -2;
    }

    if (T.size() < nT) {
        return -2;
    }

    if (N.size() < nT) {
        return -2;
    }

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Save stream flags
    std::ios::fmtflags streamFlags(m_fileHandle.flags());

    // Write header
    std::stringstream sheader;
    sheader << ASCII_SOLID_BEGIN << " " << name;

    std::string header = sheader.str();
    header = utils::string::trim(header);
    m_fileHandle << header << std::endl;

    sheader.str("");

    // Write solid
    for (std::size_t i = 0; i < nT; ++i) {
        // Facet header
        m_fileHandle << "  " << ASCII_FACET_BEGIN;

        // Facet normal
        m_fileHandle << " normal ";
        m_fileHandle << std::scientific << N[i][0] << " ";
        m_fileHandle << std::scientific << N[i][1] << " ";
        m_fileHandle << std::scientific << N[i][2];;
        m_fileHandle << std::endl;

        // Facet vertices
        m_fileHandle << "    outer loop" << std::endl;
        for (std::size_t j = 0; j < 3; ++j) {
            m_fileHandle << "      vertex ";
            m_fileHandle << std::scientific << V[T[i][j]][0] << " ";
            m_fileHandle << std::scientific << V[T[i][j]][1] << " ";
            m_fileHandle << std::scientific << V[T[i][j]][2];
            m_fileHandle << std::endl;
        }
        m_fileHandle << "    endloop"    << std::endl;

        // Facet footer
        m_fileHandle << "  " << ASCII_FACET_END << std::endl;
    }

    // Solid footer
    std::stringstream sfooter;
    sfooter << ASCII_SOLID_END << " " << name;

    std::string footer;
    footer = sfooter.str();
    footer = utils::string::trim(footer);
    m_fileHandle << footer << std::endl;

    sfooter.str("");

    // Restore stream flags
    m_fileHandle.flags(streamFlags);

    return 0;
}

/*!
    Write the specified solid data to a binary STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nV are the number of vertices of the solid
    \param nT are the number of facets of the solid
    \param V is the list of vertex coordinates
    \param N is the list of facet normals
    \param T is the facet->vertex connectivity
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
        - error = -2: input variable are not self-consistent
*/
int STLWriter::writeSolidBinary(const std::string &name, std::size_t nV, std::size_t nT,
                                const std::vector<std::array<double, 3>> &V, const std::vector<std::array<double, 3>> &N,
                                const std::vector<std::array<std::size_t, 3>> &T)
{
    BITPIT_UNUSED(name);

    // Check input variable
    if (V.size() < nV) {
        return -2;
    }

    if (T.size() < nT) {
        return -2;
    }

    if (N.size() < nT) {
        return -2;
    }

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Write header
    for (std::size_t i = 0; i < BINARY_HEADER_SIZE / sizeof(BINARY_UINT8); ++i) {
        BINARY_UINT8 headerCharacter = 0;
        m_fileHandle.write(reinterpret_cast<char*>(&headerCharacter), sizeof(BINARY_UINT8));
    }

    // Write number of facets
    BINARY_UINT32 nFacets = (BINARY_UINT32) nT;
    m_fileHandle.write(reinterpret_cast<char *>(&nFacets), sizeof(BINARY_UINT32));

    // Write facets
    for (std::size_t i = 0; i < nT; ++i) {
        // Normals
        for (std::size_t j = 0; j < 3; ++j) {
            BINARY_REAL32 N_ij = (BINARY_REAL32) N[i][j];
            m_fileHandle.write(reinterpret_cast<char *>(&N_ij), sizeof(BINARY_REAL32));
        }

        // Vertex
        std::size_t T_size = T[i].size();
        for (std::size_t j = 0; j < T_size; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                BINARY_REAL32 V_ijk = (BINARY_REAL32) V[T[i][j]][k];
                m_fileHandle.write(reinterpret_cast<char *>(&V_ijk), sizeof(BINARY_REAL32));
            }
        }

        // Attribute byte count
        //
        // In the standard format, this should be zero because most software
        // does not understand anything else.
        // Attribute byte count
        BINARY_UINT16 attributeByteCount = 0;
        m_fileHandle.write(reinterpret_cast<char *>(&attributeByteCount), sizeof(BINARY_UINT16));
    }

    return 0;
}

}
