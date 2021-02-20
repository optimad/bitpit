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

#include <cassert>
#include <cstring>

#include "bitpit_common.hpp"
#include "bitpit_operators.hpp"

#include "STL.hpp"
#include "logger.hpp"

namespace bitpit {

/*!
    \class STLBase
    \brief Base class for the STL writer and the STL reader.
*/

// Header is 80 characters long
const std::size_t STLBase::BINARY_HEADER_SIZE = 80 * sizeof(STLBase::BINARY_UINT8);

// An empty binary file contains the header plus a float_32 stating that there
// are zero trinagles in the file.
const std::size_t STLBase::BINARY_MINIMUM_SIZE = STLBase::BINARY_HEADER_SIZE + sizeof(STLBase::BINARY_UINT32);

// Each facet is defined by the following information:
//  - Normal: 3 float_32
//  - Vertices' coordinates: 3x float_32
//  - Attribute byte count: 1 unit_16
const std::size_t STLBase::BINARY_FACET_SIZE = 3 * sizeof(BINARY_REAL32) + 3 * 3 * sizeof(BINARY_REAL32) + sizeof(BINARY_UINT16);

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
    std::string word;
    std::stringstream sname;

    int inspectionError = 0;
    while (m_lineStream.readLine(m_fileHandle) >= 0) {
        // Get keyword
        if ((m_lineStream >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            // Initialize inspection info
            int solidIndex = info->nSolids;

            ++info->nSolids;
            info->solidErrors.emplace_back();
            info->solidNames.emplace_back("");
            info->solidFacetCount.emplace_back(0);
            info->solidVertexCount.emplace_back(0);

            // Get solid name
            sname.str("");
            while (m_lineStream >> word) {
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
    std::string word;
    long last_valid_pos;

    int inspectError = 0;
    while (true) {
        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        m_lineStream.readLine(m_fileHandle);
        if (!(m_lineStream >> word)) {
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
    std::string word;
    long last_valid_pos;

    last_valid_pos = m_fileHandle.tellg();
    m_lineStream.readLine(m_fileHandle);
    if ((!(m_lineStream >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
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
            if ((m_lineStream >> word) && (word.compare("normal") == 0)) {
                normal_found = true;

                int nxyz = 0;
                while (m_lineStream >> word) {
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
            while (m_lineStream >> word) {
                nxyz++;
            }

            if (nxyz != 3) {
                (*errors)[5] = true;
            }
        }

        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        m_lineStream.readLine(m_fileHandle);
        if (!(m_lineStream >> word)) {
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
        std::array<char, BINARY_FACET_SIZE> facetData;
        m_fileHandle.read(facetData.data(), BINARY_FACET_SIZE);
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
    out << "Inspection info" << "\n";

    out << "  Filename : " << getFilename() << "\n";

    Format format = getFormat();
    if (format == FormatBinary) {
        out << "  Format   : binary" << "\n";
    } else {
        out << "  Format   : ASCII" << "\n";
    }

    out << "\n";
    if (info.nSolids > 0) {
        out << "  Solid count    : " << info.nSolids << "\n";

        for (int i = 0; i < info.nSolids; ++i) {
            out << "\n";
            out << "  Solid index    : " << i << "\n";
            out << "  Solid name     : " << info.solidNames[i] << "\n";
            out << "  Solid facets   : " << info.solidFacetCount[i] << "\n";
            out << "  Solid vertices : " << info.solidVertexCount[i] << "\n";

            if (info.solidErrors[i][0]) {
                out << "    **ERROR** Unterminated solid block." << "\n";
            }
            if (info.solidErrors[i][1]) {
                out << "    **ERROR** Unterminated facet block." << "\n";
            }
            if (info.solidErrors[i][2]) {
                out << "    **ERROR** Normal data are missing." << "\n";
            }
            if (info.solidErrors[i][3]) {
                out << "    **ERROR** Wrong number of components for normal data." << "\n";
            }
            if (info.solidErrors[i][4]) {
                out << "    **ERROR** Wrong number of vertices in facet block." << "\n";
            }
            if (info.solidErrors[i][5]) {
                out << "    **ERROR** Wrong number of coordinates for vertice." << "\n";
            }
        }
    } else {
        out << "  STL file contains no solids." << "\n";
    }

    out.flush();
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
int STLReader::readSolid(std::string *name, std::size_t *nV, std::size_t *nT,
                         std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                         std::vector<std::array<std::size_t, 3>> *T)
{
    return readSolid("", name, nV, nT, V, N, T);
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
    // Read header
    std::size_t nSolidFacets;
    int headerError = readHeader(solid, name, &nSolidFacets);
    if (headerError != 0) {
        return headerError;
    }

    // Read facet data
    V->resize(*nV + 3 * nSolidFacets, {{0., 0., 0.}});
    N->resize(*nT + nSolidFacets, {{0., 0., 0.}});
    T->resize(*nT + nSolidFacets, {{0, 0, 0}});

    for (std::size_t i = 0; i < *nT; ++i) {
        // Read facet data
        std::array<double, 3> *V0 = V->data() + *nV + 3 * i;
        std::array<double, 3> *V1 = V0 + 1;
        std::array<double, 3> *V2 = V1 + 1;

        int facetError = readFacet(V0, V1, V2, N->data() + i);
        if (facetError != 0) {
            return facetError;
        }

        // Update facet->vertex connectivity
        (*T)[i][0] = *nV + 3 * i;
        (*T)[i][1] = (*T)[i][0] + 1;
        (*T)[i][2] = (*T)[i][1] + 1;
    }

    // Read footer
    int footerError = readFooter(solid);
    if (footerError != 0) {
        return footerError;
    }

    return 0;
}

/*!
    Read the header of the STL file.

    This routine assumes that the file stream is already open.

    \param[out] name on output will contain the name of the solid that has
    been actually read
    \param[out] nT on input will contain the number of factes of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: header was not found
*/
int STLReader::readHeader(std::string *name, std::size_t *nT)
{
    return readHeader("", name, nT);
}

/*!
    Read the header of the STL file.

    This routine assumes that the file stream is already open.

    \param[in] solid is the name of the solid, if the name is empty, the
    first header will be read
    \param[out] name on output will contain the name of the solid that has
    been actually read
    \param[out] nT on input will contain the number of factes of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: header was not found
*/
int STLReader::readHeader(const std::string &solid, std::string *name, std::size_t *nT)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = readHeaderASCII(solid, name, nT);
    } else {
        std::string trimmedSolid = solid;
        utils::string::trim(trimmedSolid);
        if (!trimmedSolid.empty()) {
            log::cout() << "WARNING: loading solids with a specific name is only supported for ASCII files." << std::endl;
            log::cout() << "         The reader will read the next solid." << std::endl;
        }

        error = readHeaderBinary(name, nT);
    }

    return error;
}

/*!
    Read the footer of an ASCII STL file.

    This routine assumes that the file stream is already open.

    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: footer was not found
*/
int STLReader::readFooter()
{
    return readFooter("");
}

/*!
    Read the footer of an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param solid is the name of the solid, if the name is empty, the first
    footer will be read
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: footer was not found
*/
int STLReader::readFooter(const std::string &solid)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = readFooterASCII(solid);
    } else {
        error = readFooterBinary();
    }

    return error;
}

/*!
    Read facet data from an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param[out] V0 on output will contain the coordinates of the first vertex
    \param[out] V1 on output will contain the coordinates of the second vertex
    \param[out] V2 on output will contain the coordinates of the third vertex
    \param[out] N on output will contain the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
        - error = -2: facet data section is malformed
        - error = -3: the facet has more than 3 vertices
*/
int STLReader::readFacet(std::array<double, 3> *V0, std::array<double, 3> *V1,
                         std::array<double, 3> *V2, std::array<double, 3> *N)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = readFacetASCII(V0, V1, V2, N);
    } else {
        error = readFacetBinary(V0, V1, V2, N);
    }

    return error;
}

/*!
    Read the header of an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param[in] solid is the name of the solid, if the name is empty, the
    first header will be read
    \param[out] name on output will contain the name of the solid that has
    been actually read
    \param[out] nT on input will contain the number of factes of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: header was not found
*/
int STLReader::readHeaderASCII(const std::string &solid, std::string *name, std::size_t *nT)
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
    std::string word;
    std::string line;

    long start_pos   = m_fileHandle.tellg();
    long current_pos = start_pos + 1;

    bool solidFound = false;
    bool wrapAround = solidKey.compare(ASCII_SOLID_BEGIN) != 0;
    while (!solidFound && (start_pos != current_pos)) {
        // Get current line
        m_lineStream.readLine(m_fileHandle);

        // Check end of file
        if (m_fileHandle.eof()) {
            if (wrapAround) {
                m_fileHandle.clear();
                m_fileHandle.seekg(0);
                wrapAround = false;
            } else {
                solidFound = false;
                break;
            }
        }
        current_pos = m_fileHandle.tellg();

        // Look for keyword "solid"
        if ((m_lineStream >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            m_lineStream.copyLine(&line);
            if (solidKey.compare(ASCII_SOLID_BEGIN) == 0 || line.compare(solidKey) == 0) {
                *name = line.erase(0, ASCII_SOLID_BEGIN.size());
                *name = utils::string::trim(*name);

                start_pos = current_pos;

                solidFound = true;
            }
        }
    }

    if (!solidFound) {
        return -2;
    }

    // Read number of facets
    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    std::array<bool, 6> solidErrors;
    inspectSolidASCII(nT, &solidErrors);

    m_fileHandle.clear();
    m_fileHandle.seekg(start_pos);

    return 0;
}

/*!
    Read the footer of an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param solid is the name of the solid, if the name is empty, the first
    footer will be read
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: footer was not found
*/
int STLReader::readFooterASCII(const std::string &solid)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Get solid key
    std::string solidKey = solid;
    utils::string::trim(solidKey);
    solidKey = ASCII_SOLID_END + " " + solidKey;
    utils::string::trim(solidKey);

    // Look for the end of solid section
    std::string word;
    std::string line;

    while (true) {
        // Get next line
        m_lineStream.readLine(m_fileHandle);

        // Get next word
        if (!(m_lineStream >> word)) {
            word = "";
        }

        // Handle the word
        if (word.compare(ASCII_SOLID_END) == 0) {
            m_lineStream.copyLine(&line);
            if (line.compare(solidKey) != 0) {
                log::cout() << "WARNING: end-solid key does not match the solid name." << std::endl;
                log::cout() << "         Expected end-solid key : " << solidKey << std::endl;
                log::cout() << "         Current end-solid key  : " << line << std::endl;
            }

            break;
        } else if (word.compare(ASCII_SOLID_BEGIN) == 0) {
            return -2;
        } else if (m_fileHandle.eof()) {
            return -2;
        }
    }

    return 0;
}

/*!
    Read facet data from an ASCII STL file.

    This routine assumes that the file stream is already open.

    \param[out] V0 on output will contain the coordinates of the first vertex
    \param[out] V1 on output will contain the coordinates of the second vertex
    \param[out] V2 on output will contain the coordinates of the third vertex
    \param[out] N on output will contain the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
        - error = -2: facet data section is malformed
        - error = -3: the facet has more than 3 vertices
*/
int STLReader::readFacetASCII(std::array<double, 3> *V0, std::array<double, 3> *V1,
                              std::array<double, 3> *V2, std::array<double, 3> *N)
{
    // Read facet data
    std::string word;
    std::string value;
    long last_valid_pos;

    int error = 0;
    int nFacetVertices = 0;
    std::string target = "facet";
    while (true) {
        // Get next line
        last_valid_pos = m_fileHandle.tellg();
        m_lineStream.readLine(m_fileHandle);
        if (!(m_lineStream >> word)) {
            continue;
        }

        // Handle the word
        //
        // Beginning of faacet section and normal are on the same line, after
        // finding the beginning of the facet section we have to extract the
        // next word without getting a new line.
        if (word.compare(ASCII_FACET_BEGIN) == 0) {
            if (word.compare(target) == 0) {
                target = "normal";

                if (!(m_lineStream >> word)) {
                    continue;
                }
            } else {
                error = -2;
                break;
            }
        }

        if (word.compare(ASCII_FACET_END) == 0) {
            if (word.compare(target) == 0) {
                break;
            } else {
                error = -2;
                break;
            }
        } else if (word.compare("normal") == 0) {
            if (word.compare(target) == 0) {
                for (int k = 0; k < 3; ++k) {
                    m_lineStream >> value;
                    (*N)[k] = stod(value);
                }
                target = "vertex";
            } else {
                error = -2;
                break;
            }
        } else if (word.compare("vertex") == 0) {
            if (word.compare(target) == 0) {
                std::array<double, 3> *coords;
                if (nFacetVertices == 0) {
                    coords = V0;
                } else if (nFacetVertices == 1) {
                    coords = V1;
                } else if (nFacetVertices == 2) {
                    coords = V2;
                    target = ASCII_FACET_END;
                } else {
                    error = -3;
                    break;
                }

                for (int k = 0; k < 3; ++k) {
                    m_lineStream >> value;
                    (*coords)[k] = stod(value);
                }

                nFacetVertices++;
            } else {
                error = -2;
                break;
            }
        } else if (word.compare(ASCII_SOLID_BEGIN) == 0) {
            error = -2;
            break;
        } else if (word.compare(ASCII_SOLID_END) == 0) {
            error = -2;
            break;
        } else if (m_fileHandle.eof()) {
            error = -2;
            break;
        }
    }

    // Restor cursor position
    if (error != 0) {
        m_fileHandle.clear();
        m_fileHandle.seekg(last_valid_pos);
    }

    return error;
}

/*!
    Read the header of a binary STL file.

    This routine assumes that the file stream is already open.

    \param[out] name since binary files have no information about solid names,
    on output it will contain an empty string
    \param[out] nT on input will contain the number of factes of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
        - error = -2: unable to read the header
*/
int STLReader::readHeaderBinary(std::string *name, std::size_t *nT)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Binary data has not information about solid names
    (*name) = "";

    // Skip header
    for (std::size_t i = 0; i < BINARY_HEADER_SIZE / sizeof(BINARY_UINT8); ++i) {
        BINARY_UINT8 headerCharacter;
        m_fileHandle.read(reinterpret_cast<char *>(&headerCharacter), sizeof(BINARY_UINT8));
    }

    // Read number of facets
    BINARY_UINT32 nSolidFacets;
    m_fileHandle.read(reinterpret_cast<char *>(&nSolidFacets), sizeof(BINARY_UINT32));
    *nT = nSolidFacets;

    // Check if the end of file has been reached
    if (m_fileHandle.eof()) {
        return -2;
    }

    return 0;
}

/*!
    Read the footer of a binary STL file.

    This routine assumes that the file stream is already open.

    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read data from output stream
*/
int STLReader::readFooterBinary()
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Nothing to do
    return 0;
}

/*!
    Read facet data from a binary STL file.

    This routine assumes that the file stream is already open.

    \param[out] V0 on output will contain the coordinates of the first vertex
    \param[out] V1 on output will contain the coordinates of the second vertex
    \param[out] V2 on output will contain the coordinates of the third vertex
    \param[out] N on output will contain the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to read from input stream
        - error = -2: unable to read facet data
*/
int STLReader::readFacetBinary(std::array<double, 3> *V0, std::array<double, 3> *V1,
                               std::array<double, 3> *V2, std::array<double, 3> *N)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Read facet data
    std::array<char, BINARY_FACET_SIZE> facetData;
    m_fileHandle.read(facetData.data(), BINARY_FACET_SIZE);
    int facetDataOffset = 0;

    // Store normal
    for (int k = 0; k < 3; ++k) {
        (*N)[k] = (double) *(reinterpret_cast<BINARY_REAL32 *>(facetData.data() + facetDataOffset));
        facetDataOffset += sizeof(BINARY_REAL32);
    }

    // Store vertex coordinates
    for (int k = 0; k < 3; ++k) {
        (*V0)[k] = (double) *(reinterpret_cast<BINARY_REAL32 *>(facetData.data() + facetDataOffset));
        facetDataOffset += sizeof(BINARY_REAL32);
    }

    for (int k = 0; k < 3; ++k) {
        (*V1)[k] = (double) *(reinterpret_cast<BINARY_REAL32 *>(facetData.data() + facetDataOffset));
        facetDataOffset += sizeof(BINARY_REAL32);
    }

    for (int k = 0; k < 3; ++k) {
        (*V2)[k] = (double) *(reinterpret_cast<BINARY_REAL32 *>(facetData.data() + facetDataOffset));
        facetDataOffset += sizeof(BINARY_REAL32);
    }

    // Check if the end of file has been reached
    if (m_fileHandle.eof()) {
        return -2;
    }

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
    \param partialWrite tells the writer that only part of the facets will
    be written, this can be used for writing a binary file incrementally
    \result Returns a negative number if an error occured, zero otherwise.
*/
int STLWriter::writeBegin(WriteMode writeMode, bool partialWrite)
{
    if (m_fileHandle.is_open()) {
        return -2;
    }

    Format format = getFormat();

    std::ios_base::openmode openMode;
    if (format == FormatBinary) {
        if (writeMode == WriteOverwrite) {
            openMode = std::ofstream::out | std::ofstream::binary;
        } else if (partialWrite && writeMode == WriteAppend) {
            openMode = std::ofstream::app | std::ofstream::binary;
        } else {
            throw std::runtime_error("Specified write mode is not supported for binary files.");
        }
    } else {
        if (writeMode == WriteOverwrite) {
            openMode = std::ofstream::out;
        } else if (writeMode == WriteAppend) {
            openMode = std::ofstream::app;
        } else {
            throw std::runtime_error("Specified write mode is not supported for ASCII files.");
        }
    }

    m_fileHandle.open(getFilename(), openMode);
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Set stream properties
    m_fileHandle << std::scientific;

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
    // Check input variables
    if (V.size() < nV) {
        return -2;
    }

    if (T.size() < nT) {
        return -2;
    }

    if (N.size() < nT) {
        return -2;
    }

    // Save stream flags
    std::ios::fmtflags streamFlags(m_fileHandle.flags());

    // Write header
    int headerError = writeHeader(name, nT);
    if (headerError != 0) {
        return headerError;
    }

    // Write facet data
    for (std::size_t i = 0; i < nT; ++i) {
        // Check connectivity
        for (int k = 0; k < 3; ++k) {
            if (T[i][k] >= nV) {
                return -2;
            }
        }

        // Write data
        int facetError = writeFacet(V[T[i][0]], V[T[i][1]], V[T[i][2]], N[i]);
        if (facetError != 0) {
            return facetError;
        }
    }

    // Write footer
    int footerError = writeFooter(name);
    if (footerError != 0) {
        return footerError;
    }

    // Restore stream flags
    m_fileHandle.flags(streamFlags);

    return 0;
}

/*!
    Write the header to the STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nT are the number of facets of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeHeader(const std::string &name, std::size_t nT)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = writeHeaderASCII(name, nT);
    } else {
        error = writeHeaderBinary(name, nT);
    }

    return error;
}

/*!
    Write the footer to the STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFooter(const std::string &name)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = writeFooterASCII(name);
    } else {
        error = writeFooterBinary(name);
    }

    return error;
}

/*!
    Write the specified facet data to the STL file.

    This routine assumes that the file stream is already open.

    \param V0 are the coordinates of the first vertex
    \param V1 are the coordinates of the second vertex
    \param V2 are the coordinates of the third vertex
    \param N is the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFacet(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                          const std::array<double, 3> &V2, const std::array<double, 3> &N)
{
    Format format = getFormat();

    int error;
    if (format == FormatASCII) {
        error = writeFacetASCII(V0, V1, V2, N);
    } else {
        error = writeFacetBinary(V0, V1, V2, N);
    }

    return error;
}

/*!
    Write the header of a binary STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nT are the number of facets of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeHeaderASCII(const std::string &name, std::size_t nT)
{
    BITPIT_UNUSED(nT);

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Write header
    std::stringstream sheader;
    sheader << ASCII_SOLID_BEGIN << " " << name;

    // Write number of facets
    std::string header = sheader.str();
    header = utils::string::trim(header);
    m_fileHandle << header << "\n";

    return 0;
}

/*!
    Write the footer of a binary STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFooterASCII(const std::string &name)
{
    BITPIT_UNUSED(name);

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Write footer
    std::stringstream sfooter;
    sfooter << ASCII_SOLID_END << " " << name;

    std::string footer;
    footer = sfooter.str();
    footer = utils::string::trim(footer);
    m_fileHandle << footer << "\n";

    return 0;
}

/*!
    Write the specified facet data to a binary STL file.

    This routine assumes that the file stream is already open.

    \param V0 are the coordinates of the first vertex
    \param V1 are the coordinates of the second vertex
    \param V2 are the coordinates of the third vertex
    \param N is the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFacetASCII(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                               const std::array<double, 3> &V2, const std::array<double, 3> &N)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Facet header
    m_fileHandle << "  " << ASCII_FACET_BEGIN;

    // Facet normal
    m_fileHandle << " normal ";
    m_fileHandle << N[0] << " ";
    m_fileHandle << N[1] << " ";
    m_fileHandle << N[2];
    m_fileHandle << "\n";

    // Facet vertices
    m_fileHandle << "    outer loop" << "\n";

    m_fileHandle << "      vertex ";
    m_fileHandle << V0[0] << " ";
    m_fileHandle << V0[1] << " ";
    m_fileHandle << V0[2];
    m_fileHandle << "\n";

    m_fileHandle << "      vertex ";
    m_fileHandle << V1[0] << " ";
    m_fileHandle << V1[1] << " ";
    m_fileHandle << V1[2];
    m_fileHandle << "\n";

    m_fileHandle << "      vertex ";
    m_fileHandle << V2[0] << " ";
    m_fileHandle << V2[1] << " ";
    m_fileHandle << V2[2];
    m_fileHandle << "\n";

    m_fileHandle << "    endloop"    << "\n";

    // Facet footer
    m_fileHandle << "  " << ASCII_FACET_END << "\n";

    return 0;
}

/*!
    Write the header of a binary STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \param nT are the number of facets of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeHeaderBinary(const std::string &name, std::size_t nT)
{
    BITPIT_UNUSED(name);

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

    return 0;
}

/*!
    Write the footer of a binary STL file.

    This routine assumes that the file stream is already open.

    \param name is the name of the solid
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFooterBinary(const std::string &name)
{
    BITPIT_UNUSED(name);

    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Nothing to do
    return 0;
}

/*!
    Write the specified facet data to a binary STL file.

    This routine assumes that the file stream is already open.

    \param V0 are the coordinates of the first vertex
    \param V1 are the coordinates of the second vertex
    \param V2 are the coordinates of the third vertex
    \param N is the normal
    \result Returns a negative number if an error occured, zero otherwise.
    The meaning of the error codes is the following:
        - error = -1: failed to write data to output stream
*/
int STLWriter::writeFacetBinary(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                                const std::array<double, 3> &V2, const std::array<double, 3> &N)
{
    // Check stream status
    if (!m_fileHandle.good()) {
        return -1;
    }

    // Normals
    for (int k = 0; k < 3; ++k) {
        BINARY_REAL32 N_k = (BINARY_REAL32) N[k];
        m_fileHandle.write(reinterpret_cast<char *>(&N_k), sizeof(BINARY_REAL32));
    }

    // Vertices
    for (int k = 0; k < 3; ++k) {
        BINARY_REAL32 V_k = (BINARY_REAL32) V0[k];
        m_fileHandle.write(reinterpret_cast<char *>(&V_k), sizeof(BINARY_REAL32));
    }

    for (int k = 0; k < 3; ++k) {
        BINARY_REAL32 V_k = (BINARY_REAL32) V1[k];
        m_fileHandle.write(reinterpret_cast<char *>(&V_k), sizeof(BINARY_REAL32));
    }

    for (int k = 0; k < 3; ++k) {
        BINARY_REAL32 V_k = (BINARY_REAL32) V2[k];
        m_fileHandle.write(reinterpret_cast<char *>(&V_k), sizeof(BINARY_REAL32));
    }

    // Attribute byte count
    //
    // In the standard format, this should be zero because most software
    // does not understand anything else.
    // Attribute byte count
    BINARY_UINT16 attributeByteCount = 0;
    m_fileHandle.write(reinterpret_cast<char *>(&attributeByteCount), sizeof(BINARY_UINT16));

    return 0;
}

}
