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

#include "bitpit_common.hpp"
#include "bitpit_operators.hpp"

#include "STL.hpp"
#include "logger.hpp"

namespace bitpit {

/*!
    \class STLObj
    \brief Interface to STL I/O function

    This class has been designed to allow an easy interface between end-user
    and STL I/O functions.
*/

const std::size_t STLObj::BINARY_HEADER_SIZE  = 80 * sizeof(STLObj::BINARY_UINT8);
const std::size_t STLObj::BINARY_MINIMUM_SIZE = STLObj::BINARY_HEADER_SIZE + sizeof(STLObj::BINARY_UINT32);

const std::string STLObj::ASCII_SOLID_BEGIN  = "solid";
const std::string STLObj::ASCII_SOLID_END    = "endsolid";
const std::string STLObj::ASCII_FACET_BEGIN  = "facet";
const std::string STLObj::ASCII_FACET_END    = "endfacet";
const std::string STLObj::ASCII_FILE_BEGIN   = STLObj::ASCII_SOLID_BEGIN + " ";
const std::string STLObj::ASCII_FILE_END     = STLObj::ASCII_SOLID_END;
const std::size_t STLObj::ASCII_MINIMUM_SIZE = STLObj::ASCII_FILE_BEGIN.length() + STLObj::ASCII_FILE_END.length();

/*!
    Default constructor for class STLObj.

    Initialize an empty interface to STL file.
*/
STLObj::STLObj()
{
    // General info
    stl_name = "";
    stl_type = false;

    // Error flags
    err = 0;

    // stl content
    data.n_solids = -1;
}

/*!
    Constructor #1 for class STLObj.
    Initialize an interface to STL file with name specified in filename.

    \param[in] filename STL file name
    \param[in] filetype boolean flag for ASCII (false) or binary (true) STL file
*/
STLObj::STLObj(std::string filename, bool filetype)
{
    // General info
    stl_name = utils::string::trim(filename);
    stl_type = filetype;

    // Error flags
    err = 0;

    // stl content
    data.n_solids = -1;
}

/*!
    Constructor #1 for class STLObj.
    Initialize an interface to STL file with name specified in filename.

    \param[in] filename STL file name
*/
STLObj::STLObj(std::string filename)
{
    // General info
    stl_name = utils::string::trim(filename);

    // Error flags
    err = 0;

    // STL content
    data.n_solids = -1;

    // Detect file format (ASCII or BINARY)
    Format fileFormat = detectFormat(stl_name);
    if (fileFormat == FormatUnknown) {
        throw std::runtime_error("Invalid STL.");
    }

    stl_type = (fileFormat == FormatBinary);
}

/*!
    Detects if the specified STL file is in binary format.

    \param[in] filename STL file name
*/
STLObj::Format STLObj::detectFormat(const std::string &filename)
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
    fileStream.read(reinterpret_cast<char*>(&nFacets), sizeof(BINARY_UINT32));
    fileStream.close();
    fileStream.clear();

    // Check that the size of the file is compatiblewith the number of facets
    //
    // Each facet contains the following information:
    //  - normal: 3 float_32;
    //  - vertices' coordinates: 3 float_32;
    //  - attribute byte count: 1 unit_16.
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
    Open the file associated to the interface.

    \param[in] mode opening mode ("in": input, "out": output, "app": append mode)
*/
void STLObj::open(const std::string &mode)
{
    // Open streams
    if (mode.compare("in") == 0) {
        if (!m_ifile_handle.is_open()) {
            // Close output stream
            close("out");

            // Open input stream
            m_ifile_handle.open(stl_name, std::ifstream::in | std::ifstream::binary);

            // Check stream status
            if (!m_ifile_handle.good()) {
                err = 1;
            } else {
                err = 0;
            }
        }
    } else if (mode.compare("out") == 0) {
        if (!m_ofile_handle.is_open()) {
            // Close input stream
            close("in");

            // Clear infos & flags
            clear();

            // Open output stream
            if (stl_type) {
                m_ofile_handle.open(stl_name, std::ifstream::out | std::ifstream::binary);
            } else {
                m_ofile_handle.open(stl_name, std::ifstream::out);
            }
        }

        // Check stream status
        if (!m_ofile_handle.good()) {
            err = 1;
        } else {
            err = 0;
        }
    } else if (mode.compare("app") == 0) {
        if (!m_ofile_handle.is_open()) {

            // Close input stream
            close("in");

            // Open output stream
            if (stl_type) { m_ofile_handle.open(stl_name, std::ifstream::app | std::ifstream::binary); }
            else          { m_ofile_handle.open(stl_name, std::ifstream::app); }
        }

        // Check stream status
        if (!m_ofile_handle.good()) {
            err = 1;
        } else {
            err = 0;
        }
    }
}

/*!
    Close the current stream to STL file.

    \param[in] mode opening mode used to open stream ("in": input, "out": output,
    "app": append)
*/
void STLObj::close(const std::string &mode)
{
    // Close streams
    if (mode.compare("in") == 0) {
        m_ifile_handle.close();
    } else if (mode.compare("out") == 0) {
        m_ofile_handle.close();
    } else if (mode.compare("app") == 0) {
        m_ofile_handle.close();
    } else if (mode.compare("") == 0) {
        close("in");
        close("out");
    }
}

/*!
    Clear info and error flags gathered on the STL file associated
    to the interface.
*/
void STLObj::clear(void)
{
    // Close
    close();

    // Reset error flags
    err = 0;
    stl_errors.resize(0);

    // Reset STL content
    data.n_solids = -1;
    data.solid_names.resize(0);
    data.solid_facets.resize(0);
}

/*!
    Display info gathered on the STL file associated to the interface.

    \param[in,out] out output stream
*/
void STLObj::display(std::ostream &out)
{
    // General info
    out << "STL object:" << std::endl;
    out << "name          : " << stl_name << std::endl;
    if (stl_type) {
        out << "file type         : binary" << std::endl;
    } else {
        out << "file type         : ASCII"  << std::endl;
    }
    if (m_ifile_handle.is_open()) {
        out << "input stream      : open" << std::endl;
    } else {
        out << "input stream      : closed" << std::endl;
    }
    if (m_ofile_handle.is_open()) {
        out << "output stream     : open" << std::endl;
    } else {
        out << "output stream     : closed" << std::endl;
    }

    // Error flags
    if (err == 1) {
        out << "**ERROR** STL file is missing!" << std::endl;
        return;
    }

    // STL content
    if (data.n_solids >= 0) {
        out << "# of stl solids   : " << data.n_solids << std::endl;
        if (data.n_solids > 0) {
            out << "solid's name      : " << data.solid_names << std::endl;
            out << "# of facet        : " << data.solid_facets << std::endl;
        }
    }

    // Errors
    if ((data.n_solids >= 0) && (stl_errors.size() > 0)) {
        out << "ERROR report:" << std::endl;
        for (int i = 0; i < data.n_solids; ++i) {
            out << "  stl solid: '" << data.solid_names[i] << "'" << std::endl;
            if (stl_errors[i][0]) {
                out << "    **WARNING** unterminated solid block!!" << std::endl;
            }
            if (stl_errors[i][1]) {
                out << "    **WARNING** unterminated facet block!!" << std::endl;
            }
            if (stl_errors[i][2]) {
                out << "    **WARNING** normal data are missing!!" << std::endl;
            }
            if (stl_errors[i][3]) {
                out << "    **WARNING** wrong number of components for normal data!!" << std::endl;
            }
            if (stl_errors[i][4]) {
                out << "    **WARNING** wrong number of vertices in facet block!!" << std::endl;
            }
            if (stl_errors[i][5]) {
                out << "    **WARNING** wrong number of coordinates for vertices!!" << std::endl;
            }
        }
    }
}

/*!
    Scan and gather info from the STL file associated to the interface
*/
void STLObj::scan()
{
    // Celar info
    clear();

    // Open input stream
    open("in");

    // Scan file
    if (stl_type) {
        err = STLObj::scanBINARY(m_ifile_handle, data.solid_names, data.solid_facets);
        data.n_solids = data.solid_names.size();
    } else {
        err = STLObj::scanASCII(m_ifile_handle, data.solid_names, data.solid_facets);
        data.n_solids = data.solid_names.size();
    }

    // Close inut stream
    close("in");
}

/*!
    Scan and check for format errors in the STL file associated to the interface
*/
void STLObj::check()
{
    // Scan file
    scan();

    // Open input stream
    open("in");

    // Check file
    if (stl_type) {
        STLObj::checkBINARY(m_ifile_handle, stl_errors);
    } else {
        STLObj::checkASCII(m_ifile_handle, stl_errors);
    }

    // Close input stream
    close("in");
}

/*!
    Load solid data from the STL file associated to the interface.

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
void STLObj::load(std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                  std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T)
{
    // Open stream
    open("in");
    if (err == 1) {
        return;
    }

    // Read STL data
    if (stl_type) {
        STLObj::readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        STLObj::readASCII(m_ifile_handle, nV, nT, V, N, T);
    }

    // Close stream
    close("in");
}

/*!
    Load solid data from the STL file associated to the interface. Overloading of
    member function STLObj::load() for container std::vector<array<double, 3>>

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
void STLObj::load(std::size_t &nV, std::size_t &nT,
                  std::vector<std::array<double, 3>> &V,
                  std::vector<std::array<double, 3>> &N,
                  std::vector<std::array<std::size_t, 3>> &T)
{
    // Open stream
    open("in");
    if (err == 1) {
        return;
    }

    // Read STL data
    if (stl_type) {
        STLObj::readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        STLObj::readASCII(m_ifile_handle, nV, nT, V, N, T);
    }

    // Close stream
    close("in");
}

/*!
    Load single solid data from the STL file associated to the interface.
    This routine assumes that the input stream is already open.

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
    \param[in,out] name name of the solid that will to be read, if the name is
    empty, the first solid found will be read. On output il will contain the name
    of the solid that has been actually read
*/
void STLObj::loadSolid(std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                       std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                       std::string &name)
{
    // Check input data
    if (stl_type && !utils::string::trim(name).empty()) {
        log::cout() << "WARNING: loading solids with a specific name is only supported for ASCII files" << std::endl;
    }

    // Read STL data
    if (stl_type) {
        STLObj::readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        STLObj::readSolidASCII(m_ifile_handle, false, nV, nT, V, N, T, name);
    }
}

/*!
    Load single solid solid data from the STL file associated to the interface.
    This routine assumes that the input stream is already open.
    Overloading of member function STLObj::load() for container
    std::vector<array<double, 3>>

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
    \param[in,out] name name of the solid that will to be read, if the name is
    empty, the first solid found will be read. On output il will contain the name
    of the solid that has been actually read
*/
void STLObj::loadSolid(std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                       std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T,
                       std::string &name)
{
    // Check inupt data
    if (stl_type && !utils::string::trim(name).empty()) {
        log::cout() << "WARNING: loading solids with a specific name is only supported for ASCII files" << std::endl;
    }

    // Read STL data
    if (stl_type) {
        STLObj::readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        STLObj::readSolidASCII(m_ifile_handle, false, nV, nT, V, N, T, name);
    }
}

/*!
    Save single solid data appending to the the STL file associated to the interface.
    This routine assumes that the output stream is already open, in append mode.

    \param[in]  name name associated to the solid.
    \param[in]  nV number of vertices of current solid.
    \param[in]  nT number of facets of current solid.
    \param[in]  V vertex coordinates list.
    \param[in]  N facet normals.
    \param[in]  T facet->vertex connectivity.
 */
void STLObj::saveSolid(const std::string &name, std::size_t &nV, std::size_t &nT,  std::vector<std::vector<double>> &V,
                       std::vector<std::vector<double>> &N,  std::vector<std::vector<std::size_t>> &T)
{
    if (stl_type) {
        STLObj::writeSolidBINARY(m_ofile_handle, nV, nT, V, N, T, name);
    } else {
        STLObj::writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }
}

/*!
*   Save single solid data appending to the the STL file associated to the interface.
*   Overloading for vertices and facet normals passed as std::vector< std::array< double, 3>>.
*   This routine assumes that the output stream is already open, in append mode.
*
*    \param[in]  name label associated to the solid.
*    \param[in]  nV number of vertices of current solid.
*    \param[in]  nT number of facets of current solid.
*    \param[in]  V vertex coordinates list.
*    \param[in]  N facet normals.
*    \param[in]  T facet->vertex connectivity.
*/
void STLObj::saveSolid(const std::string &name, std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                       std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T)
{
    if (stl_type) {
        STLObj::writeSolidBINARY(m_ofile_handle, nV, nT, V, N, T, name);
    } else {
        STLObj::writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }
}

/*!
    Dummy function for recursive variadic template STLObj::save().
*/
void STLObj::save()
{
    // Close output stream
    close("out");
}

/*!
    Dummy function for recursive variadic template STLObj::load().
*/
void STLObj::load()
{
    // Close input stream
    close("in");
}

/*!
    Scan ASCII STL file and retrieve infos about content.

    \param[in,out] file_handle input stream from STL file
    \param[out] solid_names list of solid names stored in the STL file
    \param[out] solid_facets number of facets for each solid found in the
    STL file (same ordering of solid_names).

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to sca STL file.
*/
unsigned int STLObj::scanASCII(std::ifstream &file_handle, std::vector<std::string> &solid_names,
                               std::vector<std::size_t> &solid_facets)
{
    // Initialize output variables
    solid_names.resize(0);
    solid_facets.resize(0);

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Scan file
    std::string line;
    std::string word;
    std::stringstream sline;
    std::stringstream sword;

    while (!file_handle.eof()) {
        // Get current line
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Get keyword
        if (sline >> word) {
            if (word.compare(ASCII_SOLID_BEGIN) == 0) {
                // Get solid name
                sword.str("");
                while (sline >> word) {
                    sword << word << " ";
                }
                word = sword.str();
                solid_names.push_back(utils::string::trim(word));

                // Get solid info
                std::size_t nF = 0;
                scanSolidASCII(file_handle, nF);
                solid_facets.push_back(nF);

            }
        }
    }

    // Restore cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Scan binary STL file and retrieve infos about content.

    \param[in,out] file_handle input stream from STL file
    \param[out] solid_names list of solid names stored in the STL file
    \param[out] solid_facets number of facets for each solid found in the
    STL file (same ordering of solid_names).

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to sca STL file.
*/
unsigned int STLObj::scanBINARY(std::ifstream &file_handle, std::vector<std::string> &solid_names,
                                std::vector<std::size_t> &solid_facets)
{
    // Resize output variables
    solid_names.resize(0);
    solid_facets.resize(0);

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Define data types
    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Get solid name
    for (int i = 0; i < 20; ++i) {
        file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
    }
    solid_names.push_back("");

    // Read number of facets
    file_handle.read(reinterpret_cast<char*>(longint4byte_ptr), 4);
    solid_facets.push_back((int)*longint4byte_ptr);

    // Reset cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Scan solid data from ASCII STL file and retrieve infos.

    \param[in,out] file_handle input stream from STL file
    \param[in,out] nT number of facets found in the solid data section

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to sca STL file.
*/
unsigned int STLObj::scanSolidASCII(std::ifstream &file_handle, std::size_t &nT)
{
    // Resize output variables
    nT = 0;

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Scan file
    std::string line;
    std::string word;
    std::stringstream sline;
    long int last_valid_pos;

    last_valid_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if (!(sline >> word)) {
        word = "";
    }

    while ((!file_handle.eof())
        && ((word.compare(ASCII_SOLID_END)   != 0)
        &&  (word.compare(ASCII_SOLID_BEGIN) != 0))) {

        // Look for keyword "facet"
        if (word.compare(ASCII_FACET_BEGIN) == 0) {
            nT++;
        }

        // Get next line
        last_valid_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }
    }

    if (word.compare(ASCII_SOLID_END) != 0) {
        file_handle.clear();
        file_handle.seekg(last_valid_pos);
    }

    return 0;
}

/*!
    Scan ASCII STL file and perform check on format error.

    \param[in] file_handle input stream from STL file.
    \param[out] err_map for each solid found in the STL file, store a bitmask
    for format errors encountered in the solid data section.
    error code:
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check STL file
*/
unsigned int STLObj::checkASCII(std::ifstream &file_handle, std::vector<std::vector<bool>> &err_map)
{
    // Resize output variables
    err_map.resize(0);

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Scan file
    std::string line;
    std::string word;
    std::stringstream sline;

    while (getline(file_handle, line)) {
        line = utils::string::trim(line);
        sline.clear(),
        sline.str(line);

        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            std::vector<bool> _map(6, false);
            checkSolidASCII(file_handle, _map);
            err_map.push_back(_map);
        }

    }

    // Reset cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Scan solid data section in ASCII STL file and perform check on format error.

    \param[in] file_handle input stream from STL file.
    \param[out] err_map bitmask for format errors encountered in the solid data section.
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check STL file
*/
unsigned int STLObj::checkSolidASCII(std::ifstream &file_handle, std::vector<bool> &err_map)
{
    // Resize output variables
    err_map.resize(6, false);

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Check file
    std::string line;
    std::string word;
    std::stringstream sline;
    long int last_valid_pos;

    while (true) {
        // Get next line
        last_valid_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }

        // Exit conditions
        if (file_handle.eof()) {
            break;
        } else if (word.compare(ASCII_SOLID_END) == 0) {
            break;
        } else if (word.compare(ASCII_SOLID_BEGIN) != 0) {
            break;
        }

        // Look for keyword "facet"
        if (word.compare(ASCII_FACET_BEGIN) == 0) {
            file_handle.seekg(last_valid_pos);
            checkFacetASCII(file_handle, err_map);
        }
    }

    // Check block temination
    if (word.compare(ASCII_SOLID_END) != 0) {
        err_map[0] = true;
        file_handle.clear();
        file_handle.seekg(last_valid_pos);
    }

    return 0;
}

/*!
    Scan facet data section in ASCII STL file and perform check on format error.

    \param[in] file_handle input stream from STL file.
    \param[out] err_map bitmask for format errors encountered in the facet data section.
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check STL file
*/
unsigned int STLObj::checkFacetASCII(std::ifstream &file_handle, std::vector<bool> &err_map)
{
    // Initialize output variables
    err_map.resize(6, false);

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Check facet data
    std::string line;
    std::string word;
    std::stringstream sline;
    long int last_valid_pos;

    last_valid_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear(),
    sline.str(line);
    if ((!(sline >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
        word = "begin";
    }

    int nV = 0;
    bool normal_found = false;
    while ((!file_handle.eof())
           && ((word.compare(ASCII_FACET_END)   != 0)
           &&  (word.compare(ASCII_FACET_BEGIN) != 0)
           &&  (word.compare(ASCII_SOLID_END)   != 0)
           &&  (word.compare(ASCII_SOLID_BEGIN) != 0))) {

        // Check facet normal or facet vertices
        if (word.compare("begin") == 0) {
            if ((sline >> word) && (word.compare("normal") == 0)) {
                normal_found = true;

                int nxyz = 0;
                while (sline >> word) {
                    nxyz++;
                }

                if (nxyz != 3) {
                    err_map[3] = true;
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
                err_map[5] = true;
            }
        }

        // Get next line
        last_valid_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear(),
        sline.str(line);
        if (!(sline >> word)) { word = ""; }
    }

    // Check if facket section is properly closed
    if (word.compare(ASCII_FACET_END) != 0) {
        err_map[1] = true;
        file_handle.clear(),
        file_handle.seekg(last_valid_pos);
    }

    // Check if normal is valid
    if (!normal_found) {
        err_map[2] = true;
    }

    // Check if number of vertices is valid
    if (nV != 3) {
        err_map[4] = true;
    }

    return 0;
}

/*!
    Scan binary STL file and perform check on format error.

    \param[in] file_handle input stream from STL file.
    \param[out] err_map for each solid found in the STL file, store a bitmask
    for format errors encountered in the solid data section.
    error code:
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check STL file
*/
unsigned int STLObj::checkBINARY(std::ifstream &file_handle, std::vector<std::vector<bool>> &err_map)
{
    // Resize output variables
    err_map.resize(1, std::vector<bool>(6, false));

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Define data types
    unsigned int  int2byte;
    unsigned int *int2byte_ptr = &int2byte;

    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Check file header
    for (int i = 0; i < 20; ++i) {
        file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
    }

    if (file_handle.eof()) {
        err_map[0][0] = true;
        return 0;
    }

    // Check number of facets
    file_handle.read(reinterpret_cast<char*>(longint4byte_ptr), 4);
    std::size_t nT = (int) *longint4byte_ptr;

    if (file_handle.eof()) {
        err_map[0][0] = true;
        return 0;
    }

    // Check facet data
    std::size_t n = 0;
    while ((!file_handle.eof()) && (n < nT)) {
        // Read normal
        for (int j = 0; j < 3; ++j) {
            file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
        }

        // Read vertex coordinates
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
            }
        }

        // Facet closing header
        file_handle.read(reinterpret_cast<char*>(int2byte_ptr), 2);

        // Update facet counter
        n++;
    }

    // Check number of facets
    if (n < nT) {
        err_map[0][1] = true;
    }

    return 0;
}

/*!
    Read solid data from ASCII STL file.

    \param[in,out] file_handle stream from STL file
    \param[in] wrapAround controls if, when the end of file is reached, the
    search for the solid to read will begin again from the beginning of the
    file.
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the solid.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the solid.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the solid. Vertices acquired from solid are appended at the end
    of V.
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the solid. Normals acquired from solid are appended at the end
    of N.
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T.
    \param[in,out] name name of the solid that will to be read, if the name is
    empty, the first solid found will be read. On output il will contain the name
    of the solid that has been actually read

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: failed to read a facet
*/
unsigned int STLObj::readSolidASCII(std::ifstream &file_handle, bool wrapAround,
                                    std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                                    std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                                    std::string &name)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Get solid key
    name = utils::string::trim(name);

    std::string name_key = ASCII_SOLID_BEGIN;
    if (!name.empty()) {
        name_key += " " + name;
        name_key = utils::string::trim(name_key);
    }

    // Scan file until solid is found
    std::string line;
    std::string word;
    std::stringstream sline;

    long int start_pos   = file_handle.tellg();
    long int current_pos = start_pos + 1;

    bool check = false;
    while (!check && (start_pos != current_pos)) {
        // Get current line
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Check end of file
        if (file_handle.eof()) {
            if (wrapAround) {
                file_handle.clear();
                file_handle.seekg(0);
            } else {
                check = false;
                break;
            }
        }
        current_pos = file_handle.tellg();

        // Look for keyword "solid"
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            if (name.empty() || line.compare(name_key) == 0) {
                name = line.erase(0, ASCII_SOLID_BEGIN.size());
                name = utils::string::trim(name);

                start_pos = current_pos;
                check = true;
            }
        }
    }

    if (!check) {
        return 0;
    }

    // Read number of facets
    std::size_t nt;

    file_handle.clear();
    file_handle.seekg(start_pos);
    scanSolidASCII(file_handle, nt);

    // Read solid data
    V.resize(nV+3*nt, std::vector<double>(3, 0.0));
    N.resize(nT+nt, std::vector<double>(3, 0.0));
    T.resize(nT+nt, std::vector<std::size_t>(3, -1));

    file_handle.clear();
    file_handle.seekg(start_pos);
    word = "begin";
    while ((!file_handle.eof())
        && (word.compare(ASCII_SOLID_END)   != 0)
        && (word.compare(ASCII_SOLID_BEGIN) != 0)) {

        // Get current line
        current_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Look for keyword "facet"
        if ((sline >> word) && (word.compare(ASCII_FACET_BEGIN) == 0)) {
            file_handle.seekg(current_pos);
            int readFacetError = readFacetASCII(file_handle, nV, nT, V, N, T);
            if (readFacetError != 0) {
                return 2;
            }
        }
    }

    if (word.compare(ASCII_SOLID_END) != 0) {
        file_handle.clear();
        file_handle.seekg(current_pos);
    }

    return 0;
}

/*!
    Read solid data from ASCII STL file. Overloading of readSolidASCII()
    for std::vector<array<double, 3>> container.

    \param[in,out] file_handle stream from STL file
    \param[in] wrapAround controls if, when the end of file is reached, the
    search for the solid to read will begin again from the beginning of the
    file.
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the solid.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the solid.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the solid. Vertices acquired from solid are appended at the end
    of V.
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the solid. Normals acquired from solid are appended at the end
    of N.
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T.
    \param[in,out] name name of the solid that will to be read, if the name is
    empty, the first solid found will be read. On output il will contain the name
    of the solid that has been actually read

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: failed to read a facet
*/
unsigned int STLObj::readSolidASCII(std::ifstream &file_handle, bool wrapAround,
                                    std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                                    std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T,
                                    std::string &name)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Get solid key
    name = utils::string::trim(name);

    std::string name_key = ASCII_SOLID_BEGIN;
    if (!name.empty()) {
        name_key += " " + name;
        name_key = utils::string::trim(name_key);
    }

    // Scan file until solid is found
    std::string line;
    std::string word;
    std::stringstream sline;

    long int start_pos   = file_handle.tellg();
    long int current_pos = start_pos + 1;

    bool check = false;
    while (!check && (start_pos != current_pos)) {
        // Get current line
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Check end of file
        if (file_handle.eof()) {
            if (wrapAround) {
                file_handle.clear();
                file_handle.seekg(0);
            } else {
                check = false;
                break;
            }
        }
        current_pos = file_handle.tellg();

        // Look for keyword "solid"
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            if (name.empty() || line.compare(name_key) == 0) {
                name = line.erase(0, ASCII_SOLID_BEGIN.size());
                name = utils::string::trim(name);

                start_pos = current_pos;
                check = true;
            }
        }
    }

    if (!check) {
        return 0;
    }

    // Read number of facets
    std::size_t nt;

    file_handle.clear();
    file_handle.seekg(start_pos);
    scanSolidASCII(file_handle, nt);

    // Read solid data
    std::array<double, 3> dummyDoubleArray;
    dummyDoubleArray.fill(0.);

    std::array<std::size_t, 3> dummyIntArray;
    dummyIntArray.fill(-1);

    V.resize(nV+3*nt, dummyDoubleArray);
    N.resize(nT+nt, dummyDoubleArray);
    T.resize(nT+nt, dummyIntArray);

    file_handle.clear();
    file_handle.seekg(start_pos);
    word = "begin";
    while ((!file_handle.eof())
        && (word.compare(ASCII_SOLID_END)   != 0)
        && (word.compare(ASCII_SOLID_BEGIN) != 0)) {

        // Get current line
        current_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);

        // Look for keyword "facet"
        if ((sline >> word) && (word.compare(ASCII_FACET_BEGIN) == 0)) {
            file_handle.seekg(current_pos);
            int readFacetError = readFacetASCII(file_handle, nV, nT, V, N, T);
            if (readFacetError != 0) {
                return 2;
            }
        }
    }

    if (word.compare(ASCII_SOLID_END) != 0) {
        file_handle.clear();
        file_handle.seekg(current_pos);
    }

    return 0;
}

/*!
    Read data from ASCII STL file.

    \param[in,out] file_handle stream from STL file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the STL file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the STL file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the STL file. Vertices acquired from the STL file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the STL file. Normals acquired from the STL file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same STL file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int STLObj::readASCII(std::ifstream &file_handle,
                               std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                               std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Read solids
    std::string line;
    std::string word;
    std::string name = "";
    std::stringstream sline;

    long int current_pos = file_handle.tellg();
    while (getline(file_handle, line)) {
        sline.clear();
        sline.str(line);
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            file_handle.clear();
            file_handle.seekg(current_pos);
            readSolidASCII(file_handle, true, nV, nT, V, N, T, name);
        }
        current_pos = file_handle.tellg();
    }

    // Restore cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Read data from ASCII STL file. Overloading of readASCII()
    for container std::vector<array<double, 3>>

    \param[in,out] file_handle stream from STL file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the STL file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the STL file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the STL file. Vertices acquired from the STL file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the STL file. Normals acquired from the STL file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same STL file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int STLObj::readASCII(std::ifstream &file_handle,
                               std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                               std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Read solids
    std::string line;
    std::string word;
    std::string name = "";
    std::stringstream sline;

    long int current_pos = file_handle.tellg();
    while (getline(file_handle, line)) {
        sline.clear();
        sline.str(line);
        if ((sline >> word) && (word.compare(ASCII_SOLID_BEGIN) == 0)) {
            file_handle.clear();
            file_handle.seekg(current_pos);
            readSolidASCII(file_handle, true, nV, nT, V, N, T, name);
        }
        current_pos = file_handle.tellg();
    }

    // Restore cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}


/*!
    Read STL facet data from ASCII STL file.

    \param[in,out] file_handle input stream from STL file
    \param[in,out] nV on input stores the number of vertices previously acquired from
    the STL file, on output stores the input value incremented by the number
    of vertices read from the facet section.
    \param[in,out] nT on input stores the number of facets previously aquired from
    the STL file, on output stores the input value incremented by 1.
    \param[in,out] V vertex coordinate list. Vertices read from facet section
    are stored in the location V[nV], V[nV+1], etc. nV is the number of vertices
    previously acquired from the STL file and passed as input to this function.
    \param[in,out] N normal list. Normal read from facet section is stored at
    position N[nT], where nT is the number of facet previously acquired from the
    STL file and passed as input to this function.
    \param[in,out] T vertex->facet connectivity. A new connectivity entry for the facet
    begin read is created at position T[nT], where nT is the number of facet previously
    acquired from the STL file and passed as input to this function.

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: one facet has more than 3 vertices
*/
unsigned int STLObj::readFacetASCII(std::ifstream &file_handle,
                                    std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                                    std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T)
{
    // Read facet data
    std::string line;
    std::string word;
    std::stringstream sline;
    long int last_valid_pos;

    last_valid_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if ((!(sline >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
        word = "begin";
    }

    int nv = 0;
    while ((!file_handle.eof())
           && ((word.compare(ASCII_FACET_END)   != 0)
           &&  (word.compare(ASCII_FACET_BEGIN) != 0)
           &&  (word.compare(ASCII_SOLID_BEGIN) != 0)
           &&  (word.compare(ASCII_SOLID_END)   != 0))) {

        // Read facet normal or facet vertices
        if (word.compare("begin") == 0) {
            if ((sline >> word) && (word.compare("normal") == 0)) {
                sline >> N[nT];
            }
        }
        else if (word.compare("vertex") == 0) {
            if (nv >= 3) {
                return(2);
            }

            sline >> V[nV+nv];
            nv++;
        }

        // Get next line
        last_valid_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }
    }

    // Restor cursor position
    if (word.compare(ASCII_FACET_END) != 0) {
        file_handle.clear();
        file_handle.seekg(last_valid_pos);
    }

    // Update triangle-vertex connectivity
    for (int i = 0; i < nv; ++i) {
        T[nT][i] = nV + i;
    }

    // Update facet/vertex counters
    nV += 3;
    nT++;

    return 0;
}

/*!
    Read STL facet data from ASCII STL file.
    Overloading of readFacetASCII() for std::vector<array<double, 3>> container.

    \param[in,out] file_handle input stream from STL file
    \param[in,out] nV on input stores the number of vertices previously acquired from
    the STL file, on output stores the input value incremented by the number
    of vertices read from the facet section.
    \param[in,out] nT on input stores the number of facets previously aquired from
    the STL file, on output stores the input value incremented by 1.
    \param[in,out] V vertex coordinate list. Vertices read from facet section
    are stored in the location V[nV], V[nV+1], etc. nV is the number of vertices
    previously acquired from the STL file and passed as input to this function.
    \param[in,out] N normal list. Normal read from facet section is stored at
    position N[nT], where nT is the number of facet previously acquired from the
    STL file and passed as input to this function.
    \param[in,out] T vertex->facet connectivity. A new connectivity entry for the facet
    begin read is created at position T[nT], where nT is the number of facet previously
    acquired from the STL file and passed as input to this function.

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: one facet has more than 3 vertices
*/
unsigned int STLObj::readFacetASCII(std::ifstream &file_handle,
                                    std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                                    std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T)
{
    // Read facet data
    std::string line;
    std::string word;
    std::stringstream sline;
    long int last_valid_pos;

    last_valid_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if ((!(sline >> word)) || (word.compare(ASCII_FACET_BEGIN) == 0)) {
        word = "begin";
    }

    int nv = 0;
    while ((!file_handle.eof())
           && ((word.compare(ASCII_FACET_END)   != 0)
           &&  (word.compare(ASCII_FACET_BEGIN) != 0)
           &&  (word.compare(ASCII_SOLID_BEGIN) != 0)
           &&  (word.compare(ASCII_SOLID_END)   != 0))) {

        // Read facet normal or facet vertices
        if (word.compare("begin") == 0) {
            if ((sline >> word) && (word.compare("normal") == 0)) {
                sline >> N[nT];
            }
        }
        else if (word.compare("vertex") == 0) {
            if (nv >= 3) {
                return(2);
            }

            sline >> V[nV+nv];
            nv++;
        }

        // Get next line
        last_valid_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) {
            word = "";
        }
    }

    // Restor cursor position
    if (word.compare(ASCII_FACET_END) != 0) {
        file_handle.clear();
        file_handle.seekg(last_valid_pos);
    }

    // Update triangle-vertex connectivity
    for (int i = 0; i < nv; ++i) {
        T[nT][i] = nV + i;
    }

    // Update facet/vertex counters
    nV += 3;
    nT++;

    return 0;
}

/*!
    Read data from binary  STL file.

    \param[in,out] file_handle stream from STL file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the STL file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the STL file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the STL file. Vertices acquired from the STL file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the STL file. Normals acquired from the STL file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same STL file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int STLObj::readBINARY(std::ifstream &file_handle,
                                std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                                std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Define data types
    unsigned int  int2byte;
    unsigned int *int2byte_ptr = &int2byte;

    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Skip header
    for (int i = 0; i < 20; ++i) {
        file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
    }

    // Read number of elements
    file_handle.read(reinterpret_cast<char*>(longint4byte_ptr), 4);
    nT = (std::size_t)*longint4byte_ptr;
    nV = 3*nT;

    // Read data
    V.resize(nV, std::vector<double>(3, 0.0));
    N.resize(nT, std::vector<double>(3, 0.0));
    T.resize(nT, std::vector<std::size_t>(3, -1));

    nV = 0;
    for (std::size_t i = 0; i < nT; ++i) {
        // Read normal
        for (int j = 0; j < 3; ++j) {
            file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
            N[i][j] = (double)*float4byte_ptr;
        }

        // Read vertex coordinates
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
                V[nV][k] = (double)*float4byte_ptr;
            }

            nV++;
        }

        // Triangle-vertex connectivity
        T[i][0] = nV - 3;
        T[i][1] = nV - 2;
        T[i][2] = nV - 1;

        // Facet closing header
        file_handle.read(reinterpret_cast<char*>(int2byte_ptr), 2);
    }

    // Restore cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Read data from binary STL file. Overloading of readBINARY() for container
    std::vector<array<double, 3>>

    \param[in,out] file_handle stream from STL file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the STL file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the STL file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the STL file. Vertices acquired from the STL file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] N facet's normals. On output stores the normal unit std::vector
    to each facet acquired from the STL file. Normals acquired from the STL file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same STL file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same STL file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int STLObj::readBINARY(std::ifstream &file_handle,
                                std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                                std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T)
{
    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Set cursor at file begin
    file_handle.clear();
    long int start_pos = file_handle.tellg();
    file_handle.seekg(0);

    // Define data types
    unsigned int  int2byte;
    unsigned int *int2byte_ptr = &int2byte;

    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Skip header
    for (int i = 0; i < 20; ++i) {
        file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
    }

    // Read number of elements
    file_handle.read(reinterpret_cast<char*>(longint4byte_ptr), 4);
    nT = (std::size_t)*longint4byte_ptr;
    nV = 3*nT;

    // Read data
    std::array<double, 3> dummyDoubleArray;
    dummyDoubleArray.fill(0.);

    std::array<std::size_t, 3> dummyIntArray;
    dummyIntArray.fill(-1);

    V.resize(nV, dummyDoubleArray);
    N.resize(nT, dummyDoubleArray);
    T.resize(nT, dummyIntArray);

    nV = 0;
    for (std::size_t i = 0; i < nT; ++i) {
        // Read normal
        for (int j = 0; j < 3; ++j) {
            file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
            N[i][j] = (double)*float4byte_ptr;
        }

        // Read vertex coordinates
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                file_handle.read(reinterpret_cast<char*>(float4byte_ptr), 4);
                V[nV][k] = (double)*float4byte_ptr;
            }
            nV++;
        }

        // Triangle-vertex connectivity
        T[i][0] = nV - 3;
        T[i][1] = nV - 2;
        T[i][2] = nV - 1;

        // Facet closing header
        file_handle.read(reinterpret_cast<char*>(int2byte_ptr), 2);
    }

    // Restore cursor position
    file_handle.clear();
    file_handle.seekg(start_pos);

    return 0;
}

/*!
    Write solid data to ASCII STL file.

    \param[in,out] file_handle stream to STL file
    \param[in,out] nV number of vertices to be written to STL file.
    \param[in,out] nT number of facet to be written to STL file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int STLObj::writeSolidASCII(std::ofstream &file_handle,
                                     std::size_t &nV, std::size_t &nT,  std::vector<std::vector<double>> &V,
                                     std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                                     const std::string &solid_name)
{
    // Check input variable
    if (V.size() < static_cast<std::size_t>(nV)) {
        return 2;
    }

    if (T.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    if (N.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Save stream flags
    std::ios::fmtflags streamFlags(file_handle.flags());

    // Write header
    std::stringstream sheader;
    sheader << ASCII_SOLID_BEGIN << " " << solid_name;

    std::string header = sheader.str();
    header = utils::string::trim(header);
    file_handle << header << std::endl;

    sheader.str("");

    // Write solid
    for (std::size_t i = 0; i < nT; ++i) {
        // Facet header
        file_handle << "  " << ASCII_FACET_BEGIN;

        // Facet normal
        std::size_t n_size = N[i].size();
        if (n_size > 0) {
            file_handle << " normal ";
            for (std::size_t j = 0; j < n_size-1; ++j) {
                file_handle << std::scientific << N[i][j] << " ";
            }
            file_handle << std::scientific << N[i][n_size-1];;
        }
        file_handle << std::endl;

        // Facet vertices
        file_handle << "    outer loop" << std::endl;
        std::size_t t_size = T[i].size();
        for (std::size_t j = 0; j < t_size; ++j) {
            std::size_t v_size = V[T[i][j]].size();
            file_handle << "      vertex ";
            if (v_size > 0) {
                for (std::size_t k = 0; k < v_size-1; ++k) {
                    file_handle << std::scientific << V[T[i][j]][k] << " ";
                }
                file_handle << std::scientific << V[T[i][j]][v_size-1];
            }
            file_handle << std::endl;
        }
        file_handle << "    endloop"    << std::endl;

        // Facet footer
        file_handle << "  " << ASCII_FACET_END << std::endl;
    }

    // Solid footer
    std::stringstream sfooter;
    sfooter << ASCII_SOLID_END << " " << solid_name;

    std::string footer;
    footer = sfooter.str();
    footer = utils::string::trim(footer);
    file_handle << footer << std::endl;

    sfooter.str("");

    // Restore stream flags
    file_handle.flags(streamFlags);

    return 0;
}

/*!
    Write solid data to ASCII STL file. Overloading of writeSolidASCII() for
    container std::vector<array<double, 3>>.

    \param[in,out] file_handle stream to STL file
    \param[in,out] nV number of vertices to be written to STL file.
    \param[in,out] nT number of facet to be written to STL file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int STLObj::writeSolidASCII(std::ofstream &file_handle,
                                     std::size_t &nV, std::size_t &nT,  std::vector<std::array<double, 3>> &V,
                                     std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T,
                                     const std::string &solid_name)
{
    // Check input variable
    if (V.size() < static_cast<std::size_t>(nV)) {
        return 2;
    }

    if (T.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    if (N.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Save stream flags
    std::ios::fmtflags streamFlags(file_handle.flags());

    // Write header
    std::stringstream sheader;
    sheader << ASCII_SOLID_BEGIN << " " << solid_name;

    std::string header = sheader.str();
    header = utils::string::trim(header);
    file_handle << header << std::endl;

    sheader.str("");

    // Write solid
    for (std::size_t i = 0; i < nT; ++i) {
        // Facet header
        file_handle << "  " << ASCII_FACET_BEGIN;

        // Facet normal
        std::size_t n_size = N[i].size();
        if (n_size > 0) {
            file_handle << " normal ";
            for (std::size_t j = 0; j < n_size-1; ++j) {
                file_handle << std::scientific << N[i][j] << " ";
            }
            file_handle << std::scientific << N[i][n_size-1];;
        }
        file_handle << std::endl;

        // Facet vertices
        file_handle << "    outer loop" << std::endl;
        std::size_t t_size = T[i].size();
        for (std::size_t j = 0; j < t_size; ++j) {
            std::size_t v_size = V[T[i][j]].size();
            file_handle << "      vertex ";
            if (v_size > 0) {
                for (std::size_t k = 0; k < v_size-1; ++k) {
                    file_handle << std::scientific << V[T[i][j]][k] << " ";
                }
                file_handle << std::scientific << V[T[i][j]][v_size-1];
            }
            file_handle << std::endl;
        }
        file_handle << "    endloop"    << std::endl;

        // Facet footer
        file_handle << "  " << ASCII_FACET_END << std::endl;
    }

    // Solid footer
    std::stringstream sfooter;
    sfooter << ASCII_SOLID_END << " " << solid_name;

    std::string footer;
    footer = sfooter.str();
    footer = utils::string::trim(footer);
    file_handle << footer << std::endl;

    sfooter.str("");

    // Restore stream flags
    file_handle.flags(streamFlags);

    return 0;
}

/*!
    Write solid data to binary STL file.

    \param[in,out] file_handle stream to STL file
    \param[in,out] nV number of vertices to be written to STL file.
    \param[in,out] nT number of facet to be written to STL file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int STLObj::writeSolidBINARY(std::ofstream &file_handle,
                                      std::size_t &nV, std::size_t &nT,  std::vector<std::vector<double>> &V,
                                      std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                                      const std::string &solid_name)
{
    BITPIT_UNUSED(solid_name);

    // Check input variable
    if (V.size() < static_cast<std::size_t>(nV)) {
        return 2;
    }

    if (T.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    if (N.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Define data types
    unsigned int  int2byte;
    unsigned int *int2byte_ptr = &int2byte;

    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Write header
    float4byte = (float)0.0;
    for (int i = 0; i < 20; ++i) {
        file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
    }

    // Write number of facets
    longint4byte = (unsigned long int) nT;
    file_handle.write(reinterpret_cast<char*>(longint4byte_ptr), 4);

    // Write facets
    for (std::size_t i = 0; i < nT; ++i) {
        // Normals
        std::size_t n_size = N[i].size();
        for (std::size_t j = 0; j < n_size; ++j) {
            float4byte = (float) N[i][j];
            file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
        }

        // Vertex
        std::size_t t_size = T[i].size();
        for (std::size_t j = 0; j < t_size; ++j) {
            std::size_t v_size = V[T[i][j]].size();
            for (std::size_t k = 0; k < v_size; ++k) {
                float4byte = (float) V[T[i][j]][k];
                file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
            }
        }

        // Attribute byte count
        //
        // In the standard format, this should be zero because most software
        // does not understand anything else.
        int2byte = 0;
        file_handle.write(reinterpret_cast<char*>(int2byte_ptr), 2);
    }

    return 0;
}

/*!
    Write solid data to binary STL file. Overloading of writeSolidBINARY()
    for container std::vector<array<double, 3>>

    \param[in,out] file_handle stream to STL file
    \param[in,out] nV number of vertices to be written to STL file.
    \param[in,out] nT number of facet to be written to STL file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int STLObj::writeSolidBINARY(std::ofstream &file_handle,
                                      std::size_t &nV, std::size_t &nT, std::vector<std::array<double, 3>> &V,
                                      std::vector<std::array<double, 3>> &N, std::vector<std::array<std::size_t, 3>> &T,
                                      const std::string &solid_name)
{
    BITPIT_UNUSED(solid_name);

    // Check input variable
    if (V.size() < static_cast<std::size_t>(nV)) {
        return 2;
    }

    if (T.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    if (N.size() < static_cast<std::size_t>(nT)) {
        return 2;
    }

    // Check stream status
    if (!file_handle.good()) {
        return 1;
    }

    // Define data types
    unsigned int  int2byte;
    unsigned int *int2byte_ptr = &int2byte;

    unsigned long int  longint4byte;
    unsigned long int *longint4byte_ptr = &longint4byte;

    float  float4byte;
    float *float4byte_ptr = &float4byte;

    // Write header
    float4byte = (float)0.0;
    for (int i = 0; i < 20; ++i) {
        file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
    }

    // Write number of facets
    longint4byte = (unsigned long int) nT;
    file_handle.write(reinterpret_cast<char*>(longint4byte_ptr), 4);

    // Write facets
    for (std::size_t i = 0; i < nT; ++i) {
        // Normals
        std::size_t n_size = N[i].size();
        for (std::size_t j = 0; j < n_size; ++j) {
            float4byte = (float) N[i][j];
            file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
        }

        // Vertex
        std::size_t t_size = T[i].size();
        for (std::size_t j = 0; j < t_size; ++j) {
            std::size_t v_size = V[T[i][j]].size();
            for (std::size_t k = 0; k < v_size; ++k) {
                float4byte = (float) V[T[i][j]][k];
                file_handle.write(reinterpret_cast<char*>(float4byte_ptr), 4);
            }
        }

        // Attribute byte count
        //
        // In the standard format, this should be zero because most software
        // does not understand anything else.
        int2byte = 0;
        file_handle.write(reinterpret_cast<char*>(int2byte_ptr), 2);
    }

    return 0;
}

}
