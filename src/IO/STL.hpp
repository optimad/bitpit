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

#ifndef __BITPIT_STL_HPP__
#define __BITPIT_STL_HPP__

#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "line_stream.hpp"

namespace bitpit {

class STLBase {

public:
    enum Format {
        FormatUnknown = -1,
        FormatASCII,
        FormatBinary
    };

    virtual ~STLBase() = default;

    const std::string & getFilename() const;
    Format getFormat() const;

protected:
    typedef uint8_t  BINARY_UINT8;
    typedef uint16_t BINARY_UINT16;
    typedef uint32_t BINARY_UINT32;
    typedef float    BINARY_REAL32;

    static const std::size_t BINARY_HEADER_SIZE;
    static const std::size_t BINARY_MINIMUM_SIZE;

    static const std::string ASCII_SOLID_BEGIN;
    static const std::string ASCII_SOLID_END;
    static const std::string ASCII_FACET_BEGIN;
    static const std::string ASCII_FACET_END;
    static const std::string ASCII_FILE_BEGIN;
    static const std::string ASCII_FILE_END;
    static const std::size_t ASCII_MINIMUM_SIZE;

    STLBase(const std::string &filename);
    STLBase(const std::string &filename, Format format);

    void setFilename(const std::string &filename);
    void setFormat(Format format);

private:
    std::string m_filename;
    Format m_format;

};

class STLReader : public STLBase {

public:
    /*!
        @struct InspectionInfo
        @ingroup STereoLithography
        @brief Structure holding inspection information
    */
    struct InspectionInfo {
        int nSolids;                                      /**< Number of solids */

        std::vector<bool> solidValid;                     /**< Errors associated with the solids */
        std::vector<std::array<bool, 6>> solidErrors;     /**< Errors associated with the solids */

        std::vector<std::string> solidNames;              /**< Name of the solids */
        std::vector<std::size_t> solidFacetCount;         /**< Number of facets of each solid */
        std::vector<std::size_t> solidVertexCount;        /**< Number of vertices of each solid */
    };

    static Format detectFormat(const std::string &filename);

    STLReader(const std::string &filename, Format format = FormatUnknown);

    int inspect(InspectionInfo *info);
    void displayInspectionInfo(const InspectionInfo &info, std::ostream &out) const;

    int readBegin();
    int readEnd();

    int readSolid(std::string *name, std::size_t *nV, std::size_t *nT,
                  std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                  std::vector<std::array<std::size_t, 3>> *T);

    int readSolid(const std::string &solid, std::string *name, std::size_t *nV, std::size_t *nT,
                  std::vector<std::array<double, 3>> *V, std::vector<std::array<double, 3>> *N,
                  std::vector<std::array<std::size_t, 3>> *T);

    int readHeader(std::string *name, std::size_t *nT);
    int readHeader(const std::string &solid, std::string *name, std::size_t *nT);

    int readFooter();
    int readFooter(const std::string &solid);

    int readFacet(std::array<double, 3> *V0, std::array<double, 3> *V1,
                  std::array<double, 3> *V2, std::array<double, 3> *N);

private:
    std::ifstream m_fileHandle;      /**< File handle */

    LineStream m_lineStream;    /**< Line stream */

    int inspectASCII(InspectionInfo *info);
    int inspectSolidASCII(std::size_t *nFactes, std::array<bool, 6> *errors);
    int inspectFacetASCII(std::array<bool, 6> *errors);

    int inspectBinary(InspectionInfo *info);

    int readHeaderASCII(const std::string &solid, std::string *name, std::size_t *nT);

    int readFooterASCII(const std::string &solid);

    int readFacetASCII(std::array<double, 3> *V0, std::array<double, 3> *V1,
                       std::array<double, 3> *V2, std::array<double, 3> *N);

    int readHeaderBinary(std::string *name, std::size_t *nT);

    int readFooterBinary();

    int readFacetBinary(std::array<double, 3> *V0, std::array<double, 3> *V1,
                        std::array<double, 3> *V2, std::array<double, 3> *N);

};

class STLWriter : public STLBase {

public:
    enum WriteMode {
        WriteOverwrite,
        WriteAppend
    };

    STLWriter(const std::string &filename, Format format);

    int writeBegin(WriteMode writeMode, bool partialWrite = false);
    int writeEnd();

    int writeSolid(const std::string &name, std::size_t nV, std::size_t nT,
                   const std::vector<std::vector<double>> &V, const std::vector<std::vector<double>> &N,
                   const std::vector<std::vector<std::size_t>> &T);

    int writeSolid(const std::string &name, std::size_t nV, std::size_t nT,
                   const std::vector<std::array<double, 3>> &V, const std::vector<std::array<double, 3>> &N,
                   const std::vector<std::array<std::size_t, 3>> &T);

    int writeHeader(const std::string &name, std::size_t nT);

    int writeFooter(const std::string &name);

    int writeFacet(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                   const std::array<double, 3> &V2, const std::array<double, 3> &N);

private:
    std::ofstream m_fileHandle;      /**< File handle */

    int writeHeaderASCII(const std::string &name, std::size_t nT);

    int writeFooterASCII(const std::string &name);

    int writeFacetASCII(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                        const std::array<double, 3> &V2, const std::array<double, 3> &N);

    int writeHeaderBinary(const std::string &name, std::size_t nT);

    int writeFooterBinary(const std::string &name);

    int writeFacetBinary(const std::array<double, 3> &V0, const std::array<double, 3> &V1,
                         const std::array<double, 3> &V2, const std::array<double, 3> &N);

};

}

#endif
