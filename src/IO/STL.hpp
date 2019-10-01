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

#ifndef __BITPIT_STL_HPP__
#define __BITPIT_STL_HPP__

#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "Operators.hpp"

namespace bitpit {

/*!
    @struct STLData
    @ingroup STereoLithography
    @brief structure holding meta-information
*/
struct STLData {
    int n_solids;                             /**< number of stl solids */
    std::vector<std::string> solid_names;     /**< solids names */
    std::vector<int> solid_facets;            /**< number of facet for each stl solid */
};

class STLObj {

public:
    enum FileFormat {
        FormatInvalid = -1,
        FormatASCII,
        FormatBinary
    };

    std::string stl_name;                         /** stl file name */
    bool stl_type;                                /** flag for binary/ASCII stl file */
    unsigned int err;                             /** general error */

    std::vector<std::vector<bool>> stl_errors;    /** error flags for each stl solid */
    STLData data;                                 /** stl data */

    STLObj();
    STLObj(std::string filename, bool filetype);
    STLObj(std::string filename);

    FileFormat detectFileFormat(const std::string &filename);

    void open(const std::string &mode);
    void close(const std::string &mode = "");
    void clear(void);

    void display(std::ostream &out);

    void scan();
    void check();

    void load(int &nV, int &nT, std::vector<std::vector<double>> &V,
              std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T);

    void load(int &nV, int &nT, std::vector<std::array<double, 3>> &V,
              std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T);

    void loadSolid(int &nV, int &nT, std::vector<std::vector<double>> &V,
                   std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T,
                   std::string &name);

    void loadSolid(int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                   std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T,
                   std::string &name);

    template <typename ... T2>
    void load(std::string name, int &nV, int &nT, std::vector<std::vector<double>> &V,
              std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T,
              T2 & ... others);

    template <typename ... T2>
    void load(std::string name, int &nV, int &nT, std::vector<std::array<double,3>> &V,
              std::vector<std::array<double,3>> &N, std::vector<std::array<int,3>> &T,
              T2 & ... others);

    void saveSolid(const std::string &name, int &nV, int &nT, std::vector<std::vector<double>> &V,
                   std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T);

    void saveSolid(const std::string &name, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                   std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T);

    template <typename ... T2>
    void save(const std::string &name, int &nV, int &nT, std::vector<std::vector<double>> &V,
              std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T,
              T2 & ... others);

    template <typename ... T2>
    void save(const std::string &name, int &nV, int &nT, std::vector<std::array<double,3>> &V,
              std::vector<std::array<double,3>> &N, std::vector<std::array<int,3>> &T,
              T2 & ... others);

    template <typename ... T2>
    void append(const std::string &name, int &nV, int &nT, std::vector<std::vector<double>> &V,
                std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T,
                T2 & ... others);

    template <typename ... T2>
    void append(const std::string &name, int &nV, int &nT, std::vector<std::array<double,3>> &V,
                std::vector<std::array<double,3>> &N, std::vector<std::array<int,3>> &T,
                T2 & ... others);

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

private:
    std::ifstream m_ifile_handle;      /**< input stream to stl file */
    std::ofstream m_ofile_handle;      /**< output stream to stl file */

    void save();
    void load();

    unsigned int scanASCII(std::ifstream &file_handle, std::vector<std::string> &solid_names,
                           std::vector<int> &solid_facets);

    unsigned int scanBINARY(std::ifstream &file_handle, std::vector<std::string> &solid_names,
                            std::vector<int> &solid_facets);

    unsigned int scanSolidASCII(std::ifstream &file_handle, int &nT);

    unsigned int checkASCII(std::ifstream &file_handle, std::vector<std::vector<bool>> &err_map);

    unsigned int checkSolidASCII(std::ifstream &file_handle, std::vector<bool> &err_map);

    unsigned int checkFacetASCII(std::ifstream &file_handle, std::vector<bool> &err_map);

    unsigned int checkBINARY(std::ifstream &file_handle, std::vector<std::vector<bool>> &err_map);

    unsigned int readSolidASCII(std::ifstream &file_handle, bool wrapAround, int &nV, int &nT,
                                std::vector<std::vector<double>> &V, std::vector<std::vector<double>> &N,
                                std::vector<std::vector<int>> &T, std::string &name);

    unsigned int readSolidASCII(std::ifstream &file_handle, bool wrapAround, int &nV, int &nT,
                                std::vector<std::array<double, 3>> &V, std::vector<std::array<double, 3>> &N,
                                std::vector<std::array<int, 3>> &T, std::string &name);

    unsigned int readASCII(std::ifstream &file_handle, int &nV, int &nT, std::vector<std::vector<double>> &V,
                           std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T);

    unsigned int readASCII(std::ifstream &file_handle, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                           std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T);

    unsigned int readFacetASCII(std::ifstream &file_handle,  int &nV, int &nT, std::vector<std::vector<double>> &V,
                                std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T);

    unsigned int readFacetASCII(std::ifstream &file_handle, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                                std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T);

    unsigned int readBINARY(std::ifstream &file_handle, int &nV, int &nT, std::vector<std::vector<double>> &V,
                            std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T);

    unsigned int readBINARY(std::ifstream &file_handle, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                            std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T);

    unsigned int writeSolidASCII(std::ofstream &file_handle, int &nV, int &nT, std::vector<std::vector<double>> &V,
                                 std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T, const std::string &solid_name = "");

    unsigned int writeSolidASCII(std::ofstream &file_handle, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                                 std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T, const std::string &solid_name = "");

    unsigned int writeSolidBINARY(std::ofstream &file_handle, int &nV, int &nT, std::vector<std::vector<double>> &V,
                                  std::vector<std::vector<double>> &N, std::vector<std::vector<int>> &T, const std::string &solid_name = "");

    unsigned int writeSolidBINARY(std::ofstream &file_handle, int &nV, int &nT, std::vector<std::array<double, 3>> &V,
                                  std::vector<std::array<double, 3>> &N, std::vector<std::array<int, 3>> &T, const std::string &solid_name = "");

};

// Template implementation
#include "STL.tpp"

}

#endif
