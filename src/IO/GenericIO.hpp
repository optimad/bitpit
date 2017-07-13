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

#ifndef __BITPIT_GENERIC_IO_HPP__
#define __BITPIT_GENERIC_IO_HPP__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"

namespace bitpit {

/*!
 * @ingroup GenericIO
 * @brief contains routines for reading/writing in ASCII/BINARY format
 */
namespace genericIO{

template<class data_T>
void  lineStream(std::fstream &str, data_T &data);

template<class data_T>
void  lineStream(std::fstream &str, std::vector<data_T> &data);

template<class data_T, size_t d>
void  lineStream(std::fstream &str, std::array<data_T,d> &data);

template<class data_T, size_t d>
void  lineStream(std::fstream &str, data_T *data, int nr);

template<class data_T>
void flushASCII(std::fstream &str, const data_T &data);

template<>
void flushASCII(std::fstream &str, const uint8_t &data);

template<class data_T>
void flushASCII(std::fstream &str, int elements_per_line, const std::vector<data_T> &data);

template<class data_T, size_t d>
void flushASCII(std::fstream &str, int elements_per_line, const std::array<data_T,d> &data);

template<class data_T>
void flushASCII(std::fstream &str, int elements_per_line, const data_T *data, int nr);

template<class data_T>
void flushASCII(std::fstream &str, int elements_per_line, const bitpit::PiercedVector<data_T> &data, bool writeIndex = false);

template<class data_T>
void flushBINARY(std::fstream &str, const data_T &data);

template<class data_T>
void flushBINARY(std::fstream &str, const std::vector<data_T> &data);

template<class data_T>
void flushBINARY(std::fstream &str, const std::vector<std::vector<data_T>> &data);

template<class data_T, size_t d>
void flushBINARY(std::fstream &str, const std::vector<std::array<data_T,d>> &data);

template<class data_T, size_t d>
void flushBINARY(std::fstream &str, const std::array<data_T,d> &data);

template<class data_T>
void flushBINARY(std::fstream &str, const data_T *data, int nr);

template<class data_T>
void flushBINARY(std::fstream &str, const bitpit::PiercedVector<data_T> &, bool writeIndex = false);

template<class data_T>
void absorbASCII(std::fstream &str, data_T &data);

template<class data_T>
void absorbASCII(std::fstream &str, std::vector<data_T> &data);

template<class data_T, size_t d>
void absorbASCII(std::fstream &str, std::array<data_T,d> &data);

template<class data_T>
void absorbASCII(std::fstream &str, data_T *data, int nr);

template<class data_T>
void absorbASCII(std::fstream &str, bitpit::PiercedVector<data_T> &data);

template<class data_T>
void absorbASCII(std::fstream &str, bitpit::PiercedVector<data_T> &data, long);

template<class data_T>
void absorbBINARY(std::fstream &str, data_T &data);

template<class data_T>
void absorbBINARY(std::fstream &str, std::vector<data_T> &data);

template<class data_T>
void absorbBINARY(std::fstream &str, std::vector<std::vector<data_T>> &data);

template<class data_T, size_t d>
void absorbBINARY(std::fstream &str, std::vector<std::array<data_T,d>> &data);

template<class data_T, size_t d>
void absorbBINARY(std::fstream &str, std::array<data_T,d> &data);

template<class data_T>
void absorbBINARY(std::fstream &str, data_T *data, int nr);

template<class data_T>
void absorbBINARY(std::fstream &str, bitpit::PiercedVector<data_T> &data);

template<class data_T>
void absorbBINARY(std::fstream &str, bitpit::PiercedVector<data_T> &data, long);

void copyUntilEOFInString(std::fstream &str, char *&buffer, int &length);

}

}

#include "GenericIO.tpp"

#endif
