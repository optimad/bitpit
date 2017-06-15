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

#ifndef __BITPIT_STRING_UTILS_HPP__
#define __BITPIT_STRING_UTILS_HPP__

/*! \file */

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace bitpit {

namespace utils {

/*!
    \ingroup common::utils::string

    \brief Namespace for string utility functions
*/
namespace string {

// Trimming operators
inline std::string &ltrim(std::string &s);
inline std::string &rtrim(std::string &s);
inline std::string &trim(std::string &s);

// Padding operators
inline std::string lfill(const int &nchars, std::string &s, char c);
inline std::string rfill(const int &nchars, std::string &s, char c);
inline std::string zeroPadNumber(int nchars, int num);

// Keyword search
bool getAfterKeyword(std::string line, std::string key, char del, std::string &result);
inline bool keywordInString(std::string line, std::string key);

// Conversion
template <class T>
void convertString(std::string input, T &output);

template <class T>
void convertString(std::string input, std::vector<T> &output);

template <class T, size_t n>
void convertString(std::string input, std::array<T,n> &output);

}

}

}

// Template implementation
#include "stringUtils.tpp"

#endif
