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
#ifndef __BITPIT_COMMON_BINARY_UTILS_HPP__
#define __BITPIT_COMMON_BINARY_UTILS_HPP__

#include <array>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

#include "commonUtils.hpp"
#include "commonUtils.tpp"

namespace bitpit {

namespace utils {

/*!
    \ingroup utils::binary
    \brief The namespace 'binary' contains routines for handling binary
    archives.
*/
namespace binary {

template<typename T, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void write(std::ostream &stream, const std::vector<T> &value);

template<>
void write(std::ostream &stream, const std::vector<bool> &value);

template<typename T, std::size_t dim, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void write(std::ostream &stream, const std::array<T, dim> &value);

template<typename T, typename std::enable_if<utils::is_iterable<T>::value>::type* = nullptr>
void write(std::ostream &stream, const T &value);

template<typename T, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void write(std::ostream &stream, const T &value);

template<typename T>
void write(std::ostream &stream, const T &value, size_t size);

template<typename T>
void write(std::ostream &stream, const T *value, size_t size);

void write(std::ostream &stream, const std::string &string);

template<typename T, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void read(std::istream &stream, std::vector<T> &value);

template<>
void read(std::istream &stream, std::vector<bool> &value);

template<typename T, std::size_t dim, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void read(std::istream &stream, std::array<T, dim> &value);

template<typename T, typename std::enable_if<utils::is_iterable<T>::value>::type* = nullptr>
void read(std::istream &stream, T &value);

template<typename T, typename std::enable_if<std::is_pod<T>::value>::type* = nullptr>
void read(std::istream &stream, T &value);

template<typename T>
void read(std::istream &stream, T &value, size_t size);

template<typename T>
void read(std::istream &stream, T *value, size_t size);

void read(std::istream &stream, std::string &string);

}

}

}

// Include templates' implementation
#include "binaryUtils.tpp"

#endif
