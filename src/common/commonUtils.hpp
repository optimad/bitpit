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

#ifndef __BITPIT_COMMON_UTILS_HPP__
#define __BITPIT_COMMON_UTILS_HPP__

/*! \file */

#include <array>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

// Stringification macro
#define BITPIT_STR2(X) #X
#define BITPIT_STR(X) BITPIT_STR2(X)

// Macros to allow using oveload in preprocessing macro
#define BITPIT_CAT(A, B) A ## B

#define BITPIT_COUNT_ARGS_MAX6(_1, _2, _3, _4, _5, _6 /* ad nauseam */, COUNT, ...) COUNT
#define BITPIT_EXPAND_ARGS_FOR_COUNT(ARGS) BITPIT_COUNT_ARGS_MAX6 ARGS
#define BITPIT_ARGS_SIZE(...) BITPIT_EXPAND_ARGS_FOR_COUNT((__VA_ARGS__, 6, 5, 4, 3, 2, 1, 0))

#define BITPIT_SELECT_OVERLOAD(NAME, NUM) BITPIT_CAT(NAME ## _, NUM)
#define BITPIT_OVERLOAD_CALL(NAME, ...) BITPIT_SELECT_OVERLOAD(NAME, BITPIT_ARGS_SIZE(__VA_ARGS__))(__VA_ARGS__)

namespace bitpit {

namespace utils {

template <typename T, typename Comparator = std::less<T> >
bool addToOrderedVector(const T &value, std::vector<T> &list, Comparator comparator = Comparator());

template <typename T, typename Comparator = std::less<T> >
typename std::vector<T>::const_iterator findInOrderedVector(const T &value, const std::vector<T> &list, Comparator comparator = Comparator());

template<typename T>
void reorderVector(std::vector<size_t>& order, std::vector<T>& v, const size_t &size);

template <class T>
void eraseValue(std::vector<T> &, const T&);

template <class T>
std::vector<T> intersectionVector(const std::vector<T>&, const std::vector<T>&);

#ifndef __BITPIT_UTILS_SRC__
extern template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
extern template bool addToOrderedVector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

extern template std::vector<long>::const_iterator findInOrderedVector<>(const long&, const std::vector<long>&, std::less<long>);
extern template std::vector<unsigned long>::const_iterator findInOrderedVector<>(const unsigned long&, const std::vector<unsigned long>&, std::less<unsigned long>);
#endif

void extractWithoutReplacement(                                               // Extract integers without replacement
    int                         ,                                             // (input) number of integers to be extracted
    int                         ,                                             // (input) upper bound of extraction interval
    std::vector<int>           &                                              // (input/output) list of extracted value
);

/*!
    Functor to compare two double precision floating point numbers.

    See: https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
*/
struct DoubleFloatingEqual
{
    /*!
        Compares the specified double precision floating point numbers
        and returns true if the numbers match.

        \param x if the first value to compare
        \param y if the second value to compare
        \result Returns true if the numbers match, false otherwise.
    */
    bool operator()(const double &x, const double &y) const
    {
        const double ABS_MAX_DIFF = 1e-14;
        const double REL_MAX_DIFF = DBL_EPSILON;

        // Check if the numbers are really close (needed when comparing
        // numbers near zero).
        double diff = std::abs(x - y);
        if (diff <= ABS_MAX_DIFF) {
            return true;
        }

        // Check if the numbers have the same sign
        if ((x < 0 && y > 0) || (x > 0 && y < 0)) {
            return false;
        }

        // Compare using a relative difference
        double abs_x   = std::abs(x);
        double abs_y   = std::abs(y);
        double largest = (abs_y > abs_x) ? abs_y : abs_x;

        return (diff <= largest * REL_MAX_DIFF);
    }
};

}

}

// Template implementation
#include "commonUtils.tpp"

#endif
