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

/**
 * \ingroup common_macro
 *
 * Create a workspace with the specified size.
 *
 * If the size is less than the specified stack size, the workspace will
 * make use of an array created on the stack, otherwise a container with
 * the specified size will be created on the heap and that container will
 * be used for the workspace.
 *
 * NOTE: this macro will always create an array on the stack, whether it
 * will be used as workspace or not depend on the requested workspace size.
 *
 * \param workspace is the name of the workspace
 * \param item_type is the type of items the workspace will contain
 * \param size is the size of the workspace
 * \param stack_size is the maximum size the workspace can have to be
 * allocated on the stack
 */
#define BITPIT_CREATE_WORKSPACE(workspace, item_type, size, stack_size) \
item_type *workspace;                                            \
std::array<item_type, stack_size> workspace##_stack;             \
std::unique_ptr<std::vector<item_type>> workspace##_heap;        \
if (size <= stack_size) {                                        \
    workspace = workspace##_stack.data();                        \
} else {                                                         \
    workspace##_heap = std::unique_ptr<std::vector<item_type>>(new std::vector<item_type>(size)); \
    workspace = workspace##_heap->data();                        \
}

namespace bitpit {

namespace utils {

template <typename T, typename Comparator = std::less<T> >
bool addToOrderedVector(const T &value, std::vector<T> &list, Comparator comparator = Comparator());

template <typename T, typename Comparator = std::less<T> >
typename std::vector<T>::const_iterator findInOrderedVector(const T &value, const std::vector<T> &list, Comparator comparator = Comparator());

template<typename T>
void reorderVector(std::vector<size_t>& order, std::vector<T>& v, std::size_t size);

template <class T>
void eraseValue(std::vector<T> &, const T&);

template <class T>
std::vector<T> intersectionVector(const std::vector<T>&, const std::vector<T>&);

#ifndef __BITPIT_COMMON_UTILS_SRC__
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

unsigned long factorial(unsigned long n);

/*!
    Functor to compare two double precision floating point numbers.

    See chapter 4.2.2 equation 35 of "The art of computer programming (vol II)"
    Donald. E. Knuth, 1998, Addison-Wesley Longman, Inc., ISBN 0-201-89684-2,
    Addison-Wesley Professional; 3rd edition.

    Also see: "https://www.boost.org/doc/libs/1_74_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html"
*/
struct DoubleFloatingEqual
{
    /*!
        Compares the specified double precision floating point numbers
        and returns true if the numbers match.

        \param x if the first value to compare
        \param y if the second value to compare
        \param relativeTolerance is the relative tolerance that will be used to
        perform the comparison
        \param absoluteTolerance is the absolute tolerance that will be used to
        perform the comparison
        \result Returns true if the numbers match, false otherwise.
    */
    bool operator()(double x, double y, double relativeTolerance = DEFAULT_REL_TOL, double absoluteTolerance = DEFAULT_ABS_TOL) const
    {
        // Check if the numbers are really close (needed when comparing
        // numbers near zero).
        double absoluteDifference = std::abs(x - y);
        if (absoluteDifference <= absoluteTolerance) {
            return true;
        }

        // Compare using a relative difference
        //
        // The check that is performed is:
        //
        //  abs(x - y) <= abs(x) * epsilon OR abs(x - y) <= abs(y) * epsilon
        //
        // The multiplication in the right side of inequalities may cause an
        // unwanted underflow condition. To prevent this, the difference is
        // scaled rather than epsilon.
        double relativeDifference = absoluteDifference / std::max(std::abs(x), std::abs(y));

        return (relativeDifference <= relativeTolerance);
    }

private:
    constexpr static const double DEFAULT_ABS_TOL = 10 * std::numeric_limits<double>::epsilon();
    constexpr static const double DEFAULT_REL_TOL = 10 * std::numeric_limits<double>::epsilon();

};

}

}

// Template implementation
#include "commonUtils.tpp"

#endif
