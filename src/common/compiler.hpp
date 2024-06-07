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

#ifndef __BITPIT_COMPILER_HPP__
#define __BITPIT_COMPILER_HPP__

/*! \file */

/**
 * \ingroup common_macro
 *
 * Unreachable macro.
 *
 * Useful for suppressing "control reaches end of non-void function" warnings.
 *
 * \param str is the error message that will be displayed it the unreachable
 * code is reached
 */
#ifdef BITPIT_ENABLE_BUILTIN_UNREACHABLE
#define BITPIT_UNREACHABLE(str)    \
do {                        \
   assert(!str);            \
   __builtin_unreachable(); \
} while (0)
#elif defined (_MSC_VER)
#define BITPIT_UNREACHABLE(str)    \
do {                        \
   assert(!str);            \
   __assume(0);             \
} while (0)
#else
#define BITPIT_UNREACHABLE(str) assert(!str)
#endif

/*!
 * \ingroup common_macro
 *
 * Unused macro.
 *
 * \param variable is the name of variable to be marked as unused
 */
#define BITPIT_UNUSED(variable)     \
do {                  \
    (void)(variable); \
} while (0)


/*!
 * \ingroup common_macro
 *
 * Deprecated macro.
 *
 * \param func id the function/method to be marked as deprecated
 */
#if defined __GNUC__
#   define BITPIT_DEPRECATED(func) func __attribute__ ((deprecated))
#   define BITPIT_DEPRECATED_FOR(func, replacement) func __attribute__ ((deprecated(" Use " #replacement)))
#elif defined __clang__
#   define BITPIT_DEPRECATED(func) func __attribute__ ((deprecated))
#   define BITPIT_DEPRECATED_FOR(func, replacement) func __attribute__ ((deprecated(" Use " #replacement)))
#elif defined(_MSC_VER)
#   define BITPIT_DEPRECATED(func) __declspec(deprecated) func
#   define BITPIT_DEPRECATED_FOR(func, replacement) __declspec(deprecated(" Use " #replacement)) func
#else
#   pragma message("WARNING: macros to declare functions as deprecated are not implemented for this compiler")
#   define BITPIT_DEPRECATED(func) func
#   define BITPIT_DEPRECATED_FOR(func, replacement) func
#endif


/*!
 * \ingroup common_macro
 *
 * Comma.
 */
#define BITPIT_COMMA ,


/*!
 * \ingroup common_misc
 *
 * Custom implementation of std::void_t for C++11 backwards compatibility.
 */
#if __cplusplus < 201703L
template <class...>
struct make_void { using type = void; };

template <typename... T>
using void_t = typename make_void<T...>::type;

template <typename... T>
using bitpit_void_t = void_t<T...>;
#else
template <typename... T>
using bitpit_void_t = std::void_t<T...>;
#endif

#endif
