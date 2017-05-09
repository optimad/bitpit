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

#ifndef __BITPIT_COMPILER_HPP__
#define __BITPIT_COMPILER_HPP__

/*! \file */

/*!
 * @ingroup macro
 * @{
 */

/**
 * Unreachable macro.
 *
 * Useful for suppressing "control reaches end of non-void function" warnings.
 *
 * @param str is the error message that will be displayed it the unreachable
 * code is reached
 */
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 5
#define BITPIT_UNREACHABLE(str)    \
do {                         \
    assert(!str);            \
    __builtin_unreachable(); \
} while (0)
#elif (defined(__clang__) && defined(__has_builtin))
# if __has_builtin(__builtin_unreachable)
#  define BITPIT_UNREACHABLE(str)  \
do {                         \
    assert(!str);            \
    __builtin_unreachable(); \
} while (0)
# endif
#endif

#ifndef BITPIT_UNREACHABLE
#define BITPIT_UNREACHABLE(str)
#endif

/*!
 * Unused macro.
 *
 * @param variable is the name of variable to be marked as unused
 */
#define BITPIT_UNUSED(variable)     \
do {                  \
    (void)(variable); \
} while (0)


/*!
 * Deprecated macro.
 *
 * @param func id the function/method to be marked as deprecated
 */
#if defined __GNUC__
#   define BITPIT_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined __clang__
#   define BITPIT_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
    #define BITPIT_DEPRECATED(func) __declspec(deprecated) func
#else
#   pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#   define BITPIT_DEPRECATED(func) func
#endif

/*!
 * @}
 */

#endif
