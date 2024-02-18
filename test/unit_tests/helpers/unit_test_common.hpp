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

#ifndef __BITPIT_TEST_HELPERS_UNIT_TEST_COMMON__
#define __BITPIT_TEST_HELPERS_UNIT_TEST_COMMON__

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE BITPIT_UNIT_TEST_SUITE_NAME
#include <boost/test/included/unit_test.hpp>

#ifdef _MSC_VER
#   ifdef min
#       undef min
#   endif
#   ifdef max
#       undef max
#   endif
#   ifdef interface
#       undef interface
#   endif
#endif

/*!
 * Display the name of the current unit test.
 *
 * \param stream is the sream that will be used to the output
 */
#define BITPIT_UNIT_TEST_DISPLAY_NAME(stream) stream << "Running test \"" << boost::unit_test::framework::current_test_case().p_name << "\"..." << std::endl;

/*!
 * \def BITPIT_UNIT_TEST_SUITE_NAME
 *
 * Defines the name of a unit test suite. It should be defined before any
 * inclusion directive to unit test helpers.
 */

#endif
