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

#ifndef __BITPIT_TEST_HELPERS_UNIT_TEST__
#define __BITPIT_TEST_HELPERS_UNIT_TEST__

/*!
 * \file unit_test.hpp
 *
 * \brief Infrastruture needed for running parallel unit tests.
 *
 * This file should be included at the beginning of a unit test suite. It
 * provides all the infrastructure needed to run the tests contained inside
 * the suite.
 *
 * Before including this file it is necessary to define the name of the
 * test suite using the macro BITPIT_UNIT_TEST_SUITE_NAME.
 */

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "helpers/unit_test_common.hpp"

#include "bitpit_IO.hpp"

/*!
 * Main program.
 *
 * \param argc is the argument count of the command line arguments
 * \param argv is the argument vector of the command line arguments
 * \return Returns 0 upon exit success.
 */
int main(int argc, char* argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    bitpit::log::manager().initialize(bitpit::log::COMBINED);

    // Run the tests
#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
    extern bool init_unit_test();

    boost::unit_test::init_unit_test_func init_func = &init_unit_test;
#else
    extern boost::unit_test::test_suite * init_unit_test_suite(int argc, char* argv[]);

    boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
#endif

    int status = boost::unit_test::unit_test_main(init_func, argc, argv);

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}

#endif
