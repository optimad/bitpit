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

#ifndef __BITPIT_TEST_HELPERS_UNIT_TEST_PARALLEL__
#define __BITPIT_TEST_HELPERS_UNIT_TEST_PARALLEL__

/*!
 * \file unit_test_parallel.hpp
 *
 * \brief Infrastruture needed for running parallel unit tests.
 *
 * This file should be included at the beginning of a parallel unit test suite.
 * It provides all the infrastructure needed to run the tests contained inside
 * the suite.
 *
 * Before including this file it is necessary to define the name of the
 * test suite using the macro #BITPIT_UNIT_TEST_SUITE_NAME.
 */

#include <mpi.h>

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
    // Initialize MPI
    assert(BITPIT_ENABLE_MPI);
    MPI_Init(&argc,&argv);

    // Initialize the logger
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bitpit::log::manager().initialize(bitpit::log::COMBINED, true, nProcs, rank);
    bitpit::log::cout().setVisibility(bitpit::log::GLOBAL);

    // Run the tests
#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
    extern bool init_unit_test();

    boost::unit_test::init_unit_test_func init_func = &init_unit_test;
#else
    extern boost::unit_test::test_suite * init_unit_test_suite(int argc, char* argv[]);

    boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
#endif

    int status = boost::unit_test::unit_test_main(init_func, argc, argv);

    // Finalize MPI
    MPI_Finalize();

    return status;
}

#endif
