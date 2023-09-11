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

#include "bitpit_containers.hpp"

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include <stdexcept>

using namespace bitpit;


/*!
* Subtest 001
*
* Testing move constructor capabilities.
*/
int subtest_001()
{
    std::cout << std::endl;
    std::cout << "Testing move constructor capabilities" << std::endl;

    // Create container
    std::cout << "Creating container..." << std::endl;

    PiercedVector<double> containerExpected;
    containerExpected.insert(21, 12.0);
    containerExpected.insert(22, 13.0);
    containerExpected.insert(23, 14.0);

    std::cout << "  Size of container ....... " << containerExpected.size() << std::endl;

    // Create a new container moving data from container with move constructor
    std::cout << "Moving container..." << std::endl;

    PiercedVector<double> container(containerExpected);
    container.checkIntegrity();
    std::cout << "  Size of original container ....... " << container.size() << std::endl;

    PiercedVector<double> movedContainer(std::move(container));
    movedContainer.checkIntegrity();
    std::cout << "  Size of moved container .......... " << movedContainer.size() << std::endl;

    if (containerExpected.size() != movedContainer.size()) {
        throw std::runtime_error("Contents of moved container doesn't match expected values");
    }

    for (auto itr = movedContainer.begin(), end = movedContainer.end(); itr != end; ++itr) {
        long id = itr.getId();
        if (containerExpected[id] != *itr) {
            throw std::runtime_error("Contents of moved container doesn't match expected values");
        }
    }

    std::cout << "Test completed." << std::endl;

    return 0;
}

/*!
* Subtest 002
*
* Testing move assignment operator capabilities.
*/
int subtest_002()
{
    std::cout << std::endl;
    std::cout << "Testing move assignment operator capabilities" << std::endl;

    // Create source containers
    std::cout << "Creating container..." << std::endl;

    PiercedVector<double> containerExpected;
    containerExpected.insert(21, 12.0);
    containerExpected.insert(22, 13.0);
    containerExpected.insert(23, 14.0);

    std::cout << "  Size of container ....... " << containerExpected.size() << std::endl;

    // Create a new container moving data from container with move assignment
    std::cout << "Moving container..." << std::endl;

    PiercedVector<double> container(containerExpected);
    container.checkIntegrity();

    PiercedVector<double> movedContainer;
    movedContainer = std::move(container);
    movedContainer.checkIntegrity();

    std::cout << "  Size of original container ....... " << container.size() << std::endl;
    std::cout << "  Size of moved container .......... " << movedContainer.size()<< std::endl;

    if (containerExpected.size() != movedContainer.size()) {
        throw std::runtime_error("Contents of moved container doesn't match expected values");
    }

    for (auto itr = movedContainer.begin(), end = movedContainer.end(); itr != end; ++itr) {
        long id = itr.getId();
        if (containerExpected[id] != *itr) {
            throw std::runtime_error("Contents of moved container doesn't match expected values");
        }
    }

    std::cout << "Test completed." << std::endl;

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Run the subtests
    std::cout << "Testing PiercedVector move constructor/operator" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }

        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        std::cout << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
