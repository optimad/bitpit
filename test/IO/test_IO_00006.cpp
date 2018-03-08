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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing auto detection of STL file type.
*/
int subtest_001()
{
    std::cout << "Reading single-solid ASCII STL" << std::endl;

    STLObj reader;
    STLObj::FileFormat format = reader.detectFileFormat("./data/cubeAscii.stl");

    if (format == STLObj::FileFormat::FormatASCII) {
        std::cout << " Test passed (ASCII format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatBinary) {
        std::cout << " Test failed (binary format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatInvalid) {
        std::cout << " Test failed (invalid file detected)" << std::endl;
    }

    return (int) (format != STLObj::FileFormat::FormatASCII);
}

/*!
* Subtest 002
*
* Testing auto detection of STL file type.
*/
int subtest_002()
{
    std::cout << "Reading single-solid binary STL" << std::endl;

    STLObj reader;
    STLObj::FileFormat format = reader.detectFileFormat("./data/cubeBinary.stl");

    if (format == STLObj::FileFormat::FormatASCII) {
        std::cout << " Test failed (ASCII format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatBinary) {
        std::cout << " Test passed (binary format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatInvalid) {
        std::cout << " Test failed (invalid file detected)" << std::endl;
    }

    return (int) (format != STLObj::FileFormat::FormatBinary);
}

/*!
* Subtest 003
*
* Testing auto detection of STL file type.
*/
int subtest_003()
{
    std::cout << "Reading multi-solid binary STL" << std::endl;

    STLObj reader;
    STLObj::FileFormat format = reader.detectFileFormat("./data/cubeAscii_MultiSolid.stl");

    if (format == STLObj::FileFormat::FormatASCII) {
        std::cout << " Test passed (ASCII format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatBinary) {
        std::cout << " Test failed (binary format detected)" << std::endl;
    } else if (format == STLObj::FileFormat::FormatInvalid) {
        std::cout << " Test failed (invalid file detected)" << std::endl;
    }

    return (int) (format != STLObj::FileFormat::FormatASCII);
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

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Testing STL binary/ascii detection" << std::endl;

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

        status = subtest_003();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
