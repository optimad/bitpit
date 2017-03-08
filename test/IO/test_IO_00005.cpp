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

int main(int argc, char *argv[]) {

    //
    // Testing OBinaryArchive
    //
    std::cout << std::endl;
    std::cout << "::: Testing OBinaryArchive... " << std::endl;
    std::cout << std::endl;

    // Open archive
    std::cout << "Creating archive... " << std::endl;

    std::stringstream headerStream;
    headerStream << "Test header" << std::endl;
    headerStream << "|01234567890123456789" << std::endl;
    headerStream << "|01234567890123456789" << std::endl;
    headerStream << "|01234567890123456789";
    std::string header(headerStream.str());

    OBinaryArchive binaryWriter("binary_archive_test.dat", 1, header);

    // Write test out_data
    std::cout << "Writing data to archive... " << std::endl;

    int out_data_int       = 1;
    long out_data_long     = 11;
    double out_data_double = 22.22;
    bool out_data_bool     = 1;

    binaryWriter.write(reinterpret_cast<const char*>(&out_data_int), sizeof(out_data_int));
    binaryWriter.write(reinterpret_cast<const char*>(&out_data_long), sizeof(out_data_long));
    binaryWriter.write(reinterpret_cast<const char*>(&out_data_double), sizeof(out_data_double));
    binaryWriter.write(reinterpret_cast<const char*>(&out_data_bool), sizeof(out_data_bool));

    // Close archive
    binaryWriter.close();

    //
    // Testing IBinaryArchive
    //
    std::cout << std::endl;
    std::cout << "::: Testing IBinaryArchive... " << std::endl;
    std::cout << std::endl;

    // Open archive
    std::cout << "Reading archive... " << std::endl;

    IBinaryArchive binaryReader("binary_archive_test.dat");

    // read test data
    std::cout << "Reading data from archive... " << std::endl;

    int in_data_int       = 1;
    long in_data_long     = 11;
    double in_data_double = 22.22;
    bool in_data_bool     = 1;

    binaryReader.read(reinterpret_cast<char*>(&in_data_int), sizeof(in_data_int));
    binaryReader.read(reinterpret_cast<char*>(&in_data_long), sizeof(in_data_long));
    binaryReader.read(reinterpret_cast<char*>(&in_data_double), sizeof(in_data_double));
    binaryReader.read(reinterpret_cast<char*>(&in_data_bool), sizeof(in_data_bool));

    std::cout << "Header: " << std::endl;
    std::cout << binaryReader.getHeader() << std::endl;

    std::cout << "Integer data: " << in_data_int << " (read) vs " << out_data_int << " (expected) " << std::endl;
    if (in_data_int != out_data_int) {
        exit(1);
    }

    std::cout << "Long integer data: " << in_data_long << " (read) vs " << out_data_long << " (expected) " << std::endl;
    if (in_data_long != out_data_long) {
        exit(1);
    }

    std::cout << "Double data: " << in_data_double << " (read) vs " << out_data_double << " (expected) " << std::endl;
    if (in_data_double != out_data_double) {
        exit(1);
    }

    std::cout << "Bool data: " << in_data_bool << " (read) vs " << out_data_bool << " (expected) " << std::endl;
    if (in_data_bool != out_data_bool) {
        exit(1);
    }

    // Close archive
    binaryReader.close();

}
