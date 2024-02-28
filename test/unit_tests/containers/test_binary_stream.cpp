/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2024 OPTIMAD engineering Srl
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

#include <array>

#define BITPIT_UNIT_TEST_SUITE_NAME unit_test_containers_binary_stream
#include "helpers/unit_test.hpp"

using namespace bitpit;

BOOST_AUTO_TEST_SUITE(binaryStream)

BOOST_AUTO_TEST_CASE(binary_stream_array_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::array<double, 3> EXPECTED_RESULT = {{1., 2., 3.}};

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::array<double, 3> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::array into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_vector_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::vector<double> EXPECTED_RESULT;

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::vector<double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::vector into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_vector_2)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::vector<double> EXPECTED_RESULT;
    EXPECTED_RESULT.push_back(1.);
    EXPECTED_RESULT.push_back(2.);
    EXPECTED_RESULT.push_back(3.);

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::vector<double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::vector into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_pair_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::pair<int, double> EXPECTED_RESULT = std::make_pair(1, 11.);

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::pair<int, double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::pair into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_map_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::map<int, double> EXPECTED_RESULT;
    EXPECTED_RESULT[0] = 11.;
    EXPECTED_RESULT[1] = 22.;
    EXPECTED_RESULT[2] = 33.;

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::map<int, double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::map into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_unordered_map_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::unordered_map<int, double> EXPECTED_RESULT;
    EXPECTED_RESULT[0] = 11.;
    EXPECTED_RESULT[1] = 22.;
    EXPECTED_RESULT[2] = 33.;

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::unordered_map<int, double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::unordered_map into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_unordered_set_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::unordered_set<double> EXPECTED_RESULT;
    EXPECTED_RESULT.insert(0.);
    EXPECTED_RESULT.insert(1.),
    EXPECTED_RESULT.insert(2.);

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::unordered_set<double> result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::unordered_set into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_double_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    double EXPECTED_RESULT = 108.;

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    double result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a double into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_long_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    long EXPECTED_RESULT = 108;

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    long result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a long into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_string_1)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::string EXPECTED_RESULT = "";

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::string result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::string into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_CASE(binary_stream_string_2)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::string EXPECTED_RESULT = "test_string";

    OBinaryStream outputStream;
    outputStream << EXPECTED_RESULT;

    std::string result;
    IBinaryStream inputStream(outputStream.data(), outputStream.getSize());
    inputStream >> result;
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Streaming a std::string into/from a binary stream failed unit test.");
    }
}

BOOST_AUTO_TEST_SUITE_END()
