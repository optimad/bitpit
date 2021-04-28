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

#include <array>
#include <vector>

#define BITPIT_UNIT_TEST_SUITE_NAME unit_test_volcartesian_cells
#include "helpers/unit_test.hpp"

#include "bitpit_common.hpp"
#include "bitpit_volcartesian.hpp"

#include "test_volcartesian_fixtures.hpp"

using namespace bitpit;

BOOST_AUTO_TEST_SUITE(normalMemoryMode)

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getInterfaceCount, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_INTERFACES;

    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    long result = patch->getInterfaceCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getInterfaceType_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::PIXEL;

    ElementType result = patch->getInterfaceType();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getInterfaceType_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::PIXEL;

    ElementType result = patch->getInterfaceType(0);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_resetInterfaces, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 0;

    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    patch->resetInterfaces();
    long result = patch->getInterfaces().size();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'resetInterfaces()' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_initializeInterfaces, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_INTERFACES;

    patch->resetInterfaces();
    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    long result = patch->getInterfaces().size();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'initializeInterfaces()' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalInterfaceArea, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = INTERFACE_AREA_X;

    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    double result = patch->evalInterfaceArea(0);
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalInterfaceArea(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalInterfaceNormal, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0., 0., 1.}};

    patch->initializeAdjacencies();
    patch->initializeInterfaces();
    std::array<double, 3> result = patch->evalInterfaceNormal(10);
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalInterfaceNormal(long id) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(lightMemoryMode)

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getInterfaceCount_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_INTERFACES;

    long result = patch->getInterfaceCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getInterfaceType_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::PIXEL;

    ElementType result = patch->getInterfaceType();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getInterfaceType_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::PIXEL;

    ElementType result = patch->getInterfaceType(0);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
    }
}

BOOST_AUTO_TEST_SUITE_END()
