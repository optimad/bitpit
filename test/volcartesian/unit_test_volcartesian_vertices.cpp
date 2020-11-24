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

#define BITPIT_UNIT_TEST_SUITE_NAME unit_test_volcartesian_cells
#include "helpers/unit_test.hpp"

#include "bitpit_common.hpp"
#include "bitpit_volcartesian.hpp"

#include "unit_test_volcartesian_fixtures.hpp"

using namespace bitpit;

BOOST_AUTO_TEST_SUITE(normalMemoryMode)

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCount_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_VERTICES;

    long result = patch->getVertexCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCount_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{N_VERTICES_X, N_VERTICES_Y, N_VERTICES_Z}};

    for (int d = 0; d < 3; ++d) {
        int result = patch->getVertexCount(d);
        if (result != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCount(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalVertexCoords_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{1 * SPACING_X, 2 * SPACING_Y, 3 * SPACING_Z}};

    std::array<double, 3> result = patch->evalVertexCoords(139);
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalVertexCoords(long id) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalVertexCoords_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{1 * SPACING_X, 2 * SPACING_Y, 3 * SPACING_Z}};

    std::array<double, 3> result = patch->evalVertexCoords({{1, 2, 3}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalVertexCoords(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCoords, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

    for (int d = 0; d < 3; ++d) {
        double result = patch->getVertexCoords(d)[1];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getVertexCoords(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexLinearId_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 139;

    long result = patch->getVertexLinearId(1, 2, 3);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexLinearId(int i, int j, int k) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexLinearId_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 139;

    long result = patch->getVertexLinearId({{1, 2, 3}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCartesianId_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId(139);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(long idx) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCartesianId_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId(101, 0);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(long cellIdx, int vertex) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getVertexCartesianId_3, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId({{1, 2, 3}}, 0);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(const std::array<int, 3> &cellIjk, int vertex) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_isVertexCartesianIdValid, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isVertexCartesianIdValid({{1, 2, 3}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'isVertexCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_locateClosestVertex, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 146;

    long result = patch->locateClosestVertex({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'locateClosestVertex(std::array<double, 3> const &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_locateClosestVertexCartesian, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{2, 3, 3}};

    std::array<int, 3> result = patch->locateClosestVertexCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'locateClosestVertexCartesian(std::array<double, 3> const &point) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractVertexSubSet_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{139, 140, 145, 146}};

    std::vector<long> result = patch->extractVertexSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractVertexSubSet_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{101, 107}};

    std::vector<long> result = patch->extractVertexSubSet(101, 107);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(int idxMin, int idxMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractVertexSubSet_3, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{140, 146}};

    const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
    std::vector<long> result = patch->extractVertexSubSet(pointMin, pointMax);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(lightMemoryMode)

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCount_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_VERTICES;

    long result = patch->getVertexCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCount_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{N_VERTICES_X, N_VERTICES_Y, N_VERTICES_Z}};

    for (int d = 0; d < 3; ++d) {
        int result = patch->getVertexCount(d);
        if (result != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCount(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_evalVertexCoords, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{1 * SPACING_X, 2 * SPACING_Y, 3 * SPACING_Z}};

    std::array<double, 3> result = patch->evalVertexCoords(139);
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalVertexCoords(long id) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCoords, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

    for (int d = 0; d < 3; ++d) {
        double result = patch->getVertexCoords(d)[1];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getVertexCoords(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexLinearId_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 139;

    long result = patch->getVertexLinearId(1, 2, 3);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexLinearId(int i, int j, int k) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexLinearId_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 139;

    long result = patch->getVertexLinearId({{1, 2, 3}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getVertexLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCartesianId_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId(139);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(long idx) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCartesianId_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId(101, 0);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(long cellIdx, int vertex) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getVertexCartesianId_3, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> result = patch->getVertexCartesianId({{1, 2, 3}}, 0);
    for (int d = 0; d < 3; ++d) {
        if (result[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getVertexCartesianId(const std::array<int, 3> &cellIjk, int vertex) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_isVertexCartesianIdValid, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isVertexCartesianIdValid({{1, 2, 3}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'isVertexCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_locateClosestVertex, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 146;

    long result = patch->locateClosestVertex({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'locateClosestVertex(std::array<double, 3> const &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_locateClosestVertexCartesian, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{2, 3, 3}};

    std::array<int, 3> result = patch->locateClosestVertexCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'locateClosestVertexCartesian(std::array<double, 3> const &point) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractVertexSubSet_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{139, 140, 145, 146}};

    std::vector<long> result = patch->extractVertexSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractVertexSubSet_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{101, 107}};

    std::vector<long> result = patch->extractVertexSubSet(101, 107);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(int idxMin, int idxMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractVertexSubSet_3, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{140, 146}};

    const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
    std::vector<long> result = patch->extractVertexSubSet(pointMin, pointMax);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
