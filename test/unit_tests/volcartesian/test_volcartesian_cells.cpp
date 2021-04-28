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

#include "test_volcartesian_fixtures.hpp"

using namespace bitpit;

BOOST_AUTO_TEST_SUITE(normalMemoryMode)

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellCount_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_CELLS;

    long result = patch->getCellCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellCount_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};

    for (int d = 0; d < 3; ++d) {
        int result = patch->getCellCount(d);
        if (result != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getCellCount(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellType_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::VOXEL;

    ElementType result = patch->getCellType();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellType() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellType_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::VOXEL;

    ElementType result = patch->getCellType(0);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellType(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellVolume_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_VOLUME;

    double result = patch->evalCellVolume(0);
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellVolume(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellVolume_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_VOLUME;

    double result = patch->evalCellVolume({{0, 0, 0}});
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellVolume(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellSize_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_SIZE;

    double result = patch->evalCellSize(0);
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellSize(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellSize_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_SIZE;

    double result = patch->evalCellSize({{0, 0, 0}});
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellSize(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellCentroid_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

    std::array<double, 3> result = patch->evalCellCentroid(0);
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalCellCentroid(long id) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_evalCellCentroid_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

    std::array<double, 3> result = patch->evalCellCentroid({{0, 0, 0}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalCellCentroid(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellCentroids, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

    for (int d = 0; d < 3; ++d) {
        double result = patch->getCellCentroids(d)[0];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getCellCentroids(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellLinearId_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 101;

    long result = patch->getCellLinearId(1, 2, 3);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellLinearId(int i, int j, int k) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellLinearId_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 101;

    long id = patch->getCellLinearId({{1, 2, 3}});
    if (id != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_getCellCartesianId, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> cartesianId = patch->getCellCartesianId(101);
    for (int d = 0; d < 3; ++d) {
        if (cartesianId[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getCellCartesianId(long idx) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_isCellCartesianIdValid, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isCellCartesianIdValid({{1, 2, 3}});
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'isCellCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_locateClosestCell, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 102;

    long result = patch->locateClosestCell({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'locateClosestCell(std::array<double, 3> const &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_locateClosestCellCartesian, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

    std::array<int, 3> result = patch->locateClosestCellCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'locateClosestCellCartesian(std::array<double, 3> const &point) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractCellSubSet_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

    std::vector<long> result = patch->extractCellSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractCellSubSet_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

    std::vector<long> result = patch->extractCellSubSet(101, 107);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(int idMin, int idMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normalMemoryMode_extractCellSubSet_3, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{102, 107}};

    const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
    std::vector<long> result = patch->extractCellSubSet(pointMin, pointMax);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(lightMemoryMode)

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellCount_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = N_CELLS;

    long result = patch->getCellCount();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellCount() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellCount_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};

    for (int d = 0; d < 3; ++d) {
        int result = patch->getCellCount(d);
        if (result != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getCellCount(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellType_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::VOXEL;

    ElementType result = patch->getCellType();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellType() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellType_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const ElementType EXPECTED_RESULT = ElementType::VOXEL;

    ElementType result = patch->getCellType(0);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellType(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_evalCellVolume, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_VOLUME;

    double result = patch->evalCellVolume(0);
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellVolume(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_evalCellSize, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = CELL_SIZE;

    double result = patch->evalCellSize(0);
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'evalCellSize(long id) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_evalCellCentroid, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

    std::array<double, 3> result = patch->evalCellCentroid(0);
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'evalCellCentroid(long id) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellCentroids, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

    for (int d = 0; d < 3; ++d) {
        double result = patch->getCellCentroids(d)[0];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getCellCentroids(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellLinearId_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 101;

    long result = patch->getCellLinearId(1, 2, 3);
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellLinearId(int i, int j, int k) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellLinearId_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 101;

    long id = patch->getCellLinearId({{1, 2, 3}});
    if (id != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getCellLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_getCellCartesianId, LightTestPatch)
{
   BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

   const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

    std::array<int, 3> cartesianId = patch->getCellCartesianId(101);
    for (int d = 0; d < 3; ++d) {
        if (cartesianId[d] != EXPECTED_RESULT[d]) {
            throw std::runtime_error("Function 'getCellCartesianId(long idx) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_isCellCartesianIdValid, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isCellCartesianIdValid({{1, 2, 3}});
    if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'isCellCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_locateClosestCell, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 102;

    long result = patch->locateClosestCell({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'locateClosestCell(std::array<double, 3> const &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_locateClosestCellCartesian, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

    std::array<int, 3> result = patch->locateClosestCellCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'locateClosestCellCartesian(std::array<double, 3> const &point) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractCellSubSet_1, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

    std::vector<long> result = patch->extractCellSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractCellSubSet_2, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

    std::vector<long> result = patch->extractCellSubSet(101, 107);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(int idMin, int idMax) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(lightMemoryMode_extractCellSubSet_3, LightTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<long, 2> EXPECTED_RESULT = {{102, 107}};

    const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
    std::vector<long> result = patch->extractCellSubSet(pointMin, pointMax);
    for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
        if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
            throw std::runtime_error("Function 'extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
