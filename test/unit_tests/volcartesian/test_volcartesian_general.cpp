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

#define BITPIT_UNIT_TEST_SUITE_NAME unit_test_volcartesian_general
#include "helpers/unit_test.hpp"

#include "bitpit_volcartesian.hpp"

#include "test_volcartesian_fixtures.hpp"

using namespace bitpit;

BOOST_AUTO_TEST_SUITE(normalMemoryMode)

BOOST_FIXTURE_TEST_CASE(method_clone, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    std::unique_ptr<PatchKernel> clonedPatch = patch->clone();
    if (!clonedPatch) {
        throw std::runtime_error("Function 'clone()' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_reset, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    patch->reset();
    if (patch->getCellCount() != 0) {
        throw std::runtime_error("Function 'reset()' failed unit test.");
    } else if (patch->getVertexCount() != 0) {
        throw std::runtime_error("Function 'reset()' failed unit test.");
    } else if (patch->getInterfaceCount() != 0) {
        throw std::runtime_error("Function 'reset()' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_getSpacing_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

    std::array<double, 3> result = patch->getSpacing();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getSpacing() const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_getSpacing_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

    for (int d = 0; d < 3; ++d) {
        double result = patch->getSpacing(d);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getSpacing(int direction) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_switchMemoryMode, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const VolCartesian::MemoryMode EXPECTED_RESULT = VolCartesian::MEMORY_LIGHT;

    patch->switchMemoryMode(VolCartesian::MEMORY_LIGHT);
    VolCartesian::MemoryMode result = patch->getMemoryMode();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'switchMemoryMode(MemoryMode mode)' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_getMemoryMode, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const VolCartesian::MemoryMode EXPECTED_RESULT = VolCartesian::MEMORY_NORMAL;

    VolCartesian::MemoryMode result = patch->getMemoryMode();
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'getMemoryMode() const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_isPointInside_1, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isPointInside({{ORIGIN_X + 0.5 * LENGTH_X, ORIGIN_Y + 0.5 * LENGTH_Y, ORIGIN_Z + 0.5 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'isPointInside(const std::array<double, 3> &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_isPointInside_2, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const bool EXPECTED_RESULT = true;

    bool result = patch->isPointInside(102, {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'isPointInside(const std::array<double, 3> &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_locatePoint, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const long EXPECTED_RESULT = 102;

    long result = patch->locatePoint({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    if (result != EXPECTED_RESULT) {
        throw std::runtime_error("Function 'locatePoint(const std::array<double, 3> &point) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_locatePointCartesian, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

    std::array<int, 3> result = patch->locatePointCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'locatePointCartesian(const std::array<double, 3> &point) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_getOrigin, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};

    std::array<double, 3> result = patch->getOrigin();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getOrigin() const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_setOrigin, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{10., 10., 10.}};

    patch->setOrigin(EXPECTED_RESULT);
    std::array<double, 3> result = patch->getOrigin();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'setOrigin(const std::array<double, 3> &origin)' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_translate, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{ORIGIN_X + 10., ORIGIN_Y + 10., ORIGIN_Z + 10.}};

    patch->translate({{10., 10., 10.}});
    std::array<double, 3> result = patch->getOrigin();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'translate(const std::array<double, 3> &translation)' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_getLengths, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};

    std::array<double, 3> result = patch->getLengths();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'getLengths() const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_setLengths, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};

    patch->setLengths(EXPECTED_RESULT);
    std::array<double, 3> result = patch->getLengths();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'setLengths(const std::array<double, 3> &lengths)' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_scale, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 3> EXPECTED_RESULT = {{10. * LENGTH_X, 10. * LENGTH_Y, 10. * LENGTH_Z}};

    patch->scale({{10., 10., 10.}}, {{10., 10., 10.}});
    std::array<double, 3> result = patch->getLengths();
    for (int d = 0; d < 3; ++d) {
        if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
            throw std::runtime_error("Function 'scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center)' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_convertToVertexData, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = 56.5;

    std::vector<double> cellData(patch->getCellCount());
    for (long n = 0; n < patch->getCellCount(); ++n) {
        cellData[n] = n;
    }

    std::vector<double> result = patch->convertToVertexData(cellData);
    if (!utils::DoubleFloatingEqual()(result[101], EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'convertToVertexData(const std::vector<double> &cellData) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_convertToCellData, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const double EXPECTED_RESULT = 163.5;

    std::vector<double> vertexData(patch->getVertexCount());
    for (long n = 0; n < patch->getVertexCount(); ++n) {
        vertexData[n] = n;
    }

    std::vector<double> result = patch->convertToCellData(vertexData);
    if (!utils::DoubleFloatingEqual()(result[101], EXPECTED_RESULT)) {
        throw std::runtime_error("Function 'convertToCellData(const std::vector<double> &vertexData) const' failed unit test.");
    }
}

BOOST_FIXTURE_TEST_CASE(method_linearCellInterpolation, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 8> EXPECTED_RESULT = {{0.07, 0.21, 0.0175, 0.0525, 0.13, 0.39, 0.0325, 0.0975}};

    std::vector<double> vertexData(patch->getVertexCount());
    for (long n = 0; n < patch->getVertexCount(); ++n) {
        vertexData[n] = n;
    }

    std::vector<int> stencil;
    std::vector<double> weights;
    std::array<double, 3> point = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    int stencilSize = patch->linearCellInterpolation(point, &stencil, &weights);
    for (int i = 0; i < stencilSize; ++i) {
        double result = weights[i];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[i])) {
            throw std::runtime_error("Function 'linearCellInterpolation(std::array<double, 3> &point, std::vector<int> &stencil, std::vector<double> &weights) const' failed unit test.");
        }
    }
}

BOOST_FIXTURE_TEST_CASE(method_linearVertexInterpolation, NormalTestPatch)
{
    BITPIT_UNIT_TEST_DISPLAY_NAME(log::cout());

    const std::array<double, 8> EXPECTED_RESULT = {{0.19125, 0.06375, 0.44625, 0.14875, 0.03375, 0.01125, 0.07875, 0.02625}};

    std::vector<double> vertexData(patch->getVertexCount());
    for (long n = 0; n < patch->getVertexCount(); ++n) {
        vertexData[n] = n;
    }

    std::vector<int> stencil;
    std::vector<double> weights;
    std::array<double, 3> point = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
    int stencilSize = patch->linearVertexInterpolation(point, &stencil, &weights);
    for (int i = 0; i < stencilSize; ++i) {
        double result = weights[i];
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[i])) {
            throw std::runtime_error("Function 'linearVertexInterpolation(std::array<double, 3> &point, std::vector<int> &stencil, std::vector<double> &weights) const' failed unit test.");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
