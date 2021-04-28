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

# ifndef __BITPIT_UNIT_TEST_VOLCARTESIAN_FIXTURES__
# define __BITPIT_UNIT_TEST_VOLCARTESIAN_FIXTURES__

#include <array>
#include <memory>

#include "bitpit_IO.hpp"
#include "bitpit_volcartesian.hpp"

/*!
 * Test patch constants for the unit tests.
 */
struct TestPatchConstants {

    const int DIMENSION = 3;

    const double ORIGIN_X = 0.;
    const double ORIGIN_Y = 0.;
    const double ORIGIN_Z = 0.;

    const double LENGTH_X = 10.;
    const double LENGTH_Y = 20.;
    const double LENGTH_Z = 30.;

    const int N_CELLS_X = 5.;
    const int N_CELLS_Y = 6.;
    const int N_CELLS_Z = 7.;
    const int N_CELLS   = N_CELLS_X * N_CELLS_Y * N_CELLS_Z;

    const int N_VERTICES_X = N_CELLS_X + 1;
    const int N_VERTICES_Y = N_CELLS_Y + 1;
    const int N_VERTICES_Z = N_CELLS_Z + 1;
    const int N_VERTICES   = N_VERTICES_X * N_VERTICES_Y * N_VERTICES_Z;

    const int N_INTERFACES_X = (N_CELLS_X + 1) * (N_CELLS_Y + 0) * (N_CELLS_Z + 0);
    const int N_INTERFACES_Y = (N_CELLS_X + 0) * (N_CELLS_Y + 1) * (N_CELLS_Z + 0);
    const int N_INTERFACES_Z = (N_CELLS_X + 0) * (N_CELLS_Y + 0) * (N_CELLS_Z + 1);
    const int N_INTERFACES   = N_INTERFACES_X + N_INTERFACES_Y + N_INTERFACES_Z;

    const double SPACING_X = (double) LENGTH_X / N_CELLS_X;
    const double SPACING_Y = (double) LENGTH_Y / N_CELLS_Y;
    const double SPACING_Z = (double) LENGTH_Z / N_CELLS_Z;

    const double CELL_VOLUME = SPACING_X * SPACING_Y * SPACING_Z;
    const double CELL_SIZE   = std::pow(CELL_VOLUME, 1. / DIMENSION);

    const double INTERFACE_AREA_X = SPACING_Y * SPACING_Z;
    const double INTERFACE_AREA_Y = SPACING_X * SPACING_Z;
    const double INTERFACE_AREA_Z = SPACING_X * SPACING_Y;

};

/*!
 * Base test patch for the unit tests.
 */
struct BaseTestPatch : TestPatchConstants {

    std::unique_ptr<bitpit::VolCartesian> patch;

protected:
    BaseTestPatch(bitpit::VolCartesian::MemoryMode memoryMode)
    {
        if (memoryMode == bitpit::VolCartesian::MEMORY_NORMAL) {
            bitpit::log::cout() << "Initializing test patch in normal memory mode..." << std::endl;
        } else if (memoryMode == bitpit::VolCartesian::MEMORY_LIGHT) {
            bitpit::log::cout() << "Initializing test patch in light memory mode..." << std::endl;
        }

        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        std::array<double, 3> lengths = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};
        std::array<int, 3> nCells = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};
        patch = std::unique_ptr<bitpit::VolCartesian>(new bitpit::VolCartesian(dimension, origin, lengths, nCells));
        patch->switchMemoryMode(memoryMode);
        if (memoryMode == bitpit::VolCartesian::MEMORY_NORMAL) {
            patch->update();
        }
    }

};

/*!
 * Test patch for unit tests with in normal memory mode.
 */
struct NormalTestPatch : BaseTestPatch {

    NormalTestPatch()
        : BaseTestPatch(bitpit::VolCartesian::MEMORY_NORMAL)
    {
    }

};

/*!
 * Test patch for unit tests with in light memory mode.
 */
struct LightTestPatch : BaseTestPatch {

    LightTestPatch()
        : BaseTestPatch(bitpit::VolCartesian::MEMORY_LIGHT)
    {
    }

};

#endif
