/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#   include <mpi.h>
#endif

#include "bitpit_IO.hpp"
#include "bitpit_LA.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing solution of linear systems.
*/
int subtest_001()
{
    int nRows = 10;
    int nCols = 10;
    int nNZ   = 10;

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    std::vector<long> rowPattern(1);
    std::vector<double> rowValues(1);

#if BITPIT_ENABLE_MPI==1
    SparseMatrix matrix(MPI_COMM_WORLD, false, nRows, nCols, nNZ);
#else
    SparseMatrix matrix(nRows, nCols, nNZ);
#endif
    for (int i = 0; i < nRows; ++i) {
        rowPattern[0] = i;
        rowValues[0]  = 1. / (double) (i + 1);

        matrix.addRow(rowPattern, rowValues);
    }
    matrix.assembly();

    // Build system
    log::cout() << "Building system..." << std::endl;

#if BITPIT_ENABLE_MPI==1
    SystemSolver system(MPI_COMM_WORLD);
#else
    SystemSolver system;
#endif
    system.initialize(matrix);

    double *rhs = system.getRHSRawPtr();
    for (int i = 0; i < nCols; ++i) {
        rhs[i] = 1.;
    }
    system.restoreRHSRawPtr(rhs);

    double *initialSolution = system.getSolutionRawPtr();
    for (int i = 0; i < nRows; ++i) {
        initialSolution[i] = 0;
    }
    system.restoreSolutionRawPtr(initialSolution);

    // Solve system
    log::cout() << "Solving system..." << std::endl;

    system.solve();

    log::cout() << std::setprecision(16) << std::scientific;

    const double *solution = system.getSolutionRawReadPtr();
    for (int i = 0; i < nRows; ++i) {
        log::cout() << "  Solution[" << i << "] = " << solution[i] << std::endl;

        double expectedSolution = i + 1;
        if (!utils::DoubleFloatingEqual()(solution[i], expectedSolution, 10)) {
            log::cout() << "  Expected solution[" << i << "] = " << expectedSolution << std::endl;
            log::cout() << "  Error[" << i << "] = " << (expectedSolution - solution[i]) << std::endl;
            throw std::runtime_error("  The solution of the system doesn't match the expected one.");
        }
    }
    system.restoreSolutionRawReadPtr(solution);

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

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Testing soltion of linear systems" << std::endl;

    int status;
    try {
        status = subtest_001();
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
