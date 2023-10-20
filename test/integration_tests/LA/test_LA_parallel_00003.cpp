/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2023 OPTIMAD engineering Srl
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
#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_LA.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing reading system solver from file.
*
* \param rank is the rank of the process
* \param nProcs is the number of processes
*/
int subtest_001(int rank, int nProcs)
{
    int nRows;
    int nCols;
    int nNZ;
    if (nProcs > 1 && rank == 0) {
        nRows = 0;
        nCols = 0;
        nNZ   = 0;
    } else {
        nRows = 10;
        nCols = 10;
        nNZ   = 10;
    }

    int blockSize = 1;

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    std::vector<long> rowPattern(1);
    std::vector<double> rowValues(1);

    SparseMatrix matrix(MPI_COMM_WORLD, true, blockSize, nRows, nCols, nNZ);

    int rowOffset = matrix.getRowGlobalOffset();
    int colOffset = matrix.getColGlobalOffset();

    for (int i = 0; i < nRows; ++i) {
        rowPattern[0] = colOffset + i;
        rowValues[0]  = 1. / (double) (rowOffset + i + 1);

        matrix.addRow(rowPattern, rowValues);
    }

    matrix.assembly();

    // Build system
    log::cout() << "Building system..." << std::endl;

    SystemSolver system;
    system.assembly(matrix);

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

    // Dump system
    log::cout() << "Dumping system..." << std::endl;

    system.dump(".", "test_parallel_00003_");

    // Restore system
    log::cout() << "Restoring system..." << std::endl;

    SystemSolver restoredSystem("", false);
    restoredSystem.restore(MPI_COMM_WORLD, true, ".", "test_parallel_00003_");

    long restoredGlobalNRows = restoredSystem.getRowGlobalCount();
    long globalNRows = system.getRowGlobalCount();

    if (restoredGlobalNRows != globalNRows) {
        throw std::runtime_error("  The restored system has not the same global number of rows of the original one.");
    }

    nRows = restoredSystem.getRowCount();
    nCols = restoredSystem.getColCount();

    // Solve restored system
    log::cout() << "Solving restored system..." << std::endl;

    restoredSystem.solve();

    log::cout() << std::setprecision(16) << std::scientific;

    if (nRows > 0) {
        const double *solution = restoredSystem.getSolutionRawReadPtr();
        for (int i = 0; i < nRows; ++i) {
            log::cout() << "  Solution[" << i << "] = " << solution[i] << std::endl;

            double expectedSolution = matrix.getRowGlobalOffset() + i + 1;
            if (!utils::DoubleFloatingEqual()(solution[i], expectedSolution, 10)) {
                log::cout() << "  Expected solution[" << i << "] = " << expectedSolution << std::endl;
                log::cout() << "  Error[" << i << "] = " << (expectedSolution - solution[i]) << std::endl;
                throw std::runtime_error("  The solution of the system doesn't match the expected one.");
            }
        }
        restoredSystem.restoreSolutionRawReadPtr(solution);
    } else {
        log::cout() << "  System matrix is empty on this process" << std::endl;
    }

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    // Initialize the logger
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing basic features of parallel octree patches" << std::endl;

    int status;
    try {
        status = subtest_001(rank, nProcs);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
