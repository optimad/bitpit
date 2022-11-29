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

#include "bitpit_IO.hpp"
#include "bitpit_LA.hpp"

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include <array>
#include <vector>

using namespace bitpit;
using namespace bitpit::linearalgebra;

/*!
* Subtest 001
*
* Testing solution of block linear system (block size equal two) assembled using
* a block assembler (same block size).
*/
int subtest_001()
{
    log::cout() << "  >> Testing solution of block linear system assembled using a block assembler." << std::endl;

    int blockSize = 2;
    int nRows     = 5;
    int nCols     = 5;
    int nNZ       = 9;

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    int nMaxNZ = 2;
    std::vector<long> rowPattern(nMaxNZ);
    std::vector<double> rowValues(nMaxNZ * blockSize * blockSize);

#if BITPIT_ENABLE_MPI==1
    SparseMatrix matrix(MPI_COMM_WORLD, false, blockSize, nRows, nCols, nNZ);
#else
    SparseMatrix matrix(blockSize, nRows, nCols, nNZ);
#endif
    for (int row = 0; row < nRows; ++row) {
        // Number of non-zero blocks
        int nRowNZ = (row != 2) ? 2 : 1;

        // Set pattern
        rowPattern.resize(nRowNZ);
        rowPattern[0] = row;
        if (nRowNZ == 2) {
            rowPattern[1] = nRows - row - 1;
        }

        // Set values
        std::fill(rowValues.begin(), rowValues.end(), 0.);
        for (int k  = 0; k < blockSize; ++k) {
            int ib0 = 0;
            int jb0 = 0;
            int i0  = k;
            int j0  = (k == 0) ? 0 : 1;
            int ij0 = linearIndexRowMajor(ib0 * blockSize + i0, jb0 * blockSize + j0, blockSize, rowPattern.size() * blockSize);
            rowValues[ij0] = (row * blockSize + k + 1);

            int ib1 = 0;
            int jb1 = (row != static_cast<int>(std::ceil(0.5 * nRows)) - 1) ? 1 : 0;
            int i1  = k;
            int j1  = (k == 0) ? 1 : 0;
            int ij1 = linearIndexRowMajor(ib1 * blockSize + i1, jb1 * blockSize + j1, blockSize, rowPattern.size() * blockSize);
            rowValues[ij1] = 11 * (row * blockSize + k + 1);
        }

        // Add row
        matrix.addRow(rowPattern, rowValues);
    }
    matrix.assembly();

    // Build system
    log::cout() << "Building system..." << std::endl;

    SystemSolver system(true, false);
    system.assembly(matrix);

    int vectorSize = system.getColElementCount();

    double *rhs = system.getRHSRawPtr();
    rhs[0] = 1101;
    rhs[1] =  895;
    rhs[2] =  713;
    rhs[3] =  555;
    rhs[4] =  421;
    rhs[5] =  311;
    rhs[6] =  225;
    rhs[7] =  163;
    rhs[8] =  125;
    rhs[9] =  111;
    system.restoreRHSRawPtr(rhs);

    double *initialSolution = system.getSolutionRawPtr();
    for (int i = 0; i < vectorSize; ++i) {
        initialSolution[i] = 0;
    }
    system.restoreSolutionRawPtr(initialSolution);

    // Solve system
    log::cout() << "Solving transposed system..." << std::endl;

    system.solve();

    log::cout() << std::setprecision(16) << std::scientific;

    const double *solution = system.getSolutionRawReadPtr();
    for (int i = 0; i < vectorSize; ++i) {
        log::cout() << "  Solution[" << i << "] = " << solution[i] << std::endl;

        double expectedSolution = i + 1;
        if (!utils::DoubleFloatingEqual()(solution[i], expectedSolution, 1e-10)) {
            log::cout() << "  Expected solution[" << i << "] = " << expectedSolution << std::endl;
            log::cout() << "  Error[" << i << "] = " << (expectedSolution - solution[i]) << std::endl;
            throw std::runtime_error("  The solution of the system doesn't match the expected one.");
        }
    }
    system.restoreSolutionRawReadPtr(solution);

    return 0;
}


/*!
* Subtest 002
*
* Testing solution of linear system (block size equal one) assembled using a
* block assembler (block size equal one).
*/
int subtest_002()
{
    log::cout() << "  >> Testing solution of linear system assembled using a block assembler." << std::endl;

    int blockSize = 2;
    int nRows     = 5;
    int nCols     = 5;
    int nNZ       = 9;

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    int nMaxNZ = 2;
    std::vector<long> rowPattern(nMaxNZ);
    std::vector<double> rowValues(nMaxNZ * blockSize * blockSize);

#if BITPIT_ENABLE_MPI==1
    SparseMatrix matrix(MPI_COMM_WORLD, false, blockSize, nRows, nCols, nNZ);
#else
    SparseMatrix matrix(blockSize, nRows, nCols, nNZ);
#endif
    for (int row = 0; row < nRows; ++row) {
        // Number of non-zero blocks
        int nRowNZ = (row != 2) ? 2 : 1;

        // Set pattern
        rowPattern.resize(nRowNZ);
        rowPattern[0] = row;
        if (nRowNZ == 2) {
            rowPattern[1] = nRows - row - 1;
        }

        // Set values
        std::fill(rowValues.begin(), rowValues.end(), 0.);
        for (int k  = 0; k < blockSize; ++k) {
            int ib0 = 0;
            int jb0 = 0;
            int i0  = k;
            int j0  = (k == 0) ? 0 : 1;
            int ij0 = linearIndexRowMajor(ib0 * blockSize + i0, jb0 * blockSize + j0, blockSize, rowPattern.size() * blockSize);
            rowValues[ij0] = (row * blockSize + k + 1);

            int ib1 = 0;
            int jb1 = (row != static_cast<int>(std::ceil(0.5 * nRows)) - 1) ? 1 : 0;
            int i1  = k;
            int j1  = (k == 0) ? 1 : 0;
            int ij1 = linearIndexRowMajor(ib1 * blockSize + i1, jb1 * blockSize + j1, blockSize, rowPattern.size() * blockSize);
            rowValues[ij1] = 11 * (row * blockSize + k + 1);
        }

        // Add row
        matrix.addRow(rowPattern, rowValues);
    }
    matrix.assembly();

    // Build system
    log::cout() << "Building system..." << std::endl;

    SystemSolver system(true, false);
    system.assembly(matrix);

    int vectorSize = system.getColElementCount();

    double *rhs = system.getRHSRawPtr();
    rhs[0] = 1101;
    rhs[1] =  895;
    rhs[2] =  713;
    rhs[3] =  555;
    rhs[4] =  421;
    rhs[5] =  311;
    rhs[6] =  225;
    rhs[7] =  163;
    rhs[8] =  125;
    rhs[9] =  111;
    system.restoreRHSRawPtr(rhs);

    double *initialSolution = system.getSolutionRawPtr();
    for (int i = 0; i < vectorSize; ++i) {
        initialSolution[i] = 0;
    }
    system.restoreSolutionRawPtr(initialSolution);

    // Solve system
    log::cout() << "Solving transposed system..." << std::endl;

    system.solve();

    log::cout() << std::setprecision(16) << std::scientific;

    const double *solution = system.getSolutionRawReadPtr();
    for (int i = 0; i < vectorSize; ++i) {
        log::cout() << "  Solution[" << i << "] = " << solution[i] << std::endl;

        double expectedSolution = i + 1;
        if (!utils::DoubleFloatingEqual()(solution[i], expectedSolution, 1e-10)) {
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
    log::cout() << "Testing solution of block linear systems..." << std::endl;

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
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
