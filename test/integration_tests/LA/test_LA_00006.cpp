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
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_LA.hpp"

using namespace bitpit;

const double SOLVER_RTOL = 1e-7;

/*!
 * Subtest 001
 *
 * Testing solution of a split block linear system and its dump.
 */
int subtest_001()
{
    log::cout() << "  >> Testing restore of a split block linear system and its solution." << std::endl;

    int nRows = 10;
    int nCols = 10;
    int nNZ   = 10;

    int blockSize = 5;

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    std::vector<long> rowPattern(1);
    std::vector<double> rowValues(blockSize * blockSize);

#if BITPIT_ENABLE_MPI==1
    SparseMatrix matrix(MPI_COMM_WORLD, true, blockSize, nRows, nCols, nNZ);

    int rowOffset = matrix.getRowGlobalOffset();
    int colOffset = matrix.getColGlobalOffset();
#else
    SparseMatrix matrix(blockSize, nRows, nCols, nNZ);

    int rowOffset = 0;
    int colOffset = 0;
#endif

    for (int i = 0; i < nRows; ++i) {
        rowPattern[0] = colOffset + i;
        for (int k = 0; k < blockSize; ++k) {
            int kk = linearalgebra::linearIndexRowMajor(k, k, blockSize, blockSize);

            rowValues[kk]  = 1. / (double) (blockSize * (rowOffset + i) + k + 1);
        }

        matrix.addRow(rowPattern, rowValues);
    }

    matrix.assembly();

    // Build solver
    log::cout() << "Building solver..." << std::endl;

    bool multigrid = false;
    bool debug     = false;
    bool transpose = false;
    SplitSystemSolver solver(transpose, multigrid, debug);

    std::vector<int> splitSizes(2);
    splitSizes[0] = 4;
    splitSizes[1] = blockSize - splitSizes[0];

    solver.assembly(matrix, SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_LOWER, splitSizes);

    double *rhs = solver.getRHSRawPtr();
    for (int i = 0; i < nCols; ++i) {
        int field = 0;
        for (std::size_t split = 0; split < splitSizes.size(); ++split) {
            for (int k = 0; k < splitSizes[split]; ++k) {
                rhs[blockSize * i + field] = std::pow(10, split);
                ++field;

            }
        }
    }
    solver.restoreRHSRawPtr(rhs);

    double *initialSolution = solver.getSolutionRawPtr();
    for (int i = 0; i < nRows; ++i) {
        initialSolution[i] = 0;
    }
    solver.restoreSolutionRawPtr(initialSolution);

    // Export system
    solver.exportMatrix("LA_example_0002_matrix.txt", SystemSolver::FILE_ASCII);
    solver.exportRHS("LA_example_0002_rhs.txt", SystemSolver::FILE_ASCII);
    solver.exportSolution("LA_example_0002_solution.txt", SystemSolver::FILE_ASCII);

    // Dump system
    log::cout() << "Dumping initialized linear system..." << std::endl;

    solver.dumpSystem("bitpit linear system", ".", "LA_example_0002_linear_system_");

    // Set KSP options
    for (int i = 0; i < solver.getSplitCount(); ++i) {
        KSPOptions &options = solver.getSplitKSPOptions(i);
        options.rtol = SOLVER_RTOL;
    }

    // Solve System
    log::cout() << "Solve linear system..." << std::endl;
    solver.solve();

    for (int split = 0; split < solver.getSplitCount(); ++split) {
        const KSPStatus &status = solver.getSplitKSPStatus(split);
        log::cout() << "  Split #" << split << " convergence reason: " << status.convergence << " in its: " << status.its << std::endl;
    }
    log::cout() << "Linear system solved." << std::endl;

    // Export solution
    log::cout() << "Export solution..." << std::endl;
    solver.exportSolution("LA_example_0002_linear_system_solution.dat", bitpit::SystemSolver::FILE_BINARY);

    // Compare computed solution with expected one
    //
    // Set comparison tolerance to the linear system relative tolerance
    log::cout() << "Comparing solutions..." << std::endl;

    log::cout() << "    Tolerance = " << SOLVER_RTOL << std::endl;

    std::size_t nRowElements = solver.getRowElementCount();

    const double *solution = solver.getSolutionRawReadPtr();
    for (std::size_t i = 0; i < nRowElements; ++i) {
        std::size_t globalDOF = rowOffset * blockSize + i;

        // Check for absolute tolerance only
        double expectedSolution;
        if (globalDOF == 0) {
            expectedSolution = 1.;
        } else if ((globalDOF == 1) || ((globalDOF + 1) % blockSize) != 0) {
            expectedSolution = static_cast<double>(globalDOF + 1);
        } else {
            expectedSolution = static_cast<double>(10 * (globalDOF + 1));
        }

        if (!utils::DoubleFloatingEqual()(expectedSolution - solution[i], 0.0, SOLVER_RTOL, SOLVER_RTOL)) {
            std::stringstream message;
            message << "Solutions do not match for global DOF #" << globalDOF << ".";
            message << " Expected solution is " << expectedSolution;
            message << ", current solution is " << solution[i] << ".";
            log::cout() << message.str() << std::endl;
            throw std::runtime_error("  The solution of the system doesn't match the expected one.");
        }
    }
    solver.restoreSolutionRawReadPtr(solution);

    return 0;
}

/*!
 * Subtest 002
 *
 * Testing restore of a split block linear system and its solution.
 */
int subtest_002()
{
    log::cout() << "  >> Testing restore of a split block linear system and its solution." << std::endl;

    // Restore solver with initial guess
    log::cout() << "Restoring linear system ..." << std::endl;
    bool multigrid = false;
    bool debug     = false;
    bool transpose = false;
    SplitSystemSolver solver(transpose, multigrid, debug);

#if BITPIT_ENABLE_MPI==1
    solver.restoreSystem(MPI_COMM_WORLD, ".", "LA_example_0002_linear_system_");
#else
    solver.restoreSystem(".", "LA_example_0002_linear_system_");
#endif
    log::cout() << "Linear system restored." << std::endl;

    // Export system
    solver.exportMatrix("LA_example_0002_restored_matrix.txt", SystemSolver::FILE_ASCII);
    solver.exportRHS("LA_example_0002_restored_rhs.txt", SystemSolver::FILE_ASCII);
    solver.exportSolution("LA_example_0002_restored_solution.txt", SystemSolver::FILE_ASCII);

    // Set KSP options
    for (int i = 0; i < solver.getSplitCount(); ++i) {
        KSPOptions &options = solver.getSplitKSPOptions(i);
        options.rtol = SOLVER_RTOL;
    }

    // Solve linear system
    log::cout() << "Solving linear system..." << std::endl;
    solver.solve();

    for (int split = 0; split < solver.getSplitCount(); ++split) {
        const KSPStatus &status = solver.getSplitKSPStatus(split);
        log::cout() << "  Split #" << split << " convergence reason: " << status.convergence << " in its: " << status.its << std::endl;
    }
    log::cout() << "Linear system solved." << std::endl;

    // Turn solution comparison on
    bool compareSolutions = true;
    if (compareSolutions) {
        // Store expected solution for comparison
        log::cout() << "Storing computed solution for comparison..." << std::endl;
        std::size_t nRowElements = solver.getRowElementCount();
        std::vector<double> expectedSolution(nRowElements, 0.0);
        double *computedSolution = solver.getSolutionRawPtr();
        for (std::size_t i = 0; i < nRowElements; ++i) {
            expectedSolution[i] = computedSolution[i];
        }
        solver.restoreSolutionRawPtr(computedSolution);

        // Import solution with solution from outside
        log::cout() << "Restoring solution..." << std::endl;
        solver.importSolution("LA_example_0002_linear_system_solution.dat");
        log::cout() << "Solution restored." << std::endl;

        // Compare computed solution with restored one
        //
        // Set comparison tolerance to the linear system relative tolerance
        log::cout() << "Comparing solutions..." << std::endl;

        log::cout() << "    Tolerance = " << SOLVER_RTOL << std::endl;

        const double *solution = solver.getSolutionRawReadPtr();
        for (std::size_t i = 0; i < nRowElements; ++i) {
            if (!utils::DoubleFloatingEqual()(expectedSolution[i] - solution[i], 0.0, SOLVER_RTOL, SOLVER_RTOL)) {
                std::stringstream message;
                message << "Solutions do not match for local DOF #" << i << ".";
                message << " Expected solution is " << expectedSolution[i];
                message << ", current solution is " << solution[i] << ".";
                log::cout() << message.str() << std::endl;
                throw std::runtime_error("  The solution of the system doesn't match the expected one.");
            }
        }
        solver.restoreSolutionRawReadPtr(solution);
    }

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
