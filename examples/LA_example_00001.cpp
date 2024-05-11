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

/*!
 * \example LA_example_00001.cpp
 *
 * \brief This example shows how to dump a SystemSolver.
 *
 * This example builds a small linear system, it dumps it to the disk,
 * it solves the linear system and it exports the solution to disk.
 * Destination of dump files is ./data/ folder and they are:
 * example_initialized_linear_system_info.dat
 * example_initialized_linear_system_A.dat.info
 * example_initialized_linear_system_A.dat
 * example_initialized_linear_system_solution.dat.info
 * example_initialized_linear_system_solution.dat
 * example_initialized_linear_system_rhs.dat.info
 * example_initialized_linear_system_rhs.dat
 * example_linear_system_solution.dat.info
 * example_linear_system_solution.dat
 *
 * <b>To run</b>: ./LA_example_00001.cpp \n
 */

/*!
 * Run the example
 *
 * Example of SystemSolver dump functionality.
 *
 * \param rank is the rank of the process
 * \param nProcs is the number of processes
 */
int run(int rank, int nProcs)
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
        rowValues[0]  = 1. / (double) (rowOffset + i + 1);

        matrix.addRow(rowPattern, rowValues);
    }

    matrix.assembly();

    // Build solver
    log::cout() << "Building solver..." << std::endl;

    SystemSolver solver;
    solver.assembly(matrix);

    double *rhs = solver.getRHSRawPtr();
    for (int i = 0; i < nCols; ++i) {
        rhs[i] = 1.;
    }
    solver.restoreRHSRawPtr(rhs);

    double *initialSolution = solver.getSolutionRawPtr();
    for (int i = 0; i < nRows; ++i) {
        initialSolution[i] = 0;
    }
    solver.restoreSolutionRawPtr(initialSolution);

    // Dump system
    log::cout() << "Dumping initialized linear system..." << std::endl;

    solver.dumpSystem("bitpit linear system", "data/", "example_initialized_linear_system_");

    // Set KSP
    KSPOptions &options = solver.getKSPOptions();
    options.atol        = 1e-50;
    options.rtol        = 1e-7;
    options.maxits      = 100;

    // Solve System
    log::cout() << "Solve linear system..." << std::endl;
    solver.solve();

    // Export solution
    log::cout() << "Export solution..." << std::endl;
    std::string solutionPath = std::string("./data/") + std::string("example_linear_system_solution.dat");
    solver.exportSolution(solutionPath, bitpit::SystemSolver::FILE_BINARY);

    return 0;
}

/*!
 * Main program.
 */
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc, &argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
    nProcs = 1;
#endif
    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Build, dump (initialization), solve linear system and dump the solution" << std::endl;

    int status;
    try {
        status = run(rank, nProcs);
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
