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
#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_LA.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing parallel solution of transposed linear systems.
*
* \param rank is the rank of the process
* \param nProcs is the number of processes
*/
int subtest_001(int rank, int nProcs)
{
    int nRows;
    int nCols;
    int nNZ;
    if (nProcs == 1) {
        nRows = 10;
        nCols = 10;
        nNZ   = 20;
    } else if (rank <= 1) {
        nRows = 5;
        nCols = 5;
        nNZ   = 10;
    } else {
        nRows = 0;
        nCols = 0;
        nNZ   = 0;
    }

    // Build matrix
    log::cout() << "Building matrix..." << std::endl;

    SparseMatrix matrix(MPI_COMM_WORLD, true, nRows, nCols, nNZ);

    std::vector<long> rowPattern(2);
    std::vector<double> rowValues(2);

    int rowOffset  = matrix.getRowGlobalOffset();
    int nGlobalRos = matrix.getRowGlobalCount();
    for (int row = 0; row < nRows; ++row) {
        int globalRow = rowOffset + row;

        rowPattern[0] = globalRow;
        rowValues[0]  = (globalRow + 1);

        rowPattern[1] = nGlobalRos - globalRow - 1;
        rowValues[1]  = 11 * (globalRow + 1);

        matrix.addRow(rowPattern, rowValues);
    }
    matrix.assembly();

    // Build system
    log::cout() << "Building system..." << std::endl;

    std::vector<double> globalRHS(10);
    globalRHS[0] = 1101;
    globalRHS[1] =  895;
    globalRHS[2] =  713;
    globalRHS[3] =  555;
    globalRHS[4] =  421;
    globalRHS[5] =  311;
    globalRHS[6] =  225;
    globalRHS[7] =  163;
    globalRHS[8] =  125;
    globalRHS[9] =  111;

    SystemSolver system(true, false);
    system.assembly(matrix);

    double *rhs = system.getRHSRawPtr();
    for (int row = 0; row < nRows; ++row) {
        int globalRow = rowOffset + row;

        rhs[row] = globalRHS[globalRow];
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

    if (nRows > 0) {
        const double *solution = system.getSolutionRawReadPtr();
        for (int i = 0; i < nRows; ++i) {
            log::cout() << "  Solution[" << i << "] = " << solution[i] << std::endl;

            double expectedSolution = matrix.getRowGlobalOffset() + i + 1;
            if (!utils::DoubleFloatingEqual()(solution[i], expectedSolution, 1e-10)) {
                log::cout() << "  Expected solution[" << i << "] = " << expectedSolution << std::endl;
                log::cout() << "  Error[" << i << "] = " << (expectedSolution - solution[i]) << std::endl;
                throw std::runtime_error("  The solution of the system doesn't match the expected one.");
            }
        }
        system.restoreSolutionRawReadPtr(solution);
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
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setDefaultVisibility(log::GLOBAL);

	// Run the subtests
    log::cout() << "Testing parallel solution of transposed linear systems." << std::endl;

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
