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
 * \example LA_example_00002.cpp
 *
 * \brief This example shows how to restore a SystemSolver.
 *
 * This example restores and solves an initialized linear system.
 * It compares the computed solution with a restored solution,
 * using the same tolerance used by the linear system solver.
 * Input data are in ./data folder and they are:
 *    example_initialized_linear_system_info.txt
 *    example_initialized_linear_system_A.txt.info
 *    example_initialized_linear_system_A.txt
 *    example_initialized_linear_system_solution.txt.info
 *    example_initialized_linear_system_solution.txt
 *    example_initialized_linear_system_rhs.txt.info
 *    example_initialized_linear_system_rhs.txt
 *    example_linear_system_solution.txt.info
 *    example_linear_system_solution.txt
 *
 * The user can rename their files using the standard names above
 * to run a custom case.
 *
 * This example contains an extension of the bare SystemSolver class
 * adding the chance to solve using algebraic multigrid (PETSc GAMG)
 * preconditioner. Only non-block matrices can be solved using GAMG.
 * A check disable algebraic multigrid for block matrices with a
 * message in the logger.
 *
 * <b>To run</b>: ./LA_example_00002.cpp \n
*/

/*!
 * \class MyLinearSolver
 *
 * \brief The class MyLinearSolver is derived form SystemSolver and it
 * allows to solve the linear system using PETSc GAMG preconditioner.
 * Only non-block matrices can be solved using GAMG.
 */
class MyLinearSolver : public SystemSolver {

public:
    /*!
     * Constructor.
     *
     * \param transpose if true the solver solves the transposed linear system.
     * \param multigrid if true the solver uses algebraic multigrid as preconditioner.
     * \param debug if true PETSc monitors on residual, reason and singular values are turned on.
     */
    MyLinearSolver(bool transpose, bool multigrid, bool debug)
        : SystemSolver(transpose, debug),
          m_multigrid(multigrid)
    {
    }

    /*!
     * Set the algebraic multigrid preconditioner.
     *
     * \param multigrid if true the algebraic multigrid preconditioner is turned on.
     */
    void setMultigridFlag(bool multigrid)
    {
        m_multigrid = multigrid;
    }

protected:
    const static PetscInt GAMG_MAX_LEVELS;
    const static PetscReal GAMG_THRESHOLD;
    const static PetscReal GAMG_THRESHOLD_SCALE;

    bool m_multigrid;

    /*!
     * Perform actions before KSP setup.
     */
    void preKSPSetupActions() override
    {
        //
        // Preconditioner configuration
        //
        PCType preconditionerType;
        if (m_multigrid) {
            preconditionerType = PCGAMG;
        } else {
#if BITPIT_ENABLE_MPI == 1
            if (this->isPartitioned()) {
                preconditionerType = PCASM;
            } else {
                preconditionerType = PCILU;
            }
#else
            preconditionerType = PCILU;
#endif
        }

        PC preconditioner;
        KSPGetPC(this->m_KSP, &preconditioner);
        PCSetType(preconditioner, preconditionerType);

        // Set preconditioner options
        if (strcmp(preconditionerType, PCASM) == 0) {
            if (this->m_KSPOptions.overlap != PETSC_DEFAULT) {
                PCASMSetOverlap(preconditioner, this->m_KSPOptions.overlap);
            }
        } else if (strcmp(preconditionerType, PCILU) == 0) {
            if (this->m_KSPOptions.levels != PETSC_DEFAULT) {
                PCFactorSetLevels(preconditioner, this->m_KSPOptions.levels);
            }
        } else if (strcmp(preconditionerType, PCGAMG) == 0) {
            std::vector<PetscReal> thresholds(GAMG_MAX_LEVELS);
            for (int i = 0; i < GAMG_MAX_LEVELS; ++i) {
                thresholds[i] = GAMG_THRESHOLD;
            }

            PCGAMGSetNSmooths(preconditioner, 0);
#if (PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 17)
#if BITPIT_ENABLE_MPI == 1
            if (this->isPartitioned()) {
                PCGAMGSetSymGraph(preconditioner, PETSC_TRUE);
            }
#endif
#endif

            PCGAMGSetNlevels(preconditioner, GAMG_MAX_LEVELS);
#if (PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 16)
            PCGAMGSetThreshold(preconditioner, thresholds.data(), GAMG_MAX_LEVELS);
            PCGAMGSetThresholdScale(preconditioner, GAMG_THRESHOLD_SCALE);
#else
            PCGAMGSetThreshold(preconditioner, GAMG_THRESHOLD);
#endif
        }

        //
        // Solver configuration
        //
        KSPSetType(this->m_KSP, KSPFGMRES);
        if (this->m_KSPOptions.restart != PETSC_DEFAULT) {
            KSPGMRESSetRestart(this->m_KSP, this->m_KSPOptions.restart);
        }
        if (this->m_KSPOptions.rtol != PETSC_DEFAULT || this->m_KSPOptions.atol != PETSC_DEFAULT || this->m_KSPOptions.maxits != PETSC_DEFAULT) {
            KSPSetTolerances(this->m_KSP, this->m_KSPOptions.rtol, this->m_KSPOptions.atol, PETSC_DEFAULT, this->m_KSPOptions.maxits);
        }
        KSPSetInitialGuessNonzero(this->m_KSP, PETSC_TRUE);
    }

    /*!
     * Perform actions after KSP setup.
     */
    void postKSPSetupActions() override
    {
        // Post-setup actions for the base class
        SystemSolver::postKSPSetupActions();

        // Post-setup actions for the preconditioner
        PC preconditioner;
        KSPGetPC(this->m_KSP, &preconditioner);

        PCType preconditionerType;
        PCGetType(preconditioner, &preconditionerType);
        if (strcmp(preconditionerType, PCGAMG) == 0) {
            KSP coarseKSP;
            PCMGGetCoarseSolve(preconditioner, &coarseKSP);

            PC coarsePreconditioner;
            KSPGetPC(coarseKSP, &coarsePreconditioner);
            PCSetType(coarsePreconditioner, PCSVD);
        }
    }
};

const PetscInt MyLinearSolver::GAMG_MAX_LEVELS       = 10;
const PetscReal MyLinearSolver::GAMG_THRESHOLD       = 0.02;
const PetscReal MyLinearSolver::GAMG_THRESHOLD_SCALE = 1.0;

/*!
 * Run the example.
 *
 * Example of SystemSolver restore functionality.
 *
 * \param rank is the rank of the process
 * \param nProcs is the number of processes
 */
int run(int rank, int nProcs)
{
    BITPIT_UNUSED(nProcs);

    // Restore solver with initial guess
    log::cout() << "Restoring linear system ..." << std::endl;
    std::string solverPath = "./data/";
    bool multigrid         = false;
    bool debug             = false;
    bool transpose         = false;
    MyLinearSolver restoredSolver(transpose, multigrid, debug);
#if BITPIT_ENABLE_MPI==1
    restoredSolver.restoreSystem(MPI_COMM_WORLD, solverPath, "example_initialized_linear_system_");
#else
    restoredSolver.restoreSystem(solverPath, "example_initialized_linear_system_");
#endif
    log::cout() << "Linear system restored." << std::endl;

    // Check if multigrid is possible. Multigrid is not possible on BAIJ matrices
    if (multigrid) {
        log::cout() << "Enable multigrid preconditioner..." << std::endl;
        int blockSize        = restoredSolver.getBlockSize();
        bool enableMultigrid = (blockSize == 1);
        restoredSolver.setMultigridFlag(enableMultigrid);
        if (enableMultigrid) {
            log::cout() << "Multigrid preconditioner enabled..." << std::endl;
        } else {
            log::cout() << "The matrix is a block matrix, no multigrid preconditioning is possible." << std::endl;
        }
    }

    // Set KSP options
    KSPOptions &options = restoredSolver.getKSPOptions();
    options.atol        = 1e-50;
    options.rtol        = 1e-7;
    options.maxits      = 100;

    // Solve linear system
    log::cout() << "Solving linear system..." << std::endl;
    restoredSolver.solve();
    const KSPStatus &status = restoredSolver.getKSPStatus();
    log::cout() << "Reason: " << status.convergence << " in its: " << status.its << std::endl;
    log::cout() << "Linear system solved." << std::endl;

    // Turn solution comparison on
    bool compareSolutions = true;
    if (compareSolutions) {
        // Store expected solution for comparison
        log::cout() << "Storing computed solution for comparison..." << std::endl;
        std::size_t vecLength = restoredSolver.getRowElementCount();
        std::vector<double> expectedSolution(vecLength, 0.0);
        double *computedSolution = restoredSolver.getSolutionRawPtr();
        for (std::size_t i = 0; i < vecLength; ++i) {
            expectedSolution[i] = computedSolution[i];
        }
        restoredSolver.restoreSolutionRawPtr(computedSolution);

        // Fill solution with solution from outside
        log::cout() << "Restoring solution..." << std::endl;
        std::string solutionPath = solverPath + "example_linear_system_solution.txt";
        restoredSolver.fillSolution(solutionPath);
        log::cout() << "Solution restored." << std::endl;

        // Set comparison tolerance to the linear system relative tolerance
        log::cout() << "Comparing solutions..." << std::endl;

        const double TOLERANCE = options.rtol;
        log::cout() << "Tolerance = " << options.rtol << std::endl;

        // Compare computed solution with restored one
        double *restoredSolution = restoredSolver.getSolutionRawPtr();
        for (std::size_t i = 0; i < vecLength; ++i) {
            // Check for absolute tolerance only
            if (!utils::DoubleFloatingEqual()(expectedSolution[i] - restoredSolution[i], 0.0, TOLERANCE, TOLERANCE)) {
                double error = expectedSolution[i] - restoredSolution[i];

                log::cout() << "    Solutions do not match (on rank #" << rank << " error for DoF # " << i << " is " << error << ")" << std::endl;
            }
        }
        restoredSolver.restoreSolutionRawPtr(restoredSolution);
        log::cout() << "Solutions compared..." << std::endl;
    }

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

    // Run the example
    log::cout() << "Linear system restore example" << std::endl;

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
