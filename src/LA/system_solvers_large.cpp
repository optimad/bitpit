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

#include <stdexcept>
#include <string>
#include <unordered_set>

#include "bitpit_common.hpp"

#include "system_solvers_large.hpp"

namespace bitpit {

/*!
 * \class SystemMatrixAssembler
 * \ingroup system_solver_large
 *
 * \brief The SystemMatrixAssembler class provides an interface for defining
 * system matrix assemblers.
 */

/*!
 * \class SystemSparseMatrixAssembler
 * \ingroup system_solver_large
 *
 * \brief The SystemSparseMatrixAssembler class defines an assembler for
 * building the system matrix form a sparse matrix.
 */

/*!
 * Constructor.
 *
 * \param matrix is the matrix
 */
SystemSparseMatrixAssembler::SystemSparseMatrixAssembler(const SparseMatrix *matrix)
    : SystemMatrixAssembler(), m_matrix(matrix)
{
}

/*!
 * Get the number of rows of the matrix.
 *
 * \result The number of rows of the matrix.
 */
long SystemSparseMatrixAssembler::getRowCount() const
{
    return m_matrix->getRowCount();
}

/*!
 * Get the number of columns of the matrix.
 *
 * \result The number of columns of the matrix.
 */
long SystemSparseMatrixAssembler::getColCount() const
{
    return m_matrix->getColCount();
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Get the global number of rows of the matrix.
 *
 * \result The global number of rows of the matrix.
 */
long SystemSparseMatrixAssembler::getRowGlobalCount() const
{
    return m_matrix->getRowGlobalCount();
}

/*!
 * Get the global number of columns of the matrix.
 *
 * \result The global number of columns of the matrix.
 */
long SystemSparseMatrixAssembler::getColGlobalCount() const
{
    return m_matrix->getColGlobalCount();
}

/*!
 * Get global row offset.
 *
 * \result The global row offset.
 */
long SystemSparseMatrixAssembler::getRowGlobalOffset() const
{
    return m_matrix->getRowGlobalOffset();
}

/*!
 * Get global column offset.
 *
 * \result The global column offset.
 */
long SystemSparseMatrixAssembler::getColGlobalOffset() const
{
    return m_matrix->getColGlobalOffset();
}
#endif

/*!
 * Get the number of non-zero elements in the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The number of non-zero elements in the specified row.
 */
long SystemSparseMatrixAssembler::getRowNZCount(long rowIndex) const
{
    return m_matrix->getRowNZCount(rowIndex);
}

/**
 * Get the maximum number of non-zero elements per row.
 *
 * \result The maximum number of non-zero elements per row.
 */
long SystemSparseMatrixAssembler::getMaxRowNZCount() const
{
    return m_matrix->getMaxRowNZCount();
}

/*!
 * Get the pattern of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the values of the specified row
 */
void SystemSparseMatrixAssembler::getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const
{
    m_matrix->getRowPattern(rowIndex, pattern);
}

/*!
 * Get the values of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param values on output will contain the values of the specified row
 */
void SystemSparseMatrixAssembler::getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    m_matrix->getRowValues(rowIndex, values);
}

/*!
 * Get the data of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the values of the specified row
 * \param values on output will contain the values of the specified row
 */
void SystemSparseMatrixAssembler::getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const
{
    m_matrix->getRowPattern(rowIndex, pattern);
    m_matrix->getRowValues(rowIndex, values);
}

/*!
 * \class PetscManager
 * \ingroup system_solver_large
 *
 * \brief The PetscManager class handles the interaction with PETSc library.
 */

/*!
 * Display the log view.
 *
 * \returns The error code.
 */
PetscErrorCode PetscManager::displayLogView()
{
    return PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
}

/*!
 * Add an initialization option.
 *
 * \param option is the option that will be added
 */
void PetscManager::addInitOption(const std::string &option)
{
    if (!areOptionsEditable()) {
        throw std::runtime_error("Initialization options can be set only before initializing the solver.");
    }

    m_options.push_back(option);
}

/*!
 * Add initialization options.
 *
 * \param argc is a non-negative value representing the number of arguments
 * passed to the program from the environment in which the program is run
 * \param argv is a pointer to the first element of an array of argc + 1
 * pointers, of which the last one is null and the previous ones, if any,
 * point to null-terminated multibyte strings that represent the arguments
 * passed to the program from the execution environment. If argv[0] is not
 * a null pointer (or, equivalently, if argc > 0), it points to a string
 * that represents the name used to invoke the program, or to an empty string.
 * The value of argv[0] is not propagated to the system solver, the solver
 * will see a dummy name.
 */
void PetscManager::addInitOptions(int argc, char **argv)
{
    if (!areOptionsEditable()) {
        throw std::runtime_error("Initialization options can be set only before initializing the solver.");
    }

    for (int i = 1; i < argc; ++i) {
        m_options.push_back(argv[i]);
    }
}

/*!
 * Add initialization options.
 *
 * \param options are the options that will be added
 */
void PetscManager::addInitOptions(const std::vector<std::string> &options)
{
    if (!areOptionsEditable()) {
        throw std::runtime_error("Initialization options can be set only before initializing the solver.");
    }

    for (const std::string &option : options) {
        m_options.push_back(option);
    }
}

/*!
 * Clear initialization options
 */
void PetscManager::clearInitOptions()
{
    m_options.clear();
}

/*!
 * Enable visualization of the log view.
 *
 * \returns The error code.
 */
void PetscManager::enableLogView()
{
    if (m_logViewEnabled) {
        return;
    }

#if (PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 7)
    PetscLogDefaultBegin();
#else
    PetscLogBegin();
#endif
    PetscRegisterFinalize(&(PetscManager::displayLogView));
    m_logViewEnabled = true;
}

/*!
 * Constructor.
 */
PetscManager::PetscManager()
    : m_externalMPIInitialization(true), m_externalPETScInitialization(true),
      m_options(std::vector<std::string>(1, "bitpit")),  m_logViewEnabled(false)
{
}

/*!
 * Destructor.
 */
PetscManager::~PetscManager()
{
    finalize(true);
}

/*!
 * Check if initialization options are editable.
 *
 * \result Returns true if initialization options are editable, false otherwise.
 */
bool PetscManager::areOptionsEditable() const
{
    PetscBool isPetscInitialized;
    PetscInitialized(&isPetscInitialized);

    return (!isPetscInitialized);
}

/*!
 * PETSc initialization.
 *
 * If PETSc is already initialized, this function is a no-op.
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * If the MPI initialization has not be called yet, it will be called by
 * this function.
 */
#endif
/*!
 * The function can be called more than once.
 *
 * \param debug if set to true, turns on logging of objects and events, once
 * the logging is enabled it cannot be disabled
 * \result Return true is PETSc hase been initialized by this function, false
 * if PETSc was already initialized.
 */
bool PetscManager::initialize(bool debug)
{
    // Early return if PETSc is already initialized
    PetscBool isPetscInitialized;
    PetscInitialized(&isPetscInitialized);
    if (isPetscInitialized) {
        return false;
    }

#if BITPIT_ENABLE_MPI==1
    // Early return if MPI is already finalized
    int isMPIFinalized;
    MPI_Finalized(&isMPIFinalized);
    if (isMPIFinalized) {
        throw std::runtime_error("PETSc finalization cannot be called after MPI finaliation.");
    }
#endif

    // Generate command line arguments
    //
    // The first argument is the executable name and it is set to a
    // dummy value.
    std::string help        = "None";
    std::string programName = "bitpit_petsc_manager";

    int argc = 1 + m_options.size();
    char **argv = new char*[argc + 1];
    argv[0] = strdup(programName.data());
    for (std::size_t i = 0; i < m_options.size(); i++) {
        argv[1 + i] = strdup(m_options[i].data());
    }
    argv[argc] = nullptr;

#if BITPIT_ENABLE_MPI==1
    // Initialize MPI
    int isMPIInitialized;
    MPI_Initialized(&isMPIInitialized);
    if (!isMPIInitialized) {
        MPI_Init(&argc, &argv);
        m_externalMPIInitialization = false;
    }
#endif

    // Initialize PETSc
    PetscInitialize(&argc, &argv, 0, help.data());
    m_externalPETScInitialization = false;

    // Enable log view
    if (debug) {
        enableLogView();
    }

    // Clean-up command line arguments
    for (int i = 0; i < argc; ++i) {
        free(argv[i]);
    }

    delete[] argv;

    // Initialization completed
    return true;
}

/*!
 * PETSc finalization.
 *
 * If PETSc finalization was already called, it will not be called again.
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * If MPI initialization has been performed by the PETSc manager and a
 * permanent finalization is requested, this function will call MPI
 * finalization.
 */
#else
/*!
 * If a permanent finalization is requested, it may not be possible to
 * re-initialize PETSc after the finalization has been performed.
 */
#endif
/*!
 * \param permanent if set to true, the function will try to finalized
 * both PETSc and its related libraries (e.g., MPI). If a permanent
 * finalization is requested, it may not be possible to re-initialize
 * PETSc after the finalization has been performed
 * \result Return true is PETSc has been finalized by this function, false
 * if PETSc was already finalized.
 */
bool PetscManager::finalize(bool permanent)
{
#if BITPIT_ENABLE_MPI==1
    // Early return if MPI is already finalized
    int isMPIFinalized;
    MPI_Finalized(&isMPIFinalized);
    if (isMPIFinalized) {
        return false;
    }
#endif

    // Early return if PETSc was initialized externally
    if (m_externalPETScInitialization) {
        return false;
    }

    // Finalize PETSc
    PetscBool isPetscFinalized;
    PetscFinalized(&isPetscFinalized);
    if (!isPetscFinalized) {
        PetscFinalize();
    }

#if BITPIT_ENABLE_MPI==1
    // Finalize MPI
    if (permanent) {
        if (!m_externalMPIInitialization) {
            MPI_Finalize();
        }
    }
#endif

    return (!isPetscFinalized);
}

/*!
 * \class SystemSolver
 * \ingroup system_solver_large
 *
 * \brief The SystemSolver class provides methods for building and solving
 * large linear systems.
 */

PetscManager SystemSolver::m_petscManager = PetscManager();

int SystemSolver::m_nInstances = 0;

/*!
 * Constructor.
 *
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool debug)
    : SystemSolver("", false, debug)
{
}

/*!
 * Constuctor
 *
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool transpose, bool debug)
    : SystemSolver("", transpose, debug)
{
}

/*!
 * Constructor.
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(const std::string &prefix, bool debug)
    : SystemSolver(prefix, false, debug)
{
}

/*!
 * Constuctor
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(const std::string &prefix, bool transpose, bool debug)
    : m_transpose(transpose),
      m_A(nullptr), m_rhs(nullptr), m_solution(nullptr),
      m_KSP(nullptr),
      m_prefix(prefix), m_assembled(false), m_setUp(false),
#if BITPIT_ENABLE_MPI==1
      m_communicator(MPI_COMM_SELF), m_partitioned(false),
#endif
      m_rowPermutation(nullptr), m_colPermutation(nullptr),
      m_forceConsistency(false)
{
    // Initialize PETSc
    if (m_nInstances == 0) {
        m_petscManager.initialize(debug);
    }

    // Set KSP debug options
    if (debug) {
#if (PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 7)
            PetscOptionsSetValue(nullptr, ("-" + m_prefix + "ksp_monitor_true_residual").c_str(), "");
            PetscOptionsSetValue(nullptr, ("-" + m_prefix + "ksp_converged_reason").c_str(), "");
            PetscOptionsSetValue(nullptr, ("-" + m_prefix + "ksp_monitor_singular_value").c_str(), "");
#else
            PetscOptionsSetValue(("-" + m_prefix + "ksp_monitor_true_residual").c_str(), "");
            PetscOptionsSetValue(("-" + m_prefix + "ksp_converged_reason").c_str(), "");
            PetscOptionsSetValue(("-" + m_prefix + "ksp_monitor_singular_value").c_str(), "");
#endif
    }

    // Increase the number of instances
    ++m_nInstances;
}

/*!
 * Destructor.
 */
SystemSolver::~SystemSolver()
{
    // Clear the solver
    clear();

    // Reset the permutations
    resetPermutations();

    // Decrease the number of instances
    --m_nInstances;

    // Finalize PETSc
    if (m_nInstances == 0) {
        m_petscManager.finalize(false);
    }
}

/*!
 * Clear the system
 */
void SystemSolver::clear()
{
    if (isSetUp()) {
        KSPDestroy(&m_KSP);
        m_KSP = nullptr;

        m_setUp = false;
    }

    if (isAssembled()) {
        MatDestroy(&m_A);
        VecDestroy(&m_rhs);
        VecDestroy(&m_solution);

#if BITPIT_ENABLE_MPI==1
        freeCommunicator();
#endif

        m_assembled = false;
    }
}

/*!
 * Set the permutations that will use internally by the solver.
 *
 * Only local permutations are supported.
 *
 * \param nRows are the rows of the matrix
 * \param rowRanks are the rank of the rows
 * \param nCols are the columns of the matrix
 * \param colRanks are the rank of the columns
 */
void SystemSolver::setPermutations(long nRows, const long *rowRanks, long nCols, const long *colRanks)
{
    // Permutation has to be set before assembling the system
    if (isAssembled()) {
        throw std::runtime_error("Unable to set the permutations. The system is already assembled.");
    }

    // Reset existing permutations
    resetPermutations();

    // Create new permutations
    PetscInt *rowPermutationsStorage;
    PetscMalloc(nRows * sizeof(PetscInt), &rowPermutationsStorage);
    for (long i = 0; i < nRows; ++i) {
        rowPermutationsStorage[i] = rowRanks[i];
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(m_communicator, nRows, rowPermutationsStorage, PETSC_OWN_POINTER, &m_rowPermutation);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nRows, rowPermutationsStorage, PETSC_OWN_POINTER, &m_rowPermutation);
#endif
    ISSetPermutation(m_rowPermutation);

    PetscInt *colPermutationsStorage;
    PetscMalloc(nCols * sizeof(PetscInt), &colPermutationsStorage);
    for (long j = 0; j < nCols; ++j) {
        colPermutationsStorage[j] = colRanks[j];
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(m_communicator, nCols, colPermutationsStorage, PETSC_OWN_POINTER, &m_colPermutation);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nCols, colPermutationsStorage, PETSC_OWN_POINTER, &m_colPermutation);
#endif
    ISSetPermutation(m_colPermutation);
}

/*!
 * Reset the permutations
 */
void SystemSolver::resetPermutations()
{
    if (m_rowPermutation) {
        ISDestroy(&m_rowPermutation);
    }

    if (m_colPermutation) {
        ISDestroy(&m_colPermutation);
    }
}

/*!
 * Assembly the system.
 *
 * \param matrix is the matrix
 */
void SystemSolver::assembly(const SparseMatrix &matrix)
{
    // Check if the matrix is assembled
    if (!matrix.isAssembled()) {
        throw std::runtime_error("Unable to assembly the system. The matrix is not yet assembled.");
    }

    // Assembly the system matrix
    SystemSparseMatrixAssembler assembler(&matrix);
#if BITPIT_ENABLE_MPI == 1
    assembly(matrix.getCommunicator(), matrix.isPartitioned(), assembler);
#else
    assembly(assembler);
#endif
}

#if BITPIT_ENABLE_MPI == 1
/*!
 * Assembly the system.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(const SystemMatrixAssembler &assembler)
{
    assembly(MPI_COMM_SELF, false, assembler);
}

/*!
 * Assembly the system.
 *
 * \param communicator is the MPI communicator
 * \param isPartitioned controls if the system is partitioned
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(MPI_Comm communicator, bool isPartitioned, const SystemMatrixAssembler &assembler)
{
#else
/*!
 * Assembly the system.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(const SystemMatrixAssembler &assembler)
{
#endif
    // Clear the system
    clear();

#if BITPIT_ENABLE_MPI == 1
    // Set the communicator
    setCommunicator(communicator);

    // Detect if the system is partitioned
    m_partitioned = isPartitioned;
#endif

    // Initialize matrix
    matrixCreate(assembler);
    matrixFill(assembler);

    // Initialize RHS and solution vectors
    vectorsCreate();

    // The system is now assembled
    m_assembled = true;
}

/*!
 * Update all the rows of the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param elements are the elements that will be used to update the rows
 */
void SystemSolver::update(const SparseMatrix &elements)
{
    update(getRowCount(), nullptr, elements);
}

/*!
 * Update the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the indices of the rows that will be updated
 * \param elements are the elements that will be used to update the rows
 */
void SystemSolver::update(long nRows, const long *rows, const SparseMatrix &elements)
{
    // Check if the element storage is assembled
    if (!elements.isAssembled()) {
        throw std::runtime_error("Unable to update the system. The element storage is not yet assembled.");
    }

    // Update matrix
    SystemSparseMatrixAssembler assembler(&elements);
    update(nRows, rows, assembler);
}

/*!
 * Update all the rows of the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::update(const SystemMatrixAssembler &assembler)
{
    update(getRowCount(), nullptr, assembler);
}

/*!
 * Update the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the indices of the rows that will be updated
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::update(long nRows, const long *rows, const SystemMatrixAssembler &assembler)
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to update the system. The system is not yet assembled.");
    }

    // Update matrix
    matrixUpdate(nRows, rows, assembler);
}

/**
* Get the number of rows of the system.
*
* \result The number of rows of the system.
*/
long SystemSolver::getRowCount() const
{
    if (!isAssembled()) {
        return 0;
    }

    PetscInt nRows;
    MatGetLocalSize(m_A, &nRows, NULL);

    return nRows;
}

/**
* Get the number of columns of the system.
*
* \result The number of columns of the system.
*/
long SystemSolver::getColCount() const
{
    if (!isAssembled()) {
        return 0;
    }

    PetscInt nCols;
    MatGetLocalSize(m_A, NULL, &nCols);

    return nCols;
}

#if BITPIT_ENABLE_MPI==1
/**
* Get the number of global rows
*
* \result The number of global rows
*/
long SystemSolver::getRowGlobalCount() const
{
    if (!isAssembled()) {
        return 0;
    }

    PetscInt nRows;
    MatGetSize(m_A, &nRows, NULL);

    return nRows;
}

/**
* Get number of global columns.
*
* \result The number of global columns.
*/
long SystemSolver::getColGlobalCount() const
{
    if (!isAssembled()) {
        return 0;
    }

    PetscInt nCols;
    MatGetSize(m_A, NULL, &nCols);

    return nCols;
}

/*!
    Checks if the matrix is partitioned.

    \result Returns true if the patch is partitioned, false otherwise.
*/
bool SystemSolver::isPartitioned() const
{
    return m_partitioned;
}
#endif

/*!
 * Check if the system is assembled.
 *
 * \return Returns true if the system is assembled, false otherwise.
 */
bool SystemSolver::isAssembled() const
{
    return m_assembled;
}

/*!
 * Check if the system is set up.
 *
 * \return Returns true if the system is set up, false otherwise.
 */
bool SystemSolver::isSetUp() const
{
    return m_setUp;
}

/*!
 * Solve the system
 */
void SystemSolver::solve()
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The system is not yet assembled.");
    }

    // Check if the system is set up
    if (!isSetUp()) {
        setUp();
    }

    // Perfrom actions before KSP solution
    preKSPSolveActions();

    // Force consistency
    if (m_forceConsistency) {
        removeNullSpaceFromRHS();
    }

    // Solve the system
    if (!m_transpose) {
        m_KSPStatus.error = KSPSolve(m_KSP, m_rhs, m_solution);
    } else {
        m_KSPStatus.error = KSPSolveTranspose(m_KSP, m_rhs, m_solution);
    }

    // Set solver info
    if (m_KSPStatus.error == 0) {
        KSPGetIterationNumber(m_KSP, &m_KSPStatus.its);
        KSPGetConvergedReason(m_KSP, &m_KSPStatus.convergence);
    } else {
        m_KSPStatus.its         = -1;
        m_KSPStatus.convergence = KSP_DIVERGED_BREAKDOWN;
    }

    // Perfrom actions after KSP solution
    postKSPSolveActions();
}

/*!
 * Solve the system
 *
 * \param rhs is the right-hand-side of the system
 * \param solution in input should contain the initial solution, on output it
 * contains the solution of the linear system
 */
void SystemSolver::solve(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Fills the vectors
    vectorsFill(rhs, solution);

    // Solve the system
    solve();

    // Export the solution
    vectorsExport(solution);
}

/*!
 * Pre-solve actions.
 */
void SystemSolver::preKSPSolveActions()
{
    // Apply permutations
    vectorsPermute(false);
}

/*!
 * Post-solve actions.
 */
void SystemSolver::postKSPSolveActions()
{
    // Invert permutations
    vectorsPermute(true);
}

/*!
 * Create the matrix.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::matrixCreate(const SystemMatrixAssembler &assembler)
{
    const PetscInt *rowRanks = nullptr;
    if (m_rowPermutation) {
        ISGetIndices(m_rowPermutation, &rowRanks);
    }

    // Set sizes
    long nRows = assembler.getRowCount();
    long nCols = assembler.getColCount();

#if BITPIT_ENABLE_MPI == 1
    long nGlobalRows = assembler.getRowGlobalCount();
    long nGlobalCols = assembler.getColGlobalCount();
#endif

    // Preallocation information
    std::vector<int> d_nnz(nRows, 0);
#if BITPIT_ENABLE_MPI == 1
    std::vector<int> o_nnz(nRows, 0);

    long firstColGlobalId = assembler.getColGlobalOffset();
    long lastColGlobalId  = firstColGlobalId + nCols - 1;

    ConstProxyVector<long> rowPattern;
#endif

    for (long row = 0; row < nRows; ++row) {
        long matrixRow = row;
        if (m_rowPermutation) {
            matrixRow = rowRanks[matrixRow];
        }

        d_nnz[row] = assembler.getRowNZCount(matrixRow);
#if BITPIT_ENABLE_MPI == 1
        if (m_partitioned) {
            assembler.getRowPattern(matrixRow, &rowPattern);

            int nRowNZ = rowPattern.size();
            for (int k = 0; k < nRowNZ; ++k) {
                long columnGlobalId = rowPattern[k];
                if (columnGlobalId < firstColGlobalId || columnGlobalId > lastColGlobalId) {
                    ++o_nnz[row];
                }
            }

            d_nnz[row] -= o_nnz[row];
        }
#endif
    }

    // Create the matrix
#if BITPIT_ENABLE_MPI == 1
    MatCreateAIJ(m_communicator, nRows, nCols, nGlobalRows, nGlobalCols, 0, d_nnz.data(), 0, o_nnz.data(), &m_A);
#else
    MatCreateSeqAIJ(PETSC_COMM_SELF, nRows, nCols, 0, d_nnz.data(), &m_A);
#endif

    // Cleanup
    if (m_rowPermutation) {
        ISRestoreIndices(m_rowPermutation, &rowRanks);
    }
}

/*!
 * Fills the matrix.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::matrixFill(const SystemMatrixAssembler &assembler)
{
    const long nRows = assembler.getRowCount();
    const long nCols = assembler.getColCount();
    const long maxRowNZ = assembler.getMaxRowNZCount();

    const PetscInt *rowRanks = nullptr;
    if (m_rowPermutation) {
        ISGetIndices(m_rowPermutation, &rowRanks);
    }

    IS invColPermutation;
    const PetscInt *colInvRanks = nullptr;
    if (m_colPermutation) {
        ISInvertPermutation(m_colPermutation, nCols, &invColPermutation);
        ISGetIndices(invColPermutation, &colInvRanks);
    }

    // Create the matrix
    if (maxRowNZ > 0) {
        std::vector<PetscInt> rowNZGlobalIds(maxRowNZ);
        std::vector<PetscScalar> rowNZValues(maxRowNZ);

        PetscInt rowGlobalOffset;
        MatGetOwnershipRangeColumn(m_A, &rowGlobalOffset, nullptr);

        PetscInt colGlobalBegin;
        PetscInt colGlobalEnd;
        MatGetOwnershipRangeColumn(m_A, &colGlobalBegin, &colGlobalEnd);

        ConstProxyVector<long> rowPattern;
        ConstProxyVector<double> rowValues;
        for (long row = 0; row < nRows; ++row) {
            long matrixRow = row;
            if (m_rowPermutation) {
                matrixRow = rowRanks[matrixRow];
            }

            assembler.getRowData(matrixRow, &rowPattern, &rowValues);

            const int nRowNZ = rowPattern.size();
            const PetscInt globalRow = rowGlobalOffset + row;
            for (int k = 0; k < nRowNZ; ++k) {
                long matrixGlobalCol = rowPattern[k];

                long globalCol = matrixGlobalCol;
                if (m_colPermutation) {
                    if (globalCol >= colGlobalBegin && globalCol < colGlobalEnd) {
                        long col = globalCol - colGlobalBegin;
                        col = colInvRanks[col];
                        globalCol = colGlobalBegin + col;
                    }
                }

                rowNZGlobalIds[k] = globalCol;
                rowNZValues[k]    = rowValues[k];
            }

            MatSetValues(m_A, 1, &globalRow, nRowNZ, rowNZGlobalIds.data(), rowNZValues.data(), INSERT_VALUES);
        }
    }

    // Let petsc build the matrix
    MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);

    // No new allocations are now allowed
    //
    // When updating the matrix it will not be possible to alter the pattern,
    // it will be possible to change only the values.
    MatSetOption(m_A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

    // Cleanup
    if (m_rowPermutation) {
        ISRestoreIndices(m_rowPermutation, &rowRanks);
    }

    if (m_colPermutation) {
        ISDestroy(&invColPermutation);
    }
}

/*!
 * Update the specified rows of the matrix.
 *
 * The contents of the specified rows will be replaced by the specified
 * elements.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the indices of the rows that will be updated, if a
 * null pointer is passed, the rows that will be updated are the rows
 * from 0 to (nRows - 1).
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::matrixUpdate(long nRows, const long *rows, const SystemMatrixAssembler &assembler)
{
    // Update element values
    PetscInt rowGlobalOffset;
    MatGetOwnershipRange(m_A, &rowGlobalOffset, nullptr);

    const long maxRowNZ = std::max(assembler.getMaxRowNZCount(), 0L);

    std::vector<PetscInt> rawRowPattern(maxRowNZ);
    std::vector<PetscScalar> rawRowValues(maxRowNZ);

    ConstProxyVector<long> rowPattern;
    ConstProxyVector<double> rowValues;
    for (long n = 0; n < nRows; ++n) {
        assembler.getRowData(n, &rowPattern, &rowValues);
        const int nRowElements = rowPattern.size();
        if (nRowElements == 0) {
            continue;
        }

        // Get global row
        long row;
        if (rows) {
            row = rows[n];
        } else {
            row = n;
        }

        const PetscInt globalRow = rowGlobalOffset + row;

        // Update values
        for (int k = 0; k < nRowElements; ++k) {
            rawRowPattern[k] = rowPattern[k];;
            rawRowValues[k]  = rowValues[k];
        }

        MatSetValues(m_A, 1, &globalRow, nRowElements, rawRowPattern.data(), rawRowValues.data(), INSERT_VALUES);
    }

    // Let petsc assembly the matrix after the update
    MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
}

/*!
 * Create RHS and solution vectors.
 */
void SystemSolver::vectorsCreate()
{
    PetscInt nRows;
    PetscInt nColumns;
    MatGetLocalSize(m_A, &nRows, &nColumns);

    PetscInt rhsSize;
    PetscInt solutionSize;
    if (!m_transpose) {
        rhsSize      = nRows;
        solutionSize = nColumns;
    } else {
        rhsSize      = nColumns;
        solutionSize = nRows;
    }

#if BITPIT_ENABLE_MPI == 1
    PetscInt nGlobalRows;
    PetscInt nGlobalColumns;
    MatGetSize(m_A, &nGlobalRows, &nGlobalColumns);

    PetscInt rhsGlobalSize;
    PetscInt solutionGlobalSize;
    if (!m_transpose) {
        rhsGlobalSize      = nGlobalRows;
        solutionGlobalSize = nGlobalColumns;
    } else {
        rhsGlobalSize      = nGlobalColumns;
        solutionGlobalSize = nGlobalRows;
    }

    PetscInt nGhosts;
    const PetscInt *ghosts;
    MatGetGhosts(m_A, &nGhosts, &ghosts);

    VecCreateGhost(m_communicator, solutionSize, solutionGlobalSize, nGhosts, ghosts, &m_solution);
    VecCreateGhost(m_communicator, rhsSize, rhsGlobalSize, nGhosts, ghosts, &m_rhs);
#else
    VecCreateSeq(PETSC_COMM_SELF, solutionSize, &m_solution);
    VecCreateSeq(PETSC_COMM_SELF, rhsSize, &m_rhs);
#endif
}

/*!
 * Apply permutations to RHS and solution vectors.
 *
 * \param invert is a flag for inverting the permutation
 */
void SystemSolver::vectorsPermute(bool invert)
{
    PetscBool petscInvert;
    if (invert) {
        petscInvert = PETSC_TRUE;
    } else {
        petscInvert = PETSC_FALSE;
    }

    if (m_colPermutation) {
        VecPermute(m_solution, m_colPermutation, petscInvert);
    }

    if (m_rowPermutation) {
        VecPermute(m_rhs, m_rowPermutation, petscInvert);
    }
}

/*!
 * Fills rhs and solution vectors.
 *
 * \param rhs is the right-hand-side of the system
 * \param solution is the solution of the linear system
 */
void SystemSolver::vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Import RHS
    int rhsSize;
    VecGetLocalSize(m_rhs, &rhsSize);

    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);
    for (int i = 0; i < rhsSize; ++i) {
        raw_rhs[i] = rhs[i];
    }
    VecRestoreArray(m_rhs, &raw_rhs);

    // Import initial solution
    int solutionSize;
    VecGetLocalSize(m_solution, &solutionSize);

    PetscScalar *raw_solution;
    VecGetArray(m_solution, &raw_solution);
    for (int i = 0; i < solutionSize; ++i) {
        raw_solution[i] = (*solution)[i];
    }
    VecRestoreArray(m_solution, &raw_solution);
}

/*!
 * Export the solution vector.
 *
 * \param solution on output it will contain the solution of the linear system
 */
void SystemSolver::vectorsExport(std::vector<double> *solution)
{
    int size;
    VecGetLocalSize(m_solution, &size);

    const PetscScalar *raw_solution;
    VecGetArrayRead(m_solution, &raw_solution);
    for (int i = 0; i < size; ++i) {
        (*solution)[i] = raw_solution[i];
    }
    VecRestoreArrayRead(m_solution, &raw_solution);
}

/*!
 * Get a raw pointer to the solution vector.
 *
 * \result A raw pointer to the solution vector.
 */
double * SystemSolver::getRHSRawPtr()
{
    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);

    return raw_rhs;
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getRHSRawPtr() const
{
    return getRHSRawReadPtr();
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getRHSRawReadPtr() const
{
    const PetscScalar *raw_rhs;
    VecGetArrayRead(m_rhs, &raw_rhs);

    return raw_rhs;
}

/*!
 * Restores the solution vector after getRHSRawPtr() has been called.
 *
 * \param raw_rhs is the location of pointer to array obtained from
 * getRHSRawPtr()
 */
void SystemSolver::restoreRHSRawPtr(double *raw_rhs)
{
    VecRestoreArray(m_rhs, &raw_rhs);
}

/*!
 * Restores the solution vector after getRHSRawReadPtr() has been called.
 *
 * \param raw_rhs is the location of pointer to array obtained from
 * getRHSRawReadPtr()
 */
void SystemSolver::restoreRHSRawReadPtr(const double *raw_rhs) const
{
    VecRestoreArrayRead(m_rhs, &raw_rhs);
}

/*!
 * Get a raw pointer to the solution vector.
 *
 * \result A raw pointer to the solution vector.
 */
double * SystemSolver::getSolutionRawPtr()
{
    PetscScalar *raw_solution;
    VecGetArray(m_solution, &raw_solution);

    return raw_solution;
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getSolutionRawPtr() const
{
    return getSolutionRawReadPtr();
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getSolutionRawReadPtr() const
{
    const PetscScalar *raw_solution;
    VecGetArrayRead(m_solution, &raw_solution);

    return raw_solution;
}

/*!
 * Restores the solution vector after getSolutionRawPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getSolutionRawPtr()
 */
void SystemSolver::restoreSolutionRawPtr(double *raw_solution)
{
    VecRestoreArray(m_solution, &raw_solution);
}

/*!
 * Restores the solution vector after getSolutionRawReadPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getSolutionRawReadPtr()
 */
void SystemSolver::restoreSolutionRawReadPtr(const double *raw_solution) const
{
    VecRestoreArrayRead(m_solution, &raw_solution);
}

/*!
 * Dump the system to file
 *
 * \param directory is the directory where the files will be saved
 * \param prefix is the prefix that will be added to the files
 * \param matrixFormat is the dump format that will be used for the matrix,
 * the ASCII format may not be able to dump large matrices
 * \param rhsFormat is the dump format that will be used for the RHS,
 * the ASCII format may not be able to dump large vectors
 * \param solutionFormat is the dump format that will be used for the solution,
 * the ASCII format may not be able to dump large vectors
 */
void SystemSolver::dump(const std::string &directory, const std::string &prefix,
                        DumpFormat matrixFormat, DumpFormat rhsFormat,
                        DumpFormat solutionFormat) const
{
    std::stringstream filePathStream;

    // Matrix
    PetscViewerType matrixViewerType;
    PetscViewerFormat matrixViewerFormat;
    if (matrixFormat == DUMP_BINARY) {
        matrixViewerType   = PETSCVIEWERBINARY;
        matrixViewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        matrixViewerType   = PETSCVIEWERASCII;
        matrixViewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer matViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &matViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &matViewer);
#endif
    PetscViewerSetType(matViewer, matrixViewerType);
    PetscViewerFileSetMode(matViewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(matViewer, matrixViewerFormat);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "A.txt";
    PetscViewerFileSetName(matViewer, filePathStream.str().c_str());
    MatView(m_A, matViewer);
    PetscViewerDestroy(&matViewer);

    // RHS
    PetscViewerType rhsViewerType;
    PetscViewerFormat rhsViewerFormat;
    if (rhsFormat == DUMP_BINARY) {
        rhsViewerType   = PETSCVIEWERBINARY;
        rhsViewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        rhsViewerType   = PETSCVIEWERASCII;
        rhsViewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer rhsViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &rhsViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &rhsViewer);
#endif
    PetscViewerSetType(rhsViewer, rhsViewerType);
    PetscViewerFileSetMode(rhsViewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(rhsViewer, rhsViewerFormat);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "rhs.txt";
    PetscViewerFileSetName(rhsViewer, filePathStream.str().c_str());
    VecView(m_rhs, rhsViewer);
    PetscViewerDestroy(&rhsViewer);

    // Solution
    PetscViewerType solutionViewerType;
    PetscViewerFormat solutionViewerFormat;
    if (solutionFormat == DUMP_BINARY) {
        solutionViewerType   = PETSCVIEWERBINARY;
        solutionViewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        solutionViewerType   = PETSCVIEWERASCII;
        solutionViewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer solutionViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &solutionViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &solutionViewer);
#endif
    PetscViewerSetType(solutionViewer, solutionViewerType);
    PetscViewerFileSetMode(solutionViewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(solutionViewer, solutionViewerFormat);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "solution.txt";
    PetscViewerFileSetName(solutionViewer, filePathStream.str().c_str());
    VecView(m_solution, solutionViewer);
    PetscViewerDestroy(&solutionViewer);
}

/*!
 * Attaches a null space to the system matrix.
 */
void SystemSolver::setNullSpace()
{
    MatNullSpace nullspace;
#if BITPIT_ENABLE_MPI==1
    MatNullSpaceCreate(m_communicator, PETSC_TRUE, 0, NULL, &nullspace);
#else
    MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_TRUE, 0, NULL, &nullspace);
#endif
    MatSetNullSpace(m_A, nullspace);
    MatNullSpaceDestroy(&nullspace);
}

/*!
 * Removes the null space from the system matrix.
 */
void SystemSolver::unsetNullSpace()
{
    MatSetNullSpace(m_A, NULL);
}

/*!
 * Setup the system.
 */
void SystemSolver::setUp()
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The system is not yet assembled.");
    }

    // Destroy existing Krylov space
    if (m_KSP) {
        KSPDestroy(&m_KSP);
    }

    // Create Krylov space
#if BITPIT_ENABLE_MPI==1
    KSPCreate(m_communicator, &m_KSP);
#else
    KSPCreate(PETSC_COMM_SELF, &m_KSP);
#endif

    // Set options prefix
    if (!m_prefix.empty()) {
        KSPSetOptionsPrefix(m_KSP, m_prefix.c_str());
    }

    // Set the matrix associated with the linear system
    KSPSetOperators(m_KSP, m_A, m_A);

    // Perform actions before KSP set up
    preKSPSetupActions();

    // Setup Krylov space
    KSPSetFromOptions(m_KSP);
    KSPSetUp(m_KSP);

    // Perform actions after KSP set up
    postKSPSetupActions();

    // Set up is now complete
    m_setUp = true;
}

/*!
 * Perform actions before KSP setup.
 */
void SystemSolver::preKSPSetupActions()
{
    // Preconditioner configuration
    PCType preconditionerType;
#if BITPIT_ENABLE_MPI == 1
    if (isPartitioned()) {
        preconditionerType = PCASM;
    } else {
        preconditionerType = PCILU;
    }
#else
    preconditionerType = PCILU;
#endif

    PC preconditioner;
    KSPGetPC(m_KSP, &preconditioner);
    PCSetType(preconditioner, preconditionerType);
    if (strcmp(preconditionerType, PCASM) == 0) {
        if (m_KSPOptions.overlap != PETSC_DEFAULT) {
            PCASMSetOverlap(preconditioner, m_KSPOptions.overlap);
        }
    } else if (strcmp(preconditionerType, PCILU) == 0) {
        if (m_KSPOptions.levels != PETSC_DEFAULT) {
            PCFactorSetLevels(preconditioner, m_KSPOptions.levels);
        }
    }

    // Solver configuration
    KSPSetType(m_KSP, KSPFGMRES);
    if (m_KSPOptions.restart != PETSC_DEFAULT) {
        KSPGMRESSetRestart(m_KSP, m_KSPOptions.restart);
    }
    if (m_KSPOptions.rtol != PETSC_DEFAULT || m_KSPOptions.atol != PETSC_DEFAULT || m_KSPOptions.maxits != PETSC_DEFAULT) {
        KSPSetTolerances(m_KSP, m_KSPOptions.rtol, m_KSPOptions.atol, PETSC_DEFAULT, m_KSPOptions.maxits);
    }
    KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);
}

/*!
 * Perform actions after KSP setup.
 */
void SystemSolver::postKSPSetupActions()
{
    // Get preconditioner information
    PC preconditioner;
    KSPGetPC(m_KSP, &preconditioner);

    PCType preconditionerType;
    PCGetType(preconditioner, &preconditionerType);

    // Set ASM sub block preconditioners
    if (strcmp(preconditionerType, PCASM) == 0) {
        KSP *subksp;
        PC subpc;
        PetscInt nlocal, first;
        PCASMGetSubKSP(preconditioner, &nlocal, &first, &subksp);
        for (PetscInt i = 0; i < nlocal; ++i) {
            KSPGetPC(subksp[i], &subpc);
            PCSetType(subpc, PCILU);
            if (m_KSPOptions.sublevels != PETSC_DEFAULT) {
                PCFactorSetLevels(subpc, m_KSPOptions.sublevels);
            }
            if (m_KSPOptions.subrtol != PETSC_DEFAULT) {
                KSPSetTolerances(subksp[i], m_KSPOptions.subrtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            }
        }
    }
}

/*!
 * Get a reference to the options associated to the Kryolov solver.
 *
 * \return A reference to the options associated to the Kryolov solver.
 */
KSPOptions & SystemSolver::getKSPOptions()
{
    return m_KSPOptions;
}

/*!
 * Get a constant reference to the options associated to the Kryolov solver.
 *
 * \return A constant reference to the options associated to the Kryolov solver.
 */
const KSPOptions & SystemSolver::getKSPOptions() const
{
    return m_KSPOptions;
}

/*!
 * Get a constant reference to the status of the Kryolov solver.
 *
 * \return A constant reference to the status of the Kryolov solver.
 */
const KSPStatus & SystemSolver::getKSPStatus() const
{
    return m_KSPStatus;
}

#if BITPIT_ENABLE_MPI==1
/*!
	Gets the MPI communicator associated to the system.

	\return The MPI communicator associated to the system.
*/
const MPI_Comm & SystemSolver::getCommunicator() const
{
	return m_communicator;
}

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void SystemSolver::setCommunicator(MPI_Comm communicator)
{
    if ((communicator != MPI_COMM_NULL) && (communicator != MPI_COMM_SELF)) {
        MPI_Comm_dup(communicator, &m_communicator);
    } else {
        m_communicator = MPI_COMM_SELF;
    }
}

/*!
	Frees the MPI communicator associated to the matrix.
*/
void SystemSolver::freeCommunicator()
{
    if (m_communicator != MPI_COMM_SELF) {
        int finalizedCalled;
        MPI_Finalized(&finalizedCalled);
        if (!finalizedCalled) {
            MPI_Comm_free(&m_communicator);
        }
    }
}
#endif

/*!
    Remove null space components from right hand side.

    The purpose of this function is to force right hand side consistency.
    This is a hack, your model should produce consistent right hand sides.
*/
void SystemSolver::removeNullSpaceFromRHS()
{
    // Get the null space from the matrix
    MatNullSpace nullspace;
    MatGetNullSpace(m_A, &nullspace);

    // Remove null space components from right hand side
    MatNullSpaceRemove(nullspace, m_rhs);
}

/*!
    Check if right hand side consistency is forced before solving the system.

    In order to force right hand side consistency, null space components
    are removed from the right hand side before solving the system.
*/
bool SystemSolver::isForceConsistencyEnabled() const
{
    return m_forceConsistency;
}

/*!
    Enable or disable forcing right hand side consistency.

    \param enable if set to true, right hand side consistency will be
    forced before solving the system
*/
void SystemSolver::enableForceConsistency(bool enable)
{
    m_forceConsistency = enable;
}

}
