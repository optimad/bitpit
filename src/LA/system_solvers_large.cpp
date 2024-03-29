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

#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "bitpit_common.hpp"

#include "system_solvers_large.hpp"

namespace bitpit {

/*!
 * \class SystemMatrixOrdering
 * \ingroup system_solver_large
 *
 * \brief The SystemMatrixOrdering class provides an interface for defining
 * classes that allows to reorder the system matrix.
 *
 * If the system is partitioned, each process can reorder only it's local
 * part of the matix.
 */

/*!
 * \class NaturalSystemMatrixOrdering
 * \ingroup system_solver_large
 *
 * \brief The NaturalSystemMatrixOrdering class defines allows to use a matrix
 * natural ordering.
 */

/*!
 * Get the rank of the specified local row.
 *
 * \param row is the local row
 * \result The rank of the specified local row.
 */
long NaturalSystemMatrixOrdering::getRowRank(long row) const
{
    return row;
}

/*!
 * Get the rank of the specified local column.
 *
 * \param col is the local column
 * \result The rank of the specified local column.
 */
long NaturalSystemMatrixOrdering::getColRank(long col) const
{
    return col;
}

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
 * Get the assembly options.
 *
 * \result The assembly options that will be used.
 */
SystemMatrixAssembler::AssemblyOptions SystemSparseMatrixAssembler::getOptions() const
{
    AssemblyOptions options;
    options.full   = true;
    options.sorted = false;

    return options;
}

/*!
 * Get the transpose flag.
 *
 * \result Returns true if the transposed system will be solved, false otherwise.
 */
bool SystemSolver::getTranspose() const
{
    return m_transpose;
}

/*!
 * Set the transpose flag.
 *
 * If the system is already assembled and the transpose flag needs to be changed,
 * both the soltion vector and the RHS one will be destroyed and re-created.
 *
 * If the transpose flag needs to be changed, the workspaces associated wit the
 * system will be cleared.
 *
 * \param transpose if set to true, transposed system will be solved
 */
void SystemSolver::setTranspose(bool transpose)
{
    // Early return if the transpose flag is already correctly set.
    if (m_transpose == transpose) {
        return;
    }

    // Clear the workspace
    clearWorkspace();

    // Re-create the vectors
    if (isAssembled()) {
        VecDestroy(&m_rhs);
        VecDestroy(&m_solution);

        vectorsCreate();
    }

    // Set the flag
    m_transpose = transpose;
}

/*!
 * Get the block size.
 *
 * \result The block size.
 */
int SystemSparseMatrixAssembler::getBlockSize() const
{
    return m_matrix->getBlockSize();
}

/*!
 * Get the number of (block) rows handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block rows, where a block row is
 * defined as a group of blockSize matrix rows.
 *
 * \result The number of (block) rows handled by the assembler.
 */
long SystemSparseMatrixAssembler::getRowCount() const
{
    return m_matrix->getRowCount();
}

/*!
 * Get the number of (block) columns handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block columns, where a block column
 * is defined as a group of blockSize matrix columns.
 *
 * \result The number of (block) columns handled by the assembler.
 */
long SystemSparseMatrixAssembler::getColCount() const
{
    return m_matrix->getColCount();
}

/*!
 * Get the number of elements in the rows handled by the assembler.
 *
 * This function will return the effective number of rows of the matrix that
 * will be assembled.
 *
 * \result The number of rows handled by the assembler.
 */
long SystemSparseMatrixAssembler::getRowElementCount() const
{
    long nRowElements = getBlockSize() * getRowCount();

    return nRowElements;
}

/*!
 * Get the number of elements in the columns handled by the assembler.
 *
 * This function will return the effective number of columns of the matrix that
 * will be assembled.
 *
 * \result The number of columns handled by the assembler.
 */
long SystemSparseMatrixAssembler::getColElementCount() const
{
    long nColElements = getBlockSize() * getColCount();

    return nColElements;
}

#if BITPIT_ENABLE_MPI==1
/**
 * Get number of global (block) rows handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block rows, where a block row
 * is defined as a group of blockSize matrix rows.
 *
 * \result The number of global rows handled by the assembler.
 */
long SystemSparseMatrixAssembler::getRowGlobalCount() const
{
    return m_matrix->getRowGlobalCount();
}

/**
 * Get number of global (block) columns handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block columns, where a block
 * column is defined as a group of blockSize matrix columns.
 *
 * \result The number of global (block) columns handled by the assembler.
 */
long SystemSparseMatrixAssembler::getColGlobalCount() const
{
    return m_matrix->getColGlobalCount();
}

/**
 * Get the number of global elements in the rows handled by the assembler.
 *
 * This function will return the effective global number of rows of the system
 * matrix.
 *
 * \result The number of global elements in the rows handled by the assembler.
 */
long SystemSparseMatrixAssembler::getRowGlobalElementCount() const
{
    long nElements = getBlockSize() * getRowGlobalCount();

    return nElements;
}

/*!
 * Get the global number of columns handled by the assembler.
 *
 * This function will return the effective global number of columns of the
 * system matrix.
 *
 * \result The global number of columns handled by the assembler.
 */
long SystemSparseMatrixAssembler::getColGlobalElementCount() const
{
    long nElements = getBlockSize() * getColGlobalCount();

    return nElements;
}

/*!
 * Get global (block) row offset.
 *
 * \result The global (block) row offset.
 */
long SystemSparseMatrixAssembler::getRowGlobalOffset() const
{
    return m_matrix->getRowGlobalOffset();
}

/*!
 * Get global (block) column offset.
 *
 * \result The global (block) column offset.
 */
long SystemSparseMatrixAssembler::getColGlobalOffset() const
{
    return m_matrix->getColGlobalOffset();
}

/*!
 * Get global offset for the elements of the row.
 *
 * This function will return the offset expressed in effective rows of the
 * system matrix.
 *
 * \result The global offset for the elements of the row.
 */
long SystemSparseMatrixAssembler::getRowGlobalElementOffset() const
{
    long offset = getBlockSize() * getRowGlobalOffset();

    return offset;
}

/*!
 * Get global offset for the elements of the column.
 *
 * This function will return the offset expressed in effective columns of the
 * system matrix.
 *
 * \result The global offset for the elements of the column.
 */
long SystemSparseMatrixAssembler::getColGlobalElementOffset() const
{
    long offset = getBlockSize() * getColGlobalOffset();

    return offset;
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
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global ids of the block columns of the row,
 * where a block column is defined as a group of blockSize matrix columns.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the pattern of the specified (block) row
 */
void SystemSparseMatrixAssembler::getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const
{
    m_matrix->getRowPattern(rowIndex, pattern);
}

/*!
 * Get the values of the specified row.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the values of all the elements of a block row,
 * where a block column is defined as a group of blockSize matrix columns. The
 * values are returned as a row-oriented logically two-dimensional array of
 * values.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
 */
void SystemSparseMatrixAssembler::getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    m_matrix->getRowValues(rowIndex, values);
}

/*!
 * Get the data of the specified row.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * the pattern defines the global ids of the block columns of the row, where a
 * block column is defined as a group of blockSize matrix columns.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * the values contain of all the elements of a block row, where a block column
 * is defined as a group of blockSize matrix columns. The values are returned
 * as a row-oriented logically two-dimensional array of values.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the pattern of the specified (block) row
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
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
    finalize();
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
 * \result Return true is PETSc has been finalized by this function, false
 * if PETSc was already finalized.
 */
bool PetscManager::finalize()
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
    if (!m_externalMPIInitialization) {
        MPI_Finalize();
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
 *
 * Rather than working with individual elements in the system matrix, it is
 * possible to employ blocks of elements. The size of the blocks can be defined
 * during assembly. When a size different that one is provided, the matrix will
 * store elements by fixed-sized dense nb × nb blocks, where nb is the size of
 * the blocks. Blocking may be advantageous when solving PDE-based simulations
 * that leads to matrices with a naturally blocked structure (with a block size
 * equal to the number of degrees of freedom per cell).
 *
 * When blocking is used, row and column indexes will count the number of blocks
 * in the row/column direction, not the number of rows/columns of the matrix.
 */

PetscManager SystemSolver::m_petscManager = PetscManager();

int SystemSolver::m_nInstances = 0;

/*!
 * Constructor.
 *
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool debug)
    : SystemSolver("", false, false, debug)
{
}

/*!
 * Constuctor
 *
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool transpose, bool debug)
    : SystemSolver("", false, transpose, debug)
{
}


/*!
 * Constuctor
 *
 * \param flatten if set to true, the system matrix will be created with a
 * unitary block size, regardless of the blocks size of the assembler
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool flatten, bool transpose, bool debug)
    : SystemSolver("", flatten, transpose, debug)
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
 * Constructor.
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(const std::string &prefix, bool transpose, bool debug)
    : SystemSolver(prefix, false, transpose, debug)
{
}

/*!
 * Constuctor
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param flatten if set to true, the system matrix will be created with a
 * unitary block size, regardless of the blocks size of the assembler
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(const std::string &prefix, bool flatten, bool transpose, bool debug)
    : m_flatten(flatten), m_transpose(transpose),
      m_A(PETSC_NULL), m_rhs(PETSC_NULL), m_solution(PETSC_NULL),
      m_KSP(PETSC_NULL),
      m_prefix(prefix), m_assembled(false), m_KSPDirty(true),
#if BITPIT_ENABLE_MPI==1
      m_communicator(MPI_COMM_SELF), m_partitioned(false),
#endif
      m_rowReordering(PETSC_NULL), m_colReordering(PETSC_NULL),
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

    // Decrease the number of instances
    --m_nInstances;
}

/*!
 * Clear the system
 */
void SystemSolver::clear()
{
    clearWorkspace();

    vectorsDestroy();
    matrixDestroy();

    clearReordering();

#if BITPIT_ENABLE_MPI==1
    freeCommunicator();
#endif

    m_assembled = false;
}

/*!
 * Clear and release the memory of all data structures needed for the solution of the system
 *
 * These data structures will be re-created the next time the system will be solved.
 */
void SystemSolver::clearWorkspace()
{
    // Clear KSP
    if (m_KSP) {
        destroyKSP();
    }
}

/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size of the matrix received in
 * input.
 *
 * \param matrix is the matrix
 */
void SystemSolver::assembly(const SparseMatrix &matrix)
{
    assembly(matrix, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size of the matrix received in
 * input.
 *
 * \param matrix is the matrix
 * \param reordering is the reordering that will be applied when assemblying the
 * system
 */
void SystemSolver::assembly(const SparseMatrix &matrix, const SystemMatrixOrdering &reordering)
{
    // Check if the matrix is assembled
    if (!matrix.isAssembled()) {
        throw std::runtime_error("Unable to assembly the system. The matrix is not yet assembled.");
    }

    // Assembly the system matrix
    SystemSparseMatrixAssembler assembler(&matrix);
#if BITPIT_ENABLE_MPI == 1
    assembly(matrix.getCommunicator(), matrix.isPartitioned(), assembler, reordering);
#else
    assembly(assembler, reordering);
#endif
}

#if BITPIT_ENABLE_MPI == 1
/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size specified by the assembler.
 *
 * \param communicator is the MPI communicator
 * \param isPartitioned controls if the system is partitioned
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(MPI_Comm communicator, bool isPartitioned, const SystemMatrixAssembler &assembler)
{
    assembly(communicator, isPartitioned, assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size specified by the assembler.
 *
 * \param communicator is the MPI communicator
 * \param isPartitioned controls if the system is partitioned
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assemblying the
 * system
 */
void SystemSolver::assembly(MPI_Comm communicator, bool isPartitioned, const SystemMatrixAssembler &assembler, const SystemMatrixOrdering &reordering)
{
#else
/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size specified by the assembler.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(const SystemMatrixAssembler &assembler)
{
    assembly(assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * If the system was created with the flatten flag set to true, the system matrix
 * will be created with a unitary block size. Otherwise, the block size of the
 * system matrix will be set equal to the block size specified by the assembler.
 *
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assemblying the
 * system
 */
void SystemSolver::assembly(const SystemMatrixAssembler &assembler, const SystemMatrixOrdering &reordering)
{
#endif
    // Clear the system
    clear();

    // Set reordering
    setReordering(assembler.getRowCount(), assembler.getColCount(), reordering);

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

/*!
 * Get the block size of the system.
 *
 * \result The block size of the system.
 */
int SystemSolver::getBlockSize() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt blockSize;
    MatGetBlockSize(m_A, &blockSize);

    return blockSize;
}

/**
 * Get the number of rows of the system.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block rows, where a block row is
 * defined as a group of blockSize matrix rows.
 *
 * \result The number of rows of the system.
 */
long SystemSolver::getRowCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nRows;
    MatGetLocalSize(m_A, &nRows, NULL);
    nRows /= getBlockSize();

    return nRows;
}

/**
 * Get the number of columns of the system.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block columns, where a block column
 * is defined as a group of blockSize matrix columns.
 *
 * \result The number of columns of the system.
 */
long SystemSolver::getColCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nCols;
    MatGetLocalSize(m_A, NULL, &nCols);
    nCols /= getBlockSize();

    return nCols;
}

/**
 * Get the number of elements in the rows of the system.
 *
 * This function will return the effective number of rows of the system matrix.
 * This value is the same as the local size used in creating the y vector for
 * the matrix-vector product y = Ax.
 *
 * \result The number of elements in the rows of the system.
 */
long SystemSolver::getRowElementCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nRows;
    MatGetLocalSize(m_A, &nRows, NULL);

    return nRows;
}

/**
 * Get the number of elements in the columns of the system.
 *
 * This function will return the effective number of columns of the system
 * matrix. This value is the same as the local size used in creating the x
 * vector for the matrix-vector product y = Ax.
 *
 * \result The number of elements in the columns of the system.
 */
long SystemSolver::getColElementCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nCols;
    MatGetLocalSize(m_A, NULL, &nCols);

    return nCols;
}

#if BITPIT_ENABLE_MPI==1
/**
 * Get number of global (block) rows.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block rows, where a block row
 * is defined as a group of blockSize matrix rows.
 *
 * \result The number of global rows
 */
long SystemSolver::getRowGlobalCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nRows;
    MatGetSize(m_A, &nRows, NULL);
    nRows /= getBlockSize();

    return nRows;
}

/**
 * Get number of global (block) columns.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block columns, where a block
 * column is defined as a group of blockSize matrix columns.
 *
 * \result The number of global (block) columns.
 */
long SystemSolver::getColGlobalCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nCols;
    MatGetSize(m_A, NULL, &nCols);
    nCols /= getBlockSize();

    return nCols;
}

/**
 * Get the number of global elements in the rows of the system.
 *
 * This function will return the effective global number of rows of the system
 * matrix.
 *
 * \result The number of global elements in the rows of the system.
 */
long SystemSolver::getRowGlobalElementCount() const
{
    if (m_A == PETSC_NULL) {
        return 0;
    }

    PetscInt nRows;
    MatGetSize(m_A, &nRows, NULL);

    return nRows;
}

/**
 * Get the number of global elements in the columns of the system.
 *
 * This function will return the effective global number of columns of the
 * system matrix.
 *
 * \result The number of global elements in the columns of the system.
 */
long SystemSolver::getColGlobalElementCount() const
{
    if (m_A == PETSC_NULL) {
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
 * This function is deprecated, it will always return true.
 *
 * \return Returns true if the system is set up, false otherwise.
 */
bool SystemSolver::isSetUp() const
{
    return true;
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

    // Prepare the KSP
    prepareKSP();

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

    // Finalize the KSP
    finalizeKSP();
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
    // Reorder vectors
    vectorsReorder(true);
}

/*!
 * Post-solve actions.
 */
void SystemSolver::postKSPSolveActions()
{
    // Reorder vectors
    vectorsReorder(false);
}

/*!
 * Create the matrix.
 *
 * Matrix will be created, but it will not be initialized.
 *
 * \param blockSize is the block size of the matrix
 */
void SystemSolver::matrixCreate(int blockSize)
{
    // Create the matrix
#if BITPIT_ENABLE_MPI == 1
    MatCreate(m_communicator, &m_A);
#else
    MatCreate(PETSC_COMM_SELF, &m_A);
#endif

    // Set matrix type
#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
        if (blockSize > 1) {
            MatSetType(m_A, MATMPIBAIJ);
        } else {
            MatSetType(m_A, MATMPIAIJ);
        }
    } else
#endif
    {
        if (blockSize > 1) {
            MatSetType(m_A, MATSEQBAIJ);
        } else {
            MatSetType(m_A, MATSEQAIJ);
        }
    }

    // Set block size
    if (blockSize > 1) {
        MatSetBlockSize(m_A, blockSize);
    }
}

/*!
 * Create the matrix.
 *
 * If the system was created with the flatten flag set to true, the matrix will
 * be created with a unitary block size. Otherwise, the block size of the system
 * matrix will be set equal to the block size specified by the assembler.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::matrixCreate(const SystemMatrixAssembler &assembler)
{
    const PetscInt *rowReordering = nullptr;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    // Set block size
    //
    // When blocks are ignored, the PETSc matrix will be created with a unitary
    // block size, regardless of the blocks size of the assembler. If this is
    // the case, preallociation information provided by the assembler will be
    // properly expanded to have preallocation information for each matrix row.
    int assemblerBlockSize = assembler.getBlockSize();

    int matrixBlockSize;
    if (m_flatten) {
        matrixBlockSize = 1;
    } else {
        matrixBlockSize = assembler.getBlockSize();
    }

    int blockExpansionSize = (matrixBlockSize != assemblerBlockSize) ? assemblerBlockSize : 1;

    // Create the matrix
    matrixCreate(matrixBlockSize);

    // Get sizes
    long nAssemblerRows = assembler.getRowCount();

    long nRowsElements = assembler.getRowElementCount();
    long nColsElements = assembler.getColElementCount();

    long nGlobalRowsElements;
    long nGlobalColsElements;
#if BITPIT_ENABLE_MPI == 1
    nGlobalRowsElements = assembler.getRowGlobalElementCount();
    nGlobalColsElements = assembler.getColGlobalElementCount();
#else
    nGlobalRowsElements = nRowsElements;
    nGlobalColsElements = nColsElements;
#endif

    MatSetSizes(m_A, nRowsElements, nColsElements, nGlobalRowsElements, nGlobalColsElements);

    // Preallocation information
    long nPreallocationRows = blockExpansionSize * nAssemblerRows;

    std::vector<int> d_nnz(nPreallocationRows, 0);
    for (long n = 0; n < nAssemblerRows; ++n) {
        long matrixRow = n;
        if (rowReordering) {
            matrixRow = rowReordering[matrixRow];
        }

        int nAssemblerRowNZ = assembler.getRowNZCount(n);

        long matrixRowOffset = matrixRow * blockExpansionSize;
        for (int i = 0; i < blockExpansionSize; ++i) {
            d_nnz[matrixRowOffset + i] = blockExpansionSize * nAssemblerRowNZ;
        }
    }


#if BITPIT_ENABLE_MPI == 1
    std::vector<int> o_nnz(nPreallocationRows, 0);
    if (m_partitioned) {
        long nAssemblerCols = assembler.getColCount();

        long assemblerDiagonalBegin = assembler.getColGlobalOffset();
        long assemblerDiagonalEnd   = assemblerDiagonalBegin + nAssemblerCols;

        ConstProxyVector<long> assemblerRowPattern(static_cast<std::size_t>(0), assembler.getMaxRowNZCount());
        for (long n = 0; n < nAssemblerRows; ++n) {
            long matrixRow = n;
            if (rowReordering) {
                matrixRow = rowReordering[matrixRow];
            }


            assembler.getRowPattern(n, &assemblerRowPattern);
            int nAssemblerRowNZ = assemblerRowPattern.size();

            long matrixRowOffset = matrixRow * blockExpansionSize;
            for (int k = 0; k < nAssemblerRowNZ; ++k) {
                long id = assemblerRowPattern[k];
                if (id < assemblerDiagonalBegin || id >= assemblerDiagonalEnd) {
                    for (int i = 0; i < blockExpansionSize; ++i) {
                        o_nnz[matrixRowOffset + i] += blockExpansionSize;
                        d_nnz[matrixRowOffset + i] -= blockExpansionSize;
                    }
                }
            }
        }
    }
#endif

#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
        if (matrixBlockSize == 1) {
            MatMPIAIJSetPreallocation(m_A, 0, d_nnz.data(), 0, o_nnz.data());
        } else {
            MatMPIBAIJSetPreallocation(m_A, matrixBlockSize, 0, d_nnz.data(), 0, o_nnz.data());
        }
    } else
#endif
    {
        if (matrixBlockSize == 1) {
            MatSeqAIJSetPreallocation(m_A, 0, d_nnz.data());
        } else {
            MatSeqBAIJSetPreallocation(m_A, matrixBlockSize, 0, d_nnz.data());
        }

    }

    // Each process will only set values for its own rows
    MatSetOption(m_A, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);

#if PETSC_VERSION_GE(3, 12, 0)
    // The first assembly will set a superset of the off-process entries
    // required for all subsequent assemblies. This avoids a rendezvous
    // step in the MatAssembly functions.
    MatSetOption(m_A, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE);
#endif

    // Cleanup
    if (m_rowReordering) {
        ISRestoreIndices(m_rowReordering, &rowReordering);
    }
}

/*!
 * Fills the matrix reading its contents from the specified assembler..
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::matrixFill(const SystemMatrixAssembler &assembler)
{
    // Fil matrix
    matrixUpdate(assembler.getRowCount(), nullptr, assembler);

    // No new allocations are now allowed
    //
    // When updating the matrix it will not be possible to alter the pattern,
    // it will be possible to change only the values.
    MatSetOption(m_A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
}

/*!
 * Fill the matrix reading its contents form the specified file.
 *
 * The input file should contain a compatible matrix stored in PETSc binary format. It's up to the
 * caller of this routine to make sure the loaded matrix is compatible with the system. If the
 * matrix file cannot be read an exception is thrown.
 *
 * \param filePath is the path of the file
 */
void SystemSolver::matrixFill(const std::string &filePath)
{
    // Check if the matrix exists
    if (!m_A) {
        throw std::runtime_error("Matrix should be created before restoring it.");
    }

    // Check if the file exists
    std::ifstream fileStream(filePath.c_str());
    if (!fileStream.good()) {
        throw std::runtime_error("The PETSc matrix file \"" + filePath + "\" doesn't exists");
    }
    fileStream.close();

    // Clear workspace
    clearWorkspace();

    // Clear reordering
    clearReordering();

    // Restore the matrix
    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, filePath.c_str());
    MatLoad(m_A, viewer);
    PetscViewerDestroy(&viewer);

    // Re-create vectors
    vectorsDestroy();
    vectorsCreate();
}

/*!
 * Destroy the matrix.
 */
void SystemSolver::matrixDestroy()
{
    if (m_A) {
        MatDestroy(&m_A);
        m_A = PETSC_NULL;
    }
}

/*!
 * Update the specified rows of the matrix.
 *
 * The contents of the specified rows will be replaced by the specified
 * elements.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the local indices of the rows that will be updated, if a
 * null pointer is passed, the rows that will be updated are the rows
 * from 0 to (nRows - 1).
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::matrixUpdate(long nRows, const long *rows, const SystemMatrixAssembler &assembler)
{
    // Updating the matrix invalidates the KSP
    m_KSPDirty = true;

    // Get block size
    const int matrixBlockSize    = getBlockSize();
    const int assemblerBlockSize = assembler.getBlockSize();
    if (assemblerBlockSize != matrixBlockSize && matrixBlockSize != 1) {
        std::string message = "Unable to update the matrix.";
        message += " The block size of the assembler is not compatible with the block size of the system matrix.";
        throw std::runtime_error(message);
    }

    // Initialize reordering
    const PetscInt *rowReordering = nullptr;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    const PetscInt *colReordering = nullptr;
    if (m_colReordering) {
        ISGetIndices(m_colReordering, &colReordering);
    }

    // Global information
    PetscInt colGlobalBegin;
    PetscInt colGlobalEnd;
    MatGetOwnershipRangeColumn(m_A, &colGlobalBegin, &colGlobalEnd);
    colGlobalBegin /= assemblerBlockSize;
    colGlobalEnd /= assemblerBlockSize;

    PetscInt rowGlobalOffset;
    MatGetOwnershipRange(m_A, &rowGlobalOffset, nullptr);
    rowGlobalOffset /= assemblerBlockSize;

    // Get the options for assembling the matrix
    SystemMatrixAssembler::AssemblyOptions assemblyOptions = assembler.getOptions();

    PetscBool matrixSortedFull = (assemblyOptions.full && assemblyOptions.sorted) ? PETSC_TRUE : PETSC_FALSE;

#if PETSC_VERSION_GE(3, 12, 0)
    // Check if it is possible to speedup insertion of values
    //
    // The option MAT_SORTED_FULL means that each process provides exactly its
    // local rows; all column indices for a given row are passed in a single call
    // to MatSetValues(), preallocation is perfect, row oriented, INSERT_VALUES
    // is used. If this options is set to PETSC_TRUE, the function MatSetValues
    // will be faster.
    //
    // This options needs at least PETSc 3.12.
    MatSetOption(m_A, MAT_SORTED_FULL, matrixSortedFull);
#endif

    // Update element values
    //
    // If the sizes of PETSc data types match the sizes of data types expected by
    // bitpit a direct update can be performed, otherwise the matrix is updated
    // using intermediate data storages.
    const long maxRowNZ = std::max(assembler.getMaxRowNZCount(), 0L);
    bool valuesDirectUpdate = (sizeof(double) == sizeof(PetscScalar));

    ConstProxyVector<long> rowPattern(static_cast<std::size_t>(0), static_cast<std::size_t>(maxRowNZ));
    std::vector<PetscInt> petscRowPattern(maxRowNZ);
    std::vector<PetscInt> petscMatrixRowPattern(assemblerBlockSize * maxRowNZ);

    ConstProxyVector<double> rowValues;
    std::vector<PetscScalar> petscRowValuesStorage;
    const PetscScalar *petscRowValues;
    if (!valuesDirectUpdate) {
        long maxRowNZElements = assemblerBlockSize * assemblerBlockSize * maxRowNZ;
        rowValues.set(ConstProxyVector<double>::INTERNAL_STORAGE, 0, maxRowNZElements);
        petscRowValuesStorage.resize(maxRowNZElements);
        petscRowValues = petscRowValuesStorage.data();
    }

    for (long n = 0; n < nRows; ++n) {
        // Get row information
        long row;
        if (rows) {
            row = rows[n];
        } else {
            row = n;
        }

        if (rowReordering) {
            row = rowReordering[row];
        }

        const PetscInt globalRow = rowGlobalOffset + row;

        // Check if it possible to perform a fast update
        //
        // A fast update allows to set all the values of a row at once (without
        // the need to get the row pattern), it can be performed if:
        //  - the system matrix has already been assembled;
        //  - the system matrix has a unitary block size;
        //  - the assembler is providing all the values of the row;
        //  - values provided by the assembler are sorted by ascending column.
        //
        // If fast update is used, row values will be set using a special PETSc
        // function (MatSetValuesRow) that allows to set all the values of a
        // row at once, without requiring the pattern of the row.
        //
        // Fast update is not related to the option MAT_SORTED_FULL, that option
        // is used to speedup the standard function MatSetValues (which still
        // requires the pattern of the row).
        bool fastUpdate = isAssembled() && (matrixBlockSize == 1) && (matrixSortedFull == PETSC_TRUE);

        // Get row data
        if (fastUpdate) {
            assembler.getRowValues(n, &rowValues);
        } else {
            assembler.getRowData(n, &rowPattern, &rowValues);
        }

        // Get values in PETSc format
        if (valuesDirectUpdate) {
            petscRowValues = reinterpret_cast<const PetscScalar *>(rowValues.data());
        } else {
            std::copy(rowValues.cbegin(), rowValues.cend(), petscRowValuesStorage.begin());
        }

        if (fastUpdate) {
            // Set values
            if (assemblerBlockSize == 1) {
                MatSetValuesRow(m_A, globalRow, petscRowValues);
            } else {
                const int nMatrixRowValues = rowValues.size() / assemblerBlockSize;

                PetscInt matrixGlobalRow = globalRow * assemblerBlockSize;
                const PetscScalar *petscMatrixRowValues = petscRowValues;
                for (int k = 0; k < assemblerBlockSize; ++k) {
                    MatSetValuesRow(m_A, matrixGlobalRow, petscMatrixRowValues);
                    ++matrixGlobalRow;
                    petscMatrixRowValues += nMatrixRowValues;
                }
            }
        } else {
            // Get the size of the pattern
            const int rowPatternSize = rowPattern.size();
            if (rowPatternSize == 0) {
                continue;
            }

            // Get pattern in PETSc format
            for (int k = 0; k < rowPatternSize; ++k) {
                long globalCol = rowPattern[k];
                if (colReordering) {
                    if (globalCol >= colGlobalBegin && globalCol < colGlobalEnd) {
                        long col = globalCol - colGlobalBegin;
                        col = colReordering[col];
                        globalCol = colGlobalBegin + col;
                    }
                }

                petscRowPattern[k] = globalCol;
            }

            // Set data
            if (matrixBlockSize > 1) {
                MatSetValuesBlocked(m_A, 1, &globalRow, rowPatternSize, petscRowPattern.data(), petscRowValues, INSERT_VALUES);
            } else {
                if (assemblerBlockSize == 1) {
                    MatSetValues(m_A, 1, &globalRow, rowPatternSize, petscRowPattern.data(), petscRowValues, INSERT_VALUES);
                } else {
                    // Expand the pattern to get a values for each element
                    const int matrixRowPatternSize = rowPatternSize * assemblerBlockSize;
                    for (int k = 0; k < rowPatternSize; ++k) {
                        int matrixPattenOffset = k * assemblerBlockSize;
                        for (int n = 0; n < assemblerBlockSize; ++n) {
                            petscMatrixRowPattern[matrixPattenOffset + n] = assemblerBlockSize * petscRowPattern[k] + n;
                        }
                    }

                    // Set data
                    PetscInt matrixGlobalRow = globalRow * assemblerBlockSize;
                    const PetscScalar *petscMatrixRowValues = petscRowValues;
                    for (int k = 0; k < assemblerBlockSize; ++k) {
                        MatSetValues(m_A, 1, &matrixGlobalRow, matrixRowPatternSize, petscMatrixRowPattern.data(), petscMatrixRowValues, INSERT_VALUES);
                        ++matrixGlobalRow;
                        petscMatrixRowValues += matrixRowPatternSize;
                    }
                }
            }
        }
    }

    // Let petsc assembly the matrix after the update
    MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);

    // Cleanup
    if (rowReordering) {
        ISRestoreIndices(m_rowReordering, &rowReordering);
    }

    if (colReordering) {
        ISRestoreIndices(m_colReordering, &colReordering);
    }
}

/*!
 * Create RHS and solution vectors.
 *
 * Vectors will be created, but they will not be initialized.
 */
void SystemSolver::vectorsCreate()
{
    if (!m_transpose) {
        MatCreateVecs(m_A, &m_rhs, &m_solution);
    } else {
        MatCreateVecs(m_A, &m_solution, &m_rhs);
    }
}

/*!
 * Destroy RHS and solution vectors.
 */
void SystemSolver::vectorsDestroy()
{
    if (m_rhs) {
        VecDestroy(&m_rhs);
        m_rhs = PETSC_NULL;
    }

    if (m_solution) {
        VecDestroy(&m_solution);
        m_solution = PETSC_NULL;
    }
}

/*!
 * Reorder RHS and solution vectors to match the order of the system matrix.
 *
 * \param invert is a flag for inverting the ordering
 */
void SystemSolver::vectorsReorder(bool invert)
{
    PetscBool petscInvert;
    if (invert) {
        petscInvert = PETSC_TRUE;
    } else {
        petscInvert = PETSC_FALSE;
    }

    if (m_colReordering) {
        VecPermute(m_rhs, m_colReordering, petscInvert);
    }

    if (m_rowReordering) {
        VecPermute(m_solution, m_rowReordering, petscInvert);
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
    // Fill RHS
    int rhsSize;
    VecGetLocalSize(m_rhs, &rhsSize);

    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);
    for (int i = 0; i < rhsSize; ++i) {
        raw_rhs[i] = rhs[i];
    }
    VecRestoreArray(m_rhs, &raw_rhs);

    // Fill initial solution
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
 * Fills rhs and solution vectors.
 *
 * \param rhs is the right-hand-side of the system
 * \param solution is the solution of the linear system
 */

/*!
 * Fill the matrix reading its contents form the specified file.
 *
 * The input file should contain a compatible matrix stored in PETSc binary format. It's up to the
 * caller of this routine to make sure the loaded matrix is compatible with the system. If the
 * matrix file cannot be read an exception is thrown.
 *
 * \param filePath is the path of the file
 */


/*!
 * Import the RHS vector from the specified file.
 *
 * The input files should contain compatible vectors stored in PETSc binary format. If the size
 * of the loaded vectors is not compatible with the matrix, an exception is thrown.
 *
 * If a file cannot be read, the corresponding vector will not be filled.
 *
 * \param rhsFilePath is the path of the file that contains the RHS vector
 * \param solutionFilePath is the path of the file that contains the solution vector
 */
void SystemSolver::vectorsFill(const std::string &rhsFilePath, const std::string &solutionFilePath)
{
    // Fill the RHS vector
    std::ifstream rhsFileStream(rhsFilePath);
    if (rhsFileStream.good()) {
        fillRHS(rhsFilePath);
    }
    rhsFileStream.close();

    // Fill the solution vector
    std::ifstream solutionFileStream(solutionFilePath);
    if (solutionFileStream.good()) {
        fillSolution(solutionFilePath);
    }
    solutionFileStream.close();
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
 * Export the matrix to file
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the matrix,
 * the ASCII format may not be able to dump large matrices
 */
void SystemSolver::exportMatrix(const std::string &filePath, FileFormat fileFormat) const
{
    PetscViewerType viewerType;
    PetscViewerFormat viewerFormat;
    if (fileFormat == FILE_BINARY) {
        viewerType   = PETSCVIEWERBINARY;
        viewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        viewerType   = PETSCVIEWERASCII;
        viewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer matViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &matViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &matViewer);
#endif
    PetscViewerSetType(matViewer, viewerType);
    PetscViewerFileSetMode(matViewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(matViewer, viewerFormat);

    PetscViewerFileSetName(matViewer, filePath.c_str());
    MatView(m_A, matViewer);
    PetscViewerDestroy(&matViewer);
}

/*!
 * Export the RHS vector to file
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the RHS,
 * the ASCII format may not be able to handle large vectors
 */
void SystemSolver::exportRHS(const std::string &filePath, FileFormat fileFormat) const
{
    PetscViewerType viewerType;
    PetscViewerFormat viewerFormat;
    if (fileFormat == FILE_BINARY) {
        viewerType   = PETSCVIEWERBINARY;
        viewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        viewerType   = PETSCVIEWERASCII;
        viewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, viewerType);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(viewer, viewerFormat);

    PetscViewerFileSetName(viewer, filePath.c_str());
    VecView(m_rhs, viewer);
    PetscViewerDestroy(&viewer);
}

/*!
 * Dump the solution vector to file
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the solution,
 * the ASCII format may not be able to handle large vectors
 */
void SystemSolver::exportSolution(const std::string &filePath, FileFormat fileFormat) const
{
    PetscViewerType viewerType;
    PetscViewerFormat viewerFormat;
    if (fileFormat == FILE_BINARY) {
        viewerType   = PETSCVIEWERBINARY;
        viewerFormat = PETSC_VIEWER_DEFAULT;
    } else {
        viewerType   = PETSCVIEWERASCII;
        viewerFormat = PETSC_VIEWER_ASCII_MATLAB;
    }

    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, viewerType);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(viewer, viewerFormat);

    PetscViewerFileSetName(viewer, filePath.c_str());
    VecView(m_solution, viewer);
    PetscViewerDestroy(&viewer);
}

/*!
 * Fill the RHS vector reading its contents from the specified file.
 *
 * The input file should contain a compatible vector stored in PETSc binary format. If the size
 * of the vector stored in the file is not compatible with the matrix, an exception is thrown.
 * An exception is also raised if the file cannot be read.
 *
 * It is possible to fill the RHS vector only after the system has been assembled.
 *
 * \param filePath is the path of the file
 */
void SystemSolver::fillRHS(const std::string &filePath)
{
    // Check if the system is assembled
    if (!m_assembled) {
        throw std::runtime_error("The RHS vector can be loaded only after assembling the system.");
    }

    // Check if the file exists
    std::ifstream fileStream(filePath);
    if (!fileStream.good()) {
        throw std::runtime_error("The file \"" + filePath + "\" cannot be read.");
    }
    fileStream.close();

    // Restore the RHS
    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, filePath.c_str());
    VecLoad(m_rhs, viewer);
    PetscViewerDestroy(&viewer);

    // Check if the imported RHS is compatible with the matrix
    PetscInt vecSize;
    VecGetLocalSize(m_rhs, &vecSize);

    PetscInt expectedVecSize;
    if (!m_transpose) {
        expectedVecSize = getRowCount();
    } else {
        expectedVecSize = getColCount();
    }
    expectedVecSize *= getBlockSize();

    if (vecSize != expectedVecSize) {
        log::cout() << "The imported RHS vector is not compatible with the matrix" << std::endl;
        log::cout() << "The size of the imported RHS vector is " << vecSize << std::endl;
        log::cout() << "The expected size of RHS vector is " << expectedVecSize << std::endl;
        throw std::runtime_error("The imported RHS vector is not compatible with the matrix");
    }
}

/*!
 * Fill the solution vector reading its contents from the specified file.
 *
 * The input file should contain a compatible vector stored in PETSc binary format. If the size
 * of the vector stored in the file is not compatible with the matrix, an exception is thrown.
 * An exception is also raised if the file cannot be read.
 *
 * It is possible to fill the solution vector only after the system has been assembled.
 *
 * \param filePath is the path of the file
 */
void SystemSolver::fillSolution(const std::string &filePath)
{
    // Check if the system is assembled
    if (!m_assembled) {
        throw std::runtime_error("The solution vector can be loaded only after assembling the system.");
    }

    // Check if the file exists
    std::ifstream fileStream(filePath);
    if (!fileStream.good()) {
        throw std::runtime_error("The file \"" + filePath + "\" cannot be read.");
    }
    fileStream.close();

    // Restore the vector
    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, filePath.c_str());
    VecLoad(m_solution, viewer);
    PetscViewerDestroy(&viewer);

    // Check if the imported solution is compatible with the matrix
    PetscInt vecSize;
    VecGetLocalSize(m_solution, &vecSize);

    PetscInt expectedVecSize;
    if (!m_transpose) {
        expectedVecSize = getColCount();
    } else {
        expectedVecSize = getRowCount();
    }
    expectedVecSize *= getBlockSize();

    if (vecSize != expectedVecSize) {
        log::cout() << "The imported solution vector is not compatible with the matrix" << std::endl;
        log::cout() << "The size of the imported solution vector is " << vecSize << std::endl;
        log::cout() << "The expected size of solution vector is " << expectedVecSize << std::endl;
        throw std::runtime_error("The imported solution vector is not compatible with the matrix");
    }
}

/*!
 * Dump the system to files.
 *
 * Only the contents of the system will be dumped, this include the matrix, the RHS vector,
 * the solution vector, and the information needed to properly restore the aforementioned
 * data structures (e.g., the block size or the transpose flag).
 *
 * \param directory is the directory where the files will be saved
 * \param prefix is the prefix that will be added to the files
 */
void SystemSolver::dumpSystem(const std::string &directory, const std::string &prefix) const
{
    // Dump system information
    int blockSize = getBlockSize();

    std::string filePath = getInfoFilePath(directory, prefix);
    std::ofstream fileStream(filePath);
    std::string fileHeader = "transpose blocksize partitioned";
    fileStream << fileHeader << std::endl;
    fileStream << m_transpose;
    fileStream << " " << blockSize;
#if BITPIT_ENABLE_MPI==1
    fileStream << " " << m_partitioned;
#endif
    fileStream << std::endl;
    fileStream.close();

    // Dump matrix
    std::string matrixFilePath = getMatrixFilePath(directory, prefix);
    exportMatrix(matrixFilePath, FILE_BINARY);

    // Dump RHS
    std::string rhsFilePath = getRHSFilePath(directory, prefix);
    exportRHS(rhsFilePath, FILE_BINARY);

    // Dump solution
    std::string solutionFilePath = getSolutionFilePath(directory, prefix);
    exportSolution(solutionFilePath, FILE_BINARY);
}

/*!
 * Restore the system from files.
 *
 * Only the contents of the system will be restored, this include the matrix, the RHS vector,
 * the solution vector, and the information needed to properly create the aforementioned data
 * structures (e.g., the block size or the transpose flag).
 *
 * Matrix file is mandatory, whereas solution and RHS files are not. If no solution or RHS
 * files are found they will be re-created from scratch (but they will not be initialized).
 */
#if BITPIT_ENABLE_MPI==1
 /*!
 * \param communicator is the MPI communicator
 * \param directory is the directory where the files will be read from
 * \param prefix is the prefix that will be was added to the files during
 */
void SystemSolver::restoreSystem(MPI_Comm communicator, const std::string &directory, const std::string &prefix)
#else
 /*!
 * \param directory is the directory where the files will be read from
 * \param prefix is the prefix that will be was added to the files during
 */
void SystemSolver::restoreSystem(const std::string &directory, const std::string &prefix)
#endif
{
    // Read system information
    std::string filePath = getInfoFilePath(directory, prefix);
    std::ifstream fileStream(filePath);
    if (!fileStream.good()) {
        throw std::runtime_error("The system solver info file \"" + filePath + "\" doesn't exists");
    }

    std::string header;
    std::getline(fileStream, header);

    bool transpose;
    fileStream >> transpose;

    int blockSize;
    fileStream >> blockSize;

#if BITPIT_ENABLE_MPI == 1
    bool isPartitioned;
    fileStream >> isPartitioned;
#endif

    fileStream.close();

    // Clear the system
    clear();

    // Set transpose flag
    setTranspose(transpose);

#if BITPIT_ENABLE_MPI == 1
    // Set the communicator
    setCommunicator(communicator);

    // Detect if the system is partitioned
    m_partitioned = isPartitioned;
#endif

    // Restore the matrix
    matrixCreate(blockSize);

    std::string matrixFilePath = getMatrixFilePath(directory, prefix);
    matrixFill(matrixFilePath);

    // Initialize RHS and solution vectors
    vectorsCreate();

    // The system is now assembled
    m_assembled = true;

    // Restore RHS and solution vectors
    std::string rhsFilePath      = getRHSFilePath(directory, prefix);
    std::string solutionFilePath = getSolutionFilePath(directory, prefix);
    vectorsFill(rhsFilePath, solutionFilePath);
}

/*!
 * Get the path of the file that will be used to dump/restore system information.
 *
 * \param directory is the directory that contains the file
 * \param prefix is the prefix that will be was added to the file
 * \result The path of the file that will be used to dump/restore system information.
 */
std::string SystemSolver::getInfoFilePath(const std::string &directory, const std::string &prefix) const
{
    std::string path = directory + "/" + prefix + "info.txt";

    return path;
}

/*!
 * Get the path of the file that will be used to dump/restore the matrix.
 *
 * \param directory is the directory that contains the file
 * \param prefix is the prefix that will be was added to the file
 * \result The path of the file that will be used to dump/restore the matrix.
 */
std::string SystemSolver::getMatrixFilePath(const std::string &directory, const std::string &prefix) const
{
    std::string path = directory + "/" + prefix + "A.txt";

    return path;
}

/*!
 * Get the path of the file that will be used to dump/restore the RHS vector.
 *
 * \param directory is the directory that contains the file
 * \param prefix is the prefix that will be was added to the file
 * \result The path of the file that will be used to dump/restore the RHS vector.
 */
std::string SystemSolver::getRHSFilePath(const std::string &directory, const std::string &prefix) const
{
    std::string path = directory + "/" + prefix + "rhs.txt";

    return path;
}

/*!
 * Get the path of the file that will be used to dump/restore the solution vector.
 *
 * \param directory is the directory that contains the file
 * \param prefix is the prefix that will be was added to the file
 * \result The path of the file that will be used to dump/restore the solution vector.
 */
std::string SystemSolver::getSolutionFilePath(const std::string &directory, const std::string &prefix) const
{
    std::string path = directory + "/" + prefix + "solution.txt";

    return path;
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
 * Set the reordering that will be applied when assemblying the matrix.
 *
 * Reordering will be applied when the system is assembled and its sole purpose
 * is to speed up the resolution of the system (e.g., reorder can be used to
 * reduce the fill-in of the LU factorization).
 *
 * Reordering is only applied internally, all public functions expects row
 * and column indices to be in natural matrix order (i.e., not reordered).
 *
 * \param nRows is the number of rows of the system matrix
 * \param nCols is the number of columns of the system matrix
 * \param reordering is the reordering that will be applied
 */
void SystemSolver::setReordering(long nRows, long nCols, const SystemMatrixOrdering &reordering)
{
    // Clear existing reordering
    clearReordering();

    // Early return if natural ordering is used
    try {
        dynamic_cast<const NaturalSystemMatrixOrdering &>(reordering);
        return;
    } catch(const std::bad_cast &exception) {
        BITPIT_UNUSED(exception);

        // A reordering other than the natural one has been passed
    }

    // Set row reordering
    PetscInt *rowReorderingStorage;
    PetscMalloc(nRows * sizeof(PetscInt), &rowReorderingStorage);
    for (long i = 0; i < nRows; ++i) {
        rowReorderingStorage[reordering.getRowRank(i)] = i;
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(m_communicator, nRows, rowReorderingStorage, PETSC_OWN_POINTER, &m_rowReordering);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nRows, rowReorderingStorage, PETSC_OWN_POINTER, &m_rowReordering);
#endif
    ISSetPermutation(m_rowReordering);

    // Create new permutations
    PetscInt *colReorderingStorage;
    PetscMalloc(nCols * sizeof(PetscInt), &colReorderingStorage);
    for (long j = 0; j < nCols; ++j) {
        colReorderingStorage[reordering.getColRank(j)] = j;
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(m_communicator, nCols, colReorderingStorage, PETSC_OWN_POINTER, &m_colReordering);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nCols, colReorderingStorage, PETSC_OWN_POINTER, &m_colReordering);
#endif
    ISSetPermutation(m_colReordering);
}

/*!
 * Clear the reordering that will be applied when assemblying the matrix.
 *
 * The function will clear any reordering preiously set. With no reordering
 * defined, the matrix will be assembled using its natural ordering.
 */
void SystemSolver::clearReordering()
{
    if (m_rowReordering) {
        ISDestroy(&m_rowReordering);
    }

    if (m_colReordering) {
        ISDestroy(&m_colReordering);
    }
}

/*!
 * Setup the system.
 *
 * This function is deprecated, there is no need to call it.
 */
void SystemSolver::setUp()
{
    // Nothing to do
}

/*!
 * Prepare the KSP for solving the system.
 */
void SystemSolver::prepareKSP()
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The system is not yet assembled.");
    }

    // Early return if the preconditioner can be reused
    if (!m_KSPDirty) {
        return;
    }

    // Create KSP
    bool setupNeeded = false;
    if (!m_KSP) {
        createKSP();
        setupNeeded = true;
    }

    // Set the matrix associated with the linear system
    KSPSetOperators(m_KSP, m_A, m_A);

    // Set up
    if (setupNeeded) {
        // Perform actions before KSP set up
        preKSPSetupActions();

        // KSP set up
        KSPSetFromOptions(m_KSP);
        KSPSetUp(m_KSP);

        // Perform actions after KSP set up
        postKSPSetupActions();
    }

    // KSP is now ready
    m_KSPDirty = false;
}

/*!
 * Finalize the KSP after the system has been solved.
 */
void SystemSolver::finalizeKSP()
{
    // Nothing to do
}

/*!
 * Create the KSP.
 */
void SystemSolver::createKSP()
{
    // Create KSP object
#if BITPIT_ENABLE_MPI==1
    KSPCreate(m_communicator, &m_KSP);
#else
    KSPCreate(PETSC_COMM_SELF, &m_KSP);
#endif

    // Set options prefix
    if (!m_prefix.empty()) {
        KSPSetOptionsPrefix(m_KSP, m_prefix.c_str());
    }
}

/*!
 * Destroy the KSP.
 */
void SystemSolver::destroyKSP()
{
    m_KSPDirty = true;
    KSPDestroy(&m_KSP);
    m_KSP = nullptr;
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
