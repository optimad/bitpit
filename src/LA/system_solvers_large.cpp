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

#include "system_solvers_large.hpp"

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_operators.hpp"

#include "petscksp.h"
#include "petscmat.h"
#include "petscsystypes.h"
#include "petscvec.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <numeric>
#include <unordered_set>

#ifndef PETSC_NULLPTR
#define PETSC_NULLPTR PETSC_NULL
#endif

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
 * Get the permutation rank of the specified local row.
 *
 * The permutation rank defines the position, in the final assembled matrix, of the
 * specified local row (i.e., the position of the specified local row after applying
 * the reordering).
 *
 * \param row is the local row
 * \result The permutation rank of the specified local row.
 */
long NaturalSystemMatrixOrdering::getRowPermutationRank(long row) const
{
    return row;
}

/*!
 * Get the permutation rank of the specified local column.
 *
 * The permutation rank defines the position, in the final assembled matrix, of the
 * specified local column (i.e., the position of the specified local column after
 * applying the reordering).
 *
 * \param col is the local column
 * \result The permutation rank of the specified local column.
 */
long NaturalSystemMatrixOrdering::getColPermutationRank(long col) const
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

#if BITPIT_ENABLE_MPI==1
/*!
 * Checks if the matrix is partitioned.
 *
 * \result Returns true if the patch is partitioned, false otherwise.
 */
bool SystemSparseMatrixAssembler::isPartitioned() const
{
    return m_matrix->isPartitioned();
}

/*!
 * Gets the MPI communicator associated to the matrix.
 *
 * \return The MPI communicator associated to the matrix.
 */
const MPI_Comm & SystemSparseMatrixAssembler::getCommunicator() const
{
    return m_matrix->getCommunicator();
}
#endif

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
 * both the solution vector and the RHS one will be destroyed and re-created.
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
 * defined as a group of block-size matrix rows.
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
 * is defined as a group of block-size matrix columns.
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
 * is defined as a group of block-size matrix rows.
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
 * column is defined as a group of block-size matrix columns.
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
 * where a block column is defined as a group of block-size matrix columns.
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
 * where a block column is defined as a group of block-size matrix columns. The
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
 * block column is defined as a group of block-size matrix columns.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * the values contain of all the elements of a block row, where a block column
 * is defined as a group of block-size matrix columns. The values are returned
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
 * \result Return true is PETSc has been initialized by this function, false
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
 *
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * The system is solved using a flexible GMRES iterative method. If the system
 * is partitioned it is preconditioned using the (restricted) additive Schwarz
 * method (ASM). On each block of the ASM preconditioner an incomplete LU
 * factorization (ILU) is used. There is one block per process. If the system
 * is not partitioned it is preconditioned using the incomplete LU factorization
 * (ILU).
 */
#else
/*!
 * The system is solved using a flexible GMRES iterative method preconditioned
 * using the incomplete LU factorization (ILU).
 */
#endif

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
 * \param flatten if set to true, block size will not be taken into account when allocating
 * the internal storage of the system matrix. Even when the internal storage is flat, block
 * size information are available and can be used by the solver to speed up the solution of
 * the system. However, since the internal storage doesn't take blocks into account, some
 * low level operations (e.g., matrix-matrix multiplications) cannot use block information.
 * Some algorithms for the solution of the system (e.g., multigrid) requires a flat storage.
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
 * \param flatten if set to true, block size will not be taken into account when allocating
 * the internal storage of the system matrix. Even when the internal storage is flat, block
 * size information are available and can be used by the solver to speed up the solution of
 * the system. However, since the internal storage doesn't take blocks into account, some
 * low level operations (e.g., matrix-matrix multiplications) cannot use block information.
 * Some algorithms for the solution of the system (e.g., multigrid) requires a flat storage.
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(const std::string &prefix, bool flatten, bool transpose, bool debug)
    : m_flatten(flatten), m_transpose(transpose),
      m_A(PETSC_NULLPTR), m_rhs(PETSC_NULLPTR), m_solution(PETSC_NULLPTR),
      m_rowReordering(PETSC_NULLPTR), m_colReordering(PETSC_NULLPTR),
      m_convergenceMonitorEnabled(debug),
      m_KSP(PETSC_NULLPTR), m_KSPDirty(true),
      m_prefix(prefix), m_assembled(false),
#if BITPIT_ENABLE_MPI==1
      m_communicator(MPI_COMM_SELF), m_partitioned(false),
#endif
      m_forceConsistency(false)
{
    // Initialize PETSc
    if (m_nInstances == 0) {
        m_petscManager.initialize(debug);
    }

    // Increase the number of instances
    ++m_nInstances;

    // Reset KSP options
    initializeKSPOptions();

    // Reset KSP status
    initializeKSPStatus();
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

    destroyKSPOptions();
    destroyKSPStatus();

    vectorsDestroy();
    matrixDestroy();

    clearReordering();

#if BITPIT_ENABLE_MPI==1
    freeCommunicator();
#endif

    m_assembled = false;
}

/*!
 * Clear and release the memory of all data structures needed for the solution of the system.
 *
 * These data structures will be re-created the next time the system will be solved.
 */
void SystemSolver::clearWorkspace()
{
    // Early return if the KSP has not been created
    if (!m_KSP) {
        return;
    }

    // Destroy the KSP
    destroyKSP();
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
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
 * After assembying th system solver, its options will be reset.
 *
 * \param matrix is the matrix
 * \param reordering is the reordering that will be applied when assembling the system
 */
void SystemSolver::assembly(const SparseMatrix &matrix, const SystemMatrixOrdering &reordering)
{
    // Check if the matrix is assembled
    if (!matrix.isAssembled()) {
        throw std::runtime_error("Unable to assembly the system. The matrix is not yet assembled.");
    }

    // Assembly the system matrix
    SystemSparseMatrixAssembler assembler(&matrix);
    assembly<SystemSolver>(static_cast<const Assembler &>(assembler), reordering);
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::assembly(const Assembler &assembler)
{
    assembly(assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
 *
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assembling the system
 */
void SystemSolver::assembly(const Assembler &assembler, const SystemMatrixOrdering &reordering)
{
    assembly<SystemSolver>(assembler, reordering);
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
    update<SystemSolver>(nRows, rows, static_cast<const Assembler &>(assembler));
}

/*!
 * Update all the rows of the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::update(const Assembler &assembler)
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
void SystemSolver::update(long nRows, const long *rows, const Assembler &assembler)
{
    update<SystemSolver>(nRows, rows, assembler);
}

/*!
 * Get the block size of the system.
 *
 * \result The block size of the system.
 */
int SystemSolver::getBlockSize() const
{
    if (m_A == PETSC_NULLPTR) {
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
 * defined as a group of block-size matrix rows.
 *
 * \result The number of rows of the system.
 */
long SystemSolver::getRowCount() const
{
    if (m_A == PETSC_NULLPTR) {
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
 * is defined as a group of block-size matrix columns.
 *
 * \result The number of columns of the system.
 */
long SystemSolver::getColCount() const
{
    if (m_A == PETSC_NULLPTR) {
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
    if (m_A == PETSC_NULLPTR) {
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
    if (m_A == PETSC_NULLPTR) {
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
 * is defined as a group of block-size matrix rows.
 *
 * \result The number of global rows
 */
long SystemSolver::getRowGlobalCount() const
{
    if (m_A == PETSC_NULLPTR) {
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
 * column is defined as a group of block-size matrix columns.
 *
 * \result The number of global (block) columns.
 */
long SystemSolver::getColGlobalCount() const
{
    if (m_A == PETSC_NULLPTR) {
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
    if (m_A == PETSC_NULLPTR) {
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
    if (m_A == PETSC_NULLPTR) {
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

    // Perform actions before KSP solution
    preKSPSolveActions();

    // Solve KSP
    solveKSP();

    // Perform actions after KSP solution
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
    vectorsFill(rhs, *solution);

    // Solve the system
    solve();

    // Export the solution
    exportVector(m_solution, solution);
}

/*!
 * Solve KSP.
 */
void SystemSolver::solveKSP()
{
    PetscErrorCode solverError;
    if (!m_transpose) {
        solverError = KSPSolve(m_KSP, m_rhs, m_solution);
    } else {
        solverError = KSPSolveTranspose(m_KSP, m_rhs, m_solution);
    }

    if (solverError) {
        const char *petscMessage = nullptr;
        PetscErrorMessage(solverError, &petscMessage, PETSC_NULLPTR);
        std::string message = "Unable to solver the system. " + std::string(petscMessage);
        throw std::runtime_error(message);
    }
}

/*!
 * Pre-solve actions.
 */
void SystemSolver::preKSPSolveActions()
{
    // Prepare KSP
    prepareKSP();

    // Reorder vectors
    vectorsReorder(true);

    // Force consistency
    if (m_forceConsistency) {
        removeNullSpaceFromRHS();
    }
}

/*!
 * Post-solve actions.
 */
void SystemSolver::postKSPSolveActions()
{
    // Fill status of KSP
    fillKSPStatus();

    // Reorder vectors
    vectorsReorder(false);

    // Finalize KSP
    finalizeKSP();
}

/*!
 *  Get the version associated with the binary dumps.
 *
 *  \result The version associated with the binary dumps.
 */
int SystemSolver::getDumpVersion() const
{
    const int DUMP_VERSION = 1;

    return DUMP_VERSION;
}

/*!
 * Assemble the matrix.
 *
 * \param assembler is the matrix assembler
 */
void SystemSolver::matrixAssembly(const Assembler &assembler)
{
    const PetscInt *rowReordering = PETSC_NULLPTR;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    // Create the matrix
    int blockSize = assembler.getBlockSize();
    createMatrix(blockSize, blockSize, &m_A);

    MatType matrixType;
    MatGetType(m_A, &matrixType);

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

    // Allocate storage
    //
    // When the internal storage of the system matrix was created without taking into account
    // block information, preallocation information should be provided for each row of each
    // block.
    int allocationExpansion;
    if (strcmp(matrixType, MATSEQAIJ) == 0) {
        allocationExpansion = blockSize;
#if BITPIT_ENABLE_MPI == 1
    } else if (strcmp(matrixType, MATMPIAIJ) == 0) {
        allocationExpansion = blockSize;
#endif
    } else {
        allocationExpansion = 1;
    }

    long nAllocatedRows = allocationExpansion * nAssemblerRows;

    std::vector<int> d_nnz(nAllocatedRows, 0);
    for (long n = 0; n < nAssemblerRows; ++n) {
        long matrixRow = n;
        if (rowReordering) {
            matrixRow = rowReordering[matrixRow];
        }

        int nAssemblerRowNZ = assembler.getRowNZCount(n);

        long matrixRowOffset = matrixRow * allocationExpansion;
        for (int i = 0; i < allocationExpansion; ++i) {
            d_nnz[matrixRowOffset + i] = allocationExpansion * nAssemblerRowNZ;
        }
    }

#if BITPIT_ENABLE_MPI == 1
    std::vector<int> o_nnz(nAllocatedRows, 0);
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

            long matrixRowOffset = matrixRow * allocationExpansion;
            for (int k = 0; k < nAssemblerRowNZ; ++k) {
                long id = assemblerRowPattern[k];
                if (id < assemblerDiagonalBegin || id >= assemblerDiagonalEnd) {
                    for (int i = 0; i < allocationExpansion; ++i) {
                        o_nnz[matrixRowOffset + i] += allocationExpansion;
                        d_nnz[matrixRowOffset + i] -= allocationExpansion;
                    }
                }
            }
        }
    }
#endif

    if (strcmp(matrixType, MATSEQAIJ) == 0) {
        MatSeqAIJSetPreallocation(m_A, 0, d_nnz.data());
    } else if (strcmp(matrixType, MATSEQBAIJ) == 0) {
        MatSeqBAIJSetPreallocation(m_A, blockSize, 0, d_nnz.data());
#if BITPIT_ENABLE_MPI == 1
    } else if (strcmp(matrixType, MATMPIAIJ) == 0) {
        MatMPIAIJSetPreallocation(m_A, 0, d_nnz.data(), 0, o_nnz.data());
    } else if (strcmp(matrixType, MATMPIBAIJ) == 0) {
        MatMPIBAIJSetPreallocation(m_A, blockSize, 0, d_nnz.data(), 0, o_nnz.data());
#endif
    } else {
        throw std::runtime_error("Matrix format not supported.");
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

    // Fill matrix
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
 * The input file should contain a compatible vector stored in PETSc binary format. It's up to
 * the caller of this routine to make sure the loaded vector is compatible with the matrix. If
 * the matrix file cannot be read an exception is thrown.
 *
 * \param filePath is the path of the file
 */
void SystemSolver::matrixFill(const std::string &filePath)
{
    // Check if the matrix exists
    if (!m_A) {
        throw std::runtime_error("Matrix should be created before filling it.");
    }

    // Fill the matrix
    fillMatrix(m_A, filePath);
}

/*!
 * Update the specified rows of the matrix.
 *
 * The contents of the specified rows will be replaced by the data provided by the given
 * assembler. If the matrix has not been assembled yet, both the pattern and the values
 * of the matrix will be updated. After the matrix has been assembled only the values
 * will be updated.
 *
 * The block size of the assembler should be equal to the block size of the matrix.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the local indices of the rows that will be updated, if a
 * null pointer is passed, the rows that will be updated are the rows
 * from 0 to (nRows - 1).
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SystemSolver::matrixUpdate(long nRows, const long *rows, const Assembler &assembler)
{
    // Updating the matrix invalidates the KSP
    m_KSPDirty = true;

    // Get block size
    const int blockSize = getBlockSize();
    if (assembler.getBlockSize() != blockSize) {
        std::string message = "Unable to update the matrix.";
        message += " The block size of the assembler is not equal to the block size of the system matrix.";
        throw std::runtime_error(message);
    }

    // Initialize reordering
    const PetscInt *rowReordering = PETSC_NULLPTR;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    const PetscInt *colReordering = PETSC_NULLPTR;
    if (m_colReordering) {
        ISGetIndices(m_colReordering, &colReordering);
    }

    // Global information
    PetscInt colGlobalBegin;
    PetscInt colGlobalEnd;
    MatGetOwnershipRangeColumn(m_A, &colGlobalBegin, &colGlobalEnd);
    colGlobalBegin /= blockSize;
    colGlobalEnd /= blockSize;

    PetscInt rowGlobalOffset;
    MatGetOwnershipRange(m_A, &rowGlobalOffset, PETSC_NULLPTR);
    rowGlobalOffset /= blockSize;

    // Get the options for assembling the matrix
    SystemMatrixAssembler::AssemblyOptions assemblyOptions = assembler.getOptions();

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
    PetscBool matrixSortedFull = (assemblyOptions.full && assemblyOptions.sorted) ? PETSC_TRUE : PETSC_FALSE;
    MatSetOption(m_A, MAT_SORTED_FULL, matrixSortedFull);
#endif

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
    bool fastUpdate = isAssembled() && (blockSize == 1) && assemblyOptions.full && assemblyOptions.sorted;

    // Update element values
    //
    // If the sizes of PETSc data types match the sizes of data types expected by
    // bitpit a direct update can be performed, otherwise the matrix is updated
    // using intermediate data storages.
    const long assemblerMaxRowNZ = std::max(assembler.getMaxRowNZCount(), 0L);

    bool patternDirectUpdate = !colReordering && (sizeof(long) == sizeof(PetscInt));
    bool valuesDirectUpdate  = (sizeof(double) == sizeof(PetscScalar));

    ConstProxyVector<long> rowPattern;
    std::vector<PetscInt> petscRowPatternStorage;
    const PetscInt *petscRowPattern;
    if (!patternDirectUpdate) {
        rowPattern.set(ConstProxyVector<long>::INTERNAL_STORAGE, 0, assemblerMaxRowNZ);
        petscRowPatternStorage.resize(assemblerMaxRowNZ);
        petscRowPattern = petscRowPatternStorage.data();
    }

    ConstProxyVector<double> rowValues;
    std::vector<PetscScalar> petscRowValuesStorage;
    const PetscScalar *petscRowValues;
    if (!valuesDirectUpdate) {
        long assemblerMaxRowNZElements = blockSize * blockSize * assemblerMaxRowNZ;
        rowValues.set(ConstProxyVector<double>::INTERNAL_STORAGE, 0, assemblerMaxRowNZElements);
        petscRowValuesStorage.resize(assemblerMaxRowNZElements);
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

        // Get row data
        if (fastUpdate) {
            assembler.getRowValues(n, &rowValues);
        } else {
            assembler.getRowData(n, &rowPattern, &rowValues);
        }

        if (rowValues.size() == 0) {
            continue;
        }

        // Get values in PETSc format
        const std::size_t rowPatternSize = rowPattern.size();

        if (valuesDirectUpdate) {
            petscRowValues = reinterpret_cast<const PetscScalar *>(rowValues.data());
        } else {
            std::copy(rowValues.cbegin(), rowValues.cend(), petscRowValuesStorage.begin());
        }

        if (fastUpdate) {
            MatSetValuesRow(m_A, globalRow, petscRowValues);
        } else {
            // Get pattern in PETSc format
            if (patternDirectUpdate) {
                petscRowPattern = reinterpret_cast<const PetscInt *>(rowPattern.data());
            } else {
                for (std::size_t k = 0; k < rowPatternSize; ++k) {
                    long globalCol = rowPattern[k];
                    if (colReordering) {
                        if (globalCol >= colGlobalBegin && globalCol < colGlobalEnd) {
                            long col = globalCol - colGlobalBegin;
                            col = colReordering[col];
                            globalCol = colGlobalBegin + col;
                        }
                    }

                    petscRowPatternStorage[k] = globalCol;
                }
            }

            // Set data
            if (blockSize > 1) {
                MatSetValuesBlocked(m_A, 1, &globalRow, rowPatternSize, petscRowPattern, petscRowValues, INSERT_VALUES);
            } else {
                MatSetValues(m_A, 1, &globalRow, rowPatternSize, petscRowPattern, petscRowValues, INSERT_VALUES);
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
 * Dump the matrix.
 *
 * \param systemStream is the stream in which system information is written
 * \param directory is the directory in which the matrix data file will be written
 * \param prefix is the prefix added to the name of the file containing matrix data
 */
void SystemSolver::matrixDump(std::ostream &systemStream, const std::string &directory,
                              const std::string &prefix) const
{
    BITPIT_UNUSED(systemStream);

    dumpMatrix(m_A, directory, prefix + "A");
}

/*!
 * Restore the matrix.
 *
 * \param systemStream is the stream from which system information is read
 * \param directory is the directory from which the matrix data file will be read
 * \param prefix is the prefix that will be was added to the files during
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * \param redistribute if set to true, the matrix will be redistributed among the available
 * processes, allowing to restore the matrix with a different number of processes than those
 * used to dump it
 */
void SystemSolver::matrixRestore(std::istream &systemStream, const std::string &directory,
                                 const std::string &prefix, bool redistribute)
#else
void SystemSolver::matrixRestore(std::istream &systemStream, const std::string &directory,
                                 const std::string &prefix)
#endif
{
    BITPIT_UNUSED(systemStream);

#if BITPIT_ENABLE_MPI==1
    restoreMatrix(directory, prefix + "A", redistribute, &m_A);
#else
    restoreMatrix(directory, prefix + "A", &m_A);
#endif
}

/*!
 * Destroy the matrix.
 */
void SystemSolver::matrixDestroy()
{
    destroyMatrix(&m_A);
}

/*!
 * Create RHS and solution vectors.
 *
 * Vectors will be created, but they will not be initialized.
 */
void SystemSolver::vectorsCreate()
{
    if (!m_transpose) {
        MatCreateVecs(m_A, &m_solution, &m_rhs);
    } else {
        MatCreateVecs(m_A, &m_rhs, &m_solution);
    }
}

/*!
 * Fill RHS and solution vectors using the given data.
 *
 * \param rhs contains the that that will be copied into the RHS vector, if the vector is
 * empty no data will be copied and function will leave the RHS vector unaltered
 * \param solution contains the that that will be copied into the solution vector, if the vector is
 * empty no data will be copied and function will leave the RHS vector unaltered
 */
void SystemSolver::vectorsFill(const std::vector<double> &rhs, const std::vector<double> &solution)
{
    if (!rhs.empty()) {
        fillVector(m_rhs, rhs);
    }

    if (!solution.empty()) {
        fillVector(m_solution, solution);
    }
}

/*!
 * Fill RHS and solution vectors reading their contents from the specified file.
 *
 * \param rhsFilePath is the file path containing the content of the RHS vector, if the path is
 * empty the function will leave the RHS vector unaltered
 * \param solutionFilePath is the file path containing the content of the solution vector, if
 * the path is empty the function will leave the solution vector unaltered
 */
void SystemSolver::vectorsFill(const std::string &rhsFilePath, const std::string &solutionFilePath)
{
    if (!rhsFilePath.empty()) {
        fillVector(m_rhs, rhsFilePath);
    }

    if (!solutionFilePath.empty()) {
        fillVector(m_solution, solutionFilePath);
    }
}

/*!
 * Reorder RHS and solution vectors to match the order of the system matrix.
 *
 * \param invert is a flag for inverting the ordering
 */
void SystemSolver::vectorsReorder(bool invert)
{
    if (!m_transpose) {
        reorderVector(m_rhs, m_colReordering, invert);
        reorderVector(m_solution, m_rowReordering, invert);
    } else {
        reorderVector(m_rhs, m_rowReordering, invert);
        reorderVector(m_solution, m_colReordering, invert);
    }
}

/*!
 * Dump RHS and solution vectors.
 *
 * \param systemStream is the stream from which system information is read
 * \param directory is the directory in which the vectors data file will be written
 * \param prefix is the prefix added to the names of the files containing vector data
 */
void SystemSolver::vectorsDump(std::ostream &systemStream, const std::string &directory,
                               const std::string &prefix) const
{
    BITPIT_UNUSED(systemStream);

    dumpVector(m_rhs, directory, prefix + "rhs");
    dumpVector(m_solution, directory, prefix + "solution");
}

/*!
 * Restore RHS and solution vectors.
 *
 * \param systemStream is the stream that contains system information
 * \param directory is the directory from which the vector data file will be read
 * \param prefix is the prefix added to the name of the file containing vectors data
 */
void SystemSolver::vectorsRestore(std::istream &systemStream, const std::string &directory,
                                  const std::string &prefix)
{
    BITPIT_UNUSED(systemStream);

    restoreVector(directory, prefix + "rhs", m_A, VECTOR_SIDE_RIGHT, &m_rhs);
    restoreVector(directory, prefix + "solution", m_A, VECTOR_SIDE_LEFT, &m_solution);
}

/*!
 * Destroy RHS and solution vectors.
 */
void SystemSolver::vectorsDestroy()
{
    destroyVector(&m_rhs);
    destroyVector(&m_solution);
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
 * Export the matrix to the specified file.
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the matrix, note that
 * the ASCII format may not be able to handle large matrices
 */
void SystemSolver::exportMatrix(const std::string &filePath, FileFormat fileFormat) const
{
    exportMatrix(m_A, filePath, fileFormat);
}

/*!
 * Import the matrix reading its contents form the specified file.
 *
 * The input file should contain a compatible matrix stored in PETSc binary format. It's up to the
 * caller of this routine to make sure the loaded matrix is compatible with the system. If the
 * matrix file cannot be read an exception is thrown.
 *
 * \param filePath is the path of the file
 */
void SystemSolver::importMatrix(const std::string &filePath)
{
    // Check if the matrix exists
    if (!m_A) {
        throw std::runtime_error("Matrix should be created before filling it.");
    }

    // Clear workspace
    clearWorkspace();

    // Clear reordering
    clearReordering();

    // Fill the matrix
    matrixFill(filePath);

    // Re-create vectors
    vectorsDestroy();
    vectorsCreate();
}

/*!
 * Export the RHS vector to the specified file.
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the RHS vector,
 * note that the ASCII format may not be able to handle large vectors
 */
void SystemSolver::exportRHS(const std::string &filePath, FileFormat fileFormat) const
{
    exportVector(m_rhs, filePath, fileFormat);
}

/*!
 * Import the RHS vector reading its contents from the specified file.
 *
 * The input file should contain a compatible vector stored in PETSc binary format. If the size
 * of the vector stored in the file is not compatible with the matrix, an exception is thrown.
 * An exception is also raised if the file cannot be read.
 *
 * It is possible to fill the RHS vector only after the system has been assembled.
 *
 * \param filePath is the path of the file that contains RHS vector data
 */
void SystemSolver::importRHS(const std::string &filePath)
{
    // Check if the system is assembled
    if (!m_assembled) {
        throw std::runtime_error("The RHS vector can be loaded only after assembling the system.");
    }

    // Fill the RHS vector
    vectorsFill(filePath, "");

    // Check if the imported RHS is compatible with the matrix
    PetscInt size;
    VecGetLocalSize(m_rhs, &size);

    PetscInt expectedSize;
    if (!m_transpose) {
        expectedSize = getColCount();
    } else {
        expectedSize = getRowCount();
    }
    expectedSize *= getBlockSize();

    if (size != expectedSize) {
        log::cout() << "The imported RHS vector is not compatible with the matrix" << std::endl;
        log::cout() << "The size of the imported RHS vector is " << size << std::endl;
        log::cout() << "The expected size of RHS vector is " << expectedSize << std::endl;
        throw std::runtime_error("The imported RHS vector is not compatible with the matrix");
    }
}

/*!
 * Export the solution vector to the specified file.
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the solution vector,
 * note that the ASCII format may not be able to handle large vectors
 */
void SystemSolver::exportSolution(const std::string &filePath, FileFormat fileFormat) const
{
    exportVector(m_solution, filePath, fileFormat);
}

/*!
 * Import the solution vector reading its contents from the specified file.
 *
 * The input file should contain a compatible vector stored in PETSc binary format. If the size
 * of the vector stored in the file is not compatible with the matrix, an exception is thrown.
 * An exception is also raised if the file cannot be read.
 *
 * It is possible to fill the solution vector only after the system has been assembled.
 *
 * \param filePath is the path of the file that contains solution vector data
 */
void SystemSolver::importSolution(const std::string &filePath)
{
    // Check if the system is assembled
    if (!m_assembled) {
        throw std::runtime_error("The solution vector can be loaded only after assembling the system.");
    }

    // Fill the solution vector
    vectorsFill("", filePath);

    // Check if the imported solution is compatible with the matrix
    PetscInt size;
    VecGetLocalSize(m_solution, &size);

    PetscInt expectedSize;
    if (!m_transpose) {
        expectedSize = getRowCount();
    } else {
        expectedSize = getColCount();
    }
    expectedSize *= getBlockSize();

    if (size != expectedSize) {
        log::cout() << "The imported solution vector is not compatible with the matrix" << std::endl;
        log::cout() << "The size of the imported solution vector is " << size << std::endl;
        log::cout() << "The expected size of solution vector is " << expectedSize << std::endl;
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
 * \param header is the header that will be written in the system information archive
 * \param directory is the directory where the files will be saved
 * \param prefix is the prefix that will be added to the files
 */
void SystemSolver::dumpSystem(const std::string &header, const std::string &directory,
                              const std::string &prefix) const
{
    int partitioningBlock = getBinaryArchiveBlock();

    // Open stream that will contain system information
    OBinaryArchive systemArchive;
    if (partitioningBlock <= 0) {
        openBinaryArchive(header, directory, prefix + "info", -1, &systemArchive);
    }
    std::ostream &systemStream = systemArchive.getStream();

    // Dump system information
    dumpInfo(systemStream);

    // Dump matrix
    matrixDump(systemStream, directory, prefix);

    // Dump vectors
    vectorsDump(systemStream, directory, prefix);

    // Open stream with system information
    closeBinaryArchive(&systemArchive);
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
 * \param redistribute if set to true, the system will be redistributed among the available
 * processes, allowing to restore the system with a different number of processes than those
 * used to dump it
 */
#endif
/*!
 * \param directory is the directory where the files will be read from
 * \param prefix is the prefix that will be was added to the files during
 */
#if BITPIT_ENABLE_MPI==1
void SystemSolver::restoreSystem(MPI_Comm communicator, bool redistribute,
                                 const std::string &directory, const std::string &prefix)
#else
void SystemSolver::restoreSystem(const std::string &directory, const std::string &prefix)
#endif
{
    // Clear the system
    clear();

#if BITPIT_ENABLE_MPI == 1
    // Set the communicator
    setCommunicator(communicator);
#endif

    // Open stream with system information
    IBinaryArchive systemArchive;
    openBinaryArchive(directory, prefix + "info", -1, &systemArchive);
    std::istream &systemStream = systemArchive.getStream();

    // Dump system information
    restoreInfo(systemStream);

    // Restore the matrix
#if BITPIT_ENABLE_MPI==1
    matrixRestore(systemStream, directory, prefix, redistribute);
#else
    matrixRestore(systemStream, directory, prefix);
#endif

    // Restore RHS and solution vectors
    vectorsRestore(systemStream, directory, prefix);

    // Close stream with system information
    closeBinaryArchive(&systemArchive);

    // The system is now assembled
    m_assembled = true;

    // Initialize KSP options
    initializeKSPOptions();

    // Initialize KSP statuses
    initializeKSPStatus();
}

/*!
 * Dump system information.
 *
 * \param systemStream is the stream in which system information is written
 */
void SystemSolver::dumpInfo(std::ostream &systemStream) const
{
    if (systemStream.good()) {
#if BITPIT_ENABLE_MPI==1
        utils::binary::write(systemStream, m_partitioned);
#endif
        utils::binary::write(systemStream, m_transpose);
    }
}

/*!
 * Restore system information.
 *
 * \param systemStream is the stream from which system information is read
 */
void SystemSolver::restoreInfo(std::istream &systemStream)
{
#if BITPIT_ENABLE_MPI == 1
    // Detect if the system is partitioned
    utils::binary::read(systemStream, m_partitioned);
#endif

    // Set transpose flag
    bool transpose;
    utils::binary::read(systemStream, transpose);
    setTranspose(transpose);

}

/*!
 * Create a matrix.
 *
 * \param rowBlockSize is the row block size of the matrix
 * \param colBlockSize is the column block size of the matrix
 * \param matrix on output will contain the newly created matrix
 */
void SystemSolver::createMatrix(int rowBlockSize, int colBlockSize, Mat *matrix) const
{
    // Create the matrix
#if BITPIT_ENABLE_MPI == 1
    MatCreate(m_communicator, matrix);
#else
    MatCreate(PETSC_COMM_SELF, matrix);
#endif

    // Set matrix type
    bool creatBlockMatrix = false;
    if (m_flatten) {
        creatBlockMatrix = false;
    } else {
        creatBlockMatrix = (rowBlockSize == colBlockSize) && (rowBlockSize != 1);
    }

#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
        if (creatBlockMatrix) {
            MatSetType(*matrix, MATMPIBAIJ);
        } else {
            MatSetType(*matrix, MATMPIAIJ);
        }
    } else
#endif
    {
        if (creatBlockMatrix) {
            MatSetType(*matrix, MATSEQBAIJ);
        } else {
            MatSetType(*matrix, MATSEQAIJ);
        }
    }

    // Set block size
    if (rowBlockSize == colBlockSize) {
        MatSetBlockSize(*matrix, rowBlockSize);
    } else if (rowBlockSize != 1 || colBlockSize != 1) {
        MatSetBlockSizes(*matrix, rowBlockSize, colBlockSize);
    }
}

/*!
 * Create a nest matrix.
 *
 * \param rowBlockSize is the row block size of the matrix
 * \param colBlockSize is the column block size of the matrix
 * \param nNestRows is the number of rows in the nest, i.e. the number of submatrices along the
 * rows of the nest
 * \param nNestCols is the number of columns in the nest, i.e. the number of submatrices along
 * the columns of the nest
 * \param subMatrices are the submatrices stored in row-major order, empty submatrices can
 * be passed using nullptr
 * \param vector on output will contain the newly created matrix
 */
void SystemSolver::createMatrix(int rowBlockSize, int colBlockSize, int nNestRows, int nNestCols,
                                Mat *subMatrices, Mat *matrix) const
{
    // Create the matrix
#if BITPIT_ENABLE_MPI == 1
    MatCreateNest(getCommunicator(), nNestRows, PETSC_NULLPTR, nNestCols, PETSC_NULLPTR, subMatrices, matrix);
#else
    MatCreateNest(PETSC_COMM_SELF, nNestRows, PETSC_NULLPTR, nNestCols, PETSC_NULLPTR, subMatrices, matrix);
#endif

    // Set block size
    if (rowBlockSize == colBlockSize) {
        MatSetBlockSize(*matrix, rowBlockSize);
    } else if (rowBlockSize != 1 || colBlockSize != 1) {
        MatSetBlockSizes(*matrix, rowBlockSize, colBlockSize);
    }
}

/*!
 * Fill the given matrix reading its contents from the specified file.
 *
 * The input file should contain a compatible matrix stored in PETSc binary format. If the
 * data file cannot be read an exception is thrown.
 *
 * \param matrix is the matrix that will be filled
 * \param filePath is the path of the file that contains matrix data
 */
void SystemSolver::fillMatrix(Mat matrix, const std::string &filePath) const
{
    // Check if the file exists
    std::ifstream fileStream(filePath.c_str());
    if (!fileStream.good()) {
        throw std::runtime_error("The PETSc matrix file \"" + filePath + "\" doesn't exists.");
    }
    fileStream.close();

    // Fill matrix content
    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, filePath.c_str());
    MatLoad(matrix, viewer);
    PetscViewerDestroy(&viewer);
}

/*!
 * Dump the specified matrix.
 *
 * In addition to the file that contains the contant of the matrix, some other files will be
 * created:
 *  - a file with the information needed for re-creating the matrix (called "<NAME>.info");
 */
#if BITPIT_ENABLE_MPI==1
/*!
 *  - a file with the information about the partitioning (called "<NAME>.partitioning"), this
 *    file is created only if the system is partitioned and bitpit MPI support is enabled.
 */
#endif
/*!
 *
 * \param matrix is the matrix that will be dumped
 * \param directory is the directory in which the matrix data file will be written
 * \param name is the name of the matrix that will be dumped
 */
void SystemSolver::dumpMatrix(Mat matrix, const std::string &directory, const std::string &name) const
{
    int partitioningBlock = getBinaryArchiveBlock();

    // Store information needed to create the matrix
    OBinaryArchive infoArchive;
    if (partitioningBlock <= 0) {
        openBinaryArchive("", directory, name + ".info", -1, &infoArchive);
        std::ostream &infoStream = infoArchive.getStream();

        bool matrixExists = matrix;
        utils::binary::write(infoStream, matrixExists);
        if (!matrixExists) {
            closeBinaryArchive(&infoArchive);
            return;
        }

        PetscInt rowBlockSize;
        PetscInt colBlockSize;
        MatGetBlockSizes(matrix, &rowBlockSize, &colBlockSize);
        utils::binary::write(infoStream, static_cast<int>(rowBlockSize));
        utils::binary::write(infoStream, static_cast<int>(colBlockSize));

        closeBinaryArchive(&infoArchive);
    }

    if (!matrix) {
        return;
    }

#if BITPIT_ENABLE_MPI==1
    // Store partitioning information
    if (isPartitioned()) {
        OBinaryArchive partitioningArchive;
        openBinaryArchive("", directory, name + ".partitioning", partitioningBlock, &partitioningArchive);
        std::ostream &partitioningStream = partitioningArchive.getStream();

        PetscInt nLocalRows;
        PetscInt nLocalCols;
        MatGetLocalSize(matrix, &nLocalRows, &nLocalCols);
        utils::binary::write(partitioningStream, static_cast<std::size_t>(nLocalRows));
        utils::binary::write(partitioningStream, static_cast<std::size_t>(nLocalCols));

        PetscInt nGlobalRows;
        PetscInt nGlobalCols;
        MatGetSize(matrix, &nGlobalRows, &nGlobalCols);
        utils::binary::write(partitioningStream, static_cast<std::size_t>(nGlobalRows));
        utils::binary::write(partitioningStream, static_cast<std::size_t>(nGlobalCols));

        closeBinaryArchive(&partitioningArchive);
    }
#endif

    // Store matrix data
    std::string filePath = getDataFilePath(directory, name);
    exportMatrix(matrix, filePath, FILE_BINARY);
}

/*!
 * Restore the specified matrix.
 *
 * \param directory is the directory from which the matrix data file will be read
 * \param name is the name of the matrix that will be dumped
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * \param redistribute if set to true, the matrix will be redistributed among the available
 * processes, allowing to restore the matrix with a different number of processes than those
 * used to dump it
 */
#endif
/*!
 * \param[out] matrix on output will contain the restored matrix
 */
#if BITPIT_ENABLE_MPI==1
void SystemSolver::restoreMatrix(const std::string &directory, const std::string &name,
                                 bool redistribute, Mat *matrix) const
{
#else
void SystemSolver::restoreMatrix(const std::string &directory, const std::string &name,
                                 Mat *matrix) const
{
#endif

    // Create matrix
    IBinaryArchive infoArchive;
    openBinaryArchive(directory, name + ".info", -1, &infoArchive);
    std::istream &infoStream = infoArchive.getStream();

    bool matrixExists;
    utils::binary::read(infoStream, matrixExists);
    if (!matrixExists) {
        *matrix = PETSC_NULLPTR;
        closeBinaryArchive(&infoArchive);
        return;
    }

    int rowBlockSize;
    int colBlockSize;
    utils::binary::read(infoStream, rowBlockSize);
    utils::binary::read(infoStream, colBlockSize);
    createMatrix(rowBlockSize, colBlockSize, matrix);

    closeBinaryArchive(&infoArchive);

#if BITPIT_ENABLE_MPI==1
    // Set partitioning information
    if (isPartitioned() && !redistribute) {
        int partitioningBlock = getBinaryArchiveBlock();

        IBinaryArchive partitioningArchive;
        openBinaryArchive(directory, name + ".partitioning", partitioningBlock, &partitioningArchive);
        std::istream &partitioningStream = partitioningArchive.getStream();

        std::size_t nLocalRows;
        std::size_t nLocalCols;
        utils::binary::read(partitioningStream, nLocalRows);
        utils::binary::read(partitioningStream, nLocalCols);
        std::size_t nGlobalRows;
        std::size_t nGlobalCols;
        utils::binary::read(partitioningStream, nGlobalRows);
        utils::binary::read(partitioningStream, nGlobalCols);
        MatSetSizes(*matrix, nLocalRows, nLocalCols, nGlobalRows, nGlobalCols);

        closeBinaryArchive(&partitioningArchive);
    }
#endif

    // Fill matrix
    std::string filePath = getDataFilePath(directory, name);
    fillMatrix(*matrix, filePath);
}

/*!
 * Export the given matrix to the specified file.
 *
 * \param matrix is the matrix that will be filled
 * \param filePath is the path of the file that will contain matrix data
 * \param fileFormat is the file format that will be used for exporting the matrix, note that
 * the ASCII format may not be able to handle large matrices
 */
void SystemSolver::exportMatrix(Mat matrix, const std::string &filePath, FileFormat fileFormat) const
{
    // Early return if the matrix doesn't exist
    bool matrixExists = matrix;
    if (!matrixExists) {
        std::ofstream dataFile(filePath);
        dataFile.close();

        std::ofstream infoFile(filePath + ".info");
        infoFile.close();

        return;
    }

    // Create the viewer
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

    // Export the matrix
    MatView(matrix, matViewer);

    // Destroy the viewer
    PetscViewerDestroy(&matViewer);
}

/*!
 * Destroy the specified matrix.
 *
 * \param matrix is the matrix that will be destroyed
 */
void SystemSolver::destroyMatrix(Mat *matrix) const
{
    if (matrix) {
        MatDestroy(matrix);
        *matrix = PETSC_NULLPTR;
    }
}

/*!
 * Create a vector.
 *
 * \param blockSize is the block size of the vector
 * \param vector on output will contain the newly created vector
 */
void SystemSolver::createVector(int blockSize, Vec *vector) const
{
    // Create the vector
#if BITPIT_ENABLE_MPI == 1
    VecCreate(m_communicator, vector);
#else
    VecCreate(PETSC_COMM_SELF, vector);
#endif

    // Set vector type
#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
        VecSetType(*vector, VECMPI);
    } else
#endif
    {
        VecSetType(*vector, VECSEQ);
    }

    // Set block size
    if (blockSize != 1) {
        VecSetBlockSize(*vector, blockSize);
    }
}

/*!
 * Create a nest vector.
 *
 * \param blockSize is the block size of the vector
 * \param nestSize is the number of subvectors that will be contained in the nest
 * \param subVectors are the subvectors, empty subvectors can be passed using nullptr
 * \param vector on output will contain the newly created vector
 */
void SystemSolver::createVector(int blockSize, int nestSize, Vec *subVectors, Vec *vector) const
{
    // Create the vector
#if BITPIT_ENABLE_MPI == 1
    VecCreateNest(getCommunicator(), nestSize, PETSC_NULLPTR, subVectors, vector);
#else
    VecCreateNest(PETSC_COMM_SELF, nestSize, PETSC_NULLPTR, subVectors, vector);
#endif

    // Set block size
    if (blockSize != 1) {
        VecSetBlockSize(*vector, blockSize);
    }
}

/*!
 * Fill the given vector reading its contents from the specified file.
 *
 * The input file should contain a compatible vector stored in PETSc binary format. If the
 * data file cannot be read an exception is thrown.
 *
 * \param vector is the vector that will be filled
 * \param filePath is the path of the file that contains vector data
 */
void SystemSolver::fillVector(Vec vector, const std::string &filePath) const
{
    // Check if the file exists
    std::ifstream fileStream(filePath);
    if (!fileStream.good()) {
        throw std::runtime_error("The file \"" + filePath + "\" cannot be read.");
    }
    fileStream.close();

    // Fill the vector
    PetscViewer viewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &viewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);
#endif
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, filePath.c_str());
    VecLoad(vector, viewer);
    PetscViewerDestroy(&viewer);
}

/*!
 * Fill the specified vector with the given data.
 *
 * \param vector is the vector that will be filled
 * \param data is the data that will be copied into the vector
 */
void SystemSolver::fillVector(Vec vector, const std::vector<double> &data) const
{
    int size;
    VecGetLocalSize(vector, &size);

    PetscScalar *petscData;
    VecGetArray(vector, &petscData);
    for (int i = 0; i < size; ++i) {
        petscData[i] = data[i];
    }
    VecRestoreArray(vector, &petscData);
}

/*!
 * Dump the specified vector.
 *
 * In addition to the file that contains the contant of the vector, some other files will be
 * created:
 *  - a file with the information needed for re-creating the vector (called "<NAME>.info");
 */
#if BITPIT_ENABLE_MPI==1
/*!
 *  - a file with the information about the partitioning (called "<NAME>.partitioning"), this
 *    file is created only if the system is partitioned and bitpit MPI support is enabled.
 */
#endif
/*!
 *
 * \param vector is the vector that will be dumped
 * \param directory is the directory in which the vector data file will be written
 * \param name is the name of the vector that will be dumped
 */
void SystemSolver::dumpVector(Vec vector, const std::string &directory, const std::string &name) const
{
    int partitioningBlock = getBinaryArchiveBlock();

    // Store information needed to create the matrix
    if (partitioningBlock <= 0) {
        OBinaryArchive infoArchive;
        openBinaryArchive("", directory, name + ".info", -1, &infoArchive);
        std::ostream &infoStream = infoArchive.getStream();

        bool vectorExists = vector;
        utils::binary::write(infoStream, vectorExists);
        if (!vectorExists) {
            closeBinaryArchive(&infoArchive);
            return;
        }

        PetscInt blockSize;
        VecGetBlockSize(vector, &blockSize);
        utils::binary::write(infoStream, static_cast<int>(blockSize));

        closeBinaryArchive(&infoArchive);
    }

    if (!vector) {
        return;
    }

#if BITPIT_ENABLE_MPI==1
    // Store partitioning information
    if (isPartitioned()) {
        int partitioningBlock = getBinaryArchiveBlock();

        OBinaryArchive partitioningArchive;
        openBinaryArchive("", directory, name + ".partitioning", partitioningBlock, &partitioningArchive);
        std::ostream &partitioningStream = partitioningArchive.getStream();

        PetscInt localSize;
        VecGetLocalSize(vector, &localSize);
        utils::binary::write(partitioningStream, static_cast<std::size_t>(localSize));

        PetscInt globalSize;
        VecGetSize(vector, &globalSize);
        utils::binary::write(partitioningStream, static_cast<std::size_t>(globalSize));

        closeBinaryArchive(&partitioningArchive);
    }
#endif

    // Store vector data
    std::string filePath = getDataFilePath(directory, name);
    exportVector(vector, filePath, FILE_BINARY);
}

/*!
 * Restore the specified vector.
 *
 * \param directory is the directory from which the vector data file will be read
 * \param name is the name of the vector that will be dumped
 * \param matrix is the matrix the vectors should be compatible with
 * \param side specifies whether the vector can be multiplied by the matrix, or whether the
 * matrix-vector product can be stored in it
 * \param[out] vector on output will contain the restored vector
 */
void SystemSolver::restoreVector(const std::string &directory, const std::string &name,
                                 Mat matrix, VectorSide side, Vec *vector) const
{
    // Get size information
    PetscInt nLocalRows;
    PetscInt nLocalCols;
    MatGetLocalSize(matrix, &nLocalRows, &nLocalCols);

    PetscInt nGlobalRows;
    PetscInt nGlobalCols;
    MatGetSize(matrix, &nGlobalRows, &nGlobalCols);

    std::size_t localSize  = std::numeric_limits<std::size_t>::max();
    std::size_t globalSize = std::numeric_limits<std::size_t>::max();
    if (!m_transpose) {
        if (side == VECTOR_SIDE_RIGHT) {
            localSize  = nLocalCols;
            globalSize = nGlobalCols;
        } else if (side == VECTOR_SIDE_LEFT) {
            localSize  = nLocalRows;
            globalSize = nGlobalRows;
        }
    } else {
        if (side == VECTOR_SIDE_RIGHT) {
            localSize  = nLocalRows;
            globalSize = nGlobalRows;
        } else if (side == VECTOR_SIDE_LEFT) {
            localSize  = nLocalCols;
            globalSize = nGlobalCols;
        }
    }

    // Restore the vector
    restoreVector(directory, name, localSize, globalSize, vector);
}

/*!
 * Restore the specified vector.
 *
 * \param directory is the directory from which the vector data file will be read
 * \param name is the name of the vector that will be dumped
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * \param redistribute if set to true, the vector will be redistributed among the available
 * processes, allowing to restore the vector with a different number of processes than those
 * used to dump it
 */
#endif
/*!
 * \param[out] vector on output will contain the restored vector
 */
#if BITPIT_ENABLE_MPI==1
void SystemSolver::restoreVector(const std::string &directory, const std::string &name,
                                 bool redistribute, Vec *vector) const
#else
void SystemSolver::restoreVector(const std::string &directory, const std::string &name,
                                 Vec *vector) const
#endif
{

    // Get size information
    std::size_t localSize  = std::numeric_limits<std::size_t>::max();
    std::size_t globalSize = std::numeric_limits<std::size_t>::max();
#if BITPIT_ENABLE_MPI==1
    if (isPartitioned()) {
        if (!redistribute) {
            int partitioningBlock = getBinaryArchiveBlock();

            IBinaryArchive partitioningArchive;
            openBinaryArchive(directory, name + ".partitioning", partitioningBlock, &partitioningArchive);
            std::istream &partitioningStream = partitioningArchive.getStream();

            utils::binary::read(partitioningStream, localSize);
            utils::binary::read(partitioningStream, globalSize);

            closeBinaryArchive(&partitioningArchive);
        }
    }
#endif

    restoreVector(directory, name, localSize, globalSize, vector);
}

/*!
 * Restore the specified vector.
 *
 * \param directory is the directory from which the vector data file will be read
 * \param name is the name of the vector that will be dumped
 * \param localSize is the local size of the restored vector
 * \param globalSize is the global size of the restored vector
 * \param[out] vector on output will contain the restored vector
 */
void SystemSolver::restoreVector(const std::string &directory, const std::string &name,
                                 std::size_t localSize, std::size_t globalSize, Vec *vector) const
{
    // Create vector
    IBinaryArchive infoArchive;
    openBinaryArchive(directory, name + ".info", -1, &infoArchive);
    std::istream &infoStream = infoArchive.getStream();

    bool vectorExists;;
    utils::binary::read(infoStream, vectorExists);
    if (!vectorExists) {
        *vector = PETSC_NULLPTR;
        closeBinaryArchive(&infoArchive);
        return;
    }

    int blockSize;
    utils::binary::read(infoStream, blockSize);
    createVector(blockSize, vector);

    closeBinaryArchive(&infoArchive);

    // Set size information
    PetscInt petscLocalSize;
    if (localSize != std::numeric_limits<std::size_t>::max()) {
        petscLocalSize = static_cast<PetscInt>(localSize);
    } else {
        petscLocalSize = PETSC_DECIDE;
    }

    PetscInt petscGlobalSize;
    if (globalSize != std::numeric_limits<std::size_t>::max()) {
        petscGlobalSize = static_cast<PetscInt>(globalSize);
    } else {
        petscGlobalSize = PETSC_DETERMINE;
    }

    if (petscLocalSize != PETSC_DECIDE || petscGlobalSize != PETSC_DETERMINE) {
        VecSetSizes(*vector, petscLocalSize, petscGlobalSize);
    }

    // Fill vector data
    std::string filePath = getDataFilePath(directory, name);
    fillVector(*vector, filePath);
}

/*!
 * Reorder the specified vector.
 *
 * \param invert is a flag for inverting the ordering
 */
void SystemSolver::reorderVector(Vec vector, IS permutations, bool invert) const
{
    if (!permutations) {
        return;
    }

    PetscBool petscInvert;
    if (invert) {
        petscInvert = PETSC_TRUE;
    } else {
        petscInvert = PETSC_FALSE;
    }

    VecPermute(vector, permutations, petscInvert);
}

/*!
 * Export the specified vector to the file given file.
 *
 * \param vector is the vector that will be exported
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the vector, note that
 * the ASCII format may not be able to handle large vectors
 */
void SystemSolver::exportVector(Vec vector, const std::string &filePath, FileFormat fileFormat) const
{
    // Early return if the matrix doesn't exist
    bool vectorExists = vector;
    if (!vectorExists) {
        std::ofstream dataFile(filePath);
        dataFile.close();

        std::ofstream infoFile(filePath + ".info");
        infoFile.close();

        return;
    }

    // Create the viewer
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

    // Export the vector
    VecView(vector, viewer);

    // Destroy the viewer
    PetscViewerDestroy(&viewer);
}

/*!
 * Export the specified vector into the given container.
 *
 * \param vector is the vector that will be exported
 * \param data on output it will contain the data of the vector of the linear system
 */
void SystemSolver::exportVector(Vec vector, std::vector<double> *data) const
{
    int size;
    VecGetLocalSize(m_solution, &size);

    const PetscScalar *vectorData;
    VecGetArrayRead(m_solution, &vectorData);
    for (int i = 0; i < size; ++i) {
        (*data)[i] = vectorData[i];
    }
    VecRestoreArrayRead(vector, &vectorData);
}

/*!
 * Destroy the specified vector.
 *
 * \param vector is the vector that will be destroyed
 **/
void SystemSolver::destroyVector(Vec *vector) const
{
    if (vector) {
        VecDestroy(vector);
        *vector = PETSC_NULLPTR;
    }
}

/*!
 * Get the path of the specified info file.
 *
 * \param directory is the directory that contains the file
 * \param name is the name, without extension, of the file
 * \result The path of the info file.
 */
std::string SystemSolver::getInfoFilePath(const std::string &directory, const std::string &name) const
{
    return getFilePath(directory, name + ".info");
}

/*!
 * Get the path of the specified data file.
 *
 * \param directory is the directory that contains the file
 * \param name is the name, without extension, of the file
 * \result The path of the data file.
 */
std::string SystemSolver::getDataFilePath(const std::string &directory, const std::string &name) const
{
    return getFilePath(directory, name + ".dat");
}

/*!
 * Get the path of the specified file.
 *
 * \param directory is the directory that contains the file
 * \param name is the name of the file
 * \result The path of the file.
 */
std::string SystemSolver::getFilePath(const std::string &directory, const std::string &name) const
{
    std::string path = directory + "/" + name;

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
 * Set the reordering that will be applied when assembling the matrix.
 *
 * Reordering will be applied when the system is assembled and its sole purpose
 * is to speed up the resolution of the system (e.g., reorder can be used to
 * reduce the fill-in of the LU factorization).
 *
 * Reordering is only applied internally, all public functions expects row
 * and column indices to be in natural matrix order (i.e., not reordered).
 *
 * If the system is partitioned, each process can reorder only it's local
 * part of the matrix.
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
        rowReorderingStorage[reordering.getRowPermutationRank(i)] = i;
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
        colReorderingStorage[reordering.getColPermutationRank(j)] = j;
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(m_communicator, nCols, colReorderingStorage, PETSC_OWN_POINTER, &m_colReordering);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nCols, colReorderingStorage, PETSC_OWN_POINTER, &m_colReordering);
#endif
    ISSetPermutation(m_colReordering);
}

/*!
 * Clear the reordering that will be applied when assembling the matrix.
 *
 * The function will clear any reordering previously set. With no reordering
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
 * Prepare the KSP before the solution of the system.
 */
void SystemSolver::prepareKSP()
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The system is not yet assembled.");
    }

    // Early return if the KSP can be reused
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
        // Perform actions before preconditioner set up
        prePreconditionerSetupActions();

        // Initialize preconditioner from options
        PC pc;
        KSPGetPC(m_KSP, &pc);
        PCSetFromOptions(pc);

        // Set up preconditioner
        setupPreconditioner();

        // Perform actions after preconditioner set up
        postPreconditionerSetupActions();

        // Perform actions before Krylov subspace method set up set up
        preKrylovSetupActions();

        // Initialize Krylov subspace from options
        KSPSetFromOptions(m_KSP);

        // Set up the Krylov subspace method
        setupKrylov();

        // Perform actions after Krylov subspace method set up
        postKrylovSetupActions();
    }

    // KSP is now ready
    m_KSPDirty = false;

    // Reset KSP status
    resetKSPStatus();
}

/*!
 * Finalize the KSP after the solution of the system.
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
    m_KSP = PETSC_NULLPTR;
}

/*!
 * Set up the preconditioner.
 */
void SystemSolver::setupPreconditioner()
{
    PC pc;
    KSPGetPC(m_KSP, &pc);
    setupPreconditioner(pc, getKSPOptions());
}

/*!
 * Set up the specified preconditioner using the given options.
 *
 * \param pc is the preconditioner to set up
 * \param options are the options that will be used to set up the preconditioner
 */
void SystemSolver::setupPreconditioner(PC pc, const KSPOptions &options) const
{
    // Set preconditioner type
    PCType pcType;
#if BITPIT_ENABLE_MPI == 1
    if (isPartitioned()) {
        pcType = PCASM;
    } else {
        pcType = PCILU;
    }
#else
    pcType = PCILU;
#endif

    PCSetType(pc, pcType);

    // Configure preconditioner
    if (strcmp(pcType, PCASM) == 0) {
        if (options.overlap != PETSC_DEFAULT) {
            PCASMSetOverlap(pc, options.overlap);
        }
    } else if (strcmp(pcType, PCILU) == 0) {
        if (options.levels != PETSC_DEFAULT) {
            PCFactorSetLevels(pc, options.levels);
        }
    }

    PCSetUp(pc);

    if (strcmp(pcType, PCASM) == 0) {
        KSP *subKSPs;
        PetscInt nSubKSPs;
        PCASMGetSubKSP(pc, &nSubKSPs, PETSC_NULLPTR, &subKSPs);

        for (PetscInt i = 0; i < nSubKSPs; ++i) {
            KSP subKSP = subKSPs[i];
            KSPSetType(subKSP, KSPPREONLY);

            PC subPC;
            KSPGetPC(subKSP, &subPC);
            PCSetType(subPC, PCILU);
            if (options.sublevels != PETSC_DEFAULT) {
                log::warning() << " Setting ASM ILU levels using the member \"sublevels\" is deprecated."
                               << " ASM ILU levels should be set using the member \"levels\"." << std::endl;

                PCFactorSetLevels(subPC, options.sublevels);
            } else if (options.levels != PETSC_DEFAULT) {
                PCFactorSetLevels(subPC, options.levels);
            }

            if (options.subrtol != PETSC_DEFAULT) {
                log::warning() << " The member \"subrtol\" is deprecated. Since ASM is only use a a preconditioner"
                               << " setting the tolerance of the Krylov subspace doens't have any effect on the"
                               << " solution of the system."
                               << std::endl;
            }
        }
    }
}

/*!
 * Perform actions before preconditioner setup.
 */
void SystemSolver::prePreconditionerSetupActions()
{
    // Nothing to do
}

/*!
 * Perform actions after preconditioner setup.
 */
void SystemSolver::postPreconditionerSetupActions()
{
    // Nothing to do
}

/*!
 * Set up the Krylov subspace method used to solve the system.
 */
void SystemSolver::setupKrylov()
{
    setupKrylov(m_KSP, getKSPOptions());
}

/*!
 * Set up the Krylov subspace method using the given options.
 *
 * This function is in charge of only setting the properties of the Krylov subspace
 * method that will be used to solve the system. There is a dedicated function to
 * set up the preconditioner.
 *
 * \param ksp is the KSP whose Krylov subspace method will be setup
 * \param options are the options that will be used to set up the KSP
 */
void SystemSolver::setupKrylov(KSP ksp, const KSPOptions &options) const
{
    KSPSetType(ksp, KSPFGMRES);
    if (options.restart != PETSC_DEFAULT) {
        KSPGMRESSetRestart(ksp, options.restart);
    }
    if (options.rtol != PETSC_DEFAULT || options.atol != PETSC_DEFAULT || options.maxits != PETSC_DEFAULT) {
        KSPSetTolerances(ksp, options.rtol, options.atol, PETSC_DEFAULT, options.maxits);
    }
    KSPSetInitialGuessNonzero(ksp, options.initial_non_zero);
}

/*!
 * Perform actions before Krylov subspace method setup.
 */
void SystemSolver::preKrylovSetupActions()
{
    // Enable convergence monitor
    if (m_convergenceMonitorEnabled) {
#if (PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 7)
        PetscOptionsSetValue(PETSC_NULLPTR, ("-" + m_prefix + "ksp_monitor_true_residual").c_str(), "");
        PetscOptionsSetValue(PETSC_NULLPTR, ("-" + m_prefix + "ksp_monitor_singular_value").c_str(), "");
        PetscOptionsSetValue(PETSC_NULLPTR, ("-" + m_prefix + "ksp_converged_reason").c_str(), "");
#else
        PetscOptionsSetValue(("-" + m_prefix + "ksp_monitor_true_residual").c_str(), "");
        PetscOptionsSetValue(("-" + m_prefix + "ksp_monitor_singular_value").c_str(), "");
        PetscOptionsSetValue(("-" + m_prefix + "ksp_converged_reason").c_str(), "");
#endif
    }
}

/*!
 * Perform actions after Krylov subspace method setup.
 */
void SystemSolver::postKrylovSetupActions()
{
}

/*!
 * Get a reference to the options associated with the Krylov solver.
 *
 * The options associated with the Krylov solver can only be accessed after assembling the system.
 *
 * \return A reference to the options associated with the Krylov solver.
 */
KSPOptions & SystemSolver::getKSPOptions()
{
    return m_KSPOptions;
}

/*!
 * Get a constant reference to the options associated with the Krylov solver.
 *
 * The options associated with the Krylov solver can only be accessed after assembling the system.
 *
 * \return A constant reference to the options associated with the Krylov solver.
 */
const KSPOptions & SystemSolver::getKSPOptions() const
{
    return m_KSPOptions;
}

/*!
 * Initialize the options associated with the KSP.
 */
void SystemSolver::initializeKSPOptions()
{
    resetKSPOptions(&m_KSPOptions);
}

/*!
 * Reset the specified KSP options.
 *
 * \param options are the options that will be reset
 */
void SystemSolver::resetKSPOptions(KSPOptions *options) const
{
    *options = KSPOptions();
}

/*!
 * Destroy the options associated with the KSP.
 */
void SystemSolver::destroyKSPOptions()
{
    resetKSPOptions(&m_KSPOptions);
}

/*!
 * Get a constant reference to the status of the Krylov solver.
 *
 * The status of the Krylov solver can only be accessed after assembling the system.
 *
 * \return A constant reference to the status of the Krylov solver.
 */
const KSPStatus & SystemSolver::getKSPStatus() const
{
    return m_KSPStatus;
}

/*!
 * Initialize the status of the KSP.
 */
void SystemSolver::initializeKSPStatus()
{
    resetKSPStatus();
}

/*!
 * Fill the status of the KSP.
 */
void SystemSolver::fillKSPStatus()
{
    fillKSPStatus(m_KSP, &m_KSPStatus);
}

/*!
 * Fill the status of the specified KSP.
 *
 * \param ksp is the KSP the status will be fetched from
 * \param status on output will contain the status of the KSP
 */
void SystemSolver::fillKSPStatus(KSP ksp, KSPStatus *status) const
{
    status->error = 0;
    KSPGetIterationNumber(ksp, &(status->its));
    KSPGetConvergedReason(ksp, &(status->convergence));
}

/*!
 * Reset the status of the KSP.
 */
void SystemSolver::resetKSPStatus()
{
    resetKSPStatus(&m_KSPStatus);
}

/*!
 * Reset the status of the specified KSP.
 *
 * \param status on output will contain the status of the KSP
 */
void SystemSolver::resetKSPStatus(KSPStatus *status) const
{
    status->error       = 0;
    status->its         = -1;
    status->convergence = KSP_CONVERGED_ITERATING;
}

/*!
 * Destroy the status of the KSP.
 */
void SystemSolver::destroyKSPStatus()
{
    resetKSPStatus();
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

/*!
 * Open an input binary archive.
 *
 * It's up to the caller of this function to close the binary archive when it is not needed
 * anymore.
 *
 * \param directory is the directory that contains the boundary archive
 * \param name is the name, without extension, of the binary archive
 * \param block is the parallel block the archive belongs to, a negative value means that the
 * archive is serial
 * \param[in, out] archive on output will contained the opened binary archive
 */
void SystemSolver::openBinaryArchive(const std::string &directory, const std::string &prefix,
                                     int block, IBinaryArchive *archive) const
{
    int expectedVersion = getDumpVersion();

    std::string path = getFilePath(directory, prefix);
    archive->open(path, "dat", block);
    bool versionsMatch = archive->checkVersion(expectedVersion);
    if (!versionsMatch) {
        std::string message = "Version stored in binary archive " + path + " does not match";
        message += " the expected version " + std::to_string(expectedVersion) + ".";

        throw std::runtime_error(message);
    }
}

/*!
 * Open an output binary archive.
 *
 * It's up to the caller of this function to close the binary archive when it is not needed
 * anymore.
 *
 * \param header is the header of the archive
 * \param directory is the directory that will contain the boundary archive
 * \param name is the name, without extension, of the binary archive
 * \param block is the parallel block the archive belongs to, a negative value means that the
 * archive is serial
 * \param[in, out] archive on output will contained the opened binary archive
 */
void SystemSolver::openBinaryArchive(const std::string &header, const std::string &directory,
                                     const std::string &prefix, int block, OBinaryArchive *archive) const
{
    int version = getDumpVersion();
    std::string path = getFilePath(directory, prefix);
    archive->open(path, "dat", version, header, block);
}

/*!
 * Close the specified binary archive.
 *
 * \param[in, out] archive is the archive that will be closed
 */
void SystemSolver::closeBinaryArchive(BinaryArchive *archive) const
{
    archive->close();
}

/*!
 * Get the block of the binary archive assigned to this process.
 *
 * \result The block of the binary archive assigned to this process.
 */
int SystemSolver::getBinaryArchiveBlock() const
{
#if BITPIT_ENABLE_MPI == 1
    if (isPartitioned()) {
        int nProcesses;
        MPI_Comm_size(getCommunicator(), &nProcesses);
        if (nProcesses <= 1) {
            return -1;
        }

        int archiveBlock;
        MPI_Comm_rank(getCommunicator(), &archiveBlock);

        return archiveBlock;
    } else {
        return -1;
    }
#else
    return -1;
#endif
}

}
