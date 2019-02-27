/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include <bitpit_IO.hpp>

#include "system_solvers_large.hpp"

namespace bitpit {

int SystemSolver::m_nInstances = 0;
std::vector<std::string> SystemSolver::m_options = std::vector<std::string>(1, "bitpit");

/*!
 * Set initialization option
 */
void SystemSolver::addInitOption(const std::string &option)
{
    if (m_nInstances != 0) {
        throw std::runtime_error("Initialization opions can be set only before initializing the solver.");
    }

    m_options.push_back(option);
}

/*!
 * Set initialization options
 */
void SystemSolver::addInitOptions(const std::vector<std::string> &options)
{
    if (m_nInstances != 0) {
        throw std::runtime_error("Initialization opions can be set only before initializing the solver.");
    }

    // The first option is the executable name and we set it to a dummy value.
    for (std::string option : options) {
        m_options.push_back(option);
    }
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Default constuctor
 *
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool debug)
    : SystemSolver(MPI_COMM_SELF, debug)
{
}

/*!
 * Constuctor
 *
 * \param communicator is the MPI communicator
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(MPI_Comm communicator, bool debug)
#else
/*!
 * Defualt constuctor
 *
 * \param debug if set to true, debug information will be printed
 */
SystemSolver::SystemSolver(bool debug)
#endif
    : m_initialized(false), m_pivotType(PIVOT_NONE),
      m_rowGlobalOffset(0), m_colGlobalOffset(0),
      m_A(nullptr), m_rhs(nullptr), m_solution(nullptr),
      m_rpivot(nullptr), m_cpivot(nullptr), m_KSP(nullptr)
{
    // Add debug options
    if (debug) {
        addInitOption("-log_view");
        addInitOption("-ksp_monitor_true_residual");
        addInitOption("-ksp_converged_reason");
    }

    // Initialize Petsc
    if (m_nInstances == 0) {
        const char help[] = "None";

        int argc = m_options.size();
        char **argv = new char*[argc];
        for (int i = 0; i < argc; ++i) {
            argv[i] = (char*) m_options[i].c_str();
        }

        PetscInitialize(&argc, &argv, 0, help);

        delete[] argv;
    }

    // Increase the number of instances
    ++m_nInstances;

#if BITPIT_ENABLE_MPI==1
    // Detect if the system is partitioned
    m_partitioned = (communicator != MPI_COMM_NULL) && (communicator != MPI_COMM_SELF);

    // Set the communicator
    if (m_partitioned) {
        MPI_Comm_dup(communicator, &m_communicator);
    } else {
        m_communicator = PETSC_COMM_SELF;
    }
#endif
}

/*!
 * Destructor
 */
SystemSolver::~SystemSolver()
{
    // Decrease the number of instances
    --m_nInstances;

#if BITPIT_ENABLE_MPI==1
    // Free the MPI communicator
    if (m_partitioned) {
        int finalizedCalled;
        MPI_Finalized(&finalizedCalled);
        if (!finalizedCalled) {
            MPI_Comm_free(&m_communicator);
        }
    }
#endif

    // Finalize petsc
    if (m_nInstances == 0) {
        PetscFinalize();
    }
}

/*!
 * Clear the system
 */
void SystemSolver::clear()
{
    if (!m_initialized) {
        return;
    }

    MatDestroy(&m_A);
    VecDestroy(&m_rhs);
    VecDestroy(&m_solution);

    KSPDestroy(&m_KSP);

    if (getPivotType() != PIVOT_NONE) {
        ISDestroy(&m_rpivot);
        ISDestroy(&m_cpivot);
    }

    m_initialized = false;
}

/*!
 * Initialize the system.
 *
 * \param matrix is the matrix
 * \param pivotType is the type of pivoting that will be used
 */
void SystemSolver::initialize(const SparseMatrix &matrix, PivotType pivotType)
{
    // Check if the matrix is finalized
    if (!matrix.isFinalized()) {
        throw std::runtime_error("Unable to initialize the system. The matrix is not yet finalized.");
    }

    // Clear the system
    clear();

    // Initialize matrix
    matrixInit(matrix);
    matrixFill(matrix);

    // Initialize pivot
    pivotInit(pivotType);
    if (getPivotType() != PIVOT_NONE) {
        matrixReorder();
    }

    // Initialize RHS and solution vectors
#if BITPIT_ENABLE_MPI == 1
    vectorsInit(matrix.extractGhostGlobalCols());
#else
    vectorsInit();
#endif

    // Initialize Krylov solver
    KSPInit();

    // The system is now initialized
    m_initialized = true;
}

/*!
 * Update the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * initialized its pattern cannot be modified.
 *
 * \param rows are the global indices of the rows that will be updated
 * \param elements are the elements that will be used to update the rows
 */
void SystemSolver::update(const std::vector<long> &rows, const SparseMatrix &elements)
{
    // Check if the element storage is finalized
    if (!elements.isFinalized()) {
        throw std::runtime_error("Unable to update the system. The element storage is not yet finalized.");
    }

    // Check if the system is initialized
    if (!m_initialized) {
        throw std::runtime_error("Unable to update the system. The system is not yet initialized.");
    }

    // Update matrix
    matrixUpdate(rows, elements);
}

/**
* Get the number of rows of the system.
*
* \result The number of rows of the system.
*/
long SystemSolver::getRowCount() const
{
    if (!isInitialized()) {
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
    if (!isInitialized()) {
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
    if (!isInitialized()) {
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
    if (!isInitialized()) {
        return 0;
    }

    PetscInt nCols;
    MatGetSize(m_A, NULL, &nCols);

    return nCols;
}
#endif

/*!
 * Check if the system is initialized.
 *
 * \return Returns true if the system is initialized, false otherwise.
 */
bool SystemSolver::isInitialized() const
{
    return m_initialized;
}

/*!
 * Solve the system
 */
void SystemSolver::solve()
{
    // Reorder the vectors
    if (getPivotType() != PIVOT_NONE) {
        vectorsReorder(PETSC_FALSE);
    }

    // Solve the system
    m_KSPStatus.error = KSPSolve(m_KSP, m_rhs, m_solution);

    // Set solver info
    if (m_KSPStatus.error == 0) {
        KSPGetIterationNumber(m_KSP, &m_KSPStatus.its);
        KSPGetConvergedReason(m_KSP, &m_KSPStatus.convergence);
    } else {
        m_KSPStatus.its         = -1;
        m_KSPStatus.convergence = KSP_DIVERGED_BREAKDOWN;
    }

    // Reorder the vectors
    if (getPivotType() != PIVOT_NONE) {
        vectorsReorder(PETSC_TRUE);
    }
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
 * Initializes the matrix.
 *
 * \param matrix is the matrix
 */
void SystemSolver::matrixInit(const SparseMatrix &matrix)
{
    long nRows = matrix.getRowCount();
    long nCols = matrix.getColCount();

    // Set row and column global offset
#if BITPIT_ENABLE_MPI == 1
    m_rowGlobalOffset = matrix.getRowGlobalOffset();
    m_colGlobalOffset = matrix.getColGlobalOffset();
#else
    m_rowGlobalOffset = 0;
    m_colGlobalOffset = 0;
#endif

#if BITPIT_ENABLE_MPI == 1
    // Evaluate the number of non-zero elements
    //
    // For each row we count the number of local non-zero elements (d_nnz) and
    // the number of non-zero elements that belong to other processors (o_nnz)
    long nGlobalCols      = matrix.getColGlobalCount();
    long nOffDiagonalCols = nGlobalCols - nCols;
    long firstColGlobalId = matrix.getColGlobalOffset();
    long lastColGlobalId  = firstColGlobalId + nCols - 1;

    std::vector<int> d_nnz(nRows, 0);
    std::vector<int> o_nnz(nRows, 0);
    if (nOffDiagonalCols > 0) {
        for (long row = 0; row < nRows; ++row) {
            ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
            int nRowNZ = rowPattern.size();
            for (int k = 0; k < nRowNZ; ++k) {
                long columnGlobalId = rowPattern[k];
                if (columnGlobalId < firstColGlobalId || columnGlobalId > lastColGlobalId) {
                    ++o_nnz[row];
                } else {
                    ++d_nnz[row];
                }
            }
        }
    } else {
        for (long row = 0; row < nRows; ++row) {
            ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
            d_nnz[row] = rowPattern.size();
        }
    }

    // Create the matrix
    MatCreateAIJ(m_communicator, nRows, nCols, PETSC_DETERMINE, PETSC_DETERMINE, 0, d_nnz.data(), 0, o_nnz.data(), &m_A);
#else
    // Evaluate the number of non-zero elements
    std::vector<int> d_nnz(nCols);

    for (long row = 0; row < nRows; ++row) {
        ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
        d_nnz[row] = rowPattern.size();
    }

    // Create the matrix
    MatCreateSeqAIJ(PETSC_COMM_SELF, nRows, nCols, 0, d_nnz.data(), &m_A);
#endif
}

/*!
 * Fills the matrix.
 *
 * \param matrix is the matrix
 */
void SystemSolver::matrixFill(const SparseMatrix &matrix)
{
    const long nRows = matrix.getRowCount();
    const long maxRowNZ = matrix.getMaxRowNZCount();

    // Create the matrix
    std::vector<PetscInt> rowNZGlobalIds(maxRowNZ);
    std::vector<PetscScalar> rowNZValues(maxRowNZ);

    for (long row = 0; row < nRows; ++row) {
        ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
        ConstProxyVector<double> rowValues = matrix.getRowValues(row);

        const int nRowNZ = rowPattern.size();
        const PetscInt globalRow = m_rowGlobalOffset + row;
        for (int k = 0; k < nRowNZ; ++k) {
            rowNZGlobalIds[k] = rowPattern[k];
            rowNZValues[k] = rowValues[k];
        }

        MatSetValues(m_A, 1, &globalRow, nRowNZ, rowNZGlobalIds.data(), rowNZValues.data(), INSERT_VALUES);
    }

    // Let petsc build the matrix
    MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
}

/*!
 * Update the specified rows of the matrix.
 *
 * The contents of the specified rows will be replaced by the specified
 * elements.
 *
 * \param rows are the global indices of the rows that will be updated
 * \param elements are the elements that will be used to update the rows
 */
void SystemSolver::matrixUpdate(const std::vector<long> &rows, const SparseMatrix &elements)
{
    const long maxRowElements = elements.getMaxRowNZCount();
    std::vector<std::size_t> elementsColsIndexes(maxRowElements);

    std::map<PetscInt, std::size_t> rowGlobalColsMap;

    std::vector<PetscScalar> rowValues(maxRowElements);

    for (std::size_t n = 0; n < rows.size(); ++n) {
        ConstProxyVector<double> elementsValues = elements.getRowValues(n);
        const int nRowElements = elementsValues.size();
        if (nRowElements == 0) {
            continue;
        }

        // Get global row
        long row = rows[n];
        const PetscInt globalRow = m_rowGlobalOffset + row;

        // Get current row elements
        PetscInt *nRowCols = nullptr;
        const PetscInt **rowGlobalCols = nullptr;
        MatGetRow(m_A, globalRow, nRowCols, rowGlobalCols, NULL);
        assert(nRowCols != nullptr);
        assert(rowGlobalCols != nullptr);

        // Get elements indices
        rowGlobalColsMap.clear();
        for (PetscInt k = 0; k < *nRowCols; ++k) {
            PetscInt globalCol = (*rowGlobalCols)[n];
            rowGlobalColsMap[globalCol] = k;
        }

        ConstProxyVector<long> elementsGlobalCols = elements.getRowPattern(row);
        for (int k = 0; k < nRowElements; ++k) {
            long elementGlobalCol = elementsGlobalCols[k];

            auto rowGlobalColsMapItr = rowGlobalColsMap.find(elementGlobalCol);
            if (rowGlobalColsMapItr == rowGlobalColsMap.end()) {
                throw std::runtime_error("The element is not in the matrix.");
            }

            elementsColsIndexes[k] = rowGlobalColsMapItr->second;
        }

        // Update values
        for (int k = 0; k < *nRowCols; ++k) {
            rowValues[k] = 0.;
        }

        for (int k = 0; k < nRowElements; ++k) {
            long elementColsIndex = elementsColsIndexes[k];
            rowValues[elementColsIndex] += elementsValues[k];
        }

        // Restore row
        MatRestoreRow(m_A, globalRow, nRowCols, rowGlobalCols, NULL);

        // update values
        MatSetValuesRow(m_A, globalRow, rowValues.data());
    }
}

/*!
 * Reorder the matrix.
 */
void SystemSolver::matrixReorder()
{
    Mat B;
    MatPermute(m_A, m_rpivot, m_cpivot, &B);
    MatDestroy(&m_A);
    MatDuplicate(B, MAT_COPY_VALUES, &m_A);
    MatDestroy(&B);
}

/*!
 * Initialize rhs and solution vectors.
 */
#if BITPIT_ENABLE_MPI == 1
/*!
 * \param ghosts is the list of global ids that are ghosts for the local
 * processor
 */
void SystemSolver::vectorsInit(const std::vector<long> &ghosts)
#else
void SystemSolver::vectorsInit()
#endif
{
    PetscInt nRows;
    PetscInt nColumns;
    MatGetLocalSize(m_A, &nRows, &nColumns);

#if BITPIT_ENABLE_MPI == 1
    PetscInt nGlobalRows;
    PetscInt nGlobalColumns;
    MatGetSize(m_A, &nGlobalRows, &nGlobalColumns);

    int nGhosts = 0;
    std::vector<PetscInt> ghostGlobalIds(ghosts.size());
    for (auto &ghostGlobalId : ghosts) {
        ghostGlobalIds[nGhosts] = ghostGlobalId;
        ++nGhosts;
    }

    VecCreateGhost(m_communicator, nColumns, nGlobalColumns, nGhosts, ghostGlobalIds.data(), &m_solution);
    VecCreateGhost(m_communicator, nRows, nGlobalRows, nGhosts, ghostGlobalIds.data(), &m_rhs);
#else
    VecCreateSeq(PETSC_COMM_SELF, nColumns, &m_solution);
    VecCreateSeq(PETSC_COMM_SELF, nRows, &m_rhs);
#endif
}

/*!
 * Reorders rhs and solution vectors.
 */
void SystemSolver::vectorsReorder(PetscBool inv)
{
    VecPermute(m_solution, m_cpivot, inv);
    VecPermute(m_rhs, m_rpivot, inv);
}

/*!
 * Fills rhs and solution vectors.
 */
void SystemSolver::vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Import RHS
    int nRows;
    VecGetLocalSize(m_rhs, &nRows);

    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);
    for (int i = 0; i < nRows; ++i) {
        raw_rhs[i] = rhs[i];
    }
    VecRestoreArray(m_rhs, &raw_rhs);

    // Import initial solution
    int nColumns;
    VecGetLocalSize(m_solution, &nColumns);

    PetscScalar *raw_solution;
    VecGetArray(m_solution, &raw_solution);
    for (int i = 0; i < nColumns; ++i) {
        raw_solution[i] = (*solution)[i];
    }
    VecRestoreArray(m_solution, &raw_solution);
}

/*!
 * Export the solution vector.
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
 * \param raw_solution is the location of pointer to array obtained from
 * getRHSRawPtr()
 */
void SystemSolver::restoreRHSRawPtr(double *raw_rhs)
{
    VecRestoreArray(m_rhs, &raw_rhs);
}

/*!
 * Restores the solution vector after getRHSRawReadPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
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
 */
void SystemSolver::dump(const std::string &directory, const std::string &prefix) const
{
    std::stringstream filePathStream;

    // Matrix
    PetscViewer matViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &matViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &matViewer);
#endif
    PetscViewerSetType(matViewer, PETSCVIEWERASCII);
    PetscViewerPushFormat(matViewer, PETSC_VIEWER_DEFAULT);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "A.txt";
    PetscViewerFileSetName(matViewer, filePathStream.str().c_str());
    MatView(m_A, matViewer);
    PetscViewerDestroy(&matViewer);

    // RHS
    PetscViewer rhsViewer;
#if BITPIT_ENABLE_MPI==1
    PetscViewerCreate(m_communicator, &rhsViewer);
#else
    PetscViewerCreate(PETSC_COMM_SELF, &rhsViewer);
#endif
    PetscViewerSetType(rhsViewer, PETSCVIEWERASCII);
    PetscViewerPushFormat(rhsViewer, PETSC_VIEWER_DEFAULT);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "rhs.txt";
    PetscViewerFileSetName(rhsViewer, filePathStream.str().c_str());
    VecView(m_rhs, rhsViewer);
    PetscViewerDestroy(&rhsViewer);
}

/*!
 * Initialize pivot
 *
 * \param pivotType is the type of pivoting that will be used
 */
void SystemSolver::pivotInit(PivotType pivotType)
{
    m_pivotType = pivotType;

    MatOrderingType petscPivotType;
    switch (m_pivotType) {

    case (PIVOT_ND):
        petscPivotType = MATORDERINGNATURAL;
        break;

    case (PIVOT_1WD):
        petscPivotType = MATORDERING1WD;
        break;

    case (PIVOT_RCM):
        petscPivotType = MATORDERINGRCM;
        break;

    case (PIVOT_MD):
        petscPivotType = MATORDERINGQMD;
        break;

    default:
        return;

    }

    MatGetOrdering(m_A, petscPivotType, &m_rpivot, &m_cpivot);
    ISSetPermutation(m_rpivot);
    ISSetPermutation(m_cpivot);
}

/*!
 * Get the pivot type.
 *
 * \result The pivot type.
 */
SystemSolver::PivotType SystemSolver::getPivotType()
{
    return m_pivotType;
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
 * Initialize the Krylov solver.
 */
void SystemSolver::KSPInit()
{
#if BITPIT_ENABLE_MPI==1
    KSPCreate(m_communicator, &m_KSP);
#else
    KSPCreate(PETSC_COMM_SELF, &m_KSP);
#endif
    KSPSetOperators(m_KSP, m_A, m_A);

    KSPSetType(m_KSP, KSPFGMRES);
    if (m_KSPOptions.restart != PETSC_DEFAULT) {
        KSPGMRESSetRestart(m_KSP, m_KSPOptions.restart);
    }
    if (m_KSPOptions.rtol != PETSC_DEFAULT || m_KSPOptions.maxits != PETSC_DEFAULT) {
        KSPSetTolerances(m_KSP, m_KSPOptions.rtol, PETSC_DEFAULT, PETSC_DEFAULT, m_KSPOptions.maxits);
    }
    KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);

    PCType preconditionerType;
#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
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
    KSPSetFromOptions(m_KSP);
    KSPSetUp(m_KSP);

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

}
