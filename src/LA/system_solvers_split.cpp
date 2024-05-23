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

#include "matrix_utilities.hpp"
#include "system_solvers_split.hpp"

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_operators.hpp"

#include "petscksp.h"
#include "petscmat.h"
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
 * \class SplitSystemMatrixAssembler
 * \ingroup system_solver_split
 *
 * \brief The SplitSystemMatrixAssembler class is the base class for defining assemblers for
 * split system solvers.
 */

/*!
 * Constructor.
 *
 * \param splitStrategy the type of split that will be applied to the system
 * \param splitSizes are the sizes of the splits
 */
SplitSystemMatrixAssembler::SplitSystemMatrixAssembler(SystemSplitStrategy splitStrategy, const std::vector<int> &splitSizes)
    : SystemMatrixAssembler(),
      m_splitStrategy(splitStrategy), m_splitSizes(splitSizes)
{
}

/*!
 * Get the split strategy will be used when solving the system.
 *
 * \result The split strategy will be used when solving the system.
 */
SystemSplitStrategy SplitSystemMatrixAssembler::getSplitStrategy() const
{
    return m_splitStrategy;
}

/*!
 * Get the number of splits that will be applied.
 *
 * \result The number of splits that will be applied.
 */
int SplitSystemMatrixAssembler::getSplitCount() const
{
    return m_splitSizes.size();
}

/*!
 * Get the size of each split.
 *
 * Each block of elements is split into pieces whose size is defined by the split sizes. The
 * sum of the split sizes must be equal to the block size. Given a matrix with a block size
 * equal to 5 and split sizes equal to [4, 1], the split matrix will look like this:
 *
 *        [                     |                     ]
 *        [  (4 * N) * (4 * M)  |  (4 * N) * (1 * M)  ]
 *        [                     |                     ]
 *        [-------------------------------------------]
 *        [                     |                     ]
 *        [  (1 * N) * (4 * M)  |  (1 * N) * (1 * M)  ]
 *        [                     |                     ]
 *
 * where M and N are the number of rows and columns respectively.
 *
 * \result The size of each split.
 */
const std::vector<int> & SplitSystemMatrixAssembler::getSplitSizes() const
{
    return m_splitSizes;
}

/*!
 * \class SplitSystemSparseMatrixAssembler
 * \ingroup system_solver_split
 *
 * \brief The SplitSystemSparseMatrixAssembler class allows to assembly a split system solver
 * from a sparse matrix.
 */

/*!
 * Constructor.
 *
 * \param matrix is the matrix
 * \param splitStrategy the type of split that will be applied to the system
 * \param splitSizes are the sizes of the splits
 */
SplitSystemSparseMatrixAssembler::SplitSystemSparseMatrixAssembler(const SparseMatrix *matrix, SystemSplitStrategy splitStrategy,
                                                                   const std::vector<int> &splitSizes)
    : SystemSparseMatrixAssembler(matrix),
      SplitSystemMatrixAssembler(splitStrategy, splitSizes)
{
}

/*!
 * \class SplitSystemSolver
 * \ingroup system_solver_split
 *
 * \brief The SplitSystemSolver class allows to solve a split linear system.
 *
 * Block matrices represent an important class of problems in numerical linear algebra and offer
 * the possibility of far more efficient iterative solvers than just treating the entire matrix
 * as black box. Following the common linear algebra definition, a block matrix is a matrix
 * divided in a small, problem-size independent (two, three or so) number of very large blocks.
 * These blocks arise naturally from the underlying physics or discretization of the problem,
 * for example, they may be associated with different variables of a physical problem. Under a
 * certain numbering of unknowns the matrix can be written as
 *
 *                          [ A00  A01  A02 ]
 *                          [ A10  A11  A12 ]
 *                          [ A20  A21  A22 ]
 *
 * where each A_ij is an entire block (see the paragraph "Solving Block Matrices" in the PETSc
 * manual).
 *
 * When assembling the matrix, a monolithic matrix should be provided. For example, assuming to
 * group the elements of the matrix in five-by-five groups (here, five is the so-called block
 * size of the matrix and usually rises when coupling variables with different meaning, for
 * example pressure, the three components of the velocity and temperature) the assembler will
 * provide the following system:
 *
 *     [           |           |          |           ]   [           ]     [           ]
 *     [  (5 x 5)  |  (5 x 5)  |    ...   |  (5 x 5)  ]   [  (5 x 5)  ]     [  (5 x 5)  ]
 *     [           |           |          |           ]   [           ]     [           ]
 *     [----------------------------------------------]   [-----------]     [-----------]
 *     [           |           |          |           ]   [           ]     [           ]
 *     [    ...    |    ...    |    ...   |  (5 x 5)  ]   [    ...    ]  =  [    ...    ]
 *     [           |           |          |           ]   [           ]     [           ]
 *     [----------------------------------------------]   [-----------]     [-----------]
 *     [           |           |          |           ]   [           ]     [           ]
 *     [  (5 x 5)  |  (5 x 5)  |    ...   |  (5 x 5)  ]   [  (5 x 5)  ]     [  (5 x 5)  ]
 *     [           |           |          |           ]   [           ]     [           ]
 *
 * Internally, the monolithic matrix will be split into blocks. For example, considering two
 * splits, the first one that group together the first four variables and the second one that
 * holds the last variable (i.e., split sizes equal to [4, 1]), the internal split system
 * will be:
 *
 *       [                     |                     ]  [           ]     [           ]
 *       [  (4 * N) x (4 * M)  |  (4 * N) x (1 * M)  ]  [  (4 x M)  ]     [  (4 x N)  ]
 *       [                     |                     ]  [           ]     [           ]
 *       [-------------------------------------------]  [-----------]  =  [-----------]
 *       [                     |                     ]  [           ]     [           ]
 *       [  (1 * N) x (4 * M)  |  (1 * N) x (1 * M)  ]  [  (1 x M)  ]     [  (1 x N)  ]
 *       [                     |                     ]  [           ]     [           ]
 *
 * where M and N are the number of rows and columns respectively.
 *
 * The PETSc PCFIELDSPLIT preconditioner is used to solve the split system. There are different
 * split strategies available:
 *
 * 1) DIAGONAL: this strategy assumes that the only non-zero blocks are the diagonal ones.
 *
 * Considering a two-by-two block block matrix
 *
 *                            [ A00    0  ]
 *                            [  0    A11 ],
 *
 * the preconditioned problem will look like
 *
 *           [ KSPSolve(A00)        0         ]  [ A00    0  ]
 *           [      0           KSPSolve(A00) ]  [  0    A11 ],
 *
 * in other words the preconditioner is:
 *
 *                                  ( [ A00    0  ] )
 *              approximate inverse (               )
 *                                  ( [  0    A11 ] ).
 *
 * The system is solved efficiently by solving each block independently from the others.
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * Blocks are solved using a flexible GMRES iterative method. If the system is partitioned
 * each block is preconditioned using the (restricted) additive Schwarz method (ASM). On
 * each block of the ASM preconditioner an incomplete LU factorization (ILU) is used. There
 * is one ASM block per process. If the system is not partitioned it is preconditioned using
 * the incomplete LU factorization (ILU).
 */
#else
/*!
 * Blocks are solved using a flexible GMRES iterative method preconditioned using the incomplete
 * LU factorization (ILU).
 */
#endif
/*!
 * 2) LOWER: this strategy assumes that the only non-zero blocks are the ones on an below
 * the diagonal.
 *
 * Considering a two-by-two block block matrix
 *
 *                            [ A00    0  ]
 *                            [ A01   A11 ],
 *
 * the preconditioner is
 *
 *                                  ( [ A00    0  ] )
 *              approximate inverse (               )
 *                                  ( [ A01   A11 ] ).
 *
 * The system is solved efficiently by first solving with A00, then applying A01 to that
 * solution, removing it from the right hand side of the second block and then solving with
 * A11, in other words
 *
 *               [ I      0   ]  [    I    0 ]  [  A00^-1  0 ]
 *               [ 0   A11^-1 ]  [ - A10   I ]  [   0      I ].
 *
 * This strategy can be seen as a "block" Gauss-Seidel with the blocks being the splits.
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * Blocks are solved using a flexible GMRES iterative method. If the system is partitioned
 * each block is preconditioned using the (restricted) additive Schwarz method (ASM). On
 * each block of the ASM preconditioner an incomplete LU factorization (ILU) is used. There
 * is one ASM block per process. If the system is not partitioned it is preconditioned using
 * the incomplete LU factorization (ILU).
 */
#else
/*!
 * Blocks are solved using a flexible GMRES iterative method preconditioned using the incomplete
 * LU factorization (ILU).
 */
#endif
/*!
 * 3) FULL: this strategy doesn't make assumptions on the structure of the blocks.
 *
 * Considering a two-by-two block block matrix
 *
 *                            [ A00   A01 ]
 *                            [ A11   A11 ],
 *
 * the preconditioned problem will look like
 *
 *           [ KSPSolve(A00)        0         ]  [ A00   A01 ]
 *           [      0           KSPSolve(A00) ]  [ A11   A11 ],
 *
 * in other words the preconditioner is:
 *
 *                                  ( [ A00    0  ] )
 *              approximate inverse (               )
 *                                  ( [  0    A11 ] )
 *
 * The preconditioner is evaluated considering only the diagonal blocks and then the full
 * system is solved.
 */
#if BITPIT_ENABLE_MPI==1
/*!
 * The system is solved using a flexible GMRES iterative method. If the system is partitioned
 * each diagonal block is preconditioned using the (restricted) additive Schwarz method (ASM).
 * On each block of the ASM preconditioner an incomplete LU factorization (ILU) is used. There
 * is one ASM block per process. If the system is not partitioned it is preconditioned using
 * the incomplete LU factorization (ILU).
 */
#else
/*!
 * The system is solved using a flexible GMRES iterative method and each diagonla block is
 * preconditioned using the incomplete LU factorization (ILU).
 */
#endif

/*!
 * Constructor.
 *
 * \param debug if set to true, debug information will be printed
 */
SplitSystemSolver::SplitSystemSolver(bool debug)
    : SplitSystemSolver("", false, false, debug)
{
}

/*!
 * Constuctor
 *
 * \param transpose if set to true, transposed system will be solved
 * \param debug if set to true, debug information will be printed
 */
SplitSystemSolver::SplitSystemSolver(bool transpose, bool debug)
    : SplitSystemSolver("", false, transpose, debug)
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
SplitSystemSolver::SplitSystemSolver(bool flatten, bool transpose, bool debug)
    : SplitSystemSolver("", flatten, transpose, debug)
{
}

/*!
 * Constructor.
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param debug if set to true, debug information will be printed
 */
SplitSystemSolver::SplitSystemSolver(const std::string &prefix, bool debug)
    : SplitSystemSolver(prefix, false, debug)
{
}

/*!
 * Constructor.
 *
 * \param prefix is the prefix string to prepend to all option requests
 * \param debug if set to true, debug information will be printed
 */
SplitSystemSolver::SplitSystemSolver(const std::string &prefix, bool transpose, bool debug)
    : SplitSystemSolver(prefix, false, transpose, debug)
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
SplitSystemSolver::SplitSystemSolver(const std::string &prefix, bool flatten, bool transpose, bool debug)
    : SystemSolver(prefix, flatten, transpose, debug),
      m_rhsSplitPermutation(PETSC_NULLPTR), m_solutionSplitPermutation(PETSC_NULLPTR)
{
}

/*!
 * Get the block size of the system.
 *
 * \result The block size of the system.
 */
int SplitSystemSolver::getBlockSize() const
{
    if (m_A == PETSC_NULLPTR) {
        return 0;
    }

    std::vector<int> splitBlockSizes = getSplitSizes();
    int blockSize = std::accumulate(splitBlockSizes.begin(), splitBlockSizes.end(), 0);

    return blockSize;
}

/*!
 * Get the split strategy will be used when solving the system.
 *
 * \result The split strategy will be used when solving the system.
 */
SystemSplitStrategy SplitSystemSolver::getSplitStrategy() const
{
    return m_splitStrategy;
}

/*!
 * Get the number of splits that will be applied.
 *
 * \result The number of splits that will be applied.
 */
int SplitSystemSolver::getSplitCount() const
{
    if (m_A == PETSC_NULLPTR) {
        return -1;
    }

    PetscInt nNestedRows;
    MatNestGetSize(m_A, &nNestedRows, PETSC_NULLPTR);

    return nNestedRows;
}

/*!
 * Get the size of each split.
 *
 * Each block of elements is split into pieces whose size is defined by the split sizes. The
 * sum of the split sizes must be equal to the block size. Given a matrix with a block size
 * equal to 5 and split sizes equal to [4, 1], the split matrix will look like this:
 *
 *        [                     |                     ]
 *        [  (4 * N) * (4 * M)  |  (4 * N) * (1 * M)  ]
 *        [                     |                     ]
 *        [-------------------------------------------]
 *        [                     |                     ]
 *        [  (1 * N) * (4 * M)  |  (1 * N) * (1 * M)  ]
 *        [                     |                     ]
 *
 * where M and N are the number of rows and columns respectively.
 *
 * \result The size of each split.
 */
std::vector<int> SplitSystemSolver::getSplitSizes() const
{
    if (m_A == PETSC_NULLPTR) {
        return std::vector<int>();
    }

    int nSplits = getSplitCount();
    std::vector<int> splitSizes(nSplits, 0);
    for (int k = 0; k < nSplits; ++k) {
        int kk = getBlockSplitLinearIndex(k, k, nSplits);

        PetscInt splitBlockSize;
        MatGetBlockSize(m_splitAs[kk], &splitBlockSize);
        splitSizes[k] = splitBlockSize;
    }

    return splitSizes;
}

/*!
 * Get the block offtest of the splits.
 *
 * \result The block offsets of the splits.
 */
std::vector<int> SplitSystemSolver::getSplitOffsets() const
{
    std::vector<int> splitBlockSizes = getSplitSizes();
    std::vector<int> splitBlockOffsets(splitBlockSizes.size());
    splitBlockOffsets[0] = 0;
    std::partial_sum(splitBlockSizes.begin(), std::prev(splitBlockSizes.end()), std::next(splitBlockOffsets.begin()));

    return splitBlockOffsets;
}

/*!
 * Post-solve actions.
 */
void SplitSystemSolver::postKSPSolveActions()
{
    // Base actions
    SystemSolver::postKSPSolveActions();

    // Fill statuses of split KSP
    fillSplitKSPStatuses();
}

/*!
 * Get a reference to the options associated to the Krylov solver of the specified split.
 *
 * The options associated with the Krylov solver can only be accessed after assembling the system.
 *
 * \param split is the split
 * \return A reference to the options associated to the Krylov solver of the specified split.
 */
KSPOptions & SplitSystemSolver::getSplitKSPOptions(int split)
{
    if (!isAssembled()) {
        throw std::runtime_error("The options associated with the Krylov solver can only be accessed after assembling the system.");
    }

    return m_splitKSPOptions[split];
}

/*!
 * Get a constant reference to the options associated to the Krylov solver of the specified split.
 *
 * The options associated with the Krylov solver can only be accessed after assembling the system.
 *
 * \param split is the split
 * \return A constant reference to the options associated to the Krylov solver of the specified
 * split.
 */
const KSPOptions & SplitSystemSolver::getSplitKSPOptions(int split) const
{
    if (!isAssembled()) {
        throw std::runtime_error("The options associated with the Krylov solver can only be accessed after assembling the system.");
    }

    return m_splitKSPOptions[split];
}

/*!
 * Initialize the options associated with the KSP.
 */
void SplitSystemSolver::initializeKSPOptions()
{
    SystemSolver::initializeKSPOptions();
    initializeSplitKSPOptions();
}

/*!
 * Initialize the options associated with the sub KSP.
 */
void SplitSystemSolver::initializeSplitKSPOptions()
{
    int nSplits = getSplitCount();
    m_splitKSPOptions.resize(nSplits);
    for (KSPOptions &splitOptions : m_splitKSPOptions) {
        resetKSPOptions(&splitOptions);
    }
}

/*!
 * Destroy the options associated with the KSP.
 */
void SplitSystemSolver::destroyKSPOptions()
{
    SystemSolver::destroyKSPOptions();
    destroySplitKSPOptions();
}

/*!
 * Destroy the options associated with the split KSPs.
 */
void SplitSystemSolver::destroySplitKSPOptions()
{
    m_splitKSPOptions.clear();
}

/*!
 * Get a constant reference to the status of the Krylov solver of the specified split.
 *
 * The status of the Krylov solver can only be accessed after assembling the system.
 *
 * \param split is the split
 * \return A constant reference to the status of the Krylov solver of the specified split.
 */
const KSPStatus & SplitSystemSolver::getSplitKSPStatus(int split) const
{
    if (!isAssembled()) {
        throw std::runtime_error("The status of the Krylov solver can only be accessed after assembling the system.");
    }

    return m_splitKSPStatuses[split];
}
/*!
 * Initialize the status of the split KSPs.
 */
void SplitSystemSolver::initializeKSPStatus()
{
    SystemSolver::initializeKSPStatus();
    initializeSplitKSPStatuses();
}

/*!
 * Initialize the options associated with the split KSPs.
 */
void SplitSystemSolver::initializeSplitKSPStatuses()
{
    int nSplits = getSplitCount();
    for (int k = static_cast<int>(m_splitKSPStatuses.size()); k < nSplits; ++k) {
        m_splitKSPStatuses.emplace_back();
        KSPStatus *splitStatus = &(m_splitKSPStatuses.back());
        resetKSPStatus(splitStatus);
    }
    m_splitKSPStatuses.resize(nSplits);
}

/*!
 * Fill the status of the KSP.
 */
void SplitSystemSolver::fillKSPStatus()
{
    SystemSolver::fillKSPStatus();
    fillSplitKSPStatuses();
}

/*!
 * Fill the status of the split KSPs.
 */
void SplitSystemSolver::fillSplitKSPStatuses()
{
    PC pc;
    KSPGetPC(m_KSP, &pc);

    PetscInt nSplits;
    KSP *splitKSPs;
    PCFieldSplitGetSubKSP(pc, &nSplits, &splitKSPs);
    for (PetscInt k = 0; k < nSplits; ++k) {
        KSPStatus *splitKSPStatus = &(m_splitKSPStatuses[k]);
        fillKSPStatus(splitKSPs[k], splitKSPStatus);
    }
}

/*!
 * Reset the status of the KSP.
 */
void SplitSystemSolver::resetKSPStatus()
{
    SystemSolver::resetKSPStatus();
    resetSplitKSPStatuses();
}

/*!
 * Reset the status of the split KSPs.
 */
void SplitSystemSolver::resetSplitKSPStatuses()
{
    for (KSPStatus &splitStatus : m_splitKSPStatuses) {
        resetKSPStatus(&splitStatus);
    }
}

/*!
 * Destroy the status of the KSP.
 */
void SplitSystemSolver::destroyKSPStatus()
{
    SystemSolver::destroyKSPStatus();
    destroySplitKSPStatuses();
}

/*!
 * Destroy the status of the split KSPs.
 */
void SplitSystemSolver::destroySplitKSPStatuses()
{
    m_splitKSPStatuses.clear();
}


/*!
 * Dump system information.
 *
 * \param systemStream is the stream in which system information is written
 */
void SplitSystemSolver::dumpInfo(std::ostream &systemStream) const
{
    SystemSolver::dumpInfo(systemStream);

    if (systemStream.good()) {
        utils::binary::write(systemStream, m_splitStrategy);
    }
}

/*!
 * Restore system information.
 *
 * \param systemStream is the stream from which system information is read
 */
void SplitSystemSolver::restoreInfo(std::istream &systemStream)
{
    SystemSolver::restoreInfo(systemStream);

    utils::binary::read(systemStream, m_splitStrategy);
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
 *
 * \param matrix is the matrix
 * \param splitStrategy the type of split that will be applied to the system
 * \param splitSizes are the sizes of the splits
 */
void SplitSystemSolver::assembly(const SparseMatrix &matrix, SystemSplitStrategy splitStrategy,
                                 const std::vector<int> &splitSizes)
{
    assembly(matrix, splitStrategy, splitSizes, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * \param matrix is the matrix
 * \param splitStrategy the type of split that will be applied to the system
 * \param splitSizes are the sizes of the splits
 * \param reordering is the reordering that will be applied when assemblying the system
 */
void SplitSystemSolver::assembly(const SparseMatrix &matrix, SystemSplitStrategy splitStrategy,
                                 const std::vector<int> &splitSizes,
                                 const SystemMatrixOrdering &reordering)
{
    // Check if the matrix is assembled
    if (!matrix.isAssembled()) {
        throw std::runtime_error("Unable to assembly the system. The matrix is not yet assembled.");
    }

    // Update matrix
    SplitSystemSparseMatrixAssembler assembler(&matrix, splitStrategy, splitSizes);
    assembly<SplitSystemSolver>(static_cast<const Assembler &>(assembler), reordering);
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
 *
 * \param assembler is the matrix assembler
 */
void SplitSystemSolver::assembly(const Assembler &assembler)
{
    assembly(assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assembling the system
 */
void SplitSystemSolver::assembly(const Assembler &assembler, const SystemMatrixOrdering &reordering)
{
    assembly<SplitSystemSolver>(assembler, reordering);
}

/*!
 * Update all the rows of the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param elements are the elements that will be used to update the rows
 */
void SplitSystemSolver::update(const SparseMatrix &elements)
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
void SplitSystemSolver::update(long nRows, const long *rows, const SparseMatrix &elements)
{
    // Check if the element storage is assembled
    if (!elements.isAssembled()) {
        throw std::runtime_error("Unable to update the system. The element storage is not yet assembled.");
    }

    // Update matrix
    SplitSystemSparseMatrixAssembler assembler(&elements, getSplitStrategy(), getSplitSizes());
    update<SplitSystemSolver>(nRows, rows, static_cast<const Assembler &>(assembler));
}

/*!
 * Update all the rows of the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param assembler is the matrix assembler for the rows that will be updated
 */
void SplitSystemSolver::update(const Assembler &assembler)
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
void SplitSystemSolver::update(long nRows, const long *rows, const Assembler &assembler)
{
    update<SplitSystemSolver>(nRows, rows, assembler);
}

/*!
 * Assemble the matrix.
 *
 * \param assembler is the matrix assembler
 */
void SplitSystemSolver::matrixAssembly(const Assembler &assembler)
{
    const PetscInt *rowReordering = PETSC_NULLPTR;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    // Matrix information
    long nRows = assembler.getRowCount();

    // Split information
    m_splitStrategy = assembler.getSplitStrategy();

    const std::vector<int> &splitBlockSizes = assembler.getSplitSizes();
    int nSplits = splitBlockSizes.size();

    // Create split matrices
    m_splitAs.assign(nSplits * nSplits, PETSC_NULLPTR);
    for (int splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (int splitCol = 0; splitCol < nSplits; ++splitCol) {
            // Create matrix
            //
            // In LOWER mode, only the lower triangular portion of the splits need to
            // be created, whereas in DIAGONAL mode only the diagonal splits need to
            // be created.
            if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_LOWER) {
                if (splitCol > splitRow) {
                    continue;
                }
            } else if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_DIAGONAL) {
                if (splitCol != splitRow) {
                    continue;
                }
            }

            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            Mat *splitMatrix = m_splitAs.data() + splitIndex;
            createMatrix(splitBlockSizes[splitRow], splitBlockSizes[splitCol], splitMatrix);

            MatType splitMatrixType;
            MatGetType(*splitMatrix, &splitMatrixType);

            // Set matrix sizes
            long nSplitRowsElements = assembler.getRowCount() * splitBlockSizes[splitRow];
            long nSplitColsElements = assembler.getColCount() * splitBlockSizes[splitCol];

            long nGlobalSplitRowsElements;
            long nGlobalSplitColsElements;
#if BITPIT_ENABLE_MPI == 1
            nGlobalSplitRowsElements = assembler.getRowGlobalCount() * splitBlockSizes[splitRow];
            nGlobalSplitColsElements = assembler.getColGlobalCount() * splitBlockSizes[splitCol];
#else
            nGlobalSplitRowsElements = nSplitRowsElements;
            nGlobalSplitColsElements = nSplitColsElements;
#endif

            MatSetSizes(*splitMatrix, nSplitRowsElements, nSplitColsElements, nGlobalSplitRowsElements, nGlobalSplitColsElements);

            // Allocate matrix storage
            //
            // When the internal storage of the system matrix was created without taking into account
            // block information, preallocation information should be provided for each row of each
            // block.
            int rowAllocationExpansion;
            int colAllocationExpansion;
            if (strcmp(splitMatrixType, MATSEQAIJ) == 0) {
                rowAllocationExpansion = splitBlockSizes[splitRow];
                colAllocationExpansion = splitBlockSizes[splitCol];
#if BITPIT_ENABLE_MPI == 1
            } else if (strcmp(splitMatrixType, MATMPIAIJ) == 0) {
                rowAllocationExpansion = splitBlockSizes[splitRow];
                colAllocationExpansion = splitBlockSizes[splitCol];
#endif
            } else {
                rowAllocationExpansion = 1;
                colAllocationExpansion = 1;
            }

            long nAllocatedRowElements = rowAllocationExpansion * nRows;

            std::vector<int> d_nnz(nAllocatedRowElements, 0);
            for (long n = 0; n < nRows; ++n) {
                long matrixRow = n;
                if (rowReordering) {
                    matrixRow = rowReordering[matrixRow];
                }

                int nAssemblerRowNZ = assembler.getRowNZCount(n);

                long matrixRowOffset = matrixRow * rowAllocationExpansion;
                for (int n = 0; n < rowAllocationExpansion; ++n) {
                    d_nnz[matrixRowOffset + n] = colAllocationExpansion * nAssemblerRowNZ;
                }
            }

#if BITPIT_ENABLE_MPI == 1
            std::vector<int> o_nnz(nAllocatedRowElements, 0);
            if (isPartitioned()) {
                long nAssemblerCols = assembler.getColCount();

                long assemblerDiagonalBegin = assembler.getColGlobalOffset();
                long assemblerDiagonalEnd   = assemblerDiagonalBegin + nAssemblerCols;

                ConstProxyVector<long> assemblerRowPattern(static_cast<std::size_t>(0), assembler.getMaxRowNZCount());
                for (long n = 0; n < nRows; ++n) {
                    long matrixRow = n;
                    if (rowReordering) {
                        matrixRow = rowReordering[matrixRow];
                    }

                    assembler.getRowPattern(n, &assemblerRowPattern);
                    int nAssemblerRowNZ = assemblerRowPattern.size();

                    long matrixRowOffset = matrixRow * rowAllocationExpansion;
                    for (int k = 0; k < nAssemblerRowNZ; ++k) {
                        long id = assemblerRowPattern[k];
                        if (id < assemblerDiagonalBegin || id >= assemblerDiagonalEnd) {
                            for (int n = 0; n < rowAllocationExpansion; ++n) {
                                o_nnz[matrixRowOffset + n] += colAllocationExpansion;
                                d_nnz[matrixRowOffset + n] -= colAllocationExpansion;
                            }
                        }
                    }
                }
            }
#endif

            if (strcmp(splitMatrixType, MATSEQAIJ) == 0) {
                MatSeqAIJSetPreallocation(*splitMatrix, 0, d_nnz.data());
            } else if (strcmp(splitMatrixType, MATSEQBAIJ) == 0) {
                MatSeqBAIJSetPreallocation(*splitMatrix, splitBlockSizes[splitRow], 0, d_nnz.data());
#if BITPIT_ENABLE_MPI == 1
            } else if (strcmp(splitMatrixType, MATMPIAIJ) == 0) {
                MatMPIAIJSetPreallocation(*splitMatrix, 0, d_nnz.data(), 0, o_nnz.data());
            } else if (strcmp(splitMatrixType, MATMPIBAIJ) == 0) {
                MatMPIBAIJSetPreallocation(*splitMatrix, splitBlockSizes[splitRow], 0, d_nnz.data(), 0, o_nnz.data());
#endif
            } else {
                throw std::runtime_error("Matrix format not supported.");
            }

            // Each process will only set values for its own rows
            MatSetOption(*splitMatrix, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);

#if PETSC_VERSION_GE(3, 12, 0)
            // The first assembly will set a superset of the off-process entries
            // required for all subsequent assemblies. This avoids a rendezvous
            // step in the MatAssembly functions.
            MatSetOption(*splitMatrix, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE);
#endif
        }
    }

    // Create nest matrix
    createMatrix(1, 1, nSplits, nSplits, m_splitAs.data(), &m_A);

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
    for (Mat splitMatrix : m_splitAs) {
        if (!splitMatrix) {
                continue;
        }

        MatSetOption(splitMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    }
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
void SplitSystemSolver::matrixFill(const std::string &filePath)
{
    // Check if the matrix exists
    if (!m_A) {
        throw std::runtime_error("Matrix should be created before filling it.");
    }

    // Fill the matrix
    int nSplits = getSplitCount();
    for (PetscInt splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (PetscInt splitCol = 0; splitCol < nSplits; ++splitCol) {
            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            std::string splitFilePath = generateSplitPath(filePath, splitRow, splitCol);
            fillMatrix(m_splitAs[splitIndex], filePath);
        }
    }
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
void SplitSystemSolver::matrixUpdate(long nRows, const long *rows, const Assembler &assembler)
{
    // Updating the matrix invalidates the KSP
    m_KSPDirty = true;

    // Get split information
    std::vector<int> splitBlockSizes = getSplitSizes();
    std::vector<int> splitBlockOffsets = getSplitOffsets();
    int nSplits = splitBlockSizes.size();

    // Get matrix information
    int blockSize = getBlockSize();

    std::vector<int> blockSplitRows(blockSize);
    for (int i = 0; i < nSplits; ++i) {
        for (int k = 0; k < splitBlockSizes[i]; ++k) {
            blockSplitRows[splitBlockOffsets[i] + k] = i;
        }
    }

    // Get assembler information
    int assemblerBlockSize = assembler.getBlockSize();
    if (assemblerBlockSize != blockSize) {
        std::string message = "Unable to update the matrix.";
        message += " The block size of the assembler is not equal to the block size of the system matrix.";
        throw std::runtime_error(message);
    }

    const long assemblerMaxRowNZ = assembler.getMaxRowNZCount();

    PetscInt colGlobalBegin;
    PetscInt colGlobalEnd;
    MatGetOwnershipRangeColumn(m_A, &colGlobalBegin, &colGlobalEnd);
    colGlobalBegin /= assemblerBlockSize;
    colGlobalEnd /= assemblerBlockSize;

    PetscInt rowGlobalOffset;
    MatGetOwnershipRange(m_A, &rowGlobalOffset, PETSC_NULLPTR);
    rowGlobalOffset /= assemblerBlockSize;

    // Initialize reordering
    const PetscInt *rowReordering = PETSC_NULLPTR;
    if (m_rowReordering) {
        ISGetIndices(m_rowReordering, &rowReordering);
    }

    const PetscInt *colReordering = PETSC_NULLPTR;
    if (m_colReordering) {
        ISGetIndices(m_colReordering, &colReordering);
    }

    // Generate the split indexes
    std::vector<std::vector<std::size_t>> scatterIndexes(nSplits);
    for (int split = 0; split < nSplits; ++split) {
        generateSplitIndexes(split, assemblerMaxRowNZ, scatterIndexes.data() + split);
    }

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

    // Insert matrix values
    //
    // If the sizes of PETSc data types match the sizes of data types expected by bitpit
    // and no column reordering is needed, a direct update of the pattern can be performed,
    // otherwise an intermediate pattern storage is required.
    bool patternDirectUpdate = !colReordering && (sizeof(long) == sizeof(PetscInt));

    ConstProxyVector<long> rowPattern;
    std::vector<PetscInt> petscRowPatternStorage;
    const PetscInt *petscRowPattern;
    if (!patternDirectUpdate) {
        rowPattern.set(ConstProxyVector<long>::INTERNAL_STORAGE, 0, assemblerMaxRowNZ);
        petscRowPatternStorage.resize(assemblerMaxRowNZ);
        petscRowPattern = petscRowPatternStorage.data();
    }

    ConstProxyVector<double> rowValues;
    std::vector<std::vector<PetscScalar>> petscSplitRowValues(m_splitAs.size());
    for (PetscInt splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (PetscInt splitCol = 0; splitCol < nSplits; ++splitCol) {
            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            if (!m_splitAs[splitIndex]) {
                continue;
            }

            int splitMatriBlockSize = splitBlockSizes[splitRow] * splitBlockSizes[splitCol];
            petscSplitRowValues[splitIndex].resize(splitMatriBlockSize * assemblerMaxRowNZ);
        }
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
        assembler.getRowData(n, &rowPattern, &rowValues);
        if (rowValues.size() == 0) {
            continue;
        }

        // Get pattern in PETSc format
        const std::size_t rowPatternSize = rowPattern.size();
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

        // Get values in PETSc format
        //
        // The assembler returns the block row values of the matrix in a row-major order.
        // We need to scatter these values to obtain the values of the block rows of the
        // split matrices.
        for (int blockRow = 0; blockRow < blockSize; ++blockRow) {
            int splitRow = blockSplitRows[blockRow];
            int splitBlockRow = blockRow - splitBlockOffsets[splitRow];

            std::size_t nRowElements = blockSize * rowPatternSize;
            std::size_t blockRowValuesOffset = nRowElements * blockRow;
            const double *blockRowValues = rowValues.data() + blockRowValuesOffset;

            for (int splitCol = 0; splitCol < nSplits; ++splitCol) {
                int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
                if (!m_splitAs[splitIndex]) {
                    continue;
                }

                const std::vector<std::size_t> &splitScatterIndexes = scatterIndexes[splitCol];
                std::size_t nSplitRowElements = splitBlockSizes[splitCol] * rowPatternSize;
                std::size_t splitBlockRowValuesOffset = nSplitRowElements * splitBlockRow;
                double *splitBlockRowValues = petscSplitRowValues[splitIndex].data() + splitBlockRowValuesOffset;
                for (std::size_t k = 0; k < nSplitRowElements; ++k) {
                    splitBlockRowValues[k] = blockRowValues[splitScatterIndexes[k]];
                }
            }
        }

        // Insert values
        for (std::size_t k = 0; k < m_splitAs.size(); ++k) {
            if (!m_splitAs[k]) {
                continue;
            }

            MatSetValuesBlocked(m_splitAs[k], 1, &globalRow, rowPatternSize, petscRowPattern,
                                petscSplitRowValues[k].data(), INSERT_VALUES);
        }
    }

    // Let PETSc assembly the matrix after the update
    for (Mat splitMatrix : m_splitAs) {
        if (!splitMatrix) {
            continue;
        }

        MatAssemblyBegin(splitMatrix, MAT_FINAL_ASSEMBLY);
    }

    for (Mat splitMatrix : m_splitAs) {
        if (!splitMatrix) {
            continue;
        }

        MatAssemblyEnd(splitMatrix, MAT_FINAL_ASSEMBLY);
    }

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
 * \param directory is the directory in which the matrix data file will be written.
 * \param prefix is the prefix added to the name of the file containing matrix data
 */
void SplitSystemSolver::matrixDump(std::ostream &systemStream, const std::string &directory,
                                   const std::string &prefix) const
{
    // Dump split matrices
    int nSplits = getSplitCount();
    if (systemStream.good()) {
        utils::binary::write(systemStream, nSplits);
    }

    for (PetscInt splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (PetscInt splitCol = 0; splitCol < nSplits; ++splitCol) {
            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            std::string splitMatrixName = generateSplitPath(prefix + "A", splitRow, splitCol);
            dumpMatrix(m_splitAs[splitIndex], directory, splitMatrixName);
        }
    }

    // Dump information for creating the main matrix
    if (systemStream.good()) {
        PetscInt rowBlockSize;
        PetscInt colBlockSize;
        MatGetBlockSizes(m_A, &rowBlockSize, &colBlockSize);
        utils::binary::write(systemStream, static_cast<int>(rowBlockSize));
        utils::binary::write(systemStream, static_cast<int>(colBlockSize));
    }
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
void SplitSystemSolver::matrixRestore(std::istream &systemStream, const std::string &directory,
                                      const std::string &prefix, bool redistribute)
#else
void SplitSystemSolver::matrixRestore(std::istream &systemStream, const std::string &directory,
                                      const std::string &prefix)
#endif
{
    // Restore split matrices
    int nSplits;
    utils::binary::read(systemStream, nSplits);

    m_splitAs.assign(nSplits * nSplits, PETSC_NULLPTR);
    for (PetscInt splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (PetscInt splitCol = 0; splitCol < nSplits; ++splitCol) {
            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            std::string splitMatrixName = generateSplitPath(prefix + "A", splitRow, splitCol);
#if BITPIT_ENABLE_MPI==1
            restoreMatrix(directory, splitMatrixName, redistribute, m_splitAs.data() + splitIndex);
#else
            restoreMatrix(directory, splitMatrixName, m_splitAs.data() + splitIndex);
#endif
        }
    }

    // Restore main matrix
    int rowBlockSize;
    utils::binary::read(systemStream, rowBlockSize);
    int colBlockSize;
    utils::binary::read(systemStream, colBlockSize);

    createMatrix(rowBlockSize, colBlockSize, nSplits, nSplits, m_splitAs.data(), &m_A);
}

/*!
 * Destroy the matrix.
 */
void SplitSystemSolver::matrixDestroy()
{
    for (Mat &splitMatrix : m_splitAs) {
        destroyMatrix(&splitMatrix);
    }

    destroyMatrix(&m_A);
}

/*!
 * Create RHS and solution vectors.
 *
 * Vectors will be created, but they will not be initialized.
 *
 * Note that MatNest does *not* depend on VecNest, and in general, Jed Brown (one of the main
 * developers of PETSc) recommend not using VecNest. For this reason, VecNest will not be used.
 */
void SplitSystemSolver::vectorsCreate()
{
    // Create vectors
    SystemSolver::vectorsCreate();

    // Create the split permutations
    vectorsCreateSplitPermutations();
}

/*!
 * Reorder RHS and solution vectors to match the order of the system matrix.
 *
 * \param invert is a flag for inverting the ordering
 */
void SplitSystemSolver::vectorsReorder(bool invert)
{
    SystemSolver::vectorsReorder(invert);

    reorderVector(m_rhs, m_solutionSplitPermutation, !invert);
    reorderVector(m_solution, m_rhsSplitPermutation, !invert);
}

/*!
 * Restore RHS and solution vectors.
 *
 * \param systemStream is the stream from which system information is read
 * \param directory is the directory from which the vector data file will be read
 * \param prefix is the prefix added to the name of the file containing vectors data
 */
void SplitSystemSolver::vectorsRestore(std::istream &systemStream, const std::string &directory,
                                       const std::string &prefix)
{
    // Create vectors
    SystemSolver::vectorsRestore(systemStream, directory, prefix);

    // Create the split permutations
    vectorsCreateSplitPermutations();
}

/*!
 * Create the split permutation for RHS and solution vectors.
 *
 * Given a vector whose elements are ordered according to the full matrix, the split permutations
 * allow to obtain a vector whose elements are ordered according to the split matrices.
 */
void SplitSystemSolver::vectorsCreateSplitPermutations()
{
    int blockSize = getBlockSize();

#if BITPIT_ENABLE_MPI == 1
    PetscInt rhsGlobalBegin;
    PetscInt rhsGlobalEnd;
    VecGetOwnershipRange(m_rhs, &rhsGlobalBegin, &rhsGlobalEnd);
    rhsGlobalBegin /= blockSize;
    rhsGlobalEnd /= blockSize;
    PetscInt rhsSize = rhsGlobalEnd - rhsGlobalBegin;

    PetscInt solutionGlobalBegin;
    PetscInt solutionGlobalEnd;
    VecGetOwnershipRange(m_solution, &solutionGlobalBegin, &solutionGlobalEnd);
    solutionGlobalBegin /= blockSize;
    solutionGlobalEnd /= blockSize;
    PetscInt solutionSize = solutionGlobalEnd - solutionGlobalBegin;

    generateSplitPermutation(rhsSize, static_cast<long>(rhsGlobalBegin), &m_rhsSplitPermutation);
    generateSplitPermutation(solutionSize, static_cast<long>(solutionGlobalBegin), &m_solutionSplitPermutation);
#else
    PetscInt rhsSize;
    VecGetSize(m_rhs, &rhsSize);
    rhsSize /= blockSize;

    PetscInt solutionSize;
    VecGetSize(m_solution, &solutionSize);
    solutionSize /= blockSize;

    generateSplitPermutation(rhsSize, &m_rhsSplitPermutation);
    generateSplitPermutation(solutionSize, &m_solutionSplitPermutation);
#endif
}

/*!
 * Destroy RHS and solution vectors.
 */
void SplitSystemSolver::vectorsDestroy()
{
    // Destroy split permutations
    if (m_rhsSplitPermutation) {
        ISDestroy(&m_rhsSplitPermutation);
    }

    if (m_solutionSplitPermutation) {
        ISDestroy(&m_solutionSplitPermutation);
    }

    // Destroy RHS and solution vectors
    SystemSolver::vectorsDestroy();
}

/*!
 * Export the matrix to the specified file.
 *
 * \param filePath is the path of the file
 * \param fileFormat is the file format that will be used for exporting the matrix, note that
 * the ASCII format may not be able to handle large matrices
 */
void SplitSystemSolver::exportMatrix(const std::string &filePath, FileFormat fileFormat) const
{
    // Export main matrix
    exportMatrix(m_A, filePath, fileFormat);

    // Export split matrices
    int nSplits = getSplitCount();
    for (int splitRow = 0; splitRow < nSplits; ++splitRow) {
        for (int splitCol = 0; splitCol < nSplits; ++splitCol) {
            int splitIndex = getBlockSplitLinearIndex(splitRow, splitCol, nSplits);
            std::string splitFilePath = generateSplitPath(filePath, splitRow, splitCol);
            exportMatrix(m_splitAs[splitIndex], splitFilePath, fileFormat);
        }
    }
}

/*!
 * Set up the preconditioner.
 */
void SplitSystemSolver::setupPreconditioner()
{
    PC pc;
    KSPGetPC(m_KSP, &pc);

    // Setup main preconditioner
    PCSetType(pc, PCFIELDSPLIT);
    if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_LOWER) {
        PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
    } else if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_DIAGONAL) {
        PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
    } else {
        PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
    }

    int nSplits = getSplitCount();
    std::vector<IS> splitIndexSets(nSplits);
    MatNestGetISs(m_A, splitIndexSets.data(), PETSC_NULLPTR);
    for (int k = 0; k < nSplits; ++k) {
        PCFieldSplitSetIS(pc, PETSC_NULLPTR, splitIndexSets[k]);
    }

    PCSetUp(pc);

    // Set preconditioners of the split
    setupSplitPreconditioners();
}

/*!
 * Set up the Krylov subspace method used to solve the system.
 */
void SplitSystemSolver::setupSplitPreconditioners()
{
    PC pc;
    KSPGetPC(m_KSP, &pc);

    PetscInt nSplits;
    KSP *splitKSPs;
    PCFieldSplitGetSubKSP(pc, &nSplits, &splitKSPs);

    for (PetscInt k = 0; k < nSplits; ++k) {
        PC splitPc;
        KSPGetPC(splitKSPs[k], &splitPc);
        const KSPOptions &splitOptions = getSplitKSPOptions(k);

        setupPreconditioner(splitPc, splitOptions);
    }
}

/*!
 * Set up the Krylov subspace method used to solve the system.
 */
void SplitSystemSolver::setupKrylov()
{
    // Setup main Krylov subspace method
    //
    // In LOWER and DIAGONAL modes, the solution of the system is performed by the Krylov subspace
    // methods of the single splits. Thus, in these modes, the main Krylov subspace method only
    // needs to apply the preconditioner.
    if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_LOWER) {
        KSPSetType(m_KSP, KSPPREONLY);
    } else if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_DIAGONAL) {
        KSPSetType(m_KSP, KSPPREONLY);
    } else {
        setupKrylov(m_KSP, getKSPOptions());
    }

    // Setup Krylov subspace methods of the splits
    setupSplitKrylovs();
}

/*!
 * Set up the Krylov subspace method used to solve the system.
 */
void SplitSystemSolver::setupSplitKrylovs()
{
    PC pc;
    KSPGetPC(m_KSP, &pc);

    PetscInt nSplits;
    KSP *splitKSPs;
    PCFieldSplitGetSubKSP(pc, &nSplits, &splitKSPs);

    for (PetscInt k = 0; k < nSplits; ++k) {
        // In FULL mode, the solution of the system is performed by the Krylov subspace
        // methods of the single splits. Thus, in these modes, the Krylov subspaces of
        // the single blocks only need to apply the preconditioner.
        const KSPOptions &splitOptions = getSplitKSPOptions(k);
        if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_LOWER) {
            setupKrylov(splitKSPs[k], splitOptions);
        } else if (m_splitStrategy == SystemSplitStrategy::SYSTEM_SPLIT_STRATEGY_DIAGONAL) {
            setupKrylov(splitKSPs[k], splitOptions);
        } else {
            KSPSetType(splitKSPs[k], KSPPREONLY);
        }
    }
}

/*!
 * Generate the split permutation.
 *
 * Given a vector whose elements are ordered according to the full matrix, the split permutations
 * allow to obtain a vector whose elements are ordered according to the split matrices.
 *
 * \param nItems is the number of items in the row/column
 */
#if BITPIT_ENABLE_MPI == 1
/*!
 * \param itemsGlobalOffset is the global offset of the items that belong to the current process
 */
#endif
/*!
 * \param split is the split
 * \param[out] permutations on output will contain the permutations
 */
#if BITPIT_ENABLE_MPI == 1
void SplitSystemSolver::generateSplitPermutation(long nItems, long itemsGlobalOffset, IS *permutation) const
#else
void SplitSystemSolver::generateSplitPermutation(long nItems, IS *permutation) const
#endif
{
   std::vector<std::size_t> workspace;

    // Get split information
    std::vector<int> splitBlockSizes = getSplitSizes();
    std::vector<int> splitBlockOffsets = getSplitOffsets();
    int nSplits = splitBlockSizes.size();

    // Matrix information
    int blockSize = getBlockSize();
    std::size_t nElements = blockSize * nItems;

    // Create permutations
    PetscInt elementsGlobalOffset;
#if BITPIT_ENABLE_MPI == 1
    elementsGlobalOffset = static_cast<PetscInt>(blockSize * itemsGlobalOffset);
#else
    elementsGlobalOffset = 0;
#endif

    PetscInt *permutationStorage;
    PetscMalloc(nElements * sizeof(PetscInt), &permutationStorage);

    PetscInt *splitPermutationStorage = permutationStorage;
    for (int split = 0; split < nSplits; ++split) {
        int splitBlockSize = splitBlockSizes[split];
        std::size_t nSplitElements = splitBlockSize * nItems;

        workspace.resize(nSplitElements);
        generateSplitIndexes(split, nItems, &workspace);
        for (std::size_t k = 0; k < nSplitElements; ++k) {
            splitPermutationStorage[k] = elementsGlobalOffset + static_cast<PetscInt>(workspace[k]);

        }
        splitPermutationStorage += nSplitElements;
    }

#if BITPIT_ENABLE_MPI == 1
    ISCreateGeneral(getCommunicator(), nItems, permutationStorage, PETSC_OWN_POINTER, permutation);
#else
    ISCreateGeneral(PETSC_COMM_SELF, nItems, permutationStorage, PETSC_OWN_POINTER, permutation);
#endif
    ISSetPermutation(*permutation);
}

/*!
 * Generate split indexes.
 *
 * Split indexes define the indexes of the matrix that are assigned to the specified split.
 * The k-th element of a split matrix row/column corresponds to the indexes[k]-th element
 * of the corresponding row/column of the full matrix.
 *
 * \param split is the split
 * \param nItems is the number of items in the row/column
 * \param indexes on output will contain the split indexes
 */
void SplitSystemSolver::generateSplitIndexes(int split, long nItems, std::vector<std::size_t> *indexes) const
{
    // Get split information
    int splitBlockSize   = getSplitSizes()[split];
    int splitBlockOffset = getSplitOffsets()[split];

    // Get matrix information
    int blockSize = getBlockSize();

    // Initialize index storage
    indexes->resize(nItems * splitBlockSize);

    // Evaluate the indexes
    for (long k = 0; k < nItems; ++k) {
        for (int n = 0; n < splitBlockSize; ++n) {
            (*indexes)[splitBlockSize * k + n] = blockSize * k + splitBlockOffset + n;
        }
    }
}

/*!
 * Generate the split path associated with the specified indexes.
 *
 * The split file path is generated adding "_split_ij" before the extension of the specified
 * path. For example, given the path "data/matrix.dat", the split path associated with the
 * index 1 is "data/matrix_split_1.dat".
 *
 * \param path is the path
 * \param i is the index associated with the split
 */
std::string SplitSystemSolver::generateSplitPath(const std::string &path, int i) const
{
    return generateSplitPath(path, std::to_string(i));
}

/*!
 * Generate the split path associated with the specified indexes.
 *
 * The split file path is generated adding "_split_ij" before the extension of the specified
 * path. For example, given the path "data/matrix.dat", the split path associated with the
 * indexes (1,0) is "data/matrix_split_10.dat".
 *
 * \param path is the path
 * \param i is the first index associated with the split
 * \param j is the second index associated with the split
 */
std::string SplitSystemSolver::generateSplitPath(const std::string &path, int i, int j) const
{
    return generateSplitPath(path, std::to_string(i) + std::to_string(j));
}

/*!
 * Generate the split path associated with the specified indexes.
 *
 * The split file path is generated adding "_split_ij" before the extension of the specified
 * path. For example, given the path "data/matrix.dat", the split path associated with the
 * index INDEX is "data/matrix_split_INDEX.dat".
 *
 * \param path is the path
 * \param index is the index associated with the split
 */
std::string SplitSystemSolver::generateSplitPath(const std::string &path, const std::string &index) const
{
    std::string splitLabel = "_split_" + index;

    std::size_t dotPosition = path.find_last_of(".");

    std::string splitPath;
    std::string extension;
    std::string basePath;
    if (dotPosition != std::string::npos) {
        extension = path.substr(dotPosition + 1);
        basePath  = path.substr(0, dotPosition);
        splitPath = basePath + splitLabel + "." + extension;
    } else {
        splitPath = path + splitLabel;
    }

    return splitPath;
}

/*!
 * Get the linear index associated with the split block with the given indices.
 *
 * The number of row/column splits is inferred from the system matrix.
 *
 * \param i is the row index of the split block
 * \param j is the column index of the split block
 * \result The linear index associated with the split block with the given indices.
 */
int SplitSystemSolver::getBlockSplitLinearIndex(int i, int j) const
{
    int nSplits = getSplitCount();

    return getBlockSplitLinearIndex(i, j, nSplits);
}

/*!
 * Get the linear index associated with the split block with the given indices.
 *
 * \param i is the row index of the split block
 * \param j is the column index of the split block
 * \param nSplits are the number of splits
 * \result The linear index associated with the split block with the given indices.
 */
int SplitSystemSolver::getBlockSplitLinearIndex(int i, int j, int nSplits) const
{
    return linearalgebra::linearIndexRowMajor(i, j, nSplits, nSplits);
}

}
