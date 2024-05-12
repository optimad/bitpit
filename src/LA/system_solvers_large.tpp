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

#ifndef __BITPIT_SYSTEM_SOLVERS_LARGE_TPP__
#define __BITPIT_SYSTEM_SOLVERS_LARGE_TPP__

namespace bitpit {

/*!
 * \class ProxySystemMatrixOrdering
 * \ingroup system_solver_large
 *
 * \brief The ProxySystemMatrixOrdering class defines allows to use a matrix
 * ordering defined in some external storage.
 *
 * If the system is partitioned, each process can reorder only it's local
 * part of the matix.
 */

/*!
 * Constructor.
 *
 * \param rowRankStorage is the storage for the row ranks
 * \param colRankStorage is the storage for the column ranks
 */
template<typename RowRankStorage, typename ColRankStorage>
ProxySystemMatrixOrdering<RowRankStorage, ColRankStorage>::ProxySystemMatrixOrdering(const RowRankStorage *rowRankStorage,
                                                                                     const ColRankStorage *colRankStorage)
    : m_rowRankStorage(rowRankStorage), m_colRankStorage(colRankStorage)
{
}

/*!
 * Get the rank of the specified local row.
 *
 * \param row is the local local row
 * \result The rank of the specified local row.
 */
template<typename RowRankStorage, typename ColRankStorage>
long ProxySystemMatrixOrdering<RowRankStorage, ColRankStorage>::getRowRank(long row) const
{
    return (*m_rowRankStorage)[row];
}

/*!
 * Get the rank of the specified local column.
 *
 * \param col is the local column
 * \result The rank of the specified local column.
 */
template<typename RowRankStorage, typename ColRankStorage>
long ProxySystemMatrixOrdering<RowRankStorage, ColRankStorage>::getColRank(long col) const
{
    return (*m_colRankStorage)[col];
}

/*!
 * Assembly the system.
 *
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assemblying the system
 */
template<typename DerivedSystemSolver>
void SystemSolver::assembly(const typename DerivedSystemSolver::Assembler &assembler,
                            const SystemMatrixOrdering &reordering)
{
    // Clear the system
    clear();

    // Set reordering
    setReordering(assembler.getRowCount(), assembler.getColCount(), reordering);

#if BITPIT_ENABLE_MPI == 1
    // Set the communicator
    setCommunicator(assembler.getCommunicator());

    // Detect if the system is partitioned
    m_partitioned = assembler.isPartitioned();
#endif

    // Initialize matrix
    static_cast<DerivedSystemSolver &>(*this).matrixCreate(assembler);
    static_cast<DerivedSystemSolver &>(*this).matrixFill(assembler);

    // Initialize RHS and solution vectors
    vectorsCreate();

    // The system is now assembled
    m_assembled = true;

    // Initialize KSP options
    initializeKSPOptions();

    // Initialize KSP statuses
    initializeKSPStatus();
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
template<typename DerivedSystemSolver>
void SystemSolver::update(long nRows, const long *rows, const typename DerivedSystemSolver::Assembler &assembler)
{
    // Check if the system is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to update the system. The system is not yet assembled.");
    }

    // Update matrix
    static_cast<DerivedSystemSolver &>(*this).matrixUpdate(nRows, rows, assembler);
}

}

#endif
