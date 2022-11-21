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

}

#endif
