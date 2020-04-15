/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include "stencil_solver.hpp"

namespace bitpit {

/*!
 * \class StencilSolverAssembler
 * \ingroup discretization
 *
 * \brief The StencilSolverAssembler class defines an assembler for
 * building the stencil solver.
 */

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
StencilSolverAssembler::StencilSolverAssembler(const std::vector<StencilScalar> *stencils)
    : StencilSolverAssembler(MPI_COMM_SELF, false, stencils)
{
}

/*!
 * Constructor.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 * \param stencils are the stencils
 */
StencilSolverAssembler::StencilSolverAssembler(MPI_Comm communicator, bool partitioned, const std::vector<StencilScalar> *stencils)
#else
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
StencilSolverAssembler::StencilSolverAssembler(const std::vector<StencilScalar> *stencils)
#endif
    : SystemMatrixAssembler(), m_stencils(stencils)
{
    // Count the DOFs
    m_nDOFs = m_stencils->size();

#if BITPIT_ENABLE_MPI==1
    m_nGlobalDOFs = m_nDOFs;
    if (partitioned) {
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalDOFs, 1, MPI_LONG, MPI_SUM, communicator);
    }
#endif

    // Count maximum non-zero elements
    m_maxRowNZ = 0;
    for (const StencilScalar &stencil : *m_stencils) {
        m_maxRowNZ = std::max(m_maxRowNZ, (long) stencil.size());
    }

#if BITPIT_ENABLE_MPI==1
    // Get offsets
    m_globalDOFOffset = 0;
    if (partitioned) {
        int nProcessors;
        MPI_Comm_size(communicator, &nProcessors);

        std::vector<long> nRankDOFs(nProcessors);
        MPI_Allgather(&m_nDOFs, 1, MPI_LONG, nRankDOFs.data(), 1, MPI_LONG, communicator);

        int rank;
        MPI_Comm_rank(communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_globalDOFOffset += nRankDOFs[i];
        }
    }
#endif
}

/*!
 * Get the number of rows of the matrix.
 *
 * \result The number of rows of the matrix.
 */
long StencilSolverAssembler::getRowCount() const
{
    return m_nDOFs;
}

/*!
 * Get the number of columns of the matrix.
 *
 * \result The number of columns of the matrix.
 */
long StencilSolverAssembler::getColCount() const
{
    return m_nDOFs;
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Get the number of global rows of the matrix.
 *
 * \result The number of global rows of the matrix.
 */
long StencilSolverAssembler::getRowGlobalCount() const
{
    return m_nGlobalDOFs;
}

/*!
 * Get the number of global columns of the matrix.
 *
 * \result The number of global columns of the matrix.
 */
long StencilSolverAssembler::getColGlobalCount() const
{
    return m_nGlobalDOFs;
}

/*!
 * Get global row offset.
 *
 * \result The global row offset.
 */
long StencilSolverAssembler::getRowGlobalOffset() const
{
    return m_globalDOFOffset;
}

/*!
 * Get global column offset.
 *
 * \result The global column offset.
 */
long StencilSolverAssembler::getColGlobalOffset() const
{
    return m_globalDOFOffset;
}
#endif

/*!
 * Get the number of non-zero elements in the specified row.
 *
 * \param row is the row of the matrix
 * \result The number of non-zero elements in the specified row.
 */
long StencilSolverAssembler::getRowNZCount(long row) const
{
    return (*m_stencils)[row].size();
}

/**
 * Get the maximum number of non-zero elements per row.
 *
 * \result The maximum number of non-zero elements per row.
 */
long StencilSolverAssembler::getMaxRowNZCount() const
{
    return m_maxRowNZ;
}

/*!
 * Get the values of the specified row.
 *
 * \param row is the row of the matrix
 * \param pattern on output will contain the values of the specified row
 */
void StencilSolverAssembler::getRowPattern(long row, ConstProxyVector<long> *pattern) const
{
    const StencilScalar &stencil = (*m_stencils)[row];
    pattern->set(stencil.patternData(), stencil.size());
}

/*!
 * Get the values of the specified row.
 *
 * \param row is the row of the matrix
 * \param pattern on output will contain the values of the specified row
 */
void StencilSolverAssembler::getRowValues(long row, ConstProxyVector<double> *values) const
{
    const StencilScalar &stencil = (*m_stencils)[row];
    values->set(stencil.weightData(), stencil.size());
}

/*!
 * Get the constant associated with the specified row.
 *
 * \param row is the row of the matrix
 * \result The constant associated with the specified row.
 */
double StencilSolverAssembler::getRowConstant(long row) const
{
    const StencilScalar &stencil = (*m_stencils)[row];

    return stencil.getConstant();
}

/*!
* \ingroup discretization
* \class StencilScalarSolver
*
* The StencilScalarSolver class handles the solution of linear systems assembled
* from scalar discretization stencils.
*/

/*!
* Constuctor
*
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
StencilScalarSolver::StencilScalarSolver(bool debug)
    : SystemSolver("", debug)
{
}

/*!
* Constuctor
*
* \param prefix is the prefix string to prepend to all option requests
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
StencilScalarSolver::StencilScalarSolver(const std::string &prefix, bool debug)
    : SystemSolver(prefix, debug)
{
}

/*!
* Clear the stencil solver
*
* \param release if it's true the memory hold by the stencil solver will be
* released, otherwise the stencil solver will be cleared but its memory will
* not be relased
 */
void StencilScalarSolver::clear(bool release)
{
    SystemSolver::clear();

    if (release) {
        std::vector<double>().swap(m_constants);
    } else {
        m_constants.clear();
    }
}

#if BITPIT_ENABLE_MPI==1
/*!
* Assembly the stencil solver.
*
* \param stencils are the stencils
*/
void StencilScalarSolver::assembly(const std::vector<StencilScalar> &stencils)
{
    assembly(MPI_COMM_SELF, false, stencils);
}

/*!
* Initialize the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param stencils are the stencils
*/
void StencilScalarSolver::assembly(MPI_Comm communicator, bool partitioned, const std::vector<StencilScalar> &stencils)
#else
/*!
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
void StencilScalarSolver::assembly(const std::vector<StencilScalar> &stencils)
#endif
{
    // Assembly system
#if BITPIT_ENABLE_MPI==1
    StencilSolverAssembler assembler(communicator, partitioned, &stencils);
    SystemSolver::assembly(communicator, partitioned, assembler);
#else
    StencilSolverAssembler assembler(&stencils);
    SystemSolver::assembly(assembler);
#endif

    // Set constants
    long nRows = assembler.getRowCount();
    m_constants.resize(nRows);
    for (long n = 0; n < nRows; ++n) {
        m_constants[n] = assembler.getRowConstant(n);
    }
}

/*!
 * Update the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param rows are the global indices of the rows that will be updated
 * \param stencils are the stencils that will be used to update the rows
 */
void StencilScalarSolver::update(const std::vector<long> &rows, const std::vector<StencilScalar> &stencils)
{
    // Assembly system
#if BITPIT_ENABLE_MPI==1
    StencilSolverAssembler assembler(getCommunicator(), isPartitioned(), &stencils);
#else
    StencilSolverAssembler assembler(&stencils);
#endif
    SystemSolver::update(rows.size(), rows.data(), assembler);

    // Set constants
    long nUpdatedRows = rows.size();
    for (long n = 0; n < nUpdatedRows; ++n) {
        long row = rows[n];
        m_constants[row] = assembler.getRowConstant(n);
    }
}

/*!
* Solve the system.
*/
void StencilScalarSolver::solve()
{
    // Check if the stencil solver is assembled
    if (!isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The stencil solver is not yet assembled.");
    }

    // Subtract constant terms to the RHS
    long nUnknowns = getRowCount();
    double *raw_rhs = getRHSRawPtr();
    for (long i = 0; i < nUnknowns; ++i) {
        raw_rhs[i] -= m_constants[i];
    }
    restoreRHSRawPtr(raw_rhs);

    // Solve the system
    SystemSolver::solve();
}

}
