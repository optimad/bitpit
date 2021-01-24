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

#include "stencil_solver.hpp"

namespace bitpit {

/*!
 * \class DiscretizationStencilSolverAssembler
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilSolverAssembler class defines an assembler
 * for building the stencil solver.
 */

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(const std::vector<stencil_t> *stencils)
    : DiscretizationStencilSolverAssembler(MPI_COMM_SELF, false, stencils)
{
}

/*!
 * Constructor.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, const std::vector<stencil_t> *stencils)
#else
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(const std::vector<stencil_t> *stencils)
#endif
    : StencilSolverAssembler(), m_stencils(stencils),
      m_blockSize(-1)
{
    // Initialize block size
    initializeBlockSize();

    // Count the DOFs
    m_nDOFs = stencils->size();

#if BITPIT_ENABLE_MPI==1
    m_nGlobalDOFs = m_nDOFs;
    if (partitioned) {
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalDOFs, 1, MPI_LONG, MPI_SUM, communicator);
    }
#endif

    // Count maximum non-zero elements
    m_maxRowNZ = 0;
    for (long n = 0; n < getRowCount(); ++n) {
        m_maxRowNZ = std::max(getRowNZCount(n), m_maxRowNZ);
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
 * Get the stencil block size.
 *
 * \result The stencil block size.
 */
template<typename stencil_t>
int DiscretizationStencilSolverAssembler<stencil_t>::getBlockSize() const
{
    return m_blockSize;
}

/*!
 * Get the number of rows of the matrix.
 *
 * \result The number of rows of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getRowCount() const
{
    return m_nDOFs;
}

/*!
 * Get the number of columns of the matrix.
 *
 * \result The number of columns of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColCount() const
{
    return m_nDOFs;
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Get the number of global rows of the matrix.
 *
 * \result The number of global rows of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getRowGlobalCount() const
{
    return m_nGlobalDOFs;
}

/*!
 * Get the number of global columns of the matrix.
 *
 * \result The number of global columns of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColGlobalCount() const
{
    return m_nGlobalDOFs;
}

/*!
 * Get global row offset.
 *
 * \result The global row offset.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getRowGlobalOffset() const
{
    return m_globalDOFOffset;
}

/*!
 * Get global column offset.
 *
 * \result The global column offset.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColGlobalOffset() const
{
    return m_globalDOFOffset;
}
#endif

/*!
 * Get the number of non-zero elements in the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The number of non-zero elements in the specified row.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getRowNZCount(long rowIndex) const
{
    const stencil_t &stencil = (*m_stencils)[rowIndex];
    std::size_t stencilSize = stencil.size();

    return (m_blockSize * stencilSize);
}

/**
 * Get the maximum number of non-zero elements per row.
 *
 * \result The maximum number of non-zero elements per row.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getMaxRowNZCount() const
{
    return m_maxRowNZ;
}

/*!
 * Get the values of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the values of the specified row
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const
{
    // Get stencil information
    const stencil_t &stencil = (*m_stencils)[rowIndex];
    std::size_t stencilSize = stencil.size();

    // Get pattern
    const long *patternData = stencil.patternData();
    if (m_blockSize == 1) {
        pattern->set(patternData, stencilSize);
    } else {
        std::size_t expandedPatternSize = m_blockSize * stencilSize;
        long *expandedPatternStorage = pattern->set(expandedPatternSize);

        std::size_t expandedPatternIdx = 0;
        for (std::size_t k = 0; k < stencilSize; ++k) {
            long patternBlockOffset = patternData[k] * m_blockSize;
            for (int i = 0; i < m_blockSize; ++i) {
                expandedPatternStorage[expandedPatternIdx++] = patternBlockOffset + i;
            }
        }
    }
}

/*!
 * Get the values of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    _getRowValues(rowIndex, values);
}

/*!
 * Get the values of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t>::_getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    // Get stencil information
    const stencil_t &stencil = (*m_stencils)[rowIndex];

    // Get values
    values->set(stencil.weightData(), stencil.size());
}

/*!
 * Get the values of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t>::_getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    // Get stencil information
    const stencil_t &stencil = (*m_stencils)[rowIndex];
    std::size_t stencilSize = stencil.size();

    // Get values
    const typename stencil_t::weight_type *weightData = stencil.weightData();

    std::size_t expandedValuesSize = m_blockSize * stencilSize;
    double *expandedValuesStorage = values->set(expandedValuesSize);

    std::size_t expandedValuesIdx = 0;
    for (std::size_t k = 0; k < stencilSize; ++k) {
        for (int i = 0; i < m_blockSize; ++i) {
            expandedValuesStorage[expandedValuesIdx++] = getRawValue(weightData[k], i);
        }
    }
}

/*!
 * Get the constant associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The constant associated with the specified row.
 */
template<typename stencil_t>
double DiscretizationStencilSolverAssembler<stencil_t>::getRowConstant(long rowIndex) const
{
    return _getRowConstant(rowIndex);
}

/*!
 * Get the constant associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The constant associated with the specified row.
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
double DiscretizationStencilSolverAssembler<stencil_t>::_getRowConstant(long rowIndex) const
{
    // Get stencil information
    const stencil_t &stencil = (*m_stencils)[rowIndex];

    // Get constant
    return getRawValue(stencil.getConstant(), 0);
}

/*!
 * Get the constant associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The constant associated with the specified row.
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
double DiscretizationStencilSolverAssembler<stencil_t>::_getRowConstant(long rowIndex) const
{
    // Get stencil information
    const stencil_t &stencil = (*m_stencils)[rowIndex];
    const typename stencil_t::weight_type &stencilConstant = stencil.getConstant();

    // Get constant
    double constant = 0.;
    for (int i = 0; i < m_blockSize; ++i) {
        constant += getRawValue(stencilConstant, i);
    }

    return constant;
}

/*!
* \ingroup discretization
* \class DiscretizationStencilSolver
*
* The DiscretizationStencilSolver class handles the solution of linear systems
* assembled from discretization stencils.
*/

/*!
* Constuctor
*
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
template<typename stencil_t>
DiscretizationStencilSolver<stencil_t>::DiscretizationStencilSolver(bool debug)
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
template<typename stencil_t>
DiscretizationStencilSolver<stencil_t>::DiscretizationStencilSolver(const std::string &prefix, bool debug)
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
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::clear(bool release)
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
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const std::vector<stencil_t> &stencils)
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
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(MPI_Comm communicator, bool partitioned, const std::vector<stencil_t> &stencils)
#else
/*!
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const std::vector<stencil_t> &stencils)
#endif
{
    // Create the assembler
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler<stencil_t> assembler(communicator, partitioned, &stencils);
#else
    DiscretizationStencilSolverAssembler<stencil_t> assembler(&stencils);
#endif

    // Assembly the system
#if BITPIT_ENABLE_MPI==1
    assembly(communicator, partitioned, assembler);
#else
    assembly(assembler);
#endif
}

#if BITPIT_ENABLE_MPI==1
/*!
* Assembly the stencil solver.
*
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const StencilSolverAssembler &assembler)
{
    assembly(MPI_COMM_SELF, false, assembler);
}

/*!
* Initialize the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(MPI_Comm communicator, bool partitioned, const StencilSolverAssembler &assembler)
#else
/*!
* Initialize the stencil solver.
*
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const StencilSolverAssembler &assembler)
#endif
{
    // Assembly system
#if BITPIT_ENABLE_MPI==1
    SystemSolver::assembly(communicator, partitioned, assembler);
#else
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
 * Update all the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param stencils are the stencils that will be used to update the rows
 */
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::update(const std::vector<stencil_t> &stencils)
{
    update(getRowCount(), nullptr, stencils);
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
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::update(const std::vector<long> &rows, const std::vector<stencil_t> &stencils)
{
    update(rows.size(), rows.data(), stencils);
}

/*!
 * Update the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param nRows is the number of stencils that will be updated
 * \param rows are the rows of the stencils that will be updated,
 * if a null pointer is passed, the stencils that will be updated are the
 * stencils from 0 to (nRows - 1).
 * \param stencils are the stencils that will be used to update the rows
 */
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::update(std::size_t nRows, const long *rows, const std::vector<stencil_t> &stencils)
{
    // Update the system
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler<stencil_t> assembler(getCommunicator(), isPartitioned(), &stencils);
#else
    DiscretizationStencilSolverAssembler<stencil_t> assembler(&stencils);
#endif
    SystemSolver::update(nRows, rows, assembler);

    // Update the constants
    for (std::size_t n = 0; n < nRows; ++n) {
        long row;
        if (rows) {
            row = rows[n];
        } else {
            row = n;
        }

        m_constants[row] = assembler.getRowConstant(n);
    }
}

/*!
* Solve the system.
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::solve()
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
