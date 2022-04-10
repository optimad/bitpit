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
 * \class DiscretizationStencilStorageInterface
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilStorageInterface class defines the interface
 * for stencil storage.
 */

/*!
 * \class DiscretizationStencilProxyBaseStorage
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilProxyBaseStorage class defines a proxy for
 * stencil storage.
 */

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 * \param stride is the stride
 */
template<typename stencil_t, typename stencil_container_t>
DiscretizationStencilProxyBaseStorage<stencil_t, stencil_container_t>::DiscretizationStencilProxyBaseStorage(const stencil_container_t *stencils, int stride)
    : DiscretizationStencilStorageInterface<stencil_t>(),
      m_stencils(stencils), m_stride(stride)
{
}

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t, typename stencil_container_t>
DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::DiscretizationStencilProxyStorage(const stencil_container_t *stencils)
    : DiscretizationStencilProxyBaseStorage<stencil_t, stencil_container_t>(stencils, 1)
{
}

/*!
 * Get the size of the container.
 *
 * \result The size of the container.
 */
template<typename stencil_t, typename stencil_container_t>
std::size_t DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::size() const
{
    return this->m_stencils->size();
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::at(long rowIndex) const
{
    return (*(this->m_stencils))[rowIndex];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowRawIndex is the raw index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::rawAt(std::size_t rowRawIndex) const
{
    return (*(this->m_stencils))[rowRawIndex];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockIndex is the index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::at(long blockIndex, int componentIdx) const
{
    return (*(this->m_stencils))[blockIndex * this->m_stride + componentIdx];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockRawIndex is the raw index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::rawAt(std::size_t blockRawIndex, int componentIdx) const
{
    return (*(this->m_stencils))[blockRawIndex * this->m_stride + componentIdx];
}

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::DiscretizationStencilProxyStorage(const PiercedStorage<stencil_t> *stencils)
    : DiscretizationStencilProxyBaseStorage<stencil_t, PiercedStorage<stencil_t>>(stencils, stencils->getFieldCount())
{
}

/*!
 * Get the size of the container.
 *
 * \result The size of the container.
 */
template<typename stencil_t>
std::size_t DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::size() const
{
    return this->m_stencils->getKernel()->size() * this->m_stride;
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::at(long rowIndex) const
{
    return this->m_stencils->at(rowIndex);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowRawIndex is the raw index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::rawAt(std::size_t rowRawIndex) const
{
    return this->m_stencils->rawAt(rowRawIndex);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockIndex is the index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::at(long blockIndex, int componentIdx) const
{
    return this->m_stencils->at(blockIndex, componentIdx);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockRawIndex is the raw index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::rawAt(std::size_t blockRawIndex, int componentIdx) const
{
    return this->m_stencils->rawAt(blockRawIndex, componentIdx);
}

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
template<typename stencil_container_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(const stencil_container_t *stencils)
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
template<typename stencil_container_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, const stencil_container_t *stencils)
    : DiscretizationStencilSolverAssembler(communicator, partitioned, std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>>(new DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>(stencils)))
{
}
#else
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
template<typename stencil_container_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(const stencil_container_t *stencils)
    : DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>>(new DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>(stencils)))
{
}
#endif

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils)
    : DiscretizationStencilSolverAssembler(MPI_COMM_SELF, false, std::move(stencils))
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
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils)
#else
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils)
#endif
    : DiscretizationStencilSolverAssembler()
{
    setStencils(std::move(stencils));
#if BITPIT_ENABLE_MPI==1
    setMatrixSizes(communicator, partitioned);
#else
    setMatrixSizes();
#endif
    setBlockSize();
    setMaximumRowNZ();
}

/*!
 * Constructor.
 */
template<typename stencil_t>
DiscretizationStencilSolverAssembler<stencil_t>::DiscretizationStencilSolverAssembler()
    : StencilSolverAssembler()
{
}

/*!
 * Set block size.
 *
 * \param blockSize is the block size
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setBlockSize(int blockSize)
{
    m_blockSize = blockSize;
}

/*!
 * Set the stencils.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setStencils(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils)
{
    m_stencils = std::move(stencils);
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Set matrix sizes.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMatrixSizes(MPI_Comm communicator, bool partitioned)
{
    setMatrixSizes(m_stencils->size(), m_stencils->size(), communicator, partitioned);
}

/*!
 * Set matrix sizes.
 *
 * \param nRows are the rows of the matrix
 * \param nCols are the columns of the matrix
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMatrixSizes(long nRows, long nCols, MPI_Comm communicator, bool partitioned)
#else
/*!
 * Set matrix sizes.
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMatrixSizes()
{
    setMatrixSizes(m_stencils->size(), m_stencils->size());
}

/*!
 * Set matrix sizes.
 *
 * \param nRows are the rows of the matrix
 * \param nCols are the columns of the matrix
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMatrixSizes(long nRows, long nCols)
#endif
{
    // Set system sizes
    m_nRows = nRows;
    m_nCols = nCols;

#if BITPIT_ENABLE_MPI==1
    // Global system sizes
    m_nGlobalRows = nRows;
    m_nGlobalCols = nCols;
    if (partitioned) {
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalRows, 1, MPI_LONG, MPI_SUM, communicator);
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalCols, 1, MPI_LONG, MPI_SUM, communicator);
    }

    // Global offsets
    m_globalRowOffset = 0;
    m_globalColOffset = 0;
    if (partitioned) {
        int nProcessors;
        MPI_Comm_size(communicator, &nProcessors);

        std::vector<long> nRankRows(nProcessors);
        MPI_Allgather(&m_nRows, 1, MPI_LONG, nRankRows.data(), 1, MPI_LONG, communicator);

        std::vector<long> nRankCols(nProcessors);
        MPI_Allgather(&m_nCols, 1, MPI_LONG, nRankCols.data(), 1, MPI_LONG, communicator);

        int rank;
        MPI_Comm_rank(communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_globalRowOffset += nRankRows[i];
            m_globalColOffset += nRankCols[i];
        }
    }
#endif
}

/*!
 * Set the maximum number of non-zero element on a single row.
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMaximumRowNZ()
{
    long maxRowNZ = 0;
    for (long n = 0; n < getRowCount(); ++n) {
        maxRowNZ = std::max(getRowNZCount(n), maxRowNZ);
    }

    setMaximumRowNZ(maxRowNZ);
}

/*!
 * Set the maximum number of non-zero element on a single row.
 *
 * \param maxRowNZ is the maximum number of non-zero element on a single row
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::setMaximumRowNZ(long maxRowNZ)
{
    m_maxRowNZ = maxRowNZ;
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
    return m_nRows;
}

/*!
 * Get the number of columns of the matrix.
 *
 * \result The number of columns of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColCount() const
{
    return m_nCols;
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
    return m_nGlobalRows;
}

/*!
 * Get the number of global columns of the matrix.
 *
 * \result The number of global columns of the matrix.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColGlobalCount() const
{
    return m_nGlobalCols;
}

/*!
 * Get global row offset.
 *
 * \result The global row offset.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getRowGlobalOffset() const
{
    return m_globalRowOffset;
}

/*!
 * Get global column offset.
 *
 * \result The global column offset.
 */
template<typename stencil_t>
long DiscretizationStencilSolverAssembler<stencil_t>::getColGlobalOffset() const
{
    return m_globalColOffset;
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
    const stencil_t &stencil = getRowStencil(rowIndex);
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
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get pattern
    getPattern(stencil, pattern);
}

/*!
 * Get the values of the specified row.
 *
 * \param stencil is the stencil
 * \param pattern on output will contain the values of the specified stencil
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::getPattern(const stencil_t &stencil, ConstProxyVector<long> *pattern) const
{
    std::size_t stencilSize = stencil.size();

    const long *patternData = stencil.patternData();
    if (m_blockSize == 1) {
        pattern->set(patternData, stencilSize);
    } else {
        std::size_t expandedPatternSize = m_blockSize * stencilSize;
        pattern->set(ConstProxyVector<long>::INTERNAL_STORAGE, expandedPatternSize);
        ConstProxyVector<long>::storage_pointer expandedPatternStorage = pattern->storedData();

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
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get values
    getValues(stencil, values);
}

/*!
 * Get the values of the specified stencil.
 *
 * \param stencil is the stencil
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t>::getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const
{
    values->set(stencil.weightData(), stencil.size());
}

/*!
 * Get the values of the specified stencil.
 *
 * \param stencil is the stencil
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t>::getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const
{
    std::size_t stencilSize = stencil.size();
    const typename stencil_t::weight_type *weightData = stencil.weightData();

    std::size_t expandedValuesSize = m_blockSize * stencilSize;
    values->set(ConstProxyVector<double>::INTERNAL_STORAGE, expandedValuesSize);
    ConstProxyVector<double>::storage_pointer expandedValuesStorage = values->storedData();

    std::size_t expandedValuesIdx = 0;
    for (std::size_t k = 0; k < stencilSize; ++k) {
        for (int i = 0; i < m_blockSize; ++i) {
            expandedValuesStorage[expandedValuesIdx++] = getRawValue(weightData[k], i);
        }
    }
}

/*!
 * Get the data of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the values of the specified row
 * \param values on output will contain the values of the specified row
 */
template<typename stencil_t>
void DiscretizationStencilSolverAssembler<stencil_t>::getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const
{
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get pattern
    getPattern(stencil, pattern);

    // Get values
    getValues(stencil, values);
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
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get constant
    return getConstant(stencil);
}

/*!
 * Get the constant associated with the specified stencil.
 *
 * \param stencil is the stencil
 * \result The constant associated with the specified stencil.
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
double DiscretizationStencilSolverAssembler<stencil_t>::getConstant(const stencil_t &stencil) const
{
    return getRawValue(stencil.getConstant(), 0);
}

/*!
 * Get the constant associated with the specified stencil.
 *
 * \param stencil is the stencil
 * \result The constant associated with the specified stencil.
 */
template<typename stencil_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
double DiscretizationStencilSolverAssembler<stencil_t>::getConstant(const stencil_t &stencil) const
{
    const typename stencil_t::weight_type &stencilConstant = stencil.getConstant();

    double constant = 0.;
    for (int i = 0; i < m_blockSize; ++i) {
        constant += getRawValue(stencilConstant, i);
    }

    return constant;
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilSolverAssembler<stencil_t>::getRowStencil(long rowIndex) const
{
    return m_stencils->at(rowIndex);
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
    : DiscretizationStencilSolver<stencil_t>("", false, debug)
{
}

/*!
* Constuctor
*
* \param transpose if set to true, transposed system will be solved
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
template<typename stencil_t>
DiscretizationStencilSolver<stencil_t>::DiscretizationStencilSolver(bool transpose, bool debug)
    : DiscretizationStencilSolver<stencil_t>("", transpose, debug)
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
    : DiscretizationStencilSolver<stencil_t>(prefix, false, debug)
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
DiscretizationStencilSolver<stencil_t>::DiscretizationStencilSolver(const std::string &prefix, bool transpose, bool debug)
    : SystemSolver(prefix, transpose, debug)
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
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const stencil_container_t &stencils)
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
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::assembly(MPI_Comm communicator, bool partitioned, const stencil_container_t &stencils)
#else
/*!
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
template<typename stencil_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const stencil_container_t &stencils)
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
void DiscretizationStencilSolver<stencil_t>::assembly(const DiscretizationStencilSolverAssembler<stencil_t> &assembler)
{
    assembly(MPI_COMM_SELF, false, static_cast<const StencilSolverAssembler &>(assembler));
}

/*!
* Assembly the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(MPI_Comm communicator, bool partitioned, const DiscretizationStencilSolverAssembler<stencil_t> &assembler)
{
    assembly(communicator, partitioned, static_cast<const StencilSolverAssembler &>(assembler));
}
#else
/*!
* Assembly the stencil solver.
*
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(const DiscretizationStencilSolverAssembler<stencil_t> &assembler)
{
    assembly(static_cast<const StencilSolverAssembler &>(assembler));
}
#endif

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
* Assembly the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param assembler is the solver assembler
*/
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::assembly(MPI_Comm communicator, bool partitioned, const StencilSolverAssembler &assembler)
#else
/*!
* Assembly the stencil solver.
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
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::update(const stencil_container_t &stencils)
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
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::update(const std::vector<long> &rows, const stencil_container_t &stencils)
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
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t>::update(std::size_t nRows, const long *rows, const stencil_container_t &stencils)
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
 * Update the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param nRows is the number of stencils that will be updated
 * \param rows are the rows of the stencils that will be updated,
 * if a null pointer is passed, the stencils that will be updated are the
 * stencils from 0 to (nRows - 1).
 * \param assembler is the solver assembler
 */
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::update(std::size_t nRows, const long *rows,
                                                    const StencilSolverAssembler &assembler)
{
    update(nRows, rows, static_cast<const StencilSolverAssembler &>(assembler));
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
 * \param assembler is the solver assembler
 */
template<typename stencil_t>
void DiscretizationStencilSolver<stencil_t>::update(std::size_t nRows, const long *rows,
                                                    const DiscretizationStencilSolverAssembler<stencil_t> &assembler)
{
    // Update the system
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
