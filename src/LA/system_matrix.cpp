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

#include "system_matrix.hpp"

namespace bitpit {

/**
* \class SparseMatrix
* \ingroup system_solver
*
* \brief Sparse matrix.
*
* m,n,M,N parameters specify the size of the matrix, and its partitioning across processors, while
* d_nz,d_nnz,o_nz,o_nnz parameters specify the approximate storage requirements for this matrix.
*
*/
#if BITPIT_ENABLE_MPI==1
/**
* The SparseMatrix class mimics the beahviour of the sparse parallel matrix
* in AIJ format of the PETSc library.
*
* The parallel matrix is partitioned across processors such that the first
* m0 rows belong to process 0, the next m1 rows belong to process 1, the
* next m2 rows belong to process 2 etc.. where m0, m1, m2,.. are the local
* number of rows defined by each processors.
*
* Each processor stores values corresponding to [m x N] submatrix, wher N
* is the global number of columns of the matrix, N is obtained summing the
* local number of columns of each processor.
*
* Although each processor stores values for the whole global row, columns are
* logically partitioned across processors, with the n0 columns belonging to
* 0th partition, the next n1 columns belonging to the next partition etc..
* where n0, n1, n2... are the local number of columns defined by each
* processor.
*
* The DIAGONAL portion of the local submatrix on any given processor is the
* submatrix corresponding to the rows and columns m,n corresponding to the
* given processor. i.e diagonal matrix on process 0 is [m0 x n0], diagonal
* matrix on process 1 is [m1 x n1] etc. The remaining portion of the local
* submatrix [m x (N-n)] constitute the OFF-DIAGONAL portion.
*
* For a square global matrix we define each processor's diagonal portion to
* be its local rows and the corresponding columns (a square submatrix); each
* processor's off-diagonal portion encompasses the remainder of the local
* matrix (a rectangular submatrix).
*/
#endif

#if BITPIT_ENABLE_MPI==1
/**
* Default constructor
*/
SparseMatrix::SparseMatrix()
    : SparseMatrix(MPI_COMM_SELF)
{
}

/**
* Constructor
*
* \param communicator is the MPI communicator
*/
SparseMatrix::SparseMatrix(MPI_Comm communicator)
    : m_communicator(communicator),
      m_nRows(0), m_nCols(0), m_nNZ(0), m_maxRowNZ(0), m_lastRow(-1),
      m_assembled(false),
      m_global_nRows(0), m_global_nCols(0), m_global_nNZ(0),
      m_global_maxRowNZ(0), m_global_rowOffset(0), m_global_colOffset(0)
{
    // Detect if the matrix is partitioned
    m_partitioned = (communicator != MPI_COMM_NULL);

    // Set the communicator
    if (m_partitioned) {
        MPI_Comm_dup(communicator, &m_communicator);
    }
}
#else
/**
* Default constructor
*/
SparseMatrix::SparseMatrix()
    : m_nRows(0), m_nCols(0), m_nNZ(0), m_maxRowNZ(0), m_lastRow(-1),
      m_assembled(false)
{
}
#endif

#if BITPIT_ENABLE_MPI==1
/**
* Constructor
*
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
SparseMatrix::SparseMatrix(long nRows, long nCols, long nNZ)
    : SparseMatrix(MPI_COMM_SELF, nRows, nCols, nNZ)
{
}

/**
* Constructor
*
* \param communicator is the MPI communicator
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
SparseMatrix::SparseMatrix(MPI_Comm communicator, long nRows, long nCols, long nNZ)
    : SparseMatrix(communicator)
{
    _initialize(nRows, nCols, nNZ);
}
#else
/**
* Constructor
*
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
SparseMatrix::SparseMatrix(long nRows, long nCols, long nNZ)
    : SparseMatrix()
{
    _initialize(nRows, nCols, nNZ);
}
#endif

/*!
 * Destructor
 */
SparseMatrix::~SparseMatrix()
{
#if ENABLE_MPI==1
    // Free the MPI communicator
    if (m_partitioned) {
        int finalizedCalled;
        MPI_Finalized(&finalizedCalled);
        if (!finalizedCalled) {
            MPI_Comm_free(&m_communicator);
        }
    }
#endif
}

/**
* Initialize the pattern.
*
*/
#if BITPIT_ENABLE_MPI==1
/*!
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
*/
#else
/*!
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
*/
#endif
/*!
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
void SparseMatrix::initialize(long nRows, long nCols, long nNZ)
{
    _initialize(nRows, nCols, nNZ);
}

/**
* Internal function to initialize the pattern.
*
*/
#if BITPIT_ENABLE_MPI==1
/*!
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
*/
#else
/*!
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
*/
#endif
/*!
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
void SparseMatrix::_initialize(long nRows, long nCols, long nNZ)
{
    assert(nRows >= 0);
    assert(nCols >= 0);
    assert(nNZ >= 0);

    clear();

    m_pattern.reserve(nRows, nNZ);
    m_values.reserve(nNZ);

    m_nRows = nRows;
    m_nCols = nCols;
#if BITPIT_ENABLE_MPI == 1
    if (m_partitioned) {
        MPI_Allreduce(&m_nRows, &m_global_nRows, 1, MPI_LONG, MPI_SUM, m_communicator);
        MPI_Allreduce(&m_nCols, &m_global_nCols, 1, MPI_LONG, MPI_SUM, m_communicator);
    }

    m_global_rowOffset = 0;
    m_global_colOffset = 0;
    if (m_partitioned) {
        int nProcessors = 1;
        MPI_Comm_size(m_communicator, &nProcessors);

        std::vector<long> nGlobalRows(nProcessors);
        MPI_Allgather(&nRows, 1, MPI_LONG, nGlobalRows.data(), 1, MPI_LONG, m_communicator);

        std::vector<long> nGlobalCols(nProcessors);
        MPI_Allgather(&nCols, 1, MPI_LONG, nGlobalCols.data(), 1, MPI_LONG, m_communicator);

        int rank;
        MPI_Comm_rank(m_communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_global_rowOffset += nGlobalRows[i];
            m_global_colOffset += nGlobalCols[i];
        }
    }
#endif
}


/**
* Clear the pattern.
*
* \param release if set to true the memory hold by the pattern will be released
*/
void SparseMatrix::clear(bool release)
{
    if (release) {
        m_pattern.clear();
        m_values.clear();
    } else {
        FlatVector2D<long>().swap(m_pattern);
        std::vector<double>().swap(m_values);
    }

    m_nRows    = -1;
    m_nCols    = -1;
    m_nNZ      = -1;
    m_maxRowNZ = -1;
    m_lastRow  = -1;

#if BITPIT_ENABLE_MPI==1
    m_global_nRows     = 0;
    m_global_nCols     = 0;
    m_global_nNZ       = 0;
    m_global_maxRowNZ  = 0;
    m_global_rowOffset = 0;
    m_global_colOffset = 0;
#endif

    m_assembled = false;
}

/*!
* Squeeze.
*
* Requests the matrix pattern to reduce its capacity to fit its size.
*/
void SparseMatrix::squeeze()
{
    m_pattern.shrinkToFit();
    m_values.shrink_to_fit();
}

/*!
* Assembly the matrix.
*
* This function should be called after adding all the rows of the matrix.
* Its purpose is to prepare the matrix for the usage.
#if BITPIT_ENABLE_MPI==1
*
* It's a collective operation, hence it has to be called by all processors.
#endif
*/
void SparseMatrix::assembly()
{
    // Early return if the matrix is already assembled.
    if (isAssembled()) {
        return;
    }

    // Assembly can be called only after adding all the rows
    if (countMissingRows() != 0) {
        throw std::runtime_error("Assembly can be called only after adding all the rows.");
    }

#if BITPIT_ENABLE_MPI==1
    // Updathe global information of the non-zero elements
    if (m_partitioned) {
        MPI_Allreduce(&m_maxRowNZ, &m_global_maxRowNZ, 1, MPI_LONG, MPI_MAX, m_communicator);
        MPI_Allreduce(&m_nNZ, &m_global_nNZ, 1, MPI_LONG, MPI_SUM, m_communicator);
    }
#endif

    // The function is now assembled
    m_assembled = true;
}

/**
* Check if the matrix is assembled and ready for use.

* \result Returns true if the matrix is assembled and ready for use.
*/
bool SparseMatrix::isAssembled() const
{
    return m_assembled;
}

/**
* Count the number of rows that needs to be added to fill the matrix.
*
* \result The number of rows that needs to be added to fill the matrix.
*/
long SparseMatrix::countMissingRows() const
{
    long nRows      = getRowCount();
    long nAddedRows = countAddedRows();

    return (nRows - nAddedRows);
}

/**
* Count the number of rows that have been added to the matrix.
*
* \result The number of rows that have been added to the matrix.
*/
long SparseMatrix::countAddedRows() const
{
    return (m_lastRow + 1);
}

/**
* Get the number of rows of the matrix.
*
* \result The number of rows of the matrix.
*/
long SparseMatrix::getRowCount() const
{
    return m_nRows;
}

/**
* Get the number of columns of the matrix.
*
* \result The number of columns of the matrix.
*/
long SparseMatrix::getColCount() const
{
    return m_nCols;
}

/**
* Get the number of non-zero elements.
*
* \result The number of non-zero elements.
*/
long SparseMatrix::getNZCount() const
{
    return m_nNZ;
}

/**
* Get the number of non-zero elements in the specified row.
*
* \result The number of non-zero elements in the specified row.
*/
long SparseMatrix::getRowNZCount(long row) const
{
    return m_pattern.getItemCount(row);
}

/**
* Get the maximum number of non-zero elements per row.
*
* \result The maximum number of non-zero elements per row.
*/
long SparseMatrix::getMaxRowNZCount() const
{
    return m_maxRowNZ;
}

#if BITPIT_ENABLE_MPI==1
/**
* Get the number of global rows
*
* \result The number of global rows
*/
long SparseMatrix::getRowGlobalCount() const
{
    return m_global_nRows;
}

/**
* Get global row offset.
*
* \result The global row offset.
*/
long SparseMatrix::getRowGlobalOffset() const
{
    return m_global_rowOffset;
}

/**
* Get number of global columns.
*
* \result The number of global columns.
*/
long SparseMatrix::getColGlobalCount() const
{
    return m_global_nCols;
}

/**
* Get global column offset.
*
* \result The global column offset.
*/
long SparseMatrix::getColGlobalOffset() const
{
    return m_global_colOffset;
}

/**
* Get the global number of non-zero elements.
*
* \result The global number of non-zero elements.
*/
long SparseMatrix::getNZGlobalCount() const
{
    return m_global_nNZ;
}

/**
* Get the global maximum number of non-zero elements per row.
*
* \result The global maximum number of non-zero elements per row.
*/
long SparseMatrix::getMaxRowNZGlobalCount() const
{
    return m_global_maxRowNZ;
}

/**
* Extract the list of local global rows.
*
* \result The the list of local global rows.
*/
std::vector<long> SparseMatrix::extractLocalGlobalRows() const
{
    const std::size_t *rowExtents = m_pattern.indices();

    std::vector<long> localGlobalRows;
    localGlobalRows.reserve(m_nRows);
    for (long i = 0; i < m_nRows; ++i) {
        std::size_t nRowNZ = rowExtents[i + 1] - rowExtents[i];
        if (nRowNZ > 0) {
            localGlobalRows.push_back(m_global_rowOffset + i);
        }
    }

    return localGlobalRows;
}

/**
* Extract the list of ghost global rows.
*
* \result The the list of ghost global rows.
*/
std::vector<long> SparseMatrix::extractGhostGlobalRows() const
{
    return std::vector<long>();
}

/**
* Extract the list of local global columns.
*
* \result The the list of local global columns.
*/
std::vector<long> SparseMatrix::extractLocalGlobalCols() const
{
    long firstGlobalCol = m_global_colOffset;
    long lastGlobalCol  = m_global_colOffset + m_nCols - 1;

    const long *globalCols = m_pattern.data();

    std::vector<long> localGlobalCols;
    localGlobalCols.reserve(m_nNZ);
    for (long k = 0; k < m_nNZ; ++k) {
        long globalCol = globalCols[k];
        if (globalCol < firstGlobalCol || globalCol >= lastGlobalCol) {
            utils::addToOrderedVector<long>(globalCol, localGlobalCols);
        }
    }

    return localGlobalCols;
}

/**
* Extract the list of ghost global columns.
*
* \result The the list of ghost global columns.
*/
std::vector<long> SparseMatrix::extractGhostGlobalCols() const
{
    long firstGlobalCol = m_global_colOffset;
    long lastGlobalCol  = m_global_colOffset + m_nCols - 1;

    const long *globalCols = m_pattern.data();

    std::vector<long> ghostGlobalCols;
    ghostGlobalCols.reserve(m_nNZ);
    for (long k = 0; k < m_nNZ; ++k) {
        long globalCol = globalCols[k];
        if (globalCol >= firstGlobalCol || globalCol < lastGlobalCol) {
            utils::addToOrderedVector<long>(globalCol, ghostGlobalCols);
        }
    }

    return ghostGlobalCols;
}
#endif

/**
* Add a row.
*
* \param rowPattern are the indexes of the non-zero columns of the matrix
* \param rowValues are the values of the non-zero columns of the matrix
*/
void SparseMatrix::addRow(const std::vector<long> &rowPattern, const std::vector<double> &rowValues)
{
    addRow(rowPattern.size(), rowPattern.data(), rowValues.data());
}

/**
* Add a row.
*
* \param nRowNZ is the number of non-zero elements in the row
* \param rowPattern are the indexes of the non-zero columns of the matrix
* \param rowValues are the values of the non-zero columns of the matrix
*/
void SparseMatrix::addRow(long nRowNZ, const long *rowPattern, const double *rowValues)
{
    if (countMissingRows() == 0) {
        throw std::runtime_error("Unable to add another row: all rows have already been defined.");
    }

    // Add the row pattern
    m_pattern.pushBack(nRowNZ, rowPattern);

    // Add row values
    m_values.insert(m_values.end(), rowValues, rowValues + nRowNZ);

    // Update the non-zero element counters
    m_maxRowNZ = std::max(nRowNZ, m_maxRowNZ);
    m_nNZ += nRowNZ;

    // Update the index of the last row
    m_lastRow++;
}

/**
* Get the pattern of the specified row.
*
* \param row is the row
* \result The pattern of the specified row.
*/
ConstProxyVector<long> SparseMatrix::getRowPattern(long row) const
{
    const std::size_t *rowExtent = m_pattern.indices(row);

    const long *rowPattern = m_pattern.data() + rowExtent[0];
    const std::size_t rowPatternSize = rowExtent[1] - rowExtent[0];

    return ConstProxyVector<long>(rowPattern, rowPatternSize);
}

/**
* Get the values of the specified row.
*
* \param row is the row
* \result The values of the specified row.
*/
ConstProxyVector<double> SparseMatrix::getRowValues(long row) const
{
    const std::size_t *rowExtent = m_pattern.indices(row);

    const double *rowValues = m_values.data() + rowExtent[0];
    const std::size_t nRowValues = rowExtent[1] - rowExtent[0];

    return ConstProxyVector<double>(rowValues, nRowValues);
}

}
