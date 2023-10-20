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

#ifndef __BITPIT_SYSTEM_MATRIX_HPP__
#define __BITPIT_SYSTEM_MATRIX_HPP__

#include <vector>

#if BITPIT_ENABLE_MPI==1
#   include <mpi.h>
#endif

#include "bitpit_containers.hpp"

namespace bitpit {

class SparseMatrix {

public:
    SparseMatrix();
    SparseMatrix(long nRows, long nCols, long nNZ);
    SparseMatrix(int blockSize, long nRows, long nCols, long nNZ);
#if BITPIT_ENABLE_MPI==1
    SparseMatrix(MPI_Comm communicator);
    SparseMatrix(MPI_Comm communicator, bool partitioned, long nRows, long nCols, long nNZ);
    SparseMatrix(MPI_Comm communicator, bool partitioned, int blockSize, long nRows, long nCols, long nNZ);
#endif

    ~SparseMatrix();

    void initialize(long nRows, long nCols, long nNZ);
    void initialize(int blockSize, long nRows, long nCols, long nNZ);
#if BITPIT_ENABLE_MPI==1
    void initialize(bool partitioned, long nRows, long nCols, long nNZ);
    void initialize(bool partitioned, int blockSize, long nRows, long nCols, long nNZ);
#endif
    void clear(bool release = false);
    void squeeze();

    virtual void display(std::ostream &stream, double negligiblity, int indent = 0) const;

    void assembly();
    bool isAssembled() const;

    long countMissingRows() const;
    long countAddedRows() const;

    int getBlockSize() const;

    long getRowCount() const;
    long getColCount() const;

    long getRowElementCount() const;
    long getColElementCount() const;

    long getNZCount() const;
    long getRowNZCount(long row) const;
    long getMaxRowNZCount() const;

    long getNZElementCount() const;
    long getRowNZElementCount(long row) const;
    long getMaxRowNZElementCount() const;

#if BITPIT_ENABLE_MPI==1
    bool isPartitioned() const;

    const MPI_Comm & getCommunicator() const;

    long getRowGlobalCount() const;
    long getRowGlobalOffset() const;

    long getRowGlobalElementCount() const;
    long getRowGlobalElementOffset() const;

    long getColGlobalCount() const;
    long getColGlobalOffset() const;

    long getColGlobalElementCount() const;
    long getColGlobalElementOffset() const;

    long getNZGlobalCount() const;
    long getMaxRowNZGlobalCount() const;

    long getNZGlobalElementCount() const;
    long getMaxRowNZGlobalElementCount() const;

    std::vector<long> extractLocalGlobalRows() const;
    std::vector<long> extractGhostGlobalRows() const;

    std::vector<long> extractLocalGlobalCols() const;
    std::vector<long> extractGhostGlobalCols() const;
#endif

    void addRow(const std::vector<long> &rowPattern, const std::vector<double> &rowValues);
    void addRow(long nRowNZ, const long *rowPattern, const double *rowValues);

    ConstProxyVector<long> getRowPattern(long row) const;
    void getRowPattern(long row, ConstProxyVector<long> *pattern) const;

    ConstProxyVector<double> getRowValues(long row) const;
    void getRowValues(long row, ConstProxyVector<double> *values) const;

    std::unique_ptr<SparseMatrix> computeTranspose() const;

protected:
    int m_blockSize;

    long m_nRows;
    long m_nCols;
    long m_nNZ;
    long m_maxRowNZ;

    long m_lastRow;

    bool m_assembled;

#if BITPIT_ENABLE_MPI==1
    bool m_partitioned;
    MPI_Comm m_communicator;

    long m_global_nRows;
    long m_global_nCols;
    long m_global_nNZ;
    long m_global_maxRowNZ;
    long m_global_rowOffset;
    long m_global_colOffset;
#endif

    FlatVector2D<long> m_pattern;
    std::vector<double> m_values;

    long * getRowPatternData(long row);
    const long * getRowPatternData(long row) const;

    double * getRowValuesData(long row);
    const double * getRowValuesData(long row) const;

    void initializePatternStorage();
    void squeezePatternStorage();
    void clearPatternStorage(bool release);

    void initializeValueStorage();
    void squeezeValueStorage();
    void clearValueStorage(bool release);

private:
    void _initialize(int blockSize, long nRows, long nCols, long nNZ);
#if BITPIT_ENABLE_MPI==1
    void _initialize(bool partitioned, int bBlocks, long nRows, long nCols, long nNZ);
#endif

#if BITPIT_ENABLE_MPI==1
    void setCommunicator(MPI_Comm communicator);
    void freeCommunicator();
#endif

    void initializePatternStorage(long capacity);
    void initializeValueStorage(long capacity);

    long getNZElementCount(long nNZ) const;

};

}

#endif
