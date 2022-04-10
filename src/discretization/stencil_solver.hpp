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

#ifndef __BITPIT_STENCIL_SOLVER_HPP__
#define __BITPIT_STENCIL_SOLVER_HPP__

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif
#include <algorithm>
#include <vector>

#include "bitpit_LA.hpp"

#include "stencil.hpp"

namespace bitpit {

class StencilSolverAssembler : public SystemMatrixAssembler {

public:
    virtual double getRowConstant(long rowIndex) const = 0;

protected:
    using SystemMatrixAssembler::SystemMatrixAssembler;

};

template<typename stencil_t>
class DiscretizationStencilStorageInterface {

public:
    virtual ~DiscretizationStencilStorageInterface() = default;

    virtual std::size_t size() const = 0;

    virtual const stencil_t & at(long rowIndex) const = 0;
    virtual const stencil_t & rawAt(std::size_t rowRawIndex) const = 0;

    virtual const stencil_t & at(long blockIndex, int componentIdx) const = 0;
    virtual const stencil_t & rawAt(std::size_t blockRawIndex, int componentIdx) const = 0;

protected:
    DiscretizationStencilStorageInterface() = default;

};

template<typename stencil_t, typename stencil_container_t>
class DiscretizationStencilProxyBaseStorage : public DiscretizationStencilStorageInterface<stencil_t> {

protected:
    const stencil_container_t *m_stencils;
    int m_stride;

    DiscretizationStencilProxyBaseStorage(const stencil_container_t *stencils, int stride);

};

template<typename stencil_t, typename stencil_container_t>
class DiscretizationStencilProxyStorage : public DiscretizationStencilProxyBaseStorage<stencil_t, stencil_container_t> {

public:
    DiscretizationStencilProxyStorage(const stencil_container_t *stencils);

    std::size_t size() const override;

    const stencil_t & at(long rowIndex) const override;
    const stencil_t & rawAt(std::size_t rowRawIndex) const override;

    const stencil_t & at(long blockIndex, int componentIdx) const override;
    const stencil_t & rawAt(std::size_t blockRawIndex, int componentIdx) const override;

};

template<typename stencil_t>
class DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>> : public DiscretizationStencilProxyBaseStorage<stencil_t, PiercedStorage<stencil_t>> {

public:
    DiscretizationStencilProxyStorage(const PiercedStorage<stencil_t> *stencils);

    std::size_t size() const override;

    const stencil_t & at(long rowIndex) const override;
    const stencil_t & rawAt(std::size_t rowRawIndex) const override;

    const stencil_t & at(long blockIndex, int componentIdx) const override;
    const stencil_t & rawAt(std::size_t blockRawIndex, int componentIdx) const override;

};

template<typename stencil_t>
class DiscretizationStencilSolverAssembler : public StencilSolverAssembler {

public:
    typedef stencil_t stencil_type;

    template<typename stencil_container_t = std::vector<stencil_t>>
    DiscretizationStencilSolverAssembler(const stencil_container_t *stencils);
#if BITPIT_ENABLE_MPI==1
    template<typename stencil_container_t = std::vector<stencil_t>>
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, const stencil_container_t *stencils);
#endif

    int getBlockSize() const;

    long getRowCount() const override;
    long getColCount() const override;

#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const override;
    long getColGlobalCount() const override;

    long getRowGlobalOffset() const override;
    long getColGlobalOffset() const override;
#endif

    long getRowNZCount(long rowIndex) const override;
    long getMaxRowNZCount() const override;

    void getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const override;
    void getRowValues(long rowIndex, ConstProxyVector<double> *values) const override;
    void getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const override;

    double getRowConstant(long rowIndex) const override;

protected:
    long m_nRows;
    long m_nCols;

#if BITPIT_ENABLE_MPI==1
    long m_nGlobalRows;
    long m_nGlobalCols;

    long m_globalRowOffset;
    long m_globalColOffset;
#endif

    std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> m_stencils;

    int m_blockSize;

    long m_maxRowNZ;

    DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils);
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils);
#endif
    DiscretizationStencilSolverAssembler();

    void setStencils(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils);

#if BITPIT_ENABLE_MPI==1
    void setMatrixSizes(MPI_Comm communicator, bool partitioned);
    void setMatrixSizes(long nRows, long nCols, MPI_Comm communicator, bool partitioned);
#else
    void setMatrixSizes();
    void setMatrixSizes(long nRows, long nCols);
#endif

    void setBlockSize();
    void setBlockSize(int blockSize);

    void setMaximumRowNZ();
    void setMaximumRowNZ(long maxRowNZ);

    virtual const stencil_t & getRowStencil(long rowIndex) const;

    void getPattern(const stencil_t &stencil, ConstProxyVector<long> *pattern) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<std::is_fundamental<U>::value>::type * = nullptr>
    void getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<!std::is_fundamental<U>::value>::type * = nullptr>
    void getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<std::is_fundamental<U>::value>::type * = nullptr>
    double getConstant(const stencil_t &stencil) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<!std::is_fundamental<U>::value>::type * = nullptr>
    double getConstant(const stencil_t &stencil) const;

    double getRawValue(const typename stencil_t::weight_type &weight, int item) const;

};

template<typename stencil_t>
class DiscretizationStencilSolver : public SystemSolver {

public:
    DiscretizationStencilSolver(bool debug = false);
    DiscretizationStencilSolver(bool transpose, bool debug);
    DiscretizationStencilSolver(const std::string &prefix, bool debug = false);
    DiscretizationStencilSolver(const std::string &prefix, bool transpose, bool debug = false);

    void clear(bool release = false);
    template<typename stencil_container_t = std::vector<stencil_t>>
    void assembly(const stencil_container_t &stencils);
    void assembly(const DiscretizationStencilSolverAssembler<stencil_t> &assembler);
    void assembly(const StencilSolverAssembler &assembler);
#if BITPIT_ENABLE_MPI==1
    template<typename stencil_container_t = std::vector<stencil_t>>
    void assembly(MPI_Comm communicator, bool partitioned, const stencil_container_t &stencils);
    void assembly(MPI_Comm communicator, bool partitioned, const DiscretizationStencilSolverAssembler<stencil_t> &assembler);
    void assembly(MPI_Comm communicator, bool partitioned, const StencilSolverAssembler &assembler);
#endif
    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(const stencil_container_t &stencils);
    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(const std::vector<long> &rows, const stencil_container_t &stencils);
    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(std::size_t nRows, const long *rows, const stencil_container_t &stencils);
    void update(std::size_t nRows, const long *rows, const StencilSolverAssembler &assembler);
    void update(std::size_t nRows, const long *rows, const DiscretizationStencilSolverAssembler<stencil_t> &assembler);

    void solve();

    std::size_t getRowStencilCount() const;

protected:
    std::size_t nStencils;

    std::vector<double> m_constants;

};

// Specializations
template<>
void DiscretizationStencilSolverAssembler<StencilScalar>::setBlockSize();

template<>
double DiscretizationStencilSolverAssembler<StencilScalar>::getRawValue(const StencilScalar::weight_type &element, int item) const;

template<>
void DiscretizationStencilSolverAssembler<StencilVector>::setBlockSize();

template<>
double DiscretizationStencilSolverAssembler<StencilVector>::getRawValue(const StencilVector::weight_type &element, int item) const;

template<>
void DiscretizationStencilSolverAssembler<StencilBlock>::setBlockSize();

template<>
double DiscretizationStencilSolverAssembler<StencilBlock>::getRawValue(const StencilBlock::weight_type &element, int item) const;

}

// Template implementation
#include "stencil_solver.tpp"

// Declaration of the typdefs
namespace bitpit {

typedef DiscretizationStencilSolverAssembler<StencilScalar> StencilScalarSolverAssembler;
typedef DiscretizationStencilSolverAssembler<StencilVector> StencilVectorSolverAssembler;
typedef DiscretizationStencilSolverAssembler<StencilBlock> StencilBlockSolverAssembler;

typedef DiscretizationStencilSolver<StencilScalar> StencilScalarSolver;
typedef DiscretizationStencilSolver<StencilVector> StencilVectorSolver;
typedef DiscretizationStencilSolver<StencilBlock> StencilBlockSolver;

}

// Explicit instantization
#ifndef __BITPIT_STENCIL_SOLVER_SRC__
namespace bitpit {

extern template class DiscretizationStencilSolverAssembler<StencilScalar>;
extern template class DiscretizationStencilSolverAssembler<StencilVector>;
extern template class DiscretizationStencilSolverAssembler<StencilBlock>;

extern template class DiscretizationStencilSolver<StencilScalar>;
extern template class DiscretizationStencilSolver<StencilVector>;
extern template class DiscretizationStencilSolver<StencilBlock>;

}
#endif

#endif
