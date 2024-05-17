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

template<typename stencil_t, typename solver_kernel_t = SystemSolver>
class DiscretizationStencilSolverAssembler : public solver_kernel_t::Assembler {

public:
    using stencil_type = stencil_t;

    using solver_kernel_type    = solver_kernel_t;
    using assembly_type         = typename solver_kernel_t::Assembler;
    using assembly_options_type = typename assembly_type::AssemblyOptions;

    template<typename stencil_container_t, typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(const stencil_container_t *stencils, AssemblerKernelArgs&&... assemblerKernelArgs);
#if BITPIT_ENABLE_MPI==1
    template<typename stencil_container_t, typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, const stencil_container_t *stencils,
                                         AssemblerKernelArgs&&... assemblerKernelArgs);
#endif

#if BITPIT_ENABLE_MPI==1
    bool isPartitioned() const override;
    const MPI_Comm & getCommunicator() const override;
#endif

    assembly_options_type getOptions() const override;

    int getBlockSize() const override;

    long getRowCount() const override;
    long getColCount() const override;

    long getRowElementCount() const override;
    long getColElementCount() const override;

#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const override;
    long getColGlobalCount() const override;

    long getRowGlobalElementCount() const override;
    long getColGlobalElementCount() const override;

    long getRowGlobalOffset() const override;
    long getColGlobalOffset() const override;

    long getRowGlobalElementOffset() const override;
    long getColGlobalElementOffset() const override;
#endif

    long getRowNZCount(long rowIndex) const override;
    long getMaxRowNZCount() const override;

    void getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const override;
    void getRowValues(long rowIndex, ConstProxyVector<double> *values) const override;
    void getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const override;

    virtual void getRowConstant(long rowIndex, bitpit::ConstProxyVector<double> *constant) const;

protected:
    using stencil_weight_type = typename stencil_type::weight_type;
    using stencil_value_type  = typename stencil_type::value_type;

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

    template<typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils,
                                         AssemblerKernelArgs&&... assemblerKernelArgs);
    template<typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(AssemblerKernelArgs&&... assemblerKernelArgs);
#if BITPIT_ENABLE_MPI==1
    template<typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned,
                                         std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils,
                                         AssemblerKernelArgs&&... assemblerKernelArgs);
    template<typename... AssemblerKernelArgs>
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned,
                                         AssemblerKernelArgs&&... assemblerKernelArgs);
#endif

    void setStencils(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils);

    void setMatrixSizes();
    void setMatrixSizes(long nRows, long nCols);

    template<typename W = stencil_weight_type, typename V = stencil_value_type, typename std::enable_if<std::is_fundamental<W>::value>::type * = nullptr>
    void setBlockSize();
    template<typename W = stencil_weight_type, typename V = stencil_value_type, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    void setBlockSize();
    template<typename W = stencil_weight_type, typename V = stencil_value_type, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    void setBlockSize();
    void setBlockSize(int blockSize);

    void setMaximumRowNZ();
    void setMaximumRowNZ(long maxRowNZ);

    virtual const stencil_t & getRowStencil(long rowIndex) const;

    void getPattern(const stencil_t &stencil, ConstProxyVector<long> *pattern) const;

    template<typename W = stencil_weight_type, typename std::enable_if<std::is_fundamental<W>::value>::type * = nullptr>
    void getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const;

    template<typename W = stencil_weight_type, typename std::enable_if<!std::is_fundamental<W>::value>::type * = nullptr>
    void getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const;

    template<typename W = stencil_weight_type, typename std::enable_if<std::is_fundamental<W>::value>::type * = nullptr>
    void getConstant(const stencil_t &stencil, bitpit::ConstProxyVector<double> *constant) const;

    template<typename W = stencil_weight_type, typename std::enable_if<!std::is_fundamental<W>::value>::type * = nullptr>
    void getConstant(const stencil_t &stencil, bitpit::ConstProxyVector<double> *constant) const;

private:
#if BITPIT_ENABLE_MPI==1
    bool m_partitioned;
    MPI_Comm m_communicator;
#endif

};

template<typename stencil_t, typename solver_kernel_t = SystemSolver>
class DiscretizationStencilSolver : public solver_kernel_t {

public:
    typedef DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t> Assembler;

    using solver_kernel_t::solver_kernel_t;

    void clear(bool release = false);

    template<typename stencil_container_t = std::vector<stencil_t>>
    void assembly(const stencil_container_t &stencils);
#if BITPIT_ENABLE_MPI==1
    template<typename stencil_container_t = std::vector<stencil_t>>
    void assembly(MPI_Comm communicator, bool partitioned, const stencil_container_t &stencils);
#endif
    void assembly(const Assembler &assembler);

    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(const stencil_container_t &stencils);
    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(const std::vector<long> &rows, const stencil_container_t &stencils);
    template<typename stencil_container_t = std::vector<stencil_t>>
    void update(std::size_t nRows, const long *rows, const stencil_container_t &stencils);
    void update(long nRows, const long *rows, const Assembler &assembler);

    void solve();

    void matrixAssembly(const Assembler &assembler);
    void matrixUpdate(long nRows, const long *rows, const Assembler &assembler);

protected:
    std::vector<double> m_constants;

    using solver_kernel_t::assembly;
    using solver_kernel_t::update;

    void assembleConstants(const Assembler &assembler);
    void updateConstants(std::size_t nRows, const long *rows, const Assembler &assembler);

};

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
