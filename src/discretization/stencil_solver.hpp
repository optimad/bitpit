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
class DiscretizationStencilSolverAssembler : public StencilSolverAssembler {

public:
    DiscretizationStencilSolverAssembler(const std::vector<stencil_t> *stencils);
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned, const std::vector<stencil_t> *stencils);
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
    double getRowConstant(long rowIndex) const override;

protected:
    const std::vector<stencil_t> *m_stencils;

    int m_blockSize;

    long m_nDOFs;
    long m_maxRowNZ;

    void initializeBlockSize();

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<std::is_fundamental<U>::value>::type * = nullptr>
    void _getRowValues(long rowIndex, ConstProxyVector<double> *values) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<!std::is_fundamental<U>::value>::type * = nullptr>
    void _getRowValues(long rowIndex, ConstProxyVector<double> *values) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<std::is_fundamental<U>::value>::type * = nullptr>
    double _getRowConstant(long rowIndex) const;

    template<typename U = typename stencil_t::weight_type, typename std::enable_if<!std::is_fundamental<U>::value>::type * = nullptr>
    double _getRowConstant(long rowIndex) const;

    double getRawValue(const typename stencil_t::weight_type &weight, int item) const;

#if BITPIT_ENABLE_MPI==1
    long m_nGlobalDOFs;
    long m_globalDOFOffset;
#endif


};

template<typename stencil_t>
class DiscretizationStencilSolver : public SystemSolver {

public:
    DiscretizationStencilSolver(bool debug = false);
    DiscretizationStencilSolver(const std::string &prefix, bool debug = false);

    void clear(bool release = false);
    void assembly(const std::vector<stencil_t> &stencils);
    void assembly(const StencilSolverAssembler &assembler);
#if BITPIT_ENABLE_MPI==1
    void assembly(MPI_Comm communicator, bool partitioned, const std::vector<stencil_t> &stencils);
    void assembly(MPI_Comm communicator, bool partitioned, const StencilSolverAssembler &assembler);
#endif
    void update(const std::vector<stencil_t> &stencils);
    void update(const std::vector<long> &rows, const std::vector<stencil_t> &stencils);
    void update(std::size_t nRows, const long *rows, const std::vector<stencil_t> &stencils);

    void solve();

    std::size_t getStencilCount() const;

protected:
    std::size_t nStencils;

    std::vector<double> m_constants;

};

// Specializations
template<>
void DiscretizationStencilSolverAssembler<StencilScalar>::initializeBlockSize();

template<>
double DiscretizationStencilSolverAssembler<StencilScalar>::getRawValue(const StencilScalar::weight_type &element, int item) const;

template<>
void DiscretizationStencilSolverAssembler<StencilVector>::initializeBlockSize();

template<>
double DiscretizationStencilSolverAssembler<StencilVector>::getRawValue(const StencilVector::weight_type &element, int item) const;

template<>
void DiscretizationStencilSolverAssembler<StencilBlock>::initializeBlockSize();

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
