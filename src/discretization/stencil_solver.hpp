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

#ifndef __BITPIT_STENCIL_SCALAR_SOLVER_HPP__
#define __BITPIT_STENCIL_SCALAR_SOLVER_HPP__

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif
#include <vector>

#include "bitpit_LA.hpp"

#include "stencil.hpp"

namespace bitpit {

class StencilSolverAssembler : public SystemMatrixAssembler {

public:
    StencilSolverAssembler(const std::vector<StencilScalar> *stencils);
#if BITPIT_ENABLE_MPI==1
    StencilSolverAssembler(MPI_Comm communicator, bool partitioned, const std::vector<StencilScalar> *stencils);
#endif

    long getRowCount() const override;
    long getColCount() const override;

#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const override;
    long getColGlobalCount() const override;

    long getRowGlobalOffset() const override;
    long getColGlobalOffset() const override;
#endif

    long getRowNZCount(long row) const override;
    long getMaxRowNZCount() const override;

    void getRowPattern(long row, ConstProxyVector<long> *pattern) const override;
    void getRowValues(long row, ConstProxyVector<double> *values) const override;

    double getRowConstant(long row) const;

protected:
    const std::vector<StencilScalar> *m_stencils;

    long m_nDOFs;
    long m_maxRowNZ;

#if BITPIT_ENABLE_MPI==1
    long m_nGlobalDOFs;
    long m_globalDOFOffset;
#endif


};

class StencilScalarSolver : public SystemSolver {

public:
    StencilScalarSolver(bool debug = false);
    StencilScalarSolver(const std::string &prefix, bool debug = false);

    void clear(bool release = false);
    void assembly(const std::vector<StencilScalar> &stencils);
#if BITPIT_ENABLE_MPI==1
    void assembly(MPI_Comm communicator, bool partitioned, const std::vector<StencilScalar> &stencils);
#endif
    void update(const std::vector<long> &rows, const std::vector<StencilScalar> &stencils);

    void solve();

protected:
    std::vector<double> m_constants;

};

}

#endif
