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
    : SystemSolver(debug)
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
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
void StencilScalarSolver::initialize(const std::vector<StencilScalar> &stencils)
{
    initialize(MPI_COMM_SELF, false, stencils);
}
/*!
* Initialize the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param stencils are the stencils
*/
void StencilScalarSolver::initialize(MPI_Comm communicator, bool partitioned, const std::vector<StencilScalar> &stencils)
#else
/*!
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
void StencilScalarSolver::initialize(const std::vector<StencilScalar> &stencils)
#endif
{
    long nDOFs = stencils.size();

    long nNZ = 0;
    for (const StencilScalar &stencil : stencils) {
        nNZ += stencil.size();
    }

#if BITPIT_ENABLE_MPI==1
    SparseMatrix matrix(communicator, partitioned, nDOFs, nDOFs, nNZ);
#else
    SparseMatrix matrix(nDOFs, nDOFs, nNZ);
#endif

    m_constants.resize(nDOFs);

    for (long n = 0; n < nDOFs; ++n) {
        const StencilScalar &stencil = stencils[n];
        matrix.addRow(stencil.size(), stencil.patternData(), stencil.weightData());
        m_constants[n] = stencil.getConstant();
    }

    matrix.assembly();

    SystemSolver::initialize(matrix);
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
    long nDOFs = stencils.size();

    long nNZ = 0;
    for (const StencilScalar &stencil : stencils) {
        nNZ += stencil.size();
    }

#if BITPIT_ENABLE_MPI==1
    SparseMatrix elements(getCommunicator(), isPartitioned(), nDOFs, nDOFs, nNZ);
#else
    SparseMatrix elements(nDOFs, nDOFs, nNZ);
#endif

    for (long n = 0; n < nDOFs; ++n) {
        long row = rows[n];

        const StencilScalar &stencil = stencils[n];
        elements.addRow(stencil.size(), stencil.patternData(), stencil.weightData());
        m_constants[row] = stencil.getConstant();
    }

    elements.assembly();

    SystemSolver::update(rows, elements);
}

/*!
* Solve the system.
*/
void StencilScalarSolver::solve()
{
    // Check if the stencil solver is initialized
    if (!isInitialized()) {
        throw std::runtime_error("Unable to solve the system. The stencil solver is not yet initialized.");
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
