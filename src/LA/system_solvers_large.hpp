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

#ifndef __BITPIT_SYSTEM_SOLVERS_LARGE_HPP__
#define __BITPIT_SYSTEM_SOLVERS_LARGE_HPP__

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <petscksp.h>

#include "system_matrix.hpp"

namespace bitpit {

struct KSPOptions {
    PetscInt restart;
    PetscInt maxits;
    PetscScalar rtol;

    PetscInt overlap;
    PetscInt levels;

    PetscInt sublevels;
    PetscScalar subrtol;

    KSPOptions()
        : restart(PETSC_DEFAULT), maxits(PETSC_DEFAULT), rtol(PETSC_DEFAULT),
          overlap(PETSC_DEFAULT), levels(PETSC_DEFAULT),
          sublevels(PETSC_DEFAULT), subrtol(PETSC_DEFAULT)
    {
    }
};

struct KSPStatus {
    PetscErrorCode error;
    PetscInt its;
    KSPConvergedReason convergence;

    KSPStatus()
        : error(0), its(-1), convergence(KSP_DIVERGED_BREAKDOWN)
    {
    }
};

class SystemSolver {

public:
    enum PivotType {
        PIVOT_NONE,  // Natural
        PIVOT_ND,    // Nested Dissection
        PIVOT_1WD,   // One-way Dissection
        PIVOT_RCM,   // Reverse Cuthill-McKee
        PIVOT_MD     // Quotient Minimum Degree
    };

    static void addInitOption(const std::string &option);
    static void addInitOptions(const std::vector<std::string> &options);

    SystemSolver(bool debug = false);
#if BITPIT_ENABLE_MPI==1
    SystemSolver(MPI_Comm communicator, bool debug = false);
#endif

    virtual ~SystemSolver();

    void clear();
    void initialize(const SparseMatrix &matrix, PivotType pivotType = PIVOT_NONE);
    bool isInitialized() const;

    void update(const std::vector<long> &rows, const SparseMatrix &elements);

    long getRowCount() const;
    long getColCount() const;
#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const;
    long getColGlobalCount() const;
#endif

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    void dump(const std::string &directory, const std::string &prefix = "") const;

    PivotType getPivotType();

    void setNullSpace();
    void unsetNullSpace();

    KSPOptions & getKSPOptions();
    const KSPOptions & getKSPOptions() const;
    const KSPStatus & getKSPStatus() const;

    double * getRHSRawPtr();
    const double * getRHSRawPtr() const;
    const double * getRHSRawReadPtr() const;
    void restoreRHSRawPtr(double *raw_rhs);
    void restoreRHSRawReadPtr(const double *raw_rhs) const;

    double * getSolutionRawPtr();
    const double * getSolutionRawPtr() const;
    const double * getSolutionRawReadPtr() const;
    void restoreSolutionRawPtr(double *raw_solution);
    void restoreSolutionRawReadPtr(const double *raw_solution) const;

protected:
    void matrixInit(const SparseMatrix &matrix);
    void matrixFill(const SparseMatrix &matrix);
    void matrixUpdate(const std::vector<long> &rows, const SparseMatrix &elements);
    void matrixReorder();

#if BITPIT_ENABLE_MPI == 1
    void vectorsInit(const std::vector<long> &ghosts);
#else
    void vectorsInit();
#endif
    void vectorsReorder(PetscBool inv);
    void vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution);
    void vectorsExport(std::vector<double> *solution);

private:
    static int m_nInstances;
    static std::vector<std::string> m_options;

    bool m_initialized;
    PivotType m_pivotType;

    long m_rowGlobalOffset;
    long m_colGlobalOffset;

#if BITPIT_ENABLE_MPI==1
    bool m_partitioned;
    MPI_Comm m_communicator;
#endif

    Mat m_A;
    Vec m_rhs;
    Vec m_solution;

    IS m_rpivot;
    IS m_cpivot;

    KSP m_KSP;
    KSPOptions m_KSPOptions;
    KSPStatus m_KSPStatus;

    void pivotInit(PivotType pivotType);

    void KSPInit();

};

}

#endif
