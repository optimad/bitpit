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
    PetscScalar atol;

    PetscInt overlap;
    PetscInt levels;

    PetscInt sublevels;
    PetscScalar subrtol;

    KSPOptions()
        : restart(PETSC_DEFAULT), maxits(PETSC_DEFAULT),
          rtol(PETSC_DEFAULT), atol(PETSC_DEFAULT),
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

class SystemMatrixAssembler {

public:
    virtual long getRowCount() const = 0;
    virtual long getColCount() const = 0;

#if BITPIT_ENABLE_MPI==1
    virtual long getRowGlobalCount() const = 0;
    virtual long getColGlobalCount() const = 0;

    virtual long getRowGlobalOffset() const = 0;
    virtual long getColGlobalOffset() const = 0;
#endif

    virtual long getRowNZCount(long rowIndex) const = 0;
    virtual long getMaxRowNZCount() const = 0;

    virtual void getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const = 0;
    virtual void getRowValues(long rowIndex, ConstProxyVector<double> *values) const = 0;

protected:
    SystemMatrixAssembler() = default;

};

class SystemSparseMatrixAssembler : public SystemMatrixAssembler {

public:
    SystemSparseMatrixAssembler(const SparseMatrix *matrix);

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

protected:
    const SparseMatrix *m_matrix;

};

class PetscManager {

public:
    PetscManager();

    virtual ~PetscManager();

    bool areOptionsEditable() const;

    bool initialize(bool debug);
    bool finalize(bool permanent);

    void addInitOption(const std::string &option);
    void addInitOptions(int argc, char **args);
    void addInitOptions(const std::vector<std::string> &options);
    void clearInitOptions();

private:
    static PetscErrorCode displayLogView();

    bool m_externalMPIInitialization;
    bool m_externalPETScInitialization;

    bool m_optionsEditable;
    std::vector<std::string> m_options;

    bool m_logViewEnabled;

    void enableLogView();

};

class SystemSolver {

public:
    enum DumpFormat {
        DUMP_BINARY,
        DUMP_ASCII
    };

    SystemSolver(bool debug = false);
    SystemSolver(bool transpose, bool debug);
    SystemSolver(const std::string &prefix, bool debug = false);
    SystemSolver(const std::string &prefix, bool transpose, bool debug);
#if BITPIT_ENABLE_MPI==1
    SystemSolver(MPI_Comm communicator, const std::string matrixPath, const std::string &prefix, bool debug);
#else
    SystemSolver(const std::string matrixFilename, const std::string &prefix, bool debug);
#endif

    virtual ~SystemSolver();

    void clear();

    void setPermutations(long nRows, const long *rowRanks, long nCols, const long *colRanks);

    void assembly(const SparseMatrix &matrix);
    void assembly(const SystemMatrixAssembler &assembler);
#if BITPIT_ENABLE_MPI==1
    void assembly(MPI_Comm communicator, bool isPartitioned, const SystemMatrixAssembler &assembler);
#endif
    bool isAssembled() const;

    void update(const SparseMatrix &elements);
    void update(long nRows, const long *rows, const SparseMatrix &elements);
    void update(const SystemMatrixAssembler &assembler);
    void update(long nRows, const long *rows, const SystemMatrixAssembler &assembler);

    void setUp();
    bool isSetUp() const;

    long getRowCount() const;
    long getColCount() const;
#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const;
    long getColGlobalCount() const;

    bool isPartitioned() const;
#endif

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    void dump(const std::string &directory, const std::string &prefix = "",
              DumpFormat matrixFormat = DUMP_BINARY, DumpFormat rhsFormat = DUMP_BINARY,
              DumpFormat solutionFormat = DUMP_BINARY) const;

    virtual void setNullSpace();
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

    bool isForceConsistencyEnabled() const;
    void enableForceConsistency(bool enable);

protected:
    bool m_transpose;

    Mat m_A;
    Vec m_rhs;
    Vec m_solution;

    KSP m_KSP;
    KSPOptions m_KSPOptions;
    KSPStatus m_KSPStatus;

    void matrixCreate(const SystemMatrixAssembler &assembler);
    void matrixRead(const std::string &matrixPath);
    void matrixFill(const SystemMatrixAssembler &assembler);
    void matrixUpdate(long nRows, const long *rows, const SystemMatrixAssembler &assembler);

    void vectorsCreate();
    void vectorsPermute(bool invert);
    void vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution);
    void vectorsExport(std::vector<double> *solution);

    virtual void preKSPSetupActions();
    virtual void postKSPSetupActions();

    virtual void preKSPSolveActions();
    virtual void postKSPSolveActions();

#if BITPIT_ENABLE_MPI==1
    const MPI_Comm & getCommunicator() const;
#endif

private:
    static PetscManager m_petscManager;

    static int m_nInstances;

    std::string m_prefix;

    bool m_assembled;
    bool m_setUp;

#if BITPIT_ENABLE_MPI==1
    MPI_Comm m_communicator;

    bool m_partitioned;

    long m_rowGlobalOffset;
    long m_colGlobalOffset;
#endif

    IS m_rowPermutation;
    IS m_colPermutation;

    bool m_forceConsistency;

#if BITPIT_ENABLE_MPI==1
    void setCommunicator(MPI_Comm communicator);
    void freeCommunicator();
#endif

    void resetPermutations();

    void removeNullSpaceFromRHS();

};

}

#endif
