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

#include "bitpit_IO.hpp"

#include "system_matrix.hpp"

namespace bitpit {

struct KSPOptions {
    PetscInt overlap; //! Overlap between a pair of subdomains for the partitioned preconditioner
    PetscInt levels; //! Number of levels of fill used by the preconditioner.

    PetscBool initial_non_zero; //! Tells the iterative solver that the initial guess is nonzero
    PetscInt restart; //! Number of iterations at which the GMRES method restarts
    PetscInt maxits; //! Maximum number of iterations
    PetscScalar rtol; //! Relative convergence tolerance, relative decrease in the preconditioned residual norm
    PetscScalar atol; //! Absolute convergence tolerance, absolute size of the preconditioned residual norm

    PetscInt sublevels; //! Deprecated, ASM ILU levels should be set using the "levels" member
    PetscScalar subrtol; //! Deprecated, it has never had any effect on the solution of the system

    KSPOptions()
        : overlap(PETSC_DEFAULT), levels(PETSC_DEFAULT),
          initial_non_zero(PETSC_TRUE), restart(PETSC_DEFAULT),
          maxits(PETSC_DEFAULT), rtol(PETSC_DEFAULT), atol(PETSC_DEFAULT),
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

class SystemMatrixOrdering
{

public:
    virtual ~SystemMatrixOrdering() = default;

    /*!
     * Get the rank of the specified local row.
     *
     * The permutation rank defines the position, in the final assembled matrix, of the
     * specified local row (i.e., the position of the specified local row after applying
     * the reordering).
     *
     * \param row is the local row
     * \result The rank of the specified local row.
     */
    virtual long getRowPermutationRank(long row) const = 0;

    /*!
     * Get the rank of the specified local column.
     *
     * The permutation rank defines the position, in the final assembled matrix, of the
     * specified local column (i.e., the position of the specified local column after
     * applying the reordering).
     *
     * \param col is the local column
     * \result The rank of the specified local column.
     */
    virtual long getColPermutationRank(long col) const = 0;

protected:
    SystemMatrixOrdering() = default;

};

class NaturalSystemMatrixOrdering : public SystemMatrixOrdering
{

public:
    NaturalSystemMatrixOrdering() = default;

    long getRowPermutationRank(long row) const override;
    long getColPermutationRank(long col) const override;

};

template<typename RowRankStorage, typename ColRankStorage>
class ProxySystemMatrixOrdering : public SystemMatrixOrdering
{

public:
    ProxySystemMatrixOrdering(const RowRankStorage *rowRankStorage, const ColRankStorage *colRankStorage);

    long getRowPermutationRank(long row) const override;
    long getColPermutationRank(long col) const override;

private:
    const RowRankStorage *m_rowRankStorage;
    const ColRankStorage *m_colRankStorage;

};

class SystemMatrixAssembler {

public:
    struct AssemblyOptions {
        bool full; //!< Controls if the assembler is providing all the non-zero values of a row
        bool sorted; //! Controls if the values provided by the assembler are sorted by ascending column
    };

    virtual ~SystemMatrixAssembler() = default;

#if BITPIT_ENABLE_MPI==1
    virtual bool isPartitioned() const = 0;
    virtual const MPI_Comm & getCommunicator() const = 0;
#endif

    virtual AssemblyOptions getOptions() const = 0;

    virtual int getBlockSize() const = 0;

    virtual long getRowCount() const = 0;
    virtual long getColCount() const = 0;

    virtual long getRowElementCount() const = 0;
    virtual long getColElementCount() const = 0;

#if BITPIT_ENABLE_MPI==1
    virtual long getRowGlobalCount() const = 0;
    virtual long getColGlobalCount() const = 0;

    virtual long getRowGlobalElementCount() const = 0;
    virtual long getColGlobalElementCount() const = 0;

    virtual long getRowGlobalOffset() const = 0;
    virtual long getColGlobalOffset() const = 0;

    virtual long getRowGlobalElementOffset() const = 0;
    virtual long getColGlobalElementOffset() const = 0;
#endif

    virtual long getRowNZCount(long rowIndex) const = 0;
    virtual long getMaxRowNZCount() const = 0;

    virtual void getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const = 0;
    virtual void getRowValues(long rowIndex, ConstProxyVector<double> *values) const = 0;
    virtual void getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const = 0;

protected:
    SystemMatrixAssembler() = default;

};

class SystemSparseMatrixAssembler : virtual public SystemMatrixAssembler {

public:
    SystemSparseMatrixAssembler(const SparseMatrix *matrix);

#if BITPIT_ENABLE_MPI==1
    bool isPartitioned() const override;
    const MPI_Comm & getCommunicator() const override;
#endif

    AssemblyOptions getOptions() const override;

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

protected:
    const SparseMatrix *m_matrix;

};

class PetscManager {

public:
    PetscManager();

    virtual ~PetscManager();

    bool areOptionsEditable() const;

    bool initialize(bool debug);
    bool finalize();

    void addInitOption(const std::string &option);
    void addInitOptions(int argc, char **args);
    void addInitOptions(const std::vector<std::string> &options);
    void clearInitOptions();

private:
    static PetscErrorCode displayLogView();

    bool m_externalMPIInitialization;
    bool m_externalPETScInitialization;

    std::vector<std::string> m_options;

    bool m_logViewEnabled;

    void enableLogView();

};

class SystemSolver {

public:
    typedef SystemMatrixAssembler Assembler;

    enum FileFormat {
        FILE_BINARY,
        FILE_ASCII
    };

    SystemSolver(bool debug = false);
    SystemSolver(bool transpose, bool debug);
    SystemSolver(bool flatten, bool transpose, bool debug);
    SystemSolver(const std::string &prefix, bool debug = false);
    SystemSolver(const std::string &prefix, bool transpose, bool debug);
    SystemSolver(const std::string &prefix, bool flatten, bool transpose, bool debug);

    virtual ~SystemSolver();

    virtual void clear();
    virtual void clearWorkspace();

    bool isAssembled() const;

    void assembly(const SparseMatrix &matrix);
    void assembly(const SparseMatrix &matrix, const SystemMatrixOrdering &reordering);
    void assembly(const Assembler &assembler);
    void assembly(const Assembler &assembler, const SystemMatrixOrdering &reordering);

    void update(const SparseMatrix &elements);
    void update(long nRows, const long *rows, const SparseMatrix &elements);
    void update(const Assembler &assembler);
    void update(long nRows, const long *rows, const Assembler &assembler);

    BITPIT_DEPRECATED(void setUp());
    BITPIT_DEPRECATED(bool isSetUp() const);

    bool getTranspose() const;
    void setTranspose(bool transpose);

    virtual int getBlockSize() const;

    long getRowCount() const;
    long getColCount() const;
    long getRowElementCount() const;
    long getColElementCount() const;
#if BITPIT_ENABLE_MPI==1
    long getRowGlobalCount() const;
    long getColGlobalCount() const;
    long getRowGlobalElementCount() const;
    long getColGlobalElementCount() const;

    bool isPartitioned() const;
#endif

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    void dumpSystem(const std::string &header, const std::string &directory, const std::string &prefix = "") const;
#if BITPIT_ENABLE_MPI==1
    void restoreSystem(MPI_Comm communicator, bool redistribute, const std::string &directory,
                       const std::string &prefix = "");
#else
    void restoreSystem(const std::string &directory, const std::string &prefix = "");
#endif

    virtual void setNullSpace();
    void unsetNullSpace();

    KSPOptions & getKSPOptions();
    const KSPOptions & getKSPOptions() const;

    const KSPStatus & getKSPStatus() const;

    virtual void exportMatrix(const std::string &filePath, FileFormat exportFormat = FILE_BINARY) const;
    virtual void importMatrix(const std::string &filePath);

    double * getRHSRawPtr();
    const double * getRHSRawPtr() const;
    const double * getRHSRawReadPtr() const;
    void restoreRHSRawPtr(double *raw_rhs);
    void restoreRHSRawReadPtr(const double *raw_rhs) const;
    virtual void exportRHS(const std::string &filePath, FileFormat exportFormat = FILE_BINARY) const;
    virtual void importRHS(const std::string &filePath);

    double * getSolutionRawPtr();
    const double * getSolutionRawPtr() const;
    const double * getSolutionRawReadPtr() const;
    void restoreSolutionRawPtr(double *raw_solution);
    void restoreSolutionRawReadPtr(const double *raw_solution) const;
    virtual void exportSolution(const std::string &filePath, FileFormat exportFormat = FILE_BINARY) const;
    virtual void importSolution(const std::string &filePath);

    bool isForceConsistencyEnabled() const;
    void enableForceConsistency(bool enable);

protected:
    enum VectorSide {
        VECTOR_SIDE_RIGHT, // Vector that the matrix can be multiplied against
        VECTOR_SIDE_LEFT, // Vector that the matrix vector product can be stored in
    };

    bool m_flatten;
    bool m_transpose;

    Mat m_A;
    Vec m_rhs;
    Vec m_solution;

    IS m_rowReordering;
    IS m_colReordering;

    bool m_convergenceMonitorEnabled;

    KSP m_KSP;
    bool m_KSPDirty;
    KSPOptions m_KSPOptions;
    KSPStatus m_KSPStatus;

    virtual int getDumpVersion() const;

    void matrixAssembly(const Assembler &assembler);
    void matrixUpdate(long nRows, const long *rows, const Assembler &assembler);
    virtual void matrixFill(const std::string &filePath);
    virtual void matrixDump(std::ostream &systemStream, const std::string &directory, const std::string &prefix) const;
#if BITPIT_ENABLE_MPI==1
    virtual void matrixRestore(std::istream &systemStream, const std::string &directory, const std::string &prefix, bool redistribute);
#else
    virtual void matrixRestore(std::istream &systemStream, const std::string &directory, const std::string &prefix);
#endif
    virtual void matrixDestroy();

    virtual void vectorsCreate();
    virtual void vectorsFill(const std::vector<double> &rhs, const std::vector<double> &solution);
    virtual void vectorsFill(const std::string &rhsFilePath, const std::string &solutionFilePath);
    virtual void vectorsReorder(bool invert);
    virtual void vectorsDump(std::ostream &systemStream, const std::string &directory, const std::string &prefix) const;
    virtual void vectorsRestore(std::istream &systemStream, const std::string &directory, const std::string &prefix);
    virtual void vectorsDestroy();

    void clearReordering();
    void setReordering(long nRows, long nCols, const SystemMatrixOrdering &reordering);

    template<typename DerivedSystemSolver>
    void assembly(const typename DerivedSystemSolver::Assembler &assembler, const SystemMatrixOrdering &reordering);

    template<typename DerivedSystemSolver>
    void update(long nRows, const long *rows, const typename DerivedSystemSolver::Assembler &assembler);

    virtual void prepareKSP();
    virtual void finalizeKSP();
    virtual void createKSP();
    virtual void destroyKSP();

    virtual void setupPreconditioner();
    virtual void setupPreconditioner(PC pc, const KSPOptions &options) const;
    virtual void prePreconditionerSetupActions();
    virtual void postPreconditionerSetupActions();

    virtual void setupKrylov();
    virtual void setupKrylov(KSP ksp, const KSPOptions &options) const;
    virtual void preKrylovSetupActions();
    virtual void postKrylovSetupActions();

    virtual void solveKSP();
    virtual void preKSPSolveActions();
    virtual void postKSPSolveActions();

    virtual void initializeKSPOptions();
    virtual void resetKSPOptions(KSPOptions *options) const;
    virtual void destroyKSPOptions();

    virtual void initializeKSPStatus();
    virtual void resetKSPStatus();
    virtual void resetKSPStatus(KSPStatus *status) const;
    virtual void fillKSPStatus();
    virtual void fillKSPStatus(KSP ksp, KSPStatus *status) const;
    virtual void destroyKSPStatus();

    virtual void dumpInfo(std::ostream &stream) const;
    virtual void restoreInfo(std::istream &stream);

    void createMatrix(int rowBlockSize, int colBlockSize, Mat *matrix) const;
    void createMatrix(int rowBlockSize, int colBlockSize, int nNestRows, int nNestCols, Mat *subMatrices, Mat *matrix) const;
    void fillMatrix(Mat matrix, const std::string &filePath) const;
    void dumpMatrix(Mat matrix, const std::string &directory, const std::string &name) const;
#if BITPIT_ENABLE_MPI==1
    void restoreMatrix(const std::string &directory, const std::string &name, bool redistribute, Mat *matrix) const;
#else
    void restoreMatrix(const std::string &directory, const std::string &name, Mat *matrix) const;
#endif
    void exportMatrix(Mat matrix, const std::string &filePath, FileFormat fileFormat) const;
    void destroyMatrix(Mat *matrix) const;

    void createVector(int blockSize, Vec *vector) const;
    void createVector(int blockSize, int nestSize, Vec *subVectors, Vec *vector) const;
    void fillVector(Vec vector, const std::string &filePath) const;
    void fillVector(Vec vector, const std::vector<double> &data) const;
    void reorderVector(Vec vector, IS permutations, bool invert) const;
    void dumpVector(Vec vector, const std::string &directory, const std::string &name) const;
    void dumpVector(Vec vector, VectorSide side, const std::string &directory, const std::string &name) const;
    void restoreVector(const std::string &directory, const std::string &name, Mat matrix, VectorSide side, Vec *vector) const;
#if BITPIT_ENABLE_MPI==1
    void restoreVector(const std::string &directory, const std::string &name, bool redistribute, Vec *vector) const;
#else
    void restoreVector(const std::string &directory, const std::string &name, Vec *vector) const;
#endif
    void restoreVector(const std::string &directory, const std::string &name, std::size_t localSize, std::size_t globalSize, Vec *vector) const;
    void exportVector(Vec vector, const std::string &filePath, FileFormat fileFormat) const;
    void exportVector(Vec vector, std::vector<double> *data) const;
    void destroyVector(Vec *vector) const;

    std::string getInfoFilePath(const std::string &directory, const std::string &name) const;
    std::string getDataFilePath(const std::string &directory, const std::string &name) const;
    std::string getFilePath(const std::string &directory, const std::string &name) const;

#if BITPIT_ENABLE_MPI==1
    const MPI_Comm & getCommunicator() const;
#endif

private:
    static PetscManager m_petscManager;

    static int m_nInstances;

    std::string m_prefix;

    bool m_assembled;

#if BITPIT_ENABLE_MPI==1
    MPI_Comm m_communicator;

    bool m_partitioned;
#endif

    bool m_forceConsistency;

#if BITPIT_ENABLE_MPI==1
    void setCommunicator(MPI_Comm communicator);
    void freeCommunicator();
#endif

    void resetPermutations();

    void removeNullSpaceFromRHS();

    void openBinaryArchive(const std::string &directory, const std::string &prefix,
                           int block, IBinaryArchive *archive) const;
    void openBinaryArchive(const std::string &header, const std::string &directory,
                           const std::string &prefix, int block, OBinaryArchive *archive) const;

    void closeBinaryArchive(BinaryArchive *archive) const;

    int getBinaryArchiveBlock() const;
};

}

// Include template implementations
#include "system_solvers_large.tpp"

#endif
