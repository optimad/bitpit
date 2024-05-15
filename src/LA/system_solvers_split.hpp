/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2024 OPTIMAD engineering Srl
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

#ifndef __BITPIT_SYSTEM_SOLVERS_SPLIT_HPP__
#define __BITPIT_SYSTEM_SOLVERS_SPLIT_HPP__

#include <string>
#include <vector>

#include "system_matrix.hpp"
#include "system_solvers_large.hpp"

namespace bitpit {

enum SystemSplitStrategy {
    SYSTEM_SPLIT_STRATEGY_DIAGONAL,
    SYSTEM_SPLIT_STRATEGY_FULL,
    SYSTEM_SPLIT_STRATEGY_LOWER,
};

class SplitSystemMatrixAssembler : virtual public SystemMatrixAssembler {

public:
    SystemSplitStrategy getSplitStrategy() const;
    int getSplitCount() const;
    const std::vector<int> & getSplitSizes() const;

protected:
    SplitSystemMatrixAssembler(SystemSplitStrategy splitStrategy, const std::vector<int> &splitSizes);

private:
    SystemSplitStrategy m_splitStrategy;
    std::vector<int> m_splitSizes;

};

class SplitSystemSparseMatrixAssembler : public SystemSparseMatrixAssembler, public SplitSystemMatrixAssembler  {

public:
    SplitSystemSparseMatrixAssembler(const SparseMatrix *matrix, SystemSplitStrategy splitStrategy,
                                     const std::vector<int> &splitSizes);

};

class SplitSystemSolver : public SystemSolver {

public:
    typedef SplitSystemMatrixAssembler Assembler;

    typedef SystemSplitStrategy splitStrategy;

    friend class SystemSolver;

    SplitSystemSolver(bool debug = false);
    SplitSystemSolver(bool transpose, bool debug);
    SplitSystemSolver(bool flatten, bool transpose, bool debug);
    SplitSystemSolver(const std::string &prefix, bool debug = false);
    SplitSystemSolver(const std::string &prefix, bool transpose, bool debug);
    SplitSystemSolver(const std::string &prefix, bool flatten, bool transpose, bool debug);

    using SystemSolver::SystemSolver;

    void assembly(const SparseMatrix &matrix, SystemSplitStrategy splitStrategy, const std::vector<int> &splitSizes);
    void assembly(const SparseMatrix &matrix, SystemSplitStrategy splitStrategy, const std::vector<int> &splitSizes,
                  const SystemMatrixOrdering &reordering);
    void assembly(const Assembler &assembler);
    void assembly(const Assembler &assembler, const SystemMatrixOrdering &reordering);

    void update(const SparseMatrix &elements);
    void update(long nRows, const long *rows, const SparseMatrix &elements);
    void update(const Assembler &assembler);
    void update(long nRows, const long *rows, const Assembler &assembler);

    int getBlockSize() const override;

    SystemSplitStrategy getSplitStrategy() const;
    int getSplitCount() const;
    std::vector<int> getSplitSizes() const;
    std::vector<int> getSplitOffsets() const;

    KSPOptions & getSplitKSPOptions(int split);
    const KSPOptions & getSplitKSPOptions(int split) const;

    const KSPStatus & getSplitKSPStatus(int split) const;

    void exportMatrix(const std::string &filePath, FileFormat exportFormat = FILE_BINARY) const override;

protected:
    SystemSplitStrategy m_splitStrategy;

    std::vector<Mat> m_splitAs;

    std::vector<KSPOptions> m_splitKSPOptions;
    std::vector<KSPStatus> m_splitKSPStatuses;

    using SystemSolver::assembly;
    using SystemSolver::update;

    void matrixAssembly(const Assembler &assembler);
    void matrixUpdate(long nRows, const long *rows, const Assembler &assembler);
    void matrixFill(const std::string &filePath) override;
    void matrixDump(std::ostream &stream, const std::string &directory, const std::string &prefix) const override;
    void matrixRestore(std::istream &stream, const std::string &directory, const std::string &prefix) override;
    void matrixDestroy() override;

    void vectorsCreate() override;
    void vectorsReorder(bool invert) override;
    void vectorsRestore(std::istream &stream, const std::string &directory, const std::string &prefix) override;
    void vectorsCreateSplitPermutations();
    void vectorsDestroy() override;

    void setupPreconditioner() override;
    using SystemSolver::setupPreconditioner;
    virtual void setupSplitPreconditioners();

    void setupKrylov() override;
    using SystemSolver::setupKrylov;
    virtual void setupSplitKrylovs();

    void postKSPSolveActions() override;

    void initializeKSPOptions() override;
    virtual void initializeSplitKSPOptions();
    using SystemSolver::resetKSPOptions;
    void destroyKSPOptions() override;
    virtual void destroySplitKSPOptions();

    void initializeKSPStatus() override;
    virtual void initializeSplitKSPStatuses();
    void fillKSPStatus() override;
    using SystemSolver::fillKSPStatus;
    virtual void fillSplitKSPStatuses();
    void resetKSPStatus() override;
    using SystemSolver::resetKSPStatus;
    virtual void resetSplitKSPStatuses();
    void destroyKSPStatus() override;
    virtual void destroySplitKSPStatuses();

    void dumpInfo(std::ostream &stream) const override;
    void restoreInfo(std::istream &stream) override;

    using SystemSolver::exportMatrix;

#if BITPIT_ENABLE_MPI == 1
    void generateSplitPermutation(long nItems, long itemGlobalOffset, IS *splitReordering) const;
#else
    void generateSplitPermutation(long nItems, IS *splitReordering) const;
#endif
    void generateSplitIndexes(int split, long nItems, std::vector<std::size_t> *indexes) const;

    std::string generateSplitPath(const std::string &path, int i) const;
    std::string generateSplitPath(const std::string &path, int i, int j) const;
    std::string generateSplitPath(const std::string &path, const std::string &index) const;

private:
    IS m_rhsSplitPermutation;
    IS m_solutionSplitPermutation;

    int getBlockSplitLinearIndex(int i, int j) const;
    int getBlockSplitLinearIndex(int i, int j, int nSplits) const;

};

}

#endif
