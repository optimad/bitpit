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

#include "stencil_solver.hpp"

#include <cmath>

namespace bitpit {

/*!
 * \class DiscretizationStencilStorageInterface
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilStorageInterface class defines the interface
 * for stencil storage.
 */

/*!
 * \class DiscretizationStencilProxyBaseStorage
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilProxyBaseStorage class defines a proxy for
 * stencil storage.
 */

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 * \param stride is the stride
 */
template<typename stencil_t, typename stencil_container_t>
DiscretizationStencilProxyBaseStorage<stencil_t, stencil_container_t>::DiscretizationStencilProxyBaseStorage(const stencil_container_t *stencils, int stride)
    : DiscretizationStencilStorageInterface<stencil_t>(),
      m_stencils(stencils), m_stride(stride)
{
}

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t, typename stencil_container_t>
DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::DiscretizationStencilProxyStorage(const stencil_container_t *stencils)
    : DiscretizationStencilProxyBaseStorage<stencil_t, stencil_container_t>(stencils, 1)
{
}

/*!
 * Get the size of the container.
 *
 * \result The size of the container.
 */
template<typename stencil_t, typename stencil_container_t>
std::size_t DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::size() const
{
    return this->m_stencils->size();
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the storage
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::at(long rowIndex) const
{
    return (*(this->m_stencils))[rowIndex];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowRawIndex is the raw index of the row in the storage
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::rawAt(std::size_t rowRawIndex) const
{
    return (*(this->m_stencils))[rowRawIndex];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockIndex is the index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::at(long blockIndex, int componentIdx) const
{
    return (*(this->m_stencils))[blockIndex * this->m_stride + componentIdx];
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockRawIndex is the raw index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename stencil_container_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>::rawAt(std::size_t blockRawIndex, int componentIdx) const
{
    return (*(this->m_stencils))[blockRawIndex * this->m_stride + componentIdx];
}

/*!
 * Constructor.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t>
DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::DiscretizationStencilProxyStorage(const PiercedStorage<stencil_t> *stencils)
    : DiscretizationStencilProxyBaseStorage<stencil_t, PiercedStorage<stencil_t>>(stencils, stencils->getFieldCount())
{
}

/*!
 * Get the size of the container.
 *
 * \result The size of the container.
 */
template<typename stencil_t>
std::size_t DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::size() const
{
    return this->m_stencils->getKernel()->size() * this->m_stride;
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the storage
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::at(long rowIndex) const
{
    return this->m_stencils->at(rowIndex);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowRawIndex is the raw index of the row in the storage
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::rawAt(std::size_t rowRawIndex) const
{
    return this->m_stencils->rawAt(rowRawIndex);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockIndex is the index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::at(long blockIndex, int componentIdx) const
{
    return this->m_stencils->at(blockIndex, componentIdx);
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param blockRawIndex is the raw index of the block
 * \param componentIdx is the index of the component inside the block
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t>
const stencil_t & DiscretizationStencilProxyStorage<stencil_t, PiercedStorage<stencil_t>>::rawAt(std::size_t blockRawIndex, int componentIdx) const
{
    return this->m_stencils->rawAt(blockRawIndex, componentIdx);
}

/*!
 * \class DiscretizationStencilSolverAssembler
 * \ingroup discretization
 *
 * \brief The DiscretizationStencilSolverAssembler class defines an assembler
 * for building the stencil solver.
 */

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t, typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(const stencil_container_t *stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(MPI_COMM_SELF, false, stencils, std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}

/*!
 * Constructor.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 * \param stencils are the stencils
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t, typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned,
                                                                                                       const stencil_container_t *stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(communicator, partitioned,
                                           std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>>(new DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>(stencils)),
                                           std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}
#else
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t, typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(const stencil_container_t *stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>>(new DiscretizationStencilProxyStorage<stencil_t, stencil_container_t>(stencils)),
                                           std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}
#endif

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param stencils are the stencils
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(MPI_COMM_SELF, false, std::move(stencils), std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}

/*!
 * Constructor.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 * \param stencils are the stencils
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned,
                                                                                                       std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(communicator, partitioned, std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
#else
/*!
 * Constructor.
 *
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 * \param stencils are the stencils
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
#endif
{
    setStencils(std::move(stencils));
    setMatrixSizes();
    setBlockSize();
    setMaximumRowNZ();
}

#if BITPIT_ENABLE_MPI==1
/*!
 * Constructor.
 *
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(AssemblerKernelArgs&&... assemblerKernelArgs)
    : DiscretizationStencilSolverAssembler(MPI_COMM_SELF, false, std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}

/*!
 * Constructor.
 *
 * \param communicator is the MPI communicator
 * \param partitioned controls if the matrix is partitioned
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(MPI_Comm communicator, bool partitioned,
                                                                                                       AssemblerKernelArgs&&... assemblerKernelArgs)
    : solver_kernel_type::Assembler(std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...),
      m_partitioned(partitioned), m_communicator(communicator)
{
}
#else
/*!
 * Constructor.
 *
 * \param assemblerKernelArgs are the arguments that will be passed to the constructor of the
 * assembler of the solver kernel
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename... AssemblerKernelArgs>
DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::DiscretizationStencilSolverAssembler(AssemblerKernelArgs&&... assemblerKernelArgs)
    : solver_kernel_type::Assembler(std::forward<AssemblerKernelArgs>(assemblerKernelArgs)...)
{
}
#endif

#if BITPIT_ENABLE_MPI==1
/*!
 * Checks if the matrix is partitioned.
 *
 * \result Returns true if the patch is partitioned, false otherwise.
 */
template<typename stencil_t, typename solver_kernel_t>
bool DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::isPartitioned() const
{
    return m_partitioned;
}

/*!
 * Gets the MPI communicator associated to the matrix.
 *
 * \return The MPI communicator associated to the matrix.
 */
template<typename stencil_t, typename solver_kernel_t>
const MPI_Comm & DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getCommunicator() const
{
    return m_communicator;
}
#endif

/*!
 * Get the assembly options.
 *
 * \result The assembly options that will be used.
 */
template<typename stencil_t, typename solver_kernel_t>
typename DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::assembly_options_type DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getOptions() const
{
    assembly_options_type options;
    options.full   = true;
    options.sorted = false;

    return options;
}

/*!
 * Set block size.
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename W, typename V, typename std::enable_if<std::is_fundamental<W>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setBlockSize()
{
    setBlockSize(1);
}

/*!
 * Set block size.
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setBlockSize()
{
    setBlockSize(sizeof(typename StencilVector::weight_type) / sizeof(typename StencilVector::weight_type::value_type));
}

/*!
 * Initialize block size.
 *
 * The block size is set equal to the square root of the weight/constant size; if the
 * square root of the weight type is not an integer number, an exception is throw.
 *
 * Block size is evaluated from the constant of the first stencil. If the size of the
 * other weights don't match the size of the evaluated size, an exception is thrown
 * only when bitpit is compiled in debug mode, otherwise the error is silently
 * ignored.
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setBlockSize()
{
    // Get block size
    //
    // Block size is evaluated from the constant of the first stencil.
    //
    // The block size is set equal to the square root of the weight size; if the square
    // root of the weight type is not an integer number, an exception is throw.
    if (getRowCount() == 0) {
        throw std::runtime_error("Unable to evaluate the block size.");
    }

    const StencilBlock &stencil = getRowStencil(0);
    std::size_t stencilConstantSize = stencil.getConstant().size();
    int blockSize = static_cast<int>(std::round(std::sqrt(stencilConstantSize)));
    if (static_cast<std::size_t>(blockSize * blockSize) != stencilConstantSize) {
        throw std::runtime_error("Weights size should be a square.");
    }
    setBlockSize(blockSize);

#ifdef DEBUG
    // Validate block size
    //
    // All weight sizes should match
    for (long i = 0; i < getRowCount(); ++i) {
        const StencilBlock &stencil = getRowStencil(i);
        const StencilBlock::weight_type *weightData = stencil.weightData();
        std::size_t stencilSize = stencil.size();

        for (std::size_t k = 0; k < stencilSize; ++k) {
            if (weightData[k].size() != m_blockSize * m_blockSize)) {
                throw std::runtime_error("All stencils weights should have the same size.");
            }
        }

        if (stencil.getConstant().size() != m_blockSize * m_blockSize) {
            throw std::runtime_error("The stencils constant should have the same size of the stencil weights.");
        }
    }
#endif
}

/*!
 * Set block size.
 *
 * \param blockSize is the block size
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setBlockSize(int blockSize)
{
    m_blockSize = blockSize;
}

/*!
 * Set the stencils.
 *
 * \param stencils are the stencils
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setStencils(std::unique_ptr<DiscretizationStencilStorageInterface<stencil_t>> &&stencils)
{
    m_stencils = std::move(stencils);
}

/*!
 * Set matrix sizes.
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setMatrixSizes()
{
    setMatrixSizes(m_stencils->size(), m_stencils->size());
}

/*!
 * Set matrix sizes.
 *
 * \param nRows are the rows of the matrix
 * \param nCols are the columns of the matrix
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setMatrixSizes(long nRows, long nCols)
{
    // Set system sizes
    m_nRows = nRows;
    m_nCols = nCols;

#if BITPIT_ENABLE_MPI==1
    // Global system sizes
    m_nGlobalRows = nRows;
    m_nGlobalCols = nCols;
    if (m_partitioned) {
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalRows, 1, MPI_LONG, MPI_SUM, m_communicator);
        MPI_Allreduce(MPI_IN_PLACE, &m_nGlobalCols, 1, MPI_LONG, MPI_SUM, m_communicator);
    }

    // Global offsets
    m_globalRowOffset = 0;
    m_globalColOffset = 0;
    if (m_partitioned) {
        int nProcessors;
        MPI_Comm_size(m_communicator, &nProcessors);

        std::vector<long> nRankRows(nProcessors);
        MPI_Allgather(&m_nRows, 1, MPI_LONG, nRankRows.data(), 1, MPI_LONG, m_communicator);

        std::vector<long> nRankCols(nProcessors);
        MPI_Allgather(&m_nCols, 1, MPI_LONG, nRankCols.data(), 1, MPI_LONG, m_communicator);

        int rank;
        MPI_Comm_rank(m_communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_globalRowOffset += nRankRows[i];
            m_globalColOffset += nRankCols[i];
        }
    }
#endif
}

/*!
 * Set the maximum number of non-zero element on a single row.
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setMaximumRowNZ()
{
    long maxRowNZ = 0;
    for (long n = 0; n < getRowCount(); ++n) {
        maxRowNZ = std::max(getRowNZCount(n), maxRowNZ);
    }

    setMaximumRowNZ(maxRowNZ);
}

/*!
 * Set the maximum number of non-zero element on a single row.
 *
 * \param maxRowNZ is the maximum number of non-zero element on a single row
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::setMaximumRowNZ(long maxRowNZ)
{
    m_maxRowNZ = maxRowNZ;
}

/*!
 * Get the stencil block size.
 *
 * \result The stencil block size.
 */
template<typename stencil_t, typename solver_kernel_t>
int DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getBlockSize() const
{
    return m_blockSize;
}

/*!
 * Get the number of (block) rows handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block rows, where a block row is
 * defined as a group of blockSize matrix rows.
 *
 * \result The number of (block) rows handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowCount() const
{
    return m_nRows;
}

/*!
 * Get the number of (block) columns handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the number of block columns, where a block column
 * is defined as a group of blockSize matrix columns.
 *
 * \result The number of (block) columns handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColCount() const
{
    return m_nCols;
}

/*!
 * Get the number of elements in the rows handled by the assembler.
 *
 * This function will return the effective number of rows of the matrix that
 * will be assembled.
 *
 * \result The number of rows handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowElementCount() const
{
    long nRowElements = getBlockSize() * getRowCount();

    return nRowElements;
}

/*!
 * Get the number of elements in the columns handled by the assembler.
 *
 * This function will return the effective number of columns of the matrix that
 * will be assembled.
 *
 * \result The number of columns handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColElementCount() const
{
    long nColElements = getBlockSize() * getColCount();

    return nColElements;
}

#if BITPIT_ENABLE_MPI==1
/**
 * Get number of global (block) rows handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block rows, where a block row
 * is defined as a group of blockSize matrix rows.
 *
 * \result The number of global rows handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowGlobalCount() const
{
    return m_nGlobalRows;
}

/**
 * Get number of global (block) columns handled by the assembler.
 *
 * If the matrix is a block matrix (i.e., the block size is greater than one),
 * this function will return the global number of block columns, where a block
 * column is defined as a group of blockSize matrix columns.
 *
 * \result The number of global (block) columns handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColGlobalCount() const
{
    return m_nGlobalCols;
}

/**
 * Get the number of global elements in the rows handled by the assembler.
 *
 * This function will return the effective global number of rows of the system
 * matrix.
 *
 * \result The number of global elements in the rows handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowGlobalElementCount() const
{
    long nElements = getBlockSize() * getRowGlobalCount();

    return nElements;
}

/*!
 * Get the global number of columns handled by the assembler.
 *
 * This function will return the effective global number of columns of the
 * system matrix.
 *
 * \result The global number of columns handled by the assembler.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColGlobalElementCount() const
{
    long nElements = getBlockSize() * getColGlobalCount();

    return nElements;
}

/*!
 * Get global (block) row offset.
 *
 * \result The global (block) row offset.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowGlobalOffset() const
{
    return m_globalRowOffset;
}

/*!
 * Get global (block) column offset.
 *
 * \result The global (block) column offset.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColGlobalOffset() const
{
    return m_globalColOffset;
}

/*!
 * Get global offset for the elements of the row.
 *
 * This function will return the offset expressed in effective rows of the
 * system matrix.
 *
 * \result The global offset for the elements of the row.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowGlobalElementOffset() const
{
    long offset = getBlockSize() * getRowGlobalOffset();

    return offset;
}

/*!
 * Get global offset for the elements of the column.
 *
 * This function will return the offset expressed in effective columns of the
 * system matrix.
 *
 * \result The global offset for the elements of the column.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getColGlobalElementOffset() const
{
    long offset = getBlockSize() * getColGlobalOffset();

    return offset;
}
#endif

/*!
 * Get the number of non-zero elements in the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The number of non-zero elements in the specified row.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowNZCount(long rowIndex) const
{
    const stencil_t &stencil = getRowStencil(rowIndex);
    std::size_t stencilSize = stencil.size();

    return stencilSize;
}

/**
 * Get the maximum number of non-zero elements per row.
 *
 * \result The maximum number of non-zero elements per row.
 */
template<typename stencil_t, typename solver_kernel_t>
long DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getMaxRowNZCount() const
{
    return m_maxRowNZ;
}

/*!
 * Get the pattern of the specified row.
 *
 * If the assembler is a block assembler (i.e., the block size is greater than
 * one), this function will return the global ids of the block columns of the
 * row, where a block column is defined as a group of blockSize assembler
 * columns.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the pattern of the specified row
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowPattern(long rowIndex, ConstProxyVector<long> *pattern) const
{
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get pattern
    getPattern(stencil, pattern);
}

/*!
 * Get the pattern of the specified stencil.
 *
 * \param stencil is the stencil
 * \param pattern on output will contain the pattern of the specified stencil
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getPattern(const stencil_t &stencil, ConstProxyVector<long> *pattern) const
{
    std::size_t stencilSize = stencil.size();

    const long *patternData = stencil.patternData();
    pattern->set(patternData, stencilSize);
}

/*!
 * Get the values of the specified (block) row.
 *
 * If the assembler is a block assembler (i.e., the block size is greater than
 * one), this function will return the values of all the elements of a block row,
 * where a block column is defined as a group of blockSize assembler columns. The
 * values are returned as a row-oriented logically two-dimensional array of
 * values.
 *
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowValues(long rowIndex, ConstProxyVector<double> *values) const
{
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get values
    getValues(stencil, values);
}

/*!
 * Get the values of the specified stencil.
 *
 * \param stencil is the stencil
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const
{
    values->set(stencil.weightData(), m_blockSize * stencil.size());
}

/*!
 * Get the values of the specified stencil.
 *
 * \param stencil is the stencil
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getValues(const stencil_t &stencil, ConstProxyVector<double> *values) const
{
    std::size_t stencilSize = stencil.size();
    const stencil_weight_type *stencilWeightData = stencil.weightData();

    int nBlockElements = m_blockSize * m_blockSize;
    std::size_t nRowValues = m_blockSize * stencilSize;
    std::size_t nValues = nBlockElements * stencilSize;
    values->set(ConstProxyVector<double>::INTERNAL_STORAGE, nValues);
    ConstProxyVector<double>::storage_pointer expandedValuesStorage = values->storedData();


    for (std::size_t k = 0; k < stencilSize; ++k) {
        const double *weightData = stencilWeightData[k].data();
        for (int i = 0; i < m_blockSize; ++i) {
            int weightOffset = linearalgebra::linearIndexRowMajor(i, 0, m_blockSize, m_blockSize);
            int valuesOffset = linearalgebra::linearIndexRowMajor(i, m_blockSize * k, m_blockSize, nRowValues);

            std::copy_n(weightData + weightOffset, m_blockSize, expandedValuesStorage + valuesOffset);
        }
    }
}

/*!
 * Get the data of the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param pattern on output will contain the values of the specified row
 * \param values on output will contain the values of the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * two-dimensional array that uses a col-major order
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowData(long rowIndex, ConstProxyVector<long> *pattern, ConstProxyVector<double> *values) const
{
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get pattern
    getPattern(stencil, pattern);

    // Get values
    getValues(stencil, values);
}

/*!
 * Get the constant associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \param constant is the constant associated with the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * one-dimensional array
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowConstant(long rowIndex, bitpit::ConstProxyVector<double> *constant) const
{
    // Get stencil information
    const stencil_t &stencil = getRowStencil(rowIndex);

    // Get constant
    getConstant(stencil, constant);
}

/*!
 * Get the constant associated with the specified stencil.
 *
 * \param stencil is the stencil
 * \param constant is the constant associated with the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * one-dimensional array
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getConstant(const stencil_t &stencil, bitpit::ConstProxyVector<double> *constant) const
{
    const stencil_weight_type &stencilConstant = stencil.getConstant();

    constant->set(ConstProxyVector<double>::INTERNAL_STORAGE, m_blockSize);
    ConstProxyVector<double>::storage_pointer constantStorage = constant->storedData();
    std::copy_n(&stencilConstant, m_blockSize, constantStorage);
}

/*!
 * Get the constant associated with the specified stencil.
 *
 * \param stencil is the stencil
 * \param constant is the constant associated with the specified (block) row.
 * If the block size is greater than one, values will be stored in a logically
 * one-dimensional array
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
void DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getConstant(const stencil_t &stencil, bitpit::ConstProxyVector<double> *constant) const
{
    const stencil_weight_type &stencilConstant = stencil.getConstant();

    constant->set(ConstProxyVector<double>::INTERNAL_STORAGE, m_blockSize);
    ConstProxyVector<double>::storage_pointer constantStorage = constant->storedData();
    for (int i = 0; i < m_blockSize; ++i) {
        int offset_i = linearalgebra::linearIndexRowMajor(i, 0, m_blockSize, m_blockSize);

        constantStorage[i] = 0;
        for (int j = 0; j < m_blockSize; ++j) {
            constantStorage[i] += stencil.getWeightManager().at(stencilConstant, offset_i + j);
        }
    }
}

/*!
 * Get the stencil associated with the specified row.
 *
 * \param rowIndex is the index of the row in the assembler
 * \result The stencil associated with the specified row.
 */
template<typename stencil_t, typename solver_kernel_t>
const stencil_t & DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t>::getRowStencil(long rowIndex) const
{
    return m_stencils->at(rowIndex);
}

/*!
* \ingroup discretization
* \class DiscretizationStencilSolver
*
* The DiscretizationStencilSolver class handles the solution of linear systems
* assembled from discretization stencils.
*/

/*!
* Clear the stencil solver
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::clear()
{
    solver_kernel_t::clear();

    std::vector<double>().swap(m_constants);
}

#if BITPIT_ENABLE_MPI==1
/*!
* Assembly the stencil solver.
*
* After assembying th system solver, its options will be reset.
*
* \param stencils are the stencils
*/
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::assembly(const stencil_container_t &stencils)
{
    assembly(MPI_COMM_SELF, false, stencils);
}

/*!
* Initialize the stencil solver.
*
* \param partitioned controls if the matrix is partitioned
* \param communicator is the MPI communicator
* \param stencils are the stencils
*/
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::assembly(MPI_Comm communicator, bool partitioned, const stencil_container_t &stencils)
#else
/*!
* Initialize the stencil solver.
*
* \param stencils are the stencils
*/
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::assembly(const stencil_container_t &stencils)
#endif
{
    // Create the assembler
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t> assembler(communicator, partitioned, &stencils);
#else
    DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t> assembler(&stencils);
#endif

    // Assembly the system
    solver_kernel_t::template assembly<DiscretizationStencilSolver<stencil_t, solver_kernel_t>>(assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the system.
 *
 * After assembying th system solver, its options will be reset.
 *
 * \param assembler is the matrix assembler
 * \param reordering is the reordering that will be applied when assembling the system
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::assembly(const Assembler &assembler)
{
    solver_kernel_t::template assembly<DiscretizationStencilSolver<stencil_t, solver_kernel_t>>(assembler, NaturalSystemMatrixOrdering());
}

/*!
 * Assembly the matrix.
 *
 * \param assembler is the matrix assembler
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::matrixAssembly(const Assembler &assembler)
{
    // Assemble matrix
    solver_kernel_t::matrixAssembly(assembler);

    // Assemble constants
    assembleConstants(assembler);
}

/*!
 * Update the specified rows of the matrix.
 *
 * The contents of the specified rows will be replaced by the data provided by the given
 * assembler. If the matrix has not been assembled yet, both the pattern and the values
 * of the matrix will be updated. After the matrix has been assembled only the values
 * will be updated.
 *
 * The block size of the assembler should be equal to the block size of the matrix.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the local indices of the rows that will be updated, if a
 * null pointer is passed, the rows that will be updated are the rows
 * from 0 to (nRows - 1).
 * \param assembler is the matrix assembler for the rows that will be updated
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::matrixUpdate(long nRows, const long *rows, const Assembler &assembler)
{
    // Update the system
    solver_kernel_t::matrixUpdate(nRows, rows, assembler);

    // Update the constants
    updateConstants(nRows, rows, assembler);
}

/*!
 * Update all the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param stencils are the stencils that will be used to update the rows
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::update(const stencil_container_t &stencils)
{
    update(this->getRowCount(), nullptr, stencils);
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
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::update(const std::vector<long> &rows, const stencil_container_t &stencils)
{
    update(rows.size(), rows.data(), stencils);
}

/*!
 * Update the stencil solver.
 *
 * Only the values of the system matrix and the values of the constants can be
 * updated, once the system is initialized its pattern cannot be modified.
 *
 * \param nRows is the number of stencils that will be updated
 * \param rows are the rows of the stencils that will be updated,
 * if a null pointer is passed, the stencils that will be updated are the
 * stencils from 0 to (nRows - 1).
 * \param stencils are the stencils that will be used to update the rows
 */
template<typename stencil_t, typename solver_kernel_t>
template<typename stencil_container_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::update(std::size_t nRows, const long *rows, const stencil_container_t &stencils)
{
#if BITPIT_ENABLE_MPI==1
    DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t> assembler(this->getCommunicator(), this->isPartitioned(), &stencils);
#else
    DiscretizationStencilSolverAssembler<stencil_t, solver_kernel_t> assembler(&stencils);
#endif

    // Update the system
    solver_kernel_t::template update<DiscretizationStencilSolver<stencil_t, solver_kernel_t>>(nRows, rows, assembler);
}


/*!
 * Update the system.
 *
 * Only the values of the system matrix can be updated, once the system is
 * assembled its pattern cannot be modified.
 *
 * \param nRows is the number of rows that will be updated
 * \param rows are the indices of the rows that will be updated
 * \param assembler is the matrix assembler for the rows that will be updated
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::update(long nRows, const long *rows, const Assembler &assembler)
{
    solver_kernel_t::template update<DiscretizationStencilSolver<stencil_t, solver_kernel_t>>(nRows, rows, assembler);
}

/*!
 * Assemble the constants associated with stencil solver.
 *
 * \param assembler is the solver assembler
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::assembleConstants(const Assembler &assembler)
{
    long nRows = assembler.getRowCount();
    int blockSize = assembler.getBlockSize();

    // Create the storage
    m_constants.resize(nRows * blockSize);

    // Update the constants
    updateConstants(nRows, nullptr, assembler);
}

/*!
 * Update the constants associated with stencil solver.
 *
 * \param nRows is the number of stencils that will be updated
 * \param rows are the rows of the stencils that will be updated,
 * if a null pointer is passed, the stencils that will be updated are the
 * stencils from 0 to (nRows - 1).
 * \param assembler is the solver assembler
 */
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::updateConstants(std::size_t nRows, const long *rows,
                                                              const Assembler &assembler)
{
    int blockSize = assembler.getBlockSize();
    ConstProxyVector<double> rowConstant(blockSize);
    for (std::size_t n = 0; n < nRows; ++n) {
        long row;
        if (rows) {
            row = rows[n];
        } else {
            row = n;
        }

        assembler.getRowConstant(n, &rowConstant);
        std::copy_n(rowConstant.data(), blockSize, m_constants.data() + row * blockSize);
    }
}

/*!
* Solve the system.
*/
template<typename stencil_t, typename solver_kernel_t>
void DiscretizationStencilSolver<stencil_t, solver_kernel_t>::solve()
{
    // Check if the stencil solver is assembled
    if (!this->isAssembled()) {
        throw std::runtime_error("Unable to solve the system. The stencil solver is not yet assembled.");
    }

    // Subtract constant terms to the RHS
    long nUnknowns = this->getBlockSize() * this->getRowCount();
    double *raw_rhs = this->getRHSRawPtr();
    for (long i = 0; i < nUnknowns; ++i) {
        raw_rhs[i] -= m_constants[i];
    }
    this->restoreRHSRawPtr(raw_rhs);

    // Solve the system
    solver_kernel_t::solve();
}

}
