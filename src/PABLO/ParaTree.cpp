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

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "bitpit_common.hpp"
#if BITPIT_ENABLE_MPI==1
#include "bitpit_communications.hpp"
#endif

#include "ParaTree.hpp"
#include <climits>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iterator>

namespace bitpit {

    // =================================================================================== //
    // NAME SPACES                                                                         //
    // =================================================================================== //
    using namespace std;

    // =================================================================================== //
    // AUXILIARY IMPLEMENTATIONS                                                           //
    // =================================================================================== //

    /*!
        \typedef ParaTree::ExchangeRanges
        \ingroup PABLO

        Defines a set of local octants' ranges that will be exchanged with
        other processes.

        A range is defined as a pair: the first entry is the local index
        referring to the first octant the will be exchanged and the second
        entry is the local index referring to the "past-the-last" octant
        that will be exchanged.

        There is a range for each each process for which an exchange will take
        place.
    */

    /*!
        \struct ParaTree::LoadBalanceRanges
        \ingroup PABLO

        Defines the range of local octants that will be exchanged during a load
        balance.
    */

    /*!
        \enum ParaTree::LoadBalanceRanges::ExchangeAction
        \ingroup PABLO

        The type of exchange action that will be performed on the octants in
        a range.
    */

    /*! Default constructor
    */
    ParaTree::LoadBalanceRanges::LoadBalanceRanges()
        : sendAction(ACTION_UNDEFINED), recvAction(ACTION_UNDEFINED)
    {
    }

    /*! Constructor
     * \param serial controls is the tree is currently serial or parallel
     * \param _sendRanges are the range of local octants that will be sent to
     * other processes
     * \param _recvRanges are the range of local octants that will be received
     * from other processes
     */
    ParaTree::LoadBalanceRanges::LoadBalanceRanges(bool serial, const ExchangeRanges &_sendRanges, const ExchangeRanges &_recvRanges)
        : sendRanges(_sendRanges), recvRanges(_recvRanges)
    {
        if (serial) {
            sendAction = ACTION_DELETE;
            recvAction = ACTION_NONE;
        } else {
            sendAction = ACTION_SEND;
            recvAction = ACTION_RECEIVE;
        }
    }

    /*! Clear the ranges
     */
    void ParaTree::LoadBalanceRanges::clear()
    {
        sendAction = ACTION_UNDEFINED;
        sendRanges.clear();

        recvAction = ACTION_UNDEFINED;
        recvRanges.clear();
    }

    // =================================================================================== //
    // CLASS IMPLEMENTATION                                                                //
    // =================================================================================== //

    const std::string	ParaTree::DEFAULT_LOG_FILE   = "PABLO";

    // =================================================================================== //
    // CONSTRUCTORS AND OPERATORS														   //
    // =================================================================================== //

    /*! Default empty constructor of ParaTree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    ParaTree::ParaTree(const std::string &logfile, MPI_Comm comm )
#else
    ParaTree::ParaTree(const std::string &logfile )
#endif
    {
#if BITPIT_ENABLE_MPI==1
        initialize(logfile, comm);
#else
        initialize(logfile);
#endif
        reset(true);
    }

    /*! Default constructor of ParaTree.
     * It builds one octant with node 0 in the Origin (0,0,0) and side of length 1.
     * \param[in] dim The space dimension of the m_octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    ParaTree::ParaTree(uint8_t dim, const std::string &logfile, MPI_Comm comm )
#else
    ParaTree::ParaTree(uint8_t dim, const std::string &logfile )
#endif
        : m_octree(dim), m_trans(dim)
    {
#if BITPIT_ENABLE_MPI==1
        initialize(dim, logfile, comm);
#else
        initialize(dim, logfile);
#endif
        reset(true);

        printHeader();
    };

    // =============================================================================== //

    /*!
     * Creates a new octree restoring the octree saved in the specified stream.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     * \param stream is the stream to read from
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    ParaTree::ParaTree(std::istream &stream, const std::string &logfile, MPI_Comm comm)
#else
    ParaTree::ParaTree(std::istream &stream, const std::string &logfile)
#endif
    {
#if BITPIT_ENABLE_MPI==1
        initialize(logfile, comm);
#else
        initialize(logfile);
#endif
        restore(stream);
    }

    // =============================================================================== //

    /*! Default Destructor of ParaTree.
     */
    ParaTree::~ParaTree(){
#if BITPIT_ENABLE_MPI==1
        freeComm();
#endif

        printFooter();
    };

    // =============================================================================== //

    /*!
     * Copy constructor of ParaTree.
     * \param[in] other class of type ParaTree
     */
    ParaTree::ParaTree(const ParaTree & other)
        : m_partitionFirstDesc(other.m_partitionFirstDesc),
          m_partitionLastDesc(other.m_partitionLastDesc),
          m_partitionRangeGlobalIdx(other.m_partitionRangeGlobalIdx),
          m_partitionRangeGlobalIdx0(other.m_partitionRangeGlobalIdx0),
          m_globalNumOctants(other.m_globalNumOctants),
          m_nproc(other.m_nproc),
          m_maxDepth(other.m_maxDepth),
          m_treeConstants(other.m_treeConstants),
          m_nofGhostLayers(other.m_nofGhostLayers),
          m_rank(other.m_rank),
          m_octree(other.m_octree),
          m_bordersPerProc(other.m_bordersPerProc),
          m_internals(other.m_internals),
          m_pborders(other.m_pborders),
          m_mapIdx(other.m_mapIdx),
          m_loadBalanceRanges(other.m_loadBalanceRanges),
          m_errorFlag(other.m_errorFlag),
          m_serial(other.m_serial),
          m_tol(other.m_tol),
          m_trans(other.m_trans),
          m_dim(other.m_dim),
          m_periodic(other.m_periodic),
          m_status(other.m_status),
          m_lastOp(other.m_lastOp),
          m_log(other.m_log)
#if BITPIT_ENABLE_MPI==1
          , m_comm(MPI_COMM_NULL)
#endif
    {
#if BITPIT_ENABLE_MPI==1
        setComm(other.m_comm);
#endif
    }

    // =================================================================================== //
    // METHODS
    // =================================================================================== //

#if BITPIT_ENABLE_MPI==1
    /*! Internal function to initialize the communications
     * \param[in] comm The MPI communicator used by the parallel octree.
     */
    void
    ParaTree::_initializeCommunications(MPI_Comm comm) {
        // Set a null communicator
        m_comm = MPI_COMM_NULL;

        //Set the number of ghost layers
        m_nofGhostLayers = 1;

        // Set the real communicator and initialize partition data
        if (comm != MPI_COMM_NULL) {
            setComm(comm);
        } else {
            _initializePartitions();
        }
    }
#endif

    /*! Internal function to initialize the partitions
     *  We always need to initialize the partitions communicator, if MPI is
     *  disabled a dummy initialization will be performed.
     */
    void
    ParaTree::_initializePartitions() {
#if BITPIT_ENABLE_MPI==1
        // Set MPI information
        if (isCommSet()) {
            MPI_Comm_size(m_comm, &m_nproc);
            MPI_Comm_rank(m_comm, &m_rank);
        } else {
            m_rank  = 0;
            m_nproc = 1;
        }
#else
        // Set dummy MPI information
        m_rank  = 0;
        m_nproc = 1;
#endif

        // Create the data structures for storing partion information
        m_partitionFirstDesc.resize(m_nproc);
        m_partitionLastDesc.resize(m_nproc);
        m_partitionRangeGlobalIdx.resize(m_nproc);
        m_partitionRangeGlobalIdx0.resize(m_nproc);
        uint64_t lastDescMorton = m_octree.getLastDescMorton();
        uint64_t firstDescMorton = m_octree.getFirstDescMorton();
        for(int p = 0; p < m_nproc; ++p){
            m_partitionRangeGlobalIdx0[p] = 0;
            m_partitionRangeGlobalIdx[p]  = m_globalNumOctants - 1;
            m_partitionLastDesc[p]  = lastDescMorton;
            m_partitionFirstDesc[p] = firstDescMorton;
        }
    }

    /*! Internal function to initialize a dummy octree
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
    void
    ParaTree::_initialize(uint8_t dim, const std::string &logfile) {
        // The octree is serial
        m_serial = true;

        // Initialize the status
        m_status = 0;

        // Initialize the logger
        initializeLogger(logfile);

        // Set the dimension to a dummy value
        setDim(dim);
    }

    /*! Initialize a dummy octree
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    void
    ParaTree::initialize(const std::string &logfile, MPI_Comm comm) {
#else
    void
    ParaTree::initialize(const std::string &logfile) {
#endif
#if BITPIT_ENABLE_MPI==1
        // Initialize communications
        _initializeCommunications(comm);
#else
        // Initialize partitions
        _initializePartitions();
#endif

        // Initialized the tree
        _initialize(0, logfile);
    }

    /*! Initialize the octree
     * \param[in] dim The space dimension of the octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
    void
#if BITPIT_ENABLE_MPI==1
    ParaTree::initialize(uint8_t dim, const std::string &logfile, MPI_Comm comm) {
#else
    ParaTree::initialize(uint8_t dim, const std::string &logfile) {
#endif
#if BITPIT_ENABLE_MPI==1
        // Initialize communications
        _initializeCommunications(comm);
#else
        // Initialize partitions
        _initializePartitions();
#endif

        // Initialized the tree
        if (dim < 2 || dim > 3) {
            throw std::runtime_error ("Invalid value for the dimension");
        }

        _initialize(dim, logfile);
    }

    /*! Re-initializes the octree
     * \param[in] dim The space dimension of the octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
    void
    ParaTree::reinitialize(uint8_t dim, const std::string &logfile) {
        // Initialize partitions
        _initializePartitions();

        // Initialized the tree
        if (dim < 2 || dim > 3) {
            throw std::runtime_error ("Invalid value for the dimension");
        }

        _initialize(dim, logfile);
    }

    /*! Initialize the logger
     */
    void
    ParaTree::initializeLogger(const std::string &logfile){
        log::manager().create(logfile, false, m_nproc, m_rank);
        m_log = &log::cout(logfile);
    }

    /*! Reset the octree
     */
    void
    ParaTree::reset(){
        reset(true);
    }

    /*! Reset the octree
     */
    void
    ParaTree::reset(bool createRoot){
        m_tol       = 1.0e-14;
        m_serial    = true;
        m_errorFlag = 0;

        m_maxDepth  = 0;

        m_octree.reset(createRoot);
        m_globalNumOctants = getNumOctants();

        m_lastOp = OP_INIT;

        m_bordersPerProc.clear();
        m_internals.clear();
        m_pborders.clear();

        m_loadBalanceRanges.clear();

        std::fill(m_periodic.begin(), m_periodic.end(), false);

        _initializePartitions();
    }

    // =============================================================================== //

    /*! Get the version associated to the binary dumps.
     *
     *  \result The version associated to the binary dumps.
     */
    int ParaTree::getDumpVersion() const
    {
        const int DUMP_VERSION = 1;

        return DUMP_VERSION;
    }

    // =============================================================================== //

    /*! Write the octree to the specified stream.
     *
     *  \param stream is the stream to write to
     *  \param full is the flag for a complete dump with mapping structureof last operation of the tree
     */
    void ParaTree::dump(std::ostream &stream, bool full)
    {
        // Version
        utils::binary::write(stream, getDumpVersion());

        // Tree data
        utils::binary::write(stream, getNproc());

        utils::binary::write(stream, getDim());

        utils::binary::write(stream, getSerial());
        utils::binary::write(stream, getNofGhostLayers());
        utils::binary::write(stream, getMaxDepth());
        utils::binary::write(stream, getStatus());
        utils::binary::write(stream, getBalanceCodimension());

        for (int i = 0; i < m_treeConstants->nFaces; i++) {
            utils::binary::write(stream, getPeriodic(i));
        }

        // Octant data
        uint32_t nOctants = getNumOctants();
        utils::binary::write(stream, nOctants);

        uint32_t nGlobalOctants = getGlobalNumOctants();
        utils::binary::write(stream, nGlobalOctants);

        for (uint32_t i = 0; i < nOctants; i++) {
            const Octant &octant = m_octree.m_octants[i];

            utils::binary::write(stream, octant.getLevel());
            utils::binary::write(stream, octant.getX());
            utils::binary::write(stream, octant.getY());
            utils::binary::write(stream, octant.getZ());
            utils::binary::write(stream, octant.getGhostLayer());

            for (size_t k = 0; k < Octant::INFO_ITEM_COUNT; ++k) {
                utils::binary::write(stream, (bool) octant.m_info[k]);
            }

            utils::binary::write(stream, octant.getBalance());
            utils::binary::write(stream, octant.getMarker());
        }

        // Information about partitioning
        for (int k = 0; k < m_nproc; ++k) {
            utils::binary::write(stream, m_partitionFirstDesc[k]);
        }

        for (int k = 0; k < m_nproc; ++k) {
            utils::binary::write(stream, m_partitionLastDesc[k]);
        }

        for (int k = 0; k < m_nproc; ++k) {
            utils::binary::write(stream, m_partitionRangeGlobalIdx[k]);
        }

        // Extended information (mapping, ...)
        utils::binary::write(stream, full);
        if (full) {
            utils::binary::write(stream, m_lastOp);
            if (m_lastOp == OP_ADAPT_MAPPED){
                for (auto idx : m_mapIdx) {
                    utils::binary::write(stream, idx);
                }

                utils::binary::write(stream, m_octree.m_lastGhostBros.size());
                for (auto lastGhostBrother : m_octree.m_lastGhostBros) {
                    utils::binary::write(stream, lastGhostBrother);
                }
            }
            else if (m_lastOp == OP_LOADBALANCE || m_lastOp == OP_LOADBALANCE_FIRST){
                for (int i = 0; i < m_nproc; ++i) {
                    utils::binary::write(stream, m_partitionRangeGlobalIdx0[i]);
                }
            }
        }

    }

    // =============================================================================== //

    /*! Restore the octree from the specified stream.
     *
     *  \param stream Stream to read from.
     */
    void ParaTree::restore(std::istream &stream)
    {
        // Version
        int version;
        utils::binary::read(stream, version);
        if (version != getDumpVersion()) {
            throw std::runtime_error ("The version of the file does not match the required version");
        }

        // Check if the number of processors matches
        int nProcs;
        utils::binary::read(stream, nProcs);
        if (nProcs != m_nproc) {
            throw std::runtime_error ("The restart was saved with a different number of processors.");
        }

        // Initialize the tree
        uint8_t dimension;
        utils::binary::read(stream, dimension);

        m_octree.initialize(dimension);
        m_trans.initialize(dimension);
        reinitialize(dimension, m_log->getName());
        reset(false);

        // Set tree properties
        utils::binary::read(stream, m_serial);
        utils::binary::read(stream, m_nofGhostLayers);
        utils::binary::read(stream, m_maxDepth);
        utils::binary::read(stream, m_status);

        bool balanceCodimension;
        utils::binary::read(stream, balanceCodimension);
        setBalanceCodimension(balanceCodimension);

        for (int i = 0; i < m_treeConstants->nFaces; i++) {
            bool periodicBorder;
            utils::binary::read(stream, periodicBorder);
            if (periodicBorder){
            	setPeriodic(i);
            }
        }

        // Restore octants
        uint32_t nOctants;
        utils::binary::read(stream, nOctants);
        m_octree.m_sizeOctants = nOctants;

        uint32_t nGlobalOctants;
        utils::binary::read(stream, nGlobalOctants);
        m_globalNumOctants = nGlobalOctants;

        m_octree.m_octants.clear();
        m_octree.m_octants.reserve(nOctants);
        for (uint32_t i = 0; i < nOctants; i++) {
            // Create octant
            uint8_t level;
            utils::binary::read(stream, level);

            uint32_t x;
            utils::binary::read(stream, x);

            uint32_t y;
            utils::binary::read(stream, y);

            uint32_t z;
            utils::binary::read(stream, z);

            Octant octant(false, m_dim, level, x, y, z);

            int ghost;
            utils::binary::read(stream, ghost);
            octant.setGhostLayer(ghost);

            // Set octant info
            for (size_t k = 0; k < Octant::INFO_ITEM_COUNT; ++k) {
                bool bit;
                utils::binary::read(stream, bit);
                octant.m_info.set(k, bit);
            }

            // Set octant 2:1 balance
            bool balance21;
            utils::binary::read(stream, balance21);
            octant.setBalance(balance21);

            // Set marker
            int8_t marker;
            utils::binary::read(stream, marker);
            octant.setMarker(marker);

            // Add octant to the list
            m_octree.m_octants.push_back(std::move(octant));
        }

        m_octree.updateLocalMaxDepth();

        // Set first last descendant
        m_partitionFirstDesc.resize(m_nproc);
        for (int k = 0; k < m_nproc; ++k) {
            uint64_t descendant;
            utils::binary::read(stream, descendant);
            m_partitionFirstDesc[k] = descendant;
        }
        m_octree.m_firstDescMorton = m_partitionFirstDesc[m_rank];

        m_partitionLastDesc.resize(m_nproc);
        for (int k = 0; k < m_nproc; ++k) {
            uint64_t descendant;
            utils::binary::read(stream, descendant);
            m_partitionLastDesc[k] = descendant;
        }
        m_octree.m_lastDescMorton = m_partitionLastDesc[m_rank];

        // Set partitions and parallel information
        m_partitionRangeGlobalIdx.resize(m_nproc);
        for (int k = 0; k < m_nproc; ++k) {
            uint64_t rangeGlobalIdx;
            utils::binary::read(stream, rangeGlobalIdx);
            m_partitionRangeGlobalIdx[k] = rangeGlobalIdx;
        }

#if BITPIT_ENABLE_MPI==1
        if (!m_serial) {
            computeGhostHalo();
        }
#endif

        // Full restore (i.e. restore with mapper of last operation)
        m_mapIdx.clear();

        for (int i = 0; i < m_nproc; ++i) {
            m_partitionRangeGlobalIdx0[i] = 0;
        }

        bool full;
        utils::binary::read(stream, full);
        if (full) {
            utils::binary::read(stream, m_lastOp);
            if (m_lastOp == OP_ADAPT_MAPPED) {
                m_mapIdx.resize(m_octree.m_octants.size());
                for (size_t i = 0; i < m_octree.m_octants.size(); ++i) {
                    utils::binary::read(stream, m_mapIdx[i]);
                }

                size_t lastGhostBrosSize;
                utils::binary::read(stream, lastGhostBrosSize);
                m_octree.m_lastGhostBros.resize(lastGhostBrosSize);
                for (size_t i = 0; i < lastGhostBrosSize; ++i) {
                    utils::binary::read(stream, m_octree.m_lastGhostBros[i]);
                }
            }
            else if (m_lastOp == OP_LOADBALANCE || m_lastOp == OP_LOADBALANCE_FIRST){
                for (int i = 0; i < m_nproc; ++i) {
                    utils::binary::read(stream, m_partitionRangeGlobalIdx0[i]);
                }
            }
        }
        else {
            m_lastOp = OP_INIT;
        }
    }

    // =============================================================================== //

    /*! Print the initial PABLO header.
     */
    void
    ParaTree::printHeader(){
        (*m_log) << log::context("PABLO");
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << "- PABLO PArallel Balanced Linear Octree -" << endl;
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " " << endl;
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << "- PABLO restart -" << endl;
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " Number of proc	:	" + to_string(static_cast<unsigned long long>(m_nproc)) << endl;
        (*m_log) << " Dimension		:	" + to_string(static_cast<unsigned long long>(m_dim)) << endl;
        (*m_log) << " Max allowed level	:	" + to_string(static_cast<unsigned long long>(TreeConstants::MAX_LEVEL)) << endl;
        (*m_log) << " Number of octants	:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " " << endl;
    }

    /*! Print the final PABLO footer.
     */
    void
    ParaTree::printFooter(){
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << "--------------- R.I.P. PABLO ----------------" << endl;
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << "---------------------------------------------" << endl;
    }

    // =================================================================================== //
    // BASIC GET/SET METHODS
    // =================================================================================== //

    /*! Get the dimension of the octree.
     * \return Dimension of the octree (2D/3D).
     */
    uint8_t
    ParaTree::getDim() const {
        return m_dim;
    };

    /*! Get the global number of octants.
     * \return Global number of octants.
     */
    uint64_t
    ParaTree::getGlobalNumOctants() const {
        return m_globalNumOctants;
    };

    /*! Get if the octree is serial.
     * \return Is the octree serial?.
     */
    bool
    ParaTree::getSerial() const {
        return m_serial;
    };

    /*! Get if the octree is parallel.
     * \return Is the octree distributed?.
     */
    bool
    ParaTree::getParallel() const {
        return (!m_serial);
    };

    /*! Get the rank of local process.
     * \return Rank of local process.
     */
    int
    ParaTree::getRank() const {
        return m_rank;
    };

    /*! Get the total number of processes.
     * \return Number of processes.
     */
    int
    ParaTree::getNproc() const {
        return m_nproc;
    };

    /*!Get the logger.
     * \return Reference to logger object.
     */
    Logger&
    ParaTree::getLog(){
        return (*m_log);
    }

    /*!Get the last operation perforfmed by the octree (initialization, adapt (mapped or unmapped), loadbalance (first or not).
     * \return Last operation performed by the octree.
     */
    ParaTree::Operation
    ParaTree::getLastOperation() const {
        return m_lastOp;
    }

#if BITPIT_ENABLE_MPI==1
    /*! Set the MPI communicator to be used for parallel communications.
     * The tree will not use the specified communicator directly, instead a
     * duplicates is ceated.
     * \param communicator is the communicator to be used for parallel
     * communications.
     */
    void
    ParaTree::setComm(MPI_Comm communicator)
    {
        // Communication can be set just once
        if (isCommSet()) {
            throw std::runtime_error ("PABLO communicator can be set just once");
        }

        // The communicator has to be valid
        if (communicator == MPI_COMM_NULL) {
            throw std::runtime_error ("PABLO communicator is not valid");
        }

        // Create a copy of the user-specified communicator
        //
        // No library routine should use MPI_COMM_WORLD as the communicator;
        // instead, a duplicate of a user-specified communicator should always
        // be used.
        MPI_Comm_dup(communicator, &m_comm);

        // Initialize partition data
        _initializePartitions();
    }

    /*! Set the MPI communicator to be used for parallel communications.
     * If the communicator is already set, it will be replaced with the new
     * one only if the two communicators are equivalent, i.e., the rank
     * of the processes have to be the same in both communicators.
     * The tree will not use the specified communicator directly, instead a
     * duplicates is ceated.
     * Previous communicator will be freed or returned depending on
     * the received arguments.
     * \param communicator is the communicator to be used for parallel
     * communications.
     * \param[in,out] previousCommunicator if a valid pointer is received,
     * on output this parameter will contain the previous communicator.
     * If this patameter is a null pointer, the previous communicator
     * will be freed.
     */
    void
    ParaTree::replaceComm(MPI_Comm communicator, MPI_Comm *previousCommunicator)
    {
        // The communicator has to be already set
        if (!isCommSet()) {
            throw std::runtime_error ("PABLO communicator is not set");
        }

        // Check if the communicator is valid
        //
        // The communicator should be equaivalent to the one currently set,
        // i.e., the rank of the processes have to be the same in both
        // communicators.
        int updatedRank;
        MPI_Comm_rank(communicator, &updatedRank);
        int currentRank = getRank();

        int isValid = (updatedRank == currentRank) ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &isValid, 1, MPI_INT, MPI_LAND, m_comm);
        if (!isValid) {
            throw std::runtime_error ("New communicator is not valid");
        }

        // Handle previous communicator
        if (previousCommunicator) {
            *previousCommunicator = m_comm;
            m_comm = MPI_COMM_NULL;
        } else {
            freeComm();
        }

        // Set the communicator
        setComm(communicator);
    }

    /*!
     * Free the MPI communicator
     */
    void
    ParaTree::freeComm() {
        if (!isCommSet()) {
            return;
        }

        // Free MPI communicator
        int finalizedCalled;
        MPI_Finalized(&finalizedCalled);
        if (finalizedCalled) {
            return;
        }

        MPI_Comm_free(&m_comm);
    }

    /*! Check if the communicator to be used for parallel communications has
     * already been set.
     * \return Returns true if the communicator has been set, false otherwise.
     */
    bool
    ParaTree::isCommSet() const {
        return (getComm() != MPI_COMM_NULL);
    }

    /*! Get thecommunicator used by octree between processes.
     * \return MPI Communicator.
     */
    MPI_Comm
    ParaTree::getComm() const {
        return m_comm;
    };
#endif

    /*! Get a constant reference to the vector containing the
     * global index of the last existing octant in each processor.
     * \return A constant reference to the vector containing the
     * global index of the last existing octant in each processor.
     */
    const std::vector<uint64_t> &
    ParaTree::getPartitionRangeGlobalIdx() const {
        return m_partitionRangeGlobalIdx;
    };

    /*! Get a constant reference to the vector containing the
     * Morton number of the first octant on each processor.
     * \return Constant reference to the vector containing the
     * Morton number of the first octant on each processor.
     */
    const std::vector<uint64_t> &
    ParaTree::getPartitionFirstDesc() const {
        return m_partitionFirstDesc;
    };

    /*! Get a constant reference to the vector containing the
     * Morton number of the last possible octant on each processor.
     * \return Constant reference to the vector containing the
     * Morton number of the last possible octant on each processor.
     */
    const std::vector<uint64_t> &
    ParaTree::getPartitionLastDesc() const {
        return m_partitionLastDesc;
    };

    /*! Get the coordinates of the origin of the octree.
     * \return Coordinates of the origin.
     */
    darray3
    ParaTree::getOrigin() const {
        return m_trans.m_origin;
    };

    /*! Get the coordinate X of the origin of the octree.
     * \return Coordinate X of the origin.
     */
    double
    ParaTree::getX0() const {
        return m_trans.m_origin[0];
    };

    /*! Get the coordinate Y of the origin of the octree.
     * \return Coordinate Y of the origin.
     */
    double
    ParaTree::getY0() const {
        return m_trans.m_origin[1];
    };

    /*! Get the coordinate Z of the origin of the octree.
     * \return Coordinate Z of the origin.
     */
    double
    ParaTree::getZ0() const {
        return m_trans.m_origin[2];
    };

    /*! Get the length of the domain.
     * \return Length of the octree.
     */
    double
    ParaTree::getL() const {
        return m_trans.m_L;
    };

    /*! Get the maximum level of refinement allowed for this octree.
     * \return Maximum refinement level for the octree.
     */
    int
    ParaTree::getMaxLevel() const {
        return TreeConstants::MAX_LEVEL;
    };

    /*! Get the length of the domain in logical domain.
     * \return Length of the octree in logical domain.
     */
    uint32_t
    ParaTree::getMaxLength() const {
        return TreeConstants::MAX_LENGTH;
    }

    /*! Get the number of nodes for each octant (4 for 2D case, 8 for 3D case)
     * \return Number of nodes for octant.
     */
    uint8_t
    ParaTree::getNnodes() const {
        return m_treeConstants->nNodes;
    }

    /*! Get the number of faces for each octant (4 for 2D case, 6 for 3D case)
     * \return Number of faces for octant.
     */
    uint8_t
    ParaTree::getNfaces() const {
        return m_treeConstants->nFaces;
    }

    /*! Get the number of edges for each octant (0 for 2D case, 12 for 3D case)
     * \return Number of edges for octant.
     */
    uint8_t
    ParaTree::getNedges() const {
        return m_treeConstants->nEdges;
    }

    /*! Get the number of possible children for each octant (4 for 2D case, 8 for 3D case)
     * \return Number of children for octant.
     */
    uint8_t
    ParaTree::getNchildren() const {
        return m_treeConstants->nChildren;
    }

    /*! Get the number of nodes for each face of an octant (2 for 2D case, 4 for 3D case)
     * \return Number of nodes for face for an octant.
     */
    uint8_t
    ParaTree::getNnodesperface() const {
        return m_treeConstants->nNodesPerFace;
    }

    /*! Get the components (in logical domain) of the 6 normals to the faces of an octant (for the 2D case consider only the first 4)
     * \param[out] normals Normals array[6][3] to the faces of an octant.
     */
    void
    ParaTree::getNormals(int8_t normals[6][3]) const {
        for (int i=0; i<6; i++){
            for (int j=0; j<3; j++){
                normals[i][j] = m_treeConstants->normals[i][j];
            }
        }
    }

    /*! Get the indices of the faces of virtual octants opposed to the 6 faces of an octant
     * (for the 2D case consider only the first 4 terms).
     * \param[out] oppface Opposed faces array[6] to the faces of an octant (oppface[i] = Index of the
     * face of an octant neighbour through the i-th face of the current octant).
     */
    void
    ParaTree::getOppface(uint8_t oppface[6]) const {
        for (int j=0; j<6; j++){
            oppface[j] = m_treeConstants->oppositeFace[j];
        }
    }

    /*! Get the face-node connectivity for 6 faces (in 2D case consider only the first 4 terms).
     * \param[out] facenode Connectivity face-node. facenode[i][0:1] = local indices of nodes
     * of the i-th face of an octant.
     */
    void
    ParaTree::getFacenode(uint8_t facenode[6][4]) const {
        for (int i=0; i<6; i++){
            for (int j=0; j<4; j++){
                facenode[i][j] = m_treeConstants->faceNode[i][j];
            }
        }
    }

    /*! Get the node-face connectivity for 8 nodes (in 2D case consider only the first 4 terms).
     * \param[out] nodeface Connectivity node-face. nodeface[i][0:1] = local indices of faces
     * sharing the i-th node of an octant.
     */
    void
    ParaTree::getNodeface(uint8_t nodeface[8][3]) const {
        for (int i=0; i<8; i++){
            for (int j=0; j<3; j++){
                nodeface[i][j] = m_treeConstants->nodeFace[i][j];
            }
        }
    }

    /*! Get the edge-face connectivity for 12 edge (in 2D case not to be considered at all).
     * \param[out] edgeface Connectivity edge-face. edgeface[i][0:1] = local indices of
     * faces sharing the i-th edge of an octant.
     */
    void
    ParaTree::getEdgeface(uint8_t edgeface[12][2]) const {
        for (int i=0; i<12; i++){
            for (int j=0; j<2; j++){
                edgeface[i][j] = m_treeConstants->edgeFace[i][j];
            }
        }
    }

    /*!Get the normals of the nodes (in 2D case consider only the first 4).
     * \param[out] nodecoeffs Components (x,y,z) of the "normals" of the nodes.
     */
    void
    ParaTree::getNodecoeffs(int8_t nodecoeffs[8][3]) const {
        for (int i=0; i<8; i++){
            nodecoeffs[i][2] = 0;
            for (int j=0; j<m_dim; j++){
                nodecoeffs[i][j] = m_treeConstants->nodeCoeffs[i][j];
            }
        }
    }

    /*!Get the normals per edge (in 2D case not to be considered at all).
     * \param[out] edgecoeffs Components (x,y,z) of the "normals" per edge.
     */
    void
    ParaTree::getEdgecoeffs(int8_t edgecoeffs[12][3]) const {
        for (int i=0; i<12; i++){
            for (int j=0; j<3; j++){
                edgecoeffs[i][j] = m_treeConstants->edgeCoeffs[i][j];
            }
        }
    }

    /*! Get the components (in logical domain) of the 6 normals to the faces of an octant (for the 2D case consider only the first 4)
     * \return Pointer to normals array[6][3] to the faces of an octant.
     */
    const int8_t
        (*ParaTree::getNormals() const) [3] {
        return m_treeConstants->normals;
    }

    /*! Get the indices of the faces of virtual octants opposed to the 6 faces of an octant
     * (for the 2D case consider only the first 4 terms).
     * \return Pointer to opposed faces array[6] to the faces of an octant (oppface[i] = Index of the
     * face of an octant neighbour through the i-th face of the current octant).
     */
    const uint8_t
        *ParaTree::getOppface() const {
        return m_treeConstants->oppositeFace;
    }

    /*! Get the face-node connectivity for 6 faces (in 2D case consider only the first 4 terms).
     * \return Pointer to connectivity face-node. facenode[i][0:1] = local indices of nodes
     * of the i-th face of an octant.
     */
    const uint8_t
        (*ParaTree::getFacenode() const)[4] {
        return m_treeConstants->faceNode;
    }

    /*! Get the node-face connectivity for 8 nodes (in 2D case consider only the first 4 terms).
     * \return Pointer to connectivity node-face. nodeface[i][0:1] = local indices of faces
     * sharing the i-th node of an octant.
     */
    const uint8_t
        (*ParaTree::getNodeface() const)[3] {
        return m_treeConstants->nodeFace;
    }

    /*! Get the edge-face connectivity for 12 edge (in 2D case not to be considered at all).
     * \return Pointer to connectivity edge-face. edgeface[i][0:1] = local indices of
     * faces sharing the i-th edge of an octant.
     */
    const uint8_t
        (*ParaTree::getEdgeface() const)[2] {
        return m_treeConstants->edgeFace;
    }

    /*!Get the normals of the nodes (in 2D case consider only the first 4).
     * \return Pointer to components (x,y,z) of the "normals" of the nodes.
     */
    const int8_t
        (*ParaTree::getNodecoeffs() const)[3] {
        return m_treeConstants->nodeCoeffs;
    };

    /*!Get the normals per edge (in 2D case not to be considered at all).
     * \return Pointer to components (x,y,z) of the "normals" per edge.
     */
    const int8_t
        (*ParaTree::getEdgecoeffs() const)[3] {
        return m_treeConstants->edgeCoeffs;
    };

    /*! Get the periodic condition of the boundaries.
     * \return Vector with the periodic conditions (true/false) of each boundary.
     */
    bvector
    ParaTree::getPeriodic() const {
        return m_periodic;
    };

    /*! Get the periodic condition of a target boundary.
     * \param[in] i Index of the target boundary face.
     * \return Boolean with the periodic conditions (true/false) of the target boundary.
     */
    bool
    ParaTree::getPeriodic(uint8_t i) const {
        return m_periodic[i];
    };

    /*!Get the tolerance used in geometric operations.
     */
    double
    ParaTree::getTol() const {
        return m_tol;
    };

    /*! Set the periodic condition of a target boundary (implicitly set the periodic face).
     * \param[in] i Index of the target boundary face (even the opp face will be set).
     */
    void
    ParaTree::setPeriodic(uint8_t i){
        m_periodic[i] = true;
        m_periodic[m_treeConstants->oppositeFace[i]] = true;
        m_octree.setPeriodic(m_periodic);
    };

    /*!Set the tolerance used in geometric operations.
     * \param[in] tol Desired tolerance.
     */
    void
    ParaTree::setTol(double tol){
        m_tol = tol;
    };

    // =================================================================================== //
    // INDEX BASED METHODS
    // =================================================================================== //

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinates X,Y,Z of node 0.
     */
    darray3
    ParaTree::getCoordinates(uint32_t idx) const {
        return m_trans.mapCoordinates(m_octree.m_octants[idx].getCoordinates());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate X of node 0.
     */
    double
    ParaTree::getX(uint32_t idx) const {
        return m_trans.mapX(m_octree.m_octants[idx].getX());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate Y of node 0.
     */
    double
    ParaTree::getY(uint32_t idx) const {
        return m_trans.mapY(m_octree.m_octants[idx].getY());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate Z of node 0.
     */
    double
    ParaTree::getZ(uint32_t idx) const {
        return m_trans.mapZ(m_octree.m_octants[idx].getZ());
    }

    /*! Get the size of an octant, i.e. the side length.
     * \param[in] idx Local index of target octant.
     * \return Size of octant.
     */
    double
    ParaTree::getSize(uint32_t idx) const {
        return m_trans.mapSize(m_octree.m_octants[idx].getSize());
    }

    /*! Get the area of an octant (for 2D case the same value of getSize).
     * \param[in] idx Local index of target octant.
     * \return Area of octant.
     */
    double
    ParaTree::getArea(uint32_t idx) const {
        return m_trans.mapArea(m_octree.m_octants[idx].getArea());
    }

    /*! Get the volume of an octant.
     * \param[in] idx Local index of target octant.
     * \return Volume of octant.
     */
    double
    ParaTree::getVolume(uint32_t idx) const {
        return m_trans.mapVolume(m_octree.m_octants[idx].getVolume());
    }

    /*! Get the coordinates of the center of an octant.
     * \param[in] idx Local index of target octant.
     * \param[out] center Coordinates of the center of octant.
     */
    void
    ParaTree::getCenter(uint32_t idx, darray3& center) const {
        darray3 center_ = m_octree.m_octants[idx].getCenter();
        m_trans.mapCenter(center_, center);
    }

    /*! Get the coordinates of the center of an octant.
     * \param[in] idx Local index of target octant.
     * \return center Coordinates of the center of octant.
     */
    darray3
    ParaTree::getCenter(uint32_t idx) const {
        darray3 center;
        darray3 center_ = m_octree.m_octants[idx].getCenter();
        m_trans.mapCenter(center_, center);
        return center;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the target face.
     * \return center Coordinates of the center of the iface-th face af octant.
     */
    darray3
    ParaTree::getFaceCenter(uint32_t idx, uint8_t iface) const {
        darray3 center;
        darray3 center_ = m_octree.m_octants[idx].getFaceCenter(iface);
        m_trans.mapCenter(center_, center);
        return center;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the target face.
     * \param[out] center Coordinates of the center of the iface-th face af octant.
     */
    void
    ParaTree::getFaceCenter(uint32_t idx, uint8_t iface, darray3& center) const {
        darray3 center_ = m_octree.m_octants[idx].getFaceCenter(iface);
        m_trans.mapCenter(center_, center);
    }

    /*! Get the coordinates of single node of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] inode Index of the target node.
     * \return Coordinates of the inode-th node of octant.
     */
    darray3
    ParaTree::getNode(uint32_t idx, uint8_t inode) const {
        darray3 node;
        u32array3 node_ = m_octree.m_octants[idx].getNode(inode);
        m_trans.mapNode(node_, node);
        return node;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] inode Index of the target node.
     * \param[out] node Coordinates of the inode-th node of octant.
     */
    void
    ParaTree::getNode(uint32_t idx, uint8_t inode, darray3& node) const {
        u32array3 node_ = m_octree.m_octants[idx].getNode(inode);
        m_trans.mapNode(node_, node);
    }

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] idx Local index of target octant.
     * \param[out] nodes Coordinates of the nodes of octant.
     */
    void
    ParaTree::getNodes(uint32_t idx, darr3vector & nodes) const {
        u32arr3vector nodes_;
        m_octree.m_octants[idx].getNodes(nodes_);
        m_trans.mapNodes(nodes_, nodes);
    }

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] idx Local index of target octant.
     * \return nodes Coordinates of the nodes of octant.
     */
    darr3vector
    ParaTree::getNodes(uint32_t idx) const{
        darr3vector nodes;
        u32arr3vector nodes_;
        m_octree.m_octants[idx].getNodes(nodes_);
        m_trans.mapNodes(nodes_, nodes);
        return nodes;
    }

    /*! Get the normal of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the face for normal computing.
     * \param[out] normal Coordinates of the normal of face.
     */
    void
    ParaTree::getNormal(uint32_t idx, uint8_t iface, darray3 & normal) const {
        i8array3 normal_;
        m_octree.m_octants[idx].getNormal(iface, normal_, m_treeConstants->normals);
        m_trans.mapNormals(normal_, normal);
    }

    /*! Get the normal of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the face for normal computing.
     * \return normal Coordinates of the normal of face.
     */
    darray3
    ParaTree::getNormal(uint32_t idx, uint8_t iface) const {
        darray3 normal;
        i8array3 normal_;
        m_octree.m_octants[idx].getNormal(iface, normal_, m_treeConstants->normals);
        m_trans.mapNormals(normal_, normal);
        return normal;
    }

    /*! Get the refinement marker of an octant.
     * \param[in] idx Local index of target octant.
     * \return Marker of octant.
     */
    int8_t
    ParaTree::getMarker(uint32_t idx) const {
        return m_octree.getMarker(idx);
    };

    /*! Get the refinement marker of an octant after a preadapt.
    * \param[in] idx Local index of target octant.
    * \return Marker of octant.
    * NOTE: if a last operation is not preadapt, it calls preadapt method.
    */
    int8_t
    ParaTree::getPreMarker(uint32_t idx){
        if (m_lastOp != OP_PRE_ADAPT) {
            throw std::runtime_error("Last operation different from preadapt, unable to call getPreMarker function");
        }

        return m_octree.getMarker(idx);
    };

    /*! Get the level of an octant.
     * \param[in] idx Local index of target octant.
     * \return Level of octant.
     */
    uint8_t
    ParaTree::getLevel(uint32_t idx) const {
        return m_octree.getLevel(idx);
    };

    /** Compute the Morton index of an octant (without level).
     * \param[in] idx Local index of target octant.
     * \return morton Morton index of the octant.
     */
    uint64_t
    ParaTree::getMorton(uint32_t idx) const {
        return m_octree.computeMorton(idx);
    };

    /** Compute the Morton index of the specified node of the idx-th octant
     * (without level).
     * \param[in] idx Local index of the target octant.
     * \param[in] inode Index of the target node.
     * \return Morton index of the node.
     */
    uint64_t
    ParaTree::getNodeMorton(uint32_t idx, uint8_t inode) const {
        return m_octree.computeNodeMorton(idx, inode);
    };

    /*! Get the balancing condition of an octant.
     * \param[in] idx Local index of target octant.
     * \return Has octant to be balanced?
     */
    bool
    ParaTree::getBalance(uint32_t idx) const {
        return m_octree.getBalance(idx);
    };

    /*! Get the bound condition of the face of the octant
     * \param[in] idx Local index of the target octant
     * \param[in] iface Index of the face
     * \return Is the face a boundary face?
     */
    bool
    ParaTree::getBound(uint32_t idx, uint8_t iface) const {
        return m_octree.m_octants[idx].getBound(iface);
    }

    /*! Get the bound condition of the face of the octant
     * \param[in] idx Local index of the target octant
     * \return Is the octant a boundary octant?
     */
    bool
    ParaTree::getBound(uint32_t idx) const {
        return m_octree.m_octants[idx].getBound();
    }

    /*! Get the partition bound condition of the face of the octant
     * \param[in] idx Local index of the target octant
     * \param[in] iface Index of the face
     * \return Is the face a partition boundary face?
     */
    bool
    ParaTree::getPbound(uint32_t idx, uint8_t iface) const {
        return m_octree.m_octants[idx].getPbound(iface);
    }

    /*! Get the partition bound condition of the face of the octant
     * \param[in] idx Local index of the target octant
     * \return Is the octant a partition boundary octant?
     */
    bool
    ParaTree::getPbound(uint32_t idx) const {
        return m_octree.m_octants[idx].getPbound();
    }

    /*! Get if the octant is new after refinement.
     * \param[in] idx Local index of target octant.
     * \return Is octant new?
     */
    bool
    ParaTree::getIsNewR(uint32_t idx) const {
        return m_octree.m_octants[idx].getIsNewR();
    };

    /*! Get if the octant is new after coarsening.
     * \param[in] idx Local index of target octant.
     * \return Is octant new?
     */
    bool
    ParaTree::getIsNewC(uint32_t idx) const {
        return m_octree.m_octants[idx].getIsNewC();
    };

    /*! Get the global index of an octant.
     * \param[in] idx Local index of target octant.
     * \return Global index of octant.
     */
    uint64_t
    ParaTree::getGlobalIdx(uint32_t idx) const {
        if (m_rank){
            return m_partitionRangeGlobalIdx[m_rank-1] + uint64_t(idx + 1);
        }
        else{
            return uint64_t(idx);
        };
    };

    /*! Get the local index of an octant.
     * \param[in] gidx Global index of target octant.
     * \param[in] rank Rank of the process owning the octant as a local one. getOwnerRank can be used to get the rank of the process owning the octant locally.
     * \return Local index of octant.
     */
    uint32_t
    ParaTree::getLocalIdx(uint64_t gidx,int rank) const {
        if (rank){
            return uint32_t(gidx - 1 - m_partitionRangeGlobalIdx[rank-1]);
        }
        else{
            return uint32_t(gidx);
        };
    };

    /*! Get the local index of an octant and the rank owning the octant.
     * \param[in] gidx Global index of target octant.
     * \param[out] rank Rank of the process owning the octant as a local one.
     * \param[out] Local index of octant.
     */
    void
    ParaTree::getLocalIdx(uint64_t gidx, uint32_t & lidx,int & rank) const {
        rank = getOwnerRank(gidx);
        lidx = getLocalIdx(gidx,rank);
    };


    /*! Get the local index of an octant.
     * \param[in] gidx Global index of target octant.
     * \return Local index of octant.
     */
    uint32_t
    ParaTree::getLocalIdx(uint64_t gidx) const {
        if (m_rank){
            return uint32_t(gidx - 1 - m_partitionRangeGlobalIdx[m_rank-1]);
        }
        else{
            return uint32_t(gidx);
        };
    };

    /*! Get the local index of a ghost octant.
     * \param[in] gidx Global index of target octant.
     * \return Local index of the ghost octant.
     */
    uint32_t
    ParaTree::getGhostLocalIdx(uint64_t gidx) const {

        uint32_t index;
        typename u64vector::const_iterator findResult;
        findResult = std::find(m_octree.m_globalIdxGhosts.begin(),m_octree.m_globalIdxGhosts.end(),gidx);
        if(findResult != m_octree.m_globalIdxGhosts.end()){
            index = std::distance(m_octree.m_globalIdxGhosts.begin(),findResult);
        }
        else{
            index = std::numeric_limits<uint32_t>::max();
        }
        return index;
    };


    /*! Get the global index of a ghost octant.
     * \param[in] idx Local index of target ghost octant.
     * \return Global index of ghost octant.
     */
    uint64_t
    ParaTree::getGhostGlobalIdx(uint32_t idx) const {
        if (idx<m_octree.m_sizeGhosts){
            return m_octree.m_globalIdxGhosts[idx];
        };
        return uint64_t(m_octree.m_sizeGhosts);
    };

    /*! Returns true if the specified global index belongs to the current
     *  process
     * \param[in] gidx Global index of target octant.
     * \return Returns true if the specified global index belongs to the
     * current process.
     */
    bool
    ParaTree::isInternal(uint64_t gidx) const {
        if (gidx > m_partitionRangeGlobalIdx[m_rank]){
            return false;
        }

        if (m_rank != 0 && gidx <= m_partitionRangeGlobalIdx[m_rank - 1]){
            return false;
        }

        return true;
    };

    /*! Get the persistent index of an octant.
     * \param[in] idx Local index of target octant.
     * \return Persistent index of octant,
     * i.e. a bitset composed by Morton index and level of octant.
     */
    bitset<72>
    ParaTree::getPersistentIdx(uint32_t idx) const {
        bitset<72> persistent = getMorton(idx);
        bitset<72> level = getLevel(idx);
        persistent <<= 8;
        persistent |= level;
        return persistent;
    };

    /*! Set the refinement marker of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
     */
    void
    ParaTree::setMarker(uint32_t idx, int8_t marker){
        if (m_lastOp == OP_PRE_ADAPT) {
            throw std::runtime_error("It is not possible to update the tree until the adaption is completed");
        }

        m_octree.setMarker(idx, marker);
    };

    /*! Set the balancing condition of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
     */
    void
    ParaTree::setBalance(uint32_t idx, bool balance){
        if (m_lastOp == OP_PRE_ADAPT) {
            throw std::runtime_error("It is not possible to update the tree until the adaption is completed");
        }

        m_octree.setBalance(idx, balance);
    };

    // =================================================================================== //
    // POINTER BASED METHODS
    // =================================================================================== //

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinates of node 0.
     */
    darray3
    ParaTree::getCoordinates(const Octant* oct) const {
        return m_trans.mapCoordinates(oct->getCoordinates());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate X of node 0.
     */
    double
    ParaTree::getX(const Octant* oct) const {
        return m_trans.mapX(oct->getX());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate Y of node 0.
     */
    double
    ParaTree::getY(const Octant* oct) const {
        return m_trans.mapY(oct->getY());
    }

    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate Z of node 0.
     */
    double
    ParaTree::getZ(const Octant* oct) const {
        return m_trans.mapZ(oct->getZ());
    }

    /*! Get the size of an octant, i.e. the side length.
     * \param[in] oct Pointer to the target octant
     * \return Size of octant.
     */
    double
    ParaTree::getSize(const Octant* oct) const {
        return m_trans.mapSize(oct->getSize());
    }

    /*! Get the area of an octant (for 2D case the same value of getSize).
     * \param[in] oct Pointer to the target octant
     * \return Area of octant.
     */
    double
    ParaTree::getArea(const Octant* oct) const {
        return m_trans.mapArea(oct->getArea());
    }

    /*! Get the volume of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Volume of octant.
     */
    double
    ParaTree::getVolume(const Octant* oct) const {
        return m_trans.mapVolume(oct->getVolume());
    }

    /*! Get the coordinates of the center of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[out] center Coordinates of the center of octant.
     */
    void
    ParaTree::getCenter(const Octant* oct, darray3& center) const {
        darray3 center_ = oct->getCenter();
        m_trans.mapCenter(center_, center);
    }

    /*! Get the coordinates of the center of an octant.
     * \param[in] oct Pointer to the target octant
     * \return center Coordinates of the center of octant.
     */
    darray3
    ParaTree::getCenter(const Octant* oct) const {
        darray3 center;
        darray3 center_ = oct->getCenter();
        m_trans.mapCenter(center_, center);
        return center;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the target face.
     * \return center Coordinates of the center of the iface-th face af octant.
     */
    darray3
    ParaTree::getFaceCenter(const Octant* oct, uint8_t iface) const {
        darray3 center;
        darray3 center_ = oct->getFaceCenter(iface);
        m_trans.mapCenter(center_, center);
        return center;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the target face.
     * \param[out] center Coordinates of the center of the iface-th face af octant.
     */
    void
    ParaTree::getFaceCenter(const Octant* oct, uint8_t iface, darray3& center) const {
        darray3 center_ = oct->getFaceCenter(iface);
        m_trans.mapCenter(center_, center);
    }

    /*! Get the coordinates of single node of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] inode Index of the target node.
     * \return Coordinates of the inode-th node of octant.
     */
    darray3
    ParaTree::getNode(const Octant* oct, uint8_t inode) const {
        darray3 node;
        u32array3 node_ = oct->getNode(inode);
        m_trans.mapNode(node_, node);
        return node;
    }

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] inode Index of the target node.
     * \param[out] node Coordinates of the inode-th node of octant.
     */
    void
    ParaTree::getNode(const Octant* oct, uint8_t inode, darray3& node) const {
        u32array3 node_ = oct->getNode(inode);
        m_trans.mapNode(node_, node);
    }

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[out] nodes Coordinates of the nodes of octant.
     */
    void
    ParaTree::getNodes(const Octant* oct, darr3vector & nodes) const {
        u32arr3vector nodes_;
        oct->getNodes(nodes_);
        m_trans.mapNodes(nodes_, nodes);
    }

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] oct Pointer to the target octant
     * \return nodes Coordinates of the nodes of octant.
     */
    darr3vector
    ParaTree::getNodes(const Octant* oct) const {
        darr3vector nodes;
        u32arr3vector nodes_;
        oct->getNodes(nodes_);
        m_trans.mapNodes(nodes_, nodes);
        return nodes;
    }

    /*! Get the normal of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face for normal computing.
     * \param[out] normal Coordinates of the normal of face.
     */
    void
    ParaTree::getNormal(const Octant* oct, uint8_t iface, darray3 & normal) const {
        i8array3 normal_;
        oct->getNormal(iface, normal_, m_treeConstants->normals);
        m_trans.mapNormals(normal_, normal);
    }

    /*! Get the normal of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face for normal computing.
     * \return normal Coordinates of the normal of face.
     */
    darray3
    ParaTree::getNormal(const Octant* oct, uint8_t iface) const {
        darray3 normal;
        i8array3 normal_;
        oct->getNormal(iface, normal_, m_treeConstants->normals);
        m_trans.mapNormals(normal_, normal);
        return normal;
    }

    /*! Get the refinement marker of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Marker of octant.
     */
    int8_t
    ParaTree::getMarker(const Octant* oct) const {
        return oct->getMarker();
    };

    /*! Get the refinement marker of an octant after a preadapt.
     * \param[in] oct Pointer to the target octant
     * \return Marker of octant.
     * NOTE: if a last operation is not preadapt, it calls preadapt method.
     */
    int8_t
    ParaTree::getPreMarker(Octant* oct){
        if (m_lastOp != OP_PRE_ADAPT) {
            throw std::runtime_error("Last operation different from preadapt, unable to call getPreMarker function");
        }

        return oct->getMarker();
    };

    /*! Get the level of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Level of octant.
     */
    uint8_t
    ParaTree::getLevel(const Octant* oct) const {
        return oct->getLevel();
    };

    /** Compute the Morton index of an octant (without level).
     * \param[in] oct Pointer to the target octant
     * \return morton Morton index of the octant.
     */
    uint64_t
    ParaTree::getMorton(const Octant* oct) const {
        return oct->computeMorton();
    };

    /** Compute the Morton index of the specified node of an octant (without
     * level).
     * \param[in] oct Pointer to the target octant
     * \param[in] inode Index of the target node.
     * \return Morton index of the node.
     */
    uint64_t
    ParaTree::getNodeMorton(const Octant* oct, uint8_t inode) const {
        return oct->computeNodeMorton(inode);
    };

    /*! Get the balancing condition of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Has octant to be balanced?
     */
    bool
    ParaTree::getBalance(const Octant* oct) const {
        return oct->getBalance();
    };

    /*! Get the bound condition of the face of the octant
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face
     * \return Is the face a boundary face?
     */
    bool
    ParaTree::getBound(const Octant* oct, uint8_t iface) const {
        return oct->getBound(iface);
    }

    /*! Get the bound condition of the octant
     * \param[in] oct Pointer to the target octant
     * \return Is the octant a boundary octant?
     */
    bool
    ParaTree::getBound(const Octant* oct) const {
        return oct->getBound();
    }

    /*! Get the partition bound condition of the face of the octant
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face
     * \return Is the face a partition boundary face?
     */
    bool
    ParaTree::getPbound(const Octant* oct, uint8_t iface) const {
        return oct->getPbound(iface);
    }

    /*! Get the partition bound condition of the face of the octant
     * \param[in] oct Pointer to the target octant
     * \return Is the octant a partition boundary octant?
     */
    bool
    ParaTree::getPbound(const Octant* oct) const {
        return oct->getPbound();
    }

    /*! Get if the octant is new after refinement.
     * \param[in] oct Pointer to the target octant
     * \return Is octant new?
     */
    bool
    ParaTree::getIsNewR(const Octant* oct) const {
        return oct->getIsNewR();
    };

    /*! Get if the octant is new after coarsening.
     * \param[in] oct Pointer to the target octant
     * \return Is octant new?
     */
    bool
    ParaTree::getIsNewC(const Octant* oct) const {
        return oct->getIsNewC();
    };

    /*! Get the local index of an octant.
     * \param[in] oct Pointer to target octant.
     * \return Local index of octant.
     */
    uint32_t
    ParaTree::getIdx(const Octant* oct) const {
#if BITPIT_ENABLE_MPI==1
        if (getIsGhost(oct)){
            return m_octree.findGhostMorton(oct->computeMorton());
        }
#endif
        return m_octree.findMorton(oct->computeMorton());
    };

    /*! Get the global index of an octant.
     * \param[in] oct Pointer to target octant.
     * \return Global index of octant.
     */
    uint64_t
    ParaTree::getGlobalIdx(const Octant* oct) const {
        uint32_t idx = getIdx(oct);
#if BITPIT_ENABLE_MPI==1
        if (getIsGhost(oct)){
            return m_octree.m_globalIdxGhosts[idx];
        }
#endif
        if (m_rank){
            return m_partitionRangeGlobalIdx[m_rank-1] + uint64_t(idx + 1);
        }
        return uint64_t(idx);
    };

    /*! Get the persistent index of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Persistent index of octant,
     * i.e. a bitset composed by Morton index and level of octant.
     */
    bitset<72>
    ParaTree::getPersistentIdx(const Octant* oct) const {
        bitset<72> persistent = getMorton(oct);
        bitset<72> level = getLevel(oct);
        persistent <<= 8;
        persistent |= level;
        return persistent;
    };

    /*! Set the refinement marker of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] marker Refinement marker of octant (n=n refinement in adapt, -n=n coarsening in adapt, default=0).
     */
    void
    ParaTree::setMarker(Octant* oct, int8_t marker){
        if (m_lastOp == OP_PRE_ADAPT) {
            throw std::runtime_error("It is not possible to update the tree until the adaption is completed");
        }

        oct->setMarker(marker);
    };

    /*! Set the balancing condition of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] balance Has octant to be 2:1 balanced in adapting procedure?
     */
    void
    ParaTree::setBalance(Octant* oct, bool balance){
        if (m_lastOp == OP_PRE_ADAPT) {
            throw std::runtime_error("It is not possible to update the tree until the adaption is completed");
        }

        oct->setBalance(balance);
    };

    // =================================================================================== //
    // LOCAL TREE GET/SET METHODS
    // =================================================================================== //

    /*! Get the status label of the octree.
     * 	\return Status label of the octree.
     */
    uint64_t
    ParaTree::getStatus()const {
        return m_status;
    }

    /*! Get the local number of octants.
     * \return Local number of octants.
     */
    uint32_t
    ParaTree::getNumOctants() const{
        return m_octree.getNumOctants();
    };

    /*! Get the local number of ghost octants.
     * \return Local number of ghost octants.
     */
    uint32_t
    ParaTree::getNumGhosts() const{
        return m_octree.getNumGhosts();
    };

    /** Get the local number of nodes.
     * \return Local total number of nodes.
     */
    uint32_t
    ParaTree::getNumNodes() const{
        return m_octree.m_nodes.size();
    }

    /*! Get the local depth of the octree.
     * \return Local depth of the octree.
     */
    uint8_t
    ParaTree::getLocalMaxDepth() const{
        return m_octree.getLocalMaxDepth();
    };

    /*! Get the local current minimum size reached by the octree.
     * \return Local current minimum size of the local partition of the octree.
     */
    double
    ParaTree::getLocalMinSize() const {
        uint32_t size = uint32_t(1)<<(TreeConstants::MAX_LEVEL-m_octree.getLocalMaxDepth());
        return m_trans.mapSize(size);
    };

    /*! Get the local current maximum size of the octree.
     * \return Local current maximum size of the local partition of the octree.
     */
    double
    ParaTree::getLocalMaxSize() const {
        uint32_t nocts = getNumOctants();
        double octSize = 0;
        double size = 0;
        for (uint32_t idx = 0; idx < nocts; idx++){
            octSize = getSize(idx);
            if (octSize > size) size = octSize;
        }
        return octSize;
    };

    /*! Get the codimension for 2:1 balancing
     * \return Maximum codimension of the entity through which the 2:1 balance is performed.
     */
    uint8_t
    ParaTree::getBalanceCodimension() const{
        return m_octree.getBalanceCodim();
    };

    /*!Get the first possible descendant with maximum refinement level of the local tree.
     * \return Constant reference to the first finest descendant of the local tree.
     */
     uint64_t ParaTree::getFirstDescMorton() const {
        return m_octree.getFirstDescMorton();
    };

    /*!Get the last possible descendant with maximum refinement level of the local tree.
     * \return Constant reference to the last finest descendant of the local tree.
     */
     uint64_t ParaTree::getLastDescMorton() const {
        return m_octree.getLastDescMorton();
    };

    /*!Get the morton index of the last possible descendant with maximum refinement level of a target octant.
     * \param[in] idx Local index of the target octant.
     * \return Morton index of the last finest descendant of the target octant.
     */
    uint64_t
    ParaTree::getLastDescMorton(uint32_t idx) const {
        return m_octree.m_octants[idx].buildLastDesc().computeMorton();
    };

    /*!Get the begin position for the iterator of the local internal octants.
     * \return Iterator begin of the local internal octants (dereferencing results in a pointer to an octant).
     */
    octantIterator
    ParaTree::getInternalOctantsBegin() {
        return m_internals.begin();
    }

    /*!Get the end position for the iterator of the local internal octants.
     * \return Iterator end of the local internal octants (dereferencing results in a pointer to an octant).
     */
    octantIterator
    ParaTree::getInternalOctantsEnd() {
        return m_internals.end();
    }

    /*!Get the begin position for the iterator of the local border of process octants.
     * \return Iterator begin of the local border of process octants (dereferencing results in a pointer to an octant).
     */
    octantIterator
    ParaTree::getPboundOctantsBegin() {
        return m_pborders.begin();
    }

    /*!Get the end position for the iterator of the local border of process octants.
     * \return Iterator end of the local border of process octants (dereferencing results in a pointer to an octant).
     */
    octantIterator
    ParaTree::getPboundOctantsEnd() {
        return m_pborders.end();
    }

    /*! Set the codimension for 2:1 balancing
     * \param[in] b21codim  Maximum codimension of the entity through which the 2:1 balance is performed (1 = 2:1 balance through edges (default); 2 = 2:1 balance through nodes and edges).
     */
    void
    ParaTree::setBalanceCodimension(uint8_t b21codim){
        if (m_lastOp == OP_PRE_ADAPT) {
            throw std::runtime_error("It is not possible to update the tree until the adaption is completed");
        }

        m_octree.setBalanceCodim(b21codim);
    };

    // =================================================================================== //
    // INTERSECTION GET/SET METHODS
    // =================================================================================== //

    /*! Get the local number of intersections.
     * \return Local number of intersections.
     */
    uint32_t
    ParaTree::getNumIntersections() const {
        return m_octree.m_intersections.size();
    }

    /*! Get a pointer to target intersection.
     * \param[in] idx Local index of intersection.
     * \return Pointer to target intersection.
     */
    Intersection*
    ParaTree::getIntersection(uint32_t idx) {
        if (idx < m_octree.m_intersections.size()){
            return &m_octree.m_intersections[idx];
        }
        return NULL;
    }

    /*! Get the level of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Level of intersection.
     */
    uint8_t
    ParaTree::getLevel(const Intersection* inter) const {
        if(inter->m_finer && inter->m_isghost)
            return m_octree.extractGhostOctant(inter->m_owners[inter->m_finer]).getLevel();
        else
            return m_octree.extractOctant(inter->m_owners[inter->m_finer]).getLevel();
    }

    /*! Get the finer owner octant of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return The finer octant of the owners of intersection (false/true = 0/1).
     */
    bool
    ParaTree::getFiner(const Intersection* inter) const {
        return inter->m_finer;
    }

    /*! Get if an intersection is a boundary domain intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Boundary or not boundary?.
     */
    bool
    ParaTree::getBound(const Intersection* inter) const {
        return inter->getBound();
    }

    /*! Get if an intersection is an intersection between an internal and a ghost element.
     * \param[in] inter Pointer to target intersection.
     * \return Ghost or not ghost?.
     */
    bool
    ParaTree::getIsGhost(const Intersection* inter) const {
        return inter->getIsGhost();
    }

    /*! Get if an intersection is a boundary intersection for a process.
     * \param[in] inter Pointer to target intersection.
     * \return Process boundary or not boundary?.
     */
    bool
    ParaTree::getPbound(const Intersection* inter) const {
        return inter->getPbound();
    }

    /*! Get the face index of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Face index of the finer octant of intersection (owners[getFiner(inter)]).
     */
    uint8_t
    ParaTree::getFace(const Intersection* inter) const {
        return inter->m_iface;
    }

    /*! Get the owner octants of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return A couple of octants owners of intersection.
     */
    u32vector
    ParaTree::getOwners(const Intersection* inter) const {
        u32vector owners(2);
        owners[0] = inter->m_owners[0];
        owners[1] = inter->m_owners[1];
        return owners;
    }

    /*! Get the owner octant of an intersection with inner normal.
     * \param[in] inter Pointer to target intersection.
     * \return Index of the octant owner with inner normal.
     */
    uint32_t
    ParaTree::getIn(const Intersection* inter) const {
        return inter->getIn();
    }

    /*! Get the owner octant of an intersection with outer normal.
     * \param[in] inter Pointer to target intersection.
     * \return Index of the octant owner with outer normal.
     */
    uint32_t
    ParaTree::getOut(const Intersection* inter) const {
        return inter->getOut();
    }

    /*! Get if the owner octant with outer normal is a ghost octant.
     * \param[in] inter Pointer to target intersection.
     * \return Is the octant owner with outer normal a ghost octant?.
     */
    bool
    ParaTree::getOutIsGhost(const Intersection* inter) const {
        return inter->getOutIsGhost();
    }

    /*! Get the size of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Size of intersection.
     */
    double
    ParaTree::getSize(const Intersection* inter) const {
        uint32_t Size;
        if(inter->m_finer && inter->m_isghost)
            Size = m_octree.extractGhostOctant(inter->m_owners[inter->m_finer]).getSize();
        else
            Size = m_octree.extractOctant(inter->m_owners[inter->m_finer]).getSize();
        return m_trans.mapSize(Size);
    }

    /*! Get the area of an intersection (for 2D case the same value of getSize).
     * \param[in] inter Pointer to target intersection.
     * \return Area of intersection.
     */
    double
    ParaTree::getArea(const Intersection* inter) const {
        uint64_t Area;
        if(inter->m_finer && inter->m_isghost)
            Area = m_octree.extractGhostOctant(inter->m_owners[1]).getArea();
        else
            Area = m_octree.extractOctant(inter->m_owners[inter->m_finer]).getArea();
        return m_trans.mapArea(Area);
    }

    /*! Get the coordinates of the center of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the center of intersection.
     */
    darray3
    ParaTree::getCenter(const Intersection* inter) const {
        darray3 center;
        Octant oct(m_dim);
        if(inter->m_finer && inter->m_isghost)
            oct = m_octree.extractGhostOctant(inter->m_owners[inter->m_finer]);
        else
            oct = m_octree.extractOctant(inter->m_owners[inter->m_finer]);
        darray3  center_ = oct.getCenter();
        int sign = ( int(2*((inter->m_iface)%2)) - 1);
        double deplace = double (sign * int(oct.getSize())) / 2;
        center_[inter->m_iface/2] = uint32_t(int(center_[inter->m_iface/2]) + deplace);
        m_trans.mapCenter(center_, center);
        return center;
    }

    /*! Get the coordinates of the nodes of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the nodes of intersection.
     */
    darr3vector
    ParaTree::getNodes(const Intersection* inter) const {
        darr3vector nodes;
        Octant oct(m_dim);
        if(inter->m_finer && inter->m_isghost)
            oct = m_octree.extractGhostOctant(inter->m_owners[inter->m_finer]);
        else
            oct = m_octree.extractOctant(inter->m_owners[inter->m_finer]);
        uint8_t iface = inter->m_iface;
        u32arr3vector nodes_all;
        oct.getNodes(nodes_all);
        u32arr3vector nodes_(m_treeConstants->nNodesPerFace);
        for (int i=0; i<m_treeConstants->nNodesPerFace; i++){
            for (int j=0; j<3; j++){
                nodes_[i][j] = nodes_all[m_treeConstants->faceNode[iface][i]][j];
            }
        }
        m_trans.mapNodesIntersection(nodes_, nodes);
        return nodes;
    }

    /*! Get the normal of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the normal of intersection.
     */
    darray3
    ParaTree::getNormal(const Intersection* inter) const {
        darray3 normal;
        Octant oct(m_dim);
        if(inter->m_finer && inter->m_isghost)
            oct = m_octree.extractGhostOctant(inter->m_owners[inter->m_finer]);
        else
            oct = m_octree.extractOctant(inter->m_owners[inter->m_finer]);
        uint8_t iface = inter->m_iface;
        i8array3 normal_;
        oct.getNormal(iface, normal_, m_treeConstants->normals);
        m_trans.mapNormals(normal_, normal);
        return normal;
    }

    // =================================================================================== //
    // OTHER GET/SET METHODS
    // =================================================================================== //

    /** Get an octant as pointer to the target octant.
     * NOTE: no checks will be performed on the octant index.
     * \param[in] idx Local index of target octant.
     * \return Pointer to target octant.
     */
    Octant*
    ParaTree::getOctant(uint32_t idx) {
        return (m_octree.m_octants.data() + idx);
    };

    /** Get an octant as constant pointer to the target octant.
     * NOTE: no checks will be performed on the octant index.
     * \param[in] idx Local index of target octant.
     * \return Constant pointer to target octant.
     */
    const Octant*
    ParaTree::getOctant(uint32_t idx) const {
        return (m_octree.m_octants.data() + idx);
    };

    /** Get a ghost octant as pointer to the target octant.
     * NOTE: no checks will be performed on the ghost octant index.
     * \param[in] idx Local index (in ghosts structure) of target ghost octant.
     * \return Pointer to target ghost octant.
     */
    Octant*
    ParaTree::getGhostOctant(uint32_t idx) {
        return (m_octree.m_ghosts.data() + idx);
    };

    /** Get a ghost octant as constant pointer to the target octant.
     * NOTE: no checks will be performed on the ghost octant index.
     * \param[in] idx Local index (in ghosts structure) of target ghost octant.
     * \return Constant pointer to target ghost octant.
     */
    const Octant*
    ParaTree::getGhostOctant(uint32_t idx) const {
        return (m_octree.m_ghosts.data() + idx);
    };

    /*! Get the nature of an octant.
     * \param[in] oct Pointer to target octant.
     * \return Is octant ghost?
     */
    bool
    ParaTree::getIsGhost(const Octant* oct) const {
        return oct->getIsGhost();
    };

    /*! Get the layer number of the ghost halo an octant belong to.
     * \param[in] oct Pointer to target octant.
     * \return the layer in the ghost halo. For internal (non-ghost) octants
     * the function will return -1.
     */
    int
    ParaTree::getGhostLayer(const Octant* oct) const {
        return oct->getGhostLayer();
    };

    /*! Get a map of elements sent to the other processes during load balance
     * \return an unordered map associating rank to sent elements given by index extremes of a chunck
     */
    const ParaTree::LoadBalanceRanges &
    ParaTree::getLoadBalanceRanges() const {
        return m_loadBalanceRanges;
    };

    /*! Get the number of ghosts layers
     * \return the number of ghosts layers
     */
    std::size_t
    ParaTree::getNofGhostLayers() const {
        return m_nofGhostLayers;
    };

    /*! Set the number of ghosts layers
     * \param[in] nofGhostLayers The number of ghost layers in the ghost halo
     */
    void
    ParaTree::setNofGhostLayers(std::size_t nofGhostLayers) {
        if (m_nofGhostLayers == 0) {
            throw std::runtime_error ("It is not possible to disable the ghost halo!");
        }

        typedef std::result_of<decltype(&Octant::getGhostLayer)(Octant)>::type layer_t;
        typedef std::make_unsigned<layer_t>::type ulayer_t;
        ulayer_t maxNofGhostLayers = std::numeric_limits<layer_t>::max() + (ulayer_t) 1;
        if (nofGhostLayers > maxNofGhostLayers) {
            throw std::runtime_error ("Halo size exceeds the maximum allowed value.");
        }

        m_nofGhostLayers = nofGhostLayers;
    };

    /*! Get a map of border octants per process
     * \return A map of border octants per process
     */
    const std::map<int, u32vector> & ParaTree::getBordersPerProc() const {
        return m_bordersPerProc;
    }

    // =================================================================================== //
    // PRIVATE GET/SET METHODS
    // =================================================================================== //

    /*! Set the dimension
     */
    void
    ParaTree::setDim(uint8_t dim){
        m_dim = dim;
        if (m_dim != 0) {
            m_treeConstants = &(TreeConstants::instance(m_dim));
            m_periodic.resize(m_treeConstants->nFaces, false);
        } else {
            m_treeConstants = nullptr;
            m_periodic.clear();
        }
    }

    /*! Set the first finer descendant of the local tree.
     */
    void
    ParaTree::setFirstDescMorton(){
        m_octree.setFirstDescMorton();
    };

    /*! Set the last finer descendant of the local tree.
     */
    void
    ParaTree::setLastDescMorton(){
        m_octree.setLastDescMorton();
    };

    // =================================================================================== //
    // OTHER METHODS												    			   //
    // =================================================================================== //

    // =================================================================================== //
    // OTHER OCTANT BASED METHODS												    			   //
    // =================================================================================== //

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant.
     * \param[in] haveIidx Boolean flag to specify if the octant is passed with a valid idx
     * \param[in] idx Index of the searching octant. Its value is not important if haveIidx is false
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).*/
    void
    ParaTree::findNeighbours(const Octant* oct, bool haveIidx, uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{

        bool	Fedge = ((codim==2) && (m_dim==3));
        bool	Fnode = (codim == m_dim);

        if (codim == 1){
            m_octree.findNeighbours(oct, haveIidx, idx, iface, neighbours, isghost, onlyinternal);
        }
        else if (Fedge){
            m_octree.findEdgeNeighbours(oct, haveIidx, idx, iface, neighbours, isghost, onlyinternal);
        }
        else if (Fnode){
            m_octree.findNodeNeighbours(oct, haveIidx, idx, iface, neighbours, isghost, onlyinternal);
        }
        else {
            neighbours.clear();
            isghost.clear();
        }

    };


    /** Finds all the neighbours of a local octant through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
    void
    ParaTree::findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const {

        const Octant* oct = &m_octree.m_octants[idx];

        findNeighbours(oct, true, idx, iface, codim, neighbours, isghost, false);

    };

    /** Finds all the neighbours of an octant through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to current octant
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs. */
    void
    ParaTree::findNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const {

        findNeighbours(oct, false, 0, iface, codim, neighbours, isghost, false);

    };

    /** Finds the internal neighbours of ghost octant through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure ( only local octants ).
     * \param[in] idx Index of current octant
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     */
    void
    ParaTree::findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours) const {

        const Octant* oct = &m_octree.m_ghosts[idx];
        bvector isghost;
        findNeighbours(oct, true, idx, iface, codim, neighbours, isghost, true);

    };

    /** Finds all the neighbours of ghost octant through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure ( only local octants ).
     * \param[in] idx Index of current octant
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const {

        const Octant* oct = &m_octree.m_ghosts[idx];
        findNeighbours(oct, true, idx, iface, codim, neighbours, isghost, false);

    };

    /** Finds all the neighbours of ghost octant through iface face/edge/node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure ( only local octants ).
     * \param[in] oct Pointer to current ghost octant
     * \param[in] iface Index of face/edge/node passed through for neighbours finding
     * \param[in] codim Codimension of the iface-th entity 1=edge, 2=node
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findGhostNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const {

        uint32_t idx = getIdx(oct);
        findNeighbours(oct, true, idx, iface, codim, neighbours, isghost, false);

    };

    /** Finds all the neighbours of a node
    * \param[in] oct Pointer to current octant
    * \param[in] inode Index of node passed through for neighbours finding
    * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
    * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
    */
    void
    ParaTree::findAllNodeNeighbours(const Octant* oct, uint32_t inode, u32vector & neighbours, bvector & isghost) const {

        u32vector neigh_edge, neigh_face;
        bvector isghost_edge, isghost_face;

        // Get vertex neighbours
        int dim = getDim();
        m_octree.findNodeNeighbours(oct, false, 0, inode, neighbours, isghost, false);

        int octantLevel = getLevel(oct);
        if (dim == 3) {
            for (int edge : m_treeConstants->nodeEdge[inode]) {
                findNeighbours(oct, edge, 2, neigh_edge, isghost_edge);
                for (std::size_t i = 0; i < neigh_edge.size(); ++i) {
                    const Octant* neighOctant;
                    if (isghost_edge[i]==0) {
                        neighOctant = &m_octree.m_octants[neigh_edge[i]];
                    }
                    else {
                        neighOctant = &m_octree.m_ghosts[neigh_edge[i]];
                    }
                    int neighOctantLevel = getLevel(neighOctant);
                    if (neighOctantLevel <= octantLevel) {
                        neighbours.push_back(neigh_edge[i]);
                        isghost.push_back(isghost_edge[i]);
                    } else if (isNodeOnOctant(oct, inode, neighOctant)) {
                        neighbours.push_back(neigh_edge[i]);
                        isghost.push_back(isghost_edge[i]);
                    }
                }
            }
        }
        for (int j = 0; j < dim; ++j) {
            int face = m_treeConstants->nodeFace[inode][j];
            findNeighbours(oct, face, 1, neigh_face, isghost_face);
            for (std::size_t i = 0; i < neigh_face.size(); ++i) {
                const Octant* neighOctant;
                if (isghost_face[i]==0) {
                    neighOctant = &m_octree.m_octants[neigh_face[i]];
                }
                else {
                    neighOctant = &m_octree.m_ghosts[neigh_face[i]];
                }
                int neighOctantLevel = getLevel(neighOctant);
                if (neighOctantLevel <= octantLevel) {
                    neighbours.push_back(neigh_face[i]);
                    isghost.push_back(isghost_face[i]);
                } else if (isNodeOnOctant(oct, inode, neighOctant)) {
                    neighbours.push_back(neigh_face[i]);
                    isghost.push_back(isghost_face[i]);
                }
            }
        }
    };

    /** Finds all the neighbours of a node
    * \param[in] idx Index of current octant
    * \param[in] inode Index of node passed through for neighbours finding
    * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
    * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
    */
    void
    ParaTree::findAllNodeNeighbours(uint32_t idx, uint32_t inode, u32vector & neighbours, bvector & isghost) {

        Octant* oct = getOctant(idx);
        findAllNodeNeighbours(oct, inode, neighbours, isghost);
    };

    /** Compute neighbours adjacencies for every internal octant, storing them in a structure provided by the user
     * \param[in] idx Index of the octant
     * \param[out] globalNeighs Vector of global neighbours indices
     */
    void
    ParaTree::findAllGlobalNeighbours(uint32_t idx, std::vector<uint64_t> &globalNeighs){

        u32vector neighs;
        bvector isghost;
        findAllCodimensionNeighbours(idx,neighs,isghost);
        size_t nofNeighs = neighs.size();

        globalNeighs.resize(nofNeighs);
        for(size_t n = 0; n < nofNeighs; ++n){
            if(isghost[n])
                globalNeighs[n] = getGhostGlobalIdx(neighs[n]);
            else
                globalNeighs[n] = getGlobalIdx(neighs[n]);
        }
    }

    /** Finds all the neighbours of an internal octant through all its boundaries of any codimension.
     * Returns a vector with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findAllCodimensionNeighbours(uint32_t idx, u32vector & neighbours, bvector & isghost){
        Octant* oct = getOctant(idx);
        findAllCodimensionNeighbours(oct,neighbours,isghost);
    }

    /** Finds all the neighbours of an internal octant through all its boundaries of any codimension.
     * Returns a vector with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree. Neighbours are not sorted by Morton.
     * \param[in] oct pointer to the current octant
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findAllCodimensionNeighbours(Octant* oct, u32vector & neighbours, bvector & isghost){
        std::array<uint8_t, 4> codimensionsIndeces;
        codimensionsIndeces[0] = 0;
        codimensionsIndeces[1] = getNfaces();
        if(m_dim == 3){
            codimensionsIndeces[2] = getNedges();
        }
        codimensionsIndeces[m_dim] = getNnodes();

        neighbours.clear();
        neighbours.reserve(26);
        isghost.clear();
        isghost.reserve(26);
        u32vector singleCodimNeighbours;
        bvector singleCodimIsGhost;
        for(uint8_t codim = 1; codim <= m_dim; ++codim){
            for(int icodim = 0; icodim < codimensionsIndeces[codim]; ++icodim){
                findNeighbours(oct,icodim,codim,singleCodimNeighbours,singleCodimIsGhost);
                for(size_t i = 0; i < singleCodimIsGhost.size(); ++i){
                    isghost.push_back(singleCodimIsGhost[i]);
                    neighbours.push_back(singleCodimNeighbours[i]);
                }
            }
        }
    }

    /** Finds all the neighbours of a ghost octant through all its boundaries of any codimension.
     * Returns a vector with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findGhostAllCodimensionNeighbours(uint32_t idx, u32vector & neighbours, bvector & isghost){
        Octant* oct = getGhostOctant(idx);
        findGhostAllCodimensionNeighbours(oct,neighbours,isghost);
    }

    /** Finds all the neighbours of a ghost octant through all its boundaries of any codimension.
     * Returns a vector with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree. Neighbours are not sorted by Morton.
     * \param[in] oct pointer to the current octant
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs.
     */
    void
    ParaTree::findGhostAllCodimensionNeighbours(Octant* oct, u32vector & neighbours, bvector & isghost){
        std::array<uint8_t, 4> codimensionsIndeces;
        codimensionsIndeces[0] = 0;
        codimensionsIndeces[1] = getNfaces();
        if(m_dim == 3){
            codimensionsIndeces[2] = getNedges();
        }
        codimensionsIndeces[m_dim] = getNnodes();

        neighbours.clear();
        neighbours.reserve(26);
        isghost.clear();
        isghost.reserve(26);
        u32vector singleCodimNeighbours;
        bvector singleCodimIsGhost;
        for(uint8_t codim = 1; codim <= m_dim; ++codim){
            for(int icodim = 0; icodim < codimensionsIndeces[codim]; ++icodim){
                findGhostNeighbours(oct,icodim,codim,singleCodimNeighbours,singleCodimIsGhost);
                for(size_t i = 0; i < singleCodimIsGhost.size(); ++i){
                    isghost.push_back(singleCodimIsGhost[i]);
                    neighbours.push_back(singleCodimNeighbours[i]);
                }
            }
        }
    }

    /** Get the internal octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Pointer to octant owner of target point (=NULL if point is outside of the domain).
     */
    Octant*
    ParaTree::getPointOwner(const dvector &point){
        uint32_t idx = getPointOwnerIdx(point);
        if(idx < numeric_limits<uint32_t>::max())
            return &m_octree.m_octants[idx];
        else
            return NULL;
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Pointer to octant owner of target point (=NULL if point is outside of the ghosted domain).
     */
    Octant*
    ParaTree::getPointOwner(const dvector &point, bool & isghost){
        uint32_t idx = getPointOwnerIdx(point,isghost);
        if(idx < numeric_limits<uint32_t>::max())
            if(isghost)
                return &m_octree.m_ghosts[idx];
            else
                return &m_octree.m_octants[idx];
        else
            return NULL;
    };


    /** Get the internal octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Pointer to octant owner of target point (=NULL if point is outside of the domain).
     */
    Octant*
    ParaTree::getPointOwner(const darray3 &point){
        uint32_t idx = getPointOwnerIdx(point);
        if(idx < numeric_limits<uint32_t>::max())
            return &m_octree.m_octants[idx];
        else
            return NULL;
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Pointer to octant owner of target point (=NULL if point is outside of the ghostd domain).
     */
    Octant*
    ParaTree::getPointOwner(const darray3 &point, bool & isghost){
        uint32_t idx = getPointOwnerIdx(point,isghost);
        if(idx < numeric_limits<uint32_t>::max())
            if(isghost)
                return &m_octree.m_ghosts[idx];
            else
                return &m_octree.m_octants[idx];
        else
            return NULL;
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const darray3 &point) const {
        return getPointOwnerIdx(point.data());
    }

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the ghosted domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const darray3 &point, bool & isghost) const {
        return getPointOwnerIdx(point.data(),isghost);
    }

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const dvector &point) const {
        assert(point.size() >= 3);
        return getPointOwnerIdx(point.data());
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the ghosted domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const dvector &point, bool & isghost) const {
        assert(point.size() >= 3);
        return getPointOwnerIdx(point.data(),isghost);
    };


    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const double * point) const {
        uint32_t noctants = m_octree.m_octants.size();
        if(noctants==0)
            return numeric_limits<uint32_t>::max();
        uint32_t idxtry = noctants/2;
        uint32_t x, y, z;
        uint64_t morton, mortontry;
        int powner = 0;
        //ParaTree works in [0,1] domain
        if (point[0] > 1+m_tol || point[1] > 1+m_tol || point[2] > 1+m_tol
            || point[0] < -m_tol || point[1] < -m_tol || point[2] < -m_tol){
            return numeric_limits<uint32_t>::max();
        }

        x = m_trans.mapX(std::min(std::max(point[0], 0.0), 1.0));
        y = m_trans.mapY(std::min(std::max(point[1], 0.0), 1.0));
        z = m_trans.mapZ(std::min(std::max(point[2], 0.0), 1.0));

        if (x == TreeConstants::MAX_LENGTH) x = x - 1;
        if (y == TreeConstants::MAX_LENGTH) y = y - 1;
        if (z == TreeConstants::MAX_LENGTH) z = z - 1;
        morton = PABLO::computeMorton(x,y,z);


        powner = 0;
        if(!m_serial) powner = findOwner(morton);

        if ((powner!=m_rank) && (!m_serial))
            return numeric_limits<uint32_t>::max();

        int32_t jump = idxtry;
        while(abs(jump) > 0){

            mortontry = m_octree.m_octants[idxtry].computeMorton();
            jump = ((mortontry<morton)-(mortontry>morton))*abs(jump)/2;
            idxtry += jump;
            if (idxtry > noctants-1){
                if (jump > 0){
                    idxtry = noctants - 1;
                    jump = 0;
                }
                else if (jump < 0){
                    idxtry = 0;
                    jump = 0;
                }
            }
        }
        if(m_octree.m_octants[idxtry].computeMorton() == morton){
            return idxtry;
        }
        else{
            // Step until the mortontry lower than morton (one idx of distance)
            {
                while(m_octree.m_octants[idxtry].computeMorton() < morton){
                    idxtry++;
                    if(idxtry > noctants-1){
                        idxtry = noctants-1;
                        break;
                    }
                }
                while(m_octree.m_octants[idxtry].computeMorton() > morton){
                    idxtry--;
                    if(idxtry > noctants-1){
                        idxtry = 0;
                        break;
                    }
                }
            }
            return idxtry;
        }
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the ghosted domain).
     */
    uint32_t
    ParaTree::getPointOwnerIdx(const double * point, bool & isghost) const {
        uint32_t noctants = m_octree.m_octants.size();
        if(noctants==0)
            return numeric_limits<uint32_t>::max();
        uint32_t idxtry = noctants/2;
        uint32_t x, y, z;
        uint64_t morton, mortontry;
        int powner = 0;
        isghost = false;
        //ParaTree works in [0,1] domain
        if (point[0] > 1+m_tol || point[1] > 1+m_tol || point[2] > 1+m_tol
            || point[0] < -m_tol || point[1] < -m_tol || point[2] < -m_tol){
            return numeric_limits<uint32_t>::max();
        }

        x = m_trans.mapX(std::min(std::max(point[0], 0.0), 1.0));
        y = m_trans.mapY(std::min(std::max(point[1], 0.0), 1.0));
        z = m_trans.mapZ(std::min(std::max(point[2], 0.0), 1.0));

        if (x == TreeConstants::MAX_LENGTH) x = x - 1;
        if (y == TreeConstants::MAX_LENGTH) y = y - 1;
        if (z == TreeConstants::MAX_LENGTH) z = z - 1;
        morton = PABLO::computeMorton(x,y,z);


        powner = 0;
        if(!m_serial) powner = findOwner(morton);

        if (powner==m_rank){

            int32_t jump = idxtry;
            while(abs(jump) > 0){
                
                mortontry = m_octree.m_octants[idxtry].computeMorton();
                jump = ((mortontry<morton)-(mortontry>morton))*abs(jump)/2;
                idxtry += jump;
                if (idxtry > noctants-1){
                    if (jump > 0){
                        idxtry = noctants - 1;
                        jump = 0;
                    }
                    else if (jump < 0){
                        idxtry = 0;
                        jump = 0;
                    }
                }
            }
            if(m_octree.m_octants[idxtry].computeMorton() == morton){
                return idxtry;
            }
            else{
                // Step until the mortontry lower than morton (one idx of distance)
                {
                    while(m_octree.m_octants[idxtry].computeMorton() < morton){
                        idxtry++;
                        if(idxtry > noctants-1){
                            idxtry = noctants-1;
                            break;
                        }
                    }
                    while(m_octree.m_octants[idxtry].computeMorton() > morton){
                        idxtry--;
                        if(idxtry > noctants-1){
                            idxtry = 0;
                            break;
                        }
                    }
                }
                return idxtry;
            }
        }
        else if((powner != m_rank) && m_serial){
            return numeric_limits<uint32_t>::max();
        }
        else{
            //GHOST SEARCH
            uint32_t nghosts = m_octree.m_ghosts.size();
            idxtry = nghosts/2;
            int32_t jump = idxtry;
            while(abs(jump) > 0){
                
                mortontry = m_octree.m_ghosts[idxtry].computeMorton();
                jump = ((mortontry<morton)-(mortontry>morton))*abs(jump)/2;
                idxtry += jump;
                if (idxtry > nghosts-1){
                    if (jump > 0){
                        idxtry = nghosts - 1;
                        jump = 0;
                    }
                    else if (jump < 0){
                        idxtry = 0;
                        jump = 0;
                    }
                }
            }
            if(m_octree.m_ghosts[idxtry].computeMorton() == morton){
                isghost = true;
                return idxtry;
            }
            else{
                // Step until the mortontry lower than morton (one idx of distance)
                {
                    while(m_octree.m_ghosts[idxtry].computeMorton() < morton){
                        idxtry++;
                        if(idxtry > nghosts-1){
                            idxtry = nghosts-1;
                            break;
                        }
                    }
                    while(m_octree.m_ghosts[idxtry].computeMorton() > morton){
                        idxtry--;
                        if(idxtry > nghosts-1){
                            idxtry = 0;
                            break;
                        }
                    }
                }

                const Octant* octtry = getGhostOctant(idxtry);
                dvector anchor_idxtry = {{getX(octtry),getY(octtry),getZ(octtry)}};
                double size_try = getSize(octtry);
                bool isInIdxtry = true;

                for(int i = 0; i < m_dim; ++i){
                    isInIdxtry *= (point[i] >= anchor_idxtry[i] && point[i] <= anchor_idxtry[i] + size_try);
                }

                if( isInIdxtry){
                    isghost = true;
                    return idxtry;
                }
                else{
                    return numeric_limits<uint32_t>::max();
                }
            }
        }///end ghosts search
    };

    /** Get the octant owner rank of an input point.
     * \param[in] point Coordinates of target point.
     * \return Owner rank of target point (negative if out of global domain).
     */
    int
    ParaTree::getPointOwnerRank(darray3 point){

        uint32_t x,y,z;
        uint64_t morton;

        if (point[0] > 1+m_tol || point[1] > 1+m_tol || point[2] > 1+m_tol
            || point[0] < -m_tol || point[1] < -m_tol || point[2] < -m_tol){
            return -1;
        }
        point[0] = min(max(point[0],0.0),1.0);
        point[1] = min(max(point[1],0.0),1.0);
        point[2] = min(max(point[2],0.0),1.0);

        x = m_trans.mapX(point[0]);
        y = m_trans.mapY(point[1]);
        z = m_trans.mapZ(point[2]);

        if ((x > TreeConstants::MAX_LENGTH) || (y > TreeConstants::MAX_LENGTH) || (z > TreeConstants::MAX_LENGTH)
            || (point[0] < m_trans.m_origin[0]) || (point[1] < m_trans.m_origin[1]) || (point[2] < m_trans.m_origin[2])){
            return -1;
        }

        if (m_serial)
            return m_rank;

        if (x == TreeConstants::MAX_LENGTH)
            x = x - 1;

        if (y == TreeConstants::MAX_LENGTH)
            y = y - 1;

        if (z == TreeConstants::MAX_LENGTH)
            z = z - 1;

        morton = PABLO::computeMorton(x, y, z);

        for (int p = 0; p < m_nproc; ++p){
            if (morton <= m_partitionLastDesc[p] && morton >= m_partitionFirstDesc[p])
                return p;
        }

        return -1;
    };

    /** Get mapping info of an octant after an adapting with tracking changes.
     * \param[in] idx Index of new octant.
     * \param[out] mapper Mapper from new octants to old octants. I.e. mapper[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
     * if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
     * if the i-th octant is new after coarsening the j-th old octant was a child of the new octant (mapper size = 4).
     * \param[out] isghost Info on ghostness of old octants.
     * I.e. isghost[i] = true/false -> the mapper[i] = j-th old octant was a local/ghost octant.
     */
    void
    ParaTree::getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost) const {

        if (idx >= m_mapIdx.size()){
            throw std::runtime_error ("Invalid value for input index in getMapping");
        }

        // Coarsening has to be handled separately, all other changes can just
        // return the value stored in the mapper.
        if (getIsNewC(idx)){
            // Count the children
            int nChildren = m_treeConstants->nChildren;

            int nInternalChildren = nChildren;
            if (idx == (getNumOctants() - 1)) {
                nInternalChildren -= m_octree.m_lastGhostBros.size();
            }

            // Fill the mapper
            mapper.resize(nChildren);
            isghost.resize(nChildren);

            for (int i = 0; i < nInternalChildren; i++){
                mapper[i]  = m_mapIdx[idx] + i;
                isghost[i] = false;
            }

            for (int i = nInternalChildren; i < nChildren; i++){
                mapper[i]  = m_octree.m_lastGhostBros[i - nInternalChildren];
                isghost[i] = true;
            }
        } else {
            mapper.resize(1);
            isghost.resize(1);

            mapper[0]  = m_mapIdx[idx];
            isghost[0] = false;
        }

    };

    /** Get mapping info of an octant after an adapting or loadbalance with tracking changes.
     * \param[in] idx Index of new octant.
     * \param[out] mapper Mapper from new octants to old octants. I.e. mapper[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
     * if the i-th octant is new after refinement or loadbalance the j-th old octant was the father of the new octant or the same octant respectively;
     * if the i-th octant is new after coarsening the j-th old octant was a child of the new octant (mapper size = 4).
     * \param[out] isghost Info on ghostness of old octants.
     * \param[out] rank Process where the octant was located before the adapt/loadbalance.
     * I.e. isghost[i] = true/false -> the mapper[i] = j-th old octant was a local/ghost octant.
     */
    void
    ParaTree::getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost, ivector & rank) const {

        if (m_lastOp == OP_ADAPT_MAPPED){
            getMapping(idx, mapper, isghost);
            int n = isghost.size();
            rank.resize(n);
            for (int i=0; i<n; i++){
                rank[i] = m_rank;
            }
        }
        else if (m_lastOp == OP_LOADBALANCE || m_lastOp == OP_LOADBALANCE_FIRST){
            mapper.resize(1);
            isghost.resize(1);
            rank.resize(1);
            uint64_t gidx = getGlobalIdx(idx);
            mapper[0] = gidx;
            for (int iproc=0; iproc<m_nproc; ++iproc){
                if (m_partitionRangeGlobalIdx0[iproc]>=gidx){
                    if (iproc > 0)
                        mapper[0] -= m_partitionRangeGlobalIdx0[iproc-1] + 1;
                    rank[0] = (m_lastOp == OP_LOADBALANCE_FIRST ? m_rank : iproc);
                    isghost[0] = false;
                    break;
                }
            }
        }
    };

    /** Get octants with marker different from zero and the related markers and ghostness info.
        * The methods has to be called after a apredapt, otherwise it calls preadapt method.
        * \param[out] idx Vector of local indices of octants with marker different from zero.
        * \param[out] markers Vector with markers related to octants in the idx list.
        * \param[out] isghost Info on ghostness of octants.
        */
    void
    ParaTree::getPreMapping(u32vector & idx, vector<int8_t> & markers, vector<bool> & isghost)
    {
        if (m_lastOp != OP_PRE_ADAPT) {
            throw std::runtime_error("Last operation different from preadapt, unable to call getPreMarker function");
        }

        idx.clear();
        markers.clear();
        isghost.clear();

        std::size_t firstGsize = m_octree.m_firstGhostBros.size();
        std::size_t lastGsize = m_octree.m_lastGhostBros.size();

        idx.reserve(getNumOctants() + firstGsize + lastGsize);
        markers.reserve(getNumOctants() + firstGsize + lastGsize);
        isghost.reserve(getNumOctants() + firstGsize + lastGsize);

        //insert first ghost brothers if present
        {
            int8_t marker;
            for (uint32_t id : m_octree.m_firstGhostBros){
                marker = m_octree.m_ghosts[id].getMarker();
                idx.push_back(id);
                markers.push_back(marker);
                isghost.push_back(true);
            }
        }

        {
            int8_t marker;
            for (uint32_t count = 0; count < m_octree.getNumOctants(); count++){
                marker = m_octree.m_octants[count].getMarker();
                if (marker != 0){
                    idx.push_back(count);
                    markers.push_back(marker);
                    isghost.push_back(false);
                }
            }
        }

        //insert last ghost brothers if present
        {
            int8_t marker;
            for (uint32_t id : m_octree.m_lastGhostBros){
                marker = m_octree.m_ghosts[id].getMarker();
                idx.push_back(id);
                markers.push_back(marker);
                isghost.push_back(true);
            }
        }
    };

    /** Check if a node lies on the specified octant.
     * \param[in] nodeOctant Pointer to the octant owning the node
     * \param[in] nodeIndex Local index of the node
     * \param[in] octant Pointer to the octant for which the check has to be
     * berformed
     */
    bool
    ParaTree::isNodeOnOctant(const Octant *nodeOctant, uint8_t nodeIndex, const Octant *octant) const {

        int dim = octant->getDim();

        // Get the coordinates of the node
        std::array<uint32_t, 3> nodeCoords = nodeOctant->getNode(nodeIndex);

        // Get minimum/maximum coordinates of the contant
        std::array<uint32_t, 3> minOctantCoords = octant->getNode(0);
        std::array<uint32_t, 3> maxOctantCoords = octant->getNode(3 + 4 * (dim - 2));

        // Check if the node intersects the bounding box octant
        //
        // NOTE: since the octants are cubes, the bounding box coincides with
        //       the octant.
        for (int i = 0; i < dim; ++i) {
            // I-th coordinate of the node
            uint32_t nodeCoord = nodeCoords[i];

            // Minimum/maximum i-th coordinate of the octant
            uint32_t minOctantBBCoord = minOctantCoords[i];
            uint32_t maxOctantBBCoord = maxOctantCoords[i];

            // Check if the node intersects the octant
            if (nodeCoord < minOctantBBCoord) {
                return false;
            } else if (nodeCoord > maxOctantBBCoord) {
                return false;
            }
        }

        return true;
    }

    /** Check if an edge lies on the specified octant.
     * \param[in] edgeOctant Pointer to the octant owning the edge
     * \param[in] edgeIndex Local index of the edge
     * \param[in] octant Pointer to the octant for which the check has to be
     * berformed
     */
    bool
    ParaTree::isEdgeOnOctant(const Octant* edgeOctant, uint8_t edgeIndex, const Octant* octant) const {

        // Edges are only defined on three-dimensional trees.
        int dim = octant->getDim();
        assert(dim == 3);

        // Get the coordinates of the edge
        const uint8_t (*edgeNodes)[2] = &m_treeConstants->edgeNode[edgeIndex];
        std::array<uint32_t, 3> minEdgeCoords = edgeOctant->getNode((*edgeNodes)[0]);
        std::array<uint32_t, 3> maxEdgeCoords = edgeOctant->getNode((*edgeNodes)[1]);

        // Get minimum/maximum coordinates of the contant
        std::array<uint32_t, 3> minOctantCoords = octant->getNode(0);
        std::array<uint32_t, 3> maxOctantCoords = octant->getNode(7);

        // Check if the edge intersects the bounding box octant
        //
        // NOTE: since the octants are cubes, the bounding box coincides with
        //       the octant.
        for (int i = 0; i < dim; ++i) {
            // Minimum/maximum i-th coordinate of the edge
            uint32_t minEdgeCoord = minEdgeCoords[i];
            uint32_t maxEdgeCoord = maxEdgeCoords[i];

            // Minimum/maximum i-th coordinate of the octant
            uint32_t minOctantBBCoord = minOctantCoords[i];
            uint32_t maxOctantBBCoord = maxOctantCoords[i];

            // Check if the edge intersects the octant
            if (minEdgeCoord < minOctantBBCoord && maxEdgeCoord < minOctantBBCoord) {
                return false;
            } else if (minEdgeCoord > maxOctantBBCoord && maxEdgeCoord > maxOctantBBCoord) {
                return false;
            }
        }

        return true;
    }

    /** Check if a face lies on the specified octant.
     * \param[in] faceOctant Pointer to the octant owning the face
     * \param[in] faceIndex Local index of the face
     * \param[in] octant Pointer to the octant for which the check has to be
     * berformed
     */
    bool
    ParaTree::isFaceOnOctant(const Octant* faceOctant, uint8_t faceIndex, const Octant* octant) const {

        int dim = octant->getDim();

        // Get minimum/maximum coordinates of the face
        const uint8_t (*faceNodes)[6][4] = &m_treeConstants->faceNode;
        std::array<uint32_t, 3> minFaceCoords = faceOctant->getNode((*faceNodes)[faceIndex][0]);
        std::array<uint32_t, 3> maxFaceCoords = faceOctant->getNode((*faceNodes)[faceIndex][2 * dim - 1]);

        // Get minimum/maximum coordinates of the contant
        std::array<uint32_t, 3> minOctantCoords = octant->getNode(0);
        std::array<uint32_t, 3> maxOctantCoords = octant->getNode(3 + 4 * (dim - 2));

        // Check if the face intersects the bounding box octant
        for (int i = 0; i < dim; ++i) {
            // Minimum/maximum i-th coordinate of the face
            uint32_t minFaceCoord = minFaceCoords[i];
            uint32_t maxFaceCoord = maxFaceCoords[i];

            // Minimum/maximum i-th coordinate of the octant
            uint32_t minOctantBBCoord = minOctantCoords[i];
            uint32_t maxOctantBBCoord = maxOctantCoords[i];

            // Check if the face intersects the octant
            if (minFaceCoord < minOctantBBCoord && maxFaceCoord < minOctantBBCoord) {
                return false;
            } else if (minFaceCoord > maxOctantBBCoord && maxFaceCoord > maxOctantBBCoord) {
                return false;
            }
        }

        return true;
    }

    // =================================================================================== //
    // OTHER PARATREE BASED METHODS												    			   //
    // =================================================================================== //


    /** Rearrange the octree markers with user setup for markers and 2:1 balancing conditions.
     *
     *  This function is analogous to pre-adapt method, but no function is mandatory after
     *  a settleMarkers call.
     *
     *  Note: the last operation tag is not changed after a settleMarkers call.
     */
    void
    ParaTree::settleMarkers(){

        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " SETTLE MARKERS " << endl;

        balance21(true, false);

        (*m_log) << " " << endl;
        (*m_log) << "---------------------------------------------" << endl;

    };

    /** Pre-adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
     *
     *  The user can call pre-adapt and then adapt or only adapt, however after the
     *  pre-adapt function has been called the adapt function is mandatory.
     *
     */
    void
    ParaTree::preadapt(){

        balance21(true, false);

        m_lastOp = OP_PRE_ADAPT;

        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " PRE-ADAPT " << endl;

        (*m_log) << " " << endl;
        (*m_log) << "---------------------------------------------" << endl;

    };

    /** Check to control if the tree has to be adapted.
     * \return Boolean true if at least one octant has marker not zero.
     */
    bool
    ParaTree::checkToAdapt(){

        bool lcheck = false;
        bool gcheck = false;
        octvector::iterator it = m_octree.m_octants.begin();
        octvector::iterator itend = m_octree.m_octants.end();
        while(!lcheck && it != itend){
            lcheck = (it->getMarker() != 0);
            ++it;
        }
        if (m_nproc == 1){
            gcheck = lcheck;
        }
        else{
#if BITPIT_ENABLE_MPI==1
        m_errorFlag = MPI_Allreduce(&lcheck,&gcheck,1,MPI_C_BOOL,MPI_LOR,m_comm);
#endif
        }
        return gcheck;
    };

    /** Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
     * \param[in] mapper_flag True/False if you want/don't want to track the changes in structure octant by a mapper.
     * \n NOTE: if mapper_flag = true the adapt method ends after a single level (refining/coarsening) adaptation.
     * \n The resulting markers will be increased/decreased by one.
     * \return Boolean true if adapt has done something.
     * \n NOTE: if  pre-adapt method is called before an adapt call the adapt method
     * do not perform pre-adapt process (no 2:1 balance check).
     */
    bool
    ParaTree::adapt(bool mapper_flag){

        bool done = false;

        done = private_adapt_mapidx(mapper_flag);
        m_status += done;
        return done;

    };

    /** Adapt the octree mesh refining all the octants by one level.
     * Optionally track the changes in structure octant by a mapper.
     * \param[in] mapper_flag True/false for tracking/not tracking the changes in structure octant .
     */
    bool
    ParaTree::adaptGlobalRefine(bool mapper_flag) {
        //TODO recoding for adapting with abs(marker) > 1
        uint32_t nocts0 = getNumOctants();
        vector<Octant>::iterator iter, iterend = m_octree.m_octants.end();

        for (iter = m_octree.m_octants.begin(); iter != iterend; ++iter){
            iter->m_info[Octant::INFO_NEW4REFINEMENT] = false;
            iter->m_info[Octant::INFO_NEW4COARSENING] = false;
            iter->m_info[Octant::INFO_AUX] = false;
        }

        // Initialize mapping
        u32vector(nocts0).swap(m_mapIdx);
        for (uint32_t i=0; i<nocts0; i++){
            m_mapIdx[i] = i;
        }

        // Update tree
        bool globalDone = false;
#if BITPIT_ENABLE_MPI==1
        if(m_serial){
#endif
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Global Refine)" << endl;
            (*m_log) << " " << endl;

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(getNumOctants())) << endl;

            // Refine
            while(m_octree.globalRefine(m_mapIdx));

            if (getNumOctants() > nocts0)
                globalDone = true;
            (*m_log) << " Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(getNumOctants())) << endl;
            nocts0 = getNumOctants();
            updateAdapt();

            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
#if BITPIT_ENABLE_MPI==1
        }
        else{
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Global Refine)" << endl;
            (*m_log) << " " << endl;

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;

            // Refine
            while(m_octree.globalRefine(m_mapIdx));

            bool localDone = false;
            if (getNumOctants() > nocts0)
                localDone = true;
            updateAdapt();
            computeGhostHalo();
            (*m_log) << " Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;
            nocts0 = getNumOctants();

            m_errorFlag = MPI_Allreduce(&localDone,&globalDone,1,MPI_C_BOOL,MPI_LOR,m_comm);
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }
#endif

        // Update last operation
        if (mapper_flag) {
            m_lastOp = OP_ADAPT_MAPPED;
        }
        else{
            m_lastOp = OP_ADAPT_UNMAPPED;
        }

        return globalDone;
    }

    /** Adapt the octree mesh coarsening all the octants by one level.
     * Optionally track the changes in structure octant by a mapper.
     * \param[in] mapper_flag True/false for tracking/not tracking the changes in structure octant .
     */
    bool
    ParaTree::adaptGlobalCoarse(bool mapper_flag) {
        //TODO recoding for adapting with abs(marker) > 1
        uint32_t nocts0 = getNumOctants();
        vector<Octant>::iterator iter, iterend = m_octree.m_octants.end();

        for (iter = m_octree.m_octants.begin(); iter != iterend; ++iter){
            iter->m_info[Octant::INFO_NEW4REFINEMENT] = false;
            iter->m_info[Octant::INFO_NEW4COARSENING] = false;
            iter->m_info[Octant::INFO_AUX] = false;
        }

        m_mapIdx.clear();
        if (mapper_flag){
            // m_mapIdx init
            m_mapIdx.resize(nocts0);
            u32vector(m_mapIdx).swap(m_mapIdx);

            for (uint32_t i=0; i<nocts0; i++){
                m_mapIdx[i] = i;
            }
        }

        bool globalDone = false;
#if BITPIT_ENABLE_MPI==1
        if(m_serial){
#endif
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Global Coarse)" << endl;
            (*m_log) << " " << endl;

            // 2:1 Balance
            balance21(true, false);

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(getNumOctants())) << endl;

            // Coarse
            while(m_octree.globalCoarse(m_mapIdx));
            updateAfterCoarse();
            balance21(false, true);
            while(m_octree.refine(m_mapIdx));
            updateAdapt();

            if (getNumOctants() < nocts0){
                globalDone = true;
            }
            nocts0 = getNumOctants();

            (*m_log) << " Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(nocts0)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
#if BITPIT_ENABLE_MPI==1
        }
        else{
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Global Coarse)" << endl;
            (*m_log) << " " << endl;

            // 2:1 Balance
            balance21(true, false);

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;

            // Coarse
            while(m_octree.globalCoarse(m_mapIdx));
            updateAfterCoarse();
            computeGhostHalo();
            balance21(false, true);
            while(m_octree.refine(m_mapIdx));
            updateAdapt();

            computeGhostHalo();
            bool localDone = false;
            if (getNumOctants() < nocts0){
                localDone = true;
            }
            nocts0 = getNumOctants();

            m_errorFlag = MPI_Allreduce(&localDone,&globalDone,1,MPI_C_BOOL,MPI_LOR,m_comm);
            (*m_log) << " Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }
#endif
        return globalDone;
    }

    /*! Get the local current maximum size of the octree.
     * \return Local current maximum size of the local partition of the octree.
     */
    uint8_t
    ParaTree::getMaxDepth() const{
        return m_maxDepth;
    };

    /** It finds the process owning the element definded by the Morton number passed as argument
     * The Morton number can be computed using the method computeMorton() of Octant.
     * \param[in] morton Morton number of the element you want find the owner of
     * \return Rank of the process owning the element
     */
    int
    ParaTree::findOwner(const uint64_t & morton) const {
        // Early return if the requested morton is on first partition
        if (morton <= m_partitionLastDesc[0]) {
            return 0;
        }
        if (morton > m_partitionLastDesc[m_nproc-1]) {
            return -1;
        }

        // Find the partition using a bisection method
        int p = -1;
        int length = m_nproc;
        int beg = 0;
        int end = m_nproc -1;
        int seed = m_nproc/2;
        while(beg != end){
            if(morton <= m_partitionLastDesc[seed]){
                end = seed;
                if(morton > m_partitionLastDesc[seed-1])
                    beg = seed;
            }
            else{
                beg = seed;
                if(morton <= m_partitionLastDesc[seed+1])
                    beg = seed + 1;
            }
            length = end - beg;
            seed = beg + length/2;
        }
        if(beg!=0){
            while(m_partitionLastDesc[beg] == m_partitionLastDesc[beg-1]){
                --beg;
                if(beg==0)
                    break;
            }
        }
        p = beg;
        return p;
    }

    /** It finds the process owning the element definded by the global index passed as argument
     * The global index can be computed using the methods getGlobalIdx or getGhostGlobalIdx.
     * \param[in] globalIndex index of the element you want find the owner of
     * \return Rank of the process owning the element
     */
    int
    ParaTree::getOwnerRank(const uint64_t & globalIndex) const {
        // Get the iterator point to the onwer rank
        std::vector<uint64_t>::const_iterator rankItr = std::lower_bound (m_partitionRangeGlobalIdx.begin(), m_partitionRangeGlobalIdx.end(), globalIndex);

        // Get the onwer rank
        int ownerRank;
        if (rankItr == m_partitionRangeGlobalIdx.end()) {
            ownerRank = -1;
        } else {
            ownerRank = std::distance(m_partitionRangeGlobalIdx.begin(), rankItr);
        }

        return ownerRank;
    }

    /** Compute the connectivity of octants and store the coordinates of nodes.
     */
    void
    ParaTree::computeConnectivity() {
        m_octree.computeConnectivity();
    }

    /** Clear the connectivity of octants.
     */
    void
    ParaTree::clearConnectivity() {
        m_octree.clearConnectivity();
    }

    /** Update the connectivity of octants.
     */
    void
    ParaTree::updateConnectivity() {
        m_octree.updateConnectivity();
    }

    /** Get the connectivity of the octants
     * \return Constant reference to the connectivity matrix of noctants*nnodes with
     * the connectivity of each octant (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector2D &
    ParaTree::getConnectivity() const {
        return m_octree.m_connectivity;
    }

    /** Get the local connectivity of an octant
     * \param[in] idx Local index of octant
     * \return Constant reference to the connectivity of the octant
     * (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector &
    ParaTree::getConnectivity(uint32_t idx) const {
        return m_octree.m_connectivity[idx];
    }

    /** Get the local connectivity of an octant
     * \param[in] oct Pointer to an octant
     * \return Constant reference to the connectivity of the octant (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector &
    ParaTree::getConnectivity(Octant* oct) const {
        return m_octree.m_connectivity[getIdx(oct)];
    }

    /** Get the logical coordinates of the nodes
     * \return Constant reference to the nodes matrix [nnodes*3] with the coordinates of the nodes.
     */
    const u32arr3vector &
    ParaTree::getNodes() const {
        return m_octree.m_nodes;
    }

    /** Get the logical coordinates of a node
     * \param[in] inode Local index of node
     * \return Constant reference to a vector containing the coordinates of the node.
     */
    const u32array3 &
    ParaTree::getNodeLogicalCoordinates(uint32_t inode) const {
        return m_octree.m_nodes[inode];
    }

    /** Get the physical coordinates of a node
     * \param[in] inode Local index of node
     * \return Vector with the coordinates of the node.
     */
    darray3
    ParaTree::getNodeCoordinates(uint32_t inode) const {
        return m_trans.mapCoordinates(m_octree.m_nodes[inode]);
    }

    /** Get the connectivity of the ghost octants
     * \return Constant reference to connectivity matrix [nghostoctants*nnodes] with
     * the connectivity of each octant (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector2D &
    ParaTree::getGhostConnectivity() const {
        return m_octree.m_ghostsConnectivity;
    }

    /** Get the local connectivity of a ghost octant
     * \param[in] idx Local index of ghost octant
     * \return Constant reference to the connectivity of the ghost octant
     * (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector &
    ParaTree::getGhostConnectivity(uint32_t idx) const {
        return m_octree.m_ghostsConnectivity[idx];
    }

    /** Get the local connectivity of a ghost octant
     * \param[in] oct Pointer to a ghost octant
     * \return Constant reference to the connectivity of the ghost octant
     * (4/8 indices of nodes for 2D/3D case).
     */
    const u32vector &
    ParaTree::getGhostConnectivity(const Octant* oct) const {
        return m_octree.m_ghostsConnectivity[getIdx(oct)];
    }



    /** Check if the grid is 2:1 balanced across intersection of balanceCodim codimension
     * \return true if the grid is 2:1 balanced
     */
    bool
    ParaTree::check21Balance(){
        
        bool balanced=true;
        uint32_t nocts = getNumOctants();
        uint8_t balanceCodim = getBalanceCodimension();
        uint8_t maxCodim = balanceCodim;
        u32vector neighs;
        bvector isghost;
        uint8_t levelDiff;
        for(uint8_t c = 1; c <= maxCodim; ++c){
            uint8_t nCodimSubelements = c==1 ? getNfaces() : (c==m_dim ? getNnodes() : getNedges());
            for(uint32_t i = 0; i < nocts; ++i){
                uint8_t level = getLevel(i);
                for(uint8_t f = 0; f < nCodimSubelements; ++f){
                    findNeighbours(i,f,c,neighs,isghost);
                    size_t nsize = neighs.size();
                    for(size_t n = 0; n < nsize; ++n){
                        const Octant* noct = isghost[n] ? getGhostOctant(neighs[n]) : getOctant(neighs[n]);
                        uint8_t nlevel = noct->getLevel();
                        levelDiff = uint8_t(abs(int(nlevel)-int(level)));
                        if(levelDiff > 1){
                            log::Visibility visi = m_log->getVisibility();
                            m_log->setVisibility(log::GLOBAL);
                            (*m_log) << "---------------------------------------------" << std::endl;
                            (*m_log) << "LOCALLY 2:1 UNBALANCED OCTREE" << std::endl;
                            (*m_log) << "I'm " << getRank() << ": element " << i << " is not 2:1 balanced across " << int(f) << " subentity of codim " << int(c) << ", relative to " << (isghost[n] ? "ghost" : "internal") << "neighbour " << neighs[n] << std::endl;
                            (*m_log) << "---------------------------------------------" << std::endl;
                            m_log->setVisibility(visi);
                            balanced = false;
                            break;
                        }
                    }
                    if(!balanced)
                        break;
                }
                if(!balanced)
                    break;
            }
            if(!balanced)
                break;
        }
        
        bool gBalanced;
#if BITPIT_ENABLE_MPI==1
        MPI_Allreduce(&balanced,&gBalanced,1,MPI_C_BOOL,MPI_LAND,getComm());
#else
        gBalanced = balanced;
#endif

        if(gBalanced){
            (*m_log) << "---------------------------------------------" << std::endl;
            (*m_log) << "CORRECTLY GLOBAL 2:1 BALANCED OCTREE" << std::endl;
            (*m_log) << "---------------------------------------------" << std::endl;
        }
        else{
            (*m_log) << "---------------------------------------------" << std::endl;
            (*m_log) << "UNCORRECTLY GLOBAL 2:1 BALANCED OCTREE" << std::endl;
            (*m_log) << "---------------------------------------------" << std::endl;
        }
        return gBalanced;
    }



#if BITPIT_ENABLE_MPI==1

    /** Distribute Load-Balancing the octants (with user defined weights) of the whole tree over
     * the processes of the job following the Morton order.
     * Until loadBalance is not called for the first time the mesh is serial.
     * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
     */
    void
    ParaTree::loadBalance(dvector* weight){

        //Write info on log
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " LOAD BALANCE " << endl;

        m_lastOp = OP_LOADBALANCE;
        if (m_nproc>1){

            uint32_t* partition = new uint32_t [m_nproc];
            if (weight == NULL)
                computePartition(partition);
            else
                computePartition(partition, weight);

            weight = NULL;

            privateLoadBalance(partition);

            delete [] partition;
            partition = NULL;

            //Write info of final partition on log
            (*m_log) << " " << endl;
            (*m_log) << " Final Parallel partition : " << endl;
            (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
            for(int ii=1; ii<m_nproc; ii++){
                (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
            }
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;

        }
        else{
            (*m_log) << " " << endl;
            (*m_log) << " Serial partition : " << endl;
            (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }

    }

    /** Distribute Load-Balanced the octants (with user defined weights) of the whole tree over
     * the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
     * The families of octants of a desired level are retained compact on the same process.
     * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
     * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
     */
    void
    ParaTree::loadBalance(uint8_t & level, dvector* weight){

        //Write info on log
        (*m_log) << "---------------------------------------------" << endl;
        (*m_log) << " LOAD BALANCE " << endl;

        m_lastOp = OP_LOADBALANCE;
        if (m_nproc>1){

            uint32_t* partition = new uint32_t [m_nproc];
            computePartition(partition, level, weight);

            privateLoadBalance(partition);

            delete [] partition;
            partition = NULL;

            //Write info of final partition on log
            (*m_log) << " " << endl;
            (*m_log) << " Final Parallel partition : " << endl;
            (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
            for(int ii=1; ii<m_nproc; ii++){
                (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
            }
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;

        }
        else{
            (*m_log) << " " << endl;
            (*m_log) << " Serial partition : " << endl;
            (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }

    };

    /**
     * Evaluate the elements of the current partition that will be exchanged
     * with other processors during the load balance.
     *
     * \param[in] weights are the weights of the local octants (if a null
     * pointer is given a uniform distribution is used)
     * \return The ranges of local ids that will be exchanged with other
     * processors.
     */
    ParaTree::LoadBalanceRanges
    ParaTree::evalLoadBalanceRanges(dvector *weights){

        // If there is only one processor no octants can be exchanged
        if (m_nproc == 1) {
            LoadBalanceRanges loadBalanceInfo;
            loadBalanceInfo.sendAction = LoadBalanceRanges::ACTION_NONE;
            loadBalanceInfo.recvAction = LoadBalanceRanges::ACTION_NONE;

            return loadBalanceInfo;
        }

        // Compute updated partition
        std::vector<uint32_t> updatedPartition(m_nproc);
        if (weights) {
            computePartition(updatedPartition.data(), weights);
        } else {
            computePartition(updatedPartition.data());
        }

        // Evaluate send ranges
        return evalLoadBalanceRanges(updatedPartition.data());
    }

    /**
     * Evaluate the elements of the current partition that will be exchanged
     * with other processors during the load balance.
     *
     * The families of octants of a desired level are retained compact on the
     * same process.
     *
     * \param[in] level is the level of the families that will be retained
     * compact on the same process
     * \param[in] weights are the weights of the local octants (if a null
     * pointer is given a uniform distribution is used)
     * \return The ranges of local ids that will be exchanged with other
     * processors.
     */
    ParaTree::LoadBalanceRanges
    ParaTree::evalLoadBalanceRanges(uint8_t level, dvector *weights){

        // If there is only one processor no octants can be sent
        if (m_nproc == 1) {
            LoadBalanceRanges loadBalanceInfo;
            loadBalanceInfo.sendAction = LoadBalanceRanges::ACTION_NONE;
            loadBalanceInfo.recvAction = LoadBalanceRanges::ACTION_NONE;

            return loadBalanceInfo;
        }

        // Compute updated partition
        std::vector<uint32_t> updatedPartition(m_nproc);
        computePartition(updatedPartition.data(), level, weights);

        // Evaluate send ranges
        return evalLoadBalanceRanges(updatedPartition.data());
    }

    /**
     * Evaluate the elements of the current partition that will be exchanged
     * with other processors during the load balance.
     *
     * \param[in] updatePartition is the pointer to the updated pattition
     * \return The ranges of local ids that will be exchanged with other
     * processors.
     */
    ParaTree::LoadBalanceRanges
    ParaTree::evalLoadBalanceRanges(const uint32_t *updatedPartition){

        // If there is only one processor no octants can be sent
        if (m_nproc == 1) {
             LoadBalanceRanges loadBalanceInfo;
            loadBalanceInfo.sendAction = LoadBalanceRanges::ACTION_NONE;
            loadBalanceInfo.recvAction = LoadBalanceRanges::ACTION_NONE;

            return loadBalanceInfo;
        }

        ExchangeRanges sendRanges = evalLoadBalanceSendRanges(updatedPartition);
        ExchangeRanges recvRanges = evalLoadBalanceRecvRanges(updatedPartition);

        return LoadBalanceRanges(m_serial, std::move(sendRanges), std::move(recvRanges));
    }

    /**
     * Evaluate the elements of the current partition that will be sent to
     * other processors after the load balance.
     *
     * \param[in] updatePartition is the pointer to the updated pattition
     * \return The range of local ids that will be sent to other processors.
     */
    ParaTree::ExchangeRanges
    ParaTree::evalLoadBalanceSendRanges(const uint32_t *updatedPartition){

        ExchangeRanges sendRanges;

        // If there is only one processor no octants can be sent
        if (m_nproc == 1) {
            return sendRanges;
        }

        // Compute current partition schema
        std::vector<uint32_t> currentPartition(m_nproc, 0);
        if (!m_serial) {
            currentPartition[0] = m_partitionRangeGlobalIdx[0] + 1;
            for (int i = 1; i < m_nproc; ++i) {
                currentPartition[i] = m_partitionRangeGlobalIdx[i] - m_partitionRangeGlobalIdx[i - 1];
            }
        } else {
            currentPartition[m_rank] = getNumOctants();
        }

        // Get the intersections
        PartitionIntersections globalIntersections = evalPartitionIntersections(currentPartition.data(), m_rank, updatedPartition);

        // Evaluate the send local indexes
        uint64_t offset = 0;
        for (int i = 0; i < m_rank; ++i) {
            offset += currentPartition[i];
        }

        for (const auto &intersectionEntry : globalIntersections) {
            int rank = intersectionEntry.first;
            const std::array<uint64_t, 2> &intersection = intersectionEntry.second;

            std::array<uint32_t, 2> &sendRange = sendRanges[rank];
            sendRange[0] = intersection[0] - offset;
            sendRange[1] = intersection[1] - offset;
        }

        return sendRanges;
    }

    /**
     * Evaluate the elements of the current partition that will be received
     * from other processors after the load balance.
     *
     * \param[in] updatePartition is the pointer to the updated pattition
     * \return The range of local ids that will be received from other
     * processors.
     */
    ParaTree::ExchangeRanges
    ParaTree::evalLoadBalanceRecvRanges(const uint32_t *updatedPartition){

        ExchangeRanges recvRanges;

        // If there is only one processor no octants can be received
        if (m_nproc == 1) {
            return recvRanges;
        }

        // Compute current partition schema
        std::vector<uint32_t> currentPartition(m_nproc, 0);
        if (!m_serial) {
            currentPartition[0] = m_partitionRangeGlobalIdx[0] + 1;
            for (int i = 1; i < m_nproc; ++i) {
                currentPartition[i] = m_partitionRangeGlobalIdx[i] - m_partitionRangeGlobalIdx[i - 1];
            }
        } else {
            currentPartition[m_rank] = getNumOctants();
        }

        // Get the intersections
        PartitionIntersections globalIntersections = evalPartitionIntersections(currentPartition.data(), m_rank, updatedPartition);

        // Evaluate the receive local indexes
        uint64_t offset = 0;
        for (int i = 0; i < m_rank; ++i) {
            offset += updatedPartition[i];
        }

        for (const auto &intersectionEntry : globalIntersections) {
            int rank = intersectionEntry.first;
            const std::array<uint64_t, 2> &intersection = intersectionEntry.second;

            std::array<uint32_t, 2> &recvRange = recvRanges[rank];
            recvRange[0] = intersection[0] - offset;
            recvRange[1] = intersection[1] - offset;
        }

        return recvRanges;
    }

    /**
     * Compute the intersections of the specified partition defined whithin
     * the partition schema A with all the partitions defined whithin the
     * partition schema B.
     *
     * Intersections are evaluated in global indexes.
     *
     * \param[in] partition_A are the number of octants contained in each
     * partition of the partition schema A
     * \param[in] rank_A is the rank associated to the partition for which the
     * intersections will be evaluated
     * \param[in] partition_B are the number of octants contained in each
     * partition of the partition schema B
     * \result The intersections of the specified partition defined whithin
     * the partition schema A with all the partitions defined whithin the
     * partition schema B.
     */
    ParaTree::PartitionIntersections
    ParaTree::evalPartitionIntersections(const uint32_t *partition_A, int rank_A, const uint32_t *partition_B){

        PartitionIntersections intersections;

        // If the partition is empty there are no intersections.
        if (partition_A[rank_A] == 0) {
            return intersections;
        }

        // Calculate partition offsets
        std::vector<uint64_t> offsets_A(m_nproc + 1, 0);
        std::vector<uint64_t> offsets_B(m_nproc + 1, 0);
        for (int i = 0; i < m_nproc; ++i) {
            offsets_A[i + 1] = offsets_A[i] + partition_A[i];
            offsets_B[i + 1] = offsets_B[i] + partition_B[i];
        }

        uint64_t beginGlobalId_A = offsets_A[m_rank];
        uint64_t endGlobalId_A   = offsets_A[m_rank + 1];

        auto firstRankItr = std::upper_bound(offsets_B.begin(), offsets_B.end(), beginGlobalId_A);
        assert(firstRankItr != offsets_B.begin());
        firstRankItr--;

        for (auto itr = firstRankItr; itr != offsets_B.end(); ++itr) {
            int rank_B = std::distance(offsets_B.begin(), itr);
            uint64_t beginGlobalId_B = offsets_B[rank_B];
            uint64_t endGlobalId_B   = offsets_B[rank_B + 1];

            std::array<uint64_t, 2> &intersection = intersections[rank_B];
            intersection[0] = std::max(beginGlobalId_A, beginGlobalId_B);
            intersection[1] = std::min(endGlobalId_A, endGlobalId_B);

            if (endGlobalId_B >= endGlobalId_A) {
                break;
            }
        }

        return intersections;
    }

    /** Distribute Load-Balancing the octants of the whole tree over
     * the processes of the job following a given partition distribution.
     * Until loadBalance is not called for the first time the mesh is serial.
     * \param[in] partition Target distribution of octants over processes.
     */
    void
    ParaTree::privateLoadBalance(uint32_t* partition){

        std::unordered_map<int, std::array<uint32_t, 2>> sendRanges = evalLoadBalanceSendRanges(partition);
        std::unordered_map<int, std::array<uint32_t, 2>> recvRanges = evalLoadBalanceRecvRanges(partition);

        m_loadBalanceRanges = LoadBalanceRanges(m_serial, sendRanges, recvRanges);

        std::array<uint32_t,4> limits = {{0,0,0,0}};

        if(m_serial)
            {
                m_lastOp = OP_LOADBALANCE_FIRST;
                (*m_log) << " " << endl;
                (*m_log) << " Initial Serial distribution : " << endl;
                for(int ii=0; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]+1)) << endl;
                }

                uint32_t stride = 0;
                for(int i = 0; i < m_rank; ++i)
                    stride += partition[i];
                LocalTree::octvector octantsCopy = m_octree.m_octants;
                LocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
                LocalTree::octvector::const_iterator last = first + partition[m_rank];

                limits[1] = stride;
                limits[2] = limits[1] + partition[m_rank];
                limits[3] = m_octree.m_octants.size();
                std::pair<int,std::array<uint32_t,4> > procLimits(m_rank,limits);

                m_octree.m_octants.assign(first, last);
                octvector(m_octree.m_octants).swap(m_octree.m_octants);
                m_octree.m_sizeOctants = m_octree.m_octants.size();

                first = octantsCopy.end();
                last = octantsCopy.end();

                //Update and ghosts here
                updateLoadBalance();
                computeGhostHalo();
            }
        else
            {
                (*m_log) << " " << endl;
                (*m_log) << " Initial Parallel partition : " << endl;
                (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
                for(int ii=1; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(ii))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
                }

                //empty ghosts
                m_octree.m_ghosts.clear();
                m_octree.m_sizeGhosts = 0;
                //compute new partition range globalidx
                assert(m_nproc > 0);
                uint64_t* newPartitionRangeGlobalidx = new uint64_t[m_nproc];
                for(int p = 0; p < m_nproc; ++p){
                    newPartitionRangeGlobalidx[p] = 0;
                    for(int pp = 0; pp <= p; ++pp)
                        newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
                    --newPartitionRangeGlobalidx[p];
                }

                //find resident octants local offset lastHead(lh) and firstTail(ft)
                int32_t lh,ft;
                if(m_rank == 0)
                    lh = -1;
                else{
                    lh = (int32_t)(newPartitionRangeGlobalidx[m_rank-1] + 1 - m_partitionRangeGlobalIdx[m_rank-1] - 1 - 1);
                }
                if(lh < 0)
                    lh = - 1;
                else if(lh > (int32_t)(getNumOctants() - 1))
                    lh = getNumOctants() - 1;

                if(m_rank == m_nproc - 1)
                    ft = getNumOctants();
                else if(m_rank == 0)
                    ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] + 1);
                else{
                    ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1]);
                }
                if(ft > (int32_t)(getNumOctants() - 1))
                    ft = getNumOctants();
                else if(ft < 0)
                    ft = 0;

                //compute size Head and size Tail
                uint32_t headSize = (uint32_t)(lh + 1);
                uint32_t tailSize = (uint32_t)(getNumOctants() - ft);
                uint32_t headOffset = headSize;
                uint32_t tailOffset = tailSize;

                //build send buffers
                DataCommunicator lbCommunicator(m_comm);

                //Compute first predecessor and first successor to send buffers to
                int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
                int64_t globalLastHead = (int64_t) lh;
                int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
                int firstPredecessor = -1;
                int firstSuccessor = m_nproc;
                if(m_rank != 0){
                    firstOctantGlobalIdx = (int64_t)(m_partitionRangeGlobalIdx[m_rank-1] + 1);
                    globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
                    globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
                    for(int pre = m_rank - 1; pre >=0; --pre){
                        if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
                            firstPredecessor = pre;
                    }
                    for(int post = m_rank + 1; post < m_nproc; ++post){
                        if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
                            firstSuccessor = post;
                    }
                }
                else if(m_rank == 0){
                    firstSuccessor = 1;
                }

                int intBuffer = 0;
                int contatore = 0;
                //build send buffers from Head
                uint32_t nofElementsFromSuccessiveToPrevious = 0;
                if(headSize != 0){
                    for(int p = firstPredecessor; p >= 0; --p){
                        if(headSize < partition[p]){
                            intBuffer = (newPartitionRangeGlobalidx[p] - partition[p] );
                            intBuffer = abs(intBuffer);
                            nofElementsFromSuccessiveToPrevious = globalLastHead - intBuffer;
                            if(nofElementsFromSuccessiveToPrevious > headSize || contatore == 1)
                                nofElementsFromSuccessiveToPrevious  = headSize;

                            std::size_t buffSize = (std::size_t)nofElementsFromSuccessiveToPrevious * (std::size_t)Octant::getBinarySize();
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
                            limits[1] = (uint32_t)lh + 1;
                            std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);

                            for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                                //WRITE octants from 0 to lh in sendBuffer[p]
                                sendBuffer << m_octree.m_octants[i];
                            }
                            if(nofElementsFromSuccessiveToPrevious == headSize)
                                break;

                            lh -= nofElementsFromSuccessiveToPrevious;
                            globalLastHead -= nofElementsFromSuccessiveToPrevious;
                            headSize = lh + 1;
                            ++contatore;
                        }
                        else{
                            nofElementsFromSuccessiveToPrevious = globalLastHead - (newPartitionRangeGlobalidx[p] - partition[p]);
                            std::size_t buffSize = (std::size_t)nofElementsFromSuccessiveToPrevious * (std::size_t)Octant::getBinarySize();
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            limits[0] = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1);
                            limits[1] = (uint32_t)lh + 1;
                            std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);

                            for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                                //WRITE octants from lh - partition[p] to lh
                                sendBuffer << m_octree.m_octants[i];
                            }
                            lh -= nofElementsFromSuccessiveToPrevious;
                            globalLastHead -= nofElementsFromSuccessiveToPrevious;
                            headSize = lh + 1;
                            if(headSize == 0)
                                break;
                        }
                    }

                }
                uint32_t nofElementsFromPreviousToSuccessive = 0;
                contatore = 0;
                //build send buffers from Tail
                if(tailSize != 0){
                    for(int p = firstSuccessor; p < m_nproc; ++p){
                        if(tailSize < partition[p]){
                            nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
                            if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
                                nofElementsFromPreviousToSuccessive = tailSize;

                            std::size_t buffSize = (std::size_t)nofElementsFromPreviousToSuccessive * (std::size_t)Octant::getBinarySize();
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            limits[0] = ft;
                            limits[1] = ft + nofElementsFromPreviousToSuccessive;
                            std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);

                            for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
                                //WRITE octants from ft to octantsSize-1
                                sendBuffer << m_octree.m_octants[i];
                            }
                            if(nofElementsFromPreviousToSuccessive == tailSize)
                                break;
                            ft += nofElementsFromPreviousToSuccessive;
                            globalFirstTail += nofElementsFromPreviousToSuccessive;
                            tailSize -= nofElementsFromPreviousToSuccessive;
                            ++contatore;
                        }
                        else{
                            nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
                            std::size_t buffSize = (std::size_t)nofElementsFromPreviousToSuccessive * (std::size_t)Octant::getBinarySize();
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
                            limits[0] = ft;
                            limits[1] = endOctants + 1;
                            std::pair<int,std::array<uint32_t,4> > procLimits(p,limits);

                            for(uint32_t i = ft; i <= endOctants; ++i ){
                                //WRITE octants from ft to ft + partition[p] -1
                                sendBuffer << m_octree.m_octants[i];
                            }
                            ft += nofElementsFromPreviousToSuccessive;
                            globalFirstTail += nofElementsFromPreviousToSuccessive;
                            tailSize -= nofElementsFromPreviousToSuccessive;
                            if(tailSize == 0)
                                break;
                        }
                    }
                }


                lbCommunicator.discoverRecvs();
                lbCommunicator.startAllRecvs();

                uint32_t nofNewHead = 0;
                uint32_t nofNewTail = 0;

                vector<int> recvRanks = lbCommunicator.getRecvRanks();
                std::sort(recvRanks.begin(),recvRanks.end());
                for(auto i : recvRanks){
                    long bufferSize = lbCommunicator.getRecvBuffer(i).getSize();
                    uint32_t nofNewPerProc = (uint32_t)(bufferSize / (uint32_t)Octant::getBinarySize());
                    if(i < m_rank)
                        nofNewHead += nofNewPerProc;
                    else if(i > m_rank)
                        nofNewTail += nofNewPerProc;
                }

                //MOVE RESIDENT TO BEGIN IN OCTANTS
                uint32_t resEnd = getNumOctants() - tailOffset;
                uint32_t nofResidents = resEnd - headOffset;
                int octCounter = 0;
                for(uint32_t i = headOffset; i < resEnd; ++i){
                    m_octree.m_octants[octCounter] = m_octree.m_octants[i];
                    ++octCounter;
                }
                uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
                m_octree.m_octants.resize(newCounter, Octant(m_dim));
                m_octree.m_sizeOctants = m_octree.m_octants.size();
                //MOVE RESIDENTS IN RIGHT POSITION
                uint32_t resCounter = nofNewHead + nofResidents - 1;
                for(uint32_t k = 0; k < nofResidents ; ++k){
                    m_octree.m_octants[resCounter - k] = m_octree.m_octants[nofResidents - k - 1];
                }

                lbCommunicator.startAllSends();

                //READ BUFFERS AND BUILD NEW OCTANTS
                newCounter = 0;
                bool jumpResident = false;
                for(int rank : recvRanks){
                    lbCommunicator.waitRecv(rank);
                    RecvBuffer & recvBuffer = lbCommunicator.getRecvBuffer(rank);
                    long bufferSize = recvBuffer.getSize();
                    uint32_t nofNewPerProc = (uint32_t)(bufferSize / (uint32_t)Octant::getBinarySize());
                    if(rank > m_rank && !jumpResident){
                        newCounter += nofResidents ;
                        jumpResident = true;
                    }
                    for(int i = nofNewPerProc - 1; i >= 0; --i){
                        recvBuffer >> m_octree.m_octants[newCounter];
                        ++newCounter;
                    }
                }
                lbCommunicator.waitAllSends();
                octvector(m_octree.m_octants).swap(m_octree.m_octants);
                m_octree.m_sizeOctants = m_octree.m_octants.size();

                delete [] newPartitionRangeGlobalidx; newPartitionRangeGlobalidx = NULL;
                //Update and ghosts here
                updateLoadBalance();
                computeGhostHalo();
            }
    };
#endif

    /*! Get the size of an octant corresponding to a target level.
     * \param[in] level Input level.
     * \return Size of an octant of input level.
     */
    double
    ParaTree::levelToSize(uint8_t & level) {
        uint32_t size = uint32_t(1)<<(TreeConstants::MAX_LEVEL-level);
        return m_trans.mapSize(size);
    }

    // =================================================================================== //
    // OTHER INTERSECTION BASED METHODS												    			   //
    // =================================================================================== //

    /** Compute the intersection between octants (local, ghost, boundary).
     */
    void
    ParaTree::computeIntersections(){
        m_octree.computeIntersections();
    }

    // =================================================================================== //
    // OTHER PRIVATE METHODS												    			   //
    // =================================================================================== //

    /*! Extract an octant from the local tree.
     * \param[in] idx Local index of target octant.
     * \return Reference to target octant.
     */
    Octant&
    ParaTree::extractOctant(uint32_t idx) {
        return m_octree.extractOctant(idx) ;
    };

    /*! Adapt the octree mesh with user setup for markers and 2:1 balancing conditions.
     * \param[in] mapflag True to track the changes in structure octant by a mapper.
     */
    bool
    ParaTree::private_adapt_mapidx(bool mapflag) {
        //TODO recoding for adapting with abs(marker) > 1

        m_loadBalanceRanges.clear();
        uint32_t nocts0 = getNumOctants();
        vector<Octant >::iterator iter, iterend = m_octree.m_octants.end();

        for (iter = m_octree.m_octants.begin(); iter != iterend; ++iter){
            iter->m_info[Octant::INFO_NEW4REFINEMENT] = false;
            iter->m_info[Octant::INFO_NEW4COARSENING] = false;
        }

        // m_mapIdx init
        u32vector().swap(m_mapIdx);
        if (mapflag) {
            m_mapIdx.resize(nocts0);
            for (uint32_t i=0; i<nocts0; i++){
                m_mapIdx[i] = i;
            }
        }

        bool globalDone = false;
#if BITPIT_ENABLE_MPI==1
        if(m_serial){
#endif
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Refine/Coarse)" << endl;
            (*m_log) << " " << endl;

            // 2:1 Balance
            if (m_lastOp != OP_PRE_ADAPT) {
                balance21(true, false);
            }

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(getNumOctants())) << endl;

            // Refine
            while(m_octree.refine(m_mapIdx));
            if (getNumOctants() > nocts0)
                globalDone = true;
            (*m_log) << " Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(getNumOctants())) << endl;
            nocts0 = getNumOctants();
            updateAdapt();

            // Coarse
            while(m_octree.coarse(m_mapIdx));
            updateAfterCoarse();
            if (getNumOctants() < nocts0){
                globalDone = true;
            }
            nocts0 = getNumOctants();

            (*m_log) << " Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(nocts0)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
#if BITPIT_ENABLE_MPI==1
        }
        else{
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " ADAPT (Refine/Coarse)" << endl;
            (*m_log) << " " << endl;

            // 2:1 Balance
            if (m_lastOp != OP_PRE_ADAPT) {
                balance21(true, false);
            }

            (*m_log) << " " << endl;
            (*m_log) << " Initial Number of octants		:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;

            // Refine
            while(m_octree.refine(m_mapIdx));
            bool localDone = false;
            if (getNumOctants() > nocts0)
                localDone = true;
            updateAdapt();
            computeGhostHalo();
            (*m_log) << " Number of octants after Refine	:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;
            nocts0 = getNumOctants();


            // Coarse
            while(m_octree.coarse(m_mapIdx));
            updateAfterCoarse();
            computeGhostHalo();
            if (getNumOctants() < nocts0){
                localDone = true;
            }
            nocts0 = getNumOctants();

            m_errorFlag = MPI_Allreduce(&localDone,&globalDone,1,MPI_C_BOOL,MPI_LOR,m_comm);
            (*m_log) << " Number of octants after Coarse	:	" + to_string(static_cast<unsigned long long>(m_globalNumOctants)) << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }
#endif

        // Update last operation
        if (mapflag) {
            m_lastOp = OP_ADAPT_MAPPED;
        }
        else{
            m_lastOp = OP_ADAPT_UNMAPPED;
        }

        return globalDone;
    }

    /*!Update the local tree after an adapt.
     */
    void
    ParaTree::updateAdapt(){
#if BITPIT_ENABLE_MPI==1
        if(m_serial)
            {
#endif
                for (int iproc=0; iproc<m_nproc; iproc++){
                    m_partitionRangeGlobalIdx0[iproc] = m_partitionRangeGlobalIdx[iproc];
                }
                m_maxDepth = m_octree.m_localMaxDepth;
                m_globalNumOctants = getNumOctants();
                for(int p = 0; p < m_nproc; ++p){
                    m_partitionRangeGlobalIdx[p] = m_globalNumOctants - 1;
                }
                m_internals.resize(getNumOctants());
                int i = 0;
                octvector::iterator itend = m_octree.m_octants.end();
                for (octvector::iterator it = m_octree.m_octants.begin(); it != itend; ++it){
                    m_internals[i] = &(*it);
                    i++;
                }
#if BITPIT_ENABLE_MPI==1
            }
        else
            {
                for (int iproc=0; iproc<m_nproc; iproc++){
                    m_partitionRangeGlobalIdx0[iproc] = m_partitionRangeGlobalIdx[iproc];
                }
                //update m_maxDepth
                m_errorFlag = MPI_Allreduce(&m_octree.m_localMaxDepth,&m_maxDepth,1,MPI_UINT8_T,MPI_MAX,m_comm);
                //update m_globalNumOctants
                uint64_t local_num_octants = getNumOctants();
                m_errorFlag = MPI_Allreduce(&local_num_octants,&m_globalNumOctants,1,MPI_UINT64_T,MPI_SUM,m_comm);
                //update m_partitionRangeGlobalIdx
                uint64_t* rbuff = new uint64_t[m_nproc];
                m_errorFlag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,rbuff,1,MPI_UINT64_T,m_comm);
                for(int p = 0; p < m_nproc; ++p){
                    m_partitionRangeGlobalIdx[p] = 0;
                    for(int pp = 0; pp <=p; ++pp)
                        m_partitionRangeGlobalIdx[p] += rbuff[pp];
                    --m_partitionRangeGlobalIdx[p];
                }
                delete [] rbuff; rbuff = NULL;
            }
#endif
    }

#if BITPIT_ENABLE_MPI==1
    /*! Compute the partition of the octree over the processes (only compute the information about
     * how distribute the mesh). This is an uniform distribution method.
     * \param[out] partition Pointer to partition information array. partition[i] = number of octants
     * to be stored on the i-th process (i-th rank).
     */
    void
    ParaTree::computePartition(uint32_t* partition){

        uint32_t division_result = 0;
        uint32_t remind = 0;

        division_result = uint32_t(m_globalNumOctants/(uint64_t)m_nproc);
        remind = (uint32_t)(m_globalNumOctants%(uint64_t)m_nproc);

        for(uint32_t i = 0; i < (uint32_t)m_nproc; ++i)
            if(i<remind)
                partition[i] = division_result + 1;
            else
                partition[i] = division_result;

    }

    /*! Compute the partition of the octree over the processes (only compute the information about
     * how distribute the mesh). This is an weighted distribution method: each process will have the same weight.
     * \param[out] partition Pointer to partition information array. partition[i] = number of octants
     * to be stored on the i-th process (i-th rank).
     * \param[in] weight Pointer to weight array. weight[i] = weight of i-th local octant.
     */
    void
    ParaTree::computePartition(uint32_t* partition, dvector* weight){
        if(m_serial){

            double division_result = 0;
            double global_weight = 0.0;
            for (unsigned int i=0; i<weight->size(); i++){
                global_weight += (*weight)[i];
            }
            division_result = global_weight/(double)m_nproc;

            //Estimate resulting weight distribution starting from proc 0 (sending tail)
            //Estimate sending weight by each proc in initial conf (sending tail)
            uint32_t i = 0, tot = 0;
            int iproc = 0;
            while (iproc < m_nproc-1){
                double partial_weight = 0.0;
                partition[iproc] = 0;
                while(partial_weight < division_result){
                    partial_weight += (*weight)[i];
                    tot++;
                    partition[iproc]++;
                    i++;
                }
                iproc++;
            }
            partition[m_nproc-1] = weight->size() - tot;
        }
        else{

            int weightSize = weight->size();
            double* gweight;
            double* lweight = new double[weightSize];

            for (unsigned int i=0; i<weight->size(); i++){
                lweight[i] = (*weight)[i];
            }

            int *oldpartition = new int[m_nproc];
            int *displays = new int[m_nproc];
            MPI_Allgather(&weightSize,1,MPI_INT,oldpartition,1,MPI_INT,m_comm);
            int globalNofOctant = 0;
            for(int i = 0; i < m_nproc; ++i){
                displays[i] = globalNofOctant;
                globalNofOctant += oldpartition[i];
            }
            gweight = new double[globalNofOctant];
            MPI_Allgatherv(lweight,weightSize,MPI_DOUBLE,gweight,oldpartition,displays,MPI_DOUBLE,m_comm);

            double division_result = 0;
            double global_weight = 0.0;
            for (int i=0; i<globalNofOctant; i++){
                global_weight += gweight[i];
            }
            division_result = global_weight/(double)m_nproc;

            //Estimate resulting weight distribution starting from proc 0 (sending tail)
            //Estimate sending weight by each proc in initial conf (sending tail)
            uint32_t i = 0, tot = 0;
            int iproc = 0;
            while (iproc < m_nproc-1){
                double partial_weight = 0.0;
                partition[iproc] = 0;
                while(partial_weight < division_result && (int32_t) i < globalNofOctant){
                    partial_weight += gweight[i];
                    tot++;
                    partition[iproc]++;
                    i++;
                }
                global_weight = 0;
                for(int j = i; j < globalNofOctant; ++j)
                    global_weight += gweight[j];
                division_result = global_weight/double(m_nproc-(iproc+1));
                iproc++;
            }
            partition[m_nproc-1] = globalNofOctant - tot;

            delete [] oldpartition;
            delete [] displays;
            delete [] lweight;
            delete [] gweight;

            //TODO CHECK OLD ALGORITHM
            //		double division_result = 0;
            //		double remind = 0;
            //		dvector local_weight(m_nproc,0.0);
            //		dvector temp_local_weight(m_nproc,0.0);
            //		dvector2D sending_weight(m_nproc, dvector(m_nproc,0.0));
            //		double* rbuff = new double[m_nproc];
            //		double global_weight = 0.0;
            //		for (int i=0; i<weight->size(); i++){
            //			local_weight[m_rank] += (*weight)[i];
            //		}
            //		m_errorFlag = MPI_Allgather(&local_weight[m_rank],1,MPI_DOUBLE,rbuff,1,MPI_DOUBLE,m_comm);
            //		for (int i=0; i<m_nproc; i++){
            //			local_weight[i] = rbuff[i];
            //			global_weight += rbuff[i];
            //		}
            //		delete [] rbuff; rbuff = NULL;
            //		division_result = global_weight/(double)m_nproc;
            //
            //		//Estimate resulting weight distribution starting from proc 0 (sending tail)
            //
            //		temp_local_weight = local_weight;
            //		//Estimate sending weight by each proc in initial conf (sending tail)
            //
            //		for (int iter = 0; iter < 1; iter++){
            //
            //			vector<double> delta(m_nproc);
            //			for (int i=0; i<m_nproc; i++){
            //				delta[i] = temp_local_weight[i] - division_result;
            //			}
            //
            //			for (int i=0; i<m_nproc-1; i++){
            //
            //				double post_weight = 0.0;
            //				for (int j=i+1; j<m_nproc; j++){
            //					post_weight += temp_local_weight[j];
            //				}
            //				if (temp_local_weight[i] > division_result){
            //
            //					delta[i] = temp_local_weight[i] - division_result;
            //					if (post_weight < division_result*(m_nproc-i-1)){
            //
            //						double post_delta =  division_result*(m_nproc-i-1) - post_weight;
            //						double delta_sending = min(local_weight[i], min(delta[i], post_delta));
            //						int jproc = i+1;
            //						double sending = 0;
            //						while (delta_sending > 0 && jproc<m_nproc){
            //							sending = min(division_result, delta_sending);
            //							sending = min(sending, (division_result-temp_local_weight[jproc]));
            //							sending = max(sending, 0.0);
            //							sending_weight[i][jproc] += sending;
            //							temp_local_weight[jproc] += sending;
            //							temp_local_weight[i] -= sending;
            //							delta_sending -= sending;
            //							delta[i] -= delta_sending;
            //							jproc++;
            //						}
            //					} //post
            //				}//weight>
            //			}//iproc
            //
            //			for (int i = m_nproc-1; i>0; i--){
            //
            //				double pre_weight = 0.0;
            //				for (int j=i-1; j>=0; j--){
            //					pre_weight += temp_local_weight[j];
            //				}
            //				if (temp_local_weight[i] > division_result){
            //
            //					delta[i] = temp_local_weight[i] - division_result;
            //					if (pre_weight < division_result*(i)){
            //
            //						double pre_delta =  division_result*(i) - pre_weight;
            //						double delta_sending = min(local_weight[i], min(temp_local_weight[i], min(delta[i], pre_delta)));
            //						int jproc = i-1;
            //						double sending = 0;
            //						while (delta_sending > 0 && jproc >=0){
            //							sending = min(division_result, delta_sending);
            //							sending = min(sending, (division_result-temp_local_weight[jproc]));
            //							sending = max(sending, 0.0);
            //							sending_weight[i][jproc] += sending;
            //							temp_local_weight[jproc] += sending;
            //							temp_local_weight[i] -= sending;
            //							delta_sending -= sending;
            //							delta[i] -= delta_sending;
            //							jproc--;
            //						}
            //					}//pre
            //				}//weight>
            //			}//iproc
            //		}//iter
            //
            //		//Update partition locally
            //		//to send
            //		u32vector sending_cell(m_nproc,0);
            //		int i = getNumOctants();;
            //		for (int jproc=m_nproc-1; jproc>m_rank; jproc--){
            //			double pack_weight = 0.0;
            //			while(pack_weight < sending_weight[m_rank][jproc] && i > 0){
            //				i--;
            //				pack_weight += (*weight)[i];
            //				sending_cell[jproc]++;
            //			}
            //		}
            //		partition[m_rank] = i;
            //		i = 0;
            //		for (int jproc=0; jproc<m_rank; jproc++){
            //			double pack_weight = 0.0;
            //			while(pack_weight < sending_weight[m_rank][jproc] && i <  getNumOctants()-1){
            //				i++;
            //				pack_weight += (*weight)[i];
            //				sending_cell[jproc]++;
            //			}
            //		}
            //		partition[m_rank] -= i;
            //
            //		//to receive
            //		u32vector rec_cell(m_nproc,0);
            //		MPI_Request* req = new MPI_Request[m_nproc*10];
            //		MPI_Status* stats = new MPI_Status[m_nproc*10];
            //		int nReq = 0;
            //		for (int iproc=0; iproc<m_nproc; iproc++){
            //			m_errorFlag = MPI_Irecv(&rec_cell[iproc],1,MPI_UINT32_T,iproc,m_rank,m_comm,&req[nReq]);
            //			++nReq;
            //		}
            //		for (int iproc=0; iproc<m_nproc; iproc++){
            //			m_errorFlag =  MPI_Isend(&sending_cell[iproc],1,MPI_UINT32_T,iproc,iproc,m_comm,&req[nReq]);
            //			++nReq;
            //		}
            //		MPI_Waitall(nReq,req,stats);
            //
            //		delete [] req; req = NULL;
            //		delete [] stats; stats = NULL;
            //
            //		i = 0;
            //		for (int jproc=0; jproc<m_nproc; jproc++){
            //			i+= rec_cell[jproc];
            //		}
            //		partition[m_rank] += i;
            //		uint32_t part = partition[m_rank];
            //		m_errorFlag = MPI_Allgather(&part,1,MPI_UINT32_T,partition,1,MPI_UINT32_T,m_comm);

        }
    };

    /*! Compute the partition of the octree over the processes (only compute the information about
     * how distribute the mesh). This is a "compact families" method: the families of octants
     * of a desired level are retained compact on the same process.
     * \param[out] partition Pointer to partition information array. partition[i] = number of octants
     * to be stored on the i-th process (i-th rank).
     * \param[in] level_ Number of level over the max depth reached in the tree at
     * which families of octants are fixed compact on the same process
     * (level=0 is uniform partition).
     * \param[in] weight Pointer to weight array. weight[i] = weight of i-th local octant.
     */
    void
    ParaTree::computePartition(uint32_t* partition, uint8_t & level_, dvector* weight) {

        uint8_t level = uint8_t(min(int(max(int(m_maxDepth) - int(level_), int(1))) , int(TreeConstants::MAX_LEVEL)));
        uint32_t* partition_temp = new uint32_t[m_nproc];
        uint8_t* boundary_proc = new uint8_t[m_nproc-1];
        uint8_t dimcomm, indcomm;
        uint8_t* glbdimcomm = new uint8_t[m_nproc];
        uint8_t* glbindcomm = new uint8_t[m_nproc];

        uint32_t Dh = uint32_t(pow(double(2),double(TreeConstants::MAX_LEVEL-level)));
        uint32_t istart, nocts, rest, forw, backw;
        uint32_t i = 0, iproc, j;
        uint64_t sum;
        int32_t* pointercomm;
        int32_t* deplace = new int32_t[m_nproc-1];

        if (weight==NULL){
            computePartition(partition_temp);
        }
        else{
            computePartition(partition_temp, weight);
        }

        j = 0;
        sum = 0;
        for (iproc=0; iproc<(uint32_t)(m_nproc-1); iproc++){
            sum += partition_temp[iproc];
            while(sum > m_partitionRangeGlobalIdx[j]){
                j++;
            }
            boundary_proc[iproc] = j;
        }
        nocts = getNumOctants();
        sum = 0;
        dimcomm = 0;
        indcomm = 0;
        for (iproc=0; iproc<(uint32_t)(m_nproc-1); iproc++){
            deplace[iproc] = 1;
            sum += partition_temp[iproc];
            if (boundary_proc[iproc] == m_rank){
                if (dimcomm == 0){
                    indcomm = iproc;
                }
                dimcomm++;
                if (m_rank!=0)
                    istart = sum - m_partitionRangeGlobalIdx[m_rank-1] - 1;
                else
                    istart = sum;

                i = istart;
                rest = m_octree.m_octants[i].getX()%Dh + m_octree.m_octants[i].getY()%Dh;
                while(rest!=0){
                    if (i==nocts){
                        i = istart + nocts;
                        break;
                    }
                    i++;
                    rest = m_octree.m_octants[i].getX()%Dh + m_octree.m_octants[i].getY()%Dh;
                }
                forw = i - istart;
                i = istart;
                rest = m_octree.m_octants[i].getX()%Dh + m_octree.m_octants[i].getY()%Dh;
                while(rest!=0){
                    if (i==0){
                        i = istart - nocts;
                        break;
                    }
                    i--;
                    rest = m_octree.m_octants[i].getX()%Dh + m_octree.m_octants[i].getY()%Dh;
                }
                backw = istart - i;
                if (forw<backw)
                    deplace[iproc] = forw;
                else
                    deplace[iproc] = -(int32_t)backw;
            }
        }

        m_errorFlag = MPI_Allgather(&dimcomm,1,MPI_UINT8_T,glbdimcomm,1,MPI_UINT8_T,m_comm);
        m_errorFlag = MPI_Allgather(&indcomm,1,MPI_UINT8_T,glbindcomm,1,MPI_UINT8_T,m_comm);
        for (iproc=0; iproc<(uint32_t)(m_nproc); iproc++){
            pointercomm = &deplace[glbindcomm[iproc]];
            m_errorFlag = MPI_Bcast(pointercomm, glbdimcomm[iproc], MPI_INT32_T, iproc, m_comm);
        }

        for (iproc=0; iproc<(uint32_t)(m_nproc); iproc++){
            if (iproc < (uint32_t)(m_nproc-1))
                partition[iproc] = partition_temp[iproc] + deplace[iproc];
            else
                partition[iproc] = partition_temp[iproc];
            if (iproc !=0)
                partition[iproc] = partition[iproc] - deplace[iproc-1];
        }

        delete [] partition_temp; partition_temp = NULL;
        delete [] boundary_proc; boundary_proc = NULL;
        delete [] glbdimcomm; glbdimcomm = NULL;
        delete [] glbindcomm; glbindcomm = NULL;
        delete [] deplace; deplace = NULL;
    }

    /*! Update the distributed octree after a LoadBalance over the processes.
     */
    void
    ParaTree::updateLoadBalance() {
        //update sizes
        m_octree.m_sizeOctants = m_octree.m_octants.size();

        m_octree.updateLocalMaxDepth();
        uint64_t* rbuff = new uint64_t[m_nproc];

        uint64_t local_num_octants = getNumOctants();
        m_errorFlag = MPI_Allgather(&local_num_octants,1,MPI_UINT64_T,rbuff,1,MPI_UINT64_T,m_comm);
        for (int iproc=0; iproc<m_nproc; iproc++){
            m_partitionRangeGlobalIdx0[iproc] = m_partitionRangeGlobalIdx[iproc];
        }
        for(int p = 0; p < m_nproc; ++p){
            m_partitionRangeGlobalIdx[p] = 0;
            for(int pp = 0; pp <=p; ++pp)
                m_partitionRangeGlobalIdx[p] += rbuff[pp];
            --m_partitionRangeGlobalIdx[p];
        }
        //update first last descendant
        if(getNumOctants()==0){
            Octant octDesc(m_dim,TreeConstants::MAX_LEVEL,pow(2,TreeConstants::MAX_LEVEL),pow(2,TreeConstants::MAX_LEVEL),(m_dim > 2 ? pow(2,TreeConstants::MAX_LEVEL) : 0));
            m_octree.m_lastDescMorton = octDesc.computeMorton();
            m_octree.m_firstDescMorton = std::numeric_limits<uint64_t>::max();
        }
        else{
            m_octree.setFirstDescMorton();
            m_octree.setLastDescMorton();
        }

        //update partition_range_position
        uint64_t lastDescMorton = m_octree.getLastDescMorton();
        m_errorFlag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,m_partitionLastDesc.data(),1,MPI_UINT64_T,m_comm);
        uint64_t firstDescMorton = m_octree.getFirstDescMorton();
        m_errorFlag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,m_partitionFirstDesc.data(),1,MPI_UINT64_T,m_comm);
        m_serial = false;

        delete [] rbuff; rbuff = NULL;
    }

    /*! Build the structure with the information about the first layer of ghost octants, partition boundary octants
     *  and parameters for communicate between processes.
     */
    void
    ParaTree::setPboundGhosts() {
        //BUILD BORDER OCTANT INDECES VECTOR (map value) TO BE SENT TO THE RIGHT PROCESS (map key)
        //find local octants to be sent as ghost to the right processes
        //it visits the local octants building virtual neighbors on each octant face
        //find the owner of these virtual neighbor and build a map (process,border octants)
        //this map contains the local octants as ghosts for neighbor processes

        // NO PBORDERS !
        //
        // TODO: provide an estimate of the border octants in order to reserve
        // the vectors that will contain them.
        m_bordersPerProc.clear();
        m_internals.resize(getNumOctants());
        m_pborders.resize(getNumOctants());

        int countpbd = 0;
        int countint = 0;
        std::set<int> neighProcs;
        uint32_t nVirtualNeighbors;
        std::vector<uint64_t> virtualNeighbors;
        for (uint32_t idx = 0; idx < m_octree.getNumOctants(); ++idx) {
            neighProcs.clear();
            Octant &octant = m_octree.m_octants[idx];

            // Virtual Face Neighbors
            for(uint8_t i = 0; i < m_treeConstants->nFaces; ++i){
                bool isFacePbound = false;
                if(octant.getBound(i) == false){
                    octant.computeFaceVirtualMortons(i, m_maxDepth, &nVirtualNeighbors, &virtualNeighbors);
                    uint32_t maxDelta = nVirtualNeighbors/2;
                    for(uint32_t j = 0; j <= maxDelta; ++j){
                        int neighProcFirst = findOwner(virtualNeighbors[j]);
                        if (neighProcFirst != m_rank) {
                            neighProcs.insert(neighProcFirst);
                            isFacePbound = true;
                        }

                        int neighProcLast = findOwner(virtualNeighbors[nVirtualNeighbors - 1 - j]);
                        if (neighProcLast != m_rank) {
                            neighProcs.insert(neighProcLast);
                            isFacePbound = true;
                        }

                        //					//TODO debug
                        //					if (abs(pBegin-pEnd) <= 1) j = maxDelta + 1;
                    }
                }
                else if(m_periodic[i]){
                    uint64_t virtualNeighbor = octant.computePeriodicMorton(i);
                    int neighProc = findOwner(virtualNeighbor);
                    if(neighProc != m_rank){
                        neighProcs.insert(neighProc);
                        isFacePbound = true;
                    }
                }

                octant.setPbound(i, isFacePbound);
            }

            // Virtual Edge Neighbors
            for(uint8_t e = 0; e < m_treeConstants->nEdges; ++e){
                octant.computeEdgeVirtualMortons(e, m_maxDepth, m_octree.m_balanceCodim, m_treeConstants->edgeFace, &nVirtualNeighbors, &virtualNeighbors);
                if(nVirtualNeighbors > 0){
                    uint32_t maxDelta = nVirtualNeighbors/2;
                    for(uint32_t ee = 0; ee <= maxDelta; ++ee){
                        int neighProcFirst = findOwner(virtualNeighbors[ee]);
                        if (neighProcFirst != m_rank) {
                            neighProcs.insert(neighProcFirst);
                        }

                        int neighProcLast = findOwner(virtualNeighbors[nVirtualNeighbors - 1- ee]);
                        if (neighProcLast != m_rank) {
                            neighProcs.insert(neighProcLast);
                        }

                        //					//TODO debug
                        //					if (abs(pBegin-pEnd) <= 1) ee = maxDelta + 1;
                    }
                }
            }

            // Virtual Corner Neighbors
            for(uint8_t c = 0; c < m_treeConstants->nNodes; ++c){
                bool hasVirtualNeighbour;
                uint64_t virtualNeighbor;
                octant.computeNodeVirtualMorton(c, m_maxDepth,m_treeConstants->nodeFace, &hasVirtualNeighbour, &virtualNeighbor);
                if(hasVirtualNeighbour){
                    int neighProc = findOwner(virtualNeighbor);
                    if (neighProc != m_rank) {
                        neighProcs.insert(neighProc);
                    }
                }
            }

            // Build list of internal and processor borders octants
            if (neighProcs.empty()){
                m_internals[countint] = &octant;
                countint++;
            } else {
                m_pborders[countpbd] = &octant;
                countpbd++;

                for (int neighProc : neighProcs) {
                    assert(neighProc != m_rank);
                    m_bordersPerProc[neighProc].push_back(idx);
                }
            }
        }
        m_pborders.shrink_to_fit();
        m_internals.shrink_to_fit();

        // Build ghosts
        buildGhostOctants(m_bordersPerProc, std::vector<AccretionData>());
    }

    /*! Build the structure with the information about the layers (from the second one) of ghost octants, partition boundary octants
     *  and parameters for communicate between processes.
     */
    void
    ParaTree::computeGhostHalo(){
        // Build first layer of ghosts
        setPboundGhosts();

        // Early return if we need to build only one layer
        if(m_nofGhostLayers <= 1){
            return;
        }

        //
        // Accrete sources
        //
        // We don't build ghost layers directly, instead we identify the
        // internal cells that are ghosts for the neighboring process and
        // we use this list to create the ghosts. We use the term "sources"
        // to identify internal cells that are ghosts for the neighboring
        // process. For each layer of ghosts, a corresponding layer of
        // sources exists.
        //
        // Sources are identified one layer at a time. The first layer is
        // already known: the processor-border octants. The neighbors of
        // processor-border octants are the second layer of sources; the
        // neighbors of the second layer of sources are the third layer,
        // and so on an so forth.
        //
        // To identify the sources, an auxiliary data structure is used. This
        // data structure is called accretion and contains the list of sources
        // currently identified (population), a list of octants to be used for
        // building the next layer of sources (seeds) and the rank on which
        // the sources gathered by the accretion will be ghosts.
        //
        // The identification of the sources starts creating one accretions for
        // each of the neighboring processors. The accretions are initialized
        // using the processor-border octants already build: those octants are
        // the first layer of sources and the seeds for the generation of the
        // second layer. Adding the internal neighbors of the internal seeds to
        // the population, accretions are grown one layer at a time. When an
        // accretion reaches a neighboring processors (i.e., when a first-layer
        // ghost enters in the list of foreign seeds), we communicate to the
        // owner of the ghost to create a new accretion and continue the search
        // for the sources. At the end of the procedure, the population of the
        // accretions on each processor will contain the desired sources.

        // Initialize cache for 1-rings of the internal octants
        std::unordered_map<uint32_t, std::vector<uint64_t>> oneRingsCache;
        oneRingsCache.reserve(getNumOctants());

        // Initialize data communicator
        DataCommunicator accretionDataCommunicator(m_comm);

        // Initialize accretions
        std::vector<AccretionData> accretions;
        initializeGhostHaloAccretions(&accretions);

        // Grow the accretions
        for(std::size_t layer = 1; layer < m_nofGhostLayers; ++layer){
            // Exchange accretions
            //
            // When a ghost is incorporated in the seeds, the accretion
            // needs to continue on the processor that owns the ghost.
            exchangeGhostHaloAccretions(&accretionDataCommunicator, &accretions);

            // Grow accretions
            growGhostHaloAccretions(&oneRingsCache, &accretions);
        }

        // To correctly identify the population of the last layer of sources,
        // we need to exchange the accretions one more time. This allows to
        // communicate the foreign seeds found during the last growth to the
        // processors that own them.
        exchangeGhostHaloAccretions(&accretionDataCommunicator, &accretions);

        //
        // Extract list of sources
        //
        // Sources are internal octants that are ghosts for other processors,
        // i.e., internal octants on processors borders (pborder octants).
        // The population of the accretions is exaclty made of the sources
        // for the target rank of the accretion.
        for(const AccretionData &accretion : accretions){
            int targetRank = accretion.targetRank;

            std::vector<uint32_t> &rankBordersPerProc = m_bordersPerProc[targetRank];
            rankBordersPerProc.resize(accretion.population.size());

            std::size_t index = 0;
            for (const auto &populationEntry : accretion.population) {
                uint64_t globalIdx = populationEntry.first;
                uint32_t localIdx = getLocalIdx(globalIdx);
                rankBordersPerProc[index] = localIdx;
                ++index;
            }

            std::sort(rankBordersPerProc.begin(), rankBordersPerProc.end());
        }

        //
        // Build the ghosts
        //
        buildGhostOctants(m_bordersPerProc, accretions);
    }

    /*! Initialize the accretions.
     *
     * Accretions are auxiliary data structures needed for generation of ghost
     * halo. An explanation of how ghost halo is generated can be found in the
     * function that builds the ghost layers.
     *
     * \param[in,out] accretions list of accretions
     */
    void
    ParaTree::initializeGhostHaloAccretions(std::vector<AccretionData> *accretions) {

        const std::size_t FIRST_LAYER = 0;

        accretions->reserve(m_bordersPerProc.size());
        for(const auto &bordersPerProcEntry : m_bordersPerProc){
            accretions->emplace_back();
            AccretionData &accretion = accretions->back();

            // Rank for which the accretion will gather sources
            int targetRank = bordersPerProcEntry.first;
            accretion.targetRank = targetRank;

            // Initialize the first layer
            //
            // The population and the seeds of the first layer are the octants
            // on processors borders.
            const std::vector<uint32_t> &rankBordersPerProc = bordersPerProcEntry.second;
            const std::size_t nRankBordersPerProc = rankBordersPerProc.size();

            accretion.population.reserve(m_nofGhostLayers * nRankBordersPerProc);
            accretion.internalSeeds.reserve(nRankBordersPerProc);
            accretion.foreignSeeds.reserve(nRankBordersPerProc);
            for(uint32_t pborderLocalIdx : rankBordersPerProc){
                uint64_t pborderGlobalIdx = getGlobalIdx(pborderLocalIdx);
                if (isInternal(pborderGlobalIdx)) {
                    accretion.population.insert({pborderGlobalIdx, FIRST_LAYER});
                    accretion.internalSeeds.insert({pborderGlobalIdx, FIRST_LAYER});
                } else {
                    accretion.foreignSeeds.insert({pborderGlobalIdx, FIRST_LAYER});
                }
            }
        }
    }

    /*! Grow the accretions adding a new layer of seeds.
     *
     * Accretions are auxiliary data structures needed for generation of ghost
     * halo. An explanation of how ghost halo is generated can be found in the
     * function that builds the ghost layers.
     *
     * \param[in,out] oneRingsCache cache for one ring evaluation
     * \param[in,out] accretions list of accretions
     */
    void
    ParaTree::growGhostHaloAccretions(std::unordered_map<uint32_t, std::vector<uint64_t>> *oneRingsCache,
                                      std::vector<AccretionData> *accretions) {

        // The neighbour of the internal seeds are the next layer of sources.
        for(AccretionData &accretion : *accretions){
            // If the accretion doesn't have internal seeds we can skip it
            std::size_t nSeeds = accretion.internalSeeds.size();
            if (nSeeds == 0) {
                continue;
            }

            // Rank for which accretion is gathering data
            const int targetRank = accretion.targetRank;

            // Seeds
            //
            // We make a copy of the seeds and then we clear the original
            // list in order to generate the seeds for the next layer.
            std::unordered_map<uint64_t, int> currentInternalSeeds;
            currentInternalSeeds.reserve(accretion.internalSeeds.size());
            accretion.internalSeeds.swap(currentInternalSeeds);

            // The next layer is obtained adding the 1-ring neighbours of
            // the internal octants of the previous layer.
            for(const auto &seedEntry : currentInternalSeeds){
                // Get seed information
                int seedLayer = seedEntry.second;
                uint64_t seedGlobalIdx = seedEntry.first;
                uint32_t seedLocalIdx = getLocalIdx(seedGlobalIdx);

                // Find the 1-ring of the source
                auto oneRingsCacheItr = oneRingsCache->find(seedLocalIdx);
                if (oneRingsCacheItr == oneRingsCache->end()) {
                    oneRingsCacheItr = oneRingsCache->insert({seedLocalIdx, std::vector<uint64_t>()}).first;
                    findAllGlobalNeighbours(seedLocalIdx, oneRingsCacheItr->second);
                    oneRingsCacheItr->second.push_back(getGlobalIdx(seedLocalIdx));
                }
                const std::vector<uint64_t> &seedOneRing = oneRingsCacheItr->second;

                // Add the 1-ring of the octant to the sources
                for(uint64_t neighGlobalIdx : seedOneRing){
                    // Discard octants already in the population
                    if (accretion.population.count(neighGlobalIdx) > 0) {
                        continue;
                    }

                    // Get neighbour information
                    bool isNeighInternal = isInternal(neighGlobalIdx);

                    int neighRank;
                    if (isNeighInternal) {
                        neighRank = getRank();
                    } else {
                        neighRank = getOwnerRank(neighGlobalIdx);
                    }

                    // Add the neighbour to the population
                    //
                    // Population should only contain internal octants.
                    if (isNeighInternal) {
                        accretion.population.insert({neighGlobalIdx, seedLayer + 1});
                    }

                    // Add the neighbour to the seeds
                    //
                    // Internal seeds contains only internal octants, foreign
                    // seeds containt ghost octants that are not owned by the
                    // target rank (if an octant is owned by the target rank,
                    // by definition it will not be a source for that rank).
                    if (isNeighInternal) {
                        accretion.internalSeeds.insert({neighGlobalIdx, seedLayer + 1});
                    } else if (neighRank != targetRank) {
                        accretion.foreignSeeds.insert({neighGlobalIdx, seedLayer + 1});
                    }
                }
            }
        }
    }

    /*! Exchange the accretions among neighbouring processors.
     *
     * Accretions are auxiliary data structures needed for generation of ghost
     * halo. An explanation of how ghost halo is generated can be found in the
     * function that builds the ghost layers.
     *
     * \param[in,out] dataCommunicator is the data communicator that will
     * handle data exchagne
     * \param[in,out] accretions list of accretions
     */
    void
    ParaTree::exchangeGhostHaloAccretions(DataCommunicator *dataCommunicator,
                                          std::vector<AccretionData> *accretions) {

        // Generate accretions that has to be sent to other processors
        //
        // When the accretion reaches a foreign process (i.e., a ghost is
        // added to the seeds), the ranks that owns the ghost seed have to
        // continue the propagation of those seeds (which will be local
        // seeds for that foreign rank) locally. Therefore we need to
        // communicate to the ranks that own ghosts seeds the information
        // about the accretion and the list of seeds.
        std::unordered_map<int, std::vector<AccretionData>> foreignAccretions;
        for(const AccretionData &accretion : *accretions){
            for(const auto &seedEntry : accretion.foreignSeeds){
                // Find the foreign accretion the ghost seed belogns to
                int seedRank = getOwnerRank(seedEntry.first);
                std::vector<AccretionData> &foreignRankAccretions = foreignAccretions[seedRank];

                auto foreignRankAccretionsItr = foreignRankAccretions.begin();
                for (; foreignRankAccretionsItr != foreignRankAccretions.end(); ++foreignRankAccretionsItr) {
                    if (foreignRankAccretionsItr->targetRank!= accretion.targetRank) {
                        continue;
                    }

                    break;
                }

                if (foreignRankAccretionsItr == foreignRankAccretions.end()) {
                    foreignRankAccretions.emplace_back();
                    foreignRankAccretionsItr = foreignRankAccretions.begin() + foreignRankAccretions.size() - 1;
                    foreignRankAccretionsItr->targetRank = accretion.targetRank;
                }

                AccretionData &foreignAccretion = *foreignRankAccretionsItr;

                // Add the seed
                foreignAccretion.internalSeeds.insert(seedEntry);
            }
        }

        // Early return if no communications are needed
        bool foreignAccrationExchangeNeeded = (foreignAccretions.size() > 0);
        MPI_Allreduce(MPI_IN_PLACE, &foreignAccrationExchangeNeeded, 1, MPI_C_BOOL, MPI_LOR, m_comm);
        if (!foreignAccrationExchangeNeeded) {
            return;
        }

        // Clear previous communications
        dataCommunicator->clearAllSends();
        dataCommunicator->clearAllRecvs();

        // Fill send buffers with accretions data
        for(const auto &foreignAccretionEntry : foreignAccretions){
            int receiverRank = foreignAccretionEntry.first;

            // Evaluate buffer size
            std::size_t buffSize = 0;
            buffSize += sizeof(std::size_t);
            for(const auto &foreignAccretion: foreignAccretionEntry.second){
                std::size_t nSeeds = foreignAccretion.internalSeeds.size();

                buffSize += sizeof(int);
                buffSize += sizeof(std::size_t);
                buffSize += nSeeds * (sizeof(uint64_t) + sizeof(int));
            }

            dataCommunicator->setSend(receiverRank, buffSize);

            // Fill buffer
            SendBuffer &sendBuffer = dataCommunicator->getSendBuffer(receiverRank);

            const std::vector<AccretionData> &foreignRankAccretions = foreignAccretionEntry.second;

            sendBuffer << foreignRankAccretions.size();
            for(const auto &foreignAccretion: foreignRankAccretions){
                sendBuffer << foreignAccretion.targetRank;
                sendBuffer << foreignAccretion.internalSeeds.size();
                for(const auto &seedEntry : foreignAccretion.internalSeeds){
                    sendBuffer << seedEntry.first;
                    sendBuffer << seedEntry.second;
                }
            }
        }

        // Start communications
        dataCommunicator->discoverRecvs();
        dataCommunicator->startAllRecvs();
        dataCommunicator->startAllSends();

        // Receive the accretiation to grow on behalf of neighbour processes
        int nCompletedRecvs = 0;
        while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
            int senderRank = dataCommunicator->waitAnyRecv();
            RecvBuffer &recvBuffer = dataCommunicator->getRecvBuffer(senderRank);

            std::size_t nForeignAccretions;
            recvBuffer >> nForeignAccretions;

            for (std::size_t k = 0; k < nForeignAccretions; ++k) {
                // Target rank
                int targetRank;
                recvBuffer >> targetRank;

                // Get the accretion to update
                //
                // If an accretions with the requeste owner and target
                // ranks doesn't exist create a new one.
                auto accretionsItr = accretions->begin();
                for(; accretionsItr != accretions->end(); ++accretionsItr){
                    if (accretionsItr->targetRank!= targetRank) {
                        continue;
                    }

                    break;
                }

                if (accretionsItr == accretions->end()) {
                    accretions->emplace_back();
                    accretionsItr = accretions->begin() + accretions->size() - 1;
                    accretionsItr->targetRank = targetRank;
                }

                AccretionData &accretion = *accretionsItr;

                // Initialize accretion seeds and population
                //
                // We are receiving only internal octants, there is no need
                // to explicitly check if and octant is internal before add
                // it to the population and to the internal seeds.
                std::size_t nSeeds;
                recvBuffer >> nSeeds;

                for(std::size_t n = 0; n < nSeeds; ++n){
                    uint64_t globalIdx;
                    recvBuffer >> globalIdx;

                    int layer;
                    recvBuffer >> layer;

                    assert(isInternal(globalIdx));
                    accretion.population.insert({globalIdx, layer});
                    accretion.internalSeeds.insert({globalIdx, layer});
                }
            }

            ++nCompletedRecvs;
        }

        // Wait until all exchanges are completed
        dataCommunicator->waitAllSends();
    }

    /*! Build ghost octants.
     *
     *  \param[in] bordersPerProc Information on partition boundary octants.
     *  \param[in] accretions Accretions used to growth the sources.
     */
    void
    ParaTree::buildGhostOctants(const std::map<int, u32vector> &bordersPerProc,
                                const std::vector<AccretionData> &accretions){

        DataCommunicator ghostDataCommunicator(m_comm);

        // Binary size of a ghost entry in the communication buffer
        const std::size_t GHOST_ENTRY_BINARY_SIZE = sizeof(uint64_t) + Octant::getBinarySize() + sizeof(int);

        // Fill the send buffers with source octants
        //
        // A source octant is an internal octants that is ghosts on other
        // process.
        for(const auto &bordersPerProcEntry : bordersPerProc){
            int rank = bordersPerProcEntry.first;
            const std::vector<uint32_t> &rankBordersPerProc = bordersPerProcEntry.second;
            std::size_t nRankBordersPerProc = rankBordersPerProc.size();

            // Get the accretion associated with this rank
            //
            // If there are no accretions it means we are building only the
            // first layer of ghosts.
            std::vector<AccretionData>::const_iterator accretionsItr;
            if (accretions.size() > 0) {
                bool accretionFound = false;
                BITPIT_UNUSED(accretionFound);
                accretionsItr = accretions.begin();
                for (; accretionsItr != accretions.end(); ++accretionsItr) {
                    if (accretionsItr->targetRank == rank) {
                        accretionFound = true;
                        break;
                    }
                }
                assert(accretionFound);
            } else {
                accretionsItr = accretions.end();
            }

            // Initialize the send
            std::size_t buffSize = GHOST_ENTRY_BINARY_SIZE * nRankBordersPerProc;
            ghostDataCommunicator.setSend(rank, buffSize);

            // Fill the buffer
            SendBuffer &sendBuffer = ghostDataCommunicator.getSendBuffer(rank);
            for(uint32_t sourceLocalIdx : rankBordersPerProc){
                // Global index
                uint64_t sourceGlobalIdx = getGlobalIdx(sourceLocalIdx);
                sendBuffer << sourceGlobalIdx;

                // Source data
                sendBuffer << m_octree.m_octants[sourceLocalIdx];

                // Layer information
                //
                // If not accretions are received we are building only the
                // first layer of ghosts.
                int layer;
                if (accretions.size() > 0) {
                    layer = accretionsItr->population.at(sourceGlobalIdx);
                } else {
                    layer = 0;
                }

                sendBuffer << layer;
            }
        }

        // Discover the receives
        ghostDataCommunicator.discoverRecvs();
        ghostDataCommunicator.startAllRecvs();
        ghostDataCommunicator.startAllSends();

        // Get the ranks from which ghosts will be received
        std::vector<int> ghostCommunicatorRecvsRanks = ghostDataCommunicator.getRecvRanks();
        std::sort(ghostCommunicatorRecvsRanks.begin(), ghostCommunicatorRecvsRanks.end());

        // Prepare ghost data structures
        m_octree.m_sizeGhosts = 0;
        for (int rank : ghostCommunicatorRecvsRanks) {
            RecvBuffer &recvBuffer = ghostDataCommunicator.getRecvBuffer(rank);
            std::size_t nRankGhosts = recvBuffer.getSize() / GHOST_ENTRY_BINARY_SIZE;
            m_octree.m_sizeGhosts += nRankGhosts;
        }

        m_octree.m_ghosts.resize(m_octree.m_sizeGhosts);
        m_octree.m_globalIdxGhosts.resize(m_octree.m_sizeGhosts);

        // Receive the ghosts
        //
        // Ghosts have to be received following the rank order.
        uint32_t ghostLocalIdx = 0;
        for (int rank : ghostCommunicatorRecvsRanks) {
            ghostDataCommunicator.waitRecv(rank);
            RecvBuffer &recvBuffer = ghostDataCommunicator.getRecvBuffer(rank);

            std::size_t nRankGhosts = recvBuffer.getSize() / GHOST_ENTRY_BINARY_SIZE;
            for(std::size_t n = 0; n < nRankGhosts; ++n){
                // Assign the global index
                uint64_t ghostGlobalIdx;
                recvBuffer >> ghostGlobalIdx;
                m_octree.m_globalIdxGhosts[ghostLocalIdx] = ghostGlobalIdx;

                // Build the ghosts
                //
                // The layer of the received octant will be overwritten with
                // the actual ghost layer.
                Octant &ghostOctant = m_octree.m_ghosts[ghostLocalIdx];
                recvBuffer >> ghostOctant;

                // Set the layer of the ghost
                int ghostLayer;
                recvBuffer >> ghostLayer;
                ghostOctant.setGhostLayer(ghostLayer);

                // Increase the ghost index
                ++ghostLocalIdx;
            }
        }

        // Wait for the communications to complete
        ghostDataCommunicator.waitAllSends();
    }

    /*! Communicate the marker of the octants and the auxiliary info[15].
     */
    void
    ParaTree::commMarker() {
        // If the tree is not partitioned, there is nothing to communicate.
        if (m_serial) {
            return;
        }

        // Binary size of a marker entry in the communication buffer
        const std::size_t MARKER_ENTRY_BINARY_SIZE = sizeof(int8_t) + sizeof(bool);

        // Fill communication buffer with level and marker
        //
        // It visits every element in m_bordersPerProc (one for every neighbor proc)
        // for every element it visits the border octants it contains and write them in the bitpit communication structure, DataCommunicator
        // this structure has a buffer for every proc containing the octants to be sent to that proc written in a char* buffer
        DataCommunicator markerCommunicator(m_comm);

        for(const auto &bordersPerProcEntry : m_bordersPerProc){
            int rank = bordersPerProcEntry.first;
            const std::vector<uint32_t> &rankBordersPerProc = m_bordersPerProc.at(rank);
            const std::size_t nRankBorders = rankBordersPerProc.size();

            std::size_t buffSize = nRankBorders * MARKER_ENTRY_BINARY_SIZE;
            markerCommunicator.setSend(rank, buffSize);

            SendBuffer sendBuffer = markerCommunicator.getSendBuffer(rank);
            for(std::size_t i = 0; i < nRankBorders; ++i){
                const Octant &octant = m_octree.m_octants[rankBordersPerProc[i]];
                sendBuffer << octant.getMarker();
                sendBuffer << octant.m_info[Octant::INFO_AUX];
            }
        }

        markerCommunicator.discoverRecvs();
        markerCommunicator.startAllRecvs();
        markerCommunicator.startAllSends();

        // Read level and marker from communication buffer
        //
        // every receive buffer is visited, and read octant by octant.
        // every ghost octant level and marker are updated
        std::vector<int> recvRanks = markerCommunicator.getRecvRanks();
        std::sort(recvRanks.begin(), recvRanks.end());

        uint32_t ghostIdx = 0;
        for(int rank : recvRanks){
            markerCommunicator.waitRecv(rank);
            RecvBuffer &recvBuffer = markerCommunicator.getRecvBuffer(rank);

            const std::size_t nRankGhosts = recvBuffer.getSize() / MARKER_ENTRY_BINARY_SIZE;
            for(std::size_t i = 0; i < nRankGhosts; ++i){
                int8_t marker;
                recvBuffer >> marker;
                m_octree.m_ghosts[ghostIdx].setMarker(marker);

                bool aux;
                recvBuffer >> aux;
                m_octree.m_ghosts[ghostIdx].m_info[Octant::INFO_AUX] = aux;

                ++ghostIdx;
            }
        }

        markerCommunicator.waitAllSends();
    }
#endif

    /*! Update the distributed octree over the processes after a coarsening procedure
     * and track the change in a mapper.
     */
    void
    ParaTree::updateAfterCoarse(){

#if BITPIT_ENABLE_MPI==1
        if(m_serial){
#endif
            updateAdapt();
#if BITPIT_ENABLE_MPI==1
        }
        else{

            updateAdapt();

            //update partition_range_position
            uint64_t lastDescMorton = m_octree.getLastDescMorton();
            m_errorFlag = MPI_Allgather(&lastDescMorton,1,MPI_UINT64_T,m_partitionLastDesc.data(),1,MPI_UINT64_T,m_comm);
            uint64_t firstDescMorton = m_octree.getFirstDescMorton();
            m_errorFlag = MPI_Allgather(&firstDescMorton,1,MPI_UINT64_T,m_partitionFirstDesc.data(),1,MPI_UINT64_T,m_comm);

            //correct first and last desc morton for empty partitions
            if (m_nproc>1){

                //Last desc
                //Attention rank = 0 can't be empty
                for (int p = 1; p < m_nproc; ++p){
                    if (m_partitionRangeGlobalIdx[p] == m_partitionRangeGlobalIdx[p-1]){
                        m_partitionLastDesc[p] = m_partitionLastDesc[p-1];
                        if (m_rank == p){
                            m_octree.m_lastDescMorton = m_partitionLastDesc[p-1];
                        }
                    }
                }

                //first desc
                //attention! Last desc of last rank (if empty partition) is the maximum uint64_t
                int pp = m_nproc-1;
                if (m_partitionRangeGlobalIdx[pp] == m_partitionRangeGlobalIdx[pp-1]){
                    m_partitionFirstDesc[pp] = std::numeric_limits<uint64_t>::max();
                    if (m_rank == pp){
                        m_octree.m_firstDescMorton = std::numeric_limits<uint64_t>::max();
                    }
                }
                for (int p = pp-1; p > 0; --p){
                    if (m_partitionRangeGlobalIdx[p] == m_partitionRangeGlobalIdx[p-1]){
                        m_partitionFirstDesc[p] = m_partitionFirstDesc[p+1];
                        if (m_rank == p){
                            m_octree.m_firstDescMorton = m_partitionFirstDesc[p+1];
                        }
                    }
                }
            }//end if nprocs>1
        }
#endif
    }

    /*!Balance 2:1 the octree.
     * \param[in] verbose If set to true output messages will be printed on
     * the logger
     * \param[in] balanceNewOctants If set to true also new octants will be
     * balanced
     */
    void
    ParaTree::balance21(bool verbose, bool balanceNewOctants){

        // Print header
        if (verbose){
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " 2:1 BALANCE (balancing Marker before Adapt)" << endl;
            (*m_log) << " " << endl;
            (*m_log) << " Iterative procedure	" << endl;
            (*m_log) << " " << endl;
        }

        // 2:1 balancing
        int iteration = 0;
        bool markersModified = true;
        while (markersModified) {
            if (verbose){
                (*m_log) << " Iteration	:	" + to_string(iteration) << endl;
            }

            // Only first iteration will process internl octants
            bool processInternals = (iteration == 0);

#if BITPIT_ENABLE_MPI==1
            // Communicate markers
            commMarker();
#endif

            // Pre-processing for 2:1 balancing
            m_octree.preBalance21(processInternals);

#if BITPIT_ENABLE_MPI==1
            // Communicate markers
            commMarker();
#endif

            // Execute local loadbalance
            markersModified = m_octree.localBalance(balanceNewOctants, processInternals);
#if BITPIT_ENABLE_MPI==1
            if (!m_serial) {
                MPI_Allreduce(MPI_IN_PLACE, &markersModified, 1, MPI_C_BOOL, MPI_LOR, m_comm);
            }
#endif

            // Increase iteration counter
            iteration++;
        }

        // Post-processing for 2:1 balancing
        //
        // The post-process is done calling the pre-process function on the
        // ghost octants.
        m_octree.preBalance21(false);

#if BITPIT_ENABLE_MPI==1
        // Communicate markers
        commMarker();
#endif

        // Print footer
        if (verbose){
            (*m_log) << " 2:1 Balancing reached " << endl;
            (*m_log) << " " << endl;
            (*m_log) << "---------------------------------------------" << endl;
        }
    };

    // =================================================================================== //
    // TESTING OUTPUT METHODS												    			   //
    // =================================================================================== //

    /** Write the physical octree mesh in .vtu format in a user-defined file.
     * If the connectivity is not stored, the method temporary computes it.
     * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
     * \param[in] filename Name of output file (PABLO will add the total number of processes p000# and the current rank s000#).
     */
    void
    ParaTree::write(const std::string &filename) {

        if (m_octree.m_connectivity.size() == 0) {
            m_octree.computeConnectivity();
        }

        stringstream name;
        name << "s" << std::setfill('0') << std::setw(4) << m_nproc << "-p" << std::setfill('0') << std::setw(4) << m_rank << "-" << filename << ".vtu";

        ofstream out(name.str().c_str());
        if(!out.is_open()){
            stringstream ss;
            ss << filename << "*.vtu cannot be opened and it won't be written." << endl;
            (*m_log) << ss.str();
            return;
        }
        int nofNodes = m_octree.m_nodes.size();
        int nofOctants = m_octree.m_connectivity.size();
        int nofGhosts = m_octree.m_ghostsConnectivity.size();
        int nofAll = nofGhosts + nofOctants;
        out << "<?xml version=\"1.0\"?>" << endl
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
            << "  <UnstructuredGrid>" << endl
            << "    <Piece NumberOfCells=\"" << m_octree.m_connectivity.size() + m_octree.m_ghostsConnectivity.size() << "\" NumberOfPoints=\"" << m_octree.m_nodes.size() << "\">" << endl;
        out << "      <Points>" << endl
            << "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        for(int i = 0; i < nofNodes; i++)
            {
                for(int j = 0; j < 3; ++j){
                    if (j==0) out << std::setprecision(6) << m_trans.mapX(m_octree.m_nodes[i][j]) << " ";
                    if (j==1) out << std::setprecision(6) << m_trans.mapY(m_octree.m_nodes[i][j]) << " ";
                    if (j==2) out << std::setprecision(6) << m_trans.mapZ(m_octree.m_nodes[i][j]) << " ";
                }
                if((i+1)%4==0 && i!=nofNodes-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Points>" << endl
            << "      <Cells>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofOctants; i++)
            {
                for(int j = 0; j < m_treeConstants->nNodes; j++)
                    {
                        int jj = j;
                        if (m_dim==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << m_octree.m_connectivity[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofOctants-1)
                    out << endl << "          ";
            }
        for(int i = 0; i < nofGhosts; i++)
            {
                for(int j = 0; j < m_treeConstants->nNodes; j++)
                    {
                        int jj = j;
                        if (m_dim==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << m_octree.m_ghostsConnectivity[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofGhosts-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                out << (i+1)*m_treeConstants->nNodes << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                int type;
                type = 5 + (m_dim*2);
                out << type << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Cells>" << endl
            << "    </Piece>" << endl
            << "  </UnstructuredGrid>" << endl
            << "</VTKFile>" << endl;


        if(m_rank == 0){
            name.str("");
            name << "s" << std::setfill('0') << std::setw(4) << m_nproc << "-" << filename << ".pvtu";
            ofstream pout(name.str().c_str());
            if(!pout.is_open()){
                stringstream ss;
                ss << filename << "*.pvtu cannot be opened and it won't be written." << endl;
                (*m_log) << ss.str();
                return;
            }

            pout << "<?xml version=\"1.0\"?>" << endl
                 << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
                 << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
                 << "    <PPointData>" << endl
                 << "    </PPointData>" << endl
                 << "    <PCellData Scalars=\"\">" << endl;
            pout << "    </PCellData>" << endl
                 << "    <PPoints>" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
                 << "    </PPoints>" << endl;
            for(int i = 0; i < m_nproc; i++)
                pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << m_nproc << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
            pout << "  </PUnstructuredGrid>" << endl
                 << "</VTKFile>";

            pout.close();

        }
#if BITPIT_ENABLE_MPI==1
        if (isCommSet()) {
            MPI_Barrier(m_comm);
        }
#endif

    }

    /** Write the physical octree mesh in .vtu format with data for test in a user-defined file.
     * If the connectivity is not stored, the method temporary computes it.
     * The method doesn't write the ghosts on file.
     * \param[in] filename Name of output file (PABLO will add the total number of processes p000# and the current rank s000#).
     * \param[in] data Vector of double with user data.
     */
    void
    ParaTree::writeTest(const std::string &filename, vector<double> data) {

        if (m_octree.m_connectivity.size() == 0) {
            m_octree.computeConnectivity();
        }

        stringstream name;
        name << "s" << std::setfill('0') << std::setw(4) << m_nproc << "-p" << std::setfill('0') << std::setw(4) << m_rank << "-" << filename << ".vtu";

        ofstream out(name.str().c_str());
        if(!out.is_open()){
            stringstream ss;
            ss << filename << "*.vtu cannot be opened and it won't be written.";
            (*m_log) << ss.str();
            return;
        }
        int nofNodes = m_octree.m_nodes.size();
        int nofOctants = m_octree.m_connectivity.size();
        int nofAll = nofOctants;
        out << "<?xml version=\"1.0\"?>" << endl
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
            << "  <UnstructuredGrid>" << endl
            << "    <Piece NumberOfCells=\"" << m_octree.m_connectivity.size() << "\" NumberOfPoints=\"" << m_octree.m_nodes.size() << "\">" << endl;
        out << "      <CellData Scalars=\"Data\">" << endl;
        out << "      <DataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        int ndata = m_octree.m_connectivity.size();
        for(int i = 0; i < ndata; i++)
            {
                out << std::setprecision(6) << data[i] << " ";
                if((i+1)%4==0 && i!=ndata-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </CellData>" << endl
            << "      <Points>" << endl
            << "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        for(int i = 0; i < nofNodes; i++)
            {
                for(int j = 0; j < 3; ++j){
                    if (j==0) out << std::setprecision(6) << m_trans.mapX(m_octree.m_nodes[i][j]) << " ";
                    if (j==1) out << std::setprecision(6) << m_trans.mapY(m_octree.m_nodes[i][j]) << " ";
                    if (j==2) out << std::setprecision(6) << m_trans.mapZ(m_octree.m_nodes[i][j]) << " ";
                }
                if((i+1)%4==0 && i!=nofNodes-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Points>" << endl
            << "      <Cells>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofOctants; i++)
            {
                for(int j = 0; j < m_treeConstants->nNodes; j++)
                    {
                        int jj = j;
                        if (m_dim==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << m_octree.m_connectivity[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofOctants-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                out << (i+1)*m_treeConstants->nNodes << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                int type;
                type = 5 + (m_dim*2);
                out << type << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Cells>" << endl
            << "    </Piece>" << endl
            << "  </UnstructuredGrid>" << endl
            << "</VTKFile>" << endl;


        if(m_rank == 0){
            name.str("");
            name << "s" << std::setfill('0') << std::setw(4) << m_nproc << "-" << filename << ".pvtu";
            ofstream pout(name.str().c_str());
            if(!pout.is_open()){
                stringstream ss;
                ss << filename << "*.pvtu cannot be opened and it won't be written." << endl;
                (*m_log) << ss.str();
                return;
            }

            pout << "<?xml version=\"1.0\"?>" << endl
                 << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
                 << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
                 << "    <PPointData>" << endl
                 << "    </PPointData>" << endl
                 << "    <PCellData Scalars=\"Data\">" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\"/>" << endl
                 << "    </PCellData>" << endl
                 << "    <PPoints>" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
                 << "    </PPoints>" << endl;
            for(int i = 0; i < m_nproc; i++)
                pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << m_nproc << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
            pout << "  </PUnstructuredGrid>" << endl
                 << "</VTKFile>";

            pout.close();

        }
#if BITPIT_ENABLE_MPI==1
        if (isCommSet()) {
            MPI_Barrier(m_comm);
        }
#endif

    }

    // =============================================================================== //

}

