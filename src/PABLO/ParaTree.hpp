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

#ifndef __BITPIT_PABLO_PARA_TREE_HPP__
#define __BITPIT_PABLO_PARA_TREE_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#include "DataLBInterface.hpp"
#include "DataCommInterface.hpp"
#include "communications.hpp"
#endif
#include "Global.hpp"
#include "Octant.hpp"
#include "LocalTree.hpp"
#include "Map.hpp"
#include "bitpit_IO.hpp"
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <bitset>
#include <algorithm>

namespace bitpit {

    // =================================================================================== //
    // TYPEDEFS																			   //
    // =================================================================================== //
    typedef std::vector<bool>				bvector;
    typedef std::vector<int>				ivector;
    typedef std::bitset<72>					octantID;
    typedef std::vector<Octant*>			ptroctvector;
    typedef ptroctvector::iterator			octantIterator;

    // =================================================================================== //
    // CLASS DEFINITION                                                                    //
    // =================================================================================== //

    /*!
     *	\ingroup		PABLO
     *	\date			17/dec/2015
     *	\authors		Marco Cisternino
     *	\authors		Edoardo Lombardi
     *
     *	\brief Para Tree is the user interface class
     *
     *	The user should (read can...) work only with this class and its methods.
     *	The sizes are intended in reference physical domain with limits [0,1]. The transformation from the logical
     *	domain to the physical domain is defined by an internal mapping.
     *
     *	The partition of the octree is performed by following the Z-curve defined by the Morton
     *	index of the octants. By default it is a balanced partition over the number of octants for each
     *	process.
     *
     *	Class ParaTree has a dimensional parameter int dim and it accepts only two
     *	 values: dim=2 and dim=3, for 2D and 3D respectively.
     */
    class ParaTree{

        // =================================================================================== //
        // MEMBERS																			   //
        // =================================================================================== //
    public:
        static const std::string	DEFAULT_LOG_FILE;			/**<Default name of logger file.*/

        typedef std::unordered_map<int, std::array<uint32_t, 2>> ExchangeRanges;

        enum Operation {
            OP_NONE,
            OP_INIT,
            OP_PRE_ADAPT,
            OP_ADAPT_MAPPED,
            OP_ADAPT_UNMAPPED,
            OP_LOADBALANCE_FIRST,
            OP_LOADBALANCE
        };

        struct LoadBalanceRanges {
            enum ExchangeAction {
                ACTION_UNDEFINED = -1,
                ACTION_NONE,
                ACTION_SEND,
                ACTION_DELETE,
                ACTION_RECEIVE
            };

            ExchangeAction sendAction;
            ExchangeRanges sendRanges;

            ExchangeAction recvAction;
            ExchangeRanges recvRanges;

            LoadBalanceRanges();
            LoadBalanceRanges(bool serial, const ExchangeRanges &_sendRanges, const ExchangeRanges &_recvRanges);

            void clear();
        };

    private:
        typedef std::unordered_map<int, std::array<uint64_t, 2>> PartitionIntersections;

        struct AccretionData {
            int targetRank;
            std::unordered_map<uint64_t, int> seeds;
            std::unordered_map<uint64_t, int> population;
        };

        //undistributed members
        std::vector<uint64_t>	m_partitionFirstDesc; 			/**<Global array containing position of the first possible octant in each processor*/
        std::vector<uint64_t>	m_partitionLastDesc; 			/**<Global array containing position of the last possible octant in each processor*/
        std::vector<uint64_t>	m_partitionRangeGlobalIdx;	 	/**<Global array containing global index of the last existing octant in each processor*/
        std::vector<uint64_t>	m_partitionRangeGlobalIdx0;	 	/**<Global array containing global index of the last existing octant in each processor before the last loadBalance (after an adapt is set equal to the actual.)*/
        uint64_t 				m_globalNumOctants;   			/**<Global number of octants in the parallel octree*/
        int 					m_nproc;						/**<Number of processes of the job*/
        uint8_t 				m_maxDepth;						/**<Global max existing level in the parallel octree*/
        Global					m_global;						/**<Global variables*/
        std::size_t 			m_nofGhostLayers;				/**<Global number of ghost layers from the process boundary expressing the depth of the ghost halo*/

        //distributed members
        int 					m_rank;							/**<Local m_rank of process*/
        LocalTree 				m_octree;						/**<Local tree in each processor*/
        std::map<int,u32vector> m_bordersPerProc;				/**<Local indices of border octants per process*/
        ptroctvector 			m_internals;					/**<Local pointers to internal octants*/
        ptroctvector 			m_pborders;						/**<Local pointers to border of process octants*/

        //distributed adpapting memebrs
        u32vector 				m_mapIdx;						/**<Local mapper for adapting. Mapper from new octants to old octants.
                                                                   m_mapIdx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
                                                                   if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
                                                                   if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
                                                                */
        //elements sent during last loadbalance operation
        LoadBalanceRanges         m_loadBalanceRanges;											  /**<Local mapper for sent elements. Each element refers to the receiver rank and collect the */

        //auxiliary members
        int 					m_errorFlag;					/**<MPI error flag*/
        bool 					m_serial;						/**<True if the octree is the same on each processor, False if the octree is distributed*/
        double					m_tol;							/**<Tolerance for geometric operations.*/

        //map members
        Map 					m_trans;						/**<Transformation map from m_logical to physical domain*/
        uint8_t					m_dim;							/**<Space dimension of the octree object (2D/3D).*/

        //boundary conditions members
        bvector 				m_periodic;						/**<Boolvector: i-th element is true if the i-th boundary face is a periodic interface.*/

        //info member
        uint64_t				m_status;						/**<Label of actual m_status of octree (incremental after an adpat
                                                                   with at least one modifyed element).*/
        Operation				m_lastOp;						/**<Last operation perforfmed by the octree (initialization, adapt (mapped or unmapped), loadbalance (first or not).*/

        //log member
        Logger* 				m_log;							/**<Log object pointer*/

        //communicator
#if BITPIT_ENABLE_MPI==1
        //TODO Duplicate communicator
        MPI_Comm 				m_comm;							/**<MPI communicator*/
#endif

        // =================================================================================== //
        // CONSTRUCTORS AND OPERATORS														   //
        // =================================================================================== //
    public:
#if BITPIT_ENABLE_MPI==1
        ParaTree(std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        ParaTree(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        ParaTree(std::istream &stream, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
#else
        ParaTree(std::string logfile = DEFAULT_LOG_FILE);
        ParaTree(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE);
        ParaTree(std::istream &stream, std::string logfile = DEFAULT_LOG_FILE);
#endif
        ~ParaTree();

        ParaTree(const ParaTree & other);

    private:
        ParaTree & operator=(const ParaTree & other) = delete;

        // =================================================================================== //
        // METHODS																			   //
        // =================================================================================== //
    public:
        virtual void	reset();

        virtual int		getDumpVersion() const;
        virtual void	dump(std::ostream &stream, bool full = true);
        virtual void	restore(std::istream &stream);

        void	printHeader();
        void	printFooter();

        // =================================================================================== //
        // BASIC GET/SET METHODS															   //
        // =================================================================================== //
        uint8_t 	getDim() const;
        uint64_t 	getGlobalNumOctants() const;
        bool		getSerial() const;
        bool		getParallel() const;
        int 		getRank() const;
        int 		getNproc() const;
        Logger& 	getLog();
        Operation	getLastOperation() const;
#if BITPIT_ENABLE_MPI==1
        void		setComm(MPI_Comm communicator);
        void		replaceComm(MPI_Comm communicator, MPI_Comm *previousCommunicator = nullptr);
        void		freeComm();
        MPI_Comm	getComm() const;
        bool		isCommSet() const;
#endif
        const std::vector<uint64_t> &getPartitionRangeGlobalIdx() const;
        const std::vector<uint64_t> &getPartitionFirstDesc() const;
        const std::vector<uint64_t> &getPartitionLastDesc() const;
        darray3		getOrigin() const;
        double		getX0() const;
        double		getY0() const;
        double		getZ0() const;
        double		getL() const;
        int 		getMaxLevel() const;
        uint32_t 	getMaxLength() const;
        uint8_t 	getNnodes() const;
        uint8_t 	getNfaces() const;
        uint8_t 	getNedges() const;
        uint8_t 	getNchildren() const;
        uint8_t 	getNnodesperface() const;
        void 		getNormals(int8_t normals[6][3]) const;
        void 		getOppface(uint8_t oppface[6]) const;
        void 		getFacenode(uint8_t facenode[6][4]) const;
        void 		getNodeface(uint8_t nodeface[8][3]) const;
        void 		getEdgeface(uint8_t edgeface[12][2]) const;
        void 		getNodecoeffs(int8_t nodecoeffs[8][3]) const;
        void 		getEdgecoeffs(int8_t edgecoeffs[12][3]) const;
        int8_t		(*getNormals())[3];
        uint8_t		(*getOppface());
        uint8_t 	(*getFacenode())[4];
        uint8_t 	(*getNodeface())[3];
        uint8_t 	(*getEdgeface())[2];
        int8_t 		(*getNodecoeffs())[3];
        int8_t 		(*getEdgecoeffs())[3];
        bvector		getPeriodic() const;
        double		getTol() const;
        bool		getPeriodic(uint8_t i) const;
        void		setPeriodic(uint8_t i);
        void		setTol(double tol = 1.0e-14);

        // =================================================================================== //
        // INDEX BASED METHODS																   //
        // =================================================================================== //
        darray3 	getCoordinates(uint32_t idx) const;
        double 		getX(uint32_t idx) const;
        double 		getY(uint32_t idx) const;
        double 		getZ(uint32_t idx) const;
        double 		getSize(uint32_t idx) const;
        double 		getArea(uint32_t idx) const;
        double 		getVolume(uint32_t idx) const;
        void 		getCenter(uint32_t idx, darray3& center) const;
        darray3 	getCenter(uint32_t idx) const;
        darray3 	getFaceCenter(uint32_t idx, uint8_t iface) const;
        void 		getFaceCenter(uint32_t idx, uint8_t iface, darray3& center) const;
        darray3 	getNode(uint32_t idx, uint8_t inode) const;
        void 		getNode(uint32_t idx, uint8_t inode, darray3& node) const;
        void 		getNodes(uint32_t idx, darr3vector & nodes) const;
        darr3vector getNodes(uint32_t idx) const;
        void 		getNormal(uint32_t idx, uint8_t iface, darray3 & normal) const;
        darray3 	getNormal(uint32_t idx, uint8_t iface) const;
        int8_t 		getMarker(uint32_t idx) const;
        uint8_t 	getLevel(uint32_t idx) const;
        uint64_t 	getMorton(uint32_t idx) const;
        uint64_t 	getNodeMorton(uint32_t idx, uint8_t inode) const;
        bool 		getBalance(uint32_t idx) const;
        bool		getBound(uint32_t idx, uint8_t iface) const;
        bool		getBound(uint32_t idx) const;
        bool		getPbound(uint32_t idx, uint8_t iface) const;
        bool		getPbound(uint32_t idx) const;
        bool 		getIsNewR(uint32_t idx) const;
        bool 		getIsNewC(uint32_t idx) const;
        uint64_t 	getGlobalIdx(uint32_t idx) const;
        uint64_t 	getGhostGlobalIdx(uint32_t idx) const;
        uint32_t    getLocalIdx(uint64_t gidx) const;
        uint32_t    getLocalIdx(uint64_t gidx,int rank) const;
        void        getLocalIdx(uint64_t gidx,uint32_t & lidx,int & rank) const;
        uint32_t    getGhostLocalIdx(uint64_t gidx) const;
        octantID	getPersistentIdx(uint32_t idx) const;
        int8_t 		getPreMarker(uint32_t idx);
        void 		setMarker(uint32_t idx, int8_t marker);
        void 		setBalance(uint32_t idx, bool balance);

        // =================================================================================== //
        // POINTER BASED METHODS															   //
        // =================================================================================== //
        darray3 	getCoordinates(const Octant* oct) const;
        double 		getX(const Octant* oct) const;
        double 		getY(const Octant* oct) const;
        double 		getZ(const Octant* oct) const;
        double 		getSize(const Octant* oct) const;
        double 		getArea(const Octant* oct) const;
        double 		getVolume(const Octant* oct) const;
        void 		getCenter(const Octant* oct, darray3& center) const;
        darray3 	getCenter(const Octant* oct) const;
        darray3 	getFaceCenter(const Octant* oct, uint8_t iface) const;
        void 		getFaceCenter(const Octant* oct, uint8_t iface, darray3& center) const;
        darray3 	getNode(const Octant* oct, uint8_t inode) const;
        void 		getNode(const Octant* oct, uint8_t inode, darray3& node) const;
        void 		getNodes(const Octant* oct, darr3vector & nodes) const;
        darr3vector getNodes(const Octant* oct) const;
        void 		getNormal(const Octant* oct, uint8_t iface, darray3 & normal) const;
        darray3 	getNormal(const Octant* oct, uint8_t iface) const;
        int8_t 		getMarker(const Octant* oct) const;
        uint8_t 	getLevel(const Octant* oct) const;
        uint64_t 	getMorton(const Octant* oct) const;
        uint64_t 	getNodeMorton(const Octant* oct, uint8_t inode) const;
        bool 		getBalance(const Octant* oct) const;
        bool		getBound(const Octant* oct, uint8_t iface) const;
        bool		getBound(const Octant* oct) const;
        bool		getPbound(const Octant* oct, uint8_t iface) const;
        bool		getPbound(const Octant* oct) const;
        bool 		getIsNewR(const Octant* oct) const;
        bool 		getIsNewC(const Octant* oct) const;
        uint32_t 	getIdx(const Octant* oct) const;
        uint64_t 	getGlobalIdx(const Octant* oct) const;
        octantID	getPersistentIdx(const Octant* oct) const;
        int8_t 		getPreMarker(Octant* oct);
        void 		setMarker(Octant* oct, int8_t marker);
        void 		setBalance(Octant* oct, bool balance);

        // =================================================================================== //
        // LOCAL TREE GET/SET METHODS														   //
        // =================================================================================== //
        uint64_t 	getStatus() const;
        uint32_t 	getNumOctants() const;
        uint32_t 	getNumGhosts() const;
        uint32_t 	getNumNodes() const;
        uint8_t 	getLocalMaxDepth() const;
        double	 	getLocalMaxSize() const;
        double	 	getLocalMinSize() const;
        uint8_t 	getBalanceCodimension() const;
        uint64_t 	getFirstDescMorton() const;
        uint64_t 	getLastDescMorton() const;
        uint64_t 	getLastDescMorton(uint32_t idx) const;
        octantIterator	getInternalOctantsBegin();
        octantIterator	getInternalOctantsEnd();
        octantIterator	getPboundOctantsBegin();
        octantIterator	getPboundOctantsEnd();
        void 		setBalanceCodimension(uint8_t b21codim);

        // =================================================================================== //
        // INTERSECTION GET/SET METHODS														   //
        // =================================================================================== //
        uint32_t 	getNumIntersections() const;
        Intersection* getIntersection(uint32_t idx);
        uint8_t 	getLevel(const Intersection* inter) const;
        bool 		getFiner(const Intersection* inter) const;
        bool 		getBound(const Intersection* inter) const;
        bool 		getIsGhost(const Intersection* inter) const;
        bool 		getPbound(const Intersection* inter) const;
        uint8_t 	getFace(const Intersection* inter) const;
        u32vector 	getOwners(const Intersection* inter) const;
        uint32_t 	getIn(const Intersection* inter) const;
        uint32_t 	getOut(const Intersection* inter) const;
        bool		getOutIsGhost(const Intersection* inter) const;
        double 		getSize(const Intersection* inter) const;
        double 		getArea(const Intersection* inter) const;
        darray3 	getCenter(const Intersection* inter) const;
        darr3vector getNodes(const Intersection* inter) const;
        darray3 	getNormal(const Intersection* inter) const;

        // =================================================================================== //
        // OTHER GET/SET METHODS															   //
        // =================================================================================== //
        Octant*	getOctant(uint32_t idx);
        const Octant* getOctant(uint32_t idx) const;
        Octant*	getGhostOctant(uint32_t idx);
        const Octant* getGhostOctant(uint32_t idx) const;
        bool getIsGhost(const Octant* oct) const;
        int getGhostLayer(const Octant* oct) const;
        const LoadBalanceRanges & getLoadBalanceRanges() const;
        std::size_t getNofGhostLayers() const;
        void setNofGhostLayers(std::size_t nofGhostLayers);
        const std::map<int, std::vector<uint32_t>> & getBordersPerProc() const;

        // =================================================================================== //
        // PRIVATE GET/SET METHODS															   //
        // =================================================================================== //
    private:
        void		setDim(uint8_t dim);
        void 		setFirstDescMorton();
        void 		setLastDescMorton();

        // =================================================================================== //
        // OTHER METHODS												    			       //
        // =================================================================================== //

#if BITPIT_ENABLE_MPI==1
        void	_initializeCommunications(MPI_Comm comm);
#endif
        void	_initializePartitions();
        void	_initialize(uint8_t dim, const std::string &logfile);

#if BITPIT_ENABLE_MPI==1
        void	initialize(const std::string &logfile, MPI_Comm comm);
        void	initialize(uint8_t dim, const std::string &logfile, MPI_Comm comm);
#else
        void	initialize(const std::string &logfile);
        void	initialize(uint8_t dim, const std::string &logfile);
#endif
        void	initializeLogger(const std::string &logfile);
        void	reinitialize(uint8_t dim, const std::string &logfile);

        void	reset(bool createRoot);

        // =================================================================================== //
        // OTHER OCTANT BASED METHODS												    	   //
        // =================================================================================== //

        void        findAllGlobalNeighbours(uint32_t idx, std::vector<uint64_t> &globalNeighs);
        void        findNeighbours(const Octant* oct, bool haveIidx, uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost, bool onlyinternals = false) const;
    public:
        void 		findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
        void 		findNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const ;
        void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours) const;
        void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
        void 		findAllNodeNeighbours(uint32_t idx, uint32_t inode, u32vector & neighbours, bvector & isghost);
        void 		findAllNodeNeighbours(const Octant* oct, uint32_t inode, u32vector & neighbours, bvector & isghost) const;
        Octant* 	getPointOwner(const dvector &point);
        Octant* 	getPointOwner(const dvector &point, bool & isghost);
        Octant* 	getPointOwner(const darray3 &point);
        Octant* 	getPointOwner(const darray3 &point, bool & isghost);
        uint32_t 	getPointOwnerIdx(const double * point) const;
        uint32_t 	getPointOwnerIdx(const double * point, bool & isghost) const;
        uint32_t 	getPointOwnerIdx(const dvector &point) const;
        uint32_t 	getPointOwnerIdx(const dvector &point, bool & isghost) const;
        uint32_t 	getPointOwnerIdx(const darray3 &point) const;
        uint32_t 	getPointOwnerIdx(const darray3 &point, bool & isghost) const;
        void 		findAllCodimensionNeighbours(uint32_t idx, u32vector & neighbours, bvector & isghost);
        void 		findAllCodimensionNeighbours(Octant* oct, u32vector & neighbours, bvector & isghost);
        void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost) const;
        void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost, ivector & rank) const;
        void 		getPreMapping(u32vector & idx, vector<int8_t> & markers, vector<bool> & isghost);
        bool 		isNodeOnOctant(const Octant* nodeOctant, uint8_t nodeIndex, const Octant* octant) const;
        bool 		isEdgeOnOctant(const Octant* edgeOctant, uint8_t edgeIndex, const Octant* octant) const;
        bool 		isFaceOnOctant(const Octant* faceOctant, uint8_t faceIndex, const Octant* octant) const;

        // =================================================================================== //
        // OTHER PARATREE BASED METHODS												    	   //
        // =================================================================================== //
        uint8_t		getMaxDepth() const;
        int 		findOwner(const uint64_t & morton) const;
        int 		getOwnerRank(const uint64_t & globalIdx) const;
        void        settleMarkers();
        void        preadapt();
        bool        checkToAdapt();
        bool        adapt(bool mapper_flag = false);
        bool 		adaptGlobalRefine(bool mapper_flag = false);
        bool 		adaptGlobalCoarse(bool mapper_flag = false);
        void 		computeConnectivity();
        void 		clearConnectivity();
        void 		updateConnectivity();
        const u32vector2D & getConnectivity() const;
        const u32vector & getConnectivity(uint32_t idx) const;
        const u32vector & getConnectivity(Octant* oct) const;
        const u32arr3vector & getNodes() const;
        const u32array3 & getNodeLogicalCoordinates(uint32_t inode) const;
        darray3 	getNodeCoordinates(uint32_t inode) const;
        const u32vector2D & getGhostConnectivity() const;
        const u32vector & getGhostConnectivity(uint32_t idx) const;
        const u32vector & getGhostConnectivity(const Octant* oct) const;
        bool        check21Balance();
#if BITPIT_ENABLE_MPI==1
        void 		loadBalance(dvector* weight = NULL);
        void 		loadBalance(uint8_t & level, dvector* weight = NULL);

        LoadBalanceRanges evalLoadBalanceRanges(dvector *weights);
        LoadBalanceRanges evalLoadBalanceRanges(uint8_t level, dvector *weights);
    private:
        void 		privateLoadBalance(uint32_t* partition);

        LoadBalanceRanges evalLoadBalanceRanges(const uint32_t *updatedPartition);

        ExchangeRanges evalLoadBalanceSendRanges(const uint32_t *updatedPartition);
        ExchangeRanges evalLoadBalanceRecvRanges(const uint32_t *updatedPartition);

        PartitionIntersections evalPartitionIntersections(const uint32_t *schema_A, int rank_A, const uint32_t *schema_B);
#endif
    public:
        double		levelToSize(uint8_t & level);

        // =================================================================================== //
        // OTHER INTERSECTION BASED METHODS										     		   //
        // =================================================================================== //
        void 		computeIntersections();

        // =================================================================================== //
        // OTHER PRIVATE METHODS												    		   //
        // =================================================================================== //
    private:
        Octant& extractOctant(uint32_t idx);
        bool 		private_adapt_mapidx(bool mapflag);
        void 		updateAdapt();
#if BITPIT_ENABLE_MPI==1
        void 		computePartition(uint32_t* partition);
        void 		computePartition(uint32_t* partition, dvector* weight);
        void 		computePartition(uint32_t* partition, uint8_t & level_, dvector* weight);
        void 		updateLoadBalance();
        void 		setPboundGhosts();
        void 		buildGhostOctants(const std::map<int, u32vector> &bordersPerProc, const std::vector<AccretionData> &accretions);

        void 		computeGhostHalo();
        void 		commMarker();
#endif
        void 		updateAfterCoarse();
        void 		balance21(bool verbose, bool balanceNewOctants);
        void		createPartitionInfo();

        bool		isInternal(uint64_t gidx) const;

        // =================================================================================== //
        // TESTING OUTPUT METHODS												    			   //
        // =================================================================================== //
    public:
        void 		write(std::string filename);
        void 		writeTest(std::string filename, dvector data);

        // =================================================================================== //
        // TEMPLATE METHODS												    			       //
        // =================================================================================== //
#if BITPIT_ENABLE_MPI==1

        /** Communicate data provided by the user between the processes.
         */
        template<class Impl>
        void
        communicate(DataCommInterface<Impl> & userData){
            //BUILD SEND BUFFERS
            DataCommunicator communicator(m_comm);
            size_t fixedDataSize = userData.fixedSize();
            std::map<int,u32vector >::iterator bitend = m_bordersPerProc.end();
            std::map<int,u32vector >::iterator bitbegin = m_bordersPerProc.begin();
            for(std::map<int,u32vector >::iterator bit = bitbegin; bit != bitend; ++bit){
                const int & key = bit->first;
                const u32vector & pborders = bit->second;
                size_t buffSize = 0;
                size_t nofPbordersPerProc = pborders.size();
                if(fixedDataSize != 0){
                    buffSize = fixedDataSize*nofPbordersPerProc;
                }
                else{
                    for(size_t i = 0; i < nofPbordersPerProc; ++i){
                        buffSize += userData.size(pborders[i]);
                    }
                }
                //enlarge buffer to store number of pborders from this proc
                buffSize += sizeof(size_t);
                //build buffer for this proc
                communicator.setSend(key,buffSize);
                SendBuffer & sendBuffer = communicator.getSendBuffer(key);
                //store number of pborders from this proc at the begining
                sendBuffer << nofPbordersPerProc;

                //WRITE SEND BUFFERS
                for(size_t j = 0; j < nofPbordersPerProc; ++j){
                    userData.gather(sendBuffer,pborders[j]);
                }
            }

            communicator.discoverRecvs();
            communicator.startAllRecvs();

            communicator.startAllSends();

            //READ RECEIVE BUFFERS
            int ghostOffset = 0;
            vector<int> recvRanks = communicator.getRecvRanks();
            std::sort(recvRanks.begin(),recvRanks.end());
            for(int rank : recvRanks){
                communicator.waitRecv(rank);
                RecvBuffer & recvBuffer = communicator.getRecvBuffer(rank);
                size_t nofGhostFromThisProc = 0;
                recvBuffer >> nofGhostFromThisProc;
                for(size_t k = 0; k < nofGhostFromThisProc; ++k){
                    userData.scatter(recvBuffer, k+ghostOffset);
                }
                ghostOffset += nofGhostFromThisProc;
            }
            communicator.waitAllSends();
        };

        /** Distribute Load-Balancing the octants (with user defined weights) of the whole tree and data provided by the user
         * over the processes of the job following the Morton order.
         * Until loadBalance is not called for the first time the mesh is serial.
         * Even distribute data provided by the user between the processes.
         * \param[in] userData User interface to distribute the data during loadBalance.
         * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
         */
        template<class Impl>
        void
        loadBalance(DataLBInterface<Impl> & userData, dvector* weight = NULL){
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

                privateLoadBalance(userData, partition);

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
                m_loadBalanceRanges.clear();

                (*m_log) << " " << endl;
                (*m_log) << " Serial partition : " << endl;
                (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
                (*m_log) << " " << endl;
                (*m_log) << "---------------------------------------------" << endl;
            }

        }

        /** Distribute Load-Balanced the octants (with user defined weights) of the whole tree and data provided by the user
         * over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
         * The families of octants of a desired level are retained compact on the same process.
         * Even distribute data provided by the user between the processes.
         * \param[in] userData User interface to distribute the data during loadBalance.
         * \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
         * \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
         */
        template<class Impl>
        void
        loadBalance(DataLBInterface<Impl> & userData, uint8_t & level, dvector* weight = NULL){

            //Write info on log
            (*m_log) << "---------------------------------------------" << endl;
            (*m_log) << " LOAD BALANCE " << endl;

            m_lastOp = OP_LOADBALANCE;
            if (m_nproc>1){

                uint32_t* partition = new uint32_t [m_nproc];
                computePartition(partition, level, weight);

                privateLoadBalance(userData, partition);

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
                m_loadBalanceRanges.clear();

                (*m_log) << " " << endl;
                (*m_log) << " Serial partition : " << endl;
                (*m_log) << " Octants for proc	"+ to_string(static_cast<unsigned long long>(0))+"	:	" + to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
                (*m_log) << " " << endl;
                (*m_log) << "---------------------------------------------" << endl;
            }

        }

        /**
        * Distribute Load-Balancing octants and user data of the whole
        * tree over the processes of the job following a given partition
        * distribution. Until loadBalance is not called for the first time
        * the mesh is serial.
        * \param[in] userData User data that will be distributed among the
        * processes.
        * \param[in] partition Target distribution of octants over processes.
        */
        template<class Impl>
        void
        privateLoadBalance(DataLBInterface<Impl> & userData,uint32_t* partition){

            if(m_serial)
            {
                m_lastOp = OP_LOADBALANCE_FIRST;
                (*m_log) << " " << endl;
                (*m_log) << " Initial Serial distribution : " << endl;
                for(int ii=0; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]+1)) << endl;
                }

                uint32_t stride = 0;
                for(int i = 0; i < m_rank; ++i)
                    stride += partition[i];
                LocalTree::octvector octantsCopy = m_octree.m_octants;
                LocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
                LocalTree::octvector::const_iterator last = first + partition[m_rank];

                m_octree.m_octants.assign(first, last);
                octvector(m_octree.m_octants).swap(m_octree.m_octants);

                first = octantsCopy.end();
                last = octantsCopy.end();

                userData.assign(stride,partition[m_rank]);

                //Update and build ghosts here
                updateLoadBalance();
                computeGhostHalo();
            }
            else
            {
                (*m_log) << " " << endl;
                (*m_log) << " Initial Parallel partition : " << endl;
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << endl;
                for(int ii=1; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << endl;
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
                else if(lh > (int64_t) m_octree.m_octants.size() - 1)
                    lh = m_octree.m_octants.size() - 1;

                if(m_rank == m_nproc - 1)
                    ft = m_octree.m_octants.size();
                else if(m_rank == 0)
                    ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] + 1);
                else{
                    ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1]);
                }
                if(ft > (int32_t)(m_octree.m_octants.size() - 1))
                    ft = m_octree.m_octants.size();
                else if(ft < 0)
                    ft = 0;

                //compute size Head and size Tail
                uint32_t headSize = (uint32_t)(lh + 1);
                uint32_t tailSize = (uint32_t)(m_octree.m_octants.size() - ft);
                uint32_t headOffset = headSize;
                uint32_t tailOffset = tailSize;

                //Communicator declaration
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
                MPI_Barrier(m_comm); //da spostare prima della prima comunicazione

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
                            //compute size of data in buffers
                            if(userData.fixedSize()){
                                buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
                            }
                            else{
                                for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                                    buffSize += userData.size(i);
                                }
                            }
                            //add room for uint32_t, number of octants in this buffer
                            buffSize += sizeof(uint32_t);
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            //store the number of octants at the beginning of the buffer
                            sendBuffer << nofElementsFromSuccessiveToPrevious;

                            for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                                sendBuffer << m_octree.m_octants[i];
                                userData.gather(sendBuffer,i);
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
                            //compute size of data in buffers
                            if(userData.fixedSize()){
                                buffSize +=  userData.fixedSize() * nofElementsFromSuccessiveToPrevious;
                            }
                            else{
                                for(int64_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
                                    buffSize += userData.size(i);
                                }
                            }
                            //add room for uint32_t, number of octants in this buffer
                            buffSize += sizeof(uint32_t);
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            //store the number of octants at the beginning of the buffer
                            sendBuffer << nofElementsFromSuccessiveToPrevious;

                            for(int64_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
                                //WRITE octants from lh - partition[p] to lh
                                sendBuffer << m_octree.m_octants[i];
                                userData.gather(sendBuffer,i);
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
                            //compute size of data in buffers
                            if(userData.fixedSize()){
                                buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
                            }
                            else{
                                for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
                                    buffSize += userData.size(i);
                                }
                            }
                            //add room for uint32_t, number of octants in this buffer
                            buffSize += sizeof(uint32_t);
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            //store the number of octants at the beginning of the buffer
                            sendBuffer << nofElementsFromPreviousToSuccessive;

                            for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
                                sendBuffer << m_octree.m_octants[i];
                                userData.gather(sendBuffer,i);
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
                            uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
                            std::size_t buffSize = (std::size_t)nofElementsFromPreviousToSuccessive * (std::size_t)Octant::getBinarySize();
                            //compute size of data in buffers
                            if(userData.fixedSize()){
                                buffSize +=  userData.fixedSize() * nofElementsFromPreviousToSuccessive;
                            }
                            else{
                                for(uint32_t i = ft; i <= endOctants; ++i){
                                    buffSize += userData.size(i);
                                }
                            }
                            //add room for uint32_t, number of octants in this buffer
                            buffSize += sizeof(uint32_t);
                            lbCommunicator.setSend(p,buffSize);
                            SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                            //store the number of octants at the beginning of the buffer
                            sendBuffer << nofElementsFromPreviousToSuccessive;

                            for(uint32_t i = ft; i <= endOctants; ++i ){
                                //WRITE octants from ft to ft + partition[p] -1
                                sendBuffer << m_octree.m_octants[i];
                                userData.gather(sendBuffer,i);
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
                lbCommunicator.startAllSends();

                uint32_t nofNewHead = 0;
                uint32_t nofNewTail = 0;

                //READ number of octants per sender
                std::map<int,uint32_t> nofNewOverProcs;
                int nCompletedRecvs = 0;
                while(nCompletedRecvs < lbCommunicator.getRecvCount()){
                    int rank = lbCommunicator.waitAnyRecv();
                    RecvBuffer & recvBuffer = lbCommunicator.getRecvBuffer(rank);
                    uint32_t nofNewPerProc;
                    recvBuffer >> nofNewPerProc;
                    nofNewOverProcs[rank] = nofNewPerProc;
                    if(rank < m_rank)
                        nofNewHead += nofNewPerProc;
                    else if(rank > m_rank)
                        nofNewTail += nofNewPerProc;
                    ++nCompletedRecvs;
                }

                lbCommunicator.waitAllSends();

                //MOVE RESIDENT TO BEGIN IN OCTANTS
                m_octree.m_sizeOctants = m_octree.m_octants.size();
                uint32_t resEnd = m_octree.m_sizeOctants - tailOffset;
                uint32_t nofResidents = resEnd - headOffset;
                uint32_t octCounter = 0;
                for(uint32_t i = headOffset; i < resEnd; ++i){
                    m_octree.m_octants[octCounter] = m_octree.m_octants[i];
                    userData.move(i,octCounter);
                    ++octCounter;
                }
                uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
                m_octree.m_octants.resize(newCounter, Octant(m_dim));
                userData.resize(newCounter);
                //MOVE RESIDENTS IN RIGHT POSITION
                uint32_t resCounter = nofNewHead + nofResidents - 1;
                for(uint32_t k = 0; k < nofResidents ; ++k){
                    m_octree.m_octants[resCounter - k] = m_octree.m_octants[nofResidents - k - 1];
                    userData.move(nofResidents - k - 1,resCounter - k);
                }

                //READ BUFFERS AND BUILD NEW OCTANTS
                newCounter = 0;
                bool jumpResident = false;

                vector<int> recvRanks = lbCommunicator.getRecvRanks();
                std::sort(recvRanks.begin(),recvRanks.end());
                for(int rank : recvRanks){
                    RecvBuffer & recvBuffer = lbCommunicator.getRecvBuffer(rank);
                    uint32_t nofNewPerProc = nofNewOverProcs[rank];
                    if(rank > m_rank && !jumpResident){
                        newCounter += nofResidents ;
                        jumpResident = true;
                    }
                    for(int i = nofNewPerProc - 1; i >= 0; --i){
                        recvBuffer >> m_octree.m_octants[newCounter];
                        userData.scatter(recvBuffer,newCounter);
                        ++newCounter;
                    }
                }
                octvector(m_octree.m_octants).swap(m_octree.m_octants);

                userData.shrink();

                delete [] newPartitionRangeGlobalidx; newPartitionRangeGlobalidx = NULL;

                //Update and ghosts here
                updateLoadBalance();
                computeGhostHalo();
                uint32_t nofGhosts = getNumGhosts();
                userData.resizeGhost(nofGhosts);

            }

            // Update load balance ranges
            std::unordered_map<int, std::array<uint32_t, 2>> sendRanges = evalLoadBalanceSendRanges(partition);
            std::unordered_map<int, std::array<uint32_t, 2>> recvRanges = evalLoadBalanceRecvRanges(partition);

            m_loadBalanceRanges = LoadBalanceRanges(m_serial, sendRanges, recvRanges);
        };
#endif

        // =============================================================================== //

    };

}

#endif /* __BITPIT_PARA_TREE_HPP__ */
