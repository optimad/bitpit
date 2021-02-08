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

#ifndef __BITPIT_PABLO_PARA_TREE_HPP__
#define __BITPIT_PABLO_PARA_TREE_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#include "DataLBInterface.hpp"
#include "DataCommInterface.hpp"
#include "bitpit_communications.hpp"
#endif
#include "tree_constants.hpp"
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
            std::unordered_map<uint64_t, int> internalSeeds;
            std::unordered_map<uint64_t, int> foreignSeeds;
            std::unordered_map<uint64_t, int> population;
        };

        //undistributed members
        std::vector<uint64_t>	m_partitionFirstDesc; 			/**<Global array containing position of the first possible octant in each process*/
        std::vector<uint64_t>	m_partitionLastDesc; 			/**<Global array containing position of the last possible octant in each process*/
        std::vector<uint64_t>	m_partitionRangeGlobalIdx;	 	/**<Global array containing global index of the last existing octant in each process*/
        std::vector<uint64_t>	m_partitionRangeGlobalIdx0;	 	/**<Global array containing global index of the last existing octant in each process before the last loadBalance (after an adapt is set equal to the actual.)*/
        uint64_t 				m_globalNumOctants;   			/**<Global number of octants in the parallel octree*/
        int 					m_nproc;						/**<Number of processes of the job*/
        int8_t 					m_maxDepth;						/**<Global max existing level in the parallel octree*/
        const TreeConstants	   *m_treeConstants;				/**<Tree constants*/
        std::size_t 			m_nofGhostLayers;				/**<Global number of ghost layers from the process boundary expressing the depth of the ghost halo*/

        //distributed members
        int 					m_rank;							/**<Local m_rank of process*/
        LocalTree 				m_octree;						/**<Local tree in each process*/
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
        bool 					m_serial;						/**<True if the octree is the same on each process, False if the octree is distributed*/
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
        ParaTree(const std::string &logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        ParaTree(uint8_t dim, const std::string &logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        ParaTree(std::istream &stream, const std::string &logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
#else
        ParaTree(const std::string &logfile = DEFAULT_LOG_FILE);
        ParaTree(uint8_t dim, const std::string &logfile = DEFAULT_LOG_FILE);
        ParaTree(std::istream &stream, const std::string &logfile = DEFAULT_LOG_FILE);
#endif
        virtual ~ParaTree();

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
        const int8_t	(*getNormals() const)[3];
        const uint8_t	*getOppface() const;
        const uint8_t	(*getFacenode() const)[4];
        const uint8_t	(*getNodeface() const)[3];
        const uint8_t	(*getEdgeface() const)[2];
        const int8_t	(*getNodecoeffs() const)[3];
        const int8_t	(*getEdgecoeffs() const)[3];
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
        uint64_t 	getLastDescMorton(const Octant* oct) const;
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
#if BITPIT_ENABLE_MPI==1
        void 		updateGlobalFirstDescMorton();
        void 		updateGlobalLasttDescMorton();
#endif

        // =================================================================================== //
        // OTHER METHODS												    			       //
        // =================================================================================== //

#if BITPIT_ENABLE_MPI==1
        void	_initializeCommunicator(MPI_Comm comm);
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
        void        findNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost, bool onlyinternals) const;
    public:
        void 		findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
        void 		findNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const ;
        void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours) const;
        void 		findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
        void 		findGhostNeighbours(const Octant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, bvector & isghost) const;
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
        void 		findGhostAllCodimensionNeighbours(uint32_t idx, u32vector & neighbours, bvector & isghost);
        void 		findGhostAllCodimensionNeighbours(Octant* oct, u32vector & neighbours, bvector & isghost);
        void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost) const;
        void 		getMapping(uint32_t & idx, u32vector & mapper, bvector & isghost, ivector & rank) const;
        void 		getPreMapping(u32vector & idx, std::vector<int8_t> & markers, std::vector<bool> & isghost);
        bool 		isNodeOnOctant(const Octant* nodeOctant, uint8_t nodeIndex, const Octant* octant) const;
        bool 		isEdgeOnOctant(const Octant* edgeOctant, uint8_t edgeIndex, const Octant* octant) const;
        bool 		isFaceOnOctant(const Octant* faceOctant, uint8_t faceIndex, const Octant* octant) const;
        int 		getPointOwnerRank(darray3 point);
        uint8_t		getFamilySplittingNode(const Octant*) const;
        void		expectedOctantAdapt(const Octant* octant, int8_t marker, octvector* result) const;

        // =================================================================================== //
        // OTHER PARATREE BASED METHODS												    	   //
        // =================================================================================== //
        int8_t		getMaxDepth() const;
        int 		findOwner(uint64_t morton) const;
        int 		getOwnerRank(uint64_t globalIdx) const;
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
        void 		computePartition(uint32_t *partition);
        void 		computePartition(const dvector *weight, uint32_t *partition);
        void 		computePartition(uint8_t level_, const dvector *weight, uint32_t *partition);
        void 		updateLoadBalance();
        void 		setPboundGhosts();
        void 		buildGhostOctants(const std::map<int, u32vector> &bordersPerProc, const std::vector<AccretionData> &accretions);
        void 		initializeGhostHaloAccretions(std::vector<AccretionData> *accretions);
        void 		growGhostHaloAccretions(std::unordered_map<uint32_t, std::vector<uint64_t>> *cachedOneRings, std::vector<AccretionData> *accretions);
        void 		exchangeGhostHaloAccretions(DataCommunicator *dataCommunicator, std::vector<AccretionData> *accretions);

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
        void 		write(const std::string &filename);
        void 		writeTest(const std::string &filename, dvector data);

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
                int  key = bit->first;
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
            std::vector<int> recvRanks = communicator.getRecvRanks();
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
            (*m_log) << "---------------------------------------------" << std::endl;
            (*m_log) << " LOAD BALANCE " << std::endl;

            m_lastOp = OP_LOADBALANCE;
            if (m_nproc>1){

                std::vector<uint32_t> partition(m_nproc);
                if (weight == NULL)
                    computePartition(partition.data());
                else
                    computePartition(weight, partition.data());

                weight = NULL;

                privateLoadBalance(partition.data(), &userData);

                //Write info of final partition on log
                (*m_log) << " " << std::endl;
                (*m_log) << " Final Parallel partition : " << std::endl;
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
                for(int ii=1; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
                }
                (*m_log) << " " << std::endl;
                (*m_log) << "---------------------------------------------" << std::endl;

            }
            else{
                m_loadBalanceRanges.clear();

                (*m_log) << " " << std::endl;
                (*m_log) << " Serial partition : " << std::endl;
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
                (*m_log) << " " << std::endl;
                (*m_log) << "---------------------------------------------" << std::endl;
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
            (*m_log) << "---------------------------------------------" << std::endl;
            (*m_log) << " LOAD BALANCE " << std::endl;

            m_lastOp = OP_LOADBALANCE;
            if (m_nproc>1){

                std::vector<uint32_t> partition(m_nproc);
                computePartition(level, weight, partition.data());

                privateLoadBalance(partition.data(), &userData);

                //Write info of final partition on log
                (*m_log) << " " << std::endl;
                (*m_log) << " Final Parallel partition : " << std::endl;
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
                for(int ii=1; ii<m_nproc; ii++){
                    (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
                }
                (*m_log) << " " << std::endl;
                (*m_log) << "---------------------------------------------" << std::endl;

            }
            else{
                m_loadBalanceRanges.clear();

                (*m_log) << " " << std::endl;
                (*m_log) << " Serial partition : " << std::endl;
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
                (*m_log) << " " << std::endl;
                (*m_log) << "---------------------------------------------" << std::endl;
            }

        }

        /**
        * Distribute Load-Balancing octants and user data of the whole
        * tree over the processes of the job following a given partition
        * distribution. Until loadBalance is not called for the first time
        * the mesh is serial.
        * \param[in] partition Target distribution of octants over processes.
        * \param[in,out] userData User data that will be distributed among the
        * processes.
        */
        template<class Impl>
        void
        privateLoadBalance(const uint32_t *partition, DataLBInterface<Impl> *userData = nullptr){

            (*m_log) << " " << std::endl;
            if (m_serial) {
                (*m_log) << " Initial Serial distribution : " << std::endl;
            } else {
                (*m_log) << " Initial Parallel partition : " << std::endl;
            }
            (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
            for(int ii=1; ii<m_nproc; ii++){
                (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
            }

            // Compute load balance ranges
            std::unordered_map<int, std::array<uint32_t, 2>> sendRanges = evalLoadBalanceSendRanges(partition);
            std::unordered_map<int, std::array<uint32_t, 2>> recvRanges = evalLoadBalanceRecvRanges(partition);

            m_loadBalanceRanges = LoadBalanceRanges(m_serial, sendRanges, recvRanges);

            // Compute information about the new partitioning
            assert(m_nproc > 0);

            uint32_t newSizeOctants = partition[m_rank];

            uint64_t newLastOctantGlobalIdx;
            uint64_t newFirstOctantGlobalIdx;
            if (newSizeOctants > 0) {
                newLastOctantGlobalIdx = partition[0];
                for (int p = 1; p < m_rank + 1; ++p) {
                    newLastOctantGlobalIdx += partition[p];
                }
                --newLastOctantGlobalIdx;

                newFirstOctantGlobalIdx = newLastOctantGlobalIdx - newSizeOctants + 1;
            } else {
                newLastOctantGlobalIdx  = 0;
                newFirstOctantGlobalIdx = 1;
            }

            // Partition internal octants
            if (m_serial) {
                m_lastOp = OP_LOADBALANCE_FIRST;

                if (newFirstOctantGlobalIdx != 0) {
                    for (uint32_t i = 0; i < newSizeOctants; ++i){
                        m_octree.m_octants[i] = m_octree.m_octants[newFirstOctantGlobalIdx + i];
                        if (userData) {
                            userData->move(newFirstOctantGlobalIdx + i, i);
                        }
                    }
                }

                m_octree.m_octants.resize(newSizeOctants);
                m_octree.m_octants.shrink_to_fit();
                m_octree.m_sizeOctants = m_octree.m_octants.size();

                if (userData) {
                    userData->resize(m_octree.m_sizeOctants);
                    userData->shrink();
                }
            } else {
                // Compute information about the current partitioning
                uint64_t lastOctantGlobalIdx  = m_partitionRangeGlobalIdx[m_rank];
                uint64_t firstOctantGlobalIdx = lastOctantGlobalIdx - m_octree.m_sizeOctants + 1;

                // Initialize communications
                DataCommunicator lbCommunicator(m_comm);

                if (!userData || userData->fixedSize()) {
                    for (const auto &entry : recvRanges) {
                        int rank = entry.first;
                        uint32_t beginRecvIdx = entry.second[0];
                        uint32_t endRecvIdx   = entry.second[1];

                        uint32_t nOctantsToReceive = endRecvIdx - beginRecvIdx;
                        std::size_t buffSize = nOctantsToReceive * Octant::getBinarySize();
                        if (userData) {
                            buffSize += nOctantsToReceive * userData->fixedSize();
                        }

                        lbCommunicator.setRecv(rank, buffSize);
                        lbCommunicator.startRecv(rank);
                    }
                }

                for (const auto &entry : sendRanges) {
                    int rank = entry.first;
                    uint32_t beginSendIdx = entry.second[0];
                    uint32_t endSendIdx   = entry.second[1];

                    uint32_t nOctantsToSend = endSendIdx - beginSendIdx;
                    std::size_t buffSize = nOctantsToSend * Octant::getBinarySize();
                    if (userData) {
                        if (userData->fixedSize()) {
                            buffSize += nOctantsToSend * userData->fixedSize();
                        }  else {
                            for (uint32_t i = beginSendIdx; i < endSendIdx; ++i) {
                                buffSize += userData->size(i);
                            }
                        }
                    }
                    lbCommunicator.setSend(rank, buffSize);

                    SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(rank);
                    for (uint32_t i = beginSendIdx; i < endSendIdx; ++i) {
                        sendBuffer << m_octree.m_octants[i];
                        userData->gather(sendBuffer, i);
                    }

                    lbCommunicator.startSend(rank);
                }

                if (userData && !userData->fixedSize()) {
                    lbCommunicator.discoverRecvs();
                    lbCommunicator.startAllRecvs();
                }

                // Move resident octants into place
                if (newSizeOctants > m_octree.m_sizeOctants) {
                    m_octree.m_octants.reserve(newSizeOctants);
                    m_octree.m_octants.resize(newSizeOctants);
                    m_octree.m_sizeOctants = m_octree.m_octants.size();

                    if (userData) {
                        userData->resize(m_octree.m_sizeOctants);
                    }
                }

                bool hasResidentOctants;
                if (newSizeOctants > 0) {
                    hasResidentOctants = true;
                    if (newFirstOctantGlobalIdx > lastOctantGlobalIdx) {
                        hasResidentOctants = false;
                    } else if (newLastOctantGlobalIdx < firstOctantGlobalIdx) {
                        hasResidentOctants = false;
                    }
                } else {
                    hasResidentOctants = false;
                }

                if (hasResidentOctants) {
                    uint32_t firstResidentGlobalOctantIdx = std::max(firstOctantGlobalIdx, newFirstOctantGlobalIdx);
                    uint32_t lastResidentGlobalOctantIdx  = std::min(lastOctantGlobalIdx, newLastOctantGlobalIdx);

                    uint32_t nofResidents = lastResidentGlobalOctantIdx - firstResidentGlobalOctantIdx + 1;
                    uint32_t newFirstResidentOffsetIdx = firstResidentGlobalOctantIdx - newFirstOctantGlobalIdx;
                    uint32_t firstResidentOffsetIdx = firstResidentGlobalOctantIdx - firstOctantGlobalIdx;

                    // If residents are moved closer to the head, we need to move
                    // them from the first to the last. Otherwise, if resident are
                    // moved closer to the tail, we need to move them in reversed
                    // order, i.e., from the last to the first (otherwise octants
                    // will be overwritten during the relocaiton).
                    if (newFirstResidentOffsetIdx != firstResidentOffsetIdx) {
                        for (uint32_t i = 0; i < nofResidents; ++i) {
                            int residentIdx;
                            if (newFirstResidentOffsetIdx < firstResidentOffsetIdx) {
                                residentIdx = i;
                            } else {
                                residentIdx = nofResidents - i - 1;
                            }

                            m_octree.m_octants[newFirstResidentOffsetIdx + residentIdx] = m_octree.m_octants[firstResidentOffsetIdx + residentIdx];
                            if (userData) {
                                userData->move(firstResidentOffsetIdx + residentIdx, newFirstResidentOffsetIdx + residentIdx);
                            }
                        }
                    }
                }

                if (newSizeOctants < m_octree.m_sizeOctants) {
                    m_octree.m_octants.resize(newSizeOctants);
                    m_octree.m_octants.shrink_to_fit();
                    m_octree.m_sizeOctants = m_octree.m_octants.size();

                    if (userData) {
                        userData->resize(m_octree.m_sizeOctants);
                        userData->shrink();
                    }
                }

                // Read buffers and build new octants
                int nCompletedRecvs = 0;
                while (nCompletedRecvs < lbCommunicator.getRecvCount()) {
                    int senderRank = lbCommunicator.waitAnyRecv();
                    RecvBuffer &recvBuffer = lbCommunicator.getRecvBuffer(senderRank);

                    const std::array<uint32_t, 2> &recvRange = recvRanges.at(senderRank);
                    uint32_t beginRecvIdx = recvRange[0];
                    uint32_t endRecvIdx   = recvRange[1];
                    assert(userData || ((endRecvIdx - beginRecvIdx) == (recvBuffer.getSize() / (Octant::getBinarySize()))));
                    assert(!userData || !userData->fixedSize() || ((endRecvIdx - beginRecvIdx) == (recvBuffer.getSize() / (Octant::getBinarySize() + userData->fixedSize()))));

                    for (uint32_t i = beginRecvIdx; i < endRecvIdx; ++i) {
                        recvBuffer >> m_octree.m_octants[i];
                        if (userData) {
                            userData->scatter(recvBuffer, i);
                        }
                    }

                    ++nCompletedRecvs;
                }

                lbCommunicator.waitAllSends();
            }

            // Update load balance information
            updateLoadBalance();

            // Update ghosts
            m_octree.m_ghosts.clear();
            m_octree.m_sizeGhosts = 0;

            computeGhostHalo();

            if (userData) {
                userData->resizeGhost(m_octree.m_sizeGhosts);
            }
        };
#endif

        // =============================================================================== //

    };

}

#endif /* __BITPIT_PARA_TREE_HPP__ */
