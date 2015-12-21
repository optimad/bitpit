#ifndef CLASSPARATREE_HPP_
#define CLASSPARATREE_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#if NOMPI==0
#include <mpi.h>
#endif
#include "classGlobal.hpp"
#include "classOctant.hpp"
#include "classLocalTree.hpp"
#include "Class_Comm_Buffer.hpp"
#include "classMap.hpp"
#include "Class_Array.hpp"
#include "Class_Data_Comm_Interface.hpp"
#include "Class_Data_LB_Interface.hpp"
#include "Class_Log.hpp"
#include <cstdint>
#include <iterator>
#include <set>
#include <algorithm>
#include <string>
#include <functional>
#include <cctype>
#include <fstream>
#include <iomanip>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// TYPEDEFS
// =================================================================================== //


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*!
 *  \ingroup        PABLO
 *  @{
 *	\date			17/dec/2015
 *	\authors		Marco Cisternino
 *	\authors		Edoardo Lombardi
 *
 *	\brief Para Tree is the user interface class
 *
 *	The user should (read can...) work only
 *	with this Class and its methods.
 *	The sizes are intended in physical domain. The transformation from the logical
 *	domain to the physical domain is defined by classMap trans.
 *
 *	The partition of the octree is performed by following the Z-curve defined by the Morton
 *	index of the octants. By default it is a balanced partition over the number of octants for each
 *	process.
 *
 *	Class classParaTree has a dimensional parameter int dim and it accepts only two values: dim=2 (Class_Para_Tree<2>)and dim=3 (Class_Para_Tree<3>), obviously for 2D and 3D respectively.
 */
class classParaTree{

	// =================================================================================== //
	// MEMBERS
	// =================================================================================== //
private:
	//undistributed members
	uint64_t* 			partition_first_desc; 			/**<Global array containing position of the first possible octant in each processor*/
	uint64_t*			partition_last_desc; 			/**<Global array containing position of the last possible octant in each processor*/
	uint64_t* 			partition_range_globalidx;	 	/**<Global array containing global index of the last existing octant in each processor*/
	uint64_t 			global_num_octants;   			/**<Global number of octants in the parallel octree*/
	map<int,u32vector> 	bordersPerProc;					/**<Local indices of border octants per process*/
	int 				nproc;							/**<Number of processes of the job*/
	uint8_t 			max_depth;						/**<Global max existing level in the parallel octree*/
	classGlobal			global;							/**<Global variables*/

	//distributed members
	int 				rank;							/**<Local rank of process*/
	classLocalTree 		octree;							/**<Local tree in each processor*/

	//distributed adpapting memebrs
	u32vector 			mapidx;							/**<Local mapper for adapting. Mapper from new octants to old octants.
														mapidx[i] = j -> the i-th octant after adapt was in the j-th position before adapt;
														if the i-th octant is new after refinement the j-th old octant was the father of the new octant;
														if the i-th octant is new after coarsening the j-th old octant was the first child of the new octant.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 */

	//auxiliary members
	int 				error_flag;						/**<MPI error flag*/
	bool 				serial;							/**<True if the octree is the same on each processor, False if the octree is distributed*/

	//map members
	classMap 			trans;							/**<Transformation map from logical to physical domain*/
	uint8_t				dim;							/**<Space dimension of the octree object (2D/3D).*/

	//info member
	uint64_t			status;							/**<Label of actual status of octree (incremental after an adpat
														with at least one modifyed element).*/

	//log member
	Class_Log 			log;							/**<Log object*/

#if NOMPI==0
	MPI_Comm 			comm;							/**<MPI communicator*/
#endif

	// =================================================================================== //
	// CONSTRUCTORS AND OPERATORS
	// =================================================================================== //
public:

#if NOMPI==0
	classParaTree(uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);// : log(logfile,comm_),comm(comm_);
	classParaTree(double X, double Y, double Z, double L, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);//:dim(2),trans(X,Y,Z,L),log(logfile,comm_),comm(comm_);
	classParaTree(double X, double Y, double Z, double L, u32vector2D & XY, u8vector & levels, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log", MPI_Comm comm_ = MPI_COMM_WORLD);//:trans(X,Y,Z,L),log(logfile,comm_),comm(comm_);
#else
	classParaTree(uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");// : log(logfile);
	classParaTree(double X, double Y, double Z, double L, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");//:dim(2),trans(X,Y,Z,L),log(logfile);
	classParaTree(double X, double Y, double Z, double L, u32vector2D & XY, u8vector & levels, uint8_t dim_ = 2, int8_t maxlevel = 20, string logfile="PABLO.log");//:trans(X,Y,Z,L),log(logfile);
#endif

	~classParaTree();

	// =================================================================================== //
	// METHODS
	// =================================================================================== //

	// =================================================================================== //
	// BASIC GET/SET METHODS
	// =================================================================================== //

	int getRank();
	int getMaxLevel();

	void setMaxLevel(int8_t maxlevel);

	// =================================================================================== //
	// INDEX BASED METHODS
	// =================================================================================== //

	double getX(uint32_t idx);
	double getY(uint32_t idx);
	double getZ(uint32_t idx);
	double getSize(uint32_t idx);
	double getArea(uint32_t idx);
	double getVolume(uint32_t idx);
	void getCenter(uint32_t idx, dvector& center);
	dvector getCenter(uint32_t idx);
	dvector getFaceCenter(uint32_t idx, uint8_t iface);
	void getFaceCenter(uint32_t idx, uint8_t iface, dvector& center);
	dvector getNode(uint32_t idx, uint8_t inode);
	void getNode(uint32_t idx, uint8_t inode, dvector& node);
	void getNodes(uint32_t idx, dvector2D & nodes);
	dvector2D getNodes(uint32_t idx);
	void getNormal(uint32_t idx, uint8_t & iface, dvector & normal);
	dvector getNormal(uint32_t idx, uint8_t & iface);
	int8_t getMarker(uint32_t idx);
	uint8_t getLevel(uint32_t idx);
	bool getBalance(uint32_t idx);
#if NOMPI==0
	bool getIsGhost(uint32_t idx);
#endif
	bool getIsNewR(uint32_t idx);
	bool getIsNewC(uint32_t idx);
	uint64_t getGlobalIdx(uint32_t idx);
	uint64_t getGhostGlobalIdx(uint32_t idx);
	void setMarker(uint32_t idx, int8_t marker);
	void setBalance(uint32_t idx, bool balance);

	// =================================================================================== //
	// POINTER BASED METHODS
	// =================================================================================== //

	double getX(classOctant* oct);
	double getY(classOctant* oct);
	double getZ(classOctant* oct);
	double getSize(classOctant* oct);
	double getArea(classOctant* oct);
	double getVolume(classOctant* oct);
	void getCenter(classOctant* oct, dvector& center);
	dvector getCenter(classOctant* oct);
	dvector getFaceCenter(classOctant* oct, uint8_t iface);
	void getFaceCenter(classOctant* oct, uint8_t iface, dvector& center);
	dvector getNode(classOctant* oct, uint8_t inode);
	void getNode(classOctant* oct, uint8_t inode, dvector& node);
	void getNodes(classOctant* oct, dvector2D & nodes);
	dvector2D getNodes(classOctant* oct);
	void getNormal(classOctant* oct, uint8_t & iface, dvector & normal);
	dvector getNormal(classOctant* oct, uint8_t & iface);
	int8_t getMarker(classOctant* oct);
	uint8_t getLevel(classOctant* oct);
	bool getBalance(classOctant* oct);
	bool getIsNewR(classOctant* oct);
	bool getIsNewC(classOctant* oct);
	void setMarker(classOctant* oct, int8_t marker);
	void setBalance(classOctant* oct, bool balance);


	// =================================================================================== //
	// LOCAL TREE GET/SET METHODS
	// =================================================================================== //

	uint64_t getStatus();
	uint32_t getNumOctants() const;
	uint32_t getNumGhosts() const;
	uint32_t getNumNodes() const;
	uint8_t getLocalMaxDepth() const;
	uint8_t getBalanceCodimension() const;
	void getBoundingBox(dvector & P0, dvector & P1);
	void getBoundingBox(darray3 & P0, darray3 & P1);
	void setBalanceCodimension(uint8_t b21codim);
	const classOctant & getFirstDesc() const;
	const classOctant & getLastDesc() const;
	uint64_t getLastDescMorton(uint32_t idx);

	// =================================================================================== //
	// INTERSECTION GET/SET METHODS
	// =================================================================================== //

	uint32_t getNumIntersections();
	classIntersection* getIntersection(uint32_t idx);
	uint8_t getLevel(classIntersection* inter);
	bool getFiner(classIntersection* inter);
	bool getBound(classIntersection* inter);
	bool getIsGhost(classIntersection* inter);
	bool getPbound(classIntersection* inter);
	uint8_t getFace(classIntersection* inter);
	u32vector getOwners(classIntersection* inter);
	uint32_t getIn(classIntersection* inter);
	uint32_t getOut(classIntersection* inter);
	double getSize(classIntersection* inter);
	double getArea(classIntersection* inter);
	dvector getCenter(classIntersection* inter);
	dvector2D getNodes(classIntersection* inter);
	dvector getNormal(classIntersection* inter);

	// =================================================================================== //
	// OTHER GET/SET METHODS
	// =================================================================================== //

	classOctant* getOctant(uint32_t idx);
	classOctant* getGhostOctant(uint32_t idx);
	uint64_t getGlobalIdx(classOctant* oct);
	uint32_t getIdx(classOctant* oct);
	uint32_t getIdx(classOctant oct);
#if NOMPI==0
	bool getIsGhost(classOctant* oct);
	bool getIsGhost(classOctant oct);
#endif

	// =================================================================================== //
	// PRIVATE GET/SET METHODS
	// =================================================================================== //
private:
	void setFirstDesc();
	void setLastDesc();

	// =================================================================================== //
	// OTHER METHODS												    			   //
	// =================================================================================== //

	// =================================================================================== //
	// OTHER OCTANT BASED METHODS												    			   //
	// =================================================================================== //
public:
	void findNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours, vector<bool> & isghost);
	void findNeighbours(classOctant* oct, uint8_t iface, uint8_t codim, u32vector & neighbours, vector<bool> & isghost);
	void findGhostNeighbours(uint32_t idx, uint8_t iface, uint8_t codim, u32vector & neighbours);
	classOctant* getPointOwner(dvector & point);
	uint32_t getPointOwnerIdx(dvector & point);
	void getMapping(uint32_t & idx, u32vector & mapper, vector<bool> & isghost);

	// =================================================================================== //
	// OTHER PARATREE BASED METHODS												    			   //
	// =================================================================================== //

	int findOwner(const uint64_t & morton);
	bool adapt(bool mapper_flag = false);
	bool adaptGlobalRefine(bool mapper_flag = false);
	bool adaptGlobalCoarse(bool mapper_flag = false);
	void computeConnectivity();
	void clearConnectivity();
	void updateConnectivity();
	const u32vector2D & getConnectivity();
	const u32vector & getConnectivity(uint32_t idx);
	const u32vector & getConnectivity(classOctant* oct);
	const u32arr3vector & getNodes();
	const u32array3 & getNodeLogicalCoordinates(uint32_t inode);
	dvector getNodeCoordinates(uint32_t inode);
	void computeGhostsConnectivity();
	void clearGhostsConnectivity();
	void updateGhostsConnectivity();
	const u32vector2D & getGhostConnectivity();
	const u32vector & getGhostConnectivity(uint32_t idx);
	const u32vector & getGhostConnectivity(classOctant* oct);
	const u32arr3vector & getGhostNodes();
	const u32array3 & getGhostNodeLogicalCoordinates(uint32_t inode);
	dvector getGhostNodeCoordinates(uint32_t inode);
	//TODO MapPablos
#if NOMPI==0
	template<class Impl>
	void communicate(Class_Data_Comm_Interface<Impl> & userData);
	void loadBalance();
	void loadBalance(uint8_t & level);
	template<class Impl>
	void loadBalance(Class_Data_LB_Interface<Impl> & userData, dvector* weight = NULL);
	template<class Impl>
	void loadBalance(Class_Data_LB_Interface<Impl> & userData, uint8_t & level);
#endif

	// =================================================================================== //
	// OTHER INTERSECTION BASED METHODS												    			   //
	// =================================================================================== //

	void computeIntersections();

	// =================================================================================== //
	// OTHER PRIVATE METHODS												    			   //
	// =================================================================================== //
private:
	classOctant& extractOctant(uint32_t idx);
	bool private_adapt();
	bool private_adapt_mapidx();
	void updateAdapt();
#if NOMPI==0
	void computePartition(uint32_t* partition);
	void computePartition(uint32_t* partition, dvector* weight);
	void computePartition(uint32_t* partition, uint8_t & level_);
	void updateLoadBalance();
	void setPboundGhosts();
	void commMarker();
#endif
	void updateAfterCoarse();
	void updateAfterCoarse(u32vector & mapidx);
	void balance21(bool const first);

	// =================================================================================== //
	// TESTING OUTPUT METHODS												    			   //
	// =================================================================================== //
public:
	void write(string filename);
	void writeTest(string filename, vector<double> data);

	// =============================================================================== //

};


/*  @}  */

#endif /* CLASSPARATREE_HPP_ */
