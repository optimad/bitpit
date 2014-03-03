/*
 * test_marco.cpp
 *
 *  Created on: 18/feb/2014
 *      Author: Marco Cisternino
 */

#include "preprocessor_defines.dat"
#include <math.h>
#include "mpi.h"
#include <iostream>
#include <string>
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree ptree;

//	ptree.octree.setMarker(0,4);
//	while(ptree.octree.refine());
//	uint64_t numoctants = ptree.octree.getNumOctants();
//
//	for (int l=0; l<4; l++){
//		for (int i=0; i<numoctants; i++){
//			double* center;
//			Class_Octant oct = ptree.octree.extractOctant(i);
//			center = oct.getCenter();
//			if ((center[0] < pow(2.0,20)*0.75) && (center[0] > pow(2.0,20)*0.25)){
//				if ((center[1] < pow(2.0,20)*0.75) && (center[1] > pow(2.0,20)*0.25)){
//					if ((center[2] < pow(2.0,20)*0.75) && (center[2] > pow(2.0,20)*0.25)){
//						ptree.octree.setMarker(i,1);
//					}
//				}
//			}
//		}
//		ptree.octree.refine();
//		numoctants = ptree.octree.getNumOctants();
//		cout << "Num Octants : " << numoctants << endl;
//	}
/*
	ptree.octree.computeConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"ciccio");
	ptree.octree.clearConnectivity();
*/
	ptree.octree.setMarker(0,1);
	ptree.octree.refine();
	ptree.updateAdapt();
//	ptree.octree.setMarker(0,1);
//	ptree.octree.refine();
//	ptree.updateAdapt();
//	ptree.octree.refine();
//	ptree.updateAdapt();
//	ptree.octree.refine();
//	ptree.updateAdapt();
//	ptree.octree.refine();
//	ptree.updateAdapt();

//	cout << "Coarse " << endl;
//	numoctants = ptree.octree.getNumOctants();
//	for (int l=0; l<4; l++){
//		for (int i=0; i<numoctants; i++){
//			double* center;
//			Class_Octant oct = ptree.octree.extractOctant(i);
//			center = oct.getCenter();
//			if ((center[0] < pow(2.0,20)*0.65) && (center[0] > pow(2.0,20)*0.35)){
//				if ((center[1] < pow(2.0,20)*0.65) && (center[1] > pow(2.0,20)*0.35)){
//					if ((center[2] < pow(2.0,20)*0.65) && (center[2] > pow(2.0,20)*0.35)){
//						ptree.octree.setMarker(i,-1);
//					}
//				}
//			}
//		}
//		cout << "Num Octants : " << numoctants << endl;
//		ptree.octree.coarse();
//		numoctants = ptree.octree.getNumOctants();
//		cout << "Num Octants : " << numoctants << endl;
//		ptree.updateAdapt();
//	}

//	cout << "I'm " << ptree.rank << " and max_depth is " << (int)ptree.max_depth << endl;
//	cout << "I'm " << ptree.rank << " and global_num_octants is " << ptree.global_num_octants << endl;
//	cout << "I'm " << ptree.rank << " and partition_range_globalidx";
//	for(int i = 0; i < ptree.nproc; ++i)
//		cout << " " << ptree.partition_range_globalidx[i];
//	cout << endl;

	ptree.loadBalance();
	ptree.updateLoadBalance();

//	cout << "I'm " << ptree.rank << " and I have " << ptree.octree.getNumOctants() << " octants" << endl;
//	cout << "I'm " << ptree.rank << " and I see ";
//	for(int i = 0; i < ptree.nproc; ++i)
//		cout << ptree.partition_range_globalidx[i] << " ";
//	cout << "as global partition" << endl;
//	cout << "I'm " << ptree.rank << " and I see ";
//	for(int i = 0; i < ptree.nproc; ++i)
//		cout << ptree.partition_last_desc[i] << " ";
//	cout << "as last descendant partition" << endl;

	ptree.setPboundGhosts();

	//TEST PARALLEL LOAD BALANCE
	if(ptree.rank == 0){
		ptree.octree.setMarker(0,1);
	}
	ptree.octree.refine();
	ptree.updateAdapt();

	ptree.loadBalance();

	ptree.octree.computeConnectivity();
	ptree.octree.computeghostsConnectivity();
	string filename = "puppa";
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,filename);

	MPI::Finalize();

	return 0;
}



