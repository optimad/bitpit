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

	ptree.octree.setMarker(0,1);
	ptree.octree.refine();
	ptree.updateRefine();
//	ptree.octree.setMarker(4,4);
//	ptree.octree.refine();
//	ptree.updateRefine();
//	ptree.octree.refine();
//	ptree.updateRefine();
//	ptree.octree.refine();
//	ptree.updateRefine();
//	ptree.octree.refine();
//	ptree.updateRefine();
	cout << "I'm " << ptree.rank << " and max_depth is " << (int)ptree.max_depth << endl;
	cout << "I'm " << ptree.rank << " and global_num_octants is " << ptree.global_num_octants << endl;
	cout << "I'm " << ptree.rank << " and partition_range_globalidx";
	for(int i = 0; i < ptree.nproc; ++i)
		cout << " " << ptree.partition_range_globalidx[i];
	cout << endl;

	ptree.loadBalance();
	ptree.updateLoadBalance();

	cout << "I'm " << ptree.rank << " and I have " << ptree.octree.getNumOctants() << " octants" << endl;
	cout << "I'm " << ptree.rank << " and I see ";
	for(int i = 0; i < ptree.nproc; ++i)
		cout << ptree.partition_range_globalidx[i] << " ";
	cout << "as global partition" << endl;
	cout << "I'm " << ptree.rank << " and I see ";
	for(int i = 0; i < ptree.nproc; ++i)
		cout << ptree.partition_last_desc[i] << " ";
	cout << "as last descendant partition" << endl;

	ptree.setPboundGhosts();

	ptree.octree.computeConnectivity();
	ptree.octree.computeghostsConnectivity();
	string filename = "puppa";
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,filename);

	MPI::Finalize();

	return 0;
}



