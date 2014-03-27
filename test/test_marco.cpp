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
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);
	{
		Class_Para_Tree ptree;

		ptree.octree.setMarker(0,1);
		ptree.octree.refine();
		ptree.updateAdapt();

		vector<double> data(ptree.octree.getNumOctants(),(double)ptree.rank);
		vector<double> gData(ptree.octree.getSizeGhost(),-1.0);

		{
			User_Data_LB<vector<double> > lbCommHandle(data);
			ptree.loadBalance(lbCommHandle);
		}
		gData.resize(ptree.octree.getSizeGhost(),0.0);

		User_Data_Comm<vector<double> > commHandle(data,gData);
		ptree.communicate(commHandle);

		if(ptree.rank == 0){
			ptree.octree.setMarker(0,1);
		}
		ptree.octree.refine();
		ptree.updateAdapt();

//
//		for(int i = 0; i < gData.size(); ++i){
//			cout << "rank: " << ptree.rank << " ghost " << i << ": " << gData[i] << endl;
//		}

		//	//TEST PARALLEL LOAD BALANCE
		//	if(ptree.rank == 3){
		//		ptree.octree.setMarker(0,1);
		//	}
		//	ptree.octree.refine();
		////	ptree.updateAdapt();
		//	if(ptree.rank == 3){
		//		ptree.octree.setMarker(0,1);
		//	}
		//	ptree.octree.refine();
		//	ptree.updateAdapt();
		//	ptree.setPboundGhosts();
		//
		//	ptree.octree.computeConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalNoGhost");
		//	ptree.octree.computeghostsConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalGhost");
		//
		//
		//	ptree.loadBalance();
		//	ptree.updateLoadBalance();
		//	ptree.setPboundGhosts();
		//	ptree.octree.clearghostsConnectivity();
		//	ptree.octree.updateConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balNoGhost");
		//
		//
		//	ptree.octree.computeghostsConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balGhost");
		//
		//	if(ptree.rank == 1){
		//		for(int i = 0; i < 6; ++i)
		//			ptree.octree.setMarker(i,-1);
		//	}
		//	if(ptree.rank == 2){
		//		for(int i = 0; i < 5; ++i)
		//			ptree.octree.setMarker(i,-1);
		//	}
		//	if(ptree.rank == 3){
		//		for(int i = 0; i < 5; ++i)
		//			ptree.octree.setMarker(i,-1);
		//	}
		//	ptree.updateAdapt();
		//	ptree.setPboundGhosts();
		//	ptree.octree.coarse();
		////	ptree.updateAdapt();
		////	ptree.setPboundGhosts();
		//	ptree.octree.clearghostsConnectivity();
		//	ptree.octree.updateConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarseunbalNoGhost");
		//	ptree.octree.computeghostsConnectivity();
		//	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarseunbalGhost");

	}
	MPI::Finalize();

	return 0;
}



