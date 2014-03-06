/*
 * test_edo.cpp
 *
 *  Created on: 5/mar/2014
 *      Author: Edoardo Lombardi
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

	ptree.octree.setBalance(0,1);
	ptree.octree.setMarker(0,1);
	while(ptree.octree.refine());
	ptree.updateAdapt();


	ptree.loadBalance();
	ptree.updateLoadBalance();
	ptree.setPboundGhosts();

	//TEST PARALLEL LOAD BALANCE
	if(ptree.rank == 1){
		ptree.octree.setMarker(1,1);
	}
	ptree.adapt();
	if(ptree.rank == 1){
		ptree.octree.setMarker(3,1);
	}
	ptree.adapt();

	ptree.octree.computeConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalNoGhost");
	ptree.octree.computeghostsConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalGhost");

	ptree.loadBalance();
	ptree.updateLoadBalance();
	ptree.setPboundGhosts();

	ptree.octree.clearghostsConnectivity();
	ptree.octree.updateConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balNoGhost");
	ptree.octree.computeghostsConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balGhost");

	if(ptree.rank == 0){
		for(int i = 0; i < ptree.octree.getNumOctants(); ++i)
			ptree.octree.setMarker(i,-1);
	}
	if(ptree.rank == 1){
		for(int i = 0; i < ptree.octree.getNumOctants(); ++i){
			ptree.octree.setMarker(i,-1);
		}
	}
	ptree.adapt();

	ptree.octree.clearghostsConnectivity();
	ptree.octree.updateConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarseunbalNoGhost");
	ptree.octree.computeghostsConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarseunbalGhost");

	ptree.loadBalance();
	ptree.updateLoadBalance();
	ptree.setPboundGhosts();

	ptree.octree.clearghostsConnectivity();
	ptree.octree.updateConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarsebalNoGhost");
	ptree.octree.computeghostsConnectivity();
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarsebalGhost");


	MPI::Finalize();

	return 0;
}



