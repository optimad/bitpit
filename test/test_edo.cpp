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

	{
		Class_Para_Tree ptree;

		//TODO TESTARE QUANDO NON BILANCIATO (setBalance true), deadlock in pboundghosts !!!
		ptree.octree.setBalance(0,false);
		ptree.octree.setMarker(0,3);
		ptree.adapt();
		ptree.loadBalance();
		uint64_t nocts = ptree.octree.getNumOctants();

		for (int l=0; l<3; l++){
			for (int i=0; i<nocts; i++){
				double* center;
				Class_Octant oct = ptree.octree.extractOctant(i);
				center = oct.getCenter();
//				if ((center[0] <  double(max_length)*0.8) && (center[0] > double(max_length)*0.2)){
//					if ((center[1] < double(max_length)*0.8) && (center[1] >  double(max_length)*0.2)){
//						if ((center[2] < double(max_length)*0.8) && (center[2] >  double(max_length)*0.2)){
				if (sqrt(pow((center[0]-double(max_length)*0.5),2.0)+pow((center[1]-double(max_length)*0.5),2.0)+pow((center[2]-double(max_length)*0.5),2.0)) <= double(max_length)*0.4){
							ptree.octree.setMarker(i,1);
						}
//					}
//				}
				delete[] center;
			}
			ptree.balance21();
			ptree.adapt();
			ptree.loadBalance();
			nocts = ptree.octree.getNumOctants();
			ptree.octree.clearghostsConnectivity();
			ptree.octree.updateConnectivity();
			writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("bbalunbalNoGhostsaa"+to_string(l)));
			ptree.octree.updateghostsConnectivity();
			writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("bbalunbalGhostsaa"+to_string(l)));
		}

//		ptree.octree.updateConnectivity();
//		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"bbalunbalNoGhostsaa");


//		for (int l=0; l<2; l++){
//			ptree.loadBalance();
//			for (int i=0; i<nocts; i++){
//				double* center;
//				Class_Octant oct = ptree.octree.extractOctant(i);
//				center = oct.getCenter();
//				if ((center[0] <  double(max_length)*0.8) && (center[0] > double(max_length)*0.2)){
//					if ((center[1] < double(max_length)*0.8) && (center[1] >  double(max_length)*0.2)){
//						if ((center[2] < double(max_length)*0.8) && (center[2] >  double(max_length)*0.2)){
//							ptree.octree.setMarker(i,1);
//						}
//					}
//				}
//				delete[] center;
//			}
//			ptree.balance21();
//			ptree.adapt();
////			cout << "loadbalance" << endl;
//			nocts = ptree.octree.getNumOctants();
//		}
//
//		ptree.octree.updateConnectivity();
//		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"bbalunbalNoGhostsbb");
//
		/*
		if(ptree.rank == 2){
			ptree.octree.setMarker(1,1);
		}
		ptree.adapt();
		if(ptree.rank == 2){
			ptree.octree.setMarker(3,1);
		}
		ptree.adapt();

		ptree.octree.computeConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalNoGhost");
		ptree.octree.computeghostsConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"unbalGhost");

		ptree.loadBalance();

		ptree.octree.clearghostsConnectivity();
		ptree.octree.updateConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balNoGhost");
		ptree.octree.computeghostsConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"balGhost");

		if(ptree.rank == 1){
			for(int i = 0; i < ptree.octree.getNumOctants(); ++i)
				ptree.octree.setMarker(i,-1);
		}
		if(ptree.rank == 2){
			for(int i = 0; i < ptree.octree.getNumOctants(); ++i){
				ptree.octree.setMarker(i,-1);
			}
		}
		if(ptree.rank == 3){
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

		ptree.octree.clearghostsConnectivity();
		ptree.octree.updateConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarsebalNoGhost");
		ptree.octree.computeghostsConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,"coarsebalGhost");
*/

	}

	MPI::Finalize();

	return 0;
}



