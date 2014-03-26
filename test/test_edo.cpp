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
		double X, Y, Z, L;
		X = 0.0; Y = 0.0; Z = 0.0; L = 100.0;
		Class_Para_Tree ptree(X, Y, Z, L);
		vector<uint32_t> mapidx;

		clock_t start = clock();
		clock_t end = clock();

		ptree.octree.setBalance(0,false);
		ptree.octree.setMarker(0,2);
		bool done = ptree.adapt(mapidx);
		ptree.loadBalance();
		uint64_t nocts = ptree.octree.getNumOctants();

		// TORUS TEST
		double R = double(max_length)*0.25;
		double r = double(max_length)*0.15;
		for (int l=0; l<2; l++){
			for (int i=0; i<nocts; i++){
				double* center;
				Class_Octant oct = ptree.octree.extractOctant(i);
				center = oct.getCenter();
				//				if (sqrt(pow((center[0]-double(max_length)*0.5),2.0)+pow((center[1]-double(max_length)*0.5),2.0)+pow((center[2]-double(max_length)*0.5),2.0)) <= double(max_length)*0.4){
				if ( pow(R -sqrt(pow((center[0]-double(max_length)*0.5),2.0)+pow((center[1]-double(max_length)*0.5),2.0)),2.0)+pow((center[2]-double(max_length)*0.5),2.0) <= pow(r,2.0) ){
					ptree.octree.setMarker(i,1);
				}
				delete[] center;
			}
			bool done = ptree.adapt(mapidx);
			ptree.loadBalance();
			nocts = ptree.octree.getNumOctants();

//		// RANDOM TEST
//		for (int l=0; l<9; l++){
//			for (int i=0; i<nocts/3; i++){
//				int j = rand() %nocts;
//				ptree.octree.setMarker(j,1);
//			}
//			ptree.balance21();
//			ptree.adapt();
//			ptree.loadBalance();
//			nocts = ptree.octree.getNumOctants();

			end = clock();
			float seconds = (float)(end - start) / CLOCKS_PER_SEC;
			writeLog(" ");
			writeLog("---------------------------------------------");
			writeLog(" ");
			writeLog(" CPU time (sec)			:	"+to_string(seconds));
			writeLog(" ");
			writeLog("---------------------------------------------");
		}

//		ptree.octree.updateConnectivity();
//		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("local_tree"));
//		ptree.octree.updateghostsConnectivity();
//		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("bbalunbalGhostsaa"));
//		ptree.octree.clearConnectivity();
//		ptree.octree.clearghostsConnectivity();

/*
		vector<double> point;
		point.resize(3);
		point[0] = 25.0;
		point[1] = 42.0;
		point[2] = 36.0;
		Class_Octant* oct = ptree.getPointOwner(point);
		if(oct!=NULL){
			cout << ptree.rank << endl;
			cout << (ptree.getX(oct)) << endl;
			cout << (ptree.getY(oct)) << endl;
			cout << (ptree.getZ(oct)) << endl;
		}
*/


		uint32_t			sizeneigh, modsize;
		vector<uint32_t> 	neigh;
		vector<bool> 		isghost;
		uint8_t 			inode;
		uint8_t 			idx;

		for (idx=0; idx<ptree.octree.getNumOctants(); idx++){
			for (inode=0; inode<8; inode++){
				ptree.octree.findNodeNeighbours(idx, inode, neigh, isghost);
				for (int i=0; i<neigh.size(); i++){
				}
			}
		}

		ptree.updateConnectivity();
		writePhysicalTree(ptree.nodes,ptree.connectivity,ptree.ghostsnodes,ptree.ghostsconnectivity,ptree,("physical_tree"));

	}

	MPI::Finalize();

	return 0;
}



