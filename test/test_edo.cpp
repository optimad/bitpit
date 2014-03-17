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

		clock_t start = clock();
		clock_t end = clock();

		ptree.octree.setBalance(0,false);
		ptree.octree.setMarker(0,3);
		ptree.adapt();
		ptree.loadBalance();
		uint64_t nocts = ptree.octree.getNumOctants();

		// TORUS TEST
		double R = double(max_length)*0.25;
		double r = double(max_length)*0.15;
		for (int l=0; l<5; l++){
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
			ptree.balance21();
			ptree.adapt();
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

		int ciccio;
		cin >> ciccio;

		ptree.octree.updateConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("bbalunbalNoGhostsaa"));
		ptree.octree.updateghostsConnectivity();
		writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("bbalunbalGhostsaa"));
		ptree.octree.clearConnectivity();
		ptree.octree.clearghostsConnectivity();

	}

	MPI::Finalize();

	return 0;
}



