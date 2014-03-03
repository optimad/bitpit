/*
 ============================================================================
 Name        : PABLO.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Compute Pi in MPI C++
 ============================================================================
 */
#include "preprocessor_defines.dat"
#include <math.h> 
#include "mpi.h" 
#include <iostream>
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree ptree;


	ptree.octree.setMarker(0,2);
	while(ptree.octree.refine());
	uint64_t numoctants = ptree.octree.getNumOctants();

	for (int l=0; l<2; l++){
		for (int i=0; i<numoctants; i++){
			double* center;
			Class_Octant oct = ptree.octree.extractOctant(i);
			center = oct.getCenter();
			if ((center[0] < pow(2.0,20)*0.75) && (center[0] > pow(2.0,20)*0.25)){
				if ((center[1] < pow(2.0,20)*0.75) && (center[1] > pow(2.0,20)*0.25)){
					if ((center[2] < pow(2.0,20)*0.75) && (center[2] > pow(2.0,20)*0.25)){
						ptree.octree.setMarker(i,1);
					}
				}
			}
			delete[] center;
		}
		ptree.octree.refine();
		numoctants = ptree.octree.getNumOctants();
		cout << "Num Octants : " << numoctants << endl;
	}
	cout << "Connectivity " << endl;
	ptree.octree.computeConnectivity();
	string filename = "test_r";
	cout << "Write " << endl;
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,filename);
	ptree.octree.clearConnectivity();


	cout << "Coarse " << endl;
	numoctants = ptree.octree.getNumOctants();
	for (int l=0; l<4; l++){
		for (int i=0; i<numoctants; i++){
			double* center;
			Class_Octant oct = ptree.octree.extractOctant(i);
			center = oct.getCenter();
			if ((center[0] < pow(2.0,20)*0.65) && (center[0] > pow(2.0,20)*0.35)){
				if ((center[1] < pow(2.0,20)*0.65) && (center[1] > pow(2.0,20)*0.35)){
					if ((center[2] < pow(2.0,20)*0.65) && (center[2] > pow(2.0,20)*0.35)){
						ptree.octree.setMarker(i,-1);
					}
				}
			}
			delete[] center;
		}
		cout << "Num Octants : " << numoctants << endl;
		ptree.octree.coarse();
		numoctants = ptree.octree.getNumOctants();
		cout << "Num Octants : " << numoctants << endl;
	}

	ptree.octree.computeConnectivity();
	filename = "test_c";
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,filename);
	ptree.octree.clearConnectivity();

	cout << "Balancing " << endl;
	bool Bdone = ptree.octree.localBalanceWithLevel();
	cout << "Bdone: " << Bdone << endl;
	cout << " refinement " << endl;
	while(ptree.octree.refine());
	cout << " refinement done " << endl;
	ptree.octree.computeConnectivity();
	filename = "test_bal";
	writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,filename);

	MPI::Finalize();

	return 0;
}

