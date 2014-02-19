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

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree ptree;


	ptree.octree.setMarker(0,2);
	cout << "Bound 0-face : " << ptree.octree.extractOctant(0).getBound(0) << endl;
	ptree.octree.refine();
	ptree.octree.setMarker(5,1);
	ptree.octree.refine();
	ptree.octree.setMarker(6,2);
	ptree.octree.refine();
	ptree.octree.setMarker(6,1);
	ptree.octree.refine();
	ptree.octree.setMarker(6,1);
	ptree.octree.refine();
	ptree.octree.setMarker(6,1);
	ptree.octree.refine();
	ptree.octree.setMarker(6,2);
	ptree.octree.refine();
	ptree.octree.refine();

	Class_Octant oct_test;
	uint8_t sizehf=0;
	for (int i=0; i<nface; i++){
		oct_test =	ptree.octree.extractOctant(6);
		uint64_t *hfneigh = oct_test.computeHalfSizeMorton(i,sizehf);
		for (int j=0; j<sizehf; j++){
			cout << "Morton half-size idx=6 iface " << i << " : " << hfneigh[j] << endl;
		}
	}

	uint64_t numoctants = ptree.octree.getNumOctants();
	cout << "Num Octants : " << numoctants << endl;

	MPI::Finalize();

	return 0;
}

