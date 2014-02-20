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
#include "Class_Para_Tree.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree ptree;
	ptree.octree.setMarker(0,1);
	cout << "Bound 0-face : " << ptree.octree.extractOctant(0).getBound(0) << endl;
	ptree.octree.refine();
	ptree.updateRefine();


	MPI::Finalize();

	return 0;
}



