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
	int num_octants = 11;
	for(int i = 0; i < num_octants; ++i)
/*	{
		ptree.octree.octants.push_back(Class_Octant(0,0,0,0));
	}*/

	MPI::Finalize();

	return 0;
}



