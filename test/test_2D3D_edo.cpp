/*
 * test_2D3D.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */

#include "preprocessor_defines.dat"
#include <mpi.h>
#include "global.hpp"
#include "Class_Para_Tree.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree<2> ptree;
	ptree.octree.setBalance(0,false);
//	ptree.octree.setMarker(0,2);
//	bool done = ptree.adapt();

//	cout << "done " << done << endl;
//	cout << "done " << done << endl;
//	cout << "done " << done << endl;

	MPI::Finalize();

}
