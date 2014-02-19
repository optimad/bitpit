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
/*
	uint8_t lev = 0;
	uint32_t x, y, z;
	x = y = z = 0;
	Class_Octant oct0(lev,x,y,z);
	vector<Class_Octant> child(8);
	oct0.buildchildren(child);
	for (int i=0; i<nchildren; i++){
		cout << "----------" << endl;
		cout << "ch : " << i << endl;
		cout << "x : " << child[i].getx() << endl;
		cout << "y : " << child[i].gety() << endl;
		cout << "z : " << child[i].getz() << endl;
		cout << "level : " << int(child[i].getlevel()) << endl;
	}
*/

	uint64_t numoctants = ptree.octree.getNumOctants();
	ptree.octree.setMarker(0,true);
	ptree.octree.refine();
	cout << ptree.octree.getNumOctants() << endl;
	cout << int(ptree.octree.getLocalMaxDepth()) << endl;
	ptree.octree.setMarker(5,true);
	ptree.octree.refine();
	cout << ptree.octree.getNumOctants() << endl;
	cout << int(ptree.octree.getLocalMaxDepth()) << endl;

	MPI::Finalize();

	return 0;
}

