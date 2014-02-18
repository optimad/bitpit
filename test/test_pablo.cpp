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
	uint8_t lev = 10;
	uint32_t x, y, z;
	x = y = z = 0;
	Class_Octant oct0(lev,x,y,z);
	Class_Octant oct1(oct0);
	cout << oct0.getsize() << endl;
	cout << oct1.getsize() << endl;
	Class_Local_Tree tree;

	MPI::Finalize();

	return 0;
}

