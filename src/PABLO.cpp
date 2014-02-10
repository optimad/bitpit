/*
 ============================================================================
 Name        : PABLO.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Compute Pi in MPI C++
 ============================================================================
 */
#include <math.h> 
#include "mpi.h" 
#include <iostream>
#include "Class_Octree.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	/*int n, rank, size, i;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, h, sum, x;

	MPI::Init(argc, argv);
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	n=1000; // number of intervals

	MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0);
	h = 1.0 / (double) n;
	sum = 0.0;
	for (i = rank + 1; i <= n; i += size) {
		x = h * ((double) i - 0.5);
		sum += (4.0 / (1.0 + x * x));
	}
	mypi = h * sum;

	MPI::COMM_WORLD.Reduce(&mypi, &pi, 1, MPI::DOUBLE, MPI::SUM, 0);
	if (rank == 0)
		cout << "pi is approximately " << pi << ", Error is "
				<< fabs(pi - PI25DT) << endl;

	MPI::Finalize();*/

	Class_Octant oct0(20);

	cout << int(oct0.level) << endl;


	return 0;
}

