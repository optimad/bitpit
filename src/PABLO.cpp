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

	uint8_t a = 10;
	uint8_t x, y, z;
	x = y = z = 0;
	Class_Octant oct0(a,x,y,z);
	int l = oct0.getlevel();
	cout << "level oct0 : " << l << endl;

	Class_Octant oct1(oct0);
	l = oct1.getlevel();
	cout << "level oct1 : " << l << endl;
	int x1 = oct1.getx();
	int y1 = oct1.gety();
	int z1 = oct1.getz();
	cout << "x oct1 : " << x1 << endl;
	cout << "y oct1 : " << y1 << endl;
	cout << "z oct1 : " << z1 << endl;
	bool balance = oct1.getbalance();
	cout << "balance oct1 : " << balance << endl;

	return 0;
}

