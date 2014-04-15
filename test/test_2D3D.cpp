/*
 * test_2D3D.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */

#include "preprocessor_defines.dat"
#include <mpi.h>
#include "global.hpp"
#include "Class_Octant.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);
	cout << "global 2D nfaces = " << (int)global2D.puppa() << endl;
	cout << "global 3D nfaces = " << (int)global3D.nfaces << endl;

	Class_Octant<3> ottante(0,1,2,3);

	cout << "get x =" << ottante.getX();

	MPI::Finalize();

}
