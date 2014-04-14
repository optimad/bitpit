/*
 * test_2D3D.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */

#include "preprocessor_defines.dat"
#include <mpi.h>
#include "global.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);


	MPI::Finalize();

}
