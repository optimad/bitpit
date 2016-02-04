//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\example cartesianpatch_example_00001.cpp

	\brief Use cartesianpatch to create a 2D Cartesian patch

	This example creates a 2D Cartesian patch on the square domain [-10,10]x[-10,10].
	The domain is discretized with a cell size of 0.5 both in x and y directions.

	<b>To run</b>: ./cartesianpatch_example_00001 \n
*/

#include <array>
#if ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_cartesianpatch.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

	std::cout << "Creating a 2D Cartesian patch" << "\n";

#if ENABLE_MPI==1
	MPI::Init(argc,argv);
#endif

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	Patch *patch_2D = new CartesianPatch(0, 2, origin, length, dh);
	patch_2D->setName("cartesian_2D_patch");
	patch_2D->update();
	patch_2D->writeMesh();

#if ENABLE_MPI==1
	MPI::Finalize();
#endif

}
