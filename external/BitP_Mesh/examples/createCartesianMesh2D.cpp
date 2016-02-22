//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\example createCartesianMesh2D.cpp

	\brief Use patchman to create a 2D Cartesian mesh

	This example creates a 2D Cartesian mesh on the square domain [-10,10]x[-10,10].
	The domain is discretized with a cell size of 0.5 both in x and y directions.

	<b>To run</b>: ./createCartesianMesh2D \n
*/

#include <array>

#include "BitP_Mesh_PATCHMAN.hpp"

int main(int argc, char *argv[]) {

	std::cout << "Creating a 2D Cartesian mesh" << "\n";

#ifndef DISABLE_MPI
	MPI::Init(argc,argv);
#endif

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	pman::Patch *patch_2D = new pman::PatchCartesian(0, 2, origin, length, dh);
	patch_2D->set_name("cartesian_2D_mesh");
	patch_2D->update();
	patch_2D->write_mesh();

#ifndef DISABLE_MPI
	MPI::Finalize();
#endif

}
