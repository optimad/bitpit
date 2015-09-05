//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <array>

#include "patchman.hpp"

int main(int argc, char *argv[]) {

	std::cout << "Testing Octree mesh" << "\n";

#ifndef DISABLE_MPI
	MPI::Init(argc,argv);
#endif

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	std::cout << "  >> 2D Mesh" << "\n";

	pman::Patch *patch_2D = new pman::PatchOctree(0, 2, origin, length, dh);
	patch_2D->set_name("backgorund_mesh");
	patch_2D->update();
	patch_2D->output_write();

	std::cout << "  >> 3D Mesh" << "\n";

	pman::Patch *patch_3D = new pman::PatchOctree(0, 3, origin, length, dh);
	patch_3D->set_name("backgorund_mesh");
	patch_3D->update();
	patch_3D->output_write();

#ifndef DISABLE_MPI
	MPI::Finalize();
#endif

}

