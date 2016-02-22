//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <array>

#include "BitP_Mesh_PATCHMAN.hpp"

int main(int argc, char *argv[]) {

	std::cout << "Testing adaption on Octree mesh" << std::endl;

#ifndef DISABLE_MPI
	MPI::Init(argc,argv);
#endif

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 1.0;

	std::cout << std::endl;
	std::cout << "  :: 2D Mesh adaption test ::" << std::endl;

	std::cout << std::endl;
	std::cout << ">> Creating the mesh" << std::endl;

	pman::Patch *patch_2D = new pman::PatchOctree(0, 2, origin, length, dh);
	patch_2D->set_name("octree_adapted_mesh_2D");
	patch_2D->update();
	patch_2D->write_mesh();

	std::cout << std::endl;
	std::cout << ">> Adapting mesh" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch_2D->get_cell_count();
		std::cout << std::endl;
		std::cout << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch_2D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_2D->extract_cell_neighs(cellId)) {
				for (auto coarseId : patch_2D->extract_cell_neighs(neighId)) {
					patch_2D->mark_cell_for_coarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch_2D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_2D->extract_cell_neighs(cellId)) {
				patch_2D->mark_cell_for_refinement(neighId);
			}
		}

		std::cout << std::endl;
		std::cout << ">> Initial number of cells... " << nCells << std::endl;

		patch_2D->update();

		nCells = patch_2D->get_cell_count();
		std::cout << ">> Final number of cells... " << nCells << std::endl;
	}
	patch_2D->write_mesh();

	std::cout << std::endl;
	std::cout << "  :: 3D Mesh adaption test ::" << std::endl;

	std::cout << std::endl;
	std::cout << ">> Creating mesh" << std::endl;

	pman::Patch *patch_3D = new pman::PatchOctree(0, 3, origin, length, dh);
	patch_3D->set_name("octree_adapted_mesh_3D");
	patch_3D->update();
	patch_3D->write_mesh();

	std::cout << std::endl;
	std::cout << ">> Adapting mesh" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch_3D->get_cell_count();
		std::cout << std::endl;
		std::cout << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch_3D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_3D->extract_cell_neighs(cellId)) {
				for (auto coarseId : patch_3D->extract_cell_neighs(neighId)) {
					patch_3D->mark_cell_for_coarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch_3D->cells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch_3D->extract_cell_neighs(cellId)) {
				patch_3D->mark_cell_for_refinement(neighId);
			}
		}

		std::cout << std::endl;
		std::cout << ">> Initial number of cells... " << nCells << std::endl;

		patch_3D->update();
		
		nCells = patch_3D->get_cell_count();
		std::cout << ">> Final number of cells... " << nCells << std::endl;
	}
	patch_3D->write_mesh();

#ifndef DISABLE_MPI
	MPI::Finalize();
#endif

}

