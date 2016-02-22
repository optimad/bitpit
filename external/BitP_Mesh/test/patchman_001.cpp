//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <array>

#include "BitP_Mesh_PATCHMAN.hpp"

int main(int argc, char *argv[]) {

#ifndef DISABLE_MPI
	MPI::Init(argc,argv);
#endif

	std::cout << "Testing Cartesian mesh" << "\n";

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	std::cout << "  >> 2D Mesh" << "\n";

	pman::Patch *patch_2D = new pman::PatchCartesian(0, 2, origin, length, dh);
	patch_2D->set_name("cartesian_backgorund_mesh_2D");
	patch_2D->update();
	patch_2D->write_mesh();

	std::cout << "\n  >> 2D Mesh Neighbour test" << "\n";

	std::vector<long> neighs_2D;

	long cellId_2D = 7;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_2D << std::endl << std::endl;

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->extract_cell_face_neighs(cellId_2D);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->extract_cell_vertex_neighs(cellId_2D, true);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face neighbours): " << std::endl;
	neighs_2D = patch_2D->extract_cell_vertex_neighs(cellId_2D, false);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_2D;

	std::cout << "  >> 3D Mesh" << "\n";

	pman::Patch *patch_3D = new pman::PatchCartesian(0, 3, origin, length, dh);
	patch_3D->set_name("cartesian_backgorund_mesh_3D");
	patch_3D->update();
	patch_3D->write_mesh();

	std::cout << "\n  >> 3D Mesh Neighbour test" << "\n";

	std::vector<long> neighs_3D;

	long cellId_3D = 13;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_3D << std::endl << std::endl;

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extract_cell_face_neighs(cellId_3D);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extract_cell_edge_neighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (excuding face neighbours): " << std::endl;
	neighs_3D = patch_3D->extract_cell_edge_neighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extract_cell_vertex_neighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face and edge neighbours): " << std::endl;
	neighs_3D = patch_3D->extract_cell_vertex_neighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_3D;

#ifndef DISABLE_MPI
	MPI::Finalize();
#endif

}
