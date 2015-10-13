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
	patch_2D->set_name("octree_backgorund_mesh_2D");
	patch_2D->update();
	patch_2D->write_mesh();

	std::cout << "\n  >> 2D Mesh Neighbour test" << "\n";

	std::vector<long> neighs_2D;

	long cellId_2D = 7;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_2D << std::endl << std::endl;

	pman::Cell &cell_2D = patch_2D->get_cell(cellId_2D);

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_2D = cell_2D.extract_face_neighs();
	for (int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_2D = cell_2D.extract_vertex_neighs(true);
	for (int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face neighbours): " << std::endl;
	neighs_2D = cell_2D.extract_vertex_neighs(false);
	for (int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_2D;

	std::cout << "  >> 3D Mesh" << "\n";

	pman::Patch *patch_3D = new pman::PatchOctree(0, 3, origin, length, dh);
	patch_3D->set_name("octree_backgorund_mesh_3D");
	patch_3D->update();
	patch_3D->write_mesh();


	std::cout << "\n  >> 3D Mesh Neighbour test" << "\n";

	std::vector<long> neighs_3D;

	long cellId_3D = 13;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_3D << std::endl << std::endl;

	pman::Cell &cell_3D = patch_3D->get_cell(cellId_3D);

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_3D = cell_3D.extract_face_neighs();
	for (int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (complete list): " << std::endl;
	neighs_3D = cell_3D.extract_edge_neighs(true);
	for (int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (excuding face neighbours): " << std::endl;
	neighs_3D = cell_3D.extract_edge_neighs(false);
	for (int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_3D = cell_3D.extract_vertex_neighs(true);
	for (int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face and edge neighbours): " << std::endl;
	neighs_3D = cell_3D.extract_vertex_neighs(false);
	for (int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_3D;

#ifndef DISABLE_MPI
	MPI::Finalize();
#endif

}

