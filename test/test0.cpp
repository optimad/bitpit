#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		/**<Instantation of a 2D para_tree object with default constructor.*/
		Class_Para_Tree<2> ptreedefault;
		/**<Write the para_tree in physical domain.*/
		ptreedefault.write("Pablo0_default");
		/**<Write the para_tree in logical domain.*/
		ptreedefault.writeLogical("Pablo0_default_logical");

		double X, Y, Z, L;
		uint8_t level0 = MAX_LEVEL_2D;
		X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
		/**<Instantation of a 2D para_tree object with custom constructor.*/
		Class_Para_Tree<2> ptreecustom(X, Y, Z, L);
		/**<Write the para_tree in physical domain.*/
		ptreecustom.write("Pablo0_custom");
		/**<Write the para_tree in logical domain.*/
		ptreecustom.writeLogical("Pablo0_custom_logical");

#if NOMPI==0
	}

	MPI::Finalize();
#endif
}

