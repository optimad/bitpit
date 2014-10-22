#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{

		Class_Para_Tree<2> ptreedefault;
		ptreedefault.write("Pablo_default");
		ptreedefault.writeLogical("Pablo_default_logical");

		double X, Y, Z, L;
		uint8_t level0 = MAX_LEVEL_2D;
		X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
		Class_Para_Tree<2> ptreecustom(X, Y, Z, L);
		ptreecustom.write("Pablo_custom");
		ptreecustom.writeLogical("Pablo_custom_logical");

	}

	MPI::Finalize();

}

