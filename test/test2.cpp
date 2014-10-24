#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{

		Class_Para_Tree<2> pablo2;
		pablo2.computeConnectivity();
		pablo2.write("Pablo2_iter0");
		for (int iter=1; iter<3; iter++){
			pablo2.adaptGlobalRefine();
			pablo2.updateConnectivity();
			pablo2.write("Pablo2_iter"+to_string(iter));
		}

		pablo2.loadBalance();

		double xc, yc;
		xc = yc = 0.5;

		double radius = 0.4;

		// Simple adapt() in upper area of domain
		for (int iter=3; iter<9; iter++){
			uint32_t nocts = pablo2.getNumOctants();
			for (int i=0; i<nocts; i++){
				Class_Octant<2> *oct = pablo2.getOctant(i);
				dvector2D nodes = pablo2.getNodes(oct);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo2.setMarker(oct, 1);
					}
				}
			}
			pablo2.adapt();
			pablo2.loadBalance();
			pablo2.updateConnectivity();
			pablo2.write("Pablo2_iter"+to_string(iter));
		}
	}

	MPI::Finalize();

}

