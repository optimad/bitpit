#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;
		Class_Para_Tree<2> pablo3;
		pablo3.computeConnectivity();
		pablo3.write("Pablo3_iter"+to_string(iter));
		for (iter=1; iter<3; iter++){
			pablo3.adaptGlobalRefine();
			pablo3.updateConnectivity();
			pablo3.write("Pablo3_iter"+to_string(iter));
		}

		pablo3.loadBalance();

		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.4;

		// Simple adapt()
		for (iter=3; iter<9; iter++){
			uint32_t nocts = pablo3.getNumOctants();
			for (int i=0; i<nocts; i++){
				pablo3.setBalance(i,false);
				dvector2D nodes = pablo3.getNodes(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo3.setMarker(i, 1);
					}
				}
			}
			pablo3.adapt();
			pablo3.loadBalance();
			pablo3.updateConnectivity();
			pablo3.write("Pablo3_iter"+to_string(iter));
		}

		pablo3.adaptGlobalCoarse();
		pablo3.updateConnectivity();
		pablo3.write("Pablo3_iter"+to_string(iter));


		xc = yc = 0.35;
		radius = 0.15;

		// Simple adapt()
		for (iter=10; iter<15; iter++){
			uint32_t nocts = pablo3.getNumOctants();
			for (int i=0; i<nocts; i++){
				pablo3.setBalance(i,false);
				dvector2D nodes = pablo3.getNodes(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo3.setMarker(i, -1);
					}
				}
			}
			pablo3.adapt();
			pablo3.loadBalance();
			pablo3.updateConnectivity();
			pablo3.write("Pablo3_iter"+to_string(iter));
		}
	}
	MPI::Finalize();

}
