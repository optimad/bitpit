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
		ptreedefault.adaptGlobalRefine();
		ptreedefault.computeConnectivity();
		ptreedefault.write("Pablo_default_iter0");

		double radius = 0.4;
		double xc, yc;
		xc = yc = 0.5;

		for (int iter=0; iter<8; iter++){
			uint32_t nocts = ptreedefault.getNumOctants();
			for (int i=0; i<nocts; i++){
				Class_Octant<2> *oct = ptreedefault.getOctant(i);
				vector<double> center = ptreedefault.getCenter(oct);
				double x = center[0];
				double y = center[1];
				if (pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)){
					ptreedefault.setMarker(oct, 1);
				}
			}
			ptreedefault.adapt();
			ptreedefault.updateConnectivity();
			ptreedefault.write("Pablo_default_iter"+to_string(iter+1));
		}
	}

	MPI::Finalize();

}

