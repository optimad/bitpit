#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {


	MPI::Init(argc, argv);

		{
			double X, Y, Z, L;
			uint8_t level0 = MAX_LEVEL_3D;
			X = 0.0; Y = 0.0; Z = 0.0; L = 1.0;
			Class_Para_Tree<3> ptree(X,Y,Z,L);
			for (int i=0; i<ptree.getNumOctants(); i++){
				Class_Octant<3>* oct = ptree.getOctant(i);
				ptree.setMarker(oct, 2);
			}
			bool done = ptree.adapt();
			ptree.loadBalance();
			for (int e=0; e<5; e++){
				for (int i=0; i<100; i++){
					int v1 = rand() % ptree.getNumOctants();
					Class_Octant<3>* oct = ptree.getOctant(v1);
					ptree.setMarker(oct, 1);
				}
				done = ptree.adapt();
				ptree.loadBalance();
				ptree.updateGhostsConnectivity();
				ptree.writeLogical("Pablo_3D"+to_string(e));
			}
		}
	MPI::Finalize();


	return 0;

}
