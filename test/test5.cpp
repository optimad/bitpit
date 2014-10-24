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
		Class_Para_Tree<2> pablo5;
		int idx = 0;
		pablo5.setBalance(idx,false);
		pablo5.computeConnectivity();
		for (iter=1; iter<5; iter++){
			pablo5.adaptGlobalRefine();
			pablo5.updateConnectivity();
		}

		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		uint32_t nocts = pablo5.getNumOctants();
		uint32_t nghosts = pablo5.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);
		for (int i=0; i<nocts; i++){
			dvector2D nodes = pablo5.getNodes(i);
			vector<double> center = pablo5.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
					if (center[0]<=xc){
						pablo5.setMarker(i,1);
					}
					else{
						pablo5.setMarker(i,-1);
					}
				}
			}
		}

		iter = 0;
		pablo5.updateConnectivity();
		pablo5.writeTest("Pablo5_iter"+to_string(iter), oct_data);

		int start = 1;
		for (iter=start; iter<start+1; iter++){
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			pablo5.adapt(mapper);
			nocts = pablo5.getNumOctants();
			oct_data_new.resize(nocts, 0.0);
			for (int i=0; i<nocts; i++){
				if (pablo5.getIsNewC(i)){
					for (int j=0; j<global2D.nchildren; j++){
						oct_data_new[i] += oct_data[mapper[i]+j]/global2D.nchildren;
					}
				}
				else{
					oct_data_new[i] += oct_data[mapper[i]];
				}
			}

			pablo5.updateConnectivity();
			pablo5.writeTest("Pablo5_iter"+to_string(iter), oct_data_new);
			oct_data = oct_data_new;
		}

	}

	MPI::Finalize();

}

