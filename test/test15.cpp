#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;
		Class_Para_Tree<2> pablo15;
		int idx = 0;
		pablo15.setBalance(idx,false);
		pablo15.computeConnectivity();
		for (iter=1; iter<6; iter++){
			pablo15.adaptGlobalRefine();
			pablo15.updateConnectivity();
		}

		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		uint32_t nocts = pablo15.getNumOctants();
		uint32_t nghosts = pablo15.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);
		for (int i=0; i<nocts; i++){
			dvector2D nodes = pablo15.getNodes(i);
			vector<double> center = pablo15.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
					if (center[0]<=xc){
						pablo15.setMarker(i,1);
					}
					else{
						pablo15.setMarker(i,-1);
					}
				}
			}
		}

		iter = 0;
		pablo15.updateConnectivity();
		pablo15.writeTest("Pablo15_iter"+to_string(iter), oct_data);

		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				dvector2D nodes = pablo15.getNodes(i);
				vector<double> center = pablo15.getCenter(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){
							pablo15.setMarker(i,1);
						}
						else{
							pablo15.setMarker(i,-1);
						}
					}
				}
			}
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			pablo15.adapt(mapper);
			nocts = pablo15.getNumOctants();
			oct_data_new.resize(nocts, 0.0);
			for (int i=0; i<nocts; i++){
				if (pablo15.getIsNewC(i)){
					for (int j=0; j<global2D.nchildren; j++){
						oct_data_new[i] += oct_data[mapper[i]+j]/global2D.nchildren;
					}
				}
				else{
					oct_data_new[i] += oct_data[mapper[i]];
				}
			}

			pablo15.updateConnectivity();
			pablo15.writeTest("Pablo15_iter"+to_string(iter), oct_data_new);
			oct_data = oct_data_new;
		}

		//TODO LOADBLANCE WITH DATA - NOT WORKING!!!!
//		User_Data_LB<vector<double> > data_lb(oct_data);
//		pablo15.loadBalance(data_lb);
		pablo15.updateConnectivity();
		pablo15.writeTest("Pablo15_iter"+to_string(iter), oct_data);


	}

	MPI::Finalize();

}



