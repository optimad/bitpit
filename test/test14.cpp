#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"
#include "User_Data_Comm.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;
		Class_Para_Tree<2> pablo14;
		pablo14.computeConnectivity();
		for (iter=1; iter<5; iter++){
			pablo14.adaptGlobalRefine();
			pablo14.updateConnectivity();
		}

		pablo14.loadBalance();

		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		uint32_t nocts = pablo14.getNumOctants();
		uint32_t nghosts = pablo14.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);
		for (int i=0; i<nocts; i++){
			dvector2D nodes = pablo14.getNodes(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = 1.0;
				}
			}
		}

		for (int i=0; i<nghosts; i++){
			Class_Octant<2> *oct = pablo14.getGhostOctant(i);
			dvector2D nodes = pablo14.getNodes(oct);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					ghost_data[i] = 1.0;
				}
			}
		}

		iter = 0;
		pablo14.updateConnectivity();
		pablo14.writeTest("Pablo14_iter"+to_string(iter), oct_data);

		int start = iter + 1;
		for (iter=start; iter<start+25; iter++){
		vector<double> oct_data_smooth(nocts, 0.0);
		vector<uint32_t> neigh, neigh_t;
		vector<bool> isghost, isghost_t;
		uint8_t iface, codim;
		for (int i=0; i<nocts; i++){
			neigh.clear();
			isghost.clear();
			codim = 1;
			for (iface=0; iface<global2D.nfaces; iface++){
				pablo14.findNeighbours(i,iface,codim,neigh_t,isghost_t);
				neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
				isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
			}
			codim = 2;
			for (iface=0; iface<global2D.nnodes; iface++){
				pablo14.findNeighbours(i,iface,codim,neigh_t,isghost_t);
				neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
				isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
			}
			oct_data_smooth[i] = oct_data[i]/(neigh.size()+1);
			for (int j=0; j<neigh.size(); j++){
				if (isghost[j]){
					oct_data_smooth[i] += ghost_data[neigh[j]]/(neigh.size()+1);
				}
				else{
					oct_data_smooth[i] += oct_data[neigh[j]]/(neigh.size()+1);
				}
			}
		}

		pablo14.updateConnectivity();
		pablo14.writeTest("Pablo14_iter"+to_string(iter), oct_data_smooth);

		User_Data_Comm<vector<double> > data_comm(oct_data_smooth, ghost_data);
		pablo14.communicate(data_comm);
		oct_data = oct_data_smooth;

	}
	}

	MPI::Finalize();

}


