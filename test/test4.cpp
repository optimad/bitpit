#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;
		int dim = 2;
		/**<Instantation of a 2D para_tree object.*/
		Class_Para_Tree<2> pablo4;

		/**<Refine globally four level and write the para_tree.*/
		for (iter=1; iter<5; iter++){
			pablo4.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo4.getNumOctants();
		uint32_t nghosts = pablo4.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			vector<vector<double> > nodes = pablo4.getNodes(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = 1.0;
				}
			}
		}

		/**<Assign a data to the ghost octants (NONE IT IS A SERIAL TEST) with at least one node inside the circle.*/
		for (int i=0; i<nghosts; i++){
			Class_Octant<2> *oct = pablo4.getGhostOctant(i);
			vector<vector<double> > nodes = pablo4.getNodes(oct);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					ghost_data[i] = 1.0;
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo4.updateConnectivity();
		pablo4.writeTest("Pablo4_iter"+to_string(iter), oct_data);

		/**<Smoothing iterations on initial data*/
		int start = 1;
		for (iter=start; iter<start+25; iter++){
			vector<double> oct_data_smooth(nocts, 0.0);
			vector<uint32_t> neigh, neigh_t;
			vector<bool> isghost, isghost_t;
			uint8_t iface, nfaces;
			int codim;
			for (int i=0; i<nocts; i++){
				neigh.clear();
				isghost.clear();

				/**<Find neighbours through edges (codim=1) and nodes (codim=2) of the octants*/
				for (codim=1; codim<dim+1; codim++){
					if (codim == 1){
						nfaces = global2D.nfaces;
					}
					else if (codim == 2){
						nfaces = global2D.nnodes;
					}
					for (iface=0; iface<nfaces; iface++){
						pablo4.findNeighbours(i,iface,codim,neigh_t,isghost_t);
						neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
						isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
					}
				}

				/**<Smoothing data with the average over the one ring neighbours of octants*/
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

			/**<Update the connectivity and write the para_tree.*/
			pablo4.updateConnectivity();
			pablo4.writeTest("Pablo4_iter"+to_string(iter), oct_data_smooth);

			oct_data = oct_data_smooth;
		}
	}

	MPI::Finalize();

}
