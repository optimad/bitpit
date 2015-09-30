#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D para_tree object.*/
		Class_Para_Tree<2> pablo;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo.adaptGlobalRefine();
		}

#if NOMPI==0
		/**<(Load)Balance the octree over the processes with communicating the data.*/
		//		User_Data_LB<vector<double> > data_lb(oct_data);
		pablo.loadBalance();
#endif


		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo.getNumOctants();
		uint32_t nghosts = pablo.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			/**<Compute the nodes of the octant.*/
			vector<vector<double> > nodes = pablo.getNodes(i);
			/**<Compute the center of the octant.*/
			vector<double> center = pablo.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
				}
			}
		}

		{
#if NOMPI==0
			/**<Communicate data and data ghost over the processes.*/
			User_Data_Comm<vector<double> > data_comm(oct_data, ghost_data);
			pablo.communicate(data_comm);
#endif
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo.updateConnectivity();
		pablo.writeTest("Pablo_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<vector<double> > nodes = pablo.getNodes(i);
				/**<Compute the center of the octant.*/
				vector<double> center = pablo.getCenter(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new, ghost_data_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo.adapt(true);
			nocts = pablo.getNumOctants();
			nghosts = pablo.getNumGhosts();
			oct_data_new.resize(nocts, 0.0);
			ghost_data_new.resize(nghosts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			cout << "iter " << iter << endl;
			for (uint32_t i=0; i<nocts; i++){
				pablo.getMapping(i, mapper, isghost);
				if (pablo.getIsNewC(i)){
					for (int j=0; j<global2D.nchildren; j++){
						if (isghost[j]){
							cout << pablo.rank << " using ghost for idx " << i << endl;
							oct_data_new[i] += ghost_data[mapper[j]]/global2D.nchildren;
						}
						else{
							oct_data_new[i] += oct_data[mapper[j]]/global2D.nchildren;
						}
					}
				}
				else if (pablo.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
				}
			}

			oct_data = oct_data_new;
			ghost_data = ghost_data_new;

#if NOMPI==0
			/**<(Load)Balance the octree over the processes with communicating the data.*/
			User_Data_LB<vector<double> > data_lb(oct_data,ghost_data);
			pablo.loadBalance(data_lb);
#endif

			nocts = pablo.getNumOctants();
			nghosts = pablo.getNumGhosts();
			ghost_data_new.resize(nghosts, 0.0);
			ghost_data = ghost_data_new;


#if NOMPI==0
			/**<Communicate data and data ghost over the processes.*/
			User_Data_Comm<vector<double> > data_comm(oct_data, ghost_data);
			pablo.communicate(data_comm);
#endif

			cout << " Out comm " << endl;

			/**<Update the connectivity and write the para_tree.*/
			pablo.updateConnectivity();
			pablo.writeTest("Pablo_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		}


		/**<Update the connectivity and write the para_tree.*/
		pablo.updateConnectivity();
		pablo.writeTest("Pablo_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);






#if NOMPI==0
	}

	MPI::Finalize();
#endif
}



