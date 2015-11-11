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
		Class_Para_Tree<2> pablo17a;
		Class_Para_Tree<2> pablo17b;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo17a.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo17a.adaptGlobalRefine();
			pablo17b.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo17a.getNumOctants();
		uint32_t nghosts = pablo17a.getNumGhosts();
		vector<double> oct_data(nocts, 0.0);
		vector<double> ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			/**<Compute the center of the octant.*/
			vector<double> center = pablo17a.getCenter(i);
			oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
			/**<Compute the nodes of the octant.*/
			vector<vector<double> > nodes = pablo17a.getNodes(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo17a.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo17a.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo17a.updateConnectivity();
		pablo17b.updateConnectivity();
		pablo17a.write("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)));
//		pablo7b.write("Pablo7b_iter"+to_string(static_cast<unsigned long long>((iter)));

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				/**<Compute the center of the octant.*/
				vector<double> center = pablo17a.getCenter(i);
				/**<Compute the nodes of the octant.*/
				vector<vector<double> > nodes = pablo17a.getNodes(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo17a.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo17a.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo17a.adapt(true);
			nocts = pablo17a.getNumOctants();
			oct_data_new.resize(nocts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (uint32_t i=0; i<nocts; i++){
				pablo17a.getMapping(i, mapper, isghost);
				if (pablo17a.getIsNewC(i)){
					for (int j=0; j<global2D.nchildren; j++){
						if (isghost[j]){
							oct_data_new[i] += ghost_data[mapper[j]]/global3D.nchildren;
						}
						else{
							oct_data_new[i] += oct_data[mapper[j]]/global3D.nchildren;
						}
					}
				}
				else if (pablo17a.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
				}
			}
			oct_data = oct_data_new;
			oct_data_new.clear();

			/**<Update the connectivity and write the para_tree.*/
			pablo17a.updateConnectivity();
			pablo17a.writeTest("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		}

#if NOMPI==0
		/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.
		 * Preserve the family compact up to 4 levels over the max deep reached in the octree.*/
		User_Data_LB<vector<double> > data_lb(oct_data,ghost_data);
		pablo17a.loadBalance(data_lb);
		pablo17b.loadBalance();

		/**<Update the connectivity and write the para_tree.*/
		pablo17a.updateConnectivity();
		pablo17a.writeTest("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

#endif


		/**<Define a mapper for two PABLO objects.*/
		vector<pair<pair<uint32_t,uint32_t>, pair<int, int> > > mapPablo;
		nocts = pablo17b.getNumOctants();
		mapPablo = pablo17b.mapPablos(pablo17a);
		vector<double> oct_data_b(nocts);
//		/**<Assign data to the new octree : process owner the first octant index in mapPablo.*/
//		for (int i=0; i<nocts; i++){
//			oct_data_b[i] = mapPablo[i].second.first;
//		}
		/**<Assign data to the new octree : first octant index in mapPablo.*/
		for (int i=0; i<nocts; i++){
			oct_data_b[i] = mapPablo[i].first.first;
		}

		/**<Update the connectivity and write the para_tree.*/
		pablo17b.updateConnectivity();
		pablo17b.writeTest("Pablo17b_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_b);

#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
