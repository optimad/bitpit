#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;
		/**<Instantation of a 2D para_tree object.*/
		Class_Para_Tree<2> pablo7a;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo7a.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo7a.adaptGlobalRefine();
		}

		/**<Instantation and copy pablo7a.*/
		Class_Para_Tree<2> pablo7b = pablo7a;

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo7a.getNumOctants();
		vector<double> oct_data(nocts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			/**<Compute the center of the octant.*/
			vector<double> center = pablo7a.getCenter(i);
			oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
			/**<Compute the nodes of the octant.*/
			vector<vector<double> > nodes = pablo7a.getNodes(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo7a.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo7a.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo7a.updateConnectivity();
		pablo7b.updateConnectivity();
		pablo7a.write("Pablo7a_iter"+to_string(iter));
//		pablo7b.write("Pablo7b_iter"+to_string(iter));

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				/**<Compute the center of the octant.*/
				vector<double> center = pablo7a.getCenter(i);
				/**<Compute the nodes of the octant.*/
				vector<vector<double> > nodes = pablo7a.getNodes(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo7a.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo7a.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			pablo7a.adapt(mapper);
			nocts = pablo7a.getNumOctants();
			oct_data_new.resize(nocts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (int i=0; i<nocts; i++){
				if (pablo7a.getIsNewC(i)){
					for (int j=0; j<global2D.nchildren; j++){
						oct_data_new[i] += oct_data[mapper[i]+j]/global2D.nchildren;
					}
				}
				else if (pablo7a.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[i]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[i]];
				}
			}
			oct_data = oct_data_new;
			oct_data_new.clear();

			/**<Update the connectivity and write the para_tree.*/
			pablo7a.updateConnectivity();
			pablo7a.writeTest("Pablo7a_iter"+to_string(iter), oct_data);

		}


			/**<Define a mapper for two PABLO objects.*/
			vector<pair<pair<uint32_t,uint32_t>, pair<int, int> > > mapPablo;
			nocts = pablo7b.getNumOctants();
			mapPablo = pablo7b.mapPablos(pablo7a);
			vector<double> oct_data_b(nocts);
			/**<Assign data to the new octree using the second octant index in mapPablo*/
			for (int i=0; i<nocts; i++){
				oct_data_b[i] = oct_data[mapPablo[i].first.second];
			}

			/**<Update the connectivity and write the para_tree.*/
			pablo7b.updateConnectivity();
			pablo7b.writeTest("Pablo7b_iter"+to_string(iter), oct_data_b);

#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
