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
		/**<Instantation of a 2D para_tree object.*/
		Class_Para_Tree<2> pablo5;
		int idx = 0;

		/**<Set NO 2:1 balance for the octree.*/
		pablo5.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo5.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo5.getNumOctants();
		uint32_t nghosts = pablo5.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			vector<vector<double> > nodes = pablo5.getNodes(i);
			vector<double> center = pablo5.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo5.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo5.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo5.updateConnectivity();
		pablo5.writeTest("Pablo5_iter"+to_string(iter), oct_data);

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				vector<vector<double> > nodes = pablo5.getNodes(i);
				vector<double> center = pablo5.getCenter(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo5.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo5.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			pablo5.adapt(mapper);
			nocts = pablo5.getNumOctants();
			oct_data_new.resize(nocts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
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

			/**<Update the connectivity and write the para_tree.*/
			pablo5.updateConnectivity();
			pablo5.writeTest("Pablo5_iter"+to_string(iter), oct_data_new);

			oct_data = oct_data_new;
		}

	}

	MPI::Finalize();

}

