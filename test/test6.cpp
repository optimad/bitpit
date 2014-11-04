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
		Class_Para_Tree<2> pablo6;
		int idx = 0;

		/**<Set NO 2:1 balance for the octree.*/
		pablo6.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo6.adaptGlobalRefine();
		}

		/**<Instantation and copy pablo6.*/
		Class_Para_Tree<2> pablo60 = pablo6;
		pablo60.write("Pablo60_iter"+to_string(iter));

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo6.getNumOctants();
		uint32_t nghosts = pablo6.getNumGhosts();

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			vector<vector<double> > nodes = pablo6.getNodes(i);
			vector<double> center = pablo6.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo6.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo6.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo6.updateConnectivity();
		pablo6.write("Pablo6_iter"+to_string(iter));

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				vector<vector<double> > nodes = pablo6.getNodes(i);
				vector<double> center = pablo6.getCenter(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo6.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo6.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			pablo6.adapt();
			nocts = pablo6.getNumOctants();

			/**<Update the connectivity and write the para_tree.*/
			pablo6.updateConnectivity();
			pablo6.write("Pablo6_iter"+to_string(iter));

			vector<pair<pair<uint32_t,uint32_t>, pair<int, int> > > mapper;
			nocts = pablo60.getNumOctants();
			mapper = pablo60.mapPablos(pablo6);
			vector<double> data(nocts);
			for (int i=0; i<nocts; i++){
				data[i] = (double) mapper[i].first.second;
				cout << data[i] << endl;
			}
			pablo60.writeTest("Pablo60_iter"+to_string(iter+1), data);
		}
	}

	MPI::Finalize();

}
