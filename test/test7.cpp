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
		Class_Para_Tree<2> pablo7;
		int idx = 0;

		/**<Set NO 2:1 balance for the octree.*/
		pablo7.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo7.adaptGlobalRefine();
		}

		/**<Instantation and copy pablo6.*/
		Class_Para_Tree<2> pablo70 = pablo7;

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo7.getNumOctants();

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			vector<vector<double> > nodes = pablo7.getNodes(i);
			vector<double> center = pablo7.getCenter(i);
			for (int j=0; j<global2D.nnodes; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo7.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo7.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo7.updateConnectivity();
		pablo70.updateConnectivity();
		pablo7.write("Pablo7_iter"+to_string(iter));
		pablo70.write("Pablo70_iter"+to_string(iter));

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				vector<vector<double> > nodes = pablo7.getNodes(i);
				vector<double> center = pablo7.getCenter(i);
				for (int j=0; j<global2D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo7.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo7.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			pablo7.adapt();
			nocts = pablo7.getNumOctants();

			/**<Update the connectivity and write the para_tree.*/
			pablo7.updateConnectivity();
			pablo7.write("Pablo7_iter"+to_string(iter));

			vector<pair<pair<uint32_t,uint32_t>, pair<int, int> > > mapper;
			nocts = pablo70.getNumOctants();
			mapper = pablo70.mapPablos(pablo7);
			vector<double> data(nocts);
			for (int i=0; i<nocts; i++){
				data[i] = (double) mapper[i].first.second;
			}
			pablo70.updateConnectivity();
			pablo70.writeTest("Pablo70_iter"+to_string(iter), data);
		}
#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
