#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;

		/**<Instantation of a 3D para_tree object.*/
		Class_Para_Tree<3> pablo120;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo120.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo120.adaptGlobalRefine();
		}

//		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
//		pablo120.loadBalance();

		/**<Define a center point and a radius.*/
		double xc, yc, zc;
		double radius = 0.15;
		double t0 = 0;
		double Dt = 0.000025;
		double omega = 2.0*3.14/0.001;

		xc = 0.25*cos(omega* Dt) + 0.5 ;
		yc = 0.25*sin(omega* Dt) + 0.5;
		zc = 100*Dt;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo120.getNumOctants();
		uint32_t nghosts = pablo120.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Adapt itend times with data injection on new octants.*/
		int itstart = 200;
		int itend = 460;

		for (iter=itstart; iter<itend; iter++){

			xc = 0.25*cos(omega* Dt* (double) iter) + 0.5 ;
			yc = 0.25*sin(omega* Dt* (double) iter) + 0.5;
			zc = 100*Dt*iter;

			for (int i=0; i<nocts; i++){
				bool inside = false;
				vector<vector<double> > nodes = pablo120.getNodes(i);
				vector<double> center = pablo120.getCenter(i);
				oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0));
				for (int j=0; j<global3D.nnodes; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					double z = nodes[j][2];
					if (pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) <= pow(radius,2.0)){
						if (pablo120.getLevel(i) < 6){
							/**<Set to refine inside the sphere.*/
							pablo120.setMarker(i,1);
						}
						else{
							pablo120.setMarker(i,0);
						}
						inside = true;
					}
				}
				if (pablo120.getLevel(i) > 5 && !inside){
					/**<Set to coarse if the octant has a level higher than 5.*/
					pablo120.setMarker(i,-1);
				}
			}

			/**<Adapt the octree.*/
			bool adapt = pablo120.adapt();

//			/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
//			pablo120.loadBalance();

			nocts = pablo120.getNumOctants();
			nghosts = pablo120.getNumGhosts();
			vector<double> oct_data_new(nocts, 0.0);

			/**<Assign to the new octant the data after an adaption.*/
			for (int i=0; i<nocts; i++){
					vector<double> center = pablo120.getCenter(i);
					oct_data_new[i] = sqrt((pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0)));
			}

			oct_data.resize(nocts);
			oct_data = oct_data_new;

			/**<Update the connectivity and write the para_tree.*/
			pablo120.updateConnectivity();
			pablo120.writeTest("Pablo120_iter"+to_string(iter), oct_data);
		}
	}

	MPI::Finalize();

}
