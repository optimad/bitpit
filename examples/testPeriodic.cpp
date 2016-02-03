#include "ParaTree.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D para_tree object.*/
		ParaTree pabloP;

		/**<Set NO 2:1 balance for the octree.*/
		uint32_t idx=0;
		pabloP.setBalance(idx,true);

		/**<Compute the connectivity and write the para_tree.*/
		pabloP.computeConnectivity();
		pabloP.updateGhostsConnectivity();
		pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Refine globally two level and write the para_tree.*/
		for (iter=1; iter<3; iter++){
			pabloP.adaptGlobalRefine();
			pabloP.updateConnectivity();
			pabloP.updateGhostsConnectivity();
			pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pabloP.loadBalance();
#endif

		pabloP.updateConnectivity();
		pabloP.updateGhostsConnectivity();
		pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.2;
		double radius = 0.2;

		/**<Simple adapt() (refine) 6 times the octants with at least one node inside the circle.*/
		for (iter=4; iter<10; iter++){
			uint32_t nocts = pabloP.getNumOctants();
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pabloP.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pabloP.setMarker(i, 1);
					}
				}
			}

			/**<Adapt octree.*/
			pabloP.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pabloP.loadBalance();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pabloP.updateConnectivity();
			pabloP.updateGhostsConnectivity();
			pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Coarse globally one level and write the para_tree.*/
		pabloP.adaptGlobalCoarse();
		pabloP.updateConnectivity();
//		pabloP.updateGhostsConnectivity();
		pabloP.clearGhostsConnectivity();
		pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));

		exit(0);

		/**<Define a center point and a radius.*/
		xc = yc = 0.1;
		radius = 0.1;

		/**<Simple adapt() 5 times the octants with at least one node inside the circle.*/
		for (iter=10; iter<15; iter++){
			uint32_t nocts = pabloP.getNumOctants();
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pabloP.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					/**<Set refinement marker=-1 (coarse it one time) for octants inside a circle.*/
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pabloP.setMarker(i, -1);
					}
				}
			}
			/**<Adapt octree.*/
			pabloP.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pabloP.loadBalance();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pabloP.updateConnectivity();
			pabloP.write("PabloPeriodic_iter"+to_string(static_cast<unsigned long long>(iter)));
		}
#if ENABLE_MPI==1
	}
	MPI::Finalize();
#endif
}
