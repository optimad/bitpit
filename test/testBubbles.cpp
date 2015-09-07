#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

/**<Declaration of a class bubble with center and radius.*/
class bubble{
public:
	double c[2];
	double r;
};


int main(int argc, char *argv[]) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D para_tree object.*/
		Class_Para_Tree<2> pabloBB;

		/**<Set 2:1 balance for the octree.*/
		pabloBB.setBalanceCodimension(1);
		int idx = 0;
		pabloBB.setBalance(idx,true);

		/**<Refine globally four level and write the para_tree.*/
		for (iter=1; iter<7; iter++){
			pabloBB.adaptGlobalRefine();
		}

#if NOMPI==0
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pabloBB.loadBalance();
#endif

		/**<Update the number of local octants.*/
		uint32_t nocts = pabloBB.getNumOctants();

		/**<Define and initialize a set of bubbles and their trajectories.*/
		time_t Time = time(NULL)
		srand(Time);
		cout << "the seed = " << Time << endl;

#if NOMPI==0
		int nb = 100;
#else
		int nb = 10;
#endif
		vector<bubble> BB;
		vector<bubble> BB0;
		vector<double> DY;
		vector<double> OM;
		vector<double> AA;

		for (int i=0; i<nb; i++){
			double randc[2];
			randc[0] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[1] = (double) (rand()) /  RAND_MAX - 0.5;
			double randr = 0.05 * (double) (rand()) / RAND_MAX + 0.02;
			double dy = 0.005 + 0.05 * (double) (rand()) / RAND_MAX;
			double omega = 0.5 * (double) (rand()) / RAND_MAX;
			double aa = 0.15 * (double) (rand()) / RAND_MAX;
			bubble bb;
			bb.c[0] = randc[0];
			bb.c[1] = randc[1];
			bb.r = randr;
			BB.push_back(bb);
			BB0.push_back(bb);
			DY.push_back(dy);
			OM.push_back(omega);
			AA.push_back(aa);
		}
		/**<Initialize time and timestep.*/
		double t0 = 0;
		double t = t0;
		double Dt = 0.5;

		/**<Adapt itend times with refinement on the interface of the bubbles.*/
		int itstart = 1;
		int itend = 200;

		/**<Set the number of refinement per time iteration.*/
		int nrefperiter = 4;

		/**<Perform time iterations.*/
		for (iter=itstart; iter<itend; iter++){
			if(pabloBB.rank==0) cout << "iter " << iter << endl;
			t += Dt;

			/**<Update bubbles position.*/
			for (int i=0; i<nb; i++){
				BB[i].c[0] = BB0[i].c[0] + AA[i]*cos(OM[i]*t);
				BB[i].c[1] = BB[i].c[1]+ Dt*DY[i];
			}

			/**<Adapting (refinement and coarsening).*/
			for (int iref=0; iref<nrefperiter; iref++){

				for (int i=0; i<nocts; i++){
					bool inside = false;
					/**<Compute the nodes of the octant.*/
					vector<vector<double> > nodes = pabloBB.getNodes(i);
					/**<Compute the center of the octant.*/
					vector<double> center = pabloBB.getCenter(i);
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double radius = BB[ib].r;
						/**<Set marker with condition on center or nodes of the octant.*/
						for (int j=0; j<global2D.nnodes; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							if ( ((!inside) &&
									(pow((x-xc),2.0)+pow((y-yc),2.0) <= 1.15*pow(radius,2.0) &&
											pow((x-xc),2.0)+pow((y-yc),2.0) >= 0.85*pow(radius,2.0)))
											|| ((!inside) && (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) <= 1.15*pow(radius,2.0) &&
													pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) >= 0.85*pow(radius,2.0)))){
								if (pabloBB.getLevel(i) < 9){
									/**<Set to refine inside a band around the interface of the bubbles.*/
									pabloBB.setMarker(i,1);
								}
								else{
									pabloBB.setMarker(i,0);
								}
								inside = true;
							}
						}
						ib++;
					}
					if (pabloBB.getLevel(i) > 6 && !inside){
						/**<Set to coarse outside the band if the octant has a level higher than 6.*/
						pabloBB.setMarker(i,-1);
					}
				}

				/**<Adapt the octree.*/
				bool adapt = pabloBB.adapt();

#if NOMPI==0
				/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
				pabloBB.loadBalance();
#endif

				/**<Update the number of local octants.*/
				nocts = pabloBB.getNumOctants();

			}

			/**<Update the connectivity and write the para_tree.*/
			pabloBB.updateConnectivity();
			pabloBB.write("PabloBubble_iter"+to_string(static_cast<unsigned long long>(iter)));
		}
#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
