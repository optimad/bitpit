#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

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

		/**<Define a set of bubbles.*/
		srand(time(NULL));
//		srand(1418143772);


		//8 proc deadlock 1418143772  at iter 51 proc

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
		double t0 = 0;
		double t = t0;
		double Dt = 0.5;
//		double Dy = 0.01;
//		double omega = 2.0*3.14/0.001;

		/**<Define vectors of data.*/
		uint32_t nocts = pabloBB.getNumOctants();
		uint32_t nghosts = pabloBB.getNumGhosts();
//		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Adapt itend times with data injection on new octants.*/
		int itstart = 1;
		int itend = 200;


		int nrefperiter = 4;

		for (iter=itstart; iter<itend; iter++){
			if(pabloBB.rank==0) cout << "iter " << iter << endl;
			t += Dt;

			for (int i=0; i<nb; i++){
				BB[i].c[0] = BB0[i].c[0] + AA[i]*cos(OM[i]*t);
				BB[i].c[1] = BB[i].c[1]+ Dt*DY[i];
			}

			for (int iref=0; iref<nrefperiter; iref++){

				for (int i=0; i<nocts; i++){
					bool inside = false;
					vector<vector<double> > nodes = pabloBB.getNodes(i);
					vector<double> center = pabloBB.getCenter(i);
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double radius = BB[ib].r;
						//oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
//						oct_data[i] = 0.0;
						for (int j=0; j<global2D.nnodes; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= 1.15*pow(radius,2.0) &&
									pow((x-xc),2.0)+pow((y-yc),2.0) >= 0.85*pow(radius,2.0)) ||
									(pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) <= 1.15*pow(radius,2.0) &&
											pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) >= 0.85*pow(radius,2.0))){
								if (pabloBB.getLevel(i) < 9){
									/**<Set to refine inside the sphere.*/
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
						/**<Set to coarse if the octant has a level higher than 5.*/
						pabloBB.setMarker(i,-1);
					}
				}

				/**<Adapt the octree.*/
				//vector<uint32_t> mapidx;
				//bool adapt = pabloBB.adapt(mapidx);
				bool adapt = pabloBB.adapt();

#if NOMPI==0
				/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
				pabloBB.loadBalance();
#endif

				nocts = pabloBB.getNumOctants();
				nghosts = pabloBB.getNumGhosts();
//				vector<double> oct_data_new(nocts, 0.0);

//				/**<Assign to the new octant the data after an adaption.*/
//				for (int i=0; i<nocts; i++){
//					vector<double> center = pabloBB.getCenter(i);
//					//oct_data_new[i] = oct_data[mapidx[i]];
//					oct_data_new[i] = 0.0;
//				}

//				oct_data.resize(nocts);
//				oct_data = oct_data_new;
//				oct_data_new.clear();

			}

			/**<Update the connectivity and write the para_tree.*/
			pabloBB.updateConnectivity();
			//pabloBB.writeTest("PabloBubble_iter"+to_string(iter), oct_data);
			pabloBB.write("PabloBubble_iter"+to_string(iter));
		}
#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
