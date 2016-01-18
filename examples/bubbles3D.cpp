#include "ParaTree.hpp"
#if ENABLE_MPI==1
#include "UserDataComm.hpp"
#include "UserDataLB.hpp"
#endif

using namespace std;

// =================================================================================== //
/*!
	\example bubbles3D.cpp

	\brief 3D dynamic adaptive mesh refinement using PABLO 

	This example creates a 2D Cartesian mesh on the square domain [-10,10]x[-10,10].
	The domain is discretized with a cell size of 0.5 both in x and y directions.

	<b>To run</b>: ./createCartesianMesh2D \n
*/

/**<Declaration of a class bubble with center and radius.*/

/*!  \cond  EXAMPLE_CLASSES */
class bubble{
public:
	double c[3];
	double r;


};
/*!  \endcond  */

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 3D para_tree object.*/
		ParaTree pabloBB(3);

		/**<Set 2:1 balance for the octree.*/
		pabloBB.setBalanceCodimension(1);
		int idx = 0;
		pabloBB.setBalance(idx,true);

		/**<Refine globally four level and write the para_tree.*/
        for (iter=1; iter<4; iter++){
			pabloBB.adaptGlobalRefine();
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pabloBB.loadBalance();
#endif

		/**<Update the number of local octants.*/
		uint32_t nocts = pabloBB.getNumOctants();

		/**<Define and initialize a set of bubbles and their trajectories.*/
		srand(time(NULL));

#if ENABLE_MPI==1
        int nb = 50;
#else
		int nb = 10;
#endif
		vector<bubble> BB;
		vector<bubble> BB0;
		vector<double> DZ;
		vector<double> OM;
		vector<double> AA;

		for (int i=0; i<nb; i++){
			double randc[3];
			randc[0] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[1] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[2] = (double) (rand()) /  RAND_MAX - 0.5;
			double randr = (0.05 * (double) (rand()) / RAND_MAX + 0.04);
			double dz = 0.005 + 0.05 * (double) (rand()) / RAND_MAX;
			double omega = 0.5 * (double) (rand()) / RAND_MAX;
			double aa = 0.15 * (double) (rand()) / RAND_MAX;
			bubble bb;
			bb.c[0] = randc[0];
			bb.c[1] = randc[1];
			bb.c[2] = randc[2];
			bb.r = randr;
			BB.push_back(bb);
			BB0.push_back(bb);
			DZ.push_back(dz);
			OM.push_back(omega);
			AA.push_back(aa);
		}
		/**<Initialize time and timestep.*/
		double t0 = 0;
		double t = t0;
		double Dt = 0.5;

		/**<Define vectors of data.*/
		vector<double> oct_data(nocts, 0.0),oct_data_ghost;
		vector<double> oct_data_new(nocts, 0.0);

		/**<Adapt itend times with refinement on the interface of the bubbles.*/
		int itstart = 1;
		int itend = 200;

		/**<Perform time iterations.*/
		for (iter=itstart; iter<itend; iter++){
			if(pabloBB.getRank()==0) cout << "iter " << iter << endl;
			t += Dt;

			/**<Update bubbles position.*/
			for (int i=0; i<nb; i++){
				BB[i].c[0] = BB0[i].c[0];// + AA[i]*cos(OM[i]*t);
				BB[i].c[1] = BB0[i].c[1];// + AA[i]*cos(OM[i]*t);
				BB[i].c[2] = BB[i].c[2]+ Dt*DZ[i];
			}

			/**<Adapting (refinement and coarsening).*/
			bool adapt = true;
			while (adapt){

				for (int i=0; i<nocts; i++){
					bool inside = false;
					/**<Compute the nodes of the octant.*/
					vector<array<double,3> > nodes = pabloBB.getNodes(i);
					/**<Compute the center of the octant.*/
					array<double,3> center = pabloBB.getCenter(i);
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double zc = BB[ib].c[2];
						double radius = BB[ib].r;
						/**<Set marker with condition on center or nodes of the octant.*/
						if (((!inside))){
							if (((pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0) <= 1.25*pow(radius,2.0) &&
											pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0) >= 0.75*pow(radius,2.0)))){
                                if (pabloBB.getLevel(i) < 7){
									/**<Set to refine inside a band around the interface of the bubbles.*/
									pabloBB.setMarker(i,1);
									oct_data[i] = 1.0;
								}
								else{
									pabloBB.setMarker(i,0);
									oct_data[i] = 1.0;
								}
								inside = true;
							}
						}
						for (int j=0; j<8; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							double z = nodes[j][2];
							if (((!inside))){
								if (((pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) <= 1.25*pow(radius,2.0) &&
										pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) >= 0.75*pow(radius,2.0)))){
									if (pabloBB.getLevel(i) < 7){
										/**<Set to refine inside a band around the interface of the bubbles.*/
										pabloBB.setMarker(i,1);
										oct_data[i] = 1.0;
									}
									else{
										pabloBB.setMarker(i,0);
										oct_data[i] = 1.0;
									}
									inside = true;
								}
							}
						}
						ib++;
					}
					if (!inside){
						oct_data[i] = 0.0;
                        if (pabloBB.getLevel(i) > 3){
							/**<Set to coarse outside the band if the octant has a level higher than 4.*/
							pabloBB.setMarker(i,3-pabloBB.getLevel(i));
						}
					}
				}

				/**<Adapt the octree.*/
				vector<uint32_t> mapidx;
				vector<bool> isghost;
				adapt = pabloBB.adapt(true);

				/**<Update the number of local octants.*/
				nocts = pabloBB.getNumOctants();

				/**<Assign to the new octant the data after an adaption.*/
				oct_data_new.resize(nocts, 0);
				for (uint32_t i=0; i<nocts; i++){
					pabloBB.getMapping(i, mapidx, isghost);
					oct_data_new[i] = oct_data[mapidx[0]];
				}
				oct_data = oct_data_new;
				vector<double>().swap(oct_data_new);

			}

#if ENABLE_MPI==1
			/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
			/**<Communicate the data of the octants and the ghost octants between the processes.*/
			UserDataLB<vector<double> > data_lb(oct_data,oct_data_ghost);
			pabloBB.loadBalance(data_lb);
			/**<Update the number of local octants.*/
			nocts = pabloBB.getNumOctants();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pabloBB.updateConnectivity();
			pabloBB.writeTest("PabloBubble3D_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);
		}

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
